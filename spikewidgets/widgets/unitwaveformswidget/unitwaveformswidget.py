import numpy as np
from matplotlib import pyplot as plt
from spikewidgets.widgets.basewidget import BaseMultiWidget


def plot_unit_waveforms(sorting=None, recording=None, channels=None, unit_ids=None, ms_before=1., ms_after=2.,
                        max_num_waveforms=100, title='', figure=None, ax=None):
    W = UnitWaveformsWidget(
        recording=recording,
        sorting=sorting,
        channels=channels,
        unit_ids=unit_ids,
        max_num_waveforms=max_num_waveforms,
        ms_before=ms_before,
        ms_after=ms_after,
        title=title,
        figure=figure,
        ax=ax
    )
    W.plot()
    return W


class UnitWaveformsWidget(BaseMultiWidget):
    def __init__(self, *, recording, sorting, channels=None, unit_ids=None, max_num_waveforms=50,
                 ms_before=1., ms_after=2., title='', figure=None, ax=None):
        BaseMultiWidget.__init__(self, figure, ax)
        self._recording = recording
        self._sorting = sorting
        self._channels = channels
        self._unit_ids = unit_ids
        self._ms_before = ms_before
        self._ms_after = ms_after
        self._title = title
        self._max_num_waveforms = max_num_waveforms

    def plot(self):
        self._do_plot()

    def _do_plot(self):
        units = self._unit_ids
        channels = self._channels
        if units is None:
            units = self._sorting.get_unit_ids()
        channel_ids = self._recording.get_channel_ids()
        M = len(channel_ids)
        channel_locations = np.zeros((M, 2))
        for ii, ch in enumerate(channel_ids):
            loc = self._recording.get_channel_property(ch, 'location')
            channel_locations[ii, :] = loc[-2:]
        if channels is None:
            channels = channel_ids
        list_spikes = []

        for unit in units:
            spiketrain = self._sorting.get_unit_spike_train(unit_id=unit)
            if spiketrain is not None:
                spikes = self._get_random_spike_waveforms(unit=unit, max_num=self._max_num_waveforms,
                                                          channels=channels)
                item = dict(
                    representative_waveforms=spikes,
                    title='Unit {}'.format(int(unit))
                )
                list_spikes.append(item)
            else:
                print(unit, ' spike train is None')
        with plt.rc_context({'axes.edgecolor': 'gray'}):
            # self._plot_spike_shapes_multi(list,channel_locations=channel_locations[np.array(channels),:])
            self._plot_spike_shapes_multi(list_spikes, channel_locations=None)

    def _get_random_spike_waveforms(self, *, unit, max_num, channels):
        st = self._sorting.get_unit_spike_train(unit_id=unit)
        num_events = len(st)
        if num_events > max_num:
            event_indices = np.random.choice(range(num_events), size=max_num, replace=False)
        else:
            event_indices = range(num_events)

        snippet_len = [int(self._ms_before * self._recording.get_sampling_frequency() / 1000),
                       int(self._ms_after * self._recording.get_sampling_frequency() / 1000)]

        spikes = self._recording.get_snippets(reference_frames=st[event_indices].astype(int),
                                              snippet_len=snippet_len, channel_ids=channels)
        if spikes.size != 0:
            spikes = np.dstack(tuple(spikes))
        else:
            spikes = np.zeros((self._recording.get_num_channels(), np.sum(snippet_len), 0))
        return spikes

    def _plot_spike_shapes_multi(self, list_spikes, *, ncols=5, **kwargs):
        if 'ylim' in kwargs:
            ylim = kwargs['ylim']
        else:
            ylim = _determine_global_ylim(list_spikes)
        if len(list_spikes) < ncols:
            ncols = len(list_spikes)
        nrows = np.ceil(len(list_spikes) / ncols)
        for i, item in enumerate(list_spikes):
            ax = self.get_tiled_ax(i, nrows, ncols)
            _plot_spike_shapes(**item, **kwargs, ylim=ylim, ax=ax)


def _plot_spike_shapes(*, ax, representative_waveforms=None, average_waveform=None, channel_locations=None,
                       ylim=None, max_representatives=None, color='blue', title=''):
    if average_waveform is None:
        if representative_waveforms is None:
            raise Exception('You must provide either average_waveform, representative waveforms, or both')
        average_waveform = np.mean(representative_waveforms, axis=2)
    M = average_waveform.shape[0]  # number of channels
    T = average_waveform.shape[1]  # number of timepoints
    if ylim is None:
        ylim = [average_waveform.min(), average_waveform.max()]
    yrange = ylim[1] - ylim[0]
    if channel_locations is None:
        channel_locations = np.zeros((M, 2))
        for m in range(M):
            channel_locations[m, :] = [0, -m]

    spacing = 1 / 0.8  # TODO: auto-determine this from the channel_locations

    xvals = np.linspace(-yrange / 2, yrange / 2, T)
    if representative_waveforms is not None:
        if max_representatives is not None:
            W0 = representative_waveforms
            if W0.shape[2] > max_representatives:
                indices = np.random.choice(range(W0.shape[2]), size=max_representatives, replace=False)
                representative_waveforms = W0[:, :, indices]
        L = representative_waveforms.shape[2]
        XX = np.zeros((T, M, L))
        YY = np.zeros((T, M, L))
        for m in range(M):
            loc = channel_locations[m, -2:] * yrange * spacing
            for j in range(L):
                XX[:, m, j] = loc[0] + xvals
                YY[:, m, j] = loc[1] + representative_waveforms[m, :, j] - representative_waveforms[m, 0, j]
        XX = XX.reshape(T, M * L)
        YY = YY.reshape(T, M * L)
        ax.plot(XX, YY, color=(0.5, 0.5, 0.5), alpha=0.4)

        XX = np.zeros((T, M))
        YY = np.zeros((T, M))
        for m in range(M):
            loc = channel_locations[m, -2:] * yrange * spacing
            XX[:, m] = loc[0] + xvals
            YY[:, m] = loc[1] + average_waveform[m, :] - average_waveform[m, 0]
        ax.plot(XX, YY, color)

    ax.set_xticks([])
    ax.set_yticks([])
    if title:
        ax.set_title(title, color='gray')


def _get_ylim_for_item(average_waveform=None, representative_waveforms=None):
    if average_waveform is None:
        if representative_waveforms is None:
            raise Exception('You must provide either average_waveform, representative waveforms, or both')
        average_waveform = np.mean(representative_waveforms, axis=2)
    return [average_waveform.min(), average_waveform.max()]


def _determine_global_ylim(list_spikes):
    ret = [0, 0]
    for item in list_spikes:
        ylim0 = _get_ylim_for_item(
            average_waveform=item.get('average_waveform', None),
            representative_waveforms=item.get('representative_waveforms', None)
        )
        ret[0] = np.minimum(ylim0[0], ret[0])
        ret[1] = np.maximum(ylim0[1], ret[1])
    return ret
