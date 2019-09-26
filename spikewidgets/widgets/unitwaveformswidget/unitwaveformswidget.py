import numpy as np
from matplotlib import pyplot as plt
import spiketoolkit as st
from spikewidgets.widgets.basewidget import BaseMultiWidget


def plot_unit_waveforms(recording, sorting, channel_ids=None, unit_ids=None, ms_before=1., ms_after=2.,
                        max_spikes_per_unit=100, channel_locs=False, figure=None, ax=None):
    """
    Plots unit waveforms.

    Parameters
    ----------
    recording: RecordingExtractor
        The recordng extractor object
    sorting: SortingExtractor
        The sorting extractor object
    channel_ids: list
        The channel ids to display.
    unit_ids: list
        List of unit ids.
    ms_before: float
        Time before peak (ms)
    ms_after: float
        Time after peak (ms)
    max_spikes_per_unit: int
        Maximum number of spikes to display per unit.
    channel_locs: bool
        If True, channel locations are used to display the waveforms.
        If False, waveforms are displayed in vertical order. (default)
    figure: matplotlib figure
        The figure to be used. If not given a figure is created
    ax: matplotlib axis
        The axis to be used. If not given an axis is created

    Returns
    -------
    W: UnitWaveformsWidget
        The output widget
    """
    W = UnitWaveformsWidget(
        recording=recording,
        sorting=sorting,
        channel_ids=channel_ids,
        unit_ids=unit_ids,
        max_spikes_per_unit=max_spikes_per_unit,
        ms_before=ms_before,
        ms_after=ms_after,
        channel_locs=channel_locs,
        figure=figure,
        ax=ax
    )
    W.plot()
    return W


class UnitWaveformsWidget(BaseMultiWidget):
    def __init__(self, *, recording, sorting, channel_ids=None, unit_ids=None, max_spikes_per_unit=50,
                 ms_before=1., ms_after=2., channel_locs=False, figure=None, ax=None):
        BaseMultiWidget.__init__(self, figure, ax)
        self._recording = recording
        self._sorting = sorting
        self._channel_ids = channel_ids
        self._unit_ids = unit_ids
        self._ms_before = ms_before
        self._ms_after = ms_after
        self._max_spikes_per_unit = max_spikes_per_unit
        self._ch_locs = channel_locs
        self.name = 'UnitWaveforms'

    def plot(self):
        self._do_plot()

    def _do_plot(self):
        unit_ids = self._unit_ids
        channel_ids = self._channel_ids
        if unit_ids is None:
            unit_ids = self._sorting.get_unit_ids()
        if channel_ids is None:
            channel_ids = self._recording.get_channel_ids()
        if 'location' in self._recording.get_shared_channel_property_names():
            all_locations = np.array(self._recording.get_channel_locations())
            channel_locations = np.array(self._recording.get_channel_locations(channel_ids=channel_ids))
        else:
            all_locations = None
            channel_locations = None
        list_spikes = []

        for unit_id in unit_ids:
            spiketrain = self._sorting.get_unit_spike_train(unit_id=unit_id)
            if spiketrain is not None:
                random_wf = st.postprocessing.get_unit_waveforms(recording=self._recording, sorting=self._sorting,
                                                                 unit_ids=[unit_id], channel_ids=channel_ids,
                                                                 ms_before=self._ms_before, ms_after=self._ms_after,
                                                                 max_spikes_per_unit=self._max_spikes_per_unit,
                                                                 save_as_features=False)
                random_wf = random_wf.swapaxes(0, 1)
                random_wf = random_wf.swapaxes(1, 2)
                spikes = random_wf
                item = dict(
                    representative_waveforms=spikes,
                    title='Unit {}'.format(int(unit_id))
                )
                list_spikes.append(item)
            else:
                print(unit_id, ' spike train is None')
        with plt.rc_context({'axes.edgecolor': 'gray'}):
            if self._ch_locs:
                self._plot_spike_shapes_multi(list_spikes, channel_locations=channel_locations,
                                              all_locations=all_locations)
            else:
                self._plot_spike_shapes_multi(list_spikes, channel_locations=None)

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
                       all_locations=None, ylim=None, max_representatives=None, color='blue', title=''):
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
        clocs = np.zeros((M, 2))
        for m in range(M):
            clocs[m, :] = [0, -m]
    else:
        if all_locations is not None:
            if all_locations.shape == channel_locations.shape:
                all_locations = channel_locations
        # normalize channel locations
        locs_xrange = np.max(all_locations[:, 0]) - np.min(all_locations[:, 0])
        locs_yrange = np.max(all_locations[:, 1]) - np.min(all_locations[:, 1])
        pitch = [np.max(np.diff(all_locations[:, 0])), np.max(np.diff(all_locations[:, 1]))]

        if locs_xrange > locs_yrange:
            ch_range = locs_xrange
            m_dim = locs_xrange / pitch[0]
        else:
            ch_range = locs_yrange
            m_dim = locs_yrange / pitch[1]

        clocs = channel_locations / ch_range * m_dim

    spacing = 1 / 0.8

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
            loc = clocs[m, -2:] * yrange * spacing
            for j in range(L):
                XX[:, m, j] = loc[0] + xvals
                YY[:, m, j] = loc[1] + representative_waveforms[m, :, j] - representative_waveforms[m, 0, j]
        XX = XX.reshape(T, M * L)
        YY = YY.reshape(T, M * L)
        ax.plot(XX, YY, color=(0.5, 0.5, 0.5), alpha=0.4)

        XX = np.zeros((T, M))
        YY = np.zeros((T, M))
        for m in range(M):
            loc = clocs[m, -2:] * yrange * spacing
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
