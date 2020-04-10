import numpy as np
from matplotlib import pyplot as plt
import spiketoolkit as st
import MEAutility as mu
from spikewidgets.widgets.basewidget import BaseMultiWidget


def plot_unit_waveforms(recording, sorting, channel_ids=None, unit_ids=None, ms_before=1., ms_after=2.,
                        max_spikes_per_unit=100, max_channels=16, channel_locs=True,
                        plot_templates=True, figure=None, ax=None):
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
    max_channels: int
        Maximum number of largest channels to plot waveform
    channel_locs: bool
        If True, channel locations are used to display the waveforms.
        If False, waveforms are displayed in vertical order. (default)
    plot_templates: bool
        If True, templates are plotted over the waveforms
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
        max_channels=max_channels,
        ms_before=ms_before,
        ms_after=ms_after,
        channel_locs=channel_locs,
        plot_templates=plot_templates,
        figure=figure,
        ax=ax
    )
    W.plot()
    return W


def plot_unit_templates(recording, sorting, channel_ids=None, unit_ids=None, ms_before=1., ms_after=2.,
                        max_spikes_per_unit=100, max_channels=16, channel_locs=True,
                        figure=None, ax=None):
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
    max_channels: int
        Maximum number of largest channels to plot waveform
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
    W = UnitTemplatesWidget(
        recording=recording,
        sorting=sorting,
        channel_ids=channel_ids,
        unit_ids=unit_ids,
        max_spikes_per_unit=max_spikes_per_unit,
        max_channels=max_channels,
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
                 max_channels=16, ms_before=1., ms_after=2., channel_locs=True, plot_templates=True,
                 figure=None, ax=None):
        BaseMultiWidget.__init__(self, figure, ax)
        self._recording = recording
        self._sorting = sorting
        self._channel_ids = channel_ids
        self._max_channels = max_channels
        self._unit_ids = unit_ids
        self._ms_before = ms_before
        self._ms_after = ms_after
        self._max_spikes_per_unit = max_spikes_per_unit
        self._ch_locs = channel_locs
        self.name = 'UnitWaveforms'
        self._plot_waveforms = True
        self._plot_templates = plot_templates

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
        max_channels_list = []
        for unit_id in unit_ids:
            spiketrain = self._sorting.get_unit_spike_train(unit_id=unit_id)
            if spiketrain is not None:
                random_wf = st.postprocessing.get_unit_waveforms(recording=self._recording, sorting=self._sorting,
                                                                 unit_ids=[unit_id], channel_ids=channel_ids,
                                                                 ms_before=self._ms_before, ms_after=self._ms_after,
                                                                 max_spikes_per_unit=self._max_spikes_per_unit,
                                                                 save_as_features=False)[0]
                if self._max_channels < self._recording.get_num_channels():
                    template = np.mean(random_wf, axis=0)
                    peak_idx = np.unravel_index(np.argmax(np.abs(template)),
                                                template.shape)[1]
                    max_channel_idxs = np.argsort(np.abs(template[:, peak_idx]))[::-1][:self._max_channels]
                    max_channels_list.append(max_channel_idxs)
                else:
                    max_channels_list.append(np.arange(len(channel_ids)))
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
                                              all_locations=all_locations, max_channels_list=max_channels_list,
                                              plot_templates=self._plot_templates, plot_waveforms=self._plot_waveforms)
            else:
                self._plot_spike_shapes_multi(list_spikes, channel_locations=None, max_channels_list=max_channels_list,
                                              plot_templates=self._plot_templates, plot_waveforms=self._plot_waveforms)

    def _plot_spike_shapes_multi(self, list_spikes, max_channels_list, *, ncols=5, **kwargs):
        vscale, ylim = _determine_global_vscale_ylim(list_spikes)
        if len(list_spikes) < ncols:
            ncols = len(list_spikes)
        nrows = np.ceil(len(list_spikes) / ncols)
        for i, item in enumerate(list_spikes):
            ax = self.get_tiled_ax(i, nrows, ncols)
            _plot_spike_shapes(**item, **kwargs, ax=ax, channels=max_channels_list[i], vscale=vscale, ylim_wf=ylim)


class UnitTemplatesWidget(UnitWaveformsWidget):
    def __init__(self, *, recording, sorting, channel_ids=None, unit_ids=None, max_spikes_per_unit=50,
                 max_channels=16, ms_before=1., ms_after=2., channel_locs=True, figure=None, ax=None):
        UnitWaveformsWidget.__init__(self, recording=recording, sorting=sorting, channel_ids=channel_ids,
                                     unit_ids=unit_ids, max_spikes_per_unit=max_spikes_per_unit,
                                     max_channels=max_channels, ms_before=ms_before, ms_after=ms_after,
                                     channel_locs=channel_locs, plot_templates=True, figure=figure, ax=ax)
        self.name = 'UnitTemplates'
        self._plot_waveforms = False


def _plot_spike_shapes(*, ax, channels, representative_waveforms=None, channel_locations=None,
                       all_locations=None, color='blue', vscale=None, plot_waveforms=True, plot_templates=True,
                       ylim_wf=None, title=''):
    if representative_waveforms is None:
        raise Exception('You must provide either average_waveform, representative waveforms, or both')
    waveforms = representative_waveforms
    average_waveform = np.mean(waveforms, axis=0)

    M = average_waveform.shape[0]  # number of channels
    if channel_locations is None:
        channel_locations = np.zeros((M, 2))
        for m in range(M):
            channel_locations[m, :] = [0, -m]
        all_locations = channel_locations
    else:
        if all_locations is not None:
            if all_locations.shape == channel_locations.shape:
                all_locations = channel_locations

    probe = mu.return_mea(info={'pos': all_locations, 'center': False})

    if vscale is None:
        vscale = 1.5 * np.max(np.abs(average_waveform))

    ylim_g = np.array([np.min(all_locations[:, 1]), np.max(all_locations[:, 1])])
    ptp = ylim_g[1] - ylim_g[0]

    # fix for linear horizontal layout
    if ptp == 0:
        ylim = ylim_wf / vscale
    else:
        ylim = [ylim_g[0] - 0.1 * ptp, ylim_g[1] + 0.1 * ptp]

    if plot_waveforms:
        ax = mu.plot_mea_recording(waveforms, probe, colors=(0.5, 0.5, 0.5), alpha=0.3, lw=0.3, ax=ax, vscale=vscale,
                                   ylim=ylim, channels=channels)
    if plot_templates:
        ax = mu.plot_mea_recording(average_waveform, probe, colors=color, ax=ax, lw=2, vscale=vscale, ylim=ylim,
                                   channels=channels)

    ax.set_xticks([])
    ax.set_yticks([])
    if title:
        ax.set_title(title, color='gray')


def _get_vscale_ylim_for_item(average_waveform=None, representative_waveforms=None):
    if average_waveform is None:
        if representative_waveforms is None:
            raise Exception('You must provide either average_waveform, representative waveforms, or both')
        average_waveform = np.mean(representative_waveforms, axis=0)
    vscale = np.max(np.abs(average_waveform))
    ylim = [np.min(average_waveform), np.max(average_waveform)]
    return vscale, ylim


def _determine_global_vscale_ylim(list_spikes):
    vscales = []
    ylims = []
    for item in list_spikes:
        vscale, ylim = _get_vscale_ylim_for_item(
            average_waveform=item.get('average_waveform', None),
            representative_waveforms=item.get('representative_waveforms', None)
        )
        vscales.append(vscale)
        ylims.append(ylim)
    vscale_global = np.median(vscales) * 1.1
    ylims = np.array(ylims)
    ylim_global = [np.min(ylims[:, 0]), np.max(ylims[:, 1])]
    ylim_global = [ylim_global[0] - 0.1 * np.ptp(ylim_global), ylim_global[1] + 0.1 * np.ptp(ylim_global)]
    return vscale_global, ylim_global
