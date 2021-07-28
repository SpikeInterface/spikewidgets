import numpy as np
from matplotlib import pyplot as plt
import spiketoolkit as st
import MEAutility as mu
from spikewidgets.widgets.basewidget import BaseMultiWidget


def plot_unit_waveforms(recording, sorting, channel_ids=None, unit_ids=None, channel_locs=True, radius=None,
                        max_channels=None, plot_templates=True, show_all_channels=True, color='k', lw=2,
                        axis_equal=False,
                        plot_channels=False, set_title=True, figure=None, ax=None, axes=None, **waveforms_kwargs):
    """
    Plots unit waveforms.

    Parameters
    ----------
    recording: RecordingExtractor
        The recording extractor object
    sorting: SortingExtractor
        The sorting extractor object
    channel_ids: list
        The channel ids to display
    unit_ids: list
        List of unit ids.
    max_channels: int
        Maximum number of largest channels to plot waveform
    channel_locs: bool
        If True, channel locations are used to display the waveforms.
        If False, waveforms are displayed in vertical order (default)
    plot_templates: bool
        If True, templates are plotted over the waveforms
    radius: float
        If not None, all channels within a circle around the peak waveform will be displayed
        Ignores max_spikes_per_unit
    set_title: bool
        Create a plot title with the unit number if True.
    plot_channels: bool
        Plot channel locations below traces, only used if channel_locs is True
    axis_equal: bool
        Equal aspext ratio for x and y axis, to visualise the array geometry to scale
    lw: float
        Line width for the traces.
    color: matplotlib color or list of colors
        Color(s) of traces.
    show_all_channels: bool
        Show the whole probe if True, or only selected channels if False
    figure: matplotlib figure
        The figure to be used. If not given a figure is created
    ax: matplotlib axis
        The axis to be used. If not given an axis is created
    axes: list of matplotlib axes
        The axes to be used for the individual plots. If not given the required axes are created. If provided, the ax
        and figure parameters are ignored
    waveforms_kwargs: keyword arguments for st.postprocessing.get_unit_waveforms()

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
        max_channels=max_channels,
        channel_locs=channel_locs,
        plot_templates=plot_templates,
        figure=figure,
        ax=ax,
        axes=axes,
        radius=radius,
        show_all_channels=show_all_channels,
        color=color,
        lw=lw,
        axis_equal=axis_equal,
        plot_channels=plot_channels,
        set_title=set_title,
        **waveforms_kwargs
    )
    W.plot()
    return W


def plot_unit_templates(recording, sorting, channel_ids=None, unit_ids=None, max_channels=None, channel_locs=True,
                        radius=None, show_all_channels=True, color='k', lw=2, axis_equal=False,
                        plot_channels=False, set_title=True, figure=None, ax=None, axes=None, **template_kwargs):
    """
    Plots unit waveforms.

    Parameters
    ----------
    recording: RecordingExtractor
        The recording extractor object
    sorting: SortingExtractor
        The sorting extractor object
    channel_ids: list
        The channel ids to display
    unit_ids: list
        List of unit ids.
    max_channels: int
        Maximum number of largest channels to plot waveform
    channel_locs: bool
        If True, channel locations are used to display the waveforms.
        If False, waveforms are displayed in vertical order. (default)
    radius: float
        If not None, all channels within a circle around the peak waveform will be displayed.
        Ignores max_spikes_per_unit
    set_title: bool
        Create a plot title with the unit number if True
    plot_channels: bool
        Plot channel locations below traces, only used if channel_locs is True
    axis_equal: bool
        Equal aspext ratio for x and y axis, to visualise the array geometry to scale
    lw: float
        Line width for the traces.
    color: matplotlib color or list of colors
        Color(s) of traces.
    show_all_channels: bool
        Show the whole probe if True, or only selected channels if False
    figure: matplotlib figure
        The figure to be used. If not given a figure is created
    ax: matplotlib axis
        The axis to be used. If not given an axis is created
    axes: list of matplotlib axes
        The axes to be used for the individual plots. If not given the required axes are created. If provided, the ax
        and figure parameters are ignored
    template_kwargs: keyword arguments for st.postprocessing.get_unit_templates()


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
        max_channels=max_channels,
        channel_locs=channel_locs,
        figure=figure,
        ax=ax,
        axes=axes,
        radius=radius,
        show_all_channels=show_all_channels,
        color=color,
        lw=lw,
        axis_equal=axis_equal,
        plot_channels=plot_channels,
        set_title=set_title,
        **template_kwargs
    )
    W.plot()
    return W


class UnitWaveformsWidget(BaseMultiWidget):
    def __init__(self, *, recording, sorting, channel_ids=None, unit_ids=None, max_channels=None,
                 channel_locs=True, plot_waveforms=True,
                 plot_templates=True, radius=None,
                 show_all_channels=True, figure=None, ax=None, axes=None, color='k', lw=2, axis_equal=False,
                 plot_channels=False, set_title=True, **kwargs):
        BaseMultiWidget.__init__(self, figure, ax, axes)
        self._recording = recording
        self._sorting = sorting
        self._channel_ids = channel_ids
        if max_channels is None:
            max_channels = recording.get_num_channels()
        self._max_channels = max_channels
        self._unit_ids = unit_ids
        self._ch_locs = channel_locs
        self.name = 'UnitWaveforms'
        self._plot_waveforms = plot_waveforms
        self._plot_templates = plot_templates
        self._radius = radius
        self._show_all_channels = show_all_channels
        self._color = color
        self._lw = lw
        self._axis_equal = axis_equal
        self._plot_channels = plot_channels
        self._set_title = set_title
        self._kwargs = kwargs

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
                if self._plot_waveforms:
                    random_wf = st.postprocessing.get_unit_waveforms(recording=self._recording, sorting=self._sorting,
                                                                     unit_ids=[unit_id], channel_ids=channel_ids,
                                                                     **self._kwargs)[0]
                else:
                    random_wf = None
                template = st.postprocessing.get_unit_templates(recording=self._recording,
                                                                sorting=self._sorting,
                                                                unit_ids=[unit_id], channel_ids=channel_ids,
                                                                **self._kwargs)[0]
                if self._radius is None:
                    if self._max_channels < template.shape[0]:
                        peak_idx = np.unravel_index(np.argmax(np.abs(template)),
                                                    template.shape)[1]
                        max_channel_idxs = np.argsort(np.abs(template[:, peak_idx]))[::-1][:self._max_channels]
                        max_channels_list.append(max_channel_idxs)
                    else:
                        if 'waveforms_channel_idxs' in self._sorting.get_unit_property_names(unit_id):
                            max_channels_list.append(self._sorting.get_unit_property(unit_id, 'waveforms_channel_idxs'))
                        else:
                            max_channels_list.append(np.arange(len(channel_ids)))
                else:
                    peak_idx = np.unravel_index(np.argmax(np.abs(template)),
                                                template.shape)[1]
                    max_channel_idx = np.argsort(np.abs(template[:, peak_idx]))[::-1][0]
                    center_location = all_locations[max_channel_idx]
                    dists = np.array([np.linalg.norm(loc - center_location) for loc in all_locations])
                    max_channels_list.append(np.where(dists <= self._radius)[0])
                spikes = random_wf
                if self._set_title:
                    item = dict(
                        representative_waveforms=spikes,
                        average_waveform=template,
                        title='Unit {}'.format(int(unit_id))
                    )
                else:
                    item = dict(
                        representative_waveforms=spikes,
                        average_waveform=template,
                        title=''
                    )
                list_spikes.append(item)
            else:
                print(unit_id, ' spike train is None')

        with plt.rc_context({'axes.edgecolor': 'gray'}):
            if self._ch_locs:
                self._plot_spike_shapes_multi(list_spikes, channel_locations=channel_locations,
                                              all_locations=all_locations, max_channels_list=max_channels_list,
                                              plot_templates=self._plot_templates, plot_waveforms=self._plot_waveforms,
                                              show_all_channels=self._show_all_channels,
                                              plot_channels=self._plot_channels)
            else:
                self._plot_spike_shapes_multi(list_spikes, channel_locations=None, max_channels_list=max_channels_list,
                                              plot_templates=self._plot_templates, plot_waveforms=self._plot_waveforms,
                                              show_all_channels=self._show_all_channels,
                                              plot_channels=self._plot_channels)

    def _plot_spike_shapes_multi(self, list_spikes, max_channels_list, *, ncols=5, **kwargs):
        vscale, ylim = _determine_global_vscale_ylim(list_spikes)
        if len(list_spikes) < ncols:
            ncols = len(list_spikes)
        nrows = np.ceil(len(list_spikes) / ncols)
        for i, item in enumerate(list_spikes):
            ax = self.get_tiled_ax(i, nrows, ncols)
            _plot_spike_shapes(**item, **kwargs, ax=ax, channels=max_channels_list[i], vscale=vscale, ylim_wf=ylim,
                               color=self._color, lw=self._lw)


class UnitTemplatesWidget(UnitWaveformsWidget):
    def __init__(self, *, recording, sorting, channel_ids=None, unit_ids=None,
                 max_channels=None, channel_locs=True, figure=None, show_all_channels=True,
                 ax=None, axes=None, radius=None, color='k', lw=2, axis_equal=False, plot_channels=False,
                 set_title=True, **template_kwargs):
        UnitWaveformsWidget.__init__(self, recording=recording, sorting=sorting, channel_ids=channel_ids,
                                     unit_ids=unit_ids, max_channels=max_channels,
                                     channel_locs=channel_locs, plot_waveforms=False,
                                     plot_templates=True, figure=figure, ax=ax, axes=axes,
                                     radius=radius, show_all_channels=show_all_channels, color=color, lw=lw,
                                     axis_equal=axis_equal, plot_channels=plot_channels, set_title=set_title,
                                     **template_kwargs)
        self.name = 'UnitTemplates'


def _plot_spike_shapes(*, ax, channels, representative_waveforms=None, average_waveform=None,
                       channel_locations=None, all_locations=None, color='blue', vscale=None, plot_waveforms=True,
                       plot_templates=True, ylim_wf=None, title='', show_all_channels=True, lw=2, axis_equal=False,
                       plot_channels=False):
    if representative_waveforms is None and average_waveform is None:
        raise Exception('You must provide either average_waveform, representative waveforms, or both')

    waveforms = representative_waveforms
    average_waveform = average_waveform
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

    if show_all_channels:
        ylim_g = np.array([np.min(all_locations[:, 1]), np.max(all_locations[:, 1])])
    else:
        ylim_g = np.array([np.min(all_locations[channels, 1]), np.max(all_locations[channels, 1])])
    ptp = ylim_g[1] - ylim_g[0]

    # fix for linear horizontal layout
    if ptp == 0:
        ylim = ylim_wf / vscale
    else:
        ylim = [ylim_g[0] - 0.2 * ptp, ylim_g[1] + 0.2 * ptp]

    if plot_channels and channel_locations is not None:
        if show_all_channels:
            mu.plot_probe(probe, type='planar', ax=ax, alpha_prb=0, alpha_elec=0.2)
        else:
            mu.plot_probe(mu.return_mea(info={'pos': channel_locations[channels], 'center': False}),
                          type='planar', ax=ax, alpha_prb=0, alpha_elec=0.2)

    if plot_waveforms:
        ax = mu.plot_mea_recording(waveforms, probe, colors=(0.5, 0.5, 0.5), alpha=0.3, lw=lw, ax=ax, vscale=vscale,
                                   ylim=ylim, channels=channels, axis_equal=axis_equal)
    if plot_templates:
        ax = mu.plot_mea_recording(average_waveform, probe, colors=color, ax=ax, lw=lw, vscale=vscale, ylim=ylim,
                                   channels=channels, axis_equal=axis_equal, hide_axis=True)

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
