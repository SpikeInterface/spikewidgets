import numpy as np
import spiketoolkit as st
import matplotlib.pylab as plt
from ..utils import LabeledRectangle
from spikewidgets.widgets.basewidget import BaseMultiWidget


def plot_unit_template_maps(recording, sorting, channel_ids=None, unit_ids=None, peak='neg', log=False, ncols=10,
                            background='on', cmap='viridis', label_color='r', figure=None, ax=None, axes=None,
                            **templates_kwargs):
    """
    Plots sorting comparison confusion matrix.

    Parameters
    ----------
    recording: RecordingExtractor
        The recordng extractor object
    sorting: SortingExtractor
        The sorting extractor object
    channel_ids: list
        The channel ids to display
    unit_ids: list
        List of unit ids.
    peak: str
        'neg', 'pos' or 'both'
    log: bool
        If True, log scale is used
    ncols: int
        Number of columns if multiple units are displayed
    background: str
        'on' or 'off'
    cmap: matplotlib colormap
        The colormap to be used (default 'viridis')
    label_color: matplotlib color
        Color to display channel name upon click
    figure: matplotlib figure
        The figure to be used. If not given a figure is created
    ax: matplotlib axis
        The axis to be used. If not given an axis is created
    axes: list of matplotlib axes
        The axes to be used for the individual plots. If not given the required axes are created. If provided, the ax
        and figure parameters are ignored
    templates_kwargs: keyword arguments for st.postprocessing.get_unit_templates()


    Returns
    -------
    W: ActivityMapWidget
        The output widget
    """
    W = UnitTemplateMapsWidget(
        recording=recording,
        sorting=sorting,
        channel_ids=channel_ids,
        unit_ids=unit_ids,
        peak=peak,
        log=log,
        ncols=ncols,
        background=background,
        cmap=cmap,
        label_color=label_color,
        figure=figure,
        ax=ax,
        axes=axes,
        **templates_kwargs
    )
    W.plot()
    return W


class UnitTemplateMapsWidget(BaseMultiWidget):
    def __init__(self,  recording, sorting, channel_ids, unit_ids, peak, log, ncols, background, cmap, label_color='r',
                 figure=None, ax=None, axes=None, **template_kwargs):
        BaseMultiWidget.__init__(self, figure, ax, axes)
        self._recording = recording
        self._sorting = sorting
        self._channel_ids = channel_ids
        self._unit_ids = unit_ids
        self._peak = peak
        self._log = log
        self._ncols = ncols
        self._bg = background
        self._cmap = cmap
        self._label_color = label_color
        self._template_kwargs = template_kwargs
        self.name = 'UnitTemplateMaps'
        assert 'location' in self._recording.get_shared_channel_property_names(), "Activity map requires 'location'" \
                                                                                  "property"

    def plot(self):
        self._do_plot()

    def _do_plot(self):
        locations = self._recording.get_channel_locations(channel_ids=self._channel_ids)
        templates = st.postprocessing.get_unit_templates(self._recording, self._sorting, channel_ids=self._channel_ids,
                                                         unit_ids=self._unit_ids,
                                                         **self._template_kwargs)
        if self._channel_ids is None:
            channel_ids = self._recording.get_channel_ids()
        else:
            channel_ids = self._channel_ids
        if self._unit_ids is None:
            unit_ids = self._sorting.get_unit_ids()
        else:
            unit_ids = self._unit_ids
        assert self._peak in ['neg', 'pos', 'both']
        if self._peak == 'min':
            fun = np.min
        elif self._peak == 'max':
            fun = np.max
        else:
            fun = np.ptp

        x = locations[:, 0]
        y = locations[:, 1]
        x_un = np.unique(x)
        y_un = np.unique(y)

        if len(y_un) == 1:
            pitch_x = np.min(np.diff(x_un))
            pitch_y = pitch_x
        elif len(x_un) == 2:
            pitch_y = np.min(np.diff(y_un))
            pitch_x = pitch_y
        else:
            pitch_x = np.min(np.diff(x_un))
            pitch_y = np.min(np.diff(y_un))

        elec_x = 0.9 * pitch_x
        elec_y = 0.9 * pitch_y

        cm = plt.get_cmap(self._cmap)

        if len(templates) <= self._ncols:
            ncols = len(templates)
            nrows = 1
        else:
            ncols = self._ncols
            nrows = np.ceil(len(templates) / ncols)

        for i, (template, unit) in enumerate(zip(templates, unit_ids)):
            ax = self.get_tiled_ax(i, nrows, ncols)
            temp_map = np.abs(fun(template, axis=1))

            if self._log:
                if np.any(temp_map < 1):
                    temp_map += (1 - np.min(temp_map))
                temp_map = np.log(temp_map)

            # normalize
            temp_map -= np.min(temp_map)
            temp_map /= np.ptp(temp_map)

            if self._bg == 'on':
                rect = plt.Rectangle((np.min(x) - pitch_x / 2, np.min(y) - pitch_y / 2),
                                     float(np.ptp(x)) + pitch_x, float(np.ptp(y)) + pitch_y,
                                     color=cm(0), edgecolor=None, alpha=0.9)
                ax.add_patch(rect)

            self._drs = []
            for (loc, tval, ch) in zip(locations, temp_map, channel_ids):
                color = cm(tval)
                rect = plt.Rectangle((loc[0] - elec_x / 2, loc[1] - elec_y / 2), elec_x, elec_y,
                                     color=color, edgecolor=None, alpha=0.9)
                ax.add_patch(rect)
                dr = LabeledRectangle(rect, ch, self._label_color)
                dr.connect()
                self._drs.append(dr)

            ax.set_title(f'Unit {unit}', color='gray')
            ax.set_xlim(np.min(x) - elec_x / 2, np.max(x) + elec_x / 2)
            ax.set_ylim(np.min(y) - elec_y / 2, np.max(y) + elec_y / 2)
            ax.axis('equal')
            ax.axis('off')
