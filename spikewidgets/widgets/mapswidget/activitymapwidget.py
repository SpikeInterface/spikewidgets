import numpy as np
import spiketoolkit as st
import matplotlib.pylab as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from ..utils import LabeledRectangle
from spikewidgets.widgets.basewidget import BaseWidget
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def plot_activity_map(recording, channel_ids=None, trange=None, cmap='viridis', background='on', label_color='r',
                      transpose=False, frame=False, colorbar=False, colorbar_bbox=None,
                      colorbar_orientation='vertical', colorbar_width=0.02, recompute_info=False,
                      ax=None, figure=None):
    """
    Plots spike rate (estimated using simple threshold detector) as 2D activity map.

    Parameters
    ----------
    recording: RecordingExtractor
        The recordng extractor object
    channel_ids: list
        The channel ids to display.
    trange: list
        List with start time and end time
    cmap: matplotlib colormap
        The colormap to be used (default 'viridis')
    background: bool
        If True, a background is added in between electrodes
    transpose: bool, optional, default: False
        Swap x and y channel coordinates if True.
    frame: bool, optional, default: False
        Draw a frame around the array if True.
    colorbar: bool
        If True, a colorbar is displayed
    colorbar_bbox: bbox
        Bbox (x,y,w,h) in figure coordinates to plot colorbar
    colorbar_orientation: str
        'vertical' or 'horizontal'
    colorbar_width: float
        Width of colorbar in figure coordinates (default 0.02)
    recompute_info: bool
        If True, spike rates are recomputed
    figure: matplotlib figure
        The figure to be used. If not given a figure is created
    ax: matplotlib axis
        The axis to be used. If not given an axis is created

    Returns
    -------
    W: ActivityMapWidget
        The output widget
    """
    W = ActivityMapWidget(
        recording=recording,
        channel_ids=channel_ids,
        trange=trange,
        background=background,
        cmap=cmap,
        label_color=label_color,
        transpose=transpose,
        frame=frame,
        figure=figure,
        ax=ax,
        colorbar=colorbar,
        colorbar_bbox=colorbar_bbox,
        colorbar_orientation=colorbar_orientation,
        colorbar_width=colorbar_width,
        recompute_info=recompute_info
    )
    W.plot()
    return W


class ActivityMapWidget(BaseWidget):
    def __init__(self, recording, channel_ids, trange, cmap, background, label_color='r', transpose=False, frame=False,
                 colorbar=False, colorbar_bbox=None, colorbar_orientation='vertical', colorbar_width=0.02,
                 recompute_info=False, figure=None, ax=None):
        BaseWidget.__init__(self, figure, ax)
        self._recording = recording
        self._channel_ids = channel_ids
        self._trange = trange
        self._transpose = transpose
        self._cmap = cmap
        self._frame = frame
        self._bg = background
        self._label_color = label_color
        self._show_colorbar = colorbar
        self._colorbar_bbox = colorbar_bbox
        self._colorbar_orient = colorbar_orientation
        self._colorbar_width = colorbar_width
        self._recompute_info = recompute_info
        self.colorbar = None
        self.name = 'ActivityMap'
        assert 'location' in self._recording.get_shared_channel_property_names(), "Activity map requires 'location'" \
                                                                                  "property"

    def plot(self):
        self._do_plot()

    def _do_plot(self):
        if self._trange is None:
            self._trange = [0, self._recording.get_num_frames()]
        else:
            assert len(self._trange) == 2, "'trange' should be a list with start and end time in seconds"
            self._trange = [int(t * self._recording.get_sampling_frequency()) for t in self._trange]

        locations = self._recording.get_channel_locations(channel_ids=self._channel_ids)
        activity = st.postprocessing.compute_channel_spiking_activity(self._recording,
                                                                      start_frame=self._trange[0],
                                                                      end_frame=self._trange[1],
                                                                      method='detection',
                                                                      align=False,
                                                                      recompute_info=self._recompute_info,
                                                                      verbose=True)
        if self._transpose:
            locations = np.roll(locations, 1, axis=1)

        x = locations[:, 0]
        y = locations[:, 1]
        x_un = np.unique(x)
        y_un = np.unique(y)

        if len(y_un) == 1:
            pitch_x = np.min(np.diff(x_un))
            pitch_y = pitch_x
        elif len(x_un) == 1:
            pitch_y = np.min(np.diff(y_un))
            pitch_x = pitch_y
        else:
            pitch_x = np.min(np.diff(x_un))
            pitch_y = np.min(np.diff(y_un))

        cm = plt.get_cmap(self._cmap)

        if self._bg == 'on':
            rect = plt.Rectangle((np.min(x) - pitch_x / 2, np.min(y) - pitch_y / 2),
                                 float(np.ptp(x)) + pitch_x, float(np.ptp(y)) + pitch_y,
                                 color=cm(0), edgecolor=None, alpha=0.9)
            self.ax.add_patch(rect)

        self._drs = []
        elec_x = 0.9 * pitch_x
        elec_y = 0.9 * pitch_y

        max_activity = np.round(np.max(activity), 2)

        # normalize
        activity -= np.min(activity)
        activity /= np.ptp(activity)

        for (loc, act, ch) in zip(locations, activity, self._recording.get_channel_ids()):
            color = cm(act)
            rect = plt.Rectangle((loc[0] - elec_x / 2, loc[1] - elec_y / 2), elec_x, elec_y,
                                 color=color, edgecolor=None, alpha=0.9)
            self.ax.add_patch(rect)
            dr = LabeledRectangle(rect, ch, self._label_color)
            dr.connect()
            self._drs.append(dr)

        if self._frame:
            rect = plt.Rectangle((np.min(x) - pitch_x, np.min(y) - pitch_y), np.max(x) - np.min(x) + 2 * pitch_x,
                                 np.max(y) - np.min(y) + 2 * pitch_y, fill=None, edgecolor='k')
            self.ax.add_patch(rect)

        self.ax.axis('equal')
        self.ax.axis('off')
        if self._show_colorbar:
            # The canvas need to be drawn to get the right transforms
            self.figure.canvas.draw()

            if self._colorbar_bbox is None:
                if self._colorbar_orient == 'vertical':
                    colorbar_width = self._colorbar_width
                    bottom_left = (np.min(x) - pitch_x, np.min(y) - pitch_y)
                    top_left = (np.min(x) - pitch_x, np.max(y) + pitch_y)
                    figure_to_data = self.figure.transFigure + self.ax.transData.inverted()
                    data_to_figure = figure_to_data.inverted()
                    corners_figure = data_to_figure.transform([bottom_left, top_left])

                    bottom_left_fig = corners_figure[0]
                    colorbar_height = (corners_figure[1] - corners_figure[0])[1]
                    bbox = (bottom_left_fig[0] - 1.5 * colorbar_width, bottom_left_fig[1],
                            colorbar_width, colorbar_height)
                else:
                    colorbar_height = self._colorbar_width
                    bottom_left = (np.min(x) - pitch_x, np.min(y) - pitch_y)
                    bottom_right = (np.max(x) + pitch_x, np.min(y) - pitch_y)
                    figure_to_data = self.figure.transFigure + self.ax.transData.inverted()
                    data_to_figure = figure_to_data.inverted()
                    corners_figure = data_to_figure.transform([bottom_left, bottom_right])

                    bottom_left_fig = corners_figure[0]
                    colorbar_width = (corners_figure[1] - corners_figure[0])[0]

                    bbox = (bottom_left_fig[0], bottom_left_fig[1] - 1.5*colorbar_height,
                            colorbar_width, colorbar_height)

            else:
                bbox = self._colorbar_bbox

            cax = self.figure.add_axes(bbox)
            scalable = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=1), cmap=self._cmap)
            self.colorbar = self.figure.colorbar(scalable, cax=cax,
                                                 orientation=self._colorbar_orient)#, shrink=0.5)
            cax.yaxis.set_ticks_position('left')
            cax.yaxis.set_label_position('left')
            self.colorbar.set_ticks((0, 1))
            self.colorbar.set_ticklabels((0, max_activity))


