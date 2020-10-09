import numpy as np
import spiketoolkit as st
import matplotlib.pylab as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from ..utils import LabeledRectangle
from spikewidgets.widgets.basewidget import BaseWidget
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def plot_activity_map(recording, channel_ids=None, trange=None, cmap='viridis', background='on', label_color='r',
                      transpose=False, frame=False, ax=None, figure=None, colorbar=False, recompute_info=False):
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
    transpose: bool, optional, default: False
        Swap x and y channel coordinates if True.
    frame: bool, optional, default: False
        Draw a frame around the array if True.
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
        recompute_info=recompute_info
    )
    W.plot()
    return W


class ActivityMapWidget(BaseWidget):

    def __init__(self, recording, channel_ids, trange, cmap, background, label_color='r', transpose=False, frame=False,
                 figure=None, ax=None, colorbar=False, recompute_info=False):
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
                                                                      normalize=False, 
                                                                      recompute_info=self._recompute_info,
                                                                      method='detection')
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
        activity = activity/(self._trange[1]-self._trange[0]) # spikes/s
        max_activity = np.ceil(np.max(activity)*100)/100 # normalise to 1 for colours
        activity = activity/max_activity
        for (loc, act, ch) in zip(locations, activity, self._recording.get_channel_ids()):
            color = cm(act)
            rect = plt.Rectangle((loc[0] - elec_x / 2, loc[1] - elec_y / 2), elec_x, elec_y,
                                 color=color, edgecolor=None, alpha=0.9)
            self.ax.add_patch(rect)
            dr = LabeledRectangle(rect, ch, self._label_color)
            dr.connect()
            self._drs.append(dr)

        self.ax.set_xlim(np.min(x) - pitch_x, np.max(x) + pitch_x)
        self.ax.set_ylim(np.min(y) - pitch_y, np.max(y) + pitch_y)
        if self._frame:
            rect = plt.Rectangle((np.min(x) - pitch_x, np.min(y) - pitch_y), np.max(x) - np.min(x) + 2 * pitch_x,
                                 np.max(y) - np.min(y) + 2 * pitch_y, fill=None, edgecolor='k')
            self.ax.add_patch(rect)
        self.ax.axis('equal')
        self.ax.axis('off')
        if self._show_colorbar:
            cax = inset_axes(self.ax, width="0.5%", height="50%", loc='upper left', bbox_to_anchor=(0.02, 0., 1, 1),bbox_transform=self.ax.transAxes)
            self.colorbar = plt.gcf().colorbar(mpl.collections.PatchCollection(self.ax.patches), cax=cax, orientation='vertical', shrink=0.2)
            cax.yaxis.set_ticks_position('left')
            cax.yaxis.set_label_position('left')
            self.colorbar.set_ticks((0,1))
            self.colorbar.set_ticklabels((0,max_activity))