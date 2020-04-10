import numpy as np
import spiketoolkit as st
import matplotlib.pylab as plt
from spikewidgets.widgets.basewidget import BaseWidget


def plot_activity_map(recording, channel_ids=None, trange=None, cmap='viridis', label_color='r',
                      ax=None, figure=None):
    """
    Plots sorting comparison confusion matrix.

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
        cmap=cmap,
        label_color=label_color,
        figure=figure,
        ax=ax,
    )
    W.plot()
    return W


class ActivityMapWidget(BaseWidget):

    def __init__(self, *, recording, channel_ids, trange, cmap, label_color='r', figure=None, ax=None):
        BaseWidget.__init__(self, figure, ax)
        self._recording = recording
        self._channel_ids = channel_ids
        self._trange = trange
        self._cmap = cmap
        self._label_color = label_color
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
                                                                      end_frame=self._trange[1])
        x = locations[:, 0]
        y = locations[:, 1]
        x_un = np.unique(x)
        y_un = np.unique(y)

        if len(y_un) == 0:
            pitch_x = np.min(np.diff(x_un))
            pitch_y = pitch_x
        elif len(x_un) == 0:
            pitch_y = np.min(np.diff(y_un))
            pitch_x = pitch_y
        else:
            pitch_x = np.min(np.diff(x_un))
            pitch_y = np.min(np.diff(y_un))

        cm = plt.get_cmap(self._cmap)

        self._drs = []
        for (loc, act, ch) in zip(locations, activity, self._recording.get_channel_ids()):
            color = cm(act)
            rect = plt.Rectangle((loc[0], loc[1]), pitch_x - 0.2 * pitch_x, pitch_y - 0.2 * pitch_y,
                                 color=color, edgecolor=None, alpha=0.9)
            self.ax.add_patch(rect)
            dr = LabeledRectangle(rect, ch, self._label_color)
            dr.connect()
            self._drs.append(dr)

        self.ax.set_xlim(np.min(x), np.max(x))
        self.ax.set_ylim(np.min(y), np.max(y))
        self.ax.axis('equal')
        self.ax.axis('off')


class LabeledRectangle:
    lock = None  # only one can be animated at a time

    def __init__(self, rect, channel, color):
        self.rect = rect
        self.press = None
        self.background = None
        self.channel_str = str(channel)
        axes = self.rect.axes
        x0, y0 = self.rect.xy
        self.text = axes.text(x0, y0, self.channel_str, color=color)
        self.text.set_visible(False)

    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.rect.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease = self.rect.figure.canvas.mpl_connect('button_release_event', self.on_release)

    def on_press(self, event):
        'on button press we will see if the mouse is over us and store some data'
        if event.inaxes != self.rect.axes:
            return
        if LabeledRectangle.lock is not None:
            return
        contains, attrd = self.rect.contains(event)
        if not contains: return
        x0, y0 = self.rect.xy
        self.press = x0, y0, event.xdata, event.ydata
        LabeledRectangle.lock = self
        self.text.set_visible(True)
        self.text.draw()

    def on_release(self, event):
        'on release we reset the press data'
        if LabeledRectangle.lock is not self:
            return
        self.press = None
        LabeledRectangle.lock = None
        self.text.set_visible(False)
        self.text.draw()

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.rect.figure.canvas.mpl_disconnect(self.cidpress)
        self.rect.figure.canvas.mpl_disconnect(self.cidrelease)

