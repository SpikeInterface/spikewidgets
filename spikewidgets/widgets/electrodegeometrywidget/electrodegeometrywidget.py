import numpy as np
from matplotlib.patches import Ellipse
from ..utils import LabeledEllipse
from spikewidgets.widgets.basewidget import BaseWidget


def plot_electrode_geometry(recording, color='C0', label_color='r', figure=None, ax=None):
    """
    Plots electrode geometry.

    Parameters
    ----------
    recording: RecordingExtractor
        The recordng extractor object
    color: matplotlib color
        The color of the electrodes
    label_color: matplotlib color
        The color of the channel label when clicking
    figure: matplotlib figure
        The figure to be used. If not given a figure is created
    ax: matplotlib axis
        The axis to be used. If not given an axis is created

    Returns
    -------
    W: UnitWaveformsWidget
        The output widget
    """
    W = ElectrodeGeometryWidget(
        recording=recording,
        color=color,
        label_color=label_color,
        figure=figure,
        ax=ax
    )
    W.plot()
    return W


class ElectrodeGeometryWidget(BaseWidget):
    def __init__(self, *, recording, color='C0', label_color='r', figure=None, ax=None):
        if 'location' not in recording.get_shared_channel_property_names():
            raise AttributeError("'location' not found as a property")
        BaseWidget.__init__(self, figure, ax)
        self._recording = recording
        self._color = color
        self._label_color = label_color
        self.name = 'ElectrodeGeometry'

    def plot(self, width=4, height=4):
        self._do_plot(width=width, height=height)

    def _do_plot(self, width, height):
        locations = np.array(self._recording.get_channel_locations())
        self.ax.axis('off')

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

        self._drs = []
        elec_x = 0.9 * pitch_x
        elec_y = 0.9 * pitch_y
        for (loc, ch) in zip(locations, self._recording.get_channel_ids()):
            ell = Ellipse((loc[0] - elec_x / 2, loc[1] - elec_y / 2), elec_x, elec_y,
                           color=self._color, alpha=0.9)
            self.ax.add_patch(ell)
            dr = LabeledEllipse(ell, ch, self._label_color)
            dr.connect()
            self._drs.append(dr)

        self.ax.axis('equal')
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.set_xlim(np.min(x) - pitch_x, np.max(x) + pitch_x)
        self.ax.set_ylim(np.min(y) - pitch_y, np.max(y) + pitch_y)

