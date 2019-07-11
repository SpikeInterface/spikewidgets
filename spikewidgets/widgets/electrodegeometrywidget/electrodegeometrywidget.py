import numpy as np
from spikewidgets.widgets.basewidget import BaseWidget


def plot_electrode_geometry(recording, markersize=20, marker='o', figure=None, ax=None):
    """
    Plots electrode geometry.

    Parameters
    ----------
    recording: RecordingExtractor
        The recordng extractor object
    markersize: int
        The size of the marker for the electrodes
    marker: str
        The matplotlib marker to use (default 'o')
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
        markersize=markersize,
        marker=marker,
        figure=figure,
        ax=ax
    )
    W.plot()
    return W


class ElectrodeGeometryWidget(BaseWidget):
    def __init__(self, *, recording, markersize=10, marker='o', figure=None, ax=None):
        BaseWidget.__init__(self, figure, ax)
        self._recording = recording
        self._ms = markersize
        self._mark = marker
        self.name = 'ElectrodeGeometry'

    def plot(self, width=4, height=4):
        self._do_plot(width=width, height=height)

    def _do_plot(self, width, height):
        R = self._recording
        geom = np.array(R.get_channel_locations())


        self.ax.axis('off')

        x = geom[:, 0]
        y = geom[:, 1]
        xmin = np.min(x)
        xmax = np.max(x)
        ymin = np.min(y)
        ymax = np.max(y)

        margin = np.maximum(xmax - xmin, ymax - ymin) * 0.2

        self.ax.scatter(x, y, marker=self._mark, s=int(self._ms))
        self.ax.axis('equal')
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.set_xlim(xmin - margin, xmax + margin)
        self.ax.set_ylim(ymin - margin, ymax + margin)
