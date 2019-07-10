import numpy as np
from matplotlib import pyplot as plt
from spikewidgets.widgets.basewidget import BaseWidget


def plot_electrode_geometry(recording, elec_size=5, figure=None, ax=None):
    W = ElectrodeGeometryWidget(
        recording=recording,
        elec_size=elec_size,
        figure=figure,
        ax=ax
    )
    W.plot()
    return W


class ElectrodeGeometryWidget(BaseWidget):
    def __init__(self, *, recording, elec_size=5, figure=None, ax=None):
        BaseWidget.__init__(self, figure, ax)
        self._recording = recording
        self._elec_size = elec_size
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

        marker_size = self._elec_size
        margin = np.maximum(xmax - xmin, ymax - ymin) * 0.2

        self.ax.scatter(x, y, marker='o', s=int(marker_size ** 2))
        self.ax.axis('equal')
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.set_xlim(xmin - margin, xmax + margin)
        self.ax.set_ylim(ymin - margin, ymax + margin)
