import numpy as np
from matplotlib import pyplot as plt


def plot_electrode_geometry(recording, elec_size=5):
    W = ElectrodeGeometryWidget(
        recording=recording,
        elec_size=elec_size,
    )
    W.plot()


class ElectrodeGeometryWidget:
    def __init__(self, *, recording, elec_size=5, ax=None):
        self._recording = recording
        self._elec_size = elec_size
        self._ax = ax

    def plot(self, width=4, height=4):
        self._do_plot(width=width, height=height)

    def _do_plot(self, width, height):
        R = self._recording
        geom = np.array(R.get_channel_locations())

        if self._ax is None:
            fig = plt.figure(figsize=(width, height))
            self._ax = fig.add_axes([0, 0, 1, 1])
        self._ax.axis('off')

        x = geom[:, 0]
        y = geom[:, 1]
        xmin = np.min(x)
        xmax = np.max(x)
        ymin = np.min(y)
        ymax = np.max(y)

        marker_size = self._elec_size
        margin = np.maximum(xmax - xmin, ymax - ymin) * 0.2

        self._ax.scatter(x, y, marker='o', s=int(marker_size ** 2))
        self._ax.axis('equal')
        self._ax.set_xticks([])
        self._ax.set_yticks([])
        self._ax.set_xlim(xmin - margin, xmax + margin)
        self._ax.set_ylim(ymin - margin, ymax + margin)
        # plt.show()