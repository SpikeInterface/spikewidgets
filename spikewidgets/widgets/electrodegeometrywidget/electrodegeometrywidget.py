import numpy as np
from matplotlib import pyplot as plt


def plot_electrode_geometry(recording, elec_size=5):
    W = ElectrodeGeometryWidget(
        recording=recording,
        elec_size=elec_size
    )
    W.plot()


class ElectrodeGeometryWidget:
    def __init__(self, *, recording, elec_size=5):
        self._recording = recording
        self._elec_size = elec_size

    def plot(self, width=1.5, height=1.5):
        self._do_plot(width=width, height=height)

    def _do_plot(self, width, height):
        R = self._recording
        geom = np.array(R.get_channel_locations())

        fig = plt.figure(figsize=(width, height))
        ax = fig.add_axes([0, 0, 1, 1])
        ax.axis('off')

        x = geom[:, 0]
        y = geom[:, 1]
        xmin = np.min(x);
        xmax = np.max(x)
        ymin = np.min(y);
        ymax = np.max(y)

        # marker_size=width*fig.dpi/6
        marker_size = self._elec_size
        margin = np.maximum(xmax - xmin, ymax - ymin) * 0.2

        plt.scatter(x, y, marker='o', s=int(marker_size ** 2))
        plt.axis('equal')
        plt.xticks([])
        plt.yticks([])
        plt.xlim(xmin - margin, xmax + margin)
        plt.ylim(ymin - margin, ymax + margin)
        # plt.show()
