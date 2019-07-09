import numpy as np
from matplotlib import pyplot as plt
from spikewidgets.widgets.basewidget import BaseWidget

def plot_rasters(sorting, sample_rate=None, unit_ids=None, color='k', figure=None, ax=None):
    if sample_rate is None:
        if sorting.get_sampling_frequency() is None:
            raise Exception("Sampling rate information is not in the SortingExtractor. "
                            "Provide the 'sample_rate' argument")
        else:
            sample_rate = sorting.get_sampling_frequency()
    W = ResterWidget(
        sorting=sorting,
        samplerate=sample_rate,
        unit_ids=unit_ids,
        color=color,
        figure=figure,
        ax=ax
    )
    W.plot()
    return W


class ResterWidget(BaseWidget):
    def __init__(self, *, sorting, samplerate, unit_ids=None, color='k', figure=None, ax=None):
        BaseWidget.__init__(self, figure, ax)
        self._SX = sorting
        self._unit_ids = unit_ids
        self._figure = None
        self._samplerate = samplerate
        self._color = color

    def plot(self):
        self._do_plot()

    def figure(self):
        return self._figure

    def _do_plot(self):
        units = self._unit_ids
        if units is None:
            units = self._SX.get_unit_ids()

        min_t = 0
        max_t = 0
        for u_i, unit in enumerate(units):
            t = self._SX.get_unit_spike_train(unit) / float(self._samplerate)

            self.ax.plot(t, u_i * np.ones_like(t), marker='|', mew=1, markersize=3,
                         ls='', color=self._color)
            if np.min(t) < min_t:
                min_t = np.min(t)
            if np.max(t) > max_t:
                max_t = np.max(t)
        self.ax.set_yticks(np.arange(len(units)))
        self.ax.set_yticklabels(units)
        self.ax.set_xlim(min_t, max_t)
