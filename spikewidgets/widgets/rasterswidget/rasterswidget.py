import numpy as np
from matplotlib import pyplot as plt
from spikewidgets.widgets.basewidget import BaseWidget


def plot_rasters(sorting, sampling_frequency=None, unit_ids=None, trange=None, color='k', figure=None, ax=None):
    """
    Plots spike train rasters.

    Parameters
    ----------
    sorting: SortingExtractor
        The sorting extractor object
    sampling_frequency: float
        The sampling frequency (if not in the sorting extractor)
    unit_ids: list
        List of unit ids
    trange: list
        List with start time and end time
    color: matplotlib color
        The color to be used
    figure: matplotlib figure
        The figure to be used. If not given a figure is created
    ax: matplotlib axis
        The axis to be used. If not given an axis is created

    Returns
    -------
    W: RasterWidget
        The output widget
    """
    if sampling_frequency is None:
        if sorting.get_sampling_frequency() is None:
            raise Exception("Sampling rate information is not in the SortingExtractor. "
                            "Provide the 'sampling_frequency' argument")
        else:
            sampling_frequency = sorting.get_sampling_frequency()
    W = RasterWidget(
        sorting=sorting,
        sampling_frequency=sampling_frequency,
        unit_ids=unit_ids,
        trange=trange,
        color=color,
        figure=figure,
        ax=ax
    )
    W.plot()
    return W


class RasterWidget(BaseWidget):
    def __init__(self, *, sorting, sampling_frequency, unit_ids=None, trange=None, color='k', figure=None, ax=None):
        BaseWidget.__init__(self, figure, ax)
        self._SX = sorting
        self._unit_ids = unit_ids
        self._figure = None
        self._sampling_frequency = sampling_frequency
        self._color = color
        self._max_frame = 0
        for unit_id in self._SX.get_unit_ids():
            spike_train = self._SX.get_unit_spike_train(unit_id)
            curr_max_frame = np.max(spike_train)
            if curr_max_frame > self._max_frame:
                self._max_frame = curr_max_frame
        self._visible_trange = trange
        if self._visible_trange is None:
            self._visible_trange = [0, self._max_frame]
        else:
            assert len(trange) == 2, "'trange' should be a list with start and end time in seconds"
            self._visible_trange = [int(t * sampling_frequency) for t in trange]
        self._visible_trange = self._fix_trange(self._visible_trange)
        self.name = 'Raster'

    def plot(self):
        self._do_plot()

    def _do_plot(self):
        units_ids = self._unit_ids
        if units_ids is None:
            units_ids = self._SX.get_unit_ids()

        with plt.rc_context({'axes.edgecolor': 'gray'}):
            for u_i, unit_id in enumerate(units_ids):
                t = self._SX.get_unit_spike_train(unit_id, start_frame=self._visible_trange[0], end_frame=self._visible_trange[1]) / float(self._sampling_frequency)
                self.ax.plot(t, u_i * np.ones_like(t), marker='|', mew=1, markersize=3,
                             ls='', color=self._color)
            visible_start_frame = self._visible_trange[0] / float(self._sampling_frequency)
            visible_end_frame = self._visible_trange[1] / float(self._sampling_frequency)
            self.ax.set_yticks(np.arange(len(units_ids)))
            self.ax.set_yticklabels(units_ids)
            self.ax.set_xlim(visible_start_frame, visible_end_frame)
            self.ax.set_xlabel('time (s)')

    def _fix_trange(self, trange):
        if trange[1] > self._max_frame:
            # trange[0] += max_t - trange[1]
            trange[1] = self._max_frame
        if trange[0] < 0:
            # trange[1] += -trange[0]
            trange[0] = 0
        # trange[0] = np.maximum(0, trange[0])
        # trange[1] = np.minimum(max_t, trange[1])
        return trange
