import numpy as np
from matplotlib import pyplot as plt
from spikewidgets.widgets.basewidget import BaseMultiWidget


def plot_isi_distribution(sorting, sampling_frequency=None, unit_ids=None, bins=10, window=1, figure=None, ax=None):
    """
    Plots spike train ISI distribution.

    Parameters
    ----------
    sorting: SortingExtractor
        The sorting extractor object
    sampling_frequency: float
        The sampling frequency (if not in the sorting extractor)
    unit_ids: list
        List of unit ids
    bins: int
        Number of bins
    window: float
        Window size in s
    figure: matplotlib figure
        The figure to be used. If not given a figure is created
    ax: matplotlib axis
        The axis to be used. If not given an axis is created

    Returns
    -------
    W: ISIDistributionWidget
        The output widget
    """
    if sampling_frequency is None:
        if sorting.get_sampling_frequency() is None:
            raise Exception("Sampling rate information is not in the SortingExtractor. "
                            "Provide the 'sampling_frequency' argument")
        else:
            sampling_frequency = sorting.get_sampling_frequency()
    W = ISIDistributionWidget(
        sorting=sorting,
        sampling_frequency=sampling_frequency,
        unit_ids=unit_ids,
        bins=bins,
        window=window,
        figure=figure,
        ax=ax
    )
    W.plot()
    return W


class ISIDistributionWidget(BaseMultiWidget):
    def __init__(self, *, sorting, sampling_frequency, unit_ids=None, bins=10, window=1, figure=None, ax=None):
        BaseMultiWidget.__init__(self, figure, ax)
        self._sorting = sorting
        self._unit_ids = unit_ids
        self._sampling_frequency = sampling_frequency
        self._bins = bins
        self._maxw = window
        self.name = 'ISIDistribution'

    def plot(self):
        self._do_plot()

    def _do_plot(self):
        units = self._unit_ids
        if units is None:
            units = self._sorting.get_unit_ids()
        list_isi = []
        for unit in units:
            times = self._sorting.get_unit_spike_train(unit_id=unit) / float(self._sampling_frequency)
            bin_counts, bin_edges = compute_isi_dist(times, bins=self._bins, maxwindow=self._maxw)
            item = dict(
                title='Unit {}'.format(int(unit)),
                bin_counts=bin_counts,
                bin_edges=bin_edges,
                xticks=[0, self._maxw / 2, self._maxw]
            )
            list_isi.append(item)
        with plt.rc_context({'axes.edgecolor': 'gray'}):
            self._plot_isi_multi(list_isi)

    def _plot_isi_multi(self, list_isi, *, ncols=5, **kwargs):
        if len(list_isi) < ncols:
            ncols = len(list_isi)
        nrows = np.ceil(len(list_isi) / ncols)
        for i, item in enumerate(list_isi):
            ax = self.get_tiled_ax(i, nrows, ncols, hspace=0.7)
            _plot_isi(**item, **kwargs, ax=ax)


def _plot_isi( *, bin_counts, bin_edges, ax, xticks=None, title=''):
    bins = bin_edges[:-1] + np.mean(np.diff(bin_edges))
    wid = np.mean(np.diff(bins))  #(bin_edges[1] - bin_edges[0]) * len(bin_edges)
    ax.bar(x=bins, height=bin_counts, width=wid, color='gray', align='edge')
    if xticks is not None:
        ax.set_xticks(xticks)
    ax.set_xlabel('dt (s)')
    ax.set_yticks([])
    if title:
        ax.set_title(title, color='gray')


def compute_isi_dist(times, *, bins, maxwindow=10):
    isi = np.diff(times)
    isi = isi[isi < maxwindow]
    # bins = np.arange(0, maxisi, binsize)
    bin_counts, bin_edges = np.histogram(isi, bins=bins, density=True)
    return bin_counts, bin_edges
