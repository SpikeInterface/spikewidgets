import numpy as np
from matplotlib import pyplot as plt


def plot_isi_distribution(sorting, sample_rate=None, unit_ids=None, bin_size=2, max_isi=100):
    if sample_rate is None:
        if sorting.get_sampling_frequency() is None:
            raise Exception("Sampling rate information is not in the SortingExtractor. "
                            "Provide the 'sample_rate' argument")
        else:
            sample_rate = sorting.get_sampling_frequency()
    W = ISIDistributionWidget(
        sorting=sorting,
        samplerate=sample_rate,
        unit_ids=unit_ids,
        bin_size=bin_size,
        max_isi=max_isi
    )
    W.plot()


class ISIDistributionWidget:
    def __init__(self, *, sorting, samplerate, unit_ids=None, bin_size=2, max_isi=100):
        self._SX = sorting
        self._unit_ids = unit_ids
        self._figure = None
        self._samplerate = samplerate
        self._binsize = bin_size
        self._maxisi = max_isi

    def plot(self):
        self._do_plot()

    def figure(self):
        return self._figure

    def _do_plot(self):
        units = self._unit_ids
        if units is None:
            units = self._SX.get_unit_ids()
        list = []
        for unit in units:
            times = self._SX.get_unit_spike_train(unit_id=unit) / float(self._samplerate) * 1000
            (bin_counts, bin_edges) = compute_isi_dist(times, binsize=self._binsize, maxisi=self._maxisi)
            item = dict(
                title=str(unit),
                bin_counts=bin_counts,
                bin_edges=bin_edges
            )
            list.append(item)
        with plt.rc_context({'axes.edgecolor': 'gray'}):
            self._plot_isi_multi(list)

    def _plot_isi(self, *, bin_counts, bin_edges, title=''):
        wid = (bin_edges[1] - bin_edges[0]) * 1000
        plt.bar(x=bin_edges[0:-1] * 1000, height=bin_counts, width=wid, color='gray', align='edge')
        plt.xlabel('dt (msec)')
        plt.gca().get_yaxis().set_ticks([])
        plt.gca().get_xaxis().set_ticks([])
        plt.gca().get_yaxis().set_ticks([])
        if title:
            plt.title(title, color='gray')

    def _plot_isi_multi(self, list, *, ncols=5, **kwargs):
        nrows = np.ceil(len(list) / ncols)
        plt.figure(figsize=(3 * ncols + 0.1, 3 * nrows + 0.1))
        for i, item in enumerate(list):
            plt.subplot(nrows, ncols, i + 1)
            self._plot_isi(**item, **kwargs)


def compute_isi_dist(times, *, binsize, maxisi=100):
    isi = np.diff(times)
    bins = np.arange(0, maxisi, binsize)
    bin_counts, bin_edges = np.histogram(isi, bins=bins)
    return bin_counts, bin_edges
