import numpy as np
from matplotlib import pyplot as plt
from spikewidgets.widgets.basewidget import BaseMultiWidget


def plot_isi_distribution(sorting, sample_rate=None, unit_ids=None, bin_size=2, max_isi=100, figure=None, ax=None):
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
        max_isi=max_isi,
        figure=figure,
        ax=ax
    )
    W.plot()
    return W


class ISIDistributionWidget(BaseMultiWidget):
    def __init__(self, *, sorting, samplerate, unit_ids=None, bin_size=2, max_isi=100, figure=None, ax=None):
        BaseMultiWidget.__init__(self, figure, ax)
        self._sorting = sorting
        self._unit_ids = unit_ids
        self._samplerate = samplerate
        self._binsize = bin_size
        self._maxisi = max_isi

    def plot(self):
        self._do_plot()

    def _do_plot(self):
        units = self._unit_ids
        if units is None:
            units = self._sorting.get_unit_ids()
        list_isi = []
        for unit in units:
            times = self._sorting.get_unit_spike_train(unit_id=unit) / float(self._samplerate)
            bin_counts, bin_edges = compute_isi_dist(times, binsize=self._binsize, maxisi=self._maxisi)
            item = dict(
                title=str(unit),
                bin_counts=bin_counts,
                bin_edges=bin_edges
            )
            list_isi.append(item)
        with plt.rc_context({'axes.edgecolor': 'gray'}):
            self._plot_isi_multi(list_isi)

    def _plot_isi_multi(self, list_isi, *, ncols=5, **kwargs):
        if len(list_isi) < ncols:
            ncols = len(list_isi)
        nrows = np.ceil(len(list_isi) / ncols)
        for i, item in enumerate(list_isi):
            ax = self.get_tiled_ax(i, nrows, ncols)
            _plot_isi(**item, **kwargs, ax=ax)


def _plot_isi( *, bin_counts, bin_edges, ax, title=''):
    wid = (bin_edges[1] - bin_edges[0]) * 1000
    ax.bar(x=bin_edges[0:-1] * 1000, height=bin_counts, width=wid, color='gray', align='edge')
    ax.set_xlabel('dt (sec)')
    ax.set_xticks([])
    ax.set_yticks([])
    if title:
        ax.set_title(title, color='gray')


def compute_isi_dist(times, *, binsize, maxisi=100):
    isi = np.diff(times)
    bins = np.arange(0, maxisi, binsize)
    bin_counts, bin_edges = np.histogram(isi, bins=bins)
    return bin_counts, bin_edges
