import numpy as np
from matplotlib import pyplot as plt


def plot_crosscorrelograms(sorting, sample_rate=None, unit_ids=None, bin_size=2, window=50):
    if sample_rate is None:
        if sorting.get_sampling_frequency() is None:
            raise Exception("Sampling rate information is not in the SortingExtractor. "
                            "Provide the 'sample_rate' argument")
        else:
            sample_rate = sorting.get_sampling_frequency()
    W = CrossCorrelogramsWidget(
        sorting=sorting,
        samplerate=sample_rate,
        unit_ids=unit_ids,
        binsize=bin_size,
        window=window
    )
    W.plot()


class CrossCorrelogramsWidget:
    def __init__(self, *, sorting, samplerate, unit_ids=None, binsize=2, window=50):
        self._SX = sorting
        self._unit_ids = unit_ids
        self._figure = None
        self._samplerate = samplerate
        self._binsize = binsize
        self._window = window

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
            times = self._SX.get_unit_spike_train(unit_id=unit)
            max_dt_msec = self._window
            bin_size_msec = self._binsize
            max_dt_tp = max_dt_msec * self._samplerate / 1000
            bin_size_tp = bin_size_msec * self._samplerate / 1000
            (bin_counts, bin_edges) = compute_autocorrelogram(times, max_dt_tp=max_dt_tp, bin_size_tp=bin_size_tp)
            item = dict(
                title=str(unit),
                bin_counts=bin_counts,
                bin_edges=bin_edges
            )
            list.append(item)
        with plt.rc_context({'axes.edgecolor': 'gray'}):
            self._plot_correlograms_multi(list)

    def _plot_correlogram(self, *, bin_counts, bin_edges, title=''):
        wid = (bin_edges[1] - bin_edges[0]) * 1000
        plt.bar(x=bin_edges[0:-1] * 1000, height=bin_counts, width=wid, color='gray', align='edge')
        plt.xlabel('dt (msec)')
        plt.gca().get_yaxis().set_ticks([])
        plt.gca().get_xaxis().set_ticks([])
        plt.gca().get_yaxis().set_ticks([])
        if title:
            plt.title(title, color='gray')

    def _plot_correlograms_multi(self, list, *, ncols=5, **kwargs):
        nrows = np.ceil(len(list) / ncols)
        plt.figure(figsize=(3 * ncols + 0.1, 3 * nrows + 0.1))
        for i, item in enumerate(list):
            plt.subplot(nrows, ncols, i + 1)
            self._plot_correlogram(**item, **kwargs)


def compute_autocorrelogram(times, *, max_dt_tp, bin_size_tp, max_samples=None):
    num_bins_left = int(max_dt_tp / bin_size_tp)  # number of bins to the left of the origin
    L = len(times)  # number of events
    times2 = np.sort(times)  # the sorted times
    step = 1  # This is the index step between an event and the next one to compare
    candidate_inds = np.arange(L)  # These are the events we are going to consider
    if max_samples is not None:
        if len(candidate_inds) > max_samples:
            candidate_inds = np.random.choice(candidate_inds, size=max_samples, replace=False)
    vals_list = []  # A list of all offsets we have accumulated
    while True:
        candidate_inds = candidate_inds[
            candidate_inds + step < L]  # we only consider events that are within workable range
        candidate_inds = candidate_inds[times2[candidate_inds + step] - times2[
            candidate_inds] <= max_dt_tp]  # we only consider event-pairs that are within max_dt_tp apart
        if len(candidate_inds) > 0:  # if we have some events to consider
            vals = times2[candidate_inds + step] - times2[candidate_inds]
            vals_list.append(vals)  # add to the autocorrelogram
            vals_list.append(-vals)  # keep it symmetric
        else:
            break  # no more to consider
        step += 1
    if len(vals_list) > 0:  # concatenate all the values
        all_vals = np.concatenate(vals_list)
    else:
        all_vals = np.array([]);
    aa = np.arange(-num_bins_left, num_bins_left + 1) * bin_size_tp
    all_vals = np.sign(all_vals) * (np.abs(
        all_vals) - bin_size_tp * 0.00001)  # a trick to make the histogram symmetric due to differences in rounding for positive and negative, i suppose
    bin_counts, bin_edges = np.histogram(all_vals, bins=aa)
    return (bin_counts, bin_edges)
