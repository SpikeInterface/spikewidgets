import numpy as np
from matplotlib import pyplot as plt
from spikewidgets.widgets.basewidget import BaseMultiWidget


def plot_amplitudes_distribution(recording, sorting, unit_ids=None, max_spikes_per_unit=100,
                                 figure=None, ax=None):
    """
    Plots waveform amplitudes distribution.

    Parameters
    ----------
    recording: RecordingExtractor
        The recording extractor object
    sorting: SortingExtractor
        The sorting extractor object
    unit_ids: list
        List of unit ids
    max_spikes_per_unit: int
        Maximum number of spikes to display per unit.
    figure: matplotlib figure
        The figure to be used. If not given a figure is created
    ax: matplotlib axis
        The axis to be used. If not given an axis is created

    Returns
    -------
    W: AmplitudeDistributionWidget
        The output widget
    """
    W = AmplitudeDistributionWidget(
        sorting=sorting,
        recording=recording,
        unit_ids=unit_ids,
        max_spikes_per_unit=max_spikes_per_unit,
        figure=figure,
        ax=ax
    )
    W.plot()
    return W


def plot_amplitudes_timeseries(recording, sorting, unit_ids=None, max_spikes_per_unit=100,
                               figure=None, ax=None):
    """
    Plots waveform amplitudes timeseries.

    Parameters
    ----------
    recording: RecordingExtractor
        The recording extractor object
    sorting: SortingExtractor
        The sorting extractor object
    unit_ids: list
        List of unit ids
    max_spikes_per_unit: int
        Maximum number of spikes to display per unit.
    figure: matplotlib figure
        The figure to be used. If not given a figure is created
    ax: matplotlib axis
        The axis to be used. If not given an axis is created

    Returns
    -------
    W: AmplitudeTimeseriesWidget
        The output widget
    """
    W = AmplitudeTimeseriesWidget(
        sorting=sorting,
        recording=recording,
        unit_ids=unit_ids,
        max_spikes_per_unit=max_spikes_per_unit,
        figure=figure,
        ax=ax
    )
    W.plot()
    return W


class AmplitudeBaseWidget(BaseMultiWidget):
    def __init__(self, recording, sorting, max_spikes_per_unit=100, figure=None, ax=None):
        BaseMultiWidget.__init__(self, figure, ax)
        self._sorting = sorting
        self._recording = recording
        self._max_spikes_per_unit = max_spikes_per_unit

    def plot(self):
        self._do_plot()

    def _do_plot(self):
        raise NotImplementedError

    def compute_amps(self, *, times):
        amps = np.zeros(len(times))
        times_ind = (times * self._recording.get_sampling_frequency()).astype(int)
        for i, t in enumerate(times_ind):
            amp = self._recording.get_traces(start_frame=t-1, end_frame=t+1)
            amps[i] = np.min(amp)
        return amps


class AmplitudeTimeseriesWidget(AmplitudeBaseWidget):
    def __init__(self, *, recording, sorting, unit_ids=None, max_spikes_per_unit=100, figure=None, ax=None):
        AmplitudeBaseWidget.__init__(self, recording, sorting, max_spikes_per_unit, figure, ax)
        self._unit_ids = unit_ids
        self.name = 'AmplitudeTimeseries'

    def _do_plot(self):
        units = self._unit_ids
        if units is None:
            units = self._sorting.get_unit_ids()
        list_amps = []
        all_amps = np.array([])
        for unit in units:
            times = self._sorting.get_unit_spike_train(unit_id=unit) / \
                    float(self._recording.get_sampling_frequency())
            if len(times) > self._max_spikes_per_unit:
                times = times[np.random.permutation(len(times))[:self._max_spikes_per_unit]]
            amps = self.compute_amps(times=times)
            item = dict(
                title='Unit {}'.format(int(unit)),
                times=times,
                amps=amps
            )
            list_amps.append(item)
            all_amps = np.concatenate((all_amps, amps))
        with plt.rc_context({'axes.edgecolor': 'gray'}):
            ylim = [np.min(all_amps) - 10, np.max(all_amps) + 10]
            self._plot_amp_time_multi(list_amps, ylim=ylim)

    def _plot_amp_time_multi(self, list_isi, *, ncols=5, ylim=None, **kwargs):
        if len(list_isi) < ncols:
            ncols = len(list_isi)
        nrows = np.ceil(len(list_isi) / ncols)
        for i, item in enumerate(list_isi):
            ax = self.get_tiled_ax(i, nrows, ncols, hspace=0.7, wspace=1)
            if i == 0 or i == ncols:
                ylab = True
            else:
                ylab = False
            _plot_amp_time(**item, **kwargs, ax=ax, ylim=ylim, ylab=ylab)


class AmplitudeDistributionWidget(AmplitudeBaseWidget):
    def __init__(self, *, recording, sorting, unit_ids=None, max_spikes_per_unit=100, figure=None, ax=None):
        AmplitudeBaseWidget.__init__(self, recording, sorting, max_spikes_per_unit, figure, ax)
        self._unit_ids = unit_ids
        self.name = 'AmplitudeDistribution'

    def _do_plot(self):
        units = self._unit_ids
        if units is None:
            units = self._sorting.get_unit_ids()
        list_amps = []
        all_amps = np.array([])
        for unit in units:
            times = self._sorting.get_unit_spike_train(unit_id=unit) / \
                    float(self._recording.get_sampling_frequency())
            if len(times) > self._max_spikes_per_unit:
                times = times[np.random.permutation(len(times))[:self._max_spikes_per_unit]]
            amps = self.compute_amps(times=times)
            item = dict(
                title='Unit {}'.format(int(unit)),
                times=times,
                amps=amps
            )
            list_amps.append(item)
            all_amps = np.concatenate((all_amps, amps))
        with plt.rc_context({'axes.edgecolor': 'gray'}):
            self._plot_amp_hist_multi(list_amps)

    def _plot_amp_hist_multi(self, list_isi, *, ncols=5, **kwargs):
        if len(list_isi) < ncols:
            ncols = len(list_isi)
        nrows = np.ceil(len(list_isi) / ncols)
        for i, item in enumerate(list_isi):
            ax = self.get_tiled_ax(i, nrows, ncols, hspace=0.7)
            _plot_amp_dist(**item, **kwargs, ax=ax)


def roundup(x, num=5):
    return int(np.ceil(x / num)) * num


def rounddown(x, num=5):
    return int(np.floor(x / num)) * num


def _plot_amp_dist(*, times, amps, ax, title=''):
    std_amps = np.std(amps)
    idxs = (amps > np.median(amps) - 4 * std_amps) & (amps < np.median(amps) + 4 * std_amps)
    amps = amps[idxs]
    h, b, _ = ax.hist(amps, density=True, color='gray', bins=10)
    ax.set_yticks([])
    ax.set_xticks([rounddown(np.mean(amps) - 2 * np.std(amps)),
                   roundup(np.mean(amps) + 2 * np.std(amps))])
    ax.set_xlabel('amp')
    if title:
        ax.set_title(title, color='gray')


def _plot_amp_time(*, times, amps, ax, title='', ylim=None, ylab=False):
    std_amps = np.std(amps)
    idxs = (amps > np.median(amps) - 4 * std_amps) & (amps < np.median(amps) + 4 * std_amps)
    amps = amps[idxs]
    times = times[idxs]
    ax.plot(times, amps, '*', color='gray', alpha=0.3)
    if ylab:
        ax.set_ylabel('amp')
    ax.set_xlabel('time (s)')
    ax.set_xticks([])
    if ylim is not None:
        ax.set_ylim(ylim)
    if title:
        ax.set_title(title, color='gray')
