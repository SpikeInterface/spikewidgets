import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from spikewidgets.widgets.basewidget import BaseWidget


def plot_timeseries(recording, channel_ids=None, trange=None, color_groups=False,
                    figure=None, ax=None):
    """
    Plots electrode geometry.

    Parameters
    ----------
    recording: RecordingExtractor
        The recordng extractor object
    channel_ids: list
        The channel ids to display.
    trange: list
        List with start time and end time
    color_groups: bool
        If True groups are plotted with different colors
    figure: matplotlib figure
        The figure to be used. If not given a figure is created
    ax: matplotlib axis
        The axis to be used. If not given an axis is created

    Returns
    -------
    W: TimeseriesWidget
        The output widget
    """
    W = TimeseriesWidget(
        recording=recording,
        channel_ids=channel_ids,
        trange=trange,
        color_groups=color_groups,
        figure=figure,
        ax=ax
    )
    W.plot()
    return W


class TimeseriesWidget(BaseWidget):
    def __init__(self, *, recording, channel_ids=None, trange=None,
                 color_groups=False, figure=None,  ax=None):
        BaseWidget.__init__(self, figure, ax)
        self._recording = recording
        self._sampling_frequency = recording.get_sampling_frequency()
        self._visible_channels = channel_ids
        if self._visible_channels is None:
            self._visible_channels = recording.get_channel_ids()
        self._visible_trange = trange
        if self._visible_trange is None:
            self._visible_trange = [0, np.minimum(10000, recording.get_num_frames())]
        else:
            assert len(trange) == 2, "'trange' should be a list with start and end time in seconds"
            self._visible_trange = [int(t * recording.get_sampling_frequency()) for t in trange]
        self._initialize_stats()
        # self._vspacing = self._mean_channel_std * 20
        self._vspacing = self._max_channel_amp * 1.5
        self._visible_trange = self._fix_trange(self._visible_trange)
        self._ax = ax
        self._color_groups = color_groups
        if color_groups:
            self._colors = []
            self._group_color_map = {}
            all_groups = recording.get_channel_groups()
            groups = np.unique(all_groups)
            N = len(groups)
            import colorsys
            HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
            self._colors = list(map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples))
            color_idx = 0
            for group in groups:
                self._group_color_map[group] = color_idx
                color_idx += 1
        self.name = 'TimeSeries'
            
    def plot(self):
        self._do_plot()

    def _do_plot(self):
        chunk0 = self._recording.get_traces(
            channel_ids=self._visible_channels,
            start_frame=self._visible_trange[0],
            end_frame=self._visible_trange[1]
        )

        self.ax.set_xlim(self._visible_trange[0] / self._sampling_frequency, self._visible_trange[1] / self._sampling_frequency)
        self.ax.set_ylim(-self._vspacing, self._vspacing * len(self._visible_channels))
        self.ax.get_xaxis().set_major_locator(MaxNLocator(prune='both'))
        self.ax.get_yaxis().set_ticks([])
        self.ax.set_xlabel('time (s)')

        self._plots = {}
        self._plot_offsets = {}
        offset0 = self._vspacing * (len(self._visible_channels) - 1)
        tt = np.arange(self._visible_trange[0], self._visible_trange[1]) / self._sampling_frequency
        for im, m in enumerate(self._visible_channels):
            self._plot_offsets[m] = offset0
            if self._color_groups:
                group = self._recording.get_channel_groups(channel_ids=[m])[0]
                group_color_idx = self._group_color_map[group]
                self._plots[m] = self.ax.plot(tt, self._plot_offsets[m] + chunk0[im, :],
                                              color=self._colors[group_color_idx])
            else:
                self._plots[m] = self.ax.plot(tt, self._plot_offsets[m] + chunk0[im, :])
            offset0 = offset0 - self._vspacing

    def _fix_trange(self, trange):
        N = self._recording.get_num_frames()
        if trange[1] > N:
            # trange[0] += N - trange[1]
            # trange[1] = N - trange[1]
            trange[1] = N
        if trange[0] < 0:
            # trange[1] += -trange[0]
            trange[0] = 0
        # trange[0] = np.maximum(0, trange[0])
        # trange[1] = np.minimum(N, trange[1])
        return trange

    def _initialize_stats(self):
        self._channel_stats = {}
        # M=self._reader.numChannels()
        # N=self._reader.numTimepoints()
        chunk0 = self._recording.get_traces(
            channel_ids=self._visible_channels,
            start_frame=self._visible_trange[0],
            end_frame=self._visible_trange[1]
        )
        # chunk0=self._reader.getChunk(channels=self._visible_channels,trange=self._visible_trange)
        M0 = chunk0.shape[0]
        for ii in range(M0):
            self._channel_stats[self._visible_channels[ii]] = _compute_channel_stats_from_data(chunk0[ii, :])
        self._mean_channel_std = np.mean([self._channel_stats[m]['std'] for m in self._visible_channels])
        self._max_channel_amp = np.max([self._channel_stats[m]['max_amp'] for m in self._visible_channels])


def _compute_channel_stats_from_data(X):
    return dict(
        mean=np.mean(X),
        std=np.std(X),
        max_amp=np.max(np.abs(X))
    )
