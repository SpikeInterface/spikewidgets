import numpy as np
from spikewidgets.widgets.basewidget import BaseWidget


def plot_spectrum(recording, channels=None, trange=None, freqrange=None, color_groups=False, color='steelblue',
                  nfft=256, figure=None, ax=None):
    """
    Plots electrode geometry.

    Parameters
    ----------
    recording: RecordingExtractor
        The recordng extractor object
    channels: list
        The channels to show
    trange: list
        List with start time and end time
    freqrange: list
        List with start frequency and end frequency
    color_groups: bool
        If True groups are plotted with different colors
    color: matplotlib color
        The color to be used
    figure: matplotlib figure
        The figure to be used. If not given a figure is created
    ax: matplotlib axis
        The axis to be used. If not given an axis is created

    Returns
    -------
    W: TimeseriesWidget
        The output widget
    """
    W = SpectrumWidget(
        recording=recording,
        channels=channels,
        trange=trange,
        freqrange=freqrange,
        color_groups=color_groups,
        color=color,
        nfft=nfft,
        figure=figure,
        ax=ax
    )
    W.plot()
    return W


def plot_spectrogram(recording, channel, trange=None, freqrange=None, cmap='viridis',
                     nfft=256, figure=None, ax=None):
    """
    Plots electrode geometry.

    Parameters
    ----------
    recording: RecordingExtractor
        The recordng extractor object
    channel: int
        The channel to plot spectrogram of
    trange: list
        List with start time and end time
    freqrange: list
        List with start frequency and end frequency
    cmap: matplotlib colorma
        The colormap to be used
    figure: matplotlib figure
        The figure to be used. If not given a figure is created
    ax: matplotlib axis
        The axis to be used. If not given an axis is created

    Returns
    -------
    W: TimeseriesWidget
        The output widget
    """
    W = SpectrogramWidget(
        recording=recording,
        channel=channel,
        trange=trange,
        freqrange=freqrange,
        cmap=cmap,
        nfft=nfft,
        figure=figure,
        ax=ax
    )
    W.plot()
    return W


class SpectrumWidget(BaseWidget):
    def __init__(self, *, recording, channels=None, trange=None, freqrange=None,
                 color_groups=False, color=None, nfft=256, figure=None, ax=None):
        BaseWidget.__init__(self, figure, ax)
        self._recording = recording
        self._samplerate = recording.get_sampling_frequency()
        self._ax = ax
        self._color_groups = color_groups
        self._trange = trange
        self._frange = freqrange
        self._color = color
        self._nfft = nfft

        if channels is None:
            self._channels = self._recording.get_channel_ids()
        else:
            self._channels = channels

        if color_groups:
            self._colors = []
            self._group_color_map = {}
            all_groups = recording.get_channel_groups()
            groups = np.unique(all_groups)
            N = len(groups)
            import colorsys
            HSV_tuples = [(x * 1.0 / N, 0.5, 0.5) for x in range(N)]
            self._colors = list(map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples))
            color_idx = 0
            for group in groups:
                self._group_color_map[group] = color_idx
                color_idx += 1
        self.name = 'Spectrum'

    def plot(self):
        if self._trange is not None:
            start_frame = int(self._trange[0] * self._recording.get_sampling_frequency())
            end_frame = int(self._trange[1] * self._recording.get_sampling_frequency())
        else:
            start_frame = 0
            end_frame = self._recording.get_num_frames()

        traces = self._recording.get_traces(channel_ids=self._channels, start_frame=start_frame, end_frame=end_frame)

        for (t, ch) in zip(traces, self._channels):
            if self._color_groups:
                group = self._recording.get_channel_groups(channel_ids=[ch])[0]
                group_color_idx = self._group_color_map[group]
                self.ax.psd(t, Fs=self._samplerate, color=self._colors[group_color_idx], NFFT=self._nfft)
            else:
                self.ax.psd(t, Fs=self._samplerate, NFFT=self._nfft)

        if self._frange is not None:
            self.ax.set_xlim(self._frange)


class SpectrogramWidget(BaseWidget):
    def __init__(self, *, recording, channel, trange=None, freqrange=None,
                 cmap=None, nfft=256, figure=None, ax=None):
        BaseWidget.__init__(self, figure, ax)
        self._recording = recording
        self._samplerate = recording.get_sampling_frequency()
        self._ax = ax
        self._trange = trange
        self._frange = freqrange
        self._cmap = cmap
        self._nfft = nfft

        assert channel is not None, "Specify which 'channel' to use"
        self._channel = channel
        self.name = 'Spectrocram'

    def plot(self):
        if self._trange is not None:
            start_frame = int(self._trange[0] * self._recording.get_sampling_frequency())
            end_frame = int(self._trange[1] * self._recording.get_sampling_frequency())
        else:
            start_frame = 0
            end_frame = self._recording.get_num_frames()

        traces = self._recording.get_traces(channel_ids=[self._channel], start_frame=start_frame, end_frame=end_frame)

        self.ax.specgram(traces[0], Fs=self._recording.get_sampling_frequency(), scale='dB', cmap=self._cmap,
                         NFFT=self._nfft)

        if self._frange is not None:
            self.ax.set_ylim(self._frange)

        if self._trange is not None:
            l_xticks = len(self.ax.get_xticks())
            self.ax.set_xticklabels(np.linspace(self._trange[0], self._trange[1], l_xticks))