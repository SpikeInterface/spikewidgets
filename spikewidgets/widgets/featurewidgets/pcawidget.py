import numpy as np
from matplotlib import pyplot as plt
from spikewidgets.widgets.basewidget import BaseMultiWidget
import spiketoolkit as st


def plot_pca_features(recording, sorting, unit_ids=None, max_spikes_per_unit=100, nproj=4, colormap=None,
                      figure=None, ax=None):
    """
    Plots unit PCA features on best projections.

    Parameters
    ----------
    recording: RecordingExtractor
        The recordng extractor object.
    sorting: SortingExtractor
        The sorting extractor object.
    unit_ids: list
        List of unit ids.
    max_spikes_per_unit: int
        Maximum number of spikes to display per unit.
    nproj: int
        Number of best projections to display.
    colormap: matplotlib colormap
        The colormap to be used. If not given default is used.
    figure: matplotlib figure
        The figure to be used. If not given a figure is created.
    ax: matplotlib axis
        The axis to be used. If not given an axis is created.

    Returns
    -------
    W: PCAWidget
        The output widget
    """
    W = PCAWidget(
        sorting=sorting,
        recording=recording,
        unit_ids=unit_ids,
        max_spikes_per_unit=max_spikes_per_unit,
        nproj=nproj,
        colormap=colormap,
        figure=figure,
        ax=ax
    )
    W.plot()
    return W


class PCAWidget(BaseMultiWidget):
    def __init__(self, *, recording, sorting, unit_ids=None, max_spikes_per_unit=100, nproj=4, colormap=None,
                 figure=None, ax=None, save_as_features=False, save_waveforms_as_features=False):
        BaseMultiWidget.__init__(self, figure, ax)
        self._sorting = sorting
        self._recording = recording
        self._unit_ids = unit_ids
        self._nproj = nproj
        self._max_spikes_per_unit = max_spikes_per_unit
        self._pca_scores = None
        self._colormap = colormap
        self._save_as_features = save_as_features
        self._save_waveforms_as_features = save_waveforms_as_features
        self.name = 'Feature'

    def _compute_pca(self):
        self._pca_scores = st.postprocessing.compute_unit_pca_scores(recording=self._recording,
                                                                     sorting=self._sorting,
                                                                     by_electrode=True,
                                                                     max_spikes_per_unit=self._max_spikes_per_unit,
                                                                     save_as_features=self._save_as_features,
                                                                     save_waveforms_as_features=self._save_waveforms_as_features)

    def plot(self):
        self._do_plot()

    def _do_plot(self):
        units = self._unit_ids
        if units is None:
            units = self._sorting.get_unit_ids()
        self._units = units

        if self._pca_scores is None:
            self._compute_pca()

        # find projections with best separation
        n_pc = self._pca_scores[0].shape[2]
        n_ch = self._pca_scores[0].shape[1]

        distances = []
        proj = []
        for ch1 in range(n_ch):
            for pc1 in range(n_pc):
                for ch2 in range(n_ch):
                    for pc2 in range(n_pc):
                        if ch1 != ch2 or pc1 != pc2:
                            dist = self.compute_cluster_average_distance(
                                pc1, ch1, pc2, ch2)
                            if [ch1, pc1, ch2, pc2] not in proj and [
                                    ch2, pc2, ch1, pc1] not in proj:
                                distances.append(dist)
                                proj.append([ch1, pc1, ch2, pc2])

        list_best_proj = np.array(
            proj)[np.argsort(distances)[::-1][:self._nproj]]
        self._plot_proj_multi(list_best_proj)

    def compute_cluster_average_distance(self, pc1, ch1, pc2, ch2):
        centroids = np.zeros((len(self._pca_scores), 2))
        for i, pcs in enumerate(self._pca_scores):
            centroids[i, 0] = np.median(pcs[:, ch1, pc1], axis=0)
            centroids[i, 1] = np.median(pcs[:, ch2, pc2], axis=0)

        dist = []
        for i, c1 in enumerate(centroids):
            for j, c2 in enumerate(centroids):
                if i > j:
                    dist.append(np.linalg.norm(c2 - c1))

        return np.mean(dist)

    def _plot_proj_multi(self, best_proj, ncols=5):
        if len(best_proj) < ncols:
            ncols = len(best_proj)
        nrows = np.ceil(len(best_proj) / ncols)
        for i, bp in enumerate(best_proj):
            ax = self.get_tiled_ax(i, nrows, ncols, hspace=0.3)
            self._plot_proj(proj=bp, ax=ax)

    def _plot_proj(self, *, proj, ax, title=''):
        ch1, pc1, ch2, pc2 = proj
        if self._colormap is not None:
            cm = plt.get_cmap(self._colormap)
            colors = [cm(i / len(self._pca_scores))
                      for i in np.arange(len(self._pca_scores))]
        else:
            colors = [None] * len(self._pca_scores)
        for i, pc in enumerate(self._pca_scores):
            if self._sorting.get_unit_ids()[i] in self._units:
                ax.plot(pc[:, ch1, pc1], pc[:, ch2, pc2],
                        '*', color=colors[i], alpha=0.3)
        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_xlabel('ch {} - pc {}'.format(ch1, pc1))
        ax.set_ylabel('ch {} - pc {}'.format(ch2, pc2))
        if title:
            ax.set_title(title, color='gray')
