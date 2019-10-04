import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from spikewidgets.widgets.basewidget import BaseWidget


def plot_multicomp_graph(multi_sorting_comparison, sorter_names=None, draw_labels=False, node_cmap='viridis',
                         edge_cmap='hot_r', alpha_edges=0.7, colorbar=False, figure=None, ax=None):
    """
    Plots multi sorting comparison graph.

    Parameters
    ----------
    multi_sorting_comparison: MultiSortingComparison
        The multi sorting comparison object
    sorter_names: list
        The names of the sorters
    draw_labels: bool
        If True unit labels are shown
    node_cmap: matplotlib colormap
        The colormap to be used for the nodes (default 'viridis')
    edge_cmap: matplotlib colormap
        The colormap to be used for the edges (default 'hot')
    alpha_edges: float
        Alpha value for edges
    colorbar: bool
        If True a colorbar for the edges is plotted
    figure: matplotlib figure
        The figure to be used. If not given a figure is created
    ax: matplotlib axis
        The axis to be used. If not given an axis is created

    Returns
    -------
    W: MultiCompGraphWidget
        The output widget
    """
    try:
        import networkx as nx
    except ImportError as e:
        raise ImportError('Install networkx to use the multi comparison widget.')
    W = MultiCompGraphWidget(
        multi_sorting_comparison=multi_sorting_comparison,
        sorternames=sorter_names,
        node_cmap=node_cmap,
        edge_cmap=edge_cmap,
        drawlabels=draw_labels,
        alpha_edges=alpha_edges,
        colorbar=colorbar,
        figure=figure,
        ax=ax
    )
    W.plot()
    return W


class MultiCompGraphWidget(BaseWidget):
    def __init__(self, *, multi_sorting_comparison, sorternames=None, drawlabels=False, node_cmap='viridis',
                 edge_cmap='hot', alpha_edges=0.5, colorbar=False, figure=None, ax=None):
        BaseWidget.__init__(self, figure, ax)
        self._msc = multi_sorting_comparison
        self._sorter_names = sorternames
        self._drawlabels = drawlabels
        self._node_cmap = node_cmap
        self._edge_cmap = edge_cmap
        self._colorbar = colorbar
        self._alpha_edges = alpha_edges
        self.name = 'MultiCompGraph'

    def plot(self):
        self._do_plot()

    def _do_plot(self):
        import networkx as nx

        g = self._msc.graph
        edge_col = []
        for e in g.edges(data=True):
            n1, n2, d = e
            edge_col.append(d['weight'])
        nodes_col = np.array([])
        for i, sort in enumerate(self._msc.get_sorting_list()):
            nodes_col = np.concatenate((nodes_col, np.array([i] * len(sort.get_unit_ids()))))
        nodes_col = nodes_col / len(self._msc.get_sorting_list())

        _ = plt.set_cmap(self._node_cmap)
        _ = nx.draw_networkx_nodes(g, pos=nx.circular_layout(sorted(g)), nodelist=sorted(g.nodes),
                                   node_color=nodes_col, node_size=20, ax=self.ax)
        _ = nx.draw_networkx_edges(g, pos=nx.circular_layout((sorted(g))), nodelist=sorted(g.nodes),
                                   edge_color=edge_col, alpha=self._alpha_edges,
                                   edge_cmap=plt.cm.get_cmap(self._edge_cmap), edge_vmin=self._msc.match_score,
                                   edge_vmax=1, ax=self.ax)
        if self._drawlabels:
            _ = nx.draw_networkx_labels(g, pos=nx.circular_layout((sorted(g))), nodelist=sorted(g.nodes), ax=self.ax)
        if self._colorbar:
            norm = matplotlib.colors.Normalize(vmin=self._msc.min_accuracy, vmax=1)
            cmap = plt.cm.get_cmap(self._edge_cmap)
            m = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
            self.figure.colorbar(m)

        self.ax.axis('off')
