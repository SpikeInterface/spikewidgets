import numpy as np
from matplotlib import pyplot as plt
from spikewidgets.widgets.basewidget import BaseWidget


def plot_multicomp_graph(multisortingcomparison, sorter_names=None, draw_labels=False, node_cmap='viridis',
                         edge_cmap='hot', title='', figure=None, ax=None):
    try:
        import networkx as nx
    except ImportError as e:
        raise ImportError('Install networkx to use the multi comparison widget.')
    W = MultiCompGraphWidget(
        multisortingcomparison=multisortingcomparison,
        sorternames=sorter_names,
        title=title,
        node_cmap=node_cmap,
        edge_cmap=edge_cmap,
        drawlabels=draw_labels,
        figure=figure,
        ax=ax
    )
    W.plot()
    return W


class MultiCompGraphWidget(BaseWidget):
    def __init__(self, *, multisortingcomparison, sorternames=None, drawlabels=False, node_cmap='viridis',
                 edge_cmap='hot', title='', figure=None, ax=None):
        BaseWidget.__init__(self, figure, ax)
        self._msc = multisortingcomparison
        self._sorter_names = sorternames
        self._drawlabels = drawlabels
        self._node_cmap = node_cmap
        self._edge_cmap = edge_cmap
        self._title = title

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
                                   node_color=nodes_col, ax=self.ax)
        _ = nx.draw_networkx_edges(g, pos=nx.circular_layout((sorted(g))), nodelist=sorted(g.nodes),
                                   edge_color=edge_col,
                                   edge_cmap=plt.cm.get_cmap(self._edge_cmap), edge_vmin=self._msc._min_accuracy,
                                   edge_vmax=1, ax=self.ax)
        if self._drawlabels:
            _ = nx.draw_networkx_labels(g, pos=nx.circular_layout((sorted(g))), nodelist=sorted(g.nodes), ax=self.ax)

        self.ax.axis('off')
