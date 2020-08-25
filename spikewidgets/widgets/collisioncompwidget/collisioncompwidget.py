import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from spikewidgets.widgets.basewidget import BaseWidget, BaseMultiWidget

import matplotlib.gridspec as gridspec

def plot_collision_gt_comparison(comp, unit_ids=None, nbins=10, figure=None, ax=None):
    """
    Plots multi sorting comparison graph.

    Parameters
    ----------
    comp: CollisionGTComparison
        The collision ground truth comparison object
    unit_ids: list
        List of considered units
    nbins: int
        Number of bins
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
    W = CollisionGTComparisonWidget(
        comp=comp,
        unit_ids=unit_ids,
        nbins=nbins,
        figure=figure,
        ax=ax
    )
    W.plot()
    return W


class CollisionGTComparisonWidget(BaseWidget):
    def __init__(self, comp, unit_ids=None, nbins=10, figure=None, ax=None):
        BaseWidget.__init__(self, figure, ax)
        if unit_ids is None:
            # take all units
            unit_ids = comp.sorting1.get_unit_ids()

        self.comp = comp
        self.unit_ids = unit_ids
        self.nbins = nbins
    
    def plot(self):
        self._do_plot()

    def _do_plot(self):
        fig = self.figure
        
        for ax in fig.axes:
            ax.remove()
        
        n = len(self.unit_ids)
        gs = gridspec.GridSpec(ncols=n, nrows=n, figure=fig)
        
        axs = np.empty((n,n), dtype=object)
        ax = None
        for r in range(n):
            for c in range(n):
                ax = fig.add_subplot(gs[r, c], sharex=ax, sharey=ax)
                if c > 0:
                    plt.setp(ax.get_yticklabels(), visible=False)
                if r < n-1:
                    plt.setp(ax.get_xticklabels(), visible=False)
                axs[r, c] = ax
        
        fs = self.comp.sorting1.get_sampling_frequency()
        
        for r in range(n):
            for c in range(r+1, n):
                
                u1 = self.unit_ids[r]
                u2 = self.unit_ids[c]
                
                bins, tp_count1, fn_count1, tp_count2, fn_count2 = self.comp.get_label_count_per_collision_bins(u1, u2, nbins=self.nbins)
                
                width = (bins[1] - bins[0]) / fs * 1000.
                lags = bins[:-1] / fs * 1000
                
                ax = axs[r, c]
                ax.bar(lags, tp_count1, width=width,  color='g')
                ax.bar(lags, fn_count1, width=width, bottom=tp_count1, color='r')
                
                ax = axs[c, r]
                ax.bar(lags, tp_count2, width=width,  color='g')
                ax.bar(lags, fn_count2, width=width, bottom=tp_count2, color='r')
        
        for r in range(n):
            ax = axs[r, 0]
            u1 = self.unit_ids[r]
            ax.set_ylabel(f'gt id{u1}')

        for c in range(n):
            ax = axs[0, c]
            u2 = self.unit_ids[c]
            ax.set_title(f'collision with \ngt id{u2}')
        
        ax = axs[0, 0]
        ax.set_xlabel('collision lag [ms]')

