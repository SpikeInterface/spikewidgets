import numpy as np
import sklearn

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec


from spikewidgets.widgets.basewidget import BaseWidget, BaseMultiWidget



def plot_comparison_collision_pair_by_pair(comp, unit_ids=None, nbins=10, figure=None, ax=None):
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
    W = ComparisonCollisionPairByPairWidget(
        comp=comp,
        unit_ids=unit_ids,
        nbins=nbins,
        figure=figure,
        ax=ax
    )
    W.plot()
    return W


def plot_comparison_collision_by_similarity(comp, templates, unit_ids=None, nbins=10, figure=None, ax=None):
    
    W = ComparisonCollisionBySimilarityWidget(
        comp=comp,
        templates=templates,
        unit_ids=unit_ids,
        nbins=nbins,
        figure=figure,
        ax=ax
    )
    W.plot()
    return W



class ComparisonCollisionPairByPairWidget(BaseWidget):
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
        
        ax = axs[-1, 0]
        ax.set_xlabel('collision lag [ms]')



class ComparisonCollisionBySimilarityWidget(BaseWidget):
    def __init__(self, comp, templates, unit_ids=None, nbins=10, figure=None, ax=None):
        BaseWidget.__init__(self, figure, ax)
        if unit_ids is None:
            # take all units
            unit_ids = comp.sorting1.get_unit_ids()

        self.comp = comp
        self.templates = templates
        self.unit_ids = unit_ids
        self.nbins = nbins
    
    def plot(self):
        self._do_plot()

    def _do_plot(self):
        #~ fig = self.figure
        ax = self.ax
        #~ for ax in fig.axes:
            #~ ax.remove()
        
        # compute similarity
        # take index of temmplate (respect unit_ids order)
        all_unit_ids = list(self.comp.sorting1.get_unit_ids())
        template_inds = [all_unit_ids.index(u) for u in self.unit_ids] 
        #~ print('template_inds', template_inds)
        templates = self.templates[template_inds, :, :].copy()
        flat_templates = templates.reshape(templates.shape[0], -1)
        similarity_matrix = sklearn.metrics.pairwise.cosine_distances(flat_templates)
        #~ print(similarity_matrix)

        n = len(self.unit_ids)
        
        
        fs = self.comp.sorting1.get_sampling_frequency()
        accuracies = []
        similarities = []
        for r in range(n):
            for c in range(r+1, n):
                
                u1 = self.unit_ids[r]
                u2 = self.unit_ids[c]
                
                bins, tp_count1, fn_count1, tp_count2, fn_count2 = self.comp.get_label_count_per_collision_bins(u1, u2, nbins=self.nbins)
                
                width = (bins[1] - bins[0]) / fs * 1000.
                lags = bins[:-1] / fs * 1000
                
                accuracy1 = tp_count1 / (tp_count1 + fn_count1)
                accuracies.append(accuracy1)
                similarities.append(similarity_matrix[r, c])

                accuracy2 = tp_count2 / (tp_count2 + fn_count1)
                accuracies.append(accuracy2)
                similarities.append(similarity_matrix[r, c])

        accuracies = np.array(accuracies)
        similarities = np.array(similarities)
        print(similarities)
        print(accuracies)
        order = np.argsort(similarities)
        
        similarities = similarities[order]
        accuracies = accuracies[order, :]
        
        self.ax.matshow(accuracies)
        
        
        
        
        

