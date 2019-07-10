from spikewidgets.widgets.basewidget import BaseWidget
import spiketoolkit as st


def plot_sorting_performance(gt_sorting_comparison, property_name, metric='accuracy', figure=None, ax=None):
    W = SortingPerformanceWidget(
        sorting_comparison=gt_sorting_comparison,
        property_name=property_name,
        figure=figure,
        metric=metric,
        ax=ax
    )
    W.plot()
    return W


class SortingPerformanceWidget(BaseWidget):
    def __init__(self, *, sorting_comparison, property_name=None, metric='accuracy', figure=None, ax=None):
        assert isinstance(sorting_comparison, st.comparison.groundtruthcomparison.GroundTruthComparison), \
            "The 'sorting_comparison' object should be a GroundTruthComparison instance"
        BaseWidget.__init__(self, figure, ax)
        self._SC = sorting_comparison
        self._property_name = property_name
        self._metric = metric
        self.name = 'SortingPerformance'

    def plot(self):
        self._do_plot()

    def _do_plot(self):
        SC = self._SC
        units = SC.sorting1.get_unit_ids()
        perf = SC.get_performance()[self._metric]
        if self._property_name is not None:
            assert self._property_name in SC.sorting1.get_unit_property_names(), "%s should be " \
                                                                                 "a property of the ground truth " \
                                                                                 "sorting extractor"
            xvals = SC.sorting1.get_units_property(unit_ids=units, property_name=self._property_name)
            self.ax.plot(xvals, perf, '.', markersize=10)
            self.ax.set_xlabel(self._property_name)
        else:
            self.ax.plot(perf, '.')
            self.ax.set_xticks([])
        self.ax.set_ylabel(self._metric)
        self.ax.set_ylim([-0.05, 1.05])
