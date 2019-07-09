from spikewidgets.widgets.basewidget import BaseWidget


def plot_sorting_accuracy(gt_sorting_comparison, property_name, figure=None, ax=None):
    W = SortingAccuracyWidget(
        sorting_comparison=gt_sorting_comparison,
        property_name=property_name,
        figure=figure,
        ax=ax
    )
    W.plot()
    return W


class SortingAccuracyWidget(BaseWidget):
    def __init__(self, *, sorting_comparison, property_name=None, figure=None, ax=None):
        BaseWidget.__init__(self, figure, ax)
        self._SC = sorting_comparison
        self._property_name = property_name

    def plot(self):
        self._do_plot()

    def _do_plot(self):
        SC = self._SC
        units = SC.sorting2.get_unit_ids()
        agreements = [SC.get_agreement_fraction(unit) for unit in units]
        if self._property_name is not None:
            assert self._property_name in SC.sorting1.get_unit_property_names(), "%s should be " \
                                                                                 "a property of the ground truth " \
                                                                                 "sorting extractor"
            xvals = SC.sorting1.get_units_property(unit_ids=units, property_name=self._property_name)
            self.ax.plot(xvals, agreements, '.')
            self.ax.set_xlabel(self._property_name)
        else:
            self.ax.plot(agreements, '.')
            self.ax.set_xticks([])
        self.ax.set_ylabel('Accuracy')
        self.ax.set_ylim([0, 1.05])
