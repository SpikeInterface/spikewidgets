from matplotlib import pyplot as plt


class SortingAccuracyWidget:
    def __init__(self, *, sorting_comparison, property_name=None):
        self._SC = sorting_comparison
        self._property_name = property_name

    def plot(self):
        self._do_plot()

    def _do_plot(self):
        SC = self._SC
        units = SC.get_sorting1().get_unit_ids()
        agreements = [SC.get_agreement_fraction(unit) for unit in units]
        if self._property_name:
            xvals = SC.get_sorting1().get_units_property(unit_ids=units, property_name=self._property_name)
            plt.plot(xvals, agreements, '.')
            plt.xlabel(self._property_name)
        else:
            plt.plot(agreements, '.')
            plt.xticks([])
        plt.ylabel('Accuracy')
        plt.show()
