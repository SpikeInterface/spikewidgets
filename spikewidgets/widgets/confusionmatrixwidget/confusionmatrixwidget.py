import numpy as np
from matplotlib import pyplot as plt
from spikewidgets.widgets.basewidget import BaseWidget


def plot_confusion_matrix(sortingcomparison, sorter_names=None, count_text=True, title='', ax=None, figure=None):
    W = ConfusionMatrixWidget(
        sortingcomparison=sortingcomparison,
        sorternames=sorter_names,
        title=title,
        count_text=count_text,
        figure=figure,
        ax=ax,
    )
    W.plot()
    return W


class ConfusionMatrixWidget(BaseWidget):
    def __init__(self, *, sortingcomparison, sorternames=None, count_text=True, title='', figure=None, ax=None):
        BaseWidget.__init__(self, figure, ax)
        self._sc = sortingcomparison
        self._sorter_names = sorternames
        self._count_text = count_text
        self._title = title

    def plot(self):
        self._do_plot()

    def _do_plot(self):
        confusion_matrix, st1_idx, st2_idx = self._sc.get_confusion_matrix()

        sorting1 = self._sc.sorting1
        sorting2 = self._sc.sorting2
        unit1_ids = sorting1.get_unit_ids()
        unit2_ids = sorting2.get_unit_ids()
        N1 = len(unit1_ids)
        N2 = len(unit2_ids)

        # Using matshow here just because it sets the ticks up nicely. imshow is faster.
        self.ax.matshow(confusion_matrix, cmap='Greens')

        if self._count_text:
            for (i, j), z in np.ndenumerate(confusion_matrix):
                if z != 0:
                    if z > np.mself.ax(confusion_matrix) / 2.:
                        self.ax.text(j, i, '{:d}'.format(z), ha='center', va='center', color='white')
                    else:
                        self.ax.text(j, i, '{:d}'.format(z), ha='center', va='center', color='black')

        self.ax.axhline(int(N1 - 1) + 0.5, color='black')
        self.ax.axvline(int(N2 - 1) + 0.5, color='black')

        # Major ticks
        self.ax.set_xticks(np.arange(0, N2 + 1))
        self.ax.set_yticks(np.arange(0, N1 + 1))
        self.ax.xaxis.tick_bottom()
        # Labels for major ticks
        self.ax.set_xticklabels(np.append(st1_idx, 'FN'), fontsize=12)
        self.ax.set_yticklabels(np.append(st2_idx, 'FP'), fontsize=12)

        if self._sorter_names is None:
            self.ax.set_xlabel(self._sc.sorting2_name, fontsize=20)
        else:
            assert len(self._sorter_names) == 2
            self.ax.set_xlabel(self._sorter_names[0], fontsize=20)

        if self._sorter_names is None:
            self.ax.set_ylabel(self._sc.sorting1_name, fontsize=20)
        else:
            assert len(self._sorter_names) == 2
            self.ax.set_xlabel(self._sorter_names[0], fontsize=20)
