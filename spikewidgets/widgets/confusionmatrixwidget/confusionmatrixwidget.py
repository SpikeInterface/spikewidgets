import numpy as np
from matplotlib import pyplot as plt


def plot_confusion_matrix(sortingcomparison, sorter_names=None, count_text=True, title=''):
    W = ConfusionMatrixWidget(
        sortingcomparison=sortingcomparison,
        sorternames=sorter_names,
        title=title,
        count_text=count_text
    )
    W.plot()


class ConfusionMatrixWidget:
    def __init__(self, *, sortingcomparison, sorternames=None, count_text=True, title=''):
        self._sc = sortingcomparison
        self._sorter_names = sorternames
        self._count_text = count_text
        self._title = title

    def plot(self):
        self._do_plot()

    def _do_plot(self):
        confusion_matrix, st1_idx, st2_idx = self._sc.get_confusion_matrix()

        sorting1 = self._sc._sorting1
        sorting2 = self._sc._sorting2
        unit1_ids = sorting1.get_unit_ids()
        unit2_ids = sorting2.get_unit_ids()
        N1 = len(unit1_ids)
        N2 = len(unit2_ids)

        fig, ax = plt.subplots()

        # Using matshow here just because it sets the ticks up nicely. imshow is faster.
        ax.matshow(confusion_matrix, cmap='Greens')

        if self._count_text:
            for (i, j), z in np.ndenumerate(confusion_matrix):
                if z != 0:
                    if z > np.max(confusion_matrix) / 2.:
                        ax.text(j, i, '{:d}'.format(z), ha='center', va='center', color='white')
                    else:
                        ax.text(j, i, '{:d}'.format(z), ha='center', va='center', color='black')

        ax.axhline(int(N1 - 1) + 0.5, color='black')
        ax.axvline(int(N2 - 1) + 0.5, color='black')

        # Major ticks
        ax.set_xticks(np.arange(0, N2 + 1))
        ax.set_yticks(np.arange(0, N1 + 1))
        ax.xaxis.tick_bottom()
        # Labels for major ticks
        ax.set_xticklabels(np.append(st1_idx, 'FN'), fontsize=12)
        ax.set_yticklabels(np.append(st2_idx, 'FP'), fontsize=12)

        if self._sorter_names is None:
            ax.set_xlabel(self._sc.sorting2_name, fontsize=20)
        else:
            assert len(self._sorter_names) == 2
            ax.set_xlabel(self._sorter_names[0], fontsize=20)

        if self._sorter_names is None:
            ax.set_ylabel(self._sc.sorting1_name, fontsize=20)
        else:
            assert len(self._sorter_names) == 2
            ax.set_xlabel(self._sorter_names[0], fontsize=20)
