import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np


class BaseWidget:
    def __init__(self, figure=None, ax=None):
        if figure is None and ax is None:
            self.figure = plt.figure()
            self.ax = self.figure.add_subplot(111)
        elif ax is None:
            self.figure = figure
            self.ax = self.figure.add_subplot(111)
        else:
            self.figure = ax.get_figure()
            self.ax = ax
        self.name = None
    
    def get_figure(self):
        return self.figure
    
    def get_ax(self):
        return self.ax
    
    def get_name(self):
        return self.name

class BaseMultiWidget:
    def __init__(self, figure=None, ax=None):
        if figure is None and ax is None:
            self.figure = plt.figure()
            self.ax = self.figure.add_subplot(111)
        elif ax is None:
            self.figure = figure
            self.ax = self.figure.add_subplot(111)
        else:
            self.figure = ax.get_figure()
            self.ax = ax
        self.axes = []
        self._gs = None
        self.ax.axis('off')
        self.name = None

    def get_tiled_ax(self, i, nrows, ncols, hspace=0.3, wspace=0.3, is_diag=False):
        if self._gs is None:
            self._gs = gridspec.GridSpecFromSubplotSpec(int(nrows), int(ncols), subplot_spec=self.ax,
                                                        hspace=hspace, wspace=wspace)
        r = int(i // ncols)
        c = int(np.mod(i, ncols))
        gs_sel = self._gs[r, c]
        ax = self.figure.add_subplot(gs_sel)
        self.axes.append(ax)
        if r == c:
            diag = True
        else:
            diag = False
        if is_diag:
            return ax, diag
        else:
            return ax
