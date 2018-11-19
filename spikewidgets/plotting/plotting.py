import matplotlib.pylab as plt
import numpy as np


def raster_plot(sorting, fs=None, ax=None, color_st=None, color=None,
                marker='|', mew=2, markersize=5, fontsize=10):
    '''

    Parameters
    ----------
    sorting
    fs
    ax
    color_st
    color
    marker
    mew
    markersize
    fontsize

    Returns
    -------

    '''
    import matplotlib.pylab as plt
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    for u_i, unit in enumerate(sorting.getUnitIds()):
        spiketrain = sorting.getUnitSpikeTrain(unit)
        if fs is not None:
            spiketrain = np.array(spiketrain) / float(fs)
        if np.any(color_st):
            import seaborn as sns
            colors = sns.color_palette("Paired", len(color_st))
            if i in color_st:
                idx = np.where(color_st == i)[0][0]
                ax.plot(spiketrain, u_i * np.ones_like(spiketrain), marker=marker,
                        mew=mew, color=colors[idx], markersize=5, ls='')
            else:
                ax.plot(spiketrain, u_i * np.ones_like(spiketrain), 'k', marker=marker,
                        mew=mew, markersize=markersize, ls='')
        elif color is not None:
            if isinstance(color, list) or isinstance(color, np.ndarray):
                ax.plot(spiketrain, u_i * np.ones_like(spiketrain), color=color[i], marker=marker,
                        mew=mew, markersize=markersize, ls='')
            else:
                ax.plot(spiketrain, u_i * np.ones_like(spiketrain), color=color, marker=marker,
                        mew=mew, markersize=markersize, ls='')
        else:
            ax.plot(spiketrain, u_i * np.ones(len(spiketrain)), 'k', marker=marker,
                    mew=mew, markersize=markersize, ls='')
    ax.axis('tight')
    ax.set_xlabel('Time (ms)', fontsize=fontsize)
    ax.set_ylabel('Spike Train Index', fontsize=fontsize)
    plt.gca().tick_params(axis='both', which='major')
    return ax


def matched_raster_plot(st, bintype=False, ax=None, overlap=False, labels=False, color_st=None, color=None, fs=10,
                        marker='|', mew=2, markersize=5):
    '''

    Parameters
    ----------
    st
    bintype
    ax

    Returns
    -------

    '''
    import matplotlib.pylab as plt
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    for i, spiketrain in enumerate(st):
        t = spiketrain.rescale(pq.s)
        if bintype:
            if spiketrain.annotations['bintype'] == 'EXCIT':
                ax.plot(t, i * np.ones_like(t), 'b', marker=marker, mew=mew, markersize=markersize, ls='')
            elif spiketrain.annotations['bintype'] == 'INHIB':
                ax.plot(t, i * np.ones_like(t), 'r', marker=marker, mew=mew, markersize=markersize, ls='')
        else:
            if not overlap and not labels:
                if np.any(color_st):
                    import seaborn as sns
                    colors = sns.color_palette("Paired", len(color_st))
                    if i in color_st:
                        idx = np.where(color_st == i)[0][0]
                        ax.plot(t, i * np.ones_like(t), marker=marker, mew=mew, color=colors[idx], markersize=5, ls='')
                    else:
                        ax.plot(t, i * np.ones_like(t), 'k', marker=marker, mew=mew, markersize=markersize, ls='')
                elif color is not None:
                    if isinstance(color, list) or isinstance(color, np.ndarray):
                        ax.plot(t, i * np.ones_like(t), color=color[i], marker=marker, mew=mew, markersize=markersize,
                                ls='')
                    else:
                        ax.plot(t, i * np.ones_like(t), color=color, marker=marker, mew=mew, markersize=markersize,
                                ls='')
                else:
                    ax.plot(t, i * np.ones_like(t), 'k', marker=marker, mew=mew, markersize=markersize, ls='')
            elif overlap:
                for j, t_sp in enumerate(spiketrain):
                    if spiketrain.annotations['overlap'][j] == 'SO':
                        ax.plot(t_sp, i, 'r', marker=marker, mew=mew, markersize=markersize, ls='')
                    elif spiketrain.annotations['overlap'][j] == 'O':
                        ax.plot(t_sp, i, 'g', marker=marker, mew=mew, markersize=markersize, ls='')
                    elif spiketrain.annotations['overlap'][j] == 'NO':
                        ax.plot(t_sp, i, 'k', marker=marker, mew=mew, markersize=markersize, ls='')
            elif labels:
                for j, t_sp in enumerate(spiketrain):
                    if 'TP' in spiketrain.annotations['labels'][j]:
                        ax.plot(t_sp, i, 'g', marker=marker, mew=mew, markersize=markersize, ls='')
                    elif 'CL' in spiketrain.annotations['labels'][j]:
                        ax.plot(t_sp, i, 'y', marker=marker, mew=mew, markersize=markersize, ls='')
                    elif 'FN' in spiketrain.annotations['labels'][j]:
                        ax.plot(t_sp, i, 'r', marker=marker, mew=mew, markersize=markersize, ls='')
                    elif 'FP' in spiketrain.annotations['labels'][j]:
                        ax.plot(t_sp, i, 'm', marker=marker, mew=mew, markersize=markersize, ls='')
                    else:
                        ax.plot(t_sp, i, 'k', marker=marker, mew=mew, markersize=markersize, ls='')

    ax.axis('tight')
    ax.set_xlim([st[0].t_start.rescale(pq.s), st[0].t_stop.rescale(pq.s)])
    ax.set_xlabel('Time (ms)', fontsize=fs)
    ax.set_ylabel('Spike Train Index', fontsize=fs)
    plt.gca().tick_params(axis='both', which='major')

    return ax
