# being a bit too dynamic
# pylint: disable=E1101
from itertools import izip
import datetime
import re

import numpy as np

from pandas.util.decorators import cache_readonly
import pandas.core.common as com
from pandas.core.index import MultiIndex
from pandas.core.series import Series, remove_na
from pandas.tseries.index import DatetimeIndex
from pandas.tseries.period import PeriodIndex
from pandas.tseries.frequencies import get_period_alias, get_base_alias
from pandas.tseries.offsets import DateOffset

try: # mpl optional
    import pandas.tseries.converter as conv
    conv.register()
except ImportError:
    pass

def _get_standard_kind(kind):
    return {'density' : 'kde'}.get(kind, kind)


def scatter_matrix(frame, alpha=0.5, figsize=None, ax=None, grid=False,
                   diagonal='hist', marker='.', **kwds):
    """
    Draw a matrix of scatter plots.

    Parameters
    ----------
    alpha : amount of transparency applied
    figsize : a tuple (width, height) in inches
    ax : Matplotlib axis object
    grid : setting this to True will show the grid
    diagonal : pick between 'kde' and 'hist' for
        either Kernel Density Estimation or Histogram
        plot in the diagonal
    kwds : other plotting keyword arguments
        To be passed to scatter function

    Examples
    --------
    >>> df = DataFrame(np.random.randn(1000, 4), columns=['A','B','C','D'])
    >>> scatter_matrix(df, alpha=0.2)
    """
    from matplotlib.artist import setp

    df = frame._get_numeric_data()
    n = df.columns.size
    fig, axes = _subplots(nrows=n, ncols=n, figsize=figsize, ax=ax,
                          squeeze=False)

    # no gaps between subplots
    fig.subplots_adjust(wspace=0, hspace=0)

    mask = com.notnull(df)

    marker = _get_marker_compat(marker)

    for i, a in zip(range(n), df.columns):
        for j, b in zip(range(n), df.columns):
            ax = axes[i, j]

            if i == j:
                values = df[a].values[mask[a].values]

                # Deal with the diagonal by drawing a histogram there.
                if diagonal == 'hist':
                    ax.hist(values)
                elif diagonal in ('kde', 'density'):
                    from scipy.stats import gaussian_kde
                    y = values
                    gkde = gaussian_kde(y)
                    ind = np.linspace(y.min(), y.max(), 1000)
                    ax.plot(ind, gkde.evaluate(ind), **kwds)
            else:
                common = (mask[a] & mask[b]).values

                ax.scatter(df[b][common], df[a][common],
                                   marker=marker, alpha=alpha, **kwds)

            ax.set_xlabel('')
            ax.set_ylabel('')

            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)

            # setup labels
            if i == 0 and j % 2 == 1:
                ax.set_xlabel(b, visible=True)
                ax.xaxis.set_visible(True)
                ax.set_xlabel(b)
                ax.xaxis.set_ticks_position('top')
                ax.xaxis.set_label_position('top')
                setp(ax.get_xticklabels(), rotation=90)
            elif i == n - 1 and j % 2 == 0:
                ax.set_xlabel(b, visible=True)
                ax.xaxis.set_visible(True)
                ax.set_xlabel(b)
                ax.xaxis.set_ticks_position('bottom')
                ax.xaxis.set_label_position('bottom')
                setp(ax.get_xticklabels(), rotation=90)
            elif j == 0 and i % 2 == 0:
                ax.set_ylabel(a, visible=True)
                ax.yaxis.set_visible(True)
                ax.set_ylabel(a)
                ax.yaxis.set_ticks_position('left')
                ax.yaxis.set_label_position('left')
            elif j == n - 1 and i % 2 == 1:
                ax.set_ylabel(a, visible=True)
                ax.yaxis.set_visible(True)
                ax.set_ylabel(a)
                ax.yaxis.set_ticks_position('right')
                ax.yaxis.set_label_position('right')

            # ax.grid(b=grid)

    axes[0, 0].yaxis.set_visible(False)
    axes[n-1, n-1].xaxis.set_visible(False)
    axes[n-1, n-1].yaxis.set_visible(False)
    axes[0, n - 1].yaxis.tick_right()

    for ax in axes.flat:
        setp(ax.get_xticklabels(), fontsize=8)
        setp(ax.get_yticklabels(), fontsize=8)

    return axes


def _gca():
    import matplotlib.pyplot as plt
    return plt.gca()

def _gcf():
    import matplotlib.pyplot as plt
    return plt.gcf()

def _get_marker_compat(marker):
    import matplotlib.lines as mlines
    import matplotlib as mpl
    if mpl.__version__ < '1.1.0' and marker == '.':
        return 'o'
    if marker not in mlines.lineMarkers:
        return 'o'
    return marker

def radviz(frame, class_column, ax=None, **kwds):
    """RadViz - a multivariate data visualization algorithm

    Parameters:
    -----------
    frame: DataFrame object
    class_column: Column name that contains information about class membership
    ax: Matplotlib axis object, optional
    kwds: Matplotlib scatter method keyword arguments, optional

    Returns:
    --------
    ax: Matplotlib axis object
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import matplotlib.text as text
    import random

    def random_color(column):
        random.seed(column)
        return [random.random() for _ in range(3)]

    def normalize(series):
        a = min(series)
        b = max(series)
        return (series - a) / (b - a)

    column_names = [column_name for column_name in frame.columns
                    if column_name != class_column]

    df = frame[column_names].apply(normalize)

    if ax is None:
        ax = plt.gca(xlim=[-1, 1], ylim=[-1, 1])

    classes = set(frame[class_column])
    to_plot = {}

    for class_ in classes:
        to_plot[class_] = [[], []]

    n = len(frame.columns) - 1
    s = np.array([(np.cos(t), np.sin(t))
                  for t in [2.0 * np.pi * (i / float(n))
                            for i in range(n)]])

    for i in range(len(frame)):
        row = df.irow(i).values
        row_ = np.repeat(np.expand_dims(row, axis=1), 2, axis=1)
        y = (s * row_).sum(axis=0) / row.sum()
        class_name = frame[class_column].iget(i)
        to_plot[class_name][0].append(y[0])
        to_plot[class_name][1].append(y[1])

    for class_ in classes:
        line = ax.scatter(to_plot[class_][0],
                          to_plot[class_][1],
                          color=random_color(class_),
                          label=com._stringify(class_), **kwds)
    ax.legend()

    ax.add_patch(patches.Circle((0.0, 0.0), radius=1.0, facecolor='none'))

    for xy, name in zip(s, column_names):

        ax.add_patch(patches.Circle(xy, radius=0.025, facecolor='gray'))

        if xy[0] < 0.0 and xy[1] < 0.0:
            ax.text(xy[0] - 0.025, xy[1] - 0.025, name,
                    ha='right', va='top', size='small')
        elif xy[0] < 0.0 and xy[1] >= 0.0:
            ax.text(xy[0] - 0.025, xy[1] + 0.025, name,
                    ha='right', va='bottom', size='small')
        elif xy[0] >= 0.0 and xy[1] < 0.0:
            ax.text(xy[0] + 0.025, xy[1] - 0.025, name,
                    ha='left', va='top', size='small')
        elif xy[0] >= 0.0 and xy[1] >= 0.0:
            ax.text(xy[0] + 0.025, xy[1] + 0.025, name,
                    ha='left', va='bottom', size='small')

    ax.axis('equal')
    return ax

def andrews_curves(data, class_column, ax=None, samples=200):
    """
    Parameters:
    data: A DataFrame containing data to be plotted, preferably
    normalized to (0.0, 1.0).
    class_column: Name of the column containing class names.
    samples: Number of points to plot in each curve.
    """
    from math import sqrt, pi, sin, cos
    import matplotlib.pyplot as plt
    import random
    def function(amplitudes):
        def f(x):
            x1 = amplitudes[0]
            result = x1 / sqrt(2.0)
            harmonic = 1.0
            for x_even, x_odd in zip(amplitudes[1::2], amplitudes[2::2]):
                result += (x_even * sin(harmonic * x) +
                            x_odd * cos(harmonic * x))
                harmonic += 1.0
            if len(amplitudes) % 2 != 0:
                result += amplitudes[-1] * sin(harmonic * x)
            return result
        return f
    def random_color(column):
        random.seed(column)
        return [random.random() for _ in range(3)]
    n = len(data)
    classes = set(data[class_column])
    class_col = data[class_column]
    columns = [data[col] for col in data.columns if (col != class_column)]
    x = [-pi + 2.0 * pi * (t / float(samples)) for t in range(samples)]
    used_legends = set([])
    if ax == None:
        ax = plt.gca(xlim=(-pi, pi))
    for i in range(n):
        row = [columns[c][i] for c in range(len(columns))]
        f = function(row)
        y = [f(t) for t in x]
        label = None
        if com._stringify(class_col[i]) not in used_legends:
            label = com._stringify(class_col[i])
            used_legends.add(label)
        ax.plot(x, y, color=random_color(class_col[i]), label=label)
    ax.legend(loc='upper right')
    ax.grid()
    return ax

def bootstrap_plot(series, fig=None, size=50, samples=500, **kwds):
    """Bootstrap plot.

    Parameters:
    -----------
    series: Time series
    fig: matplotlib figure object, optional
    size: number of data points to consider during each sampling
    samples: number of times the bootstrap procedure is performed
    kwds: optional keyword arguments for plotting commands, must be accepted by both hist and plot

    Returns:
    --------
    fig: matplotlib figure
    """
    import random
    import matplotlib
    import matplotlib.pyplot as plt
    data = series.values
    samplings = [random.sample(data, size) for _ in range(samples)]
    means = np.array([np.mean(sampling) for sampling in samplings])
    medians = np.array([np.median(sampling) for sampling in samplings])
    midranges = np.array([(min(sampling) + max(sampling)) * 0.5 for sampling in samplings])
    if fig == None:
        fig = plt.figure()
    x = range(samples)
    axes = []
    ax1 = fig.add_subplot(2, 3, 1)
    ax1.set_xlabel("Sample")
    axes.append(ax1)
    ax1.plot(x, means, **kwds)
    ax2 = fig.add_subplot(2, 3, 2)
    ax2.set_xlabel("Sample")
    axes.append(ax2)
    ax2.plot(x, medians, **kwds)
    ax3 = fig.add_subplot(2, 3, 3)
    ax3.set_xlabel("Sample")
    axes.append(ax3)
    ax3.plot(x, midranges, **kwds)
    ax4 = fig.add_subplot(2, 3, 4)
    ax4.set_xlabel("Mean")
    axes.append(ax4)
    ax4.hist(means, **kwds)
    ax5 = fig.add_subplot(2, 3, 5)
    ax5.set_xlabel("Median")
    axes.append(ax5)
    ax5.hist(medians, **kwds)
    ax6 = fig.add_subplot(2, 3, 6)
    ax6.set_xlabel("Midrange")
    axes.append(ax6)
    ax6.hist(midranges, **kwds)
    for axis in axes:
        plt.setp(axis.get_xticklabels(), fontsize=8)
        plt.setp(axis.get_yticklabels(), fontsize=8)
    return fig

def parallel_coordinates(data, class_column, cols=None, ax=None, **kwds):
    """Parallel coordinates plotting.

    Parameters:
    -----------
    data: A DataFrame containing data to be plotted
    class_column: Column name containing class names
    cols: A list of column names to use, optional
    ax: matplotlib axis object, optional
    kwds: A list of keywords for matplotlib plot method

    Returns:
    --------
    ax: matplotlib axis object
    """
    import matplotlib.pyplot as plt
    import random
    def random_color(column):
        random.seed(column)
        return [random.random() for _ in range(3)]
    n = len(data)
    classes = set(data[class_column])
    class_col = data[class_column]

    if cols is None:
        df = data.drop(class_column, axis=1)
    else:
        df = data[cols]

    used_legends = set([])

    ncols = len(df.columns)
    x = range(ncols)

    if ax == None:
        ax = plt.gca()

    for i in range(n):
        row = df.irow(i).values
        y = row
        label = None
        kls = class_col.iget_value(i)
        if com._stringify(kls) not in used_legends:
            label = com._stringify(kls)
            used_legends.add(label)
        ax.plot(x, y, color=random_color(kls), label=label, **kwds)

    for i in range(ncols):
        ax.axvline(i, linewidth=1, color='black')

    ax.set_xticks(x)
    ax.set_xticklabels(df.columns)
    ax.legend(loc='upper right')
    ax.grid()
    return ax

def lag_plot(series, ax=None, **kwds):
    """Lag plot for time series.

    Parameters:
    -----------
    series: Time series
    ax: Matplotlib axis object, optional
    kwds: Matplotlib scatter method keyword arguments, optional

    Returns:
    --------
    ax: Matplotlib axis object
    """
    import matplotlib.pyplot as plt
    data = series.values
    y1 = data[:-1]
    y2 = data[1:]
    if ax == None:
        ax = plt.gca()
    ax.set_xlabel("y(t)")
    ax.set_ylabel("y(t + 1)")
    ax.scatter(y1, y2, **kwds)
    return ax

def autocorrelation_plot(series, ax=None):
    """Autocorrelation plot for time series.

    Parameters:
    -----------
    series: Time series
    ax: Matplotlib axis object, optional

    Returns:
    -----------
    ax: Matplotlib axis object
    """
    import matplotlib.pyplot as plt
    n = len(series)
    data = np.asarray(series)
    if ax == None:
        ax = plt.gca(xlim=(1, n), ylim=(-1.0, 1.0))
    mean = np.mean(data)
    c0 = np.sum((data - mean) ** 2) / float(n)
    def r(h):
        return ((data[:n - h] - mean) * (data[h:] - mean)).sum() / float(n) / c0
    x = np.arange(n) + 1
    y = map(r, x)
    z95 = 1.959963984540054
    z99 = 2.5758293035489004
    ax.axhline(y=z99/np.sqrt(n), linestyle='--', color='grey')
    ax.axhline(y=z95/np.sqrt(n), color='grey')
    ax.axhline(y=0.0, color='black')
    ax.axhline(y=-z95/np.sqrt(n), color='grey')
    ax.axhline(y=-z99/np.sqrt(n), linestyle='--', color='grey')
    ax.set_xlabel("Lag")
    ax.set_ylabel("Autocorrelation")
    ax.plot(x, y)
    ax.grid()
    return ax

def grouped_hist(data, column=None, by=None, ax=None, bins=50, log=False,
                 figsize=None, layout=None, sharex=False, sharey=False,
                 rot=90):
    """

    Returns
    -------
    fig : matplotlib.Figure
    """
    # if isinstance(data, DataFrame):
    #     data = data[column]

    def plot_group(group, ax):
        ax.hist(group.dropna(), bins=bins)

    fig, axes = _grouped_plot(plot_group, data, column=column,
                              by=by, sharex=sharex, sharey=sharey,
                              figsize=figsize, layout=layout, rot=rot)
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.9,
                        hspace=0.3, wspace=0.2)
    return fig

class MPLPlot(object):
    """
    Base class for assembling a pandas plot using matplotlib

    Parameters
    ----------
    data :

    """
    _default_rot = 0

    _pop_attributes = ['label', 'style', 'logy', 'logx', 'loglog']
    _attr_defaults = {'logy': False, 'logx': False, 'loglog': False}

    def __init__(self, data, kind=None, by=None, subplots=False, sharex=True,
                 sharey=False, use_index=True,
                 figsize=None, grid=None, legend=True, rot=None,
                 ax=None, fig=None, title=None, xlim=None, ylim=None,
                 xticks=None, yticks=None,
                 sort_columns=False, fontsize=None,
                 secondary_y=False, **kwds):

        self.data = data
        self.by = by

        self.kind = kind

        self.sort_columns = sort_columns

        self.subplots = subplots
        self.sharex = sharex
        self.sharey = sharey
        self.figsize = figsize

        self.xticks = xticks
        self.yticks = yticks
        self.xlim = xlim
        self.ylim = ylim
        self.title = title
        self.use_index = use_index

        self.fontsize = fontsize
        self.rot = rot

        if grid is None:
            grid = False if secondary_y else True

        self.grid = grid
        self.legend = legend

        for attr in self._pop_attributes:
            value = kwds.pop(attr, self._attr_defaults.get(attr, None))
            setattr(self, attr, value)

        self.ax = ax
        self.fig = fig
        self.axes = None

        if not isinstance(secondary_y, (bool, tuple, list, np.ndarray)):
            secondary_y = [secondary_y]
        self.secondary_y = secondary_y

        self.kwds = kwds

    def _iter_data(self):
        from pandas.core.frame import DataFrame
        if isinstance(self.data, (Series, np.ndarray)):
            yield self.label, np.asarray(self.data)
        elif isinstance(self.data, DataFrame):
            df = self.data

            if self.sort_columns:
                columns = com._try_sort(df.columns)
            else:
                columns = df.columns

            for col in columns:
                # # is this right?
                # empty = df[col].count() == 0
                # values = df[col].values if not empty else np.zeros(len(df))

                values = df[col].values
                yield col, values

    @property
    def nseries(self):
        if self.data.ndim == 1:
            return 1
        else:
            return self.data.shape[1]

    def draw(self):
        self.plt.draw_if_interactive()

    def generate(self):
        self._args_adjust()
        self._compute_plot_data()
        self._setup_subplots()
        self._make_plot()
        self._post_plot_logic()
        self._adorn_subplots()

    def _args_adjust(self):
        pass

    def _maybe_right_yaxis(self, ax):
        _types = (list, tuple, np.ndarray)
        sec_true = isinstance(self.secondary_y, bool) and self.secondary_y
        list_sec = isinstance(self.secondary_y, _types)
        has_sec = list_sec and len(self.secondary_y) > 0
        all_sec = list_sec and len(self.secondary_y) == self.nseries

        if (sec_true or has_sec) and not hasattr(ax, 'right_ax'):
            orig_ax, new_ax = ax, ax.twinx()
            orig_ax.right_ax, new_ax.left_ax = new_ax, orig_ax

            if len(orig_ax.get_lines()) == 0: # no data on left y
                orig_ax.get_yaxis().set_visible(False)

            if len(new_ax.get_lines()) == 0:
                new_ax.get_yaxis().set_visible(False)

            if sec_true or all_sec:
                ax = new_ax
        else:
            ax.get_yaxis().set_visible(True)

        return ax

    def _setup_subplots(self):
        if self.subplots:
            nrows, ncols = self._get_layout()
            if self.ax is None:
                fig, axes = _subplots(nrows=nrows, ncols=ncols,
                                      sharex=self.sharex, sharey=self.sharey,
                                      figsize=self.figsize,
                                      secondary_y=self.secondary_y,
                                      data=self.data)
            else:
                fig, axes = _subplots(nrows=nrows, ncols=ncols,
                                      sharex=self.sharex, sharey=self.sharey,
                                      figsize=self.figsize, ax=self.ax,
                                      secondary_y=self.secondary_y,
                                      data=self.data)
        else:
            if self.ax is None:
                fig = self.plt.figure(figsize=self.figsize)
                ax = fig.add_subplot(111)
                ax = self._maybe_right_yaxis(ax)
            else:
                fig = self.ax.get_figure()
                ax = self._maybe_right_yaxis(self.ax)

            axes = [ax]

        self.fig = fig
        self.axes = axes

    def _get_layout(self):
        return (len(self.data.columns), 1)

    def _compute_plot_data(self):
        pass

    def _make_plot(self):
        raise NotImplementedError

    def _post_plot_logic(self):
        pass

    def _adorn_subplots(self):
        to_adorn = self.axes

        # todo: sharex, sharey handling?

        for ax in to_adorn:
            if self.yticks is not None:
                ax.set_yticks(self.yticks)

            if self.xticks is not None:
                ax.set_xticks(self.xticks)

            if self.ylim is not None:
                ax.set_ylim(self.ylim)

            if self.xlim is not None:
                ax.set_xlim(self.xlim)

            ax.grid(self.grid)

        if self.title:
            if self.subplots:
                self.fig.suptitle(self.title)
            else:
                self.axes[0].set_title(self.title)

        if self._need_to_set_index:
            labels = [_stringify(key) for key in self.data.index]
            labels = dict(zip(range(len(self.data.index)), labels))

            for ax_ in self.axes:
                # ax_.set_xticks(self.xticks)
                xticklabels = [labels.get(x, '') for x in ax_.get_xticks()]
                ax_.set_xticklabels(xticklabels, rotation=self.rot)

    @property
    def legend_title(self):
        if hasattr(self.data, 'columns'):
            if not isinstance(self.data.columns, MultiIndex):
                name = self.data.columns.name
                if name is not None:
                    name = com._stringify(name)
                return name
            else:
                stringified = map(com._stringify,
                                  self.data.columns.names)
                return ','.join(stringified)
        else:
            return None

    @cache_readonly
    def plt(self):
        import matplotlib.pyplot as plt
        return plt

    _need_to_set_index = False

    def _get_xticks(self, convert_period=False):
        index = self.data.index
        is_datetype = index.inferred_type in ('datetime', 'date',
                                              'datetime64', 'time')

        if self.use_index:
            if convert_period and isinstance(index, PeriodIndex):
                index = index.to_timestamp()
                x = index._mpl_repr()
            elif index.is_numeric() or is_datetype:
                """
                Matplotlib supports numeric values or datetime objects as
                xaxis values. Taking LBYL approach here, by the time
                matplotlib raises exception when using non numeric/datetime
                values for xaxis, several actions are already taken by plt.
                """
                x = index._mpl_repr()
            else:
                self._need_to_set_index = True
                x = range(len(index))
        else:
            x = range(len(index))

        return x

    def _get_plot_function(self):
        if self.logy:
            plotf = self.plt.Axes.semilogy
        elif self.logx:
            plotf = self.plt.Axes.semilogx
        elif self.loglog:
            plotf = self.plt.Axes.loglog
        else:
            plotf = self.plt.Axes.plot

        return plotf

    def _get_index_name(self):
        if isinstance(self.data.index, MultiIndex):
            name = self.data.index.names
            if any(x is not None for x in name):
                name = ','.join([com._stringify(x) for x in name])
            else:
                name = None
        else:
            name = self.data.index.name
            if name is not None:
                name = com._stringify(name)

        return name

    def _get_ax(self, i):
        # get the twinx ax if appropriate
        if self.subplots:
            ax = self.axes[i]
        else:
            ax = self.axes[0]

        if self.on_right(i):
            if hasattr(ax, 'right_ax'):
                ax = ax.right_ax
        elif hasattr(ax, 'left_ax'):
            ax = ax.left_ax

        ax.get_yaxis().set_visible(True)
        return ax

    def on_right(self, i):
        from pandas.core.frame import DataFrame
        if isinstance(self.secondary_y, bool):
            return self.secondary_y

        if (isinstance(self.data, DataFrame) and
            isinstance(self.secondary_y, (tuple, list, np.ndarray))):
            return self.data.columns[i] in self.secondary_y

    def _get_style(self, i, col_name):
        style = ''
        if self.subplots:
            style = 'k'

        if self.style is not None:
            if isinstance(self.style, list):
                try:
                    style = self.style[i]
                except IndexError:
                    pass
            elif isinstance(self.style, dict):
                style = self.style.get(col_name, style)
            else:
                style = self.style

        return style or None

class KdePlot(MPLPlot):
    def __init__(self, data, **kwargs):
        MPLPlot.__init__(self, data, **kwargs)

    def _make_plot(self):
        from scipy.stats import gaussian_kde
        plotf = self._get_plot_function()
        for i, (label, y) in enumerate(self._iter_data()):
            ax = self._get_ax(i)
            style = self._get_style(i, label)

            label = com._stringify(label)

            gkde = gaussian_kde(y)
            sample_range = max(y) - min(y)
            ind = np.linspace(min(y) - 0.5 * sample_range,
                max(y) + 0.5 * sample_range, 1000)
            ax.set_ylabel("Density")

            y = gkde.evaluate(ind)
            kwds = self.kwds.copy()
            kwds['label'] = label
            if style is None:
                args = (ax, ind, y)
            else:
                args = (ax, ind, y, style)

            plotf(*args, **kwds)
            ax.grid(self.grid)

    def _post_plot_logic(self):
        if self.subplots and self.legend:
            for ax in self.axes:
                ax.legend(loc='best')

class LinePlot(MPLPlot):

    def __init__(self, data, **kwargs):
        self.mark_right = kwargs.pop('mark_right', True)
        MPLPlot.__init__(self, data, **kwargs)

    def _index_freq(self):
        from pandas.core.frame import DataFrame
        if isinstance(self.data, (Series, DataFrame)):
            freq = getattr(self.data.index, 'freq', None)
            if freq is None:
                freq = getattr(self.data.index, 'inferred_freq', None)
                if freq == 'B':
                    weekdays = np.unique(self.data.index.dayofweek)
                    if (5 in weekdays) or (6 in weekdays):
                        freq = None
            return freq

    def _is_dynamic_freq(self, freq):
        if isinstance(freq, DateOffset):
            freq = freq.rule_code
        else:
            freq = get_base_alias(freq)
        freq = get_period_alias(freq)
        return freq is not None

    def _use_dynamic_x(self):
        freq = self._index_freq()

        ax = self._get_ax(0)
        ax_freq = getattr(ax, 'freq', None)
        if freq is None: # convert irregular if axes has freq info
            freq = ax_freq
        else: # do not use tsplot if irregular was plotted first
            if (ax_freq is None) and (len(ax.get_lines()) > 0):
                return False

        return (freq is not None) and self._is_dynamic_freq(freq)

    def _get_colors(self):
        import matplotlib.pyplot as plt
        cycle = ''.join(plt.rcParams.get('axes.color_cycle', list('bgrcmyk')))
        has_colors = 'colors' in self.kwds
        colors = self.kwds.pop('colors', cycle)
        return has_colors, colors

    def _make_plot(self):
        # this is slightly deceptive
        if self.use_index and self._use_dynamic_x():
            data = self._maybe_convert_index(self.data)
            self._make_ts_plot(data, **self.kwds)
        else:
            lines = []
            labels = []
            x = self._get_xticks(convert_period=True)

            has_colors, colors = self._get_colors()
            def _maybe_add_color(kwargs, style, i):
                if (not has_colors and
                    (style is None or re.match('[a-z]+', style) is None)
                    and 'color' not in kwargs):
                    kwargs['color'] = colors[i % len(colors)]

            plotf = self._get_plot_function()

            for i, (label, y) in enumerate(self._iter_data()):
                ax = self._get_ax(i)
                style = self._get_style(i, label)
                kwds = self.kwds.copy()

                _maybe_add_color(kwds, style, i)

                label = _stringify(label)

                mask = com.isnull(y)
                if mask.any():
                    y = np.ma.array(y)
                    y = np.ma.masked_where(mask, y)

                kwds['label'] = label
                if style is None:
                    args = (ax, x, y)
                else:
                    args = (ax, x, y, style)

                newline = plotf(*args, **kwds)[0]
                lines.append(newline)
                leg_label = label
                if self.mark_right and self.on_right(i):
                    leg_label += ' (right)'
                labels.append(leg_label)
                ax.grid(self.grid)

            self._make_legend(lines, labels)

    def _make_ts_plot(self, data, **kwargs):
        from pandas.tseries.plotting import tsplot
        kwargs = kwargs.copy()
        has_colors, colors = self._get_colors()

        plotf = self._get_plot_function()
        lines = []
        labels = []

        def _maybe_add_color(kwargs, style, i):
            if (not has_colors and
                (style is None or re.match('[a-z]+', style) is None)):
                kwargs['color'] = colors[i % len(colors)]

        def to_leg_label(label, i):
            if self.mark_right and self.on_right(i):
                return label + ' (right)'
            return label

        if isinstance(data, Series):
            ax = self._get_ax(0) #self.axes[0]
            style = self.style or ''
            label = com._stringify(self.label)
            kwds = kwargs.copy()
            _maybe_add_color(kwds, style, 0)

            newlines = tsplot(data, plotf, ax=ax, label=label, style=self.style,
                             **kwds)
            ax.grid(self.grid)
            lines.append(newlines[0])
            leg_label = to_leg_label(label, 0)
            labels.append(leg_label)
        else:
            for i, col in enumerate(data.columns):
                label = com._stringify(col)
                ax = self._get_ax(i)
                style = self._get_style(i, col)
                kwds = kwargs.copy()

                _maybe_add_color(kwds, style, i)

                newlines = tsplot(data[col], plotf, ax=ax, label=label,
                                  style=style, **kwds)

                lines.append(newlines[0])
                leg_label = to_leg_label(label, i)
                labels.append(leg_label)
                ax.grid(self.grid)

        self._make_legend(lines, labels)

    def _make_legend(self, lines, labels):
        ax, leg = self._get_ax_legend(self.axes[0])

        if not self.subplots:
            if leg is not None:
                ext_lines = leg.get_lines()
                ext_labels = [x.get_text() for x in leg.get_texts()]
                ext_lines.extend(lines)
                ext_labels.extend(labels)
                ax.legend(ext_lines, ext_labels, loc='best',
                          title=self.legend_title)
            elif self.legend:
                ax.legend(lines, labels, loc='best', title=self.legend_title)

    def _get_ax_legend(self, ax):
        leg = ax.get_legend()
        other_ax = (getattr(ax, 'right_ax', None) or
                    getattr(ax, 'left_ax', None))
        other_leg = None
        if other_ax is not None:
            other_leg = other_ax.get_legend()
        if leg is None and other_leg is not None:
            leg = other_leg
            ax = other_ax
        return ax, leg

    def _maybe_convert_index(self, data):
        # tsplot converts automatically, but don't want to convert index
        # over and over for DataFrames
        from pandas.core.frame import DataFrame
        if (isinstance(data.index, DatetimeIndex) and
            isinstance(data, DataFrame)):
            freq = getattr(data.index, 'freq', None)

            if freq is None:
                freq = getattr(data.index, 'inferred_freq', None)
            if isinstance(freq, DateOffset):
                freq = freq.rule_code
            freq = get_period_alias(freq)

            if freq is None:
                ax = self._get_ax(0)
                freq = getattr(ax, 'freq', None)

            if freq is None:
                raise ValueError('Could not get frequency alias for plotting')

            data = DataFrame(data.values,
                             index=data.index.to_period(freq=freq),
                             columns=data.columns)
        return data

    def _post_plot_logic(self):
        df = self.data

        condition = (not self._use_dynamic_x()
                     and df.index.is_all_dates
                     and not self.subplots
                     or (self.subplots and self.sharex))

        index_name = self._get_index_name()

        rot = 30
        if self.rot is not None:
            rot = self.rot

        for ax in self.axes:
            if condition:
                format_date_labels(ax, rot=rot)
            elif self.rot is not None:
                for l in ax.get_xticklabels():
                    l.set_rotation(self.rot)

            if index_name is not None:
                ax.set_xlabel(index_name)

        if self.subplots and self.legend:
            for ax in self.axes:
                ax.legend(loc='best')


class BarPlot(MPLPlot):

    _default_rot = {'bar' : 90, 'barh' : 0}

    def __init__(self, data, **kwargs):
        self.stacked = kwargs.pop('stacked', False)
        self.ax_pos = np.arange(len(data)) + 0.25
        MPLPlot.__init__(self, data, **kwargs)

    def _args_adjust(self):
        if self.rot is None:
            self.rot = self._default_rot[self.kind]

    @property
    def bar_f(self):
        if self.kind == 'bar':
            def f(ax, x, y, w, start=None, **kwds):
                return ax.bar(x, y, w, bottom=start, **kwds)
        elif self.kind == 'barh':
            def f(ax, x, y, w, start=None, **kwds):
                return ax.barh(x, y, w, left=start, **kwds)
        else:
            raise NotImplementedError

        return f

    def _make_plot(self):
        colors = self.kwds.get('color', 'brgyk')
        rects = []
        labels = []

        ax = self._get_ax(0) #self.axes[0]

        bar_f = self.bar_f

        pos_prior = neg_prior = np.zeros(len(self.data))

        K = self.nseries

        for i, (label, y) in enumerate(self._iter_data()):
            label = com._stringify(label)
            kwds = self.kwds.copy()
            kwds['color'] = colors[i % len(colors)]

            if self.subplots:
                ax = self._get_ax(i) #self.axes[i]
                rect = bar_f(ax, self.ax_pos, y, 0.5, start=pos_prior, **kwds)
                ax.set_title(label)
            elif self.stacked:
                mask = y > 0
                start = np.where(mask, pos_prior, neg_prior)
                rect = bar_f(ax, self.ax_pos, y, 0.5, start=start,
                             label=label, **kwds)
                pos_prior = pos_prior + np.where(mask, y, 0)
                neg_prior = neg_prior + np.where(mask, 0, y)
            else:
                rect = bar_f(ax, self.ax_pos + i * 0.75 / K, y, 0.75 / K,
                             start=pos_prior, label=label, **kwds)
            rects.append(rect)
            labels.append(label)

        if self.legend and not self.subplots:
            patches =[r[0] for r in rects]
            self.axes[0].legend(patches, labels, loc='best',
                                title=self.legend_title)

    def _post_plot_logic(self):
        for ax in self.axes:
            str_index = [_stringify(key) for key in self.data.index]

            name = self._get_index_name()
            if self.kind == 'bar':
                ax.set_xlim([self.ax_pos[0] - 0.25, self.ax_pos[-1] + 1])
                ax.set_xticks(self.ax_pos + 0.375)
                ax.set_xticklabels(str_index, rotation=self.rot,
                                   fontsize=self.fontsize)
                ax.axhline(0, color='k', linestyle='--')
                if name is not None:
                    ax.set_xlabel(name)
            else:
                # horizontal bars
                ax.set_ylim([self.ax_pos[0] - 0.25, self.ax_pos[-1] + 1])
                ax.set_yticks(self.ax_pos + 0.375)
                ax.set_yticklabels(str_index, rotation=self.rot,
                                   fontsize=self.fontsize)
                ax.axvline(0, color='k', linestyle='--')
                if name is not None:
                    ax.set_ylabel(name)

        #if self.subplots and self.legend:
        #    self.axes[0].legend(loc='best')

class BoxPlot(MPLPlot):
    pass


class HistPlot(MPLPlot):
    pass


def plot_frame(frame=None, x=None, y=None, subplots=False, sharex=True,
               sharey=False, use_index=True, figsize=None, grid=False,
               legend=True, rot=None, ax=None, style=None, title=None, xlim=None,
               ylim=None, logy=False, xticks=None, yticks=None, kind='line',
               sort_columns=False, fontsize=None, secondary_y=False, **kwds):

    """
    Make line or bar plot of DataFrame's series with the index on the x-axis
    using matplotlib / pylab.

    Parameters
    ----------
    x : label or position, default None
    y : label or position, default None
        Allows plotting of one column versus another
    subplots : boolean, default False
        Make separate subplots for each time series
    sharex : boolean, default True
        In case subplots=True, share x axis
    sharey : boolean, default False
        In case subplots=True, share y axis
    use_index : boolean, default True
        Use index as ticks for x axis
    stacked : boolean, default False
        If True, create stacked bar plot. Only valid for DataFrame input
    sort_columns: boolean, default False
        Sort column names to determine plot ordering
    title : string
        Title to use for the plot
    grid : boolean, default True
        Axis grid lines
    legend : boolean, default True
        Place legend on axis subplots

    ax : matplotlib axis object, default None
    style : list or dict
        matplotlib line style per column
    kind : {'line', 'bar', 'barh'}
        bar : vertical bar plot
        barh : horizontal bar plot
    logy : boolean, default False
        For line plots, use log scaling on y axis
    xticks : sequence
        Values to use for the xticks
    yticks : sequence
        Values to use for the yticks
    xlim : 2-tuple/list
    ylim : 2-tuple/list
    rot : int, default None
        Rotation for ticks
    secondary_y : boolean or sequence, default False
        Whether to plot on the secondary y-axis
        If dict then can select which columns to plot on secondary y-axis
    kwds : keywords
        Options to pass to matplotlib plotting method

    Returns
    -------
    ax_or_axes : matplotlib.AxesSubplot or list of them
    """
    kind = _get_standard_kind(kind.lower().strip())
    if kind == 'line':
        klass = LinePlot
    elif kind in ('bar', 'barh'):
        klass = BarPlot
    elif kind == 'kde':
        klass = KdePlot
    else:
        raise ValueError('Invalid chart type given %s' % kind)

    if x is not None:
        if com.is_integer(x) and not frame.columns.holds_integer():
            x = frame.columns[x]
        frame = frame.set_index(x)

    if y is not None:
        if com.is_integer(y) and not frame.columns.holds_integer():
            y = frame.columns[y]
        return plot_series(frame[y], label=y, kind=kind, use_index=True,
                           rot=rot, xticks=xticks, yticks=yticks,
                           xlim=xlim, ylim=ylim, ax=ax, style=style,
                           grid=grid, logy=logy, secondary_y=secondary_y,
                           **kwds)

    plot_obj = klass(frame, kind=kind, subplots=subplots, rot=rot,
                     legend=legend, ax=ax, style=style, fontsize=fontsize,
                     use_index=use_index, sharex=sharex, sharey=sharey,
                     xticks=xticks, yticks=yticks, xlim=xlim, ylim=ylim,
                     title=title, grid=grid, figsize=figsize, logy=logy,
                     sort_columns=sort_columns, secondary_y=secondary_y,
                     **kwds)
    plot_obj.generate()
    plot_obj.draw()
    if subplots:
        return plot_obj.axes
    else:
        return plot_obj.axes[0]

def plot_series(series, label=None, kind='line', use_index=True, rot=None,
                xticks=None, yticks=None, xlim=None, ylim=None,
                ax=None, style=None, grid=None, logy=False, secondary_y=False,
                **kwds):
    """
    Plot the input series with the index on the x-axis using matplotlib

    Parameters
    ----------
    label : label argument to provide to plot
    kind : {'line', 'bar'}
    rot : int, default 30
        Rotation for tick labels
    use_index : boolean, default True
        Plot index as axis tick labels
    ax : matplotlib axis object
        If not passed, uses gca()
    style : string, default matplotlib default
        matplotlib line style to use
    ax : matplotlib axis object
        If not passed, uses gca()
    kind : {'line', 'bar', 'barh'}
        bar : vertical bar plot
        barh : horizontal bar plot
    logy : boolean, default False
        For line plots, use log scaling on y axis
    xticks : sequence
        Values to use for the xticks
    yticks : sequence
        Values to use for the yticks
    xlim : 2-tuple/list
    ylim : 2-tuple/list
    rot : int, default None
        Rotation for ticks
    kwds : keywords
        Options to pass to matplotlib plotting method

    Notes
    -----
    See matplotlib documentation online for more on this subject
    """
    kind = _get_standard_kind(kind.lower().strip())
    if kind == 'line':
        klass = LinePlot
    elif kind in ('bar', 'barh'):
        klass = BarPlot
    elif kind == 'kde':
        klass = KdePlot

    if ax is None:
        ax = _gca()
        if ax.get_yaxis().get_ticks_position().strip().lower() == 'right':
            fig = _gcf()
            axes = fig.get_axes()
            for i in range(len(axes))[::-1]:
                ax = axes[i]
                ypos = ax.get_yaxis().get_ticks_position().strip().lower()
                if ypos == 'left':
                    break

    # is there harm in this?
    if label is None:
        label = series.name

    plot_obj = klass(series, kind=kind, rot=rot, logy=logy,
                     ax=ax, use_index=use_index, style=style,
                     xticks=xticks, yticks=yticks, xlim=xlim, ylim=ylim,
                     legend=False, grid=grid, label=label,
                     secondary_y=secondary_y, **kwds)

    plot_obj.generate()
    plot_obj.draw()

    return plot_obj.ax

def boxplot(data, column=None, by=None, ax=None, fontsize=None,
            rot=0, grid=True, figsize=None, **kwds):
    """
    Make a box plot from DataFrame column optionally grouped by some columns or
    other inputs

    Parameters
    ----------
    data : DataFrame or Series
    column : column name or list of names, or vector
        Can be any valid input to groupby
    by : string or sequence
        Column in the DataFrame to group by
    fontsize : int or string
    rot : label rotation angle
    kwds : other plotting keyword arguments to be passed to matplotlib boxplot
           function

    Returns
    -------
    ax : matplotlib.axes.AxesSubplot
    """
    from pandas import Series, DataFrame
    if isinstance(data, Series):
        data = DataFrame({'x' : data})
        column = 'x'

    def plot_group(grouped, ax):
        keys, values = zip(*grouped)
        keys = [_stringify(x) for x in keys]
        values = [remove_na(v) for v in values]
        ax.boxplot(values, **kwds)
        if kwds.get('vert', 1):
            ax.set_xticklabels(keys, rotation=rot, fontsize=fontsize)
        else:
            ax.set_yticklabels(keys, rotation=rot, fontsize=fontsize)

    if column == None:
        columns = None
    else:
        if isinstance(column, (list, tuple)):
            columns = column
        else:
            columns = [column]

    if by is not None:
        if not isinstance(by, (list, tuple)):
            by = [by]

        fig, axes = _grouped_plot_by_column(plot_group, data, columns=columns,
                                            by=by, grid=grid, figsize=figsize,
                                            ax=ax)

        # Return axes in multiplot case, maybe revisit later # 985
        ret = axes
    else:
        if ax is None:
            ax = _gca()
        fig = ax.get_figure()
        data = data._get_numeric_data()
        if columns:
            cols = columns
        else:
            cols = data.columns
        keys = [_stringify(x) for x in cols]

        # Return boxplot dict in single plot case

        clean_values = [remove_na(x) for x in data[cols].values.T]
        bp = ax.boxplot(clean_values, **kwds)
        if kwds.get('vert', 1):
            ax.set_xticklabels(keys, rotation=rot, fontsize=fontsize)
        else:
            ax.set_yticklabels(keys, rotation=rot, fontsize=fontsize)
        ax.grid(grid)

        ret = bp

    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.9, wspace=0.2)
    return ret


def _stringify(x):
    if isinstance(x, tuple):
        return '|'.join(com._stringify(y) for y in x)
    else:
        return com._stringify(x)


def format_date_labels(ax, rot):
    # mini version of autofmt_xdate
    try:
        for label in ax.get_xticklabels():
            label.set_ha('right')
            label.set_rotation(rot)
        fig = ax.get_figure()
        fig.subplots_adjust(bottom=0.2)
    except Exception: # pragma: no cover
        pass


def scatter_plot(data, x, y, by=None, ax=None, figsize=None, grid=False):
    """

    Returns
    -------
    fig : matplotlib.Figure
    """
    import matplotlib.pyplot as plt

    def plot_group(group, ax):
        xvals = group[x].values
        yvals = group[y].values
        ax.scatter(xvals, yvals)
        ax.grid(grid)

    if by is not None:
        fig = _grouped_plot(plot_group, data, by=by, figsize=figsize, ax=ax)
    else:
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            fig = ax.get_figure()
        plot_group(data, ax)
        ax.set_ylabel(com._stringify(y))
        ax.set_xlabel(com._stringify(x))

        ax.grid(grid)

    return fig


def hist_frame(data, grid=True, xlabelsize=None, xrot=None,
               ylabelsize=None, yrot=None, ax=None,
               sharex=False, sharey=False, **kwds):
    """
    Draw Histogram the DataFrame's series using matplotlib / pylab.

    Parameters
    ----------
    grid : boolean, default True
        Whether to show axis grid lines
    xlabelsize : int, default None
        If specified changes the x-axis label size
    xrot : float, default None
        rotation of x axis labels
    ylabelsize : int, default None
        If specified changes the y-axis label size
    yrot : float, default None
        rotation of y axis labels
    ax : matplotlib axes object, default None
    sharex : bool, if True, the X axis will be shared amongst all subplots.
    sharey : bool, if True, the Y axis will be shared amongst all subplots.
    kwds : other plotting keyword arguments
        To be passed to hist function
    """
    import matplotlib.pyplot as plt
    n = len(data.columns)
    rows, cols = 1, 1
    while rows * cols < n:
        if cols > rows:
            rows += 1
        else:
            cols += 1
    _, axes = _subplots(nrows=rows, ncols=cols, ax=ax, squeeze=False,
                        sharex=sharex, sharey=sharey)

    for i, col in enumerate(com._try_sort(data.columns)):
        ax = axes[i / cols][i % cols]
        ax.xaxis.set_visible(True)
        ax.yaxis.set_visible(True)
        ax.hist(data[col].dropna().values, **kwds)
        ax.set_title(col)
        ax.grid(grid)

        if xlabelsize is not None:
            plt.setp(ax.get_xticklabels(), fontsize=xlabelsize)
        if xrot is not None:
            plt.setp(ax.get_xticklabels(), rotation=xrot)
        if ylabelsize is not None:
            plt.setp(ax.get_yticklabels(), fontsize=ylabelsize)
        if yrot is not None:
            plt.setp(ax.get_yticklabels(), rotation=yrot)

    for j in range(i + 1, rows * cols):
        ax = axes[j / cols, j % cols]
        ax.set_visible(False)

    ax.get_figure().subplots_adjust(wspace=0.3, hspace=0.3)

    return axes

def hist_series(self, ax=None, grid=True, xlabelsize=None, xrot=None,
                ylabelsize=None, yrot=None, **kwds):
    """
    Draw histogram of the input series using matplotlib

    Parameters
    ----------
    ax : matplotlib axis object
        If not passed, uses gca()
    grid : boolean, default True
        Whether to show axis grid lines
    xlabelsize : int, default None
        If specified changes the x-axis label size
    xrot : float, default None
        rotation of x axis labels
    ylabelsize : int, default None
        If specified changes the y-axis label size
    yrot : float, default None
        rotation of y axis labels
    kwds : keywords
        To be passed to the actual plotting function

    Notes
    -----
    See matplotlib documentation online for more on this

    """
    import matplotlib.pyplot as plt

    if ax is None:
        ax = plt.gca()

    values = self.dropna().values

    ax.hist(values, **kwds)
    ax.grid(grid)

    if xlabelsize is not None:
        plt.setp(ax.get_xticklabels(), fontsize=xlabelsize)
    if xrot is not None:
        plt.setp(ax.get_xticklabels(), rotation=xrot)
    if ylabelsize is not None:
        plt.setp(ax.get_yticklabels(), fontsize=ylabelsize)
    if yrot is not None:
        plt.setp(ax.get_yticklabels(), rotation=yrot)

    return ax

def boxplot_frame_groupby(grouped, subplots=True, column=None, fontsize=None,
                          rot=0, grid=True, figsize=None, **kwds):
    """
    Make box plots from DataFrameGroupBy data.

    Parameters
    ----------
    subplots :
        * ``False`` - no subplots will be used
        * ``True`` - create a subplot for each group
    column : column name or list of names, or vector
        Can be any valid input to groupby
    fontsize : int or string
    rot : label rotation angle
    kwds : other plotting keyword arguments to be passed to matplotlib boxplot
           function

    Returns
    -------
    dict of key/value = group key/DataFrame.boxplot return value
    or DataFrame.boxplot return value in case subplots=figures=False

    Examples
    --------
    >>> import pandas
    >>> import numpy as np
    >>> import itertools
    >>>
    >>> tuples = [t for t in itertools.product(range(1000), range(4))]
    >>> index = pandas.MultiIndex.from_tuples(tuples, names=['lvl0', 'lvl1'])
    >>> data = np.random.randn(len(index),4)
    >>> df = pandas.DataFrame(data, columns=list('ABCD'), index=index)
    >>>
    >>> grouped = df.groupby(level='lvl1')
    >>> boxplot_frame_groupby(grouped)
    >>>
    >>> grouped = df.unstack(level='lvl1').groupby(level=0, axis=1)
    >>> boxplot_frame_groupby(grouped, subplots=False)
    """
    if subplots is True:
        nrows, ncols = _get_layout(len(grouped))
        _, axes = _subplots(nrows=nrows, ncols=ncols, squeeze=False,
                            sharex=False, sharey=True)
        axes = axes.reshape(-1) if len(grouped) > 1 else axes

        ret = {}
        for (key, group), ax in zip(grouped, axes):
            d = group.boxplot(ax=ax, column=column, fontsize=fontsize,
                              rot=rot, grid=grid, figsize=figsize, **kwds)
            ax.set_title(_stringify(key))
            ret[key] = d
    else:
        from pandas.tools.merge import concat
        keys, frames = zip(*grouped)
        if grouped.axis == 0:
            df = concat(frames, keys=keys, axis=1)
        else:
            if len(frames) > 1:
                df = frames[0].join(frames[1::])
            else:
                df = frames[0]
        ret = df.boxplot(column=column, fontsize=fontsize, rot=rot,
                         grid=grid, figsize=figsize, **kwds)
    return ret

def _grouped_plot(plotf, data, column=None, by=None, numeric_only=True,
                  figsize=None, sharex=True, sharey=True, layout=None,
                  rot=0, ax=None):
    from pandas.core.frame import DataFrame
    import matplotlib.pyplot as plt

    # allow to specify mpl default with 'default'
    if figsize is None or figsize == 'default':
        figsize = (10, 5)               # our default

    grouped = data.groupby(by)
    if column is not None:
        grouped = grouped[column]

    ngroups = len(grouped)

    nrows, ncols = layout or _get_layout(ngroups)

    if figsize is None:
        # our favorite default beating matplotlib's idea of the
        # default size
        figsize = (10, 5)
    fig, axes = _subplots(nrows=nrows, ncols=ncols, figsize=figsize,
                          sharex=sharex, sharey=sharey, ax=ax)

    if isinstance(axes, plt.Axes):
        ravel_axes = [axes]
    else:
        ravel_axes = []
        for row in axes:
            if isinstance(row, plt.Axes):
                ravel_axes.append(row)
            else:
                ravel_axes.extend(row)

    for i, (key, group) in enumerate(grouped):
        ax = ravel_axes[i]
        if numeric_only and isinstance(group, DataFrame):
            group = group._get_numeric_data()
        plotf(group, ax)
        ax.set_title(com._stringify(key))

    return fig, axes

def _grouped_plot_by_column(plotf, data, columns=None, by=None,
                            numeric_only=True, grid=False,
                            figsize=None, ax=None):
    import matplotlib.pyplot as plt

    grouped = data.groupby(by)
    if columns is None:
        columns = data._get_numeric_data().columns - by
    ngroups = len(columns)

    nrows, ncols = _get_layout(ngroups)
    fig, axes = _subplots(nrows=nrows, ncols=ncols,
                          sharex=True, sharey=True,
                          figsize=figsize, ax=ax)

    if isinstance(axes, plt.Axes):
        ravel_axes = [axes]
    else:
        ravel_axes = []
        for row in axes:
            if isinstance(row, plt.Axes):
                ravel_axes.append(row)
            else:
                ravel_axes.extend(row)

    for i, col in enumerate(columns):
        ax = ravel_axes[i]
        gp_col = grouped[col]
        plotf(gp_col, ax)
        ax.set_title(col)
        ax.set_xlabel(com._stringify(by))
        ax.grid(grid)

    byline = by[0] if len(by) == 1 else by
    fig.suptitle('Boxplot grouped by %s' % byline)

    return fig, axes

def _get_layout(nplots):
    if nplots == 1:
        return (1, 1)
    elif nplots == 2:
        return (1, 2)
    elif nplots < 4:
        return (2, 2)

    k = 1
    while k ** 2 < nplots:
        k += 1

    if (k - 1) * k >= nplots:
        return k, (k - 1)
    else:
        return k, k

# copied from matplotlib/pyplot.py for compatibility with matplotlib < 1.0

def _subplots(nrows=1, ncols=1, sharex=False, sharey=False, squeeze=True,
              subplot_kw=None, ax=None, secondary_y=False, data=None,
              **fig_kw):
    """Create a figure with a set of subplots already made.

    This utility wrapper makes it convenient to create common layouts of
    subplots, including the enclosing figure object, in a single call.

    Keyword arguments:

    nrows : int
      Number of rows of the subplot grid.  Defaults to 1.

    ncols : int
      Number of columns of the subplot grid.  Defaults to 1.

    sharex : bool
      If True, the X axis will be shared amongst all subplots.

    sharey : bool
      If True, the Y axis will be shared amongst all subplots.

    squeeze : bool

      If True, extra dimensions are squeezed out from the returned axis object:
        - if only one subplot is constructed (nrows=ncols=1), the resulting
        single Axis object is returned as a scalar.
        - for Nx1 or 1xN subplots, the returned object is a 1-d numpy object
        array of Axis objects are returned as numpy 1-d arrays.
        - for NxM subplots with N>1 and M>1 are returned as a 2d array.

      If False, no squeezing at all is done: the returned axis object is always
      a 2-d array containing Axis instances, even if it ends up being 1x1.

    subplot_kw : dict
      Dict with keywords passed to the add_subplot() call used to create each
      subplots.

    fig_kw : dict
      Dict with keywords passed to the figure() call.  Note that all keywords
      not recognized above will be automatically included here.

    ax : Matplotlib axis object, default None

    secondary_y : boolean or sequence of ints, default False
        If True then y-axis will be on the right

    Returns:

    fig, ax : tuple
      - fig is the Matplotlib Figure object
      - ax can be either a single axis object or an array of axis objects if
      more than one subplot was created.  The dimensions of the resulting array
      can be controlled with the squeeze keyword, see above.

    **Examples:**

    x = np.linspace(0, 2*np.pi, 400)
    y = np.sin(x**2)

    # Just a figure and one subplot
    f, ax = plt.subplots()
    ax.plot(x, y)
    ax.set_title('Simple plot')

    # Two subplots, unpack the output array immediately
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    ax1.plot(x, y)
    ax1.set_title('Sharing Y axis')
    ax2.scatter(x, y)

    # Four polar axes
    plt.subplots(2, 2, subplot_kw=dict(polar=True))
    """
    import matplotlib.pyplot as plt
    from pandas.core.frame import DataFrame

    if subplot_kw is None:
        subplot_kw = {}

    if ax is None:
        fig = plt.figure(**fig_kw)
    else:
        fig = ax.get_figure()
        fig.clear()

    # Create empty object array to hold all axes.  It's easiest to make it 1-d
    # so we can just append subplots upon creation, and then
    nplots = nrows*ncols
    axarr = np.empty(nplots, dtype=object)

    def on_right(i):
        if isinstance(secondary_y, bool):
            return secondary_y
        if isinstance(data, DataFrame):
            return data.columns[i] in secondary_y

    # Create first subplot separately, so we can share it if requested
    ax0 = fig.add_subplot(nrows, ncols, 1, **subplot_kw)
    if on_right(0):
        orig_ax = ax0
        ax0 = ax0.twinx()
        orig_ax.get_yaxis().set_visible(False)
        orig_ax.right_ax = ax0
        ax0.left_ax = orig_ax

    if sharex:
        subplot_kw['sharex'] = ax0
    if sharey:
        subplot_kw['sharey'] = ax0
    axarr[0] = ax0

    # Note off-by-one counting because add_subplot uses the MATLAB 1-based
    # convention.
    for i in range(1, nplots):
        ax = fig.add_subplot(nrows, ncols, i+1, **subplot_kw)
        if on_right(i):
            orig_ax = ax
            ax = ax.twinx()
            orig_ax.get_yaxis().set_visible(False)
        axarr[i] = ax

    if nplots > 1:
        if sharex and nrows > 1:
            for i, ax in enumerate(axarr):
                if np.ceil(float(i + 1) / ncols) < nrows: # only last row
                    [label.set_visible(False) for label in ax.get_xticklabels()]
        if sharey and ncols > 1:
            for i, ax in enumerate(axarr):
                if (i % ncols) != 0: # only first column
                    [label.set_visible(False) for label in ax.get_yticklabels()]

    if squeeze:
        # Reshape the array to have the final desired dimension (nrow,ncol),
        # though discarding unneeded dimensions that equal 1.  If we only have
        # one subplot, just return it instead of a 1-element array.
        if nplots==1:
            axes = axarr[0]
        else:
            axes = axarr.reshape(nrows, ncols).squeeze()
    else:
        # returned axis array will be always 2-d, even if nrows=ncols=1
        axes = axarr.reshape(nrows, ncols)

    return fig, axes

if __name__ == '__main__':
    # import pandas.rpy.common as com
    # sales = com.load_data('sanfrancisco.home.sales', package='nutshell')
    # top10 = sales['zip'].value_counts()[:10].index
    # sales2 = sales[sales.zip.isin(top10)]
    # _ = scatter_plot(sales2, 'squarefeet', 'price', by='zip')

    # plt.show()

    import matplotlib.pyplot as plt

    import pandas.tools.plotting as plots
    import pandas.core.frame as fr
    reload(plots)
    reload(fr)
    from pandas.core.frame import DataFrame

    data = DataFrame([[3, 6, -5], [4, 8, 2], [4, 9, -6],
                      [4, 9, -3], [2, 5, -1]],
                     columns=['A', 'B', 'C'])
    data.plot(kind='barh', stacked=True)

    plt.show()
