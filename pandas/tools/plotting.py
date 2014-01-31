# being a bit too dynamic
# pylint: disable=E1101
import datetime
import warnings
import re
from contextlib import contextmanager
from distutils.version import LooseVersion

import numpy as np

from pandas.util.decorators import cache_readonly
import pandas.core.common as com
from pandas.core.index import MultiIndex
from pandas.core.series import Series, remove_na
from pandas.tseries.index import DatetimeIndex
from pandas.tseries.period import PeriodIndex, Period
from pandas.tseries.frequencies import get_period_alias, get_base_alias
from pandas.tseries.offsets import DateOffset
from pandas.compat import range, lrange, lmap, map, zip
import pandas.compat as compat

try:  # mpl optional
    import pandas.tseries.converter as conv
    conv.register()  # needs to override so set_xlim works with str/number
except ImportError:
    pass

# Extracted from https://gist.github.com/huyng/816622
# this is the rcParams set when setting display.with_mpl_style
# to True.
mpl_stylesheet = {
    'axes.axisbelow': True,
     'axes.color_cycle': ['#348ABD',
      '#7A68A6',
      '#A60628',
      '#467821',
      '#CF4457',
      '#188487',
      '#E24A33'],
     'axes.edgecolor': '#bcbcbc',
     'axes.facecolor': '#eeeeee',
     'axes.grid': True,
     'axes.labelcolor': '#555555',
     'axes.labelsize': 'large',
     'axes.linewidth': 1.0,
     'axes.titlesize': 'x-large',
     'figure.edgecolor': 'white',
     'figure.facecolor': 'white',
     'figure.figsize': (6.0, 4.0),
     'figure.subplot.hspace': 0.5,
     'font.family': 'monospace',
     'font.monospace': ['Andale Mono',
      'Nimbus Mono L',
      'Courier New',
      'Courier',
      'Fixed',
      'Terminal',
      'monospace'],
     'font.size': 10,
     'interactive': True,
     'keymap.all_axes': ['a'],
     'keymap.back': ['left', 'c', 'backspace'],
     'keymap.forward': ['right', 'v'],
     'keymap.fullscreen': ['f'],
     'keymap.grid': ['g'],
     'keymap.home': ['h', 'r', 'home'],
     'keymap.pan': ['p'],
     'keymap.save': ['s'],
     'keymap.xscale': ['L', 'k'],
     'keymap.yscale': ['l'],
     'keymap.zoom': ['o'],
     'legend.fancybox': True,
     'lines.antialiased': True,
     'lines.linewidth': 1.0,
     'patch.antialiased': True,
     'patch.edgecolor': '#EEEEEE',
     'patch.facecolor': '#348ABD',
     'patch.linewidth': 0.5,
     'toolbar': 'toolbar2',
     'xtick.color': '#555555',
     'xtick.direction': 'in',
     'xtick.major.pad': 6.0,
     'xtick.major.size': 0.0,
     'xtick.minor.pad': 6.0,
     'xtick.minor.size': 0.0,
     'ytick.color': '#555555',
     'ytick.direction': 'in',
     'ytick.major.pad': 6.0,
     'ytick.major.size': 0.0,
     'ytick.minor.pad': 6.0,
     'ytick.minor.size': 0.0
}

def _get_standard_kind(kind):
    return {'density': 'kde'}.get(kind, kind)

def _get_standard_colors(num_colors=None, colormap=None, color_type='default',
                         color=None):
    import matplotlib.pyplot as plt

    if color is None and colormap is not None:
        if isinstance(colormap, compat.string_types):
            import matplotlib.cm as cm
            cmap = colormap
            colormap = cm.get_cmap(colormap)
            if colormap is None:
                raise ValueError("Colormap {0} is not recognized".format(cmap))
        colors = lmap(colormap, np.linspace(0, 1, num=num_colors))
    elif color is not None:
        if colormap is not None:
            warnings.warn("'color' and 'colormap' cannot be used "
                          "simultaneously. Using 'color'")
        colors = color
    else:
        if color_type == 'default':
            colors = plt.rcParams.get('axes.color_cycle', list('bgrcmyk'))
            if isinstance(colors, compat.string_types):
                colors = list(colors)
        elif color_type == 'random':
            import random
            def random_color(column):
                random.seed(column)
                return [random.random() for _ in range(3)]

            colors = lmap(random_color, lrange(num_colors))
        else:
            raise NotImplementedError

    if len(colors) != num_colors:
        multiple = num_colors//len(colors) - 1
        mod = num_colors % len(colors)

        colors += multiple * colors
        colors += colors[:mod]

    return colors

class _Options(dict):
    """
    Stores pandas plotting options.
    Allows for parameter aliasing so you can just use parameter names that are
    the same as the plot function parameters, but is stored in a canonical
    format that makes it easy to breakdown into groups later
    """

    # alias so the names are same as plotting method parameter names
    _ALIASES = {'x_compat': 'xaxis.compat'}
    _DEFAULT_KEYS = ['xaxis.compat']

    def __init__(self):
        self['xaxis.compat'] = False

    def __getitem__(self, key):
        key = self._get_canonical_key(key)
        if key not in self:
            raise ValueError('%s is not a valid pandas plotting option' % key)
        return super(_Options, self).__getitem__(key)

    def __setitem__(self, key, value):
        key = self._get_canonical_key(key)
        return super(_Options, self).__setitem__(key, value)

    def __delitem__(self, key):
        key = self._get_canonical_key(key)
        if key in self._DEFAULT_KEYS:
            raise ValueError('Cannot remove default parameter %s' % key)
        return super(_Options, self).__delitem__(key)

    def __contains__(self, key):
        key = self._get_canonical_key(key)
        return super(_Options, self).__contains__(key)

    def reset(self):
        """
        Reset the option store to its initial state

        Returns
        -------
        None
        """
        self.__init__()

    def _get_canonical_key(self, key):
        return self._ALIASES.get(key, key)

    @contextmanager
    def use(self, key, value):
        """
        Temporarily set a parameter value using the with statement.
        Aliasing allowed.
        """
        old_value = self[key]
        try:
            self[key] = value
            yield self
        finally:
            self[key] = old_value


plot_params = _Options()


def scatter_matrix(frame, alpha=0.5, figsize=None, ax=None, grid=False,
                   diagonal='hist', marker='.', density_kwds=None,
                   hist_kwds=None, range_padding=0.05, **kwds):
    """
    Draw a matrix of scatter plots.

    Parameters
    ----------
    frame : DataFrame
    alpha : float, optional
        amount of transparency applied
    figsize : (float,float), optional
        a tuple (width, height) in inches
    ax : Matplotlib axis object, optional
    grid : bool, optional
        setting this to True will show the grid
    diagonal : {'hist', 'kde'}
        pick between 'kde' and 'hist' for
        either Kernel Density Estimation or Histogram
        plot in the diagonal
    marker : str, optional
        Matplotlib marker type, default '.'    
    hist_kwds : other plotting keyword arguments
        To be passed to hist function
    density_kwds : other plotting keyword arguments
        To be passed to kernel density estimate plot
    range_padding : float, optional
        relative extension of axis range in x and y
        with respect to (x_max - x_min) or (y_max - y_min),
        default 0.05
    kwds : other plotting keyword arguments
        To be passed to scatter function

    Examples
    --------
    >>> df = DataFrame(np.random.randn(1000, 4), columns=['A','B','C','D'])
    >>> scatter_matrix(df, alpha=0.2)
    """
    import matplotlib.pyplot as plt
    from matplotlib.artist import setp

    df = frame._get_numeric_data()
    n = df.columns.size
    fig, axes = _subplots(nrows=n, ncols=n, figsize=figsize, ax=ax,
                          squeeze=False)

    # no gaps between subplots
    fig.subplots_adjust(wspace=0, hspace=0)

    mask = com.notnull(df)

    marker = _get_marker_compat(marker)

    hist_kwds = hist_kwds or {}
    density_kwds = density_kwds or {}

    # workaround because `c='b'` is hardcoded in matplotlibs scatter method
    kwds.setdefault('c', plt.rcParams['patch.facecolor'])

    boundaries_list = []
    for a in df.columns:
        values = df[a].values[mask[a].values]
        rmin_, rmax_ = np.min(values), np.max(values)
        rdelta_ext = (rmax_ - rmin_) * range_padding / 2.
        boundaries_list.append((rmin_ - rdelta_ext, rmax_+ rdelta_ext))

    for i, a in zip(lrange(n), df.columns):
        for j, b in zip(lrange(n), df.columns):
            ax = axes[i, j]

            if i == j:
                values = df[a].values[mask[a].values]

                # Deal with the diagonal by drawing a histogram there.
                if diagonal == 'hist':
                    ax.hist(values, **hist_kwds)

                elif diagonal in ('kde', 'density'):
                    from scipy.stats import gaussian_kde
                    y = values
                    gkde = gaussian_kde(y)
                    ind = np.linspace(y.min(), y.max(), 1000)
                    ax.plot(ind, gkde.evaluate(ind), **density_kwds)

                ax.set_xlim(boundaries_list[i])

            else:
                common = (mask[a] & mask[b]).values

                ax.scatter(df[b][common], df[a][common],
                           marker=marker, alpha=alpha, **kwds)

                ax.set_xlim(boundaries_list[j])
                ax.set_ylim(boundaries_list[i])

            ax.set_xlabel('')
            ax.set_ylabel('')

            _label_axis(ax, kind='x', label=b, position='bottom', rotate=True)

            _label_axis(ax, kind='y', label=a, position='left')

            if j!= 0:
                ax.yaxis.set_visible(False)
            if i != n-1:
                ax.xaxis.set_visible(False)

    for ax in axes.flat:
        setp(ax.get_xticklabels(), fontsize=8)
        setp(ax.get_yticklabels(), fontsize=8)

    return axes

def _label_axis(ax, kind='x', label='', position='top',
    ticks=True, rotate=False):

    from matplotlib.artist import setp
    if kind == 'x':
        ax.set_xlabel(label, visible=True)
        ax.xaxis.set_visible(True)
        ax.xaxis.set_ticks_position(position)
        ax.xaxis.set_label_position(position)
        if rotate:
            setp(ax.get_xticklabels(), rotation=90)
    elif kind == 'y':
        ax.yaxis.set_visible(True)
        ax.set_ylabel(label, visible=True)
        # ax.set_ylabel(a)
        ax.yaxis.set_ticks_position(position)
        ax.yaxis.set_label_position(position)
    return





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


def radviz(frame, class_column, ax=None, colormap=None, **kwds):
    """RadViz - a multivariate data visualization algorithm

    Parameters:
    -----------
    frame: DataFrame object
    class_column: Column name that contains information about class membership
    ax: Matplotlib axis object, optional
    colormap : str or matplotlib colormap object, default None
        Colormap to select colors from. If string, load colormap with that name
        from matplotlib.
    kwds: Matplotlib scatter method keyword arguments, optional

    Returns:
    --------
    ax: Matplotlib axis object
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches

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

    colors = _get_standard_colors(num_colors=len(classes), colormap=colormap,
                                  color_type='random', color=kwds.get('color'))

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

    for i, class_ in enumerate(classes):
        ax.scatter(to_plot[class_][0], to_plot[class_][1], color=colors[i],
                   label=com.pprint_thing(class_), **kwds)
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


def andrews_curves(data, class_column, ax=None, samples=200, colormap=None,
                   **kwds):
    """
    Parameters:
    -----------
    data : DataFrame
        Data to be plotted, preferably normalized to (0.0, 1.0)
    class_column : Name of the column containing class names
    ax : matplotlib axes object, default None
    samples : Number of points to plot in each curve
    colormap : str or matplotlib colormap object, default None
        Colormap to select colors from. If string, load colormap with that name
        from matplotlib.
    kwds : Optional plotting arguments to be passed to matplotlib

    Returns:
    --------
    ax: Matplotlib axis object

    """
    from math import sqrt, pi, sin, cos
    import matplotlib.pyplot as plt

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

    n = len(data)
    class_col = data[class_column]
    uniq_class = class_col.drop_duplicates()
    columns = [data[col] for col in data.columns if (col != class_column)]
    x = [-pi + 2.0 * pi * (t / float(samples)) for t in range(samples)]
    used_legends = set([])

    colors = _get_standard_colors(num_colors=len(uniq_class), colormap=colormap,
                                  color_type='random', color=kwds.get('color'))
    col_dict = dict([(klass, col) for klass, col in zip(uniq_class, colors)])
    if ax is None:
        ax = plt.gca(xlim=(-pi, pi))
    for i in range(n):
        row = [columns[c][i] for c in range(len(columns))]
        f = function(row)
        y = [f(t) for t in x]
        label = None
        if com.pprint_thing(class_col[i]) not in used_legends:
            label = com.pprint_thing(class_col[i])
            used_legends.add(label)
            ax.plot(x, y, color=col_dict[class_col[i]], label=label, **kwds)
        else:
            ax.plot(x, y, color=col_dict[class_col[i]], **kwds)

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
    kwds: optional keyword arguments for plotting commands, must be accepted
        by both hist and plot

    Returns:
    --------
    fig: matplotlib figure
    """
    import random
    import matplotlib.pyplot as plt

    # random.sample(ndarray, int) fails on python 3.3, sigh
    data = list(series.values)
    samplings = [random.sample(data, size) for _ in range(samples)]

    means = np.array([np.mean(sampling) for sampling in samplings])
    medians = np.array([np.median(sampling) for sampling in samplings])
    midranges = np.array([(min(sampling) + max(sampling)) * 0.5
                          for sampling in samplings])
    if fig is None:
        fig = plt.figure()
    x = lrange(samples)
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


def parallel_coordinates(data, class_column, cols=None, ax=None, colors=None,
                         use_columns=False, xticks=None, colormap=None, **kwds):
    """Parallel coordinates plotting.

    Parameters
    ----------
    data: DataFrame
        A DataFrame containing data to be plotted
    class_column: str
        Column name containing class names
    cols: list, optional
        A list of column names to use
    ax: matplotlib.axis, optional
        matplotlib axis object
    colors: list or tuple, optional
        Colors to use for the different classes
    use_columns: bool, optional
        If true, columns will be used as xticks
    xticks: list or tuple, optional
        A list of values to use for xticks
    colormap: str or matplotlib colormap, default None
        Colormap to use for line colors.
    kwds: list, optional
        A list of keywords for matplotlib plot method

    Returns
    -------
    ax: matplotlib axis object

    Examples
    --------
    >>> from pandas import read_csv
    >>> from pandas.tools.plotting import parallel_coordinates
    >>> from matplotlib import pyplot as plt
    >>> df = read_csv('https://raw.github.com/pydata/pandas/master/pandas/tests/data/iris.csv')
    >>> parallel_coordinates(df, 'Name', colors=('#556270', '#4ECDC4', '#C7F464'))
    >>> plt.show()
    """
    import matplotlib.pyplot as plt


    n = len(data)
    classes = set(data[class_column])
    class_col = data[class_column]

    if cols is None:
        df = data.drop(class_column, axis=1)
    else:
        df = data[cols]

    used_legends = set([])

    ncols = len(df.columns)

    # determine values to use for xticks
    if use_columns is True:
        if not np.all(np.isreal(list(df.columns))):
            raise ValueError('Columns must be numeric to be used as xticks')
        x = df.columns
    elif xticks is not None:
        if not np.all(np.isreal(xticks)):
            raise ValueError('xticks specified must be numeric')
        elif len(xticks) != ncols:
            raise ValueError('Length of xticks must match number of columns')
        x = xticks
    else:
        x = lrange(ncols)

    if ax is None:
        ax = plt.gca()

    color_values = _get_standard_colors(num_colors=len(classes),
                                        colormap=colormap, color_type='random',
                                        color=colors)

    colors = dict(zip(classes, color_values))

    for i in range(n):
        row = df.irow(i).values
        y = row
        kls = class_col.iget_value(i)
        if com.pprint_thing(kls) not in used_legends:
            label = com.pprint_thing(kls)
            used_legends.add(label)
            ax.plot(x, y, color=colors[kls],
                    label=label, **kwds)
        else:
            ax.plot(x, y, color=colors[kls], **kwds)

    for i in x:
        ax.axvline(i, linewidth=1, color='black')

    ax.set_xticks(x)
    ax.set_xticklabels(df.columns)
    ax.set_xlim(x[0], x[-1])
    ax.legend(loc='upper right')
    ax.grid()
    return ax


def lag_plot(series, lag=1, ax=None, **kwds):
    """Lag plot for time series.

    Parameters:
    -----------
    series: Time series
    lag: lag of the scatter plot, default 1
    ax: Matplotlib axis object, optional
    kwds: Matplotlib scatter method keyword arguments, optional

    Returns:
    --------
    ax: Matplotlib axis object
    """
    import matplotlib.pyplot as plt

    # workaround because `c='b'` is hardcoded in matplotlibs scatter method
    kwds.setdefault('c', plt.rcParams['patch.facecolor'])

    data = series.values
    y1 = data[:-lag]
    y2 = data[lag:]
    if ax is None:
        ax = plt.gca()
    ax.set_xlabel("y(t)")
    ax.set_ylabel("y(t + %s)" % lag)
    ax.scatter(y1, y2, **kwds)
    return ax


def autocorrelation_plot(series, ax=None, **kwds):
    """Autocorrelation plot for time series.

    Parameters:
    -----------
    series: Time series
    ax: Matplotlib axis object, optional
    kwds : keywords
        Options to pass to matplotlib plotting method

    Returns:
    -----------
    ax: Matplotlib axis object
    """
    import matplotlib.pyplot as plt
    n = len(series)
    data = np.asarray(series)
    if ax is None:
        ax = plt.gca(xlim=(1, n), ylim=(-1.0, 1.0))
    mean = np.mean(data)
    c0 = np.sum((data - mean) ** 2) / float(n)

    def r(h):
        return ((data[:n - h] - mean) * (data[h:] - mean)).sum() / float(n) / c0
    x = np.arange(n) + 1
    y = lmap(r, x)
    z95 = 1.959963984540054
    z99 = 2.5758293035489004
    ax.axhline(y=z99 / np.sqrt(n), linestyle='--', color='grey')
    ax.axhline(y=z95 / np.sqrt(n), color='grey')
    ax.axhline(y=0.0, color='black')
    ax.axhline(y=-z95 / np.sqrt(n), color='grey')
    ax.axhline(y=-z99 / np.sqrt(n), linestyle='--', color='grey')
    ax.set_xlabel("Lag")
    ax.set_ylabel("Autocorrelation")
    ax.plot(x, y, **kwds)
    if 'label' in kwds:
        ax.legend()
    ax.grid()
    return ax


def grouped_hist(data, column=None, by=None, ax=None, bins=50, figsize=None,
                 layout=None, sharex=False, sharey=False, rot=90, grid=True,
                 **kwargs):
    """
    Grouped histogram

    Parameters
    ----------
    data: Series/DataFrame
    column: object, optional
    by: object, optional
    ax: axes, optional
    bins: int, default 50
    figsize: tuple, optional
    layout: optional
    sharex: boolean, default False
    sharey: boolean, default False
    rot: int, default 90
    grid: bool, default True
    kwargs: dict, keyword arguments passed to matplotlib.Axes.hist

    Returns
    -------
    axes: collection of Matplotlib Axes
    """
    def plot_group(group, ax):
        ax.hist(group.dropna().values, bins=bins, **kwargs)

    fig, axes = _grouped_plot(plot_group, data, column=column,
                              by=by, sharex=sharex, sharey=sharey,
                              figsize=figsize, layout=layout, rot=rot)
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.9,
                        hspace=0.5, wspace=0.3)
    return axes


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
                 secondary_y=False, colormap=None, **kwds):

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

        self.colormap = colormap

        self.kwds = kwds

        self._validate_color_args()

    def _validate_color_args(self):
        from pandas import DataFrame
        if 'color' not in self.kwds and 'colors' in self.kwds:
            warnings.warn(("'colors' is being deprecated. Please use 'color'"
                           "instead of 'colors'"))
            colors = self.kwds.pop('colors')
            self.kwds['color'] = colors

        if ('color' in self.kwds and
            (isinstance(self.data, Series) or
             isinstance(self.data, DataFrame) and len(self.data.columns) == 1)):
            # support series.plot(color='green')
            self.kwds['color'] = [self.kwds['color']]

        if ('color' in self.kwds or 'colors' in self.kwds) and \
                self.colormap is not None:
            warnings.warn("'color' and 'colormap' cannot be used "
                          "simultaneously. Using 'color'")

        if 'color' in self.kwds and self.style is not None:
            # need only a single match
            if re.match('^[a-z]+?', self.style) is not None:
                raise ValueError("Cannot pass 'style' string with a color "
                                 "symbol and 'color' keyword argument. Please"
                                 " use one or the other or pass 'style' "
                                 "without a color symbol")

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
            new_ax._get_lines.color_cycle = orig_ax._get_lines.color_cycle

            orig_ax.right_ax, new_ax.left_ax = new_ax, orig_ax
            new_ax.right_ax = new_ax

            if len(orig_ax.get_lines()) == 0:  # no data on left y
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
                if self.figsize is not None:
                    fig.set_size_inches(self.figsize)
                ax = self._maybe_right_yaxis(self.ax)

            axes = [ax]

        self.fig = fig
        self.axes = axes

    def _get_layout(self):
        return (len(self.data.columns), 1)

    def _compute_plot_data(self):
        numeric_data = self.data.convert_objects()._get_numeric_data()

        try:
            is_empty = numeric_data.empty
        except AttributeError:
            is_empty = not len(numeric_data)

        # no empty frames or series allowed
        if is_empty:
            raise TypeError('Empty {0!r}: no numeric data to '
                            'plot'.format(numeric_data.__class__.__name__))

        self.data = numeric_data

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
            labels = [com.pprint_thing(key) for key in self.data.index]
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
                    name = com.pprint_thing(name)
                return name
            else:
                stringified = map(com.pprint_thing,
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
                self.data = self.data.reindex(index=index.order())
                x = self.data.index.to_timestamp()._mpl_repr()
            elif index.is_numeric():
                """
                Matplotlib supports numeric values or datetime objects as
                xaxis values. Taking LBYL approach here, by the time
                matplotlib raises exception when using non numeric/datetime
                values for xaxis, several actions are already taken by plt.
                """
                x = index._mpl_repr()
            elif is_datetype:
                self.data = self.data.sort_index()
                x = self.data.index._mpl_repr()
            else:
                self._need_to_set_index = True
                x = lrange(len(index))
        else:
            x = lrange(len(index))

        return x

    def _is_datetype(self):
        index = self.data.index
        return (isinstance(index, (PeriodIndex, DatetimeIndex)) or
                index.inferred_type in ('datetime', 'date', 'datetime64',
                                        'time'))

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
                name = ','.join([com.pprint_thing(x) for x in name])
            else:
                name = None
        else:
            name = self.data.index.name
            if name is not None:
                name = com.pprint_thing(name)

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

    def _get_colors(self):
        from pandas.core.frame import DataFrame
        if isinstance(self.data, DataFrame):
            num_colors = len(self.data.columns)
        else:
            num_colors = 1

        return _get_standard_colors(num_colors=num_colors,
                                    colormap=self.colormap,
                                    color=self.kwds.get('color'))

    def _maybe_add_color(self, colors, kwds, style, i):
        has_color = 'color' in kwds or self.colormap is not None
        if has_color and (style is None or re.match('[a-z]+', style) is None):
            kwds['color'] = colors[i % len(colors)]

    def _get_marked_label(self, label, col_num):
        if self.on_right(col_num):
            return label + ' (right)'
        else:
            return label


class KdePlot(MPLPlot):
    def __init__(self, data, bw_method=None, ind=None, **kwargs):
        MPLPlot.__init__(self, data, **kwargs)
        self.bw_method=bw_method
        self.ind=ind

    def _make_plot(self):
        from scipy.stats import gaussian_kde
        from scipy import __version__ as spv
        from distutils.version import LooseVersion
        plotf = self._get_plot_function()
        colors = self._get_colors()
        for i, (label, y) in enumerate(self._iter_data()):
            ax = self._get_ax(i)
            style = self._get_style(i, label)

            label = com.pprint_thing(label)

            if LooseVersion(spv) >= '0.11.0':
                gkde = gaussian_kde(y, bw_method=self.bw_method)
            else:
                gkde = gaussian_kde(y)
                if self.bw_method is not None:
                    msg = ('bw_method was added in Scipy 0.11.0.' +
                           ' Scipy version in use is %s.' % spv)
                    warnings.warn(msg)

            sample_range = max(y) - min(y)

            if self.ind is None:
                ind = np.linspace(min(y) - 0.5 * sample_range,
                                  max(y) + 0.5 * sample_range, 1000)
            else:
                ind = self.ind

            ax.set_ylabel("Density")

            y = gkde.evaluate(ind)
            kwds = self.kwds.copy()
            kwds['label'] = label
            self._maybe_add_color(colors, kwds, style, i)
            if style is None:
                args = (ax, ind, y)
            else:
                args = (ax, ind, y, style)

            plotf(*args, **kwds)
            ax.grid(self.grid)

    def _post_plot_logic(self):
        if self.legend:
            for ax in self.axes:
                ax.legend(loc='best')

class ScatterPlot(MPLPlot):
    def __init__(self, data, x, y, **kwargs):
        MPLPlot.__init__(self, data, **kwargs)
        self.kwds.setdefault('c', self.plt.rcParams['patch.facecolor'])
        if x is None or y is None:
            raise ValueError( 'scatter requires and x and y column')
        if com.is_integer(x) and not self.data.columns.holds_integer():
            x = self.data.columns[x]
        if com.is_integer(y) and not self.data.columns.holds_integer():
            y = self.data.columns[y]
        self.x = x
        self.y = y


    def _make_plot(self):
        x, y, data = self.x, self.y, self.data
        ax = self.axes[0]
        ax.scatter(data[x].values, data[y].values, **self.kwds)

    def _post_plot_logic(self):
        ax = self.axes[0]
        x, y = self.x, self.y
        ax.set_ylabel(com.pprint_thing(y))
        ax.set_xlabel(com.pprint_thing(x))


class LinePlot(MPLPlot):

    def __init__(self, data, **kwargs):
        self.mark_right = kwargs.pop('mark_right', True)
        MPLPlot.__init__(self, data, **kwargs)
        self.x_compat = plot_params['x_compat']
        if 'x_compat' in self.kwds:
            self.x_compat = bool(self.kwds.pop('x_compat'))

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
        return freq is not None and self._no_base(freq)

    def _no_base(self, freq):
        # hack this for 0.10.1, creating more technical debt...sigh
        from pandas.core.frame import DataFrame
        if (isinstance(self.data, (Series, DataFrame))
            and isinstance(self.data.index, DatetimeIndex)):
            import pandas.tseries.frequencies as freqmod
            base = freqmod.get_freq(freq)
            x = self.data.index
            if (base <= freqmod.FreqGroup.FR_DAY):
                return x[:1].is_normalized

            return Period(x[0], freq).to_timestamp(tz=x.tz) == x[0]
        return True

    def _use_dynamic_x(self):
        freq = self._index_freq()

        ax = self._get_ax(0)
        ax_freq = getattr(ax, 'freq', None)
        if freq is None:  # convert irregular if axes has freq info
            freq = ax_freq
        else:  # do not use tsplot if irregular was plotted first
            if (ax_freq is None) and (len(ax.get_lines()) > 0):
                return False

        return (freq is not None) and self._is_dynamic_freq(freq)

    def _make_plot(self):
        # this is slightly deceptive
        if not self.x_compat and self.use_index and self._use_dynamic_x():
            data = self._maybe_convert_index(self.data)
            self._make_ts_plot(data, **self.kwds)
        else:
            lines = []
            labels = []
            x = self._get_xticks(convert_period=True)

            plotf = self._get_plot_function()
            colors = self._get_colors()

            for i, (label, y) in enumerate(self._iter_data()):
                ax = self._get_ax(i)
                style = self._get_style(i, label)
                kwds = self.kwds.copy()
                self._maybe_add_color(colors, kwds, style, i)

                label = com.pprint_thing(label)  # .encode('utf-8')

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

                if self.mark_right:
                    labels.append(self._get_marked_label(label, i))
                else:
                    labels.append(label)

                ax.grid(self.grid)

                if self._is_datetype():
                    left, right = _get_xlim(lines)
                    ax.set_xlim(left, right)

            self._make_legend(lines, labels)

    def _make_ts_plot(self, data, **kwargs):
        from pandas.tseries.plotting import tsplot
        kwargs = kwargs.copy()
        colors = self._get_colors()

        plotf = self._get_plot_function()
        lines = []
        labels = []

        def _plot(data, col_num, ax, label, style, **kwds):
            newlines = tsplot(data, plotf, ax=ax, label=label,
                                style=style, **kwds)
            ax.grid(self.grid)
            lines.append(newlines[0])

            if self.mark_right:
                labels.append(self._get_marked_label(label, col_num))
            else:
                labels.append(label)

        if isinstance(data, Series):
            ax = self._get_ax(0)  # self.axes[0]
            style = self.style or ''
            label = com.pprint_thing(self.label)
            kwds = kwargs.copy()
            self._maybe_add_color(colors, kwds, style, 0)

            _plot(data, 0, ax, label, self.style, **kwds)
        else:
            for i, col in enumerate(data.columns):
                label = com.pprint_thing(col)
                ax = self._get_ax(i)
                style = self._get_style(i, col)
                kwds = kwargs.copy()

                self._maybe_add_color(colors, kwds, style, i)

                _plot(data[col], i, ax, label, style, **kwds)

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
            freq = get_base_alias(freq)
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

    _default_rot = {'bar': 90, 'barh': 0}

    def __init__(self, data, **kwargs):
        self.mark_right = kwargs.pop('mark_right', True)
        self.stacked = kwargs.pop('stacked', False)
        self.ax_pos = np.arange(len(data)) + 0.25
        if self.stacked:
            self.tickoffset = 0.25
        else:
            self.tickoffset = 0.375
        self.bar_width = 0.5
        self.log = kwargs.pop('log',False)
        MPLPlot.__init__(self, data, **kwargs)

    def _args_adjust(self):
        if self.rot is None:
            self.rot = self._default_rot[self.kind]

    @property
    def bar_f(self):
        if self.kind == 'bar':
            def f(ax, x, y, w, start=None, **kwds):
                return ax.bar(x, y, w, bottom=start,log=self.log, **kwds)
        elif self.kind == 'barh':
            def f(ax, x, y, w, start=None, log=self.log, **kwds):
                return ax.barh(x, y, w, left=start, **kwds)
        else:
            raise NotImplementedError

        return f

    def _make_plot(self):
        import matplotlib as mpl

        # mpl decided to make their version string unicode across all Python
        # versions for mpl >= 1.3 so we have to call str here for python 2
        mpl_le_1_2_1 = str(mpl.__version__) <= LooseVersion('1.2.1')

        colors = self._get_colors()
        ncolors = len(colors)
        rects = []
        labels = []

        bar_f = self.bar_f

        pos_prior = neg_prior = np.zeros(len(self.data))

        K = self.nseries

        for i, (label, y) in enumerate(self._iter_data()):
            ax = self._get_ax(i)
            label = com.pprint_thing(label)
            kwds = self.kwds.copy()
            kwds['color'] = colors[i % ncolors]

            start = 0
            if self.log:
                start = 1
                if any(y < 1):
                    # GH3254
                    start = 0 if mpl_le_1_2_1 else None

            if self.subplots:
                rect = bar_f(ax, self.ax_pos, y,  self.bar_width,
                             start=start, **kwds)
                ax.set_title(label)
            elif self.stacked:
                mask = y > 0
                start = np.where(mask, pos_prior, neg_prior)
                rect = bar_f(ax, self.ax_pos, y, self.bar_width, start=start,
                             label=label, **kwds)
                pos_prior = pos_prior + np.where(mask, y, 0)
                neg_prior = neg_prior + np.where(mask, 0, y)
            else:
                rect = bar_f(ax, self.ax_pos + i * 0.75 / K, y, 0.75 / K,
                             start=start, label=label, **kwds)
            rects.append(rect)
            if self.mark_right:
                labels.append(self._get_marked_label(label, i))
            else:
                labels.append(label)

        if self.legend and not self.subplots:
            patches = [r[0] for r in rects]
            self.axes[0].legend(patches, labels, loc='best',
                                title=self.legend_title)

    def _post_plot_logic(self):
        for ax in self.axes:
            if self.use_index:
                str_index = [com.pprint_thing(key) for key in self.data.index]
            else:
                str_index = [com.pprint_thing(key) for key in
                             range(self.data.shape[0])]
            name = self._get_index_name()
            if self.kind == 'bar':
                ax.set_xlim([self.ax_pos[0] - 0.25, self.ax_pos[-1] + 1])
                ax.set_xticks(self.ax_pos + self.tickoffset)
                ax.set_xticklabels(str_index, rotation=self.rot,
                                   fontsize=self.fontsize)
                if not self.log: # GH3254+
                    ax.axhline(0, color='k', linestyle='--')
                if name is not None:
                    ax.set_xlabel(name)
            else:
                # horizontal bars
                ax.set_ylim([self.ax_pos[0] - 0.25, self.ax_pos[-1] + 1])
                ax.set_yticks(self.ax_pos + self.tickoffset)
                ax.set_yticklabels(str_index, rotation=self.rot,
                                   fontsize=self.fontsize)
                ax.axvline(0, color='k', linestyle='--')
                if name is not None:
                    ax.set_ylabel(name)

        # if self.subplots and self.legend:
        #    self.axes[0].legend(loc='best')


class BoxPlot(MPLPlot):
    pass


class HistPlot(MPLPlot):
    pass


def plot_frame(frame=None, x=None, y=None, subplots=False, sharex=True,
               sharey=False, use_index=True, figsize=None, grid=None,
               legend=True, rot=None, ax=None, style=None, title=None,
               xlim=None, ylim=None, logx=False, logy=False, xticks=None,
               yticks=None, kind='line', sort_columns=False, fontsize=None,
               secondary_y=False, **kwds):

    """
    Make line, bar, or scatter plots of DataFrame series with the index on the x-axis
    using matplotlib / pylab.

    Parameters
    ----------
    frame : DataFrame
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
    grid : boolean, default None (matlab style default)
        Axis grid lines
    legend : boolean, default True
        Place legend on axis subplots

    ax : matplotlib axis object, default None
    style : list or dict
        matplotlib line style per column
    kind : {'line', 'bar', 'barh', 'kde', 'density', 'scatter'}
        bar : vertical bar plot
        barh : horizontal bar plot
        kde/density : Kernel Density Estimation plot
        scatter: scatter plot
    logx : boolean, default False
        For line plots, use log scaling on x axis
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
        If a list/tuple, which columns to plot on secondary y-axis
    mark_right: boolean, default True
        When using a secondary_y axis, should the legend label the axis of
        the various columns automatically
    colormap : str or matplotlib colormap object, default None
        Colormap to select colors from. If string, load colormap with that name
        from matplotlib.
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
    elif kind == 'scatter':
        klass = ScatterPlot
    else:
        raise ValueError('Invalid chart type given %s' % kind)

    if kind == 'scatter':
        plot_obj = klass(frame,  x=x, y=y, kind=kind, subplots=subplots,
                         rot=rot,legend=legend, ax=ax, style=style,
                         fontsize=fontsize, use_index=use_index, sharex=sharex,
                         sharey=sharey, xticks=xticks, yticks=yticks,
                         xlim=xlim, ylim=ylim, title=title, grid=grid,
                         figsize=figsize, logx=logx, logy=logy,
                         sort_columns=sort_columns, secondary_y=secondary_y,
                         **kwds)
    else:
        if x is not None:
            if com.is_integer(x) and not frame.columns.holds_integer():
                x = frame.columns[x]
            frame = frame.set_index(x)

        if y is not None:
            if com.is_integer(y) and not frame.columns.holds_integer():
                y = frame.columns[y]
            label = x if x is not None else frame.index.name
            label = kwds.pop('label', label)
            ser = frame[y]
            ser.index.name = label
            return plot_series(ser, label=label, kind=kind,
                               use_index=use_index,
                               rot=rot, xticks=xticks, yticks=yticks,
                               xlim=xlim, ylim=ylim, ax=ax, style=style,
                               grid=grid, logx=logx, logy=logy,
                               secondary_y=secondary_y, title=title,
                               figsize=figsize, fontsize=fontsize, **kwds)

        else:
            plot_obj = klass(frame, kind=kind, subplots=subplots, rot=rot,
                             legend=legend, ax=ax, style=style, fontsize=fontsize,
                             use_index=use_index, sharex=sharex, sharey=sharey,
                             xticks=xticks, yticks=yticks, xlim=xlim, ylim=ylim,
                             title=title, grid=grid, figsize=figsize, logx=logx,
                             logy=logy, sort_columns=sort_columns,
                             secondary_y=secondary_y, **kwds)

    plot_obj.generate()
    plot_obj.draw()
    if subplots:
        return plot_obj.axes
    else:
        return plot_obj.axes[0]


def plot_series(series, label=None, kind='line', use_index=True, rot=None,
                xticks=None, yticks=None, xlim=None, ylim=None,
                ax=None, style=None, grid=None, legend=False, logx=False,
                logy=False, secondary_y=False, **kwds):
    """
    Plot the input series with the index on the x-axis using matplotlib

    Parameters
    ----------
    label : label argument to provide to plot
    kind : {'line', 'bar', 'barh', 'kde', 'density'}
        bar : vertical bar plot
        barh : horizontal bar plot
        kde/density : Kernel Density Estimation plot
    use_index : boolean, default True
        Plot index as axis tick labels
    rot : int, default None
        Rotation for tick labels
    xticks : sequence
        Values to use for the xticks
    yticks : sequence
        Values to use for the yticks
    xlim : 2-tuple/list
    ylim : 2-tuple/list
    ax : matplotlib axis object
        If not passed, uses gca()
    style : string, default matplotlib default
        matplotlib line style to use
    grid : matplotlib grid
    legend: matplotlib legend
    logx : boolean, default False
        For line plots, use log scaling on x axis
    logy : boolean, default False
        For line plots, use log scaling on y axis
    secondary_y : boolean or sequence of ints, default False
        If True then y-axis will be on the right
    figsize : a tuple (width, height) in inches
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
    else:
        raise ValueError('Invalid chart type given %s' % kind)

    """
    If no axis is specified, we check whether there are existing figures.
    If so, we get the current axis and check whether yaxis ticks are on the
    right. Ticks for the plot of the series will be on the right unless
    there is at least one axis with ticks on the left.

    If we do not check for whether there are existing figures, _gca() will
    create a figure with the default figsize, causing the figsize= parameter to
    be ignored.
    """
    import matplotlib.pyplot as plt
    if ax is None and len(plt.get_fignums()) > 0:
        ax = _gca()
        if ax.get_yaxis().get_ticks_position().strip().lower() == 'right':
            fig = _gcf()
            axes = fig.get_axes()
            for i in reversed(range(len(axes))):
                ax = axes[i]
                ypos = ax.get_yaxis().get_ticks_position().strip().lower()
                if ypos == 'left':
                    break

    # is there harm in this?
    if label is None:
        label = series.name

    plot_obj = klass(series, kind=kind, rot=rot, logx=logx, logy=logy,
                     ax=ax, use_index=use_index, style=style,
                     xticks=xticks, yticks=yticks, xlim=xlim, ylim=ylim,
                     legend=legend, grid=grid, label=label,
                     secondary_y=secondary_y, **kwds)

    plot_obj.generate()
    plot_obj.draw()

    # plot_obj.ax is None if we created the first figure
    return plot_obj.axes[0]


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
    ax : Matplotlib axis object, optional
    fontsize : int or string
    rot : label rotation angle
    figsize : A tuple (width, height) in inches
    grid : Setting this to True will show the grid
    kwds : other plotting keyword arguments to be passed to matplotlib boxplot
           function

    Returns
    -------
    ax : matplotlib.axes.AxesSubplot
    """
    from pandas import Series, DataFrame
    if isinstance(data, Series):
        data = DataFrame({'x': data})
        column = 'x'


    def _get_colors():
        return _get_standard_colors(color=kwds.get('color'), num_colors=1)

    def maybe_color_bp(bp):
        if 'color' not in kwds :
            from matplotlib.artist import setp
            setp(bp['boxes'],color=colors[0],alpha=1)
            setp(bp['whiskers'],color=colors[0],alpha=1)
            setp(bp['medians'],color=colors[2],alpha=1)

    def plot_group(grouped, ax):
        keys, values = zip(*grouped)
        keys = [com.pprint_thing(x) for x in keys]
        values = [remove_na(v) for v in values]
        bp = ax.boxplot(values, **kwds)
        if kwds.get('vert', 1):
            ax.set_xticklabels(keys, rotation=rot, fontsize=fontsize)
        else:
            ax.set_yticklabels(keys, rotation=rot, fontsize=fontsize)
        maybe_color_bp(bp)

    colors = _get_colors()
    if column is None:
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
        keys = [com.pprint_thing(x) for x in cols]

        # Return boxplot dict in single plot case

        clean_values = [remove_na(x) for x in data[cols].values.T]

        bp = ax.boxplot(clean_values, **kwds)
        maybe_color_bp(bp)

        if kwds.get('vert', 1):
            ax.set_xticklabels(keys, rotation=rot, fontsize=fontsize)
        else:
            ax.set_yticklabels(keys, rotation=rot, fontsize=fontsize)
        ax.grid(grid)

        ret = bp

    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.9, wspace=0.2)
    return ret


def format_date_labels(ax, rot):
    # mini version of autofmt_xdate
    try:
        for label in ax.get_xticklabels():
            label.set_ha('right')
            label.set_rotation(rot)
        fig = ax.get_figure()
        fig.subplots_adjust(bottom=0.2)
    except Exception:  # pragma: no cover
        pass


def scatter_plot(data, x, y, by=None, ax=None, figsize=None, grid=False, **kwargs):
    """
    Make a scatter plot from two DataFrame columns

    Parameters
    ----------
    data : DataFrame
    x : Column name for the x-axis values
    y : Column name for the y-axis values
    ax : Matplotlib axis object
    figsize : A tuple (width, height) in inches
    grid : Setting this to True will show the grid
    kwargs : other plotting keyword arguments
        To be passed to scatter function

    Returns
    -------
    fig : matplotlib.Figure
    """
    import matplotlib.pyplot as plt

    # workaround because `c='b'` is hardcoded in matplotlibs scatter method
    kwargs.setdefault('c', plt.rcParams['patch.facecolor'])

    def plot_group(group, ax):
        xvals = group[x].values
        yvals = group[y].values
        ax.scatter(xvals, yvals, **kwargs)
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
        ax.set_ylabel(com.pprint_thing(y))
        ax.set_xlabel(com.pprint_thing(x))

        ax.grid(grid)

    return fig


def hist_frame(data, column=None, by=None, grid=True, xlabelsize=None,
               xrot=None, ylabelsize=None, yrot=None, ax=None, sharex=False,
               sharey=False, figsize=None, layout=None, **kwds):
    """
    Draw histogram of the DataFrame's series using matplotlib / pylab.

    Parameters
    ----------
    data : DataFrame
    column : string or sequence
        If passed, will be used to limit data to a subset of columns
    by : object, optional
        If passed, then used to form histograms for separate groups
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
    figsize : tuple
        The size of the figure to create in inches by default
    layout: (optional) a tuple (rows, columns) for the layout of the histograms
    kwds : other plotting keyword arguments
        To be passed to hist function
    """
    import matplotlib.pyplot as plt

    if column is not None:
        if not isinstance(column, (list, np.ndarray)):
            column = [column]
        data = data[column]

    if by is not None:
        axes = grouped_hist(data, by=by, ax=ax, grid=grid, figsize=figsize,
                            sharex=sharex, sharey=sharey, layout=layout,
                            **kwds)

        for ax in axes.ravel():
            if xlabelsize is not None:
                plt.setp(ax.get_xticklabels(), fontsize=xlabelsize)
            if xrot is not None:
                plt.setp(ax.get_xticklabels(), rotation=xrot)
            if ylabelsize is not None:
                plt.setp(ax.get_yticklabels(), fontsize=ylabelsize)
            if yrot is not None:
                plt.setp(ax.get_yticklabels(), rotation=yrot)

        return axes

    n = len(data.columns)

    if layout is not None:
        if not isinstance(layout, (tuple, list)) or len(layout) != 2:
            raise ValueError('Layout must be a tuple of (rows, columns)')

        rows, cols = layout
        if rows * cols < n:
            raise ValueError('Layout of %sx%s is incompatible with %s columns' % (rows, cols, n))
    else:
        rows, cols = 1, 1
        while rows * cols < n:
            if cols > rows:
                rows += 1
            else:
                cols += 1
    fig, axes = _subplots(nrows=rows, ncols=cols, ax=ax, squeeze=False,
                          sharex=sharex, sharey=sharey, figsize=figsize)

    for i, col in enumerate(com._try_sort(data.columns)):
        ax = axes[i / cols, i % cols]
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

    fig.subplots_adjust(wspace=0.3, hspace=0.3)

    return axes


def hist_series(self, by=None, ax=None, grid=True, xlabelsize=None,
                xrot=None, ylabelsize=None, yrot=None, figsize=None, **kwds):
    """
    Draw histogram of the input series using matplotlib

    Parameters
    ----------
    by : object, optional
        If passed, then used to form histograms for separate groups
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
    figsize : tuple, default None
        figure size in inches by default
    kwds : keywords
        To be passed to the actual plotting function

    Notes
    -----
    See matplotlib documentation online for more on this

    """
    import matplotlib.pyplot as plt

    if by is None:
        if kwds.get('layout', None) is not None:
            raise ValueError("The 'layout' keyword is not supported when "
                             "'by' is None")
        # hack until the plotting interface is a bit more unified
        fig = kwds.pop('figure', plt.gcf() if plt.get_fignums() else
                       plt.figure(figsize=figsize))
        if (figsize is not None and tuple(figsize) !=
            tuple(fig.get_size_inches())):
            fig.set_size_inches(*figsize, forward=True)
        if ax is None:
            ax = fig.gca()
        elif ax.get_figure() != fig:
            raise AssertionError('passed axis not bound to passed figure')
        values = self.dropna().values

        ax.hist(values, **kwds)
        ax.grid(grid)
        axes = np.array([ax])
    else:
        if 'figure' in kwds:
            raise ValueError("Cannot pass 'figure' when using the "
                             "'by' argument, since a new 'Figure' instance "
                             "will be created")
        axes = grouped_hist(self, by=by, ax=ax, grid=grid, figsize=figsize,
                            **kwds)

    for ax in axes.ravel():
        if xlabelsize is not None:
            plt.setp(ax.get_xticklabels(), fontsize=xlabelsize)
        if xrot is not None:
            plt.setp(ax.get_xticklabels(), rotation=xrot)
        if ylabelsize is not None:
            plt.setp(ax.get_yticklabels(), fontsize=ylabelsize)
        if yrot is not None:
            plt.setp(ax.get_yticklabels(), rotation=yrot)

    if axes.ndim == 1 and len(axes) == 1:
        return axes[0]
    return axes


def boxplot_frame_groupby(grouped, subplots=True, column=None, fontsize=None,
                          rot=0, grid=True, figsize=None, **kwds):
    """
    Make box plots from DataFrameGroupBy data.

    Parameters
    ----------
    grouped : Grouped DataFrame
    subplots :
        * ``False`` - no subplots will be used
        * ``True`` - create a subplot for each group
    column : column name or list of names, or vector
        Can be any valid input to groupby
    fontsize : int or string
    rot : label rotation angle
    grid : Setting this to True will show the grid
    figsize : A tuple (width, height) in inches
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
            ax.set_title(com.pprint_thing(key))
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
                  rot=0, ax=None, **kwargs):
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

    if nrows * ncols < ngroups:
        raise ValueError("Number of plots in 'layout' must greater than or "
                         "equal to the number " "of groups in 'by'")

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
        plotf(group, ax, **kwargs)
        ax.set_title(com.pprint_thing(key))

    return fig, axes


def _grouped_plot_by_column(plotf, data, columns=None, by=None,
                            numeric_only=True, grid=False,
                            figsize=None, ax=None, **kwargs):
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
        plotf(gp_col, ax, **kwargs)
        ax.set_title(col)
        ax.set_xlabel(com.pprint_thing(by))
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

    ax : Matplotlib axis object, optional

    secondary_y : boolean or sequence of ints, default False
        If True then y-axis will be on the right

    data : DataFrame, optional
        If secondary_y is a sequence, data is used to select columns.

    fig_kw : Other keyword arguments to be passed to the figure() call.
        Note that all keywords not recognized above will be
        automatically included here.


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
    nplots = nrows * ncols
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
        ax0._get_lines.color_cycle = orig_ax._get_lines.color_cycle

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
        ax = fig.add_subplot(nrows, ncols, i + 1, **subplot_kw)
        if on_right(i):
            orig_ax = ax
            ax = ax.twinx()
            ax._get_lines.color_cycle = orig_ax._get_lines.color_cycle

            orig_ax.get_yaxis().set_visible(False)
        axarr[i] = ax

    if nplots > 1:
        if sharex and nrows > 1:
            for i, ax in enumerate(axarr):
                if np.ceil(float(i + 1) / ncols) < nrows:  # only last row
                    [label.set_visible(
                        False) for label in ax.get_xticklabels()]
        if sharey and ncols > 1:
            for i, ax in enumerate(axarr):
                if (i % ncols) != 0:  # only first column
                    [label.set_visible(
                        False) for label in ax.get_yticklabels()]

    if squeeze:
        # Reshape the array to have the final desired dimension (nrow,ncol),
        # though discarding unneeded dimensions that equal 1.  If we only have
        # one subplot, just return it instead of a 1-element array.
        if nplots == 1:
            axes = axarr[0]
        else:
            axes = axarr.reshape(nrows, ncols).squeeze()
    else:
        # returned axis array will be always 2-d, even if nrows=ncols=1
        axes = axarr.reshape(nrows, ncols)

    return fig, axes


def _get_xlim(lines):
    left, right = np.inf, -np.inf
    for l in lines:
        x = l.get_xdata()
        left = min(_maybe_convert_date(x[0]), left)
        right = max(_maybe_convert_date(x[-1]), right)
    return left, right


def _maybe_convert_date(x):
    if not com.is_integer(x):
        conv_func = conv._dt_to_float_ordinal
        if isinstance(x, datetime.time):
            conv_func = conv._to_ordinalf
        x = conv_func(x)
    return x

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
