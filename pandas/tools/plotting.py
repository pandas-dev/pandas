# being a bit too dynamic
# pylint: disable=E1101
import datetime
import warnings
import re
from collections import namedtuple
from contextlib import contextmanager
from distutils.version import LooseVersion

import numpy as np

from pandas.util.decorators import cache_readonly, deprecate_kwarg
import pandas.core.common as com
from pandas.core.generic import _shared_docs, _shared_doc_kwargs
from pandas.core.index import MultiIndex
from pandas.core.series import Series, remove_na
from pandas.tseries.index import DatetimeIndex
from pandas.tseries.period import PeriodIndex, Period
from pandas.tseries.frequencies import get_period_alias, get_base_alias
from pandas.tseries.offsets import DateOffset
from pandas.compat import range, lrange, lmap, map, zip, string_types
import pandas.compat as compat
from pandas.util.decorators import Appender

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

def radviz(frame, class_column, ax=None, color=None, colormap=None, **kwds):
    """RadViz - a multivariate data visualization algorithm

    Parameters:
    -----------
    frame: DataFrame
    class_column: str
        Column name containing class names
    ax: Matplotlib axis object, optional
    color: list or tuple, optional
        Colors to use for the different classes
    colormap : str or matplotlib colormap object, default None
        Colormap to select colors from. If string, load colormap with that name
        from matplotlib.
    kwds: keywords
        Options to pass to matplotlib scatter plotting method

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

    n = len(frame)
    classes = frame[class_column].drop_duplicates()
    class_col = frame[class_column]
    df = frame.drop(class_column, axis=1).apply(normalize)

    if ax is None:
        ax = plt.gca(xlim=[-1, 1], ylim=[-1, 1])

    to_plot = {}
    colors = _get_standard_colors(num_colors=len(classes), colormap=colormap,
                                  color_type='random', color=color)

    for kls in classes:
        to_plot[kls] = [[], []]

    n = len(frame.columns) - 1
    s = np.array([(np.cos(t), np.sin(t))
                  for t in [2.0 * np.pi * (i / float(n))
                            for i in range(n)]])

    for i in range(n):
        row = df.iloc[i].values
        row_ = np.repeat(np.expand_dims(row, axis=1), 2, axis=1)
        y = (s * row_).sum(axis=0) / row.sum()
        kls = class_col.iat[i]
        to_plot[kls][0].append(y[0])
        to_plot[kls][1].append(y[1])

    for i, kls in enumerate(classes):
        ax.scatter(to_plot[kls][0], to_plot[kls][1], color=colors[i],
                   label=com.pprint_thing(kls), **kwds)
    ax.legend()

    ax.add_patch(patches.Circle((0.0, 0.0), radius=1.0, facecolor='none'))

    for xy, name in zip(s, df.columns):

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

@deprecate_kwarg(old_arg_name='data', new_arg_name='frame')
def andrews_curves(frame, class_column, ax=None, samples=200, color=None,
                   colormap=None, **kwds):
    """
    Parameters:
    -----------
    frame : DataFrame
        Data to be plotted, preferably normalized to (0.0, 1.0)
    class_column : Name of the column containing class names
    ax : matplotlib axes object, default None
    samples : Number of points to plot in each curve
    color: list or tuple, optional
        Colors to use for the different classes
    colormap : str or matplotlib colormap object, default None
        Colormap to select colors from. If string, load colormap with that name
        from matplotlib.
    kwds: keywords
        Options to pass to matplotlib plotting method

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

    n = len(frame)
    class_col = frame[class_column]
    classes = frame[class_column].drop_duplicates()
    df = frame.drop(class_column, axis=1)
    x = [-pi + 2.0 * pi * (t / float(samples)) for t in range(samples)]
    used_legends = set([])

    color_values = _get_standard_colors(num_colors=len(classes),
                                        colormap=colormap, color_type='random',
                                        color=color)
    colors = dict(zip(classes, color_values))
    if ax is None:
        ax = plt.gca(xlim=(-pi, pi))
    for i in range(n):
        row = df.iloc[i].values
        f = function(row)
        y = [f(t) for t in x]
        kls = class_col.iat[i]
        label = com.pprint_thing(kls)
        if label not in used_legends:
            used_legends.add(label)
            ax.plot(x, y, color=colors[kls], label=label, **kwds)
        else:
            ax.plot(x, y, color=colors[kls], **kwds)

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

@deprecate_kwarg(old_arg_name='colors', new_arg_name='color')
@deprecate_kwarg(old_arg_name='data', new_arg_name='frame')
def parallel_coordinates(frame, class_column, cols=None, ax=None, color=None,
                         use_columns=False, xticks=None, colormap=None,
                         **kwds):
    """Parallel coordinates plotting.

    Parameters
    ----------
    frame: DataFrame
    class_column: str
        Column name containing class names
    cols: list, optional
        A list of column names to use
    ax: matplotlib.axis, optional
        matplotlib axis object
    color: list or tuple, optional
        Colors to use for the different classes
    use_columns: bool, optional
        If true, columns will be used as xticks
    xticks: list or tuple, optional
        A list of values to use for xticks
    colormap: str or matplotlib colormap, default None
        Colormap to use for line colors.
    kwds: keywords
        Options to pass to matplotlib plotting method

    Returns
    -------
    ax: matplotlib axis object

    Examples
    --------
    >>> from pandas import read_csv
    >>> from pandas.tools.plotting import parallel_coordinates
    >>> from matplotlib import pyplot as plt
    >>> df = read_csv('https://raw.github.com/pydata/pandas/master/pandas/tests/data/iris.csv')
    >>> parallel_coordinates(df, 'Name', color=('#556270', '#4ECDC4', '#C7F464'))
    >>> plt.show()
    """
    import matplotlib.pyplot as plt

    n = len(frame)
    classes = frame[class_column].drop_duplicates()
    class_col = frame[class_column]

    if cols is None:
        df = frame.drop(class_column, axis=1)
    else:
        df = frame[cols]

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
                                        color=color)

    colors = dict(zip(classes, color_values))

    for i in range(n):
        y = df.iloc[i].values
        kls = class_col.iat[i]
        label = com.pprint_thing(kls)
        if label not in used_legends:
            used_legends.add(label)
            ax.plot(x, y, color=colors[kls], label=label, **kwds)
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


class MPLPlot(object):
    """
    Base class for assembling a pandas plot using matplotlib

    Parameters
    ----------
    data :

    """
    _default_rot = 0

    _pop_attributes = ['label', 'style', 'logy', 'logx', 'loglog',
                       'mark_right']
    _attr_defaults = {'logy': False, 'logx': False, 'loglog': False,
                      'mark_right': True}

    def __init__(self, data, kind=None, by=None, subplots=False, sharex=True,
                 sharey=False, use_index=True,
                 figsize=None, grid=None, legend=True, rot=None,
                 ax=None, fig=None, title=None, xlim=None, ylim=None,
                 xticks=None, yticks=None,
                 sort_columns=False, fontsize=None,
                 secondary_y=False, colormap=None,
                 table=False, **kwds):

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
        self.legend_handles = []
        self.legend_labels = []

        for attr in self._pop_attributes:
            value = kwds.pop(attr, self._attr_defaults.get(attr, None))
            setattr(self, attr, value)

        self.ax = ax
        self.fig = fig
        self.axes = None

        # parse errorbar input if given
        xerr = kwds.pop('xerr', None)
        yerr = kwds.pop('yerr', None)
        self.errors = {}
        for kw, err in zip(['xerr', 'yerr'], [xerr, yerr]):
            self.errors[kw] = self._parse_errorbars(kw, err)

        if not isinstance(secondary_y, (bool, tuple, list, np.ndarray)):
            secondary_y = [secondary_y]
        self.secondary_y = secondary_y

        # ugly TypeError if user passes matplotlib's `cmap` name.
        # Probably better to accept either.
        if 'cmap' in kwds and colormap:
            raise TypeError("Only specify one of `cmap` and `colormap`.")
        elif 'cmap' in kwds:
            self.colormap = kwds.pop('cmap')
        else:
            self.colormap = colormap

        self.table = table

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

    def _iter_data(self, data=None, keep_index=False):
        if data is None:
            data = self.data

        from pandas.core.frame import DataFrame
        if isinstance(data, (Series, np.ndarray)):
            if keep_index is True:
                yield self.label, data
            else:
                yield self.label, np.asarray(data)
        elif isinstance(data, DataFrame):
            if self.sort_columns:
                columns = com._try_sort(data.columns)
            else:
                columns = data.columns

            for col in columns:
                # # is this right?
                # empty = df[col].count() == 0
                # values = df[col].values if not empty else np.zeros(len(df))

                if keep_index is True:
                    yield col, data[col]
                else:
                    yield col, data[col].values

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
        self._add_table()
        self._make_legend()
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
            fig, axes = _subplots(nrows=nrows, ncols=ncols,
                                  sharex=self.sharex, sharey=self.sharey,
                                  figsize=self.figsize, ax=self.ax,
                                  secondary_y=self.secondary_y,
                                  data=self.data)
            if not com.is_list_like(axes):
                axes = np.array([axes])
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

        if self.logx or self.loglog:
            [a.set_xscale('log') for a in axes]
        if self.logy or self.loglog:
            [a.set_yscale('log') for a in axes]

        self.fig = fig
        self.axes = axes

    def _get_layout(self):
        from pandas.core.frame import DataFrame
        if isinstance(self.data, DataFrame):
            return (len(self.data.columns), 1)
        else:
            return (1, 1)

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

    def _add_table(self):
        if self.table is False:
            return
        elif self.table is True:
            from pandas.core.frame import DataFrame
            if isinstance(self.data, Series):
                data = DataFrame(self.data, columns=[self.data.name])
            elif isinstance(self.data, DataFrame):
                data = self.data
            data = data.transpose()
        else:
            data = self.table
        ax = self._get_ax(0)
        table(ax, data)

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

    def _add_legend_handle(self, handle, label, index=None):
        if not label is None:
            if self.mark_right and index is not None:
                if self.on_right(index):
                    label = label + ' (right)'
            self.legend_handles.append(handle)
            self.legend_labels.append(label)

    def _make_legend(self):
        ax, leg = self._get_ax_legend(self.axes[0])

        handles = []
        labels = []
        title = ''

        if not self.subplots:
            if not leg is None:
                title = leg.get_title().get_text()
                handles = leg.legendHandles
                labels = [x.get_text() for x in leg.get_texts()]

            if self.legend:
                if self.legend == 'reverse':
                    self.legend_handles = reversed(self.legend_handles)
                    self.legend_labels = reversed(self.legend_labels)

                handles += self.legend_handles
                labels += self.legend_labels
                if not self.legend_title is None:
                    title = self.legend_title

            if len(handles) > 0:
                ax.legend(handles, labels, loc='best', title=title)

        elif self.subplots and self.legend:
            for ax in self.axes:
                ax.legend(loc='best')


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
        '''
        Returns the matplotlib plotting function (plot or errorbar) based on
        the presence of errorbar keywords.
        '''

        if all(e is None for e in self.errors.values()):
            plotf = self.plt.Axes.plot
        else:
            plotf = self.plt.Axes.errorbar

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

    def _get_colors(self, num_colors=None, color_kwds='color'):
        from pandas.core.frame import DataFrame
        if num_colors is None:
            if isinstance(self.data, DataFrame):
                num_colors = len(self.data.columns)
            else:
                num_colors = 1

        return _get_standard_colors(num_colors=num_colors,
                                    colormap=self.colormap,
                                    color=self.kwds.get(color_kwds))

    def _maybe_add_color(self, colors, kwds, style, i):
        has_color = 'color' in kwds or self.colormap is not None
        if has_color and (style is None or re.match('[a-z]+', style) is None):
            kwds['color'] = colors[i % len(colors)]

    def _parse_errorbars(self, label, err):
        '''
        Look for error keyword arguments and return the actual errorbar data
        or return the error DataFrame/dict

        Error bars can be specified in several ways:
            Series: the user provides a pandas.Series object of the same
                    length as the data
            ndarray: provides a np.ndarray of the same length as the data
            DataFrame/dict: error values are paired with keys matching the
                    key in the plotted DataFrame
            str: the name of the column within the plotted DataFrame
        '''

        if err is None:
            return None

        from pandas import DataFrame, Series

        def match_labels(data, e):
            e = e.reindex_axis(data.index)
            return e

        # key-matched DataFrame
        if isinstance(err, DataFrame):

            err = match_labels(self.data, err)
        # key-matched dict
        elif isinstance(err, dict):
            pass

        # Series of error values
        elif isinstance(err, Series):
            # broadcast error series across data
            err = match_labels(self.data, err)
            err = np.atleast_2d(err)
            err = np.tile(err, (self.nseries, 1))

        # errors are a column in the dataframe
        elif isinstance(err, string_types):
            evalues = self.data[err].values
            self.data = self.data[self.data.columns.drop(err)]
            err = np.atleast_2d(evalues)
            err = np.tile(err, (self.nseries, 1))

        elif com.is_list_like(err):
            if com.is_iterator(err):
                err = np.atleast_2d(list(err))
            else:
                # raw error values
                err = np.atleast_2d(err)

            err_shape = err.shape

            # asymmetrical error bars
            if err.ndim == 3:
                if (err_shape[0] != self.nseries) or \
                    (err_shape[1] != 2) or \
                    (err_shape[2] != len(self.data)):
                    msg = "Asymmetrical error bars should be provided " + \
                    "with the shape (%u, 2, %u)" % \
                        (self.nseries, len(self.data))
                    raise ValueError(msg)

            # broadcast errors to each data series
            if len(err) == 1:
                err = np.tile(err, (self.nseries, 1))

        elif com.is_number(err):
            err = np.tile([err], (self.nseries, len(self.data)))

        else:
            msg = "No valid %s detected" % label
            raise ValueError(msg)

        return err

    def _get_errorbars(self, label=None, index=None, xerr=True, yerr=True):
        from pandas import DataFrame
        errors = {}

        for kw, flag in zip(['xerr', 'yerr'], [xerr, yerr]):
            if flag:
                err = self.errors[kw]
                # user provided label-matched dataframe of errors
                if isinstance(err, (DataFrame, dict)):
                    if label is not None and label in err.keys():
                        err = err[label]
                    else:
                        err = None
                elif index is not None and err is not None:
                    err = err[index]

                if err is not None:
                    errors[kw] = err
        return errors


class KdePlot(MPLPlot):
    def __init__(self, data, bw_method=None, ind=None, **kwargs):
        MPLPlot.__init__(self, data, **kwargs)
        self.bw_method=bw_method
        self.ind=ind

    def _make_plot(self):
        from scipy.stats import gaussian_kde
        from scipy import __version__ as spv
        from distutils.version import LooseVersion
        plotf = self.plt.Axes.plot
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

            newlines = plotf(*args, **kwds)
            self._add_legend_handle(newlines[0], label)


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

    def _get_layout(self):
        return (1, 1)

    def _make_plot(self):
        x, y, data = self.x, self.y, self.data
        ax = self.axes[0]

        if self.legend and hasattr(self, 'label'):
            label = self.label
        else:
            label = None
        scatter = ax.scatter(data[x].values, data[y].values, label=label,
                             **self.kwds)

        self._add_legend_handle(scatter, label)

        errors_x = self._get_errorbars(label=x, index=0, yerr=False)
        errors_y = self._get_errorbars(label=y, index=1, xerr=False)
        if len(errors_x) > 0 or len(errors_y) > 0:
            err_kwds = dict(errors_x, **errors_y)
            if 'color' in self.kwds:
                err_kwds['color'] = self.kwds['color']
            ax.errorbar(data[x].values, data[y].values, linestyle='none', **err_kwds)

    def _post_plot_logic(self):
        ax = self.axes[0]
        x, y = self.x, self.y
        ax.set_ylabel(com.pprint_thing(y))
        ax.set_xlabel(com.pprint_thing(x))


class HexBinPlot(MPLPlot):
    def __init__(self, data, x, y, C=None, **kwargs):
        MPLPlot.__init__(self, data, **kwargs)

        if x is None or y is None:
            raise ValueError('hexbin requires and x and y column')
        if com.is_integer(x) and not self.data.columns.holds_integer():
            x = self.data.columns[x]
        if com.is_integer(y) and not self.data.columns.holds_integer():
            y = self.data.columns[y]

        if com.is_integer(C) and not self.data.columns.holds_integer():
            C = self.data.columns[C]

        self.x = x
        self.y = y
        self.C = C

    def _get_layout(self):
        return (1, 1)

    def _make_plot(self):
        import matplotlib.pyplot as plt

        x, y, data, C = self.x, self.y, self.data, self.C
        ax = self.axes[0]
        # pandas uses colormap, matplotlib uses cmap.
        cmap = self.colormap or 'BuGn'
        cmap = plt.cm.get_cmap(cmap)
        cb = self.kwds.pop('colorbar', True)

        if C is None:
            c_values = None
        else:
            c_values = data[C].values

        ax.hexbin(data[x].values, data[y].values, C=c_values, cmap=cmap,
                  **self.kwds)
        if cb:
            img = ax.collections[0]
            self.fig.colorbar(img, ax=ax)

    def _post_plot_logic(self):
        ax = self.axes[0]
        x, y = self.x, self.y
        ax.set_ylabel(com.pprint_thing(y))
        ax.set_xlabel(com.pprint_thing(x))


class LinePlot(MPLPlot):

    def __init__(self, data, **kwargs):
        self.stacked = kwargs.pop('stacked', False)
        if self.stacked:
            data = data.fillna(value=0)

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

    def _is_ts_plot(self):
        # this is slightly deceptive
        return not self.x_compat and self.use_index and self._use_dynamic_x()

    def _make_plot(self):
        self._pos_prior = np.zeros(len(self.data))
        self._neg_prior = np.zeros(len(self.data))

        if self._is_ts_plot():
            data = self._maybe_convert_index(self.data)
            self._make_ts_plot(data)
        else:
            from pandas.core.frame import DataFrame
            lines = []
            x = self._get_xticks(convert_period=True)

            plotf = self._get_plot_function()
            colors = self._get_colors()

            for i, (label, y) in enumerate(self._iter_data()):
                ax = self._get_ax(i)
                style = self._get_style(i, label)
                kwds = self.kwds.copy()
                self._maybe_add_color(colors, kwds, style, i)

                errors = self._get_errorbars(label=label, index=i)
                kwds = dict(kwds, **errors)

                label = com.pprint_thing(label)  # .encode('utf-8')
                kwds['label'] = label

                y_values = self._get_stacked_values(y, label)

                if not self.stacked:
                    mask = com.isnull(y_values)
                    if mask.any():
                        y_values = np.ma.array(y_values)
                        y_values = np.ma.masked_where(mask, y_values)

                # prevent style kwarg from going to errorbar, where it is unsupported
                if style is not None and plotf.__name__ != 'errorbar':
                    args = (ax, x, y_values, style)
                else:
                    args = (ax, x, y_values)

                newlines = plotf(*args, **kwds)
                self._add_legend_handle(newlines[0], label, index=i)

                lines.append(newlines[0])

                if self.stacked and not self.subplots:
                    if (y >= 0).all():
                        self._pos_prior += y
                    elif (y <= 0).all():
                        self._neg_prior += y

            if self._is_datetype():
                left, right = _get_xlim(lines)
                ax.set_xlim(left, right)

    def _get_stacked_values(self, y, label):
        if self.stacked:
            if (y >= 0).all():
                return self._pos_prior + y
            elif (y <= 0).all():
                return self._neg_prior + y
            else:
                raise ValueError('When stacked is True, each column must be either all positive or negative.'
                                 '{0} contains both positive and negative values'.format(label))
        else:
            return y

    def _get_ts_plot_function(self):
        from pandas.tseries.plotting import tsplot
        plotf = self._get_plot_function()

        def _plot(data, ax, label, style, **kwds):
            # errorbar function does not support style argument
            if plotf.__name__ == 'errorbar':
                lines = tsplot(data, plotf, ax=ax, label=label,
                               **kwds)
                return lines
            else:
                lines = tsplot(data, plotf, ax=ax, label=label,
                               style=style, **kwds)
                return lines
        return _plot

    def _make_ts_plot(self, data, **kwargs):
        colors = self._get_colors()
        plotf = self._get_ts_plot_function()

        it = self._iter_data(data=data, keep_index=True)
        for i, (label, y) in enumerate(it):
            ax = self._get_ax(i)
            style = self._get_style(i, label)
            kwds = self.kwds.copy()

            self._maybe_add_color(colors, kwds, style, i)

            errors = self._get_errorbars(label=label, index=i, xerr=False)
            kwds = dict(kwds, **errors)

            label = com.pprint_thing(label)

            y_values = self._get_stacked_values(y, label)

            newlines = plotf(y_values, ax, label, style, **kwds)
            self._add_legend_handle(newlines[0], label, index=i)

            if self.stacked and not self.subplots:
                if (y >= 0).all():
                    self._pos_prior += y
                elif (y <= 0).all():
                    self._neg_prior += y

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


class AreaPlot(LinePlot):

    def __init__(self, data, **kwargs):
        kwargs.setdefault('stacked', True)
        data = data.fillna(value=0)
        LinePlot.__init__(self, data, **kwargs)

        if not self.stacked:
            # use smaller alpha to distinguish overlap
            self.kwds.setdefault('alpha', 0.5)

    def _get_plot_function(self):
        if self.logy or self.loglog:
            raise ValueError("Log-y scales are not supported in area plot")
        else:
            f = LinePlot._get_plot_function(self)

            def plotf(*args, **kwds):
                lines = f(*args, **kwds)

                # insert fill_between starting point
                y = args[2]
                if (y >= 0).all():
                    start = self._pos_prior
                elif (y <= 0).all():
                    start = self._neg_prior
                else:
                    start = np.zeros(len(y))

                # get x data from the line
                # to retrieve x coodinates of tsplot
                xdata = lines[0].get_data()[0]
                # remove style
                args = (args[0], xdata, start, y)

                if not 'color' in kwds:
                    kwds['color'] = lines[0].get_color()

                self.plt.Axes.fill_between(*args, **kwds)
                return lines

        return plotf

    def _add_legend_handle(self, handle, label, index=None):
        from matplotlib.patches import Rectangle
        # Because fill_between isn't supported in legend,
        # specifically add Rectangle handle here
        alpha = self.kwds.get('alpha', 0.5)
        handle = Rectangle((0, 0), 1, 1, fc=handle.get_color(), alpha=alpha)
        LinePlot._add_legend_handle(self, handle, label, index=index)

    def _post_plot_logic(self):
        LinePlot._post_plot_logic(self)

        if self._is_ts_plot():
            pass
        else:
            if self.xlim is None:
                for ax in self.axes:
                    ax.set_xlim(0, len(self.data)-1)

        if self.ylim is None:
            if (self.data >= 0).all().all():
                for ax in self.axes:
                    ax.set_ylim(0, None)
            elif (self.data <= 0).all().all():
                for ax in self.axes:
                    ax.set_ylim(None, 0)


class BarPlot(MPLPlot):

    _default_rot = {'bar': 90, 'barh': 0}

    def __init__(self, data, **kwargs):
        self.stacked = kwargs.pop('stacked', False)

        self.bar_width = kwargs.pop('width', 0.5)
        pos = kwargs.pop('position', 0.5)

        kwargs['align'] = kwargs.pop('align', 'center')
        self.tick_pos = np.arange(len(data))

        self.bottom = kwargs.pop('bottom', None)
        self.left = kwargs.pop('left', None)

        self.log = kwargs.pop('log',False)
        MPLPlot.__init__(self, data, **kwargs)

        if self.stacked or self.subplots:
            self.tickoffset = self.bar_width * pos
        elif kwargs['align'] == 'edge':
            K = self.nseries
            w = self.bar_width / K
            self.tickoffset = self.bar_width * (pos - 0.5) + w * 1.5
        else:
            K = self.nseries
            w = self.bar_width / K
            self.tickoffset = self.bar_width * pos + w
        self.ax_pos = self.tick_pos - self.tickoffset

    def _args_adjust(self):
        if self.rot is None:
            self.rot = self._default_rot[self.kind]

        if com.is_list_like(self.bottom):
            self.bottom = np.array(self.bottom)
        if com.is_list_like(self.left):
            self.left = np.array(self.left)

    def _get_plot_function(self):
        if self.kind == 'bar':
            def f(ax, x, y, w, start=None, **kwds):
                if self.bottom is not None:
                    start = start + self.bottom
                return ax.bar(x, y, w, bottom=start,log=self.log, **kwds)
        elif self.kind == 'barh':
            def f(ax, x, y, w, start=None, log=self.log, **kwds):
                if self.left is not None:
                    start = start + self.left
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

        bar_f = self._get_plot_function()
        pos_prior = neg_prior = np.zeros(len(self.data))
        K = self.nseries

        for i, (label, y) in enumerate(self._iter_data()):
            ax = self._get_ax(i)
            kwds = self.kwds.copy()
            kwds['color'] = colors[i % ncolors]

            errors = self._get_errorbars(label=label, index=i)
            kwds = dict(kwds, **errors)

            label = com.pprint_thing(label)

            if (('yerr' in kwds) or ('xerr' in kwds)) \
                and (kwds.get('ecolor') is None):
                kwds['ecolor'] = mpl.rcParams['xtick.color']

            start = 0
            if self.log:
                start = 1
                if any(y < 1):
                    # GH3254
                    start = 0 if mpl_le_1_2_1 else None

            if self.subplots:
                w = self.bar_width / 2
                rect = bar_f(ax, self.ax_pos + w, y, self.bar_width,
                             start=start, label=label, **kwds)
                ax.set_title(label)
            elif self.stacked:
                mask = y > 0
                start = np.where(mask, pos_prior, neg_prior)
                w = self.bar_width / 2
                rect = bar_f(ax, self.ax_pos + w, y, self.bar_width,
                             start=start, label=label, **kwds)
                pos_prior = pos_prior + np.where(mask, y, 0)
                neg_prior = neg_prior + np.where(mask, 0, y)
            else:
                w = self.bar_width / K
                rect = bar_f(ax, self.ax_pos + (i + 1.5) * w, y, w,
                             start=start, label=label, **kwds)

            self._add_legend_handle(rect, label, index=i)

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
                ax.set_xticks(self.tick_pos)
                ax.set_xticklabels(str_index, rotation=self.rot,
                                   fontsize=self.fontsize)
                if not self.log: # GH3254+
                    ax.axhline(0, color='k', linestyle='--')
                if name is not None:
                    ax.set_xlabel(name)
            elif self.kind == 'barh':
                # horizontal bars
                ax.set_ylim([self.ax_pos[0] - 0.25, self.ax_pos[-1] + 1])
                ax.set_yticks(self.tick_pos)
                ax.set_yticklabels(str_index, rotation=self.rot,
                                   fontsize=self.fontsize)
                ax.axvline(0, color='k', linestyle='--')
                if name is not None:
                    ax.set_ylabel(name)
            else:
                raise NotImplementedError(self.kind)

        # if self.subplots and self.legend:
        #    self.axes[0].legend(loc='best')


class PiePlot(MPLPlot):

    def __init__(self, data, kind=None, **kwargs):
        data = data.fillna(value=0)
        if (data < 0).any().any():
            raise ValueError("{0} doesn't allow negative values".format(kind))
        MPLPlot.__init__(self, data, kind=kind, **kwargs)

    def _args_adjust(self):
        self.grid = False
        self.logy = False
        self.logx = False
        self.loglog = False

    def _get_layout(self):
        from pandas import DataFrame
        if isinstance(self.data, DataFrame):
            return (1, len(self.data.columns))
        else:
            return (1, 1)

    def _validate_color_args(self):
        pass

    def _make_plot(self):
        self.kwds.setdefault('colors', self._get_colors(num_colors=len(self.data),
                                                        color_kwds='colors'))

        for i, (label, y) in enumerate(self._iter_data()):
            ax = self._get_ax(i)
            if label is not None:
                label = com.pprint_thing(label)
                ax.set_ylabel(label)

            kwds = self.kwds.copy()

            idx = [com.pprint_thing(v) for v in self.data.index]
            labels = kwds.pop('labels', idx)
            # labels is used for each wedge's labels
            results = ax.pie(y, labels=labels, **kwds)

            if kwds.get('autopct', None) is not None:
                patches, texts, autotexts = results
            else:
                patches, texts = results
                autotexts = []

            if self.fontsize is not None:
                for t in texts + autotexts:
                    t.set_fontsize(self.fontsize)

            # leglabels is used for legend labels
            leglabels = labels if labels is not None else idx
            for p, l in zip(patches, leglabels):
                self._add_legend_handle(p, l)


class BoxPlot(MPLPlot):
    pass


class HistPlot(MPLPlot):
    pass

# kinds supported by both dataframe and series
_common_kinds = ['line', 'bar', 'barh', 'kde', 'density', 'area']
# kinds supported by dataframe
_dataframe_kinds = ['scatter', 'hexbin']
# kinds supported only by series or dataframe single column
_series_kinds = ['pie']
_all_kinds = _common_kinds + _dataframe_kinds + _series_kinds

_plot_klass = {'line': LinePlot, 'bar': BarPlot, 'barh': BarPlot,
               'kde': KdePlot,
               'scatter': ScatterPlot, 'hexbin': HexBinPlot,
               'area': AreaPlot, 'pie': PiePlot}


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
    yerr : DataFrame (with matching labels), Series, list-type (tuple, list,
        ndarray), or str of column name containing y error values
    xerr : similar functionality as yerr, but for x error values
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
    legend : False/True/'reverse'
        Place legend on axis subplots

    ax : matplotlib axis object, default None
    style : list or dict
        matplotlib line style per column
    kind : {'line', 'bar', 'barh', 'kde', 'density', 'area', scatter', 'hexbin'}
        line : line plot
        bar : vertical bar plot
        barh : horizontal bar plot
        kde/density : Kernel Density Estimation plot
        area : area plot
        scatter : scatter plot
        hexbin : hexbin plot
    logx : boolean, default False
        Use log scaling on x axis
    logy : boolean, default False
        Use log scaling on y axis
    loglog : boolean, default False
        Use log scaling on both x and y axes
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
    position : float
        Specify relative alignments for bar plot layout.
        From 0 (left/bottom-end) to 1 (right/top-end). Default is 0.5 (center)
    kwds : keywords
        Options to pass to matplotlib plotting method

    Returns
    -------
    ax_or_axes : matplotlib.AxesSubplot or list of them

    Notes
    -----

    If `kind`='hexbin', you can control the size of the bins with the
    `gridsize` argument. By default, a histogram of the counts around each
    `(x, y)` point is computed. You can specify alternative aggregations
    by passing values to the `C` and `reduce_C_function` arguments.
    `C` specifies the value at each `(x, y)` point and `reduce_C_function`
    is a function of one argument that reduces all the values in a bin to
    a single number (e.g. `mean`, `max`, `sum`, `std`).
    """

    kind = _get_standard_kind(kind.lower().strip())
    if kind in _all_kinds:
        klass = _plot_klass[kind]
    else:
        raise ValueError('Invalid chart type given %s' % kind)

    if kind in _dataframe_kinds:
        plot_obj = klass(frame,  x=x, y=y, kind=kind, subplots=subplots,
                         rot=rot,legend=legend, ax=ax, style=style,
                         fontsize=fontsize, use_index=use_index, sharex=sharex,
                         sharey=sharey, xticks=xticks, yticks=yticks,
                         xlim=xlim, ylim=ylim, title=title, grid=grid,
                         figsize=figsize, logx=logx, logy=logy,
                         sort_columns=sort_columns, secondary_y=secondary_y,
                         **kwds)
    elif kind in _series_kinds:
        if y is None and subplots is False:
            msg = "{0} requires either y column or 'subplots=True'"
            raise ValueError(msg.format(kind))
        elif y is not None:
            if com.is_integer(y) and not frame.columns.holds_integer():
                y = frame.columns[y]
            frame = frame[y]  # converted to series actually
            frame.index.name = y

        plot_obj = klass(frame,  kind=kind, subplots=subplots,
                         rot=rot,legend=legend, ax=ax, style=style,
                         fontsize=fontsize, use_index=use_index, sharex=sharex,
                         sharey=sharey, xticks=xticks, yticks=yticks,
                         xlim=xlim, ylim=ylim, title=title, grid=grid,
                         figsize=figsize,
                         sort_columns=sort_columns,
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

            for kw in ['xerr', 'yerr']:
                if (kw in kwds) and \
                (isinstance(kwds[kw], string_types) or com.is_integer(kwds[kw])):
                    try:
                        kwds[kw] = frame[kwds[kw]]
                    except (IndexError, KeyError, TypeError):
                        pass

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
    kind : {'line', 'bar', 'barh', 'kde', 'density', 'area'}
        line : line plot
        bar : vertical bar plot
        barh : horizontal bar plot
        kde/density : Kernel Density Estimation plot
        area : area plot
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
        Use log scaling on x axis
    logy : boolean, default False
        Use log scaling on y axis
    loglog : boolean, default False
        Use log scaling on both x and y axes
    secondary_y : boolean or sequence of ints, default False
        If True then y-axis will be on the right
    figsize : a tuple (width, height) in inches
    position : float
        Specify relative alignments for bar plot layout.
        From 0 (left/bottom-end) to 1 (right/top-end). Default is 0.5 (center)
    kwds : keywords
        Options to pass to matplotlib plotting method

    Notes
    -----
    See matplotlib documentation online for more on this subject
    """

    kind = _get_standard_kind(kind.lower().strip())
    if kind in _common_kinds or kind in _series_kinds:
        klass = _plot_klass[kind]
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


_shared_docs['boxplot'] = """
    Make a box plot from DataFrame column optionally grouped by some columns or
    other inputs

    Parameters
    ----------
    data : the pandas object holding the data
    column : column name or list of names, or vector
        Can be any valid input to groupby
    by : string or sequence
        Column in the DataFrame to group by
    ax : Matplotlib axes object, optional
    fontsize : int or string
    rot : label rotation angle
    figsize : A tuple (width, height) in inches
    grid : Setting this to True will show the grid
    layout : tuple (optional)
        (rows, columns) for the layout of the plot
    return_type : {'axes', 'dict', 'both'}, default 'dict'
        The kind of object to return. 'dict' returns a dictionary
        whose values are the matplotlib Lines of the boxplot;
        'axes' returns the matplotlib axes the boxplot is drawn on;
        'both' returns a namedtuple with the axes and dict.

        When grouping with ``by``, a dict mapping columns to ``return_type``
        is returned.

    kwds : other plotting keyword arguments to be passed to matplotlib boxplot
           function

    Returns
    -------
    lines : dict
    ax : matplotlib Axes
    (ax, lines): namedtuple

    Notes
    -----
    Use ``return_type='dict'`` when you want to tweak the appearance
    of the lines after plotting. In this case a dict containing the Lines
    making up the boxes, caps, fliers, medians, and whiskers is returned.
    """


@Appender(_shared_docs['boxplot'] % _shared_doc_kwargs)
def boxplot(data, column=None, by=None, ax=None, fontsize=None,
            rot=0, grid=True, figsize=None, layout=None, return_type=None,
            **kwds):

    # validate return_type:
    valid_types = (None, 'axes', 'dict', 'both')
    if return_type not in valid_types:
        raise ValueError("return_type")


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
        return bp

    colors = _get_colors()
    if column is None:
        columns = None
    else:
        if isinstance(column, (list, tuple)):
            columns = column
        else:
            columns = [column]

    BP = namedtuple("Boxplot", ['ax', 'lines'])  # namedtuple to hold results

    if by is not None:
        fig, axes, d = _grouped_plot_by_column(plot_group, data, columns=columns,
                                               by=by, grid=grid, figsize=figsize,
                                               ax=ax, layout=layout)

        # Return axes in multiplot case, maybe revisit later # 985
        if return_type is None:
            ret = axes
        if return_type == 'axes':
            ret = compat.OrderedDict()
            axes = _flatten(axes)[:len(d)]
            for k, ax in zip(d.keys(), axes):
                ret[k] = ax
        elif return_type == 'dict':
            ret = d
        elif return_type == 'both':
            ret = compat.OrderedDict()
            axes = _flatten(axes)[:len(d)]
            for (k, line), ax in zip(d.items(), axes):
                ret[k] = BP(ax=ax, lines=line)
    else:
        if layout is not None:
            raise ValueError("The 'layout' keyword is not supported when "
                             "'by' is None")
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

        ret = ax

        if return_type is None:
            msg = ("\nThe default value for 'return_type' will change to "
                   "'axes' in a future release.\n To use the future behavior "
                   "now, set return_type='axes'.\n To keep the previous "
                   "behavior and silence this warning, set "
                   "return_type='dict'.")
            warnings.warn(msg, FutureWarning)
            return_type = 'dict'
        if return_type == 'dict':
            ret = bp
        elif return_type == 'both':
            ret = BP(ax=ret, lines=bp)

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
               sharey=False, figsize=None, layout=None, bins=10, **kwds):
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
    bins: integer, default 10
        Number of histogram bins to be used
    kwds : other plotting keyword arguments
        To be passed to hist function
    """
    import matplotlib.pyplot as plt

    if by is not None:
        axes = grouped_hist(data, column=column, by=by, ax=ax, grid=grid, figsize=figsize,
                            sharex=sharex, sharey=sharey, layout=layout, bins=bins,
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

    if column is not None:
        if not isinstance(column, (list, np.ndarray)):
            column = [column]
        data = data[column]
    data = data._get_numeric_data()
    naxes = len(data.columns)

    nrows, ncols = _get_layout(naxes, layout=layout)
    fig, axes = _subplots(nrows=nrows, ncols=ncols, naxes=naxes, ax=ax, squeeze=False,
                          sharex=sharex, sharey=sharey, figsize=figsize)

    for i, col in enumerate(com._try_sort(data.columns)):
        ax = axes[i // ncols, i % ncols]
        ax.xaxis.set_visible(True)
        ax.yaxis.set_visible(True)
        ax.hist(data[col].dropna().values, bins=bins, **kwds)
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

    fig.subplots_adjust(wspace=0.3, hspace=0.3)

    return axes


def hist_series(self, by=None, ax=None, grid=True, xlabelsize=None,
                xrot=None, ylabelsize=None, yrot=None, figsize=None, bins=10, **kwds):
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
    bins: integer, default 10
        Number of histogram bins to be used
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

        ax.hist(values, bins=bins, **kwds)
        ax.grid(grid)
        axes = np.array([ax])
    else:
        if 'figure' in kwds:
            raise ValueError("Cannot pass 'figure' when using the "
                             "'by' argument, since a new 'Figure' instance "
                             "will be created")
        axes = grouped_hist(self, by=by, ax=ax, grid=grid, figsize=figsize,
                            bins=bins, **kwds)

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


def boxplot_frame_groupby(grouped, subplots=True, column=None, fontsize=None,
                          rot=0, grid=True, figsize=None, layout=None, **kwds):
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
    layout : tuple (optional)
        (rows, columns) for the layout of the plot
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
        naxes = len(grouped)
        nrows, ncols = _get_layout(naxes, layout=layout)
        _, axes = _subplots(nrows=nrows, ncols=ncols, naxes=naxes, squeeze=False,
                            sharex=False, sharey=True)
        axes = _flatten(axes)

        ret = compat.OrderedDict()
        for (key, group), ax in zip(grouped, axes):
            d = group.boxplot(ax=ax, column=column, fontsize=fontsize,
                              rot=rot, grid=grid, **kwds)
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
                         grid=grid, figsize=figsize, layout=layout, **kwds)
    return ret


def _grouped_plot(plotf, data, column=None, by=None, numeric_only=True,
                  figsize=None, sharex=True, sharey=True, layout=None,
                  rot=0, ax=None, **kwargs):
    from pandas import DataFrame

    # allow to specify mpl default with 'default'
    if figsize is None or figsize == 'default':
        figsize = (10, 5)               # our default

    grouped = data.groupby(by)
    if column is not None:
        grouped = grouped[column]

    naxes = len(grouped)
    nrows, ncols = _get_layout(naxes, layout=layout)
    if figsize is None:
        # our favorite default beating matplotlib's idea of the
        # default size
        figsize = (10, 5)
    fig, axes = _subplots(nrows=nrows, ncols=ncols, naxes=naxes,
                          figsize=figsize, sharex=sharex, sharey=sharey, ax=ax)

    ravel_axes = _flatten(axes)

    for i, (key, group) in enumerate(grouped):
        ax = ravel_axes[i]
        if numeric_only and isinstance(group, DataFrame):
            group = group._get_numeric_data()
        plotf(group, ax, **kwargs)
        ax.set_title(com.pprint_thing(key))

    return fig, axes


def _grouped_plot_by_column(plotf, data, columns=None, by=None,
                            numeric_only=True, grid=False,
                            figsize=None, ax=None, layout=None, **kwargs):
    from pandas.core.frame import DataFrame

    grouped = data.groupby(by)
    if columns is None:
        if not isinstance(by, (list, tuple)):
            by = [by]
        columns = data._get_numeric_data().columns - by
    naxes = len(columns)

    if ax is None:
        nrows, ncols = _get_layout(naxes, layout=layout)
        fig, axes = _subplots(nrows=nrows, ncols=ncols, naxes=naxes,
                              sharex=True, sharey=True,
                              figsize=figsize, ax=ax)
    else:
        if naxes > 1:
            raise ValueError("Using an existing axis is not supported when plotting multiple columns.")
        fig = ax.get_figure()
        axes = ax.get_axes()

    ravel_axes = _flatten(axes)

    out_dict = compat.OrderedDict()
    for i, col in enumerate(columns):
        ax = ravel_axes[i]
        gp_col = grouped[col]
        re_plotf = plotf(gp_col, ax, **kwargs)
        ax.set_title(col)
        ax.set_xlabel(com.pprint_thing(by))
        ax.grid(grid)
        out_dict[col] = re_plotf

    byline = by[0] if len(by) == 1 else by
    fig.suptitle('Boxplot grouped by %s' % byline)

    return fig, axes, out_dict


def table(ax, data, rowLabels=None, colLabels=None,
          **kwargs):

    """
    Helper function to convert DataFrame and Series to matplotlib.table

    Parameters
    ----------
    `ax`: Matplotlib axes object
    `data`: DataFrame or Series
        data for table contents
    `kwargs`: keywords, optional
        keyword arguments which passed to matplotlib.table.table.
        If `rowLabels` or `colLabels` is not specified, data index or column name will be used.

    Returns
    -------
    matplotlib table object
    """
    from pandas import DataFrame
    if isinstance(data, Series):
        data = DataFrame(data, columns=[data.name])
    elif isinstance(data, DataFrame):
        pass
    else:
        raise ValueError('Input data must be dataframe or series')

    if rowLabels is None:
        rowLabels = data.index

    if colLabels is None:
        colLabels = data.columns

    cellText = data.values

    import matplotlib.table
    table = matplotlib.table.table(ax, cellText=cellText,
        rowLabels=rowLabels, colLabels=colLabels, **kwargs)
    return table


def _get_layout(nplots, layout=None):
    if layout is not None:
        if not isinstance(layout, (tuple, list)) or len(layout) != 2:
            raise ValueError('Layout must be a tuple of (rows, columns)')

        nrows, ncols = layout
        if nrows * ncols < nplots:
            raise ValueError('Layout of %sx%s must be larger than required size %s' %
                (nrows, ncols, nplots))

        return layout

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


def _subplots(nrows=1, ncols=1, naxes=None, sharex=False, sharey=False, squeeze=True,
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

    naxes : int
      Number of required axes. Exceeded axes are set invisible. Default is nrows * ncols.

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

    if naxes is None:
        naxes = nrows * ncols
    elif nplots < naxes:
        raise ValueError("naxes {0} is larger than layour size defined by nrows * ncols".format(naxes))

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

    if naxes != nplots:
        for ax in axarr[naxes:]:
            ax.set_visible(False)

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


def _flatten(axes):
    if not com.is_list_like(axes):
        axes = [axes]
    elif isinstance(axes, np.ndarray):
        axes = axes.ravel()
    return axes


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
