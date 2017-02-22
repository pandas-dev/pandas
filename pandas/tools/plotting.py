# being a bit too dynamic
# pylint: disable=E1101
from __future__ import division

import warnings
import re
from math import ceil
from collections import namedtuple
from contextlib import contextmanager
from distutils.version import LooseVersion

import numpy as np

from pandas.types.common import (is_list_like,
                                 is_integer,
                                 is_number,
                                 is_hashable,
                                 is_iterator)
from pandas.types.missing import isnull, notnull

from pandas.util.decorators import cache_readonly, deprecate_kwarg
from pandas.core.base import PandasObject

from pandas.core.common import AbstractMethodError, _try_sort
from pandas.core.generic import _shared_docs, _shared_doc_kwargs
from pandas.core.index import Index, MultiIndex
from pandas.core.series import Series, remove_na
from pandas.tseries.period import PeriodIndex
from pandas.compat import range, lrange, lmap, map, zip, string_types
import pandas.compat as compat
from pandas.formats.printing import pprint_thing
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


def _mpl_le_1_2_1():
    try:
        import matplotlib as mpl
        return (str(mpl.__version__) <= LooseVersion('1.2.1') and
                str(mpl.__version__)[0] != '0')
    except ImportError:
        return False


def _mpl_ge_1_3_1():
    try:
        import matplotlib
        # The or v[0] == '0' is because their versioneer is
        # messed up on dev
        return (matplotlib.__version__ >= LooseVersion('1.3.1') or
                matplotlib.__version__[0] == '0')
    except ImportError:
        return False


def _mpl_ge_1_4_0():
    try:
        import matplotlib
        return (matplotlib.__version__ >= LooseVersion('1.4') or
                matplotlib.__version__[0] == '0')
    except ImportError:
        return False


def _mpl_ge_1_5_0():
    try:
        import matplotlib
        return (matplotlib.__version__ >= LooseVersion('1.5') or
                matplotlib.__version__[0] == '0')
    except ImportError:
        return False


def _mpl_ge_2_0_0():
    try:
        import matplotlib
        return matplotlib.__version__ >= LooseVersion('2.0')
    except ImportError:
        return False


if _mpl_ge_1_5_0():
    # Compat with mp 1.5, which uses cycler.
    import cycler
    colors = mpl_stylesheet.pop('axes.color_cycle')
    mpl_stylesheet['axes.prop_cycle'] = cycler.cycler('color', colors)


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
        colors = list(color) if is_list_like(color) else color
    else:
        if color_type == 'default':
            # need to call list() on the result to copy so we don't
            # modify the global rcParams below
            try:
                colors = [c['color']
                          for c in list(plt.rcParams['axes.prop_cycle'])]
            except KeyError:
                colors = list(plt.rcParams.get('axes.color_cycle',
                                               list('bgrcmyk')))
            if isinstance(colors, compat.string_types):
                colors = list(colors)
        elif color_type == 'random':
            import random

            def random_color(column):
                random.seed(column)
                return [random.random() for _ in range(3)]

            colors = lmap(random_color, lrange(num_colors))
        else:
            raise ValueError("color_type must be either 'default' or 'random'")

    if isinstance(colors, compat.string_types):
        import matplotlib.colors
        conv = matplotlib.colors.ColorConverter()

        def _maybe_valid_colors(colors):
            try:
                [conv.to_rgba(c) for c in colors]
                return True
            except ValueError:
                return False

        # check whether the string can be convertable to single color
        maybe_single_color = _maybe_valid_colors([colors])
        # check whether each character can be convertable to colors
        maybe_color_cycle = _maybe_valid_colors(list(colors))
        if maybe_single_color and maybe_color_cycle and len(colors) > 1:
            msg = ("'{0}' can be parsed as both single color and "
                   "color cycle. Specify each color using a list "
                   "like ['{0}'] or {1}")
            raise ValueError(msg.format(colors, list(colors)))
        elif maybe_single_color:
            colors = [colors]
        else:
            # ``colors`` is regarded as color cycle.
            # mpl will raise error any of them is invalid
            pass

    if len(colors) != num_colors:
        multiple = num_colors // len(colors) - 1
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

    df = frame._get_numeric_data()
    n = df.columns.size
    naxes = n * n
    fig, axes = _subplots(naxes=naxes, figsize=figsize, ax=ax,
                          squeeze=False)

    # no gaps between subplots
    fig.subplots_adjust(wspace=0, hspace=0)

    mask = notnull(df)

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
        boundaries_list.append((rmin_ - rdelta_ext, rmax_ + rdelta_ext))

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

            ax.set_xlabel(b)
            ax.set_ylabel(a)

            if j != 0:
                ax.yaxis.set_visible(False)
            if i != n - 1:
                ax.xaxis.set_visible(False)

    if len(df.columns) > 1:
        lim1 = boundaries_list[0]
        locs = axes[0][1].yaxis.get_majorticklocs()
        locs = locs[(lim1[0] <= locs) & (locs <= lim1[1])]
        adj = (locs - lim1[0]) / (lim1[1] - lim1[0])

        lim0 = axes[0][0].get_ylim()
        adj = adj * (lim0[1] - lim0[0]) + lim0[0]
        axes[0][0].yaxis.set_ticks(adj)

        if np.all(locs == locs.astype(int)):
            # if all ticks are int
            locs = locs.astype(int)
        axes[0][0].yaxis.set_ticklabels(locs)

    _set_ticks_props(axes, xlabelsize=8, xrot=90, ylabelsize=8, yrot=0)

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

    m = len(frame.columns) - 1
    s = np.array([(np.cos(t), np.sin(t))
                  for t in [2.0 * np.pi * (i / float(m))
                            for i in range(m)]])

    for i in range(n):
        row = df.iloc[i].values
        row_ = np.repeat(np.expand_dims(row, axis=1), 2, axis=1)
        y = (s * row_).sum(axis=0) / row.sum()
        kls = class_col.iat[i]
        to_plot[kls][0].append(y[0])
        to_plot[kls][1].append(y[1])

    for i, kls in enumerate(classes):
        ax.scatter(to_plot[kls][0], to_plot[kls][1], color=colors[i],
                   label=pprint_thing(kls), **kwds)
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
    Generates a matplotlib plot of Andrews curves, for visualising clusters of
    multivariate data.

    Andrews curves have the functional form:

    f(t) = x_1/sqrt(2) + x_2 sin(t) + x_3 cos(t) +
           x_4 sin(2t) + x_5 cos(2t) + ...

    Where x coefficients correspond to the values of each dimension and t is
    linearly spaced between -pi and +pi. Each row of frame then corresponds to
    a single curve.

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
    from math import sqrt, pi
    import matplotlib.pyplot as plt

    def function(amplitudes):
        def f(t):
            x1 = amplitudes[0]
            result = x1 / sqrt(2.0)

            # Take the rest of the coefficients and resize them
            # appropriately. Take a copy of amplitudes as otherwise numpy
            # deletes the element from amplitudes itself.
            coeffs = np.delete(np.copy(amplitudes), 0)
            coeffs.resize(int((coeffs.size + 1) / 2), 2)

            # Generate the harmonics and arguments for the sin and cos
            # functions.
            harmonics = np.arange(0, coeffs.shape[0]) + 1
            trig_args = np.outer(harmonics, t)

            result += np.sum(coeffs[:, 0, np.newaxis] * np.sin(trig_args) +
                             coeffs[:, 1, np.newaxis] * np.cos(trig_args),
                             axis=0)
            return result
        return f

    n = len(frame)
    class_col = frame[class_column]
    classes = frame[class_column].drop_duplicates()
    df = frame.drop(class_column, axis=1)
    t = np.linspace(-pi, pi, samples)
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
        y = f(t)
        kls = class_col.iat[i]
        label = pprint_thing(kls)
        if label not in used_legends:
            used_legends.add(label)
            ax.plot(t, y, color=colors[kls], label=label, **kwds)
        else:
            ax.plot(t, y, color=colors[kls], **kwds)

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
@deprecate_kwarg(old_arg_name='data', new_arg_name='frame', stacklevel=3)
def parallel_coordinates(frame, class_column, cols=None, ax=None, color=None,
                         use_columns=False, xticks=None, colormap=None,
                         axvlines=True, axvlines_kwds=None, **kwds):
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
    axvlines: bool, optional
        If true, vertical lines will be added at each xtick
    axvlines_kwds: keywords, optional
        Options to be passed to axvline method for vertical lines
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
    >>> df = read_csv('https://raw.github.com/pandas-dev/pandas/master'
                      '/pandas/tests/data/iris.csv')
    >>> parallel_coordinates(df, 'Name', color=('#556270',
                             '#4ECDC4', '#C7F464'))
    >>> plt.show()
    """
    if axvlines_kwds is None:
        axvlines_kwds = {'linewidth': 1, 'color': 'black'}
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
        label = pprint_thing(kls)
        if label not in used_legends:
            used_legends.add(label)
            ax.plot(x, y, color=colors[kls], label=label, **kwds)
        else:
            ax.plot(x, y, color=colors[kls], **kwds)

    if axvlines:
        for i in x:
            ax.axvline(i, **axvlines_kwds)

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
        return ((data[:n - h] - mean) *
                (data[h:] - mean)).sum() / float(n) / c0
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

    @property
    def _kind(self):
        """Specify kind str. Must be overridden in child class"""
        raise NotImplementedError

    _layout_type = 'vertical'
    _default_rot = 0
    orientation = None
    _pop_attributes = ['label', 'style', 'logy', 'logx', 'loglog',
                       'mark_right', 'stacked']
    _attr_defaults = {'logy': False, 'logx': False, 'loglog': False,
                      'mark_right': True, 'stacked': False}

    def __init__(self, data, kind=None, by=None, subplots=False, sharex=None,
                 sharey=False, use_index=True,
                 figsize=None, grid=None, legend=True, rot=None,
                 ax=None, fig=None, title=None, xlim=None, ylim=None,
                 xticks=None, yticks=None,
                 sort_columns=False, fontsize=None,
                 secondary_y=False, colormap=None,
                 table=False, layout=None, **kwds):

        self.data = data
        self.by = by

        self.kind = kind

        self.sort_columns = sort_columns

        self.subplots = subplots

        if sharex is None:
            if ax is None:
                self.sharex = True
            else:
                # if we get an axis, the users should do the visibility
                # setting...
                self.sharex = False
        else:
            self.sharex = sharex

        self.sharey = sharey
        self.figsize = figsize
        self.layout = layout

        self.xticks = xticks
        self.yticks = yticks
        self.xlim = xlim
        self.ylim = ylim
        self.title = title
        self.use_index = use_index

        self.fontsize = fontsize

        if rot is not None:
            self.rot = rot
            # need to know for format_date_labels since it's rotated to 30 by
            # default
            self._rot_set = True
        else:
            self._rot_set = False
            self.rot = self._default_rot

        if grid is None:
            grid = False if secondary_y else self.plt.rcParams['axes.grid']

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

        if not isinstance(secondary_y, (bool, tuple, list, np.ndarray, Index)):
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
        if 'color' not in self.kwds and 'colors' in self.kwds:
            warnings.warn(("'colors' is being deprecated. Please use 'color'"
                           "instead of 'colors'"))
            colors = self.kwds.pop('colors')
            self.kwds['color'] = colors

        if ('color' in self.kwds and self.nseries == 1):
            # support series.plot(color='green')
            self.kwds['color'] = [self.kwds['color']]

        if ('color' in self.kwds or 'colors' in self.kwds) and \
                self.colormap is not None:
            warnings.warn("'color' and 'colormap' cannot be used "
                          "simultaneously. Using 'color'")

        if 'color' in self.kwds and self.style is not None:
            if is_list_like(self.style):
                styles = self.style
            else:
                styles = [self.style]
            # need only a single match
            for s in styles:
                if re.match('^[a-z]+?', s) is not None:
                    raise ValueError(
                        "Cannot pass 'style' string with a color "
                        "symbol and 'color' keyword argument. Please"
                        " use one or the other or pass 'style' "
                        "without a color symbol")

    def _iter_data(self, data=None, keep_index=False, fillna=None):
        if data is None:
            data = self.data
        if fillna is not None:
            data = data.fillna(fillna)

        # TODO: unused?
        # if self.sort_columns:
        #     columns = _try_sort(data.columns)
        # else:
        #     columns = data.columns

        for col, values in data.iteritems():
            if keep_index is True:
                yield col, values
            else:
                yield col, values.values

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
        self._adorn_subplots()

        for ax in self.axes:
            self._post_plot_logic_common(ax, self.data)
            self._post_plot_logic(ax, self.data)

    def _args_adjust(self):
        pass

    def _has_plotted_object(self, ax):
        """check whether ax has data"""
        return (len(ax.lines) != 0 or
                len(ax.artists) != 0 or
                len(ax.containers) != 0)

    def _maybe_right_yaxis(self, ax, axes_num):
        if not self.on_right(axes_num):
            # secondary axes may be passed via ax kw
            return self._get_ax_layer(ax)

        if hasattr(ax, 'right_ax'):
            # if it has right_ax proparty, ``ax`` must be left axes
            return ax.right_ax
        elif hasattr(ax, 'left_ax'):
            # if it has left_ax proparty, ``ax`` must be right axes
            return ax
        else:
            # otherwise, create twin axes
            orig_ax, new_ax = ax, ax.twinx()
            # TODO: use Matplotlib public API when available
            new_ax._get_lines = orig_ax._get_lines
            new_ax._get_patches_for_fill = orig_ax._get_patches_for_fill
            orig_ax.right_ax, new_ax.left_ax = new_ax, orig_ax

            if not self._has_plotted_object(orig_ax):  # no data on left y
                orig_ax.get_yaxis().set_visible(False)
            return new_ax

    def _setup_subplots(self):
        if self.subplots:
            fig, axes = _subplots(naxes=self.nseries,
                                  sharex=self.sharex, sharey=self.sharey,
                                  figsize=self.figsize, ax=self.ax,
                                  layout=self.layout,
                                  layout_type=self._layout_type)
        else:
            if self.ax is None:
                fig = self.plt.figure(figsize=self.figsize)
                axes = fig.add_subplot(111)
            else:
                fig = self.ax.get_figure()
                if self.figsize is not None:
                    fig.set_size_inches(self.figsize)
                axes = self.ax

        axes = _flatten(axes)

        if self.logx or self.loglog:
            [a.set_xscale('log') for a in axes]
        if self.logy or self.loglog:
            [a.set_yscale('log') for a in axes]

        self.fig = fig
        self.axes = axes

    @property
    def result(self):
        """
        Return result axes
        """
        if self.subplots:
            if self.layout is not None and not is_list_like(self.ax):
                return self.axes.reshape(*self.layout)
            else:
                return self.axes
        else:
            sec_true = isinstance(self.secondary_y, bool) and self.secondary_y
            all_sec = (is_list_like(self.secondary_y) and
                       len(self.secondary_y) == self.nseries)
            if (sec_true or all_sec):
                # if all data is plotted on secondary, return right axes
                return self._get_ax_layer(self.axes[0], primary=False)
            else:
                return self.axes[0]

    def _compute_plot_data(self):
        data = self.data

        if isinstance(data, Series):
            label = self.label
            if label is None and data.name is None:
                label = 'None'
            data = data.to_frame(name=label)

        numeric_data = data._convert(datetime=True)._get_numeric_data()

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
        raise AbstractMethodError(self)

    def _add_table(self):
        if self.table is False:
            return
        elif self.table is True:
            data = self.data.transpose()
        else:
            data = self.table
        ax = self._get_ax(0)
        table(ax, data)

    def _post_plot_logic_common(self, ax, data):
        """Common post process for each axes"""
        labels = [pprint_thing(key) for key in data.index]
        labels = dict(zip(range(len(data.index)), labels))

        if self.orientation == 'vertical' or self.orientation is None:
            if self._need_to_set_index:
                xticklabels = [labels.get(x, '') for x in ax.get_xticks()]
                ax.set_xticklabels(xticklabels)
            self._apply_axis_properties(ax.xaxis, rot=self.rot,
                                        fontsize=self.fontsize)
            self._apply_axis_properties(ax.yaxis, fontsize=self.fontsize)
        elif self.orientation == 'horizontal':
            if self._need_to_set_index:
                yticklabels = [labels.get(y, '') for y in ax.get_yticks()]
                ax.set_yticklabels(yticklabels)
            self._apply_axis_properties(ax.yaxis, rot=self.rot,
                                        fontsize=self.fontsize)
            self._apply_axis_properties(ax.xaxis, fontsize=self.fontsize)
        else:  # pragma no cover
            raise ValueError

    def _post_plot_logic(self, ax, data):
        """Post process for each axes. Overridden in child classes"""
        pass

    def _adorn_subplots(self):
        """Common post process unrelated to data"""
        if len(self.axes) > 0:
            all_axes = self._get_subplots()
            nrows, ncols = self._get_axes_layout()
            _handle_shared_axes(axarr=all_axes, nplots=len(all_axes),
                                naxes=nrows * ncols, nrows=nrows,
                                ncols=ncols, sharex=self.sharex,
                                sharey=self.sharey)

        for ax in self.axes:
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
                if is_list_like(self.title):
                    if len(self.title) != self.nseries:
                        msg = ('The length of `title` must equal the number '
                               'of columns if using `title` of type `list` '
                               'and `subplots=True`.\n'
                               'length of title = {}\n'
                               'number of columns = {}').format(
                            len(self.title), self.nseries)
                        raise ValueError(msg)

                    for (ax, title) in zip(self.axes, self.title):
                        ax.set_title(title)
                else:
                    self.fig.suptitle(self.title)
            else:
                if is_list_like(self.title):
                    msg = ('Using `title` of type `list` is not supported '
                           'unless `subplots=True` is passed')
                    raise ValueError(msg)
                self.axes[0].set_title(self.title)

    def _apply_axis_properties(self, axis, rot=None, fontsize=None):
        labels = axis.get_majorticklabels() + axis.get_minorticklabels()
        for label in labels:
            if rot is not None:
                label.set_rotation(rot)
            if fontsize is not None:
                label.set_fontsize(fontsize)

    @property
    def legend_title(self):
        if not isinstance(self.data.columns, MultiIndex):
            name = self.data.columns.name
            if name is not None:
                name = pprint_thing(name)
            return name
        else:
            stringified = map(pprint_thing,
                              self.data.columns.names)
            return ','.join(stringified)

    def _add_legend_handle(self, handle, label, index=None):
        if label is not None:
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
            if leg is not None:
                title = leg.get_title().get_text()
                handles = leg.legendHandles
                labels = [x.get_text() for x in leg.get_texts()]

            if self.legend:
                if self.legend == 'reverse':
                    self.legend_handles = reversed(self.legend_handles)
                    self.legend_labels = reversed(self.legend_labels)

                handles += self.legend_handles
                labels += self.legend_labels
                if self.legend_title is not None:
                    title = self.legend_title

            if len(handles) > 0:
                ax.legend(handles, labels, loc='best', title=title)

        elif self.subplots and self.legend:
            for ax in self.axes:
                if ax.get_visible():
                    ax.legend(loc='best')

    def _get_ax_legend(self, ax):
        leg = ax.get_legend()
        other_ax = (getattr(ax, 'left_ax', None) or
                    getattr(ax, 'right_ax', None))
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

    @staticmethod
    def mpl_ge_1_3_1():
        return _mpl_ge_1_3_1()

    @staticmethod
    def mpl_ge_1_5_0():
        return _mpl_ge_1_5_0()

    _need_to_set_index = False

    def _get_xticks(self, convert_period=False):
        index = self.data.index
        is_datetype = index.inferred_type in ('datetime', 'date',
                                              'datetime64', 'time')

        if self.use_index:
            if convert_period and isinstance(index, PeriodIndex):
                self.data = self.data.reindex(index=index.sort_values())
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

    @classmethod
    def _plot(cls, ax, x, y, style=None, is_errorbar=False, **kwds):
        mask = isnull(y)
        if mask.any():
            y = np.ma.array(y)
            y = np.ma.masked_where(mask, y)

        if isinstance(x, Index):
            x = x._mpl_repr()

        if is_errorbar:
            if 'xerr' in kwds:
                kwds['xerr'] = np.array(kwds.get('xerr'))
            if 'yerr' in kwds:
                kwds['yerr'] = np.array(kwds.get('yerr'))
            return ax.errorbar(x, y, **kwds)
        else:
            # prevent style kwarg from going to errorbar, where it is
            # unsupported
            if style is not None:
                args = (x, y, style)
            else:
                args = (x, y)
            return ax.plot(*args, **kwds)

    def _get_index_name(self):
        if isinstance(self.data.index, MultiIndex):
            name = self.data.index.names
            if any(x is not None for x in name):
                name = ','.join([pprint_thing(x) for x in name])
            else:
                name = None
        else:
            name = self.data.index.name
            if name is not None:
                name = pprint_thing(name)

        return name

    @classmethod
    def _get_ax_layer(cls, ax, primary=True):
        """get left (primary) or right (secondary) axes"""
        if primary:
            return getattr(ax, 'left_ax', ax)
        else:
            return getattr(ax, 'right_ax', ax)

    def _get_ax(self, i):
        # get the twinx ax if appropriate
        if self.subplots:
            ax = self.axes[i]
            ax = self._maybe_right_yaxis(ax, i)
            self.axes[i] = ax
        else:
            ax = self.axes[0]
            ax = self._maybe_right_yaxis(ax, i)

        ax.get_yaxis().set_visible(True)
        return ax

    def on_right(self, i):
        if isinstance(self.secondary_y, bool):
            return self.secondary_y

        if isinstance(self.secondary_y, (tuple, list, np.ndarray, Index)):
            return self.data.columns[i] in self.secondary_y

    def _apply_style_colors(self, colors, kwds, col_num, label):
        """
        Manage style and color based on column number and its label.
        Returns tuple of appropriate style and kwds which "color" may be added.
        """
        style = None
        if self.style is not None:
            if isinstance(self.style, list):
                try:
                    style = self.style[col_num]
                except IndexError:
                    pass
            elif isinstance(self.style, dict):
                style = self.style.get(label, style)
            else:
                style = self.style

        has_color = 'color' in kwds or self.colormap is not None
        nocolor_style = style is None or re.match('[a-z]+', style) is None
        if (has_color or self.subplots) and nocolor_style:
            kwds['color'] = colors[col_num % len(colors)]
        return style, kwds

    def _get_colors(self, num_colors=None, color_kwds='color'):
        if num_colors is None:
            num_colors = self.nseries

        return _get_standard_colors(num_colors=num_colors,
                                    colormap=self.colormap,
                                    color=self.kwds.get(color_kwds))

    def _parse_errorbars(self, label, err):
        """
        Look for error keyword arguments and return the actual errorbar data
        or return the error DataFrame/dict

        Error bars can be specified in several ways:
            Series: the user provides a pandas.Series object of the same
                    length as the data
            ndarray: provides a np.ndarray of the same length as the data
            DataFrame/dict: error values are paired with keys matching the
                    key in the plotted DataFrame
            str: the name of the column within the plotted DataFrame
        """

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

        elif is_list_like(err):
            if is_iterator(err):
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

        elif is_number(err):
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

    def _get_subplots(self):
        from matplotlib.axes import Subplot
        return [ax for ax in self.axes[0].get_figure().get_axes()
                if isinstance(ax, Subplot)]

    def _get_axes_layout(self):
        axes = self._get_subplots()
        x_set = set()
        y_set = set()
        for ax in axes:
            # check axes coordinates to estimate layout
            points = ax.get_position().get_points()
            x_set.add(points[0][0])
            y_set.add(points[0][1])
        return (len(y_set), len(x_set))


class PlanePlot(MPLPlot):
    """
    Abstract class for plotting on plane, currently scatter and hexbin.
    """

    _layout_type = 'single'

    def __init__(self, data, x, y, **kwargs):
        MPLPlot.__init__(self, data, **kwargs)
        if x is None or y is None:
            raise ValueError(self._kind + ' requires and x and y column')
        if is_integer(x) and not self.data.columns.holds_integer():
            x = self.data.columns[x]
        if is_integer(y) and not self.data.columns.holds_integer():
            y = self.data.columns[y]
        self.x = x
        self.y = y

    @property
    def nseries(self):
        return 1

    def _post_plot_logic(self, ax, data):
        x, y = self.x, self.y
        ax.set_ylabel(pprint_thing(y))
        ax.set_xlabel(pprint_thing(x))


class ScatterPlot(PlanePlot):
    _kind = 'scatter'

    def __init__(self, data, x, y, s=None, c=None, **kwargs):
        if s is None:
            # hide the matplotlib default for size, in case we want to change
            # the handling of this argument later
            s = 20
        super(ScatterPlot, self).__init__(data, x, y, s=s, **kwargs)
        if is_integer(c) and not self.data.columns.holds_integer():
            c = self.data.columns[c]
        self.c = c

    def _make_plot(self):
        x, y, c, data = self.x, self.y, self.c, self.data
        ax = self.axes[0]

        c_is_column = is_hashable(c) and c in self.data.columns

        # plot a colorbar only if a colormap is provided or necessary
        cb = self.kwds.pop('colorbar', self.colormap or c_is_column)

        # pandas uses colormap, matplotlib uses cmap.
        cmap = self.colormap or 'Greys'
        cmap = self.plt.cm.get_cmap(cmap)
        color = self.kwds.pop("color", None)
        if c is not None and color is not None:
            raise TypeError('Specify exactly one of `c` and `color`')
        elif c is None and color is None:
            c_values = self.plt.rcParams['patch.facecolor']
        elif color is not None:
            c_values = color
        elif c_is_column:
            c_values = self.data[c].values
        else:
            c_values = c

        if self.legend and hasattr(self, 'label'):
            label = self.label
        else:
            label = None
        scatter = ax.scatter(data[x].values, data[y].values, c=c_values,
                             label=label, cmap=cmap, **self.kwds)
        if cb:
            img = ax.collections[0]
            kws = dict(ax=ax)
            if self.mpl_ge_1_3_1():
                kws['label'] = c if c_is_column else ''
            self.fig.colorbar(img, **kws)

        if label is not None:
            self._add_legend_handle(scatter, label)
        else:
            self.legend = False

        errors_x = self._get_errorbars(label=x, index=0, yerr=False)
        errors_y = self._get_errorbars(label=y, index=0, xerr=False)
        if len(errors_x) > 0 or len(errors_y) > 0:
            err_kwds = dict(errors_x, **errors_y)
            err_kwds['ecolor'] = scatter.get_facecolor()[0]
            ax.errorbar(data[x].values, data[y].values,
                        linestyle='none', **err_kwds)


class HexBinPlot(PlanePlot):
    _kind = 'hexbin'

    def __init__(self, data, x, y, C=None, **kwargs):
        super(HexBinPlot, self).__init__(data, x, y, **kwargs)
        if is_integer(C) and not self.data.columns.holds_integer():
            C = self.data.columns[C]
        self.C = C

    def _make_plot(self):
        x, y, data, C = self.x, self.y, self.data, self.C
        ax = self.axes[0]
        # pandas uses colormap, matplotlib uses cmap.
        cmap = self.colormap or 'BuGn'
        cmap = self.plt.cm.get_cmap(cmap)
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

    def _make_legend(self):
        pass


class LinePlot(MPLPlot):
    _kind = 'line'
    _default_rot = 0
    orientation = 'vertical'

    def __init__(self, data, **kwargs):
        MPLPlot.__init__(self, data, **kwargs)
        if self.stacked:
            self.data = self.data.fillna(value=0)
        self.x_compat = plot_params['x_compat']
        if 'x_compat' in self.kwds:
            self.x_compat = bool(self.kwds.pop('x_compat'))

    def _is_ts_plot(self):
        # this is slightly deceptive
        return not self.x_compat and self.use_index and self._use_dynamic_x()

    def _use_dynamic_x(self):
        from pandas.tseries.plotting import _use_dynamic_x
        return _use_dynamic_x(self._get_ax(0), self.data)

    def _make_plot(self):
        if self._is_ts_plot():
            from pandas.tseries.plotting import _maybe_convert_index
            data = _maybe_convert_index(self._get_ax(0), self.data)

            x = data.index      # dummy, not used
            plotf = self._ts_plot
            it = self._iter_data(data=data, keep_index=True)
        else:
            x = self._get_xticks(convert_period=True)
            plotf = self._plot
            it = self._iter_data()

        stacking_id = self._get_stacking_id()
        is_errorbar = any(e is not None for e in self.errors.values())

        colors = self._get_colors()
        for i, (label, y) in enumerate(it):
            ax = self._get_ax(i)
            kwds = self.kwds.copy()
            style, kwds = self._apply_style_colors(colors, kwds, i, label)

            errors = self._get_errorbars(label=label, index=i)
            kwds = dict(kwds, **errors)

            label = pprint_thing(label)  # .encode('utf-8')
            kwds['label'] = label

            newlines = plotf(ax, x, y, style=style, column_num=i,
                             stacking_id=stacking_id,
                             is_errorbar=is_errorbar,
                             **kwds)
            self._add_legend_handle(newlines[0], label, index=i)

            lines = _get_all_lines(ax)
            left, right = _get_xlim(lines)
            ax.set_xlim(left, right)

    @classmethod
    def _plot(cls, ax, x, y, style=None, column_num=None,
              stacking_id=None, **kwds):
        # column_num is used to get the target column from protf in line and
        # area plots
        if column_num == 0:
            cls._initialize_stacker(ax, stacking_id, len(y))
        y_values = cls._get_stacked_values(ax, stacking_id, y, kwds['label'])
        lines = MPLPlot._plot(ax, x, y_values, style=style, **kwds)
        cls._update_stacker(ax, stacking_id, y)
        return lines

    @classmethod
    def _ts_plot(cls, ax, x, data, style=None, **kwds):
        from pandas.tseries.plotting import (_maybe_resample,
                                             _decorate_axes,
                                             format_dateaxis)
        # accept x to be consistent with normal plot func,
        # x is not passed to tsplot as it uses data.index as x coordinate
        # column_num must be in kwds for stacking purpose
        freq, data = _maybe_resample(data, ax, kwds)

        # Set ax with freq info
        _decorate_axes(ax, freq, kwds)
        # digging deeper
        if hasattr(ax, 'left_ax'):
            _decorate_axes(ax.left_ax, freq, kwds)
        if hasattr(ax, 'right_ax'):
            _decorate_axes(ax.right_ax, freq, kwds)
        ax._plot_data.append((data, cls._kind, kwds))

        lines = cls._plot(ax, data.index, data.values, style=style, **kwds)
        # set date formatter, locators and rescale limits
        format_dateaxis(ax, ax.freq, data.index)
        return lines

    def _get_stacking_id(self):
        if self.stacked:
            return id(self.data)
        else:
            return None

    @classmethod
    def _initialize_stacker(cls, ax, stacking_id, n):
        if stacking_id is None:
            return
        if not hasattr(ax, '_stacker_pos_prior'):
            ax._stacker_pos_prior = {}
        if not hasattr(ax, '_stacker_neg_prior'):
            ax._stacker_neg_prior = {}
        ax._stacker_pos_prior[stacking_id] = np.zeros(n)
        ax._stacker_neg_prior[stacking_id] = np.zeros(n)

    @classmethod
    def _get_stacked_values(cls, ax, stacking_id, values, label):
        if stacking_id is None:
            return values
        if not hasattr(ax, '_stacker_pos_prior'):
            # stacker may not be initialized for subplots
            cls._initialize_stacker(ax, stacking_id, len(values))

        if (values >= 0).all():
            return ax._stacker_pos_prior[stacking_id] + values
        elif (values <= 0).all():
            return ax._stacker_neg_prior[stacking_id] + values

        raise ValueError('When stacked is True, each column must be either '
                         'all positive or negative.'
                         '{0} contains both positive and negative values'
                         .format(label))

    @classmethod
    def _update_stacker(cls, ax, stacking_id, values):
        if stacking_id is None:
            return
        if (values >= 0).all():
            ax._stacker_pos_prior[stacking_id] += values
        elif (values <= 0).all():
            ax._stacker_neg_prior[stacking_id] += values

    def _post_plot_logic(self, ax, data):
        condition = (not self._use_dynamic_x() and
                     data.index.is_all_dates and
                     not self.subplots or
                     (self.subplots and self.sharex))

        index_name = self._get_index_name()

        if condition:
            # irregular TS rotated 30 deg. by default
            # probably a better place to check / set this.
            if not self._rot_set:
                self.rot = 30
            format_date_labels(ax, rot=self.rot)

        if index_name is not None and self.use_index:
            ax.set_xlabel(index_name)


class AreaPlot(LinePlot):
    _kind = 'area'

    def __init__(self, data, **kwargs):
        kwargs.setdefault('stacked', True)
        data = data.fillna(value=0)
        LinePlot.__init__(self, data, **kwargs)

        if not self.stacked:
            # use smaller alpha to distinguish overlap
            self.kwds.setdefault('alpha', 0.5)

        if self.logy or self.loglog:
            raise ValueError("Log-y scales are not supported in area plot")

    @classmethod
    def _plot(cls, ax, x, y, style=None, column_num=None,
              stacking_id=None, is_errorbar=False, **kwds):

        if column_num == 0:
            cls._initialize_stacker(ax, stacking_id, len(y))
        y_values = cls._get_stacked_values(ax, stacking_id, y, kwds['label'])

        # need to remove label, because subplots uses mpl legend as it is
        line_kwds = kwds.copy()
        if cls.mpl_ge_1_5_0():
            line_kwds.pop('label')
        lines = MPLPlot._plot(ax, x, y_values, style=style, **line_kwds)

        # get data from the line to get coordinates for fill_between
        xdata, y_values = lines[0].get_data(orig=False)

        # unable to use ``_get_stacked_values`` here to get starting point
        if stacking_id is None:
            start = np.zeros(len(y))
        elif (y >= 0).all():
            start = ax._stacker_pos_prior[stacking_id]
        elif (y <= 0).all():
            start = ax._stacker_neg_prior[stacking_id]
        else:
            start = np.zeros(len(y))

        if 'color' not in kwds:
            kwds['color'] = lines[0].get_color()

        rect = ax.fill_between(xdata, start, y_values, **kwds)
        cls._update_stacker(ax, stacking_id, y)

        # LinePlot expects list of artists
        res = [rect] if cls.mpl_ge_1_5_0() else lines
        return res

    def _add_legend_handle(self, handle, label, index=None):
        if not self.mpl_ge_1_5_0():
            from matplotlib.patches import Rectangle
            # Because fill_between isn't supported in legend,
            # specifically add Rectangle handle here
            alpha = self.kwds.get('alpha', None)
            handle = Rectangle((0, 0), 1, 1, fc=handle.get_color(),
                               alpha=alpha)
        LinePlot._add_legend_handle(self, handle, label, index=index)

    def _post_plot_logic(self, ax, data):
        LinePlot._post_plot_logic(self, ax, data)

        if self.ylim is None:
            if (data >= 0).all().all():
                ax.set_ylim(0, None)
            elif (data <= 0).all().all():
                ax.set_ylim(None, 0)


class BarPlot(MPLPlot):
    _kind = 'bar'
    _default_rot = 90
    orientation = 'vertical'

    def __init__(self, data, **kwargs):
        self.bar_width = kwargs.pop('width', 0.5)
        pos = kwargs.pop('position', 0.5)
        kwargs.setdefault('align', 'center')
        self.tick_pos = np.arange(len(data))

        self.bottom = kwargs.pop('bottom', 0)
        self.left = kwargs.pop('left', 0)

        self.log = kwargs.pop('log', False)
        MPLPlot.__init__(self, data, **kwargs)

        if self.stacked or self.subplots:
            self.tickoffset = self.bar_width * pos
            if kwargs['align'] == 'edge':
                self.lim_offset = self.bar_width / 2
            else:
                self.lim_offset = 0
        else:
            if kwargs['align'] == 'edge':
                w = self.bar_width / self.nseries
                self.tickoffset = self.bar_width * (pos - 0.5) + w * 0.5
                self.lim_offset = w * 0.5
            else:
                self.tickoffset = self.bar_width * pos
                self.lim_offset = 0

        self.ax_pos = self.tick_pos - self.tickoffset

    def _args_adjust(self):
        if is_list_like(self.bottom):
            self.bottom = np.array(self.bottom)
        if is_list_like(self.left):
            self.left = np.array(self.left)

    @classmethod
    def _plot(cls, ax, x, y, w, start=0, log=False, **kwds):
        return ax.bar(x, y, w, bottom=start, log=log, **kwds)

    @property
    def _start_base(self):
        return self.bottom

    def _make_plot(self):
        import matplotlib as mpl

        colors = self._get_colors()
        ncolors = len(colors)

        pos_prior = neg_prior = np.zeros(len(self.data))
        K = self.nseries

        for i, (label, y) in enumerate(self._iter_data(fillna=0)):
            ax = self._get_ax(i)
            kwds = self.kwds.copy()
            kwds['color'] = colors[i % ncolors]

            errors = self._get_errorbars(label=label, index=i)
            kwds = dict(kwds, **errors)

            label = pprint_thing(label)

            if (('yerr' in kwds) or ('xerr' in kwds)) \
                    and (kwds.get('ecolor') is None):
                kwds['ecolor'] = mpl.rcParams['xtick.color']

            start = 0
            if self.log and (y >= 1).all():
                start = 1
            start = start + self._start_base

            if self.subplots:
                w = self.bar_width / 2
                rect = self._plot(ax, self.ax_pos + w, y, self.bar_width,
                                  start=start, label=label,
                                  log=self.log, **kwds)
                ax.set_title(label)
            elif self.stacked:
                mask = y > 0
                start = np.where(mask, pos_prior, neg_prior) + self._start_base
                w = self.bar_width / 2
                rect = self._plot(ax, self.ax_pos + w, y, self.bar_width,
                                  start=start, label=label,
                                  log=self.log, **kwds)
                pos_prior = pos_prior + np.where(mask, y, 0)
                neg_prior = neg_prior + np.where(mask, 0, y)
            else:
                w = self.bar_width / K
                rect = self._plot(ax, self.ax_pos + (i + 0.5) * w, y, w,
                                  start=start, label=label,
                                  log=self.log, **kwds)
            self._add_legend_handle(rect, label, index=i)

    def _post_plot_logic(self, ax, data):
        if self.use_index:
            str_index = [pprint_thing(key) for key in data.index]
        else:
            str_index = [pprint_thing(key) for key in range(data.shape[0])]
        name = self._get_index_name()

        s_edge = self.ax_pos[0] - 0.25 + self.lim_offset
        e_edge = self.ax_pos[-1] + 0.25 + self.bar_width + self.lim_offset

        self._decorate_ticks(ax, name, str_index, s_edge, e_edge)

    def _decorate_ticks(self, ax, name, ticklabels, start_edge, end_edge):
        ax.set_xlim((start_edge, end_edge))
        ax.set_xticks(self.tick_pos)
        ax.set_xticklabels(ticklabels)
        if name is not None and self.use_index:
            ax.set_xlabel(name)


class BarhPlot(BarPlot):
    _kind = 'barh'
    _default_rot = 0
    orientation = 'horizontal'

    @property
    def _start_base(self):
        return self.left

    @classmethod
    def _plot(cls, ax, x, y, w, start=0, log=False, **kwds):
        return ax.barh(x, y, w, left=start, log=log, **kwds)

    def _decorate_ticks(self, ax, name, ticklabels, start_edge, end_edge):
        # horizontal bars
        ax.set_ylim((start_edge, end_edge))
        ax.set_yticks(self.tick_pos)
        ax.set_yticklabels(ticklabels)
        if name is not None and self.use_index:
            ax.set_ylabel(name)


class HistPlot(LinePlot):
    _kind = 'hist'

    def __init__(self, data, bins=10, bottom=0, **kwargs):
        self.bins = bins        # use mpl default
        self.bottom = bottom
        # Do not call LinePlot.__init__ which may fill nan
        MPLPlot.__init__(self, data, **kwargs)

    def _args_adjust(self):
        if is_integer(self.bins):
            # create common bin edge
            values = (self.data._convert(datetime=True)._get_numeric_data())
            values = np.ravel(values)
            values = values[~isnull(values)]

            hist, self.bins = np.histogram(
                values, bins=self.bins,
                range=self.kwds.get('range', None),
                weights=self.kwds.get('weights', None))

        if is_list_like(self.bottom):
            self.bottom = np.array(self.bottom)

    @classmethod
    def _plot(cls, ax, y, style=None, bins=None, bottom=0, column_num=0,
              stacking_id=None, **kwds):
        if column_num == 0:
            cls._initialize_stacker(ax, stacking_id, len(bins) - 1)
        y = y[~isnull(y)]

        base = np.zeros(len(bins) - 1)
        bottom = bottom + \
            cls._get_stacked_values(ax, stacking_id, base, kwds['label'])
        # ignore style
        n, bins, patches = ax.hist(y, bins=bins, bottom=bottom, **kwds)
        cls._update_stacker(ax, stacking_id, n)
        return patches

    def _make_plot(self):
        colors = self._get_colors()
        stacking_id = self._get_stacking_id()

        for i, (label, y) in enumerate(self._iter_data()):
            ax = self._get_ax(i)

            kwds = self.kwds.copy()

            label = pprint_thing(label)
            kwds['label'] = label

            style, kwds = self._apply_style_colors(colors, kwds, i, label)
            if style is not None:
                kwds['style'] = style

            kwds = self._make_plot_keywords(kwds, y)
            artists = self._plot(ax, y, column_num=i,
                                 stacking_id=stacking_id, **kwds)
            self._add_legend_handle(artists[0], label, index=i)

    def _make_plot_keywords(self, kwds, y):
        """merge BoxPlot/KdePlot properties to passed kwds"""
        # y is required for KdePlot
        kwds['bottom'] = self.bottom
        kwds['bins'] = self.bins
        return kwds

    def _post_plot_logic(self, ax, data):
        if self.orientation == 'horizontal':
            ax.set_xlabel('Frequency')
        else:
            ax.set_ylabel('Frequency')

    @property
    def orientation(self):
        if self.kwds.get('orientation', None) == 'horizontal':
            return 'horizontal'
        else:
            return 'vertical'


class KdePlot(HistPlot):
    _kind = 'kde'
    orientation = 'vertical'

    def __init__(self, data, bw_method=None, ind=None, **kwargs):
        MPLPlot.__init__(self, data, **kwargs)
        self.bw_method = bw_method
        self.ind = ind

    def _args_adjust(self):
        pass

    def _get_ind(self, y):
        if self.ind is None:
            # np.nanmax() and np.nanmin() ignores the missing values
            sample_range = np.nanmax(y) - np.nanmin(y)
            ind = np.linspace(np.nanmin(y) - 0.5 * sample_range,
                              np.nanmax(y) + 0.5 * sample_range, 1000)
        else:
            ind = self.ind
        return ind

    @classmethod
    def _plot(cls, ax, y, style=None, bw_method=None, ind=None,
              column_num=None, stacking_id=None, **kwds):
        from scipy.stats import gaussian_kde
        from scipy import __version__ as spv

        y = remove_na(y)

        if LooseVersion(spv) >= '0.11.0':
            gkde = gaussian_kde(y, bw_method=bw_method)
        else:
            gkde = gaussian_kde(y)
            if bw_method is not None:
                msg = ('bw_method was added in Scipy 0.11.0.' +
                       ' Scipy version in use is %s.' % spv)
                warnings.warn(msg)

        y = gkde.evaluate(ind)
        lines = MPLPlot._plot(ax, ind, y, style=style, **kwds)
        return lines

    def _make_plot_keywords(self, kwds, y):
        kwds['bw_method'] = self.bw_method
        kwds['ind'] = self._get_ind(y)
        return kwds

    def _post_plot_logic(self, ax, data):
        ax.set_ylabel('Density')


class PiePlot(MPLPlot):
    _kind = 'pie'
    _layout_type = 'horizontal'

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

    def _validate_color_args(self):
        pass

    def _make_plot(self):
        colors = self._get_colors(
            num_colors=len(self.data), color_kwds='colors')
        self.kwds.setdefault('colors', colors)

        for i, (label, y) in enumerate(self._iter_data()):
            ax = self._get_ax(i)
            if label is not None:
                label = pprint_thing(label)
                ax.set_ylabel(label)

            kwds = self.kwds.copy()

            def blank_labeler(label, value):
                if value == 0:
                    return ''
                else:
                    return label

            idx = [pprint_thing(v) for v in self.data.index]
            labels = kwds.pop('labels', idx)
            # labels is used for each wedge's labels
            # Blank out labels for values of 0 so they don't overlap
            # with nonzero wedges
            if labels is not None:
                blabels = [blank_labeler(l, value) for
                           l, value in zip(labels, y)]
            else:
                blabels = None
            results = ax.pie(y, labels=blabels, **kwds)

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


class BoxPlot(LinePlot):
    _kind = 'box'
    _layout_type = 'horizontal'

    _valid_return_types = (None, 'axes', 'dict', 'both')
    # namedtuple to hold results
    BP = namedtuple("Boxplot", ['ax', 'lines'])

    def __init__(self, data, return_type='axes', **kwargs):
        # Do not call LinePlot.__init__ which may fill nan
        if return_type not in self._valid_return_types:
            raise ValueError(
                "return_type must be {None, 'axes', 'dict', 'both'}")

        self.return_type = return_type
        MPLPlot.__init__(self, data, **kwargs)

    def _args_adjust(self):
        if self.subplots:
            # Disable label ax sharing. Otherwise, all subplots shows last
            # column label
            if self.orientation == 'vertical':
                self.sharex = False
            else:
                self.sharey = False

    @classmethod
    def _plot(cls, ax, y, column_num=None, return_type='axes', **kwds):
        if y.ndim == 2:
            y = [remove_na(v) for v in y]
            # Boxplot fails with empty arrays, so need to add a NaN
            #   if any cols are empty
            # GH 8181
            y = [v if v.size > 0 else np.array([np.nan]) for v in y]
        else:
            y = remove_na(y)
        bp = ax.boxplot(y, **kwds)

        if return_type == 'dict':
            return bp, bp
        elif return_type == 'both':
            return cls.BP(ax=ax, lines=bp), bp
        else:
            return ax, bp

    def _validate_color_args(self):
        if 'color' in self.kwds:
            if self.colormap is not None:
                warnings.warn("'color' and 'colormap' cannot be used "
                              "simultaneously. Using 'color'")
            self.color = self.kwds.pop('color')

            if isinstance(self.color, dict):
                valid_keys = ['boxes', 'whiskers', 'medians', 'caps']
                for key, values in compat.iteritems(self.color):
                    if key not in valid_keys:
                        raise ValueError("color dict contains invalid "
                                         "key '{0}' "
                                         "The key must be either {1}"
                                         .format(key, valid_keys))
        else:
            self.color = None

        # get standard colors for default
        colors = _get_standard_colors(num_colors=3,
                                      colormap=self.colormap,
                                      color=None)
        # use 2 colors by default, for box/whisker and median
        # flier colors isn't needed here
        # because it can be specified by ``sym`` kw
        self._boxes_c = colors[0]
        self._whiskers_c = colors[0]
        self._medians_c = colors[2]
        self._caps_c = 'k'          # mpl default

    def _get_colors(self, num_colors=None, color_kwds='color'):
        pass

    def maybe_color_bp(self, bp):
        if isinstance(self.color, dict):
            boxes = self.color.get('boxes', self._boxes_c)
            whiskers = self.color.get('whiskers', self._whiskers_c)
            medians = self.color.get('medians', self._medians_c)
            caps = self.color.get('caps', self._caps_c)
        else:
            # Other types are forwarded to matplotlib
            # If None, use default colors
            boxes = self.color or self._boxes_c
            whiskers = self.color or self._whiskers_c
            medians = self.color or self._medians_c
            caps = self.color or self._caps_c

        from matplotlib.artist import setp
        setp(bp['boxes'], color=boxes, alpha=1)
        setp(bp['whiskers'], color=whiskers, alpha=1)
        setp(bp['medians'], color=medians, alpha=1)
        setp(bp['caps'], color=caps, alpha=1)

    def _make_plot(self):
        if self.subplots:
            self._return_obj = Series()

            for i, (label, y) in enumerate(self._iter_data()):
                ax = self._get_ax(i)
                kwds = self.kwds.copy()

                ret, bp = self._plot(ax, y, column_num=i,
                                     return_type=self.return_type, **kwds)
                self.maybe_color_bp(bp)
                self._return_obj[label] = ret

                label = [pprint_thing(label)]
                self._set_ticklabels(ax, label)
        else:
            y = self.data.values.T
            ax = self._get_ax(0)
            kwds = self.kwds.copy()

            ret, bp = self._plot(ax, y, column_num=0,
                                 return_type=self.return_type, **kwds)
            self.maybe_color_bp(bp)
            self._return_obj = ret

            labels = [l for l, _ in self._iter_data()]
            labels = [pprint_thing(l) for l in labels]
            if not self.use_index:
                labels = [pprint_thing(key) for key in range(len(labels))]
            self._set_ticklabels(ax, labels)

    def _set_ticklabels(self, ax, labels):
        if self.orientation == 'vertical':
            ax.set_xticklabels(labels)
        else:
            ax.set_yticklabels(labels)

    def _make_legend(self):
        pass

    def _post_plot_logic(self, ax, data):
        pass

    @property
    def orientation(self):
        if self.kwds.get('vert', True):
            return 'vertical'
        else:
            return 'horizontal'

    @property
    def result(self):
        if self.return_type is None:
            return super(BoxPlot, self).result
        else:
            return self._return_obj


# kinds supported by both dataframe and series
_common_kinds = ['line', 'bar', 'barh',
                 'kde', 'density', 'area', 'hist', 'box']
# kinds supported by dataframe
_dataframe_kinds = ['scatter', 'hexbin']
# kinds supported only by series or dataframe single column
_series_kinds = ['pie']
_all_kinds = _common_kinds + _dataframe_kinds + _series_kinds

_klasses = [LinePlot, BarPlot, BarhPlot, KdePlot, HistPlot, BoxPlot,
            ScatterPlot, HexBinPlot, AreaPlot, PiePlot]

_plot_klass = {}
for klass in _klasses:
    _plot_klass[klass._kind] = klass


def _plot(data, x=None, y=None, subplots=False,
          ax=None, kind='line', **kwds):
    kind = _get_standard_kind(kind.lower().strip())
    if kind in _all_kinds:
        klass = _plot_klass[kind]
    else:
        raise ValueError("%r is not a valid plot kind" % kind)

    from pandas import DataFrame
    if kind in _dataframe_kinds:
        if isinstance(data, DataFrame):
            plot_obj = klass(data, x=x, y=y, subplots=subplots, ax=ax,
                             kind=kind, **kwds)
        else:
            raise ValueError("plot kind %r can only be used for data frames"
                             % kind)

    elif kind in _series_kinds:
        if isinstance(data, DataFrame):
            if y is None and subplots is False:
                msg = "{0} requires either y column or 'subplots=True'"
                raise ValueError(msg.format(kind))
            elif y is not None:
                if is_integer(y) and not data.columns.holds_integer():
                    y = data.columns[y]
                # converted to series actually. copy to not modify
                data = data[y].copy()
                data.index.name = y
        plot_obj = klass(data, subplots=subplots, ax=ax, kind=kind, **kwds)
    else:
        if isinstance(data, DataFrame):
            if x is not None:
                if is_integer(x) and not data.columns.holds_integer():
                    x = data.columns[x]
                data = data.set_index(x)

            if y is not None:
                if is_integer(y) and not data.columns.holds_integer():
                    y = data.columns[y]
                label = kwds['label'] if 'label' in kwds else y
                series = data[y].copy()  # Don't modify
                series.name = label

                for kw in ['xerr', 'yerr']:
                    if (kw in kwds) and \
                        (isinstance(kwds[kw], string_types) or
                            is_integer(kwds[kw])):
                        try:
                            kwds[kw] = data[kwds[kw]]
                        except (IndexError, KeyError, TypeError):
                            pass
                data = series
        plot_obj = klass(data, subplots=subplots, ax=ax, kind=kind, **kwds)

    plot_obj.generate()
    plot_obj.draw()
    return plot_obj.result


df_kind = """- 'scatter' : scatter plot
        - 'hexbin' : hexbin plot"""
series_kind = ""

df_coord = """x : label or position, default None
    y : label or position, default None
        Allows plotting of one column versus another"""
series_coord = ""

df_unique = """stacked : boolean, default False in line and
        bar plots, and True in area plot. If True, create stacked plot.
    sort_columns : boolean, default False
        Sort column names to determine plot ordering
    secondary_y : boolean or sequence, default False
        Whether to plot on the secondary y-axis
        If a list/tuple, which columns to plot on secondary y-axis"""
series_unique = """label : label argument to provide to plot
    secondary_y : boolean or sequence of ints, default False
        If True then y-axis will be on the right"""

df_ax = """ax : matplotlib axes object, default None
    subplots : boolean, default False
        Make separate subplots for each column
    sharex : boolean, default True if ax is None else False
        In case subplots=True, share x axis and set some x axis labels to
        invisible; defaults to True if ax is None otherwise False if an ax
        is passed in; Be aware, that passing in both an ax and sharex=True
        will alter all x axis labels for all axis in a figure!
    sharey : boolean, default False
        In case subplots=True, share y axis and set some y axis labels to
        invisible
    layout : tuple (optional)
        (rows, columns) for the layout of subplots"""
series_ax = """ax : matplotlib axes object
        If not passed, uses gca()"""

df_note = """- If `kind` = 'scatter' and the argument `c` is the name of a dataframe
      column, the values of that column are used to color each point.
    - If `kind` = 'hexbin', you can control the size of the bins with the
      `gridsize` argument. By default, a histogram of the counts around each
      `(x, y)` point is computed. You can specify alternative aggregations
      by passing values to the `C` and `reduce_C_function` arguments.
      `C` specifies the value at each `(x, y)` point and `reduce_C_function`
      is a function of one argument that reduces all the values in a bin to
      a single number (e.g. `mean`, `max`, `sum`, `std`)."""
series_note = ""

_shared_doc_df_kwargs = dict(klass='DataFrame', klass_obj='df',
                             klass_kind=df_kind, klass_coord=df_coord,
                             klass_ax=df_ax, klass_unique=df_unique,
                             klass_note=df_note)
_shared_doc_series_kwargs = dict(klass='Series', klass_obj='s',
                                 klass_kind=series_kind,
                                 klass_coord=series_coord, klass_ax=series_ax,
                                 klass_unique=series_unique,
                                 klass_note=series_note)

_shared_docs['plot'] = """
    Make plots of %(klass)s using matplotlib / pylab.

    *New in version 0.17.0:* Each plot kind has a corresponding method on the
    ``%(klass)s.plot`` accessor:
    ``%(klass_obj)s.plot(kind='line')`` is equivalent to
    ``%(klass_obj)s.plot.line()``.

    Parameters
    ----------
    data : %(klass)s
    %(klass_coord)s
    kind : str
        - 'line' : line plot (default)
        - 'bar' : vertical bar plot
        - 'barh' : horizontal bar plot
        - 'hist' : histogram
        - 'box' : boxplot
        - 'kde' : Kernel Density Estimation plot
        - 'density' : same as 'kde'
        - 'area' : area plot
        - 'pie' : pie plot
        %(klass_kind)s
    %(klass_ax)s
    figsize : a tuple (width, height) in inches
    use_index : boolean, default True
        Use index as ticks for x axis
    title : string or list
        Title to use for the plot. If a string is passed, print the string at
        the top of the figure. If a list is passed and `subplots` is True,
        print each item in the list above the corresponding subplot.
    grid : boolean, default None (matlab style default)
        Axis grid lines
    legend : False/True/'reverse'
        Place legend on axis subplots
    style : list or dict
        matplotlib line style per column
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
        Rotation for ticks (xticks for vertical, yticks for horizontal plots)
    fontsize : int, default None
        Font size for xticks and yticks
    colormap : str or matplotlib colormap object, default None
        Colormap to select colors from. If string, load colormap with that name
        from matplotlib.
    colorbar : boolean, optional
        If True, plot colorbar (only relevant for 'scatter' and 'hexbin' plots)
    position : float
        Specify relative alignments for bar plot layout.
        From 0 (left/bottom-end) to 1 (right/top-end). Default is 0.5 (center)
    layout : tuple (optional)
        (rows, columns) for the layout of the plot
    table : boolean, Series or DataFrame, default False
        If True, draw a table using the data in the DataFrame and the data will
        be transposed to meet matplotlib's default layout.
        If a Series or DataFrame is passed, use passed data to draw a table.
    yerr : DataFrame, Series, array-like, dict and str
        See :ref:`Plotting with Error Bars <visualization.errorbars>` for
        detail.
    xerr : same types as yerr.
    %(klass_unique)s
    mark_right : boolean, default True
        When using a secondary_y axis, automatically mark the column
        labels with "(right)" in the legend
    kwds : keywords
        Options to pass to matplotlib plotting method

    Returns
    -------
    axes : matplotlib.AxesSubplot or np.array of them

    Notes
    -----

    - See matplotlib documentation online for more on this subject
    - If `kind` = 'bar' or 'barh', you can specify relative alignments
      for bar plot layout by `position` keyword.
      From 0 (left/bottom-end) to 1 (right/top-end). Default is 0.5 (center)
    %(klass_note)s

    """


@Appender(_shared_docs['plot'] % _shared_doc_df_kwargs)
def plot_frame(data, x=None, y=None, kind='line', ax=None,
               subplots=False, sharex=None, sharey=False, layout=None,
               figsize=None, use_index=True, title=None, grid=None,
               legend=True, style=None, logx=False, logy=False, loglog=False,
               xticks=None, yticks=None, xlim=None, ylim=None,
               rot=None, fontsize=None, colormap=None, table=False,
               yerr=None, xerr=None,
               secondary_y=False, sort_columns=False,
               **kwds):
    return _plot(data, kind=kind, x=x, y=y, ax=ax,
                 subplots=subplots, sharex=sharex, sharey=sharey,
                 layout=layout, figsize=figsize, use_index=use_index,
                 title=title, grid=grid, legend=legend,
                 style=style, logx=logx, logy=logy, loglog=loglog,
                 xticks=xticks, yticks=yticks, xlim=xlim, ylim=ylim,
                 rot=rot, fontsize=fontsize, colormap=colormap, table=table,
                 yerr=yerr, xerr=xerr,
                 secondary_y=secondary_y, sort_columns=sort_columns,
                 **kwds)


@Appender(_shared_docs['plot'] % _shared_doc_series_kwargs)
def plot_series(data, kind='line', ax=None,                    # Series unique
                figsize=None, use_index=True, title=None, grid=None,
                legend=False, style=None, logx=False, logy=False, loglog=False,
                xticks=None, yticks=None, xlim=None, ylim=None,
                rot=None, fontsize=None, colormap=None, table=False,
                yerr=None, xerr=None,
                label=None, secondary_y=False,                 # Series unique
                **kwds):

    import matplotlib.pyplot as plt
    """
    If no axes is specified, check whether there are existing figures
    If there is no existing figures, _gca() will
    create a figure with the default figsize, causing the figsize=parameter to
    be ignored.
    """
    if ax is None and len(plt.get_fignums()) > 0:
        ax = _gca()
        ax = MPLPlot._get_ax_layer(ax)
    return _plot(data, kind=kind, ax=ax,
                 figsize=figsize, use_index=use_index, title=title,
                 grid=grid, legend=legend,
                 style=style, logx=logx, logy=logy, loglog=loglog,
                 xticks=xticks, yticks=yticks, xlim=xlim, ylim=ylim,
                 rot=rot, fontsize=fontsize, colormap=colormap, table=table,
                 yerr=yerr, xerr=xerr,
                 label=label, secondary_y=secondary_y,
                 **kwds)


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
    return_type : {None, 'axes', 'dict', 'both'}, default None
        The kind of object to return. The default is ``axes``
        'axes' returns the matplotlib axes the boxplot is drawn on;
        'dict' returns a dictionary  whose values are the matplotlib
        Lines of the boxplot;
        'both' returns a namedtuple with the axes and dict.

        When grouping with ``by``, a Series mapping columns to ``return_type``
        is returned, unless ``return_type`` is None, in which case a NumPy
        array of axes is returned with the same shape as ``layout``.
        See the prose documentation for more.

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
    if return_type not in BoxPlot._valid_return_types:
        raise ValueError("return_type must be {'axes', 'dict', 'both'}")

    from pandas import Series, DataFrame
    if isinstance(data, Series):
        data = DataFrame({'x': data})
        column = 'x'

    def _get_colors():
        return _get_standard_colors(color=kwds.get('color'), num_colors=1)

    def maybe_color_bp(bp):
        if 'color' not in kwds:
            from matplotlib.artist import setp
            setp(bp['boxes'], color=colors[0], alpha=1)
            setp(bp['whiskers'], color=colors[0], alpha=1)
            setp(bp['medians'], color=colors[2], alpha=1)

    def plot_group(keys, values, ax):
        keys = [pprint_thing(x) for x in keys]
        values = [remove_na(v) for v in values]
        bp = ax.boxplot(values, **kwds)
        if fontsize is not None:
            ax.tick_params(axis='both', labelsize=fontsize)
        if kwds.get('vert', 1):
            ax.set_xticklabels(keys, rotation=rot)
        else:
            ax.set_yticklabels(keys, rotation=rot)
        maybe_color_bp(bp)

        # Return axes in multiplot case, maybe revisit later # 985
        if return_type == 'dict':
            return bp
        elif return_type == 'both':
            return BoxPlot.BP(ax=ax, lines=bp)
        else:
            return ax

    colors = _get_colors()
    if column is None:
        columns = None
    else:
        if isinstance(column, (list, tuple)):
            columns = column
        else:
            columns = [column]

    if by is not None:
        # Prefer array return type for 2-D plots to match the subplot layout
        # https://github.com/pandas-dev/pandas/pull/12216#issuecomment-241175580
        result = _grouped_plot_by_column(plot_group, data, columns=columns,
                                         by=by, grid=grid, figsize=figsize,
                                         ax=ax, layout=layout,
                                         return_type=return_type)
    else:
        if return_type is None:
            return_type = 'axes'
        if layout is not None:
            raise ValueError("The 'layout' keyword is not supported when "
                             "'by' is None")

        if ax is None:
            ax = _gca()
        data = data._get_numeric_data()
        if columns is None:
            columns = data.columns
        else:
            data = data[columns]

        result = plot_group(columns, data.values.T, ax)
        ax.grid(grid)

    return result


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


def scatter_plot(data, x, y, by=None, ax=None, figsize=None, grid=False,
                 **kwargs):
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
        ax.set_ylabel(pprint_thing(y))
        ax.set_xlabel(pprint_thing(x))

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
    sharex : boolean, default True if ax is None else False
        In case subplots=True, share x axis and set some x axis labels to
        invisible; defaults to True if ax is None otherwise False if an ax
        is passed in; Be aware, that passing in both an ax and sharex=True
        will alter all x axis labels for all subplots in a figure!
    sharey : boolean, default False
        In case subplots=True, share y axis and set some y axis labels to
        invisible
    figsize : tuple
        The size of the figure to create in inches by default
    layout : tuple, optional
        Tuple of (rows, columns) for the layout of the histograms
    bins : integer, default 10
        Number of histogram bins to be used
    kwds : other plotting keyword arguments
        To be passed to hist function
    """

    if by is not None:
        axes = grouped_hist(data, column=column, by=by, ax=ax, grid=grid,
                            figsize=figsize, sharex=sharex, sharey=sharey,
                            layout=layout, bins=bins, xlabelsize=xlabelsize,
                            xrot=xrot, ylabelsize=ylabelsize,
                            yrot=yrot, **kwds)
        return axes

    if column is not None:
        if not isinstance(column, (list, np.ndarray, Index)):
            column = [column]
        data = data[column]
    data = data._get_numeric_data()
    naxes = len(data.columns)

    fig, axes = _subplots(naxes=naxes, ax=ax, squeeze=False,
                          sharex=sharex, sharey=sharey, figsize=figsize,
                          layout=layout)
    _axes = _flatten(axes)

    for i, col in enumerate(_try_sort(data.columns)):
        ax = _axes[i]
        ax.hist(data[col].dropna().values, bins=bins, **kwds)
        ax.set_title(col)
        ax.grid(grid)

    _set_ticks_props(axes, xlabelsize=xlabelsize, xrot=xrot,
                     ylabelsize=ylabelsize, yrot=yrot)
    fig.subplots_adjust(wspace=0.3, hspace=0.3)

    return axes


def hist_series(self, by=None, ax=None, grid=True, xlabelsize=None,
                xrot=None, ylabelsize=None, yrot=None, figsize=None,
                bins=10, **kwds):
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

        _set_ticks_props(axes, xlabelsize=xlabelsize, xrot=xrot,
                         ylabelsize=ylabelsize, yrot=yrot)

    else:
        if 'figure' in kwds:
            raise ValueError("Cannot pass 'figure' when using the "
                             "'by' argument, since a new 'Figure' instance "
                             "will be created")
        axes = grouped_hist(self, by=by, ax=ax, grid=grid, figsize=figsize,
                            bins=bins, xlabelsize=xlabelsize, xrot=xrot,
                            ylabelsize=ylabelsize, yrot=yrot, **kwds)

    if hasattr(axes, 'ndim'):
        if axes.ndim == 1 and len(axes) == 1:
            return axes[0]
    return axes


def grouped_hist(data, column=None, by=None, ax=None, bins=50, figsize=None,
                 layout=None, sharex=False, sharey=False, rot=90, grid=True,
                 xlabelsize=None, xrot=None, ylabelsize=None, yrot=None,
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

    xrot = xrot or rot

    fig, axes = _grouped_plot(plot_group, data, column=column,
                              by=by, sharex=sharex, sharey=sharey, ax=ax,
                              figsize=figsize, layout=layout, rot=rot)

    _set_ticks_props(axes, xlabelsize=xlabelsize, xrot=xrot,
                     ylabelsize=ylabelsize, yrot=yrot)

    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.9,
                        hspace=0.5, wspace=0.3)
    return axes


def boxplot_frame_groupby(grouped, subplots=True, column=None, fontsize=None,
                          rot=0, grid=True, ax=None, figsize=None,
                          layout=None, **kwds):
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
    ax : Matplotlib axis object, default None
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
        fig, axes = _subplots(naxes=naxes, squeeze=False,
                              ax=ax, sharex=False, sharey=True,
                              figsize=figsize, layout=layout)
        axes = _flatten(axes)

        ret = Series()
        for (key, group), ax in zip(grouped, axes):
            d = group.boxplot(ax=ax, column=column, fontsize=fontsize,
                              rot=rot, grid=grid, **kwds)
            ax.set_title(pprint_thing(key))
            ret.loc[key] = d
        fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1,
                            right=0.9, wspace=0.2)
    else:
        from pandas.tools.concat import concat
        keys, frames = zip(*grouped)
        if grouped.axis == 0:
            df = concat(frames, keys=keys, axis=1)
        else:
            if len(frames) > 1:
                df = frames[0].join(frames[1::])
            else:
                df = frames[0]
        ret = df.boxplot(column=column, fontsize=fontsize, rot=rot,
                         grid=grid, ax=ax, figsize=figsize,
                         layout=layout, **kwds)
    return ret


def _grouped_plot(plotf, data, column=None, by=None, numeric_only=True,
                  figsize=None, sharex=True, sharey=True, layout=None,
                  rot=0, ax=None, **kwargs):
    from pandas import DataFrame

    if figsize == 'default':
        # allowed to specify mpl default with 'default'
        warnings.warn("figsize='default' is deprecated. Specify figure"
                      "size by tuple instead", FutureWarning, stacklevel=4)
        figsize = None

    grouped = data.groupby(by)
    if column is not None:
        grouped = grouped[column]

    naxes = len(grouped)
    fig, axes = _subplots(naxes=naxes, figsize=figsize,
                          sharex=sharex, sharey=sharey, ax=ax,
                          layout=layout)

    _axes = _flatten(axes)

    for i, (key, group) in enumerate(grouped):
        ax = _axes[i]
        if numeric_only and isinstance(group, DataFrame):
            group = group._get_numeric_data()
        plotf(group, ax, **kwargs)
        ax.set_title(pprint_thing(key))

    return fig, axes


def _grouped_plot_by_column(plotf, data, columns=None, by=None,
                            numeric_only=True, grid=False,
                            figsize=None, ax=None, layout=None,
                            return_type=None, **kwargs):
    grouped = data.groupby(by)
    if columns is None:
        if not isinstance(by, (list, tuple)):
            by = [by]
        columns = data._get_numeric_data().columns.difference(by)
    naxes = len(columns)
    fig, axes = _subplots(naxes=naxes, sharex=True, sharey=True,
                          figsize=figsize, ax=ax, layout=layout)

    _axes = _flatten(axes)

    result = Series()
    ax_values = []

    for i, col in enumerate(columns):
        ax = _axes[i]
        gp_col = grouped[col]
        keys, values = zip(*gp_col)
        re_plotf = plotf(keys, values, ax, **kwargs)
        ax.set_title(col)
        ax.set_xlabel(pprint_thing(by))
        ax_values.append(re_plotf)
        ax.grid(grid)

    result = Series(ax_values, index=columns)

    # Return axes in multiplot case, maybe revisit later # 985
    if return_type is None:
        result = axes

    byline = by[0] if len(by) == 1 else by
    fig.suptitle('Boxplot grouped by %s' % byline)
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.9, wspace=0.2)

    return result


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
        If `rowLabels` or `colLabels` is not specified, data index or column
        name will be used.

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
        raise ValueError('Input data must be DataFrame or Series')

    if rowLabels is None:
        rowLabels = data.index

    if colLabels is None:
        colLabels = data.columns

    cellText = data.values

    import matplotlib.table
    table = matplotlib.table.table(ax, cellText=cellText,
                                   rowLabels=rowLabels,
                                   colLabels=colLabels, **kwargs)
    return table


def _get_layout(nplots, layout=None, layout_type='box'):
    if layout is not None:
        if not isinstance(layout, (tuple, list)) or len(layout) != 2:
            raise ValueError('Layout must be a tuple of (rows, columns)')

        nrows, ncols = layout

        # Python 2 compat
        ceil_ = lambda x: int(ceil(x))
        if nrows == -1 and ncols > 0:
            layout = nrows, ncols = (ceil_(float(nplots) / ncols), ncols)
        elif ncols == -1 and nrows > 0:
            layout = nrows, ncols = (nrows, ceil_(float(nplots) / nrows))
        elif ncols <= 0 and nrows <= 0:
            msg = "At least one dimension of layout must be positive"
            raise ValueError(msg)

        if nrows * ncols < nplots:
            raise ValueError('Layout of %sx%s must be larger than '
                             'required size %s' % (nrows, ncols, nplots))

        return layout

    if layout_type == 'single':
        return (1, 1)
    elif layout_type == 'horizontal':
        return (1, nplots)
    elif layout_type == 'vertical':
        return (nplots, 1)

    layouts = {1: (1, 1), 2: (1, 2), 3: (2, 2), 4: (2, 2)}
    try:
        return layouts[nplots]
    except KeyError:
        k = 1
        while k ** 2 < nplots:
            k += 1

        if (k - 1) * k >= nplots:
            return k, (k - 1)
        else:
            return k, k

# copied from matplotlib/pyplot.py and modified for pandas.plotting


def _subplots(naxes=None, sharex=False, sharey=False, squeeze=True,
              subplot_kw=None, ax=None, layout=None, layout_type='box',
              **fig_kw):
    """Create a figure with a set of subplots already made.

    This utility wrapper makes it convenient to create common layouts of
    subplots, including the enclosing figure object, in a single call.

    Keyword arguments:

    naxes : int
      Number of required axes. Exceeded axes are set invisible. Default is
      nrows * ncols.

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

    layout : tuple
      Number of rows and columns of the subplot grid.
      If not specified, calculated from naxes and layout_type

    layout_type : {'box', 'horziontal', 'vertical'}, default 'box'
      Specify how to layout the subplot grid.

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

    if subplot_kw is None:
        subplot_kw = {}

    if ax is None:
        fig = plt.figure(**fig_kw)
    else:
        if is_list_like(ax):
            ax = _flatten(ax)
            if layout is not None:
                warnings.warn("When passing multiple axes, layout keyword is "
                              "ignored", UserWarning)
            if sharex or sharey:
                warnings.warn("When passing multiple axes, sharex and sharey "
                              "are ignored. These settings must be specified "
                              "when creating axes", UserWarning,
                              stacklevel=4)
            if len(ax) == naxes:
                fig = ax[0].get_figure()
                return fig, ax
            else:
                raise ValueError("The number of passed axes must be {0}, the "
                                 "same as the output plot".format(naxes))

        fig = ax.get_figure()
        # if ax is passed and a number of subplots is 1, return ax as it is
        if naxes == 1:
            if squeeze:
                return fig, ax
            else:
                return fig, _flatten(ax)
        else:
            warnings.warn("To output multiple subplots, the figure containing "
                          "the passed axes is being cleared", UserWarning,
                          stacklevel=4)
            fig.clear()

    nrows, ncols = _get_layout(naxes, layout=layout, layout_type=layout_type)
    nplots = nrows * ncols

    # Create empty object array to hold all axes.  It's easiest to make it 1-d
    # so we can just append subplots upon creation, and then
    axarr = np.empty(nplots, dtype=object)

    # Create first subplot separately, so we can share it if requested
    ax0 = fig.add_subplot(nrows, ncols, 1, **subplot_kw)

    if sharex:
        subplot_kw['sharex'] = ax0
    if sharey:
        subplot_kw['sharey'] = ax0
    axarr[0] = ax0

    # Note off-by-one counting because add_subplot uses the MATLAB 1-based
    # convention.
    for i in range(1, nplots):
        kwds = subplot_kw.copy()
        # Set sharex and sharey to None for blank/dummy axes, these can
        # interfere with proper axis limits on the visible axes if
        # they share axes e.g. issue #7528
        if i >= naxes:
            kwds['sharex'] = None
            kwds['sharey'] = None
        ax = fig.add_subplot(nrows, ncols, i + 1, **kwds)
        axarr[i] = ax

    if naxes != nplots:
        for ax in axarr[naxes:]:
            ax.set_visible(False)

    _handle_shared_axes(axarr, nplots, naxes, nrows, ncols, sharex, sharey)

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


def _remove_labels_from_axis(axis):
    for t in axis.get_majorticklabels():
        t.set_visible(False)

    try:
        # set_visible will not be effective if
        # minor axis has NullLocator and NullFormattor (default)
        import matplotlib.ticker as ticker
        if isinstance(axis.get_minor_locator(), ticker.NullLocator):
            axis.set_minor_locator(ticker.AutoLocator())
        if isinstance(axis.get_minor_formatter(), ticker.NullFormatter):
            axis.set_minor_formatter(ticker.FormatStrFormatter(''))
        for t in axis.get_minorticklabels():
            t.set_visible(False)
    except Exception:   # pragma no cover
        raise
    axis.get_label().set_visible(False)


def _handle_shared_axes(axarr, nplots, naxes, nrows, ncols, sharex, sharey):
    if nplots > 1:

        if nrows > 1:
            try:
                # first find out the ax layout,
                # so that we can correctly handle 'gaps"
                layout = np.zeros((nrows + 1, ncols + 1), dtype=np.bool)
                for ax in axarr:
                    layout[ax.rowNum, ax.colNum] = ax.get_visible()

                for ax in axarr:
                    # only the last row of subplots should get x labels -> all
                    # other off layout handles the case that the subplot is
                    # the last in the column, because below is no subplot/gap.
                    if not layout[ax.rowNum + 1, ax.colNum]:
                        continue
                    if sharex or len(ax.get_shared_x_axes()
                                     .get_siblings(ax)) > 1:
                        _remove_labels_from_axis(ax.xaxis)

            except IndexError:
                # if gridspec is used, ax.rowNum and ax.colNum may different
                # from layout shape. in this case, use last_row logic
                for ax in axarr:
                    if ax.is_last_row():
                        continue
                    if sharex or len(ax.get_shared_x_axes()
                                     .get_siblings(ax)) > 1:
                        _remove_labels_from_axis(ax.xaxis)

        if ncols > 1:
            for ax in axarr:
                # only the first column should get y labels -> set all other to
                # off as we only have labels in teh first column and we always
                # have a subplot there, we can skip the layout test
                if ax.is_first_col():
                    continue
                if sharey or len(ax.get_shared_y_axes().get_siblings(ax)) > 1:
                    _remove_labels_from_axis(ax.yaxis)


def _flatten(axes):
    if not is_list_like(axes):
        return np.array([axes])
    elif isinstance(axes, (np.ndarray, Index)):
        return axes.ravel()
    return np.array(axes)


def _get_all_lines(ax):
    lines = ax.get_lines()

    if hasattr(ax, 'right_ax'):
        lines += ax.right_ax.get_lines()

    if hasattr(ax, 'left_ax'):
        lines += ax.left_ax.get_lines()

    return lines


def _get_xlim(lines):
    left, right = np.inf, -np.inf
    for l in lines:
        x = l.get_xdata(orig=False)
        left = min(x[0], left)
        right = max(x[-1], right)
    return left, right


def _set_ticks_props(axes, xlabelsize=None, xrot=None,
                     ylabelsize=None, yrot=None):
    import matplotlib.pyplot as plt

    for ax in _flatten(axes):
        if xlabelsize is not None:
            plt.setp(ax.get_xticklabels(), fontsize=xlabelsize)
        if xrot is not None:
            plt.setp(ax.get_xticklabels(), rotation=xrot)
        if ylabelsize is not None:
            plt.setp(ax.get_yticklabels(), fontsize=ylabelsize)
        if yrot is not None:
            plt.setp(ax.get_yticklabels(), rotation=yrot)
    return axes


class BasePlotMethods(PandasObject):

    def __init__(self, data):
        self._data = data

    def __call__(self, *args, **kwargs):
        raise NotImplementedError


class SeriesPlotMethods(BasePlotMethods):
    """Series plotting accessor and method

    Examples
    --------
    >>> s.plot.line()
    >>> s.plot.bar()
    >>> s.plot.hist()

    Plotting methods can also be accessed by calling the accessor as a method
    with the ``kind`` argument:
    ``s.plot(kind='line')`` is equivalent to ``s.plot.line()``
    """

    def __call__(self, kind='line', ax=None,
                 figsize=None, use_index=True, title=None, grid=None,
                 legend=False, style=None, logx=False, logy=False,
                 loglog=False, xticks=None, yticks=None,
                 xlim=None, ylim=None,
                 rot=None, fontsize=None, colormap=None, table=False,
                 yerr=None, xerr=None,
                 label=None, secondary_y=False, **kwds):
        return plot_series(self._data, kind=kind, ax=ax, figsize=figsize,
                           use_index=use_index, title=title, grid=grid,
                           legend=legend, style=style, logx=logx, logy=logy,
                           loglog=loglog, xticks=xticks, yticks=yticks,
                           xlim=xlim, ylim=ylim, rot=rot, fontsize=fontsize,
                           colormap=colormap, table=table, yerr=yerr,
                           xerr=xerr, label=label, secondary_y=secondary_y,
                           **kwds)
    __call__.__doc__ = plot_series.__doc__

    def line(self, **kwds):
        """
        Line plot

        .. versionadded:: 0.17.0

        Parameters
        ----------
        **kwds : optional
            Keyword arguments to pass on to :py:meth:`pandas.Series.plot`.

        Returns
        -------
        axes : matplotlib.AxesSubplot or np.array of them
        """
        return self(kind='line', **kwds)

    def bar(self, **kwds):
        """
        Vertical bar plot

        .. versionadded:: 0.17.0

        Parameters
        ----------
        **kwds : optional
            Keyword arguments to pass on to :py:meth:`pandas.Series.plot`.

        Returns
        -------
        axes : matplotlib.AxesSubplot or np.array of them
        """
        return self(kind='bar', **kwds)

    def barh(self, **kwds):
        """
        Horizontal bar plot

        .. versionadded:: 0.17.0

        Parameters
        ----------
        **kwds : optional
            Keyword arguments to pass on to :py:meth:`pandas.Series.plot`.

        Returns
        -------
        axes : matplotlib.AxesSubplot or np.array of them
        """
        return self(kind='barh', **kwds)

    def box(self, **kwds):
        """
        Boxplot

        .. versionadded:: 0.17.0

        Parameters
        ----------
        **kwds : optional
            Keyword arguments to pass on to :py:meth:`pandas.Series.plot`.

        Returns
        -------
        axes : matplotlib.AxesSubplot or np.array of them
        """
        return self(kind='box', **kwds)

    def hist(self, bins=10, **kwds):
        """
        Histogram

        .. versionadded:: 0.17.0

        Parameters
        ----------
        bins: integer, default 10
            Number of histogram bins to be used
        **kwds : optional
            Keyword arguments to pass on to :py:meth:`pandas.Series.plot`.

        Returns
        -------
        axes : matplotlib.AxesSubplot or np.array of them
        """
        return self(kind='hist', bins=bins, **kwds)

    def kde(self, **kwds):
        """
        Kernel Density Estimate plot

        .. versionadded:: 0.17.0

        Parameters
        ----------
        **kwds : optional
            Keyword arguments to pass on to :py:meth:`pandas.Series.plot`.

        Returns
        -------
        axes : matplotlib.AxesSubplot or np.array of them
        """
        return self(kind='kde', **kwds)

    density = kde

    def area(self, **kwds):
        """
        Area plot

        .. versionadded:: 0.17.0

        Parameters
        ----------
        **kwds : optional
            Keyword arguments to pass on to :py:meth:`pandas.Series.plot`.

        Returns
        -------
        axes : matplotlib.AxesSubplot or np.array of them
        """
        return self(kind='area', **kwds)

    def pie(self, **kwds):
        """
        Pie chart

        .. versionadded:: 0.17.0

        Parameters
        ----------
        **kwds : optional
            Keyword arguments to pass on to :py:meth:`pandas.Series.plot`.

        Returns
        -------
        axes : matplotlib.AxesSubplot or np.array of them
        """
        return self(kind='pie', **kwds)


class FramePlotMethods(BasePlotMethods):
    """DataFrame plotting accessor and method

    Examples
    --------
    >>> df.plot.line()
    >>> df.plot.scatter('x', 'y')
    >>> df.plot.hexbin()

    These plotting methods can also be accessed by calling the accessor as a
    method with the ``kind`` argument:
    ``df.plot(kind='line')`` is equivalent to ``df.plot.line()``
    """

    def __call__(self, x=None, y=None, kind='line', ax=None,
                 subplots=False, sharex=None, sharey=False, layout=None,
                 figsize=None, use_index=True, title=None, grid=None,
                 legend=True, style=None, logx=False, logy=False, loglog=False,
                 xticks=None, yticks=None, xlim=None, ylim=None,
                 rot=None, fontsize=None, colormap=None, table=False,
                 yerr=None, xerr=None,
                 secondary_y=False, sort_columns=False, **kwds):
        return plot_frame(self._data, kind=kind, x=x, y=y, ax=ax,
                          subplots=subplots, sharex=sharex, sharey=sharey,
                          layout=layout, figsize=figsize, use_index=use_index,
                          title=title, grid=grid, legend=legend, style=style,
                          logx=logx, logy=logy, loglog=loglog, xticks=xticks,
                          yticks=yticks, xlim=xlim, ylim=ylim, rot=rot,
                          fontsize=fontsize, colormap=colormap, table=table,
                          yerr=yerr, xerr=xerr, secondary_y=secondary_y,
                          sort_columns=sort_columns, **kwds)
    __call__.__doc__ = plot_frame.__doc__

    def line(self, x=None, y=None, **kwds):
        """
        Line plot

        .. versionadded:: 0.17.0

        Parameters
        ----------
        x, y : label or position, optional
            Coordinates for each point.
        **kwds : optional
            Keyword arguments to pass on to :py:meth:`pandas.DataFrame.plot`.

        Returns
        -------
        axes : matplotlib.AxesSubplot or np.array of them
        """
        return self(kind='line', x=x, y=y, **kwds)

    def bar(self, x=None, y=None, **kwds):
        """
        Vertical bar plot

        .. versionadded:: 0.17.0

        Parameters
        ----------
        x, y : label or position, optional
            Coordinates for each point.
        **kwds : optional
            Keyword arguments to pass on to :py:meth:`pandas.DataFrame.plot`.

        Returns
        -------
        axes : matplotlib.AxesSubplot or np.array of them
        """
        return self(kind='bar', x=x, y=y, **kwds)

    def barh(self, x=None, y=None, **kwds):
        """
        Horizontal bar plot

        .. versionadded:: 0.17.0

        Parameters
        ----------
        x, y : label or position, optional
            Coordinates for each point.
        **kwds : optional
            Keyword arguments to pass on to :py:meth:`pandas.DataFrame.plot`.

        Returns
        -------
        axes : matplotlib.AxesSubplot or np.array of them
        """
        return self(kind='barh', x=x, y=y, **kwds)

    def box(self, by=None, **kwds):
        """
        Boxplot

        .. versionadded:: 0.17.0

        Parameters
        ----------
        by : string or sequence
            Column in the DataFrame to group by.
        \*\*kwds : optional
            Keyword arguments to pass on to :py:meth:`pandas.DataFrame.plot`.

        Returns
        -------
        axes : matplotlib.AxesSubplot or np.array of them
        """
        return self(kind='box', by=by, **kwds)

    def hist(self, by=None, bins=10, **kwds):
        """
        Histogram

        .. versionadded:: 0.17.0

        Parameters
        ----------
        by : string or sequence
            Column in the DataFrame to group by.
        bins: integer, default 10
            Number of histogram bins to be used
        **kwds : optional
            Keyword arguments to pass on to :py:meth:`pandas.DataFrame.plot`.

        Returns
        -------
        axes : matplotlib.AxesSubplot or np.array of them
        """
        return self(kind='hist', by=by, bins=bins, **kwds)

    def kde(self, **kwds):
        """
        Kernel Density Estimate plot

        .. versionadded:: 0.17.0

        Parameters
        ----------
        **kwds : optional
            Keyword arguments to pass on to :py:meth:`pandas.DataFrame.plot`.

        Returns
        -------
        axes : matplotlib.AxesSubplot or np.array of them
        """
        return self(kind='kde', **kwds)

    density = kde

    def area(self, x=None, y=None, **kwds):
        """
        Area plot

        .. versionadded:: 0.17.0

        Parameters
        ----------
        x, y : label or position, optional
            Coordinates for each point.
        **kwds : optional
            Keyword arguments to pass on to :py:meth:`pandas.DataFrame.plot`.

        Returns
        -------
        axes : matplotlib.AxesSubplot or np.array of them
        """
        return self(kind='area', x=x, y=y, **kwds)

    def pie(self, y=None, **kwds):
        """
        Pie chart

        .. versionadded:: 0.17.0

        Parameters
        ----------
        y : label or position, optional
            Column to plot.
        **kwds : optional
            Keyword arguments to pass on to :py:meth:`pandas.DataFrame.plot`.

        Returns
        -------
        axes : matplotlib.AxesSubplot or np.array of them
        """
        return self(kind='pie', y=y, **kwds)

    def scatter(self, x, y, s=None, c=None, **kwds):
        """
        Scatter plot

        .. versionadded:: 0.17.0

        Parameters
        ----------
        x, y : label or position, optional
            Coordinates for each point.
        s : scalar or array_like, optional
            Size of each point.
        c : label or position, optional
            Color of each point.
        **kwds : optional
            Keyword arguments to pass on to :py:meth:`pandas.DataFrame.plot`.

        Returns
        -------
        axes : matplotlib.AxesSubplot or np.array of them
        """
        return self(kind='scatter', x=x, y=y, c=c, s=s, **kwds)

    def hexbin(self, x, y, C=None, reduce_C_function=None, gridsize=None,
               **kwds):
        """
        Hexbin plot

        .. versionadded:: 0.17.0

        Parameters
        ----------
        x, y : label or position, optional
            Coordinates for each point.
        C : label or position, optional
            The value at each `(x, y)` point.
        reduce_C_function : callable, optional
            Function of one argument that reduces all the values in a bin to
            a single number (e.g. `mean`, `max`, `sum`, `std`).
        gridsize : int, optional
            Number of bins.
        **kwds : optional
            Keyword arguments to pass on to :py:meth:`pandas.DataFrame.plot`.

        Returns
        -------
        axes : matplotlib.AxesSubplot or np.array of them
        """
        if reduce_C_function is not None:
            kwds['reduce_C_function'] = reduce_C_function
        if gridsize is not None:
            kwds['gridsize'] = gridsize
        return self(kind='hexbin', x=x, y=y, C=C, **kwds)
