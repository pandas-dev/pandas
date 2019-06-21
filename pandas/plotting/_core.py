import importlib
from typing import List, Type  # noqa

from pandas.util._decorators import Appender

from pandas.core.dtypes.common import is_integer, is_list_like
from pandas.core.dtypes.generic import ABCDataFrame, ABCSeries

import pandas
from pandas.core.base import PandasObject
from pandas.core.generic import _shared_doc_kwargs, _shared_docs

# Trigger matplotlib import, which implicitly registers our
# converts. Implicit registration is deprecated, and when enforced
# we can lazily import matplotlib.
try:
    import pandas.plotting._matplotlib  # noqa
except ImportError:
    pass

df_kind = """- 'scatter' : scatter plot
        - 'hexbin' : hexbin plot"""
series_kind = ""

df_coord = """x : label or position, default None
    y : label, position or list of label, positions, default None
        Allows plotting of one column versus another"""
series_coord = ""

df_unique = """stacked : bool, default False in line and
        bar plots, and True in area plot. If True, create stacked plot.
    sort_columns : bool, default False
        Sort column names to determine plot ordering
    secondary_y : bool or sequence, default False
        Whether to plot on the secondary y-axis
        If a list/tuple, which columns to plot on secondary y-axis"""
series_unique = """label : label argument to provide to plot
    secondary_y : bool or sequence of ints, default False
        If True then y-axis will be on the right"""

df_ax = """ax : matplotlib axes object, default None
    subplots : bool, default False
        Make separate subplots for each column
    sharex : bool, default True if ax is None else False
        In case subplots=True, share x axis and set some x axis labels to
        invisible; defaults to True if ax is None otherwise False if an ax
        is passed in; Be aware, that passing in both an ax and sharex=True
        will alter all x axis labels for all axis in a figure!
    sharey : bool, default False
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
    use_index : bool, default True
        Use index as ticks for x axis
    title : string or list
        Title to use for the plot. If a string is passed, print the string at
        the top of the figure. If a list is passed and `subplots` is True,
        print each item in the list above the corresponding subplot.
    grid : bool, default None (matlab style default)
        Axis grid lines
    legend : False/True/'reverse'
        Place legend on axis subplots
    style : list or dict
        matplotlib line style per column
    logx : bool or 'sym', default False
        Use log scaling or symlog scaling on x axis
        .. versionchanged:: 0.25.0

    logy : bool or 'sym' default False
        Use log scaling or symlog scaling on y axis
        .. versionchanged:: 0.25.0

    loglog : bool or 'sym', default False
        Use log scaling or symlog scaling on both x and y axes
        .. versionchanged:: 0.25.0

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
    colorbar : bool, optional
        If True, plot colorbar (only relevant for 'scatter' and 'hexbin' plots)
    position : float
        Specify relative alignments for bar plot layout.
        From 0 (left/bottom-end) to 1 (right/top-end). Default is 0.5 (center)
    table : bool, Series or DataFrame, default False
        If True, draw a table using the data in the DataFrame and the data will
        be transposed to meet matplotlib's default layout.
        If a Series or DataFrame is passed, use passed data to draw a table.
    yerr : DataFrame, Series, array-like, dict and str
        See :ref:`Plotting with Error Bars <visualization.errorbars>` for
        detail.
    xerr : same types as yerr.
    %(klass_unique)s
    mark_right : bool, default True
        When using a secondary_y axis, automatically mark the column
        labels with "(right)" in the legend
    `**kwds` : keywords
        Options to pass to matplotlib plotting method

    Returns
    -------
    :class:`matplotlib.axes.Axes` or numpy.ndarray of them

    Notes
    -----

    - See matplotlib documentation online for more on this subject
    - If `kind` = 'bar' or 'barh', you can specify relative alignments
      for bar plot layout by `position` keyword.
      From 0 (left/bottom-end) to 1 (right/top-end). Default is 0.5 (center)
    %(klass_note)s
    """

_shared_docs['boxplot'] = """
    Make a box plot from DataFrame columns.

    Make a box-and-whisker plot from DataFrame columns, optionally grouped
    by some other columns. A box plot is a method for graphically depicting
    groups of numerical data through their quartiles.
    The box extends from the Q1 to Q3 quartile values of the data,
    with a line at the median (Q2). The whiskers extend from the edges
    of box to show the range of the data. The position of the whiskers
    is set by default to `1.5 * IQR (IQR = Q3 - Q1)` from the edges of the box.
    Outlier points are those past the end of the whiskers.

    For further details see
    Wikipedia's entry for `boxplot <https://en.wikipedia.org/wiki/Box_plot>`_.

    Parameters
    ----------
    column : str or list of str, optional
        Column name or list of names, or vector.
        Can be any valid input to :meth:`pandas.DataFrame.groupby`.
    by : str or array-like, optional
        Column in the DataFrame to :meth:`pandas.DataFrame.groupby`.
        One box-plot will be done per value of columns in `by`.
    ax : object of class matplotlib.axes.Axes, optional
        The matplotlib axes to be used by boxplot.
    fontsize : float or str
        Tick label font size in points or as a string (e.g., `large`).
    rot : int or float, default 0
        The rotation angle of labels (in degrees)
        with respect to the screen coordinate system.
    grid : bool, default True
        Setting this to True will show the grid.
    figsize : A tuple (width, height) in inches
        The size of the figure to create in matplotlib.
    layout : tuple (rows, columns), optional
        For example, (3, 5) will display the subplots
        using 3 columns and 5 rows, starting from the top-left.
    return_type : {'axes', 'dict', 'both'} or None, default 'axes'
        The kind of object to return. The default is ``axes``.

        * 'axes' returns the matplotlib axes the boxplot is drawn on.
        * 'dict' returns a dictionary whose values are the matplotlib
          Lines of the boxplot.
        * 'both' returns a namedtuple with the axes and dict.
        * when grouping with ``by``, a Series mapping columns to
          ``return_type`` is returned.

          If ``return_type`` is `None`, a NumPy array
          of axes with the same shape as ``layout`` is returned.
    **kwds
        All other plotting keyword arguments to be passed to
        :func:`matplotlib.pyplot.boxplot`.

    Returns
    -------
    result
        See Notes.

    See Also
    --------
    Series.plot.hist: Make a histogram.
    matplotlib.pyplot.boxplot : Matplotlib equivalent plot.

    Notes
    -----
    The return type depends on the `return_type` parameter:

    * 'axes' : object of class matplotlib.axes.Axes
    * 'dict' : dict of matplotlib.lines.Line2D objects
    * 'both' : a namedtuple with structure (ax, lines)

    For data grouped with ``by``, return a Series of the above or a numpy
    array:

    * :class:`~pandas.Series`
    * :class:`~numpy.array` (for ``return_type = None``)

    Use ``return_type='dict'`` when you want to tweak the appearance
    of the lines after plotting. In this case a dict containing the Lines
    making up the boxes, caps, fliers, medians, and whiskers is returned.

    Examples
    --------

    Boxplots can be created for every column in the dataframe
    by ``df.boxplot()`` or indicating the columns to be used:

    .. plot::
        :context: close-figs

        >>> np.random.seed(1234)
        >>> df = pd.DataFrame(np.random.randn(10,4),
        ...                   columns=['Col1', 'Col2', 'Col3', 'Col4'])
        >>> boxplot = df.boxplot(column=['Col1', 'Col2', 'Col3'])

    Boxplots of variables distributions grouped by the values of a third
    variable can be created using the option ``by``. For instance:

    .. plot::
        :context: close-figs

        >>> df = pd.DataFrame(np.random.randn(10, 2),
        ...                   columns=['Col1', 'Col2'])
        >>> df['X'] = pd.Series(['A', 'A', 'A', 'A', 'A',
        ...                      'B', 'B', 'B', 'B', 'B'])
        >>> boxplot = df.boxplot(by='X')

    A list of strings (i.e. ``['X', 'Y']``) can be passed to boxplot
    in order to group the data by combination of the variables in the x-axis:

    .. plot::
        :context: close-figs

        >>> df = pd.DataFrame(np.random.randn(10,3),
        ...                   columns=['Col1', 'Col2', 'Col3'])
        >>> df['X'] = pd.Series(['A', 'A', 'A', 'A', 'A',
        ...                      'B', 'B', 'B', 'B', 'B'])
        >>> df['Y'] = pd.Series(['A', 'B', 'A', 'B', 'A',
        ...                      'B', 'A', 'B', 'A', 'B'])
        >>> boxplot = df.boxplot(column=['Col1', 'Col2'], by=['X', 'Y'])

    The layout of boxplot can be adjusted giving a tuple to ``layout``:

    .. plot::
        :context: close-figs

        >>> boxplot = df.boxplot(column=['Col1', 'Col2'], by='X',
        ...                      layout=(2, 1))

    Additional formatting can be done to the boxplot, like suppressing the grid
    (``grid=False``), rotating the labels in the x-axis (i.e. ``rot=45``)
    or changing the fontsize (i.e. ``fontsize=15``):

    .. plot::
        :context: close-figs

        >>> boxplot = df.boxplot(grid=False, rot=45, fontsize=15)

    The parameter ``return_type`` can be used to select the type of element
    returned by `boxplot`.  When ``return_type='axes'`` is selected,
    the matplotlib axes on which the boxplot is drawn are returned:

        >>> boxplot = df.boxplot(column=['Col1','Col2'], return_type='axes')
        >>> type(boxplot)
        <class 'matplotlib.axes._subplots.AxesSubplot'>

    When grouping with ``by``, a Series mapping columns to ``return_type``
    is returned:

        >>> boxplot = df.boxplot(column=['Col1', 'Col2'], by='X',
        ...                      return_type='axes')
        >>> type(boxplot)
        <class 'pandas.core.series.Series'>

    If ``return_type`` is `None`, a NumPy array of axes with the same shape
    as ``layout`` is returned:

        >>> boxplot =  df.boxplot(column=['Col1', 'Col2'], by='X',
        ...                       return_type=None)
        >>> type(boxplot)
        <class 'numpy.ndarray'>
    """

_shared_docs['kde'] = """
        Generate Kernel Density Estimate plot using Gaussian kernels.

        In statistics, `kernel density estimation`_ (KDE) is a non-parametric
        way to estimate the probability density function (PDF) of a random
        variable. This function uses Gaussian kernels and includes automatic
        bandwidth determination.

        .. _kernel density estimation:
            https://en.wikipedia.org/wiki/Kernel_density_estimation

        Parameters
        ----------
        bw_method : str, scalar or callable, optional
            The method used to calculate the estimator bandwidth. This can be
            'scott', 'silverman', a scalar constant or a callable.
            If None (default), 'scott' is used.
            See :class:`scipy.stats.gaussian_kde` for more information.
        ind : NumPy array or integer, optional
            Evaluation points for the estimated PDF. If None (default),
            1000 equally spaced points are used. If `ind` is a NumPy array, the
            KDE is evaluated at the points passed. If `ind` is an integer,
            `ind` number of equally spaced points are used.
        **kwds : optional
            Additional keyword arguments are documented in
            :meth:`pandas.%(this-datatype)s.plot`.

        Returns
        -------
        matplotlib.axes.Axes or numpy.ndarray of them

        See Also
        --------
        scipy.stats.gaussian_kde : Representation of a kernel-density
            estimate using Gaussian kernels. This is the function used
            internally to estimate the PDF.
        %(sibling-datatype)s.plot.kde : Generate a KDE plot for a
            %(sibling-datatype)s.

        Examples
        --------
        %(examples)s
        """


def hist_series(self, by=None, ax=None, grid=True, xlabelsize=None,
                xrot=None, ylabelsize=None, yrot=None, figsize=None,
                bins=10, **kwds):
    """
    Draw histogram of the input series using matplotlib.

    Parameters
    ----------
    by : object, optional
        If passed, then used to form histograms for separate groups
    ax : matplotlib axis object
        If not passed, uses gca()
    grid : bool, default True
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
    bins : integer or sequence, default 10
        Number of histogram bins to be used. If an integer is given, bins + 1
        bin edges are calculated and returned. If bins is a sequence, gives
        bin edges, including left edge of first bin and right edge of last
        bin. In this case, bins is returned unmodified.
    `**kwds` : keywords
        To be passed to the actual plotting function

    Returns
    -------
    matplotlib.AxesSubplot
        A histogram plot.

    See Also
    --------
    matplotlib.axes.Axes.hist : Plot a histogram using matplotlib.
    """
    plot_backend = _get_plot_backend()
    return plot_backend.hist_series(self, by=by, ax=ax, grid=grid,
                                    xlabelsize=xlabelsize, xrot=xrot,
                                    ylabelsize=ylabelsize, yrot=yrot,
                                    figsize=figsize, bins=bins, **kwds)


def hist_frame(data, column=None, by=None, grid=True, xlabelsize=None,
               xrot=None, ylabelsize=None, yrot=None, ax=None, sharex=False,
               sharey=False, figsize=None, layout=None, bins=10, **kwds):
    """
    Make a histogram of the DataFrame's.

    A `histogram`_ is a representation of the distribution of data.
    This function calls :meth:`matplotlib.pyplot.hist`, on each series in
    the DataFrame, resulting in one histogram per column.

    .. _histogram: https://en.wikipedia.org/wiki/Histogram

    Parameters
    ----------
    data : DataFrame
        The pandas object holding the data.
    column : string or sequence
        If passed, will be used to limit data to a subset of columns.
    by : object, optional
        If passed, then used to form histograms for separate groups.
    grid : bool, default True
        Whether to show axis grid lines.
    xlabelsize : int, default None
        If specified changes the x-axis label size.
    xrot : float, default None
        Rotation of x axis labels. For example, a value of 90 displays the
        x labels rotated 90 degrees clockwise.
    ylabelsize : int, default None
        If specified changes the y-axis label size.
    yrot : float, default None
        Rotation of y axis labels. For example, a value of 90 displays the
        y labels rotated 90 degrees clockwise.
    ax : Matplotlib axes object, default None
        The axes to plot the histogram on.
    sharex : bool, default True if ax is None else False
        In case subplots=True, share x axis and set some x axis labels to
        invisible; defaults to True if ax is None otherwise False if an ax
        is passed in.
        Note that passing in both an ax and sharex=True will alter all x axis
        labels for all subplots in a figure.
    sharey : bool, default False
        In case subplots=True, share y axis and set some y axis labels to
        invisible.
    figsize : tuple
        The size in inches of the figure to create. Uses the value in
        `matplotlib.rcParams` by default.
    layout : tuple, optional
        Tuple of (rows, columns) for the layout of the histograms.
    bins : integer or sequence, default 10
        Number of histogram bins to be used. If an integer is given, bins + 1
        bin edges are calculated and returned. If bins is a sequence, gives
        bin edges, including left edge of first bin and right edge of last
        bin. In this case, bins is returned unmodified.
    **kwds
        All other plotting keyword arguments to be passed to
        :meth:`matplotlib.pyplot.hist`.

    Returns
    -------
    matplotlib.AxesSubplot or numpy.ndarray of them

    See Also
    --------
    matplotlib.pyplot.hist : Plot a histogram using matplotlib.

    Examples
    --------

    .. plot::
        :context: close-figs

        This example draws a histogram based on the length and width of
        some animals, displayed in three bins

        >>> df = pd.DataFrame({
        ...     'length': [1.5, 0.5, 1.2, 0.9, 3],
        ...     'width': [0.7, 0.2, 0.15, 0.2, 1.1]
        ...     }, index= ['pig', 'rabbit', 'duck', 'chicken', 'horse'])
        >>> hist = df.hist(bins=3)
    """
    plot_backend = _get_plot_backend()
    return plot_backend.hist_frame(data, column=column, by=by, grid=grid,
                                   xlabelsize=xlabelsize, xrot=xrot,
                                   ylabelsize=ylabelsize, yrot=yrot,
                                   ax=ax, sharex=sharex, sharey=sharey,
                                   figsize=figsize, layout=layout, bins=bins,
                                   **kwds)


@Appender(_shared_docs['boxplot'] % _shared_doc_kwargs)
def boxplot(data, column=None, by=None, ax=None, fontsize=None,
            rot=0, grid=True, figsize=None, layout=None, return_type=None,
            **kwds):
    plot_backend = _get_plot_backend()
    return plot_backend.boxplot(data, column=column, by=by, ax=ax,
                                fontsize=fontsize, rot=rot, grid=grid,
                                figsize=figsize, layout=layout,
                                return_type=return_type, **kwds)


@Appender(_shared_docs['boxplot'] % _shared_doc_kwargs)
def boxplot_frame(self, column=None, by=None, ax=None, fontsize=None, rot=0,
                  grid=True, figsize=None, layout=None,
                  return_type=None, **kwds):
    plot_backend = _get_plot_backend()
    return plot_backend.boxplot_frame(self, column=column, by=by, ax=ax,
                                      fontsize=fontsize, rot=rot, grid=grid,
                                      figsize=figsize, layout=layout,
                                      return_type=return_type, **kwds)


def boxplot_frame_groupby(grouped, subplots=True, column=None, fontsize=None,
                          rot=0, grid=True, ax=None, figsize=None,
                          layout=None, sharex=False, sharey=True, **kwds):
    """
    Make box plots from DataFrameGroupBy data.

    Parameters
    ----------
    grouped : Grouped DataFrame
    subplots : bool
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
    sharex : bool, default False
        Whether x-axes will be shared among subplots

        .. versionadded:: 0.23.1
    sharey : bool, default True
        Whether y-axes will be shared among subplots

        .. versionadded:: 0.23.1
    `**kwds` : Keyword Arguments
        All other plotting keyword arguments to be passed to
        matplotlib's boxplot function

    Returns
    -------
    dict of key/value = group key/DataFrame.boxplot return value
    or DataFrame.boxplot return value in case subplots=figures=False

    Examples
    --------
    >>> import itertools
    >>> tuples = [t for t in itertools.product(range(1000), range(4))]
    >>> index = pd.MultiIndex.from_tuples(tuples, names=['lvl0', 'lvl1'])
    >>> data = np.random.randn(len(index),4)
    >>> df = pd.DataFrame(data, columns=list('ABCD'), index=index)
    >>>
    >>> grouped = df.groupby(level='lvl1')
    >>> boxplot_frame_groupby(grouped)
    >>>
    >>> grouped = df.unstack(level='lvl1').groupby(level=0, axis=1)
    >>> boxplot_frame_groupby(grouped, subplots=False)
    """
    plot_backend = _get_plot_backend()
    return plot_backend.boxplot_frame_groupby(
        grouped, subplots=subplots, column=column, fontsize=fontsize, rot=rot,
        grid=grid, ax=ax, figsize=figsize, layout=layout, sharex=sharex,
        sharey=sharey, **kwds)


# kinds supported by both dataframe and series
_common_kinds = ['line', 'bar', 'barh',
                 'kde', 'density', 'area', 'hist', 'box']
# kinds supported by dataframe
_dataframe_kinds = ['scatter', 'hexbin']
# kinds supported only by series or dataframe single column
_series_kinds = ['pie']
_all_kinds = _common_kinds + _dataframe_kinds + _series_kinds


def _get_standard_kind(kind):
    return {'density': 'kde'}.get(kind, kind)


def _get_plot_backend():
    """
    Return the plotting backend to use (e.g. `pandas.plotting._matplotlib`).

    The plotting system of pandas has been using matplotlib, but the idea here
    is that it can also work with other third-party backends. In the future,
    this function will return the backend from a pandas option, and all the
    rest of the code in this file will use the backend specified there for the
    plotting.

    The backend is imported lazily, as matplotlib is a soft dependency, and
    pandas can be used without it being installed.
    """
    backend_str = pandas.get_option('plotting.backend')
    if backend_str == 'matplotlib':
        backend_str = 'pandas.plotting._matplotlib'
    return importlib.import_module(backend_str)


def _plot_classes():
    plot_backend = _get_plot_backend()
    # TODO restore type annotations if we create a base class for plot classes
    # (a parent of MPLPlot, and classes of other backends)
    classes = [plot_backend.LinePlot, plot_backend.BarPlot,
               plot_backend.BarhPlot, plot_backend.AreaPlot,
               plot_backend.HistPlot, plot_backend.BoxPlot,
               plot_backend.ScatterPlot, plot_backend.HexBinPlot,
               plot_backend.KdePlot, plot_backend.PiePlot]
    return {class_._kind: class_ for class_ in classes}


def _plot(data, x=None, y=None, subplots=False,
          ax=None, kind='line', **kwds):
    kind = _get_standard_kind(kind.lower().strip())
    if kind in _all_kinds:
        klass = _plot_classes()[kind]
    else:
        raise ValueError("%r is not a valid plot kind" % kind)

    if kind in _dataframe_kinds:
        if isinstance(data, ABCDataFrame):
            plot_obj = klass(data, x=x, y=y, subplots=subplots, ax=ax,
                             kind=kind, **kwds)
        else:
            raise ValueError("plot kind %r can only be used for data frames"
                             % kind)

    elif kind in _series_kinds:
        if isinstance(data, ABCDataFrame):
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
        if isinstance(data, ABCDataFrame):
            data_cols = data.columns
            if x is not None:
                if is_integer(x) and not data.columns.holds_integer():
                    x = data_cols[x]
                elif not isinstance(data[x], ABCSeries):
                    raise ValueError("x must be a label or position")
                data = data.set_index(x)

            if y is not None:
                # check if we have y as int or list of ints
                int_ylist = is_list_like(y) and all(is_integer(c) for c in y)
                int_y_arg = is_integer(y) or int_ylist
                if int_y_arg and not data.columns.holds_integer():
                    y = data_cols[y]

                label_kw = kwds['label'] if 'label' in kwds else False
                for kw in ['xerr', 'yerr']:
                    if (kw in kwds) and \
                        (isinstance(kwds[kw], str) or
                            is_integer(kwds[kw])):
                        try:
                            kwds[kw] = data[kwds[kw]]
                        except (IndexError, KeyError, TypeError):
                            pass

                # don't overwrite
                data = data[y].copy()

                if isinstance(data, ABCSeries):
                    label_name = label_kw or y
                    data.name = label_name
                else:
                    match = is_list_like(label_kw) and len(label_kw) == len(y)
                    if label_kw and not match:
                        raise ValueError(
                            "label should be list-like and same length as y"
                        )
                    label_name = label_kw or data.columns
                    data.columns = label_name
        plot_obj = klass(data, subplots=subplots, ax=ax, kind=kind, **kwds)

    plot_obj.generate()
    plot_obj.draw()
    return plot_obj.result


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

    # FIXME move this into _matplotlib
    import matplotlib.pyplot as plt
    if ax is None and len(plt.get_fignums()) > 0:
        with plt.rc_context():
            ax = plt.gca()
        ax = getattr(ax, 'left_ax', ax)

    return _plot(data, kind=kind, ax=ax,
                 figsize=figsize, use_index=use_index, title=title,
                 grid=grid, legend=legend,
                 style=style, logx=logx, logy=logy, loglog=loglog,
                 xticks=xticks, yticks=yticks, xlim=xlim, ylim=ylim,
                 rot=rot, fontsize=fontsize, colormap=colormap, table=table,
                 yerr=yerr, xerr=xerr,
                 label=label, secondary_y=secondary_y,
                 **kwds)


class BasePlotMethods(PandasObject):

    def __init__(self, data):
        self._parent = data  # can be Series or DataFrame

    def __call__(self, *args, **kwargs):
        raise NotImplementedError


class SeriesPlotMethods(BasePlotMethods):
    """
    Series plotting accessor and method.

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
        return plot_series(self._parent, kind=kind, ax=ax, figsize=figsize,
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
        Line plot.

        Parameters
        ----------
        `**kwds` : optional
            Additional keyword arguments are documented in
            :meth:`pandas.Series.plot`.

        Returns
        -------
        :class:`matplotlib.axes.Axes` or numpy.ndarray of them

        Examples
        --------

        .. plot::
            :context: close-figs

            >>> s = pd.Series([1, 3, 2])
            >>> s.plot.line()
        """
        return self(kind='line', **kwds)

    def bar(self, **kwds):
        """
        Vertical bar plot.

        Parameters
        ----------
        `**kwds` : optional
            Additional keyword arguments are documented in
            :meth:`pandas.Series.plot`.

        Returns
        -------
        :class:`matplotlib.axes.Axes` or numpy.ndarray of them
        """
        return self(kind='bar', **kwds)

    def barh(self, **kwds):
        """
        Horizontal bar plot.

        Parameters
        ----------
        `**kwds` : optional
            Additional keyword arguments are documented in
            :meth:`pandas.Series.plot`.

        Returns
        -------
        :class:`matplotlib.axes.Axes` or numpy.ndarray of them
        """
        return self(kind='barh', **kwds)

    def box(self, **kwds):
        """
        Boxplot.

        Parameters
        ----------
        `**kwds` : optional
            Additional keyword arguments are documented in
            :meth:`pandas.Series.plot`.

        Returns
        -------
        :class:`matplotlib.axes.Axes` or numpy.ndarray of them
        """
        return self(kind='box', **kwds)

    def hist(self, bins=10, **kwds):
        """
        Histogram.

        Parameters
        ----------
        bins : integer, default 10
            Number of histogram bins to be used
        `**kwds` : optional
            Additional keyword arguments are documented in
            :meth:`pandas.Series.plot`.

        Returns
        -------
        :class:`matplotlib.axes.Axes` or numpy.ndarray of them
        """
        return self(kind='hist', bins=bins, **kwds)

    @Appender(_shared_docs['kde'] % {
        'this-datatype': 'Series',
        'sibling-datatype': 'DataFrame',
        'examples': """
        Given a Series of points randomly sampled from an unknown
        distribution, estimate its PDF using KDE with automatic
        bandwidth determination and plot the results, evaluating them at
        1000 equally spaced points (default):

        .. plot::
            :context: close-figs

            >>> s = pd.Series([1, 2, 2.5, 3, 3.5, 4, 5])
            >>> ax = s.plot.kde()

        A scalar bandwidth can be specified. Using a small bandwidth value can
        lead to over-fitting, while using a large bandwidth value may result
        in under-fitting:

        .. plot::
            :context: close-figs

            >>> ax = s.plot.kde(bw_method=0.3)

        .. plot::
            :context: close-figs

            >>> ax = s.plot.kde(bw_method=3)

        Finally, the `ind` parameter determines the evaluation points for the
        plot of the estimated PDF:

        .. plot::
            :context: close-figs

            >>> ax = s.plot.kde(ind=[1, 2, 3, 4, 5])
        """.strip()
    })
    def kde(self, bw_method=None, ind=None, **kwds):
        return self(kind='kde', bw_method=bw_method, ind=ind, **kwds)

    density = kde

    def area(self, **kwds):
        """
        Area plot.

        Parameters
        ----------
        `**kwds` : optional
            Additional keyword arguments are documented in
            :meth:`pandas.Series.plot`.

        Returns
        -------
        :class:`matplotlib.axes.Axes` or numpy.ndarray of them
        """
        return self(kind='area', **kwds)

    def pie(self, **kwds):
        """
        Pie chart.

        Parameters
        ----------
        `**kwds` : optional
            Additional keyword arguments are documented in
            :meth:`pandas.Series.plot`.

        Returns
        -------
        :class:`matplotlib.axes.Axes` or numpy.ndarray of them
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
        return plot_frame(self._parent, kind=kind, x=x, y=y, ax=ax,
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
        Plot DataFrame columns as lines.

        This function is useful to plot lines using DataFrame's values
        as coordinates.

        Parameters
        ----------
        x : int or str, optional
            Columns to use for the horizontal axis.
            Either the location or the label of the columns to be used.
            By default, it will use the DataFrame indices.
        y : int, str, or list of them, optional
            The values to be plotted.
            Either the location or the label of the columns to be used.
            By default, it will use the remaining DataFrame numeric columns.
        **kwds
            Keyword arguments to pass on to :meth:`DataFrame.plot`.

        Returns
        -------
        :class:`matplotlib.axes.Axes` or :class:`numpy.ndarray`
            Return an ndarray when ``subplots=True``.

        See Also
        --------
        matplotlib.pyplot.plot : Plot y versus x as lines and/or markers.

        Examples
        --------

        .. plot::
            :context: close-figs

            The following example shows the populations for some animals
            over the years.

            >>> df = pd.DataFrame({
            ...    'pig': [20, 18, 489, 675, 1776],
            ...    'horse': [4, 25, 281, 600, 1900]
            ...    }, index=[1990, 1997, 2003, 2009, 2014])
            >>> lines = df.plot.line()

        .. plot::
           :context: close-figs

           An example with subplots, so an array of axes is returned.

           >>> axes = df.plot.line(subplots=True)
           >>> type(axes)
           <class 'numpy.ndarray'>

        .. plot::
            :context: close-figs

            The following example shows the relationship between both
            populations.

            >>> lines = df.plot.line(x='pig', y='horse')
        """
        return self(kind='line', x=x, y=y, **kwds)

    def bar(self, x=None, y=None, **kwds):
        """
        Vertical bar plot.

        A bar plot is a plot that presents categorical data with
        rectangular bars with lengths proportional to the values that they
        represent. A bar plot shows comparisons among discrete categories. One
        axis of the plot shows the specific categories being compared, and the
        other axis represents a measured value.

        Parameters
        ----------
        x : label or position, optional
            Allows plotting of one column versus another. If not specified,
            the index of the DataFrame is used.
        y : label or position, optional
            Allows plotting of one column versus another. If not specified,
            all numerical columns are used.
        **kwds
            Additional keyword arguments are documented in
            :meth:`DataFrame.plot`.

        Returns
        -------
        matplotlib.axes.Axes or np.ndarray of them
            An ndarray is returned with one :class:`matplotlib.axes.Axes`
            per column when ``subplots=True``.

        See Also
        --------
        DataFrame.plot.barh : Horizontal bar plot.
        DataFrame.plot : Make plots of a DataFrame.
        matplotlib.pyplot.bar : Make a bar plot with matplotlib.

        Examples
        --------
        Basic plot.

        .. plot::
            :context: close-figs

            >>> df = pd.DataFrame({'lab':['A', 'B', 'C'], 'val':[10, 30, 20]})
            >>> ax = df.plot.bar(x='lab', y='val', rot=0)

        Plot a whole dataframe to a bar plot. Each column is assigned a
        distinct color, and each row is nested in a group along the
        horizontal axis.

        .. plot::
            :context: close-figs

            >>> speed = [0.1, 17.5, 40, 48, 52, 69, 88]
            >>> lifespan = [2, 8, 70, 1.5, 25, 12, 28]
            >>> index = ['snail', 'pig', 'elephant',
            ...          'rabbit', 'giraffe', 'coyote', 'horse']
            >>> df = pd.DataFrame({'speed': speed,
            ...                    'lifespan': lifespan}, index=index)
            >>> ax = df.plot.bar(rot=0)

        Instead of nesting, the figure can be split by column with
        ``subplots=True``. In this case, a :class:`numpy.ndarray` of
        :class:`matplotlib.axes.Axes` are returned.

        .. plot::
            :context: close-figs

            >>> axes = df.plot.bar(rot=0, subplots=True)
            >>> axes[1].legend(loc=2)  # doctest: +SKIP

        Plot a single column.

        .. plot::
            :context: close-figs

            >>> ax = df.plot.bar(y='speed', rot=0)

        Plot only selected categories for the DataFrame.

        .. plot::
            :context: close-figs

            >>> ax = df.plot.bar(x='lifespan', rot=0)
        """
        return self(kind='bar', x=x, y=y, **kwds)

    def barh(self, x=None, y=None, **kwds):
        """
        Make a horizontal bar plot.

        A horizontal bar plot is a plot that presents quantitative data with
        rectangular bars with lengths proportional to the values that they
        represent. A bar plot shows comparisons among discrete categories. One
        axis of the plot shows the specific categories being compared, and the
        other axis represents a measured value.

        Parameters
        ----------
        x : label or position, default DataFrame.index
            Column to be used for categories.
        y : label or position, default All numeric columns in dataframe
            Columns to be plotted from the DataFrame.
        **kwds
            Keyword arguments to pass on to :meth:`DataFrame.plot`.

        Returns
        -------
        :class:`matplotlib.axes.Axes` or numpy.ndarray of them

        See Also
        --------
        DataFrame.plot.bar: Vertical bar plot.
        DataFrame.plot : Make plots of DataFrame using matplotlib.
        matplotlib.axes.Axes.bar : Plot a vertical bar plot using matplotlib.

        Examples
        --------
        Basic example

        .. plot::
            :context: close-figs

            >>> df = pd.DataFrame({'lab':['A', 'B', 'C'], 'val':[10, 30, 20]})
            >>> ax = df.plot.barh(x='lab', y='val')

        Plot a whole DataFrame to a horizontal bar plot

        .. plot::
            :context: close-figs

            >>> speed = [0.1, 17.5, 40, 48, 52, 69, 88]
            >>> lifespan = [2, 8, 70, 1.5, 25, 12, 28]
            >>> index = ['snail', 'pig', 'elephant',
            ...          'rabbit', 'giraffe', 'coyote', 'horse']
            >>> df = pd.DataFrame({'speed': speed,
            ...                    'lifespan': lifespan}, index=index)
            >>> ax = df.plot.barh()

        Plot a column of the DataFrame to a horizontal bar plot

        .. plot::
            :context: close-figs

            >>> speed = [0.1, 17.5, 40, 48, 52, 69, 88]
            >>> lifespan = [2, 8, 70, 1.5, 25, 12, 28]
            >>> index = ['snail', 'pig', 'elephant',
            ...          'rabbit', 'giraffe', 'coyote', 'horse']
            >>> df = pd.DataFrame({'speed': speed,
            ...                    'lifespan': lifespan}, index=index)
            >>> ax = df.plot.barh(y='speed')

        Plot DataFrame versus the desired column

        .. plot::
            :context: close-figs

            >>> speed = [0.1, 17.5, 40, 48, 52, 69, 88]
            >>> lifespan = [2, 8, 70, 1.5, 25, 12, 28]
            >>> index = ['snail', 'pig', 'elephant',
            ...          'rabbit', 'giraffe', 'coyote', 'horse']
            >>> df = pd.DataFrame({'speed': speed,
            ...                    'lifespan': lifespan}, index=index)
            >>> ax = df.plot.barh(x='lifespan')
        """
        return self(kind='barh', x=x, y=y, **kwds)

    def box(self, by=None, **kwds):
        r"""
        Make a box plot of the DataFrame columns.

        A box plot is a method for graphically depicting groups of numerical
        data through their quartiles.
        The box extends from the Q1 to Q3 quartile values of the data,
        with a line at the median (Q2). The whiskers extend from the edges
        of box to show the range of the data. The position of the whiskers
        is set by default to 1.5*IQR (IQR = Q3 - Q1) from the edges of the
        box. Outlier points are those past the end of the whiskers.

        For further details see Wikipedia's
        entry for `boxplot <https://en.wikipedia.org/wiki/Box_plot>`__.

        A consideration when using this chart is that the box and the whiskers
        can overlap, which is very common when plotting small sets of data.

        Parameters
        ----------
        by : string or sequence
            Column in the DataFrame to group by.
        **kwds : optional
            Additional keywords are documented in
            :meth:`DataFrame.plot`.

        Returns
        -------
        :class:`matplotlib.axes.Axes` or numpy.ndarray of them

        See Also
        --------
        DataFrame.boxplot: Another method to draw a box plot.
        Series.plot.box: Draw a box plot from a Series object.
        matplotlib.pyplot.boxplot: Draw a box plot in matplotlib.

        Examples
        --------
        Draw a box plot from a DataFrame with four columns of randomly
        generated data.

        .. plot::
            :context: close-figs

            >>> data = np.random.randn(25, 4)
            >>> df = pd.DataFrame(data, columns=list('ABCD'))
            >>> ax = df.plot.box()
        """
        return self(kind='box', by=by, **kwds)

    def hist(self, by=None, bins=10, **kwds):
        """
        Draw one histogram of the DataFrame's columns.

        A histogram is a representation of the distribution of data.
        This function groups the values of all given Series in the DataFrame
        into bins and draws all bins in one :class:`matplotlib.axes.Axes`.
        This is useful when the DataFrame's Series are in a similar scale.

        Parameters
        ----------
        by : str or sequence, optional
            Column in the DataFrame to group by.
        bins : int, default 10
            Number of histogram bins to be used.
        **kwds
            Additional keyword arguments are documented in
            :meth:`DataFrame.plot`.

        Returns
        -------
        class:`matplotlib.AxesSubplot`
            Return a histogram plot.

        See Also
        --------
        DataFrame.hist : Draw histograms per DataFrame's Series.
        Series.hist : Draw a histogram with Series' data.

        Examples
        --------
        When we draw a dice 6000 times, we expect to get each value around 1000
        times. But when we draw two dices and sum the result, the distribution
        is going to be quite different. A histogram illustrates those
        distributions.

        .. plot::
            :context: close-figs

            >>> df = pd.DataFrame(
            ...     np.random.randint(1, 7, 6000),
            ...     columns = ['one'])
            >>> df['two'] = df['one'] + np.random.randint(1, 7, 6000)
            >>> ax = df.plot.hist(bins=12, alpha=0.5)
        """
        return self(kind='hist', by=by, bins=bins, **kwds)

    @Appender(_shared_docs['kde'] % {
        'this-datatype': 'DataFrame',
        'sibling-datatype': 'Series',
        'examples': """
        Given several Series of points randomly sampled from unknown
        distributions, estimate their PDFs using KDE with automatic
        bandwidth determination and plot the results, evaluating them at
        1000 equally spaced points (default):

        .. plot::
            :context: close-figs

            >>> df = pd.DataFrame({
            ...     'x': [1, 2, 2.5, 3, 3.5, 4, 5],
            ...     'y': [4, 4, 4.5, 5, 5.5, 6, 6],
            ... })
            >>> ax = df.plot.kde()

        A scalar bandwidth can be specified. Using a small bandwidth value can
        lead to over-fitting, while using a large bandwidth value may result
        in under-fitting:

        .. plot::
            :context: close-figs

            >>> ax = df.plot.kde(bw_method=0.3)

        .. plot::
            :context: close-figs

            >>> ax = df.plot.kde(bw_method=3)

        Finally, the `ind` parameter determines the evaluation points for the
        plot of the estimated PDF:

        .. plot::
            :context: close-figs

            >>> ax = df.plot.kde(ind=[1, 2, 3, 4, 5, 6])
        """.strip()
    })
    def kde(self, bw_method=None, ind=None, **kwds):
        return self(kind='kde', bw_method=bw_method, ind=ind, **kwds)

    density = kde

    def area(self, x=None, y=None, **kwds):
        """
        Draw a stacked area plot.

        An area plot displays quantitative data visually.
        This function wraps the matplotlib area function.

        Parameters
        ----------
        x : label or position, optional
            Coordinates for the X axis. By default uses the index.
        y : label or position, optional
            Column to plot. By default uses all columns.
        stacked : bool, default True
            Area plots are stacked by default. Set to False to create a
            unstacked plot.
        **kwds : optional
            Additional keyword arguments are documented in
            :meth:`DataFrame.plot`.

        Returns
        -------
        matplotlib.axes.Axes or numpy.ndarray
            Area plot, or array of area plots if subplots is True.

        See Also
        --------
        DataFrame.plot : Make plots of DataFrame using matplotlib / pylab.

        Examples
        --------
        Draw an area plot based on basic business metrics:

        .. plot::
            :context: close-figs

            >>> df = pd.DataFrame({
            ...     'sales': [3, 2, 3, 9, 10, 6],
            ...     'signups': [5, 5, 6, 12, 14, 13],
            ...     'visits': [20, 42, 28, 62, 81, 50],
            ... }, index=pd.date_range(start='2018/01/01', end='2018/07/01',
            ...                        freq='M'))
            >>> ax = df.plot.area()

        Area plots are stacked by default. To produce an unstacked plot,
        pass ``stacked=False``:

        .. plot::
            :context: close-figs

            >>> ax = df.plot.area(stacked=False)

        Draw an area plot for a single column:

        .. plot::
            :context: close-figs

            >>> ax = df.plot.area(y='sales')

        Draw with a different `x`:

        .. plot::
            :context: close-figs

            >>> df = pd.DataFrame({
            ...     'sales': [3, 2, 3],
            ...     'visits': [20, 42, 28],
            ...     'day': [1, 2, 3],
            ... })
            >>> ax = df.plot.area(x='day')
        """
        return self(kind='area', x=x, y=y, **kwds)

    def pie(self, y=None, **kwds):
        """
        Generate a pie plot.

        A pie plot is a proportional representation of the numerical data in a
        column. This function wraps :meth:`matplotlib.pyplot.pie` for the
        specified column. If no column reference is passed and
        ``subplots=True`` a pie plot is drawn for each numerical column
        independently.

        Parameters
        ----------
        y : int or label, optional
            Label or position of the column to plot.
            If not provided, ``subplots=True`` argument must be passed.
        **kwds
            Keyword arguments to pass on to :meth:`DataFrame.plot`.

        Returns
        -------
        matplotlib.axes.Axes or np.ndarray of them
            A NumPy array is returned when `subplots` is True.

        See Also
        --------
        Series.plot.pie : Generate a pie plot for a Series.
        DataFrame.plot : Make plots of a DataFrame.

        Examples
        --------
        In the example below we have a DataFrame with the information about
        planet's mass and radius. We pass the the 'mass' column to the
        pie function to get a pie plot.

        .. plot::
            :context: close-figs

            >>> df = pd.DataFrame({'mass': [0.330, 4.87 , 5.97],
            ...                    'radius': [2439.7, 6051.8, 6378.1]},
            ...                   index=['Mercury', 'Venus', 'Earth'])
            >>> plot = df.plot.pie(y='mass', figsize=(5, 5))

        .. plot::
            :context: close-figs

            >>> plot = df.plot.pie(subplots=True, figsize=(6, 3))
        """
        return self(kind='pie', y=y, **kwds)

    def scatter(self, x, y, s=None, c=None, **kwds):
        """
        Create a scatter plot with varying marker point size and color.

        The coordinates of each point are defined by two dataframe columns and
        filled circles are used to represent each point. This kind of plot is
        useful to see complex correlations between two variables. Points could
        be for instance natural 2D coordinates like longitude and latitude in
        a map or, in general, any pair of metrics that can be plotted against
        each other.

        Parameters
        ----------
        x : int or str
            The column name or column position to be used as horizontal
            coordinates for each point.
        y : int or str
            The column name or column position to be used as vertical
            coordinates for each point.
        s : scalar or array_like, optional
            The size of each point. Possible values are:

            - A single scalar so all points have the same size.

            - A sequence of scalars, which will be used for each point's size
              recursively. For instance, when passing [2,14] all points size
              will be either 2 or 14, alternatively.

        c : str, int or array_like, optional
            The color of each point. Possible values are:

            - A single color string referred to by name, RGB or RGBA code,
              for instance 'red' or '#a98d19'.

            - A sequence of color strings referred to by name, RGB or RGBA
              code, which will be used for each point's color recursively. For
              instance ['green','yellow'] all points will be filled in green or
              yellow, alternatively.

            - A column name or position whose values will be used to color the
              marker points according to a colormap.

        **kwds
            Keyword arguments to pass on to :meth:`DataFrame.plot`.

        Returns
        -------
        :class:`matplotlib.axes.Axes` or numpy.ndarray of them

        See Also
        --------
        matplotlib.pyplot.scatter : Scatter plot using multiple input data
            formats.

        Examples
        --------
        Let's see how to draw a scatter plot using coordinates from the values
        in a DataFrame's columns.

        .. plot::
            :context: close-figs

            >>> df = pd.DataFrame([[5.1, 3.5, 0], [4.9, 3.0, 0], [7.0, 3.2, 1],
            ...                    [6.4, 3.2, 1], [5.9, 3.0, 2]],
            ...                   columns=['length', 'width', 'species'])
            >>> ax1 = df.plot.scatter(x='length',
            ...                       y='width',
            ...                       c='DarkBlue')

        And now with the color determined by a column as well.

        .. plot::
            :context: close-figs

            >>> ax2 = df.plot.scatter(x='length',
            ...                       y='width',
            ...                       c='species',
            ...                       colormap='viridis')
        """
        return self(kind='scatter', x=x, y=y, c=c, s=s, **kwds)

    def hexbin(self, x, y, C=None, reduce_C_function=None, gridsize=None,
               **kwds):
        """
        Generate a hexagonal binning plot.

        Generate a hexagonal binning plot of `x` versus `y`. If `C` is `None`
        (the default), this is a histogram of the number of occurrences
        of the observations at ``(x[i], y[i])``.

        If `C` is specified, specifies values at given coordinates
        ``(x[i], y[i])``. These values are accumulated for each hexagonal
        bin and then reduced according to `reduce_C_function`,
        having as default the NumPy's mean function (:meth:`numpy.mean`).
        (If `C` is specified, it must also be a 1-D sequence
        of the same length as `x` and `y`, or a column label.)

        Parameters
        ----------
        x : int or str
            The column label or position for x points.
        y : int or str
            The column label or position for y points.
        C : int or str, optional
            The column label or position for the value of `(x, y)` point.
        reduce_C_function : callable, default `np.mean`
            Function of one argument that reduces all the values in a bin to
            a single number (e.g. `np.mean`, `np.max`, `np.sum`, `np.std`).
        gridsize : int or tuple of (int, int), default 100
            The number of hexagons in the x-direction.
            The corresponding number of hexagons in the y-direction is
            chosen in a way that the hexagons are approximately regular.
            Alternatively, gridsize can be a tuple with two elements
            specifying the number of hexagons in the x-direction and the
            y-direction.
        **kwds
            Additional keyword arguments are documented in
            :meth:`DataFrame.plot`.

        Returns
        -------
        matplotlib.AxesSubplot
            The matplotlib ``Axes`` on which the hexbin is plotted.

        See Also
        --------
        DataFrame.plot : Make plots of a DataFrame.
        matplotlib.pyplot.hexbin : Hexagonal binning plot using matplotlib,
            the matplotlib function that is used under the hood.

        Examples
        --------
        The following examples are generated with random data from
        a normal distribution.

        .. plot::
            :context: close-figs

            >>> n = 10000
            >>> df = pd.DataFrame({'x': np.random.randn(n),
            ...                    'y': np.random.randn(n)})
            >>> ax = df.plot.hexbin(x='x', y='y', gridsize=20)

        The next example uses `C` and `np.sum` as `reduce_C_function`.
        Note that `'observations'` values ranges from 1 to 5 but the result
        plot shows values up to more than 25. This is because of the
        `reduce_C_function`.

        .. plot::
            :context: close-figs

            >>> n = 500
            >>> df = pd.DataFrame({
            ...     'coord_x': np.random.uniform(-3, 3, size=n),
            ...     'coord_y': np.random.uniform(30, 50, size=n),
            ...     'observations': np.random.randint(1,5, size=n)
            ...     })
            >>> ax = df.plot.hexbin(x='coord_x',
            ...                     y='coord_y',
            ...                     C='observations',
            ...                     reduce_C_function=np.sum,
            ...                     gridsize=10,
            ...                     cmap="viridis")
        """
        if reduce_C_function is not None:
            kwds['reduce_C_function'] = reduce_C_function
        if gridsize is not None:
            kwds['gridsize'] = gridsize
        return self(kind='hexbin', x=x, y=y, C=C, **kwds)
