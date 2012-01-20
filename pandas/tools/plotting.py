import numpy as np

def scatter_matrix(data):
    pass

def _gca():
    import matplotlib.pyplot as plt
    return plt.gca()

def _gcf():
    import matplotlib.pyplot as plt
    return plt.gcf()

def hist(data, column, by=None, ax=None, fontsize=None):
    keys, values = zip(*data.groupby(by)[column])
    if ax is None:
        ax = _gca()
    ax.boxplot(values)
    ax.set_xticklabels(keys, rotation=0, fontsize=fontsize)
    return ax

def grouped_hist(data, column, by=None, ax=None, bins=50, log=False,
                 figsize=None):
    """

    Returns
    -------
    fig : matplotlib.Figure
    """
    def plot_group(group, ax):
        ax.hist(group[column].dropna(), bins=bins)
    fig = _grouped_plot(plot_group, data, by=by, sharex=False,
                        sharey=False, figsize=figsize)
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.9,
                        hspace=0.3, wspace=0.2)
    return fig


def boxplot(data, column=None, by=None, ax=None, fontsize=None,
            rot=0, grid=True, figsize=None):
    """
    Make a box plot from DataFrame column optionally grouped by some columns or
    other inputs

    Parameters
    ----------
    data : DataFrame
    column : column name or list of names, or vector
        Can be any valid input to groupby
    by : string or sequence
        Column in the DataFrame to group by
    fontsize : int or string

    Returns
    -------
    ax : matplotlib.axes.AxesSubplot
    """
    def plot_group(grouped, ax):
        keys, values = zip(*grouped)
        keys = [_stringify(x) for x in keys]
        ax.boxplot(values)
        ax.set_xticklabels(keys, rotation=rot, fontsize=fontsize)

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
                                            by=by, grid=grid, figsize=figsize)
        ax = axes
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
        ax.boxplot(list(data[cols].values.T))
        ax.set_xticklabels(keys, rotation=rot, fontsize=fontsize)
        ax.grid(grid)

    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.9, wspace=0.2)
    return ax

def _stringify(x):
    if isinstance(x, tuple):
        return '|'.join(str(y) for y in x)
    else:
        return str(x)

def format_date_labels(ax):
    # mini version of autofmt_xdate
    try:
        for label in ax.get_xticklabels():
            label.set_ha('right')
            label.set_rotation(30)
        fig = ax.get_figure()
        fig.subplots_adjust(bottom=0.2)
    except Exception: # pragma: no cover
        pass


def scatter_plot(data, x, y, by=None, ax=None, figsize=None):
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

    if by is not None:
        fig = _grouped_plot(plot_group, data, by=by, figsize=figsize)
    else:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plot_group(data, ax)
        ax.set_ylabel(str(y))
        ax.set_xlabel(str(x))

    return fig

def _grouped_plot(plotf, data, by=None, numeric_only=True, figsize=None,
                  sharex=True, sharey=True):
    import matplotlib.pyplot as plt

    # allow to specify mpl default with 'default'
    if not (isinstance(figsize, str) and figsize == 'default'):
        figsize = (10, 5)               # our default

    grouped = data.groupby(by)
    ngroups = len(grouped)

    nrows, ncols = _get_layout(ngroups)
    if figsize is None:
        # our favorite default beating matplotlib's idea of the
        # default size
        figsize = (10, 5)
    fig, axes = subplots(nrows=nrows, ncols=ncols, figsize=figsize,
                         sharex=sharex, sharey=sharey)

    ravel_axes = []
    for row in axes:
        ravel_axes.extend(row)

    for i, (key, group) in enumerate(grouped):
        ax = ravel_axes[i]
        if numeric_only:
            group = group._get_numeric_data()
        plotf(group, ax)
        ax.set_title(str(key))

    return fig, axes

def _grouped_plot_by_column(plotf, data, columns=None, by=None,
                            numeric_only=True, grid=False,
                            figsize=None):
    import matplotlib.pyplot as plt

    grouped = data.groupby(by)
    if columns is None:
        columns = data._get_numeric_data().columns - by
    ngroups = len(columns)

    nrows, ncols = _get_layout(ngroups)
    fig, axes = subplots(nrows=nrows, ncols=ncols,
                         sharex=True, sharey=True,
                         figsize=figsize)

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
        ax.set_xlabel(str(by))
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

def subplots(nrows=1, ncols=1, sharex=False, sharey=False, squeeze=True,
              subplot_kw=None, **fig_kw):
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

    sharex : bool
      If True, the Y axis will be shared amongst all subplots.

    squeeze : bool

      If True, extra dimensions are squeezed out from the returned axis object:
        - if only one subplot is constructed (nrows=ncols=1), the resulting
        single Axis object is returned as a scalar.
        - for Nx1 or 1xN subplots, the returned object is a 1-d numpy object
        array of Axis objects are returned as numpy 1-d arrays.
        - for NxM subplots with N>1 and M>1 are returned as a 2d array.

      If False, no squeezing at all is done: the returned axis object is always
      a 2-d array contaning Axis instances, even if it ends up being 1x1.

    subplot_kw : dict
      Dict with keywords passed to the add_subplot() call used to create each
      subplots.

    fig_kw : dict
      Dict with keywords passed to the figure() call.  Note that all keywords
      not recognized above will be automatically included here.

    Returns:

    fig, ax : tuple
      - fig is the Matplotlib Figure object
      - ax can be either a single axis object or an array of axis objects if
      more than one supblot was created.  The dimensions of the resulting array
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

    fig = plt.figure(**fig_kw)

    # Create empty object array to hold all axes.  It's easiest to make it 1-d
    # so we can just append subplots upon creation, and then
    nplots = nrows*ncols
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
        axarr[i] = fig.add_subplot(nrows, ncols, i+1, **subplot_kw)

    if squeeze:
        # Reshape the array to have the final desired dimension (nrow,ncol),
        # though discarding unneeded dimensions that equal 1.  If we only have
        # one subplot, just return it instead of a 1-element array.
        if nplots==1:
            return fig, axarr[0]
        else:
            return fig, axarr.reshape(nrows, ncols).squeeze()
    else:
        # returned axis array will be always 2-d, even if nrows=ncols=1
        return fig, axarr.reshape(nrows, ncols)

if __name__ == '__main__':
    import pandas.rpy.common as com
    sales = com.load_data('sanfrancisco.home.sales', package='nutshell')
    top10 = sales['zip'].value_counts()[:10].index
    sales2 = sales[sales.zip.isin(top10)]

    fig = scatter_plot(sales2, 'squarefeet', 'price', by='zip')

    # plt.show()
