import matplotlib.pyplot as plt

def scatter_matrix(data):
    pass

def hist(data, column, by=None, ax=None, fontsize=None):
    keys, values = zip(*data.groupby(by)[column])

    if ax is None:
        ax = plt.gca()
    ax.boxplot(values)
    ax.set_xticklabels(keys, rotation=0, fontsize=fontsize)
    return ax

def boxplot(data, column=None, by=None, ax=None, fontsize=None,
            rot=0, grid=True):
    """
    Make a box plot from DataFrame column optionally grouped by some columns or
    other inputs

    Parameters
    ----------
    data : DataFrame
    column : column names or list of names, or vector
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

    if by is not None:
        if not isinstance(by, (list, tuple)):
            by = [by]

        columns = None if column is None else [column]
        fig, axes = _grouped_plot_by_column(plot_group, data, columns=columns,
                                            by=by)
        ax = axes
    else:
        if ax is None:
            ax = plt.gca()

        data = data._get_numeric_data()
        keys = [_stringify(x) for x in data.columns]
        ax.boxplot(list(data.values.T))
        ax.set_xticklabels(keys, rotation=rot, fontsize=fontsize)

    plt.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.9, wspace=0.1)
    return ax

def _stringify(x):
    if isinstance(x, tuple):
        return '|'.join(str(y) for y in x)
    else:
        return str(x)

def scatter_plot(data, x, y, by=None, ax=None):
    """

    Returns
    -------
    fig : matplotlib.Figure
    """
    def plot_group(group, ax):
        xvals = group[x].values
        yvals = group[y].values
        ax.scatter(xvals, yvals)

    if by is not None:
        fig = _grouped_plot(plot_group, data, by=by)
    else:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plot_group(data, ax)
        ax.set_ylabel(str(y))
        ax.set_xlabel(str(x))

    return fig

def _grouped_plot(plotf, data, by=None, numeric_only=True):
    grouped = data.groupby(by)
    ngroups = len(grouped)

    nrows, ncols = _get_layout(ngroups)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                             sharex=True, sharey=True)

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
                            numeric_only=True):
    grouped = data.groupby(by)
    if columns is None:
        columns = data._get_numeric_data().columns - by
    ngroups = len(columns)

    nrows, ncols = _get_layout(ngroups)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                             sharex=True, sharey=True)

    if isinstance(axes, plt.Axes):
        ravel_axes = [axes]
    else:
        ravel_axes = []
        for row in axes:
            ravel_axes.extend(row)

    for i, col in enumerate(columns):
        ax = ravel_axes[i]
        gp_col = grouped[col]
        plotf(gp_col, ax)
        ax.set_title(col)

    fig.suptitle('Boxplot grouped by %s' % by)

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

if __name__ == '__main__':
    import pandas.rpy.common as com
    sales = com.load_data('sanfrancisco.home.sales', package='nutshell')
    top10 = sales['zip'].value_counts()[:10].index
    sales2 = sales[sales.zip.isin(top10)]

    fig = scatter_plot(sales2, 'squarefeet', 'price', by='zip')

    plt.show()
