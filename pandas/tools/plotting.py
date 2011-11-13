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

def boxplot(data, column, by=None, ax=None, fontsize=None, rot=0):
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
    keys, values = zip(*data.groupby(by)[column])

    if ax is None:
        ax = plt.gca()
    ax.boxplot(values)
    ax.set_xticklabels(keys, rotation=rot, fontsize=fontsize)

    ax.set_xlabel(str(by))
    ax.set_ylabel(str(column))

    plt.subplots_adjust(bottom=0.15)
    return ax

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

def _grouped_plot(plotf, data, by=None):
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
        plotf(group, ax)
        ax.set_title(str(key))

    return fig

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
