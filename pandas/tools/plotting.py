from pandas.util.decorators import cache_readonly
import pandas.core.common as com

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

    def __init__(self, data, kind=None, by=None, subplots=False, sharex=True,
                 sharey=False, use_index=True,
                 figsize=None, grid=True, legend=True, rot=None,
                 ax=None, fig=None, title=None,
                 xlim=None, ylim=None,
                 xticks=None, yticks=None,
                 sort_columns=True, fontsize=None, **kwds):

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

        self.grid = grid
        self.legend = legend

        self.ax = ax
        self.fig = fig
        self.axes = None

        self.kwds = kwds

        self._args_adjust()
        self._compute_plot_data()
        self._setup_subplots()
        self._make_plot()
        self._post_plot_logic()
        self._adorn_subplots()

    def draw(self):
        self.plt.draw_if_interactive()

    def _args_adjust(self):
        pass

    def _setup_subplots(self):
        if self.subplots:
            nrows, ncols = self._get_layout()
            fig, axes = _subplots(nrows=nrows, ncols=ncols,
                                  sharex=self.sharex, sharey=self.sharey,
                                  figsize=self.figsize)
        else:
            if self.ax is None:
                fig = self.plt.figure(figsize=self.figsize)
                self.ax = fig.add_subplot(111)
            else:
                fig = self.ax.get_figure()

            axes = [self.ax]

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
        if self.subplots:
            to_adorn = self.axes
        else:
            to_adorn = [self.ax]

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

        if self.legend and not self.subplots:
            self.ax.legend(loc='best')

        if self.title:
            if self.subplots:
                self.fig.suptitle(self.title)
            else:
                self.ax.set_title(self.title)

    @cache_readonly
    def plt(self):
        import matplotlib.pyplot as plt
        return plt


class LinePlot(MPLPlot):

    def _make_plot(self):
        df = self.data

        if self.use_index:
            if df.index.is_numeric() or df.index.is_datetype():
                """
                Matplotlib supports numeric values or datetime objects as
                xaxis values. Taking LBYL approach here, by the time
                matplotlib raises exception when using non numeric/datetime
                values for xaxis, several actions are already taken by plt.
                """
                need_to_set_xticklabels = False
                x = df.index
            else:
                need_to_set_xticklabels = True
                x = range(len(df))
        else:
            need_to_set_xticklabels = False
            x = range(len(df))

        if self.sort_columns:
            columns = com._try_sort(df.columns)
        else:
            columns = df.columns

        for i, col in enumerate(columns):
            empty = df[col].count() == 0
            y = df[col].values if not empty else np.zeros(x.shape)

            if self.subplots:
                ax = self.axes[i]

                # kind of a hack
                ax.plot(x, y, 'k', label=str(col), **self.kwds)
                ax.legend(loc='best')
            else:
                ax = self.ax
                ax.plot(x, y, label=str(col), **self.kwds)

            ax.grid(self.grid)

        if need_to_set_xticklabels:
            xticklabels = [_stringify(key) for key in df.index]
            for ax_ in self.axes:
                ax_.set_xticks(x)
                ax_.set_xticklabels(xticklabels, rotation=self.rot)

    def _post_plot_logic(self):
        df = self.data

        condition = (df.index.is_all_dates
                     and not self.subplots
                     or (self.subplots and self.sharex))
        if condition:
            try:
                self.fig.autofmt_xdate()
            except Exception:  # pragma: no cover
                pass


class BarPlot(MPLPlot):
    _default_rot = {'bar' : 90, 'barh' : 0}

    def __init__(self, data, **kwargs):
        self.stacked = kwargs.pop('stacked', False)
        self.ax_pos = np.arange(len(data)) + 0.25
        MPLPlot.__init__(self, data, **kwargs)

    def _args_adjust(self):
        if self.rot is None:
            self.rot = self._default_rot[self.kind]

        if self.fontsize is None:
            if len(self.data) < 10:
                self.fontsize = 12
            else:
                self.fontsize = 10

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
        df = self.data

        N, K = df.shape

        colors = 'rgbyk'
        rects = []
        labels = []

        ax = self.axes[0]

        bar_f = self.bar_f

        pos_prior = neg_prior = np.zeros(N)
        for i, col in enumerate(df.columns):
            empty = df[col].count() == 0
            y = df[col].values if not empty else np.zeros(len(df))

            if self.subplots:
                ax = self.axes[i]
                rect = bar_f(ax, self.ax_pos, y, 0.5, start=pos_prior,
                             linewidth=1, **self.kwds)
                ax.set_title(col)
            elif self.stacked:
                mask = y > 0
                start = np.where(mask, pos_prior, neg_prior)
                rect = bar_f(ax, self.ax_pos, y, 0.5, start=start,
                             color=colors[i % len(colors)],
                             label=str(col), linewidth=1,
                             **self.kwds)
                pos_prior = pos_prior + np.where(mask, y, 0)
                neg_prior = neg_prior + np.where(mask, 0, y)
            else:
                rect = bar_f(ax, self.ax_pos + i * 0.75 / K, y, 0.75 / K,
                             start=np.zeros(N), label=str(col),
                             color=colors[i % len(colors)],
                             **self.kwds)
            rects.append(rect)
            labels.append(col)

        if self.legend and not self.subplots:
            patches =[r[0] for r in rects]

            # Legend to the right of the plot
            # ax.legend(patches, labels, bbox_to_anchor=(1.05, 1),
            #           loc=2, borderaxespad=0.)
            # self.fig.subplots_adjust(right=0.80)

            ax.legend(patches, labels, loc='best')


        self.fig.subplots_adjust(top=0.8)

    def _post_plot_logic(self):
        for ax in self.axes:
            str_index = [_stringify(key) for key in self.data.index]
            if self.kind == 'bar':
                ax.set_xlim([self.ax_pos[0] - 0.25, self.ax_pos[-1] + 1])
                ax.set_xticks(self.ax_pos + 0.375)
                ax.set_xticklabels(str_index, rotation=self.rot,
                                   fontsize=self.fontsize)
                ax.axhline(0, color='k', linestyle='--')
            else:
                # horizontal bars
                ax.set_ylim([self.ax_pos[0] - 0.25, self.ax_pos[-1] + 1])
                ax.set_yticks(self.ax_pos + 0.375)
                ax.set_yticklabels(str_index, rotation=self.rot,
                                   fontsize=self.fontsize)
                ax.axvline(0, color='k', linestyle='--')

class BoxPlot(MPLPlot):
    pass


class HistPlot(MPLPlot):
    pass


def plot_frame(frame=None, subplots=False, sharex=True, sharey=False,
               use_index=True,
               figsize=None, grid=True, legend=True, rot=None,
               ax=None, title=None,
               xlim=None, ylim=None,
               xticks=None, yticks=None,
               kind='line',
               sort_columns=True, fontsize=None, **kwds):
    """
    Make line plot of DataFrame's series with the index on the x-axis using
    matplotlib / pylab.

    Parameters
    ----------
    subplots : boolean, default False
        Make separate subplots for each time series
    sharex : boolean, default True
        In case subplots=True, share x axis
    sharey : boolean, default False
        In case subplots=True, share y axis
    use_index : boolean, default True
        Use index as ticks for x axis
    kind : {'line', 'bar'}
    sort_columns: boolean, default True
        Sort column names to determine plot ordering
    kwds : keywords
        Options to pass to Axis.plot

    Notes
    -----
    This method doesn't make much sense for cross-sections,
    and will error.
    """
    if kind == 'line':
        klass = LinePlot
    elif kind in ('bar', 'barh'):
        klass = BarPlot

    plot_obj = klass(frame, kind=kind, subplots=subplots, rot=rot,
                     legend=legend, ax=ax, fontsize=fontsize,
                     use_index=use_index, sharex=sharex, sharey=sharey,
                     xticks=xticks, yticks=yticks, xlim=xlim, ylim=ylim,
                     title=title, grid=grid, figsize=figsize,
                     sort_columns=sort_columns, **kwds)
    plot_obj.draw()
    if subplots:
        return plot_obj.axes
    else:
        return plot_obj.axes[0]


def plot_series(series, label=None, kind='line', use_index=True, rot=30,
                ax=None, style='-', grid=True, logy=False, **kwds):
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
    style : string, default '-'
        matplotlib line style to use
    kwds : keywords
        To be passed to the actual plotting function

    Notes
    -----
    See matplotlib documentation online for more on this subject
    Intended to be used in ipython --pylab mode
    """
    import matplotlib.pyplot as plt

    if label is not None:
        kwds = kwds.copy()
        kwds['label'] = label

    N = len(series)

    if ax is None:
        ax = plt.gca()

    if kind == 'line':
        if use_index:
            if series.index.is_numeric() or series.index.is_datetype():
                """
                Matplotlib supports numeric values or datetime objects as
                xaxis values. Taking LBYL approach here, by the time
                matplotlib raises exception when using non numeric/datetime
                values for xaxis, several actions are already taken by plt.
                """
                need_to_set_xticklabels = False
                x = np.asarray(series.index)
            else:
                need_to_set_xticklabels = True
                x = range(len(series))
        else:
            need_to_set_xticklabels = False
            x = range(len(series))

        if logy:
            ax.semilogy(x, series.values.astype(float), style, **kwds)
        else:
            ax.plot(x, series.values.astype(float), style, **kwds)
        format_date_labels(ax)

        if need_to_set_xticklabels:
            ax.set_xticks(x)
            ax.set_xticklabels([_stringify(key) for key in series.index],
                               rotation=rot)
    elif kind == 'bar':
        xinds = np.arange(N) + 0.25
        ax.bar(xinds, series.values.astype(float), 0.5,
               bottom=np.zeros(N), linewidth=1, **kwds)

        if N < 10:
            fontsize = 12
        else:
            fontsize = 10

        ax.set_xticks(xinds + 0.25)
        ax.set_xticklabels([_stringify(key) for key in series.index],
                           rotation=rot,
                           fontsize=fontsize)
    elif kind == 'barh':
        yinds = np.arange(N) + 0.25
        ax.barh(yinds, series.values.astype(float), 0.5,
                left=np.zeros(N), linewidth=1, **kwds)

        if N < 10:
            fontsize = 12
        else:
            fontsize = 10

        ax.set_yticks(yinds + 0.25)
        ax.set_yticklabels([_stringify(key) for key in series.index],
                           rotation=rot,
                           fontsize=fontsize)

    ax.grid(grid)
    plt.draw_if_interactive()

    return ax


def boxplot(data, column=None, by=None, ax=None, fontsize=None,
            rot=0, grid=True, figsize=None):
    """
    Make a box plot from DataFrame column optionally grouped b ysome columns or
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


def hist_frame(data, grid=True, **kwds):
    """
    Draw Histogram the DataFrame's series using matplotlib / pylab.

    Parameters
    ----------
    kwds : other plotting keyword arguments
        To be passed to hist function
    """
    n = len(data.columns)
    k = 1
    while k ** 2 < n:
        k += 1
    _, axes = _subplots(nrows=k, ncols=k)

    for i, col in enumerate(com._try_sort(data.columns)):
        ax = axes[i / k][i % k]
        ax.hist(data[col].dropna().values, **kwds)
        ax.set_title(col)
        ax.grid(grid)

    return axes


def hist_series(self, ax=None, grid=True, **kwds):
    """
    Draw histogram of the input series using matplotlib

    Parameters
    ----------
    ax : matplotlib axis object
        If not passed, uses gca()
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

    return ax


def _grouped_plot(plotf, data, column=None, by=None, numeric_only=True,
                  figsize=None, sharex=True, sharey=True, layout=None,
                  rot=0):
    from pandas.core.frame import DataFrame

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
                          sharex=sharex, sharey=sharey)

    ravel_axes = []
    for row in axes:
        ravel_axes.extend(row)

    for i, (key, group) in enumerate(grouped):
        ax = ravel_axes[i]
        if numeric_only and isinstance(group, DataFrame):
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
    fig, axes = _subplots(nrows=nrows, ncols=ncols,
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

def _subplots(nrows=1, ncols=1, sharex=False, sharey=False, squeeze=True,
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
