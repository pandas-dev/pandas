import numpy as np

from pandas.core.dtypes.common import is_integer, is_list_like
from pandas.core.dtypes.generic import ABCDataFrame, ABCIndexClass
from pandas.core.dtypes.missing import isna, remove_na_arraylike

import pandas.core.common as com

from pandas.io.formats.printing import pprint_thing
from pandas.plotting._matplotlib.core import LinePlot, MPLPlot
from pandas.plotting._matplotlib.tools import _flatten, _set_ticks_props, _subplots


class HistPlot(LinePlot):
    _kind = "hist"

    def __init__(self, data, bins=10, bottom=0, **kwargs):
        self.bins = bins  # use mpl default
        self.bottom = bottom

        # Do not call LinePlot.__init__ which may fill nan
        MPLPlot.__init__(self, data, **kwargs)

    def _args_adjust(self):

        # calculate bin number separately in different subplots
        # where subplots are created based on by argument
        if is_integer(self.bins):
            if self.by is None:
                self.bins = self._caculcate_bins(self.data)

            else:
                grouped = self.data.groupby(self.by)[self.column]
                bins_list = []
                for key, group in grouped:
                    bins_list.append(self._caculcate_bins(group))
                self.bins = bins_list

        if is_list_like(self.bottom):
            self.bottom = np.array(self.bottom)

    def _caculcate_bins(self, data):
        """Calculate bins given data"""

        values = data._convert(datetime=True)._get_numeric_data()
        values = np.ravel(values)
        values = values[~isna(values)]

        hist, bins = np.histogram(
            values,
            bins=self.bins,
            range=self.kwds.get("range", None),
            weights=self.kwds.get("weights", None),
        )
        return bins

    @classmethod
    def _plot(
        cls,
        ax,
        y,
        style=None,
        bins=None,
        bottom=0,
        column_num=0,
        stacking_id=None,
        **kwds,
    ):
        if column_num == 0:
            cls._initialize_stacker(ax, stacking_id, len(bins) - 1)
        y = y[~isna(y)]

        base = np.zeros(len(bins) - 1)
        bottom = bottom + cls._get_stacked_values(ax, stacking_id, base, kwds["label"])
        # ignore style
        n, bins, patches = ax.hist(y, bins=bins, bottom=bottom, **kwds)
        cls._update_stacker(ax, stacking_id, n)
        return patches

    @classmethod
    def _group_plot(
        cls,
        axes,
        data,
        fig,
        labels,
        bins=None,
        rot=90,
        xrot=None,
        xlabelsize=None,
        ylabelsize=None,
        yrot=None,
        **kwds,
    ):
        if "figure" in kwds:
            raise ValueError(
                "Cannot pass 'figure' when using the "
                "'by' argument, since a new 'Figure' instance "
                "will be created"
            )

        xrot = xrot or rot

        for i, (label, y) in enumerate(data):
            ax = axes[i]
            if len(y.shape) > 1:
                notna = [col[~isna(col)] for col in y.T]
                y_notna = np.array(np.array(notna).T)
            else:
                y_notna = y[~isna(y)]
            ax.hist(y_notna, bins[i], label=labels, **kwds)
            ax.set_title(pprint_thing(label))

        _set_ticks_props(
            axes, xlabelsize=xlabelsize, xrot=xrot, ylabelsize=ylabelsize, yrot=yrot
        )

        fig.subplots_adjust(
            bottom=0.15, top=0.9, left=0.1, right=0.9, hspace=0.8, wspace=0.3
        )
        return axes

    def _make_plot(self):
        colors = self._get_colors()
        stacking_id = self._get_stacking_id()
        if self.by is None:
            for i, (label, y) in enumerate(self._iter_data()):
                ax = self._get_ax(i)

                kwds = self.kwds.copy()
                label = pprint_thing(label)
                kwds["label"] = label

                style, kwds = self._apply_style_colors(colors, kwds, i, label)
                if style is not None:
                    kwds["style"] = style

                kwds = self._make_plot_keywords(kwds, y)
                artists = self._plot(
                    ax, y, column_num=i, stacking_id=stacking_id, **kwds
                )
                self._add_legend_handle(artists[0], label, index=i)

        else:
            kwds = self.kwds.copy()
            kwds = self._make_plot_keywords(kwds, None)
            data = self._iter_data()
            self._group_plot(self.axes, data, self.fig, self.column, **kwds)

    def _make_plot_keywords(self, kwds, y):
        """merge BoxPlot/KdePlot properties to passed kwds"""
        # y is required for KdePlot
        kwds["bottom"] = self.bottom
        kwds["bins"] = self.bins
        return kwds

    def _post_plot_logic(self, ax, data):
        if self.orientation == "horizontal":
            ax.set_xlabel("Frequency")
        else:
            ax.set_ylabel("Frequency")

    @property
    def orientation(self):
        if self.kwds.get("orientation", None) == "horizontal":
            return "horizontal"
        else:
            return "vertical"


class KdePlot(HistPlot):
    _kind = "kde"
    orientation = "vertical"

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
            ind = np.linspace(
                np.nanmin(y) - 0.5 * sample_range,
                np.nanmax(y) + 0.5 * sample_range,
                1000,
            )
        elif is_integer(self.ind):
            sample_range = np.nanmax(y) - np.nanmin(y)
            ind = np.linspace(
                np.nanmin(y) - 0.5 * sample_range,
                np.nanmax(y) + 0.5 * sample_range,
                self.ind,
            )
        else:
            ind = self.ind
        return ind

    @classmethod
    def _plot(
        cls,
        ax,
        y,
        style=None,
        bw_method=None,
        ind=None,
        column_num=None,
        stacking_id=None,
        **kwds,
    ):
        from scipy.stats import gaussian_kde

        y = remove_na_arraylike(y)
        gkde = gaussian_kde(y, bw_method=bw_method)

        y = gkde.evaluate(ind)
        lines = MPLPlot._plot(ax, ind, y, style=style, **kwds)
        return lines

    def _make_plot_keywords(self, kwds, y):
        kwds["bw_method"] = self.bw_method
        kwds["ind"] = self._get_ind(y)
        return kwds

    def _post_plot_logic(self, ax, data):
        ax.set_ylabel("Density")


def _grouped_plot(
    plotf,
    data,
    column=None,
    by=None,
    numeric_only=True,
    figsize=None,
    sharex=True,
    sharey=True,
    layout=None,
    rot=0,
    ax=None,
    **kwargs,
):

    if figsize == "default":
        # allowed to specify mpl default with 'default'
        raise ValueError(
            "figsize='default' is no longer supported. "
            "Specify figure size by tuple instead"
        )

    grouped = data.groupby(by)
    if column is not None:
        grouped = grouped[column]

    naxes = len(grouped)
    fig, axes = _subplots(
        naxes=naxes, figsize=figsize, sharex=sharex, sharey=sharey, ax=ax, layout=layout
    )

    _axes = _flatten(axes)

    for i, (key, group) in enumerate(grouped):
        ax = _axes[i]
        if numeric_only and isinstance(group, ABCDataFrame):
            group = group._get_numeric_data()
        plotf(group, ax, **kwargs)
        ax.set_title(pprint_thing(key))

    return fig, axes


def _grouped_hist(
    data,
    column=None,
    by=None,
    ax=None,
    bins=50,
    figsize=None,
    layout=None,
    sharex=False,
    sharey=False,
    rot=90,
    grid=True,
    xlabelsize=None,
    xrot=None,
    ylabelsize=None,
    yrot=None,
    **kwargs,
):
    """
    Grouped histogram

    Parameters
    ----------
    data : Series/DataFrame
    column : object, optional
    by : object, optional
    ax : axes, optional
    bins : int, default 50
    figsize : tuple, optional
    layout : optional
    sharex : bool, default False
    sharey : bool, default False
    rot : int, default 90
    grid : bool, default True
    kwargs : dict, keyword arguments passed to matplotlib.Axes.hist

    Returns
    -------
    collection of Matplotlib Axes
    """

    def plot_group(group, ax):
        ax.hist(group.dropna().values, bins=bins, **kwargs)

    if xrot is None:
        xrot = rot

    fig, axes = _grouped_plot(
        plot_group,
        data,
        column=column,
        by=by,
        sharex=sharex,
        sharey=sharey,
        ax=ax,
        figsize=figsize,
        layout=layout,
        rot=rot,
    )

    _set_ticks_props(
        axes, xlabelsize=xlabelsize, xrot=xrot, ylabelsize=ylabelsize, yrot=yrot
    )

    fig.subplots_adjust(
        bottom=0.15, top=0.9, left=0.1, right=0.9, hspace=0.5, wspace=0.3
    )
    return axes


def hist_series(
    self,
    by=None,
    ax=None,
    grid=True,
    xlabelsize=None,
    xrot=None,
    ylabelsize=None,
    yrot=None,
    figsize=None,
    bins=10,
    **kwds,
):
    import matplotlib.pyplot as plt

    if by is None:
        if kwds.get("layout", None) is not None:
            raise ValueError("The 'layout' keyword is not supported when 'by' is None")
        # hack until the plotting interface is a bit more unified
        fig = kwds.pop(
            "figure", plt.gcf() if plt.get_fignums() else plt.figure(figsize=figsize)
        )
        if figsize is not None and tuple(figsize) != tuple(fig.get_size_inches()):
            fig.set_size_inches(*figsize, forward=True)
        if ax is None:
            ax = fig.gca()
        elif ax.get_figure() != fig:
            raise AssertionError("passed axis not bound to passed figure")
        values = self.dropna().values

        ax.hist(values, bins=bins, **kwds)
        ax.grid(grid)
        axes = np.array([ax])

        _set_ticks_props(
            axes, xlabelsize=xlabelsize, xrot=xrot, ylabelsize=ylabelsize, yrot=yrot
        )

    else:
        if "figure" in kwds:
            raise ValueError(
                "Cannot pass 'figure' when using the "
                "'by' argument, since a new 'Figure' instance "
                "will be created"
            )
        axes = _grouped_hist(
            self,
            by=by,
            ax=ax,
            grid=grid,
            figsize=figsize,
            bins=bins,
            xlabelsize=xlabelsize,
            xrot=xrot,
            ylabelsize=ylabelsize,
            yrot=yrot,
            **kwds,
        )

    if hasattr(axes, "ndim"):
        if axes.ndim == 1 and len(axes) == 1:
            return axes[0]
    return axes


def hist_frame(
    data,
    column=None,
    by=None,
    grid=True,
    xlabelsize=None,
    xrot=None,
    ylabelsize=None,
    yrot=None,
    ax=None,
    sharex=False,
    sharey=False,
    figsize=None,
    layout=None,
    bins=10,
    **kwds,
):
    if by is not None:
        axes = _grouped_hist(
            data,
            column=column,
            by=by,
            ax=ax,
            grid=grid,
            figsize=figsize,
            sharex=sharex,
            sharey=sharey,
            layout=layout,
            bins=bins,
            xlabelsize=xlabelsize,
            xrot=xrot,
            ylabelsize=ylabelsize,
            yrot=yrot,
            **kwds,
        )
        return axes

    if column is not None:
        if not isinstance(column, (list, np.ndarray, ABCIndexClass)):
            column = [column]
        data = data[column]
    data = data._get_numeric_data()
    naxes = len(data.columns)

    if naxes == 0:
        raise ValueError("hist method requires numerical columns, nothing to plot.")

    fig, axes = _subplots(
        naxes=naxes,
        ax=ax,
        squeeze=False,
        sharex=sharex,
        sharey=sharey,
        figsize=figsize,
        layout=layout,
    )
    _axes = _flatten(axes)

    for i, col in enumerate(com.try_sort(data.columns)):
        ax = _axes[i]
        ax.hist(data[col].dropna().values, bins=bins, **kwds)
        ax.set_title(col)
        ax.grid(grid)

    _set_ticks_props(
        axes, xlabelsize=xlabelsize, xrot=xrot, ylabelsize=ylabelsize, yrot=yrot
    )
    fig.subplots_adjust(wspace=0.3, hspace=0.3)

    return axes
