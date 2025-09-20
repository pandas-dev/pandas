from __future__ import annotations

import functools
import warnings
from collections.abc import Callable, Hashable, Iterable, MutableMapping
from typing import TYPE_CHECKING, Any, Literal, Union, cast, overload

import numpy as np
import pandas as pd

from xarray.core.utils import attempt_import
from xarray.plot.facetgrid import _easy_facetgrid
from xarray.plot.utils import (
    _LINEWIDTH_RANGE,
    _MARKERSIZE_RANGE,
    _add_colorbar,
    _add_legend,
    _assert_valid_xy,
    _determine_guide,
    _ensure_plottable,
    _guess_coords_to_plot,
    _infer_interval_breaks,
    _infer_xy_labels,
    _Normalize,
    _process_cmap_cbar_kwargs,
    _rescale_imshow_rgb,
    _resolve_intervals_1dplot,
    _resolve_intervals_2dplot,
    _set_concise_date,
    _update_axes,
    get_axis,
    label_from_attrs,
)
from xarray.structure.alignment import broadcast
from xarray.structure.concat import concat

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.collections import PathCollection, QuadMesh
    from matplotlib.colors import Colormap, Normalize
    from matplotlib.container import BarContainer
    from matplotlib.contour import QuadContourSet
    from matplotlib.image import AxesImage
    from matplotlib.patches import Polygon
    from mpl_toolkits.mplot3d.art3d import Line3D, Poly3DCollection
    from numpy.typing import ArrayLike

    from xarray.core.dataarray import DataArray
    from xarray.core.types import (
        AspectOptions,
        ExtendOptions,
        HueStyleOptions,
        ScaleOptions,
        T_DataArray,
    )
    from xarray.plot.facetgrid import FacetGrid

_styles: dict[str, Any] = {
    # Add a white border to make it easier seeing overlapping markers:
    "scatter.edgecolors": "w",
}


def _infer_line_data(
    darray: DataArray, x: Hashable | None, y: Hashable | None, hue: Hashable | None
) -> tuple[DataArray, DataArray, DataArray | None, str]:
    ndims = len(darray.dims)

    if x is not None and y is not None:
        raise ValueError("Cannot specify both x and y kwargs for line plots.")

    if x is not None:
        _assert_valid_xy(darray, x, "x")

    if y is not None:
        _assert_valid_xy(darray, y, "y")

    if ndims == 1:
        huename = None
        hueplt = None
        huelabel = ""

        if x is not None:
            xplt = darray[x]
            yplt = darray

        elif y is not None:
            xplt = darray
            yplt = darray[y]

        else:  # Both x & y are None
            dim = darray.dims[0]
            xplt = darray[dim]
            yplt = darray

    else:
        if x is None and y is None and hue is None:
            raise ValueError("For 2D inputs, please specify either hue, x or y.")

        if y is None:
            if hue is not None:
                _assert_valid_xy(darray, hue, "hue")
            xname, huename = _infer_xy_labels(darray=darray, x=x, y=hue)
            xplt = darray[xname]
            if xplt.ndim > 1:
                if huename in darray.dims:
                    otherindex = 1 if darray.dims.index(huename) == 0 else 0
                    otherdim = darray.dims[otherindex]
                    yplt = darray.transpose(otherdim, huename, transpose_coords=False)
                    xplt = xplt.transpose(otherdim, huename, transpose_coords=False)
                else:
                    raise ValueError(
                        "For 2D inputs, hue must be a dimension"
                        " i.e. one of " + repr(darray.dims)
                    )

            else:
                (xdim,) = darray[xname].dims
                (huedim,) = darray[huename].dims
                yplt = darray.transpose(xdim, huedim)

        else:
            yname, huename = _infer_xy_labels(darray=darray, x=y, y=hue)
            yplt = darray[yname]
            if yplt.ndim > 1:
                if huename in darray.dims:
                    otherindex = 1 if darray.dims.index(huename) == 0 else 0
                    otherdim = darray.dims[otherindex]
                    xplt = darray.transpose(otherdim, huename, transpose_coords=False)
                    yplt = yplt.transpose(otherdim, huename, transpose_coords=False)
                else:
                    raise ValueError(
                        "For 2D inputs, hue must be a dimension"
                        " i.e. one of " + repr(darray.dims)
                    )

            else:
                (ydim,) = darray[yname].dims
                (huedim,) = darray[huename].dims
                xplt = darray.transpose(ydim, huedim)

        huelabel = label_from_attrs(darray[huename])
        hueplt = darray[huename]

    return xplt, yplt, hueplt, huelabel


def _prepare_plot1d_data(
    darray: T_DataArray,
    coords_to_plot: MutableMapping[str, Hashable],
    plotfunc_name: str | None = None,
    _is_facetgrid: bool = False,
) -> dict[str, T_DataArray]:
    """
    Prepare data for usage with plt.scatter.

    Parameters
    ----------
    darray : T_DataArray
        Base DataArray.
    coords_to_plot : MutableMapping[str, Hashable]
        Coords that will be plotted.
    plotfunc_name : str | None
        Name of the plotting function that will be used.

    Returns
    -------
    plts : dict[str, T_DataArray]
        Dict of DataArrays that will be sent to matplotlib.

    Examples
    --------
    >>> # Make sure int coords are plotted:
    >>> a = xr.DataArray(
    ...     data=[1, 2],
    ...     coords={1: ("x", [0, 1], {"units": "s"})},
    ...     dims=("x",),
    ...     name="a",
    ... )
    >>> plts = xr.plot.dataarray_plot._prepare_plot1d_data(
    ...     a, coords_to_plot={"x": 1, "z": None, "hue": None, "size": None}
    ... )
    >>> # Check which coords to plot:
    >>> print({k: v.name for k, v in plts.items()})
    {'y': 'a', 'x': 1}
    """
    # If there are more than 1 dimension in the array than stack all the
    # dimensions so the plotter can plot anything:
    if darray.ndim > 1:
        # When stacking dims the lines will continue connecting. For floats
        # this can be solved by adding a nan element in between the flattening
        # points:
        dims_T = []
        if np.issubdtype(darray.dtype, np.floating):
            for v in ["z", "x"]:
                dim = coords_to_plot.get(v, None)
                if (dim is not None) and (dim in darray.dims):
                    darray_nan = np.nan * darray.isel({dim: -1})
                    darray = concat(
                        [darray, darray_nan],
                        dim=dim,
                        coords="minimal",
                        compat="override",
                        join="exact",
                    )
                    dims_T.append(coords_to_plot[v])

        # Lines should never connect to the same coordinate when stacked,
        # transpose to avoid this as much as possible:
        darray = darray.transpose(..., *dims_T)

        # Array is now ready to be stacked:
        darray = darray.stack(_stacked_dim=darray.dims)

    # Broadcast together all the chosen variables:
    plts = dict(y=darray)
    plts.update(
        {k: darray.coords[v] for k, v in coords_to_plot.items() if v is not None}
    )
    plts = dict(zip(plts.keys(), broadcast(*(plts.values())), strict=True))

    return plts


# return type is Any due to the many different possibilities
def plot(
    darray: DataArray,
    *,
    row: Hashable | None = None,
    col: Hashable | None = None,
    col_wrap: int | None = None,
    ax: Axes | None = None,
    hue: Hashable | None = None,
    subplot_kws: dict[str, Any] | None = None,
    **kwargs: Any,
) -> Any:
    """
    Default plot of DataArray using :py:mod:`matplotlib:matplotlib.pyplot`.

    Calls xarray plotting function based on the dimensions of
    the squeezed DataArray.

    =============== ===========================
    Dimensions      Plotting function
    =============== ===========================
    1               :py:func:`xarray.plot.line`
    2               :py:func:`xarray.plot.pcolormesh`
    Anything else   :py:func:`xarray.plot.hist`
    =============== ===========================

    Parameters
    ----------
    darray : DataArray
    row : Hashable or None, optional
        If passed, make row faceted plots on this dimension name.
    col : Hashable or None, optional
        If passed, make column faceted plots on this dimension name.
    col_wrap : int or None, optional
        Use together with ``col`` to wrap faceted plots.
    ax : matplotlib axes object, optional
        Axes on which to plot. By default, use the current axes.
        Mutually exclusive with ``size``, ``figsize`` and facets.
    hue : Hashable or None, optional
        If passed, make faceted line plots with hue on this dimension name.
    subplot_kws : dict, optional
        Dictionary of keyword arguments for Matplotlib subplots
        (see :py:meth:`matplotlib:matplotlib.figure.Figure.add_subplot`).
    **kwargs : optional
        Additional keyword arguments for Matplotlib.

    See Also
    --------
    xarray.DataArray.squeeze
    """
    darray = darray.squeeze(
        d for d, s in darray.sizes.items() if s == 1 and d not in (row, col, hue)
    ).compute()

    plot_dims = set(darray.dims)
    plot_dims.discard(row)
    plot_dims.discard(col)
    plot_dims.discard(hue)

    ndims = len(plot_dims)

    plotfunc: Callable

    if ndims == 0 or darray.size == 0:
        raise TypeError("No numeric data to plot.")
    if ndims in (1, 2):
        if row or col:
            kwargs["subplot_kws"] = subplot_kws
            kwargs["row"] = row
            kwargs["col"] = col
            kwargs["col_wrap"] = col_wrap
        if ndims == 1:
            plotfunc = line
            kwargs["hue"] = hue
        elif ndims == 2:
            if hue:
                plotfunc = line
                kwargs["hue"] = hue
            else:
                plotfunc = pcolormesh
                kwargs["subplot_kws"] = subplot_kws
    else:
        if row or col or hue:
            raise ValueError(
                "Only 1d and 2d plots are supported for facets in xarray. "
                "See the package `Seaborn` for more options."
            )
        plotfunc = hist

    kwargs["ax"] = ax

    return plotfunc(darray, **kwargs)


@overload
def line(  # type: ignore[misc,unused-ignore]  # None is hashable :(
    darray: DataArray,
    *args: Any,
    row: None = None,  # no wrap -> primitive
    col: None = None,  # no wrap -> primitive
    figsize: Iterable[float] | None = None,
    aspect: AspectOptions = None,
    size: float | None = None,
    ax: Axes | None = None,
    hue: Hashable | None = None,
    x: Hashable | None = None,
    y: Hashable | None = None,
    xincrease: bool | None = None,
    yincrease: bool | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    add_legend: bool = True,
    _labels: bool = True,
    **kwargs: Any,
) -> list[Line3D]: ...


@overload
def line(
    darray: T_DataArray,
    *args: Any,
    row: Hashable,  # wrap -> FacetGrid
    col: Hashable | None = None,
    figsize: Iterable[float] | None = None,
    aspect: AspectOptions = None,
    size: float | None = None,
    ax: Axes | None = None,
    hue: Hashable | None = None,
    x: Hashable | None = None,
    y: Hashable | None = None,
    xincrease: bool | None = None,
    yincrease: bool | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    add_legend: bool = True,
    _labels: bool = True,
    **kwargs: Any,
) -> FacetGrid[T_DataArray]: ...


@overload
def line(
    darray: T_DataArray,
    *args: Any,
    row: Hashable | None = None,
    col: Hashable,  # wrap -> FacetGrid
    figsize: Iterable[float] | None = None,
    aspect: AspectOptions = None,
    size: float | None = None,
    ax: Axes | None = None,
    hue: Hashable | None = None,
    x: Hashable | None = None,
    y: Hashable | None = None,
    xincrease: bool | None = None,
    yincrease: bool | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    add_legend: bool = True,
    _labels: bool = True,
    **kwargs: Any,
) -> FacetGrid[T_DataArray]: ...


# This function signature should not change so that it can use
# matplotlib format strings
def line(
    darray: T_DataArray,
    *args: Any,
    row: Hashable | None = None,
    col: Hashable | None = None,
    figsize: Iterable[float] | None = None,
    aspect: AspectOptions = None,
    size: float | None = None,
    ax: Axes | None = None,
    hue: Hashable | None = None,
    x: Hashable | None = None,
    y: Hashable | None = None,
    xincrease: bool | None = None,
    yincrease: bool | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    add_legend: bool = True,
    _labels: bool = True,
    **kwargs: Any,
) -> list[Line3D] | FacetGrid[T_DataArray]:
    """
    Line plot of DataArray values.

    Wraps :py:func:`matplotlib:matplotlib.pyplot.plot`.

    Parameters
    ----------
    darray : DataArray
        Either 1D or 2D. If 2D, one of ``hue``, ``x`` or ``y`` must be provided.
    row : Hashable, optional
        If passed, make row faceted plots on this dimension name.
    col : Hashable, optional
        If passed, make column faceted plots on this dimension name.
    figsize : tuple, optional
        A tuple (width, height) of the figure in inches.
        Mutually exclusive with ``size`` and ``ax``.
    aspect : "auto", "equal", scalar or None, optional
        Aspect ratio of plot, so that ``aspect * size`` gives the *width* in
        inches. Only used if a ``size`` is provided.
    size : scalar, optional
        If provided, create a new figure for the plot with the given size:
        *height* (in inches) of each plot. See also: ``aspect``.
    ax : matplotlib axes object, optional
        Axes on which to plot. By default, the current is used.
        Mutually exclusive with ``size`` and ``figsize``.
    hue : Hashable, optional
        Dimension or coordinate for which you want multiple lines plotted.
        If plotting against a 2D coordinate, ``hue`` must be a dimension.
    x, y : Hashable, optional
        Dimension, coordinate or multi-index level for *x*, *y* axis.
        Only one of these may be specified.
        The other will be used for values from the DataArray on which this
        plot method is called.
    xincrease : bool or None, optional
        Should the values on the *x* axis be increasing from left to right?
        if ``None``, use the default for the Matplotlib function.
    yincrease : bool or None, optional
        Should the values on the *y* axis be increasing from top to bottom?
        if ``None``, use the default for the Matplotlib function.
    xscale, yscale : {'linear', 'symlog', 'log', 'logit'}, optional
        Specifies scaling for the *x*- and *y*-axis, respectively.
    xticks, yticks : array-like, optional
        Specify tick locations for *x*- and *y*-axis.
    xlim, ylim : tuple[float, float], optional
        Specify *x*- and *y*-axis limits.
    add_legend : bool, default: True
        Add legend with *y* axis coordinates (2D inputs only).
    *args, **kwargs : optional
        Additional arguments to :py:func:`matplotlib:matplotlib.pyplot.plot`.

    Returns
    -------
    primitive : list of Line3D or FacetGrid
        When either col or row is given, returns a FacetGrid, otherwise
        a list of matplotlib Line3D objects.
    """
    # Handle facetgrids first
    if row or col:
        allargs = locals().copy()
        allargs.update(allargs.pop("kwargs"))
        allargs.pop("darray")
        return _easy_facetgrid(darray, line, kind="line", **allargs)

    ndims = len(darray.dims)
    if ndims == 0 or darray.size == 0:
        # TypeError to be consistent with pandas
        raise TypeError("No numeric data to plot.")
    if ndims > 2:
        raise ValueError(
            "Line plots are for 1- or 2-dimensional DataArrays. "
            f"Passed DataArray has {ndims} "
            "dimensions"
        )

    # The allargs dict passed to _easy_facetgrid above contains args
    if args == ():
        args = kwargs.pop("args", ())
    else:
        assert "args" not in kwargs

    ax = get_axis(figsize, size, aspect, ax)
    xplt, yplt, hueplt, hue_label = _infer_line_data(darray, x, y, hue)

    # Remove pd.Intervals if contained in xplt.values and/or yplt.values.
    xplt_val, yplt_val, x_suffix, y_suffix, kwargs = _resolve_intervals_1dplot(
        xplt.to_numpy(), yplt.to_numpy(), kwargs
    )
    xlabel = label_from_attrs(xplt, extra=x_suffix)
    ylabel = label_from_attrs(yplt, extra=y_suffix)

    _ensure_plottable(xplt_val, yplt_val)

    primitive = ax.plot(xplt_val, yplt_val, *args, **kwargs)

    if _labels:
        if xlabel is not None:
            ax.set_xlabel(xlabel)

        if ylabel is not None:
            ax.set_ylabel(ylabel)

        ax.set_title(darray._title_for_slice())

    if darray.ndim == 2 and add_legend:
        assert hueplt is not None
        ax.legend(handles=primitive, labels=list(hueplt.to_numpy()), title=hue_label)

    if isinstance(xplt.dtype, np.dtype) and np.issubdtype(xplt.dtype, np.datetime64):  # type: ignore[redundant-expr]
        _set_concise_date(ax, axis="x")

    _update_axes(ax, xincrease, yincrease, xscale, yscale, xticks, yticks, xlim, ylim)

    return primitive


@overload
def step(  # type: ignore[misc,unused-ignore]  # None is hashable :(
    darray: DataArray,
    *args: Any,
    where: Literal["pre", "post", "mid"] = "pre",
    drawstyle: str | None = None,
    ds: str | None = None,
    row: None = None,  # no wrap -> primitive
    col: None = None,  # no wrap -> primitive
    **kwargs: Any,
) -> list[Line3D]: ...


@overload
def step(
    darray: DataArray,
    *args: Any,
    where: Literal["pre", "post", "mid"] = "pre",
    drawstyle: str | None = None,
    ds: str | None = None,
    row: Hashable,  # wrap -> FacetGrid
    col: Hashable | None = None,
    **kwargs: Any,
) -> FacetGrid[DataArray]: ...


@overload
def step(
    darray: DataArray,
    *args: Any,
    where: Literal["pre", "post", "mid"] = "pre",
    drawstyle: str | None = None,
    ds: str | None = None,
    row: Hashable | None = None,
    col: Hashable,  # wrap -> FacetGrid
    **kwargs: Any,
) -> FacetGrid[DataArray]: ...


def step(
    darray: DataArray,
    *args: Any,
    where: Literal["pre", "post", "mid"] = "pre",
    drawstyle: str | None = None,
    ds: str | None = None,
    row: Hashable | None = None,
    col: Hashable | None = None,
    **kwargs: Any,
) -> list[Line3D] | FacetGrid[DataArray]:
    """
    Step plot of DataArray values.

    Similar to :py:func:`matplotlib:matplotlib.pyplot.step`.

    Parameters
    ----------
    where : {'pre', 'post', 'mid'}, default: 'pre'
        Define where the steps should be placed:

        - ``'pre'``: The y value is continued constantly to the left from
          every *x* position, i.e. the interval ``(x[i-1], x[i]]`` has the
          value ``y[i]``.
        - ``'post'``: The y value is continued constantly to the right from
          every *x* position, i.e. the interval ``[x[i], x[i+1])`` has the
          value ``y[i]``.
        - ``'mid'``: Steps occur half-way between the *x* positions.

        Note that this parameter is ignored if one coordinate consists of
        :py:class:`pandas.Interval` values, e.g. as a result of
        :py:func:`xarray.Dataset.groupby_bins`. In this case, the actual
        boundaries of the interval are used.
    drawstyle, ds : str or None, optional
        Additional drawstyle. Only use one of drawstyle and ds.
    row : Hashable, optional
        If passed, make row faceted plots on this dimension name.
    col : Hashable, optional
        If passed, make column faceted plots on this dimension name.
    *args, **kwargs : optional
        Additional arguments for :py:func:`xarray.plot.line`.

    Returns
    -------
    primitive : list of Line3D or FacetGrid
        When either col or row is given, returns a FacetGrid, otherwise
        a list of matplotlib Line3D objects.
    """
    if where not in {"pre", "post", "mid"}:
        raise ValueError("'where' argument to step must be 'pre', 'post' or 'mid'")

    if ds is not None:
        if drawstyle is None:
            drawstyle = ds
        else:
            raise TypeError("ds and drawstyle are mutually exclusive")
    if drawstyle is None:
        drawstyle = ""
    drawstyle = "steps-" + where + drawstyle

    return line(darray, *args, drawstyle=drawstyle, col=col, row=row, **kwargs)


def hist(
    darray: DataArray,
    *args: Any,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    aspect: AspectOptions = None,
    ax: Axes | None = None,
    xincrease: bool | None = None,
    yincrease: bool | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    **kwargs: Any,
) -> tuple[np.ndarray, np.ndarray, BarContainer | Polygon]:
    """
    Histogram of DataArray.

    Wraps :py:func:`matplotlib:matplotlib.pyplot.hist`.

    Plots *N*-dimensional arrays by first flattening the array.

    Parameters
    ----------
    darray : DataArray
        Can have any number of dimensions.
    figsize : Iterable of float, optional
        A tuple (width, height) of the figure in inches.
        Mutually exclusive with ``size`` and ``ax``.
    aspect : "auto", "equal", scalar or None, optional
        Aspect ratio of plot, so that ``aspect * size`` gives the *width* in
        inches. Only used if a ``size`` is provided.
    size : scalar, optional
        If provided, create a new figure for the plot with the given size:
        *height* (in inches) of each plot. See also: ``aspect``.
    ax : matplotlib axes object, optional
        Axes on which to plot. By default, use the current axes.
        Mutually exclusive with ``size`` and ``figsize``.
    xincrease : bool or None, optional
        Should the values on the *x* axis be increasing from left to right?
        if ``None``, use the default for the Matplotlib function.
    yincrease : bool or None, optional
        Should the values on the *y* axis be increasing from top to bottom?
        if ``None``, use the default for the Matplotlib function.
    xscale, yscale : {'linear', 'symlog', 'log', 'logit'}, optional
        Specifies scaling for the *x*- and *y*-axis, respectively.
    xticks, yticks : array-like, optional
        Specify tick locations for *x*- and *y*-axis.
    xlim, ylim : tuple[float, float], optional
        Specify *x*- and *y*-axis limits.
    **kwargs : optional
        Additional keyword arguments to :py:func:`matplotlib:matplotlib.pyplot.hist`.

    """
    assert len(args) == 0

    if darray.ndim == 0 or darray.size == 0:
        # TypeError to be consistent with pandas
        raise TypeError("No numeric data to plot.")

    ax = get_axis(figsize, size, aspect, ax)

    no_nan_arr = np.ravel(darray.to_numpy())
    no_nan = no_nan_arr[pd.notnull(no_nan_arr)]

    n, bins, patches = cast(
        tuple[np.ndarray, np.ndarray, Union["BarContainer", "Polygon"]],
        ax.hist(no_nan, **kwargs),
    )

    ax.set_title(darray._title_for_slice())
    ax.set_xlabel(label_from_attrs(darray))

    _update_axes(ax, xincrease, yincrease, xscale, yscale, xticks, yticks, xlim, ylim)

    return n, bins, patches


def _plot1d(plotfunc):
    """Decorator for common 1d plotting logic."""
    commondoc = """
Parameters
----------
darray : DataArray
    Must be 2 dimensional, unless creating faceted plots.
x : Hashable or None, optional
    Coordinate for x axis. If None use darray.dims[1].
y : Hashable or None, optional
    Coordinate for y axis. If None use darray.dims[0].
z : Hashable or None, optional
    If specified plot 3D and use this coordinate for *z* axis.
hue : Hashable or None, optional
    Dimension or coordinate for which you want multiple lines plotted.
markersize: Hashable or None, optional
    scatter only. Variable by which to vary size of scattered points.
linewidth: Hashable or None, optional
    Variable by which to vary linewidth.
row : Hashable, optional
    If passed, make row faceted plots on this dimension name.
col : Hashable, optional
    If passed, make column faceted plots on this dimension name.
col_wrap : int, optional
    Use together with ``col`` to wrap faceted plots
ax : matplotlib axes object, optional
    If None, uses the current axis. Not applicable when using facets.
figsize : Iterable[float] or None, optional
    A tuple (width, height) of the figure in inches.
    Mutually exclusive with ``size`` and ``ax``.
size : scalar, optional
    If provided, create a new figure for the plot with the given size.
    Height (in inches) of each plot. See also: ``aspect``.
aspect : "auto", "equal", scalar or None, optional
    Aspect ratio of plot, so that ``aspect * size`` gives the width in
    inches. Only used if a ``size`` is provided.
xincrease : bool or None, default: True
    Should the values on the x axes be increasing from left to right?
    if None, use the default for the matplotlib function.
yincrease : bool or None, default: True
    Should the values on the y axes be increasing from top to bottom?
    if None, use the default for the matplotlib function.
add_legend : bool or None, optional
    If True use xarray metadata to add a legend.
add_colorbar : bool or None, optional
    If True add a colorbar.
add_labels : bool or None, optional
    If True use xarray metadata to label axes
add_title : bool or None, optional
    If True use xarray metadata to add a title
subplot_kws : dict, optional
    Dictionary of keyword arguments for matplotlib subplots. Only applies
    to FacetGrid plotting.
xscale : {'linear', 'symlog', 'log', 'logit'} or None, optional
    Specifies scaling for the x-axes.
yscale : {'linear', 'symlog', 'log', 'logit'} or None, optional
    Specifies scaling for the y-axes.
xticks : ArrayLike or None, optional
    Specify tick locations for x-axes.
yticks : ArrayLike or None, optional
    Specify tick locations for y-axes.
xlim : tuple[float, float] or None, optional
    Specify x-axes limits.
ylim : tuple[float, float] or None, optional
    Specify y-axes limits.
cmap : matplotlib colormap name or colormap, optional
    The mapping from data values to color space. Either a
    Matplotlib colormap name or object. If not provided, this will
    be either ``'viridis'`` (if the function infers a sequential
    dataset) or ``'RdBu_r'`` (if the function infers a diverging
    dataset).
    See :doc:`Choosing Colormaps in Matplotlib <matplotlib:users/explain/colors/colormaps>`
    for more information.

    If *seaborn* is installed, ``cmap`` may also be a
    `seaborn color palette <https://seaborn.pydata.org/tutorial/color_palettes.html>`_.
    Note: if ``cmap`` is a seaborn color palette,
    ``levels`` must also be specified.
vmin : float or None, optional
    Lower value to anchor the colormap, otherwise it is inferred from the
    data and other keyword arguments. When a diverging dataset is inferred,
    setting `vmin` or `vmax` will fix the other by symmetry around
    ``center``. Setting both values prevents use of a diverging colormap.
    If discrete levels are provided as an explicit list, both of these
    values are ignored.
vmax : float or None, optional
    Upper value to anchor the colormap, otherwise it is inferred from the
    data and other keyword arguments. When a diverging dataset is inferred,
    setting `vmin` or `vmax` will fix the other by symmetry around
    ``center``. Setting both values prevents use of a diverging colormap.
    If discrete levels are provided as an explicit list, both of these
    values are ignored.
norm : matplotlib.colors.Normalize, optional
    If ``norm`` has ``vmin`` or ``vmax`` specified, the corresponding
    kwarg must be ``None``.
extend : {'neither', 'both', 'min', 'max'}, optional
    How to draw arrows extending the colorbar beyond its limits. If not
    provided, ``extend`` is inferred from ``vmin``, ``vmax`` and the data limits.
levels : int or array-like, optional
    Split the colormap (``cmap``) into discrete color intervals. If an integer
    is provided, "nice" levels are chosen based on the data range: this can
    imply that the final number of levels is not exactly the expected one.
    Setting ``vmin`` and/or ``vmax`` with ``levels=N`` is equivalent to
    setting ``levels=np.linspace(vmin, vmax, N)``.
**kwargs : optional
    Additional arguments to wrapped matplotlib function

Returns
-------
artist :
    The same type of primitive artist that the wrapped matplotlib
    function returns
    """

    # Build on the original docstring
    plotfunc.__doc__ = f"{plotfunc.__doc__}\n{commondoc}"

    @functools.wraps(
        plotfunc, assigned=("__module__", "__name__", "__qualname__", "__doc__")
    )
    def newplotfunc(
        darray: DataArray,
        *args: Any,
        x: Hashable | None = None,
        y: Hashable | None = None,
        z: Hashable | None = None,
        hue: Hashable | None = None,
        hue_style: HueStyleOptions = None,
        markersize: Hashable | None = None,
        linewidth: Hashable | None = None,
        row: Hashable | None = None,
        col: Hashable | None = None,
        col_wrap: int | None = None,
        ax: Axes | None = None,
        figsize: Iterable[float] | None = None,
        size: float | None = None,
        aspect: float | None = None,
        xincrease: bool | None = True,
        yincrease: bool | None = True,
        add_legend: bool | None = None,
        add_colorbar: bool | None = None,
        add_labels: bool | Iterable[bool] = True,
        add_title: bool = True,
        subplot_kws: dict[str, Any] | None = None,
        xscale: ScaleOptions = None,
        yscale: ScaleOptions = None,
        xticks: ArrayLike | None = None,
        yticks: ArrayLike | None = None,
        xlim: tuple[float, float] | None = None,
        ylim: tuple[float, float] | None = None,
        cmap: str | Colormap | None = None,
        vmin: float | None = None,
        vmax: float | None = None,
        norm: Normalize | None = None,
        extend: ExtendOptions = None,
        levels: ArrayLike | None = None,
        **kwargs,
    ) -> Any:
        # All 1d plots in xarray share this function signature.
        # Method signature below should be consistent.

        if TYPE_CHECKING:
            import matplotlib.pyplot as plt
        else:
            plt = attempt_import("matplotlib.pyplot")

        if subplot_kws is None:
            subplot_kws = dict()

        # Handle facetgrids first
        if row or col:
            if z is not None:
                subplot_kws.update(projection="3d")

            allargs = locals().copy()
            allargs.update(allargs.pop("kwargs"))
            allargs.pop("darray")
            allargs.pop("plt")
            allargs["plotfunc"] = globals()[plotfunc.__name__]

            return _easy_facetgrid(darray, kind="plot1d", **allargs)

        if darray.ndim == 0 or darray.size == 0:
            # TypeError to be consistent with pandas
            raise TypeError("No numeric data to plot.")

        # The allargs dict passed to _easy_facetgrid above contains args
        if args == ():
            args = kwargs.pop("args", ())

        if args:
            assert "args" not in kwargs
            # TODO: Deprecated since 2022.10:
            msg = "Using positional arguments is deprecated for plot methods, use keyword arguments instead."
            assert x is None
            x = args[0]
            if len(args) > 1:
                assert y is None
                y = args[1]
            if len(args) > 2:
                assert z is None
                z = args[2]
            if len(args) > 3:
                assert hue is None
                hue = args[3]
            if len(args) > 4:
                raise ValueError(msg)
            else:
                warnings.warn(msg, DeprecationWarning, stacklevel=2)
        del args

        if hue_style is not None:
            # TODO: Not used since 2022.10. Deprecated since 2023.07.
            warnings.warn(
                (
                    "hue_style is no longer used for plot1d plots "
                    "and the argument will eventually be removed. "
                    "Convert numbers to string for a discrete hue "
                    "and use add_legend or add_colorbar to control which guide to display."
                ),
                DeprecationWarning,
                stacklevel=2,
            )

        _is_facetgrid = kwargs.pop("_is_facetgrid", False)

        if plotfunc.__name__ == "scatter":
            size_ = kwargs.pop("_size", markersize)
            size_r = _MARKERSIZE_RANGE

            # Remove any nulls, .where(m, drop=True) doesn't work when m is
            # a dask array, so load the array to memory.
            # It will have to be loaded to memory at some point anyway:
            darray = darray.compute()
            darray = darray.where(darray.notnull(), drop=True)
        else:
            size_ = kwargs.pop("_size", linewidth)
            size_r = _LINEWIDTH_RANGE

        # Get data to plot:
        coords_to_plot: MutableMapping[str, Hashable | None] = dict(
            x=x, z=z, hue=hue, size=size_
        )
        if not _is_facetgrid:
            # Guess what coords to use if some of the values in coords_to_plot are None:
            coords_to_plot = _guess_coords_to_plot(darray, coords_to_plot, kwargs)
        plts = _prepare_plot1d_data(darray, coords_to_plot, plotfunc.__name__)
        xplt = plts.pop("x", None)
        yplt = plts.pop("y", None)
        zplt = plts.pop("z", None)
        kwargs.update(zplt=zplt)
        hueplt = plts.pop("hue", None)
        sizeplt = plts.pop("size", None)

        # Handle size and hue:
        hueplt_norm = _Normalize(data=hueplt)
        kwargs.update(hueplt=hueplt_norm.values)
        sizeplt_norm = _Normalize(
            data=sizeplt, width=size_r, _is_facetgrid=_is_facetgrid
        )
        kwargs.update(sizeplt=sizeplt_norm.values)
        cmap_params_subset = kwargs.pop("cmap_params_subset", {})
        cbar_kwargs = kwargs.pop("cbar_kwargs", {})

        if hueplt_norm.data is not None:
            if not hueplt_norm.data_is_numeric:
                # Map hue values back to its original value:
                cbar_kwargs.update(format=hueplt_norm.format, ticks=hueplt_norm.ticks)
                levels = kwargs.get("levels", hueplt_norm.levels)

            cmap_params, cbar_kwargs = _process_cmap_cbar_kwargs(
                plotfunc,
                cast("DataArray", hueplt_norm.values).data,
                **locals(),
            )

            # subset that can be passed to scatter, hist2d
            if not cmap_params_subset:
                ckw = {vv: cmap_params[vv] for vv in ("vmin", "vmax", "norm", "cmap")}
                cmap_params_subset.update(**ckw)

        with plt.rc_context(_styles):
            if z is not None:
                import mpl_toolkits

                if ax is None:
                    subplot_kws.update(projection="3d")
                ax = get_axis(figsize, size, aspect, ax, **subplot_kws)
                assert isinstance(ax, mpl_toolkits.mplot3d.axes3d.Axes3D)

                # Using 30, 30 minimizes rotation of the plot. Making it easier to
                # build on your intuition from 2D plots:
                ax.view_init(azim=30, elev=30, vertical_axis="y")
            else:
                ax = get_axis(figsize, size, aspect, ax, **subplot_kws)

            primitive = plotfunc(
                xplt,
                yplt,
                ax=ax,
                add_labels=add_labels,
                **cmap_params_subset,
                **kwargs,
            )

        if np.any(np.asarray(add_labels)) and add_title:
            ax.set_title(darray._title_for_slice())

        add_colorbar_, add_legend_ = _determine_guide(
            hueplt_norm,
            sizeplt_norm,
            add_colorbar,
            add_legend,
            plotfunc_name=plotfunc.__name__,
        )

        if add_colorbar_:
            if "label" not in cbar_kwargs:
                cbar_kwargs["label"] = label_from_attrs(hueplt_norm.data)

            _add_colorbar(
                primitive, ax, kwargs.get("cbar_ax"), cbar_kwargs, cmap_params
            )

        if add_legend_:
            if plotfunc.__name__ in ["scatter", "line"]:
                _add_legend(
                    (
                        hueplt_norm
                        if add_legend or not add_colorbar_
                        else _Normalize(None)
                    ),
                    sizeplt_norm,
                    primitive,
                    legend_ax=ax,
                    plotfunc=plotfunc.__name__,
                )
            else:
                hueplt_norm_values: list[np.ndarray | None]
                if hueplt_norm.data is not None:
                    hueplt_norm_values = list(hueplt_norm.data.to_numpy())
                else:
                    hueplt_norm_values = [hueplt_norm.data]

                if plotfunc.__name__ == "hist":
                    ax.legend(
                        handles=primitive[-1],
                        labels=hueplt_norm_values,
                        title=label_from_attrs(hueplt_norm.data),
                    )
                else:
                    ax.legend(
                        handles=primitive,
                        labels=hueplt_norm_values,
                        title=label_from_attrs(hueplt_norm.data),
                    )

        _update_axes(
            ax, xincrease, yincrease, xscale, yscale, xticks, yticks, xlim, ylim
        )

        return primitive

    # we want to actually expose the signature of newplotfunc
    # and not the copied **kwargs from the plotfunc which
    # functools.wraps adds, so delete the wrapped attr
    del newplotfunc.__wrapped__

    return newplotfunc


def _add_labels(
    add_labels: bool | Iterable[bool],
    darrays: Iterable[DataArray | None],
    suffixes: Iterable[str],
    ax: Axes,
) -> None:
    """Set x, y, z labels."""
    add_labels = [add_labels] * 3 if isinstance(add_labels, bool) else add_labels
    axes: tuple[Literal["x", "y", "z"], ...] = ("x", "y", "z")
    for axis, add_label, darray, suffix in zip(
        axes, add_labels, darrays, suffixes, strict=True
    ):
        if darray is None:
            continue

        if add_label:
            label = label_from_attrs(darray, extra=suffix)
            if label is not None:
                getattr(ax, f"set_{axis}label")(label)

        if np.issubdtype(darray.dtype, np.datetime64):
            _set_concise_date(ax, axis=axis)


@overload
def scatter(  # type: ignore[misc,unused-ignore]  # None is hashable :(
    darray: DataArray,
    *args: Any,
    x: Hashable | None = None,
    y: Hashable | None = None,
    z: Hashable | None = None,
    hue: Hashable | None = None,
    hue_style: HueStyleOptions = None,
    markersize: Hashable | None = None,
    linewidth: Hashable | None = None,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    aspect: float | None = None,
    ax: Axes | None = None,
    row: None = None,  # no wrap -> primitive
    col: None = None,  # no wrap -> primitive
    col_wrap: int | None = None,
    xincrease: bool | None = True,
    yincrease: bool | None = True,
    add_legend: bool | None = None,
    add_colorbar: bool | None = None,
    add_labels: bool | Iterable[bool] = True,
    add_title: bool = True,
    subplot_kws: dict[str, Any] | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: ArrayLike | None = None,
    ylim: ArrayLike | None = None,
    cmap: str | Colormap | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
    norm: Normalize | None = None,
    extend: ExtendOptions = None,
    levels: ArrayLike | None = None,
    **kwargs,
) -> PathCollection: ...


@overload
def scatter(
    darray: T_DataArray,
    *args: Any,
    x: Hashable | None = None,
    y: Hashable | None = None,
    z: Hashable | None = None,
    hue: Hashable | None = None,
    hue_style: HueStyleOptions = None,
    markersize: Hashable | None = None,
    linewidth: Hashable | None = None,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    aspect: float | None = None,
    ax: Axes | None = None,
    row: Hashable | None = None,
    col: Hashable,  # wrap -> FacetGrid
    col_wrap: int | None = None,
    xincrease: bool | None = True,
    yincrease: bool | None = True,
    add_legend: bool | None = None,
    add_colorbar: bool | None = None,
    add_labels: bool | Iterable[bool] = True,
    add_title: bool = True,
    subplot_kws: dict[str, Any] | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: ArrayLike | None = None,
    ylim: ArrayLike | None = None,
    cmap: str | Colormap | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
    norm: Normalize | None = None,
    extend: ExtendOptions = None,
    levels: ArrayLike | None = None,
    **kwargs,
) -> FacetGrid[T_DataArray]: ...


@overload
def scatter(
    darray: T_DataArray,
    *args: Any,
    x: Hashable | None = None,
    y: Hashable | None = None,
    z: Hashable | None = None,
    hue: Hashable | None = None,
    hue_style: HueStyleOptions = None,
    markersize: Hashable | None = None,
    linewidth: Hashable | None = None,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    aspect: float | None = None,
    ax: Axes | None = None,
    row: Hashable,  # wrap -> FacetGrid
    col: Hashable | None = None,
    col_wrap: int | None = None,
    xincrease: bool | None = True,
    yincrease: bool | None = True,
    add_legend: bool | None = None,
    add_colorbar: bool | None = None,
    add_labels: bool | Iterable[bool] = True,
    add_title: bool = True,
    subplot_kws: dict[str, Any] | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: ArrayLike | None = None,
    ylim: ArrayLike | None = None,
    cmap: str | Colormap | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
    norm: Normalize | None = None,
    extend: ExtendOptions = None,
    levels: ArrayLike | None = None,
    **kwargs,
) -> FacetGrid[T_DataArray]: ...


@_plot1d
def scatter(
    xplt: DataArray | None,
    yplt: DataArray | None,
    ax: Axes,
    add_labels: bool | Iterable[bool] = True,
    **kwargs,
) -> PathCollection:
    """Scatter variables against each other.

    Wraps :py:func:`matplotlib:matplotlib.pyplot.scatter`.
    """
    if "u" in kwargs or "v" in kwargs:
        raise ValueError("u, v are not allowed in scatter plots.")

    zplt: DataArray | None = kwargs.pop("zplt", None)
    hueplt: DataArray | None = kwargs.pop("hueplt", None)
    sizeplt: DataArray | None = kwargs.pop("sizeplt", None)

    if hueplt is not None:
        kwargs.update(c=hueplt.to_numpy().ravel())

    if sizeplt is not None:
        kwargs.update(s=sizeplt.to_numpy().ravel())

    plts_or_none = (xplt, yplt, zplt)
    _add_labels(add_labels, plts_or_none, ("", "", ""), ax)

    xplt_np = None if xplt is None else xplt.to_numpy().ravel()
    yplt_np = None if yplt is None else yplt.to_numpy().ravel()
    zplt_np = None if zplt is None else zplt.to_numpy().ravel()
    plts_np = tuple(p for p in (xplt_np, yplt_np, zplt_np) if p is not None)

    if len(plts_np) == 3:
        import mpl_toolkits

        assert isinstance(ax, mpl_toolkits.mplot3d.axes3d.Axes3D)
        return ax.scatter(xplt_np, yplt_np, zplt_np, **kwargs)

    if len(plts_np) == 2:
        return ax.scatter(plts_np[0], plts_np[1], **kwargs)

    raise ValueError("At least two variables required for a scatter plot.")


def _plot2d(plotfunc):
    """Decorator for common 2d plotting logic."""
    commondoc = """
Parameters
----------
darray : DataArray
    Must be two-dimensional, unless creating faceted plots.
x : Hashable or None, optional
    Coordinate for *x* axis. If ``None``, use ``darray.dims[1]``.
y : Hashable or None, optional
    Coordinate for *y* axis. If ``None``, use ``darray.dims[0]``.
figsize : Iterable or float or None, optional
    A tuple (width, height) of the figure in inches.
    Mutually exclusive with ``size`` and ``ax``.
size : scalar, optional
    If provided, create a new figure for the plot with the given size:
    *height* (in inches) of each plot. See also: ``aspect``.
aspect : "auto", "equal", scalar or None, optional
    Aspect ratio of plot, so that ``aspect * size`` gives the *width* in
    inches. Only used if a ``size`` is provided.
ax : matplotlib axes object, optional
    Axes on which to plot. By default, use the current axes.
    Mutually exclusive with ``size`` and ``figsize``.
row : Hashable or None, optional
    If passed, make row faceted plots on this dimension name.
col : Hashable or None, optional
    If passed, make column faceted plots on this dimension name.
col_wrap : int, optional
    Use together with ``col`` to wrap faceted plots.
xincrease : None, True, or False, optional
    Should the values on the *x* axis be increasing from left to right?
    If ``None``, use the default for the Matplotlib function.
yincrease : None, True, or False, optional
    Should the values on the *y* axis be increasing from top to bottom?
    If ``None``, use the default for the Matplotlib function.
add_colorbar : bool, optional
    Add colorbar to axes.
add_labels : bool, optional
    Use xarray metadata to label axes.
vmin : float or None, optional
    Lower value to anchor the colormap, otherwise it is inferred from the
    data and other keyword arguments. When a diverging dataset is inferred,
    setting `vmin` or `vmax` will fix the other by symmetry around
    ``center``. Setting both values prevents use of a diverging colormap.
    If discrete levels are provided as an explicit list, both of these
    values are ignored.
vmax : float or None, optional
    Upper value to anchor the colormap, otherwise it is inferred from the
    data and other keyword arguments. When a diverging dataset is inferred,
    setting `vmin` or `vmax` will fix the other by symmetry around
    ``center``. Setting both values prevents use of a diverging colormap.
    If discrete levels are provided as an explicit list, both of these
    values are ignored.
cmap : matplotlib colormap name or colormap, optional
    The mapping from data values to color space. If not provided, this
    will be either be ``'viridis'`` (if the function infers a sequential
    dataset) or ``'RdBu_r'`` (if the function infers a diverging dataset).
    See :doc:`Choosing Colormaps in Matplotlib <matplotlib:users/explain/colors/colormaps>`
    for more information.
    If *seaborn* is installed, ``cmap`` may also be a
    `seaborn color palette <https://seaborn.pydata.org/tutorial/color_palettes.html>`_.
    Note: if ``cmap`` is a seaborn color palette and the plot type
    is not ``'contour'`` or ``'contourf'``, ``levels`` must also be specified.
center : float or False, optional
    The value at which to center the colormap. Passing this value implies
    use of a diverging colormap. Setting it to ``False`` prevents use of a
    diverging colormap.
robust : bool, optional
    If ``True`` and ``vmin`` or ``vmax`` are absent, the colormap range is
    computed with 2nd and 98th percentiles instead of the extreme values.
extend : {'neither', 'both', 'min', 'max'}, optional
    How to draw arrows extending the colorbar beyond its limits. If not
    provided, ``extend`` is inferred from ``vmin``, ``vmax`` and the data limits.
levels : int or array-like, optional
    Split the colormap (``cmap``) into discrete color intervals. If an integer
    is provided, "nice" levels are chosen based on the data range: this can
    imply that the final number of levels is not exactly the expected one.
    Setting ``vmin`` and/or ``vmax`` with ``levels=N`` is equivalent to
    setting ``levels=np.linspace(vmin, vmax, N)``.
infer_intervals : bool, optional
    Only applies to pcolormesh. If ``True``, the coordinate intervals are
    passed to pcolormesh. If ``False``, the original coordinates are used
    (this can be useful for certain map projections). The default is to
    always infer intervals, unless the mesh is irregular and plotted on
    a map projection.
colors : str or array-like of color-like, optional
    A single color or a sequence of colors. If the plot type is not ``'contour'``
    or ``'contourf'``, the ``levels`` argument is required.
subplot_kws : dict, optional
    Dictionary of keyword arguments for Matplotlib subplots. Only used
    for 2D and faceted plots.
    (see :py:meth:`matplotlib:matplotlib.figure.Figure.add_subplot`).
cbar_ax : matplotlib axes object, optional
    Axes in which to draw the colorbar.
cbar_kwargs : dict, optional
    Dictionary of keyword arguments to pass to the colorbar
    (see :meth:`matplotlib:matplotlib.figure.Figure.colorbar`).
xscale : {'linear', 'symlog', 'log', 'logit'} or None, optional
    Specifies scaling for the x-axes.
yscale : {'linear', 'symlog', 'log', 'logit'} or None, optional
    Specifies scaling for the y-axes.
xticks : ArrayLike or None, optional
    Specify tick locations for x-axes.
yticks : ArrayLike or None, optional
    Specify tick locations for y-axes.
xlim : tuple[float, float] or None, optional
    Specify x-axes limits.
ylim : tuple[float, float] or None, optional
    Specify y-axes limits.
norm : matplotlib.colors.Normalize, optional
    If ``norm`` has ``vmin`` or ``vmax`` specified, the corresponding
    kwarg must be ``None``.
**kwargs : optional
    Additional keyword arguments to wrapped Matplotlib function.

Returns
-------
artist :
    The same type of primitive artist that the wrapped Matplotlib
    function returns.
    """

    # Build on the original docstring
    plotfunc.__doc__ = f"{plotfunc.__doc__}\n{commondoc}"

    @functools.wraps(
        plotfunc, assigned=("__module__", "__name__", "__qualname__", "__doc__")
    )
    def newplotfunc(
        darray: DataArray,
        *args: Any,
        x: Hashable | None = None,
        y: Hashable | None = None,
        figsize: Iterable[float] | None = None,
        size: float | None = None,
        aspect: float | None = None,
        ax: Axes | None = None,
        row: Hashable | None = None,
        col: Hashable | None = None,
        col_wrap: int | None = None,
        xincrease: bool | None = True,
        yincrease: bool | None = True,
        add_colorbar: bool | None = None,
        add_labels: bool = True,
        vmin: float | None = None,
        vmax: float | None = None,
        cmap: str | Colormap | None = None,
        center: float | Literal[False] | None = None,
        robust: bool = False,
        extend: ExtendOptions = None,
        levels: ArrayLike | None = None,
        infer_intervals=None,
        colors: str | ArrayLike | None = None,
        subplot_kws: dict[str, Any] | None = None,
        cbar_ax: Axes | None = None,
        cbar_kwargs: dict[str, Any] | None = None,
        xscale: ScaleOptions = None,
        yscale: ScaleOptions = None,
        xticks: ArrayLike | None = None,
        yticks: ArrayLike | None = None,
        xlim: tuple[float, float] | None = None,
        ylim: tuple[float, float] | None = None,
        norm: Normalize | None = None,
        **kwargs: Any,
    ) -> Any:
        # All 2d plots in xarray share this function signature.

        if args:
            # TODO: Deprecated since 2022.10:
            msg = "Using positional arguments is deprecated for plot methods, use keyword arguments instead."
            assert x is None
            x = args[0]
            if len(args) > 1:
                assert y is None
                y = args[1]
            if len(args) > 2:
                raise ValueError(msg)
            else:
                warnings.warn(msg, DeprecationWarning, stacklevel=2)
        del args

        # Decide on a default for the colorbar before facetgrids
        if add_colorbar is None:
            add_colorbar = True
            if plotfunc.__name__ == "contour" or (
                plotfunc.__name__ == "surface" and cmap is None
            ):
                add_colorbar = False
        imshow_rgb = plotfunc.__name__ == "imshow" and darray.ndim == (
            3 + (row is not None) + (col is not None)
        )
        if imshow_rgb:
            # Don't add a colorbar when showing an image with explicit colors
            add_colorbar = False
            # Matplotlib does not support normalising RGB data, so do it here.
            # See eg. https://github.com/matplotlib/matplotlib/pull/10220
            if robust or vmax is not None or vmin is not None:
                darray = _rescale_imshow_rgb(darray.as_numpy(), vmin, vmax, robust)
                vmin, vmax, robust = None, None, False

        if subplot_kws is None:
            subplot_kws = dict()

        if plotfunc.__name__ == "surface" and not kwargs.get("_is_facetgrid"):
            if ax is None:
                # TODO: Importing Axes3D is no longer necessary in matplotlib >= 3.2.
                # Remove when minimum requirement of matplotlib is 3.2:
                from mpl_toolkits.mplot3d import Axes3D

                # delete so it does not end up in locals()
                del Axes3D

                # Need to create a "3d" Axes instance for surface plots
                subplot_kws["projection"] = "3d"

            # In facet grids, shared axis labels don't make sense for surface plots
            sharex = False
            sharey = False

        # Handle facetgrids first
        if row or col:
            allargs = locals().copy()
            del allargs["darray"]
            del allargs["imshow_rgb"]
            allargs.update(allargs.pop("kwargs"))
            # Need the decorated plotting function
            allargs["plotfunc"] = globals()[plotfunc.__name__]
            return _easy_facetgrid(darray, kind="dataarray", **allargs)

        if darray.ndim == 0 or darray.size == 0:
            # TypeError to be consistent with pandas
            raise TypeError("No numeric data to plot.")

        if (
            plotfunc.__name__ == "surface"
            and not kwargs.get("_is_facetgrid")
            and ax is not None
        ):
            import mpl_toolkits

            if not isinstance(ax, mpl_toolkits.mplot3d.Axes3D):
                raise ValueError(
                    "If ax is passed to surface(), it must be created with "
                    'projection="3d"'
                )

        rgb = kwargs.pop("rgb", None)
        if rgb is not None and plotfunc.__name__ != "imshow":
            raise ValueError('The "rgb" keyword is only valid for imshow()')
        elif rgb is not None and not imshow_rgb:
            raise ValueError(
                'The "rgb" keyword is only valid for imshow()'
                "with a three-dimensional array (per facet)"
            )

        xlab, ylab = _infer_xy_labels(
            darray=darray, x=x, y=y, imshow=imshow_rgb, rgb=rgb
        )

        xval = darray[xlab]
        yval = darray[ylab]

        if xval.ndim > 1 or yval.ndim > 1 or plotfunc.__name__ == "surface":
            # Passing 2d coordinate values, need to ensure they are transposed the same
            # way as darray.
            # Also surface plots always need 2d coordinates
            xval = xval.broadcast_like(darray)
            yval = yval.broadcast_like(darray)
            dims = darray.dims
        else:
            dims = (yval.dims[0], xval.dims[0])

        # May need to transpose for correct x, y labels
        # xlab may be the name of a coord, we have to check for dim names
        if imshow_rgb:
            # For RGB[A] images, matplotlib requires the color dimension
            # to be last.  In Xarray the order should be unimportant, so
            # we transpose to (y, x, color) to make this work.
            yx_dims = (ylab, xlab)
            dims = yx_dims + tuple(d for d in darray.dims if d not in yx_dims)

        if dims != darray.dims:
            darray = darray.transpose(*dims, transpose_coords=True)

        # better to pass the ndarrays directly to plotting functions
        xvalnp = xval.to_numpy()
        yvalnp = yval.to_numpy()

        # Pass the data as a masked ndarray too
        zval = darray.to_masked_array(copy=False)

        # Replace pd.Intervals if contained in xval or yval.
        xplt, xlab_extra = _resolve_intervals_2dplot(xvalnp, plotfunc.__name__)
        yplt, ylab_extra = _resolve_intervals_2dplot(yvalnp, plotfunc.__name__)

        _ensure_plottable(xplt, yplt, zval)

        cmap_params, cbar_kwargs = _process_cmap_cbar_kwargs(
            plotfunc,
            zval.data,
            **locals(),
            _is_facetgrid=kwargs.pop("_is_facetgrid", False),
        )

        if "contour" in plotfunc.__name__:
            # extend is a keyword argument only for contour and contourf, but
            # passing it to the colorbar is sufficient for imshow and
            # pcolormesh
            kwargs["extend"] = cmap_params["extend"]
            kwargs["levels"] = cmap_params["levels"]
            # if colors == a single color, matplotlib draws dashed negative
            # contours. we lose this feature if we pass cmap and not colors
            if colors is not None:
                cmap_params["cmap"] = None
                kwargs["colors"] = colors

        if "pcolormesh" == plotfunc.__name__:
            kwargs["infer_intervals"] = infer_intervals
            kwargs["xscale"] = xscale
            kwargs["yscale"] = yscale

        if "imshow" == plotfunc.__name__ and isinstance(aspect, str):
            # forbid usage of mpl strings
            raise ValueError("plt.imshow's `aspect` kwarg is not available in xarray")

        ax = get_axis(figsize, size, aspect, ax, **subplot_kws)

        primitive = plotfunc(
            xplt,
            yplt,
            zval,
            ax=ax,
            cmap=cmap_params["cmap"],
            vmin=cmap_params["vmin"],
            vmax=cmap_params["vmax"],
            norm=cmap_params["norm"],
            **kwargs,
        )

        # Label the plot with metadata
        if add_labels:
            ax.set_xlabel(label_from_attrs(darray[xlab], xlab_extra))
            ax.set_ylabel(label_from_attrs(darray[ylab], ylab_extra))
            ax.set_title(darray._title_for_slice())
            if plotfunc.__name__ == "surface":
                import mpl_toolkits

                assert isinstance(ax, mpl_toolkits.mplot3d.axes3d.Axes3D)
                ax.set_zlabel(label_from_attrs(darray))

        if add_colorbar:
            if add_labels and "label" not in cbar_kwargs:
                cbar_kwargs["label"] = label_from_attrs(darray)
            cbar = _add_colorbar(primitive, ax, cbar_ax, cbar_kwargs, cmap_params)
        elif cbar_ax is not None or cbar_kwargs:
            # inform the user about keywords which aren't used
            raise ValueError(
                "cbar_ax and cbar_kwargs can't be used with add_colorbar=False."
            )

        # origin kwarg overrides yincrease
        if "origin" in kwargs:
            yincrease = None

        _update_axes(
            ax, xincrease, yincrease, xscale, yscale, xticks, yticks, xlim, ylim
        )

        if np.issubdtype(xplt.dtype, np.datetime64):
            _set_concise_date(ax, "x")

        return primitive

    # we want to actually expose the signature of newplotfunc
    # and not the copied **kwargs from the plotfunc which
    # functools.wraps adds, so delete the wrapped attr
    del newplotfunc.__wrapped__

    return newplotfunc


@overload
def imshow(  # type: ignore[misc,unused-ignore]  # None is hashable :(
    darray: DataArray,
    x: Hashable | None = None,
    y: Hashable | None = None,
    *,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    aspect: float | None = None,
    ax: Axes | None = None,
    row: None = None,  # no wrap -> primitive
    col: None = None,  # no wrap -> primitive
    col_wrap: int | None = None,
    xincrease: bool | None = True,
    yincrease: bool | None = True,
    add_colorbar: bool | None = None,
    add_labels: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str | Colormap | None = None,
    center: float | Literal[False] | None = None,
    robust: bool = False,
    extend: ExtendOptions = None,
    levels: ArrayLike | None = None,
    infer_intervals=None,
    colors: str | ArrayLike | None = None,
    subplot_kws: dict[str, Any] | None = None,
    cbar_ax: Axes | None = None,
    cbar_kwargs: dict[str, Any] | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: ArrayLike | None = None,
    ylim: ArrayLike | None = None,
    norm: Normalize | None = None,
    **kwargs: Any,
) -> AxesImage: ...


@overload
def imshow(
    darray: T_DataArray,
    x: Hashable | None = None,
    y: Hashable | None = None,
    *,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    aspect: AspectOptions = None,
    ax: Axes | None = None,
    row: Hashable | None = None,
    col: Hashable,  # wrap -> FacetGrid
    col_wrap: int | None = None,
    xincrease: bool | None = True,
    yincrease: bool | None = True,
    add_colorbar: bool | None = None,
    add_labels: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str | Colormap | None = None,
    center: float | Literal[False] | None = None,
    robust: bool = False,
    extend: ExtendOptions = None,
    levels: ArrayLike | None = None,
    infer_intervals=None,
    colors: str | ArrayLike | None = None,
    subplot_kws: dict[str, Any] | None = None,
    cbar_ax: Axes | None = None,
    cbar_kwargs: dict[str, Any] | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: ArrayLike | None = None,
    ylim: ArrayLike | None = None,
    norm: Normalize | None = None,
    **kwargs: Any,
) -> FacetGrid[T_DataArray]: ...


@overload
def imshow(
    darray: T_DataArray,
    x: Hashable | None = None,
    y: Hashable | None = None,
    *,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    aspect: AspectOptions = None,
    ax: Axes | None = None,
    row: Hashable,  # wrap -> FacetGrid
    col: Hashable | None = None,
    col_wrap: int | None = None,
    xincrease: bool | None = True,
    yincrease: bool | None = True,
    add_colorbar: bool | None = None,
    add_labels: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str | Colormap | None = None,
    center: float | Literal[False] | None = None,
    robust: bool = False,
    extend: ExtendOptions = None,
    levels: ArrayLike | None = None,
    infer_intervals=None,
    colors: str | ArrayLike | None = None,
    subplot_kws: dict[str, Any] | None = None,
    cbar_ax: Axes | None = None,
    cbar_kwargs: dict[str, Any] | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: ArrayLike | None = None,
    ylim: ArrayLike | None = None,
    norm: Normalize | None = None,
    **kwargs: Any,
) -> FacetGrid[T_DataArray]: ...


@_plot2d
def imshow(
    x: np.ndarray, y: np.ndarray, z: np.ma.core.MaskedArray, ax: Axes, **kwargs: Any
) -> AxesImage:
    """
    Image plot of 2D DataArray.

    Wraps :py:func:`matplotlib:matplotlib.pyplot.imshow`.

    While other plot methods require the DataArray to be strictly
    two-dimensional, ``imshow`` also accepts a 3D array where some
    dimension can be interpreted as RGB or RGBA color channels and
    allows this dimension to be specified via the kwarg ``rgb=``.

    Unlike :py:func:`matplotlib:matplotlib.pyplot.imshow`, which ignores ``vmin``/``vmax``
    for RGB(A) data,
    xarray *will* use ``vmin`` and ``vmax`` for RGB(A) data
    by applying a single scaling factor and offset to all bands.
    Passing  ``robust=True`` infers ``vmin`` and ``vmax``
    :ref:`in the usual way <robust-plotting>`.
    Additionally the y-axis is not inverted by default, you can
    restore the matplotlib behavior by setting `yincrease=False`.

    .. note::
        This function needs uniformly spaced coordinates to
        properly label the axes. Call :py:meth:`DataArray.plot` to check.

    The pixels are centered on the coordinates. For example, if the coordinate
    value is 3.2, then the pixels for those coordinates will be centered on 3.2.
    """

    if x.ndim != 1 or y.ndim != 1:
        raise ValueError(
            "imshow requires 1D coordinates, try using pcolormesh or contour(f)"
        )

    def _center_pixels(x):
        """Center the pixels on the coordinates."""
        if np.issubdtype(x.dtype, str):
            # When using strings as inputs imshow converts it to
            # integers. Choose extent values which puts the indices in
            # in the center of the pixels:
            return 0 - 0.5, len(x) - 0.5

        try:
            # Center the pixels assuming uniform spacing:
            xstep = 0.5 * (x[1] - x[0])
        except IndexError:
            # Arbitrary default value, similar to matplotlib behaviour:
            xstep = 0.1

        return x[0] - xstep, x[-1] + xstep

    # Center the pixels:
    left, right = _center_pixels(x)
    top, bottom = _center_pixels(y)

    defaults: dict[str, Any] = {"origin": "upper", "interpolation": "nearest"}

    if not hasattr(ax, "projection"):
        # not for cartopy geoaxes
        defaults["aspect"] = "auto"

    # Allow user to override these defaults
    defaults.update(kwargs)

    if defaults["origin"] == "upper":
        defaults["extent"] = [left, right, bottom, top]
    else:
        defaults["extent"] = [left, right, top, bottom]

    if z.ndim == 3:
        # matplotlib imshow uses black for missing data, but Xarray makes
        # missing data transparent.  We therefore add an alpha channel if
        # there isn't one, and set it to transparent where data is masked.
        if z.shape[-1] == 3:
            safe_dtype = np.promote_types(z.dtype, np.uint8)
            alpha = np.ma.ones(z.shape[:2] + (1,), dtype=safe_dtype)
            if np.issubdtype(z.dtype, np.integer):
                alpha[:] = 255
            z = np.ma.concatenate((z, alpha), axis=2)
        else:
            z = z.copy()
        z[np.any(z.mask, axis=-1), -1] = 0

    primitive = ax.imshow(z, **defaults)

    # If x or y are strings the ticklabels have been replaced with
    # integer indices. Replace them back to strings:
    for axis, v in [("x", x), ("y", y)]:
        if np.issubdtype(v.dtype, str):
            getattr(ax, f"set_{axis}ticks")(np.arange(len(v)))
            getattr(ax, f"set_{axis}ticklabels")(v)

    return primitive


@overload
def contour(  # type: ignore[misc,unused-ignore]  # None is hashable :(
    darray: DataArray,
    x: Hashable | None = None,
    y: Hashable | None = None,
    *,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    aspect: float | None = None,
    ax: Axes | None = None,
    row: None = None,  # no wrap -> primitive
    col: None = None,  # no wrap -> primitive
    col_wrap: int | None = None,
    xincrease: bool | None = True,
    yincrease: bool | None = True,
    add_colorbar: bool | None = None,
    add_labels: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str | Colormap | None = None,
    center: float | Literal[False] | None = None,
    robust: bool = False,
    extend: ExtendOptions = None,
    levels: ArrayLike | None = None,
    infer_intervals=None,
    colors: str | ArrayLike | None = None,
    subplot_kws: dict[str, Any] | None = None,
    cbar_ax: Axes | None = None,
    cbar_kwargs: dict[str, Any] | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: ArrayLike | None = None,
    ylim: ArrayLike | None = None,
    norm: Normalize | None = None,
    **kwargs: Any,
) -> QuadContourSet: ...


@overload
def contour(
    darray: T_DataArray,
    x: Hashable | None = None,
    y: Hashable | None = None,
    *,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    aspect: AspectOptions = None,
    ax: Axes | None = None,
    row: Hashable | None = None,
    col: Hashable,  # wrap -> FacetGrid
    col_wrap: int | None = None,
    xincrease: bool | None = True,
    yincrease: bool | None = True,
    add_colorbar: bool | None = None,
    add_labels: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str | Colormap | None = None,
    center: float | Literal[False] | None = None,
    robust: bool = False,
    extend: ExtendOptions = None,
    levels: ArrayLike | None = None,
    infer_intervals=None,
    colors: str | ArrayLike | None = None,
    subplot_kws: dict[str, Any] | None = None,
    cbar_ax: Axes | None = None,
    cbar_kwargs: dict[str, Any] | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: ArrayLike | None = None,
    ylim: ArrayLike | None = None,
    norm: Normalize | None = None,
    **kwargs: Any,
) -> FacetGrid[T_DataArray]: ...


@overload
def contour(
    darray: T_DataArray,
    x: Hashable | None = None,
    y: Hashable | None = None,
    *,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    aspect: AspectOptions = None,
    ax: Axes | None = None,
    row: Hashable,  # wrap -> FacetGrid
    col: Hashable | None = None,
    col_wrap: int | None = None,
    xincrease: bool | None = True,
    yincrease: bool | None = True,
    add_colorbar: bool | None = None,
    add_labels: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str | Colormap | None = None,
    center: float | Literal[False] | None = None,
    robust: bool = False,
    extend: ExtendOptions = None,
    levels: ArrayLike | None = None,
    infer_intervals=None,
    colors: str | ArrayLike | None = None,
    subplot_kws: dict[str, Any] | None = None,
    cbar_ax: Axes | None = None,
    cbar_kwargs: dict[str, Any] | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: ArrayLike | None = None,
    ylim: ArrayLike | None = None,
    norm: Normalize | None = None,
    **kwargs: Any,
) -> FacetGrid[T_DataArray]: ...


@_plot2d
def contour(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, ax: Axes, **kwargs: Any
) -> QuadContourSet:
    """
    Contour plot of 2D DataArray.

    Wraps :py:func:`matplotlib:matplotlib.pyplot.contour`.
    """
    primitive = ax.contour(x, y, z, **kwargs)
    return primitive


@overload
def contourf(  # type: ignore[misc,unused-ignore]  # None is hashable :(
    darray: DataArray,
    x: Hashable | None = None,
    y: Hashable | None = None,
    *,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    aspect: AspectOptions = None,
    ax: Axes | None = None,
    row: None = None,  # no wrap -> primitive
    col: None = None,  # no wrap -> primitive
    col_wrap: int | None = None,
    xincrease: bool | None = True,
    yincrease: bool | None = True,
    add_colorbar: bool | None = None,
    add_labels: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str | Colormap | None = None,
    center: float | Literal[False] | None = None,
    robust: bool = False,
    extend: ExtendOptions = None,
    levels: ArrayLike | None = None,
    infer_intervals=None,
    colors: str | ArrayLike | None = None,
    subplot_kws: dict[str, Any] | None = None,
    cbar_ax: Axes | None = None,
    cbar_kwargs: dict[str, Any] | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: ArrayLike | None = None,
    ylim: ArrayLike | None = None,
    norm: Normalize | None = None,
    **kwargs: Any,
) -> QuadContourSet: ...


@overload
def contourf(
    darray: T_DataArray,
    x: Hashable | None = None,
    y: Hashable | None = None,
    *,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    aspect: AspectOptions = None,
    ax: Axes | None = None,
    row: Hashable | None = None,
    col: Hashable,  # wrap -> FacetGrid
    col_wrap: int | None = None,
    xincrease: bool | None = True,
    yincrease: bool | None = True,
    add_colorbar: bool | None = None,
    add_labels: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str | Colormap | None = None,
    center: float | Literal[False] | None = None,
    robust: bool = False,
    extend: ExtendOptions = None,
    levels: ArrayLike | None = None,
    infer_intervals=None,
    colors: str | ArrayLike | None = None,
    subplot_kws: dict[str, Any] | None = None,
    cbar_ax: Axes | None = None,
    cbar_kwargs: dict[str, Any] | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: ArrayLike | None = None,
    ylim: ArrayLike | None = None,
    norm: Normalize | None = None,
    **kwargs: Any,
) -> FacetGrid[T_DataArray]: ...


@overload
def contourf(
    darray: T_DataArray,
    x: Hashable | None = None,
    y: Hashable | None = None,
    *,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    aspect: AspectOptions = None,
    ax: Axes | None = None,
    row: Hashable,  # wrap -> FacetGrid
    col: Hashable | None = None,
    col_wrap: int | None = None,
    xincrease: bool | None = True,
    yincrease: bool | None = True,
    add_colorbar: bool | None = None,
    add_labels: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str | Colormap | None = None,
    center: float | Literal[False] | None = None,
    robust: bool = False,
    extend: ExtendOptions = None,
    levels: ArrayLike | None = None,
    infer_intervals=None,
    colors: str | ArrayLike | None = None,
    subplot_kws: dict[str, Any] | None = None,
    cbar_ax: Axes | None = None,
    cbar_kwargs: dict[str, Any] | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: ArrayLike | None = None,
    ylim: ArrayLike | None = None,
    norm: Normalize | None = None,
    **kwargs: Any,
) -> FacetGrid[T_DataArray]: ...


@_plot2d
def contourf(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, ax: Axes, **kwargs: Any
) -> QuadContourSet:
    """
    Filled contour plot of 2D DataArray.

    Wraps :py:func:`matplotlib:matplotlib.pyplot.contourf`.
    """
    primitive = ax.contourf(x, y, z, **kwargs)
    return primitive


@overload
def pcolormesh(  # type: ignore[misc,unused-ignore]  # None is hashable :(
    darray: DataArray,
    x: Hashable | None = None,
    y: Hashable | None = None,
    *,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    aspect: AspectOptions = None,
    ax: Axes | None = None,
    row: None = None,  # no wrap -> primitive
    col: None = None,  # no wrap -> primitive
    col_wrap: int | None = None,
    xincrease: bool | None = True,
    yincrease: bool | None = True,
    add_colorbar: bool | None = None,
    add_labels: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str | Colormap | None = None,
    center: float | Literal[False] | None = None,
    robust: bool = False,
    extend: ExtendOptions = None,
    levels: ArrayLike | None = None,
    infer_intervals=None,
    colors: str | ArrayLike | None = None,
    subplot_kws: dict[str, Any] | None = None,
    cbar_ax: Axes | None = None,
    cbar_kwargs: dict[str, Any] | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: ArrayLike | None = None,
    ylim: ArrayLike | None = None,
    norm: Normalize | None = None,
    **kwargs: Any,
) -> QuadMesh: ...


@overload
def pcolormesh(
    darray: T_DataArray,
    x: Hashable | None = None,
    y: Hashable | None = None,
    *,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    aspect: AspectOptions = None,
    ax: Axes | None = None,
    row: Hashable | None = None,
    col: Hashable,  # wrap -> FacetGrid
    col_wrap: int | None = None,
    xincrease: bool | None = True,
    yincrease: bool | None = True,
    add_colorbar: bool | None = None,
    add_labels: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str | Colormap | None = None,
    center: float | Literal[False] | None = None,
    robust: bool = False,
    extend: ExtendOptions = None,
    levels: ArrayLike | None = None,
    infer_intervals=None,
    colors: str | ArrayLike | None = None,
    subplot_kws: dict[str, Any] | None = None,
    cbar_ax: Axes | None = None,
    cbar_kwargs: dict[str, Any] | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: ArrayLike | None = None,
    ylim: ArrayLike | None = None,
    norm: Normalize | None = None,
    **kwargs: Any,
) -> FacetGrid[T_DataArray]: ...


@overload
def pcolormesh(
    darray: T_DataArray,
    x: Hashable | None = None,
    y: Hashable | None = None,
    *,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    aspect: AspectOptions = None,
    ax: Axes | None = None,
    row: Hashable,  # wrap -> FacetGrid
    col: Hashable | None = None,
    col_wrap: int | None = None,
    xincrease: bool | None = True,
    yincrease: bool | None = True,
    add_colorbar: bool | None = None,
    add_labels: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str | Colormap | None = None,
    center: float | Literal[False] | None = None,
    robust: bool = False,
    extend: ExtendOptions = None,
    levels: ArrayLike | None = None,
    infer_intervals=None,
    colors: str | ArrayLike | None = None,
    subplot_kws: dict[str, Any] | None = None,
    cbar_ax: Axes | None = None,
    cbar_kwargs: dict[str, Any] | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: ArrayLike | None = None,
    ylim: ArrayLike | None = None,
    norm: Normalize | None = None,
    **kwargs: Any,
) -> FacetGrid[T_DataArray]: ...


@_plot2d
def pcolormesh(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    ax: Axes,
    xscale: ScaleOptions | None = None,
    yscale: ScaleOptions | None = None,
    infer_intervals=None,
    **kwargs: Any,
) -> QuadMesh:
    """
    Pseudocolor plot of 2D DataArray.

    Wraps :py:func:`matplotlib:matplotlib.pyplot.pcolormesh`.
    """

    # decide on a default for infer_intervals (GH781)
    x = np.asarray(x)
    if infer_intervals is None:
        if hasattr(ax, "projection"):
            if len(x.shape) == 1:
                infer_intervals = True
            else:
                infer_intervals = False
        else:
            infer_intervals = True

    if any(np.issubdtype(k.dtype, str) for k in (x, y)):
        # do not infer intervals if any axis contains str ticks, see #6775
        infer_intervals = False

    if infer_intervals and (
        (np.shape(x)[0] == np.shape(z)[1])
        or ((x.ndim > 1) and (np.shape(x)[1] == np.shape(z)[1]))
    ):
        if x.ndim == 1:
            x = _infer_interval_breaks(x, check_monotonic=True, scale=xscale)
        else:
            # we have to infer the intervals on both axes
            x = _infer_interval_breaks(x, axis=1, scale=xscale)
            x = _infer_interval_breaks(x, axis=0, scale=xscale)

    if infer_intervals and (np.shape(y)[0] == np.shape(z)[0]):
        if y.ndim == 1:
            y = _infer_interval_breaks(y, check_monotonic=True, scale=yscale)
        else:
            # we have to infer the intervals on both axes
            y = _infer_interval_breaks(y, axis=1, scale=yscale)
            y = _infer_interval_breaks(y, axis=0, scale=yscale)

    ax.grid(False)
    primitive = ax.pcolormesh(x, y, z, **kwargs)

    # by default, pcolormesh picks "round" values for bounds
    # this results in ugly looking plots with lots of surrounding whitespace
    if not hasattr(ax, "projection") and x.ndim == 1 and y.ndim == 1:
        # not a cartopy geoaxis
        ax.set_xlim(x[0], x[-1])
        ax.set_ylim(y[0], y[-1])

    return primitive


@overload
def surface(
    darray: DataArray,
    x: Hashable | None = None,
    y: Hashable | None = None,
    *,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    aspect: AspectOptions = None,
    ax: Axes | None = None,
    row: None = None,  # no wrap -> primitive
    col: None = None,  # no wrap -> primitive
    col_wrap: int | None = None,
    xincrease: bool | None = True,
    yincrease: bool | None = True,
    add_colorbar: bool | None = None,
    add_labels: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str | Colormap | None = None,
    center: float | Literal[False] | None = None,
    robust: bool = False,
    extend: ExtendOptions = None,
    levels: ArrayLike | None = None,
    infer_intervals=None,
    colors: str | ArrayLike | None = None,
    subplot_kws: dict[str, Any] | None = None,
    cbar_ax: Axes | None = None,
    cbar_kwargs: dict[str, Any] | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: ArrayLike | None = None,
    ylim: ArrayLike | None = None,
    norm: Normalize | None = None,
    **kwargs: Any,
) -> Poly3DCollection: ...


@overload
def surface(
    darray: T_DataArray,
    x: Hashable | None = None,
    y: Hashable | None = None,
    *,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    aspect: AspectOptions = None,
    ax: Axes | None = None,
    row: Hashable | None = None,
    col: Hashable,  # wrap -> FacetGrid
    col_wrap: int | None = None,
    xincrease: bool | None = True,
    yincrease: bool | None = True,
    add_colorbar: bool | None = None,
    add_labels: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str | Colormap | None = None,
    center: float | Literal[False] | None = None,
    robust: bool = False,
    extend: ExtendOptions = None,
    levels: ArrayLike | None = None,
    infer_intervals=None,
    colors: str | ArrayLike | None = None,
    subplot_kws: dict[str, Any] | None = None,
    cbar_ax: Axes | None = None,
    cbar_kwargs: dict[str, Any] | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: ArrayLike | None = None,
    ylim: ArrayLike | None = None,
    norm: Normalize | None = None,
    **kwargs: Any,
) -> FacetGrid[T_DataArray]: ...


@overload
def surface(
    darray: T_DataArray,
    x: Hashable | None = None,
    y: Hashable | None = None,
    *,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    aspect: AspectOptions = None,
    ax: Axes | None = None,
    row: Hashable,  # wrap -> FacetGrid
    col: Hashable | None = None,
    col_wrap: int | None = None,
    xincrease: bool | None = True,
    yincrease: bool | None = True,
    add_colorbar: bool | None = None,
    add_labels: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str | Colormap | None = None,
    center: float | Literal[False] | None = None,
    robust: bool = False,
    extend: ExtendOptions = None,
    levels: ArrayLike | None = None,
    infer_intervals=None,
    colors: str | ArrayLike | None = None,
    subplot_kws: dict[str, Any] | None = None,
    cbar_ax: Axes | None = None,
    cbar_kwargs: dict[str, Any] | None = None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: ArrayLike | None = None,
    ylim: ArrayLike | None = None,
    norm: Normalize | None = None,
    **kwargs: Any,
) -> FacetGrid[T_DataArray]: ...


@_plot2d
def surface(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, ax: Axes, **kwargs: Any
) -> Poly3DCollection:
    """
    Surface plot of 2D DataArray.

    Wraps :py:meth:`matplotlib:mpl_toolkits.mplot3d.axes3d.Axes3D.plot_surface`.
    """
    import mpl_toolkits

    assert isinstance(ax, mpl_toolkits.mplot3d.axes3d.Axes3D)
    primitive = ax.plot_surface(x, y, z, **kwargs)
    return primitive
