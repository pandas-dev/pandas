from __future__ import annotations

import functools
import inspect
import warnings
from collections.abc import Callable, Hashable, Iterable
from typing import TYPE_CHECKING, Any, TypeVar, overload

from xarray.plot import dataarray_plot
from xarray.plot.facetgrid import _easy_facetgrid
from xarray.plot.utils import (
    _add_colorbar,
    _get_nice_quiver_magnitude,
    _infer_meta_data,
    _process_cmap_cbar_kwargs,
    get_axis,
)
from xarray.structure.alignment import broadcast

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.collections import LineCollection, PathCollection
    from matplotlib.colors import Colormap, Normalize
    from matplotlib.quiver import Quiver
    from numpy.typing import ArrayLike

    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset
    from xarray.core.types import (
        AspectOptions,
        ExtendOptions,
        HueStyleOptions,
        ScaleOptions,
    )
    from xarray.plot.facetgrid import FacetGrid


def _dsplot(plotfunc):
    commondoc = """
Parameters
----------

ds : Dataset
x : Hashable or None, optional
    Variable name for x-axis.
y : Hashable or None, optional
    Variable name for y-axis.
u : Hashable or None, optional
    Variable name for the *u* velocity (in *x* direction).
    quiver/streamplot plots only.
v : Hashable or None, optional
    Variable name for the *v* velocity (in *y* direction).
    quiver/streamplot plots only.
hue: Hashable or None, optional
    Variable by which to color scatter points or arrows.
hue_style: {'continuous', 'discrete'} or None, optional
    How to use the ``hue`` variable:

    - ``'continuous'`` -- continuous color scale
        (default for numeric ``hue`` variables)
    - ``'discrete'`` -- a color for each unique value, using the default color cycle
        (default for non-numeric ``hue`` variables)

row : Hashable or None, optional
    If passed, make row faceted plots on this dimension name.
col : Hashable or None, optional
    If passed, make column faceted plots on this dimension name.
col_wrap : int, optional
    Use together with ``col`` to wrap faceted plots.
ax : matplotlib axes object or None, optional
    If ``None``, use the current axes. Not applicable when using facets.
figsize : Iterable[float] or None, optional
    A tuple (width, height) of the figure in inches.
    Mutually exclusive with ``size`` and ``ax``.
size : scalar, optional
    If provided, create a new figure for the plot with the given size.
    Height (in inches) of each plot. See also: ``aspect``.
aspect : "auto", "equal", scalar or None, optional
    Aspect ratio of plot, so that ``aspect * size`` gives the width in
    inches. Only used if a ``size`` is provided.
sharex : bool or None, optional
    If True all subplots share the same x-axis.
sharey : bool or None, optional
    If True all subplots share the same y-axis.
add_guide: bool or None, optional
    Add a guide that depends on ``hue_style``:

    - ``'continuous'`` -- build a colorbar
    - ``'discrete'`` -- build a legend

subplot_kws : dict or None, optional
    Dictionary of keyword arguments for Matplotlib subplots
    (see :py:meth:`matplotlib:matplotlib.figure.Figure.add_subplot`).
    Only applies to FacetGrid plotting.
cbar_kwargs : dict, optional
    Dictionary of keyword arguments to pass to the colorbar
    (see :meth:`matplotlib:matplotlib.figure.Figure.colorbar`).
cbar_ax : matplotlib axes object, optional
    Axes in which to draw the colorbar.
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
infer_intervals: bool | None
    If True the intervals are inferred.
center : float, optional
    The value at which to center the colormap. Passing this value implies
    use of a diverging colormap. Setting it to ``False`` prevents use of a
    diverging colormap.
robust : bool, optional
    If ``True`` and ``vmin`` or ``vmax`` are absent, the colormap range is
    computed with 2nd and 98th percentiles instead of the extreme values.
colors : str or array-like of color-like, optional
    A single color or a list of colors. The ``levels`` argument
    is required.
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
    Additional keyword arguments to wrapped Matplotlib function.
    """

    # Build on the original docstring
    plotfunc.__doc__ = f"{plotfunc.__doc__}\n{commondoc}"

    @functools.wraps(
        plotfunc, assigned=("__module__", "__name__", "__qualname__", "__doc__")
    )
    def newplotfunc(
        ds: Dataset,
        *args: Any,
        x: Hashable | None = None,
        y: Hashable | None = None,
        u: Hashable | None = None,
        v: Hashable | None = None,
        hue: Hashable | None = None,
        hue_style: HueStyleOptions = None,
        row: Hashable | None = None,
        col: Hashable | None = None,
        col_wrap: int | None = None,
        ax: Axes | None = None,
        figsize: Iterable[float] | None = None,
        size: float | None = None,
        aspect: AspectOptions = None,
        sharex: bool = True,
        sharey: bool = True,
        add_guide: bool | None = None,
        subplot_kws: dict[str, Any] | None = None,
        cbar_kwargs: dict[str, Any] | None = None,
        cbar_ax: Axes | None = None,
        cmap: str | Colormap | None = None,
        vmin: float | None = None,
        vmax: float | None = None,
        norm: Normalize | None = None,
        infer_intervals: bool | None = None,
        center: float | None = None,
        robust: bool | None = None,
        colors: str | ArrayLike | None = None,
        extend: ExtendOptions = None,
        levels: ArrayLike | None = None,
        **kwargs: Any,
    ) -> Any:
        if args:
            # TODO: Deprecated since 2022.10:
            msg = "Using positional arguments is deprecated for plot methods, use keyword arguments instead."
            assert x is None
            x = args[0]
            if len(args) > 1:
                assert y is None
                y = args[1]
            if len(args) > 2:
                assert u is None
                u = args[2]
            if len(args) > 3:
                assert v is None
                v = args[3]
            if len(args) > 4:
                assert hue is None
                hue = args[4]
            if len(args) > 5:
                raise ValueError(msg)
            else:
                warnings.warn(msg, DeprecationWarning, stacklevel=2)
        del args

        _is_facetgrid = kwargs.pop("_is_facetgrid", False)
        if _is_facetgrid:  # facetgrid call
            meta_data = kwargs.pop("meta_data")
        else:
            meta_data = _infer_meta_data(
                ds, x, y, hue, hue_style, add_guide, funcname=plotfunc.__name__
            )

        hue_style = meta_data["hue_style"]

        # handle facetgrids first
        if col or row:
            allargs = locals().copy()
            allargs["plotfunc"] = globals()[plotfunc.__name__]
            allargs["data"] = ds
            # remove kwargs to avoid passing the information twice
            for arg in ["meta_data", "kwargs", "ds"]:
                del allargs[arg]

            return _easy_facetgrid(kind="dataset", **allargs, **kwargs)

        figsize = kwargs.pop("figsize", None)
        ax = get_axis(figsize, size, aspect, ax)

        if hue_style == "continuous" and hue is not None:
            if _is_facetgrid:
                cbar_kwargs = meta_data["cbar_kwargs"]
                cmap_params = meta_data["cmap_params"]
            else:
                cmap_params, cbar_kwargs = _process_cmap_cbar_kwargs(
                    plotfunc, ds[hue].values, **locals()
                )

            # subset that can be passed to scatter, hist2d
            cmap_params_subset = {
                vv: cmap_params[vv] for vv in ["vmin", "vmax", "norm", "cmap"]
            }

        else:
            cmap_params_subset = {}

        if (u is not None or v is not None) and plotfunc.__name__ not in (
            "quiver",
            "streamplot",
        ):
            raise ValueError("u, v are only allowed for quiver or streamplot plots.")

        primitive = plotfunc(
            ds=ds,
            x=x,
            y=y,
            ax=ax,
            u=u,
            v=v,
            hue=hue,
            hue_style=hue_style,
            cmap_params=cmap_params_subset,
            **kwargs,
        )

        if _is_facetgrid:  # if this was called from Facetgrid.map_dataset,
            return primitive  # finish here. Else, make labels

        if meta_data.get("xlabel", None):
            ax.set_xlabel(meta_data.get("xlabel"))
        if meta_data.get("ylabel", None):
            ax.set_ylabel(meta_data.get("ylabel"))

        if meta_data["add_legend"]:
            ax.legend(handles=primitive, title=meta_data.get("hue_label", None))
        if meta_data["add_colorbar"]:
            cbar_kwargs = {} if cbar_kwargs is None else cbar_kwargs
            if "label" not in cbar_kwargs:
                cbar_kwargs["label"] = meta_data.get("hue_label", None)
            _add_colorbar(primitive, ax, cbar_ax, cbar_kwargs, cmap_params)

        if meta_data["add_quiverkey"]:
            magnitude = _get_nice_quiver_magnitude(ds[u], ds[v])
            units = ds[u].attrs.get("units", "")
            ax.quiverkey(
                primitive,
                X=0.85,
                Y=0.9,
                U=magnitude,
                label=f"{magnitude}\n{units}",
                labelpos="E",
                coordinates="figure",
            )

        if plotfunc.__name__ in ("quiver", "streamplot"):
            title = ds[u]._title_for_slice()
        else:
            title = ds[x]._title_for_slice()
        ax.set_title(title)

        return primitive

    # we want to actually expose the signature of newplotfunc
    # and not the copied **kwargs from the plotfunc which
    # functools.wraps adds, so delete the wrapped attr
    del newplotfunc.__wrapped__

    return newplotfunc


@overload
def quiver(  # type: ignore[misc,unused-ignore]  # None is hashable :(
    ds: Dataset,
    *args: Any,
    x: Hashable | None = None,
    y: Hashable | None = None,
    u: Hashable | None = None,
    v: Hashable | None = None,
    hue: Hashable | None = None,
    hue_style: HueStyleOptions = None,
    col: None = None,  # no wrap -> primitive
    row: None = None,  # no wrap -> primitive
    ax: Axes | None = None,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    col_wrap: int | None = None,
    sharex: bool = True,
    sharey: bool = True,
    aspect: AspectOptions = None,
    subplot_kws: dict[str, Any] | None = None,
    add_guide: bool | None = None,
    cbar_kwargs: dict[str, Any] | None = None,
    cbar_ax: Axes | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
    norm: Normalize | None = None,
    infer_intervals: bool | None = None,
    center: float | None = None,
    levels: ArrayLike | None = None,
    robust: bool | None = None,
    colors: str | ArrayLike | None = None,
    extend: ExtendOptions = None,
    cmap: str | Colormap | None = None,
    **kwargs: Any,
) -> Quiver: ...


@overload
def quiver(
    ds: Dataset,
    *args: Any,
    x: Hashable | None = None,
    y: Hashable | None = None,
    u: Hashable | None = None,
    v: Hashable | None = None,
    hue: Hashable | None = None,
    hue_style: HueStyleOptions = None,
    col: Hashable,  # wrap -> FacetGrid
    row: Hashable | None = None,
    ax: Axes | None = None,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    col_wrap: int | None = None,
    sharex: bool = True,
    sharey: bool = True,
    aspect: AspectOptions = None,
    subplot_kws: dict[str, Any] | None = None,
    add_guide: bool | None = None,
    cbar_kwargs: dict[str, Any] | None = None,
    cbar_ax: Axes | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
    norm: Normalize | None = None,
    infer_intervals: bool | None = None,
    center: float | None = None,
    levels: ArrayLike | None = None,
    robust: bool | None = None,
    colors: str | ArrayLike | None = None,
    extend: ExtendOptions = None,
    cmap: str | Colormap | None = None,
    **kwargs: Any,
) -> FacetGrid[Dataset]: ...


@overload
def quiver(
    ds: Dataset,
    *args: Any,
    x: Hashable | None = None,
    y: Hashable | None = None,
    u: Hashable | None = None,
    v: Hashable | None = None,
    hue: Hashable | None = None,
    hue_style: HueStyleOptions = None,
    col: Hashable | None = None,
    row: Hashable,  # wrap -> FacetGrid
    ax: Axes | None = None,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    col_wrap: int | None = None,
    sharex: bool = True,
    sharey: bool = True,
    aspect: AspectOptions = None,
    subplot_kws: dict[str, Any] | None = None,
    add_guide: bool | None = None,
    cbar_kwargs: dict[str, Any] | None = None,
    cbar_ax: Axes | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
    norm: Normalize | None = None,
    infer_intervals: bool | None = None,
    center: float | None = None,
    levels: ArrayLike | None = None,
    robust: bool | None = None,
    colors: str | ArrayLike | None = None,
    extend: ExtendOptions = None,
    cmap: str | Colormap | None = None,
    **kwargs: Any,
) -> FacetGrid[Dataset]: ...


@_dsplot
def quiver(
    ds: Dataset,
    x: Hashable,
    y: Hashable,
    ax: Axes,
    u: Hashable,
    v: Hashable,
    **kwargs: Any,
) -> Quiver:
    """Quiver plot of Dataset variables.

    Wraps :py:func:`matplotlib:matplotlib.pyplot.quiver`.
    """
    import matplotlib as mpl

    if x is None or y is None or u is None or v is None:
        raise ValueError("Must specify x, y, u, v for quiver plots.")

    dx, dy, du, dv = broadcast(ds[x], ds[y], ds[u], ds[v])

    args = [dx.values, dy.values, du.values, dv.values]
    hue = kwargs.pop("hue")
    cmap_params = kwargs.pop("cmap_params")

    if hue:
        args.append(ds[hue].values)

        # TODO: Fix this by always returning a norm with vmin, vmax in cmap_params
        if not cmap_params["norm"]:
            cmap_params["norm"] = mpl.colors.Normalize(
                cmap_params.pop("vmin"), cmap_params.pop("vmax")
            )

    kwargs.pop("hue_style")
    kwargs.setdefault("pivot", "middle")
    hdl = ax.quiver(*args, **kwargs, **cmap_params)
    return hdl


@overload
def streamplot(  # type: ignore[misc,unused-ignore]  # None is hashable :(
    ds: Dataset,
    *args: Any,
    x: Hashable | None = None,
    y: Hashable | None = None,
    u: Hashable | None = None,
    v: Hashable | None = None,
    hue: Hashable | None = None,
    hue_style: HueStyleOptions = None,
    col: None = None,  # no wrap -> primitive
    row: None = None,  # no wrap -> primitive
    ax: Axes | None = None,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    col_wrap: int | None = None,
    sharex: bool = True,
    sharey: bool = True,
    aspect: AspectOptions = None,
    subplot_kws: dict[str, Any] | None = None,
    add_guide: bool | None = None,
    cbar_kwargs: dict[str, Any] | None = None,
    cbar_ax: Axes | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
    norm: Normalize | None = None,
    infer_intervals: bool | None = None,
    center: float | None = None,
    levels: ArrayLike | None = None,
    robust: bool | None = None,
    colors: str | ArrayLike | None = None,
    extend: ExtendOptions = None,
    cmap: str | Colormap | None = None,
    **kwargs: Any,
) -> LineCollection: ...


@overload
def streamplot(
    ds: Dataset,
    *args: Any,
    x: Hashable | None = None,
    y: Hashable | None = None,
    u: Hashable | None = None,
    v: Hashable | None = None,
    hue: Hashable | None = None,
    hue_style: HueStyleOptions = None,
    col: Hashable,  # wrap -> FacetGrid
    row: Hashable | None = None,
    ax: Axes | None = None,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    col_wrap: int | None = None,
    sharex: bool = True,
    sharey: bool = True,
    aspect: AspectOptions = None,
    subplot_kws: dict[str, Any] | None = None,
    add_guide: bool | None = None,
    cbar_kwargs: dict[str, Any] | None = None,
    cbar_ax: Axes | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
    norm: Normalize | None = None,
    infer_intervals: bool | None = None,
    center: float | None = None,
    levels: ArrayLike | None = None,
    robust: bool | None = None,
    colors: str | ArrayLike | None = None,
    extend: ExtendOptions = None,
    cmap: str | Colormap | None = None,
    **kwargs: Any,
) -> FacetGrid[Dataset]: ...


@overload
def streamplot(
    ds: Dataset,
    *args: Any,
    x: Hashable | None = None,
    y: Hashable | None = None,
    u: Hashable | None = None,
    v: Hashable | None = None,
    hue: Hashable | None = None,
    hue_style: HueStyleOptions = None,
    col: Hashable | None = None,
    row: Hashable,  # wrap -> FacetGrid
    ax: Axes | None = None,
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    col_wrap: int | None = None,
    sharex: bool = True,
    sharey: bool = True,
    aspect: AspectOptions = None,
    subplot_kws: dict[str, Any] | None = None,
    add_guide: bool | None = None,
    cbar_kwargs: dict[str, Any] | None = None,
    cbar_ax: Axes | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
    norm: Normalize | None = None,
    infer_intervals: bool | None = None,
    center: float | None = None,
    levels: ArrayLike | None = None,
    robust: bool | None = None,
    colors: str | ArrayLike | None = None,
    extend: ExtendOptions = None,
    cmap: str | Colormap | None = None,
    **kwargs: Any,
) -> FacetGrid[Dataset]: ...


@_dsplot
def streamplot(
    ds: Dataset,
    x: Hashable,
    y: Hashable,
    ax: Axes,
    u: Hashable,
    v: Hashable,
    **kwargs: Any,
) -> LineCollection:
    """Plot streamlines of Dataset variables.

    Wraps :py:func:`matplotlib:matplotlib.pyplot.streamplot`.
    """
    import matplotlib as mpl

    if x is None or y is None or u is None or v is None:
        raise ValueError("Must specify x, y, u, v for streamplot plots.")

    # Matplotlib's streamplot has strong restrictions on what x and y can be, so need to
    # get arrays transposed the 'right' way around. 'x' cannot vary within 'rows', so
    # the dimension of x must be the second dimension. 'y' cannot vary with 'columns' so
    # the dimension of y must be the first dimension. If x and y are both 2d, assume the
    # user has got them right already.
    xdim = ds[x].dims[0] if len(ds[x].dims) == 1 else None
    ydim = ds[y].dims[0] if len(ds[y].dims) == 1 else None
    if xdim is not None and ydim is None:
        ydims = set(ds[y].dims) - {xdim}
        if len(ydims) == 1:
            ydim = next(iter(ydims))
    if ydim is not None and xdim is None:
        xdims = set(ds[x].dims) - {ydim}
        if len(xdims) == 1:
            xdim = next(iter(xdims))

    dx, dy, du, dv = broadcast(ds[x], ds[y], ds[u], ds[v])

    if xdim is not None and ydim is not None:
        # Need to ensure the arrays are transposed correctly
        dx = dx.transpose(ydim, xdim)
        dy = dy.transpose(ydim, xdim)
        du = du.transpose(ydim, xdim)
        dv = dv.transpose(ydim, xdim)

    hue = kwargs.pop("hue")
    cmap_params = kwargs.pop("cmap_params")

    if hue:
        if xdim is not None and ydim is not None:
            ds[hue] = ds[hue].transpose(ydim, xdim)
        kwargs["color"] = ds[hue].values

        # TODO: Fix this by always returning a norm with vmin, vmax in cmap_params
        if not cmap_params["norm"]:
            cmap_params["norm"] = mpl.colors.Normalize(
                cmap_params.pop("vmin"), cmap_params.pop("vmax")
            )

    kwargs.pop("hue_style")
    hdl = ax.streamplot(
        dx.values, dy.values, du.values, dv.values, **kwargs, **cmap_params
    )

    # Return .lines so colorbar creation works properly
    return hdl.lines


F = TypeVar("F", bound=Callable)


def _update_doc_to_dataset(dataarray_plotfunc: Callable) -> Callable[[F], F]:
    """
    Add a common docstring by reusing the DataArray one.

    TODO: Reduce code duplication.

    * The goal is to reduce code duplication by moving all Dataset
      specific plots to the DataArray side and use this thin wrapper to
      handle the conversion between Dataset and DataArray.
    * Improve docstring handling, maybe reword the DataArray versions to
      explain Datasets better.

    Parameters
    ----------
    dataarray_plotfunc : Callable
        Function that returns a finished plot primitive.
    """

    # Build on the original docstring
    da_doc = dataarray_plotfunc.__doc__
    if da_doc is None:
        raise NotImplementedError("DataArray plot method requires a docstring")

    da_str = """
    Parameters
    ----------
    darray : DataArray
    """
    ds_str = """

    The `y` DataArray will be used as base, any other variables are added as coords.

    Parameters
    ----------
    ds : Dataset
    """
    # TODO: improve this?
    if da_str in da_doc:
        ds_doc = da_doc.replace(da_str, ds_str).replace("darray", "ds")
    else:
        ds_doc = da_doc

    @functools.wraps(dataarray_plotfunc)
    def wrapper(dataset_plotfunc: F) -> F:
        dataset_plotfunc.__doc__ = ds_doc
        return dataset_plotfunc

    return wrapper  # type: ignore[return-value]


def _normalize_args(
    plotmethod: str, args: tuple[Any, ...], kwargs: dict[str, Any]
) -> dict[str, Any]:
    from xarray.core.dataarray import DataArray

    # Determine positional arguments keyword by inspecting the
    # signature of the plotmethod:
    locals_ = dict(
        inspect.signature(getattr(DataArray().plot, plotmethod))
        .bind(*args, **kwargs)
        .arguments.items()
    )
    locals_.update(locals_.pop("kwargs", {}))

    return locals_


def _temp_dataarray(ds: Dataset, y: Hashable, locals_: dict[str, Any]) -> DataArray:
    """Create a temporary datarray with extra coords."""
    from xarray.core.dataarray import DataArray

    coords = dict(ds[y].coords)
    dims = set(ds[y].dims)

    # Add extra coords to the DataArray from valid kwargs, if using all
    # kwargs there is a risk that we add unnecessary dataarrays as
    # coords straining RAM further for example:
    # ds.both and extend="both" would add ds.both to the coords:
    valid_coord_kwargs = {"x", "z", "markersize", "hue", "row", "col", "u", "v"}
    coord_kwargs = locals_.keys() & valid_coord_kwargs
    for k in coord_kwargs:
        key = locals_[k]
        darray = ds.get(key)
        if darray is not None:
            coords[key] = darray
            dims.update(darray.dims)

    # Trim dataset from unnecessary dims:
    ds_trimmed = ds.drop_dims(ds.sizes.keys() - dims)  # TODO: Use ds.dims in the future

    # The dataarray has to include all the dims. Broadcast to that shape
    # and add the additional coords:
    _y = ds[y].broadcast_like(ds_trimmed)

    return DataArray(_y, coords=coords)


@overload
def scatter(  # type: ignore[misc,unused-ignore]  # None is hashable :(
    ds: Dataset,
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
    **kwargs: Any,
) -> PathCollection: ...


@overload
def scatter(
    ds: Dataset,
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
    **kwargs: Any,
) -> FacetGrid[DataArray]: ...


@overload
def scatter(
    ds: Dataset,
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
    **kwargs: Any,
) -> FacetGrid[DataArray]: ...


@_update_doc_to_dataset(dataarray_plot.scatter)
def scatter(
    ds: Dataset,
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
    **kwargs: Any,
) -> PathCollection | FacetGrid[DataArray]:
    """Scatter plot Dataset data variables against each other."""
    locals_ = locals()
    del locals_["ds"]
    locals_.update(locals_.pop("kwargs", {}))
    da = _temp_dataarray(ds, y, locals_)

    return da.plot.scatter(*locals_.pop("args", ()), **locals_)
