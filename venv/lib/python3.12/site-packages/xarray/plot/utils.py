from __future__ import annotations

import itertools
import textwrap
import warnings
from collections.abc import (
    Callable,
    Hashable,
    Iterable,
    Mapping,
    MutableMapping,
    Sequence,
)
from datetime import date, datetime
from inspect import getfullargspec
from typing import TYPE_CHECKING, Any, Literal, cast, overload

import numpy as np
import pandas as pd

from xarray.core.indexes import PandasMultiIndex
from xarray.core.options import OPTIONS
from xarray.core.utils import (
    attempt_import,
    is_scalar,
    module_available,
)
from xarray.namedarray.pycompat import DuckArrayModule

nc_time_axis_available = module_available("nc_time_axis")


try:
    import cftime
except ImportError:
    cftime = None


if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.colors import Normalize
    from matplotlib.ticker import FuncFormatter
    from numpy.typing import ArrayLike

    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset
    from xarray.core.types import AspectOptions, ScaleOptions

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        plt: Any = None  # type: ignore[no-redef]

ROBUST_PERCENTILE = 2.0

# copied from seaborn
_MARKERSIZE_RANGE = (18.0, 36.0, 72.0)
_LINEWIDTH_RANGE = (1.5, 1.5, 6.0)


def _determine_extend(calc_data, vmin, vmax):
    extend_min = calc_data.min() < vmin
    extend_max = calc_data.max() > vmax
    if extend_min and extend_max:
        return "both"
    elif extend_min:
        return "min"
    elif extend_max:
        return "max"
    else:
        return "neither"


def _build_discrete_cmap(cmap, levels, extend, filled):
    """
    Build a discrete colormap and normalization of the data.
    """
    import matplotlib as mpl

    if len(levels) == 1:
        levels = [levels[0], levels[0]]

    if not filled:
        # non-filled contour plots
        extend = "max"

    if extend == "both":
        ext_n = 2
    elif extend in ["min", "max"]:
        ext_n = 1
    else:
        ext_n = 0

    n_colors = len(levels) + ext_n - 1
    pal = _color_palette(cmap, n_colors)

    new_cmap, cnorm = mpl.colors.from_levels_and_colors(levels, pal, extend=extend)
    # copy the old cmap name, for easier testing
    new_cmap.name = getattr(cmap, "name", cmap)

    # copy colors to use for bad, under, and over values in case they have been
    # set to non-default values
    try:
        # matplotlib<3.2 only uses bad color for masked values
        bad = cmap(np.ma.masked_invalid([np.nan]))[0]
    except TypeError:
        # cmap was a str or list rather than a color-map object, so there are
        # no bad, under or over values to check or copy
        pass
    else:
        under = cmap(-np.inf)
        over = cmap(np.inf)

        new_cmap.set_bad(bad)

        # Only update under and over if they were explicitly changed by the user
        # (i.e. are different from the lowest or highest values in cmap). Otherwise
        # leave unchanged so new_cmap uses its default values (its own lowest and
        # highest values).
        if under != cmap(0):
            new_cmap.set_under(under)
        if over != cmap(cmap.N - 1):
            new_cmap.set_over(over)

    return new_cmap, cnorm


def _color_palette(cmap, n_colors):
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    colors_i = np.linspace(0, 1.0, n_colors)
    if isinstance(cmap, list | tuple):
        # expand or truncate the list of colors to n_colors
        cmap = list(itertools.islice(itertools.cycle(cmap), n_colors))
        cmap = ListedColormap(cmap)
        pal = cmap(colors_i)
    elif isinstance(cmap, str):
        # we have some sort of named palette
        try:
            # is this a matplotlib cmap?
            cmap = plt.get_cmap(cmap)
            pal = cmap(colors_i)
        except ValueError:
            # ValueError happens when mpl doesn't like a colormap, try seaborn
            try:
                from seaborn import color_palette

                pal = color_palette(cmap, n_colors=n_colors)
            except (ValueError, ImportError):
                # or maybe we just got a single color as a string
                cmap = ListedColormap([cmap] * n_colors)
                pal = cmap(colors_i)
    else:
        # cmap better be a LinearSegmentedColormap (e.g. viridis)
        pal = cmap(colors_i)

    return pal


# _determine_cmap_params is adapted from Seaborn:
# https://github.com/mwaskom/seaborn/blob/v0.6/seaborn/matrix.py#L158
# Used under the terms of Seaborn's license, see licenses/SEABORN_LICENSE.


def _determine_cmap_params(
    plot_data,
    vmin=None,
    vmax=None,
    cmap=None,
    center=None,
    robust=False,
    extend=None,
    levels=None,
    filled=True,
    norm=None,
    _is_facetgrid=False,
):
    """
    Use some heuristics to set good defaults for colorbar and range.

    Parameters
    ----------
    plot_data : Numpy array
        Doesn't handle xarray objects

    Returns
    -------
    cmap_params : dict
        Use depends on the type of the plotting function
    """
    if TYPE_CHECKING:
        import matplotlib as mpl
    else:
        mpl = attempt_import("matplotlib")

    if isinstance(levels, Iterable):
        levels = sorted(levels)

    calc_data = np.ravel(plot_data[np.isfinite(plot_data)])

    # Handle all-NaN input data gracefully
    if calc_data.size == 0:
        # Arbitrary default for when all values are NaN
        calc_data = np.array(0.0)

    # Setting center=False prevents a divergent cmap
    possibly_divergent = center is not False

    # Set center to 0 so math below makes sense but remember its state
    center_is_none = False
    if center is None:
        center = 0
        center_is_none = True

    # Setting both vmin and vmax prevents a divergent cmap
    if (vmin is not None) and (vmax is not None):
        possibly_divergent = False

    # Setting vmin or vmax implies linspaced levels
    user_minmax = (vmin is not None) or (vmax is not None)

    # vlim might be computed below
    vlim = None

    # save state; needed later
    vmin_was_none = vmin is None
    vmax_was_none = vmax is None

    if vmin is None:
        if robust:
            vmin = np.percentile(calc_data, ROBUST_PERCENTILE)
        else:
            vmin = calc_data.min()
    elif possibly_divergent:
        vlim = abs(vmin - center)

    if vmax is None:
        if robust:
            vmax = np.percentile(calc_data, 100 - ROBUST_PERCENTILE)
        else:
            vmax = calc_data.max()
    elif possibly_divergent:
        vlim = abs(vmax - center)

    if possibly_divergent:
        levels_are_divergent = (
            isinstance(levels, Iterable) and levels[0] * levels[-1] < 0
        )
        # kwargs not specific about divergent or not: infer defaults from data
        divergent = (vmin < 0 < vmax) or not center_is_none or levels_are_divergent
    else:
        divergent = False

    # A divergent map should be symmetric around the center value
    if divergent:
        if vlim is None:
            vlim = max(abs(vmin - center), abs(vmax - center))
        vmin, vmax = -vlim, vlim

    # Now add in the centering value and set the limits
    vmin += center
    vmax += center

    # now check norm and harmonize with vmin, vmax
    if norm is not None:
        if norm.vmin is None:
            norm.vmin = vmin
        else:
            if not vmin_was_none and vmin != norm.vmin:
                raise ValueError("Cannot supply vmin and a norm with a different vmin.")
            vmin = norm.vmin

        if norm.vmax is None:
            norm.vmax = vmax
        else:
            if not vmax_was_none and vmax != norm.vmax:
                raise ValueError("Cannot supply vmax and a norm with a different vmax.")
            vmax = norm.vmax

    # if BoundaryNorm, then set levels
    if isinstance(norm, mpl.colors.BoundaryNorm):
        levels = norm.boundaries

    # Choose default colormaps if not provided
    if cmap is None:
        if divergent:
            cmap = OPTIONS["cmap_divergent"]
        else:
            cmap = OPTIONS["cmap_sequential"]

    # Handle discrete levels
    if levels is not None:
        if is_scalar(levels):
            if user_minmax:
                levels = np.linspace(vmin, vmax, levels)
            elif levels == 1:
                levels = np.asarray([(vmin + vmax) / 2])
            else:
                # N in MaxNLocator refers to bins, not ticks
                ticker = mpl.ticker.MaxNLocator(levels - 1)
                levels = ticker.tick_values(vmin, vmax)
        vmin, vmax = levels[0], levels[-1]

    # GH3734
    if vmin == vmax:
        vmin, vmax = mpl.ticker.LinearLocator(2).tick_values(vmin, vmax)

    if extend is None:
        extend = _determine_extend(calc_data, vmin, vmax)

    if (levels is not None) and (not isinstance(norm, mpl.colors.BoundaryNorm)):
        cmap, newnorm = _build_discrete_cmap(cmap, levels, extend, filled)
        norm = newnorm if norm is None else norm

    # vmin & vmax needs to be None if norm is passed
    # TODO: always return a norm with vmin and vmax
    if norm is not None:
        vmin = None
        vmax = None

    return dict(
        vmin=vmin, vmax=vmax, cmap=cmap, extend=extend, levels=levels, norm=norm
    )


def _infer_xy_labels_3d(
    darray: DataArray | Dataset,
    x: Hashable | None,
    y: Hashable | None,
    rgb: Hashable | None,
) -> tuple[Hashable, Hashable]:
    """
    Determine x and y labels for showing RGB images.

    Attempts to infer which dimension is RGB/RGBA by size and order of dims.

    """
    assert rgb is None or rgb != x
    assert rgb is None or rgb != y
    # Start by detecting and reporting invalid combinations of arguments
    assert darray.ndim == 3
    not_none = [a for a in (x, y, rgb) if a is not None]
    if len(set(not_none)) < len(not_none):
        raise ValueError(
            "Dimension names must be None or unique strings, but imshow was "
            f"passed x={x!r}, y={y!r}, and rgb={rgb!r}."
        )
    for label in not_none:
        if label not in darray.dims:
            raise ValueError(f"{label!r} is not a dimension")

    # Then calculate rgb dimension if certain and check validity
    could_be_color = [
        label
        for label in darray.dims
        if darray[label].size in (3, 4) and label not in (x, y)
    ]
    if rgb is None and not could_be_color:
        raise ValueError(
            "A 3-dimensional array was passed to imshow(), but there is no "
            "dimension that could be color.  At least one dimension must be "
            "of size 3 (RGB) or 4 (RGBA), and not given as x or y."
        )
    if rgb is None and len(could_be_color) == 1:
        rgb = could_be_color[0]
    if rgb is not None and darray[rgb].size not in (3, 4):
        raise ValueError(
            f"Cannot interpret dim {rgb!r} of size {darray[rgb].size} as RGB or RGBA."
        )

    # If rgb dimension is still unknown, there must be two or three dimensions
    # in could_be_color.  We therefore warn, and use a heuristic to break ties.
    if rgb is None:
        assert len(could_be_color) in (2, 3)
        rgb = could_be_color[-1]
        warnings.warn(
            "Several dimensions of this array could be colors.  Xarray "
            f"will use the last possible dimension ({rgb!r}) to match "
            "matplotlib.pyplot.imshow.  You can pass names of x, y, "
            "and/or rgb dimensions to override this guess.",
            stacklevel=2,
        )
    assert rgb is not None

    # Finally, we pick out the red slice and delegate to the 2D version:
    return _infer_xy_labels(darray.isel({rgb: 0}), x, y)


def _infer_xy_labels(
    darray: DataArray | Dataset,
    x: Hashable | None,
    y: Hashable | None,
    imshow: bool = False,
    rgb: Hashable | None = None,
) -> tuple[Hashable, Hashable]:
    """
    Determine x and y labels. For use in _plot2d

    darray must be a 2 dimensional data array, or 3d for imshow only.
    """
    if (x is not None) and (x == y):
        raise ValueError("x and y cannot be equal.")

    if imshow and darray.ndim == 3:
        return _infer_xy_labels_3d(darray, x, y, rgb)

    if x is None and y is None:
        if darray.ndim != 2:
            raise ValueError("DataArray must be 2d")
        y, x = darray.dims
    elif x is None:
        _assert_valid_xy(darray, y, "y")
        x = darray.dims[0] if y == darray.dims[1] else darray.dims[1]
    elif y is None:
        _assert_valid_xy(darray, x, "x")
        y = darray.dims[0] if x == darray.dims[1] else darray.dims[1]
    else:
        _assert_valid_xy(darray, x, "x")
        _assert_valid_xy(darray, y, "y")

        if darray._indexes.get(x, 1) is darray._indexes.get(y, 2) and isinstance(
            darray._indexes[x], PandasMultiIndex
        ):
            raise ValueError("x and y cannot be levels of the same MultiIndex")

    return x, y


# TODO: Can by used to more than x or y, rename?
def _assert_valid_xy(
    darray: DataArray | Dataset, xy: Hashable | None, name: str
) -> None:
    """
    make sure x and y passed to plotting functions are valid
    """

    # MultiIndex cannot be plotted; no point in allowing them here
    multiindex_dims = {
        idx.dim
        for idx in darray.xindexes.get_unique()
        if isinstance(idx, PandasMultiIndex)
    }

    valid_xy = (set(darray.dims) | set(darray.coords)) - multiindex_dims

    if (xy is not None) and (xy not in valid_xy):
        valid_xy_str = "', '".join(sorted(str(v) for v in valid_xy))
        raise ValueError(
            f"{name} must be one of None, '{valid_xy_str}'. Received '{xy}' instead."
        )


def get_axis(
    figsize: Iterable[float] | None = None,
    size: float | None = None,
    aspect: AspectOptions = None,
    ax: Axes | None = None,
    **subplot_kws: Any,
) -> Axes:
    if TYPE_CHECKING:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
    else:
        mpl = attempt_import("matplotlib")
        plt = attempt_import("matplotlib.pyplot")

    if figsize is not None:
        if ax is not None:
            raise ValueError("cannot provide both `figsize` and `ax` arguments")
        if size is not None:
            raise ValueError("cannot provide both `figsize` and `size` arguments")
        _, ax = plt.subplots(figsize=figsize, subplot_kw=subplot_kws)
        return ax

    if size is not None:
        if ax is not None:
            raise ValueError("cannot provide both `size` and `ax` arguments")
        if aspect is None or aspect == "auto":
            width, height = mpl.rcParams["figure.figsize"]
            faspect = width / height
        elif aspect == "equal":
            faspect = 1
        else:
            faspect = aspect
        figsize = (size * faspect, size)
        _, ax = plt.subplots(figsize=figsize, subplot_kw=subplot_kws)
        return ax

    if aspect is not None:
        raise ValueError("cannot provide `aspect` argument without `size`")

    if subplot_kws and ax is not None:
        raise ValueError("cannot use subplot_kws with existing ax")

    if ax is None:
        ax = _maybe_gca(**subplot_kws)

    return ax


def _maybe_gca(**subplot_kws: Any) -> Axes:
    import matplotlib.pyplot as plt

    # can call gcf unconditionally: either it exists or would be created by plt.axes
    f = plt.gcf()

    # only call gca if an active axes exists
    if f.axes:
        # can not pass kwargs to active axes
        return plt.gca()

    return plt.axes(**subplot_kws)


def _get_units_from_attrs(da: DataArray) -> str:
    """Extracts and formats the unit/units from a attributes."""
    pint_array_type = DuckArrayModule("pint").type
    units = " [{}]"
    if isinstance(da.data, pint_array_type):
        return units.format(str(da.data.units))
    if "units" in da.attrs:
        return units.format(da.attrs["units"])
    if "unit" in da.attrs:
        return units.format(da.attrs["unit"])
    return ""


def label_from_attrs(da: DataArray | None, extra: str = "") -> str:
    """Makes informative labels if variable metadata (attrs) follows
    CF conventions."""
    if da is None:
        return ""

    name: str = "{}"
    if "long_name" in da.attrs:
        name = name.format(da.attrs["long_name"])
    elif "standard_name" in da.attrs:
        name = name.format(da.attrs["standard_name"])
    elif da.name is not None:
        name = name.format(da.name)
    else:
        name = ""

    units = _get_units_from_attrs(da)

    # Treat `name` differently if it's a latex sequence
    if name.startswith("$") and (name.count("$") % 2 == 0):
        return "$\n$".join(
            textwrap.wrap(name + extra + units, 60, break_long_words=False)
        )
    else:
        return "\n".join(textwrap.wrap(name + extra + units, 30))


def _interval_to_mid_points(array: Iterable[pd.Interval]) -> np.ndarray:
    """
    Helper function which returns an array
    with the Intervals' mid points.
    """

    return np.array([x.mid for x in array])


def _interval_to_bound_points(array: Sequence[pd.Interval]) -> np.ndarray:
    """
    Helper function which returns an array
    with the Intervals' boundaries.
    """

    array_boundaries = np.array([x.left for x in array])
    array_boundaries = np.concatenate((array_boundaries, np.array([array[-1].right])))

    return array_boundaries


def _interval_to_double_bound_points(
    xarray: Iterable[pd.Interval], yarray: Iterable
) -> tuple[np.ndarray, np.ndarray]:
    """
    Helper function to deal with a xarray consisting of pd.Intervals. Each
    interval is replaced with both boundaries. I.e. the length of xarray
    doubles. yarray is modified so it matches the new shape of xarray.
    """

    xarray1 = np.array([x.left for x in xarray])
    xarray2 = np.array([x.right for x in xarray])

    xarray_out = np.array(
        list(itertools.chain.from_iterable(zip(xarray1, xarray2, strict=True)))
    )
    yarray_out = np.array(
        list(itertools.chain.from_iterable(zip(yarray, yarray, strict=True)))
    )

    return xarray_out, yarray_out


def _resolve_intervals_1dplot(
    xval: np.ndarray, yval: np.ndarray, kwargs: dict
) -> tuple[np.ndarray, np.ndarray, str, str, dict]:
    """
    Helper function to replace the values of x and/or y coordinate arrays
    containing pd.Interval with their mid-points or - for step plots - double
    points which double the length.
    """
    x_suffix = ""
    y_suffix = ""

    # Is it a step plot? (see matplotlib.Axes.step)
    if kwargs.get("drawstyle", "").startswith("steps-"):
        remove_drawstyle = False

        # Convert intervals to double points
        x_is_interval = _valid_other_type(xval, pd.Interval)
        y_is_interval = _valid_other_type(yval, pd.Interval)
        if x_is_interval and y_is_interval:
            raise TypeError("Can't step plot intervals against intervals.")
        elif x_is_interval:
            xval, yval = _interval_to_double_bound_points(xval, yval)
            remove_drawstyle = True
        elif y_is_interval:
            yval, xval = _interval_to_double_bound_points(yval, xval)
            remove_drawstyle = True

        # Remove steps-* to be sure that matplotlib is not confused
        if remove_drawstyle:
            del kwargs["drawstyle"]

    # Is it another kind of plot?
    else:
        # Convert intervals to mid points and adjust labels
        if _valid_other_type(xval, pd.Interval):
            xval = _interval_to_mid_points(xval)
            x_suffix = "_center"
        if _valid_other_type(yval, pd.Interval):
            yval = _interval_to_mid_points(yval)
            y_suffix = "_center"

    # return converted arguments
    return xval, yval, x_suffix, y_suffix, kwargs


def _resolve_intervals_2dplot(val, func_name):
    """
    Helper function to replace the values of a coordinate array containing
    pd.Interval with their mid-points or - for pcolormesh - boundaries which
    increases length by 1.
    """
    label_extra = ""
    if _valid_other_type(val, pd.Interval):
        if func_name == "pcolormesh":
            val = _interval_to_bound_points(val)
        else:
            val = _interval_to_mid_points(val)
            label_extra = "_center"

    return val, label_extra


def _valid_other_type(
    x: ArrayLike, types: type[object] | tuple[type[object], ...]
) -> bool:
    """
    Do all elements of x have a type from types?
    """
    return all(isinstance(el, types) for el in np.ravel(x))


def _valid_numpy_subdtype(x, numpy_types):
    """
    Is any dtype from numpy_types superior to the dtype of x?
    """
    # If any of the types given in numpy_types is understood as numpy.generic,
    # all possible x will be considered valid.  This is probably unwanted.
    for t in numpy_types:
        assert not np.issubdtype(np.generic, t)

    return any(np.issubdtype(x.dtype, t) for t in numpy_types)


def _ensure_plottable(*args) -> None:
    """
    Raise exception if there is anything in args that can't be plotted on an
    axis by matplotlib.
    """
    numpy_types: tuple[type[object], ...] = (
        np.floating,
        np.integer,
        np.timedelta64,
        np.datetime64,
        np.bool_,
        np.str_,
    )
    other_types: tuple[type[object], ...] = (datetime, date)
    cftime_datetime_types: tuple[type[object], ...] = (
        () if cftime is None else (cftime.datetime,)
    )
    other_types += cftime_datetime_types

    for x in args:
        if not (
            _valid_numpy_subdtype(np.asarray(x), numpy_types)
            or _valid_other_type(np.asarray(x), other_types)
        ):
            raise TypeError(
                "Plotting requires coordinates to be numeric, boolean, "
                "or dates of type numpy.datetime64, "
                "datetime.datetime, cftime.datetime or "
                f"pandas.Interval. Received data of type {np.asarray(x).dtype} instead."
            )
        if _valid_other_type(np.asarray(x), cftime_datetime_types):
            if nc_time_axis_available:
                # Register cftime datetypes to matplotlib.units.registry,
                # otherwise matplotlib will raise an error:
                import nc_time_axis  # noqa: F401
            else:
                raise ImportError(
                    "Plotting of arrays of cftime.datetime "
                    "objects or arrays indexed by "
                    "cftime.datetime objects requires the "
                    "optional `nc-time-axis` (v1.2.0 or later) "
                    "package."
                )


def _is_numeric(arr):
    numpy_types = [np.floating, np.integer]
    return _valid_numpy_subdtype(arr, numpy_types)


def _add_colorbar(primitive, ax, cbar_ax, cbar_kwargs, cmap_params):
    cbar_kwargs.setdefault("extend", cmap_params["extend"])
    if cbar_ax is None:
        cbar_kwargs.setdefault("ax", ax)
    else:
        cbar_kwargs.setdefault("cax", cbar_ax)

    # dont pass extend as kwarg if it is in the mappable
    if hasattr(primitive, "extend"):
        cbar_kwargs.pop("extend")

    fig = ax.get_figure()
    cbar = fig.colorbar(primitive, **cbar_kwargs)

    return cbar


def _rescale_imshow_rgb(darray, vmin, vmax, robust):
    assert robust or vmin is not None or vmax is not None

    # Calculate vmin and vmax automatically for `robust=True`
    if robust:
        if vmax is None:
            vmax = np.nanpercentile(darray, 100 - ROBUST_PERCENTILE)
        if vmin is None:
            vmin = np.nanpercentile(darray, ROBUST_PERCENTILE)
    # If not robust and one bound is None, calculate the default other bound
    # and check that an interval between them exists.
    elif vmax is None:
        vmax = 255 if np.issubdtype(darray.dtype, np.integer) else 1
        if vmax < vmin:
            raise ValueError(
                f"vmin={vmin!r} is less than the default vmax ({vmax!r}) - you must supply "
                "a vmax > vmin in this case."
            )
    elif vmin is None:
        vmin = 0
        if vmin > vmax:
            raise ValueError(
                f"vmax={vmax!r} is less than the default vmin (0) - you must supply "
                "a vmin < vmax in this case."
            )
    # Scale interval [vmin .. vmax] to [0 .. 1], with darray as 64-bit float
    # to avoid precision loss, integer over/underflow, etc with extreme inputs.
    # After scaling, downcast to 32-bit float.  This substantially reduces
    # memory usage after we hand `darray` off to matplotlib.
    darray = ((darray.astype("f8") - vmin) / (vmax - vmin)).astype("f4")
    return np.minimum(np.maximum(darray, 0), 1)


def _update_axes(
    ax: Axes,
    xincrease: bool | None,
    yincrease: bool | None,
    xscale: ScaleOptions = None,
    yscale: ScaleOptions = None,
    xticks: ArrayLike | None = None,
    yticks: ArrayLike | None = None,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
) -> None:
    """
    Update axes with provided parameters
    """
    if xincrease is None:
        pass
    elif (xincrease and ax.xaxis_inverted()) or (
        not xincrease and not ax.xaxis_inverted()
    ):
        ax.invert_xaxis()

    if yincrease is None:
        pass
    elif (yincrease and ax.yaxis_inverted()) or (
        not yincrease and not ax.yaxis_inverted()
    ):
        ax.invert_yaxis()

    # The default xscale, yscale needs to be None.
    # If we set a scale it resets the axes formatters,
    # This means that set_xscale('linear') on a datetime axis
    # will remove the date labels. So only set the scale when explicitly
    # asked to. https://github.com/matplotlib/matplotlib/issues/8740
    if xscale is not None:
        ax.set_xscale(xscale)
    if yscale is not None:
        ax.set_yscale(yscale)

    if xticks is not None:
        ax.set_xticks(xticks)
    if yticks is not None:
        ax.set_yticks(yticks)

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)


def _is_monotonic(coord, axis=0):
    """
    >>> _is_monotonic(np.array([0, 1, 2]))
    np.True_
    >>> _is_monotonic(np.array([2, 1, 0]))
    np.True_
    >>> _is_monotonic(np.array([0, 2, 1]))
    np.False_
    """
    if coord.shape[axis] < 3:
        return True
    else:
        n = coord.shape[axis]
        delta_pos = coord.take(np.arange(1, n), axis=axis) >= coord.take(
            np.arange(0, n - 1), axis=axis
        )
        delta_neg = coord.take(np.arange(1, n), axis=axis) <= coord.take(
            np.arange(0, n - 1), axis=axis
        )
        return np.all(delta_pos) or np.all(delta_neg)


def _infer_interval_breaks(coord, axis=0, scale=None, check_monotonic=False):
    """
    >>> _infer_interval_breaks(np.arange(5))
    array([-0.5,  0.5,  1.5,  2.5,  3.5,  4.5])
    >>> _infer_interval_breaks([[0, 1], [3, 4]], axis=1)
    array([[-0.5,  0.5,  1.5],
           [ 2.5,  3.5,  4.5]])
    >>> _infer_interval_breaks(np.logspace(-2, 2, 5), scale="log")
    array([3.16227766e-03, 3.16227766e-02, 3.16227766e-01, 3.16227766e+00,
           3.16227766e+01, 3.16227766e+02])
    """
    coord = np.asarray(coord)

    if check_monotonic and not _is_monotonic(coord, axis=axis):
        raise ValueError(
            "The input coordinate is not sorted in increasing "
            f"order along axis {axis}. This can lead to unexpected "
            "results. Consider calling the `sortby` method on "
            "the input DataArray. To plot data with categorical "
            "axes, consider using the `heatmap` function from "
            "the `seaborn` statistical plotting library."
        )

    # If logscale, compute the intervals in the logarithmic space
    if scale == "log":
        if (coord <= 0).any():
            raise ValueError(
                "Found negative or zero value in coordinates. "
                "Coordinates must be positive on logscale plots."
            )
        coord = np.log10(coord)

    deltas = 0.5 * np.diff(coord, axis=axis)
    if deltas.size == 0:
        deltas = np.array(0.0)
    first = np.take(coord, [0], axis=axis) - np.take(deltas, [0], axis=axis)
    last = np.take(coord, [-1], axis=axis) + np.take(deltas, [-1], axis=axis)
    trim_last = tuple(
        slice(None, -1) if n == axis else slice(None) for n in range(coord.ndim)
    )
    interval_breaks = np.concatenate(
        [first, coord[trim_last] + deltas, last], axis=axis
    )
    if scale == "log":
        # Recovert the intervals into the linear space
        return np.power(10, interval_breaks)
    return interval_breaks


def _process_cmap_cbar_kwargs(
    func,
    data,
    cmap=None,
    colors=None,
    cbar_kwargs: Iterable[tuple[str, Any]] | Mapping[str, Any] | None = None,
    levels=None,
    _is_facetgrid=False,
    **kwargs,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """
    Parameters
    ----------
    func : plotting function
    data : ndarray,
        Data values

    Returns
    -------
    cmap_params : dict
    cbar_kwargs : dict
    """
    if func.__name__ == "surface":
        # Leave user to specify cmap settings for surface plots
        kwargs["cmap"] = cmap
        return {
            k: kwargs.get(k)
            for k in ["vmin", "vmax", "cmap", "extend", "levels", "norm"]
        }, {}

    cbar_kwargs = {} if cbar_kwargs is None else dict(cbar_kwargs)

    # colors is mutually exclusive with cmap
    if cmap and colors:
        raise ValueError("Can't specify both cmap and colors.")

    # colors is only valid when levels is supplied or the plot is of type
    # contour or contourf
    if colors and (("contour" not in func.__name__) and (levels is None)):
        raise ValueError("Can only specify colors with contour or levels")

    # we should not be getting a list of colors in cmap anymore
    # is there a better way to do this test?
    if isinstance(cmap, list | tuple):
        raise ValueError(
            "Specifying a list of colors in cmap is deprecated. "
            "Use colors keyword instead."
        )

    cmap_kwargs = {
        "plot_data": data,
        "levels": levels,
        "cmap": colors or cmap,
        "filled": func.__name__ != "contour",
    }

    cmap_args = getfullargspec(_determine_cmap_params).args
    cmap_kwargs.update((a, kwargs[a]) for a in cmap_args if a in kwargs)
    if not _is_facetgrid:
        cmap_params = _determine_cmap_params(**cmap_kwargs)
    else:
        cmap_params = {
            k: cmap_kwargs[k]
            for k in ["vmin", "vmax", "cmap", "extend", "levels", "norm"]
        }

    return cmap_params, cbar_kwargs


def _get_nice_quiver_magnitude(u, v):
    import matplotlib as mpl

    ticker = mpl.ticker.MaxNLocator(3)
    mean = np.mean(np.hypot(u.to_numpy(), v.to_numpy()))
    magnitude = ticker.tick_values(0, mean)[-2]
    return magnitude


# Copied from matplotlib, tweaked so func can return strings.
# https://github.com/matplotlib/matplotlib/issues/19555
def legend_elements(
    self, prop="colors", num="auto", fmt=None, func=lambda x: x, **kwargs
):
    """
    Create legend handles and labels for a PathCollection.

    Each legend handle is a `.Line2D` representing the Path that was drawn,
    and each label is a string what each Path represents.

    This is useful for obtaining a legend for a `~.Axes.scatter` plot;
    e.g.::

        scatter = plt.scatter([1, 2, 3], [4, 5, 6], c=[7, 2, 3])
        plt.legend(*scatter.legend_elements())

    creates three legend elements, one for each color with the numerical
    values passed to *c* as the labels.

    Also see the :ref:`automatedlegendcreation` example.


    Parameters
    ----------
    prop : {"colors", "sizes"}, default: "colors"
        If "colors", the legend handles will show the different colors of
        the collection. If "sizes", the legend will show the different
        sizes. To set both, use *kwargs* to directly edit the `.Line2D`
        properties.
    num : int, None, "auto" (default), array-like, or `~.ticker.Locator`
        Target number of elements to create.
        If None, use all unique elements of the mappable array. If an
        integer, target to use *num* elements in the normed range.
        If *"auto"*, try to determine which option better suits the nature
        of the data.
        The number of created elements may slightly deviate from *num* due
        to a `~.ticker.Locator` being used to find useful locations.
        If a list or array, use exactly those elements for the legend.
        Finally, a `~.ticker.Locator` can be provided.
    fmt : str, `~matplotlib.ticker.Formatter`, or None (default)
        The format or formatter to use for the labels. If a string must be
        a valid input for a `~.StrMethodFormatter`. If None (the default),
        use a `~.ScalarFormatter`.
    func : function, default: ``lambda x: x``
        Function to calculate the labels.  Often the size (or color)
        argument to `~.Axes.scatter` will have been pre-processed by the
        user using a function ``s = f(x)`` to make the markers visible;
        e.g. ``size = np.log10(x)``.  Providing the inverse of this
        function here allows that pre-processing to be inverted, so that
        the legend labels have the correct values; e.g. ``func = lambda
        x: 10**x``.
    **kwargs
        Allowed keyword arguments are *color* and *size*. E.g. it may be
        useful to set the color of the markers if *prop="sizes"* is used;
        similarly to set the size of the markers if *prop="colors"* is
        used. Any further parameters are passed onto the `.Line2D`
        instance. This may be useful to e.g. specify a different
        *markeredgecolor* or *alpha* for the legend handles.

    Returns
    -------
    handles : list of `.Line2D`
        Visual representation of each element of the legend.
    labels : list of str
        The string labels for elements of the legend.
    """
    import matplotlib as mpl

    mlines = mpl.lines

    handles = []
    labels = []

    if prop == "colors":
        arr = self.get_array()
        if arr is None:
            warnings.warn(
                "Collection without array used. Make sure to "
                "specify the values to be colormapped via the "
                "`c` argument.",
                stacklevel=2,
            )
            return handles, labels
        _size = kwargs.pop("size", mpl.rcParams["lines.markersize"])

        def _get_color_and_size(value):
            return self.cmap(self.norm(value)), _size

    elif prop == "sizes":
        if isinstance(self, mpl.collections.LineCollection):
            arr = self.get_linewidths()
        else:
            arr = self.get_sizes()
        _color = kwargs.pop("color", "k")

        def _get_color_and_size(value):
            return _color, np.sqrt(value)

    else:
        raise ValueError(
            "Valid values for `prop` are 'colors' or "
            f"'sizes'. You supplied '{prop}' instead."
        )

    # Get the unique values and their labels:
    values = np.unique(arr)
    label_values = np.asarray(func(values))
    label_values_are_numeric = np.issubdtype(label_values.dtype, np.number)

    # Handle the label format:
    if fmt is None and label_values_are_numeric:
        fmt = mpl.ticker.ScalarFormatter(useOffset=False, useMathText=True)
    elif fmt is None and not label_values_are_numeric:
        fmt = mpl.ticker.StrMethodFormatter("{x}")
    elif isinstance(fmt, str):
        fmt = mpl.ticker.StrMethodFormatter(fmt)
    fmt.create_dummy_axis()

    if num == "auto":
        num = 9
        if len(values) <= num:
            num = None

    if label_values_are_numeric:
        label_values_min = label_values.min()
        label_values_max = label_values.max()
        fmt.axis.set_view_interval(label_values_min, label_values_max)
        fmt.axis.set_data_interval(label_values_min, label_values_max)

        if num is not None:
            # Labels are numerical but larger than the target
            # number of elements, reduce to target using matplotlibs
            # ticker classes:
            if isinstance(num, mpl.ticker.Locator):
                loc = num
            elif np.iterable(num):
                loc = mpl.ticker.FixedLocator(num)
            else:
                num = int(num)
                loc = mpl.ticker.MaxNLocator(
                    nbins=num, min_n_ticks=num - 1, steps=[1, 2, 2.5, 3, 5, 6, 8, 10]
                )

            # Get nicely spaced label_values:
            label_values = loc.tick_values(label_values_min, label_values_max)

            # Remove extrapolated label_values:
            cond = (label_values >= label_values_min) & (
                label_values <= label_values_max
            )
            label_values = label_values[cond]

            # Get the corresponding values by creating a linear interpolant
            # with small step size:
            values_interp = np.linspace(values.min(), values.max(), 256)
            label_values_interp = func(values_interp)
            ix = np.argsort(label_values_interp)
            values = np.interp(label_values, label_values_interp[ix], values_interp[ix])
    elif num is not None and not label_values_are_numeric:
        # Labels are not numerical so modifying label_values is not
        # possible, instead filter the array with nicely distributed
        # indexes:
        if type(num) is int:
            loc = mpl.ticker.LinearLocator(num)
        else:
            raise ValueError("`num` only supports integers for non-numeric labels.")

        ind = loc.tick_values(0, len(label_values) - 1).astype(int)
        label_values = label_values[ind]
        values = values[ind]

    # Some formatters requires set_locs:
    if hasattr(fmt, "set_locs"):
        fmt.set_locs(label_values)

    # Default settings for handles, add or override with kwargs:
    kw = dict(markeredgewidth=self.get_linewidths()[0], alpha=self.get_alpha())
    kw.update(kwargs)

    for val, lab in zip(values, label_values, strict=True):
        color, size = _get_color_and_size(val)

        if isinstance(self, mpl.collections.PathCollection):
            kw.update(linestyle="", marker=self.get_paths()[0], markersize=size)
        elif isinstance(self, mpl.collections.LineCollection):
            kw.update(linestyle=self.get_linestyle()[0], linewidth=size)

        h = mlines.Line2D([0], [0], color=color, **kw)

        handles.append(h)
        labels.append(fmt(lab))

    return handles, labels


def _legend_add_subtitle(handles, labels, text):
    """Add a subtitle to legend handles."""
    import matplotlib.pyplot as plt

    if text and len(handles) > 1:
        # Create a blank handle that's not visible, the
        # invisibility will be used to discern which are subtitles
        # or not:
        blank_handle = plt.Line2D([], [], label=text)
        blank_handle.set_visible(False)

        # Subtitles are shown first:
        handles = [blank_handle] + handles
        labels = [text] + labels

    return handles, labels


def _adjust_legend_subtitles(legend):
    """Make invisible-handle "subtitles" entries look more like titles."""
    import matplotlib.pyplot as plt

    # Legend title not in rcParams until 3.0
    font_size = plt.rcParams.get("legend.title_fontsize", None)
    hpackers = legend.findobj(plt.matplotlib.offsetbox.VPacker)[0].get_children()
    hpackers = [v for v in hpackers if isinstance(v, plt.matplotlib.offsetbox.HPacker)]
    for hpack in hpackers:
        areas = hpack.get_children()
        if len(areas) < 2:
            continue
        draw_area, text_area = areas

        handles = draw_area.get_children()

        # Assume that all artists that are not visible are
        # subtitles:
        if not all(artist.get_visible() for artist in handles):
            # Remove the dummy marker which will bring the text
            # more to the center:
            draw_area.set_width(0)
            for text in text_area.get_children():
                if font_size is not None:
                    # The sutbtitles should have the same font size
                    # as normal legend titles:
                    text.set_size(font_size)


def _infer_meta_data(ds, x, y, hue, hue_style, add_guide, funcname):
    dvars = set(ds.variables.keys())
    error_msg = f" must be one of ({', '.join(sorted(str(v) for v in dvars))})"

    if x not in dvars:
        raise ValueError(f"Expected 'x' {error_msg}. Received {x} instead.")

    if y not in dvars:
        raise ValueError(f"Expected 'y' {error_msg}. Received {y} instead.")

    if hue is not None and hue not in dvars:
        raise ValueError(f"Expected 'hue' {error_msg}. Received {hue} instead.")

    if hue:
        hue_is_numeric = _is_numeric(ds[hue].values)

        if hue_style is None:
            hue_style = "continuous" if hue_is_numeric else "discrete"

        if not hue_is_numeric and (hue_style == "continuous"):
            raise ValueError(
                f"Cannot create a colorbar for a non numeric coordinate: {hue}"
            )

        if add_guide is None or add_guide is True:
            add_colorbar = hue_style == "continuous"
            add_legend = hue_style == "discrete"
        else:
            add_colorbar = False
            add_legend = False
    else:
        if add_guide is True and funcname not in ("quiver", "streamplot"):
            raise ValueError("Cannot set add_guide when hue is None.")
        add_legend = False
        add_colorbar = False

    if (add_guide or add_guide is None) and funcname == "quiver":
        add_quiverkey = True
        if hue:
            add_colorbar = True
            if not hue_style:
                hue_style = "continuous"
            elif hue_style != "continuous":
                raise ValueError(
                    "hue_style must be 'continuous' or None for .plot.quiver or "
                    ".plot.streamplot"
                )
    else:
        add_quiverkey = False

    if (add_guide or add_guide is None) and funcname == "streamplot" and hue:
        add_colorbar = True
        if not hue_style:
            hue_style = "continuous"
        elif hue_style != "continuous":
            raise ValueError(
                "hue_style must be 'continuous' or None for .plot.quiver or "
                ".plot.streamplot"
            )

    if hue_style is not None and hue_style not in ["discrete", "continuous"]:
        raise ValueError("hue_style must be either None, 'discrete' or 'continuous'.")

    if hue:
        hue_label = label_from_attrs(ds[hue])
        hue = ds[hue]
    else:
        hue_label = None
        hue = None

    return {
        "add_colorbar": add_colorbar,
        "add_legend": add_legend,
        "add_quiverkey": add_quiverkey,
        "hue_label": hue_label,
        "hue_style": hue_style,
        "xlabel": label_from_attrs(ds[x]),
        "ylabel": label_from_attrs(ds[y]),
        "hue": hue,
    }


@overload
def _parse_size(
    data: None,
    norm: tuple[float | None, float | None, bool] | Normalize | None,
) -> None: ...


@overload
def _parse_size(
    data: DataArray,
    norm: tuple[float | None, float | None, bool] | Normalize | None,
) -> pd.Series: ...


# copied from seaborn
def _parse_size(
    data: DataArray | None,
    norm: tuple[float | None, float | None, bool] | Normalize | None,
) -> pd.Series | None:
    import matplotlib as mpl

    if data is None:
        return None

    flatdata = data.values.flatten()

    if not _is_numeric(flatdata):
        levels = np.unique(flatdata)
        numbers = np.arange(1, 1 + len(levels))[::-1]
    else:
        levels = numbers = np.sort(np.unique(flatdata))

    min_width, default_width, max_width = _MARKERSIZE_RANGE
    # width_range = min_width, max_width

    if norm is None:
        norm = mpl.colors.Normalize()
    elif isinstance(norm, tuple):
        norm = mpl.colors.Normalize(*norm)
    elif not isinstance(norm, mpl.colors.Normalize):
        err = "``size_norm`` must be None, tuple, or Normalize object."
        raise ValueError(err)
    assert isinstance(norm, mpl.colors.Normalize)

    norm.clip = True
    if not norm.scaled():
        norm(np.asarray(numbers))
    # limits = norm.vmin, norm.vmax

    scl = norm(numbers)
    widths = np.asarray(min_width + scl * (max_width - min_width))
    if scl.mask.any():
        widths[scl.mask] = 0
    sizes = dict(zip(levels, widths, strict=True))

    return pd.Series(sizes)


class _Normalize(Sequence):
    """
    Normalize numerical or categorical values to numerical values.

    The class includes helper methods that simplifies transforming to
    and from normalized values.

    Parameters
    ----------
    data : DataArray
        DataArray to normalize.
    width : Sequence of three numbers, optional
        Normalize the data to these (min, default, max) values.
        The default is None.
    """

    _data: DataArray | None
    _data_unique: np.ndarray
    _data_unique_index: np.ndarray
    _data_unique_inverse: np.ndarray
    _data_is_numeric: bool
    _width: tuple[float, float, float] | None

    __slots__ = (
        "_data",
        "_data_is_numeric",
        "_data_unique",
        "_data_unique_index",
        "_data_unique_inverse",
        "_width",
    )

    def __init__(
        self,
        data: DataArray | None,
        width: tuple[float, float, float] | None = None,
        _is_facetgrid: bool = False,
    ) -> None:
        self._data = data
        self._width = width if not _is_facetgrid else None

        pint_array_type = DuckArrayModule("pint").type
        to_unique = (
            data.to_numpy()  # type: ignore[union-attr]
            if isinstance(data if data is None else data.data, pint_array_type)
            else data
        )
        data_unique, data_unique_inverse = np.unique(to_unique, return_inverse=True)  # type: ignore[call-overload]
        self._data_unique = data_unique
        self._data_unique_index = np.arange(0, data_unique.size)
        self._data_unique_inverse = data_unique_inverse
        self._data_is_numeric = False if data is None else _is_numeric(data)

    def __repr__(self) -> str:
        with np.printoptions(precision=4, suppress=True, threshold=5):
            return (
                f"<_Normalize(data, width={self._width})>\n"
                f"{self._data_unique} -> {self._values_unique}"
            )

    def __len__(self) -> int:
        return len(self._data_unique)

    def __getitem__(self, key):
        return self._data_unique[key]

    @property
    def data(self) -> DataArray | None:
        return self._data

    @property
    def data_is_numeric(self) -> bool:
        """
        Check if data is numeric.

        Examples
        --------
        >>> a = xr.DataArray(["b", "a", "a", "b", "c"])
        >>> _Normalize(a).data_is_numeric
        False

        >>> a = xr.DataArray([0.5, 0, 0, 0.5, 2, 3])
        >>> _Normalize(a).data_is_numeric
        True

        >>> # TODO: Datetime should be numeric right?
        >>> a = xr.DataArray(pd.date_range("2000-1-1", periods=4))
        >>> _Normalize(a).data_is_numeric
        False

        # TODO: Timedelta should be numeric right?
        >>> a = xr.DataArray(pd.timedelta_range("-1D", periods=4, freq="D"))
        >>> _Normalize(a).data_is_numeric
        True
        """
        return self._data_is_numeric

    @overload
    def _calc_widths(self, y: np.ndarray) -> np.ndarray: ...

    @overload
    def _calc_widths(self, y: DataArray) -> DataArray: ...

    def _calc_widths(self, y: np.ndarray | DataArray) -> np.ndarray | DataArray:
        """
        Normalize the values so they're in between self._width.
        """
        if self._width is None:
            return y

        xmin, xdefault, xmax = self._width

        diff_maxy_miny = np.max(y) - np.min(y)
        if diff_maxy_miny == 0:
            # Use default with if y is constant:
            widths = xdefault + 0 * y
        else:
            # Normalize in between xmin and xmax:
            k = (y - np.min(y)) / diff_maxy_miny
            widths = xmin + k * (xmax - xmin)
        return widths

    @overload
    def _indexes_centered(self, x: np.ndarray) -> np.ndarray: ...

    @overload
    def _indexes_centered(self, x: DataArray) -> DataArray: ...

    def _indexes_centered(self, x: np.ndarray | DataArray) -> np.ndarray | DataArray:
        """
        Offset indexes to make sure being in the center of self.levels.
        ["a", "b", "c"] -> [1, 3, 5]
        """
        return x * 2 + 1

    @property
    def values(self) -> DataArray | None:
        """
        Return a normalized number array for the unique levels.

        Examples
        --------
        >>> a = xr.DataArray(["b", "a", "a", "b", "c"])
        >>> _Normalize(a).values
        <xarray.DataArray (dim_0: 5)> Size: 40B
        array([3, 1, 1, 3, 5])
        Dimensions without coordinates: dim_0

        >>> _Normalize(a, width=(18, 36, 72)).values
        <xarray.DataArray (dim_0: 5)> Size: 40B
        array([45., 18., 18., 45., 72.])
        Dimensions without coordinates: dim_0

        >>> a = xr.DataArray([0.5, 0, 0, 0.5, 2, 3])
        >>> _Normalize(a).values
        <xarray.DataArray (dim_0: 6)> Size: 48B
        array([0.5, 0. , 0. , 0.5, 2. , 3. ])
        Dimensions without coordinates: dim_0

        >>> _Normalize(a, width=(18, 36, 72)).values
        <xarray.DataArray (dim_0: 6)> Size: 48B
        array([27., 18., 18., 27., 54., 72.])
        Dimensions without coordinates: dim_0

        >>> _Normalize(a * 0, width=(18, 36, 72)).values
        <xarray.DataArray (dim_0: 6)> Size: 48B
        array([36., 36., 36., 36., 36., 36.])
        Dimensions without coordinates: dim_0

        """
        if self.data is None:
            return None

        val: DataArray
        if self.data_is_numeric:
            val = self.data
        else:
            arr = self._indexes_centered(self._data_unique_inverse)
            val = self.data.copy(data=arr.reshape(self.data.shape))

        return self._calc_widths(val)

    @property
    def _values_unique(self) -> np.ndarray | None:
        """
        Return unique values.

        Examples
        --------
        >>> a = xr.DataArray(["b", "a", "a", "b", "c"])
        >>> _Normalize(a)._values_unique
        array([1, 3, 5])

        >>> _Normalize(a, width=(18, 36, 72))._values_unique
        array([18., 45., 72.])

        >>> a = xr.DataArray([0.5, 0, 0, 0.5, 2, 3])
        >>> _Normalize(a)._values_unique
        array([0. , 0.5, 2. , 3. ])

        >>> _Normalize(a, width=(18, 36, 72))._values_unique
        array([18., 27., 54., 72.])
        """
        if self.data is None:
            return None

        val: np.ndarray
        if self.data_is_numeric:
            val = self._data_unique
        else:
            val = self._indexes_centered(self._data_unique_index)

        return self._calc_widths(val)

    @property
    def ticks(self) -> np.ndarray | None:
        """
        Return ticks for plt.colorbar if the data is not numeric.

        Examples
        --------
        >>> a = xr.DataArray(["b", "a", "a", "b", "c"])
        >>> _Normalize(a).ticks
        array([1, 3, 5])
        """
        val: np.ndarray | None
        if self.data_is_numeric:
            val = None
        else:
            val = self._indexes_centered(self._data_unique_index)

        return val

    @property
    def levels(self) -> np.ndarray:
        """
        Return discrete levels that will evenly bound self.values.
        ["a", "b", "c"] -> [0, 2, 4, 6]

        Examples
        --------
        >>> a = xr.DataArray(["b", "a", "a", "b", "c"])
        >>> _Normalize(a).levels
        array([0, 2, 4, 6])
        """
        return (
            np.append(self._data_unique_index, np.max(self._data_unique_index) + 1) * 2
        )

    @property
    def _lookup(self) -> pd.Series:
        if self._values_unique is None:
            raise ValueError("self.data can't be None.")

        return pd.Series(dict(zip(self._values_unique, self._data_unique, strict=True)))

    def _lookup_arr(self, x) -> np.ndarray:
        # Use reindex to be less sensitive to float errors. reindex only
        # works with sorted index.
        # Return as numpy array since legend_elements
        # seems to require that:
        return self._lookup.sort_index().reindex(x, method="nearest").to_numpy()

    @property
    def format(self) -> FuncFormatter:
        """
        Return a FuncFormatter that maps self.values elements back to
        the original value as a string. Useful with plt.colorbar.

        Examples
        --------
        >>> a = xr.DataArray([0.5, 0, 0, 0.5, 2, 3])
        >>> aa = _Normalize(a, width=(0, 0.5, 1))
        >>> aa._lookup
        0.000000    0.0
        0.166667    0.5
        0.666667    2.0
        1.000000    3.0
        dtype: float64
        >>> aa.format(1)
        '3.0'
        """
        import matplotlib.pyplot as plt

        def _func(x: Any, pos: Any | None = None):
            return f"{self._lookup_arr([x])[0]}"

        return plt.FuncFormatter(_func)

    @property
    def func(self) -> Callable[[Any, Any | None], Any]:
        """
        Return a lambda function that maps self.values elements back to
        the original value as a numpy array. Useful with ax.legend_elements.

        Examples
        --------
        >>> a = xr.DataArray([0.5, 0, 0, 0.5, 2, 3])
        >>> aa = _Normalize(a, width=(0, 0.5, 1))
        >>> aa._lookup
        0.000000    0.0
        0.166667    0.5
        0.666667    2.0
        1.000000    3.0
        dtype: float64
        >>> aa.func([0.16, 1])
        array([0.5, 3. ])
        """

        def _func(x: Any, pos: Any | None = None):
            return self._lookup_arr(x)

        return _func


def _determine_guide(
    hueplt_norm: _Normalize,
    sizeplt_norm: _Normalize,
    add_colorbar: bool | None = None,
    add_legend: bool | None = None,
    plotfunc_name: str | None = None,
) -> tuple[bool, bool]:
    if plotfunc_name == "hist":
        return False, False

    if (add_colorbar) and hueplt_norm.data is None:
        raise KeyError("Cannot create a colorbar when hue is None.")
    if add_colorbar is None:
        if hueplt_norm.data is not None:
            add_colorbar = True
        else:
            add_colorbar = False

    if add_legend and hueplt_norm.data is None and sizeplt_norm.data is None:
        raise KeyError("Cannot create a legend when hue and markersize is None.")
    if add_legend is None:
        if (
            not add_colorbar
            and (hueplt_norm.data is not None and hueplt_norm.data_is_numeric is False)
        ) or sizeplt_norm.data is not None:
            add_legend = True
        else:
            add_legend = False

    return add_colorbar, add_legend


def _add_legend(
    hueplt_norm: _Normalize,
    sizeplt_norm: _Normalize,
    primitive,
    legend_ax,
    plotfunc: str,
):
    primitive = primitive if isinstance(primitive, list) else [primitive]

    handles, labels = [], []
    for huesizeplt, prop in [
        (hueplt_norm, "colors"),
        (sizeplt_norm, "sizes"),
    ]:
        if huesizeplt.data is not None:
            # Get legend handles and labels that displays the
            # values correctly. Order might be different because
            # legend_elements uses np.unique instead of pd.unique,
            # FacetGrid.add_legend might have troubles with this:
            hdl, lbl = [], []
            for p in primitive:
                hdl_, lbl_ = legend_elements(p, prop, num="auto", func=huesizeplt.func)
                hdl += hdl_
                lbl += lbl_

            # Only save unique values:
            u, ind = np.unique(lbl, return_index=True)
            ind = np.argsort(ind)
            lbl = cast(list, u[ind].tolist())
            hdl = cast(list, np.array(hdl)[ind].tolist())

            # Add a subtitle:
            hdl, lbl = _legend_add_subtitle(hdl, lbl, label_from_attrs(huesizeplt.data))
            handles += hdl
            labels += lbl
    legend = legend_ax.legend(handles, labels, framealpha=0.5)
    _adjust_legend_subtitles(legend)

    return legend


def _guess_coords_to_plot(
    darray: DataArray,
    coords_to_plot: MutableMapping[str, Hashable | None],
    kwargs: dict,
    default_guess: tuple[str, ...] = ("x",),
    # TODO: Can this be normalized, plt.cbook.normalize_kwargs?
    ignore_guess_kwargs: tuple[tuple[str, ...], ...] = ((),),
) -> MutableMapping[str, Hashable]:
    """
    Guess what coords to plot if some of the values in coords_to_plot are None which
    happens when the user has not defined all available ways of visualizing
    the data.

    Parameters
    ----------
    darray : DataArray
        The DataArray to check for available coords.
    coords_to_plot : MutableMapping[str, Hashable]
        Coords defined by the user to plot.
    kwargs : dict
        Extra kwargs that will be sent to matplotlib.
    default_guess : Iterable[str], optional
        Default values and order to retrieve dims if values in dims_plot is
        missing, default: ("x", "hue", "size").
    ignore_guess_kwargs : tuple[tuple[str, ...], ...]
        Matplotlib arguments to ignore.

    Examples
    --------
    >>> ds = xr.tutorial.scatter_example_dataset(seed=42)
    >>> # Only guess x by default:
    >>> xr.plot.utils._guess_coords_to_plot(
    ...     ds.A,
    ...     coords_to_plot={"x": None, "z": None, "hue": None, "size": None},
    ...     kwargs={},
    ... )
    {'x': 'x', 'z': None, 'hue': None, 'size': None}

    >>> # Guess all plot dims with other default values:
    >>> xr.plot.utils._guess_coords_to_plot(
    ...     ds.A,
    ...     coords_to_plot={"x": None, "z": None, "hue": None, "size": None},
    ...     kwargs={},
    ...     default_guess=("x", "hue", "size"),
    ...     ignore_guess_kwargs=((), ("c", "color"), ("s",)),
    ... )
    {'x': 'x', 'z': None, 'hue': 'y', 'size': 'z'}

    >>> # Don't guess ´size´, since the matplotlib kwarg ´s´ has been defined:
    >>> xr.plot.utils._guess_coords_to_plot(
    ...     ds.A,
    ...     coords_to_plot={"x": None, "z": None, "hue": None, "size": None},
    ...     kwargs={"s": 5},
    ...     default_guess=("x", "hue", "size"),
    ...     ignore_guess_kwargs=((), ("c", "color"), ("s",)),
    ... )
    {'x': 'x', 'z': None, 'hue': 'y', 'size': None}

    >>> # Prioritize ´size´ over ´s´:
    >>> xr.plot.utils._guess_coords_to_plot(
    ...     ds.A,
    ...     coords_to_plot={"x": None, "z": None, "hue": None, "size": "x"},
    ...     kwargs={"s": 5},
    ...     default_guess=("x", "hue", "size"),
    ...     ignore_guess_kwargs=((), ("c", "color"), ("s",)),
    ... )
    {'x': 'y', 'z': None, 'hue': 'z', 'size': 'x'}
    """
    coords_to_plot_exist = {k: v for k, v in coords_to_plot.items() if v is not None}
    available_coords = tuple(
        k for k in darray.coords.keys() if k not in coords_to_plot_exist.values()
    )

    # If dims_plot[k] isn't defined then fill with one of the available dims, unless
    # one of related mpl kwargs has been used. This should have similar behaviour as
    # * plt.plot(x, y) -> Multiple lines with different colors if y is 2d.
    # * plt.plot(x, y, color="red") -> Multiple red lines if y is 2d.
    for k, dim, ign_kws in zip(
        default_guess, available_coords, ignore_guess_kwargs, strict=False
    ):
        if coords_to_plot.get(k, None) is None and all(
            kwargs.get(ign_kw) is None for ign_kw in ign_kws
        ):
            coords_to_plot[k] = dim

    for k, dim in coords_to_plot.items():
        _assert_valid_xy(darray, dim, k)

    return coords_to_plot


def _set_concise_date(ax: Axes, axis: Literal["x", "y", "z"] = "x") -> None:
    """
    Use ConciseDateFormatter which is meant to improve the
    strings chosen for the ticklabels, and to minimize the
    strings used in those tick labels as much as possible.

    https://matplotlib.org/stable/gallery/ticks/date_concise_formatter.html

    Parameters
    ----------
    ax : Axes
        Figure axes.
    axis : Literal["x", "y", "z"], optional
        Which axis to make concise. The default is "x".
    """
    import matplotlib.dates as mdates

    locator = mdates.AutoDateLocator()
    formatter = mdates.ConciseDateFormatter(locator)
    _axis = getattr(ax, f"{axis}axis")
    _axis.set_major_locator(locator)
    _axis.set_major_formatter(formatter)
