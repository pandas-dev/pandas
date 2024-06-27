from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Literal, TypedDict

from xarray.core.utils import FrozenDict

if TYPE_CHECKING:
    from matplotlib.colors import Colormap

    Options = Literal[
        "arithmetic_join",
        "cmap_divergent",
        "cmap_sequential",
        "display_max_rows",
        "display_values_threshold",
        "display_style",
        "display_width",
        "display_expand_attrs",
        "display_expand_coords",
        "display_expand_data_vars",
        "display_expand_data",
        "display_expand_groups",
        "display_expand_indexes",
        "display_default_indexes",
        "enable_cftimeindex",
        "file_cache_maxsize",
        "keep_attrs",
        "warn_for_unclosed_files",
        "use_bottleneck",
        "use_numbagg",
        "use_opt_einsum",
        "use_flox",
    ]

    class T_Options(TypedDict):
        arithmetic_broadcast: bool
        arithmetic_join: Literal["inner", "outer", "left", "right", "exact"]
        cmap_divergent: str | Colormap
        cmap_sequential: str | Colormap
        display_max_rows: int
        display_values_threshold: int
        display_style: Literal["text", "html"]
        display_width: int
        display_expand_attrs: Literal["default", True, False]
        display_expand_coords: Literal["default", True, False]
        display_expand_data_vars: Literal["default", True, False]
        display_expand_data: Literal["default", True, False]
        display_expand_groups: Literal["default", True, False]
        display_expand_indexes: Literal["default", True, False]
        display_default_indexes: Literal["default", True, False]
        enable_cftimeindex: bool
        file_cache_maxsize: int
        keep_attrs: Literal["default", True, False]
        warn_for_unclosed_files: bool
        use_bottleneck: bool
        use_flox: bool
        use_numbagg: bool
        use_opt_einsum: bool


OPTIONS: T_Options = {
    "arithmetic_broadcast": True,
    "arithmetic_join": "inner",
    "cmap_divergent": "RdBu_r",
    "cmap_sequential": "viridis",
    "display_max_rows": 12,
    "display_values_threshold": 200,
    "display_style": "html",
    "display_width": 80,
    "display_expand_attrs": "default",
    "display_expand_coords": "default",
    "display_expand_data_vars": "default",
    "display_expand_data": "default",
    "display_expand_groups": "default",
    "display_expand_indexes": "default",
    "display_default_indexes": False,
    "enable_cftimeindex": True,
    "file_cache_maxsize": 128,
    "keep_attrs": "default",
    "warn_for_unclosed_files": False,
    "use_bottleneck": True,
    "use_flox": True,
    "use_numbagg": True,
    "use_opt_einsum": True,
}

_JOIN_OPTIONS = frozenset(["inner", "outer", "left", "right", "exact"])
_DISPLAY_OPTIONS = frozenset(["text", "html"])


def _positive_integer(value: int) -> bool:
    return isinstance(value, int) and value > 0


_VALIDATORS = {
    "arithmetic_broadcast": lambda value: isinstance(value, bool),
    "arithmetic_join": _JOIN_OPTIONS.__contains__,
    "display_max_rows": _positive_integer,
    "display_values_threshold": _positive_integer,
    "display_style": _DISPLAY_OPTIONS.__contains__,
    "display_width": _positive_integer,
    "display_expand_attrs": lambda choice: choice in [True, False, "default"],
    "display_expand_coords": lambda choice: choice in [True, False, "default"],
    "display_expand_data_vars": lambda choice: choice in [True, False, "default"],
    "display_expand_data": lambda choice: choice in [True, False, "default"],
    "display_expand_indexes": lambda choice: choice in [True, False, "default"],
    "display_default_indexes": lambda choice: choice in [True, False, "default"],
    "enable_cftimeindex": lambda value: isinstance(value, bool),
    "file_cache_maxsize": _positive_integer,
    "keep_attrs": lambda choice: choice in [True, False, "default"],
    "use_bottleneck": lambda value: isinstance(value, bool),
    "use_numbagg": lambda value: isinstance(value, bool),
    "use_opt_einsum": lambda value: isinstance(value, bool),
    "use_flox": lambda value: isinstance(value, bool),
    "warn_for_unclosed_files": lambda value: isinstance(value, bool),
}


def _set_file_cache_maxsize(value) -> None:
    from xarray.backends.file_manager import FILE_CACHE

    FILE_CACHE.maxsize = value


def _warn_on_setting_enable_cftimeindex(enable_cftimeindex):
    warnings.warn(
        "The enable_cftimeindex option is now a no-op "
        "and will be removed in a future version of xarray.",
        FutureWarning,
    )


_SETTERS = {
    "enable_cftimeindex": _warn_on_setting_enable_cftimeindex,
    "file_cache_maxsize": _set_file_cache_maxsize,
}


def _get_boolean_with_default(option: Options, default: bool) -> bool:
    global_choice = OPTIONS[option]

    if global_choice == "default":
        return default
    elif isinstance(global_choice, bool):
        return global_choice
    else:
        raise ValueError(
            f"The global option {option} must be one of True, False or 'default'."
        )


def _get_keep_attrs(default: bool) -> bool:
    return _get_boolean_with_default("keep_attrs", default)


class set_options:
    """
    Set options for xarray in a controlled context.

    Parameters
    ----------
    arithmetic_join : {"inner", "outer", "left", "right", "exact"}, default: "inner"
        DataArray/Dataset alignment in binary operations:

        - "outer": use the union of object indexes
        - "inner": use the intersection of object indexes
        - "left": use indexes from the first object with each dimension
        - "right": use indexes from the last object with each dimension
        - "exact": instead of aligning, raise `ValueError` when indexes to be
          aligned are not equal
        - "override": if indexes are of same size, rewrite indexes to be
          those of the first object with that dimension. Indexes for the same
          dimension must have the same size in all objects.

    cmap_divergent : str or matplotlib.colors.Colormap, default: "RdBu_r"
        Colormap to use for divergent data plots. If string, must be
        matplotlib built-in colormap. Can also be a Colormap object
        (e.g. mpl.colormaps["magma"])
    cmap_sequential : str or matplotlib.colors.Colormap, default: "viridis"
        Colormap to use for nondivergent data plots. If string, must be
        matplotlib built-in colormap. Can also be a Colormap object
        (e.g. mpl.colormaps["magma"])
    display_expand_attrs : {"default", True, False}
        Whether to expand the attributes section for display of
        ``DataArray`` or ``Dataset`` objects. Can be

        * ``True`` : to always expand attrs
        * ``False`` : to always collapse attrs
        * ``default`` : to expand unless over a pre-defined limit
    display_expand_coords : {"default", True, False}
        Whether to expand the coordinates section for display of
        ``DataArray`` or ``Dataset`` objects. Can be

        * ``True`` : to always expand coordinates
        * ``False`` : to always collapse coordinates
        * ``default`` : to expand unless over a pre-defined limit
    display_expand_data : {"default", True, False}
        Whether to expand the data section for display of ``DataArray``
        objects. Can be

        * ``True`` : to always expand data
        * ``False`` : to always collapse data
        * ``default`` : to expand unless over a pre-defined limit
    display_expand_data_vars : {"default", True, False}
        Whether to expand the data variables section for display of
        ``Dataset`` objects. Can be

        * ``True`` : to always expand data variables
        * ``False`` : to always collapse data variables
        * ``default`` : to expand unless over a pre-defined limit
    display_expand_indexes : {"default", True, False}
        Whether to expand the indexes section for display of
        ``DataArray`` or ``Dataset``. Can be

        * ``True`` : to always expand indexes
        * ``False`` : to always collapse indexes
        * ``default`` : to expand unless over a pre-defined limit (always collapse for html style)
    display_max_rows : int, default: 12
        Maximum display rows.
    display_values_threshold : int, default: 200
        Total number of array elements which trigger summarization rather
        than full repr for variable data views (numpy arrays).
    display_style : {"text", "html"}, default: "html"
        Display style to use in jupyter for xarray objects.
    display_width : int, default: 80
        Maximum display width for ``repr`` on xarray objects.
    file_cache_maxsize : int, default: 128
        Maximum number of open files to hold in xarray's
        global least-recently-usage cached. This should be smaller than
        your system's per-process file descriptor limit, e.g.,
        ``ulimit -n`` on Linux.
    keep_attrs : {"default", True, False}
        Whether to keep attributes on xarray Datasets/dataarrays after
        operations. Can be

        * ``True`` : to always keep attrs
        * ``False`` : to always discard attrs
        * ``default`` : to use original logic that attrs should only
          be kept in unambiguous circumstances
    use_bottleneck : bool, default: True
        Whether to use ``bottleneck`` to accelerate 1D reductions and
        1D rolling reduction operations.
    use_flox : bool, default: True
        Whether to use ``numpy_groupies`` and `flox`` to
        accelerate groupby and resampling reductions.
    use_numbagg : bool, default: True
        Whether to use ``numbagg`` to accelerate reductions.
        Takes precedence over ``use_bottleneck`` when both are True.
    use_opt_einsum : bool, default: True
        Whether to use ``opt_einsum`` to accelerate dot products.
    warn_for_unclosed_files : bool, default: False
        Whether or not to issue a warning when unclosed files are
        deallocated. This is mostly useful for debugging.

    Examples
    --------
    It is possible to use ``set_options`` either as a context manager:

    >>> ds = xr.Dataset({"x": np.arange(1000)})
    >>> with xr.set_options(display_width=40):
    ...     print(ds)
    ...
    <xarray.Dataset> Size: 8kB
    Dimensions:  (x: 1000)
    Coordinates:
      * x        (x) int64 8kB 0 1 ... 999
    Data variables:
        *empty*

    Or to set global options:

    >>> xr.set_options(display_width=80)  # doctest: +ELLIPSIS
    <xarray.core.options.set_options object at 0x...>
    """

    def __init__(self, **kwargs):
        self.old = {}
        for k, v in kwargs.items():
            if k not in OPTIONS:
                raise ValueError(
                    f"argument name {k!r} is not in the set of valid options {set(OPTIONS)!r}"
                )
            if k in _VALIDATORS and not _VALIDATORS[k](v):
                if k == "arithmetic_join":
                    expected = f"Expected one of {_JOIN_OPTIONS!r}"
                elif k == "display_style":
                    expected = f"Expected one of {_DISPLAY_OPTIONS!r}"
                else:
                    expected = ""
                raise ValueError(
                    f"option {k!r} given an invalid value: {v!r}. " + expected
                )
            self.old[k] = OPTIONS[k]
        self._apply_update(kwargs)

    def _apply_update(self, options_dict):
        for k, v in options_dict.items():
            if k in _SETTERS:
                _SETTERS[k](v)
        OPTIONS.update(options_dict)

    def __enter__(self):
        return

    def __exit__(self, type, value, traceback):
        self._apply_update(self.old)


def get_options():
    """
    Get options for xarray.

    See Also
    ----------
    set_options

    """
    return FrozenDict(OPTIONS)
