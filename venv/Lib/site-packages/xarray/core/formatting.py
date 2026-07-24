"""String formatting routines for __repr__."""

from __future__ import annotations

import contextlib
import functools
import math
from collections import ChainMap, defaultdict
from collections.abc import Collection, Hashable, Mapping, Sequence
from datetime import datetime, timedelta
from itertools import chain, zip_longest
from reprlib import recursive_repr
from textwrap import indent
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
from pandas.errors import OutOfBoundsDatetime

from xarray.core.datatree_render import RenderDataTree
from xarray.core.duck_array_ops import array_all, array_any, array_equiv, astype, ravel
from xarray.core.extension_array import PandasExtensionArray
from xarray.core.indexing import (
    BasicIndexer,
    ExplicitlyIndexed,
    MemoryCachedArray,
)
from xarray.core.options import OPTIONS, _get_boolean_with_default
from xarray.core.treenode import group_subtrees
from xarray.core.utils import is_duck_array
from xarray.namedarray.pycompat import array_type, to_duck_array

if TYPE_CHECKING:
    from xarray.core.coordinates import AbstractCoordinates
    from xarray.core.datatree import DataTree
    from xarray.core.variable import Variable

UNITS = ("B", "kB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")


def pretty_print(x, numchars: int):
    """Given an object `x`, call `str(x)` and format the returned string so
    that it is numchars long, padding with trailing spaces or truncating with
    ellipses as necessary
    """
    s = maybe_truncate(x, numchars)
    return s + " " * max(numchars - len(s), 0)


def maybe_truncate(obj, maxlen=500):
    s = str(obj)
    if len(s) > maxlen:
        s = s[: (maxlen - 3)] + "..."
    return s


def wrap_indent(text, start="", length=None):
    if length is None:
        length = len(start)
    indent = "\n" + " " * length
    return start + indent.join(x for x in text.splitlines())


def _get_indexer_at_least_n_items(shape, n_desired, from_end):
    assert 0 < n_desired <= math.prod(shape)
    cum_items = np.cumprod(shape[::-1])
    n_steps = np.argmax(cum_items >= n_desired)
    stop = math.ceil(float(n_desired) / np.r_[1, cum_items][n_steps])
    indexer = (
        ((-1 if from_end else 0),) * (len(shape) - 1 - n_steps)
        + ((slice(-stop, None) if from_end else slice(stop)),)
        + (slice(None),) * n_steps
    )
    return indexer


def first_n_items(array, n_desired):
    """Returns the first n_desired items of an array"""
    # Unfortunately, we can't just do array.flat[:n_desired] here because it
    # might not be a numpy.ndarray. Moreover, access to elements of the array
    # could be very expensive (e.g. if it's only available over DAP), so go out
    # of our way to get them in a single call to __getitem__ using only slices.
    from xarray.core.variable import Variable

    if n_desired < 1:
        raise ValueError("must request at least one item")

    if array.size == 0:
        # work around for https://github.com/numpy/numpy/issues/5195
        return []

    if n_desired < array.size:
        indexer = _get_indexer_at_least_n_items(array.shape, n_desired, from_end=False)
        if isinstance(array, ExplicitlyIndexed):
            indexer = BasicIndexer(indexer)
        array = array[indexer]

    # We pass variable objects in to handle indexing
    # with indexer above. It would not work with our
    # lazy indexing classes at the moment, so we cannot
    # pass Variable._data
    if isinstance(array, Variable):
        array = array._data
    return ravel(to_duck_array(array))[:n_desired]


def last_n_items(array, n_desired):
    """Returns the last n_desired items of an array"""
    # Unfortunately, we can't just do array.flat[-n_desired:] here because it
    # might not be a numpy.ndarray. Moreover, access to elements of the array
    # could be very expensive (e.g. if it's only available over DAP), so go out
    # of our way to get them in a single call to __getitem__ using only slices.
    from xarray.core.variable import Variable

    if (n_desired == 0) or (array.size == 0):
        return []

    if n_desired < array.size:
        indexer = _get_indexer_at_least_n_items(array.shape, n_desired, from_end=True)
        if isinstance(array, ExplicitlyIndexed):
            indexer = BasicIndexer(indexer)
        array = array[indexer]

    # We pass variable objects in to handle indexing
    # with indexer above. It would not work with our
    # lazy indexing classes at the moment, so we cannot
    # pass Variable._data
    if isinstance(array, Variable):
        array = array._data
    return ravel(to_duck_array(array))[-n_desired:]


def last_item(array):
    """Returns the last item of an array."""
    indexer = (slice(-1, None),) * array.ndim
    return ravel(to_duck_array(array[indexer]))


def calc_max_rows_first(max_rows: int) -> int:
    """Calculate the first rows to maintain the max number of rows."""
    return max_rows // 2 + max_rows % 2


def calc_max_rows_last(max_rows: int) -> int:
    """Calculate the last rows to maintain the max number of rows."""
    return max_rows // 2


def format_timestamp(t):
    """Cast given object to a Timestamp and return a nicely formatted string"""
    try:
        timestamp = pd.Timestamp(t)
        datetime_str = timestamp.isoformat(sep=" ")
    except OutOfBoundsDatetime:
        datetime_str = str(t)

    try:
        date_str, time_str = datetime_str.split()
    except ValueError:
        # catch NaT and others that don't split nicely
        return datetime_str
    else:
        if time_str == "00:00:00":
            return date_str
        else:
            return f"{date_str}T{time_str}"


def format_timedelta(t, timedelta_format=None):
    """Cast given object to a Timestamp and return a nicely formatted string"""
    timedelta_str = str(pd.Timedelta(t))
    try:
        days_str, time_str = timedelta_str.split(" days ")
    except ValueError:
        # catch NaT and others that don't split nicely
        return timedelta_str
    else:
        if timedelta_format == "date":
            return days_str + " days"
        elif timedelta_format == "time":
            return time_str
        else:
            return timedelta_str


def format_item(x, timedelta_format=None, quote_strings=True):
    """Returns a succinct summary of an object as a string"""
    if isinstance(x, PandasExtensionArray):
        # We want to bypass PandasExtensionArray's repr here
        # because its __repr__ is PandasExtensionArray(array=[...])
        # and this function is only for single elements.
        return str(x.array[0])
    if isinstance(x, np.datetime64 | datetime):
        return format_timestamp(x)
    if isinstance(x, np.timedelta64 | timedelta):
        return format_timedelta(x, timedelta_format=timedelta_format)
    elif isinstance(x, str | bytes):
        if hasattr(x, "dtype"):
            x = x.item()
        return repr(x) if quote_strings else x
    elif hasattr(x, "dtype") and np.issubdtype(x.dtype, np.floating) and x.shape == ():
        return f"{x.item():.4}"
    else:
        return str(x)


def format_items(x):
    """Returns a succinct summaries of all items in a sequence as strings"""
    x = to_duck_array(x)
    timedelta_format = "datetime"
    if not isinstance(x, PandasExtensionArray) and np.issubdtype(
        x.dtype, np.timedelta64
    ):
        x = astype(x, dtype="timedelta64[ns]")
        day_part = x[~pd.isnull(x)].astype("timedelta64[D]").astype("timedelta64[ns]")
        time_needed = x[~pd.isnull(x)] != day_part
        day_needed = day_part != np.timedelta64(0, "ns")
        if array_all(np.logical_not(day_needed)):
            timedelta_format = "time"
        elif array_all(np.logical_not(time_needed)):
            timedelta_format = "date"

    formatted = [format_item(xi, timedelta_format) for xi in x]
    return formatted


def format_array_flat(array, max_width: int):
    """Return a formatted string for as many items in the flattened version of
    array that will fit within max_width characters.
    """
    # every item will take up at least two characters, but we always want to
    # print at least first and last items
    max_possibly_relevant = min(max(array.size, 1), max(math.ceil(max_width / 2.0), 2))
    relevant_front_items = format_items(
        first_n_items(array, (max_possibly_relevant + 1) // 2)
    )
    relevant_back_items = format_items(last_n_items(array, max_possibly_relevant // 2))
    # interleave relevant front and back items:
    #     [a, b, c] and [y, z] -> [a, z, b, y, c]
    relevant_items = sum(
        zip_longest(relevant_front_items, reversed(relevant_back_items)), ()
    )[:max_possibly_relevant]

    cum_len = np.cumsum([len(s) + 1 for s in relevant_items]) - 1
    if (array.size > 2) and (
        (max_possibly_relevant < array.size) or array_any(cum_len > max_width)
    ):
        padding = " ... "
        max_len = max(int(np.argmax(cum_len + len(padding) - 1 > max_width)), 2)
        count = min(array.size, max_len)
    else:
        count = array.size
        padding = "" if (count <= 1) else " "

    num_front = (count + 1) // 2
    num_back = count - num_front
    # note that num_back is 0 <--> array.size is 0 or 1
    #                         <--> relevant_back_items is []
    pprint_str = "".join(
        [
            " ".join(relevant_front_items[:num_front]),
            padding,
            " ".join(relevant_back_items[-num_back:]),
        ]
    )

    # As a final check, if it's still too long even with the limit in values,
    # replace the end with an ellipsis
    # NB: this will still returns a full 3-character ellipsis when max_width < 3
    if len(pprint_str) > max_width:
        pprint_str = pprint_str[: max(max_width - 3, 0)] + "..."

    return pprint_str


# mapping of tuple[modulename, classname] to repr
_KNOWN_TYPE_REPRS = {
    ("numpy", "ndarray"): "np.ndarray",
    ("sparse._coo.core", "COO"): "sparse.COO",
}


def inline_dask_repr(array):
    """Similar to dask.array.DataArray.__repr__, but without
    redundant information that's already printed by the repr
    function of the xarray wrapper.
    """
    assert isinstance(array, array_type("dask")), array

    chunksize = tuple(c[0] for c in array.chunks)

    if hasattr(array, "_meta"):
        meta = array._meta
        identifier = (type(meta).__module__, type(meta).__name__)
        meta_repr = _KNOWN_TYPE_REPRS.get(identifier, ".".join(identifier))
        meta_string = f", meta={meta_repr}"
    else:
        meta_string = ""

    return f"dask.array<chunksize={chunksize}{meta_string}>"


def inline_sparse_repr(array):
    """Similar to sparse.COO.__repr__, but without the redundant shape/dtype."""
    sparse_array_type = array_type("sparse")
    assert isinstance(array, sparse_array_type), array
    return f"<{type(array).__name__}: nnz={array.nnz:d}, fill_value={array.fill_value}>"


def inline_variable_array_repr(var, max_width):
    """Build a one-line summary of a variable's data."""
    if hasattr(var._data, "_repr_inline_"):
        return var._data._repr_inline_(max_width)
    if getattr(var, "_in_memory", False):
        return format_array_flat(var, max_width)
    dask_array_type = array_type("dask")
    if isinstance(var._data, dask_array_type):
        return inline_dask_repr(var.data)
    sparse_array_type = array_type("sparse")
    if isinstance(var._data, sparse_array_type):
        return inline_sparse_repr(var.data)
    if hasattr(var._data, "__array_function__"):
        return maybe_truncate(repr(var._data).replace("\n", " "), max_width)
    # internal xarray array type
    return "..."


def summarize_variable(
    name: Hashable,
    var: Variable,
    col_width: int | None = None,
    max_width: int | None = None,
    is_index: bool = False,
):
    """Summarize a variable in one line, e.g., for the Dataset.__repr__."""
    variable = getattr(var, "variable", var)

    if max_width is None:
        max_width_options = OPTIONS["display_width"]
        if not isinstance(max_width_options, int):
            raise TypeError(f"`max_width` value of `{max_width}` is not a valid int")
        else:
            max_width = max_width_options

    marker = "*" if is_index else " "
    first_col = f"  {marker} {name} "
    if col_width is not None:
        first_col = pretty_print(first_col, col_width)

    if variable.dims:
        dims_str = ", ".join(map(str, variable.dims))
        dims_str = f"({dims_str}) "
    else:
        dims_str = ""

    front_str = f"{first_col}{dims_str}{variable.dtype} {render_human_readable_nbytes(variable.nbytes)} "

    values_width = max_width - len(front_str)
    values_str = inline_variable_array_repr(variable, values_width)

    return f"{front_str}{values_str}"


def summarize_attr(key, value, col_width=None):
    """Summary for __repr__ - use ``X.attrs[key]`` for full value."""
    # Indent key and add ':', then right-pad if col_width is not None
    k_str = f"    {key}:"
    if col_width is not None:
        k_str = pretty_print(k_str, col_width)
    # Replace tabs and newlines, so we print on one line in known width
    v_str = str(value).replace("\t", "\\t").replace("\n", "\\n")
    # Finally, truncate to the desired display width
    return maybe_truncate(f"{k_str} {v_str}", OPTIONS["display_width"])


EMPTY_REPR = "    *empty*"


def _calculate_col_width(col_items):
    max_name_length = max((len(str(s)) for s in col_items), default=0)
    col_width = max(max_name_length, 7) + 6
    return col_width


def _mapping_repr(
    mapping,
    title,
    summarizer,
    expand_option_name,
    col_width=None,
    max_rows=None,
    indexes=None,
):
    if col_width is None:
        col_width = _calculate_col_width(mapping)

    summarizer_kwargs = defaultdict(dict)
    if indexes is not None:
        summarizer_kwargs = {k: {"is_index": k in indexes} for k in mapping}

    summary = [f"{title}:"]
    if mapping:
        len_mapping = len(mapping)
        if not _get_boolean_with_default(expand_option_name, default=True):
            summary = [f"{summary[0]} ({len_mapping})"]
        elif max_rows is not None and len_mapping > max_rows:
            summary = [f"{summary[0]} ({max_rows}/{len_mapping})"]
            first_rows = calc_max_rows_first(max_rows)
            keys = list(mapping.keys())
            summary += [
                summarizer(k, mapping[k], col_width, **summarizer_kwargs[k])
                for k in keys[:first_rows]
            ]
            if max_rows > 1:
                last_rows = calc_max_rows_last(max_rows)
                summary += [pretty_print("    ...", col_width) + " ..."]
                summary += [
                    summarizer(k, mapping[k], col_width, **summarizer_kwargs[k])
                    for k in keys[-last_rows:]
                ]
        else:
            summary += [
                summarizer(k, v, col_width, **summarizer_kwargs[k])
                for k, v in mapping.items()
            ]
    else:
        summary += [EMPTY_REPR]
    return "\n".join(summary)


data_vars_repr = functools.partial(
    _mapping_repr,
    title="Data variables",
    summarizer=summarize_variable,
    expand_option_name="display_expand_data_vars",
)

attrs_repr = functools.partial(
    _mapping_repr,
    title="Attributes",
    summarizer=summarize_attr,
    expand_option_name="display_expand_attrs",
)


def _coord_sort_key(coord, dims):
    """Sort key for coordinate ordering.

        Orders by:
        1. Primary: index of the matching dimension in dataset dims.
        2. Secondary: dimension coordinates (name == dim) come before non-dimension coordinates

        This groups non-dimension coordinates right after their associated dimension
    coordinate.
    """
    name, var = coord

    # Dimension coordinates come first within their dim section
    if name in dims:
        return (dims.index(name), 0)

    # Non-dimension coordinates come second within their dim section
    # Check the var.dims list in backwards order to put (x, y) after (x) and (y)
    for d in var.dims[::-1]:
        if d in dims:
            return (dims.index(d), 1)

    # Scalar coords or coords with dims not in dataset dims go at the end
    return (len(dims), 1)


def coords_repr(coords: AbstractCoordinates, col_width=None, max_rows=None):
    if col_width is None:
        col_width = _calculate_col_width(coords)
    dims = tuple(coords._data.dims)
    dim_ordered_coords = sorted(
        coords.items(), key=functools.partial(_coord_sort_key, dims=dims)
    )
    return _mapping_repr(
        dict(dim_ordered_coords),
        title="Coordinates",
        summarizer=summarize_variable,
        expand_option_name="display_expand_coords",
        col_width=col_width,
        indexes=coords.xindexes,
        max_rows=max_rows,
    )


def inherited_coords_repr(node: DataTree, col_width=None, max_rows=None):
    coords = inherited_vars(node._coord_variables)
    if col_width is None:
        col_width = _calculate_col_width(coords)
    return _mapping_repr(
        coords,
        title="Inherited coordinates",
        summarizer=summarize_variable,
        expand_option_name="display_expand_coords",
        col_width=col_width,
        indexes=node._indexes,
        max_rows=max_rows,
    )


def inline_index_repr(index: pd.Index, max_width: int) -> str:
    if hasattr(index, "_repr_inline_"):
        repr_ = index._repr_inline_(max_width=max_width)
    else:
        # fallback for the `pandas.Index` subclasses from
        # `Indexes.get_pandas_indexes` / `xr_obj.indexes`
        repr_ = repr(index)

    return repr_


def summarize_index(
    names: tuple[Hashable, ...],
    index,
    col_width: int,
    max_width: int | None = None,
) -> str:
    if max_width is None:
        max_width = OPTIONS["display_width"]

    def prefixes(length: int) -> list[str]:
        if length in (0, 1):
            return [" "]

        return ["┌"] + ["│"] * max(length - 2, 0) + ["└"]

    preformatted = [
        pretty_print(f"  {prefix} {name}", col_width)
        for prefix, name in zip(prefixes(len(names)), names, strict=True)
    ]

    head, *tail = preformatted
    index_width = max_width - len(head)
    repr_ = inline_index_repr(index, max_width=index_width)
    return "\n".join([head + repr_] + [line.rstrip() for line in tail])


def filter_nondefault_indexes(indexes, filter_indexes: bool):
    from xarray.core.indexes import PandasIndex, PandasMultiIndex

    if not filter_indexes:
        return indexes

    default_indexes = (PandasIndex, PandasMultiIndex)

    return {
        key: index
        for key, index in indexes.items()
        if not isinstance(index, default_indexes)
    }


def indexes_repr(indexes, max_rows: int | None = None, title: str = "Indexes") -> str:
    col_width = _calculate_col_width(chain.from_iterable(indexes))

    return _mapping_repr(
        indexes,
        title,
        summarize_index,
        "display_expand_indexes",
        col_width=col_width,
        max_rows=max_rows,
    )


def dim_summary(obj):
    elements = [f"{k}: {v}" for k, v in obj.sizes.items()]
    return ", ".join(elements)


def _element_formatter(
    elements: Collection[Hashable],
    col_width: int,
    max_rows: int | None = None,
    delimiter: str = ", ",
) -> str:
    """
    Formats elements for better readability.

    Once it becomes wider than the display width it will create a newline and
    continue indented to col_width.
    Once there are more rows than the maximum displayed rows it will start
    removing rows.

    Parameters
    ----------
    elements : Collection of hashable
        Elements to join together.
    col_width : int
        The width to indent to if a newline has been made.
    max_rows : int, optional
        The maximum number of allowed rows. The default is None.
    delimiter : str, optional
        Delimiter to use between each element. The default is ", ".
    """
    elements_len = len(elements)
    out = [""]
    length_row = 0
    for i, v in enumerate(elements):
        delim = delimiter if i < elements_len - 1 else ""
        v_delim = f"{v}{delim}"
        length_element = len(v_delim)
        length_row += length_element

        # Create a new row if the next elements makes the print wider than
        # the maximum display width:
        if col_width + length_row > OPTIONS["display_width"]:
            out[-1] = out[-1].rstrip()  # Remove trailing whitespace.
            out.append("\n" + pretty_print("", col_width) + v_delim)
            length_row = length_element
        else:
            out[-1] += v_delim

    # If there are too many rows of dimensions trim some away:
    if max_rows and (len(out) > max_rows):
        first_rows = calc_max_rows_first(max_rows)
        last_rows = calc_max_rows_last(max_rows)
        out = (
            out[:first_rows]
            + ["\n" + pretty_print("", col_width) + "..."]
            + (out[-last_rows:] if max_rows > 1 else [])
        )
    return "".join(out)


def dim_summary_limited(
    sizes: Mapping[Any, int], col_width: int, max_rows: int | None = None
) -> str:
    elements = [f"{k}: {v}" for k, v in sizes.items()]
    return _element_formatter(elements, col_width, max_rows)


def unindexed_dims_repr(dims, coords, max_rows: int | None = None):
    unindexed_dims = [d for d in dims if d not in coords]
    if unindexed_dims:
        dims_start = "Dimensions without coordinates: "
        dims_str = _element_formatter(
            unindexed_dims, col_width=len(dims_start), max_rows=max_rows
        )
        return dims_start + dims_str
    else:
        return None


@contextlib.contextmanager
def set_numpy_options(*args, **kwargs):
    original = np.get_printoptions()
    np.set_printoptions(*args, **kwargs)
    try:
        yield
    finally:
        np.set_printoptions(**original)


def limit_lines(string: str, *, limit: int):
    """
    If the string is more lines than the limit,
    this returns the middle lines replaced by an ellipsis
    """
    lines = string.splitlines()
    if len(lines) > limit:
        string = "\n".join(chain(lines[: limit // 2], ["..."], lines[-limit // 2 :]))
    return string


def short_array_repr(array):
    from xarray.core.common import AbstractArray

    if isinstance(array, AbstractArray):
        array = array.data
    if isinstance(array, pd.api.extensions.ExtensionArray):
        return repr(array)
    array = to_duck_array(array)

    # default to lower precision so a full (abbreviated) line can fit on
    # one line with the default display_width
    options = {
        "precision": 6,
        "linewidth": OPTIONS["display_width"],
        "threshold": OPTIONS["display_values_threshold"],
    }
    if array.ndim < 3:
        edgeitems = 3
    elif array.ndim == 3:
        edgeitems = 2
    else:
        edgeitems = 1
    options["edgeitems"] = edgeitems
    with set_numpy_options(**options):
        return repr(array)


def short_data_repr(array):
    """Format "data" for DataArray and Variable."""
    internal_data = getattr(array, "variable", array)._data

    if isinstance(array, np.ndarray):
        return short_array_repr(array)
    elif is_duck_array(internal_data):
        return limit_lines(repr(array.data), limit=40)
    elif getattr(array, "_in_memory", None):
        return short_array_repr(array)
    else:
        # internal xarray array type
        return f"[{array.size} values with dtype={array.dtype}]"


def _get_indexes_dict(indexes):
    return {
        tuple(index_vars.keys()): idx for idx, index_vars in indexes.group_by_index()
    }


@recursive_repr("<recursive array>")
def array_repr(arr):
    from xarray.core.variable import Variable

    max_rows = OPTIONS["display_max_rows"]

    # used for DataArray, Variable and IndexVariable
    if hasattr(arr, "name") and arr.name is not None:
        name_str = f"{arr.name!r} "
    else:
        name_str = ""

    if (
        isinstance(arr, Variable)
        or _get_boolean_with_default("display_expand_data", default=True)
        or isinstance(arr.variable._data, MemoryCachedArray)
    ):
        data_repr = short_data_repr(arr)
    else:
        data_repr = inline_variable_array_repr(arr.variable, OPTIONS["display_width"])

    start = f"<xarray.{type(arr).__name__} {name_str}"
    dims = dim_summary_limited(arr.sizes, col_width=len(start) + 1, max_rows=max_rows)
    nbytes_str = render_human_readable_nbytes(arr.nbytes)
    summary = [
        f"{start}({dims})> Size: {nbytes_str}",
        data_repr,
    ]
    if hasattr(arr, "coords"):
        if arr.coords:
            col_width = _calculate_col_width(arr.coords)
            summary.append(
                coords_repr(arr.coords, col_width=col_width, max_rows=max_rows)
            )

        unindexed_dims_str = unindexed_dims_repr(
            arr.dims, arr.coords, max_rows=max_rows
        )
        if unindexed_dims_str:
            summary.append(unindexed_dims_str)

        display_default_indexes = _get_boolean_with_default(
            "display_default_indexes", False
        )

        xindexes = filter_nondefault_indexes(
            _get_indexes_dict(arr.xindexes), not display_default_indexes
        )

        if xindexes:
            summary.append(indexes_repr(xindexes, max_rows=max_rows))

    if arr.attrs:
        summary.append(attrs_repr(arr.attrs, max_rows=max_rows))

    return "\n".join(summary)


@recursive_repr("<recursive Dataset>")
def dataset_repr(ds):
    nbytes_str = render_human_readable_nbytes(ds.nbytes)
    summary = [f"<xarray.{type(ds).__name__}> Size: {nbytes_str}"]

    col_width = _calculate_col_width(ds.variables)
    max_rows = OPTIONS["display_max_rows"]

    dims_start = pretty_print("Dimensions:", col_width)
    dims_values = dim_summary_limited(
        ds.sizes, col_width=col_width + 1, max_rows=max_rows
    )
    summary.append(f"{dims_start}({dims_values})")

    if ds.coords:
        summary.append(coords_repr(ds.coords, col_width=col_width, max_rows=max_rows))

    unindexed_dims_str = unindexed_dims_repr(ds.dims, ds.coords, max_rows=max_rows)
    if unindexed_dims_str:
        summary.append(unindexed_dims_str)

    summary.append(data_vars_repr(ds.data_vars, col_width=col_width, max_rows=max_rows))

    display_default_indexes = _get_boolean_with_default(
        "display_default_indexes", False
    )
    xindexes = filter_nondefault_indexes(
        _get_indexes_dict(ds.xindexes), not display_default_indexes
    )
    if xindexes:
        summary.append(indexes_repr(xindexes, max_rows=max_rows))

    if ds.attrs:
        summary.append(attrs_repr(ds.attrs, max_rows=max_rows))

    return "\n".join(summary)


def dims_and_coords_repr(ds) -> str:
    """Partial Dataset repr for use inside DataTree inheritance errors."""
    summary = []

    col_width = _calculate_col_width(ds.coords)
    max_rows = OPTIONS["display_max_rows"]

    dims_start = pretty_print("Dimensions:", col_width)
    dims_values = dim_summary_limited(
        ds.sizes, col_width=col_width + 1, max_rows=max_rows
    )
    summary.append(f"{dims_start}({dims_values})")

    if ds.coords:
        summary.append(coords_repr(ds.coords, col_width=col_width, max_rows=max_rows))

    unindexed_dims_str = unindexed_dims_repr(ds.dims, ds.coords, max_rows=max_rows)
    if unindexed_dims_str:
        summary.append(unindexed_dims_str)

    return "\n".join(summary)


def diff_name_summary(a, b) -> str:
    if a.name != b.name:
        return f"Differing names:\n    {a.name!r} != {b.name!r}"
    else:
        return ""


def diff_dim_summary(a, b) -> str:
    if a.sizes != b.sizes:
        return f"Differing dimensions:\n    ({dim_summary(a)}) != ({dim_summary(b)})"
    else:
        return ""


def _diff_mapping_repr(
    a_mapping,
    b_mapping,
    compat,
    title,
    summarizer,
    col_width=None,
    a_indexes=None,
    b_indexes=None,
):
    def compare_attr(a, b):
        if is_duck_array(a) or is_duck_array(b):
            return array_equiv(a, b)
        else:
            return a == b

    def extra_items_repr(extra_keys, mapping, ab_side, kwargs):
        extra_repr = [
            summarizer(k, mapping[k], col_width, **kwargs[k]) for k in extra_keys
        ]
        if extra_repr:
            header = f"{title} only on the {ab_side} object:"
            return [header] + extra_repr
        else:
            return []

    a_keys = set(a_mapping)
    b_keys = set(b_mapping)

    summary = []

    diff_items = []

    a_summarizer_kwargs = defaultdict(dict)
    if a_indexes is not None:
        a_summarizer_kwargs = {k: {"is_index": k in a_indexes} for k in a_mapping}
    b_summarizer_kwargs = defaultdict(dict)
    if b_indexes is not None:
        b_summarizer_kwargs = {k: {"is_index": k in b_indexes} for k in b_mapping}

    for k in a_keys & b_keys:
        try:
            # compare xarray variable
            if not callable(compat):
                compatible = getattr(a_mapping[k].variable, compat)(
                    b_mapping[k].variable
                )
            else:
                compatible = compat(a_mapping[k].variable, b_mapping[k].variable)
            is_variable = True
        except AttributeError:
            # compare attribute value
            compatible = compare_attr(a_mapping[k], b_mapping[k])
            is_variable = False

        if not compatible:
            temp = [
                summarizer(k, a_mapping[k], col_width, **a_summarizer_kwargs[k]),
                summarizer(k, b_mapping[k], col_width, **b_summarizer_kwargs[k]),
            ]

            if compat == "identical" and is_variable:
                attrs_summary = []
                a_attrs = a_mapping[k].attrs
                b_attrs = b_mapping[k].attrs

                attrs_to_print = set(a_attrs) ^ set(b_attrs)
                attrs_to_print.update(
                    {
                        k
                        for k in set(a_attrs) & set(b_attrs)
                        if not compare_attr(a_attrs[k], b_attrs[k])
                    }
                )
                for m in (a_mapping, b_mapping):
                    attr_s = "\n".join(
                        "    " + summarize_attr(ak, av)
                        for ak, av in m[k].attrs.items()
                        if ak in attrs_to_print
                    )
                    if attr_s:
                        attr_s = "    Differing variable attributes:\n" + attr_s
                    attrs_summary.append(attr_s)

                temp = [
                    f"{var_s}\n{attr_s}" if attr_s else var_s
                    for var_s, attr_s in zip(temp, attrs_summary, strict=True)
                ]

                # TODO: It should be possible recursively use _diff_mapping_repr
                #       instead of explicitly handling variable attrs specially.
                #       That would require some refactoring.
                # newdiff = _diff_mapping_repr(
                #     {k: v for k,v in a_attrs.items() if k in attrs_to_print},
                #     {k: v for k,v in b_attrs.items() if k in attrs_to_print},
                #     compat=compat,
                #     summarizer=summarize_attr,
                #     title="Variable Attributes"
                # )
                # temp += [newdiff]

            diff_items += [
                ab_side + s[1:] for ab_side, s in zip(("L", "R"), temp, strict=True)
            ]

    if diff_items:
        summary += [f"Differing {title.lower()}:"] + diff_items

    summary += extra_items_repr(a_keys - b_keys, a_mapping, "left", a_summarizer_kwargs)
    summary += extra_items_repr(
        b_keys - a_keys, b_mapping, "right", b_summarizer_kwargs
    )

    return "\n".join(summary)


def diff_coords_repr(a, b, compat, col_width=None):
    return _diff_mapping_repr(
        a,
        b,
        compat,
        "Coordinates",
        summarize_variable,
        col_width=col_width,
        a_indexes=a.xindexes,
        b_indexes=b.xindexes,
    )


diff_data_vars_repr = functools.partial(
    _diff_mapping_repr, title="Data variables", summarizer=summarize_variable
)


diff_attrs_repr = functools.partial(
    _diff_mapping_repr, title="Attributes", summarizer=summarize_attr
)


def _compat_to_str(compat):
    if callable(compat):
        compat = compat.__name__

    if compat == "equals":
        return "equal"
    elif compat == "allclose":
        return "close"
    else:
        return compat


def diff_indexes_repr(a_indexes, b_indexes, col_width: int = 20) -> str:
    """Generate diff representation for indexes."""
    a_keys = set(a_indexes.keys())
    b_keys = set(b_indexes.keys())

    summary = []

    if only_a := a_keys - b_keys:
        summary.append(f"Indexes only on the left object:  {sorted(only_a)}")

    if only_b := b_keys - a_keys:
        summary.append(f"Indexes only on the right object: {sorted(only_b)}")

    # Check for indexes on the same coordinates but with different types or values
    common_keys = a_keys & b_keys
    diff_items = []

    for key in sorted(common_keys):
        a_idx = a_indexes[key]
        b_idx = b_indexes[key]

        # Check if indexes differ
        indexes_equal = False
        if type(a_idx) is type(b_idx):
            try:
                indexes_equal = a_idx.equals(b_idx)
            except NotImplementedError:
                # Fall back to variable comparison
                a_var = a_indexes.variables[key]
                b_var = b_indexes.variables[key]
                indexes_equal = a_var.equals(b_var)

        if not indexes_equal:
            # Format the index values similar to variable diff
            try:
                a_repr = inline_index_repr(
                    a_indexes.to_pandas_indexes()[key], max_width=70
                )
                b_repr = inline_index_repr(
                    b_indexes.to_pandas_indexes()[key], max_width=70
                )
            except TypeError:
                # Custom indexes may not support to_pandas_index()
                a_repr = repr(a_idx)
                b_repr = repr(b_idx)
            diff_items.append(f"L   {key!s:<{col_width}} {a_repr}")
            diff_items.append(f"R   {key!s:<{col_width}} {b_repr}")

    if diff_items:
        summary.append("Differing indexes:\n" + "\n".join(diff_items))

    return "\n".join(summary)


def diff_array_repr(a, b, compat):
    # used for DataArray, Variable and IndexVariable
    summary = [
        f"Left and right {type(a).__name__} objects are not {_compat_to_str(compat)}"
    ]

    if dims_diff := diff_dim_summary(a, b):
        summary.append(dims_diff)

    if callable(compat):
        equiv = compat
    else:
        equiv = array_equiv

    if not equiv(a.data, b.data):
        temp = [wrap_indent(short_array_repr(obj), start="    ") for obj in (a, b)]
        diff_data_repr = [
            ab_side + "\n" + ab_data_repr
            for ab_side, ab_data_repr in zip(("L", "R"), temp, strict=True)
        ]
        summary += ["Differing values:"] + diff_data_repr

    if hasattr(a, "coords"):
        col_width = _calculate_col_width(set(a.coords) | set(b.coords))
        if coords_diff := diff_coords_repr(
            a.coords, b.coords, compat, col_width=col_width
        ):
            summary.append(coords_diff)

    if compat == "identical":
        if hasattr(a, "xindexes") and (
            indexes_diff := diff_indexes_repr(a.xindexes, b.xindexes)
        ):
            summary.append(indexes_diff)

        if attrs_diff := diff_attrs_repr(a.attrs, b.attrs, compat):
            summary.append(attrs_diff)

    return "\n".join(summary)


def diff_treestructure(a: DataTree, b: DataTree) -> str | None:
    """
    Return a summary of why two trees are not isomorphic.
    If they are isomorphic return None.
    """
    # .group_subtrees walks nodes in breadth-first-order, in order to produce as
    # shallow of a diff as possible
    for path, (node_a, node_b) in group_subtrees(a, b):
        if node_a.children.keys() != node_b.children.keys():
            path_str = "root node" if path == "." else f"node {path!r}"
            child_summary = f"{list(node_a.children)} vs {list(node_b.children)}"
            diff = f"Children at {path_str} do not match: {child_summary}"
            return diff

    return None


def diff_dataset_repr(a, b, compat):
    summary = [
        f"Left and right {type(a).__name__} objects are not {_compat_to_str(compat)}"
    ]

    col_width = _calculate_col_width(set(list(a.variables) + list(b.variables)))

    if dims_diff := diff_dim_summary(a, b):
        summary.append(dims_diff)
    if coords_diff := diff_coords_repr(a.coords, b.coords, compat, col_width=col_width):
        summary.append(coords_diff)
    if data_diff := diff_data_vars_repr(
        a.data_vars, b.data_vars, compat, col_width=col_width
    ):
        summary.append(data_diff)

    if compat == "identical":
        if indexes_diff := diff_indexes_repr(a.xindexes, b.xindexes):
            summary.append(indexes_diff)

        if attrs_diff := diff_attrs_repr(a.attrs, b.attrs, compat):
            summary.append(attrs_diff)

    return "\n".join(summary)


def diff_nodewise_summary(a: DataTree, b: DataTree, compat):
    """Iterates over all corresponding nodes, recording differences between data at each location."""

    summary = []
    for path, (node_a, node_b) in group_subtrees(a, b):
        a_ds, b_ds = node_a.dataset, node_b.dataset

        if not a_ds._all_compat(b_ds, compat):
            path_str = "root node" if path == "." else f"node {path!r}"
            dataset_diff = diff_dataset_repr(a_ds, b_ds, compat)
            data_diff = indent(
                "\n".join(dataset_diff.split("\n", 1)[1:]), prefix="    "
            )
            nodediff = f"Data at {path_str} does not match:\n{data_diff}"
            summary.append(nodediff)

    return "\n\n".join(summary)


def diff_datatree_repr(a: DataTree, b: DataTree, compat):
    summary = [
        f"Left and right {type(a).__name__} objects are not {_compat_to_str(compat)}"
    ]

    if compat == "identical" and (diff_name := diff_name_summary(a, b)):
        summary.append(diff_name)

    treestructure_diff = diff_treestructure(a, b)

    # If the trees structures are different there is no point comparing each node,
    # and doing so would raise an error.
    # TODO we could show any differences in nodes up to the first place that structure differs?
    if treestructure_diff is not None:
        summary.append(treestructure_diff)
    elif compat != "isomorphic":
        nodewise_diff = diff_nodewise_summary(a, b, compat)
        summary.append(nodewise_diff)

    return "\n\n".join(summary)


def inherited_vars(mapping: ChainMap) -> dict:
    return {k: v for k, v in mapping.parents.items() if k not in mapping.maps[0]}


def _datatree_node_repr(node: DataTree, root: bool) -> str:
    summary = [f"Group: {node.path}"]

    col_width = _calculate_col_width(node.variables)
    max_rows = OPTIONS["display_max_rows"]

    inherited_coords = inherited_vars(node._coord_variables)

    # Only show dimensions if also showing a variable or coordinates section.
    show_dims = (
        node._node_coord_variables
        or (root and inherited_coords)
        or node._data_variables
    )

    dim_sizes = node.sizes if root else node._node_dims

    if show_dims:
        # Includes inherited dimensions.
        dims_start = pretty_print("Dimensions:", col_width)
        dims_values = dim_summary_limited(
            dim_sizes, col_width=col_width + 1, max_rows=max_rows
        )
        summary.append(f"{dims_start}({dims_values})")

    if node._node_coord_variables:
        node_coords = node.to_dataset(inherit=False).coords
        summary.append(coords_repr(node_coords, col_width=col_width, max_rows=max_rows))

    if root and inherited_coords:
        summary.append(
            inherited_coords_repr(node, col_width=col_width, max_rows=max_rows)
        )

    if show_dims:
        unindexed_dims_str = unindexed_dims_repr(
            dim_sizes, node.coords, max_rows=max_rows
        )
        if unindexed_dims_str:
            summary.append(unindexed_dims_str)

    if node._data_variables:
        summary.append(
            data_vars_repr(node._data_variables, col_width=col_width, max_rows=max_rows)
        )

    # TODO: only show indexes defined at this node, with a separate section for
    # inherited indexes (if root=True)
    display_default_indexes = _get_boolean_with_default(
        "display_default_indexes", False
    )
    xindexes = filter_nondefault_indexes(
        _get_indexes_dict(node.xindexes), not display_default_indexes
    )
    if xindexes:
        summary.append(indexes_repr(xindexes, max_rows=max_rows))

    if node.attrs:
        summary.append(attrs_repr(node.attrs, max_rows=max_rows))

    return "\n".join(summary)


def datatree_repr(dt: DataTree) -> str:
    """A printable representation of the structure of this entire tree."""
    max_children = OPTIONS["display_max_children"]

    renderer = RenderDataTree(dt, maxchildren=max_children)

    name_info = "" if dt.name is None else f" {dt.name!r}"
    header = f"<xarray.DataTree{name_info}>"

    lines = [header]
    root = True

    for pre, fill, node in renderer:
        if isinstance(node, str):
            lines.append(f"{fill}{node}")
            continue

        node_repr = _datatree_node_repr(node, root=root)
        root = False  # only the first node is the root

        # TODO: figure out if we can restructure this logic to move child groups
        # up higher in the repr, directly below the <xarray.DataTree> header.
        # This would be more consistent with the HTML repr.
        raw_repr_lines = node_repr.splitlines()

        node_line = f"{pre}{raw_repr_lines[0]}"
        lines.append(node_line)

        for line in raw_repr_lines[1:]:
            if len(node.children) > 0:
                lines.append(f"{fill}{renderer.style.vertical}{line}")
            else:
                lines.append(f"{fill}{' ' * len(renderer.style.vertical)}{line}")

    return "\n".join(lines)


def shorten_list_repr(items: Sequence, max_items: int) -> str:
    if len(items) <= max_items:
        return repr(items)
    else:
        first_half = repr(items[: max_items // 2])[
            1:-1
        ]  # Convert to string and remove brackets
        second_half = repr(items[-max_items // 2 :])[
            1:-1
        ]  # Convert to string and remove brackets
        return f"[{first_half}, ..., {second_half}]"


def render_human_readable_nbytes(
    nbytes: int,
    /,
    *,
    attempt_constant_width: bool = False,
) -> str:
    """Renders simple human-readable byte count representation

    This is only a quick representation that should not be relied upon for precise needs.

    To get the exact byte count, please use the ``nbytes`` attribute directly.

    Parameters
    ----------
    nbytes
        Byte count
    attempt_constant_width
        For reasonable nbytes sizes, tries to render a fixed-width representation.

    Returns
    -------
        Human-readable representation of the byte count
    """
    dividend = float(nbytes)
    divisor = 1000.0
    last_unit_available = UNITS[-1]

    for unit in UNITS:
        if dividend < divisor or unit == last_unit_available:
            break
        dividend /= divisor

    dividend_str = f"{dividend:.0f}"
    unit_str = f"{unit}"

    if attempt_constant_width:
        dividend_str = dividend_str.rjust(3)
        unit_str = unit_str.ljust(2)

    string = f"{dividend_str}{unit_str}"
    return string
