from __future__ import annotations

import re
import textwrap

from xarray.core.dataset import Dataset
from xarray.core.datatree_mapping import map_over_subtree

"""
Module which specifies the subset of xarray.Dataset's API which we wish to copy onto DataTree.

Structured to mirror the way xarray defines Dataset's various operations internally, but does not actually import from
xarray's internals directly, only the public-facing xarray.Dataset class.
"""


_MAPPED_DOCSTRING_ADDENDUM = (
    "This method was copied from xarray.Dataset, but has been altered to "
    "call the method on the Datasets stored in every node of the subtree. "
    "See the `map_over_subtree` function for more details."
)

# TODO equals, broadcast_equals etc.
# TODO do dask-related private methods need to be exposed?
_DATASET_DASK_METHODS_TO_MAP = [
    "load",
    "compute",
    "persist",
    "unify_chunks",
    "chunk",
    "map_blocks",
]
_DATASET_METHODS_TO_MAP = [
    "as_numpy",
    "set_coords",
    "reset_coords",
    "info",
    "isel",
    "sel",
    "head",
    "tail",
    "thin",
    "broadcast_like",
    "reindex_like",
    "reindex",
    "interp",
    "interp_like",
    "rename",
    "rename_dims",
    "rename_vars",
    "swap_dims",
    "expand_dims",
    "set_index",
    "reset_index",
    "reorder_levels",
    "stack",
    "unstack",
    "merge",
    "drop_vars",
    "drop_sel",
    "drop_isel",
    "drop_dims",
    "transpose",
    "dropna",
    "fillna",
    "interpolate_na",
    "ffill",
    "bfill",
    "combine_first",
    "reduce",
    "map",
    "diff",
    "shift",
    "roll",
    "sortby",
    "quantile",
    "rank",
    "differentiate",
    "integrate",
    "cumulative_integrate",
    "filter_by_attrs",
    "polyfit",
    "pad",
    "idxmin",
    "idxmax",
    "argmin",
    "argmax",
    "query",
    "curvefit",
]
_ALL_DATASET_METHODS_TO_MAP = _DATASET_DASK_METHODS_TO_MAP + _DATASET_METHODS_TO_MAP

_DATA_WITH_COORDS_METHODS_TO_MAP = [
    "squeeze",
    "clip",
    "assign_coords",
    "where",
    "close",
    "isnull",
    "notnull",
    "isin",
    "astype",
]

REDUCE_METHODS = ["all", "any"]
NAN_REDUCE_METHODS = [
    "max",
    "min",
    "mean",
    "prod",
    "sum",
    "std",
    "var",
    "median",
]
NAN_CUM_METHODS = ["cumsum", "cumprod"]
_TYPED_DATASET_OPS_TO_MAP = [
    "__add__",
    "__sub__",
    "__mul__",
    "__pow__",
    "__truediv__",
    "__floordiv__",
    "__mod__",
    "__and__",
    "__xor__",
    "__or__",
    "__lt__",
    "__le__",
    "__gt__",
    "__ge__",
    "__eq__",
    "__ne__",
    "__radd__",
    "__rsub__",
    "__rmul__",
    "__rpow__",
    "__rtruediv__",
    "__rfloordiv__",
    "__rmod__",
    "__rand__",
    "__rxor__",
    "__ror__",
    "__iadd__",
    "__isub__",
    "__imul__",
    "__ipow__",
    "__itruediv__",
    "__ifloordiv__",
    "__imod__",
    "__iand__",
    "__ixor__",
    "__ior__",
    "__neg__",
    "__pos__",
    "__abs__",
    "__invert__",
    "round",
    "argsort",
    "conj",
    "conjugate",
]
# TODO NUM_BINARY_OPS apparently aren't defined on DatasetArithmetic, and don't appear to be injected anywhere...
_ARITHMETIC_METHODS_TO_MAP = (
    REDUCE_METHODS
    + NAN_REDUCE_METHODS
    + NAN_CUM_METHODS
    + _TYPED_DATASET_OPS_TO_MAP
    + ["__array_ufunc__"]
)


def _wrap_then_attach_to_cls(
    target_cls_dict, source_cls, methods_to_set, wrap_func=None
):
    """
    Attach given methods on a class, and optionally wrap each method first. (i.e. with map_over_subtree).

    Result is like having written this in the classes' definition:
    ```
    @wrap_func
    def method_name(self, *args, **kwargs):
        return self.method(*args, **kwargs)
    ```

    Every method attached here needs to have a return value of Dataset or DataArray in order to construct a new tree.

    Parameters
    ----------
    target_cls_dict : MappingProxy
        The __dict__ attribute of the class which we want the methods to be added to. (The __dict__ attribute can also
        be accessed by calling vars() from within that classes' definition.) This will be updated by this function.
    source_cls : class
        Class object from which we want to copy methods (and optionally wrap them). Should be the actual class object
        (or instance), not just the __dict__.
    methods_to_set : Iterable[Tuple[str, callable]]
        The method names and definitions supplied as a list of (method_name_string, method) pairs.
        This format matches the output of inspect.getmembers().
    wrap_func : callable, optional
        Function to decorate each method with. Must have the same return type as the method.
    """
    for method_name in methods_to_set:
        orig_method = getattr(source_cls, method_name)
        wrapped_method = (
            wrap_func(orig_method) if wrap_func is not None else orig_method
        )
        target_cls_dict[method_name] = wrapped_method

        if wrap_func is map_over_subtree:
            # Add a paragraph to the method's docstring explaining how it's been mapped
            orig_method_docstring = orig_method.__doc__

            if orig_method_docstring is not None:
                new_method_docstring = insert_doc_addendum(
                    orig_method_docstring, _MAPPED_DOCSTRING_ADDENDUM
                )
                setattr(target_cls_dict[method_name], "__doc__", new_method_docstring)


def insert_doc_addendum(docstring: str | None, addendum: str) -> str | None:
    """Insert addendum after first paragraph or at the end of the docstring.

    There are a number of Dataset's functions that are wrapped. These come from
    Dataset directly as well as the mixins: DataWithCoords, DatasetAggregations, and DatasetOpsMixin.

    The majority of the docstrings fall into a parseable pattern. Those that
    don't, just have the addendum appeneded after. None values are returned.

    """
    if docstring is None:
        return None

    pattern = re.compile(
        r"^(?P<start>(\S+)?(.*?))(?P<paragraph_break>\n\s*\n)(?P<whitespace>[ ]*)(?P<rest>.*)",
        re.DOTALL,
    )
    capture = re.match(pattern, docstring)
    if capture is None:
        ### single line docstring.
        return (
            docstring
            + "\n\n"
            + textwrap.fill(
                addendum,
                subsequent_indent="    ",
                width=79,
            )
        )

    if len(capture.groups()) == 6:
        return (
            capture["start"]
            + capture["paragraph_break"]
            + capture["whitespace"]
            + ".. note::\n"
            + textwrap.fill(
                addendum,
                initial_indent=capture["whitespace"] + "    ",
                subsequent_indent=capture["whitespace"] + "    ",
                width=79,
            )
            + capture["paragraph_break"]
            + capture["whitespace"]
            + capture["rest"]
        )
    else:
        return docstring


class MappedDatasetMethodsMixin:
    """
    Mixin to add methods defined specifically on the Dataset class such as .query(), but wrapped to map over all nodes
    in the subtree.
    """

    _wrap_then_attach_to_cls(
        target_cls_dict=vars(),
        source_cls=Dataset,
        methods_to_set=_ALL_DATASET_METHODS_TO_MAP,
        wrap_func=map_over_subtree,
    )


class MappedDataWithCoords:
    """
    Mixin to add coordinate-aware Dataset methods such as .where(), but wrapped to map over all nodes in the subtree.
    """

    # TODO add mapped versions of groupby, weighted, rolling, rolling_exp, coarsen, resample
    _wrap_then_attach_to_cls(
        target_cls_dict=vars(),
        source_cls=Dataset,
        methods_to_set=_DATA_WITH_COORDS_METHODS_TO_MAP,
        wrap_func=map_over_subtree,
    )


class DataTreeArithmeticMixin:
    """
    Mixin to add Dataset arithmetic operations such as __add__, reduction methods such as .mean(), and enable numpy
    ufuncs such as np.sin(), but wrapped to map over all nodes in the subtree.
    """

    _wrap_then_attach_to_cls(
        target_cls_dict=vars(),
        source_cls=Dataset,
        methods_to_set=_ARITHMETIC_METHODS_TO_MAP,
        wrap_func=map_over_subtree,
    )
