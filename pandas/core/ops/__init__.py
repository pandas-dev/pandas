"""
Arithmetic operations for PandasObjects

This is not a public API.
"""
import operator
from typing import TYPE_CHECKING, Optional, Set, Type

import numpy as np

from pandas._libs import lib
from pandas._libs.ops_dispatch import maybe_dispatch_ufunc_to_dunder_op  # noqa:F401
from pandas._typing import Level
from pandas.util._decorators import Appender

from pandas.core.dtypes.common import is_list_like
from pandas.core.dtypes.generic import ABCDataFrame, ABCIndexClass, ABCSeries
from pandas.core.dtypes.missing import isna

from pandas.core import algorithms
from pandas.core.construction import extract_array
from pandas.core.ops.array_ops import (
    arithmetic_op,
    comparison_op,
    get_array_op,
    logical_op,
)
from pandas.core.ops.array_ops import comp_method_OBJECT_ARRAY  # noqa:F401
from pandas.core.ops.common import unpack_zerodim_and_defer
from pandas.core.ops.docstrings import (
    _arith_doc_FRAME,
    _flex_comp_doc_FRAME,
    _make_flex_doc,
    _op_descriptions,
)
from pandas.core.ops.invalid import invalid_comparison  # noqa:F401
from pandas.core.ops.mask_ops import kleene_and, kleene_or, kleene_xor  # noqa: F401
from pandas.core.ops.methods import (  # noqa:F401
    add_flex_arithmetic_methods,
    add_special_arithmetic_methods,
)
from pandas.core.ops.roperator import (  # noqa:F401
    radd,
    rand_,
    rdiv,
    rdivmod,
    rfloordiv,
    rmod,
    rmul,
    ror_,
    rpow,
    rsub,
    rtruediv,
    rxor,
)

if TYPE_CHECKING:
    from pandas import DataFrame, Series  # noqa:F401

# -----------------------------------------------------------------------------
# constants
ARITHMETIC_BINOPS: Set[str] = {
    "add",
    "sub",
    "mul",
    "pow",
    "mod",
    "floordiv",
    "truediv",
    "divmod",
    "radd",
    "rsub",
    "rmul",
    "rpow",
    "rmod",
    "rfloordiv",
    "rtruediv",
    "rdivmod",
}


COMPARISON_BINOPS: Set[str] = {"eq", "ne", "lt", "gt", "le", "ge"}

# -----------------------------------------------------------------------------
# Ops Wrapping Utilities


def get_op_result_name(left, right):
    """
    Find the appropriate name to pin to an operation result.  This result
    should always be either an Index or a Series.

    Parameters
    ----------
    left : {Series, Index}
    right : object

    Returns
    -------
    name : object
        Usually a string
    """
    # `left` is always a Series when called from within ops
    if isinstance(right, (ABCSeries, ABCIndexClass)):
        name = _maybe_match_name(left, right)
    else:
        name = left.name
    return name


def _maybe_match_name(a, b):
    """
    Try to find a name to attach to the result of an operation between
    a and b.  If only one of these has a `name` attribute, return that
    name.  Otherwise return a consensus name if they match of None if
    they have different names.

    Parameters
    ----------
    a : object
    b : object

    Returns
    -------
    name : str or None

    See Also
    --------
    pandas.core.common.consensus_name_attr
    """
    a_has = hasattr(a, "name")
    b_has = hasattr(b, "name")
    if a_has and b_has:
        if a.name == b.name:
            return a.name
        else:
            # TODO: what if they both have np.nan for their names?
            return None
    elif a_has:
        return a.name
    elif b_has:
        return b.name
    return None


# -----------------------------------------------------------------------------


def _get_frame_op_default_axis(name: str) -> Optional[str]:
    """
    Only DataFrame cares about default_axis, specifically:
    special methods have default_axis=None and flex methods
    have default_axis='columns'.

    Parameters
    ----------
    name : str

    Returns
    -------
    default_axis: str or None
    """
    if name.replace("__r", "__") in ["__and__", "__or__", "__xor__"]:
        # bool methods
        return "columns"
    elif name.startswith("__"):
        # __add__, __mul__, ...
        return None
    else:
        # add, mul, ...
        return "columns"


def _get_op_name(op, special: bool) -> str:
    """
    Find the name to attach to this method according to conventions
    for special and non-special methods.

    Parameters
    ----------
    op : binary operator
    special : bool

    Returns
    -------
    op_name : str
    """
    opname = op.__name__.strip("_")
    if special:
        opname = f"__{opname}__"
    return opname


# -----------------------------------------------------------------------------
# Masking NA values and fallbacks for operations numpy does not support


def fill_binop(left, right, fill_value):
    """
    If a non-None fill_value is given, replace null entries in left and right
    with this value, but only in positions where _one_ of left/right is null,
    not both.

    Parameters
    ----------
    left : array-like
    right : array-like
    fill_value : object

    Returns
    -------
    left : array-like
    right : array-like

    Notes
    -----
    Makes copies if fill_value is not None and NAs are present.
    """
    if fill_value is not None:
        left_mask = isna(left)
        right_mask = isna(right)

        # one but not both
        mask = left_mask ^ right_mask

        if left_mask.any():
            # Avoid making a copy if we can
            left = left.copy()
            left[left_mask & mask] = fill_value

        if right_mask.any():
            # Avoid making a copy if we can
            right = right.copy()
            right[right_mask & mask] = fill_value

    return left, right


# -----------------------------------------------------------------------------
# Dispatch logic


def dispatch_to_series(left, right, func, axis: Optional[int] = None):
    """
    Evaluate the frame operation func(left, right) by evaluating
    column-by-column, dispatching to the Series implementation.

    Parameters
    ----------
    left : DataFrame
    right : scalar, Series, or DataFrame
    func : arithmetic or comparison operator
    axis : {None, 0, 1}

    Returns
    -------
    DataFrame
    """
    # Get the appropriate array-op to apply to each column/block's values.
    array_op = get_array_op(func)

    right = lib.item_from_zerodim(right)
    if not is_list_like(right):
        # i.e. scalar, faster than checking np.ndim(right) == 0
        bm = left._mgr.apply(array_op, right=right)
        return type(left)(bm)

    elif isinstance(right, ABCDataFrame):
        assert left.index.equals(right.index)
        assert left.columns.equals(right.columns)
        # TODO: The previous assertion `assert right._indexed_same(left)`
        #  fails in cases with empty columns reached via
        #  _frame_arith_method_with_reindex

        bm = left._mgr.operate_blockwise(right._mgr, array_op)
        return type(left)(bm)

    elif isinstance(right, ABCSeries) and axis == 1:
        # axis=1 means we want to operate row-by-row
        assert right.index.equals(left.columns)

        right = right._values
        # maybe_align_as_frame ensures we do not have an ndarray here
        assert not isinstance(right, np.ndarray)

        arrays = [array_op(l, r) for l, r in zip(left._iter_column_arrays(), right)]

    elif isinstance(right, ABCSeries):
        assert right.index.equals(left.index)  # Handle other cases later
        right = right._values

        arrays = [array_op(l, right) for l in left._iter_column_arrays()]

    else:
        # Remaining cases have less-obvious dispatch rules
        raise NotImplementedError(right)

    return type(left)._from_arrays(
        arrays, left.columns, left.index, verify_integrity=False
    )


# -----------------------------------------------------------------------------
# Series


def _align_method_SERIES(left: "Series", right, align_asobject: bool = False):
    """ align lhs and rhs Series """
    # ToDo: Different from _align_method_FRAME, list, tuple and ndarray
    # are not coerced here
    # because Series has inconsistencies described in #13637

    if isinstance(right, ABCSeries):
        # avoid repeated alignment
        if not left.index.equals(right.index):

            if align_asobject:
                # to keep original value's dtype for bool ops
                left = left.astype(object)
                right = right.astype(object)

            left, right = left.align(right, copy=False)

    return left, right


def _arith_method_SERIES(cls, op, special):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    assert special  # non-special uses _flex_method_SERIES
    op_name = _get_op_name(op, special)

    @unpack_zerodim_and_defer(op_name)
    def wrapper(left, right):

        left, right = _align_method_SERIES(left, right)
        res_name = get_op_result_name(left, right)

        lvalues = extract_array(left, extract_numpy=True)
        rvalues = extract_array(right, extract_numpy=True)
        result = arithmetic_op(lvalues, rvalues, op)

        return left._construct_result(result, name=res_name)

    wrapper.__name__ = op_name
    return wrapper


def _comp_method_SERIES(cls, op, special):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    assert special  # non-special uses _flex_method_SERIES
    op_name = _get_op_name(op, special)

    @unpack_zerodim_and_defer(op_name)
    def wrapper(self, other):

        res_name = get_op_result_name(self, other)

        if isinstance(other, ABCSeries) and not self._indexed_same(other):
            raise ValueError("Can only compare identically-labeled Series objects")

        lvalues = extract_array(self, extract_numpy=True)
        rvalues = extract_array(other, extract_numpy=True)

        res_values = comparison_op(lvalues, rvalues, op)

        return self._construct_result(res_values, name=res_name)

    wrapper.__name__ = op_name
    return wrapper


def _bool_method_SERIES(cls, op, special):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    assert special  # non-special uses _flex_method_SERIES
    op_name = _get_op_name(op, special)

    @unpack_zerodim_and_defer(op_name)
    def wrapper(self, other):
        self, other = _align_method_SERIES(self, other, align_asobject=True)
        res_name = get_op_result_name(self, other)

        lvalues = extract_array(self, extract_numpy=True)
        rvalues = extract_array(other, extract_numpy=True)

        res_values = logical_op(lvalues, rvalues, op)
        return self._construct_result(res_values, name=res_name)

    wrapper.__name__ = op_name
    return wrapper


def _flex_method_SERIES(cls, op, special):
    assert not special  # "special" also means "not flex"
    name = _get_op_name(op, special)
    doc = _make_flex_doc(name, "series")

    @Appender(doc)
    def flex_wrapper(self, other, level=None, fill_value=None, axis=0):
        # validate axis
        if axis is not None:
            self._get_axis_number(axis)

        if isinstance(other, ABCSeries):
            return self._binop(other, op, level=level, fill_value=fill_value)
        elif isinstance(other, (np.ndarray, list, tuple)):
            if len(other) != len(self):
                raise ValueError("Lengths must be equal")
            other = self._constructor(other, self.index)
            return self._binop(other, op, level=level, fill_value=fill_value)
        else:
            if fill_value is not None:
                self = self.fillna(fill_value)

            return op(self, other)

    flex_wrapper.__name__ = name
    return flex_wrapper


# -----------------------------------------------------------------------------
# DataFrame


def _align_method_FRAME(
    left, right, axis, flex: Optional[bool] = False, level: Level = None
):
    """
    Convert rhs to meet lhs dims if input is list, tuple or np.ndarray.

    Parameters
    ----------
    left : DataFrame
    right : Any
    axis: int, str, or None
    flex: bool or None, default False
        Whether this is a flex op, in which case we reindex.
        None indicates not to check for alignment.
    level : int or level name, default None

    Returns
    -------
    left : DataFrame
    right : Any
    """

    def to_series(right):
        msg = "Unable to coerce to Series, length must be {req_len}: given {given_len}"
        if axis is not None and left._get_axis_name(axis) == "index":
            if len(left.index) != len(right):
                raise ValueError(
                    msg.format(req_len=len(left.index), given_len=len(right))
                )
            right = left._constructor_sliced(right, index=left.index)
        else:
            if len(left.columns) != len(right):
                raise ValueError(
                    msg.format(req_len=len(left.columns), given_len=len(right))
                )
            right = left._constructor_sliced(right, index=left.columns)
        return right

    if isinstance(right, np.ndarray):

        if right.ndim == 1:
            right = to_series(right)

        elif right.ndim == 2:
            if right.shape == left.shape:
                right = left._constructor(right, index=left.index, columns=left.columns)

            elif right.shape[0] == left.shape[0] and right.shape[1] == 1:
                # Broadcast across columns
                right = np.broadcast_to(right, left.shape)
                right = left._constructor(right, index=left.index, columns=left.columns)

            elif right.shape[1] == left.shape[1] and right.shape[0] == 1:
                # Broadcast along rows
                right = to_series(right[0, :])

            else:
                raise ValueError(
                    "Unable to coerce to DataFrame, shape "
                    f"must be {left.shape}: given {right.shape}"
                )

        elif right.ndim > 2:
            raise ValueError(
                "Unable to coerce to Series/DataFrame, "
                f"dimension must be <= 2: {right.shape}"
            )

    elif is_list_like(right) and not isinstance(right, (ABCSeries, ABCDataFrame)):
        # GH17901
        right = to_series(right)

    if flex is not None and isinstance(right, ABCDataFrame):
        if not left._indexed_same(right):
            if flex:
                left, right = left.align(right, join="outer", level=level, copy=False)
            else:
                raise ValueError(
                    "Can only compare identically-labeled DataFrame objects"
                )
    elif isinstance(right, ABCSeries):
        # axis=1 is default for DataFrame-with-Series op
        axis = left._get_axis_number(axis) if axis is not None else 1
        left, right = left.align(
            right, join="outer", axis=axis, level=level, copy=False
        )
        right = _maybe_align_series_as_frame(left, right, axis)

    return left, right


def _should_reindex_frame_op(
    left: "DataFrame", right, op, axis, default_axis, fill_value, level
) -> bool:
    """
    Check if this is an operation between DataFrames that will need to reindex.
    """
    assert isinstance(left, ABCDataFrame)

    if op is operator.pow or op is rpow:
        # GH#32685 pow has special semantics for operating with null values
        return False

    if not isinstance(right, ABCDataFrame):
        return False

    if fill_value is None and level is None and axis is default_axis:
        # TODO: any other cases we should handle here?
        cols = left.columns.intersection(right.columns)

        # Intersection is always unique so we have to check the unique columns
        left_uniques = left.columns.unique()
        right_uniques = right.columns.unique()
        if not (cols.equals(left_uniques) and cols.equals(right_uniques)):
            return True

    return False


def _frame_arith_method_with_reindex(
    left: "DataFrame", right: "DataFrame", op
) -> "DataFrame":
    """
    For DataFrame-with-DataFrame operations that require reindexing,
    operate only on shared columns, then reindex.

    Parameters
    ----------
    left : DataFrame
    right : DataFrame
    op : binary operator

    Returns
    -------
    DataFrame
    """
    # GH#31623, only operate on shared columns
    cols, lcols, rcols = left.columns.join(
        right.columns, how="inner", level=None, return_indexers=True
    )

    new_left = left.iloc[:, lcols]
    new_right = right.iloc[:, rcols]
    result = op(new_left, new_right)

    # Do the join on the columns instead of using _align_method_FRAME
    #  to avoid constructing two potentially large/sparse DataFrames
    join_columns, _, _ = left.columns.join(
        right.columns, how="outer", level=None, return_indexers=True
    )

    if result.columns.has_duplicates:
        # Avoid reindexing with a duplicate axis.
        # https://github.com/pandas-dev/pandas/issues/35194
        indexer, _ = result.columns.get_indexer_non_unique(join_columns)
        indexer = algorithms.unique1d(indexer)
        result = result._reindex_with_indexers(
            {1: [join_columns, indexer]}, allow_dups=True
        )
    else:
        result = result.reindex(join_columns, axis=1)

    return result


def _maybe_align_series_as_frame(frame: "DataFrame", series: "Series", axis: int):
    """
    If the Series operand is not EA-dtype, we can broadcast to 2D and operate
    blockwise.
    """
    rvalues = series._values
    if not isinstance(rvalues, np.ndarray):
        # TODO(EA2D): no need to special-case with 2D EAs
        if rvalues.dtype == "datetime64[ns]" or rvalues.dtype == "timedelta64[ns]":
            # We can losslessly+cheaply cast to ndarray
            rvalues = np.asarray(rvalues)
        else:
            return series

    if axis == 0:
        rvalues = rvalues.reshape(-1, 1)
    else:
        rvalues = rvalues.reshape(1, -1)

    rvalues = np.broadcast_to(rvalues, frame.shape)
    return type(frame)(rvalues, index=frame.index, columns=frame.columns)


def _arith_method_FRAME(cls: Type["DataFrame"], op, special: bool):
    # This is the only function where `special` can be either True or False
    op_name = _get_op_name(op, special)
    default_axis = _get_frame_op_default_axis(op_name)

    na_op = get_array_op(op)

    if op_name in _op_descriptions:
        # i.e. include "add" but not "__add__"
        doc = _make_flex_doc(op_name, "dataframe")
    else:
        doc = _arith_doc_FRAME % op_name

    @Appender(doc)
    def f(self, other, axis=default_axis, level=None, fill_value=None):

        if _should_reindex_frame_op(
            self, other, op, axis, default_axis, fill_value, level
        ):
            return _frame_arith_method_with_reindex(self, other, op)

        if isinstance(other, ABCSeries) and fill_value is not None:
            # TODO: We could allow this in cases where we end up going
            #  through the DataFrame path
            raise NotImplementedError(f"fill_value {fill_value} not supported.")

        axis = self._get_axis_number(axis) if axis is not None else 1

        # TODO: why are we passing flex=True instead of flex=not special?
        #  15 tests fail if we pass flex=not special instead
        self, other = _align_method_FRAME(self, other, axis, flex=True, level=level)

        if isinstance(other, ABCDataFrame):
            # Another DataFrame
            new_data = self._combine_frame(other, na_op, fill_value)

        elif isinstance(other, ABCSeries):
            new_data = dispatch_to_series(self, other, op, axis=axis)
        else:
            # in this case we always have `np.ndim(other) == 0`
            if fill_value is not None:
                self = self.fillna(fill_value)

            new_data = dispatch_to_series(self, other, op)

        return self._construct_result(new_data)

    f.__name__ = op_name

    return f


def _flex_comp_method_FRAME(cls: Type["DataFrame"], op, special: bool):
    assert not special  # "special" also means "not flex"
    op_name = _get_op_name(op, special)
    default_axis = _get_frame_op_default_axis(op_name)
    assert default_axis == "columns", default_axis  # because we are not "special"

    doc = _flex_comp_doc_FRAME.format(
        op_name=op_name, desc=_op_descriptions[op_name]["desc"]
    )

    @Appender(doc)
    def f(self, other, axis=default_axis, level=None):
        axis = self._get_axis_number(axis) if axis is not None else 1

        self, other = _align_method_FRAME(self, other, axis, flex=True, level=level)

        new_data = dispatch_to_series(self, other, op, axis=axis)
        return self._construct_result(new_data)

    f.__name__ = op_name

    return f


def _comp_method_FRAME(cls: Type["DataFrame"], op, special: bool):
    assert special  # "special" also means "not flex"
    op_name = _get_op_name(op, special)

    @Appender(f"Wrapper for comparison method {op_name}")
    def f(self, other):
        axis = 1  # only relevant for Series other case

        self, other = _align_method_FRAME(self, other, axis, level=None, flex=False)

        # See GH#4537 for discussion of scalar op behavior
        new_data = dispatch_to_series(self, other, op, axis=axis)
        return self._construct_result(new_data)

    f.__name__ = op_name

    return f
