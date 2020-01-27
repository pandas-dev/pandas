"""
Arithmetic operations for PandasObjects

This is not a public API.
"""
import datetime
import operator
from typing import Set, Tuple, Union

import numpy as np

from pandas._libs import Timedelta, Timestamp, lib
from pandas._libs.ops_dispatch import maybe_dispatch_ufunc_to_dunder_op  # noqa:F401
from pandas.util._decorators import Appender

from pandas.core.dtypes.common import is_list_like, is_timedelta64_dtype
from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCExtensionArray,
    ABCIndexClass,
    ABCSeries,
)
from pandas.core.dtypes.missing import isna

from pandas.core.construction import extract_array
from pandas.core.ops.array_ops import (
    arithmetic_op,
    comparison_op,
    define_na_arithmetic_op,
    get_array_op,
    logical_op,
)
from pandas.core.ops.array_ops import comp_method_OBJECT_ARRAY  # noqa:F401
from pandas.core.ops.common import unpack_zerodim_and_defer
from pandas.core.ops.dispatch import should_series_dispatch
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


COMPARISON_BINOPS: Set[str] = {
    "eq",
    "ne",
    "lt",
    "gt",
    "le",
    "ge",
}

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


def maybe_upcast_for_op(obj, shape: Tuple[int, ...]):
    """
    Cast non-pandas objects to pandas types to unify behavior of arithmetic
    and comparison operations.

    Parameters
    ----------
    obj: object
    shape : tuple[int]

    Returns
    -------
    out : object

    Notes
    -----
    Be careful to call this *after* determining the `name` attribute to be
    attached to the result of the arithmetic operation.
    """
    from pandas.core.arrays import DatetimeArray, TimedeltaArray

    if type(obj) is datetime.timedelta:
        # GH#22390  cast up to Timedelta to rely on Timedelta
        # implementation; otherwise operation against numeric-dtype
        # raises TypeError
        return Timedelta(obj)
    elif isinstance(obj, np.datetime64):
        # GH#28080 numpy casts integer-dtype to datetime64 when doing
        #  array[int] + datetime64, which we do not allow
        if isna(obj):
            # Avoid possible ambiguities with pd.NaT
            obj = obj.astype("datetime64[ns]")
            right = np.broadcast_to(obj, shape)
            return DatetimeArray(right)

        return Timestamp(obj)

    elif isinstance(obj, np.timedelta64):
        if isna(obj):
            # wrapping timedelta64("NaT") in Timedelta returns NaT,
            #  which would incorrectly be treated as a datetime-NaT, so
            #  we broadcast and wrap in a TimedeltaArray
            obj = obj.astype("timedelta64[ns]")
            right = np.broadcast_to(obj, shape)
            return TimedeltaArray(right)

        # In particular non-nanosecond timedelta64 needs to be cast to
        #  nanoseconds, or else we get undesired behavior like
        #  np.timedelta64(3, 'D') / 2 == np.timedelta64(1, 'D')
        return Timedelta(obj)

    elif isinstance(obj, np.ndarray) and is_timedelta64_dtype(obj.dtype):
        # GH#22390 Unfortunately we need to special-case right-hand
        # timedelta64 dtypes because numpy casts integer dtypes to
        # timedelta64 when operating with timedelta64
        return TimedeltaArray._from_sequence(obj)
    return obj


# -----------------------------------------------------------------------------


def _get_frame_op_default_axis(name):
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


def _get_opstr(op):
    """
    Find the operation string, if any, to pass to numexpr for this
    operation.

    Parameters
    ----------
    op : binary operator

    Returns
    -------
    op_str : string or None
    """

    return {
        operator.add: "+",
        radd: "+",
        operator.mul: "*",
        rmul: "*",
        operator.sub: "-",
        rsub: "-",
        operator.truediv: "/",
        rtruediv: "/",
        operator.floordiv: "//",
        rfloordiv: "//",
        operator.mod: None,  # TODO: Why None for mod but '%' for rmod?
        rmod: "%",
        operator.pow: "**",
        rpow: "**",
        operator.eq: "==",
        operator.ne: "!=",
        operator.le: "<=",
        operator.lt: "<",
        operator.ge: ">=",
        operator.gt: ">",
        operator.and_: "&",
        rand_: "&",
        operator.or_: "|",
        ror_: "|",
        operator.xor: "^",
        rxor: "^",
        divmod: None,
        rdivmod: None,
    }[op]


def _get_op_name(op, special):
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


def dispatch_to_series(left, right, func, str_rep=None, axis=None):
    """
    Evaluate the frame operation func(left, right) by evaluating
    column-by-column, dispatching to the Series implementation.

    Parameters
    ----------
    left : DataFrame
    right : scalar or DataFrame
    func : arithmetic or comparison operator
    str_rep : str or None, default None
    axis : {None, 0, 1, "index", "columns"}

    Returns
    -------
    DataFrame
    """
    # Note: we use iloc to access columns for compat with cases
    #       with non-unique columns.
    import pandas.core.computation.expressions as expressions

    right = lib.item_from_zerodim(right)
    if lib.is_scalar(right) or np.ndim(right) == 0:

        # Get the appropriate array-op to apply to each block's values.
        array_op = get_array_op(func, str_rep=str_rep)
        bm = left._data.apply(array_op, right=right)
        return type(left)(bm)

    elif isinstance(right, ABCDataFrame):
        assert right._indexed_same(left)

        def column_op(a, b):
            return {i: func(a.iloc[:, i], b.iloc[:, i]) for i in range(len(a.columns))}

    elif isinstance(right, ABCSeries) and axis == "columns":
        # We only get here if called via _combine_series_frame,
        # in which case we specifically want to operate row-by-row
        assert right.index.equals(left.columns)

        if right.dtype == "timedelta64[ns]":
            # ensure we treat NaT values as the correct dtype
            # Note: we do not do this unconditionally as it may be lossy or
            #  expensive for EA dtypes.
            right = np.asarray(right)

            def column_op(a, b):
                return {i: func(a.iloc[:, i], b[i]) for i in range(len(a.columns))}

        else:

            def column_op(a, b):
                return {i: func(a.iloc[:, i], b.iloc[i]) for i in range(len(a.columns))}

    elif isinstance(right, ABCSeries):
        assert right.index.equals(left.index)  # Handle other cases later

        def column_op(a, b):
            return {i: func(a.iloc[:, i], b) for i in range(len(a.columns))}

    else:
        # Remaining cases have less-obvious dispatch rules
        raise NotImplementedError(right)

    new_data = expressions.evaluate(column_op, str_rep, left, right)
    return new_data


# -----------------------------------------------------------------------------
# Series


def _align_method_SERIES(left, right, align_asobject=False):
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


def _construct_result(
    left: ABCSeries,
    result: Union[np.ndarray, ABCExtensionArray],
    index: ABCIndexClass,
    name,
):
    """
    Construct an appropriately-labelled Series from the result of an op.

    Parameters
    ----------
    left : Series
    result : ndarray or ExtensionArray
    index : Index
    name : object

    Returns
    -------
    Series
        In the case of __divmod__ or __rdivmod__, a 2-tuple of Series.
    """
    if isinstance(result, tuple):
        # produced by divmod or rdivmod
        return (
            _construct_result(left, result[0], index=index, name=name),
            _construct_result(left, result[1], index=index, name=name),
        )

    # We do not pass dtype to ensure that the Series constructor
    #  does inference in the case where `result` has object-dtype.
    out = left._constructor(result, index=index)
    out = out.__finalize__(left)

    # Set the result's name after __finalize__ is called because __finalize__
    #  would set it back to self.name
    out.name = name
    return out


def _arith_method_SERIES(cls, op, special):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    str_rep = _get_opstr(op)
    op_name = _get_op_name(op, special)

    @unpack_zerodim_and_defer(op_name)
    def wrapper(left, right):

        left, right = _align_method_SERIES(left, right)
        res_name = get_op_result_name(left, right)

        lvalues = extract_array(left, extract_numpy=True)
        rvalues = extract_array(right, extract_numpy=True)
        result = arithmetic_op(lvalues, rvalues, op, str_rep)

        return _construct_result(left, result, index=left.index, name=res_name)

    wrapper.__name__ = op_name
    return wrapper


def _comp_method_SERIES(cls, op, special):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    op_name = _get_op_name(op, special)

    @unpack_zerodim_and_defer(op_name)
    def wrapper(self, other):

        res_name = get_op_result_name(self, other)

        if isinstance(other, ABCSeries) and not self._indexed_same(other):
            raise ValueError("Can only compare identically-labeled Series objects")

        lvalues = extract_array(self, extract_numpy=True)
        rvalues = extract_array(other, extract_numpy=True)

        res_values = comparison_op(lvalues, rvalues, op)

        return _construct_result(self, res_values, index=self.index, name=res_name)

    wrapper.__name__ = op_name
    return wrapper


def _bool_method_SERIES(cls, op, special):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    op_name = _get_op_name(op, special)

    @unpack_zerodim_and_defer(op_name)
    def wrapper(self, other):
        self, other = _align_method_SERIES(self, other, align_asobject=True)
        res_name = get_op_result_name(self, other)

        lvalues = extract_array(self, extract_numpy=True)
        rvalues = extract_array(other, extract_numpy=True)

        res_values = logical_op(lvalues, rvalues, op)
        return _construct_result(self, res_values, index=self.index, name=res_name)

    wrapper.__name__ = op_name
    return wrapper


def _flex_method_SERIES(cls, op, special):
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


def _combine_series_frame(left, right, func, axis: int):
    """
    Apply binary operator `func` to self, other using alignment and fill
    conventions determined by the axis argument.

    Parameters
    ----------
    left : DataFrame
    right : Series
    func : binary operator
    axis : {0, 1}

    Returns
    -------
    result : DataFrame
    """
    # We assume that self.align(other, ...) has already been called
    if axis == 0:
        new_data = left._combine_match_index(right, func)
    else:
        new_data = dispatch_to_series(left, right, func, axis="columns")

    return left._construct_result(new_data)


def _align_method_FRAME(left, right, axis):
    """ convert rhs to meet lhs dims if input is list, tuple or np.ndarray """

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
                f"Unable to coerce to Series/DataFrame, dim must be <= 2: {right.shape}"
            )

    elif is_list_like(right) and not isinstance(right, (ABCSeries, ABCDataFrame)):
        # GH17901
        right = to_series(right)

    return right


def _arith_method_FRAME(cls, op, special):
    str_rep = _get_opstr(op)
    op_name = _get_op_name(op, special)
    default_axis = _get_frame_op_default_axis(op_name)

    na_op = define_na_arithmetic_op(op, str_rep)
    is_logical = str_rep in ["&", "|", "^"]

    if op_name in _op_descriptions:
        # i.e. include "add" but not "__add__"
        doc = _make_flex_doc(op_name, "dataframe")
    else:
        doc = _arith_doc_FRAME % op_name

    @Appender(doc)
    def f(self, other, axis=default_axis, level=None, fill_value=None):

        other = _align_method_FRAME(self, other, axis)

        if isinstance(other, ABCDataFrame):
            # Another DataFrame
            pass_op = op if should_series_dispatch(self, other, op) else na_op
            pass_op = pass_op if not is_logical else op

            left, right = self.align(other, join="outer", level=level, copy=False)
            new_data = left._combine_frame(right, pass_op, fill_value)
            return left._construct_result(new_data)

        elif isinstance(other, ABCSeries):
            # For these values of `axis`, we end up dispatching to Series op,
            # so do not want the masked op.
            pass_op = op if axis in [0, "columns", None] else na_op
            pass_op = pass_op if not is_logical else op

            if fill_value is not None:
                raise NotImplementedError(f"fill_value {fill_value} not supported.")

            axis = self._get_axis_number(axis) if axis is not None else 1
            self, other = self.align(
                other, join="outer", axis=axis, level=level, copy=False
            )
            return _combine_series_frame(self, other, pass_op, axis=axis)
        else:
            # in this case we always have `np.ndim(other) == 0`
            if fill_value is not None:
                self = self.fillna(fill_value)

            new_data = dispatch_to_series(self, other, op, str_rep)
            return self._construct_result(new_data)

    f.__name__ = op_name

    return f


def _flex_comp_method_FRAME(cls, op, special):
    str_rep = _get_opstr(op)
    op_name = _get_op_name(op, special)
    default_axis = _get_frame_op_default_axis(op_name)

    doc = _flex_comp_doc_FRAME.format(
        op_name=op_name, desc=_op_descriptions[op_name]["desc"]
    )

    @Appender(doc)
    def f(self, other, axis=default_axis, level=None):

        other = _align_method_FRAME(self, other, axis)

        if isinstance(other, ABCDataFrame):
            # Another DataFrame
            if not self._indexed_same(other):
                self, other = self.align(other, "outer", level=level, copy=False)
            new_data = dispatch_to_series(self, other, op, str_rep)
            return self._construct_result(new_data)

        elif isinstance(other, ABCSeries):
            axis = self._get_axis_number(axis) if axis is not None else 1
            self, other = self.align(
                other, join="outer", axis=axis, level=level, copy=False
            )
            return _combine_series_frame(self, other, op, axis=axis)
        else:
            # in this case we always have `np.ndim(other) == 0`
            new_data = dispatch_to_series(self, other, op)
            return self._construct_result(new_data)

    f.__name__ = op_name

    return f


def _comp_method_FRAME(cls, op, special):
    str_rep = _get_opstr(op)
    op_name = _get_op_name(op, special)

    @Appender(f"Wrapper for comparison method {op_name}")
    def f(self, other):

        other = _align_method_FRAME(self, other, axis=None)

        if isinstance(other, ABCDataFrame):
            # Another DataFrame
            if not self._indexed_same(other):
                raise ValueError(
                    "Can only compare identically-labeled DataFrame objects"
                )
            new_data = dispatch_to_series(self, other, op, str_rep)

        elif isinstance(other, ABCSeries):
            # axis=1 is default for DataFrame-with-Series op
            self, other = self.align(
                other, join="outer", axis=1, level=None, copy=False
            )
            new_data = dispatch_to_series(self, other, op, axis="columns")

        else:

            # straight boolean comparisons we want to allow all columns
            # (regardless of dtype to pass thru) See #4537 for discussion.
            new_data = dispatch_to_series(self, other, op)

        return self._construct_result(new_data)

    f.__name__ = op_name

    return f
