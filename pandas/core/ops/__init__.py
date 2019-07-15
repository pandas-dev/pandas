"""
Arithmetic operations for PandasObjects

This is not a public API.
"""
import datetime
import operator
import textwrap
from typing import Any, Callable
import warnings

import numpy as np

from pandas._libs import Timedelta, Timestamp, lib, ops as libops
from pandas.errors import NullFrequencyError
from pandas.util._decorators import Appender

from pandas.core.dtypes.cast import (
    construct_1d_object_array_from_listlike,
    find_common_type,
    maybe_upcast_putmask,
)
from pandas.core.dtypes.common import (
    ensure_object,
    is_bool_dtype,
    is_categorical_dtype,
    is_datetime64_dtype,
    is_datetime64tz_dtype,
    is_datetimelike_v_numeric,
    is_extension_array_dtype,
    is_integer_dtype,
    is_list_like,
    is_object_dtype,
    is_period_dtype,
    is_scalar,
    is_timedelta64_dtype,
    needs_i8_conversion,
)
from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCDatetimeArray,
    ABCIndex,
    ABCIndexClass,
    ABCSeries,
    ABCSparseArray,
    ABCSparseSeries,
    ABCTimedeltaArray,
)
from pandas.core.dtypes.missing import isna, notna

import pandas as pd
from pandas._typing import ArrayLike
import pandas.core.common as com

from . import missing
from .docstrings import (
    _arith_doc_FRAME,
    _flex_comp_doc_FRAME,
    _make_flex_doc,
    _op_descriptions,
)
from .roperator import (  # noqa:F401
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
    # `left` is always a pd.Series when called from within ops
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


def maybe_upcast_for_op(obj):
    """
    Cast non-pandas objects to pandas types to unify behavior of arithmetic
    and comparison operations.

    Parameters
    ----------
    obj: object

    Returns
    -------
    out : object

    Notes
    -----
    Be careful to call this *after* determining the `name` attribute to be
    attached to the result of the arithmetic operation.
    """
    if type(obj) is datetime.timedelta:
        # GH#22390  cast up to Timedelta to rely on Timedelta
        # implementation; otherwise operation against numeric-dtype
        # raises TypeError
        return Timedelta(obj)
    elif isinstance(obj, np.timedelta64) and not isna(obj):
        # In particular non-nanosecond timedelta64 needs to be cast to
        #  nanoseconds, or else we get undesired behavior like
        #  np.timedelta64(3, 'D') / 2 == np.timedelta64(1, 'D')
        # The isna check is to avoid casting timedelta64("NaT"), which would
        #  return NaT and incorrectly be treated as a datetime-NaT.
        return Timedelta(obj)
    elif isinstance(obj, np.ndarray) and is_timedelta64_dtype(obj):
        # GH#22390 Unfortunately we need to special-case right-hand
        # timedelta64 dtypes because numpy casts integer dtypes to
        # timedelta64 when operating with timedelta64
        return pd.TimedeltaIndex(obj)
    return obj


# -----------------------------------------------------------------------------


def make_invalid_op(name):
    """
    Return a binary method that always raises a TypeError.

    Parameters
    ----------
    name : str

    Returns
    -------
    invalid_op : function
    """

    def invalid_op(self, other=None):
        raise TypeError(
            "cannot perform {name} with this index type: "
            "{typ}".format(name=name, typ=type(self).__name__)
        )

    invalid_op.__name__ = name
    return invalid_op


def _gen_eval_kwargs(name):
    """
    Find the keyword arguments to pass to numexpr for the given operation.

    Parameters
    ----------
    name : str

    Returns
    -------
    eval_kwargs : dict

    Examples
    --------
    >>> _gen_eval_kwargs("__add__")
    {}

    >>> _gen_eval_kwargs("rtruediv")
    {'reversed': True, 'truediv': True}
    """
    kwargs = {}

    # Series appear to only pass __add__, __radd__, ...
    # but DataFrame gets both these dunder names _and_ non-dunder names
    # add, radd, ...
    name = name.replace("__", "")

    if name.startswith("r"):
        if name not in ["radd", "rand", "ror", "rxor"]:
            # Exclude commutative operations
            kwargs["reversed"] = True

    if name in ["truediv", "rtruediv"]:
        kwargs["truediv"] = True

    if name in ["ne"]:
        kwargs["masker"] = True

    return kwargs


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


def _get_opstr(op, cls):
    """
    Find the operation string, if any, to pass to numexpr for this
    operation.

    Parameters
    ----------
    op : binary operator
    cls : class

    Returns
    -------
    op_str : string or None
    """
    # numexpr is available for non-sparse classes
    subtyp = getattr(cls, "_subtyp", "")
    use_numexpr = "sparse" not in subtyp

    if not use_numexpr:
        # if we're not using numexpr, then don't pass a str_rep
        return None

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
        opname = "__{opname}__".format(opname=opname)
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
    Makes copies if fill_value is not None
    """
    # TODO: can we make a no-copy implementation?
    if fill_value is not None:
        left_mask = isna(left)
        right_mask = isna(right)
        left = left.copy()
        right = right.copy()

        # one but not both
        mask = left_mask ^ right_mask
        left[left_mask & mask] = fill_value
        right[right_mask & mask] = fill_value
    return left, right


def mask_cmp_op(x, y, op):
    """
    Apply the function `op` to only non-null points in x and y.

    Parameters
    ----------
    x : array-like
    y : array-like
    op : binary operation

    Returns
    -------
    result : ndarray[bool]
    """
    xrav = x.ravel()
    result = np.empty(x.size, dtype=bool)
    if isinstance(y, (np.ndarray, ABCSeries)):
        yrav = y.ravel()
        mask = notna(xrav) & notna(yrav)
        result[mask] = op(np.array(list(xrav[mask])), np.array(list(yrav[mask])))
    else:
        mask = notna(xrav)
        result[mask] = op(np.array(list(xrav[mask])), y)

    if op == operator.ne:  # pragma: no cover
        np.putmask(result, ~mask, True)
    else:
        np.putmask(result, ~mask, False)
    result = result.reshape(x.shape)
    return result


def masked_arith_op(x, y, op):
    """
    If the given arithmetic operation fails, attempt it again on
    only the non-null elements of the input array(s).

    Parameters
    ----------
    x : np.ndarray
    y : np.ndarray, Series, Index
    op : binary operator
    """
    # For Series `x` is 1D so ravel() is a no-op; calling it anyway makes
    # the logic valid for both Series and DataFrame ops.
    xrav = x.ravel()
    assert isinstance(x, (np.ndarray, ABCSeries)), type(x)
    if isinstance(y, (np.ndarray, ABCSeries, ABCIndexClass)):
        dtype = find_common_type([x.dtype, y.dtype])
        result = np.empty(x.size, dtype=dtype)

        # PeriodIndex.ravel() returns int64 dtype, so we have
        # to work around that case.  See GH#19956
        yrav = y if is_period_dtype(y) else y.ravel()
        mask = notna(xrav) & notna(yrav)

        if yrav.shape != mask.shape:
            # FIXME: GH#5284, GH#5035, GH#19448
            # Without specifically raising here we get mismatched
            # errors in Py3 (TypeError) vs Py2 (ValueError)
            # Note: Only = an issue in DataFrame case
            raise ValueError("Cannot broadcast operands together.")

        if mask.any():
            with np.errstate(all="ignore"):
                result[mask] = op(xrav[mask], com.values_from_object(yrav[mask]))

    else:
        assert is_scalar(y), type(y)
        assert isinstance(x, np.ndarray), type(x)
        # mask is only meaningful for x
        result = np.empty(x.size, dtype=x.dtype)
        mask = notna(xrav)

        # 1 ** np.nan is 1. So we have to unmask those.
        if op == pow:
            mask = np.where(x == 1, False, mask)
        elif op == rpow:
            mask = np.where(y == 1, False, mask)

        if mask.any():
            with np.errstate(all="ignore"):
                result[mask] = op(xrav[mask], y)

    result, changed = maybe_upcast_putmask(result, ~mask, np.nan)
    result = result.reshape(x.shape)  # 2D compat
    return result


def invalid_comparison(left, right, op):
    """
    If a comparison has mismatched types and is not necessarily meaningful,
    follow python3 conventions by:

        - returning all-False for equality
        - returning all-True for inequality
        - raising TypeError otherwise

    Parameters
    ----------
    left : array-like
    right : scalar, array-like
    op : operator.{eq, ne, lt, le, gt}

    Raises
    ------
    TypeError : on inequality comparisons
    """
    if op is operator.eq:
        res_values = np.zeros(left.shape, dtype=bool)
    elif op is operator.ne:
        res_values = np.ones(left.shape, dtype=bool)
    else:
        raise TypeError(
            "Invalid comparison between dtype={dtype} and {typ}".format(
                dtype=left.dtype, typ=type(right).__name__
            )
        )
    return res_values


# -----------------------------------------------------------------------------
# Dispatch logic


def should_series_dispatch(left, right, op):
    """
    Identify cases where a DataFrame operation should dispatch to its
    Series counterpart.

    Parameters
    ----------
    left : DataFrame
    right : DataFrame
    op : binary operator

    Returns
    -------
    override : bool
    """
    if left._is_mixed_type or right._is_mixed_type:
        return True

    if not len(left.columns) or not len(right.columns):
        # ensure obj.dtypes[0] exists for each obj
        return False

    ldtype = left.dtypes.iloc[0]
    rdtype = right.dtypes.iloc[0]

    if (is_timedelta64_dtype(ldtype) and is_integer_dtype(rdtype)) or (
        is_timedelta64_dtype(rdtype) and is_integer_dtype(ldtype)
    ):
        # numpy integer dtypes as timedelta64 dtypes in this scenario
        return True

    if is_datetime64_dtype(ldtype) and is_object_dtype(rdtype):
        # in particular case where right is an array of DateOffsets
        return True

    return False


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

        def column_op(a, b):
            return {i: func(a.iloc[:, i], b) for i in range(len(a.columns))}

    elif isinstance(right, ABCDataFrame):
        assert right._indexed_same(left)

        def column_op(a, b):
            return {i: func(a.iloc[:, i], b.iloc[:, i]) for i in range(len(a.columns))}

    elif isinstance(right, ABCSeries) and axis == "columns":
        # We only get here if called via left._combine_match_columns,
        # in which case we specifically want to operate row-by-row
        assert right.index.equals(left.columns)

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

    result = left._constructor(new_data, index=left.index, copy=False)
    # Pin columns instead of passing to constructor for compat with
    # non-unique columns case
    result.columns = left.columns
    return result


def dispatch_to_index_op(op, left, right, index_class):
    """
    Wrap Series left in the given index_class to delegate the operation op
    to the index implementation.  DatetimeIndex and TimedeltaIndex perform
    type checking, timezone handling, overflow checks, etc.

    Parameters
    ----------
    op : binary operator (operator.add, operator.sub, ...)
    left : Series
    right : object
    index_class : DatetimeIndex or TimedeltaIndex

    Returns
    -------
    result : object, usually DatetimeIndex, TimedeltaIndex, or Series
    """
    left_idx = index_class(left)

    # avoid accidentally allowing integer add/sub.  For datetime64[tz] dtypes,
    # left_idx may inherit a freq from a cached DatetimeIndex.
    # See discussion in GH#19147.
    if getattr(left_idx, "freq", None) is not None:
        left_idx = left_idx._shallow_copy(freq=None)
    try:
        result = op(left_idx, right)
    except NullFrequencyError:
        # DatetimeIndex and TimedeltaIndex with freq == None raise ValueError
        # on add/sub of integers (or int-like).  We re-raise as a TypeError.
        raise TypeError(
            "incompatible type for a datetime/timedelta "
            "operation [{name}]".format(name=op.__name__)
        )
    return result


def dispatch_to_extension_op(op, left, right):
    """
    Assume that left or right is a Series backed by an ExtensionArray,
    apply the operator defined by op.
    """

    # The op calls will raise TypeError if the op is not defined
    # on the ExtensionArray

    # unbox Series and Index to arrays
    if isinstance(left, (ABCSeries, ABCIndexClass)):
        new_left = left._values
    else:
        new_left = left

    if isinstance(right, (ABCSeries, ABCIndexClass)):
        new_right = right._values
    else:
        new_right = right

    res_values = op(new_left, new_right)
    res_name = get_op_result_name(left, right)

    if op.__name__ in ["divmod", "rdivmod"]:
        return _construct_divmod_result(left, res_values, left.index, res_name)

    return _construct_result(left, res_values, left.index, res_name)


# -----------------------------------------------------------------------------
# Functions that add arithmetic methods to objects, given arithmetic factory
# methods


def _get_method_wrappers(cls):
    """
    Find the appropriate operation-wrappers to use when defining flex/special
    arithmetic, boolean, and comparison operations with the given class.

    Parameters
    ----------
    cls : class

    Returns
    -------
    arith_flex : function or None
    comp_flex : function or None
    arith_special : function
    comp_special : function
    bool_special : function

    Notes
    -----
    None is only returned for SparseArray
    """
    if issubclass(cls, ABCSparseSeries):
        # Be sure to catch this before ABCSeries and ABCSparseArray,
        # as they will both come see SparseSeries as a subclass
        arith_flex = _flex_method_SERIES
        comp_flex = _flex_method_SERIES
        arith_special = _arith_method_SPARSE_SERIES
        comp_special = _arith_method_SPARSE_SERIES
        bool_special = _bool_method_SERIES
        # TODO: I don't think the functions defined by bool_method are tested
    elif issubclass(cls, ABCSeries):
        # Just Series; SparseSeries is caught above
        arith_flex = _flex_method_SERIES
        comp_flex = _flex_method_SERIES
        arith_special = _arith_method_SERIES
        comp_special = _comp_method_SERIES
        bool_special = _bool_method_SERIES
    elif issubclass(cls, ABCDataFrame):
        # Same for DataFrame and SparseDataFrame
        arith_flex = _arith_method_FRAME
        comp_flex = _flex_comp_method_FRAME
        arith_special = _arith_method_FRAME
        comp_special = _comp_method_FRAME
        bool_special = _arith_method_FRAME
    return arith_flex, comp_flex, arith_special, comp_special, bool_special


def _create_methods(cls, arith_method, comp_method, bool_method, special):
    # creates actual methods based upon arithmetic, comp and bool method
    # constructors.

    have_divmod = issubclass(cls, ABCSeries)
    # divmod is available for Series and SparseSeries

    # yapf: disable
    new_methods = dict(
        add=arith_method(cls, operator.add, special),
        radd=arith_method(cls, radd, special),
        sub=arith_method(cls, operator.sub, special),
        mul=arith_method(cls, operator.mul, special),
        truediv=arith_method(cls, operator.truediv, special),
        floordiv=arith_method(cls, operator.floordiv, special),
        # Causes a floating point exception in the tests when numexpr enabled,
        # so for now no speedup
        mod=arith_method(cls, operator.mod, special),
        pow=arith_method(cls, operator.pow, special),
        # not entirely sure why this is necessary, but previously was included
        # so it's here to maintain compatibility
        rmul=arith_method(cls, rmul, special),
        rsub=arith_method(cls, rsub, special),
        rtruediv=arith_method(cls, rtruediv, special),
        rfloordiv=arith_method(cls, rfloordiv, special),
        rpow=arith_method(cls, rpow, special),
        rmod=arith_method(cls, rmod, special))
    # yapf: enable
    new_methods["div"] = new_methods["truediv"]
    new_methods["rdiv"] = new_methods["rtruediv"]
    if have_divmod:
        # divmod doesn't have an op that is supported by numexpr
        new_methods["divmod"] = arith_method(cls, divmod, special)
        new_methods["rdivmod"] = arith_method(cls, rdivmod, special)

    new_methods.update(
        dict(
            eq=comp_method(cls, operator.eq, special),
            ne=comp_method(cls, operator.ne, special),
            lt=comp_method(cls, operator.lt, special),
            gt=comp_method(cls, operator.gt, special),
            le=comp_method(cls, operator.le, special),
            ge=comp_method(cls, operator.ge, special),
        )
    )

    if bool_method:
        new_methods.update(
            dict(
                and_=bool_method(cls, operator.and_, special),
                or_=bool_method(cls, operator.or_, special),
                # For some reason ``^`` wasn't used in original.
                xor=bool_method(cls, operator.xor, special),
                rand_=bool_method(cls, rand_, special),
                ror_=bool_method(cls, ror_, special),
                rxor=bool_method(cls, rxor, special),
            )
        )

    if special:
        dunderize = lambda x: "__{name}__".format(name=x.strip("_"))
    else:
        dunderize = lambda x: x
    new_methods = {dunderize(k): v for k, v in new_methods.items()}
    return new_methods


def add_methods(cls, new_methods):
    for name, method in new_methods.items():
        # For most methods, if we find that the class already has a method
        # of the same name, it is OK to over-write it.  The exception is
        # inplace methods (__iadd__, __isub__, ...) for SparseArray, which
        # retain the np.ndarray versions.
        force = not (issubclass(cls, ABCSparseArray) and name.startswith("__i"))
        if force or name not in cls.__dict__:
            setattr(cls, name, method)


# ----------------------------------------------------------------------
# Arithmetic
def add_special_arithmetic_methods(cls):
    """
    Adds the full suite of special arithmetic methods (``__add__``,
    ``__sub__``, etc.) to the class.

    Parameters
    ----------
    cls : class
        special methods will be defined and pinned to this class
    """
    _, _, arith_method, comp_method, bool_method = _get_method_wrappers(cls)
    new_methods = _create_methods(
        cls, arith_method, comp_method, bool_method, special=True
    )
    # inplace operators (I feel like these should get passed an `inplace=True`
    # or just be removed

    def _wrap_inplace_method(method):
        """
        return an inplace wrapper for this method
        """

        def f(self, other):
            result = method(self, other)

            # this makes sure that we are aligned like the input
            # we are updating inplace so we want to ignore is_copy
            self._update_inplace(
                result.reindex_like(self, copy=False)._data, verify_is_copy=False
            )

            return self

        f.__name__ = "__i{name}__".format(name=method.__name__.strip("__"))
        return f

    new_methods.update(
        dict(
            __iadd__=_wrap_inplace_method(new_methods["__add__"]),
            __isub__=_wrap_inplace_method(new_methods["__sub__"]),
            __imul__=_wrap_inplace_method(new_methods["__mul__"]),
            __itruediv__=_wrap_inplace_method(new_methods["__truediv__"]),
            __ifloordiv__=_wrap_inplace_method(new_methods["__floordiv__"]),
            __imod__=_wrap_inplace_method(new_methods["__mod__"]),
            __ipow__=_wrap_inplace_method(new_methods["__pow__"]),
        )
    )

    new_methods.update(
        dict(
            __iand__=_wrap_inplace_method(new_methods["__and__"]),
            __ior__=_wrap_inplace_method(new_methods["__or__"]),
            __ixor__=_wrap_inplace_method(new_methods["__xor__"]),
        )
    )

    add_methods(cls, new_methods=new_methods)


def add_flex_arithmetic_methods(cls):
    """
    Adds the full suite of flex arithmetic methods (``pow``, ``mul``, ``add``)
    to the class.

    Parameters
    ----------
    cls : class
        flex methods will be defined and pinned to this class
    """
    flex_arith_method, flex_comp_method, _, _, _ = _get_method_wrappers(cls)
    new_methods = _create_methods(
        cls, flex_arith_method, flex_comp_method, bool_method=None, special=False
    )
    new_methods.update(
        dict(
            multiply=new_methods["mul"],
            subtract=new_methods["sub"],
            divide=new_methods["div"],
        )
    )
    # opt out of bool flex methods for now
    assert not any(kname in new_methods for kname in ("ror_", "rxor", "rand_"))

    add_methods(cls, new_methods=new_methods)


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


def _construct_result(left, result, index, name, dtype=None):
    """
    If the raw op result has a non-None name (e.g. it is an Index object) and
    the name argument is None, then passing name to the constructor will
    not be enough; we still need to override the name attribute.
    """
    out = left._constructor(result, index=index, dtype=dtype)
    out = out.__finalize__(left)
    out.name = name
    return out


def _construct_divmod_result(left, result, index, name, dtype=None):
    """divmod returns a tuple of like indexed series instead of a single series.
    """
    return (
        _construct_result(left, result[0], index=index, name=name, dtype=dtype),
        _construct_result(left, result[1], index=index, name=name, dtype=dtype),
    )


def _arith_method_SERIES(cls, op, special):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    str_rep = _get_opstr(op, cls)
    op_name = _get_op_name(op, special)
    eval_kwargs = _gen_eval_kwargs(op_name)
    construct_result = (
        _construct_divmod_result if op in [divmod, rdivmod] else _construct_result
    )

    def na_op(x, y):
        """
        Return the result of evaluating op on the passed in values.

        If native types are not compatible, try coersion to object dtype.

        Parameters
        ----------
        x : array-like
        y : array-like or scalar

        Returns
        -------
        array-like

        Raises
        ------
        TypeError : invalid operation
        """
        import pandas.core.computation.expressions as expressions

        try:
            result = expressions.evaluate(op, str_rep, x, y, **eval_kwargs)
        except TypeError:
            result = masked_arith_op(x, y, op)

        return missing.dispatch_fill_zeros(op, x, y, result)

    def wrapper(left, right):
        if isinstance(right, ABCDataFrame):
            return NotImplemented

        left, right = _align_method_SERIES(left, right)
        res_name = get_op_result_name(left, right)
        right = maybe_upcast_for_op(right)

        if is_categorical_dtype(left):
            raise TypeError(
                "{typ} cannot perform the operation "
                "{op}".format(typ=type(left).__name__, op=str_rep)
            )

        elif is_datetime64_dtype(left) or is_datetime64tz_dtype(left):
            # Give dispatch_to_index_op a chance for tests like
            # test_dt64_series_add_intlike, which the index dispatching handles
            # specifically.
            result = dispatch_to_index_op(op, left, right, pd.DatetimeIndex)
            return construct_result(
                left, result, index=left.index, name=res_name, dtype=result.dtype
            )

        elif is_extension_array_dtype(left) or (
            is_extension_array_dtype(right) and not is_scalar(right)
        ):
            # GH#22378 disallow scalar to exclude e.g. "category", "Int64"
            return dispatch_to_extension_op(op, left, right)

        elif is_timedelta64_dtype(left):
            result = dispatch_to_index_op(op, left, right, pd.TimedeltaIndex)
            return construct_result(left, result, index=left.index, name=res_name)

        elif is_timedelta64_dtype(right):
            # We should only get here with non-scalar or timedelta64('NaT')
            #  values for right
            # Note: we cannot use dispatch_to_index_op because
            #  that may incorrectly raise TypeError when we
            #  should get NullFrequencyError
            orig_right = right
            if is_scalar(right):
                # broadcast and wrap in a TimedeltaIndex
                assert np.isnat(right)
                right = np.broadcast_to(right, left.shape)
                right = pd.TimedeltaIndex(right)

            assert isinstance(right, (pd.TimedeltaIndex, ABCTimedeltaArray, ABCSeries))
            try:
                result = op(left._values, right)
            except NullFrequencyError:
                if orig_right is not right:
                    # i.e. scalar timedelta64('NaT')
                    #  We get a NullFrequencyError because we broadcast to
                    #  TimedeltaIndex, but this should be TypeError.
                    raise TypeError(
                        "incompatible type for a datetime/timedelta "
                        "operation [{name}]".format(name=op.__name__)
                    )
                raise

            # We do not pass dtype to ensure that the Series constructor
            #  does inference in the case where `result` has object-dtype.
            return construct_result(left, result, index=left.index, name=res_name)

        elif isinstance(right, (ABCDatetimeArray, pd.DatetimeIndex)):
            result = op(left._values, right)
            return construct_result(left, result, index=left.index, name=res_name)

        lvalues = left.values
        rvalues = right
        if isinstance(rvalues, (ABCSeries, ABCIndexClass)):
            rvalues = rvalues._values

        with np.errstate(all="ignore"):
            result = na_op(lvalues, rvalues)
        return construct_result(
            left, result, index=left.index, name=res_name, dtype=None
        )

    wrapper.__name__ = op_name
    return wrapper


def _comp_method_OBJECT_ARRAY(op, x, y):
    if isinstance(y, list):
        y = construct_1d_object_array_from_listlike(y)
    if isinstance(y, (np.ndarray, ABCSeries, ABCIndex)):
        if not is_object_dtype(y.dtype):
            y = y.astype(np.object_)

        if isinstance(y, (ABCSeries, ABCIndex)):
            y = y.values

        result = libops.vec_compare(x, y, op)
    else:
        result = libops.scalar_compare(x, y, op)
    return result


def _comp_method_SERIES(cls, op, special):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    op_name = _get_op_name(op, special)
    masker = _gen_eval_kwargs(op_name).get("masker", False)

    def na_op(x, y):
        # TODO:
        # should have guarantess on what x, y can be type-wise
        # Extension Dtypes are not called here

        # Checking that cases that were once handled here are no longer
        # reachable.
        assert not (is_categorical_dtype(y) and not is_scalar(y))

        if is_object_dtype(x.dtype):
            result = _comp_method_OBJECT_ARRAY(op, x, y)

        elif is_datetimelike_v_numeric(x, y):
            return invalid_comparison(x, y, op)

        else:

            # we want to compare like types
            # we only want to convert to integer like if
            # we are not NotImplemented, otherwise
            # we would allow datetime64 (but viewed as i8) against
            # integer comparisons

            # we have a datetime/timedelta and may need to convert
            assert not needs_i8_conversion(x)
            mask = None
            if not is_scalar(y) and needs_i8_conversion(y):
                mask = isna(x) | isna(y)
                y = y.view("i8")
                x = x.view("i8")

            method = getattr(x, op_name, None)
            if method is not None:
                with np.errstate(all="ignore"):
                    result = method(y)
                if result is NotImplemented:
                    return invalid_comparison(x, y, op)
            else:
                result = op(x, y)

            if mask is not None and mask.any():
                result[mask] = masker

        return result

    def wrapper(self, other, axis=None):
        # Validate the axis parameter
        if axis is not None:
            self._get_axis_number(axis)

        res_name = get_op_result_name(self, other)

        if isinstance(other, list):
            # TODO: same for tuples?
            other = np.asarray(other)

        if isinstance(other, ABCDataFrame):  # pragma: no cover
            # Defer to DataFrame implementation; fail early
            return NotImplemented

        elif isinstance(other, ABCSeries) and not self._indexed_same(other):
            raise ValueError("Can only compare identically-labeled " "Series objects")

        elif is_categorical_dtype(self):
            # Dispatch to Categorical implementation; pd.CategoricalIndex
            # behavior is non-canonical GH#19513
            res_values = dispatch_to_index_op(op, self, other, pd.Categorical)
            return self._constructor(res_values, index=self.index, name=res_name)

        elif is_datetime64_dtype(self) or is_datetime64tz_dtype(self):
            # Dispatch to DatetimeIndex to ensure identical
            # Series/Index behavior
            if isinstance(other, datetime.date) and not isinstance(
                other, datetime.datetime
            ):
                # https://github.com/pandas-dev/pandas/issues/21152
                # Compatibility for difference between Series comparison w/
                # datetime and date
                msg = (
                    "Comparing Series of datetimes with 'datetime.date'.  "
                    "Currently, the 'datetime.date' is coerced to a "
                    "datetime. In the future pandas will not coerce, "
                    "and {future}. "
                    "To retain the current behavior, "
                    "convert the 'datetime.date' to a datetime with "
                    "'pd.Timestamp'."
                )

                if op in {operator.lt, operator.le, operator.gt, operator.ge}:
                    future = "a TypeError will be raised"
                else:
                    future = (
                        "'the values will not compare equal to the " "'datetime.date'"
                    )
                msg = "\n".join(textwrap.wrap(msg.format(future=future)))
                warnings.warn(msg, FutureWarning, stacklevel=2)
                other = Timestamp(other)

            res_values = dispatch_to_index_op(op, self, other, pd.DatetimeIndex)

            return self._constructor(res_values, index=self.index, name=res_name)

        elif is_timedelta64_dtype(self):
            res_values = dispatch_to_index_op(op, self, other, pd.TimedeltaIndex)
            return self._constructor(res_values, index=self.index, name=res_name)

        elif is_extension_array_dtype(self) or (
            is_extension_array_dtype(other) and not is_scalar(other)
        ):
            # Note: the `not is_scalar(other)` condition rules out
            # e.g. other == "category"
            return dispatch_to_extension_op(op, self, other)

        elif isinstance(other, ABCSeries):
            # By this point we have checked that self._indexed_same(other)
            res_values = na_op(self.values, other.values)
            # rename is needed in case res_name is None and res_values.name
            # is not.
            return self._constructor(
                res_values, index=self.index, name=res_name
            ).rename(res_name)

        elif isinstance(other, (np.ndarray, ABCIndexClass)):
            # do not check length of zerodim array
            # as it will broadcast
            if other.ndim != 0 and len(self) != len(other):
                raise ValueError("Lengths must match to compare")

            res_values = na_op(self.values, np.asarray(other))
            result = self._constructor(res_values, index=self.index)
            # rename is needed in case res_name is None and self.name
            # is not.
            return result.__finalize__(self).rename(res_name)

        elif is_scalar(other) and isna(other):
            # numpy does not like comparisons vs None
            if op is operator.ne:
                res_values = np.ones(len(self), dtype=bool)
            else:
                res_values = np.zeros(len(self), dtype=bool)
            return self._constructor(
                res_values, index=self.index, name=res_name, dtype="bool"
            )

        else:
            values = self.to_numpy()

            with np.errstate(all="ignore"):
                res = na_op(values, other)
            if is_scalar(res):
                raise TypeError(
                    "Could not compare {typ} type with Series".format(typ=type(other))
                )

            # always return a full value series here
            res_values = com.values_from_object(res)
            return self._constructor(
                res_values, index=self.index, name=res_name, dtype="bool"
            )

    wrapper.__name__ = op_name
    return wrapper


def _bool_method_SERIES(cls, op, special):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    op_name = _get_op_name(op, special)

    def na_op(x, y):
        try:
            result = op(x, y)
        except TypeError:
            assert not isinstance(y, (list, ABCSeries, ABCIndexClass))
            if isinstance(y, np.ndarray):
                # bool-bool dtype operations should be OK, should not get here
                assert not (is_bool_dtype(x) and is_bool_dtype(y))
                x = ensure_object(x)
                y = ensure_object(y)
                result = libops.vec_binop(x, y, op)
            else:
                # let null fall thru
                assert lib.is_scalar(y)
                if not isna(y):
                    y = bool(y)
                try:
                    result = libops.scalar_binop(x, y, op)
                except (
                    TypeError,
                    ValueError,
                    AttributeError,
                    OverflowError,
                    NotImplementedError,
                ):
                    raise TypeError(
                        "cannot compare a dtyped [{dtype}] array "
                        "with a scalar of type [{typ}]".format(
                            dtype=x.dtype, typ=type(y).__name__
                        )
                    )

        return result

    fill_int = lambda x: x.fillna(0)
    fill_bool = lambda x: x.fillna(False).astype(bool)

    def wrapper(self, other):
        is_self_int_dtype = is_integer_dtype(self.dtype)

        self, other = _align_method_SERIES(self, other, align_asobject=True)
        res_name = get_op_result_name(self, other)

        if isinstance(other, ABCDataFrame):
            # Defer to DataFrame implementation; fail early
            return NotImplemented

        elif isinstance(other, (ABCSeries, ABCIndexClass)):
            is_other_int_dtype = is_integer_dtype(other.dtype)
            other = fill_int(other) if is_other_int_dtype else fill_bool(other)

            ovalues = other.values
            finalizer = lambda x: x

        else:
            # scalars, list, tuple, np.array
            is_other_int_dtype = is_integer_dtype(np.asarray(other))
            if is_list_like(other) and not isinstance(other, np.ndarray):
                # TODO: Can we do this before the is_integer_dtype check?
                # could the is_integer_dtype check be checking the wrong
                # thing?  e.g. other = [[0, 1], [2, 3], [4, 5]]?
                other = construct_1d_object_array_from_listlike(other)

            ovalues = other
            finalizer = lambda x: x.__finalize__(self)

        # For int vs int `^`, `|`, `&` are bitwise operators and return
        #   integer dtypes.  Otherwise these are boolean ops
        filler = fill_int if is_self_int_dtype and is_other_int_dtype else fill_bool
        res_values = na_op(self.values, ovalues)
        unfilled = self._constructor(res_values, index=self.index, name=res_name)
        filled = filler(unfilled)
        return finalizer(filled)

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

            return self._constructor(op(self, other), self.index).__finalize__(self)

    flex_wrapper.__name__ = name
    return flex_wrapper


# -----------------------------------------------------------------------------
# DataFrame


def _combine_series_frame(self, other, func, fill_value=None, axis=None, level=None):
    """
    Apply binary operator `func` to self, other using alignment and fill
    conventions determined by the fill_value, axis, and level kwargs.

    Parameters
    ----------
    self : DataFrame
    other : Series
    func : binary operator
    fill_value : object, default None
    axis : {0, 1, 'columns', 'index', None}, default None
    level : int or None, default None

    Returns
    -------
    result : DataFrame
    """
    if fill_value is not None:
        raise NotImplementedError(
            "fill_value {fill} not supported.".format(fill=fill_value)
        )

    if axis is not None:
        axis = self._get_axis_number(axis)
        if axis == 0:
            return self._combine_match_index(other, func, level=level)
        else:
            return self._combine_match_columns(other, func, level=level)
    else:
        if not len(other):
            return self * np.nan

        if not len(self):
            # Ambiguous case, use _series so works with DataFrame
            return self._constructor(
                data=self._series, index=self.index, columns=self.columns
            )

        # default axis is columns
        return self._combine_match_columns(other, func, level=level)


def _align_method_FRAME(left, right, axis):
    """ convert rhs to meet lhs dims if input is list, tuple or np.ndarray """

    def to_series(right):
        msg = (
            "Unable to coerce to Series, length must be {req_len}: " "given {given_len}"
        )
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
                    "must be {req_shape}: given {given_shape}".format(
                        req_shape=left.shape, given_shape=right.shape
                    )
                )

        elif right.ndim > 2:
            raise ValueError(
                "Unable to coerce to Series/DataFrame, dim "
                "must be <= 2: {dim}".format(dim=right.shape)
            )

    elif is_list_like(right) and not isinstance(right, (ABCSeries, ABCDataFrame)):
        # GH17901
        right = to_series(right)

    return right


def _arith_method_FRAME(cls, op, special):
    str_rep = _get_opstr(op, cls)
    op_name = _get_op_name(op, special)
    eval_kwargs = _gen_eval_kwargs(op_name)
    default_axis = _get_frame_op_default_axis(op_name)

    def na_op(x, y):
        import pandas.core.computation.expressions as expressions

        try:
            result = expressions.evaluate(op, str_rep, x, y, **eval_kwargs)
        except TypeError:
            result = masked_arith_op(x, y, op)

        return missing.dispatch_fill_zeros(op, x, y, result)

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
            return self._combine_frame(other, pass_op, fill_value, level)
        elif isinstance(other, ABCSeries):
            # For these values of `axis`, we end up dispatching to Series op,
            # so do not want the masked op.
            pass_op = op if axis in [0, "columns", None] else na_op
            return _combine_series_frame(
                self, other, pass_op, fill_value=fill_value, axis=axis, level=level
            )
        else:
            if fill_value is not None:
                self = self.fillna(fill_value)

            assert np.ndim(other) == 0
            return self._combine_const(other, op)

    f.__name__ = op_name

    return f


def _flex_comp_method_FRAME(cls, op, special):
    str_rep = _get_opstr(op, cls)
    op_name = _get_op_name(op, special)
    default_axis = _get_frame_op_default_axis(op_name)

    def na_op(x, y):
        try:
            with np.errstate(invalid="ignore"):
                result = op(x, y)
        except TypeError:
            result = mask_cmp_op(x, y, op)
        return result

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
            return dispatch_to_series(self, other, na_op, str_rep)

        elif isinstance(other, ABCSeries):
            return _combine_series_frame(
                self, other, na_op, fill_value=None, axis=axis, level=level
            )
        else:
            assert np.ndim(other) == 0, other
            return self._combine_const(other, na_op)

    f.__name__ = op_name

    return f


def _comp_method_FRAME(cls, func, special):
    str_rep = _get_opstr(func, cls)
    op_name = _get_op_name(func, special)

    @Appender("Wrapper for comparison method {name}".format(name=op_name))
    def f(self, other):

        other = _align_method_FRAME(self, other, axis=None)

        if isinstance(other, ABCDataFrame):
            # Another DataFrame
            if not self._indexed_same(other):
                raise ValueError(
                    "Can only compare identically-labeled " "DataFrame objects"
                )
            return dispatch_to_series(self, other, func, str_rep)

        elif isinstance(other, ABCSeries):
            return _combine_series_frame(
                self, other, func, fill_value=None, axis=None, level=None
            )
        else:

            # straight boolean comparisons we want to allow all columns
            # (regardless of dtype to pass thru) See #4537 for discussion.
            res = self._combine_const(other, func)
            return res.fillna(True).astype(bool)

    f.__name__ = op_name

    return f


# -----------------------------------------------------------------------------
# Sparse


def _cast_sparse_series_op(left, right, opname):
    """
    For SparseSeries operation, coerce to float64 if the result is expected
    to have NaN or inf values

    Parameters
    ----------
    left : SparseArray
    right : SparseArray
    opname : str

    Returns
    -------
    left : SparseArray
    right : SparseArray
    """
    from pandas.core.sparse.api import SparseDtype

    opname = opname.strip("_")

    # TODO: This should be moved to the array?
    if is_integer_dtype(left) and is_integer_dtype(right):
        # series coerces to float64 if result should have NaN/inf
        if opname in ("floordiv", "mod") and (right.to_dense() == 0).any():
            left = left.astype(SparseDtype(np.float64, left.fill_value))
            right = right.astype(SparseDtype(np.float64, right.fill_value))
        elif opname in ("rfloordiv", "rmod") and (left.to_dense() == 0).any():
            left = left.astype(SparseDtype(np.float64, left.fill_value))
            right = right.astype(SparseDtype(np.float64, right.fill_value))

    return left, right


def _arith_method_SPARSE_SERIES(cls, op, special):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    op_name = _get_op_name(op, special)

    def wrapper(self, other):
        if isinstance(other, ABCDataFrame):
            return NotImplemented
        elif isinstance(other, ABCSeries):
            if not isinstance(other, ABCSparseSeries):
                other = other.to_sparse(fill_value=self.fill_value)
            return _sparse_series_op(self, other, op, op_name)
        elif is_scalar(other):
            with np.errstate(all="ignore"):
                new_values = op(self.values, other)
            return self._constructor(new_values, index=self.index, name=self.name)
        else:  # pragma: no cover
            raise TypeError(
                "operation with {other} not supported".format(other=type(other))
            )

    wrapper.__name__ = op_name
    return wrapper


def _sparse_series_op(left, right, op, name):
    left, right = left.align(right, join="outer", copy=False)
    new_index = left.index
    new_name = get_op_result_name(left, right)

    from pandas.core.arrays.sparse import _sparse_array_op

    lvalues, rvalues = _cast_sparse_series_op(left.values, right.values, name)
    result = _sparse_array_op(lvalues, rvalues, op, name)
    return left._constructor(result, index=new_index, name=new_name)


def maybe_dispatch_ufunc_to_dunder_op(
    self: ArrayLike, ufunc: Callable, method: str, *inputs: ArrayLike, **kwargs: Any
):
    """
    Dispatch a ufunc to the equivalent dunder method.

    Parameters
    ----------
    self : ArrayLike
        The array whose dunder method we dispatch to
    ufunc : Callable
        A NumPy ufunc
    method : {'reduce', 'accumulate', 'reduceat', 'outer', 'at', '__call__'}
    inputs : ArrayLike
        The input arrays.
    kwargs : Any
        The additional keyword arguments, e.g. ``out``.

    Returns
    -------
    result : Any
        The result of applying the ufunc
    """
    # special has the ufuncs we dispatch to the dunder op on
    special = {
        "add",
        "sub",
        "mul",
        "pow",
        "mod",
        "floordiv",
        "truediv",
        "divmod",
        "eq",
        "ne",
        "lt",
        "gt",
        "le",
        "ge",
        "remainder",
        "matmul",
    }
    aliases = {
        "subtract": "sub",
        "multiply": "mul",
        "floor_divide": "floordiv",
        "true_divide": "truediv",
        "power": "pow",
        "remainder": "mod",
        "divide": "div",
        "equal": "eq",
        "not_equal": "ne",
        "less": "lt",
        "less_equal": "le",
        "greater": "gt",
        "greater_equal": "ge",
    }

    # For op(., Array) -> Array.__r{op}__
    flipped = {
        "lt": "__gt__",
        "le": "__ge__",
        "gt": "__lt__",
        "ge": "__le__",
        "eq": "__eq__",
        "ne": "__ne__",
    }

    op_name = ufunc.__name__
    op_name = aliases.get(op_name, op_name)

    def not_implemented(*args, **kwargs):
        return NotImplemented

    if method == "__call__" and op_name in special and kwargs.get("out") is None:
        if isinstance(inputs[0], type(self)):
            name = "__{}__".format(op_name)
            return getattr(self, name, not_implemented)(inputs[1])
        else:
            name = flipped.get(op_name, "__r{}__".format(op_name))
            return getattr(self, name, not_implemented)(inputs[0])
    else:
        return NotImplemented
