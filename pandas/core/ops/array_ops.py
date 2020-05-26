"""
Functions for arithmetic and comparison operations on NumPy arrays and
ExtensionArrays.
"""
from datetime import timedelta
from functools import partial
import operator
from typing import Any, Tuple
import warnings

import numpy as np

from pandas._libs import Timedelta, Timestamp, lib, ops as libops
from pandas._typing import ArrayLike

from pandas.core.dtypes.cast import (
    construct_1d_object_array_from_listlike,
    find_common_type,
    maybe_upcast_putmask,
)
from pandas.core.dtypes.common import (
    ensure_object,
    is_bool_dtype,
    is_integer_dtype,
    is_list_like,
    is_object_dtype,
    is_scalar,
)
from pandas.core.dtypes.generic import ABCExtensionArray, ABCIndex, ABCSeries
from pandas.core.dtypes.missing import isna, notna

from pandas.core.ops import missing
from pandas.core.ops.dispatch import should_extension_dispatch
from pandas.core.ops.invalid import invalid_comparison
from pandas.core.ops.roperator import rpow


def comp_method_OBJECT_ARRAY(op, x, y):
    if isinstance(y, list):
        y = construct_1d_object_array_from_listlike(y)

    if isinstance(y, (np.ndarray, ABCSeries, ABCIndex)):
        # Note: these checks can be for ABCIndex and not ABCIndexClass
        #  because that is the only object-dtype class.
        if not is_object_dtype(y.dtype):
            y = y.astype(np.object_)

        if isinstance(y, (ABCSeries, ABCIndex)):
            y = y._values

        if x.shape != y.shape:
            raise ValueError("Shapes must match", x.shape, y.shape)
        result = libops.vec_compare(x.ravel(), y.ravel(), op)
    else:
        result = libops.scalar_compare(x.ravel(), y, op)
    return result.reshape(x.shape)


def masked_arith_op(x: np.ndarray, y, op):
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
    assert isinstance(x, np.ndarray), type(x)
    if isinstance(y, np.ndarray):
        dtype = find_common_type([x.dtype, y.dtype])
        result = np.empty(x.size, dtype=dtype)

        if len(x) != len(y):
            raise ValueError(x.shape, y.shape)
        else:
            ymask = notna(y)

        # NB: ravel() is only safe since y is ndarray; for e.g. PeriodIndex
        #  we would get int64 dtype, see GH#19956
        yrav = y.ravel()
        mask = notna(xrav) & ymask.ravel()

        # See GH#5284, GH#5035, GH#19448 for historical reference
        if mask.any():
            with np.errstate(all="ignore"):
                result[mask] = op(xrav[mask], yrav[mask])

    else:
        if not is_scalar(y):
            raise TypeError(
                f"Cannot broadcast np.ndarray with operand of type { type(y) }"
            )

        # mask is only meaningful for x
        result = np.empty(x.size, dtype=x.dtype)
        mask = notna(xrav)

        # 1 ** np.nan is 1. So we have to unmask those.
        if op is pow:
            mask = np.where(x == 1, False, mask)
        elif op is rpow:
            mask = np.where(y == 1, False, mask)

        if mask.any():
            with np.errstate(all="ignore"):
                result[mask] = op(xrav[mask], y)

    result, _ = maybe_upcast_putmask(result, ~mask, np.nan)
    result = result.reshape(x.shape)  # 2D compat
    return result


def na_arithmetic_op(left, right, op, is_cmp: bool = False):
    """
    Return the result of evaluating op on the passed in values.

    If native types are not compatible, try coercion to object dtype.

    Parameters
    ----------
    left : np.ndarray
    right : np.ndarray or scalar
    is_cmp : bool, default False
        If this a comparison operation.

    Returns
    -------
    array-like

    Raises
    ------
    TypeError : invalid operation
    """
    import pandas.core.computation.expressions as expressions

    try:
        result = expressions.evaluate(op, left, right)
    except TypeError:
        if is_cmp:
            # numexpr failed on comparison op, e.g. ndarray[float] > datetime
            #  In this case we do not fall back to the masked op, as that
            #  will handle complex numbers incorrectly, see GH#32047
            raise
        result = masked_arith_op(left, right, op)

    if is_cmp and (is_scalar(result) or result is NotImplemented):
        # numpy returned a scalar instead of operating element-wise
        # e.g. numeric array vs str
        return invalid_comparison(left, right, op)

    return missing.dispatch_fill_zeros(op, left, right, result)


def arithmetic_op(left: ArrayLike, right: Any, op):
    """
    Evaluate an arithmetic operation `+`, `-`, `*`, `/`, `//`, `%`, `**`, ...

    Parameters
    ----------
    left : np.ndarray or ExtensionArray
    right : object
        Cannot be a DataFrame or Index.  Series is *not* excluded.
    op : {operator.add, operator.sub, ...}
        Or one of the reversed variants from roperator.

    Returns
    -------
    ndarray or ExtensionArray
        Or a 2-tuple of these in the case of divmod or rdivmod.
    """

    # NB: We assume that extract_array has already been called
    #  on `left` and `right`.
    lvalues = maybe_upcast_datetimelike_array(left)
    rvalues = maybe_upcast_datetimelike_array(right)
    rvalues = maybe_upcast_for_op(rvalues, lvalues.shape)

    if should_extension_dispatch(lvalues, rvalues) or isinstance(rvalues, Timedelta):
        # Timedelta is included because numexpr will fail on it, see GH#31457
        res_values = op(lvalues, rvalues)

    else:
        with np.errstate(all="ignore"):
            res_values = na_arithmetic_op(lvalues, rvalues, op)

    return res_values


def comparison_op(left: ArrayLike, right: Any, op) -> ArrayLike:
    """
    Evaluate a comparison operation `=`, `!=`, `>=`, `>`, `<=`, or `<`.

    Parameters
    ----------
    left : np.ndarray or ExtensionArray
    right : object
        Cannot be a DataFrame, Series, or Index.
    op : {operator.eq, operator.ne, operator.gt, operator.ge, operator.lt, operator.le}

    Returns
    -------
    ndarray or ExtensionArray
    """
    # NB: We assume extract_array has already been called on left and right
    lvalues = maybe_upcast_datetimelike_array(left)
    rvalues = right

    rvalues = lib.item_from_zerodim(rvalues)
    if isinstance(rvalues, list):
        # TODO: same for tuples?
        rvalues = np.asarray(rvalues)

    if isinstance(rvalues, (np.ndarray, ABCExtensionArray)):
        # TODO: make this treatment consistent across ops and classes.
        #  We are not catching all listlikes here (e.g. frozenset, tuple)
        #  The ambiguous case is object-dtype.  See GH#27803
        if len(lvalues) != len(rvalues):
            raise ValueError(
                "Lengths must match to compare", lvalues.shape, rvalues.shape
            )

    if should_extension_dispatch(lvalues, rvalues):
        # Call the method on lvalues
        res_values = op(lvalues, rvalues)

    elif is_scalar(rvalues) and isna(rvalues):
        # numpy does not like comparisons vs None
        if op is operator.ne:
            res_values = np.ones(lvalues.shape, dtype=bool)
        else:
            res_values = np.zeros(lvalues.shape, dtype=bool)

    elif is_object_dtype(lvalues.dtype):
        res_values = comp_method_OBJECT_ARRAY(op, lvalues, rvalues)

    else:
        with warnings.catch_warnings():
            # suppress warnings from numpy about element-wise comparison
            warnings.simplefilter("ignore", DeprecationWarning)
            with np.errstate(all="ignore"):
                res_values = na_arithmetic_op(lvalues, rvalues, op, is_cmp=True)

    return res_values


def na_logical_op(x: np.ndarray, y, op):
    try:
        # For exposition, write:
        #  yarr = isinstance(y, np.ndarray)
        #  yint = is_integer(y) or (yarr and y.dtype.kind == "i")
        #  ybool = is_bool(y) or (yarr and y.dtype.kind == "b")
        #  xint = x.dtype.kind == "i"
        #  xbool = x.dtype.kind == "b"
        # Then Cases where this goes through without raising include:
        #  (xint or xbool) and (yint or bool)
        result = op(x, y)
    except TypeError:
        if isinstance(y, np.ndarray):
            # bool-bool dtype operations should be OK, should not get here
            assert not (is_bool_dtype(x.dtype) and is_bool_dtype(y.dtype))
            x = ensure_object(x)
            y = ensure_object(y)
            result = libops.vec_binop(x.ravel(), y.ravel(), op)
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
            ) as err:
                typ = type(y).__name__
                raise TypeError(
                    f"Cannot perform '{op.__name__}' with a dtyped [{x.dtype}] array "
                    f"and scalar of type [{typ}]"
                ) from err

    return result.reshape(x.shape)


def logical_op(left: ArrayLike, right: Any, op) -> ArrayLike:
    """
    Evaluate a logical operation `|`, `&`, or `^`.

    Parameters
    ----------
    left : np.ndarray or ExtensionArray
    right : object
        Cannot be a DataFrame, Series, or Index.
    op : {operator.and_, operator.or_, operator.xor}
        Or one of the reversed variants from roperator.

    Returns
    -------
    ndarray or ExtensionArray
    """
    fill_int = lambda x: x

    def fill_bool(x, left=None):
        # if `left` is specifically not-boolean, we do not cast to bool
        if x.dtype.kind in ["c", "f", "O"]:
            # dtypes that can hold NA
            mask = isna(x)
            if mask.any():
                x = x.astype(object)
                x[mask] = False

        if left is None or is_bool_dtype(left.dtype):
            x = x.astype(bool)
        return x

    is_self_int_dtype = is_integer_dtype(left.dtype)

    right = lib.item_from_zerodim(right)
    if is_list_like(right) and not hasattr(right, "dtype"):
        # e.g. list, tuple
        right = construct_1d_object_array_from_listlike(right)

    # NB: We assume extract_array has already been called on left and right
    lvalues = maybe_upcast_datetimelike_array(left)
    rvalues = right

    if should_extension_dispatch(lvalues, rvalues):
        # Call the method on lvalues
        res_values = op(lvalues, rvalues)

    else:
        if isinstance(rvalues, np.ndarray):
            is_other_int_dtype = is_integer_dtype(rvalues.dtype)
            rvalues = rvalues if is_other_int_dtype else fill_bool(rvalues, lvalues)

        else:
            # i.e. scalar
            is_other_int_dtype = lib.is_integer(rvalues)

        # For int vs int `^`, `|`, `&` are bitwise operators and return
        #   integer dtypes.  Otherwise these are boolean ops
        filler = fill_int if is_self_int_dtype and is_other_int_dtype else fill_bool

        res_values = na_logical_op(lvalues, rvalues, op)
        res_values = filler(res_values)  # type: ignore

    return res_values


def get_array_op(op):
    """
    Return a binary array operation corresponding to the given operator op.

    Parameters
    ----------
    op : function
        Binary operator from operator or roperator module.

    Returns
    -------
    functools.partial
    """
    if isinstance(op, partial):
        # We get here via dispatch_to_series in DataFrame case
        # TODO: avoid getting here
        return op

    op_name = op.__name__.strip("_").lstrip("r")
    if op_name == "arith_op":
        # Reached via DataFrame._combine_frame
        return op

    if op_name in {"eq", "ne", "lt", "le", "gt", "ge"}:
        return partial(comparison_op, op=op)
    elif op_name in {"and", "or", "xor", "rand", "ror", "rxor"}:
        return partial(logical_op, op=op)
    elif op_name in {
        "add",
        "sub",
        "mul",
        "truediv",
        "floordiv",
        "mod",
        "divmod",
        "pow",
    }:
        return partial(arithmetic_op, op=op)
    else:
        raise NotImplementedError(op_name)


def maybe_upcast_datetimelike_array(obj: ArrayLike) -> ArrayLike:
    """
    If we have an ndarray that is either datetime64 or timedelta64, wrap in EA.

    Parameters
    ----------
    obj : ndarray or ExtensionArray

    Returns
    -------
    ndarray or ExtensionArray
    """
    if isinstance(obj, np.ndarray):
        if obj.dtype.kind == "m":
            from pandas.core.arrays import TimedeltaArray

            return TimedeltaArray._from_sequence(obj)
        if obj.dtype.kind == "M":
            from pandas.core.arrays import DatetimeArray

            return DatetimeArray._from_sequence(obj)

    return obj


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

    if type(obj) is timedelta:
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

    elif isinstance(obj, np.ndarray) and obj.dtype.kind == "m":
        # GH#22390 Unfortunately we need to special-case right-hand
        # timedelta64 dtypes because numpy casts integer dtypes to
        # timedelta64 when operating with timedelta64
        return TimedeltaArray._from_sequence(obj)
    return obj
