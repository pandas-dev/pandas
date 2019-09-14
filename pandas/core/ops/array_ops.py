"""
Functions for arithmetic and comparison operations on NumPy arrays and
ExtensionArrays.
"""
import operator

import numpy as np

from pandas._libs import Timestamp, lib, ops as libops

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
from pandas.core.dtypes.generic import (
    ABCDatetimeArray,
    ABCDatetimeIndex,
    ABCExtensionArray,
    ABCIndex,
    ABCIndexClass,
    ABCSeries,
    ABCTimedeltaArray,
    ABCTimedeltaIndex,
)
from pandas.core.dtypes.missing import isna, notna

from pandas.core.construction import extract_array
from pandas.core.ops import missing
from pandas.core.ops.invalid import invalid_comparison
from pandas.core.ops.roperator import rpow


def comp_method_OBJECT_ARRAY(op, x, y):
    if isinstance(y, list):
        y = construct_1d_object_array_from_listlike(y)

    # TODO: Should the checks below be ABCIndexClass?
    if isinstance(y, (np.ndarray, ABCSeries, ABCIndex)):
        # TODO: should this be ABCIndexClass??
        if not is_object_dtype(y.dtype):
            y = y.astype(np.object_)

        if isinstance(y, (ABCSeries, ABCIndex)):
            y = y.values

        result = libops.vec_compare(x, y, op)
    else:
        result = libops.scalar_compare(x, y, op)
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
    assert isinstance(x, np.ndarray), type(x)
    if isinstance(y, np.ndarray):
        dtype = find_common_type([x.dtype, y.dtype])
        result = np.empty(x.size, dtype=dtype)

        # NB: ravel() is only safe since y is ndarray; for e.g. PeriodIndex
        #  we would get int64 dtype, see GH#19956
        yrav = y.ravel()
        mask = notna(xrav) & notna(yrav)

        if yrav.shape != mask.shape:
            # FIXME: GH#5284, GH#5035, GH#19448
            # Without specifically raising here we get mismatched
            # errors in Py3 (TypeError) vs Py2 (ValueError)
            # Note: Only = an issue in DataFrame case
            raise ValueError("Cannot broadcast operands together.")

        if mask.any():
            with np.errstate(all="ignore"):
                result[mask] = op(xrav[mask], yrav[mask])

    else:
        if not is_scalar(y):
            raise TypeError(type(y))

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

    result, changed = maybe_upcast_putmask(result, ~mask, np.nan)
    result = result.reshape(x.shape)  # 2D compat
    return result


def define_na_arithmetic_op(op, str_rep, eval_kwargs):
    def na_op(x, y):
        return na_arithmetic_op(x, y, op, str_rep, eval_kwargs)

    return na_op


def na_arithmetic_op(left, right, op, str_rep, eval_kwargs):
    """
    Return the result of evaluating op on the passed in values.

    If native types are not compatible, try coersion to object dtype.

    Parameters
    ----------
    left : np.ndarray
    right : np.ndarray or scalar
    str_rep : str or None
    eval_kwargs : kwargs to pass to expressions

    Returns
    -------
    array-like

    Raises
    ------
    TypeError : invalid operation
    """
    import pandas.core.computation.expressions as expressions

    try:
        result = expressions.evaluate(op, str_rep, left, right, **eval_kwargs)
    except TypeError:
        result = masked_arith_op(left, right, op)

    return missing.dispatch_fill_zeros(op, left, right, result)


def arithmetic_op(left, right, op, str_rep, eval_kwargs):

    from pandas.core.ops import (
        maybe_upcast_for_op,
        should_extension_dispatch,
        dispatch_to_extension_op,
    )

    keep_null_freq = isinstance(
        right,
        (
            ABCDatetimeIndex,
            ABCDatetimeArray,
            ABCTimedeltaIndex,
            ABCTimedeltaArray,
            Timestamp,
        ),
    )

    # NB: We assume that extract_array has already been called on `left`, but
    #  cannot make the same assumption about `right`.  This is because we need
    #  to define `keep_null_freq` before calling extract_array on it.
    lvalues = left
    rvalues = extract_array(right, extract_numpy=True)

    rvalues = maybe_upcast_for_op(rvalues, lvalues.shape)

    if should_extension_dispatch(left, rvalues) or isinstance(
        rvalues, (ABCTimedeltaArray, ABCDatetimeArray, Timestamp)
    ):
        # TimedeltaArray, DatetimeArray, and Timestamp are included here
        #  because they have `freq` attribute which is handled correctly
        #  by dispatch_to_extension_op.
        res_values = dispatch_to_extension_op(op, lvalues, rvalues, keep_null_freq)

    else:
        with np.errstate(all="ignore"):
            res_values = na_arithmetic_op(lvalues, rvalues, op, str_rep, eval_kwargs)

    return res_values


def comparison_op(left, right, op):
    from pandas.core.ops import should_extension_dispatch, dispatch_to_extension_op

    # NB: We assume extract_array has already been called on left and right
    lvalues = left
    rvalues = right

    rvalues = lib.item_from_zerodim(rvalues)
    if isinstance(rvalues, list):
        # TODO: same for tuples?
        rvalues = np.asarray(rvalues)

    if isinstance(rvalues, (np.ndarray, ABCExtensionArray, ABCIndexClass)):
        # TODO: make this treatment consistent across ops and classes.
        #  We are not catching all listlikes here (e.g. frozenset, tuple)
        #  The ambiguous case is object-dtype.  See GH#27803
        if len(lvalues) != len(rvalues):
            raise ValueError("Lengths must match to compare")

    if should_extension_dispatch(lvalues, rvalues):
        res_values = dispatch_to_extension_op(op, lvalues, rvalues)

    elif is_scalar(rvalues) and isna(rvalues):
        # numpy does not like comparisons vs None
        if op is operator.ne:
            res_values = np.ones(len(lvalues), dtype=bool)
        else:
            res_values = np.zeros(len(lvalues), dtype=bool)

    elif is_object_dtype(lvalues.dtype):
        res_values = comp_method_OBJECT_ARRAY(op, lvalues, rvalues)

    else:
        op_name = "__{op}__".format(op=op.__name__)
        method = getattr(lvalues, op_name)
        with np.errstate(all="ignore"):
            res_values = method(rvalues)

        if res_values is NotImplemented:
            res_values = invalid_comparison(lvalues, rvalues, op)
        if is_scalar(res_values):
            raise TypeError(
                "Could not compare {typ} type with Series".format(typ=type(rvalues))
            )

    return res_values


def logical_op(left, right, op):
    from pandas.core.ops import should_extension_dispatch, dispatch_to_extension_op

    def na_op(x, y):
        try:
            result = op(x, y)
        except TypeError:
            if isinstance(y, np.ndarray):
                # bool-bool dtype operations should be OK, should not get here
                assert not (is_bool_dtype(x.dtype) and is_bool_dtype(y.dtype))
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
    lvalues = left
    rvalues = right

    if should_extension_dispatch(lvalues, rvalues):
        res_values = dispatch_to_extension_op(op, lvalues, rvalues)

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

        res_values = na_op(lvalues, rvalues)
        res_values = filler(res_values)

    return res_values
