"""
Functions for arithmetic and comparison operations on NumPy arrays and
ExtensionArrays.
"""
import numpy as np

from pandas._libs import ops as libops

from pandas.core.dtypes.cast import (
    construct_1d_object_array_from_listlike,
    find_common_type,
    maybe_upcast_putmask,
)
from pandas.core.dtypes.common import is_object_dtype, is_scalar
from pandas.core.dtypes.generic import ABCIndex, ABCSeries
from pandas.core.dtypes.missing import notna

from pandas.core.ops import missing
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

    return na_op
