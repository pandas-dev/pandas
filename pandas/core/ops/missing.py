"""
Missing data handling for arithmetic operations.

In particular, pandas conventions regarding divison by zero differ
from numpy in the following ways:
    - FIXME: fill this in
"""
import operator

import numpy as np

from pandas.core.dtypes.common import is_float_dtype, is_integer_dtype, is_scalar


def fill_zeros(result, x, y, name, fill):
    """
    If this is a reversed op, then flip x,y

    If we have an integer value (or array in y)
    and we have 0's, fill them with the fill,
    return the result.

    Mask the nan's from x.
    """
    if fill is None or is_float_dtype(result):
        return result

    if name.startswith(("r", "__r")):
        x, y = y, x

    is_variable_type = hasattr(y, "dtype") or hasattr(y, "type")
    is_scalar_type = is_scalar(y)

    if not is_variable_type and not is_scalar_type:
        return result

    if is_scalar_type:
        y = np.array(y)

    if is_integer_dtype(y):

        if (y == 0).any():

            # GH 7325, mask and nans must be broadcastable (also: PR 9308)
            # Raveling and then reshaping makes np.putmask faster
            mask = ((y == 0) & ~np.isnan(result)).ravel()

            shape = result.shape
            result = result.astype("float64", copy=False).ravel()

            np.putmask(result, mask, fill)

            # if we have a fill of inf, then sign it correctly
            # (GH 6178 and PR 9308)
            if np.isinf(fill):
                signs = y if name.startswith(("r", "__r")) else x
                signs = np.sign(signs.astype("float", copy=False))
                negative_inf_mask = (signs.ravel() < 0) & mask
                np.putmask(result, negative_inf_mask, -fill)

            if "floordiv" in name:  # (PR 9308)
                nan_mask = ((y == 0) & (x == 0)).ravel()
                np.putmask(result, nan_mask, np.nan)

            result = result.reshape(shape)

    return result


def mask_zero_div_zero(x, y, result, copy=False):
    """
    Set results of 0 / 0 or 0 // 0 to np.nan, regardless of the dtypes
    of the numerator or the denominator.

    Parameters
    ----------
    x : ndarray
    y : ndarray
    result : ndarray
    copy : bool (default False)
        Whether to always create a new array or try to fill in the existing
        array if possible.

    Returns
    -------
    filled_result : ndarray

    Examples
    --------
    >>> x = np.array([1, 0, -1], dtype=np.int64)
    >>> y = 0       # int 0; numpy behavior is different with float
    >>> result = x / y
    >>> result      # raw numpy result does not fill division by zero
    array([0, 0, 0])
    >>> mask_zero_div_zero(x, y, result)
    array([ inf,  nan, -inf])
    """
    if is_scalar(y):
        y = np.array(y)

    zmask = y == 0
    if zmask.any():
        shape = result.shape

        nan_mask = (zmask & (x == 0)).ravel()
        neginf_mask = (zmask & (x < 0)).ravel()
        posinf_mask = (zmask & (x > 0)).ravel()

        if nan_mask.any() or neginf_mask.any() or posinf_mask.any():
            # Fill negative/0 with -inf, positive/0 with +inf, 0/0 with NaN
            result = result.astype("float64", copy=copy).ravel()

            np.putmask(result, nan_mask, np.nan)
            np.putmask(result, posinf_mask, np.inf)
            np.putmask(result, neginf_mask, -np.inf)

            result = result.reshape(shape)

    return result


def dispatch_missing(op, left, right, result):
    """
    Fill nulls caused by division by zero, casting to a different dtype
    if necessary.

    Parameters
    ----------
    op : function (operator.add, operator.div, ...)
    left : object (Index for non-reversed ops)
    right : object (Index fof reversed ops)
    result : ndarray

    Returns
    -------
    result : ndarray
    """
    opstr = "__{opname}__".format(opname=op.__name__).replace("____", "__")
    if op in [operator.truediv, operator.floordiv, getattr(operator, "div", None)]:
        result = mask_zero_div_zero(left, right, result)
    elif op is operator.mod:
        result = fill_zeros(result, left, right, opstr, np.nan)
    elif op is divmod:
        res0 = mask_zero_div_zero(left, right, result[0])
        res1 = fill_zeros(result[1], left, right, opstr, np.nan)
        result = (res0, res1)
    return result
