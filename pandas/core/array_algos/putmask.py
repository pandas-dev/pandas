"""
EA-compatible analogue to to np.putmask
"""
from typing import Any
import warnings

import numpy as np

from pandas._libs import lib
from pandas._typing import ArrayLike

from pandas.core.dtypes.cast import convert_scalar_for_putitemlike, maybe_promote
from pandas.core.dtypes.common import is_float_dtype, is_integer_dtype, is_list_like
from pandas.core.dtypes.missing import isna_compat


def putmask_inplace(values: ArrayLike, mask: np.ndarray, value: Any) -> None:
    """
    ExtensionArray-compatible implementation of np.putmask.  The main
    difference is we do not handle repeating or truncating like numpy.

    Parameters
    ----------
    mask : np.ndarray[bool]
        We assume _extract_bool_array has already been called.
    value : Any
    """

    if lib.is_scalar(value) and isinstance(values, np.ndarray):
        value = convert_scalar_for_putitemlike(value, values.dtype)

    if not isinstance(values, np.ndarray) or (
        values.dtype == object and not lib.is_scalar(value)
    ):
        # GH#19266 using np.putmask gives unexpected results with listlike value
        if is_list_like(value) and len(value) == len(values):
            values[mask] = value[mask]
        else:
            values[mask] = value
    else:
        # GH#37833 np.putmask is more performant than __setitem__
        np.putmask(values, mask, value)


def putmask_smart(values: np.ndarray, mask: np.ndarray, new) -> np.ndarray:
    """
    Return a new ndarray, try to preserve dtype if possible.

    Parameters
    ----------
    values : np.ndarray
        `values`, updated in-place.
    mask : np.ndarray[bool]
        Applies to both sides (array like).
    new : `new values` either scalar or an array like aligned with `values`

    Returns
    -------
    values : ndarray with updated values
        this *may* be a copy of the original

    See Also
    --------
    ndarray.putmask
    """
    # we cannot use np.asarray() here as we cannot have conversions
    # that numpy does when numeric are mixed with strings

    # n should be the length of the mask or a scalar here
    if not is_list_like(new):
        new = np.repeat(new, len(mask))

    # see if we are only masking values that if putted
    # will work in the current dtype
    try:
        nn = new[mask]
    except TypeError:
        # TypeError: only integer scalar arrays can be converted to a scalar index
        pass
    else:
        # make sure that we have a nullable type if we have nulls
        if not isna_compat(values, nn[0]):
            pass
        elif not (is_float_dtype(nn.dtype) or is_integer_dtype(nn.dtype)):
            # only compare integers/floats
            pass
        elif not (is_float_dtype(values.dtype) or is_integer_dtype(values.dtype)):
            # only compare integers/floats
            pass
        else:

            # we ignore ComplexWarning here
            with warnings.catch_warnings(record=True):
                warnings.simplefilter("ignore", np.ComplexWarning)
                nn_at = nn.astype(values.dtype)

            comp = nn == nn_at
            if is_list_like(comp) and comp.all():
                nv = values.copy()
                nv[mask] = nn_at
                return nv

    new = np.asarray(new)

    if values.dtype.kind == new.dtype.kind:
        # preserves dtype if possible
        return _putmask_preserve(values, new, mask)

    # change the dtype if needed
    dtype, _ = maybe_promote(new.dtype)

    values = values.astype(dtype)

    return _putmask_preserve(values, new, mask)


def _putmask_preserve(new_values: np.ndarray, new, mask: np.ndarray):
    try:
        new_values[mask] = new[mask]
    except (IndexError, ValueError):
        new_values[mask] = new
    return new_values
