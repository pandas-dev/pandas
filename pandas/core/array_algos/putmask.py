"""
EA-compatible analogue to to np.putmask
"""
from __future__ import annotations

from typing import Any
import warnings

import numpy as np

from pandas._libs import lib
from pandas._typing import (
    ArrayLike,
    npt,
)

from pandas.core.dtypes.cast import (
    can_hold_element,
    convert_scalar_for_putitemlike,
    find_common_type,
    infer_dtype_from,
)
from pandas.core.dtypes.common import (
    is_float_dtype,
    is_integer_dtype,
    is_interval_dtype,
    is_list_like,
)
from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCIndex,
    ABCSeries,
)
from pandas.core.dtypes.missing import (
    is_valid_na_for_dtype,
    isna_compat,
    na_value_for_dtype,
)

from pandas.core.arrays import ExtensionArray
from pandas.core.arrays._mixins import NDArrayBackedExtensionArray


def putmask_inplace(values: ArrayLike, mask: npt.NDArray[np.bool_], value: Any) -> None:
    """
    ExtensionArray-compatible implementation of np.putmask.  The main
    difference is we do not handle repeating or truncating like numpy.

    Parameters
    ----------
    values: np.ndarray or ExtensionArray
    mask : np.ndarray[bool]
        We assume extract_bool_array has already been called.
    value : Any
    """

    if lib.is_scalar(value) and isinstance(values, np.ndarray):
        value = convert_scalar_for_putitemlike(value, values.dtype)

    if (
        not isinstance(values, np.ndarray)
        or (values.dtype == object and not lib.is_scalar(value))
        # GH#43424: np.putmask raises TypeError if we cannot cast between types with
        # rule = "safe", a stricter guarantee we may not have here
        or (
            isinstance(value, np.ndarray) and not np.can_cast(value.dtype, values.dtype)
        )
    ):
        # GH#19266 using np.putmask gives unexpected results with listlike value
        #  along with object dtype
        if is_list_like(value) and len(value) == len(values):
            values[mask] = value[mask]
        else:
            values[mask] = value
    else:
        # GH#37833 np.putmask is more performant than __setitem__
        np.putmask(values, mask, value)


def putmask_smart(values: np.ndarray, mask: npt.NDArray[np.bool_], new) -> np.ndarray:
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
    np.putmask
    """
    # we cannot use np.asarray() here as we cannot have conversions
    # that numpy does when numeric are mixed with strings

    if not is_list_like(new):
        new = np.broadcast_to(new, mask.shape)

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

    # new = np.asarray(new)

    if values.dtype.kind == new.dtype.kind:
        # preserves dtype if possible
        np.putmask(values, mask, new)
        return values

    dtype = find_common_type([values.dtype, new.dtype])
    # error: Argument 1 to "astype" of "_ArrayOrScalarCommon" has incompatible type
    # "Union[dtype[Any], ExtensionDtype]"; expected "Union[dtype[Any], None, type,
    # _SupportsDType, str, Union[Tuple[Any, int], Tuple[Any, Union[int, Sequence[int]]],
    # List[Any], _DTypeDict, Tuple[Any, Any]]]"
    values = values.astype(dtype)  # type: ignore[arg-type]

    np.putmask(values, mask, new)
    return values


def putmask_without_repeat(
    values: np.ndarray, mask: npt.NDArray[np.bool_], new: Any
) -> None:
    """
    np.putmask will truncate or repeat if `new` is a listlike with
    len(new) != len(values).  We require an exact match.

    Parameters
    ----------
    values : np.ndarray
    mask : np.ndarray[bool]
    new : Any
    """
    if getattr(new, "ndim", 0) >= 1:
        new = new.astype(values.dtype, copy=False)

    # TODO: this prob needs some better checking for 2D cases
    nlocs = mask.sum()
    if nlocs > 0 and is_list_like(new) and getattr(new, "ndim", 1) == 1:
        if nlocs == len(new):
            # GH#30567
            # If length of ``new`` is less than the length of ``values``,
            # `np.putmask` would first repeat the ``new`` array and then
            # assign the masked values hence produces incorrect result.
            # `np.place` on the other hand uses the ``new`` values at it is
            # to place in the masked locations of ``values``
            np.place(values, mask, new)
            # i.e. values[mask] = new
        elif mask.shape[-1] == len(new) or len(new) == 1:
            np.putmask(values, mask, new)
        else:
            raise ValueError("cannot assign mismatch length to masked array")
    else:
        np.putmask(values, mask, new)


def validate_putmask(
    values: ArrayLike, mask: np.ndarray
) -> tuple[npt.NDArray[np.bool_], bool]:
    """
    Validate mask and check if this putmask operation is a no-op.
    """
    mask = extract_bool_array(mask)
    if mask.shape != values.shape:
        raise ValueError("putmask: mask and data must be the same size")

    noop = not mask.any()
    return mask, noop


def extract_bool_array(mask: ArrayLike) -> npt.NDArray[np.bool_]:
    """
    If we have a SparseArray or BooleanArray, convert it to ndarray[bool].
    """
    if isinstance(mask, ExtensionArray):
        # We could have BooleanArray, Sparse[bool], ...
        #  Except for BooleanArray, this is equivalent to just
        #  np.asarray(mask, dtype=bool)
        mask = mask.to_numpy(dtype=bool, na_value=False)

    mask = np.asarray(mask, dtype=bool)
    return mask


def setitem_datetimelike_compat(values: np.ndarray, num_set: int, other):
    """
    Parameters
    ----------
    values : np.ndarray
    num_set : int
        For putmask, this is mask.sum()
    other : Any
    """
    if values.dtype == object:
        dtype, _ = infer_dtype_from(other, pandas_dtype=True)

        if isinstance(dtype, np.dtype) and dtype.kind in ["m", "M"]:
            # https://github.com/numpy/numpy/issues/12550
            #  timedelta64 will incorrectly cast to int
            if not is_list_like(other):
                other = [other] * num_set
            else:
                other = list(other)

    return other


def putmask_flexible_ndarray(array: np.ndarray, mask, new):
    """
    Putmask implementation for ArrayManager putmask for ndarray.

    Flexible version that will upcast if needed.
    """
    mask, noop = validate_putmask(array, mask)
    assert not isinstance(new, (ABCIndex, ABCSeries, ABCDataFrame))

    # if we are passed a scalar None, convert it here
    if not array.dtype == "object" and is_valid_na_for_dtype(new, array.dtype):
        new = na_value_for_dtype(array.dtype, compat=False)

    if can_hold_element(array, new):
        # error: Argument 1 to "putmask_without_repeat" has incompatible type
        # "Union[ndarray, ExtensionArray]"; expected "ndarray"
        putmask_without_repeat(array, mask, new)  # type: ignore[arg-type]
        return array

    elif noop:
        return array

    dtype, _ = infer_dtype_from(new)
    if dtype.kind in ["m", "M"]:
        array = array.astype(object)
        # convert to list to avoid numpy coercing datetimelikes to integers
        new = setitem_datetimelike_compat(
            array, mask.sum(), new  # type: ignore[arg-type]
        )
        # putmask_smart below converts it back to array
        np.putmask(array, mask, new)
        return array

    new_values = putmask_smart(array, mask, new)
    return new_values


def _coerce_to_target_dtype(array, new):
    dtype, _ = infer_dtype_from(new, pandas_dtype=True)
    new_dtype = find_common_type([array.dtype, dtype])
    return array.astype(new_dtype, copy=False)


def putmask_flexible_ea(array: ExtensionArray, mask, new):
    """
    Putmask implementation for ArrayManager putmask for EA.

    Flexible version that will upcast if needed.
    """
    mask = extract_bool_array(mask)

    if isinstance(array, NDArrayBackedExtensionArray):

        if not can_hold_element(array, new):
            array = _coerce_to_target_dtype(array, new)
            return putmask_flexible_ndarray(array, mask, new)

        array.putmask(mask, new)
        return array

    if isinstance(new, (np.ndarray, ExtensionArray)) and len(new) == len(mask):
        new = new[mask]

    try:
        array[mask] = new
    except TypeError:
        if not is_interval_dtype(array.dtype):
            # Discussion about what we want to support in the general
            #  case GH#39584
            raise

        array = _coerce_to_target_dtype(array, new)
        if array.dtype == np.dtype("object"):
            # For now at least, only support casting e.g.
            #  Interval[int64]->Interval[float64],
            raise
        return putmask_flexible_ea(array, mask, new)

    return array
