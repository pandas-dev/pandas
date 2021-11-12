"""
Wrappers around array_algos with internals-specific logic
"""
from __future__ import annotations

import numpy as np

from pandas.core.dtypes.cast import (
    can_hold_element,
    find_common_type,
    infer_dtype_from,
)
from pandas.core.dtypes.common import is_interval_dtype
from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCIndex,
    ABCSeries,
)
from pandas.core.dtypes.missing import (
    is_valid_na_for_dtype,
    na_value_for_dtype,
)

from pandas.core.array_algos.putmask import (
    extract_bool_array,
    putmask_smart,
    putmask_without_repeat,
    setitem_datetimelike_compat,
    validate_putmask,
)
from pandas.core.arrays import ExtensionArray
from pandas.core.arrays._mixins import NDArrayBackedExtensionArray


def putmask_flexible(array: np.ndarray | ExtensionArray, mask, new):
    """
    Putmask implementation for ArrayManager putmask for ndarray.

    Flexible version that will upcast if needed.
    """
    if isinstance(array, np.ndarray):
        return putmask_flexible_ndarray(array, mask=mask, new=new)
    else:
        return putmask_flexible_ea(array, mask=mask, new=new)


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
        putmask_without_repeat(array, mask, new)
        return array

    elif noop:
        return array

    dtype, _ = infer_dtype_from(new)
    if dtype.kind in ["m", "M"]:
        array = array.astype(object)
        # convert to list to avoid numpy coercing datetimelikes to integers
        new = setitem_datetimelike_compat(array, mask.sum(), new)
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
            return putmask_flexible(array, mask, new)

    try:
        array._putmask(mask, new)
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
        return putmask_flexible(array, mask, new)

    return array
