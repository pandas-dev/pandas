"""
Implementation of DataFrame.filter and Series.filter with a boolean mask.
"""

from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Literal,
)

import numpy as np

from pandas._libs import lib

from pandas.core.dtypes.common import is_bool_dtype
from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCMultiIndex,
    ABCSeries,
)
from pandas.core.dtypes.missing import isna

from pandas.core.construction import (
    array as pd_array,
    extract_array,
)
from pandas.core.indexing import check_bool_indexer

if TYPE_CHECKING:
    from pandas._typing import (
        AnyArrayLike,
        AxisInt,
        NDFrameT,
    )

    from pandas import Series


def is_mask(key: object) -> bool:
    """
    Whether ``key`` is a boolean mask, possibly holding missing values.

    A boolean dtype is a mask. Otherwise ``key`` must be a list or an
    object-dtype array-like that is one-dimensional with every non-missing
    element a bool. A tuple is never a mask.
    """
    if isinstance(key, ABCDataFrame):
        # A boolean DataFrame (e.g. df > 1) is an attempted mask, so classify
        # it as one to get filter_mask's "must be one-dimensional" error
        # rather than the generic error a non-mask gets.
        return all(is_bool_dtype(dtype) for dtype in key.dtypes)
    if isinstance(key, list):
        values = np.asarray(key, dtype=object)
    else:
        dtype = getattr(key, "dtype", None)
        if dtype is None or isinstance(key, ABCMultiIndex):
            return False
        if is_bool_dtype(dtype):
            return True
        if dtype != np.object_:
            return False
        values = np.asarray(key, dtype=object)
    if values.ndim != 1:
        # lib.is_bool_array iterates a 2-D array flat, so a list of bool tuples
        # (labels for a MultiIndex) would otherwise be mistaken for a mask.
        return False
    return lib.is_bool_array(values, skipna=True)


def resembles_mask(key: object) -> bool:
    """
    Whether ``key``, passed positionally, was likely intended as a mask.

    This only decides whether to warn; ``key`` selects labels regardless. At
    least one boolean is required so that a list of missing labels such as
    ``[np.nan]`` does not warn. A DataFrame is excluded because the label path
    raises for it anyway, and a warning before an error is just noise.
    """
    if isinstance(key, ABCDataFrame) or not is_mask(key):
        return False
    if not isinstance(key, list) and is_bool_dtype(key.dtype):
        return True
    return not isna(np.asarray(key, dtype=object)).all()


def filter_mask(
    obj: NDFrameT,
    mask: list | AnyArrayLike,
    axis: AxisInt,
    na: Literal["raise"] | bool,
) -> NDFrameT:
    """
    Select the entries of ``obj`` along ``axis`` where ``mask`` is True.
    """
    labels = obj._get_axis(axis)

    if getattr(mask, "ndim", 1) != 1:
        raise ValueError(
            f"The mask passed to {type(obj).__name__}.filter must be one-dimensional"
        )

    values = extract_array(mask, extract_numpy=True)
    if isinstance(values, list):
        values = np.asarray(values, dtype=object)
    if isinstance(values, np.ndarray) and values.dtype == np.bool_:
        # A NumPy bool array cannot hold missing values, so skip the copy
        # through BooleanArray that the NA handling below requires.
        np_mask = values
    else:
        values = pd_array(values, dtype="boolean")
        if values.isna().any():
            if na == "raise":
                raise ValueError("The mask contains missing values")
            values = values.fillna(na)
        np_mask = values.to_numpy(dtype=bool)
    key: Series | np.ndarray
    if isinstance(mask, ABCSeries):
        key = mask._constructor(np_mask, index=mask.index, copy=False)
    else:
        key = np_mask
    indexer = check_bool_indexer(labels, key)
    return obj.loc(axis=axis)[indexer]
