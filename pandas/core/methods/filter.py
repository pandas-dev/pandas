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
    ABCSeries,
)

from pandas.core.construction import array as pd_array

if TYPE_CHECKING:
    from pandas._typing import (
        AnyArrayLike,
        AxisInt,
        NDFrameT,
    )

    from pandas import (
        Index,
        Series,
    )


def is_mask(key: object) -> bool:
    """
    Whether ``key`` is a boolean mask, possibly holding missing values.
    """
    if isinstance(key, list):
        if len(key) == 0:
            return False
        return lib.is_bool_array(np.asarray(key, dtype=object), skipna=True)
    if isinstance(key, ABCDataFrame):
        # rejected as not one-dimensional by filter_mask
        return True
    dtype = getattr(key, "dtype", None)
    if dtype is None:
        return False
    if dtype == np.object_:
        return lib.is_bool_array(np.asarray(key), skipna=True)
    return is_bool_dtype(dtype)


def has_bool_labels(labels: Index) -> bool:
    """
    Whether ``labels`` contains the values True or False.
    """
    if labels.dtype == bool:
        return True
    if labels.dtype != np.object_:
        return False
    return any(isinstance(label, (bool, np.bool_)) for label in labels)


def filter_mask(
    obj: NDFrameT,
    mask: list | AnyArrayLike,
    axis: AxisInt,
    na: Literal["raise"] | bool,
) -> NDFrameT:
    """
    Select the entries of ``obj`` along ``axis`` where ``mask`` is True.
    """
    from pandas.core.indexing import check_bool_indexer

    labels = obj._get_axis(axis)

    if getattr(mask, "ndim", 1) != 1:
        raise ValueError(
            f"The mask passed to {type(obj).__name__}.filter must be one-dimensional"
        )

    values = mask._values if isinstance(mask, ABCSeries) else mask
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
