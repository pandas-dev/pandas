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

    from pandas import (
        Index,
        Series,
    )


def is_mask(key: object, labels: Index) -> bool:
    """
    Whether ``key`` is a boolean mask rather than a list-like of labels.

    A boolean dtype is always a mask. Otherwise ``key`` is a mask when it is
    one-dimensional and every non-missing element is a bool. A tuple is always
    a sequence of labels. ``labels`` is the axis that label-based selection
    would use; it only matters when ``key`` consists entirely of missing
    values, which is a list of labels when ``labels`` contains a missing value
    and an uninformative mask otherwise.
    """
    if isinstance(key, ABCDataFrame):
        # a boolean DataFrame is rejected as not one-dimensional by filter_mask
        return all(is_bool_dtype(dtype) for dtype in key.dtypes)
    if isinstance(key, list):
        if len(key) == 0:
            return False
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
        # e.g. a list of tuples selecting labels from a MultiIndex
        return False
    if not lib.is_bool_array(values, skipna=True):
        return False
    if isna(values).all():
        return not labels.hasnans
    return True


def has_bool_labels(labels: Index) -> bool:
    """
    Whether ``labels`` contains the values True or False.
    """
    if isinstance(labels, ABCMultiIndex):
        # labels are tuples
        return False
    if is_bool_dtype(labels.dtype):
        return True
    if labels.dtype != np.object_:
        return False
    # inferred_type is cached on the Index, avoiding the loop below on
    # repeated calls with a homogeneous object-dtype axis
    inferred = labels.inferred_type
    if inferred == "boolean":
        return True
    if not inferred.startswith("mixed"):
        return False
    return any(lib.is_bool(label) for label in labels)


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
        # fast path: no missing values are possible
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
