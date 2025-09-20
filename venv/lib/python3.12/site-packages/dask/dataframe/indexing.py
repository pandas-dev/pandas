from __future__ import annotations

import bisect
from collections import defaultdict
from datetime import datetime

import numpy as np
import pandas as pd

from dask.dataframe import methods
from dask.dataframe._compat import PANDAS_GE_300
from dask.dataframe.utils import is_index_like


def _partition_of_index_value(divisions, val):
    """In which partition does this value lie?

    >>> _partition_of_index_value([0, 5, 10], 3)
    0
    >>> _partition_of_index_value([0, 5, 10], 8)
    1
    >>> _partition_of_index_value([0, 5, 10], 100)
    1
    >>> _partition_of_index_value([0, 5, 10], 5)  # left-inclusive divisions
    1
    """
    if divisions[0] is None:
        msg = "Can not use loc on DataFrame without known divisions"
        raise ValueError(msg)
    val = _coerce_loc_index(divisions, val)
    i = bisect.bisect_right(divisions, val)
    return min(len(divisions) - 2, max(0, i - 1))


def _partitions_of_index_values(divisions, values):
    """Return defaultdict of division and values pairs
    Each key corresponds to the division which values are index values belong
    to the division.

    >>> sorted(_partitions_of_index_values([0, 5, 10], [3]).items())
    [(0, [3])]
    >>> sorted(_partitions_of_index_values([0, 5, 10], [3, 8, 5]).items())
    [(0, [3]), (1, [8, 5])]
    """
    if divisions[0] is None:
        msg = "Can not use loc on DataFrame without known divisions"
        raise ValueError(msg)

    results = defaultdict(list)

    # values here can be lots of things, including arrays, lists,
    # and extension arrays like Categorical. Not all of these
    # are supported by tolist, so we ignore the ones we can't handle.
    try:
        tolist = methods.tolist_dispatch.dispatch(type(values))
    except TypeError:
        pass
    else:
        values = tolist(values)

    for val in values:
        i = bisect.bisect_right(divisions, val)
        div = min(len(divisions) - 2, max(0, i - 1))
        results[div].append(val)
    return results


def _coerce_loc_index(divisions, o):
    """Transform values to be comparable against divisions

    This is particularly valuable to use with pandas datetimes
    """
    if divisions and isinstance(divisions[0], datetime):
        return pd.Timestamp(o)
    if divisions and isinstance(divisions[0], np.datetime64):
        return np.datetime64(o).astype(divisions[0].dtype)
    return o


def _maybe_partial_time_string(index, indexer, unit="ns"):
    """
    Convert indexer for partial string selection
    if data has DatetimeIndex/PeriodIndex
    """
    # do not pass dd.Index
    assert is_index_like(index)
    unit = unit or "ns"

    if not isinstance(index, (pd.DatetimeIndex, pd.PeriodIndex)):
        return indexer

    if isinstance(indexer, slice):
        if isinstance(indexer.start, str):
            start = index._maybe_cast_slice_bound(indexer.start, "left")
        else:
            start = indexer.start

        if isinstance(indexer.stop, str):
            stop = index._maybe_cast_slice_bound(indexer.stop, "right")
        else:
            stop = indexer.stop
        if PANDAS_GE_300 and hasattr(start, "as_unit"):
            start = None if start is None else start.as_unit(unit)
            stop = None if stop is None else stop.as_unit(unit)
        return slice(start, stop)

    elif isinstance(indexer, str):
        start = index._maybe_cast_slice_bound(indexer, "left")
        stop = index._maybe_cast_slice_bound(indexer, "right")
        if PANDAS_GE_300 and hasattr(start, "as_unit"):
            start = None if start is None else start.as_unit(unit)
            stop = None if stop is None else stop.as_unit(unit)
        return slice(min(start, stop), max(start, stop))

    return indexer
