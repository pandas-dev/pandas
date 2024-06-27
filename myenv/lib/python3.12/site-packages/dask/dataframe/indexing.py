from __future__ import annotations

import bisect
from collections import defaultdict
from datetime import datetime

import numpy as np
import pandas as pd
from pandas.api.types import is_bool_dtype

from dask.array.core import Array
from dask.base import tokenize
from dask.dataframe import methods
from dask.dataframe._compat import IndexingError
from dask.dataframe.core import Series, new_dd_object
from dask.dataframe.utils import is_index_like, is_series_like, meta_nonempty
from dask.highlevelgraph import HighLevelGraph
from dask.utils import is_arraylike


class _IndexerBase:
    def __init__(self, obj):
        self.obj = obj

    @property
    def _name(self):
        return self.obj._name

    @property
    def _meta_indexer(self):
        raise NotImplementedError

    def _make_meta(self, iindexer, cindexer):
        """
        get metadata
        """
        if cindexer is None:
            return self.obj
        else:
            return self._meta_indexer[:, cindexer]

    def __dask_tokenize__(self):
        return type(self).__name__, tokenize(self.obj)


class _iLocIndexer(_IndexerBase):
    @property
    def _meta_indexer(self):
        return self.obj._meta.iloc

    def __getitem__(self, key):
        # dataframe
        msg = (
            "'DataFrame.iloc' only supports selecting columns. "
            "It must be used like 'df.iloc[:, column_indexer]'."
        )
        if not isinstance(key, tuple):
            raise NotImplementedError(msg)

        if len(key) > 2:
            raise ValueError("Too many indexers")

        iindexer, cindexer = key

        if iindexer != slice(None):
            raise NotImplementedError(msg)

        if not self.obj.columns.is_unique:
            # if there are any duplicate column names, do an iloc
            return self._iloc(iindexer, cindexer)
        else:
            # otherwise dispatch to dask.dataframe.core.DataFrame.__getitem__
            col_names = self.obj.columns[cindexer]
            return self.obj.__getitem__(col_names)

    def _iloc(self, iindexer, cindexer):
        assert iindexer == slice(None)
        meta = self._make_meta(iindexer, cindexer)

        return self.obj.map_partitions(methods.iloc, cindexer, meta=meta)


class _LocIndexer(_IndexerBase):
    """Helper class for the .loc accessor"""

    @property
    def _meta_indexer(self):
        return self.obj._meta.loc

    def __getitem__(self, key):
        if isinstance(key, tuple):
            # multi-dimensional selection
            if len(key) > self.obj.ndim:
                # raise from pandas
                msg = "Too many indexers"
                raise IndexingError(msg)

            iindexer = key[0]
            cindexer = key[1]
        else:
            # if self.obj is Series, cindexer is always None
            iindexer = key
            cindexer = None
        return self._loc(iindexer, cindexer)

    def _loc(self, iindexer, cindexer):
        """Helper function for the .loc accessor"""
        if isinstance(iindexer, Series):
            return self._loc_series(iindexer, cindexer)
        elif isinstance(iindexer, Array):
            return self._loc_array(iindexer, cindexer)
        elif callable(iindexer):
            return self._loc(iindexer(self.obj), cindexer)

        if self.obj.known_divisions:
            iindexer = self._maybe_partial_time_string(iindexer)

            if isinstance(iindexer, slice):
                return self._loc_slice(iindexer, cindexer)
            elif is_series_like(iindexer) and not is_bool_dtype(iindexer.dtype):
                return self._loc_list(iindexer.values, cindexer)
            elif isinstance(iindexer, list) or is_arraylike(iindexer):
                return self._loc_list(iindexer, cindexer)
            else:
                # element should raise KeyError
                return self._loc_element(iindexer, cindexer)
        else:
            if isinstance(iindexer, (list, np.ndarray)) or (
                is_series_like(iindexer) and not is_bool_dtype(iindexer.dtype)
            ):
                # applying map_partitions to each partition
                # results in duplicated NaN rows
                msg = (
                    "Cannot index with list against unknown division. "
                    "Try setting divisions using ``ddf.set_index``"
                )
                raise KeyError(msg)
            elif not isinstance(iindexer, slice):
                iindexer = slice(iindexer, iindexer)

            meta = self._make_meta(iindexer, cindexer)
            return self.obj.map_partitions(
                methods.try_loc, iindexer, cindexer, meta=meta
            )

    def _maybe_partial_time_string(self, iindexer):
        """
        Convert index-indexer for partial time string slicing
        if obj.index is DatetimeIndex / PeriodIndex
        """
        idx = meta_nonempty(self.obj._meta.index)
        iindexer = _maybe_partial_time_string(idx, iindexer)
        return iindexer

    def _loc_series(self, iindexer, cindexer):
        if not is_bool_dtype(iindexer.dtype):
            raise KeyError(
                "Cannot index with non-boolean dask Series. Try passing computed "
                "values instead (e.g. ``ddf.loc[iindexer.compute()]``)"
            )
        meta = self._make_meta(iindexer, cindexer)
        return self.obj.map_partitions(
            methods.loc, iindexer, cindexer, token="loc-series", meta=meta
        )

    def _loc_array(self, iindexer, cindexer):
        iindexer_series = iindexer.to_dask_dataframe("_", self.obj.index)
        return self._loc_series(iindexer_series, cindexer)

    def _loc_list(self, iindexer, cindexer):
        name = "loc-%s" % tokenize(iindexer, self.obj)
        parts = self._get_partitions(iindexer)
        meta = self._make_meta(iindexer, cindexer)

        if len(iindexer):
            dsk = {}
            divisions = []
            items = sorted(parts.items())
            for i, (div, indexer) in enumerate(items):
                dsk[name, i] = (methods.loc, (self._name, div), indexer, cindexer)
                # append minimum value as division
                divisions.append(sorted(indexer)[0])
            # append maximum value of the last division
            divisions.append(sorted(items[-1][1])[-1])
            graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self.obj])
        else:
            divisions = [None, None]
            dsk = {(name, 0): meta.head(0)}
            graph = HighLevelGraph.from_collections(name, dsk)
        return new_dd_object(graph, name, meta=meta, divisions=divisions)

    def _loc_element(self, iindexer, cindexer):
        name = "loc-%s" % tokenize(iindexer, self.obj)
        part = self._get_partitions(iindexer)

        if iindexer < self.obj.divisions[0] or iindexer > self.obj.divisions[-1]:
            raise KeyError("the label [%s] is not in the index" % str(iindexer))

        dsk = {
            (name, 0): (
                methods.loc,
                (self._name, part),
                slice(iindexer, iindexer),
                cindexer,
            )
        }

        meta = self._make_meta(iindexer, cindexer)
        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self.obj])
        return new_dd_object(graph, name, meta=meta, divisions=[iindexer, iindexer])

    def _get_partitions(self, keys):
        if isinstance(keys, list) or is_arraylike(keys):
            return _partitions_of_index_values(self.obj.divisions, keys)
        else:
            # element
            return _partition_of_index_value(self.obj.divisions, keys)

    def _coerce_loc_index(self, key):
        return _coerce_loc_index(self.obj.divisions, key)

    def _loc_slice(self, iindexer, cindexer):
        name = "loc-%s" % tokenize(iindexer, cindexer, self)

        assert isinstance(iindexer, slice)
        assert iindexer.step in (None, 1)

        if iindexer.start is not None:
            start = self._get_partitions(iindexer.start)
        else:
            start = 0
        if iindexer.stop is not None:
            stop = self._get_partitions(iindexer.stop)
        else:
            stop = self.obj.npartitions - 1

        if iindexer.start is None and self.obj.known_divisions:
            istart = (
                self.obj.divisions[0]
                if iindexer.stop is None
                else min(self.obj.divisions[0], iindexer.stop)
            )
        else:
            istart = self._coerce_loc_index(iindexer.start)
        if iindexer.stop is None and self.obj.known_divisions:
            istop = (
                self.obj.divisions[-1]
                if iindexer.start is None
                else max(self.obj.divisions[-1], iindexer.start)
            )
        else:
            istop = self._coerce_loc_index(iindexer.stop)

        if stop == start:
            dsk = {
                (name, 0): (
                    methods.loc,
                    (self._name, start),
                    slice(iindexer.start, iindexer.stop),
                    cindexer,
                )
            }
            divisions = [istart, istop]
        else:
            dsk = {
                (name, 0): (
                    methods.loc,
                    (self._name, start),
                    slice(iindexer.start, None),
                    cindexer,
                )
            }
            for i in range(1, stop - start):
                if cindexer is None:
                    dsk[name, i] = (self._name, start + i)
                else:
                    dsk[name, i] = (
                        methods.loc,
                        (self._name, start + i),
                        slice(None, None),
                        cindexer,
                    )

            dsk[name, stop - start] = (
                methods.loc,
                (self._name, stop),
                slice(None, iindexer.stop),
                cindexer,
            )

            if iindexer.start is None:
                div_start = self.obj.divisions[0]
            else:
                div_start = max(istart, self.obj.divisions[start])

            if iindexer.stop is None:
                div_stop = self.obj.divisions[-1]
            else:
                div_stop = min(istop, self.obj.divisions[stop + 1])

            divisions = (
                (div_start,) + self.obj.divisions[start + 1 : stop + 1] + (div_stop,)
            )

        assert len(divisions) == len(dsk) + 1

        meta = self._make_meta(iindexer, cindexer)
        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self.obj])
        return new_dd_object(graph, name, meta=meta, divisions=divisions)


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


def _maybe_partial_time_string(index, indexer):
    """
    Convert indexer for partial string selection
    if data has DatetimeIndex/PeriodIndex
    """
    # do not pass dd.Index
    assert is_index_like(index)

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
        return slice(start, stop)

    elif isinstance(indexer, str):
        start = index._maybe_cast_slice_bound(indexer, "left")
        stop = index._maybe_cast_slice_bound(indexer, "right")
        return slice(min(start, stop), max(start, stop))

    return indexer
