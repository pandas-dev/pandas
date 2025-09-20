from __future__ import annotations

import functools
from collections.abc import Callable

import numpy as np
import pandas as pd
from pandas.api.types import is_bool_dtype
from pandas.errors import IndexingError

from dask._task_spec import Alias, DataNode, Task, TaskRef, convert_legacy_graph
from dask.array import Array
from dask.dataframe import methods
from dask.dataframe.dask_expr import from_dask_array
from dask.dataframe.dask_expr._collection import Series, new_collection
from dask.dataframe.dask_expr._expr import (
    Blockwise,
    Expr,
    MaybeAlignPartitions,
    Partitions,
    Projection,
    are_co_aligned,
    plain_column_projection,
)
from dask.dataframe.dispatch import meta_nonempty
from dask.dataframe.indexing import (
    _coerce_loc_index,
    _maybe_partial_time_string,
    _partition_of_index_value,
    _partitions_of_index_values,
)
from dask.dataframe.utils import is_scalar
from dask.typing import Key
from dask.utils import is_arraylike, is_series_like


class Indexer:
    def __init__(self, obj):
        self.obj = obj


class ILocIndexer(Indexer):
    def __getitem__(self, key):
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

        if len(self.obj.columns) == len(set(self.obj.columns)):
            col_names = self.obj.columns[cindexer]
            if not is_scalar(col_names):
                col_names = list(col_names)
            return new_collection(Projection(self.obj, col_names))
        else:
            raise NotImplementedError


class LocIndexer(Indexer):
    def __getitem__(self, key):
        if isinstance(key, tuple):
            if len(key) > self.obj.ndim:
                raise IndexingError("Too many indexers")

            iindexer = key[0]
            cindexer = key[1]
            pd_loc = self.obj._meta.loc[:, cindexer]
            if isinstance(pd_loc, pd.DataFrame):
                cindexer = list(pd_loc.columns)
        else:
            iindexer = key
            cindexer = None

        if isinstance(cindexer, np.generic):
            cindexer = cindexer.item()

        return self._loc(iindexer, cindexer)

    def _loc(self, iindexer, cindexer):
        if iindexer is None or isinstance(iindexer, slice) and iindexer == slice(None):
            if not isinstance(cindexer, Callable):
                return new_collection(Projection(self.obj, cindexer))
        if isinstance(iindexer, Series):
            return self._loc_series(iindexer, cindexer)
        elif isinstance(iindexer, Array):
            return self._loc_array(iindexer, cindexer)
        elif callable(iindexer):
            return self._loc(iindexer(self.obj), cindexer)

        if self.obj.known_divisions:
            idx = self.obj.index._meta
            unit = idx.unit if hasattr(idx, "unit") else None
            iindexer = self._maybe_partial_time_string(iindexer, unit=unit)

            if isinstance(iindexer, slice):
                return self._loc_slice(iindexer, cindexer)
            elif is_series_like(iindexer) and not is_bool_dtype(iindexer.dtype):
                return new_collection(LocList(self.obj, iindexer.values, cindexer))
            elif isinstance(iindexer, list) or is_arraylike(iindexer):
                if len(iindexer) == 0:
                    return new_collection(LocEmpty(self.obj._meta, cindexer))
                return new_collection(LocList(self.obj, iindexer, cindexer))
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

            return new_collection(LocUnknown(self.obj, iindexer, cindexer))

    def _loc_series(self, iindexer, cindexer, check_alignment=True):
        if not is_bool_dtype(iindexer.dtype):
            raise KeyError(
                "Cannot index with non-boolean dask Series. Try passing computed "
                "values instead (e.g. ``ddf.loc[iindexer.compute()]``)"
            )
        frame = self.obj.expr
        if cindexer is not None:
            frame = Projection(frame, cindexer)
        if check_alignment and not are_co_aligned(frame, iindexer.expr):
            return new_collection(LocAlign(frame, iindexer))
        return new_collection(Loc(frame, iindexer))

    def _loc_array(self, iindexer, cindexer):
        iindexer_series = from_dask_array(iindexer, columns="_", index=self.obj.index)
        return self._loc_series(iindexer_series, cindexer, check_alignment=False)

    def _maybe_partial_time_string(self, iindexer, unit):
        """
        Convert index-indexer for partial time string slicing
        if obj.index is DatetimeIndex / PeriodIndex
        """
        idx = meta_nonempty(self.obj._meta.index)
        iindexer = _maybe_partial_time_string(idx, iindexer, unit)
        return iindexer

    def _loc_slice(self, iindexer, cindexer):
        assert isinstance(iindexer, slice)
        if iindexer.step is not None and iindexer.step < 0:
            obj = ReverseDataFrame(self.obj)
            iindexer = slice(iindexer.start, iindexer.stop, -iindexer.step)
        else:
            obj = self.obj
        assert iindexer.step in (None, 1)
        return new_collection(LocSlice(obj, iindexer, cindexer))

    def _loc_element(self, iindexer, cindexer):
        if iindexer < self.obj.divisions[0] or iindexer > self.obj.divisions[-1]:
            raise KeyError("the label [%s] is not in the index" % str(iindexer))
        return new_collection(LocElement(self.obj, iindexer, cindexer))


def _reverse_partition(df):
    return df.loc[::-1]


class ReverseDataFrame(Expr):
    _parameters = ["frame"]

    @functools.cached_property
    def _meta(self):
        return self.frame._meta

    def _divisions(self):
        return (None,) * (self.frame.npartitions + 1)

    @functools.cached_property
    def npartitions(self):
        return self.frame.npartitions

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Projection):
            return plain_column_projection(self, parent, dependents)

    def _layer(self) -> dict:
        npartitions = self.npartitions
        return {
            (self._name, i): Task(
                TaskRef((self._name, i)),
                _reverse_partition,
                TaskRef((self.frame._name, npartitions - i - 1)),
            )
            for i in range(self.npartitions)
        }


class LocBase(Blockwise):
    _parameters = ["frame", "iindexer", "cindexer"]
    operation = staticmethod(methods.loc)

    @functools.cached_property
    def _meta(self):
        if self.cindexer is None:
            return self.frame._meta
        else:
            return self.frame._meta.loc[:, self.cindexer]

    @functools.cached_property
    def _layer_cache(self):
        return convert_legacy_graph(self._layer())

    def _task(self, name: Key, index: int) -> Task:
        t = self._layer_cache[(self._name, index)]
        if isinstance(t, Alias):
            return Alias(name, t.target)  # type: ignore
        elif t.key != name:
            return Task(name, lambda x: x, t)
        return t


class LocUnknown(Blockwise):
    _parameters = ["frame", "iindexer", "cindexer"]
    operation = staticmethod(methods.try_loc)


class LocElement(LocBase):
    def _divisions(self):
        return (self.iindexer, self.iindexer)

    def _lower(self):
        if self.frame.npartitions == 1:
            return

        part = _get_partitions(self.frame, self.iindexer)
        return type(self)(Partitions(self.frame, [part]), self.iindexer, self.cindexer)

    def _layer(self) -> dict:
        part = _get_partitions(self.frame, self.iindexer)
        return {
            (self._name, 0): Task(
                (self._name, 0),
                methods.loc,
                TaskRef((self.frame._name, part)),
                slice(self.iindexer, self.iindexer),
                self.cindexer,
            )
        }


class LocList(LocBase):
    def _lower(self):
        parts = _get_partitions(self.frame, self.iindexer)
        parts = sorted(parts.keys())
        if len(parts) == 0:
            parts = [0]
        if self.frame.npartitions == len(parts):
            return
        return type(self)(Partitions(self.frame, parts), self.iindexer, self.cindexer)

    @functools.cached_property
    def _layer_information(self):
        dsk = {}
        parts = _get_partitions(self.frame, self.iindexer)
        if len(self.iindexer):
            divisions = []
            items = sorted(parts.items())
            for i, (div, indexer) in enumerate(items):
                dsk[self._name, i] = Task(
                    (self._name, i),
                    methods.loc,
                    TaskRef((self.frame._name, div)),
                    indexer,
                    self.cindexer,
                )
                divisions.append(sorted(indexer)[0])
            divisions.append(sorted(items[-1][1])[-1])
            return dsk, divisions
        else:
            divisions = [None, None]
            dsk = {(self._name, 0): DataNode((self._name, 0), self._meta)}
            return dsk, divisions

    def _divisions(self):
        return self._layer_information[1]

    def _layer(self) -> dict:
        return self._layer_information[0]


class LocEmpty(LocList):
    _parameters = ["meta", "cindexer"]

    def _lower(self):
        return None

    @functools.cached_property
    def _meta(self):
        if self.cindexer is None:
            return self.operand("meta")
        else:
            return self.operand("meta").loc[:, self.cindexer]

    @functools.cached_property
    def _layer_information(self):
        divisions = [None, None]
        dsk = {(self._name, 0): DataNode((self._name, 0), self._meta)}
        return dsk, divisions


class LocSlice(LocBase):
    _projection_passthrough = True

    @functools.cached_property
    def start(self):
        if self.iindexer.start is not None:
            start = _get_partitions(self.frame, self.iindexer.start)
        else:
            start = 0
        return start

    @functools.cached_property
    def stop(self):
        if self.iindexer.stop is not None:
            stop = _get_partitions(self.frame, self.iindexer.stop)
        else:
            stop = self.frame.npartitions - 1
        return stop

    @functools.cached_property
    def istart(self):
        if self.iindexer.start is None and self.frame.known_divisions:
            istart = (
                self.frame.divisions[0]
                if self.iindexer.stop is None
                else min(self.frame.divisions[0], self.iindexer.stop)
            )
        else:
            istart = coerce_loc_index(self.frame, self.iindexer.start)
        return istart

    @functools.cached_property
    def istop(self):
        if self.iindexer.stop is None and self.frame.known_divisions:
            istop = (
                self.frame.divisions[-1]
                if self.iindexer.start is None
                else max(self.frame.divisions[-1], self.iindexer.start)
            )
        else:
            istop = coerce_loc_index(self.frame, self.iindexer.stop)
        return istop

    def _divisions(self):
        if self.stop == self.start:
            return (self.istart, self.istop)

        if self.iindexer.start is None:
            div_start = self.frame.divisions[0]
        else:
            div_start = max(self.istart, self.frame.divisions[self.start])

        if self.iindexer.stop is None:
            div_stop = self.frame.divisions[-1]
        else:
            div_stop = min(self.istop, self.frame.divisions[self.stop + 1])

        return (
            (div_start,)
            + self.frame.divisions[self.start + 1 : self.stop + 1]
            + (div_stop,)
        )

    def _lower(self):
        parts = list(range(self.start, self.stop + 1))
        if self.frame.npartitions == len(parts):
            return
        return type(self)(Partitions(self.frame, parts), self.iindexer, self.cindexer)

    def _layer(self) -> dict:
        if self.stop == self.start:
            return {
                (self._name, 0): Task(
                    (self._name, 0),
                    methods.loc,
                    TaskRef((self.frame._name, self.start)),
                    slice(self.iindexer.start, self.iindexer.stop),
                    self.cindexer,
                )
            }

        dsk = {
            (self._name, 0): Task(
                (self._name, 0),
                methods.loc,
                TaskRef((self.frame._name, self.start)),
                slice(self.iindexer.start, None),
                self.cindexer,
            )
        }
        for i in range(1, self.stop - self.start):
            if self.cindexer is None:
                dsk[self._name, i] = Alias((self.frame._name, self.start + i))  # type: ignore
            else:
                dsk[self._name, i] = Task(
                    (self._name, i),
                    methods.loc,
                    TaskRef((self.frame._name, self.start + i)),
                    slice(None, None),
                    self.cindexer,
                )

        dsk[self._name, self.stop - self.start] = Task(
            (self._name, self.stop - self.start),
            methods.loc,
            TaskRef((self.frame._name, self.stop)),
            slice(None, self.iindexer.stop),
            self.cindexer,
        )
        return dsk


class Loc(Blockwise):
    _parameters = ["frame", "iindexer", "cindexer"]
    _defaults = {"cindexer": None}
    operation = staticmethod(methods.loc)


class LocAlign(MaybeAlignPartitions):
    _parameters = ["frame", "iindexer", "cindexer"]
    _defaults = {"cindexer": None}
    _expr_cls = Loc


def coerce_loc_index(obj, key):
    return _coerce_loc_index(obj.divisions, key)


def _get_partitions(obj, keys):
    if isinstance(keys, list) or is_arraylike(keys):
        return _partitions_of_index_values(obj.divisions, keys)
    else:
        # element
        return _partition_of_index_value(obj.divisions, keys)
