from __future__ import annotations

import functools
import math
import operator

import numpy as np
import pyarrow as pa

from dask import delayed
from dask._task_spec import DataNode, List, Task
from dask.dataframe import methods
from dask.dataframe._pyarrow import to_pyarrow_string
from dask.dataframe.core import apply_and_enforce, is_dataframe_like
from dask.dataframe.dask_expr._expr import (
    Blockwise,
    Expr,
    Lengths,
    Literal,
    PartitionsFiltered,
    Projection,
    determine_column_projection,
    no_default,
)
from dask.dataframe.dask_expr._reductions import Len
from dask.dataframe.dask_expr._util import _BackendData, _convert_to_list
from dask.dataframe.dispatch import make_meta
from dask.dataframe.io.io import _meta_from_array, sorted_division_locations
from dask.tokenize import _tokenize_deterministic
from dask.typing import Key
from dask.utils import funcname, is_series_like


class IO(Expr):
    def __str__(self):
        return f"{type(self).__name__}({self._name[-7:]})"


class FromGraph(IO):
    """A DataFrame created from an opaque Dask task graph

    This is used in persist, for example, and would also be used in any
    conversion from legacy dataframes.
    """

    _parameters = ["layer", "_meta", "divisions", "keys", "name_prefix"]

    @property
    def _meta(self):
        return self.operand("_meta")

    def _divisions(self):
        return self.operand("divisions")

    @functools.cached_property
    def _name(self):
        return (
            self.operand("name_prefix") + "-" + _tokenize_deterministic(*self.operands)
        )

    def _layer(self):
        dsk = dict(self.operand("layer"))
        # The name may not actually match the layers name therefore rewrite this
        # using an alias
        for part, k in enumerate(self.operand("keys")):
            dsk[(self._name, part)] = k
        return dsk


class BlockwiseIO(Blockwise, IO):
    _absorb_projections = False

    @functools.cached_property
    def _fusion_compression_factor(self):
        return 1

    def _simplify_up(self, parent, dependents):
        if (
            self._absorb_projections
            and isinstance(parent, Projection)
            and is_dataframe_like(self._meta)
        ):
            # Column projection
            parent_columns = parent.operand("columns")
            proposed_columns = determine_column_projection(self, parent, dependents)
            proposed_columns = _convert_to_list(proposed_columns)
            proposed_columns = [col for col in self.columns if col in proposed_columns]
            if set(proposed_columns) == set(self.columns):
                # Already projected or nothing to do
                return
            substitutions = {"columns": proposed_columns}
            result = self.substitute_parameters(substitutions)
            if result.columns != parent_columns:
                result = result[parent_columns]
            return result

    def _tune_up(self, parent):
        if self._fusion_compression_factor >= 1:
            return
        if isinstance(parent, FusedIO):
            return
        return parent.substitute(self, FusedIO(self))


class FusedIO(BlockwiseIO):
    _parameters = ["_expr"]

    @functools.cached_property
    def _name(self):
        return self.operand("_expr")._funcname + "-fused-" + self.deterministic_token

    @functools.cached_property
    def _meta(self):
        return self.operand("_expr")._meta

    def dependencies(self):
        return []

    @functools.cached_property
    def npartitions(self):
        return len(self._fusion_buckets)

    def _divisions(self):
        divisions = self.operand("_expr")._divisions()
        new_divisions = [divisions[b[0]] for b in self._fusion_buckets]
        if new_divisions[0] is None:
            new_divisions.append(None)
        else:
            new_divisions.append(divisions[-1])
        return tuple(new_divisions)

    def _task(self, name: Key, index: int) -> Task:
        expr = self.operand("_expr")
        bucket = self._fusion_buckets[index]
        # FIXME: This will likely require a wrapper
        return Task(
            name,
            methods.concat,
            List(*(expr._filtered_task(name, i) for i in bucket)),
            _data_producer=True,
        )

    @functools.cached_property
    def _fusion_buckets(self):
        partitions = self.operand("_expr")._partitions
        npartitions = len(partitions)

        step = math.ceil(1 / self.operand("_expr")._fusion_compression_factor)
        step = min(step, math.ceil(math.sqrt(npartitions)), 100)

        buckets = [partitions[i : i + step] for i in range(0, npartitions, step)]
        return buckets

    def _tune_up(self, parent):
        return


class FusedParquetIO(FusedIO):
    _parameters = ["_expr"]

    @functools.cached_property
    def _name(self):
        return (
            funcname(type(self.operand("_expr"))).lower()
            + "-fused-parq-"
            + self.deterministic_token
        )

    @staticmethod
    def _load_multiple_files(
        frag_filters,
        columns,
        schema,
        **to_pandas_kwargs,
    ):
        from dask.dataframe.dask_expr.io.parquet import ReadParquetPyarrowFS

        tables = (
            ReadParquetPyarrowFS._fragment_to_table(
                frag,
                filter,
                columns,
                schema,
            )
            for frag, filter in frag_filters
        )
        table = pa.concat_tables(tables, promote_options="permissive")
        return ReadParquetPyarrowFS._table_to_pandas(table, **to_pandas_kwargs)

    def _task(self, name: str, index: int) -> Task:  # type: ignore
        expr = self.operand("_expr")
        bucket = self._fusion_buckets[index]
        fragments_filters = []
        assert bucket
        for i in bucket:
            subtask = expr._filtered_task(name, i)
            # This is unique / same for all tasks
            to_pandas_kwargs = subtask.kwargs
            assert len(subtask.args) == 1
            frag_to_table_task = subtask.args[0]
            fragments_filters.append(
                (
                    frag_to_table_task.kwargs["fragment_wrapper"],
                    frag_to_table_task.kwargs["filters"],
                )
            )
            columns = frag_to_table_task.kwargs["columns"]
            schema = frag_to_table_task.kwargs["schema"]
        return Task(
            name,
            self._load_multiple_files,
            fragments_filters,
            columns,
            schema,
            **to_pandas_kwargs,
            _data_producer=True,
        )


class FromMap(PartitionsFiltered, BlockwiseIO):
    _parameters = [
        "func",
        "iterables",
        "args",
        "kwargs",
        "user_meta",
        "enforce_metadata",
        "user_divisions",
        "label",
        "_partitions",
    ]
    _defaults = {
        "user_meta": no_default,
        "enforce_metadata": False,
        "user_divisions": None,
        "label": None,
        "_partitions": None,
    }
    _absorb_projections = False

    @functools.cached_property
    def _name(self):
        if self.label is None:
            return funcname(self.func).lower() + "-" + self.deterministic_token
        else:
            return self.label + "-" + self.deterministic_token

    @functools.cached_property
    def _meta(self):
        if self.operand("user_meta") is not no_default:
            meta = self.operand("user_meta")
            return make_meta(meta)
        else:
            vals = [v[0] for v in self.iterables]
            meta = delayed(self.func)(*vals, *self.args, **self.kwargs)
            return delayed(make_meta)(meta).compute()

    def _divisions(self):
        if self.operand("user_divisions"):
            return self.operand("user_divisions")
        else:
            npartitions = len(self.iterables[0])
            return (None,) * (npartitions + 1)

    @property
    def apply_func(self):
        if self.enforce_metadata:
            return apply_and_enforce
        return self.func

    @functools.cached_property
    def apply_kwargs(self):
        kwargs = self.kwargs
        if self.enforce_metadata:
            kwargs = kwargs.copy()
            kwargs.update(
                {
                    "_func": self.func,
                    "_meta": self._meta,
                }
            )
        return kwargs

    def _filtered_task(self, name: Key, index: int) -> Task:
        vals = [v[index] for v in self.iterables]
        if self.enforce_metadata:
            return Task(
                name,
                apply_and_enforce,
                *vals,
                *self.args,
                **self.apply_kwargs,
                _data_producer=True,
            )
        return Task(
            name, self.func, *vals, *self.args, **self.apply_kwargs, _data_producer=True
        )


class FromMapProjectable(FromMap):
    _parameters = [
        "func",
        "iterables",
        "columns",
        "args",
        "kwargs",
        "columns_arg_required",
        "user_meta",
        "enforce_metadata",
        "user_divisions",
        "label",
        "_partitions",
        "_series",
    ]
    _defaults = {
        "user_meta": no_default,
        "enforce_metadata": False,
        "user_divisions": None,
        "label": None,
        "_partitions": None,
        "_series": False,
    }
    _absorb_projections = True

    @functools.cached_property
    def columns_operand(self):
        return _convert_to_list(self.operand("columns"))

    @property
    def columns(self):
        if self.columns_operand is None:
            return list(self.frame_meta.columns)
        else:
            return self.columns_operand

    @functools.cached_property
    def _series(self):
        # Only need to convert to _series if func
        # doesn't produce a Series already
        return self.operand("_series") and self.frame_meta.ndim > 1

    @functools.cached_property
    def kwargs(self):
        options = self.operand("kwargs")
        if self.columns_arg_required or self.columns_operand:
            options = options.copy()
            options["columns"] = self.columns
        return options

    @functools.cached_property
    def apply_kwargs(self):
        kwargs = self.kwargs
        if self.enforce_metadata:
            kwargs = kwargs.copy()
            kwargs.update(
                {
                    "_func": self.func,
                    "_meta": self.frame_meta,
                }
            )
        return kwargs

    @functools.cached_property
    def frame_meta(self):
        # This is our `_meta` result before possibly
        # converting to a Series
        meta = super()._meta
        if meta.ndim > 1 and self.columns_operand is not None:
            return meta[self.columns_operand]
        return meta

    @property
    def _meta(self):
        # This is our final `_meta` result
        # (may need to be a Series)
        meta = self.frame_meta
        if self._series:
            assert len(self.columns_operand) > 0
            return meta[self.columns_operand[0]]
        return meta

    def _filtered_task(self, name: Key, index: int) -> Task:
        tsk = super()._filtered_task(name, index)
        if self._series:
            return Task(name, operator.getitem, tsk, self.columns[0])
        return tsk


class FromPandas(PartitionsFiltered, BlockwiseIO):
    """The only way today to get a real dataframe"""

    _parameters = [
        "frame",
        "npartitions",
        "sort",
        "chunksize",
        "columns",
        "pyarrow_strings_enabled",
        "_partitions",
        "_series",
        "_pd_length_stats",
    ]
    _defaults = {
        "npartitions": None,
        "sort": True,
        "columns": None,
        "_partitions": None,
        "_series": False,
        "chunksize": None,
        "pyarrow_strings_enabled": True,
        "_pd_length_stats": None,
    }
    _pd_length_stats: tuple | None
    _absorb_projections = True

    @functools.cached_property
    def frame(self):
        frame = self.operand("frame")._data
        if self.sort and not frame.index.is_monotonic_increasing:
            frame = frame.sort_index()
            return _BackendData(frame)
        return self.operand("frame")

    @functools.cached_property
    def _meta(self):
        if self.pyarrow_strings_enabled:
            meta = make_meta(to_pyarrow_string(self.frame.head(1)))
        else:
            meta = self.frame.head(0)

        if self.operand("columns") is not None:
            return meta[self.columns[0]] if self._series else meta[self.columns]
        return meta

    @functools.cached_property
    def columns(self):
        columns_operand = self.operand("columns")
        if columns_operand is None:
            try:
                return list(self.frame.columns)
            except AttributeError:
                if self.ndim == 1:
                    return [self.name]
                return []
        else:
            return _convert_to_list(columns_operand)

    @functools.cached_property
    def _divisions_and_locations(self):
        assert isinstance(self.frame, _BackendData)
        npartitions = self.operand("npartitions")
        sort = self.sort
        key = (npartitions, sort)
        _division_info_cache = self.frame._division_info
        if key not in _division_info_cache:
            data = self.frame._data
            nrows = len(data)
            if nrows == 0:
                npartitions = 1 if not npartitions else npartitions
                locations = [0] * (npartitions + 1)
                divisions = (None,) * len(locations)
            elif sort or self.frame._data.index.is_monotonic_increasing:
                divisions, locations = sorted_division_locations(
                    data.index,
                    npartitions=npartitions,
                    chunksize=self.operand("chunksize"),
                )
            else:
                if npartitions is None:
                    chunksize = self.operand("chunksize")
                else:
                    chunksize = int(math.ceil(nrows / npartitions))
                locations = list(range(0, nrows, chunksize)) + [len(data)]
                divisions = (None,) * len(locations)
            _division_info_cache[key] = divisions, locations
        return _division_info_cache[key]

    def _get_lengths(self) -> tuple | None:
        if self._pd_length_stats is None:
            locations = self._locations()
            self._pd_length_stats = tuple(
                offset - locations[i]
                for i, offset in enumerate(locations[1:])
                if not self._filtered or i in self._partitions
            )
        return self._pd_length_stats

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Lengths):
            _lengths = self._get_lengths()
            if _lengths:
                return Literal(_lengths)

        if isinstance(parent, Len):
            _lengths = self._get_lengths()
            if _lengths:
                return Literal(sum(_lengths))

        if isinstance(parent, Projection):
            return super()._simplify_up(parent, dependents)

    def _divisions(self):
        return self._divisions_and_locations[0]

    @functools.cached_property
    def npartitions(self):
        if self._filtered:
            return super().npartitions
        return len(self._divisions_and_locations[0]) - 1

    def _locations(self):
        return self._divisions_and_locations[1]

    def _filtered_task(self, name: Key, index: int) -> DataNode:  # type: ignore
        start, stop = self._locations()[index : index + 2]
        part = self.frame.iloc[start:stop]
        if self.pyarrow_strings_enabled:
            part = to_pyarrow_string(part)
        if self.operand("columns") is not None:
            return DataNode(
                name, part[self.columns[0]] if self._series else part[self.columns]
            )
        return DataNode(name, part)

    def __str__(self):
        if self._absorb_projections and self.operand("columns"):
            if self._series:
                return f"df[{self.columns[0]}]"
            return f"df[{self.columns}]"
        return "df"

    __repr__ = __str__


class FromPandasDivisions(FromPandas):
    _parameters = [
        "frame",
        "divisions",
        "columns",
        "pyarrow_strings_enabled",
        "_partitions",
        "_series",
        "_pd_length_stats",
    ]
    _defaults = {
        "columns": None,
        "_partitions": None,
        "_series": False,
        "_pd_length_stats": None,
    }
    sort = True

    @functools.cached_property
    def _name(self):
        return "from_pd_divs" + "-" + self.deterministic_token

    @property
    def _divisions_and_locations(self):
        assert isinstance(self.frame, _BackendData)
        key = tuple(self.operand("divisions"))
        _division_info_cache = self.frame._division_info
        if key not in _division_info_cache:
            data = self.frame._data
            if data.index.is_unique:
                indexer = data.index.get_indexer(key, method="bfill")
            else:
                # get_indexer for doesn't support method
                indexer = np.searchsorted(data.index.values, key, side="left")
            indexer[-1] = len(data)
            _division_info_cache[key] = key, indexer
        return _division_info_cache[key]


class FromScalars(IO):
    _parameters = ["meta", "names"]

    @property
    def _scalars(self):
        return self.dependencies()

    def _divisions(self):
        return (min(self.names), max(self.names))

    @functools.cached_property
    def _meta(self):
        return type(self.meta)(
            [s._meta for s in self._scalars], index=self.names, name=self.meta.name
        )

    def _layer(self) -> dict:
        return {
            (self._name, 0): (
                type(self.meta),
                [(s._name, 0) for s in self._scalars],
                self.names,
                None,
                self.meta.name,
            )
        }

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Projection):
            if sorted(parent.columns) == sorted(self.names):
                return
            new_names, new_scalars = [], []
            for n, s in zip(self.names, self._scalars):
                if n in parent.columns:
                    new_names.append(n)
                    new_scalars.append(s)
            return type(parent)(
                type(self)(self.meta, new_names, *new_scalars), *parent.operands[1:]
            )


class FromArray(PartitionsFiltered, BlockwiseIO):
    _parameters = [
        "frame",
        "chunksize",
        "original_columns",
        "meta",
        "columns",
        "_partitions",
    ]
    _defaults = {
        "chunksize": 50_000,
        "original_columns": None,
        "meta": None,
        "columns": None,
        "_partitions": None,
    }
    _pd_length_stats = None
    _absorb_projections = True

    @functools.cached_property
    def _meta(self):
        meta = _meta_from_array(
            self.frame, self.operand("original_columns"), self.operand("meta")
        )
        if self.operand("columns") is not None:
            return meta[self.operand("columns")]
        return meta

    @functools.cached_property
    def original_columns(self):
        if self.operand("original_columns") is None:
            if is_series_like(self._meta):
                return [0]
            return list(range(len(self._meta.columns)))
        return self.operand("original_columns")

    @functools.cached_property
    def _column_indices(self):
        if self.operand("columns") is None:
            return slice(0, len(self.original_columns))
        return [
            i
            for i, col in enumerate(self.original_columns)
            if col in self.operand("columns")
        ]

    def _divisions(self):
        divisions = tuple(range(0, len(self.frame), self.chunksize))
        divisions = divisions + (len(self.frame) - 1,)
        return divisions

    @functools.cached_property
    def unfiltered_divisions(self):
        return self._divisions()

    def _filtered_task(self, name: Key, index: int) -> Task:
        data = self.frame[slice(index * self.chunksize, (index + 1) * self.chunksize)]
        if index == len(self.unfiltered_divisions) - 2:
            idx = range(
                self.unfiltered_divisions[index],
                self.unfiltered_divisions[index + 1] + 1,
            )
        else:
            idx = range(
                self.unfiltered_divisions[index], self.unfiltered_divisions[index + 1]
            )

        if is_series_like(self._meta):
            return Task(
                name,
                type(self._meta),
                data,
                idx,
                self._meta.dtype,
                self._meta.name,
                _data_producer=True,
            )
        else:
            if data.ndim == 2:
                data = data[:, self._column_indices]
            return Task(
                name,
                type(self._meta),
                data,
                idx,
                self._meta.columns,
                _data_producer=True,
            )
