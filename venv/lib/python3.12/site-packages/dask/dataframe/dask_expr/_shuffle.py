from __future__ import annotations

import functools
import math
import operator
import uuid

import numpy as np
import pandas as pd
import tlz as toolz
from pandas import CategoricalDtype

from dask import compute
from dask._task_spec import Task, TaskRef
from dask.dataframe.core import _concat
from dask.dataframe.dask_expr._expr import (
    Assign,
    Blockwise,
    Expr,
    Filter,
    PartitionsFiltered,
    Projection,
    ToSeriesIndex,
    determine_column_projection,
    is_filter_pushdown_available,
)
from dask.dataframe.dask_expr._reductions import (
    All,
    Any,
    Count,
    DropDuplicates,
    Len,
    Max,
    Mean,
    MemoryUsage,
    Min,
    Mode,
    NBytes,
    NFirst,
    NLargest,
    NLast,
    NSmallest,
    Prod,
    Size,
    Sum,
    Unique,
    ValueCounts,
)
from dask.dataframe.dask_expr._repartition import Repartition, RepartitionToFewer
from dask.dataframe.dask_expr._util import LRU, _convert_to_list
from dask.dataframe.dispatch import is_categorical_dtype, make_meta
from dask.dataframe.shuffle import (
    barrier,
    collect,
    ensure_cleanup_on_exception,
    maybe_buffered_partd,
    partitioning_index,
    set_partitions_pre,
    shuffle_group,
    shuffle_group_2,
    shuffle_group_get,
)
from dask.utils import (
    M,
    digit,
    get_default_shuffle_method,
    insert,
    is_index_like,
    is_series_like,
)


class ShuffleBase(Expr):
    _parameters = [
        "frame",
        "partitioning_index",
        "npartitions_out",
        "ignore_index",
        "method",
        "options",
        "index_shuffle",
        "original_partitioning_index",
    ]
    _defaults = {
        "ignore_index": False,
        "method": None,
        "options": None,
        "index_shuffle": None,
        "original_partitioning_index": None,
    }
    _is_length_preserving = True
    _filter_passthrough = True

    def __str__(self):
        return f"Shuffle({self._name[-7:]})"

    def _node_label_args(self):
        return [self.frame, self.partitioning_index]

    @functools.cached_property
    def _partitioning_index(self):
        partitioning_index = self.partitioning_index
        if isinstance(partitioning_index, (str, int)):
            partitioning_index = [partitioning_index]
        return partitioning_index

    @functools.cached_property
    def unique_partition_mapping_columns_from_shuffle(self):
        idx = self.original_partitioning_index or self._partitioning_index
        return {tuple(idx)} if isinstance(idx, list) else set()

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Filter) and self._filter_passthrough_available(
            parent, dependents
        ):
            return self._filter_simplification(parent)
        if isinstance(parent, Projection):
            # Move the column projection to come
            # before the abstract Shuffle
            projection = _convert_to_list(
                determine_column_projection(self, parent, dependents)
            )
            partitioning_index = self._partitioning_index

            target = self.frame
            new_projection = [
                col
                for col in target.columns
                if (col in partitioning_index or col in projection)
            ]
            if set(new_projection) < set(target.columns):
                return type(self)(target[new_projection], *self.operands[1:])[
                    parent.operand("columns")
                ]

        if isinstance(
            parent,
            (
                Unique,
                DropDuplicates,
                Sum,
                Prod,
                Max,
                Any,
                All,
                Min,
                Len,
                Size,
                NBytes,
                Mean,
                Count,
                Mode,
                NLargest,
                NSmallest,
                ValueCounts,
                MemoryUsage,
            ),
        ):
            return type(parent)(self.frame, *parent.operands[1:])

    def _layer(self):
        raise NotImplementedError(
            f"{self} is abstract! Please call `simplify`"
            f"before generating a task graph."
        )

    @functools.cached_property
    def _meta(self):
        meta = self.frame._meta
        if self.ignore_index and self.method == "tasks":
            meta = meta.reset_index(drop=True)
        return meta

    def _divisions(self):
        return (None,) * (self.npartitions_out + 1)


class Shuffle(ShuffleBase):
    """Abstract shuffle class

    Parameters
    ----------
    frame: Expr
        The DataFrame-like expression to shuffle.
    partitioning_index: str, list
        Column and/or index names to hash and partition by.
    npartitions: int
        Number of output partitions.
    ignore_index: bool
        Whether to ignore the index during this shuffle operation.
    method: str or Callable
        Label or callback function to convert a shuffle operation
        to its necessary components.
    options: dict
        Algorithm-specific options.
    index_shuffle : bool
        Whether to perform the shuffle on the index.
    """

    def _lower(self):
        # Use `method` to decide how to compose a
        # shuffle operation from concerete expressions

        # Reduce partition count if necessary
        frame = self.frame
        npartitions_out = self.npartitions_out
        method = self.method or get_default_shuffle_method()

        if npartitions_out < frame.npartitions and method != "p2p":
            frame = Repartition(frame, new_partitions=npartitions_out)

        ops = [
            self.partitioning_index,
            self.npartitions_out,
            self.ignore_index,
            self.options,
            self.original_partitioning_index,
        ]
        if method == "p2p":
            return P2PShuffle(frame, *ops)
        elif method == "disk":
            return DiskShuffle(frame, *ops)
        elif method == "simple":
            return SimpleShuffle(frame, *ops)
        elif method == "tasks":
            return TaskShuffle(frame, *ops)
        else:
            raise ValueError(f"{method} not supported")


def _is_numeric_cast_type(dtype):
    return (
        pd.api.types.is_numeric_dtype(dtype)
        or isinstance(dtype, CategoricalDtype)
        and pd.api.types.is_numeric_dtype(dtype.categories)
    )


class RearrangeByColumn(ShuffleBase):
    def _lower(self):
        frame = self.frame
        partitioning_index = self.partitioning_index
        npartitions_out = self.npartitions_out
        ignore_index = self.ignore_index
        options = self.options
        index_shuffle = self.index_shuffle

        # Normalize partitioning_index

        if isinstance(partitioning_index, str):
            partitioning_index = [partitioning_index]
        if index_shuffle:
            pass
        elif not isinstance(partitioning_index, (list, Expr)):
            raise ValueError(
                f"{type(partitioning_index)} not a supported type for partitioning_index"
            )

        if not isinstance(partitioning_index, Expr) and not index_shuffle:
            cs = [col for col in partitioning_index if col not in frame.columns]
            if len(cs) == 1:
                frame = Assign(frame, "_partitions_0", frame.index)
                partitioning_index = partitioning_index.copy()
                partitioning_index[partitioning_index.index(cs[0])] = "_partitions_0"

        # Assign new "_partitions" column
        index_added = AssignPartitioningIndex(
            frame,
            partitioning_index,
            "_partitions",
            npartitions_out,
            frame._meta,
            index_shuffle,
        )

        # Apply shuffle
        shuffled = Shuffle(
            index_added,
            "_partitions",
            npartitions_out,
            ignore_index,
            self.method,
            options,
            original_partitioning_index=self._partitioning_index,
        )
        if frame.ndim == 1:
            # Reduce back to series
            return shuffled[index_added.columns[0]]

        # Drop "_partitions" column and return
        return shuffled[
            [c for c in shuffled.columns if c not in ["_partitions", "_partitions_0"]]
        ]


class SimpleShuffle(PartitionsFiltered, Shuffle):
    _parameters = [
        "frame",
        "partitioning_index",
        "npartitions_out",
        "ignore_index",
        "options",
        "original_partitioning_index",
        "_partitions",
    ]

    _defaults = {
        "_partitions": None,
        "original_partitioning_index": None,
        "partitioning_index": "_partitions",
        "ignore_index": False,
        "options": None,
    }

    @functools.cached_property
    def _meta(self):
        return self.frame._meta

    @staticmethod
    def _shuffle_group(df, _filter, *args):
        """Filter the output of `shuffle_group`"""
        if _filter is None:
            return shuffle_group(df, *args)
        return {k: v for k, v in shuffle_group(df, *args).items() if k in _filter}

    def _layer(self):
        """Construct graph for a simple shuffle operation."""
        shuffle_group_name = "group-" + self._name
        split_name = "split-" + self._name
        npartitions = self.npartitions_out

        dsk = {}
        _filter = self._partitions if self._filtered else None
        for global_part, part_out in enumerate(self._partitions):
            _concat_list = [
                (split_name, part_out, part_in)
                for part_in in range(self.frame.npartitions)
            ]
            dsk[(self._name, global_part)] = (
                _concat,
                _concat_list,
                self.ignore_index,
            )
            for _, _part_out, _part_in in _concat_list:
                dsk[(split_name, _part_out, _part_in)] = (
                    operator.getitem,
                    (shuffle_group_name, _part_in),
                    _part_out,
                )
                if (shuffle_group_name, _part_in) not in dsk:
                    dsk[(shuffle_group_name, _part_in)] = (
                        self._shuffle_group,
                        (self.frame._name, _part_in),
                        _filter,
                        self.partitioning_index,
                        0,
                        npartitions,
                        npartitions,
                        self.ignore_index,
                        npartitions,
                    )

        return dsk

    def _lower(self):
        return None


class TaskShuffle(SimpleShuffle):
    """Staged task-based shuffle implementation"""

    @functools.cached_property
    def _meta(self):
        meta = self.frame._meta
        if self.ignore_index:
            meta = meta.reset_index(drop=True)
        return meta

    def _layer(self):
        max_branch = (self.options or {}).get("max_branch", None) or 32
        npartitions_input = self.frame.npartitions
        if len(self._partitions) <= max_branch or npartitions_input <= max_branch:
            # We are creating a small number of output partitions,
            # or starting with a small number of input partitions.
            # No need for staged shuffling. Staged shuffling will
            # sometimes require extra work/communication in this case.
            return super()._layer()

        # Calculate number of stages and splits per stage
        npartitions = self.npartitions_out
        stages = int(math.ceil(math.log(npartitions_input) / math.log(max_branch)))
        if stages > 1:
            nsplits = int(math.ceil(npartitions_input ** (1 / stages)))
        else:
            nsplits = npartitions_input

        # Construct global data-movement plan
        inputs = [
            tuple(digit(i, j, nsplits) for j in range(stages))
            for i in range(nsplits**stages)
        ]
        inp_part_map = {inp: i for i, inp in enumerate(inputs)}
        parts_out = range(len(inputs))

        # Build graph
        dsk = {}
        name = self.frame._name
        meta_input = make_meta(self.frame._meta)
        for stage in range(stages):
            # Define names
            name_input = name
            if stage == (stages - 1) and npartitions == npartitions_input:
                name = self._name
                parts_out = self._partitions
                _filter = parts_out if self._filtered else None
            else:
                name = f"stage-{stage}-{self._name}"
                _filter = None

            shuffle_group_name = "group-" + name
            split_name = "split-" + name

            for global_part, part in enumerate(parts_out):
                out = inputs[part]

                _concat_list = []  # get_item tasks to concat for this output partition
                for i in range(nsplits):
                    # Get out each individual dataframe piece from the dicts
                    _inp = insert(out, stage, i)
                    _idx = out[stage]
                    _concat_list.append((split_name, _idx, _inp))

                # concatenate those pieces together, with their friends
                dsk[(name, global_part)] = (
                    _concat,
                    _concat_list,
                    self.ignore_index,
                )

                for _, _idx, _inp in _concat_list:
                    dsk[(split_name, _idx, _inp)] = (
                        operator.getitem,
                        (shuffle_group_name, _inp),
                        _idx,
                    )

                    if (shuffle_group_name, _inp) not in dsk:
                        # Initial partitions (output of previous stage)
                        _part = inp_part_map[_inp]
                        if stage == 0:
                            if _part < npartitions_input:
                                input_key = (name_input, _part)
                            else:
                                # In order to make sure that to_serialize() serialize the
                                # empty dataframe input, we add it as a key.
                                input_key = (shuffle_group_name, _inp, "empty")
                                dsk[input_key] = meta_input
                        else:
                            input_key = (name_input, _part)

                        # Convert partition into dict of dataframe pieces
                        dsk[(shuffle_group_name, _inp)] = (
                            self._shuffle_group,
                            input_key,
                            _filter,
                            self.partitioning_index,
                            stage,
                            nsplits,
                            npartitions_input,
                            self.ignore_index,
                            npartitions,
                        )

        if npartitions != npartitions_input:
            repartition_group_name = "repartition-group-" + name

            dsk2 = {
                (repartition_group_name, i): (
                    shuffle_group_2,
                    (name, i),
                    self.partitioning_index,
                    self.ignore_index,
                    npartitions,
                )
                for i in range(npartitions_input)
            }

            for i, p in enumerate(self._partitions):
                dsk2[(self._name, i)] = (
                    shuffle_group_get,
                    (repartition_group_name, p % npartitions_input),
                    p,
                )

            dsk.update(dsk2)
        return dsk


class DiskShuffle(SimpleShuffle):
    """Disk-based shuffle implementation"""

    @staticmethod
    def _shuffle_group(df, col, _filter, p):
        with ensure_cleanup_on_exception(p):
            g = df.groupby(col)
            d = {i: g.get_group(i) for i in g.groups if i in _filter}
            p.append(d, fsync=True)

    def _layer(self):
        from dask.dataframe.dispatch import partd_encode_dispatch

        column = self.partitioning_index
        df = self.frame

        always_new_token = uuid.uuid1().hex

        p = ("zpartd-" + always_new_token,)
        encode_cls = partd_encode_dispatch(df._meta)
        dsk1 = {p: (maybe_buffered_partd(encode_cls=encode_cls),)}

        # Partition data on disk
        name = "shuffle-partition-" + always_new_token
        dsk2 = {
            (name, i): (self._shuffle_group, key, column, self._partitions, p)
            for i, key in enumerate(df.__dask_keys__())
        }

        # Barrier
        barrier_token = ("barrier-" + always_new_token,)
        dsk3 = {barrier_token: (barrier, list(dsk2))}

        # Collect groups
        dsk4 = {
            (self._name, j): (collect, p, k, df._meta, barrier_token)
            for j, k in enumerate(self._partitions)
        }

        return toolz.merge(dsk1, dsk2, dsk3, dsk4)


def _shuffle_transfer(
    input: pd.DataFrame,
    id,
    input_partition: int,
) -> int:
    from distributed.shuffle._shuffle import shuffle_transfer

    return shuffle_transfer(
        input,
        id,
        input_partition,
    )


class P2PShuffle(SimpleShuffle):
    """P2P worker-based shuffle implementation"""

    @functools.cached_property
    def _meta(self):
        return self.frame._meta.drop(columns=self.partitioning_index)

    def _layer(self):
        from distributed.shuffle._core import (
            P2PBarrierTask,
            ShuffleId,
            barrier_key,
            p2p_barrier,
        )
        from distributed.shuffle._shuffle import DataFrameShuffleSpec, shuffle_unpack

        dsk = {}
        token = self._name.split("-")[-1]
        shuffle_id = ShuffleId(token)
        _barrier_key = barrier_key(shuffle_id)
        name = "shuffle-transfer-" + token

        parts_out = (
            self._partitions if self._filtered else list(range(self.npartitions_out))
        )
        # Avoid embedding a materialized list unless necessary
        parts_out_arg = (
            tuple(self._partitions) if self._filtered else self.npartitions_out
        )

        transfer_keys = list()
        for i in range(self.frame.npartitions):
            t = Task(
                (name, i),
                _shuffle_transfer,
                TaskRef((self.frame._name, i)),
                token,
                i,
            )
            dsk[t.key] = t
            transfer_keys.append(t.ref())

        barrier = P2PBarrierTask(
            _barrier_key,
            p2p_barrier,
            token,
            *transfer_keys,
            spec=DataFrameShuffleSpec(
                id=shuffle_id,
                npartitions=self.npartitions_out,
                column=self.partitioning_index,
                meta=self.frame._meta,
                parts_out=parts_out_arg,
                disk=True,
                drop_column=True,
            ),
        )
        dsk[barrier.key] = barrier

        # TODO: Decompose p2p Into transfer/barrier + unpack
        name = self._name
        for i, part_out in enumerate(parts_out):
            t = Task(
                (name, i),
                shuffle_unpack,
                token,
                part_out,
                barrier.ref(),
            )
            dsk[t.key] = t
        return dsk


#
# Helper logic
#


def _select_columns_or_index(df, columns_or_index):
    """
    Make a column selection that may include the index

    Parameters
    ----------
    columns_or_index
        Column or index name, or a list of these
    """

    # Ensure columns_or_index is a list
    columns_or_index = (
        columns_or_index if isinstance(columns_or_index, list) else [columns_or_index]
    )

    column_names = [n for n in columns_or_index if _is_column_label_reference(df, n)]

    selected_df = df[column_names]
    if _contains_index_name(df, columns_or_index):
        # Index name was included
        selected_df = selected_df.assign(_index=df.index)

    return selected_df


def _is_column_label_reference(df, key):
    """
    Test whether a key is a column label reference

    To be considered a column label reference, `key` must match the name of at
    least one column.
    """
    return (
        not isinstance(key, Expr)
        and (np.isscalar(key) or pd.api.types.is_scalar(key) or isinstance(key, tuple))
        and key in df.columns
    )


def _contains_index_name(df, columns_or_index):
    """
    Test whether the input contains a reference to the index of the df
    """
    if isinstance(columns_or_index, list):
        return any(_is_index_level_reference(df, n) for n in columns_or_index)
    else:
        return _is_index_level_reference(df, columns_or_index)


def _is_index_level_reference(df, key):
    """
    Test whether a key is an index level reference

    To be considered an index level reference, `key` must match the index name
    and must NOT match the name of any column.
    """
    index_name = df.index._meta.name if isinstance(df, Expr) else df.index.name
    return (
        index_name is not None
        and not isinstance(key, Expr)
        and (np.isscalar(key) or pd.api.types.is_scalar(key) or isinstance(key, tuple))
        and key == index_name
        and key not in getattr(df, "columns", ())
    )


class AssignPartitioningIndex(Blockwise):
    """Assign a partitioning index

    This class is used to construct a hash-based
    partitioning index for shuffling.

    Parameters
    ----------
    frame: Expr
        Frame-like expression being partitioned.
    partitioning_index: Expr or list
        Index-like expression or list of columns to construct
        the partitioning-index from.
    index_name: str
        New column name to assign.
    npartitions_out: int
        Number of partitions after repartitioning is finished.
    index_shuffle : bool, default False
        Whether we are using solely the index for the shuffle
    """

    _parameters = [
        "frame",
        "partitioning_index",
        "index_name",
        "npartitions_out",
        "meta",
        "index_shuffle",
    ]

    _defaults = {
        "cast_dtype": None,
        "index_shuffle": False,
        "partitioning_index": "_partitions",
        "index_name": "_partitions",
    }
    _preserves_partitioning_information = True

    @staticmethod
    def operation(df, index, name: str, npartitions: int, meta, index_shuffle: bool):  # type: ignore
        """Construct a hash-based partitioning index"""

        def _get_index(idx, obj):
            if hasattr(idx, "ndim"):
                if idx.ndim == 1:
                    idx = idx.to_frame()
            elif index_shuffle:
                # set default index name, otherwise we will end up with 0
                name = {"name": "_index"} if idx == ["_index"] else {}
                idx = obj.index.to_frame(**name)
            else:
                idx = _select_columns_or_index(obj, idx)
            return idx

        # meta and df dtypes can deviate, this is why we do the cast here
        meta_index = _get_index(index, meta)
        index = _get_index(index, df)

        dtypes = {}
        for col, dtype in meta_index.dtypes.items():
            if _is_numeric_cast_type(dtype):
                dtypes[col] = np.float64
        if dtypes:
            index = index.astype(dtypes, errors="ignore")

        index = partitioning_index(index, npartitions)
        if df.ndim == 1:
            df = df.to_frame()
        return df.assign(**{name: index})


class BaseSetIndexSortValues(Expr):
    _is_length_preserving = True

    def _divisions(self):
        if "user_divisions" in self._parameters and self.user_divisions is not None:
            return self.user_divisions
        if self._npartitions_input == 1:
            return (None, None)

        if (
            is_index_like(self._divisions_column._meta)
            and self._divisions_column.known_divisions
            and self._divisions_column.npartitions == self.frame.npartitions
        ):
            return self.other.divisions

        divisions, mins, maxes, presorted = _get_divisions(
            self.frame,
            self._divisions_column,
            self._npartitions_input,
            self.ascending,
            upsample=self.upsample,
        )
        if presorted and len(mins) == self._npartitions_input:
            divisions = mins.copy() + [maxes[-1]]
        return divisions

    @property
    def _npartitions_input(self):
        return self.operand("npartitions") or self.frame.npartitions

    @property
    def npartitions(self):
        return self.operand("npartitions") or len(self._divisions()) - 1


class SetIndex(BaseSetIndexSortValues):
    """Abstract ``set_index`` class.

    Simplifies (later lowers) either to Blockwise ops if we are already sorted
    or to ``SetPartition`` which handles shuffling.

    Parameters
    ----------
    frame: Expr
        Frame-like expression where the index is set.
    _other: Expr | Scalar
        Either a Series-like expression to use as Index or a scalar defining the column.
    drop: bool
        Whether we drop the old column.
    sorted: str
        No need for shuffling if we are already sorted.
    user_divisions: int
        Divisions as passed by the user.
    upsample: float
        Used to increase the number of samples for quantiles.
    """

    _parameters = [
        "frame",
        "_other",
        "drop",
        "user_divisions",
        "partition_size",
        "ascending",
        "npartitions",
        "upsample",
        "shuffle_method",
        "append",
        "options",  # Options for the chosen shuffle method
    ]
    _defaults = {
        "drop": True,
        "user_divisions": None,
        "partition_size": 128e6,
        "ascending": True,
        "npartitions": None,
        "upsample": 1.0,
        "shuffle_method": None,
        "options": None,
        "append": False,
    }
    _filter_passthrough = True

    @property
    def _projection_columns(self):
        return self.columns + (
            [self._other] if not isinstance(self._other, Expr) else []
        )

    @functools.cached_property
    def _meta(self):
        if isinstance(self._other, Expr):
            other = self._other._meta
        else:
            other = self._other
        return self.frame._meta.set_index(other, drop=self.drop)

    @property
    def _divisions_column(self):
        return self.other

    @property
    def other(self):
        if isinstance(self._other, Expr):
            return self._other
        return self.frame[self._other]

    def _lower(self):
        if (
            self.operand("npartitions") == 1
            or self.frame.npartitions == 1
            and (self.user_divisions is None or len(self.user_divisions) == 2)
            and self.operand("npartitions") is None
        ):
            expr = self.frame
            if self.frame.npartitions > 1:
                expr = RepartitionToFewer(expr, 1)

            index_set = SetIndexBlockwise(expr, self._other, self.drop, None)
            return SortIndexBlockwise(index_set)

        if self.user_divisions is None:
            divisions = self._divisions()
            if (
                is_index_like(self._divisions_column._meta)
                and self.other.divisions == divisions
            ):
                presorted = True
            else:
                presorted = _get_divisions(
                    self.frame,
                    self.other,
                    self._npartitions_input,
                    self.ascending,
                    upsample=self.upsample,
                )[3]

            if presorted and self.npartitions == self.frame.npartitions:
                index_set = SetIndexBlockwise(
                    self.frame, self._other, self.drop, divisions, self.append
                )
                return SortIndexBlockwise(index_set)

        return SetPartition(
            self.frame,
            self._other,
            self.drop,
            self._npartitions_input if self.user_divisions is None else None,
            self.ascending,
            self.upsample,
            self.user_divisions,
            self.shuffle_method,
            self.options,
        )

    def _simplify_up(self, parent, dependents):
        from dask.dataframe.dask_expr._expr import Filter, Head, Tail

        # TODO, handle setting index with other frame
        if (
            isinstance(parent, Head)
            and isinstance(self._other, (int, str))
            and self._other in self.frame.columns
        ):
            head = NFirst(self.frame, n=parent.n, _columns=self._other, ascending=True)
            return SetIndex(head, _other=self._other)

        if (
            isinstance(parent, Tail)
            and isinstance(self._other, (int, str))
            and self._other in self.frame.columns
        ):
            tail = NLast(self.frame, n=parent.n, _columns=self._other, ascending=True)
            return SetIndex(tail, _other=self._other)

        if isinstance(parent, Projection):
            addition_columns = (
                [self._other] if not isinstance(self._other, Expr) else []
            )
            columns = determine_column_projection(
                self, parent, dependents, additional_columns=addition_columns
            )
            columns = _convert_to_list(columns)
            columns = [c for c in self.frame.columns if c in columns]
            if self.frame.columns == columns:
                return
            return type(parent)(
                type(self)(self.frame[columns], *self.operands[1:]),
                parent.operand("columns"),
            )
        if isinstance(parent, Filter) and self._filter_passthrough_available(
            parent, dependents
        ):
            return self._filter_simplification(parent)

    def _filter_passthrough_available(self, parent, dependents):
        if is_filter_pushdown_available(self, parent, dependents):
            from dask.dataframe.dask_expr._expr import Index

            p = parent.predicate
            return not any(isinstance(x, Index) for x in p.walk())
        return False


class SortValues(BaseSetIndexSortValues):
    _parameters = [
        "frame",
        "by",
        "ascending",
        "na_position",
        "npartitions",
        "partition_size",
        "sort_function",
        "sort_function_kwargs",
        "upsample",
        "ignore_index",
        "shuffle_method",
        "options",  # Options for the chosen shuffle method
    ]
    _defaults = {
        "partition_size": 128e6,
        "ascending": True,
        "npartitions": None,
        "na_position": "last",
        "sort_function": None,
        "sort_function_kwargs": None,
        "upsample": 1.0,
        "ignore_index": False,
        "shuffle_method": None,
    }
    _filter_passthrough = True

    def _divisions(self):
        if self.frame.npartitions == 1:
            # Protect against triggering calculations when we only have one division
            return (None, None)

        divisions, mins, maxes, presorted = _get_divisions(
            self.frame,
            self.frame[self.by[0]],
            self._npartitions_input,
            self._divisions_ascending,
            upsample=self.upsample,
        )
        if presorted:
            return self.frame.divisions
        return (None,) * len(divisions)

    @property
    def _divisions_ascending(self) -> bool:
        divisions_ascending = self.ascending
        if not isinstance(divisions_ascending, bool):
            divisions_ascending = divisions_ascending[0]
        assert isinstance(divisions_ascending, bool)
        return divisions_ascending

    @property
    def sort_function(self):
        if self.operand("sort_function") is not None:
            return self.operand("sort_function")
        return M.sort_values

    @property
    def sort_function_kwargs(self):
        sort_kwargs = {
            "by": self.by,
            "ascending": self.ascending,
            "na_position": self.na_position,
            "ignore_index": self.ignore_index,
        }
        if self.operand("sort_function_kwargs") is not None:
            sort_kwargs.update(self.operand("sort_function_kwargs"))
        return sort_kwargs

    @functools.cached_property
    def _meta(self):
        return self.frame._meta

    @functools.cached_property
    def _meta_by_dtype(self):
        dtype = self._meta.dtypes[self.by]
        if is_series_like(dtype):
            dtype = dtype.iloc[0]
        return dtype

    def _lower(self):
        if self.frame.npartitions == 1:
            return SortValuesBlockwise(
                self.frame, self.sort_function, self.sort_function_kwargs
            )

        _divisions_by = self.frame[self.by[0]]
        divisions, _, _, presorted = _get_divisions(
            self.frame,
            _divisions_by,
            self._npartitions_input,
            self._divisions_ascending,
            upsample=self.upsample,
        )
        if presorted and self.npartitions == self.frame.npartitions:
            return SortValuesBlockwise(
                self.frame, self.sort_function, self.sort_function_kwargs
            )

        partitions = _SetPartitionsPreSetIndex(
            _divisions_by,
            _divisions_by._meta._constructor(divisions).sort_values(),
            ascending=self._divisions_ascending,
        )
        assigned = Assign(self.frame, "_partitions", partitions)
        shuffled = Shuffle(
            assigned,
            "_partitions",
            npartitions_out=len(divisions) - 1,
            ignore_index=self.ignore_index,
            method=self.shuffle_method,
            options=self.options,
        )
        shuffled = Projection(shuffled, self.frame.columns)
        return SortValuesBlockwise(
            shuffled, self.sort_function, self.sort_function_kwargs
        )

    def _simplify_up(self, parent, dependents):
        from dask.dataframe.dask_expr._expr import Filter, Head, Tail

        if isinstance(parent, Head):
            return NFirst(
                self.frame, n=parent.n, _columns=self.by, ascending=self.ascending
            )

        if isinstance(parent, Tail):
            return NLast(
                self.frame, n=parent.n, _columns=self.by, ascending=self.ascending
            )

        if isinstance(parent, Filter) and self._filter_passthrough_available(
            parent, dependents
        ):
            return self._filter_simplification(parent)
        if isinstance(parent, Projection):
            columns = determine_column_projection(
                self, parent, dependents, additional_columns=self.by
            )
            columns = _convert_to_list(columns)
            columns = [col for col in self.frame.columns if col in columns]
            if self.frame.columns == columns:
                return
            return type(parent)(
                type(self)(self.frame[columns], *self.operands[1:]),
                parent.operand("columns"),
            )
        if (
            isinstance(parent, Repartition)
            and parent.operand("new_partitions") is not None
        ):
            return type(self)(
                type(parent)(self.frame, *parent.operands[1:]), *self.operands[1:]
            )


class SetPartition(SetIndex):
    """Shuffles the DataFrame according to its new divisions.

    Simplifies the Expression to blockwise pre-processing, shuffle and
    blockwise post-processing expressions.

    Parameters
    ----------
    frame: Expr
        Frame-like expression where the index is set.
    _other: Expr | Scalar
        Either a Series-like expression to use as Index or a scalar defining the column.
    drop: bool
        Whether to drop the old column.
    new_divisions: int
        Divisions of the resulting expression.
    """

    _parameters = [
        "frame",
        "_other",
        "drop",
        "npartitions",
        "ascending",
        "upsample",
        "user_divisions",
        "shuffle_method",
        "options",  # Shuffle method options
    ]

    def _lower(self):
        divisions = self.other._meta._constructor(self._divisions())
        partitions = _SetPartitionsPreSetIndex(self.other, divisions)
        assigned = Assign(self.frame, "_partitions", partitions)
        if isinstance(self._other, Expr):
            assigned = Assign(assigned, "_index", self._other)
        shuffled = Shuffle(
            assigned,
            "_partitions",
            npartitions_out=len(self._divisions()) - 1,
            ignore_index=True,
            method=self.shuffle_method,
            options=self.options,
        )
        shuffled = Projection(
            shuffled, [c for c in assigned.columns if c != "_partitions"]
        )

        if isinstance(self._other, Expr):
            drop, set_name = True, "_index"
        else:
            drop, set_name = self.drop, self.other._meta.name
        lru_key = (
            self.other._name,
            self._npartitions_input,
            self.ascending,
            128e6,
            self.upsample,
        )
        computed_divisions = divisions_lru.get(lru_key)
        index_set = _SetIndexPost(
            shuffled,
            self.other._meta.name,
            drop,
            set_name,
            self.frame._meta.columns.dtype,
            computed_divisions,
            self.user_divisions,
        )
        return SortIndexBlockwise(index_set)


class _SetPartitionsPreSetIndex(Blockwise):
    _parameters = ["frame", "new_divisions", "ascending", "na_position"]
    _defaults = {"ascending": True, "na_position": "last"}
    operation = staticmethod(set_partitions_pre)
    _is_length_preserving = True

    @functools.cached_property
    def _meta(self):
        return make_meta(self.frame._meta._constructor([0]))


class _SetIndexPost(Blockwise):
    _parameters = [
        "frame",
        "index_name",
        "drop",
        "set_name",
        "column_dtype",
        "computed_divisions",
        "user_divisions",
    ]
    _is_length_preserving = True

    @property
    def _args(self) -> list:
        return self.operands[:5]

    @staticmethod
    def operation(df, index_name, drop, set_name, column_dtype):
        df = df.set_index(set_name, drop=drop)
        df.index.name = index_name
        df.columns = df.columns.astype(column_dtype)
        return df

    def _get_culled_divisions(self, divisions):
        if self.frame.npartitions < len(divisions) - 1:
            part_filter = list(self.frame.find_operations(PartitionsFiltered))
            if len(part_filter) > 0:
                return tuple(
                    [
                        div
                        for i, div in enumerate(divisions)
                        if i in part_filter[0]._partitions
                    ]
                    + [divisions[-1]]
                )
            else:
                return self.frame.divisions

        return divisions

    def _divisions(self):
        if self.operand("user_divisions") is not None:
            return self._get_culled_divisions(self.operand("user_divisions"))
        assert self.computed_divisions is not None
        return self._get_culled_divisions(self.computed_divisions[0])


class SortIndexBlockwise(Blockwise):
    _projection_passthrough = True
    _parameters = ["frame"]
    operation = M.sort_index
    _is_length_preserving = True


class SortValuesBlockwise(Blockwise):
    _projection_passthrough = False
    _parameters = ["frame", "sort_function", "sort_kwargs"]
    _keyword_only = ["sort_function", "sort_kwargs"]
    _is_length_preserving = True

    @staticmethod
    def operation(*args, **kwargs):
        sort_func = kwargs.pop("sort_function")
        sort_kwargs = kwargs.pop("sort_kwargs")
        return sort_func(*args, **kwargs, **sort_kwargs)

    @functools.cached_property
    def _meta(self):
        return self.frame._meta


class SetIndexBlockwise(Blockwise):
    _parameters = ["frame", "other", "drop", "new_divisions", "append"]
    _defaults = {"append": False, "new_divisions": None, "drop": True}
    _keyword_only = ["drop", "new_divisions", "append"]
    _is_length_preserving = True
    _preserves_partitioning_information = True

    @staticmethod
    def operation(df, *args, new_divisions, **kwargs):
        return df.set_index(*args, **kwargs)

    def _divisions(self):
        if self.new_divisions is None:
            return (None,) * (self.frame.npartitions + 1)
        return tuple(self.new_divisions)

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Projection):
            columns = determine_column_projection(
                self,
                parent,
                dependents,
                additional_columns=_convert_to_list(self.other),
            )
            if self.frame.columns == columns:
                return
            columns = [col for col in self.frame.columns if col in columns]
            return type(parent)(
                type(self)(self.frame[columns], *self.operands[1:]),
                parent.operand("columns"),
            )


divisions_lru = LRU(10)  # type: ignore


def _get_divisions(
    frame,
    other,
    npartitions: int,
    ascending: bool = True,
    partition_size: float = 128e6,
    upsample: float = 1.0,
):
    key = (other._name, npartitions, ascending, partition_size, upsample)
    if key in divisions_lru:
        return divisions_lru[key]
    result = _calculate_divisions(
        frame, other, npartitions, ascending, partition_size, upsample
    )
    divisions_lru[key] = result
    return result


def _calculate_divisions(
    frame,
    other,
    npartitions: int,
    ascending: bool = True,
    partition_size: float = 128e6,
    upsample: float = 1.0,
):
    from dask.dataframe.dask_expr import RepartitionQuantiles, new_collection

    if is_index_like(other._meta):
        other = ToSeriesIndex(other)

    if is_categorical_dtype(other._meta.dtype):
        other = new_collection(other).cat.as_ordered()._expr

    try:
        divisions, mins, maxes = compute(
            new_collection(RepartitionQuantiles(other, npartitions, upsample=upsample)),
            new_collection(other).map_partitions(M.min),
            new_collection(other).map_partitions(M.max),
        )
    except TypeError as e:
        # When there are nulls and a column is non-numeric, a TypeError is sometimes raised as a result of
        # 1) computing mins/maxes above, 2) every null being switched to NaN, and 3) NaN being a float.
        # Also, Pandas ExtensionDtypes may cause TypeErrors when dealing with special nulls such as pd.NaT or pd.NA.
        # If this happens, we hint the user about eliminating nulls beforehand.
        if not pd.api.types.is_numeric_dtype(other._meta.dtype):
            obj, suggested_method = (
                ("column", f"`.dropna(subset=['{other.name}'])`")
                if any(other._name == frame[c]._name for c in frame.columns)
                else ("series", "`.loc[series[~series.isna()]]`")
            )
            raise NotImplementedError(
                f"Divisions calculation failed for non-numeric {obj} '{other.name}'.\n"
                f"This is probably due to the presence of nulls, which Dask does not entirely support in the index.\n"
                f"We suggest you try with {suggested_method}."
            ) from e
        # For numeric types there shouldn't be problems with nulls, so we raise as-it-is this particular TypeError
        else:
            raise e

    sizes = []  # type: ignore

    empty_dataframe_detected = pd.isna(divisions).all()
    if empty_dataframe_detected:
        total = sum(sizes)
        npartitions = max(math.ceil(total / partition_size), 1)
        npartitions = min(npartitions, frame.npartitions)
        n = divisions.size
        try:
            divisions = np.interp(
                x=np.linspace(0, n - 1, npartitions + 1),
                xp=np.linspace(0, n - 1, n),
                fp=divisions.tolist(),
            ).tolist()
        except (TypeError, ValueError):  # str type
            indexes = np.linspace(0, n - 1, npartitions + 1).astype(int)
            divisions = divisions.iloc[indexes].tolist()
    else:
        # Drop duplicate divisions returned by partition quantiles
        n = divisions.size
        divisions = (
            list(divisions.iloc[: n - 1].unique()) + divisions.iloc[n - 1 :].tolist()
        )

    mins = mins.bfill()
    maxes = maxes.bfill()
    if isinstance(other._meta.dtype, pd.CategoricalDtype):
        dtype = other._meta.dtype
        mins = mins.astype(dtype)
        maxes = maxes.astype(dtype)

    if mins.isna().any() or maxes.isna().any():
        presorted = False
    else:
        n = mins.size
        maxes2 = (maxes.iloc[: n - 1] if ascending else maxes.iloc[1:]).reset_index(
            drop=True
        )
        mins2 = (mins.iloc[1:] if ascending else mins.iloc[: n - 1]).reset_index(
            drop=True
        )
        presorted = (
            mins.tolist() == mins.sort_values(ascending=ascending).tolist()
            and maxes.tolist() == maxes.sort_values(ascending=ascending).tolist()
            and (maxes2 < mins2).all()
        )
    return divisions, mins.tolist(), maxes.tolist(), presorted
