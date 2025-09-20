from __future__ import annotations

import functools
import warnings

import pandas as pd
from toolz import merge_sorted, unique

from dask.core import flatten
from dask.dataframe import methods
from dask.dataframe.dask_expr._expr import (
    AsType,
    Blockwise,
    Expr,
    Projection,
    ToFrame,
    are_co_aligned,
    determine_column_projection,
)
from dask.dataframe.dask_expr._util import _convert_to_list
from dask.dataframe.dispatch import make_meta, meta_nonempty
from dask.dataframe.multi import concat_and_check
from dask.dataframe.utils import check_meta, strip_unknown_categories
from dask.utils import apply, is_dataframe_like, is_series_like


class Concat(Expr):
    _parameters = [
        "join",
        "ignore_order",
        "axis",
        "ignore_unknown_divisions",
        "interleave_partitions",
        "_kwargs",
    ]
    _defaults = {
        "join": "outer",
        "ignore_order": False,
        "_kwargs": {},
        "axis": 0,
        "ignore_unknown_divisions": False,
        "interleave_partitions": False,
    }

    def __str__(self):
        s = (
            "frames="
            + str(self.dependencies())
            + ", "
            + ", ".join(self._operands_for_repr())
        )
        return f"{type(self).__name__}({s})"

    @property
    def _frames(self):
        return self.dependencies()

    @functools.cached_property
    def _meta(self):
        # ignore DataFrame without columns to avoid dtype upcasting
        return make_meta(
            methods.concat(
                [
                    meta_nonempty(df._meta)
                    for df in self._frames
                    if df.ndim < 2 or len(df._meta.columns) > 0
                ],
                join=self.join,
                filter_warning=False,
                axis=self.axis,
                ignore_order=self.ignore_order,
                **self._kwargs,
            )
        )

    def _divisions(self):
        dfs = self._frames

        if self.axis == 1:
            if self._are_co_alinged_or_single_partition:
                if {df.npartitions for df in self._frames} == {1}:
                    divisions = set(
                        flatten([e.divisions for e in dfs], container=tuple)
                    )
                    return min(divisions), max(divisions)
                return dfs[0].divisions
            elif self._all_known_divisions:
                divisions = list(unique(merge_sorted(*[df.divisions for df in dfs])))
                if len(divisions) == 1:  # single value for index
                    divisions = (divisions[0], divisions[0])
                return divisions
            return (None,) * (max(df.npartitions for df in dfs) + 1)

        if self._monotonic_divisions:
            divisions = []
            for df in dfs[:-1]:
                # remove last to concatenate with next
                divisions += df.divisions[:-1]
            divisions += dfs[-1].divisions
            return divisions
        elif self._all_known_divisions and self.interleave_partitions:
            divisions = list(unique(merge_sorted(*[df.divisions for df in dfs])))
            if len(divisions) == 1:  # single value for index
                divisions = (divisions[0], divisions[0])
            return divisions
        return [None] * (sum(df.npartitions for df in dfs) + 1)

    @functools.cached_property
    def interleave_partitions(self):
        if "interleave_partitions" in self._parameters:
            return self.operand("interleave_partitions")
        return False

    @functools.cached_property
    def _all_known_divisions(self):
        dfs = self._frames
        return all(df.known_divisions for df in dfs)

    @functools.cached_property
    def _monotonic_divisions(self):
        dfs = self._frames

        if self._all_known_divisions:
            # each DataFrame's division must be greater than previous one
            return all(
                dfs[i].divisions[-1] < dfs[i + 1].divisions[0]
                for i in range(len(dfs) - 1)
            )
        return False

    @functools.cached_property
    def _are_co_alinged_or_single_partition(self):
        return are_co_aligned(*self._frames) or {
            df.npartitions for df in self._frames
        } == {1}

    def _lower(self):
        dfs = self._frames
        if self.axis == 1:
            if self._are_co_alinged_or_single_partition:
                return ConcatIndexed(
                    self.ignore_order, self._kwargs, self.axis, self.join, *dfs
                )

            elif (
                all(not df.known_divisions for df in dfs)
                and len({df.npartitions for df in dfs}) == 1
            ):
                if not self.ignore_unknown_divisions:
                    warnings.warn(
                        "Concatenating dataframes with unknown divisions.\n"
                        "We're assuming that the indices of each dataframes"
                        " are \n aligned. This assumption is not generally "
                        "safe."
                    )
                return ConcatUnindexed(
                    self.ignore_order, self._kwargs, self.axis, self.join, *dfs
                )
            elif self._all_known_divisions:
                from dask.dataframe.dask_expr._repartition import Repartition

                divs = self._divisions()
                cast_dfs = [
                    Repartition(df, new_divisions=divs, force=True) for df in dfs
                ]
                return StackPartitionInterleaved(
                    self.join,
                    self.ignore_order,
                    self.axis,
                    self.ignore_unknown_divisions,
                    self.interleave_partitions,
                    self._kwargs,
                    *cast_dfs,
                )

            else:
                raise ValueError(
                    "Unable to concatenate DataFrame with unknown "
                    "division specifying axis=1"
                )

        cast_dfs = []
        for df in dfs:
            # dtypes of all dfs need to be coherent
            # refer to https://github.com/dask/dask/issues/4685
            # and https://github.com/dask/dask/issues/5968.
            if is_dataframe_like(df._meta):
                shared_columns = list(set(df.columns).intersection(self._meta.columns))
                needs_astype = {
                    col: self._meta[col].dtype
                    for col in shared_columns
                    if df._meta[col].dtype != self._meta[col].dtype
                    and not isinstance(df[col]._meta.dtype, pd.CategoricalDtype)
                }

                if needs_astype:
                    cast_dfs.append(AsType(df, dtypes=needs_astype))
                else:
                    cast_dfs.append(df)
            elif is_series_like(df) and is_series_like(self._meta):
                if df.dtype != self._meta.dtype and not isinstance(
                    df.dtype, pd.CategoricalDtype
                ):
                    cast_dfs.append(AsType(df, dtypes=self._meta.dtype))
                else:
                    cast_dfs.append(df)
            else:
                cast_dfs.append(df)

        if (
            not self._monotonic_divisions
            and self._all_known_divisions
            and self.interleave_partitions
        ):
            from dask.dataframe.dask_expr._repartition import Repartition

            divs = self._divisions()
            cast_dfs = [
                Repartition(df, new_divisions=divs, force=True) for df in cast_dfs
            ]
            return StackPartitionInterleaved(
                self.join,
                self.ignore_order,
                self.axis,
                self.ignore_unknown_divisions,
                self.interleave_partitions,
                self._kwargs,
                *cast_dfs,
            )

        return StackPartition(
            self.join,
            self.ignore_order,
            self.axis,
            self.ignore_unknown_divisions,
            self.interleave_partitions,
            self._kwargs,
            *cast_dfs,
        )

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Projection):

            def get_columns_or_name(e: Expr):
                return e.columns if e.ndim == 2 else [e.name]

            columns = determine_column_projection(self, parent, dependents)
            columns = _convert_to_list(columns)
            columns_frame = [
                [col for col in get_columns_or_name(frame) if col in columns]
                for frame in self._frames
            ]
            if all(
                set(cols) == set(get_columns_or_name(frame))
                and len(cols) == len(get_columns_or_name(frame))
                for frame, cols in zip(self._frames, columns_frame)
            ):
                return

            frames = [
                (
                    frame[cols]
                    if set(cols) != set(get_columns_or_name(frame))
                    or len(cols) != len(get_columns_or_name(frame))
                    else frame
                )
                for frame, cols in zip(self._frames, columns_frame)
                if len(cols) > 0
            ]
            result = type(self)(
                self.join,
                self.ignore_order,
                self.axis,
                self.ignore_unknown_divisions,
                self.interleave_partitions,
                self._kwargs,
                *frames,
            )

            if result.columns == _convert_to_list(parent.operand("columns")):
                if result.ndim == parent.ndim:
                    return result
                elif result.ndim < parent.ndim:
                    return ToFrame(result)

            return type(parent)(result, *parent.operands[1:])


class StackPartition(Concat):

    def _layer(self):
        dsk, i = {}, 0
        kwargs = self._kwargs.copy()
        kwargs["ignore_order"] = self.ignore_order
        ctr = 0
        meta = strip_unknown_categories(self._meta)
        for df in self._frames:
            try:
                check_meta(df._meta, self._meta)
                match = True
            except (ValueError, TypeError):
                match = False

            for i in range(df.npartitions):
                if match:
                    dsk[(self._name, ctr)] = df._name, i
                else:
                    dsk[(self._name, ctr)] = (
                        apply,
                        methods.concat,
                        [
                            [meta, (df._name, i)],
                            self.axis,
                            self.join,
                            False,
                            True,
                        ],
                        kwargs,
                    )
                ctr += 1
        return dsk

    def _lower(self):
        return


class StackPartitionInterleaved(StackPartition):
    def _divisions(self):
        return self._frames[0].divisions

    def _layer(self):
        dsk, i = {}, 0
        kwargs = self._kwargs.copy()
        kwargs["ignore_order"] = self.ignore_order

        dfs = self._frames
        for i in range(self.npartitions):
            dsk[(self._name, i)] = (
                apply,
                methods.concat,
                [
                    [(df._name, i) for df in dfs],
                    self.axis,
                    self.join,
                    False,
                    True,
                ],
                kwargs,
            )
        return dsk


class ConcatUnindexed(Blockwise):
    _parameters = ["ignore_order", "_kwargs", "axis", "join"]
    _defaults = {"ignore_order": False, "_kwargs": {}, "axis": 1, "join": "outer"}
    _keyword_only = ["ignore_order", "_kwargs", "axis", "join"]

    @functools.cached_property
    def _meta(self):
        return methods.concat(
            [df._meta for df in self.dependencies()],
            ignore_order=self.ignore_order,
            axis=self.axis,
            join=self.join,
            **self.operand("_kwargs"),
        )

    @staticmethod
    def operation(*args, ignore_order, _kwargs, axis, join):
        return concat_and_check(args, ignore_order=ignore_order)


class ConcatIndexed(ConcatUnindexed):
    @staticmethod
    def operation(*args, ignore_order, _kwargs, axis, join):
        return methods.concat(args, ignore_order=ignore_order, axis=axis, join=join)

    def _broadcast_dep(self, dep: Expr):
        return dep.npartitions == 1
