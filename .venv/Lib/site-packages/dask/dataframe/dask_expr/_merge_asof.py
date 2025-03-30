from __future__ import annotations

import functools

import pandas as pd

from dask.dataframe import methods
from dask.dataframe.dask_expr import SetIndexBlockwise, new_collection
from dask.dataframe.dask_expr._expr import MapPartitions, RenameAxis, ResetIndex
from dask.dataframe.dask_expr._merge import Merge
from dask.dataframe.dask_expr._util import _BackendData
from dask.dataframe.dask_expr.io import FromPandas
from dask.dataframe.dispatch import make_meta, meta_nonempty
from dask.dataframe.multi import merge_asof_padded, pair_partitions
from dask.dataframe.utils import pyarrow_strings_enabled
from dask.utils import apply


class MergeAsof(Merge):
    _parameters = [
        "left",
        "right",
        "left_on",
        "right_on",
        "left_index",
        "right_index",
        "left_by",
        "right_by",
        "suffixes",
        "tolerance",
        "allow_exact_matches",
        "direction",
    ]
    _defaults = {
        "left_on": None,
        "right_on": None,
        "left_index": False,
        "right_index": False,
        "left_by": None,
        "right_by": None,
        "suffixes": ("_x", "_y"),
        "tolerance": None,
        "allow_exact_matches": True,
        "direction": "backward",
    }

    @functools.cached_property
    def _kwargs(self):
        return {
            "left_on": self.left_on,
            "right_on": self.right_on,
            "left_index": self.left_index,
            "right_index": self.right_index,
            "left_by": self.left_by,
            "right_by": self.right_by,
            "suffixes": self.suffixes,
            "tolerance": self.tolerance,
            "allow_exact_matches": self.allow_exact_matches,
            "direction": self.direction,
        }

    @functools.cached_property
    def _left(self):
        left = self.left
        if self.left_on is not None:
            if self.right_index:
                left = ResetIndex(left)
            return new_collection(left).set_index(self.left_on, sorted=True)
        return left

    def _divisions(self):
        if (self.left_on or self.right_on) and (
            not self.right_index or not self.left.known_divisions
        ):
            return (None,) * (self.left.npartitions + 1)
        elif self.left_on or self.right_on:
            return self.left.divisions
        return self._left.divisions

    @functools.cached_property
    def _meta(self):
        return make_meta(
            pd.merge_asof(
                meta_nonempty(self.left._meta),
                meta_nonempty(self.right._meta),
                **self._kwargs,
            )
        )

    def _lower(self):
        left = self._left
        right = self.right
        left_on = self.left_on
        right_on = self.right_on
        right_index = self.right_index
        left_by = self.left_by
        right_by = self.right_by

        ixname = ixcol = divs = None
        if left_on is not None:
            if right_index:
                divs = self.left.divisions if self.left.known_divisions else None
                ixname = self.left.index.name
                ixcol = left.columns[0]

        if right_on is not None:
            right = (
                new_collection(right)
                .set_index(right_on, drop=(left_on == right_on), sorted=True)
                .expr
            )

        if (
            not left.known_divisions
            and left.npartitions > 1
            or not right.known_divisions
            and right.npartitions > 1
        ):
            raise ValueError("merge_asof input must be sorted!")

        left_index, right_index = True, True

        if left.npartitions == right.npartitions == 1:
            return MapPartitions(
                self.left,
                pd.merge_asof,
                self._meta,
                True,
                True,
                False,
                True,
                None,
                None,
                self._kwargs,
                self.right,
            )

        if all(map(pd.isnull, left.divisions)):
            return FromPandas(
                _BackendData(self._meta),
                npartitions=left.npartitions,
                pyarrow_strings_enabled=pyarrow_strings_enabled(),
            )

        if all(map(pd.isnull, right.divisions)):
            return MapPartitions(
                left,
                pd.merge_asof,
                self._meta,
                True,
                True,
                False,
                True,
                None,
                None,
                {"left_index": True, "right_index": True},
                right,
            )

        result = MergeAsofIndexed(
            left,
            right,
            left_index,
            right_index,
            left_by,
            right_by,
            self.suffixes,
            self.tolerance,
            self.allow_exact_matches,
            self.direction,
        )

        if left_on or right_on:
            result = ResetIndex(result)
            if ixcol is not None:
                if divs is not None:
                    result = SetIndexBlockwise(result, ixcol, new_divisions=divs)
                else:
                    result = SetIndexBlockwise(result, ixcol)
                result = RenameAxis(result, ixname)

        return result


class MergeAsofIndexed(MergeAsof):
    _parameters = [
        "left",
        "right",
        "left_index",
        "right_index",
        "left_by",
        "right_by",
        "suffixes",
        "tolerance",
        "allow_exact_matches",
        "direction",
    ]

    def _divisions(self):
        return self.left.divisions

    @functools.cached_property
    def left_on(self):
        # For optimisation
        return self.left_by

    @functools.cached_property
    def right_on(self):
        # For optimisation
        return self.right_by

    @functools.cached_property
    def _kwargs(self):
        return {
            "left_index": self.left_index,
            "right_index": self.right_index,
            "suffixes": self.suffixes,
            "left_by": self.left_by,
            "right_by": self.right_by,
            "tolerance": self.tolerance,
            "allow_exact_matches": self.allow_exact_matches,
            "direction": self.direction,
        }

    @functools.cached_property
    def _meta(self):
        return make_meta(
            pd.merge_asof(
                meta_nonempty(self.left._meta),
                meta_nonempty(self.right._meta),
                **self._kwargs,
            )
        )

    def _lower(self):
        return None

    def _layer(self) -> dict:
        dsk = dict()
        tails = heads = None
        tails_name = "prefix-reduction-" + self._name
        heads_name = "suffix_reduction-" + self._name
        if self.direction in ["backward", "nearest"]:
            tails = compute_tails(self.right, tails_name, by=self.right_by)
            dsk.update(tails)
        if self.direction in ["forward", "nearest"]:
            heads = compute_heads(self.right, heads_name, by=self.right_by)
            dsk.update(heads)

        for i, J in enumerate(
            pair_partitions(self.left.divisions, self.right.divisions)
        ):
            frames = []
            for j, lower, upper in J:
                slice = (
                    methods.boundary_slice,
                    (self.left._name, i),
                    lower,
                    upper,
                    False,
                )
                tail = (tails_name, j) if tails is not None else None
                head = (heads_name, j) if heads is not None else None
                frames.append(
                    (
                        apply,
                        merge_asof_padded,
                        [slice, (self.right._name, j), tail, head],
                        self._kwargs,
                    )
                )
            dsk[(self._name, i)] = (methods.concat, frames)
        return dsk


def most_recent_tail(left, right):
    if len(right.index) == 0:
        return left
    return right.tail(1)


def most_recent_tail_summary(left, right, by=None):
    return pd.concat([left, right]).drop_duplicates(subset=by, keep="last")


def compute_tails(ddf, name, by=None) -> dict:
    """For each partition, returns the last row of the most recent nonempty
    partition.
    """
    empty = ddf._meta.iloc[0:0]

    if by is None:
        return prefix_reduction(most_recent_tail, ddf, empty, name)
    else:
        kwargs = {"by": by}
        return prefix_reduction(most_recent_tail_summary, ddf, empty, name, **kwargs)


def prefix_reduction(f, ddf, identity, name, **kwargs):
    """Computes the prefix sums of f on df

    If df has partitions [P1, P2, ..., Pn], then returns the DataFrame with
    partitions [f(identity, P1),
                f(f(identity, P1), P2),
                f(f(f(identity, P1), P2), P3),
                ...]

    Parameters
    ----------
    f : callable
        an associative function f
    ddf : dd.DataFrame
    identity : pd.DataFrame
        an identity element of f, that is f(identity, df) = f(df, identity) = df
    """
    dsk = dict()
    n = len(ddf.divisions) - 1

    N = 1
    while N < n:
        N *= 2
    for i in range(n):
        dsk[(name, i, 1, 0)] = (apply, f, [(ddf._name, i), identity], kwargs)
    for i in range(n, N):
        dsk[(name, i, 1, 0)] = identity

    d = 1
    while d < N:
        for i in range(0, N, 2 * d):
            dsk[(name, i + 2 * d - 1, 2 * d, 0)] = (
                apply,
                f,
                [(name, i + d - 1, d, 0), (name, i + 2 * d - 1, d, 0)],
                kwargs,
            )
        d *= 2

    dsk[(name, N - 1, N, 1)] = identity

    while d > 1:
        d //= 2
        for i in range(0, N, 2 * d):
            dsk[(name, i + d - 1, d, 1)] = (name, i + 2 * d - 1, 2 * d, 1)
            dsk[(name, i + 2 * d - 1, d, 1)] = (
                apply,
                f,
                [(name, i + 2 * d - 1, 2 * d, 1), (name, i + d - 1, d, 0)],
                kwargs,
            )

    for i in range(n):
        dsk[(name, i)] = (apply, f, [(name, i, 1, 1), identity], kwargs)

    return dsk


def most_recent_head(left, right):
    if len(left.index) == 0:
        return right
    return left.head(1)


def most_recent_head_summary(left, right, by=None):
    return pd.concat([left, right]).drop_duplicates(subset=by, keep="first")


def compute_heads(ddf, name, by=None) -> dict:
    """For each partition, returns the first row of the next nonempty
    partition.
    """
    empty = ddf._meta.iloc[0:0]

    if by is None:
        return suffix_reduction(most_recent_head, ddf, empty, name)
    else:
        kwargs = {"by": by}
        return suffix_reduction(most_recent_head_summary, ddf, empty, name, **kwargs)


def suffix_reduction(f, ddf, identity, name, **kwargs):
    """Computes the suffix sums of f on df

    If df has partitions [P1, P2, ..., Pn], then returns the DataFrame with
    partitions [f(P1, f(P2, ...f(Pn, identity)...)),
                f(P2, ...f(Pn, identity)...),
                ...f(Pn, identity)...,
                ...]

    Parameters
    ----------
    f : callable
        an associative function f
    ddf : dd.DataFrame
    identity : pd.DataFrame
        an identity element of f, that is f(identity, df) = f(df, identity) = df
    kwargs : ??
        keyword arguments of f ??
    """
    dsk = dict()
    n = len(ddf.divisions) - 1

    N = 1
    while N < n:
        N *= 2
    for i in range(n):
        dsk[(name, i, 1, 0)] = (apply, f, [(ddf._name, n - 1 - i), identity], kwargs)
    for i in range(n, N):
        dsk[(name, i, 1, 0)] = identity

    d = 1
    while d < N:
        for i in range(0, N, 2 * d):
            dsk[(name, i + 2 * d - 1, 2 * d, 0)] = (
                apply,
                f,
                [(name, i + 2 * d - 1, d, 0), (name, i + d - 1, d, 0)],
                kwargs,
            )
        d *= 2

    dsk[(name, N - 1, N, 1)] = identity

    while d > 1:
        d //= 2
        for i in range(0, N, 2 * d):
            dsk[(name, i + d - 1, d, 1)] = (name, i + 2 * d - 1, 2 * d, 1)
            dsk[(name, i + 2 * d - 1, d, 1)] = (
                apply,
                f,
                [(name, i + d - 1, d, 0), (name, i + 2 * d - 1, 2 * d, 1)],
                kwargs,
            )

    for i in range(n):
        dsk[(name, i)] = (apply, f, [(name, n - 1 - i, 1, 1), identity], kwargs)

    return dsk
