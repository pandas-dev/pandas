from __future__ import annotations

import functools
from collections import namedtuple
from numbers import Integral

import pandas as pd
from pandas.core.window import Rolling as pd_Rolling

from dask.dataframe.dask_expr._collection import new_collection
from dask.dataframe.dask_expr._expr import (
    Blockwise,
    Expr,
    MapOverlap,
    Projection,
    determine_column_projection,
    make_meta,
)
from dask.utils import derived_from

BlockwiseDep = namedtuple(typename="BlockwiseDep", field_names=["iterable"])  # type: ignore


def _rolling_agg(
    frame,
    window,
    kwargs,
    how,
    how_args,
    how_kwargs,
    groupby_kwargs=None,
    groupby_slice=None,
):
    if groupby_kwargs is not None:
        frame = frame.groupby(**groupby_kwargs)
        if groupby_slice:
            frame = frame[groupby_slice]
    rolling = frame.rolling(window, **kwargs)
    result = getattr(rolling, how)(*how_args, **(how_kwargs or {}))
    if groupby_kwargs is not None:
        return result.sort_index(level=-1)
    return result


class RollingReduction(Expr):
    _parameters = [
        "frame",
        "window",
        "kwargs",
        "how_args",
        "how_kwargs",
        "groupby_kwargs",
        "groupby_slice",
    ]
    _defaults = {
        "kwargs": None,
        "how_args": (),
        "how_kwargs": None,
        "groupby_kwargs": None,
        "groupby_slice": None,
    }
    how: str | None = None

    @functools.cached_property
    def npartitions(self):
        return self.frame.npartitions

    def _divisions(self):
        return self.frame.divisions

    @functools.cached_property
    def _meta(self):
        meta = _rolling_agg(
            self.frame._meta,
            window=self.window,
            kwargs=self.kwargs,
            how=self.how,
            how_args=self.how_args,
            how_kwargs=self.how_kwargs,
            groupby_kwargs=self.groupby_kwargs,
            groupby_slice=self.groupby_slice,
        )
        return make_meta(meta)

    @functools.cached_property
    def kwargs(self):
        return {} if self.operand("kwargs") is None else self.operand("kwargs")

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Projection):
            by = self.groupby_kwargs.get("by", []) if self.groupby_kwargs else []
            by_columns = by if not isinstance(by, Expr) else []
            columns = determine_column_projection(self, parent, dependents, by_columns)
            columns = [col for col in self.frame.columns if col in columns]
            if columns == self.frame.columns:
                return
            if self.groupby_kwargs is not None:
                return type(parent)(
                    type(self)(self.frame[columns], *self.operands[1:]),
                    *parent.operands[1:],
                )
            if len(columns) == 1:
                columns = columns[0]
            return type(self)(self.frame[columns], *self.operands[1:])

    @property
    def _is_blockwise_op(self):
        return (
            self.kwargs.get("axis") in (1, "columns")
            or (isinstance(self.window, Integral) and self.window <= 1)
            or self.frame.npartitions == 1
        )

    def _lower(self):
        if self._is_blockwise_op:
            return RollingAggregation(
                self.frame,
                self.window,
                self.kwargs,
                self.how,
                list(self.how_args),
                self.how_kwargs,
                groupby_kwargs=self.groupby_kwargs,
                groupby_slice=self.groupby_slice,
            )

        if self.kwargs.get("center"):
            before = self.window // 2
            after = self.window - before - 1
        elif not isinstance(self.window, int):
            before = pd.Timedelta(self.window)
            after = 0
        else:
            before = self.window - 1
            after = 0

        return MapOverlap(
            frame=self.frame,
            func=_rolling_agg,
            before=before,
            after=after,
            meta=self._meta,
            enforce_metadata=True,
            kwargs=dict(
                window=self.window,
                kwargs=self.kwargs,
                how=self.how,
                how_args=self.how_args,
                how_kwargs=self.how_kwargs,
                groupby_kwargs=self.groupby_kwargs,
                groupby_slice=self.groupby_slice,
            ),
        )


class RollingAggregation(Blockwise):
    _parameters = [
        "frame",
        "window",
        "kwargs",
        "how",
        "how_args",
        "how_kwargs",
        "groupby_kwargs",
        "groupby_slice",
    ]

    operation = staticmethod(_rolling_agg)

    @functools.cached_property
    def _meta(self):
        return self.frame._meta


class RollingCount(RollingReduction):
    how = "count"


class RollingSum(RollingReduction):
    how = "sum"


class RollingMean(RollingReduction):
    how = "mean"


class RollingMin(RollingReduction):
    how = "min"


class RollingMax(RollingReduction):
    how = "max"


class RollingVar(RollingReduction):
    how = "var"


class RollingStd(RollingReduction):
    how = "std"


class RollingMedian(RollingReduction):
    how = "median"


class RollingQuantile(RollingReduction):
    how = "quantile"


class RollingSkew(RollingReduction):
    how = "skew"


class RollingKurt(RollingReduction):
    how = "kurt"


class RollingAgg(RollingReduction):
    how = "agg"

    def _simplify_up(self, parent, dependents):
        # Disable optimization in `agg`; function may access other columns
        return


class RollingApply(RollingReduction):
    how = "apply"


class RollingCov(RollingReduction):
    how = "cov"


class Rolling:
    """Aggregate using one or more operations

    The purpose of this class is to expose an API similar
    to Pandas' `Rolling` for dask-expr
    """

    def __init__(
        self,
        obj,
        window,
        groupby_kwargs=None,
        groupby_slice=None,
        min_periods=None,
        center=False,
        win_type=None,
    ):
        if obj.divisions[0] is None and len(obj.divisions) > 2:
            msg = (
                "Can only rolling dataframes with known divisions\n"
                "See https://docs.dask.org/en/latest/dataframe-design.html#partitions\n"
                "for more information."
            )
            raise ValueError(msg)
        self.obj = obj
        self.window = window
        self.groupby_kwargs = groupby_kwargs
        self.groupby_slice = groupby_slice
        self.min_periods = min_periods
        self.center = center
        self.win_type = win_type

        # Allow pandas to raise if appropriate
        obj._meta.rolling(window, **self.kwargs)

    @functools.cached_property
    def kwargs(self):
        return dict(
            min_periods=self.min_periods, center=self.center, win_type=self.win_type
        )

    def _single_agg(self, expr_cls, how_args=(), how_kwargs=None):
        return new_collection(
            expr_cls(
                self.obj,
                self.window,
                kwargs=self.kwargs,
                how_args=how_args,
                how_kwargs=how_kwargs,
                groupby_kwargs=self.groupby_kwargs,
                groupby_slice=self.groupby_slice,
            )
        )

    @derived_from(pd_Rolling)
    def cov(self):
        return self._single_agg(RollingCov)

    @derived_from(pd_Rolling)
    def apply(self, func, *args, **kwargs):
        return self._single_agg(RollingApply, how_args=(func, *args), how_kwargs=kwargs)

    @derived_from(pd_Rolling)
    def count(self):
        return self._single_agg(RollingCount)

    @derived_from(pd_Rolling)
    def sum(self):
        return self._single_agg(RollingSum)

    @derived_from(pd_Rolling)
    def mean(self):
        return self._single_agg(RollingMean)

    @derived_from(pd_Rolling)
    def min(self):
        return self._single_agg(RollingMin)

    @derived_from(pd_Rolling)
    def max(self):
        return self._single_agg(RollingMax)

    @derived_from(pd_Rolling)
    def var(self):
        return self._single_agg(RollingVar)

    @derived_from(pd_Rolling)
    def std(self):
        return self._single_agg(RollingStd)

    @derived_from(pd_Rolling)
    def median(self):
        return self._single_agg(RollingMedian)

    @derived_from(pd_Rolling)
    def quantile(self, q):
        return self._single_agg(RollingQuantile, how_args=(q,))

    @derived_from(pd_Rolling)
    def skew(self):
        return self._single_agg(RollingSkew)

    @derived_from(pd_Rolling)
    def kurt(self):
        return self._single_agg(RollingKurt)

    @derived_from(pd_Rolling)
    def agg(self, func, *args, **kwargs):
        return self._single_agg(RollingAgg, how_args=(func, *args), how_kwargs=kwargs)
