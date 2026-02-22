from __future__ import annotations

import functools

import numpy as np
from pandas.core.dtypes.common import is_datetime64_any_dtype, is_timedelta64_dtype

from dask.dataframe._compat import PANDAS_GE_300
from dask.dataframe.dask_expr._expr import (
    Blockwise,
    DropnaSeries,
    Filter,
    Head,
    Sqrt,
    ToNumeric,
)
from dask.dataframe.dask_expr._quantile import SeriesQuantile
from dask.dataframe.dask_expr._reductions import Reduction, Size, ValueCounts
from dask.dataframe.dispatch import make_meta, meta_nonempty
from dask.dataframe.methods import (
    describe_nonnumeric_aggregate,
    describe_numeric_aggregate,
)


class DescribeNumeric(Reduction):
    _parameters = ["frame", "split_every", "percentiles", "percentile_method"]
    _defaults = {
        "percentiles": None,
        "split_every": None,
        "percentile_method": "default",
    }

    @functools.cached_property
    def _meta(self):
        return make_meta(meta_nonempty(self.frame._meta).describe())

    def _divisions(self):
        return (None, None)

    def _lower(self):
        frame = self.frame
        if self.percentiles is None:
            percentiles = self.percentiles or [0.25, 0.5, 0.75]
        else:
            percentiles = np.array(self.percentiles)
            if not PANDAS_GE_300:
                percentiles = np.append(percentiles, 0.5)
            percentiles = np.unique(percentiles)
            percentiles = list(percentiles)

        is_td_col = is_timedelta64_dtype(frame._meta.dtype)
        is_dt_col = is_datetime64_any_dtype(frame._meta.dtype)
        if is_td_col or is_dt_col:
            frame = ToNumeric(DropnaSeries(frame))

        stats = [
            frame.count(split_every=self.split_every),
            frame.mean(split_every=self.split_every),
            Sqrt(frame.var(split_every=self.split_every)),
            frame.min(split_every=self.split_every),
            SeriesQuantile(frame, q=percentiles, method=self.percentile_method),
            frame.max(split_every=self.split_every),
        ]
        try:
            unit = getattr(self.frame._meta.array, "unit", None)
        except AttributeError:
            # cudf Series has no array attribute
            unit = None
        return DescribeNumericAggregate(
            self.frame._meta.name,
            is_td_col,
            is_dt_col,
            unit,
            *stats,
        )


class DescribeNumericAggregate(Blockwise):
    _parameters = ["name", "is_timedelta_col", "is_datetime_col", "unit"]
    _defaults = {"is_timedelta_col": False, "is_datetime_col": False}

    def _broadcast_dep(self, dep):
        return dep.npartitions == 1

    @staticmethod
    def operation(name, is_timedelta_col, is_datetime_col, unit, *stats):
        return describe_numeric_aggregate(
            stats, name, is_timedelta_col, is_datetime_col, unit
        )


class DescribeNonNumeric(DescribeNumeric):
    _parameters = ["frame", "split_every"]

    def _lower(self):
        frame = self.frame
        vcounts = ValueCounts(frame, split_every=self.split_every, sort=True)
        count_unique = Size(Filter(vcounts, vcounts > 0))
        stats = [
            count_unique,
            frame.count(split_every=self.split_every),
            Head(vcounts, n=1),
        ]
        return DescribeNonNumericAggregate(frame._meta.name, *stats)


class DescribeNonNumericAggregate(DescribeNumericAggregate):
    _parameters = ["name"]
    _defaults = {}

    @staticmethod
    def operation(name, *stats):
        return describe_nonnumeric_aggregate(stats, name)
