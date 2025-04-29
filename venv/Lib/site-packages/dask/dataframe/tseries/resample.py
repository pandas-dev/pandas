from __future__ import annotations

import functools
from collections import namedtuple

import numpy as np
import pandas as pd
from pandas.core.resample import Resampler as pd_Resampler

from dask.dataframe import methods
from dask.dataframe.dask_expr._collection import new_collection
from dask.dataframe.dask_expr._expr import (
    Blockwise,
    Expr,
    Projection,
    make_meta,
    plain_column_projection,
)
from dask.dataframe.dask_expr._repartition import Repartition
from dask.dataframe.dispatch import meta_nonempty
from dask.utils import derived_from

BlockwiseDep = namedtuple(typename="BlockwiseDep", field_names=["iterable"])  # type: ignore


def _resample_series(
    series,
    start,
    end,
    reindex_closed,
    rule,
    resample_kwargs,
    how,
    fill_value,
    how_args,
    how_kwargs,
):
    out = getattr(series.resample(rule, **resample_kwargs), how)(
        *how_args, **how_kwargs
    )
    if reindex_closed is None:
        inclusive = "both"
    else:
        inclusive = reindex_closed
    closed_kwargs = {"inclusive": inclusive}

    new_index = pd.date_range(
        start.tz_localize(None),
        end.tz_localize(None),
        freq=rule,
        **closed_kwargs,
        name=out.index.name,
        unit=out.index.unit,
    ).tz_localize(start.tz, nonexistent="shift_forward")

    if not out.index.isin(new_index).all():
        raise ValueError(
            "Index is not contained within new index. This can often be "
            "resolved by using larger partitions, or unambiguous "
            "frequencies: 'Q', 'A'..."
        )

    return out.reindex(new_index, fill_value=fill_value)


def _resample_bin_and_out_divs(divisions, rule, closed="left", label="left"):
    rule = pd.tseries.frequencies.to_offset(rule)
    g = pd.Grouper(freq=rule, how="count", closed=closed, label=label)

    # Determine bins to apply `how` to. Disregard labeling scheme.
    divs = pd.Series(range(len(divisions)), index=divisions)
    temp = divs.resample(rule, closed=closed, label="left").count()
    tempdivs = temp.loc[temp > 0].index

    # Cleanup closed == 'right' and label == 'right'
    res = pd.offsets.Nano() if isinstance(rule, pd.offsets.Tick) else pd.offsets.Day()
    if g.closed == "right":
        newdivs = tempdivs + res
    else:
        newdivs = tempdivs
    if g.label == "right":
        outdivs = tempdivs + rule
    else:
        outdivs = tempdivs

    newdivs = methods.tolist(newdivs)
    outdivs = methods.tolist(outdivs)

    # Adjust ends
    if newdivs[0] < divisions[0]:
        newdivs[0] = divisions[0]
    if newdivs[-1] < divisions[-1]:
        if len(newdivs) < len(divs):
            setter = lambda a, val: a.append(val)
        else:
            setter = lambda a, val: a.__setitem__(-1, val)
        setter(newdivs, divisions[-1] + res)
        if outdivs[-1] > divisions[-1]:
            setter(outdivs, outdivs[-1])
        elif outdivs[-1] < divisions[-1]:
            setter(outdivs, temp.index[-1])

    return tuple(map(pd.Timestamp, newdivs)), tuple(map(pd.Timestamp, outdivs))


class ResampleReduction(Expr):
    _parameters = [
        "frame",
        "rule",
        "kwargs",
        "how_args",
        "how_kwargs",
    ]
    _defaults = {
        "closed": None,
        "label": None,
        "kwargs": None,
        "how_args": (),
        "how_kwargs": None,
    }
    how: str | None = None
    fill_value = np.nan

    @functools.cached_property
    def npartitions(self):
        return self.frame.npartitions

    def _divisions(self):
        return self._resample_divisions[1]

    @functools.cached_property
    def _meta(self):
        resample = meta_nonempty(self.frame._meta).resample(self.rule, **self.kwargs)
        meta = getattr(resample, self.how)(*self.how_args, **self.how_kwargs or {})
        return make_meta(meta)

    @functools.cached_property
    def kwargs(self):
        return {} if self.operand("kwargs") is None else self.operand("kwargs")

    @functools.cached_property
    def how_kwargs(self):
        return {} if self.operand("how_kwargs") is None else self.operand("how_kwargs")

    @functools.cached_property
    def _resample_divisions(self):
        return _resample_bin_and_out_divs(
            self.frame.divisions, self.rule, **self.kwargs or {}
        )

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Projection):
            return plain_column_projection(self, parent, dependents)

    def _lower(self):
        partitioned = Repartition(
            self.frame, new_divisions=self._resample_divisions[0], force=True
        )
        output_divisions = self._resample_divisions[1]
        return ResampleAggregation(
            partitioned,
            BlockwiseDep(output_divisions[:-1]),
            BlockwiseDep(output_divisions[1:]),
            BlockwiseDep(["left"] * (len(output_divisions[1:]) - 1) + [None]),
            self.rule,
            self.kwargs,
            self.how,
            self.fill_value,
            list(self.how_args),
            self.how_kwargs,
        )


class ResampleAggregation(Blockwise):
    _parameters = [
        "frame",
        "divisions_left",
        "divisions_right",
        "closed",
        "rule",
        "kwargs",
        "how",
        "fill_value",
        "how_args",
        "how_kwargs",
    ]
    operation = staticmethod(_resample_series)

    @functools.cached_property
    def _meta(self):
        return self.frame._meta

    def _divisions(self):
        return list(self.divisions_left.iterable) + [self.divisions_right.iterable[-1]]

    def _blockwise_arg(self, arg, i):
        if isinstance(arg, BlockwiseDep):
            return arg.iterable[i]
        return super()._blockwise_arg(arg, i)


class ResampleCount(ResampleReduction):
    how = "count"
    fill_value = 0


class ResampleSum(ResampleReduction):
    how = "sum"
    fill_value = 0


class ResampleProd(ResampleReduction):
    how = "prod"


class ResampleMean(ResampleReduction):
    how = "mean"


class ResampleMin(ResampleReduction):
    how = "min"


class ResampleMax(ResampleReduction):
    how = "max"


class ResampleFirst(ResampleReduction):
    how = "first"


class ResampleLast(ResampleReduction):
    how = "last"


class ResampleVar(ResampleReduction):
    how = "var"


class ResampleStd(ResampleReduction):
    how = "std"


class ResampleSize(ResampleReduction):
    how = "size"
    fill_value = 0


class ResampleNUnique(ResampleReduction):
    how = "nunique"
    fill_value = 0


class ResampleMedian(ResampleReduction):
    how = "median"


class ResampleQuantile(ResampleReduction):
    how = "quantile"


class ResampleOhlc(ResampleReduction):
    how = "ohlc"


class ResampleSem(ResampleReduction):
    how = "sem"


class ResampleAgg(ResampleReduction):
    how = "agg"

    def _simplify_up(self, parent, dependents):
        # Disable optimization in `agg`; function may access other columns
        return


class Resampler:
    """Aggregate using one or more operations

    The purpose of this class is to expose an API similar
    to Pandas' `Resampler` for dask-expr
    """

    def __init__(self, obj, rule, **kwargs):
        if obj.divisions[0] is None:
            msg = (
                "Can only resample dataframes with known divisions\n"
                "See https://docs.dask.org/en/latest/dataframe-design.html#partitions\n"
                "for more information."
            )
            raise ValueError(msg)
        self.obj = obj
        self.rule = rule
        self.kwargs = kwargs

    def _single_agg(self, expr_cls, how_args=(), how_kwargs=None):
        return new_collection(
            expr_cls(
                self.obj,
                self.rule,
                self.kwargs,
                how_args=how_args,
                how_kwargs=how_kwargs,
            )
        )

    @derived_from(pd_Resampler)
    def count(self):
        return self._single_agg(ResampleCount)

    @derived_from(pd_Resampler)
    def sum(self):
        return self._single_agg(ResampleSum)

    @derived_from(pd_Resampler)
    def prod(self):
        return self._single_agg(ResampleProd)

    @derived_from(pd_Resampler)
    def mean(self):
        return self._single_agg(ResampleMean)

    @derived_from(pd_Resampler)
    def min(self):
        return self._single_agg(ResampleMin)

    @derived_from(pd_Resampler)
    def max(self):
        return self._single_agg(ResampleMax)

    @derived_from(pd_Resampler)
    def first(self):
        return self._single_agg(ResampleFirst)

    @derived_from(pd_Resampler)
    def last(self):
        return self._single_agg(ResampleLast)

    @derived_from(pd_Resampler)
    def var(self):
        return self._single_agg(ResampleVar)

    @derived_from(pd_Resampler)
    def std(self):
        return self._single_agg(ResampleStd)

    @derived_from(pd_Resampler)
    def size(self):
        return self._single_agg(ResampleSize)

    @derived_from(pd_Resampler)
    def nunique(self):
        return self._single_agg(ResampleNUnique)

    @derived_from(pd_Resampler)
    def median(self):
        return self._single_agg(ResampleMedian)

    @derived_from(pd_Resampler)
    def quantile(self):
        return self._single_agg(ResampleQuantile)

    @derived_from(pd_Resampler)
    def ohlc(self):
        return self._single_agg(ResampleOhlc)

    @derived_from(pd_Resampler)
    def sem(self):
        return self._single_agg(ResampleSem)

    @derived_from(pd_Resampler)
    def agg(self, func, *args, **kwargs):
        return self._single_agg(ResampleAgg, how_args=(func, *args), how_kwargs=kwargs)
