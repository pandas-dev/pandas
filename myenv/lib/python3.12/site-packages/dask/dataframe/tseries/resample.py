from __future__ import annotations

import numpy as np
import pandas as pd
from pandas.core.resample import Resampler as pd_Resampler

from dask.base import tokenize
from dask.dataframe import methods
from dask.dataframe._compat import PANDAS_GE_140
from dask.dataframe.core import DataFrame, Series
from dask.highlevelgraph import HighLevelGraph
from dask.utils import derived_from


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

    if PANDAS_GE_140:
        if reindex_closed is None:
            inclusive = "both"
        else:
            inclusive = reindex_closed
        closed_kwargs = {"inclusive": inclusive}
    else:
        closed_kwargs = {"closed": reindex_closed}

    new_index = pd.date_range(
        start.tz_localize(None),
        end.tz_localize(None),
        freq=rule,
        **closed_kwargs,
        name=out.index.name,
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


class Resampler:
    """Class for resampling timeseries data.

    This class is commonly encountered when using ``obj.resample(...)`` which
    return ``Resampler`` objects.

    Parameters
    ----------
    obj : Dask DataFrame or Series
        Data to be resampled.
    rule : str, tuple, datetime.timedelta, DateOffset or None
        The offset string or object representing the target conversion.
    kwargs : optional
        Keyword arguments passed to underlying pandas resampling function.

    Returns
    -------
    Resampler instance of the appropriate type
    """

    def __init__(self, obj, rule, **kwargs):
        if not obj.known_divisions:
            msg = (
                "Can only resample dataframes with known divisions\n"
                "See https://docs.dask.org/en/latest/dataframe-design.html#partitions\n"
                "for more information."
            )
            raise ValueError(msg)
        self.obj = obj
        self._rule = pd.tseries.frequencies.to_offset(rule)
        self._kwargs = kwargs

    def _agg(
        self,
        how,
        meta=None,
        fill_value=np.nan,
        how_args=(),
        how_kwargs=None,
    ):
        """Aggregate using one or more operations

        Parameters
        ----------
        how : str
            Name of aggregation operation
        fill_value : scalar, optional
            Value to use for missing values, applied during upsampling.
            Default is NaN.
        how_args : optional
            Positional arguments for aggregation operation.
        how_kwargs : optional
            Keyword arguments for aggregation operation.

        Returns
        -------
        Dask DataFrame or Series
        """
        if how_kwargs is None:
            how_kwargs = {}

        rule = self._rule
        kwargs = self._kwargs
        name = "resample-" + tokenize(
            self.obj, rule, kwargs, how, *how_args, **how_kwargs
        )

        # Create a grouper to determine closed and label conventions
        newdivs, outdivs = _resample_bin_and_out_divs(
            self.obj.divisions, rule, **kwargs
        )

        # Repartition divs into bins. These won't match labels after mapping
        partitioned = self.obj.repartition(newdivs, force=True)

        keys = partitioned.__dask_keys__()
        dsk = {}

        args = zip(keys, outdivs, outdivs[1:], ["left"] * (len(keys) - 1) + [None])
        for i, (k, s, e, c) in enumerate(args):
            dsk[(name, i)] = (
                _resample_series,
                k,
                s,
                e,
                c,
                rule,
                kwargs,
                how,
                fill_value,
                list(how_args),
                how_kwargs,
            )

        # Infer output metadata
        meta_r = self.obj._meta_nonempty.resample(self._rule, **self._kwargs)
        meta = getattr(meta_r, how)(*how_args, **how_kwargs)

        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[partitioned])
        if isinstance(meta, pd.DataFrame):
            return DataFrame(graph, name, meta, outdivs)
        return Series(graph, name, meta, outdivs)

    @derived_from(pd_Resampler)
    def agg(self, agg_funcs, *args, **kwargs):
        return self._agg("agg", how_args=(agg_funcs,) + args, how_kwargs=kwargs)

    @derived_from(pd_Resampler)
    def count(self):
        return self._agg("count", fill_value=0)

    @derived_from(pd_Resampler)
    def first(self):
        return self._agg("first")

    @derived_from(pd_Resampler)
    def last(self):
        return self._agg("last")

    @derived_from(pd_Resampler)
    def mean(self):
        return self._agg("mean")

    @derived_from(pd_Resampler)
    def min(self):
        return self._agg("min")

    @derived_from(pd_Resampler)
    def median(self):
        return self._agg("median")

    @derived_from(pd_Resampler)
    def max(self):
        return self._agg("max")

    @derived_from(pd_Resampler)
    def nunique(self):
        return self._agg("nunique", fill_value=0)

    @derived_from(pd_Resampler)
    def ohlc(self):
        return self._agg("ohlc")

    @derived_from(pd_Resampler)
    def prod(self):
        return self._agg("prod")

    @derived_from(pd_Resampler)
    def sem(self):
        return self._agg("sem")

    @derived_from(pd_Resampler)
    def std(self):
        return self._agg("std")

    @derived_from(pd_Resampler)
    def size(self):
        return self._agg("size", fill_value=0)

    @derived_from(pd_Resampler)
    def sum(self):
        return self._agg("sum", fill_value=0)

    @derived_from(pd_Resampler)
    def var(self):
        return self._agg("var")

    @derived_from(pd_Resampler)
    def quantile(self):
        return self._agg("quantile")
