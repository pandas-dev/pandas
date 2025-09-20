from __future__ import annotations

import functools
import operator

import numpy as np
import pandas as pd

from dask._task_spec import Task
from dask.dataframe.dask_expr._collection import new_collection
from dask.dataframe.dask_expr._expr import ArrowStringConversion
from dask.dataframe.dask_expr.io import BlockwiseIO, PartitionsFiltered
from dask.dataframe.utils import pyarrow_strings_enabled
from dask.typing import Key

__all__ = ["timeseries"]


class Timeseries(PartitionsFiltered, BlockwiseIO):
    _absorb_projections = True

    _parameters = [
        "start",
        "end",
        "dtypes",
        "freq",
        "partition_freq",
        "seed",
        "kwargs",
        "columns",
        "_partitions",
        "_series",
    ]
    _defaults = {
        "start": "2000-01-01",
        "end": "2000-01-31",
        "dtypes": {"name": "string", "id": int, "x": float, "y": float},
        "freq": "1s",
        "partition_freq": "1d",
        "seed": None,
        "kwargs": {},
        "_partitions": None,
        "_series": False,
    }

    @functools.cached_property
    def _meta(self):
        result = self._make_timeseries_part("2000", "2000", 0).iloc[:0]
        if self._series:
            return result[result.columns[0]]
        return result

    def _divisions(self):
        return pd.date_range(start=self.start, end=self.end, freq=self.partition_freq)

    @property
    def _dtypes(self):
        dtypes = self.operand("dtypes")
        return {col: dtypes[col] for col in self.operand("columns")}

    @functools.cached_property
    def random_state(self):
        npartitions = len(self._divisions()) - 1
        ndtypes = max(len(self.operand("dtypes")), 1)
        random_state = np.random.RandomState(self.seed)
        n = npartitions * ndtypes
        random_data = random_state.bytes(n * 4)  # `n` 32-bit integers
        l = list(np.frombuffer(random_data, dtype=np.uint32).reshape((n,)))
        assert len(l) == n
        return l

    @functools.cached_property
    def _make_timeseries_part(self):
        return MakeTimeseriesPart(
            self.operand("dtypes"),
            list(self._dtypes.keys()),
            self.freq,
            self.kwargs,
        )

    def _filtered_task(self, name: Key, index: int) -> Task:
        full_divisions = self._divisions()
        ndtypes = max(len(self.operand("dtypes")), 1)
        task = Task(
            name,
            self._make_timeseries_part,
            full_divisions[index],
            full_divisions[index + 1],
            self.random_state[index * ndtypes],
        )
        if self._series:
            return Task(name, operator.getitem, task, self.operand("columns")[0])
        return task


names = [
    "Alice",
    "Bob",
    "Charlie",
    "Dan",
    "Edith",
    "Frank",
    "George",
    "Hannah",
    "Ingrid",
    "Jerry",
    "Kevin",
    "Laura",
    "Michael",
    "Norbert",
    "Oliver",
    "Patricia",
    "Quinn",
    "Ray",
    "Sarah",
    "Tim",
    "Ursula",
    "Victor",
    "Wendy",
    "Xavier",
    "Yvonne",
    "Zelda",
]


def make_string(n, rstate):
    return rstate.choice(names, size=n)


def make_categorical(n, rstate):
    return pd.Categorical.from_codes(rstate.randint(0, len(names), size=n), names)


def make_float(n, rstate):
    return rstate.rand(n) * 2 - 1


def make_int(n, rstate, lam=1000):
    return rstate.poisson(lam, size=n)


make = {
    float: make_float,
    int: make_int,
    str: make_string,
    object: make_string,
    "string": make_string,
    "category": make_categorical,
}


class MakeTimeseriesPart:
    def __init__(self, dtypes, columns, freq, kwargs):
        self.dtypes = dtypes
        self.columns = columns
        self.freq = freq
        self.kwargs = kwargs

    def __call__(self, start, end, state_data):
        dtypes = self.dtypes
        columns = self.columns
        freq = self.freq
        kwargs = self.kwargs
        state = np.random.RandomState(state_data)
        index = pd.date_range(start=start, end=end, freq=freq, name="timestamp")
        data = {}
        for k, dt in dtypes.items():
            kws = {
                kk.rsplit("_", 1)[1]: v
                for kk, v in kwargs.items()
                if kk.rsplit("_", 1)[0] == k
            }
            # Note: we compute data for all dtypes in order, not just those in the output
            # columns. This ensures the same output given the same state_data, regardless
            # of whether there is any column projection.
            # cf. https://github.com/dask/dask/pull/9538#issuecomment-1267461887
            result = make[dt](len(index), state, **kws)
            if k in columns:
                data[k] = result
        df = pd.DataFrame(data, index=index, columns=columns)
        if df.index[-1] == end:
            df = df.iloc[:-1]
        return df


def timeseries(
    start="2000-01-01",
    end="2000-01-31",
    freq="1s",
    partition_freq="1D",
    dtypes=None,
    seed=None,
    **kwargs,
):
    """Create timeseries dataframe with random data

    Parameters
    ----------
    start: datetime (or datetime-like string)
        Start of time series
    end: datetime (or datetime-like string)
        End of time series
    dtypes: dict (optional)
        Mapping of column names to types.
        Valid types include {float, int, str, 'category'}
    freq: string
        String like '2s' or '1H' or '12W' for the time series frequency
    partition_freq: string
        String like '1M' or '2Y' to divide the dataframe into partitions
    seed: int (optional)
        Randomstate seed
    kwargs:
        Keywords to pass down to individual column creation functions.
        Keywords should be prefixed by the column name and then an underscore.

    Examples
    --------
    >>> from dask.dataframe.dask_expr.datasets import timeseries
    >>> df = timeseries(
    ...     start='2000', end='2010',
    ...     dtypes={'value': float, 'name': str, 'id': int},
    ...     freq='2h', partition_freq='1D', seed=1
    ... )
    >>> df.head()  # doctest: +SKIP
                           id      name     value
    2000-01-01 00:00:00   969     Jerry -0.309014
    2000-01-01 02:00:00  1010       Ray -0.760675
    2000-01-01 04:00:00  1016  Patricia -0.063261
    2000-01-01 06:00:00   960   Charlie  0.788245
    2000-01-01 08:00:00  1031     Kevin  0.466002
    """
    if dtypes is None:
        dtypes = {"name": "string", "id": int, "x": float, "y": float}

    if seed is None:
        seed = np.random.randint(2e9)

    expr = Timeseries(
        start,
        end,
        dtypes,
        freq,
        partition_freq,
        seed,
        kwargs,
        columns=list(dtypes.keys()),
    )
    if pyarrow_strings_enabled():
        return new_collection(ArrowStringConversion(expr))
    return new_collection(expr)
