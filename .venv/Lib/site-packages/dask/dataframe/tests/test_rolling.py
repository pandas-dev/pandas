from __future__ import annotations

import datetime

import numpy as np
import pandas as pd
import pytest

import dask.dataframe as dd
from dask.dataframe._compat import PANDAS_GE_210
from dask.dataframe.utils import assert_eq

N = 40
df = pd.DataFrame(
    {
        "a": np.random.randn(N).cumsum(),
        "b": np.random.randint(100, size=(N,)),
        "c": np.random.randint(100, size=(N,)),
        "d": np.random.randint(100, size=(N,)),
        "e": np.random.randint(100, size=(N,)),
    }
)
ddf = dd.from_pandas(df, 3)

idx = (
    pd.date_range("2016-01-01", freq="3s", periods=100).union(
        pd.date_range("2016-01-01", freq="5s", periods=100)
    )
)[:N]

idx_constant_freq = (pd.date_range("2016-01-01", freq="1s", periods=100))[:N]

ts_data = {
    "a": np.random.randn(N).cumsum(),
    "b": np.random.randint(100, size=(N,)),
    "c": np.random.randint(100, size=(N,)),
    "d": np.random.randint(100, size=(N,)),
    "e": np.random.randint(100, size=(N,)),
}

ts = pd.DataFrame(
    ts_data,
    index=idx,
)

ts_constant_freq = pd.DataFrame(ts_data, index=idx_constant_freq)

dts = dd.from_pandas(ts, 3)


def shifted_sum(df, before, after, c=0):
    a = df.shift(before)
    b = df.shift(-after)
    return df + a + b + c


@pytest.mark.parametrize("npartitions", [1, 4])
@pytest.mark.parametrize("use_dask_input", [True, False])
def test_map_overlap(npartitions, use_dask_input):
    ddf = df
    if use_dask_input:
        ddf = dd.from_pandas(df, npartitions)

    for before, after in [(0, 3), (3, 0), (3, 3), (0, 0)]:
        # DataFrame
        res = dd.map_overlap(shifted_sum, ddf, before, after, before, after, c=2)
        sol = shifted_sum(df, before, after, c=2)
        assert_eq(res, sol)

        # Series
        res = dd.map_overlap(shifted_sum, ddf.b, before, after, before, after, c=2)
        sol = shifted_sum(df.b, before, after, c=2)
        assert_eq(res, sol)


@pytest.mark.parametrize("use_dask_input", [True, False])
@pytest.mark.parametrize("npartitions", [1, 4])
@pytest.mark.parametrize("enforce_metadata", [True, False])
@pytest.mark.parametrize("transform_divisions", [True, False])
@pytest.mark.parametrize("align_dataframes", [True, False])
@pytest.mark.parametrize(
    "overlap_setup",
    [
        (df, 0, 3),
        (df, 3, 0),
        (df, 3, 3),
        (df, 0, 0),
        (
            ts_constant_freq,
            datetime.timedelta(seconds=3),
            datetime.timedelta(seconds=3),
        ),
        (ts_constant_freq, datetime.timedelta(seconds=3), 0),
    ],
)
def test_map_overlap_multiple_dataframes(
    use_dask_input,
    npartitions,
    enforce_metadata,
    transform_divisions,
    align_dataframes,
    overlap_setup,
):
    dataframe, before, after = overlap_setup

    ddf = dataframe
    ddf2 = dataframe * 2
    if use_dask_input:
        ddf = dd.from_pandas(ddf, npartitions)
        ddf2 = dd.from_pandas(ddf2, 2 if align_dataframes else npartitions)

    def get_shifted_sum_arg(overlap):
        return (
            overlap.seconds - 1 if isinstance(overlap, datetime.timedelta) else overlap
        )

    before_shifted_sum, after_shifted_sum = get_shifted_sum_arg(
        before
    ), get_shifted_sum_arg(after)

    # DataFrame
    res = dd.map_overlap(
        shifted_sum,
        ddf,
        before,
        after,
        before_shifted_sum,
        after_shifted_sum,
        ddf2,
        align_dataframes=align_dataframes,
        transform_divisions=transform_divisions,
        enforce_metadata=enforce_metadata,
    )
    sol = shifted_sum(dataframe, before_shifted_sum, after_shifted_sum, dataframe * 2)
    assert_eq(res, sol)

    # Series
    res = dd.map_overlap(
        shifted_sum,
        ddf.b,
        before,
        after,
        before_shifted_sum,
        after_shifted_sum,
        ddf2.b,
        align_dataframes=align_dataframes,
        transform_divisions=transform_divisions,
        enforce_metadata=enforce_metadata,
    )
    sol = shifted_sum(
        dataframe.b, before_shifted_sum, after_shifted_sum, dataframe.b * 2
    )
    assert_eq(res, sol)


@pytest.mark.parametrize("npartitions", [1, 4])
@pytest.mark.parametrize("enforce_metadata", [True, False])
@pytest.mark.parametrize("transform_divisions", [True, False])
@pytest.mark.parametrize("align_dataframes", [True, False])
def test_map_overlap_names(
    npartitions, enforce_metadata, transform_divisions, align_dataframes
):
    ddf = dd.from_pandas(df, npartitions)

    res = ddf.map_overlap(
        shifted_sum,
        0,
        3,
        0,
        3,
        c=2,
        align_dataframes=align_dataframes,
        transform_divisions=transform_divisions,
        enforce_metadata=enforce_metadata,
    )
    res2 = ddf.map_overlap(
        shifted_sum,
        0,
        3,
        0,
        3,
        c=2,
        align_dataframes=align_dataframes,
        transform_divisions=transform_divisions,
        enforce_metadata=enforce_metadata,
    )
    assert set(res.dask) == set(res2.dask)

    res3 = ddf.map_overlap(
        shifted_sum,
        0,
        3,
        0,
        3,
        c=3,
        align_dataframes=align_dataframes,
        transform_divisions=transform_divisions,
        enforce_metadata=enforce_metadata,
    )
    assert res3._name != res._name
    # Difference is just the final map
    diff = res3.dask.keys() - res.dask.keys()
    assert len(diff) == npartitions

    res4 = ddf.map_overlap(shifted_sum, 3, 0, 0, 3, c=2)
    assert res4._name != res._name


def test_map_overlap_errors():
    # Non-integer
    with pytest.raises(ValueError):
        ddf.map_overlap(shifted_sum, 0.5, 3, 0, 2, c=2)

    # Negative
    with pytest.raises(ValueError):
        ddf.map_overlap(shifted_sum, 0, -5, 0, 2, c=2)

    # Partition size < window size
    with pytest.raises(NotImplementedError):
        ddf.map_overlap(shifted_sum, 0, 100, 0, 100, c=2).compute()

    # Timedelta offset with non-datetime
    with pytest.raises(TypeError):
        ddf.map_overlap(shifted_sum, pd.Timedelta("1s"), pd.Timedelta("1s"), 0, 2, c=2)

    # String timedelta offset with non-datetime
    with pytest.raises(TypeError):
        ddf.map_overlap(shifted_sum, "1s", "1s", 0, 2, c=2)


def test_map_overlap_provide_meta():
    df = pd.DataFrame(
        {"x": [1, 2, 4, 7, 11], "y": [1.0, 2.0, 3.0, 4.0, 5.0]}
    ).rename_axis("myindex")
    ddf = dd.from_pandas(df, npartitions=2)

    # Provide meta spec, but not full metadata
    res = ddf.map_overlap(
        lambda df: df.rolling(2).sum(), 2, 0, meta={"x": "i8", "y": "i8"}
    )
    sol = df.rolling(2).sum()
    assert_eq(res, sol)


def mad(x):
    return np.fabs(x - x.mean()).mean()


rolling_method_args_check_less_precise = [
    ("count", (), False),
    ("sum", (), False),
    ("mean", (), False),
    ("median", (), False),
    ("min", (), False),
    ("max", (), False),
    ("std", (), True),
    ("var", (), True),
    ("skew", (), True),  # here and elsewhere, results for kurt and skew are
    ("kurt", (), True),  # checked with check_less_precise=True so that we are
    # only looking at 3ish decimal places for the equality check
    # rather than 5ish. I have encountered a case where a test
    # seems to have failed due to numerical problems with kurt.
    # So far, I am only weakening the check for kurt and skew,
    # as they involve third degree powers and higher
    ("quantile", (0.38,), False),
    ("apply", (mad,), False),
]


@pytest.mark.parametrize(
    "method,args,check_less_precise", rolling_method_args_check_less_precise
)
@pytest.mark.parametrize("window", [1, 2, 4, 5])
@pytest.mark.parametrize("center", [True, False])
def test_rolling_methods(method, args, window, center, check_less_precise):
    if check_less_precise:
        check_less_precise = {"atol": 1e-3, "rtol": 1e-3}
    else:
        check_less_precise = {}
    if method == "count":
        min_periods = 0
    else:
        min_periods = None
    # DataFrame
    prolling = df.rolling(window, center=center, min_periods=min_periods)
    drolling = ddf.rolling(window, center=center, min_periods=min_periods)
    if method == "apply":
        kwargs = {"raw": False}
    else:
        kwargs = {}

    assert_eq(
        getattr(prolling, method)(*args, **kwargs),
        getattr(drolling, method)(*args, **kwargs),
        **check_less_precise,
    )

    # Series
    prolling = df.a.rolling(window, center=center, min_periods=min_periods)
    drolling = ddf.a.rolling(window, center=center, min_periods=min_periods)
    assert_eq(
        getattr(prolling, method)(*args, **kwargs),
        getattr(drolling, method)(*args, **kwargs),
        **check_less_precise,
    )


@pytest.mark.parametrize("window", [1, 2, 4, 5])
@pytest.mark.parametrize("center", [True, False])
def test_rolling_cov(window, center):
    # DataFrame
    prolling = df.drop("a", axis=1).rolling(window, center=center)
    drolling = ddf.drop("a", axis=1).rolling(window, center=center)
    assert_eq(prolling.cov(), drolling.cov())

    # Series
    prolling = df.b.rolling(window, center=center)
    drolling = ddf.b.rolling(window, center=center)
    assert_eq(prolling.cov(), drolling.cov())


def test_rolling_names():
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    a = dd.from_pandas(df, npartitions=2)
    assert sorted(a.rolling(2).sum().dask) == sorted(a.rolling(2).sum().dask)


def test_rolling_partition_size():
    df = pd.DataFrame(np.random.randn(50, 2))
    ddf = dd.from_pandas(df, npartitions=5)

    for obj, dobj in [(df, ddf), (df[0], ddf[0])]:
        assert_eq(obj.rolling(10).mean(), dobj.rolling(10).mean())
        assert_eq(obj.rolling(11).mean(), dobj.rolling(11).mean())
        with pytest.raises(NotImplementedError):
            dobj.rolling(12).mean().compute()


def test_time_rolling_constructor():
    result = dts.rolling("4s")
    assert result.window == "4s"
    assert result.min_periods is None
    assert result.win_type is None


@pytest.mark.parametrize(
    "method,args,check_less_precise", rolling_method_args_check_less_precise
)
@pytest.mark.parametrize("window", ["1s", "2s", "3s", pd.offsets.Second(5)])
def test_time_rolling_methods(method, args, window, check_less_precise):
    if check_less_precise:
        check_less_precise = {"atol": 1e-3, "rtol": 1e-3}
    else:
        check_less_precise = {}

    # DataFrame
    if method == "apply":
        kwargs = {"raw": False}
    else:
        kwargs = {}
    prolling = ts.rolling(window)
    drolling = dts.rolling(window)
    assert_eq(
        getattr(prolling, method)(*args, **kwargs),
        getattr(drolling, method)(*args, **kwargs),
        **check_less_precise,
    )

    # Series
    prolling = ts.a.rolling(window)
    drolling = dts.a.rolling(window)
    assert_eq(
        getattr(prolling, method)(*args, **kwargs),
        getattr(drolling, method)(*args, **kwargs),
        **check_less_precise,
    )


@pytest.mark.parametrize("window", ["1s", "2s", "3s", pd.offsets.Second(5)])
def test_time_rolling_cov(window):
    # DataFrame
    prolling = ts.drop("a", axis=1).rolling(window)
    drolling = dts.drop("a", axis=1).rolling(window)
    assert_eq(prolling.cov(), drolling.cov())

    # Series
    prolling = ts.b.rolling(window)
    drolling = dts.b.rolling(window)
    assert_eq(prolling.cov(), drolling.cov())


@pytest.mark.parametrize(
    "window,N",
    [("1s", 10), ("2s", 10), ("10s", 10), ("10h", 10), ("10s", 100), ("10h", 100)],
)
def test_time_rolling_large_window_fixed_chunks(window, N):
    df = pd.DataFrame(
        {
            "a": pd.date_range("2016-01-01 00:00:00", periods=N, freq="1s"),
            "b": np.random.randint(100, size=(N,)),
        }
    )
    df = df.set_index("a")
    ddf = dd.from_pandas(df, 5)
    assert_eq(ddf.rolling(window).sum(), df.rolling(window).sum())
    assert_eq(ddf.rolling(window).count(), df.rolling(window).count())
    assert_eq(ddf.rolling(window).mean(), df.rolling(window).mean())


@pytest.mark.parametrize("window", ["2s", "5s", "20s", "10h"])
def test_time_rolling_large_window_variable_chunks(window):
    df = pd.DataFrame(
        {
            "a": pd.date_range("2016-01-01 00:00:00", periods=100, freq="1s"),
            "b": np.random.randint(100, size=(100,)),
        }
    )
    ddf = dd.from_pandas(df, 5)
    ddf = ddf.repartition(divisions=[0, 5, 20, 28, 33, 54, 79, 80, 82, 99])
    df = df.set_index("a")
    ddf = ddf.set_index("a")
    assert_eq(ddf.rolling(window).sum(), df.rolling(window).sum())
    assert_eq(ddf.rolling(window).count(), df.rolling(window).count())
    assert_eq(ddf.rolling(window).mean(), df.rolling(window).mean())


@pytest.mark.parametrize("before, after", [("6s", "6s"), ("2s", "2s"), ("6s", "2s")])
def test_time_rolling(before, after):
    window = before
    expected = dts.compute().rolling(window).count()

    # String timedelta
    result = dts.map_overlap(lambda x: x.rolling(window).count(), before, after)
    assert_eq(result, expected)

    # Timedelta
    before = pd.Timedelta(before)
    after = pd.Timedelta(after)
    result = dts.map_overlap(lambda x: x.rolling(window).count(), before, after)
    assert_eq(result, expected)


def test_rolling_agg_aggregate():
    df = pd.DataFrame({"A": range(5), "B": range(0, 10, 2)})
    ddf = dd.from_pandas(df, npartitions=3)

    assert_eq(
        df.rolling(window=3).agg(["mean", "std"]),
        ddf.rolling(window=3).agg(["mean", "std"]),
    )

    assert_eq(
        df.rolling(window=3).agg({"A": "sum", "B": lambda x: np.std(x, ddof=1)}),
        ddf.rolling(window=3).agg({"A": "sum", "B": lambda x: np.std(x, ddof=1)}),
    )

    assert_eq(
        df.rolling(window=3).agg(["sum", "mean"]),
        ddf.rolling(window=3).agg(["sum", "mean"]),
    )

    assert_eq(
        df.rolling(window=3).agg({"A": ["sum", "mean"]}),
        ddf.rolling(window=3).agg({"A": ["sum", "mean"]}),
    )

    kwargs = {"raw": True}
    assert_eq(
        df.rolling(window=3).apply(lambda x: np.std(x, ddof=1), **kwargs),
        ddf.rolling(window=3).apply(lambda x: np.std(x, ddof=1), **kwargs),
    )


@pytest.mark.skipif(not PANDAS_GE_210, reason="buggy pandas implementation")
def test_rolling_numba_engine():
    pytest.importorskip("numba")
    df = pd.DataFrame({"A": range(5), "B": range(0, 10, 2)})
    ddf = dd.from_pandas(df, npartitions=3)

    def f(x):
        return np.sum(x) + 5

    assert_eq(
        df.rolling(3).apply(f, engine="numba", raw=True),
        ddf.rolling(3).apply(f, engine="numba", raw=True),
    )


def test_groupby_rolling():
    df = pd.DataFrame(
        {
            "column1": range(600),
            "group1": 5 * ["g" + str(i) for i in range(120)],
        },
        index=pd.date_range("20190101", periods=60).repeat(10),
    )

    ddf = dd.from_pandas(df, npartitions=8)

    expected = df.groupby("group1").rolling("15D").sum()
    actual = ddf.groupby("group1").rolling("15D").sum()

    assert_eq(expected, actual, check_divisions=False)

    expected = df.groupby("group1").column1.rolling("15D").mean()
    actual = ddf.groupby("group1").column1.rolling("15D").mean()

    assert_eq(expected, actual, check_divisions=False)
