from __future__ import annotations

from datetime import datetime

import numpy as np
import pytest

from dask.dataframe.dask_expr import from_pandas
from dask.dataframe.dask_expr.tests._util import _backend_library, assert_eq, xfail_gpu
from dask.utils import M

# Set DataFrame backend for this module
pd = _backend_library()


@pytest.fixture
def pdf():
    pdf = pd.DataFrame({"x": range(100)})
    pdf["y"] = pdf.x // 7  # Not unique; duplicates span different partitions
    yield pdf


@pytest.fixture
def df(pdf):
    yield from_pandas(pdf, npartitions=10)


def test_median(pdf, df):
    assert_eq(df.repartition(npartitions=1).median(axis=0), pdf.median(axis=0))
    assert_eq(df.repartition(npartitions=1).x.median(), pdf.x.median())
    assert_eq(df.x.median_approximate(), 49.0, atol=1)
    assert_eq(df.median_approximate(), pdf.median(), atol=1)

    # Ensure `median` redirects to `median_approximate` appropriately
    for axis in [None, 0, "rows"]:
        with pytest.raises(
            NotImplementedError, match="See the `median_approximate` method instead"
        ):
            df.median(axis=axis)

    with pytest.raises(
        NotImplementedError, match="See the `median_approximate` method instead"
    ):
        df.x.median()


def test_min_dt(pdf):
    pdf["dt"] = "a"
    df = from_pandas(pdf, npartitions=10)
    assert_eq(df.min(numeric_only=True), pdf.min(numeric_only=True))
    assert_eq(df.max(numeric_only=True), pdf.max(numeric_only=True))
    assert_eq(df.count(numeric_only=True), pdf.count(numeric_only=True))
    assert_eq(df.mean(numeric_only=True), pdf.mean(numeric_only=True))


@pytest.mark.parametrize(
    "series",
    [
        [0, 1, 0, 0, 1, 0],  # not monotonic
        [0, 1, 1, 2, 2, 3],  # monotonic increasing
        [0, 1, 2, 3, 4, 5],  # strictly monotonic increasing
        [0, 0, 0, 0, 0, 0],  # both monotonic increasing and monotonic decreasing
        [0, 1, 2, 0, 1, 2],  # Partitions are individually monotonic; whole series isn't
    ],
)
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("cls", ["Series", "Index"])
def test_monotonic(series, reverse, cls):
    if reverse:
        series = series[::-1]
    pds = pd.Series(series, index=series)
    ds = from_pandas(pds, 2, sort=False)
    if cls == "Index":
        pds = pds.index
        ds = ds.index

    # assert_eq fails due to numpy.bool vs. plain bool mismatch
    # See https://github.com/dask/dask/pull/10671
    assert ds.is_monotonic_increasing.compute() == pds.is_monotonic_increasing
    assert ds.is_monotonic_decreasing.compute() == pds.is_monotonic_decreasing


@pytest.mark.parametrize("split_every", [False, None, 5])
@pytest.mark.parametrize("split_out", [1, True])
def test_drop_duplicates(pdf, df, split_every, split_out):
    assert_eq(
        df.drop_duplicates(split_every=split_every, split_out=split_out),
        pdf.drop_duplicates(),
        check_index=split_out is not True,
    )
    assert_eq(
        df.x.drop_duplicates(split_every=split_every, split_out=split_out),
        pdf.x.drop_duplicates(),
        check_index=split_out is not True,
    )


def test_series_reduction_name():
    ser = from_pandas(pd.Series(range(10)), npartitions=2)
    df = ser.drop_duplicates().to_frame()
    assert_eq(df, df)


@pytest.mark.parametrize("split_every", [False, None, 5])
@pytest.mark.parametrize("split_out", [1, True])
def test_value_counts(pdf, df, split_every, split_out):
    assert_eq(
        df.x.value_counts(split_every=split_every, split_out=split_out),
        pdf.x.value_counts(),
        check_index=split_out is not True,
    )


def test_value_counts_sort(pdf):
    pdf = pdf.iloc[4:]
    df = from_pandas(pdf, npartitions=10)
    pd.testing.assert_series_equal(
        df.y.value_counts(split_out=1).compute(),
        pdf.y.value_counts(),
    )


@pytest.mark.parametrize("split_every", [None, 5])
@pytest.mark.parametrize("split_out", [1, True])
def test_unique(pdf, df, split_every, split_out):
    assert_eq(
        df.x.unique(split_every=split_every, split_out=split_out),
        pd.Series(pdf.x.unique(), name="x"),
        check_index=split_out is not True,
    )


@pytest.mark.parametrize(
    "reduction", ["sum", "prod", "min", "max", "any", "all", "count"]
)
@pytest.mark.parametrize(
    "split_every,expect_tasks", [(False, 21), (None, 23), (5, 23), (2, 31)]
)
def test_dataframe_split_every(pdf, df, split_every, expect_tasks, reduction):
    assert_eq(
        getattr(df, reduction)(split_every=split_every),
        getattr(pdf, reduction)(),
    )
    q = getattr(df, reduction)(split_every=split_every).optimize(fuse=False)
    assert len(q.__dask_graph__()) == expect_tasks


@pytest.mark.parametrize(
    "split_every, expect_tasks", [(False, 53), (None, 57), (5, 57), (2, 73)]
)
def test_dataframe_mode_split_every(pdf, df, split_every, expect_tasks):
    assert_eq(df.mode(split_every=split_every), pdf.mode())
    q = df.mode(split_every=split_every).optimize(fuse=False)
    assert len(q.__dask_graph__()) == expect_tasks


@pytest.mark.parametrize(
    "reduction", ["sum", "prod", "min", "max", "any", "all", "mode", "count"]
)
@pytest.mark.parametrize(
    "split_every,expect_tasks", [(False, 31), (None, 33), (5, 33), (2, 41)]
)
def test_series_split_every(pdf, df, split_every, expect_tasks, reduction):
    assert_eq(
        getattr(df.x, reduction)(split_every=split_every),
        getattr(pdf.x, reduction)(),
    )
    q = getattr(df.x, reduction)(split_every=split_every).optimize(fuse=False)
    assert len(q.__dask_graph__()) == expect_tasks


def test_reduction_with_strings():
    pdf = pd.DataFrame({"x": ["A", "B", "C"]})
    df = from_pandas(pdf, npartitions=2)
    assert_eq(
        df["x"].unique(), pd.Series(pdf["x"].unique(), name="x"), check_index=False
    )


@pytest.mark.parametrize("split_every", [-1, 0, 1])
@pytest.mark.parametrize(
    "reduction",
    [
        "sum",
        "prod",
        "min",
        "max",
        "any",
        "all",
        "mode",
        "count",
        "nunique_approx",
    ],
)
def test_split_every_lt2(df, reduction, split_every):
    with pytest.raises(ValueError, match="split_every must be greater than 1 or False"):
        # TODO validate parameters before graph materialization
        getattr(df.x, reduction)(split_every=split_every).__dask_graph__()


@pytest.mark.parametrize("split_every", [-1, 0, 1])
@pytest.mark.parametrize("reduction", ["drop_duplicates", "unique", "value_counts"])
def test_split_every_lt2_split_out(df, reduction, split_every):
    """split_out=True ignores split_every; force split_out=1"""
    with pytest.raises(ValueError, match="split_every must be greater than 1 or False"):
        # TODO validate parameters before graph materialization
        getattr(df.x, reduction)(split_out=1, split_every=split_every).__dask_graph__()


@pytest.mark.parametrize("split_every", [None, False, 2, 10])
def test_nunique_approx(split_every, pdf, df):
    approx = df.nunique_approx(split_every=split_every).compute()
    exact = len(df.drop_duplicates())
    assert abs(approx - exact) <= 2 or abs(approx - exact) / exact < 0.05


def test_unique_base(df, pdf):
    with pytest.raises(
        AttributeError, match="'DataFrame' object has no attribute 'unique'"
    ):
        df.unique()

    # pandas returns a numpy array while we return a Series/Index
    assert_eq(df.x.unique(), pd.Series(pdf.x.unique(), name="x"), check_index=False)
    assert_eq(df.index.unique(split_out=1), pd.Index(pdf.index.unique()))
    np.testing.assert_array_equal(
        df.index.unique().compute().sort_values().values,
        pd.Index(pdf.index.unique()).values,
    )


def test_value_counts_split_out_normalize(df, pdf):
    result = df.x.value_counts(split_out=2, normalize=True)
    expected = pdf.x.value_counts(normalize=True)
    assert_eq(result, expected)


@pytest.mark.parametrize("method", ["sum", "prod", "product"])
@pytest.mark.parametrize("min_count", [0, 9])
def test_series_agg_with_min_count(method, min_count):
    df = pd.DataFrame([[1]], columns=["a"])
    ddf = from_pandas(df, npartitions=1)
    func = getattr(ddf["a"], method)
    result = func(min_count=min_count).compute()
    if min_count == 0:
        assert result == 1
    else:
        # TODO: dtype is wrong
        assert_eq(result, np.nan, check_dtype=False)

    assert_eq(
        getattr(ddf, method)(min_count=min_count),
        getattr(df, method)(min_count=min_count),
    )


@pytest.mark.parametrize(
    "func",
    [
        M.max,
        M.min,
        M.any,
        M.all,
        M.sum,
        M.prod,
        M.product,
        M.count,
        M.mean,
        M.std,
        M.var,
        M.sem,
        pytest.param(
            M.idxmin, marks=xfail_gpu("https://github.com/rapidsai/cudf/issues/9602")
        ),
        pytest.param(
            M.idxmax, marks=xfail_gpu("https://github.com/rapidsai/cudf/issues/9602")
        ),
        pytest.param(
            lambda df: df.size,
            marks=pytest.mark.xfail(reason="scalars don't work yet"),
        ),
    ],
)
def test_reductions(func, pdf, df):
    result = func(df)
    assert result.known_divisions
    assert_eq(result, func(pdf))
    result = func(df.x)
    assert not result.known_divisions
    assert_eq(result, func(pdf.x))
    # check_dtype False because sub-selection of columns that is pushed through
    # is not reflected in the meta calculation
    assert_eq(func(df)["x"], func(pdf)["x"], check_dtype=False)

    result = func(df, axis=1)
    assert result.known_divisions
    assert_eq(result, func(pdf, axis=1))


def test_skew_kurt():
    pdf = pd.DataFrame(
        {
            "a": [1, 2, 6, 4, 4, 6, 4, 3, 7] * 1000,
            "b": [4, 2, 7, 3, 3, 1, 1, 1, 2] * 1000,
        },
    )
    df = from_pandas(pdf, npartitions=2)

    assert_eq(df.kurtosis().round(2), pdf.kurtosis().round(2))
    assert_eq(df.skew().round(2), pdf.skew().round(2))


@pytest.mark.parametrize(
    "func",
    [
        M.max,
        M.min,
        M.any,
        M.all,
        lambda idx: idx.size,
    ],
)
def test_index_reductions(func, pdf, df):
    result = func(df.index)
    assert not result.known_divisions
    assert_eq(result, func(pdf.index))


@pytest.mark.parametrize(
    "func",
    [
        lambda idx: idx.index,
        M.sum,
        M.prod,
        M.mean,
        M.std,
        M.var,
        M.idxmin,
        M.idxmax,
    ],
)
def test_unimplemented_on_index(func, pdf, df):
    # Methods/properties of Series that don't exist on Index
    with pytest.raises(AttributeError):
        func(pdf.index)
    with pytest.raises(AttributeError, match="'Index' object has no attribute '"):
        func(df.index)


def test_cov_corr(df, pdf):
    assert_eq(df.cov(), pdf.cov())
    assert_eq(df.corr(), pdf.corr())

    assert_eq(df.x.cov(df.y), pdf.x.cov(pdf.y))
    assert_eq(df.x.corr(df.y), pdf.x.corr(pdf.y))
    assert_eq(df.x.autocorr(lag=1), pdf.x.autocorr(lag=1))


def test_reduction_on_empty_df():
    pdf = pd.DataFrame()
    df = from_pandas(pdf)
    assert_eq(df.sum(), pdf.sum())


@pytest.mark.parametrize("axis", [0, 1])
@pytest.mark.parametrize(
    "skipna",
    [
        True,
        pytest.param(
            False, marks=xfail_gpu("cudf requires skipna=True when nulls are present.")
        ),
    ],
)
@pytest.mark.parametrize("ddof", [1, 2])
def test_std_kwargs(axis, skipna, ddof):
    pdf = pd.DataFrame(
        {"x": range(30), "y": [1, 2, None] * 10, "z": ["dog", "cat"] * 15}
    )
    df = from_pandas(pdf, npartitions=3)
    assert_eq(
        pdf.std(axis=axis, skipna=skipna, ddof=ddof, numeric_only=True),
        df.std(axis=axis, skipna=skipna, ddof=ddof, numeric_only=True),
    )


def test_mean_series_axis_none(df, pdf):
    assert_eq(df.x.mean(axis=None), pdf.x.mean(axis=None))


def test_mode_numeric_only():
    df = pd.DataFrame(
        {
            "int": [1, 2, 3, 4, 5, 6, 7, 8],
            "float": [1.0, 2.0, 3.0, 4.0, np.nan, 6.0, 7.0, 8.0],
            "dt": [pd.NaT] + [datetime(2010, i, 1) for i in range(1, 8)],
            "timedelta": pd.to_timedelta([1, 2, 3, 4, 5, 6, 7, np.nan]),
        }
    )
    ddf = from_pandas(df, npartitions=2)

    assert_eq(ddf.mode(numeric_only=False), df.mode(numeric_only=False))
    assert_eq(ddf.mode(), df.mode())
    assert_eq(ddf.mode(numeric_only=True), df.mode(numeric_only=True))


def test_divmod():
    df1 = pd.Series(np.random.rand(10))
    df2 = pd.Series(np.random.rand(10))

    ddf1 = from_pandas(df1, npartitions=3)
    ddf2 = from_pandas(df2, npartitions=3)

    result = divmod(ddf1, 2.0)
    expected = divmod(df1, 2.0)
    assert_eq(result[0], expected[0])
    assert_eq(result[1], expected[1])

    result = divmod(ddf1, ddf2)
    expected = divmod(df1, df2)
    assert_eq(result[0], expected[0])
    assert_eq(result[1], expected[1])


def test_value_counts_with_normalize():
    df = pd.DataFrame({"x": [1, 2, 1, 3, 3, 1, 4]})
    ddf = from_pandas(df, npartitions=3)
    result = ddf.x.value_counts(normalize=True)
    expected = df.x.value_counts(normalize=True)
    assert_eq(result, expected)

    result2 = ddf.x.value_counts(split_every=2, normalize=True)
    assert_eq(result2, expected)
    assert result._name != result2._name

    result3 = ddf.x.value_counts(split_out=2, normalize=True)
    assert_eq(result3, expected)
    assert result._name != result3._name


def test_reduction_method():
    df = pd.DataFrame({"x": range(50), "y": range(50, 100)})
    ddf = from_pandas(df, npartitions=4)

    chunk = lambda x, val=0: (x >= val).sum()
    agg = lambda x: x.sum()

    # Output of chunk is a scalar
    res = ddf.x.reduction(chunk, aggregate=agg)
    assert_eq(res, df.x.count())

    # Output of chunk is a series
    res = ddf.reduction(chunk, aggregate=agg)
    assert res._name == ddf.reduction(chunk, aggregate=agg)._name
    assert_eq(res, df.count())

    # Test with keywords
    res2 = ddf.reduction(chunk, aggregate=agg, chunk_kwargs={"val": 25})
    assert (
        res2._name
        == ddf.reduction(chunk, aggregate=agg, chunk_kwargs={"val": 25})._name
    )
    assert res2._name != res._name
    assert_eq(res2, (df >= 25).sum())

    # Output of chunk is a dataframe
    def sum_and_count(x):
        return pd.DataFrame({"sum": x.sum(), "count": x.count()})

    res = ddf.reduction(sum_and_count, aggregate=lambda x: x.groupby(level=0).sum())

    assert_eq(res, pd.DataFrame({"sum": df.sum(), "count": df.count()}))


def test_reduction_method_split_every():
    df = pd.Series([1] * 60)
    ddf = from_pandas(df, npartitions=15)

    def chunk(x, constant=0):
        return x.sum() + constant

    def combine(x, constant=0):
        return x.sum() + constant + 1

    def agg(x, constant=0):
        return x.sum() + constant + 2

    f = lambda n: ddf.reduction(
        chunk,
        aggregate=agg,
        combine=combine,
        chunk_kwargs=dict(constant=1.0),
        combine_kwargs=dict(constant=2.0),
        aggregate_kwargs=dict(constant=3.0),
        split_every=n,
    )

    # Keywords are different for each step
    assert f(3).compute() == 60 + 15 + 7 * (2 + 1) + (3 + 2)
    # Keywords are same for each step
    res = ddf.reduction(
        chunk, aggregate=agg, combine=combine, constant=3.0, split_every=3
    )
    assert res.compute() == 60 + 15 * 3 + 7 * (3 + 1) + (3 + 2)
    # No combine provided, combine is agg
    res = ddf.reduction(chunk, aggregate=agg, constant=3.0, split_every=3)
    assert res.compute() == 60 + 15 * 3 + 8 * (3 + 2)

    # split_every must be >= 2
    with pytest.raises(ValueError):
        f(1)

    # combine_kwargs with no combine provided
    with pytest.raises(ValueError):
        ddf.reduction(
            chunk,
            aggregate=agg,
            split_every=3,
            chunk_kwargs=dict(constant=1.0),
            combine_kwargs=dict(constant=2.0),
            aggregate_kwargs=dict(constant=3.0),
        )


def test_reduction_split_every_false():
    pdf = pd.DataFrame({"a": [1]})
    df = from_pandas(pdf, npartitions=1)
    result = df.reduction(chunk=lambda x: x, split_every=False)
    assert_eq(result, pdf)


@pytest.mark.parametrize("key", [0, 1, 2])
def test_unique_numerical_columns(key):
    pdf = pd.DataFrame({key: [0, 3, 0, 1, 2, 0, 3, 1, 1]})
    df = from_pandas(pdf, 3)
    assert_eq(
        df[key].unique(), pd.Series(pdf[key].unique(), name=key), check_index=False
    )


def test_cat_value_counts_large_unknown_categories():
    pdf = pd.DataFrame({"x": np.random.randint(1, 1_000_000, (250_000,))})
    df = from_pandas(pdf, npartitions=50)
    df["x"] = df["x"].astype("category")
    result = df.x.value_counts()
    assert result.npartitions == 50  # unknown
    pdf["x"] = pdf["x"].astype("category")
    expected = pdf.x.value_counts()
    assert_eq(result, expected, check_index=False, check_dtype=False)

    df = from_pandas(pdf, npartitions=50)
    result = df.x.value_counts()
    assert result.npartitions == 3  # known but large
    assert_eq(result, expected, check_index=False, check_dtype=False)


def test_reductions_timestamps_display():
    data = pd.to_datetime(["2024-10-02 12:00:00", "2024-10-02 14:00:00"])
    df = from_pandas(pd.DataFrame({"valid_time": data}))
    assert df["valid_time"].min().__repr__()


def test_value_counts_shuffle_properly():
    pdf = pd.DataFrame(
        {
            "A": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            "B": [11, 12, 13, 4, 15, 6, 17, 8, 19, 10],
        }
    )
    df = from_pandas(pdf, npartitions=2)
    result = (df["A"] == df["B"]).value_counts()
    expected = (pdf["A"] == pdf["B"]).value_counts()
    assert_eq(result, expected)
