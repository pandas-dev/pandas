from datetime import datetime
import warnings

import numpy as np
from numpy.random import randn
import pytest

import pandas.util._test_decorators as td

import pandas as pd
from pandas import DataFrame, DatetimeIndex, Index, Series
import pandas._testing as tm
from pandas.core.window.common import _flex_binary_moment
from pandas.tests.window.common import (
    check_pairwise_moment,
    moments_consistency_cov_data,
    moments_consistency_is_constant,
    moments_consistency_mock_mean,
    moments_consistency_series_data,
    moments_consistency_std_data,
    moments_consistency_var_data,
    moments_consistency_var_debiasing_factors,
)


def _rolling_consistency_cases():
    for window in [1, 2, 3, 10, 20]:
        for min_periods in {0, 1, 2, 3, 4, window}:
            if min_periods and (min_periods > window):
                continue
            for center in [False, True]:
                yield window, min_periods, center


# binary moments
def test_rolling_cov(series):
    A = series
    B = A + randn(len(A))

    result = A.rolling(window=50, min_periods=25).cov(B)
    tm.assert_almost_equal(result[-1], np.cov(A[-50:], B[-50:])[0, 1])


def test_rolling_corr(series):
    A = series
    B = A + randn(len(A))

    result = A.rolling(window=50, min_periods=25).corr(B)
    tm.assert_almost_equal(result[-1], np.corrcoef(A[-50:], B[-50:])[0, 1])

    # test for correct bias correction
    a = tm.makeTimeSeries()
    b = tm.makeTimeSeries()
    a[:5] = np.nan
    b[:10] = np.nan

    result = a.rolling(window=len(a), min_periods=1).corr(b)
    tm.assert_almost_equal(result[-1], a.corr(b))


@pytest.mark.parametrize("func", ["cov", "corr"])
def test_rolling_pairwise_cov_corr(func, frame):
    check_pairwise_moment(frame, "rolling", func, window=10, min_periods=5)


@pytest.mark.parametrize("method", ["corr", "cov"])
def test_flex_binary_frame(method, frame):
    series = frame[1]

    res = getattr(series.rolling(window=10), method)(frame)
    res2 = getattr(frame.rolling(window=10), method)(series)
    exp = frame.apply(lambda x: getattr(series.rolling(window=10), method)(x))

    tm.assert_frame_equal(res, exp)
    tm.assert_frame_equal(res2, exp)

    frame2 = frame.copy()
    frame2.values[:] = np.random.randn(*frame2.shape)

    res3 = getattr(frame.rolling(window=10), method)(frame2)
    exp = DataFrame(
        {k: getattr(frame[k].rolling(window=10), method)(frame2[k]) for k in frame}
    )
    tm.assert_frame_equal(res3, exp)


@pytest.mark.slow
@pytest.mark.parametrize(
    "window,min_periods,center", list(_rolling_consistency_cases())
)
def test_rolling_apply_consistency(
    consistency_data, base_functions, no_nan_functions, window, min_periods, center
):
    x, is_constant, no_nans = consistency_data

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message=".*(empty slice|0 for slice).*", category=RuntimeWarning,
        )
        # test consistency between rolling_xyz() and either (a)
        # rolling_apply of Series.xyz(), or (b) rolling_apply of
        # np.nanxyz()
        functions = base_functions

        # GH 8269
        if no_nans:
            functions = no_nan_functions + base_functions
        for (f, require_min_periods, name) in functions:
            rolling_f = getattr(
                x.rolling(window=window, center=center, min_periods=min_periods), name,
            )

            if (
                require_min_periods
                and (min_periods is not None)
                and (min_periods < require_min_periods)
            ):
                continue

            if name == "count":
                rolling_f_result = rolling_f()
                rolling_apply_f_result = x.rolling(
                    window=window, min_periods=min_periods, center=center
                ).apply(func=f, raw=True)
            else:
                if name in ["cov", "corr"]:
                    rolling_f_result = rolling_f(pairwise=False)
                else:
                    rolling_f_result = rolling_f()
                rolling_apply_f_result = x.rolling(
                    window=window, min_periods=min_periods, center=center
                ).apply(func=f, raw=True)

            # GH 9422
            if name in ["sum", "prod"]:
                tm.assert_equal(rolling_f_result, rolling_apply_f_result)


@pytest.mark.parametrize("window", range(7))
def test_rolling_corr_with_zero_variance(window):
    # GH 18430
    s = pd.Series(np.zeros(20))
    other = pd.Series(np.arange(20))

    assert s.rolling(window=window).corr(other=other).isna().all()


def test_flex_binary_moment():
    # GH3155
    # don't blow the stack
    msg = "arguments to moment function must be of type np.ndarray/Series/DataFrame"
    with pytest.raises(TypeError, match=msg):
        _flex_binary_moment(5, 6, None)


def test_corr_sanity():
    # GH 3155
    df = DataFrame(
        np.array(
            [
                [0.87024726, 0.18505595],
                [0.64355431, 0.3091617],
                [0.92372966, 0.50552513],
                [0.00203756, 0.04520709],
                [0.84780328, 0.33394331],
                [0.78369152, 0.63919667],
            ]
        )
    )

    res = df[0].rolling(5, center=True).corr(df[1])
    assert all(np.abs(np.nan_to_num(x)) <= 1 for x in res)

    # and some fuzzing
    for _ in range(10):
        df = DataFrame(np.random.rand(30, 2))
        res = df[0].rolling(5, center=True).corr(df[1])
        try:
            assert all(np.abs(np.nan_to_num(x)) <= 1 for x in res)
        except AssertionError:
            print(res)


def test_rolling_cov_diff_length():
    # GH 7512
    s1 = Series([1, 2, 3], index=[0, 1, 2])
    s2 = Series([1, 3], index=[0, 2])
    result = s1.rolling(window=3, min_periods=2).cov(s2)
    expected = Series([None, None, 2.0])
    tm.assert_series_equal(result, expected)

    s2a = Series([1, None, 3], index=[0, 1, 2])
    result = s1.rolling(window=3, min_periods=2).cov(s2a)
    tm.assert_series_equal(result, expected)


def test_rolling_corr_diff_length():
    # GH 7512
    s1 = Series([1, 2, 3], index=[0, 1, 2])
    s2 = Series([1, 3], index=[0, 2])
    result = s1.rolling(window=3, min_periods=2).corr(s2)
    expected = Series([None, None, 1.0])
    tm.assert_series_equal(result, expected)

    s2a = Series([1, None, 3], index=[0, 1, 2])
    result = s1.rolling(window=3, min_periods=2).corr(s2a)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "f",
    [
        lambda x: x.rolling(window=10, min_periods=5).cov(x, pairwise=False),
        lambda x: x.rolling(window=10, min_periods=5).corr(x, pairwise=False),
        lambda x: x.rolling(window=10, min_periods=5).max(),
        lambda x: x.rolling(window=10, min_periods=5).min(),
        lambda x: x.rolling(window=10, min_periods=5).sum(),
        lambda x: x.rolling(window=10, min_periods=5).mean(),
        lambda x: x.rolling(window=10, min_periods=5).std(),
        lambda x: x.rolling(window=10, min_periods=5).var(),
        lambda x: x.rolling(window=10, min_periods=5).skew(),
        lambda x: x.rolling(window=10, min_periods=5).kurt(),
        lambda x: x.rolling(window=10, min_periods=5).quantile(quantile=0.5),
        lambda x: x.rolling(window=10, min_periods=5).median(),
        lambda x: x.rolling(window=10, min_periods=5).apply(sum, raw=False),
        lambda x: x.rolling(window=10, min_periods=5).apply(sum, raw=True),
        lambda x: x.rolling(win_type="boxcar", window=10, min_periods=5).mean(),
    ],
)
@td.skip_if_no_scipy
def test_rolling_functions_window_non_shrinkage(f):
    # GH 7764
    s = Series(range(4))
    s_expected = Series(np.nan, index=s.index)
    df = DataFrame([[1, 5], [3, 2], [3, 9], [-1, 0]], columns=["A", "B"])
    df_expected = DataFrame(np.nan, index=df.index, columns=df.columns)

    s_result = f(s)
    tm.assert_series_equal(s_result, s_expected)

    df_result = f(df)
    tm.assert_frame_equal(df_result, df_expected)


def test_rolling_functions_window_non_shrinkage_binary():

    # corr/cov return a MI DataFrame
    df = DataFrame(
        [[1, 5], [3, 2], [3, 9], [-1, 0]],
        columns=Index(["A", "B"], name="foo"),
        index=Index(range(4), name="bar"),
    )
    df_expected = DataFrame(
        columns=Index(["A", "B"], name="foo"),
        index=pd.MultiIndex.from_product([df.index, df.columns], names=["bar", "foo"]),
        dtype="float64",
    )
    functions = [
        lambda x: (x.rolling(window=10, min_periods=5).cov(x, pairwise=True)),
        lambda x: (x.rolling(window=10, min_periods=5).corr(x, pairwise=True)),
    ]
    for f in functions:
        df_result = f(df)
        tm.assert_frame_equal(df_result, df_expected)


def test_rolling_skew_edge_cases():

    all_nan = Series([np.NaN] * 5)

    # yields all NaN (0 variance)
    d = Series([1] * 5)
    x = d.rolling(window=5).skew()
    tm.assert_series_equal(all_nan, x)

    # yields all NaN (window too small)
    d = Series(np.random.randn(5))
    x = d.rolling(window=2).skew()
    tm.assert_series_equal(all_nan, x)

    # yields [NaN, NaN, NaN, 0.177994, 1.548824]
    d = Series([-1.50837035, -0.1297039, 0.19501095, 1.73508164, 0.41941401])
    expected = Series([np.NaN, np.NaN, np.NaN, 0.177994, 1.548824])
    x = d.rolling(window=4).skew()
    tm.assert_series_equal(expected, x)


def test_rolling_kurt_edge_cases():

    all_nan = Series([np.NaN] * 5)

    # yields all NaN (0 variance)
    d = Series([1] * 5)
    x = d.rolling(window=5).kurt()
    tm.assert_series_equal(all_nan, x)

    # yields all NaN (window too small)
    d = Series(np.random.randn(5))
    x = d.rolling(window=3).kurt()
    tm.assert_series_equal(all_nan, x)

    # yields [NaN, NaN, NaN, 1.224307, 2.671499]
    d = Series([-1.50837035, -0.1297039, 0.19501095, 1.73508164, 0.41941401])
    expected = Series([np.NaN, np.NaN, np.NaN, 1.224307, 2.671499])
    x = d.rolling(window=4).kurt()
    tm.assert_series_equal(expected, x)


def test_rolling_skew_eq_value_fperr():
    # #18804 all rolling skew for all equal values should return Nan
    a = Series([1.1] * 15).rolling(window=10).skew()
    assert np.isnan(a).all()


def test_rolling_kurt_eq_value_fperr():
    # #18804 all rolling kurt for all equal values should return Nan
    a = Series([1.1] * 15).rolling(window=10).kurt()
    assert np.isnan(a).all()


def test_rolling_max_gh6297():
    """Replicate result expected in GH #6297"""
    indices = [datetime(1975, 1, i) for i in range(1, 6)]
    # So that we can have 2 datapoints on one of the days
    indices.append(datetime(1975, 1, 3, 6, 0))
    series = Series(range(1, 7), index=indices)
    # Use floats instead of ints as values
    series = series.map(lambda x: float(x))
    # Sort chronologically
    series = series.sort_index()

    expected = Series(
        [1.0, 2.0, 6.0, 4.0, 5.0],
        index=DatetimeIndex([datetime(1975, 1, i, 0) for i in range(1, 6)], freq="D"),
    )
    x = series.resample("D").max().rolling(window=1).max()
    tm.assert_series_equal(expected, x)


def test_rolling_max_resample():

    indices = [datetime(1975, 1, i) for i in range(1, 6)]
    # So that we can have 3 datapoints on last day (4, 10, and 20)
    indices.append(datetime(1975, 1, 5, 1))
    indices.append(datetime(1975, 1, 5, 2))
    series = Series(list(range(0, 5)) + [10, 20], index=indices)
    # Use floats instead of ints as values
    series = series.map(lambda x: float(x))
    # Sort chronologically
    series = series.sort_index()

    # Default how should be max
    expected = Series(
        [0.0, 1.0, 2.0, 3.0, 20.0],
        index=DatetimeIndex([datetime(1975, 1, i, 0) for i in range(1, 6)], freq="D"),
    )
    x = series.resample("D").max().rolling(window=1).max()
    tm.assert_series_equal(expected, x)

    # Now specify median (10.0)
    expected = Series(
        [0.0, 1.0, 2.0, 3.0, 10.0],
        index=DatetimeIndex([datetime(1975, 1, i, 0) for i in range(1, 6)], freq="D"),
    )
    x = series.resample("D").median().rolling(window=1).max()
    tm.assert_series_equal(expected, x)

    # Now specify mean (4+10+20)/3
    v = (4.0 + 10.0 + 20.0) / 3.0
    expected = Series(
        [0.0, 1.0, 2.0, 3.0, v],
        index=DatetimeIndex([datetime(1975, 1, i, 0) for i in range(1, 6)], freq="D"),
    )
    x = series.resample("D").mean().rolling(window=1).max()
    tm.assert_series_equal(expected, x)


def test_rolling_min_resample():

    indices = [datetime(1975, 1, i) for i in range(1, 6)]
    # So that we can have 3 datapoints on last day (4, 10, and 20)
    indices.append(datetime(1975, 1, 5, 1))
    indices.append(datetime(1975, 1, 5, 2))
    series = Series(list(range(0, 5)) + [10, 20], index=indices)
    # Use floats instead of ints as values
    series = series.map(lambda x: float(x))
    # Sort chronologically
    series = series.sort_index()

    # Default how should be min
    expected = Series(
        [0.0, 1.0, 2.0, 3.0, 4.0],
        index=DatetimeIndex([datetime(1975, 1, i, 0) for i in range(1, 6)], freq="D"),
    )
    r = series.resample("D").min().rolling(window=1)
    tm.assert_series_equal(expected, r.min())


def test_rolling_median_resample():

    indices = [datetime(1975, 1, i) for i in range(1, 6)]
    # So that we can have 3 datapoints on last day (4, 10, and 20)
    indices.append(datetime(1975, 1, 5, 1))
    indices.append(datetime(1975, 1, 5, 2))
    series = Series(list(range(0, 5)) + [10, 20], index=indices)
    # Use floats instead of ints as values
    series = series.map(lambda x: float(x))
    # Sort chronologically
    series = series.sort_index()

    # Default how should be median
    expected = Series(
        [0.0, 1.0, 2.0, 3.0, 10],
        index=DatetimeIndex([datetime(1975, 1, i, 0) for i in range(1, 6)], freq="D"),
    )
    x = series.resample("D").median().rolling(window=1).median()
    tm.assert_series_equal(expected, x)


def test_rolling_median_memory_error():
    # GH11722
    n = 20000
    Series(np.random.randn(n)).rolling(window=2, center=False).median()
    Series(np.random.randn(n)).rolling(window=2, center=False).median()


def test_rolling_min_max_numeric_types():

    # GH12373
    types_test = [np.dtype(f"f{width}") for width in [4, 8]]
    types_test.extend(
        [np.dtype(f"{sign}{width}") for width in [1, 2, 4, 8] for sign in "ui"]
    )
    for data_type in types_test:
        # Just testing that these don't throw exceptions and that
        # the return type is float64. Other tests will cover quantitative
        # correctness
        result = DataFrame(np.arange(20, dtype=data_type)).rolling(window=5).max()
        assert result.dtypes[0] == np.dtype("f8")
        result = DataFrame(np.arange(20, dtype=data_type)).rolling(window=5).min()
        assert result.dtypes[0] == np.dtype("f8")


def test_moment_functions_zero_length():
    # GH 8056
    s = Series(dtype=np.float64)
    s_expected = s
    df1 = DataFrame()
    df1_expected = df1
    df2 = DataFrame(columns=["a"])
    df2["a"] = df2["a"].astype("float64")
    df2_expected = df2

    functions = [
        lambda x: x.rolling(window=10).count(),
        lambda x: x.rolling(window=10, min_periods=5).cov(x, pairwise=False),
        lambda x: x.rolling(window=10, min_periods=5).corr(x, pairwise=False),
        lambda x: x.rolling(window=10, min_periods=5).max(),
        lambda x: x.rolling(window=10, min_periods=5).min(),
        lambda x: x.rolling(window=10, min_periods=5).sum(),
        lambda x: x.rolling(window=10, min_periods=5).mean(),
        lambda x: x.rolling(window=10, min_periods=5).std(),
        lambda x: x.rolling(window=10, min_periods=5).var(),
        lambda x: x.rolling(window=10, min_periods=5).skew(),
        lambda x: x.rolling(window=10, min_periods=5).kurt(),
        lambda x: x.rolling(window=10, min_periods=5).quantile(0.5),
        lambda x: x.rolling(window=10, min_periods=5).median(),
        lambda x: x.rolling(window=10, min_periods=5).apply(sum, raw=False),
        lambda x: x.rolling(window=10, min_periods=5).apply(sum, raw=True),
        lambda x: x.rolling(win_type="boxcar", window=10, min_periods=5).mean(),
    ]
    for f in functions:
        try:
            s_result = f(s)
            tm.assert_series_equal(s_result, s_expected)

            df1_result = f(df1)
            tm.assert_frame_equal(df1_result, df1_expected)

            df2_result = f(df2)
            tm.assert_frame_equal(df2_result, df2_expected)
        except (ImportError):

            # scipy needed for rolling_window
            continue


def test_moment_functions_zero_length_pairwise():

    df1 = DataFrame()
    df2 = DataFrame(columns=Index(["a"], name="foo"), index=Index([], name="bar"))
    df2["a"] = df2["a"].astype("float64")

    df1_expected = DataFrame(
        index=pd.MultiIndex.from_product([df1.index, df1.columns]), columns=Index([]),
    )
    df2_expected = DataFrame(
        index=pd.MultiIndex.from_product(
            [df2.index, df2.columns], names=["bar", "foo"]
        ),
        columns=Index(["a"], name="foo"),
        dtype="float64",
    )

    functions = [
        lambda x: (x.rolling(window=10, min_periods=5).cov(x, pairwise=True)),
        lambda x: (x.rolling(window=10, min_periods=5).corr(x, pairwise=True)),
    ]

    for f in functions:
        df1_result = f(df1)
        tm.assert_frame_equal(df1_result, df1_expected)

        df2_result = f(df2)
        tm.assert_frame_equal(df2_result, df2_expected)


@pytest.mark.slow
@pytest.mark.parametrize(
    "window,min_periods,center", list(_rolling_consistency_cases())
)
def test_rolling_consistency_var(consistency_data, window, min_periods, center):
    x, is_constant, no_nans = consistency_data
    moments_consistency_var_data(
        x=x,
        is_constant=is_constant,
        min_periods=min_periods,
        count=lambda x: (
            x.rolling(window=window, min_periods=min_periods, center=center).count()
        ),
        mean=lambda x: (
            x.rolling(window=window, min_periods=min_periods, center=center).mean()
        ),
        var_unbiased=lambda x: (
            x.rolling(window=window, min_periods=min_periods, center=center).var()
        ),
        var_biased=lambda x: (
            x.rolling(window=window, min_periods=min_periods, center=center).var(ddof=0)
        ),
    )


@pytest.mark.slow
@pytest.mark.parametrize(
    "window,min_periods,center", list(_rolling_consistency_cases())
)
def test_rolling_consistency_std(consistency_data, window, min_periods, center):
    x, is_constant, no_nans = consistency_data
    moments_consistency_std_data(
        x=x,
        var_unbiased=lambda x: (
            x.rolling(window=window, min_periods=min_periods, center=center).var()
        ),
        std_unbiased=lambda x: (
            x.rolling(window=window, min_periods=min_periods, center=center).std()
        ),
        var_biased=lambda x: (
            x.rolling(window=window, min_periods=min_periods, center=center).var(ddof=0)
        ),
        std_biased=lambda x: (
            x.rolling(window=window, min_periods=min_periods, center=center).std(ddof=0)
        ),
    )


@pytest.mark.slow
@pytest.mark.parametrize(
    "window,min_periods,center", list(_rolling_consistency_cases())
)
def test_rolling_consistency_cov(consistency_data, window, min_periods, center):
    x, is_constant, no_nans = consistency_data
    moments_consistency_cov_data(
        x=x,
        var_unbiased=lambda x: (
            x.rolling(window=window, min_periods=min_periods, center=center).var()
        ),
        cov_unbiased=lambda x, y: (
            x.rolling(window=window, min_periods=min_periods, center=center).cov(y)
        ),
        var_biased=lambda x: (
            x.rolling(window=window, min_periods=min_periods, center=center).var(ddof=0)
        ),
        cov_biased=lambda x, y: (
            x.rolling(window=window, min_periods=min_periods, center=center).cov(
                y, ddof=0
            )
        ),
    )


@pytest.mark.slow
@pytest.mark.parametrize(
    "window,min_periods,center", list(_rolling_consistency_cases())
)
def test_rolling_consistency_series(consistency_data, window, min_periods, center):
    x, is_constant, no_nans = consistency_data
    moments_consistency_series_data(
        x=x,
        mean=lambda x: (
            x.rolling(window=window, min_periods=min_periods, center=center).mean()
        ),
        corr=lambda x, y: (
            x.rolling(window=window, min_periods=min_periods, center=center).corr(y)
        ),
        var_unbiased=lambda x: (
            x.rolling(window=window, min_periods=min_periods, center=center).var()
        ),
        std_unbiased=lambda x: (
            x.rolling(window=window, min_periods=min_periods, center=center).std()
        ),
        cov_unbiased=lambda x, y: (
            x.rolling(window=window, min_periods=min_periods, center=center).cov(y)
        ),
        var_biased=lambda x: (
            x.rolling(window=window, min_periods=min_periods, center=center).var(ddof=0)
        ),
        std_biased=lambda x: (
            x.rolling(window=window, min_periods=min_periods, center=center).std(ddof=0)
        ),
        cov_biased=lambda x, y: (
            x.rolling(window=window, min_periods=min_periods, center=center).cov(
                y, ddof=0
            )
        ),
    )


@pytest.mark.slow
@pytest.mark.parametrize(
    "window,min_periods,center", list(_rolling_consistency_cases())
)
def test_rolling_consistency(consistency_data, window, min_periods, center):
    x, is_constant, no_nans = consistency_data
    # suppress warnings about empty slices, as we are deliberately testing
    # with empty/0-length Series/DataFrames
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message=".*(empty slice|0 for slice).*", category=RuntimeWarning,
        )

        # test consistency between different rolling_* moments
        moments_consistency_mock_mean(
            x=x,
            mean=lambda x: (
                x.rolling(window=window, min_periods=min_periods, center=center).mean()
            ),
            mock_mean=lambda x: (
                x.rolling(window=window, min_periods=min_periods, center=center)
                .sum()
                .divide(
                    x.rolling(
                        window=window, min_periods=min_periods, center=center
                    ).count()
                )
            ),
        )

        moments_consistency_is_constant(
            x=x,
            is_constant=is_constant,
            min_periods=min_periods,
            count=lambda x: (
                x.rolling(window=window, min_periods=min_periods, center=center).count()
            ),
            mean=lambda x: (
                x.rolling(window=window, min_periods=min_periods, center=center).mean()
            ),
            corr=lambda x, y: (
                x.rolling(window=window, min_periods=min_periods, center=center).corr(y)
            ),
        )

        moments_consistency_var_debiasing_factors(
            x=x,
            var_unbiased=lambda x: (
                x.rolling(window=window, min_periods=min_periods, center=center).var()
            ),
            var_biased=lambda x: (
                x.rolling(window=window, min_periods=min_periods, center=center).var(
                    ddof=0
                )
            ),
            var_debiasing_factors=lambda x: (
                x.rolling(window=window, min_periods=min_periods, center=center)
                .count()
                .divide(
                    (
                        x.rolling(
                            window=window, min_periods=min_periods, center=center
                        ).count()
                        - 1.0
                    ).replace(0.0, np.nan)
                )
            ),
        )
