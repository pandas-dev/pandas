import warnings

import numpy as np
from numpy.random import randn
import pytest

from pandas import DataFrame, Index, MultiIndex, Series, isna, notna
import pandas._testing as tm
from pandas.tests.window.common import (
    moments_consistency_cov_data,
    moments_consistency_is_constant,
    moments_consistency_mock_mean,
    moments_consistency_series_data,
    moments_consistency_std_data,
    moments_consistency_var_data,
    moments_consistency_var_debiasing_factors,
)


def _check_expanding(
    func, static_comp, preserve_nan=True, series=None, frame=None, nan_locs=None
):

    series_result = func(series)
    assert isinstance(series_result, Series)
    frame_result = func(frame)
    assert isinstance(frame_result, DataFrame)

    result = func(series)
    tm.assert_almost_equal(result[10], static_comp(series[:11]))

    if preserve_nan:
        assert result.iloc[nan_locs].isna().all()


def _check_expanding_has_min_periods(func, static_comp, has_min_periods):
    ser = Series(randn(50))

    if has_min_periods:
        result = func(ser, min_periods=30)
        assert result[:29].isna().all()
        tm.assert_almost_equal(result.iloc[-1], static_comp(ser[:50]))

        # min_periods is working correctly
        result = func(ser, min_periods=15)
        assert isna(result.iloc[13])
        assert notna(result.iloc[14])

        ser2 = Series(randn(20))
        result = func(ser2, min_periods=5)
        assert isna(result[3])
        assert notna(result[4])

        # min_periods=0
        result0 = func(ser, min_periods=0)
        result1 = func(ser, min_periods=1)
        tm.assert_almost_equal(result0, result1)
    else:
        result = func(ser)
        tm.assert_almost_equal(result.iloc[-1], static_comp(ser[:50]))


def test_expanding_corr(series):
    A = series.dropna()
    B = (A + randn(len(A)))[:-5]

    result = A.expanding().corr(B)

    rolling_result = A.rolling(window=len(A), min_periods=1).corr(B)

    tm.assert_almost_equal(rolling_result, result)


def test_expanding_count(series):
    result = series.expanding(min_periods=0).count()
    tm.assert_almost_equal(
        result, series.rolling(window=len(series), min_periods=0).count()
    )


def test_expanding_quantile(series):
    result = series.expanding().quantile(0.5)

    rolling_result = series.rolling(window=len(series), min_periods=1).quantile(0.5)

    tm.assert_almost_equal(result, rolling_result)


def test_expanding_cov(series):
    A = series
    B = (A + randn(len(A)))[:-5]

    result = A.expanding().cov(B)

    rolling_result = A.rolling(window=len(A), min_periods=1).cov(B)

    tm.assert_almost_equal(rolling_result, result)


def test_expanding_cov_pairwise(frame):
    result = frame.expanding().cov()

    rolling_result = frame.rolling(window=len(frame), min_periods=1).cov()

    tm.assert_frame_equal(result, rolling_result)


def test_expanding_corr_pairwise(frame):
    result = frame.expanding().corr()

    rolling_result = frame.rolling(window=len(frame), min_periods=1).corr()
    tm.assert_frame_equal(result, rolling_result)


@pytest.mark.parametrize("has_min_periods", [True, False])
@pytest.mark.parametrize(
    "func,static_comp",
    [("sum", np.sum), ("mean", np.mean), ("max", np.max), ("min", np.min)],
    ids=["sum", "mean", "max", "min"],
)
def test_expanding_func(func, static_comp, has_min_periods, series, frame, nan_locs):
    def expanding_func(x, min_periods=1, center=False, axis=0):
        exp = x.expanding(min_periods=min_periods, center=center, axis=axis)
        return getattr(exp, func)()

    _check_expanding(
        expanding_func,
        static_comp,
        preserve_nan=False,
        series=series,
        frame=frame,
        nan_locs=nan_locs,
    )
    _check_expanding_has_min_periods(expanding_func, static_comp, has_min_periods)


@pytest.mark.parametrize("has_min_periods", [True, False])
def test_expanding_apply(engine_and_raw, has_min_periods, series, frame, nan_locs):

    engine, raw = engine_and_raw

    def expanding_mean(x, min_periods=1):

        exp = x.expanding(min_periods=min_periods)
        result = exp.apply(lambda x: x.mean(), raw=raw, engine=engine)
        return result

    # TODO(jreback), needed to add preserve_nan=False
    # here to make this pass
    _check_expanding(
        expanding_mean,
        np.mean,
        preserve_nan=False,
        series=series,
        frame=frame,
        nan_locs=nan_locs,
    )
    _check_expanding_has_min_periods(expanding_mean, np.mean, has_min_periods)


@pytest.mark.parametrize("min_periods", [0, 1, 2, 3, 4])
def test_expanding_apply_consistency(
    consistency_data, base_functions, no_nan_functions, min_periods
):
    x, is_constant, no_nans = consistency_data

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message=".*(empty slice|0 for slice).*", category=RuntimeWarning,
        )
        # test consistency between expanding_xyz() and either (a)
        # expanding_apply of Series.xyz(), or (b) expanding_apply of
        # np.nanxyz()
        functions = base_functions

        # GH 8269
        if no_nans:
            functions = base_functions + no_nan_functions
        for (f, require_min_periods, name) in functions:
            expanding_f = getattr(x.expanding(min_periods=min_periods), name)

            if (
                require_min_periods
                and (min_periods is not None)
                and (min_periods < require_min_periods)
            ):
                continue

            if name == "count":
                expanding_f_result = expanding_f()
                expanding_apply_f_result = x.expanding(min_periods=0).apply(
                    func=f, raw=True
                )
            else:
                if name in ["cov", "corr"]:
                    expanding_f_result = expanding_f(pairwise=False)
                else:
                    expanding_f_result = expanding_f()
                expanding_apply_f_result = x.expanding(min_periods=min_periods).apply(
                    func=f, raw=True
                )

            # GH 9422
            if name in ["sum", "prod"]:
                tm.assert_equal(expanding_f_result, expanding_apply_f_result)


@pytest.mark.parametrize("min_periods", [0, 1, 2, 3, 4])
def test_moments_consistency_var(consistency_data, min_periods):
    x, is_constant, no_nans = consistency_data
    moments_consistency_var_data(
        x=x,
        is_constant=is_constant,
        min_periods=min_periods,
        count=lambda x: x.expanding(min_periods=min_periods).count(),
        mean=lambda x: x.expanding(min_periods=min_periods).mean(),
        var_unbiased=lambda x: x.expanding(min_periods=min_periods).var(),
        var_biased=lambda x: x.expanding(min_periods=min_periods).var(ddof=0),
    )


@pytest.mark.parametrize("min_periods", [0, 1, 2, 3, 4])
def test_expanding_consistency_std(consistency_data, min_periods):
    x, is_constant, no_nans = consistency_data
    moments_consistency_std_data(
        x=x,
        var_unbiased=lambda x: x.expanding(min_periods=min_periods).var(),
        std_unbiased=lambda x: x.expanding(min_periods=min_periods).std(),
        var_biased=lambda x: x.expanding(min_periods=min_periods).var(ddof=0),
        std_biased=lambda x: x.expanding(min_periods=min_periods).std(ddof=0),
    )


@pytest.mark.parametrize("min_periods", [0, 1, 2, 3, 4])
def test_expanding_consistency_cov(consistency_data, min_periods):
    x, is_constant, no_nans = consistency_data
    moments_consistency_cov_data(
        x=x,
        var_unbiased=lambda x: x.expanding(min_periods=min_periods).var(),
        cov_unbiased=lambda x, y: x.expanding(min_periods=min_periods).cov(y),
        var_biased=lambda x: x.expanding(min_periods=min_periods).var(ddof=0),
        cov_biased=lambda x, y: x.expanding(min_periods=min_periods).cov(y, ddof=0),
    )


@pytest.mark.parametrize("min_periods", [0, 1, 2, 3, 4])
def test_expanding_consistency_series(consistency_data, min_periods):
    x, is_constant, no_nans = consistency_data
    moments_consistency_series_data(
        x=x,
        mean=lambda x: x.expanding(min_periods=min_periods).mean(),
        corr=lambda x, y: x.expanding(min_periods=min_periods).corr(y),
        var_unbiased=lambda x: x.expanding(min_periods=min_periods).var(),
        std_unbiased=lambda x: x.expanding(min_periods=min_periods).std(),
        cov_unbiased=lambda x, y: x.expanding(min_periods=min_periods).cov(y),
        var_biased=lambda x: x.expanding(min_periods=min_periods).var(ddof=0),
        std_biased=lambda x: x.expanding(min_periods=min_periods).std(ddof=0),
        cov_biased=lambda x, y: x.expanding(min_periods=min_periods).cov(y, ddof=0),
    )


@pytest.mark.slow
@pytest.mark.parametrize("min_periods", [0, 1, 2, 3, 4])
def test_expanding_consistency(consistency_data, min_periods):
    x, is_constant, no_nans = consistency_data
    # suppress warnings about empty slices, as we are deliberately testing
    # with empty/0-length Series/DataFrames
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message=".*(empty slice|0 for slice).*", category=RuntimeWarning,
        )

        # test consistency between different expanding_* moments
        moments_consistency_mock_mean(
            x=x,
            mean=lambda x: x.expanding(min_periods=min_periods).mean(),
            mock_mean=lambda x: x.expanding(min_periods=min_periods).sum()
            / x.expanding().count(),
        )

        moments_consistency_is_constant(
            x=x,
            is_constant=is_constant,
            min_periods=min_periods,
            count=lambda x: x.expanding().count(),
            mean=lambda x: x.expanding(min_periods=min_periods).mean(),
            corr=lambda x, y: x.expanding(min_periods=min_periods).corr(y),
        )

        moments_consistency_var_debiasing_factors(
            x=x,
            var_unbiased=lambda x: x.expanding(min_periods=min_periods).var(),
            var_biased=lambda x: x.expanding(min_periods=min_periods).var(ddof=0),
            var_debiasing_factors=lambda x: (
                x.expanding().count()
                / (x.expanding().count() - 1.0).replace(0.0, np.nan)
            ),
        )


@pytest.mark.parametrize(
    "f",
    [
        lambda x: (x.expanding(min_periods=5).cov(x, pairwise=True)),
        lambda x: (x.expanding(min_periods=5).corr(x, pairwise=True)),
    ],
)
def test_moment_functions_zero_length_pairwise(f):

    df1 = DataFrame()
    df2 = DataFrame(columns=Index(["a"], name="foo"), index=Index([], name="bar"))
    df2["a"] = df2["a"].astype("float64")

    df1_expected = DataFrame(
        index=MultiIndex.from_product([df1.index, df1.columns]), columns=Index([])
    )
    df2_expected = DataFrame(
        index=MultiIndex.from_product([df2.index, df2.columns], names=["bar", "foo"]),
        columns=Index(["a"], name="foo"),
        dtype="float64",
    )

    df1_result = f(df1)
    tm.assert_frame_equal(df1_result, df1_expected)

    df2_result = f(df2)
    tm.assert_frame_equal(df2_result, df2_expected)


@pytest.mark.parametrize(
    "f",
    [
        lambda x: x.expanding().count(),
        lambda x: x.expanding(min_periods=5).cov(x, pairwise=False),
        lambda x: x.expanding(min_periods=5).corr(x, pairwise=False),
        lambda x: x.expanding(min_periods=5).max(),
        lambda x: x.expanding(min_periods=5).min(),
        lambda x: x.expanding(min_periods=5).sum(),
        lambda x: x.expanding(min_periods=5).mean(),
        lambda x: x.expanding(min_periods=5).std(),
        lambda x: x.expanding(min_periods=5).var(),
        lambda x: x.expanding(min_periods=5).skew(),
        lambda x: x.expanding(min_periods=5).kurt(),
        lambda x: x.expanding(min_periods=5).quantile(0.5),
        lambda x: x.expanding(min_periods=5).median(),
        lambda x: x.expanding(min_periods=5).apply(sum, raw=False),
        lambda x: x.expanding(min_periods=5).apply(sum, raw=True),
    ],
)
def test_moment_functions_zero_length(f):
    # GH 8056
    s = Series(dtype=np.float64)
    s_expected = s
    df1 = DataFrame()
    df1_expected = df1
    df2 = DataFrame(columns=["a"])
    df2["a"] = df2["a"].astype("float64")
    df2_expected = df2

    s_result = f(s)
    tm.assert_series_equal(s_result, s_expected)

    df1_result = f(df1)
    tm.assert_frame_equal(df1_result, df1_expected)

    df2_result = f(df2)
    tm.assert_frame_equal(df2_result, df2_expected)


def test_expanding_apply_empty_series(engine_and_raw):
    engine, raw = engine_and_raw
    ser = Series([], dtype=np.float64)
    tm.assert_series_equal(
        ser, ser.expanding().apply(lambda x: x.mean(), raw=raw, engine=engine)
    )


def test_expanding_apply_min_periods_0(engine_and_raw):
    # GH 8080
    engine, raw = engine_and_raw
    s = Series([None, None, None])
    result = s.expanding(min_periods=0).apply(lambda x: len(x), raw=raw, engine=engine)
    expected = Series([1.0, 2.0, 3.0])
    tm.assert_series_equal(result, expected)


def test_expanding_cov_diff_index():
    # GH 7512
    s1 = Series([1, 2, 3], index=[0, 1, 2])
    s2 = Series([1, 3], index=[0, 2])
    result = s1.expanding().cov(s2)
    expected = Series([None, None, 2.0])
    tm.assert_series_equal(result, expected)

    s2a = Series([1, None, 3], index=[0, 1, 2])
    result = s1.expanding().cov(s2a)
    tm.assert_series_equal(result, expected)

    s1 = Series([7, 8, 10], index=[0, 1, 3])
    s2 = Series([7, 9, 10], index=[0, 2, 3])
    result = s1.expanding().cov(s2)
    expected = Series([None, None, None, 4.5])
    tm.assert_series_equal(result, expected)


def test_expanding_corr_diff_index():
    # GH 7512
    s1 = Series([1, 2, 3], index=[0, 1, 2])
    s2 = Series([1, 3], index=[0, 2])
    result = s1.expanding().corr(s2)
    expected = Series([None, None, 1.0])
    tm.assert_series_equal(result, expected)

    s2a = Series([1, None, 3], index=[0, 1, 2])
    result = s1.expanding().corr(s2a)
    tm.assert_series_equal(result, expected)

    s1 = Series([7, 8, 10], index=[0, 1, 3])
    s2 = Series([7, 9, 10], index=[0, 2, 3])
    result = s1.expanding().corr(s2)
    expected = Series([None, None, None, 1.0])
    tm.assert_series_equal(result, expected)


def test_expanding_cov_pairwise_diff_length():
    # GH 7512
    df1 = DataFrame([[1, 5], [3, 2], [3, 9]], columns=Index(["A", "B"], name="foo"))
    df1a = DataFrame(
        [[1, 5], [3, 9]], index=[0, 2], columns=Index(["A", "B"], name="foo")
    )
    df2 = DataFrame(
        [[5, 6], [None, None], [2, 1]], columns=Index(["X", "Y"], name="foo")
    )
    df2a = DataFrame(
        [[5, 6], [2, 1]], index=[0, 2], columns=Index(["X", "Y"], name="foo")
    )
    # TODO: xref gh-15826
    # .loc is not preserving the names
    result1 = df1.expanding().cov(df2, pairwise=True).loc[2]
    result2 = df1.expanding().cov(df2a, pairwise=True).loc[2]
    result3 = df1a.expanding().cov(df2, pairwise=True).loc[2]
    result4 = df1a.expanding().cov(df2a, pairwise=True).loc[2]
    expected = DataFrame(
        [[-3.0, -6.0], [-5.0, -10.0]],
        columns=Index(["A", "B"], name="foo"),
        index=Index(["X", "Y"], name="foo"),
    )
    tm.assert_frame_equal(result1, expected)
    tm.assert_frame_equal(result2, expected)
    tm.assert_frame_equal(result3, expected)
    tm.assert_frame_equal(result4, expected)


def test_expanding_corr_pairwise_diff_length():
    # GH 7512
    df1 = DataFrame(
        [[1, 2], [3, 2], [3, 4]], columns=["A", "B"], index=Index(range(3), name="bar"),
    )
    df1a = DataFrame(
        [[1, 2], [3, 4]], index=Index([0, 2], name="bar"), columns=["A", "B"]
    )
    df2 = DataFrame(
        [[5, 6], [None, None], [2, 1]],
        columns=["X", "Y"],
        index=Index(range(3), name="bar"),
    )
    df2a = DataFrame(
        [[5, 6], [2, 1]], index=Index([0, 2], name="bar"), columns=["X", "Y"]
    )
    result1 = df1.expanding().corr(df2, pairwise=True).loc[2]
    result2 = df1.expanding().corr(df2a, pairwise=True).loc[2]
    result3 = df1a.expanding().corr(df2, pairwise=True).loc[2]
    result4 = df1a.expanding().corr(df2a, pairwise=True).loc[2]
    expected = DataFrame(
        [[-1.0, -1.0], [-1.0, -1.0]], columns=["A", "B"], index=Index(["X", "Y"])
    )
    tm.assert_frame_equal(result1, expected)
    tm.assert_frame_equal(result2, expected)
    tm.assert_frame_equal(result3, expected)
    tm.assert_frame_equal(result4, expected)


def test_expanding_apply_args_kwargs(engine_and_raw):
    def mean_w_arg(x, const):
        return np.mean(x) + const

    engine, raw = engine_and_raw

    df = DataFrame(np.random.rand(20, 3))

    expected = df.expanding().apply(np.mean, engine=engine, raw=raw) + 20.0

    result = df.expanding().apply(mean_w_arg, engine=engine, raw=raw, args=(20,))
    tm.assert_frame_equal(result, expected)

    result = df.expanding().apply(mean_w_arg, raw=raw, kwargs={"const": 20})
    tm.assert_frame_equal(result, expected)
