import numpy as np
import pytest

from pandas import Series
import pandas._testing as tm


def _rolling_consistency_cases():
    for window in [1, 2, 3, 10, 20]:
        for min_periods in {0, 1, 2, 3, 4, window}:
            if min_periods and (min_periods > window):
                continue
            for center in [False, True]:
                yield window, min_periods, center


@pytest.mark.parametrize(
    "window,min_periods,center", list(_rolling_consistency_cases())
)
@pytest.mark.parametrize("f", [lambda v: Series(v).sum(), np.nansum])
def test_rolling_apply_consistency_sum_nans(
    consistency_data, window, min_periods, center, f
):
    x, is_constant, no_nans = consistency_data

    if f is np.nansum and min_periods == 0:
        pass
    else:
        rolling_f_result = x.rolling(
            window=window, min_periods=min_periods, center=center
        ).sum()
        rolling_apply_f_result = x.rolling(
            window=window, min_periods=min_periods, center=center
        ).apply(func=f, raw=True)
        tm.assert_equal(rolling_f_result, rolling_apply_f_result)


@pytest.mark.parametrize(
    "window,min_periods,center", list(_rolling_consistency_cases())
)
@pytest.mark.parametrize("f", [lambda v: Series(v).sum(), np.nansum, np.sum])
def test_rolling_apply_consistency_sum_no_nans(
    consistency_data, window, min_periods, center, f
):
    x, is_constant, no_nans = consistency_data

    if no_nans:
        if f is np.nansum and min_periods == 0:
            pass
        else:
            rolling_f_result = x.rolling(
                window=window, min_periods=min_periods, center=center
            ).sum()
            rolling_apply_f_result = x.rolling(
                window=window, min_periods=min_periods, center=center
            ).apply(func=f, raw=True)
            tm.assert_equal(rolling_f_result, rolling_apply_f_result)


@pytest.mark.parametrize(
    "window,min_periods,center", list(_rolling_consistency_cases())
)
@pytest.mark.parametrize("ddof", [0, 1])
def test_moments_consistency_var(consistency_data, window, min_periods, center, ddof):
    x, is_constant, no_nans = consistency_data

    mean_x = x.rolling(window=window, min_periods=min_periods, center=center).mean()
    var_x = x.rolling(window=window, min_periods=min_periods, center=center).var(
        ddof=ddof
    )
    assert not (var_x < 0).any().any()

    if ddof == 0:
        # check that biased var(x) == mean(x^2) - mean(x)^2
        mean_x2 = (
            (x * x)
            .rolling(window=window, min_periods=min_periods, center=center)
            .mean()
        )
        tm.assert_equal(var_x, mean_x2 - (mean_x * mean_x))


@pytest.mark.parametrize(
    "window,min_periods,center", list(_rolling_consistency_cases())
)
@pytest.mark.parametrize("ddof", [0, 1])
def test_moments_consistency_var_constant(
    consistency_data, window, min_periods, center, ddof
):
    x, is_constant, no_nans = consistency_data

    if is_constant:
        count_x = x.rolling(
            window=window, min_periods=min_periods, center=center
        ).count()
        var_x = x.rolling(window=window, min_periods=min_periods, center=center).var(
            ddof=ddof
        )

        # check that variance of constant series is identically 0
        assert not (var_x > 0).any().any()
        expected = x * np.nan
        expected[count_x >= max(min_periods, 1)] = 0.0
        if ddof == 1:
            expected[count_x < 2] = np.nan
        tm.assert_equal(var_x, expected)


@pytest.mark.parametrize(
    "window,min_periods,center", list(_rolling_consistency_cases())
)
@pytest.mark.parametrize("ddof", [0, 1])
def test_rolling_consistency_std(consistency_data, window, min_periods, center, ddof):
    x, is_constant, no_nans = consistency_data

    var_x = x.rolling(window=window, min_periods=min_periods, center=center).var(
        ddof=ddof
    )
    std_x = x.rolling(window=window, min_periods=min_periods, center=center).std(
        ddof=ddof
    )
    assert not (var_x < 0).any().any()
    assert not (std_x < 0).any().any()

    # check that var(x) == std(x)^2
    tm.assert_equal(var_x, std_x * std_x)


@pytest.mark.parametrize(
    "window,min_periods,center", list(_rolling_consistency_cases())
)
@pytest.mark.parametrize("ddof", [0, 1])
def test_rolling_consistency_cov(consistency_data, window, min_periods, center, ddof):
    x, is_constant, no_nans = consistency_data
    var_x = x.rolling(window=window, min_periods=min_periods, center=center).var(
        ddof=ddof
    )
    assert not (var_x < 0).any().any()

    cov_x_x = x.rolling(window=window, min_periods=min_periods, center=center).cov(
        x, ddof=ddof
    )
    assert not (cov_x_x < 0).any().any()

    # check that var(x) == cov(x, x)
    tm.assert_equal(var_x, cov_x_x)


@pytest.mark.parametrize(
    "window,min_periods,center", list(_rolling_consistency_cases())
)
@pytest.mark.parametrize("ddof", [0, 1])
def test_rolling_consistency_series_cov_corr(
    consistency_data, window, min_periods, center, ddof
):
    x, is_constant, no_nans = consistency_data

    if isinstance(x, Series):
        var_x_plus_y = (
            (x + x)
            .rolling(window=window, min_periods=min_periods, center=center)
            .var(ddof=ddof)
        )
        var_x = x.rolling(window=window, min_periods=min_periods, center=center).var(
            ddof=ddof
        )
        var_y = x.rolling(window=window, min_periods=min_periods, center=center).var(
            ddof=ddof
        )
        cov_x_y = x.rolling(window=window, min_periods=min_periods, center=center).cov(
            x, ddof=ddof
        )
        # check that cov(x, y) == (var(x+y) - var(x) -
        # var(y)) / 2
        tm.assert_equal(cov_x_y, 0.5 * (var_x_plus_y - var_x - var_y))

        # check that corr(x, y) == cov(x, y) / (std(x) *
        # std(y))
        corr_x_y = x.rolling(
            window=window, min_periods=min_periods, center=center
        ).corr(x)
        std_x = x.rolling(window=window, min_periods=min_periods, center=center).std(
            ddof=ddof
        )
        std_y = x.rolling(window=window, min_periods=min_periods, center=center).std(
            ddof=ddof
        )
        tm.assert_equal(corr_x_y, cov_x_y / (std_x * std_y))

        if ddof == 0:
            # check that biased cov(x, y) == mean(x*y) -
            # mean(x)*mean(y)
            mean_x = x.rolling(
                window=window, min_periods=min_periods, center=center
            ).mean()
            mean_y = x.rolling(
                window=window, min_periods=min_periods, center=center
            ).mean()
            mean_x_times_y = (
                (x * x)
                .rolling(window=window, min_periods=min_periods, center=center)
                .mean()
            )
            tm.assert_equal(cov_x_y, mean_x_times_y - (mean_x * mean_y))


@pytest.mark.parametrize(
    "window,min_periods,center", list(_rolling_consistency_cases())
)
def test_rolling_consistency_mean(consistency_data, window, min_periods, center):
    x, is_constant, no_nans = consistency_data

    result = x.rolling(window=window, min_periods=min_periods, center=center).mean()
    expected = (
        x.rolling(window=window, min_periods=min_periods, center=center)
        .sum()
        .divide(
            x.rolling(window=window, min_periods=min_periods, center=center).count()
        )
    )
    tm.assert_equal(result, expected.astype("float64"))


@pytest.mark.parametrize(
    "window,min_periods,center", list(_rolling_consistency_cases())
)
def test_rolling_consistency_constant(consistency_data, window, min_periods, center):
    x, is_constant, no_nans = consistency_data

    if is_constant:
        count_x = x.rolling(
            window=window, min_periods=min_periods, center=center
        ).count()
        mean_x = x.rolling(window=window, min_periods=min_periods, center=center).mean()
        # check that correlation of a series with itself is either 1 or NaN
        corr_x_x = x.rolling(
            window=window, min_periods=min_periods, center=center
        ).corr(x)

        exp = x.max() if isinstance(x, Series) else x.max().max()

        # check mean of constant series
        expected = x * np.nan
        expected[count_x >= max(min_periods, 1)] = exp
        tm.assert_equal(mean_x, expected)

        # check correlation of constant series with itself is NaN
        expected[:] = np.nan
        tm.assert_equal(corr_x_x, expected)


@pytest.mark.parametrize(
    "window,min_periods,center", list(_rolling_consistency_cases())
)
def test_rolling_consistency_var_debiasing_factors(
    consistency_data, window, min_periods, center
):
    x, is_constant, no_nans = consistency_data

    # check variance debiasing factors
    var_unbiased_x = x.rolling(
        window=window, min_periods=min_periods, center=center
    ).var()
    var_biased_x = x.rolling(window=window, min_periods=min_periods, center=center).var(
        ddof=0
    )
    var_debiasing_factors_x = (
        x.rolling(window=window, min_periods=min_periods, center=center)
        .count()
        .divide(
            (
                x.rolling(window=window, min_periods=min_periods, center=center).count()
                - 1.0
            ).replace(0.0, np.nan)
        )
    )
    tm.assert_equal(var_unbiased_x, var_biased_x * var_debiasing_factors_x)
