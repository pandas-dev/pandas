import numpy as np

from pandas import Series
import pandas._testing as tm


def check_pairwise_moment(frame, dispatch, name, **kwargs):
    def get_result(obj, obj2=None):
        return getattr(getattr(obj, dispatch)(**kwargs), name)(obj2)

    result = get_result(frame)
    result = result.loc[(slice(None), 1), 5]
    result.index = result.index.droplevel(1)
    expected = get_result(frame[1], frame[5])
    expected.index = expected.index._with_freq(None)
    tm.assert_series_equal(result, expected, check_names=False)


def ew_func(A, B, com, name, **kwargs):
    return getattr(A.ewm(com, **kwargs), name)(B)


def check_binary_ew(name, A, B):

    result = ew_func(A=A, B=B, com=20, name=name, min_periods=5)
    assert np.isnan(result.values[:14]).all()
    assert not np.isnan(result.values[14:]).any()


def check_binary_ew_min_periods(name, min_periods, A, B):
    # GH 7898
    result = ew_func(A, B, 20, name=name, min_periods=min_periods)
    # binary functions (ewmcov, ewmcorr) with bias=False require at
    # least two values
    assert np.isnan(result.values[:11]).all()
    assert not np.isnan(result.values[11:]).any()

    # check series of length 0
    empty = Series([], dtype=np.float64)
    result = ew_func(empty, empty, 50, name=name, min_periods=min_periods)
    tm.assert_series_equal(result, empty)

    # check series of length 1
    result = ew_func(
        Series([1.0]), Series([1.0]), 50, name=name, min_periods=min_periods
    )
    tm.assert_series_equal(result, Series([np.NaN]))


def moments_consistency_mock_mean(x, mean, mock_mean):
    mean_x = mean(x)
    # check that correlation of a series with itself is either 1 or NaN

    if mock_mean:
        # check that mean equals mock_mean
        expected = mock_mean(x)
        tm.assert_equal(mean_x, expected.astype("float64"))


def moments_consistency_is_constant(x, is_constant, min_periods, count, mean, corr):
    count_x = count(x)
    mean_x = mean(x)
    # check that correlation of a series with itself is either 1 or NaN
    corr_x_x = corr(x, x)

    if is_constant:
        exp = x.max() if isinstance(x, Series) else x.max().max()

        # check mean of constant series
        expected = x * np.nan
        expected[count_x >= max(min_periods, 1)] = exp
        tm.assert_equal(mean_x, expected)

        # check correlation of constant series with itself is NaN
        expected[:] = np.nan
        tm.assert_equal(corr_x_x, expected)


def moments_consistency_var_debiasing_factors(
    x, var_biased, var_unbiased, var_debiasing_factors
):
    if var_unbiased and var_biased and var_debiasing_factors:
        # check variance debiasing factors
        var_unbiased_x = var_unbiased(x)
        var_biased_x = var_biased(x)
        var_debiasing_factors_x = var_debiasing_factors(x)
        tm.assert_equal(var_unbiased_x, var_biased_x * var_debiasing_factors_x)


def moments_consistency_var_data(
    x, is_constant, min_periods, count, mean, var_unbiased, var_biased
):
    count_x = count(x)
    mean_x = mean(x)
    for var in [var_biased, var_unbiased]:
        var_x = var(x)
        assert not (var_x < 0).any().any()

        if var is var_biased:
            # check that biased var(x) == mean(x^2) - mean(x)^2
            mean_x2 = mean(x * x)
            tm.assert_equal(var_x, mean_x2 - (mean_x * mean_x))

        if is_constant:
            # check that variance of constant series is identically 0
            assert not (var_x > 0).any().any()
            expected = x * np.nan
            expected[count_x >= max(min_periods, 1)] = 0.0
            if var is var_unbiased:
                expected[count_x < 2] = np.nan
            tm.assert_equal(var_x, expected)


def moments_consistency_std_data(x, std_unbiased, var_unbiased, std_biased, var_biased):
    for (std, var) in [(std_biased, var_biased), (std_unbiased, var_unbiased)]:
        var_x = var(x)
        std_x = std(x)
        assert not (var_x < 0).any().any()
        assert not (std_x < 0).any().any()

        # check that var(x) == std(x)^2
        tm.assert_equal(var_x, std_x * std_x)


def moments_consistency_cov_data(x, cov_unbiased, var_unbiased, cov_biased, var_biased):
    for (cov, var) in [(cov_biased, var_biased), (cov_unbiased, var_unbiased)]:
        var_x = var(x)
        assert not (var_x < 0).any().any()
        if cov:
            cov_x_x = cov(x, x)
            assert not (cov_x_x < 0).any().any()

            # check that var(x) == cov(x, x)
            tm.assert_equal(var_x, cov_x_x)


def moments_consistency_series_data(
    x,
    corr,
    mean,
    std_biased,
    std_unbiased,
    cov_unbiased,
    var_unbiased,
    var_biased,
    cov_biased,
):
    if isinstance(x, Series):
        y = x
        mean_x = mean(x)
        if not x.isna().equals(y.isna()):
            # can only easily test two Series with similar
            # structure
            pass

        # check that cor(x, y) is symmetric
        corr_x_y = corr(x, y)
        corr_y_x = corr(y, x)
        tm.assert_equal(corr_x_y, corr_y_x)

        for (std, var, cov) in [
            (std_biased, var_biased, cov_biased),
            (std_unbiased, var_unbiased, cov_unbiased),
        ]:
            var_x = var(x)
            std_x = std(x)

            if cov:
                # check that cov(x, y) is symmetric
                cov_x_y = cov(x, y)
                cov_y_x = cov(y, x)
                tm.assert_equal(cov_x_y, cov_y_x)

                # check that cov(x, y) == (var(x+y) - var(x) -
                # var(y)) / 2
                var_x_plus_y = var(x + y)
                var_y = var(y)
                tm.assert_equal(cov_x_y, 0.5 * (var_x_plus_y - var_x - var_y))

                # check that corr(x, y) == cov(x, y) / (std(x) *
                # std(y))
                std_y = std(y)
                tm.assert_equal(corr_x_y, cov_x_y / (std_x * std_y))

                if cov is cov_biased:
                    # check that biased cov(x, y) == mean(x*y) -
                    # mean(x)*mean(y)
                    mean_y = mean(y)
                    mean_x_times_y = mean(x * y)
                    tm.assert_equal(cov_x_y, mean_x_times_y - (mean_x * mean_y))
