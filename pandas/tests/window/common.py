from datetime import datetime

import numpy as np
from numpy.random import randn

from pandas import DataFrame, Series, bdate_range, notna
import pandas.util.testing as tm

N, K = 100, 10


class Base:

    _nan_locs = np.arange(20, 40)
    _inf_locs = np.array([])

    def _create_data(self):
        arr = randn(N)
        arr[self._nan_locs] = np.NaN

        self.arr = arr
        self.rng = bdate_range(datetime(2009, 1, 1), periods=N)
        self.series = Series(arr.copy(), index=self.rng)
        self.frame = DataFrame(randn(N, K), index=self.rng, columns=np.arange(K))


# create the data only once as we are not setting it
def _create_consistency_data():
    def create_series():
        return [
            Series(dtype=object),
            Series([np.nan]),
            Series([np.nan, np.nan]),
            Series([3.0]),
            Series([np.nan, 3.0]),
            Series([3.0, np.nan]),
            Series([1.0, 3.0]),
            Series([2.0, 2.0]),
            Series([3.0, 1.0]),
            Series(
                [5.0, 5.0, 5.0, 5.0, np.nan, np.nan, np.nan, 5.0, 5.0, np.nan, np.nan]
            ),
            Series(
                [
                    np.nan,
                    5.0,
                    5.0,
                    5.0,
                    np.nan,
                    np.nan,
                    np.nan,
                    5.0,
                    5.0,
                    np.nan,
                    np.nan,
                ]
            ),
            Series(
                [
                    np.nan,
                    np.nan,
                    5.0,
                    5.0,
                    np.nan,
                    np.nan,
                    np.nan,
                    5.0,
                    5.0,
                    np.nan,
                    np.nan,
                ]
            ),
            Series(
                [
                    np.nan,
                    3.0,
                    np.nan,
                    3.0,
                    4.0,
                    5.0,
                    6.0,
                    np.nan,
                    np.nan,
                    7.0,
                    12.0,
                    13.0,
                    14.0,
                    15.0,
                ]
            ),
            Series(
                [
                    np.nan,
                    5.0,
                    np.nan,
                    2.0,
                    4.0,
                    0.0,
                    9.0,
                    np.nan,
                    np.nan,
                    3.0,
                    12.0,
                    13.0,
                    14.0,
                    15.0,
                ]
            ),
            Series(
                [
                    2.0,
                    3.0,
                    np.nan,
                    3.0,
                    4.0,
                    5.0,
                    6.0,
                    np.nan,
                    np.nan,
                    7.0,
                    12.0,
                    13.0,
                    14.0,
                    15.0,
                ]
            ),
            Series(
                [
                    2.0,
                    5.0,
                    np.nan,
                    2.0,
                    4.0,
                    0.0,
                    9.0,
                    np.nan,
                    np.nan,
                    3.0,
                    12.0,
                    13.0,
                    14.0,
                    15.0,
                ]
            ),
            Series(range(10)),
            Series(range(20, 0, -2)),
        ]

    def create_dataframes():
        return [
            DataFrame(),
            DataFrame(columns=["a"]),
            DataFrame(columns=["a", "a"]),
            DataFrame(columns=["a", "b"]),
            DataFrame(np.arange(10).reshape((5, 2))),
            DataFrame(np.arange(25).reshape((5, 5))),
            DataFrame(np.arange(25).reshape((5, 5)), columns=["a", "b", 99, "d", "d"]),
        ] + [DataFrame(s) for s in create_series()]

    def is_constant(x):
        values = x.values.ravel()
        return len(set(values[notna(values)])) == 1

    def no_nans(x):
        return x.notna().all().all()

    # data is a tuple(object, is_constant, no_nans)
    data = create_series() + create_dataframes()

    return [(x, is_constant(x), no_nans(x)) for x in data]


_consistency_data = _create_consistency_data()


class ConsistencyBase(Base):
    base_functions = [
        (lambda v: Series(v).count(), None, "count"),
        (lambda v: Series(v).max(), None, "max"),
        (lambda v: Series(v).min(), None, "min"),
        (lambda v: Series(v).sum(), None, "sum"),
        (lambda v: Series(v).mean(), None, "mean"),
        (lambda v: Series(v).std(), 1, "std"),
        (lambda v: Series(v).cov(Series(v)), None, "cov"),
        (lambda v: Series(v).corr(Series(v)), None, "corr"),
        (lambda v: Series(v).var(), 1, "var"),
        # restore once GH 8086 is fixed
        # lambda v: Series(v).skew(), 3, 'skew'),
        # (lambda v: Series(v).kurt(), 4, 'kurt'),
        # restore once GH 8084 is fixed
        # lambda v: Series(v).quantile(0.3), None, 'quantile'),
        (lambda v: Series(v).median(), None, "median"),
        (np.nanmax, 1, "max"),
        (np.nanmin, 1, "min"),
        (np.nansum, 1, "sum"),
        (np.nanmean, 1, "mean"),
        (lambda v: np.nanstd(v, ddof=1), 1, "std"),
        (lambda v: np.nanvar(v, ddof=1), 1, "var"),
        (np.nanmedian, 1, "median"),
    ]
    no_nan_functions = [
        (np.max, None, "max"),
        (np.min, None, "min"),
        (np.sum, None, "sum"),
        (np.mean, None, "mean"),
        (lambda v: np.std(v, ddof=1), 1, "std"),
        (lambda v: np.var(v, ddof=1), 1, "var"),
        (np.median, None, "median"),
    ]

    def _create_data(self):
        super()._create_data()
        self.data = _consistency_data

    def _test_moments_consistency_mock_mean(self, mean, mock_mean):
        for (x, is_constant, no_nans) in self.data:
            mean_x = mean(x)
            # check that correlation of a series with itself is either 1 or NaN

            if mock_mean:
                # check that mean equals mock_mean
                expected = mock_mean(x)
                tm.assert_equal(mean_x, expected.astype("float64"))

    def _test_moments_consistency_is_constant(self, min_periods, count, mean, corr):
        for (x, is_constant, no_nans) in self.data:
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

    def _test_moments_consistency_var_debiasing_factors(
        self, var_biased=None, var_unbiased=None, var_debiasing_factors=None
    ):
        for (x, is_constant, no_nans) in self.data:
            if var_unbiased and var_biased and var_debiasing_factors:
                # check variance debiasing factors
                var_unbiased_x = var_unbiased(x)
                var_biased_x = var_biased(x)
                var_debiasing_factors_x = var_debiasing_factors(x)
                tm.assert_equal(var_unbiased_x, var_biased_x * var_debiasing_factors_x)

    def _test_moments_consistency(
        self,
        min_periods,
        count,
        mean,
        corr,
        var_unbiased=None,
        std_unbiased=None,
        cov_unbiased=None,
        var_biased=None,
        std_biased=None,
        cov_biased=None,
    ):

        for (x, is_constant, no_nans) in self.data:
            count_x = count(x)
            mean_x = mean(x)

            for (std, var, cov) in [
                (std_biased, var_biased, cov_biased),
                (std_unbiased, var_unbiased, cov_unbiased),
            ]:

                # check that var(x), std(x), and cov(x) are all >= 0
                var_x = var(x)
                std_x = std(x)
                assert not (var_x < 0).any().any()
                assert not (std_x < 0).any().any()
                if cov:
                    cov_x_x = cov(x, x)
                    assert not (cov_x_x < 0).any().any()

                    # check that var(x) == cov(x, x)
                    tm.assert_equal(var_x, cov_x_x)

                # check that var(x) == std(x)^2
                tm.assert_equal(var_x, std_x * std_x)

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

                if isinstance(x, Series):
                    for (y, is_constant, no_nans) in self.data:
                        if not x.isna().equals(y.isna()):
                            # can only easily test two Series with similar
                            # structure
                            continue

                        # check that cor(x, y) is symmetric
                        corr_x_y = corr(x, y)
                        corr_y_x = corr(y, x)
                        tm.assert_equal(corr_x_y, corr_y_x)

                        if cov:
                            # check that cov(x, y) is symmetric
                            cov_x_y = cov(x, y)
                            cov_y_x = cov(y, x)
                            tm.assert_equal(cov_x_y, cov_y_x)

                            # check that cov(x, y) == (var(x+y) - var(x) -
                            # var(y)) / 2
                            var_x_plus_y = var(x + y)
                            var_y = var(y)
                            tm.assert_equal(
                                cov_x_y, 0.5 * (var_x_plus_y - var_x - var_y)
                            )

                            # check that corr(x, y) == cov(x, y) / (std(x) *
                            # std(y))
                            std_y = std(y)
                            tm.assert_equal(corr_x_y, cov_x_y / (std_x * std_y))

                            if cov is cov_biased:
                                # check that biased cov(x, y) == mean(x*y) -
                                # mean(x)*mean(y)
                                mean_y = mean(y)
                                mean_x_times_y = mean(x * y)
                                tm.assert_equal(
                                    cov_x_y, mean_x_times_y - (mean_x * mean_y)
                                )

    def _check_pairwise_moment(self, dispatch, name, **kwargs):
        def get_result(obj, obj2=None):
            return getattr(getattr(obj, dispatch)(**kwargs), name)(obj2)

        result = get_result(self.frame)
        result = result.loc[(slice(None), 1), 5]
        result.index = result.index.droplevel(1)
        expected = get_result(self.frame[1], self.frame[5])
        tm.assert_series_equal(result, expected, check_names=False)
