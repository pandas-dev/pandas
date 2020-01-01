import numpy as np
from numpy.random import randn
import pytest

import pandas as pd
from pandas import DataFrame, Series, concat
from pandas.tests.window.common import Base, ConsistencyBase
import pandas.util.testing as tm


@pytest.mark.filterwarnings("ignore:can't resolve package:ImportWarning")
class TestMoments(Base):
    def setup_method(self, method):
        self._create_data()

    def test_ewma(self):
        self._check_ew(name="mean")

        vals = pd.Series(np.zeros(1000))
        vals[5] = 1
        result = vals.ewm(span=100, adjust=False).mean().sum()
        assert np.abs(result - 1) < 1e-2

    @pytest.mark.parametrize("adjust", [True, False])
    @pytest.mark.parametrize("ignore_na", [True, False])
    def test_ewma_cases(self, adjust, ignore_na):
        # try adjust/ignore_na args matrix

        s = Series([1.0, 2.0, 4.0, 8.0])

        if adjust:
            expected = Series([1.0, 1.6, 2.736842, 4.923077])
        else:
            expected = Series([1.0, 1.333333, 2.222222, 4.148148])

        result = s.ewm(com=2.0, adjust=adjust, ignore_na=ignore_na).mean()
        tm.assert_series_equal(result, expected)

    def test_ewma_nan_handling(self):
        s = Series([1.0] + [np.nan] * 5 + [1.0])
        result = s.ewm(com=5).mean()
        tm.assert_series_equal(result, Series([1.0] * len(s)))

        s = Series([np.nan] * 2 + [1.0] + [np.nan] * 2 + [1.0])
        result = s.ewm(com=5).mean()
        tm.assert_series_equal(result, Series([np.nan] * 2 + [1.0] * 4))

        # GH 7603
        s0 = Series([np.nan, 1.0, 101.0])
        s1 = Series([1.0, np.nan, 101.0])
        s2 = Series([np.nan, 1.0, np.nan, np.nan, 101.0, np.nan])
        s3 = Series([1.0, np.nan, 101.0, 50.0])
        com = 2.0
        alpha = 1.0 / (1.0 + com)

        def simple_wma(s, w):
            return (s.multiply(w).cumsum() / w.cumsum()).fillna(method="ffill")

        for (s, adjust, ignore_na, w) in [
            (s0, True, False, [np.nan, (1.0 - alpha), 1.0]),
            (s0, True, True, [np.nan, (1.0 - alpha), 1.0]),
            (s0, False, False, [np.nan, (1.0 - alpha), alpha]),
            (s0, False, True, [np.nan, (1.0 - alpha), alpha]),
            (s1, True, False, [(1.0 - alpha) ** 2, np.nan, 1.0]),
            (s1, True, True, [(1.0 - alpha), np.nan, 1.0]),
            (s1, False, False, [(1.0 - alpha) ** 2, np.nan, alpha]),
            (s1, False, True, [(1.0 - alpha), np.nan, alpha]),
            (
                s2,
                True,
                False,
                [np.nan, (1.0 - alpha) ** 3, np.nan, np.nan, 1.0, np.nan],
            ),
            (s2, True, True, [np.nan, (1.0 - alpha), np.nan, np.nan, 1.0, np.nan]),
            (
                s2,
                False,
                False,
                [np.nan, (1.0 - alpha) ** 3, np.nan, np.nan, alpha, np.nan],
            ),
            (s2, False, True, [np.nan, (1.0 - alpha), np.nan, np.nan, alpha, np.nan]),
            (s3, True, False, [(1.0 - alpha) ** 3, np.nan, (1.0 - alpha), 1.0]),
            (s3, True, True, [(1.0 - alpha) ** 2, np.nan, (1.0 - alpha), 1.0]),
            (
                s3,
                False,
                False,
                [
                    (1.0 - alpha) ** 3,
                    np.nan,
                    (1.0 - alpha) * alpha,
                    alpha * ((1.0 - alpha) ** 2 + alpha),
                ],
            ),
            (
                s3,
                False,
                True,
                [(1.0 - alpha) ** 2, np.nan, (1.0 - alpha) * alpha, alpha],
            ),
        ]:
            expected = simple_wma(s, Series(w))
            result = s.ewm(com=com, adjust=adjust, ignore_na=ignore_na).mean()

            tm.assert_series_equal(result, expected)
            if ignore_na is False:
                # check that ignore_na defaults to False
                result = s.ewm(com=com, adjust=adjust).mean()
                tm.assert_series_equal(result, expected)

    def test_ewmvar(self):
        self._check_ew(name="var")

    def test_ewmvol(self):
        self._check_ew(name="vol")

    def test_ewma_span_com_args(self):
        A = self.series.ewm(com=9.5).mean()
        B = self.series.ewm(span=20).mean()
        tm.assert_almost_equal(A, B)

        with pytest.raises(ValueError):
            self.series.ewm(com=9.5, span=20)
        with pytest.raises(ValueError):
            self.series.ewm().mean()

    def test_ewma_halflife_arg(self):
        A = self.series.ewm(com=13.932726172912965).mean()
        B = self.series.ewm(halflife=10.0).mean()
        tm.assert_almost_equal(A, B)

        with pytest.raises(ValueError):
            self.series.ewm(span=20, halflife=50)
        with pytest.raises(ValueError):
            self.series.ewm(com=9.5, halflife=50)
        with pytest.raises(ValueError):
            self.series.ewm(com=9.5, span=20, halflife=50)
        with pytest.raises(ValueError):
            self.series.ewm()

    def test_ewm_alpha(self):
        # GH 10789
        s = Series(self.arr)
        a = s.ewm(alpha=0.61722699889169674).mean()
        b = s.ewm(com=0.62014947789973052).mean()
        c = s.ewm(span=2.240298955799461).mean()
        d = s.ewm(halflife=0.721792864318).mean()
        tm.assert_series_equal(a, b)
        tm.assert_series_equal(a, c)
        tm.assert_series_equal(a, d)

    def test_ewm_alpha_arg(self):
        # GH 10789
        s = self.series
        with pytest.raises(ValueError):
            s.ewm()
        with pytest.raises(ValueError):
            s.ewm(com=10.0, alpha=0.5)
        with pytest.raises(ValueError):
            s.ewm(span=10.0, alpha=0.5)
        with pytest.raises(ValueError):
            s.ewm(halflife=10.0, alpha=0.5)

    def test_ewm_domain_checks(self):
        # GH 12492
        s = Series(self.arr)
        msg = "comass must satisfy: comass >= 0"
        with pytest.raises(ValueError, match=msg):
            s.ewm(com=-0.1)
        s.ewm(com=0.0)
        s.ewm(com=0.1)

        msg = "span must satisfy: span >= 1"
        with pytest.raises(ValueError, match=msg):
            s.ewm(span=-0.1)
        with pytest.raises(ValueError, match=msg):
            s.ewm(span=0.0)
        with pytest.raises(ValueError, match=msg):
            s.ewm(span=0.9)
        s.ewm(span=1.0)
        s.ewm(span=1.1)

        msg = "halflife must satisfy: halflife > 0"
        with pytest.raises(ValueError, match=msg):
            s.ewm(halflife=-0.1)
        with pytest.raises(ValueError, match=msg):
            s.ewm(halflife=0.0)
        s.ewm(halflife=0.1)

        msg = "alpha must satisfy: 0 < alpha <= 1"
        with pytest.raises(ValueError, match=msg):
            s.ewm(alpha=-0.1)
        with pytest.raises(ValueError, match=msg):
            s.ewm(alpha=0.0)
        s.ewm(alpha=0.1)
        s.ewm(alpha=1.0)
        with pytest.raises(ValueError, match=msg):
            s.ewm(alpha=1.1)

    @pytest.mark.parametrize("method", ["mean", "vol", "var"])
    def test_ew_empty_series(self, method):
        vals = pd.Series([], dtype=np.float64)

        ewm = vals.ewm(3)
        result = getattr(ewm, method)()
        tm.assert_almost_equal(result, vals)

    def _check_ew(self, name=None, preserve_nan=False):
        series_result = getattr(self.series.ewm(com=10), name)()
        assert isinstance(series_result, Series)

        frame_result = getattr(self.frame.ewm(com=10), name)()
        assert type(frame_result) == DataFrame

        result = getattr(self.series.ewm(com=10), name)()
        if preserve_nan:
            assert result[self._nan_locs].isna().all()

        # excluding NaNs correctly
        arr = randn(50)
        arr[:10] = np.NaN
        arr[-10:] = np.NaN
        s = Series(arr)

        # check min_periods
        # GH 7898
        result = getattr(s.ewm(com=50, min_periods=2), name)()
        assert result[:11].isna().all()
        assert not result[11:].isna().any()

        for min_periods in (0, 1):
            result = getattr(s.ewm(com=50, min_periods=min_periods), name)()
            if name == "mean":
                assert result[:10].isna().all()
                assert not result[10:].isna().any()
            else:
                # ewm.std, ewm.vol, ewm.var (with bias=False) require at least
                # two values
                assert result[:11].isna().all()
                assert not result[11:].isna().any()

            # check series of length 0
            result = getattr(
                Series(dtype=object).ewm(com=50, min_periods=min_periods), name
            )()
            tm.assert_series_equal(result, Series(dtype="float64"))

            # check series of length 1
            result = getattr(Series([1.0]).ewm(50, min_periods=min_periods), name)()
            if name == "mean":
                tm.assert_series_equal(result, Series([1.0]))
            else:
                # ewm.std, ewm.vol, ewm.var with bias=False require at least
                # two values
                tm.assert_series_equal(result, Series([np.NaN]))

        # pass in ints
        result2 = getattr(Series(np.arange(50)).ewm(span=10), name)()
        assert result2.dtype == np.float_


class TestEwmMomentsConsistency(ConsistencyBase):
    def setup_method(self, method):
        self._create_data()

    def test_ewmcov(self):
        self._check_binary_ew("cov")

    def test_ewmcov_pairwise(self):
        self._check_pairwise_moment("ewm", "cov", span=10, min_periods=5)

    def test_ewmcorr(self):
        self._check_binary_ew("corr")

    def test_ewmcorr_pairwise(self):
        self._check_pairwise_moment("ewm", "corr", span=10, min_periods=5)

    def _check_binary_ew(self, name):
        def func(A, B, com, **kwargs):
            return getattr(A.ewm(com, **kwargs), name)(B)

        A = Series(randn(50), index=np.arange(50))
        B = A[2:] + randn(48)

        A[:10] = np.NaN
        B[-10:] = np.NaN

        result = func(A, B, 20, min_periods=5)
        assert np.isnan(result.values[:14]).all()
        assert not np.isnan(result.values[14:]).any()

        # GH 7898
        for min_periods in (0, 1, 2):
            result = func(A, B, 20, min_periods=min_periods)
            # binary functions (ewmcov, ewmcorr) with bias=False require at
            # least two values
            assert np.isnan(result.values[:11]).all()
            assert not np.isnan(result.values[11:]).any()

            # check series of length 0
            empty = Series([], dtype=np.float64)
            result = func(empty, empty, 50, min_periods=min_periods)
            tm.assert_series_equal(result, empty)

            # check series of length 1
            result = func(Series([1.0]), Series([1.0]), 50, min_periods=min_periods)
            tm.assert_series_equal(result, Series([np.NaN]))

        msg = "Input arrays must be of the same type!"
        # exception raised is Exception
        with pytest.raises(Exception, match=msg):
            func(A, randn(50), 20, min_periods=5)

    @pytest.mark.slow
    @pytest.mark.parametrize("min_periods", [0, 1, 2, 3, 4])
    @pytest.mark.parametrize("adjust", [True, False])
    @pytest.mark.parametrize("ignore_na", [True, False])
    def test_ewm_consistency(self, min_periods, adjust, ignore_na):
        def _weights(s, com, adjust, ignore_na):
            if isinstance(s, DataFrame):
                if not len(s.columns):
                    return DataFrame(index=s.index, columns=s.columns)
                w = concat(
                    [
                        _weights(
                            s.iloc[:, i], com=com, adjust=adjust, ignore_na=ignore_na
                        )
                        for i, _ in enumerate(s.columns)
                    ],
                    axis=1,
                )
                w.index = s.index
                w.columns = s.columns
                return w

            w = Series(np.nan, index=s.index)
            alpha = 1.0 / (1.0 + com)
            if ignore_na:
                w[s.notna()] = _weights(
                    s[s.notna()], com=com, adjust=adjust, ignore_na=False
                )
            elif adjust:
                for i in range(len(s)):
                    if s.iat[i] == s.iat[i]:
                        w.iat[i] = pow(1.0 / (1.0 - alpha), i)
            else:
                sum_wts = 0.0
                prev_i = -1
                for i in range(len(s)):
                    if s.iat[i] == s.iat[i]:
                        if prev_i == -1:
                            w.iat[i] = 1.0
                        else:
                            w.iat[i] = alpha * sum_wts / pow(1.0 - alpha, i - prev_i)
                        sum_wts += w.iat[i]
                        prev_i = i
            return w

        def _variance_debiasing_factors(s, com, adjust, ignore_na):
            weights = _weights(s, com=com, adjust=adjust, ignore_na=ignore_na)
            cum_sum = weights.cumsum().fillna(method="ffill")
            cum_sum_sq = (weights * weights).cumsum().fillna(method="ffill")
            numerator = cum_sum * cum_sum
            denominator = numerator - cum_sum_sq
            denominator[denominator <= 0.0] = np.nan
            return numerator / denominator

        def _ewma(s, com, min_periods, adjust, ignore_na):
            weights = _weights(s, com=com, adjust=adjust, ignore_na=ignore_na)
            result = (
                s.multiply(weights)
                .cumsum()
                .divide(weights.cumsum())
                .fillna(method="ffill")
            )
            result[
                s.expanding().count() < (max(min_periods, 1) if min_periods else 1)
            ] = np.nan
            return result

        com = 3.0
        self._test_moments_consistency_mock_mean(
            mean=lambda x: x.ewm(
                com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na
            ).mean(),
            mock_mean=lambda x: _ewma(
                x, com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na
            ),
        )

        self._test_moments_consistency_is_constant(
            min_periods=min_periods,
            count=lambda x: x.expanding().count(),
            mean=lambda x: x.ewm(
                com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na
            ).mean(),
            corr=lambda x, y: x.ewm(
                com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na
            ).corr(y),
        )

        self._test_moments_consistency_var_debiasing_factors(
            var_unbiased=lambda x: (
                x.ewm(
                    com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na
                ).var(bias=False)
            ),
            var_biased=lambda x: (
                x.ewm(
                    com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na
                ).var(bias=True)
            ),
            var_debiasing_factors=lambda x: (
                _variance_debiasing_factors(
                    x, com=com, adjust=adjust, ignore_na=ignore_na
                )
            ),
        )
        # test consistency between different ewm* moments
        self._test_moments_consistency(
            min_periods=min_periods,
            count=lambda x: x.expanding().count(),
            mean=lambda x: x.ewm(
                com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na
            ).mean(),
            corr=lambda x, y: x.ewm(
                com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na
            ).corr(y),
            var_unbiased=lambda x: (
                x.ewm(
                    com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na
                ).var(bias=False)
            ),
            std_unbiased=lambda x: (
                x.ewm(
                    com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na
                ).std(bias=False)
            ),
            cov_unbiased=lambda x, y: (
                x.ewm(
                    com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na
                ).cov(y, bias=False)
            ),
            var_biased=lambda x: (
                x.ewm(
                    com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na
                ).var(bias=True)
            ),
            std_biased=lambda x: x.ewm(
                com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na
            ).std(bias=True),
            cov_biased=lambda x, y: (
                x.ewm(
                    com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na
                ).cov(y, bias=True)
            ),
        )
