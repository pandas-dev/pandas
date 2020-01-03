import copy
from datetime import datetime
import warnings

import numpy as np
from numpy.random import randn
import pytest

import pandas.util._test_decorators as td

import pandas as pd
from pandas import DataFrame, Index, Series, isna, notna
from pandas.core.window.common import _flex_binary_moment
from pandas.tests.window.common import Base, ConsistencyBase
import pandas.util.testing as tm

import pandas.tseries.offsets as offsets


@pytest.mark.filterwarnings("ignore:can't resolve package:ImportWarning")
class TestMoments(Base):
    def setup_method(self, method):
        self._create_data()

    def test_centered_axis_validation(self):

        # ok
        Series(np.ones(10)).rolling(window=3, center=True, axis=0).mean()

        # bad axis
        with pytest.raises(ValueError):
            Series(np.ones(10)).rolling(window=3, center=True, axis=1).mean()

        # ok ok
        DataFrame(np.ones((10, 10))).rolling(window=3, center=True, axis=0).mean()
        DataFrame(np.ones((10, 10))).rolling(window=3, center=True, axis=1).mean()

        # bad axis
        with pytest.raises(ValueError):
            (DataFrame(np.ones((10, 10))).rolling(window=3, center=True, axis=2).mean())

    def test_rolling_sum(self, raw):
        self._check_moment_func(
            np.nansum, name="sum", zero_min_periods_equal=False, raw=raw
        )

    def test_rolling_count(self, raw):
        counter = lambda x: np.isfinite(x).astype(float).sum()
        self._check_moment_func(
            counter, name="count", has_min_periods=False, fill_value=0, raw=raw
        )

    def test_rolling_mean(self, raw):
        self._check_moment_func(np.mean, name="mean", raw=raw)

    @td.skip_if_no_scipy
    def test_cmov_mean(self):
        # GH 8238
        vals = np.array(
            [6.95, 15.21, 4.72, 9.12, 13.81, 13.49, 16.68, 9.48, 10.63, 14.48]
        )
        result = Series(vals).rolling(5, center=True).mean()
        expected = Series(
            [
                np.nan,
                np.nan,
                9.962,
                11.27,
                11.564,
                12.516,
                12.818,
                12.952,
                np.nan,
                np.nan,
            ]
        )
        tm.assert_series_equal(expected, result)

    @td.skip_if_no_scipy
    def test_cmov_window(self):
        # GH 8238
        vals = np.array(
            [6.95, 15.21, 4.72, 9.12, 13.81, 13.49, 16.68, 9.48, 10.63, 14.48]
        )
        result = Series(vals).rolling(5, win_type="boxcar", center=True).mean()
        expected = Series(
            [
                np.nan,
                np.nan,
                9.962,
                11.27,
                11.564,
                12.516,
                12.818,
                12.952,
                np.nan,
                np.nan,
            ]
        )
        tm.assert_series_equal(expected, result)

    @td.skip_if_no_scipy
    def test_cmov_window_corner(self):
        # GH 8238
        # all nan
        vals = pd.Series([np.nan] * 10)
        result = vals.rolling(5, center=True, win_type="boxcar").mean()
        assert np.isnan(result).all()

        # empty
        vals = pd.Series([], dtype=object)
        result = vals.rolling(5, center=True, win_type="boxcar").mean()
        assert len(result) == 0

        # shorter than window
        vals = pd.Series(np.random.randn(5))
        result = vals.rolling(10, win_type="boxcar").mean()
        assert np.isnan(result).all()
        assert len(result) == 5

    @td.skip_if_no_scipy
    @pytest.mark.parametrize(
        "f,xp",
        [
            (
                "mean",
                [
                    [np.nan, np.nan],
                    [np.nan, np.nan],
                    [9.252, 9.392],
                    [8.644, 9.906],
                    [8.87, 10.208],
                    [6.81, 8.588],
                    [7.792, 8.644],
                    [9.05, 7.824],
                    [np.nan, np.nan],
                    [np.nan, np.nan],
                ],
            ),
            (
                "std",
                [
                    [np.nan, np.nan],
                    [np.nan, np.nan],
                    [3.789706, 4.068313],
                    [3.429232, 3.237411],
                    [3.589269, 3.220810],
                    [3.405195, 2.380655],
                    [3.281839, 2.369869],
                    [3.676846, 1.801799],
                    [np.nan, np.nan],
                    [np.nan, np.nan],
                ],
            ),
            (
                "var",
                [
                    [np.nan, np.nan],
                    [np.nan, np.nan],
                    [14.36187, 16.55117],
                    [11.75963, 10.48083],
                    [12.88285, 10.37362],
                    [11.59535, 5.66752],
                    [10.77047, 5.61628],
                    [13.51920, 3.24648],
                    [np.nan, np.nan],
                    [np.nan, np.nan],
                ],
            ),
            (
                "sum",
                [
                    [np.nan, np.nan],
                    [np.nan, np.nan],
                    [46.26, 46.96],
                    [43.22, 49.53],
                    [44.35, 51.04],
                    [34.05, 42.94],
                    [38.96, 43.22],
                    [45.25, 39.12],
                    [np.nan, np.nan],
                    [np.nan, np.nan],
                ],
            ),
        ],
    )
    def test_cmov_window_frame(self, f, xp):
        # Gh 8238
        df = DataFrame(
            np.array(
                [
                    [12.18, 3.64],
                    [10.18, 9.16],
                    [13.24, 14.61],
                    [4.51, 8.11],
                    [6.15, 11.44],
                    [9.14, 6.21],
                    [11.31, 10.67],
                    [2.94, 6.51],
                    [9.42, 8.39],
                    [12.44, 7.34],
                ]
            )
        )
        xp = DataFrame(np.array(xp))

        roll = df.rolling(5, win_type="boxcar", center=True)
        rs = getattr(roll, f)()

        tm.assert_frame_equal(xp, rs)

    @td.skip_if_no_scipy
    def test_cmov_window_na_min_periods(self):
        # min_periods
        vals = Series(np.random.randn(10))
        vals[4] = np.nan
        vals[8] = np.nan

        xp = vals.rolling(5, min_periods=4, center=True).mean()
        rs = vals.rolling(5, win_type="boxcar", min_periods=4, center=True).mean()
        tm.assert_series_equal(xp, rs)

    @td.skip_if_no_scipy
    def test_cmov_window_regular(self, win_types):
        # GH 8238
        vals = np.array(
            [6.95, 15.21, 4.72, 9.12, 13.81, 13.49, 16.68, 9.48, 10.63, 14.48]
        )
        xps = {
            "hamming": [
                np.nan,
                np.nan,
                8.71384,
                9.56348,
                12.38009,
                14.03687,
                13.8567,
                11.81473,
                np.nan,
                np.nan,
            ],
            "triang": [
                np.nan,
                np.nan,
                9.28667,
                10.34667,
                12.00556,
                13.33889,
                13.38,
                12.33667,
                np.nan,
                np.nan,
            ],
            "barthann": [
                np.nan,
                np.nan,
                8.4425,
                9.1925,
                12.5575,
                14.3675,
                14.0825,
                11.5675,
                np.nan,
                np.nan,
            ],
            "bohman": [
                np.nan,
                np.nan,
                7.61599,
                9.1764,
                12.83559,
                14.17267,
                14.65923,
                11.10401,
                np.nan,
                np.nan,
            ],
            "blackmanharris": [
                np.nan,
                np.nan,
                6.97691,
                9.16438,
                13.05052,
                14.02156,
                15.10512,
                10.74574,
                np.nan,
                np.nan,
            ],
            "nuttall": [
                np.nan,
                np.nan,
                7.04618,
                9.16786,
                13.02671,
                14.03559,
                15.05657,
                10.78514,
                np.nan,
                np.nan,
            ],
            "blackman": [
                np.nan,
                np.nan,
                7.73345,
                9.17869,
                12.79607,
                14.20036,
                14.57726,
                11.16988,
                np.nan,
                np.nan,
            ],
            "bartlett": [
                np.nan,
                np.nan,
                8.4425,
                9.1925,
                12.5575,
                14.3675,
                14.0825,
                11.5675,
                np.nan,
                np.nan,
            ],
        }

        xp = Series(xps[win_types])
        rs = Series(vals).rolling(5, win_type=win_types, center=True).mean()
        tm.assert_series_equal(xp, rs)

    @td.skip_if_no_scipy
    def test_cmov_window_regular_linear_range(self, win_types):
        # GH 8238
        vals = np.array(range(10), dtype=np.float)
        xp = vals.copy()
        xp[:2] = np.nan
        xp[-2:] = np.nan
        xp = Series(xp)

        rs = Series(vals).rolling(5, win_type=win_types, center=True).mean()
        tm.assert_series_equal(xp, rs)

    @td.skip_if_no_scipy
    def test_cmov_window_regular_missing_data(self, win_types):
        # GH 8238
        vals = np.array(
            [6.95, 15.21, 4.72, 9.12, 13.81, 13.49, 16.68, np.nan, 10.63, 14.48]
        )
        xps = {
            "bartlett": [
                np.nan,
                np.nan,
                9.70333,
                10.5225,
                8.4425,
                9.1925,
                12.5575,
                14.3675,
                15.61667,
                13.655,
            ],
            "blackman": [
                np.nan,
                np.nan,
                9.04582,
                11.41536,
                7.73345,
                9.17869,
                12.79607,
                14.20036,
                15.8706,
                13.655,
            ],
            "barthann": [
                np.nan,
                np.nan,
                9.70333,
                10.5225,
                8.4425,
                9.1925,
                12.5575,
                14.3675,
                15.61667,
                13.655,
            ],
            "bohman": [
                np.nan,
                np.nan,
                8.9444,
                11.56327,
                7.61599,
                9.1764,
                12.83559,
                14.17267,
                15.90976,
                13.655,
            ],
            "hamming": [
                np.nan,
                np.nan,
                9.59321,
                10.29694,
                8.71384,
                9.56348,
                12.38009,
                14.20565,
                15.24694,
                13.69758,
            ],
            "nuttall": [
                np.nan,
                np.nan,
                8.47693,
                12.2821,
                7.04618,
                9.16786,
                13.02671,
                14.03673,
                16.08759,
                13.65553,
            ],
            "triang": [
                np.nan,
                np.nan,
                9.33167,
                9.76125,
                9.28667,
                10.34667,
                12.00556,
                13.82125,
                14.49429,
                13.765,
            ],
            "blackmanharris": [
                np.nan,
                np.nan,
                8.42526,
                12.36824,
                6.97691,
                9.16438,
                13.05052,
                14.02175,
                16.1098,
                13.65509,
            ],
        }

        xp = Series(xps[win_types])
        rs = Series(vals).rolling(5, win_type=win_types, min_periods=3).mean()
        tm.assert_series_equal(xp, rs)

    @td.skip_if_no_scipy
    def test_cmov_window_special(self, win_types_special):
        # GH 8238
        kwds = {
            "kaiser": {"beta": 1.0},
            "gaussian": {"std": 1.0},
            "general_gaussian": {"power": 2.0, "width": 2.0},
            "exponential": {"tau": 10},
        }

        vals = np.array(
            [6.95, 15.21, 4.72, 9.12, 13.81, 13.49, 16.68, 9.48, 10.63, 14.48]
        )

        xps = {
            "gaussian": [
                np.nan,
                np.nan,
                8.97297,
                9.76077,
                12.24763,
                13.89053,
                13.65671,
                12.01002,
                np.nan,
                np.nan,
            ],
            "general_gaussian": [
                np.nan,
                np.nan,
                9.85011,
                10.71589,
                11.73161,
                13.08516,
                12.95111,
                12.74577,
                np.nan,
                np.nan,
            ],
            "kaiser": [
                np.nan,
                np.nan,
                9.86851,
                11.02969,
                11.65161,
                12.75129,
                12.90702,
                12.83757,
                np.nan,
                np.nan,
            ],
            "exponential": [
                np.nan,
                np.nan,
                9.83364,
                11.10472,
                11.64551,
                12.66138,
                12.92379,
                12.83770,
                np.nan,
                np.nan,
            ],
        }

        xp = Series(xps[win_types_special])
        rs = (
            Series(vals)
            .rolling(5, win_type=win_types_special, center=True)
            .mean(**kwds[win_types_special])
        )
        tm.assert_series_equal(xp, rs)

    @td.skip_if_no_scipy
    def test_cmov_window_special_linear_range(self, win_types_special):
        # GH 8238
        kwds = {
            "kaiser": {"beta": 1.0},
            "gaussian": {"std": 1.0},
            "general_gaussian": {"power": 2.0, "width": 2.0},
            "slepian": {"width": 0.5},
            "exponential": {"tau": 10},
        }

        vals = np.array(range(10), dtype=np.float)
        xp = vals.copy()
        xp[:2] = np.nan
        xp[-2:] = np.nan
        xp = Series(xp)

        rs = (
            Series(vals)
            .rolling(5, win_type=win_types_special, center=True)
            .mean(**kwds[win_types_special])
        )
        tm.assert_series_equal(xp, rs)

    def test_rolling_median(self, raw):
        self._check_moment_func(np.median, name="median", raw=raw)

    def test_rolling_min(self, raw):
        self._check_moment_func(np.min, name="min", raw=raw)

        a = pd.Series([1, 2, 3, 4, 5])
        result = a.rolling(window=100, min_periods=1).min()
        expected = pd.Series(np.ones(len(a)))
        tm.assert_series_equal(result, expected)

        with pytest.raises(ValueError):
            pd.Series([1, 2, 3]).rolling(window=3, min_periods=5).min()

    def test_rolling_max(self, raw):
        self._check_moment_func(np.max, name="max", raw=raw)

        a = pd.Series([1, 2, 3, 4, 5], dtype=np.float64)
        b = a.rolling(window=100, min_periods=1).max()
        tm.assert_almost_equal(a, b)

        with pytest.raises(ValueError):
            pd.Series([1, 2, 3]).rolling(window=3, min_periods=5).max()

    @pytest.mark.parametrize("q", [0.0, 0.1, 0.5, 0.9, 1.0])
    def test_rolling_quantile(self, q, raw):
        def scoreatpercentile(a, per):
            values = np.sort(a, axis=0)

            idx = int(per / 1.0 * (values.shape[0] - 1))

            if idx == values.shape[0] - 1:
                retval = values[-1]

            else:
                qlow = float(idx) / float(values.shape[0] - 1)
                qhig = float(idx + 1) / float(values.shape[0] - 1)
                vlow = values[idx]
                vhig = values[idx + 1]
                retval = vlow + (vhig - vlow) * (per - qlow) / (qhig - qlow)

            return retval

        def quantile_func(x):
            return scoreatpercentile(x, q)

        self._check_moment_func(quantile_func, name="quantile", quantile=q, raw=raw)

    def test_rolling_quantile_np_percentile(self):
        # #9413: Tests that rolling window's quantile default behavior
        # is analogous to Numpy's percentile
        row = 10
        col = 5
        idx = pd.date_range("20100101", periods=row, freq="B")
        df = DataFrame(np.random.rand(row * col).reshape((row, -1)), index=idx)

        df_quantile = df.quantile([0.25, 0.5, 0.75], axis=0)
        np_percentile = np.percentile(df, [25, 50, 75], axis=0)

        tm.assert_almost_equal(df_quantile.values, np.array(np_percentile))

    @pytest.mark.parametrize("quantile", [0.0, 0.1, 0.45, 0.5, 1])
    @pytest.mark.parametrize(
        "interpolation", ["linear", "lower", "higher", "nearest", "midpoint"]
    )
    @pytest.mark.parametrize(
        "data",
        [
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
            [8.0, 1.0, 3.0, 4.0, 5.0, 2.0, 6.0, 7.0],
            [0.0, np.nan, 0.2, np.nan, 0.4],
            [np.nan, np.nan, np.nan, np.nan],
            [np.nan, 0.1, np.nan, 0.3, 0.4, 0.5],
            [0.5],
            [np.nan, 0.7, 0.6],
        ],
    )
    def test_rolling_quantile_interpolation_options(
        self, quantile, interpolation, data
    ):
        # Tests that rolling window's quantile behavior is analogous to
        # Series' quantile for each interpolation option
        s = Series(data)

        q1 = s.quantile(quantile, interpolation)
        q2 = s.expanding(min_periods=1).quantile(quantile, interpolation).iloc[-1]

        if np.isnan(q1):
            assert np.isnan(q2)
        else:
            assert q1 == q2

    def test_invalid_quantile_value(self):
        data = np.arange(5)
        s = Series(data)

        msg = "Interpolation 'invalid' is not supported"
        with pytest.raises(ValueError, match=msg):
            s.rolling(len(data), min_periods=1).quantile(0.5, interpolation="invalid")

    def test_rolling_quantile_param(self):
        ser = Series([0.0, 0.1, 0.5, 0.9, 1.0])

        with pytest.raises(ValueError):
            ser.rolling(3).quantile(-0.1)

        with pytest.raises(ValueError):
            ser.rolling(3).quantile(10.0)

        with pytest.raises(TypeError):
            ser.rolling(3).quantile("foo")

    def test_rolling_apply(self, raw):
        # suppress warnings about empty slices, as we are deliberately testing
        # with a 0-length Series

        def f(x):
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message=".*(empty slice|0 for slice).*",
                    category=RuntimeWarning,
                )
                return x[np.isfinite(x)].mean()

        self._check_moment_func(np.mean, name="apply", func=f, raw=raw)

    def test_rolling_std(self, raw):
        self._check_moment_func(lambda x: np.std(x, ddof=1), name="std", raw=raw)
        self._check_moment_func(
            lambda x: np.std(x, ddof=0), name="std", ddof=0, raw=raw
        )

    def test_rolling_std_1obs(self):
        vals = pd.Series([1.0, 2.0, 3.0, 4.0, 5.0])

        result = vals.rolling(1, min_periods=1).std()
        expected = pd.Series([np.nan] * 5)
        tm.assert_series_equal(result, expected)

        result = vals.rolling(1, min_periods=1).std(ddof=0)
        expected = pd.Series([0.0] * 5)
        tm.assert_series_equal(result, expected)

        result = pd.Series([np.nan, np.nan, 3, 4, 5]).rolling(3, min_periods=2).std()
        assert np.isnan(result[2])

    def test_rolling_std_neg_sqrt(self):
        # unit test from Bottleneck

        # Test move_nanstd for neg sqrt.

        a = pd.Series(
            [
                0.0011448196318903589,
                0.00028718669878572767,
                0.00028718669878572767,
                0.00028718669878572767,
                0.00028718669878572767,
            ]
        )
        b = a.rolling(window=3).std()
        assert np.isfinite(b[2:]).all()

        b = a.ewm(span=3).std()
        assert np.isfinite(b[2:]).all()

    def test_rolling_var(self, raw):
        self._check_moment_func(lambda x: np.var(x, ddof=1), name="var", raw=raw)
        self._check_moment_func(
            lambda x: np.var(x, ddof=0), name="var", ddof=0, raw=raw
        )

    @td.skip_if_no_scipy
    def test_rolling_skew(self, raw):
        from scipy.stats import skew

        self._check_moment_func(lambda x: skew(x, bias=False), name="skew", raw=raw)

    @td.skip_if_no_scipy
    def test_rolling_kurt(self, raw):
        from scipy.stats import kurtosis

        self._check_moment_func(lambda x: kurtosis(x, bias=False), name="kurt", raw=raw)

    def _check_moment_func(
        self,
        static_comp,
        name,
        raw,
        has_min_periods=True,
        has_center=True,
        has_time_rule=True,
        fill_value=None,
        zero_min_periods_equal=True,
        **kwargs,
    ):

        # inject raw
        if name == "apply":
            kwargs = copy.copy(kwargs)
            kwargs["raw"] = raw

        def get_result(obj, window, min_periods=None, center=False):
            r = obj.rolling(window=window, min_periods=min_periods, center=center)
            return getattr(r, name)(**kwargs)

        series_result = get_result(self.series, window=50)
        assert isinstance(series_result, Series)
        tm.assert_almost_equal(series_result.iloc[-1], static_comp(self.series[-50:]))

        frame_result = get_result(self.frame, window=50)
        assert isinstance(frame_result, DataFrame)
        tm.assert_series_equal(
            frame_result.iloc[-1, :],
            self.frame.iloc[-50:, :].apply(static_comp, axis=0, raw=raw),
            check_names=False,
        )

        # check time_rule works
        if has_time_rule:
            win = 25
            minp = 10
            series = self.series[::2].resample("B").mean()
            frame = self.frame[::2].resample("B").mean()

            if has_min_periods:
                series_result = get_result(series, window=win, min_periods=minp)
                frame_result = get_result(frame, window=win, min_periods=minp)
            else:
                series_result = get_result(series, window=win)
                frame_result = get_result(frame, window=win)

            last_date = series_result.index[-1]
            prev_date = last_date - 24 * offsets.BDay()

            trunc_series = self.series[::2].truncate(prev_date, last_date)
            trunc_frame = self.frame[::2].truncate(prev_date, last_date)

            tm.assert_almost_equal(series_result[-1], static_comp(trunc_series))

            tm.assert_series_equal(
                frame_result.xs(last_date),
                trunc_frame.apply(static_comp, raw=raw),
                check_names=False,
            )

        # excluding NaNs correctly
        obj = Series(randn(50))
        obj[:10] = np.NaN
        obj[-10:] = np.NaN
        if has_min_periods:
            result = get_result(obj, 50, min_periods=30)
            tm.assert_almost_equal(result.iloc[-1], static_comp(obj[10:-10]))

            # min_periods is working correctly
            result = get_result(obj, 20, min_periods=15)
            assert isna(result.iloc[23])
            assert not isna(result.iloc[24])

            assert not isna(result.iloc[-6])
            assert isna(result.iloc[-5])

            obj2 = Series(randn(20))
            result = get_result(obj2, 10, min_periods=5)
            assert isna(result.iloc[3])
            assert notna(result.iloc[4])

            if zero_min_periods_equal:
                # min_periods=0 may be equivalent to min_periods=1
                result0 = get_result(obj, 20, min_periods=0)
                result1 = get_result(obj, 20, min_periods=1)
                tm.assert_almost_equal(result0, result1)
        else:
            result = get_result(obj, 50)
            tm.assert_almost_equal(result.iloc[-1], static_comp(obj[10:-10]))

        # window larger than series length (#7297)
        if has_min_periods:
            for minp in (0, len(self.series) - 1, len(self.series)):
                result = get_result(self.series, len(self.series) + 1, min_periods=minp)
                expected = get_result(self.series, len(self.series), min_periods=minp)
                nan_mask = isna(result)
                tm.assert_series_equal(nan_mask, isna(expected))

                nan_mask = ~nan_mask
                tm.assert_almost_equal(result[nan_mask], expected[nan_mask])
        else:
            result = get_result(self.series, len(self.series) + 1)
            expected = get_result(self.series, len(self.series))
            nan_mask = isna(result)
            tm.assert_series_equal(nan_mask, isna(expected))

            nan_mask = ~nan_mask
            tm.assert_almost_equal(result[nan_mask], expected[nan_mask])

        # check center=True
        if has_center:
            if has_min_periods:
                result = get_result(obj, 20, min_periods=15, center=True)
                expected = get_result(
                    pd.concat([obj, Series([np.NaN] * 9)]), 20, min_periods=15
                )[9:].reset_index(drop=True)
            else:
                result = get_result(obj, 20, center=True)
                expected = get_result(pd.concat([obj, Series([np.NaN] * 9)]), 20)[
                    9:
                ].reset_index(drop=True)

            tm.assert_series_equal(result, expected)

            # shifter index
            s = ["x{x:d}".format(x=x) for x in range(12)]

            if has_min_periods:
                minp = 10

                series_xp = (
                    get_result(
                        self.series.reindex(list(self.series.index) + s),
                        window=25,
                        min_periods=minp,
                    )
                    .shift(-12)
                    .reindex(self.series.index)
                )
                frame_xp = (
                    get_result(
                        self.frame.reindex(list(self.frame.index) + s),
                        window=25,
                        min_periods=minp,
                    )
                    .shift(-12)
                    .reindex(self.frame.index)
                )

                series_rs = get_result(
                    self.series, window=25, min_periods=minp, center=True
                )
                frame_rs = get_result(
                    self.frame, window=25, min_periods=minp, center=True
                )

            else:
                series_xp = (
                    get_result(
                        self.series.reindex(list(self.series.index) + s), window=25
                    )
                    .shift(-12)
                    .reindex(self.series.index)
                )
                frame_xp = (
                    get_result(
                        self.frame.reindex(list(self.frame.index) + s), window=25
                    )
                    .shift(-12)
                    .reindex(self.frame.index)
                )

                series_rs = get_result(self.series, window=25, center=True)
                frame_rs = get_result(self.frame, window=25, center=True)

            if fill_value is not None:
                series_xp = series_xp.fillna(fill_value)
                frame_xp = frame_xp.fillna(fill_value)
            tm.assert_series_equal(series_xp, series_rs)
            tm.assert_frame_equal(frame_xp, frame_rs)


def _rolling_consistency_cases():
    for window in [1, 2, 3, 10, 20]:
        for min_periods in {0, 1, 2, 3, 4, window}:
            if min_periods and (min_periods > window):
                continue
            for center in [False, True]:
                yield window, min_periods, center


class TestRollingMomentsConsistency(ConsistencyBase):
    def setup_method(self, method):
        self._create_data()

    @pytest.mark.slow
    @pytest.mark.parametrize(
        "window,min_periods,center", list(_rolling_consistency_cases())
    )
    def test_rolling_consistency(self, window, min_periods, center):

        # suppress warnings about empty slices, as we are deliberately testing
        # with empty/0-length Series/DataFrames
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=".*(empty slice|0 for slice).*",
                category=RuntimeWarning,
            )

            # test consistency between different rolling_* moments
            self._test_moments_consistency_mock_mean(
                mean=lambda x: (
                    x.rolling(
                        window=window, min_periods=min_periods, center=center
                    ).mean()
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

            self._test_moments_consistency_is_constant(
                min_periods=min_periods,
                count=lambda x: (x.rolling(window=window, center=center).count()),
                mean=lambda x: (
                    x.rolling(
                        window=window, min_periods=min_periods, center=center
                    ).mean()
                ),
                corr=lambda x, y: (
                    x.rolling(
                        window=window, min_periods=min_periods, center=center
                    ).corr(y)
                ),
            )

            self._test_moments_consistency_var_debiasing_factors(
                var_unbiased=lambda x: (
                    x.rolling(
                        window=window, min_periods=min_periods, center=center
                    ).var()
                ),
                var_biased=lambda x: (
                    x.rolling(
                        window=window, min_periods=min_periods, center=center
                    ).var(ddof=0)
                ),
                var_debiasing_factors=lambda x: (
                    x.rolling(window=window, center=center)
                    .count()
                    .divide(
                        (x.rolling(window=window, center=center).count() - 1.0).replace(
                            0.0, np.nan
                        )
                    )
                ),
            )

            self._test_moments_consistency(
                min_periods=min_periods,
                count=lambda x: (x.rolling(window=window, center=center).count()),
                mean=lambda x: (
                    x.rolling(
                        window=window, min_periods=min_periods, center=center
                    ).mean()
                ),
                corr=lambda x, y: (
                    x.rolling(
                        window=window, min_periods=min_periods, center=center
                    ).corr(y)
                ),
                var_unbiased=lambda x: (
                    x.rolling(
                        window=window, min_periods=min_periods, center=center
                    ).var()
                ),
                std_unbiased=lambda x: (
                    x.rolling(
                        window=window, min_periods=min_periods, center=center
                    ).std()
                ),
                cov_unbiased=lambda x, y: (
                    x.rolling(
                        window=window, min_periods=min_periods, center=center
                    ).cov(y)
                ),
                var_biased=lambda x: (
                    x.rolling(
                        window=window, min_periods=min_periods, center=center
                    ).var(ddof=0)
                ),
                std_biased=lambda x: (
                    x.rolling(
                        window=window, min_periods=min_periods, center=center
                    ).std(ddof=0)
                ),
                cov_biased=lambda x, y: (
                    x.rolling(
                        window=window, min_periods=min_periods, center=center
                    ).cov(y, ddof=0)
                ),
            )

            # test consistency between rolling_xyz() and either (a)
            # rolling_apply of Series.xyz(), or (b) rolling_apply of
            # np.nanxyz()
            for (x, is_constant, no_nans) in self.data:
                functions = self.base_functions

                # GH 8269
                if no_nans:
                    functions = self.base_functions + self.no_nan_functions
                for (f, require_min_periods, name) in functions:
                    rolling_f = getattr(
                        x.rolling(
                            window=window, center=center, min_periods=min_periods
                        ),
                        name,
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
                            window=window, min_periods=0, center=center
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

    # binary moments
    def test_rolling_cov(self):
        A = self.series
        B = A + randn(len(A))

        result = A.rolling(window=50, min_periods=25).cov(B)
        tm.assert_almost_equal(result[-1], np.cov(A[-50:], B[-50:])[0, 1])

    def test_rolling_cov_pairwise(self):
        self._check_pairwise_moment("rolling", "cov", window=10, min_periods=5)

    def test_rolling_corr(self):
        A = self.series
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

    def test_rolling_corr_pairwise(self):
        self._check_pairwise_moment("rolling", "corr", window=10, min_periods=5)

    @pytest.mark.parametrize("window", range(7))
    def test_rolling_corr_with_zero_variance(self, window):
        # GH 18430
        s = pd.Series(np.zeros(20))
        other = pd.Series(np.arange(20))

        assert s.rolling(window=window).corr(other=other).isna().all()

    def test_flex_binary_moment(self):
        # GH3155
        # don't blow the stack
        msg = (
            "arguments to moment function must be of type"
            " np.ndarray/Series/DataFrame"
        )
        with pytest.raises(TypeError, match=msg):
            _flex_binary_moment(5, 6, None)

    def test_corr_sanity(self):
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

    @pytest.mark.parametrize("method", ["corr", "cov"])
    def test_flex_binary_frame(self, method):
        series = self.frame[1]

        res = getattr(series.rolling(window=10), method)(self.frame)
        res2 = getattr(self.frame.rolling(window=10), method)(series)
        exp = self.frame.apply(lambda x: getattr(series.rolling(window=10), method)(x))

        tm.assert_frame_equal(res, exp)
        tm.assert_frame_equal(res2, exp)

        frame2 = self.frame.copy()
        frame2.values[:] = np.random.randn(*frame2.shape)

        res3 = getattr(self.frame.rolling(window=10), method)(frame2)
        exp = DataFrame(
            {
                k: getattr(self.frame[k].rolling(window=10), method)(frame2[k])
                for k in self.frame
            }
        )
        tm.assert_frame_equal(res3, exp)

    def test_rolling_cov_diff_length(self):
        # GH 7512
        s1 = Series([1, 2, 3], index=[0, 1, 2])
        s2 = Series([1, 3], index=[0, 2])
        result = s1.rolling(window=3, min_periods=2).cov(s2)
        expected = Series([None, None, 2.0])
        tm.assert_series_equal(result, expected)

        s2a = Series([1, None, 3], index=[0, 1, 2])
        result = s1.rolling(window=3, min_periods=2).cov(s2a)
        tm.assert_series_equal(result, expected)

    def test_rolling_corr_diff_length(self):
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
    def test_rolling_functions_window_non_shrinkage(self, f):
        # GH 7764
        s = Series(range(4))
        s_expected = Series(np.nan, index=s.index)
        df = DataFrame([[1, 5], [3, 2], [3, 9], [-1, 0]], columns=["A", "B"])
        df_expected = DataFrame(np.nan, index=df.index, columns=df.columns)

        s_result = f(s)
        tm.assert_series_equal(s_result, s_expected)

        df_result = f(df)
        tm.assert_frame_equal(df_result, df_expected)

    def test_rolling_functions_window_non_shrinkage_binary(self):

        # corr/cov return a MI DataFrame
        df = DataFrame(
            [[1, 5], [3, 2], [3, 9], [-1, 0]],
            columns=Index(["A", "B"], name="foo"),
            index=Index(range(4), name="bar"),
        )
        df_expected = DataFrame(
            columns=Index(["A", "B"], name="foo"),
            index=pd.MultiIndex.from_product(
                [df.index, df.columns], names=["bar", "foo"]
            ),
            dtype="float64",
        )
        functions = [
            lambda x: (x.rolling(window=10, min_periods=5).cov(x, pairwise=True)),
            lambda x: (x.rolling(window=10, min_periods=5).corr(x, pairwise=True)),
        ]
        for f in functions:
            df_result = f(df)
            tm.assert_frame_equal(df_result, df_expected)

    def test_rolling_skew_edge_cases(self):

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

    def test_rolling_kurt_edge_cases(self):

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

    def test_rolling_skew_eq_value_fperr(self):
        # #18804 all rolling skew for all equal values should return Nan
        a = Series([1.1] * 15).rolling(window=10).skew()
        assert np.isnan(a).all()

    def test_rolling_kurt_eq_value_fperr(self):
        # #18804 all rolling kurt for all equal values should return Nan
        a = Series([1.1] * 15).rolling(window=10).kurt()
        assert np.isnan(a).all()

    def test_rolling_max_gh6297(self):
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
            index=[datetime(1975, 1, i, 0) for i in range(1, 6)],
        )
        x = series.resample("D").max().rolling(window=1).max()
        tm.assert_series_equal(expected, x)

    def test_rolling_max_resample(self):

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
            index=[datetime(1975, 1, i, 0) for i in range(1, 6)],
        )
        x = series.resample("D").max().rolling(window=1).max()
        tm.assert_series_equal(expected, x)

        # Now specify median (10.0)
        expected = Series(
            [0.0, 1.0, 2.0, 3.0, 10.0],
            index=[datetime(1975, 1, i, 0) for i in range(1, 6)],
        )
        x = series.resample("D").median().rolling(window=1).max()
        tm.assert_series_equal(expected, x)

        # Now specify mean (4+10+20)/3
        v = (4.0 + 10.0 + 20.0) / 3.0
        expected = Series(
            [0.0, 1.0, 2.0, 3.0, v],
            index=[datetime(1975, 1, i, 0) for i in range(1, 6)],
        )
        x = series.resample("D").mean().rolling(window=1).max()
        tm.assert_series_equal(expected, x)

    def test_rolling_min_resample(self):

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
            index=[datetime(1975, 1, i, 0) for i in range(1, 6)],
        )
        r = series.resample("D").min().rolling(window=1)
        tm.assert_series_equal(expected, r.min())

    def test_rolling_median_resample(self):

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
            index=[datetime(1975, 1, i, 0) for i in range(1, 6)],
        )
        x = series.resample("D").median().rolling(window=1).median()
        tm.assert_series_equal(expected, x)

    def test_rolling_median_memory_error(self):
        # GH11722
        n = 20000
        Series(np.random.randn(n)).rolling(window=2, center=False).median()
        Series(np.random.randn(n)).rolling(window=2, center=False).median()

    def test_rolling_min_max_numeric_types(self):

        # GH12373
        types_test = [np.dtype("f{}".format(width)) for width in [4, 8]]
        types_test.extend(
            [
                np.dtype("{}{}".format(sign, width))
                for width in [1, 2, 4, 8]
                for sign in "ui"
            ]
        )
        for data_type in types_test:
            # Just testing that these don't throw exceptions and that
            # the return type is float64. Other tests will cover quantitative
            # correctness
            result = DataFrame(np.arange(20, dtype=data_type)).rolling(window=5).max()
            assert result.dtypes[0] == np.dtype("f8")
            result = DataFrame(np.arange(20, dtype=data_type)).rolling(window=5).min()
            assert result.dtypes[0] == np.dtype("f8")

    def test_moment_functions_zero_length(self):
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

    def test_moment_functions_zero_length_pairwise(self):

        df1 = DataFrame()
        df2 = DataFrame(columns=Index(["a"], name="foo"), index=Index([], name="bar"))
        df2["a"] = df2["a"].astype("float64")

        df1_expected = DataFrame(
            index=pd.MultiIndex.from_product([df1.index, df1.columns]),
            columns=Index([]),
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
