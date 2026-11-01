from fractions import Fraction

import numpy as np
import pytest

from pandas.core.dtypes.common import is_integer

import pandas as pd
import pandas._testing as tm
from pandas.core.indexes.datetimes import Timestamp


class TestSeriesQuantile:
    def test_quantile(self, datetime_series):
        q = datetime_series.quantile(0.1)
        assert q == np.percentile(datetime_series.dropna(), 10)

        q = datetime_series.quantile(0.9)
        assert q == np.percentile(datetime_series.dropna(), 90)

        # object dtype
        q = pd.Series(datetime_series, dtype=object).quantile(0.9)
        assert q == np.percentile(datetime_series.dropna(), 90)

        # datetime64[ns] dtype
        dts = datetime_series.index.to_series()
        q = dts.quantile(0.2)
        assert q == Timestamp("2000-01-10 19:12:00")

        # timedelta64[ns] dtype
        tds = dts.diff()
        q = tds.quantile(0.25)
        assert q == pd.to_timedelta("24:00:00")

        # GH7661
        result = pd.Series([np.timedelta64("NaT", "ns")]).sum()
        assert result == pd.Timedelta(0)

        msg = "percentiles should all be in the interval \\[0, 1\\]"
        for invalid in [-1, 2, [0.5, -1], [0.5, 2]]:
            with pytest.raises(ValueError, match=msg):
                datetime_series.quantile(invalid)

        s = pd.Series(np.random.default_rng(2).standard_normal(100))
        percentile_array = [-0.5, 0.25, 1.5]
        with pytest.raises(ValueError, match=msg):
            s.quantile(percentile_array)

    def test_quantile_multi(self, datetime_series, unit):
        datetime_series.index = datetime_series.index.as_unit(unit)
        qs = [0.1, 0.9]
        result = datetime_series.quantile(qs)
        expected = pd.Series(
            [
                np.percentile(datetime_series.dropna(), 10),
                np.percentile(datetime_series.dropna(), 90),
            ],
            index=qs,
            name=datetime_series.name,
        )
        tm.assert_series_equal(result, expected)

        dts = datetime_series.index.to_series()
        dts.name = "xxx"
        result = dts.quantile((0.2, 0.2))
        expected = pd.Series(
            [Timestamp("2000-01-10 19:12:00"), Timestamp("2000-01-10 19:12:00")],
            index=[0.2, 0.2],
            name="xxx",
            dtype=f"M8[{unit}]",
        )
        tm.assert_series_equal(result, expected)

        result = datetime_series.quantile([])
        expected = pd.Series(
            [],
            name=datetime_series.name,
            index=pd.Index([], dtype=float),
            dtype="float64",
        )
        tm.assert_series_equal(result, expected)

    def test_quantile_interpolation(self, datetime_series):
        # see gh-10174

        # interpolation = linear (default case)
        q = datetime_series.quantile(0.1, interpolation="linear")
        assert q == np.percentile(datetime_series.dropna(), 10)
        q1 = datetime_series.quantile(0.1)
        assert q1 == np.percentile(datetime_series.dropna(), 10)

        # test with and without interpolation keyword
        assert q == q1

    def test_quantile_interpolation_dtype(self):
        # GH #10174

        # interpolation = linear (default case)
        q = pd.Series([1, 3, 4]).quantile(0.5, interpolation="lower")
        assert q == np.percentile(np.array([1, 3, 4]), 50)
        assert is_integer(q)

        q = pd.Series([1, 3, 4]).quantile(0.5, interpolation="higher")
        assert q == np.percentile(np.array([1, 3, 4]), 50)
        assert is_integer(q)

    def test_quantile_object_dtype_preserves_numpy_scalars(self):
        # GH#64266
        ser = pd.Series([np.int8(1), np.int8(3)], dtype=object)
        result = ser.quantile(0.5, interpolation="lower")
        assert result == 1
        assert type(result) is np.int8

    def test_quantile_nan(self):
        # GH 13098
        ser = pd.Series([1, 2, 3, 4, np.nan])
        result = ser.quantile(0.5)
        expected = 2.5
        assert result == expected

        # all nan/empty
        s1 = pd.Series([], dtype=object)
        cases = [s1, pd.Series([np.nan, np.nan])]

        for ser in cases:
            res = ser.quantile(0.5)
            assert np.isnan(res)

            res = ser.quantile([0.5])
            tm.assert_series_equal(res, pd.Series([np.nan], index=[0.5]))

            res = ser.quantile([0.2, 0.3])
            tm.assert_series_equal(res, pd.Series([np.nan, np.nan], index=[0.2, 0.3]))

    @pytest.mark.parametrize(
        "case",
        [
            [
                Timestamp("2011-01-01"),
                Timestamp("2011-01-02"),
                Timestamp("2011-01-03"),
            ],
            [
                Timestamp("2011-01-01", tz="US/Eastern"),
                Timestamp("2011-01-02", tz="US/Eastern"),
                Timestamp("2011-01-03", tz="US/Eastern"),
            ],
            [pd.Timedelta("1 days"), pd.Timedelta("2 days"), pd.Timedelta("3 days")],
            # NaT
            [
                Timestamp("2011-01-01"),
                Timestamp("2011-01-02"),
                Timestamp("2011-01-03"),
                pd.NaT,
            ],
            [
                Timestamp("2011-01-01", tz="US/Eastern"),
                Timestamp("2011-01-02", tz="US/Eastern"),
                Timestamp("2011-01-03", tz="US/Eastern"),
                pd.NaT,
            ],
            [
                pd.Timedelta("1 days"),
                pd.Timedelta("2 days"),
                pd.Timedelta("3 days"),
                pd.NaT,
            ],
        ],
    )
    def test_quantile_box(self, case):
        ser = pd.Series(case, name="XXX")
        res = ser.quantile(0.5)
        assert res == case[1]

        res = ser.quantile([0.5])
        exp = pd.Series([case[1]], index=[0.5], name="XXX")
        tm.assert_series_equal(res, exp)

    def test_datetime_timedelta_quantiles(self):
        # covers #9694
        assert pd.isna(pd.Series([], dtype="M8[ns]").quantile(0.5))
        assert pd.isna(pd.Series([], dtype="m8[ns]").quantile(0.5))

    def test_quantile_nat(self):
        res = pd.Series([pd.NaT, pd.NaT]).quantile(0.5)
        assert res is pd.NaT

        res = pd.Series([pd.NaT, pd.NaT]).quantile([0.5])
        tm.assert_series_equal(res, pd.Series([pd.NaT], index=[0.5]))

    @pytest.mark.parametrize(
        "values, dtype",
        [([0, 0, 0, 1, 2, 3], "Sparse[int]"), ([0.0, None, 1.0, 2.0], "Sparse[float]")],
    )
    def test_quantile_sparse(self, values, dtype):
        ser = pd.Series(values, dtype=dtype)
        result = ser.quantile([0.5])
        expected = pd.Series(np.asarray(ser)).quantile([0.5]).astype("Sparse[float]")
        tm.assert_series_equal(result, expected)

    def test_quantile_empty_float64(self):
        # floats
        ser = pd.Series([], dtype="float64")

        res = ser.quantile(0.5)
        assert np.isnan(res)

        res = ser.quantile([0.5])
        exp = pd.Series([np.nan], index=[0.5])
        tm.assert_series_equal(res, exp)

    def test_quantile_empty_int64(self):
        # int
        ser = pd.Series([], dtype="int64")

        res = ser.quantile(0.5)
        assert np.isnan(res)

        res = ser.quantile([0.5])
        exp = pd.Series([np.nan], index=[0.5])
        tm.assert_series_equal(res, exp)

    def test_quantile_empty_dt64(self):
        # datetime
        ser = pd.Series([], dtype="datetime64[ns]")

        res = ser.quantile(0.5)
        assert res is pd.NaT

        res = ser.quantile([0.5])
        exp = pd.Series([pd.NaT], index=[0.5], dtype=ser.dtype)
        tm.assert_series_equal(res, exp)

    @pytest.mark.parametrize("dtype", [int, float, "Int64"])
    def test_quantile_dtypes(self, dtype):
        result = pd.Series([1, 2, 3], dtype=dtype).quantile(np.arange(0, 1, 0.25))
        expected = pd.Series(np.arange(1, 3, 0.5), index=np.arange(0, 1, 0.25))
        if dtype == "Int64":
            expected = expected.astype("Float64")
        tm.assert_series_equal(result, expected)

    def test_quantile_all_na(self, any_int_ea_dtype):
        # GH#50681
        ser = pd.Series([pd.NA, pd.NA], dtype=any_int_ea_dtype)
        with tm.assert_produces_warning(None):
            result = ser.quantile([0.1, 0.5])
        expected = pd.Series([pd.NA, pd.NA], dtype=any_int_ea_dtype, index=[0.1, 0.5])
        tm.assert_series_equal(result, expected)

    def test_quantile_dtype_size(self, any_int_ea_dtype):
        # GH#50681
        ser = pd.Series([pd.NA, pd.NA, 1], dtype=any_int_ea_dtype)
        result = ser.quantile([0.1, 0.5])
        expected = pd.Series([1, 1], dtype=any_int_ea_dtype, index=[0.1, 0.5])
        tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("typ", ["datetime64", "timedelta64"])
def test_quantile_empty_datetimelike(typ, unit):
    ser = pd.Series([], dtype=f"{typ}[{unit}]")
    result = ser.quantile()
    assert result is pd.NaT


@pytest.mark.parametrize("typ", ["datetime64", "timedelta64"])
@pytest.mark.parametrize(
    "interpolation", ["linear", "lower", "higher", "midpoint", "nearest"]
)
def test_quantile_datetimelike_endpoints_exact(typ, interpolation):
    # GH#49110
    values = np.array([1567703790155681000, 1567703793155681000], dtype="i8")
    ser = pd.Series(values.view(f"{typ}[ns]"))
    result = ser.quantile([0, 1], interpolation=interpolation)
    expected = pd.Series([ser.min(), ser.max()], index=[0.0, 1.0])
    tm.assert_series_equal(result, expected)


def test_quantile_dt64_linear_exact():
    # GH#49110
    ser = pd.Series(
        np.array([1567703790155681000, 1567703793155681001], dtype="i8").view("M8[ns]")
    )
    result = ser.quantile(0.5)
    assert result.value == 1567703791655681000


def test_quantile_dt64_fraction_truncates_toward_zero():
    # GH#49110 fractional nanoseconds truncate toward zero
    ser = pd.Series(np.array([0, 1], dtype="i8").view("M8[ns]"))
    assert ser.quantile(0.5).value == 0
    ser = pd.Series(np.array([-1, 0], dtype="i8").view("M8[ns]"))
    assert ser.quantile(0.5).value == 0
    ser = pd.Series(np.array([-3, -2], dtype="i8").view("M8[ns]"))
    assert ser.quantile(0.5).value == -2


def test_quantile_dt64_with_nat_exact():
    # GH#49110
    values = np.array([1567703790155681000, 1567703793155681000], dtype="i8")
    ser = pd.Series([*values.view("M8[ns]"), pd.NaT], dtype="M8[ns]")
    result = ser.quantile([0, 1])
    expected = pd.Series(values.view("M8[ns]"), index=[0.0, 1.0])
    tm.assert_series_equal(result, expected)


def test_quantile_dt64_full_int64_range_no_overflow():
    # GH#49110 span wider than int64 must not overflow
    lo = -(2**62) - 2**61
    hi = 2**62 + 2**61
    ser = pd.Series(np.array([lo, hi], dtype="i8").view("M8[ns]"))
    result = ser.quantile([0, 0.5, 1])
    expected = pd.Series(
        np.array([lo, 0, hi], dtype="i8").view("M8[ns]"), index=[0.0, 0.5, 1.0]
    )
    tm.assert_series_equal(result, expected)


def test_quantile_dt64_far_apart_neighbors():
    # GH#49110 neighbors more than 2**53 apart: result stays between them,
    #  within the documented error of gap * 2**-52
    lo, hi = -(2**62), 2**62 + 12345
    ser = pd.Series(np.array([lo, hi], dtype="i8").view("M8[ns]"))
    qs = [0.1, 0.5, 0.9, np.nextafter(1.0, 0.0)]
    result = ser.quantile(qs).values.view("i8")
    for q, res in zip(qs, result, strict=True):
        exact = lo + int((hi - lo) * Fraction(q))
        assert lo <= res <= hi
        assert abs(int(res) - exact) <= (hi - lo) >> 52


def test_quantile_dt64tz_endpoints_exact():
    # GH#49110
    dti = pd.DatetimeIndex(
        ["2019-09-05 17:16:30.155681", "2019-09-05 17:16:33.155681"], tz="UTC"
    ).as_unit("ns")
    ser = pd.Series(dti)
    result = ser.quantile([0, 1])
    expected = pd.Series(dti, index=[0.0, 1.0])
    tm.assert_series_equal(result, expected)
