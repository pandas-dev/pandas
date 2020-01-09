import numpy as np
import pytest

from pandas.errors import NullFrequencyError

import pandas as pd
from pandas import (
    DatetimeIndex,
    Index,
    NaT,
    Series,
    TimedeltaIndex,
    date_range,
    offsets,
)
import pandas._testing as tm

from pandas.tseries.offsets import BDay


class TestShift:
    def test_shift(self, datetime_series):
        shifted = datetime_series.shift(1)
        unshifted = shifted.shift(-1)

        tm.assert_index_equal(shifted.index, datetime_series.index)
        tm.assert_index_equal(unshifted.index, datetime_series.index)
        tm.assert_numpy_array_equal(
            unshifted.dropna().values, datetime_series.values[:-1]
        )

        offset = BDay()
        shifted = datetime_series.shift(1, freq=offset)
        unshifted = shifted.shift(-1, freq=offset)

        tm.assert_series_equal(unshifted, datetime_series)

        unshifted = datetime_series.shift(0, freq=offset)
        tm.assert_series_equal(unshifted, datetime_series)

        shifted = datetime_series.shift(1, freq="B")
        unshifted = shifted.shift(-1, freq="B")

        tm.assert_series_equal(unshifted, datetime_series)

        # corner case
        unshifted = datetime_series.shift(0)
        tm.assert_series_equal(unshifted, datetime_series)

        # Shifting with PeriodIndex
        ps = tm.makePeriodSeries()
        shifted = ps.shift(1)
        unshifted = shifted.shift(-1)
        tm.assert_index_equal(shifted.index, ps.index)
        tm.assert_index_equal(unshifted.index, ps.index)
        tm.assert_numpy_array_equal(unshifted.dropna().values, ps.values[:-1])

        shifted2 = ps.shift(1, "B")
        shifted3 = ps.shift(1, BDay())
        tm.assert_series_equal(shifted2, shifted3)
        tm.assert_series_equal(ps, shifted2.shift(-1, "B"))

        msg = "Given freq D does not match PeriodIndex freq B"
        with pytest.raises(ValueError, match=msg):
            ps.shift(freq="D")

        # legacy support
        shifted4 = ps.shift(1, freq="B")
        tm.assert_series_equal(shifted2, shifted4)

        shifted5 = ps.shift(1, freq=BDay())
        tm.assert_series_equal(shifted5, shifted4)

        # 32-bit taking
        # GH#8129
        index = date_range("2000-01-01", periods=5)
        for dtype in ["int32", "int64"]:
            s1 = Series(np.arange(5, dtype=dtype), index=index)
            p = s1.iloc[1]
            result = s1.shift(periods=p)
            expected = Series([np.nan, 0, 1, 2, 3], index=index)
            tm.assert_series_equal(result, expected)

        # GH#8260
        # with tz
        s = Series(
            date_range("2000-01-01 09:00:00", periods=5, tz="US/Eastern"), name="foo"
        )
        result = s - s.shift()

        exp = Series(TimedeltaIndex(["NaT"] + ["1 days"] * 4), name="foo")
        tm.assert_series_equal(result, exp)

        # incompat tz
        s2 = Series(date_range("2000-01-01 09:00:00", periods=5, tz="CET"), name="foo")
        msg = "DatetimeArray subtraction must have the same timezones or no timezones"
        with pytest.raises(TypeError, match=msg):
            s - s2

    def test_shift2(self):
        ts = Series(
            np.random.randn(5), index=date_range("1/1/2000", periods=5, freq="H")
        )

        result = ts.shift(1, freq="5T")
        exp_index = ts.index.shift(1, freq="5T")
        tm.assert_index_equal(result.index, exp_index)

        # GH#1063, multiple of same base
        result = ts.shift(1, freq="4H")
        exp_index = ts.index + offsets.Hour(4)
        tm.assert_index_equal(result.index, exp_index)

        idx = DatetimeIndex(["2000-01-01", "2000-01-02", "2000-01-04"])
        msg = "Cannot shift with no freq"
        with pytest.raises(NullFrequencyError, match=msg):
            idx.shift(1)

    def test_shift_fill_value(self):
        # GH#24128
        ts = Series(
            [1.0, 2.0, 3.0, 4.0, 5.0], index=date_range("1/1/2000", periods=5, freq="H")
        )

        exp = Series(
            [0.0, 1.0, 2.0, 3.0, 4.0], index=date_range("1/1/2000", periods=5, freq="H")
        )
        # check that fill value works
        result = ts.shift(1, fill_value=0.0)
        tm.assert_series_equal(result, exp)

        exp = Series(
            [0.0, 0.0, 1.0, 2.0, 3.0], index=date_range("1/1/2000", periods=5, freq="H")
        )
        result = ts.shift(2, fill_value=0.0)
        tm.assert_series_equal(result, exp)

        ts = pd.Series([1, 2, 3])
        res = ts.shift(2, fill_value=0)
        assert res.dtype == ts.dtype

    def test_shift_categorical_fill_value(self):
        ts = pd.Series(["a", "b", "c", "d"], dtype="category")
        res = ts.shift(1, fill_value="a")
        expected = pd.Series(
            pd.Categorical(
                ["a", "a", "b", "c"], categories=["a", "b", "c", "d"], ordered=False
            )
        )
        tm.assert_equal(res, expected)

        # check for incorrect fill_value
        msg = "'fill_value=f' is not present in this Categorical's categories"
        with pytest.raises(ValueError, match=msg):
            ts.shift(1, fill_value="f")

    def test_shift_dst(self):
        # GH#13926
        dates = date_range("2016-11-06", freq="H", periods=10, tz="US/Eastern")
        s = Series(dates)

        res = s.shift(0)
        tm.assert_series_equal(res, s)
        assert res.dtype == "datetime64[ns, US/Eastern]"

        res = s.shift(1)
        exp_vals = [NaT] + dates.astype(object).values.tolist()[:9]
        exp = Series(exp_vals)
        tm.assert_series_equal(res, exp)
        assert res.dtype == "datetime64[ns, US/Eastern]"

        res = s.shift(-2)
        exp_vals = dates.astype(object).values.tolist()[2:] + [NaT, NaT]
        exp = Series(exp_vals)
        tm.assert_series_equal(res, exp)
        assert res.dtype == "datetime64[ns, US/Eastern]"

        for ex in [10, -10, 20, -20]:
            res = s.shift(ex)
            exp = Series([NaT] * 10, dtype="datetime64[ns, US/Eastern]")
            tm.assert_series_equal(res, exp)
            assert res.dtype == "datetime64[ns, US/Eastern]"

    def test_tshift(self, datetime_series):
        # PeriodIndex
        ps = tm.makePeriodSeries()
        shifted = ps.tshift(1)
        unshifted = shifted.tshift(-1)

        tm.assert_series_equal(unshifted, ps)

        shifted2 = ps.tshift(freq="B")
        tm.assert_series_equal(shifted, shifted2)

        shifted3 = ps.tshift(freq=BDay())
        tm.assert_series_equal(shifted, shifted3)

        msg = "Given freq M does not match PeriodIndex freq B"
        with pytest.raises(ValueError, match=msg):
            ps.tshift(freq="M")

        # DatetimeIndex
        shifted = datetime_series.tshift(1)
        unshifted = shifted.tshift(-1)

        tm.assert_series_equal(datetime_series, unshifted)

        shifted2 = datetime_series.tshift(freq=datetime_series.index.freq)
        tm.assert_series_equal(shifted, shifted2)

        inferred_ts = Series(
            datetime_series.values, Index(np.asarray(datetime_series.index)), name="ts"
        )
        shifted = inferred_ts.tshift(1)
        unshifted = shifted.tshift(-1)
        tm.assert_series_equal(shifted, datetime_series.tshift(1))
        tm.assert_series_equal(unshifted, inferred_ts)

        no_freq = datetime_series[[0, 5, 7]]
        msg = "Freq was not given and was not set in the index"
        with pytest.raises(ValueError, match=msg):
            no_freq.tshift()

    def test_shift_int(self, datetime_series):
        ts = datetime_series.astype(int)
        shifted = ts.shift(1)
        expected = ts.astype(float).shift(1)
        tm.assert_series_equal(shifted, expected)

    def test_shift_object_non_scalar_fill(self):
        # shift requires scalar fill_value except for object dtype
        ser = Series(range(3))
        with pytest.raises(ValueError, match="fill_value must be a scalar"):
            ser.shift(1, fill_value=[])

        df = ser.to_frame()
        with pytest.raises(ValueError, match="fill_value must be a scalar"):
            df.shift(1, fill_value=np.arange(3))

        obj_ser = ser.astype(object)
        result = obj_ser.shift(1, fill_value={})
        assert result[0] == {}

        obj_df = obj_ser.to_frame()
        result = obj_df.shift(1, fill_value={})
        assert result.iloc[0, 0] == {}

    def test_shift_categorical(self):
        # GH#9416
        s = pd.Series(["a", "b", "c", "d"], dtype="category")

        tm.assert_series_equal(s.iloc[:-1], s.shift(1).shift(-1).dropna())

        sp1 = s.shift(1)
        tm.assert_index_equal(s.index, sp1.index)
        assert np.all(sp1.values.codes[:1] == -1)
        assert np.all(s.values.codes[:-1] == sp1.values.codes[1:])

        sn2 = s.shift(-2)
        tm.assert_index_equal(s.index, sn2.index)
        assert np.all(sn2.values.codes[-2:] == -1)
        assert np.all(s.values.codes[2:] == sn2.values.codes[:-2])

        tm.assert_index_equal(s.values.categories, sp1.values.categories)
        tm.assert_index_equal(s.values.categories, sn2.values.categories)
