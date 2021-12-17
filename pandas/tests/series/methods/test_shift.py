import numpy as np
import pytest

from pandas.errors import NullFrequencyError

from pandas import (
    DatetimeIndex,
    Series,
    TimedeltaIndex,
    date_range,
    offsets,
)
import pandas._testing as tm


class TestShift:
    @pytest.mark.parametrize(
        "ser",
        [
            Series([np.arange(5)]),
            date_range("1/1/2011", periods=24, freq="H"),
            Series(range(5), index=date_range("2017", periods=5)),
        ],
    )
    @pytest.mark.parametrize("shift_size", [0, 1, 2])
    def test_shift_always_copy(self, ser, shift_size):
        # GH22397
        assert ser.shift(shift_size) is not ser

    def test_shift(self, datetime_series):
        shifted = datetime_series.shift(1)
        unshifted = shifted.shift(-1)

        tm.assert_index_equal(shifted.index, datetime_series.index)
        tm.assert_index_equal(unshifted.index, datetime_series.index)
        tm.assert_numpy_array_equal(
            unshifted.dropna().values, datetime_series.values[:-1]
        )

    def test_shift_with_tz(self):
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
        s = Series(["a", "b", "c", "d"], dtype="category")

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

    @pytest.mark.parametrize(
        "input_data, output_data",
        [(np.empty(shape=(0,)), []), (np.ones(shape=(2,)), [np.nan, 1.0])],
    )
    def test_shift_non_writable_array(self, input_data, output_data):
        # GH21049 Verify whether non writable numpy array is shiftable
        input_data.setflags(write=False)

        result = Series(input_data).shift(1)
        expected = Series(output_data, dtype="float64")

        tm.assert_series_equal(result, expected)
