from datetime import datetime

import numpy as np
import pytest

from pandas.errors import OutOfBoundsDatetime

import pandas as pd
import pandas._testing as tm


class TestSeriesClip:
    def test_clip(self, datetime_series):
        val = datetime_series.median()

        assert datetime_series.clip(lower=val).min() == val
        assert datetime_series.clip(upper=val).max() == val

        result = datetime_series.clip(-0.5, 0.5)
        expected = np.clip(datetime_series, -0.5, 0.5)
        tm.assert_series_equal(result, expected)
        assert isinstance(expected, pd.Series)

    def test_clip_types_and_nulls(self):
        sers = [
            pd.Series([np.nan, 1.0, 2.0, 3.0]),
            pd.Series([None, "a", "b", "c"]),
            pd.Series(pd.to_datetime([np.nan, 1, 2, 3], unit="D")),
        ]

        for s in sers:
            thresh = s[2]
            lower = s.clip(lower=thresh)
            upper = s.clip(upper=thresh)
            assert lower[pd.notna(lower)].min() == thresh
            assert upper[pd.notna(upper)].max() == thresh
            assert list(pd.isna(s)) == list(pd.isna(lower))
            assert list(pd.isna(s)) == list(pd.isna(upper))

    def test_series_clipping_with_na_values(self, any_numeric_ea_dtype):
        # Ensure that clipping method can handle NA values with out failing
        # GH#40581

        ser = pd.Series([pd.NA, 1.0, 3.0], dtype=any_numeric_ea_dtype)
        s_clipped_upper = ser.clip(upper=2.0)
        s_clipped_lower = ser.clip(lower=2.0)

        expected_upper = pd.Series([pd.NA, 1.0, 2.0], dtype=any_numeric_ea_dtype)
        expected_lower = pd.Series([pd.NA, 2.0, 3.0], dtype=any_numeric_ea_dtype)

        tm.assert_series_equal(s_clipped_upper, expected_upper)
        tm.assert_series_equal(s_clipped_lower, expected_lower)

    def test_clip_with_na_args(self):
        """Should process np.nan argument as None"""
        # GH#17276
        s = pd.Series([1, 2, 3])

        tm.assert_series_equal(s.clip(np.nan), pd.Series([1, 2, 3]))
        tm.assert_series_equal(s.clip(upper=np.nan, lower=np.nan), pd.Series([1, 2, 3]))

        # GH#19992

        res = s.clip(lower=[0, 4, np.nan])
        tm.assert_series_equal(res, pd.Series([1, 4, 3.0]))
        res = s.clip(upper=[1, np.nan, 1])
        tm.assert_series_equal(res, pd.Series([1, 2, 1.0]))

        # GH#40420
        s = pd.Series([1, 2, 3])
        result = s.clip(0, [np.nan, np.nan, np.nan])
        tm.assert_series_equal(s, result)

    def test_clip_against_series(self):
        # GH#6966

        s = pd.Series([1.0, 1.0, 4.0])

        lower = pd.Series([1.0, 2.0, 3.0])
        upper = pd.Series([1.5, 2.5, 3.5])

        tm.assert_series_equal(s.clip(lower, upper), pd.Series([1.0, 2.0, 3.5]))
        tm.assert_series_equal(s.clip(1.5, upper), pd.Series([1.5, 1.5, 3.5]))

    @pytest.mark.parametrize("inplace", [True, False])
    @pytest.mark.parametrize("upper", [[1, 2, 3], np.asarray([1, 2, 3])])
    def test_clip_against_list_like(self, inplace, upper):
        # GH#15390
        original = pd.Series([5, 6, 7])
        result = original.clip(upper=upper, inplace=inplace)
        expected = pd.Series([1, 2, 3])

        if inplace:
            assert result is original
        tm.assert_series_equal(result, expected, check_exact=True)

    def test_clip_with_datetimes(self):
        # GH#11838
        # naive and tz-aware datetimes

        t = pd.Timestamp("2015-12-01 09:30:30")
        s = pd.Series(
            [pd.Timestamp("2015-12-01 09:30:00"), pd.Timestamp("2015-12-01 09:31:00")]
        )
        result = s.clip(upper=t)
        expected = pd.Series(
            [pd.Timestamp("2015-12-01 09:30:00"), pd.Timestamp("2015-12-01 09:30:30")]
        )
        tm.assert_series_equal(result, expected)

        t = pd.Timestamp("2015-12-01 09:30:30", tz="US/Eastern")
        s = pd.Series(
            [
                pd.Timestamp("2015-12-01 09:30:00", tz="US/Eastern"),
                pd.Timestamp("2015-12-01 09:31:00", tz="US/Eastern"),
            ]
        )
        result = s.clip(upper=t)
        expected = pd.Series(
            [
                pd.Timestamp("2015-12-01 09:30:00", tz="US/Eastern"),
                pd.Timestamp("2015-12-01 09:30:30", tz="US/Eastern"),
            ]
        )
        tm.assert_series_equal(result, expected)

    def test_clip_with_timestamps_and_oob_datetimes_object(self):
        # GH-42794
        ser = pd.Series([datetime(1, 1, 1), datetime(9999, 9, 9)], dtype=object)

        result = ser.clip(lower=pd.Timestamp.min, upper=pd.Timestamp.max)
        expected = pd.Series([pd.Timestamp.min, pd.Timestamp.max], dtype=object)

        tm.assert_series_equal(result, expected)

    def test_clip_with_timestamps_and_oob_datetimes_non_nano(self):
        # GH#56410
        dtype = "M8[us]"
        ser = pd.Series([datetime(1, 1, 1), datetime(9999, 9, 9)], dtype=dtype)

        msg = (
            r"Incompatible \(high-resolution\) value for dtype='datetime64\[us\]'. "
            "Explicitly cast before operating"
        )
        with pytest.raises(OutOfBoundsDatetime, match=msg):
            ser.clip(lower=pd.Timestamp.min, upper=pd.Timestamp.max)

        lower = pd.Timestamp.min.as_unit("us")
        upper = pd.Timestamp.max.as_unit("us")
        result = ser.clip(lower=lower, upper=upper)
        expected = pd.Series([lower, upper], dtype=dtype)

        tm.assert_series_equal(result, expected)

    def test_clip_with_scalar_array_lower(self):
        # GH#59053
        ser = pd.Series([-1, 2, 3])
        lower = np.array(0)
        result = ser.clip(lower=lower)
        expected = pd.Series([0, 2, 3])
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("bound", [1e-8, 1e-15])
    def test_clip_int_with_tiny_float_bound(self, bound):
        # GH#25066 clipping an int Series with a tiny float bound (1e-8 lands
        # exactly on the old downcast atol) must keep the float values, not
        # silently round them back to int
        ser = pd.Series([0, 1, 0, 1])
        result = ser.clip(bound, 1)
        expected = pd.Series([bound, 1.0, bound, 1.0])
        tm.assert_series_equal(result, expected)
        # the np.clip ufunc dispatches to Series.clip
        tm.assert_series_equal(np.clip(ser, bound, 1), expected)
