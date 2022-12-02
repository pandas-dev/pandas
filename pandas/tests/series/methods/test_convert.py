from datetime import datetime

import pytest

from pandas import (
    Series,
    Timestamp,
)
import pandas._testing as tm


class TestConvert:
    def test_convert(self):
        # GH#10265
        dt = datetime(2001, 1, 1, 0, 0)
        td = dt - datetime(2000, 1, 1, 0, 0)

        # Test coercion with mixed types
        ser = Series(["a", "3.1415", dt, td])

        # Test standard conversion returns original
        results = ser._convert(datetime=True)
        tm.assert_series_equal(results, ser)

        results = ser._convert(timedelta=True)
        tm.assert_series_equal(results, ser)

    def test_convert_numeric_strings_with_other_true_args(self):
        # test pass-through and non-conversion when other types selected
        ser = Series(["1.0", "2.0", "3.0"])
        results = ser._convert(datetime=True, timedelta=True)
        tm.assert_series_equal(results, ser)

    def test_convert_datetime_objects(self):
        ser = Series(
            [datetime(2001, 1, 1, 0, 0), datetime(2001, 1, 1, 0, 0)], dtype="O"
        )
        results = ser._convert(datetime=True, timedelta=True)
        expected = Series([datetime(2001, 1, 1, 0, 0), datetime(2001, 1, 1, 0, 0)])
        tm.assert_series_equal(results, expected)
        results = ser._convert(datetime=False, timedelta=True)
        tm.assert_series_equal(results, ser)

    def test_convert_datetime64(self):
        # no-op if already dt64 dtype
        ser = Series(
            [
                datetime(2001, 1, 1, 0, 0),
                datetime(2001, 1, 2, 0, 0),
                datetime(2001, 1, 3, 0, 0),
            ]
        )

        result = ser._convert(datetime=True)
        expected = Series(
            [Timestamp("20010101"), Timestamp("20010102"), Timestamp("20010103")],
            dtype="M8[ns]",
        )
        tm.assert_series_equal(result, expected)

        result = ser._convert(datetime=True)
        tm.assert_series_equal(result, expected)

    def test_convert_timedeltas(self):
        td = datetime(2001, 1, 1, 0, 0) - datetime(2000, 1, 1, 0, 0)
        ser = Series([td, td], dtype="O")
        results = ser._convert(datetime=True, timedelta=True)
        expected = Series([td, td])
        tm.assert_series_equal(results, expected)
        results = ser._convert(datetime=True, timedelta=False)
        tm.assert_series_equal(results, ser)

    def test_convert_preserve_non_object(self):
        # preserve if non-object
        ser = Series([1], dtype="float32")
        result = ser._convert(datetime=True)
        tm.assert_series_equal(result, ser)

    def test_convert_no_arg_error(self):
        ser = Series(["1.0", "2"])
        msg = r"At least one of datetime or timedelta must be True\."
        with pytest.raises(ValueError, match=msg):
            ser._convert()

    def test_convert_preserve_bool(self):
        ser = Series([1, True, 3, 5], dtype=object)
        res = ser._convert(datetime=True)
        tm.assert_series_equal(res, ser)

    def test_convert_preserve_all_bool(self):
        ser = Series([False, True, False, False], dtype=object)
        res = ser._convert(datetime=True)
        expected = Series([False, True, False, False], dtype=bool)
        tm.assert_series_equal(res, expected)
