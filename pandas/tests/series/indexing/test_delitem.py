import pytest

import pandas as pd
import pandas._testing as tm


class TestSeriesDelItem:
    def test_delitem(self):
        # GH#5542
        # should delete the item inplace
        s = pd.Series(range(5))
        del s[0]

        expected = pd.Series(range(1, 5), index=range(1, 5))
        tm.assert_series_equal(s, expected)

        del s[1]
        expected = pd.Series(range(2, 5), index=range(2, 5))
        tm.assert_series_equal(s, expected)

        # only 1 left, del, add, del
        s = pd.Series(1)
        del s[0]
        tm.assert_series_equal(
            s, pd.Series(dtype="int64", index=pd.Index([], dtype="int64"))
        )
        s[0] = 1
        tm.assert_series_equal(s, pd.Series(1))
        del s[0]
        tm.assert_series_equal(
            s, pd.Series(dtype="int64", index=pd.Index([], dtype="int64"))
        )

    def test_delitem_object_index(self):
        # Index(dtype=object)
        s = pd.Series(1, index=pd.Index(["a"], dtype="str"))
        del s["a"]
        tm.assert_series_equal(
            s, pd.Series(dtype="int64", index=pd.Index([], dtype="str"))
        )
        s["a"] = 1
        tm.assert_series_equal(s, pd.Series(1, index=pd.Index(["a"], dtype="str")))
        del s["a"]
        tm.assert_series_equal(
            s, pd.Series(dtype="int64", index=pd.Index([], dtype="str"))
        )

    def test_delitem_missing_key(self):
        # empty
        s = pd.Series(dtype=object)

        with pytest.raises(KeyError, match=r"^0$"):
            del s[0]

    def test_delitem_extension_dtype(self):
        # GH#40386
        # DatetimeTZDtype
        dti = pd.date_range("2016-01-01", periods=3, tz="US/Pacific")
        ser = pd.Series(dti)

        expected = ser[[0, 2]]
        del ser[1]
        assert ser.dtype == dti.dtype
        tm.assert_series_equal(ser, expected)

        # PeriodDtype
        pi = dti.tz_localize(None).to_period("D")
        ser = pd.Series(pi)

        expected = ser[:2]
        del ser[2]
        assert ser.dtype == pi.dtype
        tm.assert_series_equal(ser, expected)
