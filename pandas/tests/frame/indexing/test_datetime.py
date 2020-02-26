import pandas as pd
from pandas import DataFrame, Index, Series, date_range, notna
import pandas._testing as tm


class TestDataFrameIndexingDatetimeWithTZ:
    def test_setitem(self, timezone_frame):

        df = timezone_frame
        idx = df["B"].rename("foo")

        # setitem
        df["C"] = idx
        tm.assert_series_equal(df["C"], Series(idx, name="C"))

        df["D"] = "foo"
        df["D"] = idx
        tm.assert_series_equal(df["D"], Series(idx, name="D"))
        del df["D"]

        # assert that A & C are not sharing the same base (e.g. they
        # are copies)
        b1 = df._data.blocks[1]
        b2 = df._data.blocks[2]
        tm.assert_extension_array_equal(b1.values, b2.values)
        assert id(b1.values._data.base) != id(b2.values._data.base)

        # with nan
        df2 = df.copy()
        df2.iloc[1, 1] = pd.NaT
        df2.iloc[1, 2] = pd.NaT
        result = df2["B"]
        tm.assert_series_equal(notna(result), Series([True, False, True], name="B"))
        tm.assert_series_equal(df2.dtypes, df.dtypes)

    def test_set_reset(self):

        idx = Index(date_range("20130101", periods=3, tz="US/Eastern"), name="foo")

        # set/reset
        df = DataFrame({"A": [0, 1, 2]}, index=idx)
        result = df.reset_index()
        assert result["foo"].dtype, "M8[ns, US/Eastern"

        df = result.set_index("foo")
        tm.assert_index_equal(df.index, idx)

    def test_scalar_assignment(self):
        # issue #19843
        df = pd.DataFrame(index=(0, 1, 2))
        df["now"] = pd.Timestamp("20130101", tz="UTC")
        expected = pd.DataFrame(
            {"now": pd.Timestamp("20130101", tz="UTC")}, index=[0, 1, 2]
        )
        tm.assert_frame_equal(df, expected)
