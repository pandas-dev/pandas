import pandas as pd
import pandas._testing as tm


class TestDataFrameReductions:
    def test_min_max_dt64_with_NaT(self):
        # Both NaT and Timestamp are in DataFrame.
        df = pd.DataFrame({"foo": [pd.NaT, pd.NaT, pd.Timestamp("2012-05-01")]})

        res = df.min()
        exp = pd.Series([pd.Timestamp("2012-05-01")], index=["foo"])
        tm.assert_series_equal(res, exp)

        res = df.max()
        exp = pd.Series([pd.Timestamp("2012-05-01")], index=["foo"])
        tm.assert_series_equal(res, exp)

        # GH12941, only NaTs are in DataFrame.
        df = pd.DataFrame({"foo": [pd.NaT, pd.NaT]})

        res = df.min()
        exp = pd.Series([pd.NaT], index=["foo"])
        tm.assert_series_equal(res, exp)

        res = df.max()
        exp = pd.Series([pd.NaT], index=["foo"])
        tm.assert_series_equal(res, exp)
