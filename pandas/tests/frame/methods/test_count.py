import pandas as pd
import pandas._testing as tm


class TestDataFrameCount:
    def test_count(self):
        # corner case
        frame = pd.DataFrame()
        ct1 = frame.count(1)
        assert isinstance(ct1, pd.Series)

        ct2 = frame.count(0)
        assert isinstance(ct2, pd.Series)

        # GH#423
        df = pd.DataFrame(index=range(10))
        result = df.count(1)
        expected = pd.Series(0, index=df.index)
        tm.assert_series_equal(result, expected)

        df = pd.DataFrame(columns=range(10))
        result = df.count(0)
        expected = pd.Series(0, index=df.columns)
        tm.assert_series_equal(result, expected)

        df = pd.DataFrame()
        result = df.count()
        expected = pd.Series(dtype="int64")
        tm.assert_series_equal(result, expected)

    def test_count_objects(self, float_string_frame):
        dm = pd.DataFrame(float_string_frame._series)
        df = pd.DataFrame(float_string_frame._series)

        tm.assert_series_equal(dm.count(), df.count())
        tm.assert_series_equal(dm.count(1), df.count(1))
