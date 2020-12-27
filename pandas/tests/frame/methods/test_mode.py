from pandas import DataFrame, Index
import pandas._testing as tm


def test_empty_df_mode():
    df = DataFrame([], columns=["a", "b"])
    result = df.mode()
    expected = DataFrame([], columns=["a", "b"], index=Index([], dtype="int64"))
    tm.assert_frame_equal(result, expected)
