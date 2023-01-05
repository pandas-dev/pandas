from pandas import (
    DataFrame,
    Series,
)
import pandas._testing as tm


class TestDataFrameSetItem:
    def test_isetitem_ea_df(self):
        # GH#49922
        df = DataFrame([[1, 2, 3], [4, 5, 6]])
        rhs = DataFrame([[11, 12], [13, 14]], dtype="Int64")

        df.isetitem([0, 1], rhs)
        expected = DataFrame(
            {
                0: Series([11, 13], dtype="Int64"),
                1: Series([12, 14], dtype="Int64"),
                2: [3, 6],
            }
        )
        tm.assert_frame_equal(df, expected)

    def test_isetitem_ea_df_scalar_indexer(self):
        # GH#49922
        df = DataFrame([[1, 2, 3], [4, 5, 6]])
        rhs = DataFrame([[11], [13]], dtype="Int64")

        df.isetitem(2, rhs)
        expected = DataFrame(
            {
                0: [1, 4],
                1: [2, 5],
                2: Series([11, 13], dtype="Int64"),
            }
        )
        tm.assert_frame_equal(df, expected)
