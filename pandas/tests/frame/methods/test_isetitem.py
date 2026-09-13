import pytest

import pandas as pd
import pandas._testing as tm


class TestDataFrameSetItem:
    def test_isetitem_ea_df(self):
        # GH#49922
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        rhs = pd.DataFrame([[11, 12], [13, 14]], dtype="Int64")

        df.isetitem([0, 1], rhs)
        expected = pd.DataFrame(
            {
                0: pd.Series([11, 13], dtype="Int64"),
                1: pd.Series([12, 14], dtype="Int64"),
                2: [3, 6],
            }
        )
        tm.assert_frame_equal(df, expected)

    def test_isetitem_ea_df_scalar_indexer(self):
        # GH#49922
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        rhs = pd.DataFrame([[11], [13]], dtype="Int64")

        df.isetitem(2, rhs)
        expected = pd.DataFrame(
            {
                0: [1, 4],
                1: [2, 5],
                2: pd.Series([11, 13], dtype="Int64"),
            }
        )
        tm.assert_frame_equal(df, expected)

    def test_isetitem_dimension_mismatch(self):
        # GH#51701
        df = pd.DataFrame({"a": [1, 2], "b": [3, 4], "c": [5, 6]})
        value = df.copy()
        with pytest.raises(ValueError, match="Got 2 positions but value has 3 columns"):
            df.isetitem([1, 2], value)

        value = df.copy()
        with pytest.raises(ValueError, match="Got 2 positions but value has 1 columns"):
            df.isetitem([1, 2], value[["a"]])
