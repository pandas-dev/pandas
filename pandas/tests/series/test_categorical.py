import pytest
import pandas as pd
import pandas._testing as tm

from pandas import Categorical, Index


class TestCategoricalSeries:
    def test_loc_new_category_series_raises(self):
        ser = pd.Series(Categorical(["a", "b", "c"]))
        msg = "Cannot setitem on a Categorical with a new category"
        with pytest.raises(ValueError, match=msg):
            ser.loc[3] = "d"

    def test_unused_category_retention(self):
        # Init case
        exp_cats = Index(["a", "b", "c", "d"])
        ser = pd.Series(Categorical(["a", "b", "c"], categories=exp_cats))
        tm.assert_index_equal(ser.cat.categories, exp_cats)

        # Modify case
        ser.loc[0] = "b"
        expected = pd.Series(Categorical(["b", "b", "c"], categories=exp_cats))
        tm.assert_index_equal(ser.cat.categories, exp_cats)
        tm.assert_series_equal(ser, expected)

    def test_loc_new_category_row_raises(self):
        data = {
            "int": [0, 1, 2],
            "cat": Categorical(["a", "b", "c"], categories=["a", "b", "c"]),
        }
        df = pd.DataFrame(data)
        msg = "Cannot setitem on a Categorical with a new category"
        with pytest.raises(ValueError, match=msg):
            df.loc[3] = [3, "d"]

    def test_loc_new_row_category_dtype_retention(self):
        df_data = {
            "int": [0, 1, 2],
            "cat": pd.Categorical(["a", "b", "c"], categories=["a", "b", "c"]),
        }
        df = pd.DataFrame(df_data)
        df.loc[3] = [3, "c"]

        expected_data = {
            "int": [0, 1, 2, 3],
            "cat": pd.Categorical(["a", "b", "c", "c"], categories=["a", "b", "c"]),
        }
        expected = pd.DataFrame(expected_data)
        tm.assert_frame_equal(df, expected)
