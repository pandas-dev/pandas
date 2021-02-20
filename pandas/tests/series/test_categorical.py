import pytest

import pandas as pd
from pandas import Categorical
import pandas._testing as tm


class TestCategoricalSeries:
    def test_setitem_undefined_category_raises(self):
        ser = pd.Series(Categorical(["a", "b", "c"]))
        msg = (
            "Cannot setitem on a Categorical with a new category, "
            "set the categories first"
        )
        with pytest.raises(ValueError, match=msg):
            ser.loc[2] = "d"

    def test_concat_undefined_category_raises(self):
        ser = pd.Series(Categorical(["a", "b", "c"]))
        msg = (
            "Cannot setitem on a Categorical with a new category, "
            "set the categories first"
        )
        with pytest.raises(ValueError, match=msg):
            ser.loc[3] = "d"

    def test_loc_category_dtype_retention(self):
        # Case 1
        df = pd.DataFrame(
            {
                "int": [0, 1, 2],
                "cat": Categorical(["a", "b", "c"], categories=["a", "b", "c"]),
            }
        )
        df.loc[3] = [3, "c"]
        expected = pd.DataFrame(
            {
                "int": [0, 1, 2, 3],
                "cat": Categorical(["a", "b", "c", "c"], categories=["a", "b", "c"]),
            }
        )
        tm.assert_frame_equal(df, expected)

        # Case 2
        ser = pd.Series(Categorical(["a", "b", "c"]))
        ser.loc[3] = "c"
        expected = pd.Series(Categorical(["a", "b", "c", "c"]))
        tm.assert_series_equal(ser, expected)

        # Case 3
        ser = pd.Series(Categorical([1, 2, 3]))
        ser.loc[3] = 3
        expected = pd.Series(Categorical([1, 2, 3, 3]))
        tm.assert_series_equal(ser, expected)

        # Case 4
        ser = pd.Series(Categorical([1, 2, 3]))
        ser.loc[3] = pd.NA
        expected = pd.Series(Categorical([1, 2, 3, pd.NA]))
        tm.assert_series_equal(ser, expected)
