import pytest

from pandas.core.dtypes.concat import _can_cast_to_categorical

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
            "Cannot concat on a Categorical with a new category, "
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

    def test_can_cast_to_categorical(self):
        # Case 1:
        # Series of identical categorical dtype should
        # be able to concat to categorical
        ser1 = pd.Series(Categorical(["a", "b", "c"]))
        ser2 = pd.Series(Categorical(["a", "b", "c"]))
        arr = [ser1, ser2]
        assert _can_cast_to_categorical(arr) is True

        # Case 2:
        # Series of non-identical categorical dtype should
        # not be able to concat to categoorical
        ser1 = pd.Series(Categorical(["a", "b", "c"]))
        ser2 = pd.Series(Categorical(["a", "b", "d"]))
        arr = [ser1, ser2]
        assert _can_cast_to_categorical(arr) is False

        # Concat of a categorical series with a series
        # containing only values identical to the
        # categorical values should be possible

        # Case 3: For string categorical values
        ser1 = pd.Series(Categorical(["a", "b", "c"]))
        ser2 = pd.Series(["a", "a", "b"])
        arr = [ser1, ser2]
        assert _can_cast_to_categorical(arr) is True

        # Case 4: For int categorical values
        ser1 = pd.Series(Categorical([1, 2, 3]))
        ser2 = pd.Series([1, 2])
        arr = [ser1, ser2]
        assert _can_cast_to_categorical(arr) is True

        # The rest should raise because not all values
        # are present in the categorical.

        # Case 5
        ser1 = pd.Series(Categorical([1, 2, 3]))
        ser2 = pd.Series([3, 4])
        arr = [ser1, ser2]
        msg = (
            "Cannot concat on a Categorical with a new category, "
            "set the categories first"
        )
        with pytest.raises(ValueError, match=msg):
            _can_cast_to_categorical(arr)

        # Case 6
        ser1 = pd.Series(Categorical(["a", "b", "c"]))
        ser2 = pd.Series(["d", "e"])
        arr = [ser1, ser2]
        with pytest.raises(ValueError, match=msg):
            _can_cast_to_categorical(arr)
