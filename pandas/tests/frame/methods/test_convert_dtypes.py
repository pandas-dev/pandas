import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class TestConvertDtypes:
    @pytest.mark.parametrize(
        "convert_integer, expected", [(False, np.dtype("int32")), (True, "Int32")]
    )
    def test_convert_dtypes(self, convert_integer, expected):
        # Specific types are tested in tests/series/test_dtypes.py
        # Just check that it works for DataFrame here
        df = pd.DataFrame(
            {
                "a": pd.Series([1, 2, 3], dtype=np.dtype("int32")),
                "b": pd.Series(["x", "y", "z"], dtype=np.dtype("O")),
            }
        )
        result = df.convert_dtypes(True, True, convert_integer, False)
        expected = pd.DataFrame(
            {
                "a": pd.Series([1, 2, 3], dtype=expected),
                "b": pd.Series(["x", "y", "z"], dtype="string"),
            }
        )
        tm.assert_frame_equal(result, expected)

    def test_convert_empty(self):
        # Empty DataFrame can pass convert_dtypes, see GH#40393
        empty_df = pd.DataFrame()
        tm.assert_frame_equal(empty_df, empty_df.convert_dtypes())
