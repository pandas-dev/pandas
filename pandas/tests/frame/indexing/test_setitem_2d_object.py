import numpy as np
import pandas as pd
import pandas._testing as tm
import pytest


class TestSetItem2DObjectArray:
    # GH#61026
    def test_setitem_2d_object_array_single_column_unravel(self):
        df = pd.DataFrame({"c1": [1, 2, 3]})
        t = np.array([["A"], ["B"], ["C"]], dtype=object)

        df["c1"] = t

        expected = pd.Series(["A", "B", "C"], name="c1")
        tm.assert_series_equal(df["c1"], expected)

    # GH#61026
    def test_setitem_2d_object_array_wrong_shape_raises(self):
        df = pd.DataFrame({"c1": [1, 2, 3]})
        t = np.array([["A", "B"], ["C", "D"], ["E", "F"]], dtype=object)

        with pytest.raises(ValueError, match="requires shape"):
            df["c1"] = t

    # GH#61026
    def test_setitem_3d_object_array_raises(self):
        df = pd.DataFrame({"c1": [1, 2, 3]})
        t = np.array([[["A"]], [["B"]], [["C"]]], dtype=object)

        with pytest.raises(ValueError, match="ndim 3"):
            df["c1"] = t

    # GH#61026
    def test_setitem_2d_string_array_regression(self):
        df = pd.DataFrame({"c1": [1, 2, 3]})
        t = np.array([["A"], ["B"], ["C"]])  # dtype '<U1'

        df["c1"] = t

        # just checking values; internal mask shape isn't part of the contract.
        assert list(df["c1"]) == ["A", "B", "C"]

    # GH#61026
    def test_setitem_multicolumn_object_array_still_ok(self):
        df = pd.DataFrame({"c1": [1, 2, 3], "c2": [4, 5, 6]})
        t = np.array([["A", "X"], ["B", "Y"], ["C", "Z"]], dtype=object)

        df[["c1", "c2"]] = t

        expected = pd.DataFrame({"c1": ["A", "B", "C"], "c2": ["X", "Y", "Z"]})
        tm.assert_frame_equal(df, expected)
