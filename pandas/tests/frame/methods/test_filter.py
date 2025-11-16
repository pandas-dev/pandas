import numpy as np
import pytest

from pandas.errors import Pandas4Warning

import pandas as pd
from pandas import DataFrame
import pandas._testing as tm


class TestDataFrameSelect:
    def test_select(self, float_frame, float_string_frame):
        # Items
        selected = float_frame.select(["A", "B", "E"])
        assert len(selected.columns) == 2
        assert "E" not in selected

        selected = float_frame.select(["A", "B", "E"], axis="columns")
        assert len(selected.columns) == 2
        assert "E" not in selected

        # Other axis
        idx = float_frame.index[0:4]
        selected = float_frame.select(idx, axis="index")
        expected = float_frame.reindex(index=idx)
        tm.assert_frame_equal(selected, expected)

        # like
        fcopy = float_frame.copy()
        fcopy["AA"] = 1

        selected = fcopy.select(like="A")
        assert len(selected.columns) == 2
        assert "AA" in selected

        # like with ints in column names
        df = DataFrame(0.0, index=[0, 1, 2], columns=[0, 1, "_A", "_B"])
        selected = df.select(like="_")
        assert len(selected.columns) == 2

        # regex with ints in column names
        # from PR #10384
        df = DataFrame(0.0, index=[0, 1, 2], columns=["A1", 1, "B", 2, "C"])
        expected = DataFrame(
            0.0, index=[0, 1, 2], columns=pd.Index([1, 2], dtype=object)
        )
        selected = df.select(regex="^[0-9]+$")
        tm.assert_frame_equal(selected, expected)

        expected = DataFrame(0.0, index=[0, 1, 2], columns=[0, "0", 1, "1"])
        # shouldn't remove anything
        selected = expected.select(regex="^[0-9]+$")
        tm.assert_frame_equal(selected, expected)

        # pass in None
        with pytest.raises(TypeError, match="Must pass"):
            float_frame.select()
        with pytest.raises(TypeError, match="Must pass"):
            float_frame.select(items=None)
        with pytest.raises(TypeError, match="Must pass"):
            float_frame.select(axis=1)

        # test mutually exclusive arguments
        with pytest.raises(TypeError, match="mutually exclusive"):
            float_frame.select(items=["one", "three"], regex="e$", like="bbi")
        with pytest.raises(TypeError, match="mutually exclusive"):
            float_frame.select(items=["one", "three"], regex="e$", axis=1)
        with pytest.raises(TypeError, match="mutually exclusive"):
            float_frame.select(items=["one", "three"], regex="e$")
        with pytest.raises(TypeError, match="mutually exclusive"):
            float_frame.select(items=["one", "three"], like="bbi", axis=0)
        with pytest.raises(TypeError, match="mutually exclusive"):
            float_frame.select(items=["one", "three"], like="bbi")

        # objects
        selected = float_string_frame.select(like="foo")
        assert "foo" in selected

        # unicode columns, won't ascii-encode
        df = float_frame.rename(columns={"B": "\u2202"})
        selected = df.select(like="C")
        assert "C" in selected

    def test_select_regex_search(self, float_frame):
        fcopy = float_frame.copy()
        fcopy["AA"] = 1

        # regex
        selected = fcopy.select(regex="[A]+")
        assert len(selected.columns) == 2
        assert "AA" in selected

        # doesn't have to be at beginning
        df = DataFrame(
            {"aBBa": [1, 2], "BBaBB": [1, 2], "aCCa": [1, 2], "aCCaBB": [1, 2]}
        )

        result = df.select(regex="BB")
        exp = df[[x for x in df.columns if "BB" in x]]
        tm.assert_frame_equal(result, exp)

    @pytest.mark.parametrize(
        "name,expected_data",
        [
            ("a", {"a": [1, 2]}),
            ("あ", {"あ": [3, 4]}),
        ],
    )
    def test_select_unicode(self, name, expected_data):
        # GH13101
        df = DataFrame({"a": [1, 2], "あ": [3, 4]})
        expected = DataFrame(expected_data)

        tm.assert_frame_equal(df.select(like=name), expected)
        tm.assert_frame_equal(df.select(regex=name), expected)

    def test_select_bytestring(self):
        # GH13101
        name = "a"
        df = DataFrame({b"a": [1, 2], b"b": [3, 4]})
        expected = DataFrame({b"a": [1, 2]})

        tm.assert_frame_equal(df.select(like=name), expected)
        tm.assert_frame_equal(df.select(regex=name), expected)

    def test_select_corner(self):
        empty = DataFrame()

        result = empty.select([])
        tm.assert_frame_equal(result, empty)

        result = empty.select(like="foo")
        tm.assert_frame_equal(result, empty)

    def test_select_regex_non_string(self):
        # GH#5798 trying to select on non-string columns should drop,
        #  not raise
        df = DataFrame(np.random.default_rng(2).random((3, 2)), columns=["STRING", 123])
        result = df.select(regex="STRING")
        expected = df[["STRING"]]
        tm.assert_frame_equal(result, expected)

    def test_select_keep_order(self):
        # GH#54980
        df = DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
        result = df.select(items=["B", "A"])
        expected = df[["B", "A"]]
        tm.assert_frame_equal(result, expected)

    def test_select_different_dtype(self):
        # GH#54980
        df = DataFrame({1: [1, 2, 3], 2: [4, 5, 6]})
        result = df.select(items=["B", "A"])
        expected = df[[]]
        tm.assert_frame_equal(result, expected)

    def test_filter_deprecated(self):
        # GH#26642
        df = DataFrame({1: [1, 2, 3], 2: [4, 5, 6]})
        msg = "DataFrame.filter is deprecated"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            df.filter(items=["B", "A"])

        ser = df[1]
        msg = "Series.filter is deprecated"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            ser.filter([0, 1])
