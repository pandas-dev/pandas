import io

import numpy as np
import pytest

from pandas._config import option_context

import pandas as pd
from pandas import (
    DataFrame,
    Series,
)
import pandas._testing as tm
from pandas.core.indexes.api import (
    Float64Index,
    Index,
    Int64Index,
    NoIndex,
    RangeIndex,
)
from pandas.tests.indexes.ranges.test_range import TestRangeIndex

# aliases to make some tests easier to read
NI = NoIndex
RI = RangeIndex
I64 = Int64Index
F64 = Float64Index
OI = Index


class TestNoIndex(TestRangeIndex):
    _index_cls = NoIndex

    @pytest.fixture
    def simple_index(self) -> Index:
        return self._index_cls(10)

    def test_constructor_name_unhashable(self, simple_index):
        # GH#29069 check that name is hashable
        # See also same-named test in tests.series.test_constructors
        idx = simple_index
        with pytest.raises(TypeError, match="Can't set name of NoIndex!"):
            type(idx)(idx, name="Joe")

    def test_create_index_existing_name(self, simple_index):
        # GH11193, when an existing index is passed, and a new name is not
        # specified, the new index should inherit the previous object name
        expected = simple_index
        with pytest.raises(ValueError, match="Can't set name of NoIndex!"):
            expected.name = "foo"

    def test_repeat(self, simple_index):
        rep = 2
        idx = simple_index.copy()
        expected = NoIndex(20)
        tm.assert_index_equal(idx.repeat(rep), expected)

        idx = simple_index
        rep = np.arange(len(idx))
        expected = NoIndex(45)
        tm.assert_index_equal(idx.repeat(rep), expected)


class TestCommon:
    @pytest.fixture
    def ser1(self):
        with option_context("mode.no_default_index", True):
            res = Series([1, 2, 3])
        return res

    @pytest.fixture
    def ser2(self):
        with option_context("mode.no_default_index", True):
            res = Series([4, 5, 6])
        return res

    @pytest.fixture
    def df1(self):
        with option_context("mode.no_default_index", True):
            res = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
        return res

    @pytest.fixture
    def df2(self):
        with option_context("mode.no_default_index", True):
            res = DataFrame({"a": [6, 5, 4], "b": [1, 4, 2]})
        return res

    def test_boolean_mask(self, ser1):
        mask = ser1 > 1
        result = ser1[mask]
        expected = Series([2, 3], index=NoIndex(2))
        tm.assert_series_equal(result, expected)

        ser1[mask] = 4
        expected = Series([1, 4, 4], index=NoIndex(3))
        tm.assert_series_equal(ser1, expected)

    def test_join(self, ser1, ser2, df1, df2):
        result = df1.join(df2, lsuffix="_df1")
        expected = DataFrame(
            {
                "a_df1": [1, 2, 3],
                "b_df1": [4, 5, 6],
                "a": [6, 5, 4],
                "b": [1, 4, 2],
            },
            index=NoIndex(3),
        )
        tm.assert_frame_equal(result, expected)

    def test_loc(self, df1):
        with pytest.raises(
            IndexError, match="Cannot use label-based indexing on NoIndex!"
        ):
            df1.loc[0, "a"]
        with pytest.raises(
            IndexError, match="Cannot use label-based indexing on NoIndex!"
        ):
            df1.loc[0]

        result = df1.loc[:, "a"]
        expected = Series([1, 2, 3], index=NoIndex(3), name="a")
        tm.assert_series_equal(result, expected)

        mask = df1["a"] > 2
        result = df1.loc[mask]
        expected = DataFrame({"a": [3], "b": [6]}, index=NoIndex(1))
        tm.assert_frame_equal(result, expected)

        result = df1.loc[df1["a"] > 2, "a"]
        expected = Series([3], index=NoIndex(1), name="a")
        tm.assert_series_equal(result, expected)

        result = df1.iloc[1:]
        expected = DataFrame(
            {
                "a": [
                    2,
                    3,
                ],
                "b": [5, 6],
            }
        )
        tm.assert_frame_equal(result, expected)

    def test_alignment(self, df1):
        with pytest.raises(TypeError, match="Can't join NoIndex of different lengths"):
            result = df1 + df1.iloc[1:]
        result = df1 + df1
        expected = DataFrame({"a": [2, 4, 6], "b": [8, 10, 12]}, index=NoIndex(3))
        tm.assert_frame_equal(result, expected)

    def test_reader(self):
        with option_context("mode.no_default_index", True):
            result = pd.read_csv(io.StringIO("data\n1\n"))
        expected = DataFrame({"data": [1]}, index=NoIndex(1))
        tm.assert_frame_equal(result, expected)

    def test_repr(self):
        with option_context("mode.no_default_index", True):
            df = DataFrame({"a": [1, 2, 3] * 50})
        result = repr(df)
        expected = (
            " a\n"
            " 1\n"
            " 2\n"
            " 3\n"
            " 1\n"
            " 2\n"
            "..\n"
            " 2\n"
            " 3\n"
            " 1\n"
            " 2\n"
            " 3\n"
            "\n[150 rows x 1 columns]"
        )
        assert result == expected

        result = repr(df["a"])
        expected = (
            "1 \n"
            "2 \n"
            "3 \n"
            "1 \n"
            "2 \n"
            "..\n"
            "2 \n"
            "3 \n"
            "1 \n"
            "2 \n"
            "3 \n"
            "Name: a, Length: 150, dtype: int64"
        )
        assert result == expected

    def test_concat(self, df1):
        result = pd.concat([df1, df1])
        expected = DataFrame(
            {
                "a": [1, 2, 3, 1, 2, 3],
                "b": [4, 5, 6, 4, 5, 6],
            },
            index=NoIndex(6),
        )
        tm.assert_frame_equal(result, expected)

    def test_merge(self, df1):
        result = df1.merge(df1, on="a")
        expected = DataFrame(
            {
                "a": [1, 2, 3],
                "b_x": [4, 5, 6],
                "b_y": [4, 5, 6],
            },
            index=NoIndex(3),
        )
        tm.assert_frame_equal(result, expected)
