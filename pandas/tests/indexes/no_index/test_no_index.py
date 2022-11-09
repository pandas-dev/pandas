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
            type(idx)(10, name="Joe")

    def test_create_index_existing_name(self, simple_index):
        # GH11193, when an existing index is passed, and a new name is not
        # specified, the new index should inherit the previous object name
        expected = simple_index
        with pytest.raises(TypeError, match="Can't set name of NoIndex!"):
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

    def test_equals_op(self, simple_index):
        # GH9947, GH10637
        index_a = simple_index

        n = len(index_a)
        index_b = index_a[0:-1]
        index_c = index_a[0:-1].append(index_a[-2:-1])
        index_d = index_a[0:1]

        msg = "Lengths must match|could not be broadcast"
        with pytest.raises(ValueError, match=msg):
            index_a == index_b
        expected = np.array([True] * n)
        tm.assert_numpy_array_equal(index_a == index_a, expected)
        tm.assert_numpy_array_equal(index_a == index_c, expected)

        # test comparisons with numpy arrays
        array_a = np.array(index_a)
        array_b = np.array(index_a[0:-1])
        array_c = np.array(index_a[0:-1].append(index_a[-2:-1]))
        array_d = np.array(index_a[0:1])
        with pytest.raises(ValueError, match=msg):
            index_a == array_b
        tm.assert_numpy_array_equal(index_a == array_a, expected)
        tm.assert_numpy_array_equal(index_a == array_c, expected)

        # test comparisons with Series
        series_a = Series(array_a)
        series_b = Series(array_b)
        series_c = Series(array_c)
        series_d = Series(array_d)
        with pytest.raises(ValueError, match=msg):
            index_a == series_b

        tm.assert_numpy_array_equal(index_a == series_a, expected)
        tm.assert_numpy_array_equal(index_a == series_c, expected)

        # cases where length is 1 for one of them
        with pytest.raises(ValueError, match="Lengths must match"):
            index_a == index_d
        with pytest.raises(ValueError, match="Lengths must match"):
            index_a == series_d
        with pytest.raises(ValueError, match="Lengths must match"):
            index_a == array_d
        msg = "Can only compare identically-labeled Series objects"
        with pytest.raises(ValueError, match=msg):
            series_a == series_d
        with pytest.raises(ValueError, match="Lengths must match"):
            series_a == array_d

        # comparing with a scalar should broadcast; note that we are excluding
        # MultiIndex because in this case each item in the index is a tuple of
        # length 2, and therefore is considered an array of length 2 in the
        # comparison instead of a scalar
        expected3 = np.array([False] * (len(index_a) - 2) + [True, False])
        # assuming the 2nd to last item is unique in the data
        item = index_a[-2]
        tm.assert_numpy_array_equal(index_a == item, expected3)
        # For RangeIndex we can convert to Int64Index
        tm.assert_series_equal(series_a == item, Series(expected3))

    def test_is_unique(self, simple_index):
        assert simple_index.is_unique

    def test_arithmetic_explicit_conversions(self):
        # GH 8608
        # add/sub are overridden explicitly for Float/Int Index
        index_cls = self._index_cls
        idx = index_cls(5)

        # float conversions
        with pytest.raises(NotImplementedError, match=None):
            idx * 3.2
        with pytest.raises(NotImplementedError, match=None):
            3.2 * idx

        # interops with numpy arrays
        a = np.zeros(5, dtype="float64")
        with pytest.raises(NotImplementedError, match=None):
            idx - a
        with pytest.raises(NotImplementedError, match=None):
            a - idx

    def test_invalid_dtype(self, invalid_dtype):
        # GH 29539
        dtype = invalid_dtype
        self._index_cls(3, dtype=dtype)

    def test_constructor_unwraps_index(self, dtype):
        result = self._index_cls(1)
        expected = np.array([0], dtype=dtype)
        tm.assert_numpy_array_equal(result._data, expected)

    def test_slice_specialised(self, simple_index):
        index = simple_index

        # scalar indexing
        res = index[1]
        expected = 1
        assert res == expected

        res = index[-1]
        expected = 9
        assert res == expected

        # slicing
        # slice value completion
        index_slice = index[:]
        expected = index
        tm.assert_index_equal(index_slice, expected)

        # positive slice values
        index_slice = index[7:10:2]
        expected = Index(np.array([0, 1]))
        tm.assert_index_equal(index_slice, expected, exact="equiv")

        # negative slice values
        index_slice = index[-1:-5:-2]
        expected = Index(np.array([0, 1]))
        tm.assert_index_equal(index_slice, expected, exact="equiv")

        # stop overshoot
        index_slice = index[2:100:4]
        expected = Index(np.array([0, 1]))
        tm.assert_index_equal(index_slice, expected, exact="equiv")

        # reverse
        index_slice = index[::-1]
        expected = Index(index.values)
        tm.assert_index_equal(index_slice, expected, exact="equiv")

        index_slice = index[-8::-1]
        expected = Index(np.array([0, 1, 2]))
        tm.assert_index_equal(index_slice, expected, exact="equiv")

        index_slice = index[-40::-1]
        expected = Index(np.array([], dtype=np.int64))
        tm.assert_index_equal(index_slice, expected, exact="equiv")

        index_slice = index[40::-1]
        expected = Index(range(10))
        tm.assert_index_equal(index_slice, expected, exact="equiv")

        index_slice = index[10::-1]
        expected = Index(range(10))
        tm.assert_index_equal(index_slice, expected, exact="equiv")


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
