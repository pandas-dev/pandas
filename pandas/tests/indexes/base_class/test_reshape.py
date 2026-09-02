"""
Tests for ndarray-like method on the base Index class
"""

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class TestReshape:
    def test_repeat(self):
        repeats = 2
        index = pd.Index([1, 2, 3])
        expected = pd.Index([1, 1, 2, 2, 3, 3])

        result = index.repeat(repeats)
        tm.assert_index_equal(result, expected)

    def test_insert(self):
        # GH 7256
        # validate neg/pos inserts
        result = pd.Index(["b", "c", "d"])

        # test 0th element
        tm.assert_index_equal(pd.Index(["a", "b", "c", "d"]), result.insert(0, "a"))

        # test Nth element that follows Python list behavior
        tm.assert_index_equal(pd.Index(["b", "c", "e", "d"]), result.insert(-1, "e"))

        # test loc +/- neq (0, -1)
        tm.assert_index_equal(result.insert(1, "z"), result.insert(-2, "z"))

        # test empty
        null_index = pd.Index([])
        tm.assert_index_equal(pd.Index(["a"]), null_index.insert(0, "a"))

    @pytest.mark.parametrize("val", [True, False])
    def test_insert_bool_into_numeric_ea(self, val):
        # GH#61709 - bool should not be silently cast to numeric
        float_idx = pd.Index([1.0, 2.0, 3.0], dtype="Float64")
        result = float_idx.insert(1, val)
        expected = pd.Index([1.0, val, 2.0, 3.0], dtype=object)
        tm.assert_index_equal(result, expected)

        int_idx = pd.Index([1, 2, 3], dtype="Int64")
        result = int_idx.insert(1, val)
        expected = pd.Index([1, val, 2, 3], dtype=object)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("val", [0, 1])
    def test_insert_int_into_boolean_ea(self, val):
        # GH#61709 - int should not be silently cast to bool
        bool_idx = pd.Index([True, False, True], dtype="boolean")
        result = bool_idx.insert(1, val)
        expected = pd.Index([True, val, False, True], dtype=object)
        tm.assert_index_equal(result, expected)

    def test_insert_missing(self, nulls_fixture, using_infer_string):
        # GH#22295
        # test there is no mangling of NA values
        expected = pd.Index(["a", nulls_fixture, "b", "c"], dtype=object)
        result = pd.Index(list("abc"), dtype=object).insert(
            1, pd.Index([nulls_fixture], dtype=object)
        )
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize(
        "val", [(1, 2), np.datetime64("2019-12-31"), np.timedelta64(1, "D")]
    )
    @pytest.mark.parametrize("loc", [-1, 2])
    def test_insert_datetime_into_object(self, loc, val):
        # GH#44509
        idx = pd.Index(["1", "2", "3"])
        result = idx.insert(loc, val)
        expected = pd.Index(["1", "2", val, "3"])
        tm.assert_index_equal(result, expected)
        assert type(expected[2]) is type(val)

    def test_insert_none_into_string_numpy(self, string_dtype_no_object):
        # GH#55365
        index = pd.Index(["a", "b", "c"], dtype=string_dtype_no_object)
        result = index.insert(-1, None)
        expected = pd.Index(["a", "b", None, "c"], dtype=string_dtype_no_object)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize(
        "pos,expected",
        [
            (0, pd.Index(["b", "c", "d"], name="index")),
            (-1, pd.Index(["a", "b", "c"], name="index")),
        ],
    )
    def test_delete(self, pos, expected):
        index = pd.Index(["a", "b", "c", "d"], name="index")
        result = index.delete(pos)
        tm.assert_index_equal(result, expected)
        assert result.name == expected.name

    def test_delete_raises(self):
        index = pd.Index(["a", "b", "c", "d"], name="index")
        msg = "index 5 is out of bounds for axis 0 with size 4"
        with pytest.raises(IndexError, match=msg):
            index.delete(5)

    def test_append_multiple(self):
        index = pd.Index(["a", "b", "c", "d", "e", "f"])

        foos = [index[:2], index[2:4], index[4:]]
        result = foos[0].append(foos[1:])
        tm.assert_index_equal(result, index)

        # empty
        result = index.append([])
        tm.assert_index_equal(result, index)
