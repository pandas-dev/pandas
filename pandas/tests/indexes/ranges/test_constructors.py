from datetime import datetime

import numpy as np
import pytest

from pandas import (
    Index,
    RangeIndex,
    Series,
)
import pandas._testing as tm
from pandas.core.indexes.base import maybe_sequence_to_range


class TestRangeIndexConstructors:
    @pytest.mark.parametrize("name", [None, "foo"])
    @pytest.mark.parametrize(
        "args, kwargs, start, stop, step",
        [
            ((5,), {}, 0, 5, 1),
            ((1, 5), {}, 1, 5, 1),
            ((1, 5, 2), {}, 1, 5, 2),
            ((0,), {}, 0, 0, 1),
            ((0, 0), {}, 0, 0, 1),
            ((), {"start": 0}, 0, 0, 1),
            ((), {"stop": 0}, 0, 0, 1),
        ],
    )
    def test_constructor(self, args, kwargs, start, stop, step, name):
        result = RangeIndex(*args, name=name, **kwargs)
        expected = Index(np.arange(start, stop, step, dtype=np.int64), name=name)
        assert isinstance(result, RangeIndex)
        assert result.name is name
        assert result._range == range(start, stop, step)
        tm.assert_index_equal(result, expected, exact="equiv")

    def test_constructor_invalid_args(self):
        msg = "RangeIndex\\(\\.\\.\\.\\) must be called with integers"
        with pytest.raises(TypeError, match=msg):
            RangeIndex()

        with pytest.raises(TypeError, match=msg):
            RangeIndex(name="Foo")

        # we don't allow on a bare Index
        msg = (
            r"Index\(\.\.\.\) must be called with a collection of some "
            r"kind, 0 was passed"
        )
        with pytest.raises(TypeError, match=msg):
            Index(0)

    @pytest.mark.parametrize(
        "args",
        [
            Index(["a", "b"]),
            Series(["a", "b"]),
            np.array(["a", "b"]),
            [],
            np.arange(0, 10),
            np.array([1]),
            [1],
        ],
    )
    def test_constructor_additional_invalid_args(self, args):
        msg = f"Value needs to be a scalar value, was type {type(args).__name__}"
        with pytest.raises(TypeError, match=msg):
            RangeIndex(args)

    @pytest.mark.parametrize("args", ["foo", datetime(2000, 1, 1, 0, 0)])
    def test_constructor_invalid_args_wrong_type(self, args):
        msg = f"Wrong type {type(args)} for value {args}"
        with pytest.raises(TypeError, match=msg):
            RangeIndex(args)

    def test_constructor_same(self):
        # pass thru w and w/o copy
        index = RangeIndex(1, 5, 2)
        result = RangeIndex(index, copy=False)
        assert result.identical(index)

        result = RangeIndex(index, copy=True)
        tm.assert_index_equal(result, index, exact=True)

        result = RangeIndex(index)
        tm.assert_index_equal(result, index, exact=True)

        with pytest.raises(
            ValueError,
            match="Incorrect `dtype` passed: expected signed integer, received float64",
        ):
            RangeIndex(index, dtype="float64")

    def test_constructor_range_object(self):
        result = RangeIndex(range(1, 5, 2))
        expected = RangeIndex(1, 5, 2)
        tm.assert_index_equal(result, expected, exact=True)

    def test_constructor_range(self):
        result = RangeIndex.from_range(range(1, 5, 2))
        expected = RangeIndex(1, 5, 2)
        tm.assert_index_equal(result, expected, exact=True)

        result = RangeIndex.from_range(range(5, 6))
        expected = RangeIndex(5, 6, 1)
        tm.assert_index_equal(result, expected, exact=True)

        # an invalid range
        result = RangeIndex.from_range(range(5, 1))
        expected = RangeIndex(0, 0, 1)
        tm.assert_index_equal(result, expected, exact=True)

        result = RangeIndex.from_range(range(5))
        expected = RangeIndex(0, 5, 1)
        tm.assert_index_equal(result, expected, exact=True)

        result = Index(range(1, 5, 2))
        expected = RangeIndex(1, 5, 2)
        tm.assert_index_equal(result, expected, exact=True)

        msg = (
            r"(RangeIndex.)?from_range\(\) got an unexpected keyword argument( 'copy')?"
        )
        with pytest.raises(TypeError, match=msg):
            RangeIndex.from_range(range(10), copy=True)

    def test_constructor_name(self):
        # GH#12288
        orig = RangeIndex(10)
        orig.name = "original"

        copy = RangeIndex(orig)
        copy.name = "copy"

        assert orig.name == "original"
        assert copy.name == "copy"

        new = Index(copy)
        assert new.name == "copy"

        new.name = "new"
        assert orig.name == "original"
        assert copy.name == "copy"
        assert new.name == "new"

    def test_constructor_corner(self):
        arr = np.array([1, 2, 3, 4], dtype=object)
        index = RangeIndex(1, 5)
        assert index.values.dtype == np.int64
        expected = Index(arr).astype("int64")

        tm.assert_index_equal(index, expected, exact="equiv")

        # non-int raise Exception
        with pytest.raises(TypeError, match=r"Wrong type \<class 'str'\>"):
            RangeIndex("1", "10", "1")
        with pytest.raises(TypeError, match=r"Wrong type \<class 'float'\>"):
            RangeIndex(1.1, 10.2, 1.3)

        # invalid passed type
        with pytest.raises(
            ValueError,
            match="Incorrect `dtype` passed: expected signed integer, received float64",
        ):
            RangeIndex(1, 5, dtype="float64")


@pytest.mark.parametrize(
    "seq",
    [
        # uint64 above INT64_MAX would silently wrap to a negative range
        np.array([2**63, 2**63 + 1, 2**63 + 2], dtype=np.uint64),
        # a two-element key whose step overflows int64
        np.array([np.iinfo(np.int64).min, np.iinfo(np.int64).max], dtype=np.int64),
        # a step that overflows int64 from one extreme end
        np.array([0, np.iinfo(np.int64).max], dtype=np.int64),
        # not a range, but its elements satisfy is_sequence_range's modular
        # (uint64) check: the true value at index 2 is 2**63, which wraps to the
        # bit pattern of INT64_MIN -- must still be rejected
        np.array([0, 2**62, np.iinfo(np.int64).min], dtype=np.int64),
    ],
)
def test_maybe_sequence_to_range_extreme_ints_not_converted(seq):
    # GH#64148 don't build an incorrect range for integer sequences whose
    # conversion to a range would overflow/wrap int64
    result = maybe_sequence_to_range(seq)
    tm.assert_numpy_array_equal(result, seq)


def test_maybe_sequence_to_range_still_converts_representable():
    # GH#64148 sequences that fit in the int64 domain are still converted
    assert maybe_sequence_to_range(np.array([10, 11, 12], dtype=np.int64)) == range(
        10, 13
    )
    # a genuinely negative arithmetic progression is still converted
    assert maybe_sequence_to_range([-5, -4, -3]) == range(-5, -2)
    # uint64 values that fit in int64 are still converted
    assert maybe_sequence_to_range(np.array([1, 2, 3], dtype=np.uint64)) == range(1, 4)
    # a genuine range spanning more than INT64_MAX: the intermediate
    # ``2 * step`` overflows int64, but the range is still representable and
    # must be converted exactly (GH#64148 UBSAN signed-overflow regression)
    span = np.array([-(2**62) - 1, -1, 2**62 - 1], dtype=np.int64)
    result = maybe_sequence_to_range(span)
    assert result == range(-(2**62) - 1, 2**63 - 1, 2**62)
    assert list(result) == span.tolist()
