import itertools
import operator

import numpy as np
import pytest

from pandas.errors import Pandas4Warning

from pandas.core.dtypes.missing import isna

import pandas as pd
import pandas._testing as tm
from pandas.core.ops.array_ops import (
    comparison_op,
    na_logical_op,
)


def test_na_logical_op_2d():
    left = np.arange(8).reshape(4, 2)
    right = left.astype(object)
    right[0, 0] = np.nan

    # Check that we fall back to the vec_binop branch
    with pytest.raises(TypeError, match="unsupported operand type"):
        operator.or_(left, right)

    result = na_logical_op(left, right, operator.or_)
    expected = right
    tm.assert_numpy_array_equal(result, expected)


def test_object_comparison_2d():
    left = np.arange(9).reshape(3, 3).astype(object)
    right = left.T

    result = comparison_op(left, right, operator.eq)
    expected = np.eye(3).astype(bool)
    tm.assert_numpy_array_equal(result, expected)

    # Ensure that cython doesn't raise on non-writeable arg, which
    #  we can get from np.broadcast_to
    right.flags.writeable = False
    result = comparison_op(left, right, operator.ne)
    tm.assert_numpy_array_equal(result, ~expected)


@pytest.mark.parametrize("rvalues", [1, [1, 1, 1], np.nan, None])
@pytest.mark.parametrize(
    "op", [operator.eq, operator.ne, operator.lt, operator.le, operator.gt, operator.ge]
)
def test_comparison_for_subclasses(rvalues, op):
    # GH#63205 Ensure subclasses of ndarray are correctly handled in comparison_op
    # Define a custom ndarray subclass
    class TestArray(np.ndarray):
        def __new__(cls, input_array):
            return np.asarray(input_array).view(cls)

        def __array_finalize__(self, obj) -> None:
            self._is_test_array = True

    def expected_with_na_handling(lvalues, rvalues, op):
        # Similar to comparison_op, handle zerodim arrays with na value separately
        if (rvalues.ndim == 0) and isna(rvalues.item()):
            # numpy does not like comparisons vs None
            if op is operator.ne:
                return np.ones(lvalues.shape, dtype=bool)
            else:
                return np.zeros(lvalues.shape, dtype=bool)
        return op(lvalues, rvalues)

    # Define test data
    lvalues = [1, 2, 3]

    # Test with both ndarray and TestArray
    result = comparison_op(np.array(lvalues), np.array(rvalues), op)
    expected = expected_with_na_handling(np.array(lvalues), np.array(rvalues), op)
    tm.assert_numpy_array_equal(result, expected)

    result = comparison_op(TestArray(lvalues), TestArray(rvalues), op)
    expected = expected_with_na_handling(TestArray(lvalues), TestArray(rvalues), op)
    tm.assert_numpy_array_equal(result, expected)


# the last three yield 1, 2, 3 -- the values several of the arrays below hold --
#  so element-wise treatment would give all-True there, scalar-like all-False
ITERATOR_BOXES = [
    pytest.param(lambda: (num for num in range(10)), id="generator"),
    pytest.param(lambda: iter([1, 2, 3]), id="iterator"),
    pytest.param(lambda: map(int, "123"), id="map"),
    pytest.param(lambda: reversed([3, 2, 1]), id="reversed"),
]

ARRAYS = [
    pd.array([1, 2, 3], dtype="Int64"),
    pd.array([1.0, 2.0, 3.0], dtype="Float64"),
    pd.array([True, False, True]),
    pd.Categorical([1, 2, 3]),
    pd.array(pd.date_range("2020", periods=3)),
    pd.array(pd.to_timedelta([1, 2, 3], unit="D")),
    pd.array(pd.period_range("2020", periods=3, freq="D")),
    pd.arrays.IntervalArray.from_breaks([1, 2, 3, 4]),
    pd.arrays.SparseArray([1, 2, 3]),
    pd.array(["a", "b", "c"], dtype=pd.StringDtype("python")),
    pd.Series([1, 2, 3]),
    pd.Index([1, 2, 3]),
    np.array([1, 2, 3]),
]


@pytest.mark.parametrize("box", ITERATOR_BOXES)
@pytest.mark.parametrize("arr", ARRAYS)
def test_cmp_iterator_treated_as_scalar(arr, box):
    # GH#31646 an iterator has no length, so we cannot compare element-wise;
    #  it is treated as scalar-like, matching ndarray behavior
    # the shape matters as much as the values: a scalar False would satisfy
    #  .any() too, and that is the failure mode the scalar path can produce
    expected = np.zeros(len(arr), dtype=bool)
    with tm.assert_produces_warning(None):
        result = arr == box()
    tm.assert_numpy_array_equal(np.asarray(result), expected)

    with tm.assert_produces_warning(None):
        result = arr != box()
    tm.assert_numpy_array_equal(np.asarray(result), ~expected)


@pytest.mark.parametrize("box", ITERATOR_BOXES)
@pytest.mark.parametrize(
    "arr",
    [
        pd.array(pd.date_range("2020", periods=3)),
        pd.array(pd.to_timedelta([1, 2, 3], unit="D")),
        pd.arrays.IntervalArray.from_breaks([1, 2, 3, 4]),
        pd.arrays.SparseArray([1, 2, 3]),
        pd.array([1, 2, 3], dtype="Int64"),
        pd.array([1.0, 2.0, 3.0], dtype="Float64"),
    ],
)
def test_arith_iterator_raises(arr, box):
    # GH#31646 scalar-like treatment means the operation is simply invalid.
    #  Only sparse and the masked dtypes change here: the datetimelike and
    #  interval boxes already raised on main, from Python rather than pandas
    with pytest.raises(TypeError, match="unsupported operand type"):
        arr + box()


TIMEDELTA_OPS = [
    pytest.param(operator.mul, id="mul"),
    pytest.param(operator.truediv, id="truediv"),
    pytest.param(operator.floordiv, id="floordiv"),
    pytest.param(operator.mod, id="mod"),
    pytest.param(divmod, id="divmod"),
    pytest.param(lambda arr, other: other / arr, id="rtruediv"),
    pytest.param(lambda arr, other: other // arr, id="rfloordiv"),
]


@pytest.mark.parametrize("box", ITERATOR_BOXES)
@pytest.mark.parametrize("op", TIMEDELTA_OPS)
def test_timedelta_arith_iterator_raises(op, box):
    # GH#31646 the multiplicative timedelta ops gated on lib.is_scalar, so an
    #  iterator fell through to np.array(other) and then len() on the 0-d result
    arr = pd.array(pd.to_timedelta([1, 2, 3], unit="D"))
    msg = "|".join(["cannot use operands", "unsupported operand", "Cannot divide"])
    with pytest.raises(TypeError, match=msg):
        op(arr, box())


@pytest.mark.parametrize("box", ITERATOR_BOXES)
def test_frame_cmp_iterator_treated_as_scalar(box):
    # GH#31646 the DataFrame alignment path must not consume the iterator
    #  looking for a length either
    df = pd.DataFrame({"a": [1, 2, 3], "b": pd.array([4, 5, 6], dtype="Int64")})
    with tm.assert_produces_warning(None):
        result = df == box()
    assert not result.to_numpy().any()

    with tm.assert_produces_warning(None):
        result = df.ne(box())
    assert result.to_numpy().all()


class ReIterable:
    """``__iter__`` but no ``__next__`` and no ``__len__``: not an iterator."""

    def __iter__(self):
        return iter([1, 2, 3])


def test_reiterable_is_not_treated_as_an_iterator():
    # GH#31646 only iterators became scalar-like; a re-readable list-like with
    #  no length is still on the GH#62423 deprecation path
    with tm.assert_produces_warning(Pandas4Warning, match="is deprecated"):
        result = pd.arrays.SparseArray([1, 2, 3]) == ReIterable()
    assert np.asarray(result).all()


class SizedSequence:
    """``__len__`` + ``__getitem__``, no ``__iter__``: ``is_list_like`` is False."""

    def __init__(self, values):
        self.values = list(values)

    def __len__(self) -> int:
        return len(self.values)

    def __getitem__(self, key):
        return self.values[key]


class ArrayCastable:
    """``__array__`` but no ``__len__`` and no ``__iter__``."""

    def __array__(self, dtype=None, copy=None):
        return np.array([1, 2, 3])


class ArrayInterface:
    """``__array_interface__`` only -- e.g. a ``PIL.Image``."""

    def __init__(self) -> None:
        self._values = np.array([1, 2, 3])
        self.__array_interface__ = self._values.__array_interface__


@pytest.mark.parametrize("cls", [ArrayCastable, ArrayInterface])
def test_array_castable_operand_still_elementwise(cls):
    # GH#31646 has_castable_attr exists to keep these element-wise; an operand
    #  with no __len__ must not be mistaken for a scalar
    with tm.assert_produces_warning(None):
        result = pd.arrays.SparseArray([1, 2, 3]) + cls()
    # check_dtype: the sparse op upcasts to int64 on 32-bit platforms
    tm.assert_numpy_array_equal(
        np.asarray(result), np.array([2, 4, 6]), check_dtype=False
    )

    arr = pd.array(pd.to_timedelta([1, 2, 3], unit="D"))
    result = arr * cls()
    tm.assert_equal(result, pd.array(pd.to_timedelta([1, 4, 9], unit="D")))

    # a SparseArray whose fill value participates: the scalar path would build
    #  the fill from the whole operand and then fail to broadcast
    sparse = pd.arrays.SparseArray([0, 1, 2]) == cls()
    tm.assert_numpy_array_equal(np.asarray(sparse), np.zeros(3, dtype=bool))
    assert sparse.fill_value


def test_sized_sequence_still_elementwise():
    # GH#31646 the op sites that used to gate on lib.is_scalar must keep sending
    #  a sized sequence down the element-wise path; NumPy coerces it anyway
    with tm.assert_produces_warning(Pandas4Warning, match="is deprecated"):
        result = pd.arrays.SparseArray([1, 2, 3]) + SizedSequence([1, 2, 3])
    tm.assert_numpy_array_equal(
        np.asarray(result), np.array([2, 4, 6]), check_dtype=False
    )

    result = pd.arrays.SparseArray([1, 2, 3]) == SizedSequence([1, 2, 3])
    assert np.asarray(result).all()

    result = pd.array([True, False, True]) & SizedSequence([True, True, True])
    tm.assert_extension_array_equal(result, pd.array([True, False, True]))

    arr = pd.array(pd.to_timedelta([1, 2, 3], unit="D"))
    result = arr * SizedSequence([1, 2, 3])
    tm.assert_equal(result, pd.array(pd.to_timedelta([1, 4, 9], unit="D")))


@pytest.mark.parametrize("box", ITERATOR_BOXES)
@pytest.mark.parametrize("dtype", ["int64[pyarrow]", "string[pyarrow]"])
def test_arrow_arith_does_not_consume_iterator(dtype, box):
    # GH#31646 _evaluate_op_method boxed the iterator as an array, draining it --
    #  and never returning for an endless one -- instead of treating it as scalar
    pa = pytest.importorskip("pyarrow")
    data = ["a", "b", "c"] if dtype.startswith("string") else [1, 2, 3]
    arr = pd.array(data, dtype=dtype)
    other = box()
    with pytest.raises((pa.ArrowInvalid, TypeError)):
        arr + other
    assert list(other), "the iterator was consumed instead of treated as a scalar"


def test_arrow_logical_op_iterator_raises_arrow_invalid():
    # GH#31646 an iterator now reaches _evaluate_op_method rather than
    #  logical_op's dtype-less-sequence gate, so the arrow scalar path reports
    #  it the way it reports any other unboxable scalar
    pa = pytest.importorskip("pyarrow")
    arr = pd.array([True, False, True], dtype="bool[pyarrow]")
    with pytest.raises(pa.ArrowInvalid, match="Could not convert"):
        arr & (num for num in range(3))


@pytest.mark.parametrize("box", ITERATOR_BOXES)
@pytest.mark.parametrize("dtype", ["int64[pyarrow]", "string[pyarrow]"])
def test_arrow_cmp_iterator_treated_as_scalar(dtype, box):
    # GH#31646 _cmp_method gated on lib.is_scalar, so an iterator matched no arm
    #  and raised NotImplementedError while every other dtype compared False
    pytest.importorskip("pyarrow")
    data = ["a", "b", "c"] if dtype.startswith("string") else [1, 2, 3]
    arr = pd.array(data, dtype=dtype)
    other = box()
    tm.assert_numpy_array_equal(np.asarray(arr == other), np.zeros(3, dtype=bool))
    assert list(other), "the iterator was consumed instead of treated as a scalar"

    # an ordering comparison against a scalar it cannot compare to still raises
    with pytest.raises(TypeError, match="Invalid comparison"):
        arr < box()


def test_arrow_arith_endless_iterator_raises():
    # GH#31646 an endless iterator must raise rather than being consumed forever
    pa = pytest.importorskip("pyarrow")
    arr = pd.array([1, 2, 3], dtype="int64[pyarrow]")
    with pytest.raises(pa.ArrowInvalid, match="Could not convert"):
        arr + itertools.count()


@pytest.mark.parametrize("box", ITERATOR_BOXES)
def test_logical_op_iterator_reaches_the_operand_message(box):
    # GH#31646 logical_op gated on bare is_list_like, so an iterator was called a
    #  "dtype-less sequence" and never reached the scalar path the rest of the
    #  operators put it on
    msg = "'other' should be pandas.NA or a bool"
    with pytest.raises(TypeError, match=msg):
        pd.array([True, False, True]) & box()
    with pytest.raises(TypeError, match=msg):
        pd.Series([True, False, True], dtype="boolean") & box()

    # numpy-backed reports it as the scalar it now is, not as a sequence
    with pytest.raises(TypeError, match="Cannot perform 'and_'"):
        pd.Series([True, False, True]) & box()

    # object dtype reaches na_logical_op's non-ndarray branch, where bool(other)
    #  would have answered True instead of raising
    with pytest.raises(TypeError, match=r"scalar of type \[\w+\]"):
        pd.Series([True, False, True], dtype=object) | box()


@pytest.mark.parametrize("box", ITERATOR_BOXES)
def test_sparse_logical_and_ordering_iterator_raises(box):
    # GH#31646 sparse computed elementwise against a length-matching iterator
    #  for these two as well, and `|` did not even return booleans
    with pytest.raises(TypeError, match="unsupported operand type"):
        pd.arrays.SparseArray([True, False, True]) | box()

    with pytest.raises(TypeError, match="not supported between instances"):
        pd.arrays.SparseArray([1, 2, 3]) > box()


def test_sparse_cmp_unrecognized_scalar():
    # GH#31646 an object that is neither list-like nor a recognized scalar was
    #  routed through the list-like branch and compared as a length-1 operand
    arr = pd.arrays.SparseArray([1, 2, 3])
    result = arr == object()
    assert not np.asarray(result).any()

    with pytest.raises(TypeError, match="unsupported operand type"):
        arr + object()


def test_boolean_logical_unrecognized_scalar():
    # GH#31646 lib.is_scalar is narrower than "not list-like", so an ordinary
    #  object fell through to the element-wise branch and hit len()
    arr = pd.array([True, False, True])
    with pytest.raises(TypeError, match="'other' should be pandas.NA or a bool"):
        arr & object()


def test_object_dtype_logical_unrecognized_scalar():
    # GH#31646 an unrecognized object reached bool(other) in na_logical_op and
    #  silently computed against True; it raised a bare AssertionError before
    ser = pd.Series([True, False, True], dtype=object)
    with pytest.raises(TypeError, match=r"scalar of type \[object\]"):
        ser & object()

    # the message names the operand, not the bool a recognized scalar is cast to
    with pytest.raises(TypeError, match=r"scalar of type \[str\]"):
        pd.Series([1, 2]) | "x"
