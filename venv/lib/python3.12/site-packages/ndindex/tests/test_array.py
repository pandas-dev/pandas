from numpy import intp, array, int16

from pytest import raises

from ..array import ArrayIndex
from ..integerarray import IntegerArray

from .helpers import assert_equal

# Everything else is tested in the subclasses

def test_ArrayIndex():
    raises(TypeError, lambda: ArrayIndex([]))

def test_attributes():
    a = array([[0, 1], [1, 0]])

    idx = IntegerArray(a)
    assert_equal(idx.array, array(a, dtype=intp))
    assert idx.dtype == intp
    assert idx.ndim == a.ndim == 2
    assert idx.shape == a.shape == (2, 2)
    assert idx.size == a.size == 4

def test_cast_raises():
    with raises(TypeError):
        a = array([[0, 1], [1, 0]])
        idx = IntegerArray(a)
        assert array(idx) == idx

def test_copy():
    idx = IntegerArray([1, 2])
    idx2 = IntegerArray(idx.raw)
    assert idx.raw is not idx2.raw
    idx3 = IntegerArray(idx.raw, _copy=False)
    assert idx.raw is idx3.raw
    raises(ValueError, lambda: IntegerArray([], _copy=False))
    raises(ValueError, lambda: IntegerArray(array([1], dtype=int16), _copy=False))
