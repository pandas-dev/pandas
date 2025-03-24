from numpy import prod, arange, array, bool_, empty, full, __version__ as np_version

NP1 = np_version.startswith('1')

from hypothesis import given, example
from hypothesis.strategies import one_of, integers

from pytest import raises

from .helpers import boolean_arrays, short_shapes, check_same, assert_equal, reduce_kwargs

from ..booleanarray import BooleanArray

def test_booleanarray_constructor():
    raises(ValueError, lambda: BooleanArray([False], shape=(1,)))
    raises(ValueError, lambda: BooleanArray([], shape=(1,)))
    raises(TypeError, lambda: BooleanArray([0]))
    raises(TypeError, lambda: BooleanArray(array(0.0)))
    raises(TypeError, lambda: BooleanArray((True,)))
    idx = BooleanArray(array([True], dtype=bool_))
    assert_equal(idx.array, array([True], dtype=bool_))

    idx = BooleanArray([], shape=(0, 1))
    assert_equal(idx.array, empty((0, 1), dtype=bool_))

    # Make sure the underlying array is immutable
    idx = BooleanArray([True])
    with raises(ValueError):
        idx.array[0] = False
    assert_equal(idx.array, array([True], dtype=bool_))

    # Make sure the underlying array is copied
    a = array([True, False])
    idx = BooleanArray(a)
    a[0] = False
    assert idx == BooleanArray([True, False])

@given(boolean_arrays, short_shapes)
def test_booleanarray_hypothesis(idx, shape):
    a = arange(prod(shape)).reshape(shape)
    check_same(a, idx)

@given(boolean_arrays, one_of(short_shapes, integers(0, 10)), reduce_kwargs)
def test_booleanarray_reduce_no_shape_hypothesis(idx, shape, kwargs):
    if isinstance(shape, int):
        a = arange(shape)
    else:
        a = arange(prod(shape)).reshape(shape)

    index = BooleanArray(idx)

    check_same(a, index.raw, ndindex_func=lambda a, x: a[x.reduce(**kwargs).raw])

@example(full((1, 9), True), 3, {})
@example(full((1, 9), True), (3, 3), {})
@example(full((1, 9), False), (3, 3), {})
@given(boolean_arrays, one_of(short_shapes, integers(0, 10)), reduce_kwargs)
def test_booleanarray_reduce_hypothesis(idx, shape, kwargs):
    if isinstance(shape, int):
        a = arange(shape)
    else:
        a = arange(prod(shape)).reshape(shape)

    index = BooleanArray(idx)

    same_exception = not NP1
    check_same(a, index.raw, ndindex_func=lambda a, x: a[x.reduce(shape, **kwargs).raw],
               same_exception=same_exception)

    try:
        reduced = index.reduce(shape, **kwargs)
    except IndexError:
        pass
    else:
        # At present, reduce() always returns the same index if it doesn't
        # give an IndexError
        assert reduced == index

        # Idempotency
        assert reduced.reduce(**kwargs) == reduced
        assert reduced.reduce(shape, **kwargs) == reduced

@given(boolean_arrays, one_of(short_shapes, integers(0, 10)))
def test_booleanarray_isempty_hypothesis(idx, shape):
    if isinstance(shape, int):
        a = arange(shape)
    else:
        a = arange(prod(shape)).reshape(shape)

    index = BooleanArray(idx)

    def raw_func(a, idx):
        return a[idx].size == 0

    def ndindex_func(a, index):
        return index.isempty(), index.isempty(shape)

    def assert_equal(raw_empty, ndindex_empty):
        isempty, isempty_shape = ndindex_empty

        # If isempty is True then a[t] should be empty
        if isempty:
            assert raw_empty, (index, shape)
        # We cannot test the converse with hypothesis. isempty may be False
        # but a[idx] could still be empty for this specific a (e.g., if a is
        # already itself empty).

        # If isempty is true with no shape it should be true for a specific
        # shape. The converse is not true because the indexed array could be
        # empty.
        if isempty:
            assert isempty_shape, (index, shape)

        # isempty() should always give the correct result for a specific
        # array after reduction
        assert isempty_shape == raw_empty, (index, shape)

    check_same(a, idx, raw_func=raw_func, ndindex_func=ndindex_func,
               assert_equal=assert_equal, same_exception=False)
