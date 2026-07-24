from itertools import product

from numpy import arange, array, intp, empty, all as np_all

from hypothesis import given, example
from hypothesis.strategies import integers, one_of

from pytest import raises

from ..ndindex import ndindex
from ..tuple import Tuple
from ..booleanarray import BooleanArray
from ..integer import Integer
from ..integerarray import IntegerArray
from .helpers import (check_same, Tuples, prod, short_shapes, iterslice,
                      reduce_kwargs, assert_equal_allow_scalar_0d)

def test_tuple_constructor():
    # Nested tuples are not allowed
    raises(ValueError, lambda: Tuple((1, 2, 3)))
    raises(ValueError, lambda: Tuple(0, (1, 2, 3)))
    raises(ValueError, lambda: Tuple(Tuple(1, 2, 3)))
    raises(ValueError, lambda: Tuple(0, Tuple(1, 2, 3)))

    # Multiple ellipses in a tuple are not allowed
    raises(IndexError, lambda: Tuple(..., 0, ...))
    raises(IndexError, lambda: Tuple(0, ..., ...))
    raises(IndexError, lambda: Tuple(..., ...))

    # Test NotImplementedError behavior for Tuples with arrays split up by
    # slices, ellipses, and newaxes.
    raises(NotImplementedError, lambda: Tuple(0, slice(None), [0]))
    raises(NotImplementedError, lambda: Tuple([0], slice(None), [0]))
    raises(NotImplementedError, lambda: Tuple([0], slice(None), [0]))
    raises(NotImplementedError, lambda: Tuple(0, ..., [0]))
    raises(NotImplementedError, lambda: Tuple([0], ..., [0]))
    raises(NotImplementedError, lambda: Tuple([0], ..., [0]))
    raises(NotImplementedError, lambda: Tuple(0, None, [0]))
    raises(NotImplementedError, lambda: Tuple([0], None, [0]))
    raises(NotImplementedError, lambda: Tuple([0], None, [0]))

    # Test NotImplementedError for boolean scalars mixed with other arrays
    Tuple(0, True, 0) # Doesn't raise
    raises(NotImplementedError, lambda: Tuple([0], True, 0))
    raises(NotImplementedError, lambda: Tuple(False, [0]))
    raises(NotImplementedError, lambda: Tuple(False, IntegerArray([0]), 0))
    raises(NotImplementedError, lambda: Tuple(False, 0, IntegerArray([0])))
    raises(NotImplementedError, lambda: Tuple(True, BooleanArray([True])))

    # Make sure this doesn't raise
    Tuple(0, slice(None), 0)
    Tuple(0, ..., 0)
    Tuple(0, None, 0)

def test_tuple_exhaustive():
    # Exhaustive tests here have to be very limited because of combinatorial
    # explosion.
    a = arange(2*2*2).reshape((2, 2, 2))
    types = {
        slice: lambda: iterslice((-1, 1), (-1, 1), (-1, 1), one_two_args=False),
        # slice: _iterslice,
        int: lambda: ((i,) for i in range(-3, 3)),
        type(...): lambda: ()
    }

    for t1, t2, t3 in product(types, repeat=3):
        for t1_args in types[t1]():
            for t2_args in types[t2]():
                for t3_args in types[t3]():
                    idx1 = t1(*t1_args)
                    idx2 = t2(*t2_args)
                    idx3 = t3(*t3_args)

                    idx = (idx1, idx2, idx3)
                    # Disable the same exception check because there could be
                    # multiple invalid indices in the tuple, and for instance
                    # numpy may give an IndexError but we would give a
                    # TypeError because we check the type first.
                    check_same(a, idx, same_exception=False)
                    try:
                        index = Tuple(*idx)
                    except (IndexError, ValueError):
                        pass
                    else:
                        assert index.has_ellipsis == (type(...) in (t1, t2, t3))

@given(Tuples, short_shapes)
def test_tuples_hypothesis(t, shape):
    a = arange(prod(shape)).reshape(shape)
    check_same(a, t, same_exception=False)

@given(Tuples, short_shapes)
def test_ellipsis_index(t, shape):
    a = arange(prod(shape)).reshape(shape)
    # Don't know if there is a better way to test ellipsis_idx
    def ndindex_func(a, index):
        return a[ndindex((*index.raw[:index.ellipsis_index], ...,
                          *index.raw[index.ellipsis_index+1:])).raw]

    check_same(a, t, ndindex_func=ndindex_func, assert_equal=assert_equal_allow_scalar_0d)

@example((True, 0, False), 1, {})
@example((..., None), (), {})
@given(Tuples, one_of(short_shapes, integers(0, 10)), reduce_kwargs)
def test_tuple_reduce_no_shape_hypothesis(t, shape, kwargs):
    if isinstance(shape, int):
        a = arange(shape)
    else:
        a = arange(prod(shape)).reshape(shape)

    index = Tuple(*t)

    check_same(a, index.raw, ndindex_func=lambda a, x:
               a[x.reduce(**kwargs).raw], same_exception=False,
               assert_equal=assert_equal_allow_scalar_0d)

    reduced = index.reduce(**kwargs)
    if isinstance(reduced, Tuple):
        assert len(reduced.args) != 1
        assert reduced == () or reduced.args[-1] != ...

    # Idempotency
    assert reduced.reduce(**kwargs) == reduced

@example((..., empty((1, 0), dtype=intp)), (1, 0), {})
@example((1, -1, [1, -1]), (3, 3, 3), {'negative_int': True})
@example((..., None), (), {})
@example((..., empty((0, 0), dtype=bool)), (0, 0), {})
@example((empty((0, 0), dtype=bool), 0), (0, 0, 1), {})
@example((array([], dtype=intp), 0), (0, 0), {})
@example((array([], dtype=intp), array(0)), (0, 0), {})
@example((array([], dtype=intp), [0]), (0, 0), {})
@example((0, 1, ..., 2, 3), (2, 3, 4, 5, 6, 7), {})
@example((0, slice(None), ..., slice(None), 3), (2, 3, 4, 5, 6, 7), {})
@example((0, ..., slice(None)), (2, 3, 4, 5, 6, 7), {})
@example((slice(None, None, -1),), (2,), {})
@example((..., slice(None, None, -1),), (2, 3, 4), {})
@example((..., False, slice(None)), 0, {})
@given(Tuples, one_of(short_shapes, integers(0, 10)), reduce_kwargs)
def test_tuple_reduce_hypothesis(t, shape, kwargs):
    if isinstance(shape, int):
        a = arange(shape)
    else:
        a = arange(prod(shape)).reshape(shape)

    index = Tuple(*t)

    check_same(a, index.raw, ndindex_func=lambda a, x: a[x.reduce(shape,
                                                                  **kwargs).raw], same_exception=False,
               assert_equal=assert_equal_allow_scalar_0d)

    negative_int = kwargs.get('negative_int', False)

    try:
        reduced = index.reduce(shape, **kwargs)
    except IndexError:
        pass
    else:
        if isinstance(reduced, Tuple):
            assert len(reduced.args) != 1
            assert reduced == () or reduced.args[-1] != ...
        # TODO: Check the other properties from the Tuple.reduce docstring.

        # Idempotency
        assert reduced.reduce(**kwargs) == reduced
        # This is currently not implemented, for example, (..., False, :)
        # takes two steps to remove the redundant slice.
        # assert reduced.reduce(shape) == reduced

        for arg in reduced.args:
            if isinstance(arg, Integer):
                if negative_int:
                    assert arg.raw < 0
                else:
                    assert arg.raw >= 0
            elif isinstance(arg, IntegerArray):
                if negative_int:
                    assert np_all(arg.raw < 0)
                else:
                    assert np_all(arg.raw >= 0)

def test_tuple_reduce_explicit():
    # Some aspects of Tuple.reduce are hard to test as properties, so include
    # some explicit tests here.

    # (Before Index, shape): After index
    tests = {
        # Make sure redundant slices are removed
        (Tuple(0, ..., slice(0, 3)), (5, 3)): Integer(0),
        (Tuple(slice(0, 5), ..., 0), (5, 3)): Tuple(..., Integer(0)),
        # Ellipsis is removed if unnecessary
        (Tuple(0, ...), (2, 3)): Integer(0),
        (Tuple(0, 1, ...), (2, 3)): Tuple(Integer(0), Integer(1)),
        (Tuple(..., 0, 1), (2, 3)): Tuple(Integer(0), Integer(1)),
    }

    for (before, shape), after in tests.items():
        reduced = before.reduce(shape)
        assert reduced == after, (before, shape)

        a = arange(prod(shape)).reshape(shape)
        check_same(a, before.raw, ndindex_func=lambda a, x:
                   a[x.reduce(shape).raw],
                   assert_equal=assert_equal_allow_scalar_0d)

        # Idempotency
        assert reduced.reduce() == reduced
        assert reduced.reduce(shape) == reduced

@example((slice(0, 0),), 2)
@example((0, slice(0, 0)), (1, 2))
@given(Tuples, one_of(short_shapes, integers(0, 10)))
def test_tuple_isempty_hypothesis(t, shape):
    if isinstance(shape, int):
        a = arange(shape)
    else:
        a = arange(prod(shape)).reshape(shape)

    T = Tuple(*t)

    try:
        ndindex(t).isempty(shape)
    except NotImplementedError:
        return
    except IndexError:
        pass

    def raw_func(a, t):
        return a[t].size == 0

    def ndindex_func(a, T):
        return T.isempty(), T.isempty(shape)

    def assert_equal(raw_empty, ndindex_empty):
        isempty, isempty_shape = ndindex_empty

        # If isempty is True then a[t] should be empty
        if isempty:
            assert raw_empty, (T, shape)
        # We cannot test the converse with hypothesis. isempty may be False
        # but a[t] could still be empty for this specific a (e.g., if a is
        # already itself empty).

        # If isempty is true with no shape it should be true for a specific
        # shape. The converse is not true because the indexed array could be
        # empty.
        if isempty:
            assert isempty_shape, (T, shape)

        # isempty() should always give the correct result for a specific
        # array after reduction
        assert isempty_shape == raw_empty, (T, shape)

    check_same(a, t, raw_func=raw_func, ndindex_func=ndindex_func,
               assert_equal=assert_equal, same_exception=False)
