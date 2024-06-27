from pytest import raises

from numpy import arange, array, full

from hypothesis import given, example
from hypothesis.strategies import integers, one_of

from ..ndindex import ndindex
from ..tuple import Tuple
from ..integer import Integer
from .helpers import ndindices, check_same, short_shapes, prod

@example(..., 0)
@example((True,), ())
@example(([[True, False], [True, False]], [True, True], slice(0, 2)), ((2, 2, 2, 3, 3)))
@example((array([], dtype=bool),), (0, 0))
@example((False, False), ())
@example(array([], dtype=bool), 0)
@example((array([], dtype=bool),), 0)
@example(array([[[True], [False]]]), (1, 1, 2))
@example(full((1, 9), False), (3, 3))
@example(([0, 1], 0), (2, 2))
@example(([0, 0, 0], [0, 0]), (2, 2))
@example((0, None, 0, ..., 0, None, 0), (2, 2, 2, 2, 2, 2, 2))
@example((0, slice(None), ..., slice(None), 3), (2, 3, 4, 5, 6, 7))
@given(ndindices, one_of(short_shapes, integers(0, 10)))
def test_newshape_hypothesis(idx, shape):
    if isinstance(shape, int):
        a = arange(shape)
    else:
        a = arange(prod(shape)).reshape(shape)

    try:
        index = ndindex(idx)
    except IndexError:
        pass
    else:
        # Make sure ndindex input gives an error
        raises(TypeError, lambda: index.newshape(Tuple(2, 1)))
        raises(TypeError, lambda: index.newshape(Integer(2)))

    def raw_func(a, idx):
        return a[idx].shape

    def ndindex_func(a, index):
        return index.newshape(shape)

    def assert_equal(raw_shape, newshape):
        assert raw_shape == newshape

    check_same(a, idx, raw_func=raw_func, ndindex_func=ndindex_func,
               assert_equal=assert_equal, same_exception=False)
