from hypothesis import given, example
from hypothesis.strategies import one_of, integers

from numpy import arange

from .helpers import ndindices, shapes, MAX_ARRAY_SIZE, check_same, prod

@example([0], (1,))
@example(..., (1, 2, 3))
@example(slice(0, 1), ())
@example(slice(0, 1), (1,))
@example((0, 1), (2, 2))
@example((0,), ())
@example([[1]], (0, 0, 1))
@example(None, ())
@given(ndindices, one_of(shapes, integers(0, MAX_ARRAY_SIZE)))
def test_isvalid_hypothesis(idx, shape):
    if isinstance(shape, int):
        a = arange(shape)
    else:
        a = arange(prod(shape)).reshape(shape)

    def raw_func(a, idx):
        try:
            a[idx]
            return True
        except Warning as w:
            # check_same unconditionally turns this warning into raise
            # IndexError, so we have to handle it separately here.
            if "Out of bound index found. This was previously ignored when the indexing result contained no elements. In the future the index error will be raised. This error occurs either due to an empty slice, or if an array has zero elements even before indexing." in w.args[0]:
                return False
            raise # pragma: no cover
        except IndexError:
            return False

    def ndindex_func(a, index):
        return index.isvalid(a.shape)

    def assert_equal(x, y):
        assert x == y

    check_same(a, idx, raw_func=raw_func, ndindex_func=ndindex_func,
               assert_equal=assert_equal)
