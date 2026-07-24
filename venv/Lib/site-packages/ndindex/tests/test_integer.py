from numpy import arange, int64, isin, bool_

from pytest import raises

from hypothesis import given, example
from hypothesis.strategies import integers, one_of

from ..integer import Integer
from ..slice import Slice
from .helpers import check_same, ints, prod, shapes, iterslice, assert_equal, reduce_kwargs

def test_integer_args():
    zero = Integer(0)
    assert zero.raw == 0
    idx = Integer(int64(0))
    assert idx == zero
    assert idx.raw == 0
    assert isinstance(idx.raw, int)
    assert Integer(zero) == zero

    raises(TypeError, lambda: Integer(1.0))
    # See the docstring of operator_index()
    raises(TypeError, lambda: Integer(True))
    raises(TypeError, lambda: Integer(bool_(True)))


    class HasIndex:
        def __init__(self, x):
            self.x = x

        def __index__(self):
            return self.x

    idx = Integer(HasIndex(0))
    assert idx.args == (0,)
    assert idx.raw == 0
    assert type(idx.args[0]) is int
    assert type(idx.raw) is int

    class HasInt:
        def __init__(self, x):
            self.x = x

        def __int__(self):
            return self.x # pragma: no cover

    raises(TypeError, lambda: Integer(HasInt(0)))

def test_integer_exhaustive():
    a = arange(10)
    for i in range(-12, 12):
        check_same(a, i)


@given(ints(), integers(5, 100))
def test_integer_hypothesis(i, size):
    a = arange(size)
    check_same(a, i)


def test_integer_len_exhaustive():
    for i in range(-12, 12):
        idx = Integer(i)
        assert len(idx) == 1


@given(ints())
def test_integer_len_hypothesis(i):
    idx = Integer(i)
    assert len(idx) == 1

def test_integer_reduce_exhaustive():
    a = arange(10)
    for i in range(-12, 12):
        for kwargs in [{'negative_int': False}, {'negative_int': True}, {}]:
            check_same(a, i, ndindex_func=lambda a, x: a[x.reduce((10,), **kwargs).raw])

            negative_int = kwargs.get('negative_int', False)

            try:
                reduced = Integer(i).reduce(10, **kwargs)
            except IndexError:
                pass
            else:
                if negative_int:
                    assert reduced.raw < 0
                else:
                    assert reduced.raw >= 0

                # Idempotency
                assert reduced.reduce(**kwargs) == reduced
                assert reduced.reduce(10, **kwargs) == reduced

@example(0, (1,), {'negative_int': True})
@given(ints(), shapes, reduce_kwargs)
def test_integer_reduce_hypothesis(i, shape, kwargs):
    a = arange(prod(shape)).reshape(shape)
    # The axis argument is tested implicitly in the Tuple.reduce test. It is
    # difficult to test here because we would have to pass in a Tuple to
    # check_same.
    check_same(a, i, ndindex_func=lambda a, x: a[x.reduce(shape, **kwargs).raw])

    negative_int = kwargs.get('negative_int', False)

    try:
        reduced = Integer(i).reduce(shape, **kwargs)
    except IndexError:
        pass
    else:
        if negative_int:
            assert reduced.raw < 0
        else:
            assert reduced.raw >= 0

        # Idempotency
        assert reduced.reduce(**kwargs) == reduced
        assert reduced.reduce(shape, **kwargs) == reduced

def test_integer_reduce_no_shape_exhaustive():
    a = arange(10)
    for i in range(-12, 12):
        check_same(a, i, ndindex_func=lambda a, x: a[x.reduce().raw])

@given(ints(), shapes, reduce_kwargs)
def test_integer_reduce_no_shape_hypothesis(i, shape, kwargs):
    a = arange(prod(shape)).reshape(shape)
    check_same(a, i, ndindex_func=lambda a, x: a[x.reduce(**kwargs).raw])

@given(ints())
def test_integer_reduce_no_shape_unchanged(i):
    idx = Integer(i)
    assert idx.reduce() == idx.reduce(negative_int=False) == idx.reduce(negative_int=True) == i

def test_integer_newshape_exhaustive():
    shape = 5
    a = arange(shape)

    def raw_func(a, idx):
        return a[idx].shape

    def ndindex_func(a, index):
        return index.newshape(shape)

    def assert_equal(raw_shape, newshape):
        assert raw_shape == newshape

    for i in range(-10, 10):
        check_same(a, i, raw_func=raw_func, ndindex_func=ndindex_func,
                   assert_equal=assert_equal)

def test_integer_as_subindex_slice_exhaustive():
    for n in range(10):
        a = arange(n)
        for i in range(-10, 10):
            try:
                a[i]
            except IndexError:
                continue

            for indexargs in iterslice():
                idx = Integer(i)

                try:
                    Index = Slice(*indexargs)
                except ValueError:
                    continue

                empty = False
                try:
                    Subindex = idx.as_subindex(Index)
                except NotImplementedError:
                    continue
                except ValueError as e:
                    assert "do not intersect" in e.args[0]
                    empty = True

                aidx = a[idx.raw]
                aindex = a[Index.raw]
                if empty:
                    assert not isin(aidx, aindex).any()
                    assert not isin(aindex, aidx).any()
                    with raises(ValueError, match="do not intersect"):
                        Index.as_subindex(idx)
                else:
                    asubindex = aindex[Subindex.raw]

                    assert_equal(asubindex.flatten(), aidx[isin(aidx, aindex)])

                    subindex2 = Index.as_subindex(idx)
                    asubindex2 = aidx[subindex2.raw]
                    assert_equal(asubindex2, asubindex)

def test_integer_isempty_exhaustive():
    for i in range(-10, 10):
        idx = Integer(i)

        isempty = idx.isempty()

        for n in range(30):
            a = arange(n)

            exception = False
            try:
                aidx = a[idx.raw]
            except IndexError:
                exception = True
            else:
                if aidx.size != 0:
                    # If a[i] doesn't give an index error, it should always be nonempty
                    assert not isempty
            # isempty() should always give the correct result for a specific
            # array shape
            try:
                isemptyn = idx.isempty(n)
            except IndexError:
                if not exception:
                    raise AssertionError(f"idx.isempty(n) raised but a[idx] did not (idx = {idx}, n = {n}).")
            else:
                if exception:
                    raise AssertionError(f"a[idx] raised but idx.isempty(n) did not (idx = {idx}, n = {n}).")
                assert isemptyn == (aidx.size == 0)

@example(1, (2, 0))
@given(ints(), one_of(shapes, integers(0, 10)))
def test_integer_isempty_hypothesis(i, shape):
    if isinstance(shape, int):
        a = arange(shape)
    else:
        a = arange(prod(shape)).reshape(shape)

    index = Integer(i)

    def raw_func(a, idx):
        return a[idx].size == 0

    def ndindex_func(a, index):
        return index.isempty(), index.isempty(shape)

    def assert_equal(raw_empty, ndindex_empty):
        isempty, isempty_shape = ndindex_empty

        # Since i is an integer, it should never be unconditionally empty
        assert not isempty
        # We cannot test the converse with hypothesis. isempty may be False
        # but a[i] could still be empty for this specific a (e.g., if a is
        # already itself empty).

        # isempty() should always give the correct result for a specific
        # array after reduction
        assert isempty_shape == raw_empty, (index, shape)

    check_same(a, i, raw_func=raw_func, ndindex_func=ndindex_func,
               assert_equal=assert_equal, same_exception=False)
