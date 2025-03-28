import numpy as np
try:
    from numpy import AxisError as np_AxisError
except ImportError: # pragma: no cover
    from numpy.exceptions import AxisError as np_AxisError

from hypothesis import assume, given, example
from hypothesis.strategies import (one_of, integers, tuples as
                                   hypothesis_tuples, just, lists, shared,
                                   )

from pytest import raises

from ..ndindex import ndindex
from ..shapetools import (asshape, iter_indices, ncycles, BroadcastError,
                          AxisError, broadcast_shapes, remove_indices,
                          unremove_indices, associated_axis,
                          normalize_skip_axes)
from ..integer import Integer
from ..tuple import Tuple
from .helpers import (prod, mutually_broadcastable_shapes_with_skipped_axes,
                      skip_axes_st, mutually_broadcastable_shapes, tuples,
                      shapes, assert_equal, cross_shapes, cross_skip_axes,
                      cross_arrays_st, matmul_shapes, matmul_skip_axes,
                      matmul_arrays_st)

@example([[(1, 1), (1, 1)], (1,)], (0,))
@example([[(0,), (0,)], ()], (0,))
@example([[(1, 2), (2, 1)], (2,)], 1)
@given(mutually_broadcastable_shapes_with_skipped_axes, skip_axes_st)
def test_iter_indices(broadcastable_shapes, skip_axes):
    # broadcasted_shape will contain None on the skip_axes, as those axes
    # might not be broadcast compatible
    shapes, broadcasted_shape = broadcastable_shapes
    # We need no more than 31 dimensions so that the np.stack call below
    # doesn't fail.
    assume(len(broadcasted_shape) < 32)

    # 1. Normalize inputs
    _skip_axes = normalize_skip_axes(shapes, skip_axes)
    _skip_axes_kwarg_default = [()]*len(shapes)

    # Skipped axes may not be broadcast compatible. Since the index for a
    # skipped axis should always be a slice(None), the result should be the
    # same if the skipped axes are all moved to the end of the shape.
    canonical_shapes = []
    for s, sk in zip(shapes, _skip_axes):
        c = remove_indices(s, sk)
        canonical_shapes.append(c)

    non_skip_shapes = [remove_indices(shape, sk) for shape, sk in zip(shapes, _skip_axes)]
    assert np.broadcast_shapes(*non_skip_shapes) == broadcasted_shape

    nitems = prod(broadcasted_shape)

    if skip_axes == (): # kwarg default
        res = iter_indices(*shapes)
    else:
        res = iter_indices(*shapes, skip_axes=skip_axes)
    broadcasted_res = iter_indices(broadcasted_shape)

    sizes = [prod(shape) for shape in shapes]
    arrays = [np.arange(size).reshape(shape) for size, shape in zip(sizes, shapes)]
    canonical_sizes = [prod(shape) for shape in canonical_shapes]
    canonical_arrays = [np.arange(size).reshape(shape) for size, shape in zip(canonical_sizes, canonical_shapes)]
    canonical_broadcasted_array = np.arange(nitems).reshape(broadcasted_shape)

    # 2. Check that iter_indices is the same whether or not the shapes are
    # broadcasted together first. Also check that every iterated index is the
    # expected type and there are as many as expected.
    vals = []
    bvals = []
    n = -1

    def _remove_slices(idx):
        assert isinstance(idx, Tuple)
        idx2 = [i for i in idx.args if i != slice(None)]
        return Tuple(*idx2)

    for n, (idxes, bidxes) in enumerate(zip(res, broadcasted_res)):
        assert len(idxes) == len(shapes)
        assert len(bidxes) == 1
        for idx, shape, sk in zip(idxes, shapes, _skip_axes):
            assert isinstance(idx, Tuple)
            assert len(idx.args) == len(shape)

            for i in range(-1, -len(idx.args) - 1, -1):
                if i in sk:
                    assert idx.args[i] == slice(None)
                else:
                    assert isinstance(idx.args[i], Integer)

        canonical_idxes = [_remove_slices(idx) for idx in idxes]
        a_indexed = tuple([a[idx.raw] for a, idx in zip(arrays, idxes)])
        canonical_a_indexed = tuple([a[idx.raw] for a, idx in
                                  zip(canonical_arrays, canonical_idxes)])
        canonical_b_indexed = canonical_broadcasted_array[bidxes[0].raw]

        for c_indexed in canonical_a_indexed:
            assert c_indexed.shape == ()
        assert canonical_b_indexed.shape == ()

        if _skip_axes != _skip_axes_kwarg_default:
            vals.append(tuple(canonical_a_indexed))
        else:
            vals.append(a_indexed)

        bvals.append(canonical_b_indexed)

    # assert both iterators have the same length
    raises(StopIteration, lambda: next(res))
    raises(StopIteration, lambda: next(broadcasted_res))

    # Check that the correct number of items are iterated
    assert n == nitems - 1
    assert len(set(vals)) == len(vals) == nitems

    # 3. Check that every element of the (broadcasted) arrays is represented
    # by an iterated index.

    # The indices should correspond to the values that would be matched up
    # if the arrays were broadcasted together.
    if not arrays:
        assert vals == [()]
    else:
        correct_vals = list(zip(*[x.flat for x in np.broadcast_arrays(*canonical_arrays)]))
        # Also test that the indices are produced in a lexicographic order
        # (even though this isn't strictly guaranteed by the iter_indices
        # docstring) in the case when there are no skip axes. The order when
        # there are skip axes is more complicated because the skipped axes are
        # iterated together.
        if _skip_axes == _skip_axes_kwarg_default:
            assert vals == correct_vals
        else:
            assert set(vals) == set(correct_vals)
        assert bvals == list(canonical_broadcasted_array.flat)

@given(cross_arrays_st(), cross_shapes, cross_skip_axes)
def test_iter_indices_cross(cross_arrays, broadcastable_shapes, _skip_axes):
    # Test iter_indices behavior against np.cross, which effectively skips the
    # crossed axis. Note that we don't test against cross products of size 2
    # because a 2 x 2 cross product just returns the z-axis (i.e., it doesn't
    # actually skip an axis in the result shape), and also that behavior is
    # going to be removed in NumPy 2.0.
    a, b = cross_arrays
    shapes, broadcasted_shape = broadcastable_shapes

    # Sanity check
    skip_axes = normalize_skip_axes([*shapes, broadcasted_shape], _skip_axes)
    for sh, sk in zip([*shapes, broadcasted_shape], skip_axes):
        assert len(sk) == 1
        assert sh[sk[0]] == 3

    res = np.cross(a, b, axisa=skip_axes[0][0], axisb=skip_axes[1][0], axisc=skip_axes[2][0])
    assert res.shape == broadcasted_shape

    for idx1, idx2, idx3 in iter_indices(*shapes, broadcasted_shape, skip_axes=_skip_axes):
        assert a[idx1.raw].shape == (3,)
        assert b[idx2.raw].shape == (3,)
        assert_equal(np.cross(
            a[idx1.raw],
            b[idx2.raw]),
                     res[idx3.raw])

@given(matmul_arrays_st(), matmul_shapes, matmul_skip_axes)
def test_iter_indices_matmul(matmul_arrays, broadcastable_shapes, skip_axes):
    # Test iter_indices behavior against np.matmul, which effectively skips the
    # contracted axis (they aren't broadcasted together, even when they are
    # broadcast compatible).
    a, b = matmul_arrays
    shapes, broadcasted_shape = broadcastable_shapes

    # Note, we don't use normalize_skip_axes here because it sorts the skip
    # axes

    ax1, ax2 = skip_axes[0]
    ax3 = skip_axes[1][1]
    n, m, k = shapes[0][ax1], shapes[0][ax2], shapes[1][ax3]

    # Sanity check
    sk0, sk1, sk2 = skip_axes
    shape1, shape2 = shapes
    assert a.shape == shape1
    assert b.shape == shape2
    assert shape1[sk0[0]] == n
    assert shape1[sk0[1]] == m
    assert shape2[sk1[0]] == m
    assert shape2[sk1[1]] == k
    assert broadcasted_shape[sk2[0]] == n
    assert broadcasted_shape[sk2[1]] == k

    res = np.matmul(a, b, axes=skip_axes)
    assert res.shape == broadcasted_shape

    is_ordered = lambda sk, shape: (Integer(sk[0]).reduce(len(shape)).raw <= Integer(sk[1]).reduce(len(shape)).raw)
    orders = [
        is_ordered(sk0, shapes[0]),
        is_ordered(sk1, shapes[1]),
        is_ordered(sk2, broadcasted_shape),
    ]

    for idx1, idx2, idx3 in iter_indices(*shapes, broadcasted_shape, skip_axes=skip_axes):
        assert a[idx1.raw].shape == (n, m) if orders[0] else (m, n)
        assert b[idx2.raw].shape == (m, k) if orders[1] else (k, m)
        sub_res_axes = [
            (0, 1) if orders[0] else (1, 0),
            (0, 1) if orders[1] else (1, 0),
            (0, 1) if orders[2] else (1, 0),
        ]
        sub_res = np.matmul(a[idx1.raw], b[idx2.raw], axes=sub_res_axes)
        assert_equal(sub_res, res[idx3.raw])

def test_iter_indices_errors():
    try:
        list(iter_indices((10,), skip_axes=(2,)))
    except AxisError as e:
        ndindex_msg = str(e)
    else:
        raise RuntimeError("iter_indices did not raise AxisError") # pragma: no cover

    # Check that the message is the same one used by NumPy
    try:
        np.sum(np.arange(10), axis=2)
    except np_AxisError as e:
        np_msg = str(e)
    else:
        raise RuntimeError("np.sum() did not raise AxisError") # pragma: no cover

    assert ndindex_msg == np_msg

    try:
        list(iter_indices((2, 3), (3, 2)))
    except BroadcastError as e:
        ndindex_msg = str(e)
    else:
        raise RuntimeError("iter_indices did not raise BroadcastError") # pragma: no cover

    try:
        np.broadcast_shapes((2, 3), (3, 2))
    except ValueError as e:
        np_msg = str(e)
    else:
        raise RuntimeError("np.broadcast_shapes() did not raise ValueError") # pragma: no cover


    if 'Mismatch' in str(np_msg): # pragma: no cover
        # Older versions of NumPy do not have the more helpful error message
        assert ndindex_msg == np_msg

    with raises(ValueError, match=r"not unique"):
        list(iter_indices((1, 2), skip_axes=(0, 1, 0)))

    raises(AxisError, lambda: list(iter_indices((0,), skip_axes=(3,))))
    raises(ValueError, lambda: list(iter_indices(skip_axes=(0,))))
    raises(TypeError, lambda: list(iter_indices(1, 2)))
    raises(TypeError, lambda: list(iter_indices(1, 2, (2, 2))))
    raises(TypeError, lambda: list(iter_indices([(1, 2), (2, 2)])))

@example(1, 1, 1)
@given(integers(0, 100), integers(0, 100), integers(0, 100))
def test_ncycles(i, n, m):
    N = ncycles(range(i), n)
    if n == 1:
        assert N == range(i)
    else:
        assert isinstance(N, ncycles)
        assert N.iterable == range(i)
        assert N.n == n
        assert f"range(0, {i})" in repr(N)
        assert str(n) in repr(N)

    L = list(N)
    assert len(L) == i*n
    for j in range(i*n):
        assert L[j] == j % i

    M = ncycles(N, m)
    if n*m == 1:
        assert M == range(i)
    else:
        assert isinstance(M, ncycles)
        assert M.iterable == range(i)
        assert M.n == n*m

@given(one_of(mutually_broadcastable_shapes,
              hypothesis_tuples(tuples(shapes), just(None))))
def test_broadcast_shapes(broadcastable_shapes):
    shapes, broadcasted_shape = broadcastable_shapes
    if broadcasted_shape is not None:
        assert broadcast_shapes(*shapes) == broadcasted_shape

    arrays = [np.empty(shape) for shape in shapes]
    broadcastable = True
    try:
        broadcasted_shape = np.broadcast(*arrays).shape
    except ValueError:
        broadcastable = False

    if broadcastable:
        assert broadcast_shapes(*shapes) == broadcasted_shape
    else:
        raises(BroadcastError, lambda: broadcast_shapes(*shapes))


@given(lists(shapes, max_size=32))
def test_broadcast_shapes_errors(shapes):
    error = True
    try:
        broadcast_shapes(*shapes)
    except BroadcastError as exc:
        e = exc
    else:
        error = False

    # The ndindex and numpy errors won't match in general, because
    # ndindex.broadcast_shapes gives an error with the first two shapes that
    # aren't broadcast compatible, but numpy doesn't always, due to different
    # implementation algorithms (e.g., the message from
    # np.broadcast_shapes((0,), (0, 2), (2, 0)) mentions the last two shapes
    # whereas ndindex.broadcast_shapes mentions the first two).

    # Instead, just confirm that the error message is correct as stated, and
    # check against the numpy error message when just broadcasting the two
    # reportedly bad shapes.

    if not error:
        try:
            np.broadcast_shapes(*shapes)
        except: # pragma: no cover
            raise RuntimeError("ndindex.broadcast_shapes raised but np.broadcast_shapes did not")
        return

    assert shapes[e.arg1] == e.shape1
    assert shapes[e.arg2] == e.shape2

    try:
        np.broadcast_shapes(e.shape1, e.shape2)
    except ValueError as np_exc:
        # Check that they do in fact not broadcast, and the error messages are
        # the same modulo the different arg positions.
        if 'Mismatch' in str(np_exc): # pragma: no cover
            # Older versions of NumPy do not have the more helpful error message
            assert str(BroadcastError(0, e.shape1, 1, e.shape2)) == str(np_exc)
    else: # pragma: no cover
        raise RuntimeError("ndindex.broadcast_shapes raised but np.broadcast_shapes did not")

    raises(TypeError, lambda: broadcast_shapes(1, 2))
    raises(TypeError, lambda: broadcast_shapes(1, 2, (2, 2)))
    raises(TypeError, lambda: broadcast_shapes([(1, 2), (2, 2)]))

@given(mutually_broadcastable_shapes_with_skipped_axes, skip_axes_st)
def test_broadcast_shapes_skip_axes(broadcastable_shapes, skip_axes):
    shapes, broadcasted_shape = broadcastable_shapes
    assert broadcast_shapes(*shapes, skip_axes=skip_axes) == broadcasted_shape

@example([[], ()], (0,))
@example([[(0, 1)], (0, 1)], (2,))
@example([[(0, 1)], (0, 1)], (0, -1))
@example([[(0, 1, 0, 0, 0), (2, 0, 0, 0)], (0, 2, 0, 0, 0)], [1])
@given(mutually_broadcastable_shapes,
       one_of(
           integers(-20, 20),
           tuples(integers(-20, 20), max_size=20),
           lists(tuples(integers(-20, 20), max_size=20), max_size=32)))
def test_broadcast_shapes_skip_axes_errors(broadcastable_shapes, skip_axes):
    shapes, broadcasted_shape = broadcastable_shapes

    # All errors should come from normalize_skip_axes, which is tested
    # separately below.
    try:
        normalize_skip_axes(shapes, skip_axes)
    except (TypeError, ValueError, IndexError) as e:
        raises(type(e), lambda: broadcast_shapes(*shapes,
                                                   skip_axes=skip_axes))
        return

    try:
        broadcast_shapes(*shapes, skip_axes=skip_axes)
    except IndexError:
        raise RuntimeError("broadcast_shapes raised but should not have") # pragma: no cover
    except BroadcastError:
        # Broadcastable shapes can become unbroadcastable after skipping axes
        # (see the @example above).
        pass

remove_indices_n = shared(integers(0, 100))

@given(remove_indices_n,
       remove_indices_n.flatmap(lambda n: lists(integers(-n, n), unique=True)))
def test_remove_indices(n, idxes):
    if idxes:
        assume(max(idxes) < n)
        assume(min(idxes) >= -n)
    a = tuple(range(n))
    b = remove_indices(a, idxes)
    if len(idxes) == 1:
        assert remove_indices(a, idxes[0]) == b

    A = list(a)
    for i in idxes:
        A[i] = None

    assert set(A) - set(b) == ({None} if idxes else set())
    assert set(b) - set(A) == set()

    # Check the order is correct
    j = 0
    for i in range(n):
        val = A[i]
        if val == None:
            assert val not in b
        else:
            assert b[j] == val
            j += 1

    # Test that unremove_indices is the inverse
    if all(i >= 0 for i in idxes) or all(i < 0 for i in idxes):
        assert unremove_indices(b, idxes) == tuple(A)
    else:
        raises(NotImplementedError, lambda: unremove_indices(b, idxes))

# Meta-test for the hypothesis strategy
@given(mutually_broadcastable_shapes_with_skipped_axes, skip_axes_st)
def test_mutually_broadcastable_shapes_with_skipped_axes(broadcastable_shapes,
                                                         skip_axes): # pragma: no cover
    shapes, broadcasted_shape = broadcastable_shapes
    _skip_axes = normalize_skip_axes(shapes, skip_axes)

    assert len(_skip_axes) == len(shapes)

    for shape in shapes:
        assert None not in shape
    assert None not in broadcasted_shape

    _shapes = [remove_indices(shape, sk) for shape, sk in zip(shapes, _skip_axes)]

    assert broadcast_shapes(*_shapes) == broadcasted_shape

@example([[(2, 10, 3, 4), (10, 3, 4)], (2, 3, 4)], (-3,))
@example([[(0, 10, 2, 3, 10, 4), (1, 10, 1, 0, 10, 2, 3, 4)],
          (1, 1, 0, 2, 3, 4)], (1, 4))
@example([[(2, 0, 3, 4)], (2, 3, 4)], (1,))
@example([[(0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0)], (0, 0, 0, 0)], (1, 2))
@given(mutually_broadcastable_shapes_with_skipped_axes, skip_axes_st)
def test_associated_axis(broadcastable_shapes, skip_axes):
    shapes, broadcasted_shape = broadcastable_shapes
    _skip_axes = normalize_skip_axes(shapes, skip_axes)

    for shape, sk in zip(shapes, _skip_axes):
        n = len(shape)
        for i in range(-len(shape), 0):
            val = shape[i]

            bval = associated_axis(broadcasted_shape, i, sk)
            if bval is None:
                assert ndindex(i).reduce(n, negative_int=True) in sk, (shape, i)
            else:
                assert val == 1 or bval == val, (shape, i)


    sk = max(_skip_axes, key=len, default=())
    for i in range(-len(broadcasted_shape)-len(sk)-10, -len(broadcasted_shape)-len(sk)):
        assert associated_axis(broadcasted_shape, i, sk) is None

# TODO: add a hypothesis test for asshape
def test_asshape():
    assert asshape(1) == (1,)
    assert asshape(np.int64(2)) == (2,)
    assert type(asshape(np.int64(2))[0]) == int
    assert asshape((1, 2)) == (1, 2)
    assert asshape([1, 2]) == (1, 2)
    assert asshape((1, 2), allow_int=False) == (1, 2)
    assert asshape([1, 2], allow_int=False) == (1, 2)
    assert asshape((np.int64(1), np.int64(2))) == (1, 2)
    assert type(asshape((np.int64(1), np.int64(2)))[0]) == int
    assert type(asshape((np.int64(1), np.int64(2)))[1]) == int
    assert asshape((-1, -2), allow_negative=True) == (-1, -2)
    assert asshape(-2, allow_negative=True) == (-2,)


    raises(TypeError, lambda: asshape(1.0))
    raises(TypeError, lambda: asshape((1.0,)))
    raises(ValueError, lambda: asshape(-1))
    raises(ValueError, lambda: asshape((1, -1)))
    raises(ValueError, lambda: asshape((1, None)))
    raises(TypeError, lambda: asshape(...))
    raises(TypeError, lambda: asshape(Integer(1)))
    raises(TypeError, lambda: asshape(Tuple(1, 2)))
    raises(TypeError, lambda: asshape((True,)))
    raises(TypeError, lambda: asshape({1, 2}))
    raises(TypeError, lambda: asshape({1: 2}))
    raises(TypeError, lambda: asshape('1'))
    raises(TypeError, lambda: asshape(1, allow_int=False))
    raises(TypeError, lambda: asshape(-1, allow_int=False))
    raises(TypeError, lambda: asshape(-1, allow_negative=True, allow_int=False))
    raises(TypeError, lambda: asshape(np.int64(1), allow_int=False))
    raises(IndexError, lambda: asshape((2, 3), 3))

@example([(0,), ()], (0, 0))
@example([(0, 1), (0,), ()], [(-1,), (0,), ()])
@example([(5,)], (10,))
@example([], [])
@example([()], [])
@example([(0, 1)], 0)
@example([(2, 3), (2, 3, 4)], [(3,), (0,)])
@example([(0, 1)], 0)
@example([(2, 3)], (0, -2))
@example([(2, 4), (2, 3, 4)], [(0,), (-3,)])
@example([(1, 2)], [(0,), (1,)])
@given(lists(tuples(integers(0))),
       one_of(integers(), tuples(integers()), lists(tuples(integers()))))
def test_normalize_skip_axes(shapes, skip_axes):
    if not shapes:
        if skip_axes in [(), []]:
            assert normalize_skip_axes(shapes, skip_axes) == []
        else:
            raises(ValueError, lambda: normalize_skip_axes(shapes, skip_axes))
        return

    min_dim = min(len(shape) for shape in shapes)

    if isinstance(skip_axes, int):
        if not (-min_dim <= skip_axes < min_dim):
            raises(AxisError, lambda: normalize_skip_axes(shapes, skip_axes))
            return
        _skip_axes = [(skip_axes,)]*len(shapes)
    elif isinstance(skip_axes, tuple):
        if not all(-min_dim <= s < min_dim for s in skip_axes):
            raises(AxisError, lambda: normalize_skip_axes(shapes, skip_axes))
            return
        _skip_axes = [skip_axes]*len(shapes)
    elif not skip_axes:
        # empty list will be interpreted as a single skip_axes tuple
        assert normalize_skip_axes(shapes, skip_axes) == [()]*len(shapes)
        return
    else:
        if len(shapes) != len(skip_axes):
            raises(ValueError, lambda: normalize_skip_axes(shapes, skip_axes))
            return
        _skip_axes = skip_axes

    try:
        res = normalize_skip_axes(shapes, skip_axes)
    except AxisError as e:
        axis, ndim = e.args
        assert any(axis in s for s in _skip_axes)
        assert any(ndim == len(shape) for shape in shapes)
        assert axis < -ndim or axis >= ndim
        return
    except ValueError as e:
        if 'not unique' in str(e):
            bad_skip_axes, bad_shape = e.skip_axes, e.shape
            assert str(bad_skip_axes) in str(e)
            assert str(bad_shape) in str(e)
            assert bad_skip_axes in _skip_axes
            assert bad_shape in shapes
            indexed = [bad_shape[i] for i in bad_skip_axes]
            assert len(indexed) != len(set(indexed))
            return
        else: # pragma: no cover
            raise

    assert isinstance(res, list)
    assert all(isinstance(x, tuple) for x in res)
    assert all(isinstance(i, int) for x in res for i in x)

    assert len(res) == len(shapes)
    for shape, new_skip_axes in zip(shapes, res):
        assert len(new_skip_axes) == len(set(new_skip_axes))
        assert new_skip_axes == tuple(sorted(new_skip_axes))
        for i in new_skip_axes:
            assert i < 0
            assert ndindex(i).reduce(len(shape), negative_int=True) == i

    # TODO: Assert the order is maintained (doesn't actually matter for now
    # but could for future applications)

def test_normalize_skip_axes_errors():
    raises(TypeError, lambda: normalize_skip_axes([(1,)], {0: 1}))
    raises(TypeError, lambda: normalize_skip_axes([(1,)], {0}))
    raises(TypeError, lambda: normalize_skip_axes([(1,)], [(0,), 0]))
    raises(TypeError, lambda: normalize_skip_axes([(1,)], [0, (0,)]))

@example(10, 5)
@given(integers(), integers())
def test_axiserror(axis, ndim):
    if ndim == 0 and axis in [0, -1]:
        # NumPy allows axis=0 or -1 for 0-d arrays
        AxisError(axis, ndim)
        return

    try:
        if ndim >= 0:
            range(ndim)[axis]
    except IndexError:
        e = AxisError(axis, ndim)
    else:
        raises(ValueError, lambda: AxisError(axis, ndim))
        return

    try:
        raise e
    except AxisError as e2:
        assert e2.args == (axis, ndim)
        if ndim <= 32 and -1000 < axis < 1000:
            a = np.empty((0,)*ndim)
            try:
                np.sum(a, axis=axis)
            except np_AxisError as e3:
                assert str(e2) == str(e3)
            else:
                raise RuntimeError("numpy didn't raise AxisError") # pragma: no cover
