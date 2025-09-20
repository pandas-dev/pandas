import numbers
import itertools
from collections.abc import Sequence
from ._crt import prod

from .ndindex import ndindex, operator_index

class BroadcastError(ValueError):
    """
    Exception raised by :func:`iter_indices()` and
    :func:`broadcast_shapes()` when the input shapes are not broadcast
    compatible.

    """
    def __init__(self, arg1, shape1, arg2, shape2):
        self.arg1 = arg1
        self.shape1 = shape1
        self.arg2 = arg2
        self.shape2 = shape2

    def __str__(self):
        arg1, shape1, arg2, shape2 = self.args
        return f"shape mismatch: objects cannot be broadcast to a single shape.  Mismatch is between arg {arg1} with shape {shape1} and arg {arg2} with shape {shape2}."

class AxisError(ValueError, IndexError):
    """
    Exception raised by :func:`iter_indices()` and
    :func:`broadcast_shapes()` when the `skip_axes` argument is out of bounds.

    This is used instead of the NumPy exception of the same name so that
    `iter_indices` does not need to depend on NumPy.

    """
    def __init__(self, axis, ndim):
        # NumPy allows axis=-1 for 0-d arrays
        if (ndim < 0 or -ndim <= axis < ndim) and not (ndim == 0 and axis == -1):
            raise ValueError(f"Invalid AxisError ({axis}, {ndim})")
        self.axis = axis
        self.ndim = ndim

    def __str__(self):
        return f"axis {self.axis} is out of bounds for array of dimension {self.ndim}"

def broadcast_shapes(*shapes, skip_axes=()):
    """
    Broadcast the input shapes `shapes` to a single shape.

    This is the same as :external+numpy:py:func:`np.broadcast_shapes()
    <numpy.broadcast_shapes>`, except is also supports skipping axes in the
    shape with `skip_axes`.

    `skip_axes` can be a tuple of integers which apply to all shapes, or a
    list of tuples of integers, one for each shape, which apply to each
    respective shape. The `skip_axes` argument works the same as in
    :func:`iter_indices`. See its docstring for more details.

    If the shapes are not broadcast compatible (excluding `skip_axes`),
    :class:`BroadcastError` is raised.

    >>> from ndindex import broadcast_shapes
    >>> broadcast_shapes((2, 3), (3,), (4, 2, 1))
    (4, 2, 3)
    >>> broadcast_shapes((2, 3), (5,), (4, 2, 1))
    Traceback (most recent call last):
    ...
    ndindex.shapetools.BroadcastError: shape mismatch: objects cannot be broadcast to a single shape.  Mismatch is between arg 0 with shape (2, 3) and arg 1 with shape (5,).

    Axes in `skip_axes` apply to each shape *before* being broadcasted. Each
    shape will be broadcasted together with these axes removed. The dimensions
    in `skip_axes` do not need to be equal or broadcast compatible with one
    another. The final broadcasted shape be the result of broadcasting all the
    non-skip axes.

    >>> broadcast_shapes((10, 3, 2), (2, 20), skip_axes=[(0,), (1,)])
    (3, 2)

    """
    shapes = [asshape(shape, allow_int=False) for shape in shapes]
    skip_axes = normalize_skip_axes(shapes, skip_axes)

    if not shapes:
        return ()

    non_skip_shapes = [remove_indices(shape, skip_axis) for shape, skip_axis in zip(shapes, skip_axes)]
    dims = [len(shape) for shape in non_skip_shapes]
    N = max(dims)

    broadcasted_shape = [1]*N

    arg = None
    for i in range(-1, -N-1, -1):
        for j in range(len(shapes)):
            if dims[j] < -i:
                continue
            shape = non_skip_shapes[j]
            broadcasted_side = broadcasted_shape[i]
            shape_side = shape[i]
            if shape_side == 1:
                continue
            elif broadcasted_side == 1:
                broadcasted_side = shape_side
                arg = j
            elif shape_side != broadcasted_side:
                raise BroadcastError(arg, shapes[arg], j, shapes[j])
            broadcasted_shape[i] = broadcasted_side

    return tuple(broadcasted_shape)

def iter_indices(*shapes, skip_axes=(), _debug=False):
    """
    Iterate indices for every element of an arrays of shape `shapes`.

    Each shape in `shapes` should be a shape tuple, which are broadcast
    compatible along the non-skipped axes. Each iteration step will produce a
    tuple of indices, one for each shape, which would correspond to the same
    elements if the arrays of the given shapes were first broadcast together.

    This is a generalization of the NumPy :external+numpy:py:class:`np.ndindex()
    <numpy.ndindex>` function (which otherwise has no relation).
    `np.ndindex()` only iterates indices for a single shape, whereas
    `iter_indices()` supports generating indices for multiple broadcast
    compatible shapes at once. This is equivalent to first broadcasting the
    arrays then generating indices for the single broadcasted shape.

    Additionally, this function supports the ability to skip axes of the
    shapes using `skip_axes`. These axes will be fully sliced in each index.
    The remaining axes will be indexed one element at a time with integer
    indices.

    `skip_axes` should be a tuple of axes to skip or a list of tuples of axes
    to skip. If it is a single tuple, it applies to all shapes. Otherwise,
    each tuple applies to each shape respectively. It can use negative
    integers, e.g., `skip_axes=(-1,)` will skip the last axis. The order of
    the axes in `skip_axes` does not matter. Mixing negative and nonnegative
    skip axes is supported, but the skip axes must refer to unique dimensions
    for each shape.

    The axes in `skip_axes` refer to the shapes *before* broadcasting (if you
    want to refer to the axes after broadcasting, either broadcast the shapes
    and arrays first, or refer to the axes using negative integers). For
    example, `iter_indices((10, 2), (20, 1, 2), skip_axes=(0,))` will skip the
    size `10` axis of `(10, 2)` and the size `20` axis of `(20, 1, 2)`. The
    result is two sets of indices, one for each element of the non-skipped
    dimensions:

    >>> from ndindex import iter_indices
    >>> for idx1, idx2 in iter_indices((10, 2), (20, 1, 2), skip_axes=(0,)):
    ...     print(idx1, idx2)
    Tuple(slice(None, None, None), 0) Tuple(slice(None, None, None), 0, 0)
    Tuple(slice(None, None, None), 1) Tuple(slice(None, None, None), 0, 1)

    The skipped axes do not themselves need to be broadcast compatible, but
    the shapes with all the skipped axes removed should be broadcast
    compatible.

    For example, suppose `a` is an array with shape `(3, 2, 4, 4)`, which we
    wish to think of as a `(3, 2)` stack of 4 x 4 matrices. We can generate an
    iterator for each matrix in the "stack" with `iter_indices((3, 2, 4, 4),
    skip_axes=(-1, -2))`:

    >>> for idx in iter_indices((3, 2, 4, 4), skip_axes=(-1, -2)):
    ...     print(idx)
    (Tuple(0, 0, slice(None, None, None), slice(None, None, None)),)
    (Tuple(0, 1, slice(None, None, None), slice(None, None, None)),)
    (Tuple(1, 0, slice(None, None, None), slice(None, None, None)),)
    (Tuple(1, 1, slice(None, None, None), slice(None, None, None)),)
    (Tuple(2, 0, slice(None, None, None), slice(None, None, None)),)
    (Tuple(2, 1, slice(None, None, None), slice(None, None, None)),)

    .. note::

       The iterates of `iter_indices` are always a tuple, even if only a
       single shape is provided (one could instead use `for idx, in
       iter_indices(...)` above).

    As another example, say `a` is shape `(1, 3)` and `b` is shape `(2, 1)`,
    and we want to generate indices for every value of the broadcasted
    operation `a + b`. We can do this by using `a[idx1.raw] + b[idx2.raw]` for every
    `idx1` and `idx2` as below:

    >>> import numpy as np
    >>> a = np.arange(3).reshape((1, 3))
    >>> b = np.arange(100, 111, 10).reshape((2, 1))
    >>> a
    array([[0, 1, 2]])
    >>> b
    array([[100],
           [110]])
    >>> for idx1, idx2 in iter_indices((1, 3), (2, 1)):
    ...     print(f"{idx1 = }; {idx2 = }; {(a[idx1.raw], b[idx2.raw]) = }") # doctest: +SKIPNP1
    idx1 = Tuple(0, 0); idx2 = Tuple(0, 0); (a[idx1.raw], b[idx2.raw]) = (np.int64(0), np.int64(100))
    idx1 = Tuple(0, 1); idx2 = Tuple(0, 0); (a[idx1.raw], b[idx2.raw]) = (np.int64(1), np.int64(100))
    idx1 = Tuple(0, 2); idx2 = Tuple(0, 0); (a[idx1.raw], b[idx2.raw]) = (np.int64(2), np.int64(100))
    idx1 = Tuple(0, 0); idx2 = Tuple(1, 0); (a[idx1.raw], b[idx2.raw]) = (np.int64(0), np.int64(110))
    idx1 = Tuple(0, 1); idx2 = Tuple(1, 0); (a[idx1.raw], b[idx2.raw]) = (np.int64(1), np.int64(110))
    idx1 = Tuple(0, 2); idx2 = Tuple(1, 0); (a[idx1.raw], b[idx2.raw]) = (np.int64(2), np.int64(110))
    >>> a + b
    array([[100, 101, 102],
           [110, 111, 112]])

    To include an index into the final broadcasted array, you can simply
    include the final broadcasted shape as one of the shapes (the function
    :func:`broadcast_shapes` is useful here).

    >>> np.broadcast_shapes((1, 3), (2, 1))
    (2, 3)
    >>> for idx1, idx2, broadcasted_idx in iter_indices((1, 3), (2, 1), (2, 3)):
    ...     print(broadcasted_idx)
    Tuple(0, 0)
    Tuple(0, 1)
    Tuple(0, 2)
    Tuple(1, 0)
    Tuple(1, 1)
    Tuple(1, 2)

    """
    skip_axes = normalize_skip_axes(shapes, skip_axes)
    shapes = [asshape(shape, allow_int=False) for shape in shapes]

    if not shapes:
        yield ()
        return

    shapes = [asshape(shape) for shape in shapes]
    S = len(shapes)

    iters = [[] for i in range(S)]
    broadcasted_shape = broadcast_shapes(*shapes, skip_axes=skip_axes)

    idxes = [-1]*S

    while any(i is not None for i in idxes):
        for s, it, shape, sk in zip(range(S), iters, shapes, skip_axes):
            i = idxes[s]
            if i is None:
                continue
            if -i > len(shape):
                if not shape:
                    pass
                elif len(shape) == len(sk):
                    # The whole shape is skipped. Just repeat the most recent slice
                    it[0] = ncycles(it[0], prod(broadcasted_shape))
                else:
                    # Find the first non-skipped axis and repeat by however
                    # many implicit axes are left in the broadcasted shape
                    for j in range(-len(shape), 0):
                        if j not in sk:
                            break
                    it[j] = ncycles(it[j], prod(broadcasted_shape[:len(sk)-len(shape)+len(broadcasted_shape)]))

                idxes[s] = None
                continue

            val = associated_axis(broadcasted_shape, i, sk)
            if i in sk:
                it.insert(0, [slice(None)])
            else:
                if val == 0:
                    return
                elif val != 1 and shape[i] == 1:
                    it.insert(0, ncycles(range(shape[i]), val))
                else:
                    it.insert(0, range(shape[i]))
            idxes[s] -= 1

    if _debug: # pragma: no cover
        print(f"{iters = }")
    for idxes in itertools.zip_longest(*[itertools.product(*i) for i in
                                         iters], fillvalue=()):
        yield tuple(ndindex(idx) for idx in idxes)

#### Internal helpers


def asshape(shape, axis=None, *, allow_int=True, allow_negative=False):
    """
    Cast `shape` as a valid NumPy shape.

    The input can be an integer `n` (if `allow_int=True`), which is equivalent
    to `(n,)`, or a tuple of integers.

    If the `axis` argument is provided, an `IndexError` is raised if it is out
    of bounds for the shape.

    The resulting shape is always a tuple of nonnegative integers. If
    `allow_negative=True`, negative integers are also allowed.

    All ndindex functions that take a shape input should use::

        shape = asshape(shape)

    or::

        shape = asshape(shape, axis=axis)

    """
    from .integer import Integer
    from .tuple import Tuple
    if isinstance(shape, (Tuple, Integer)):
        raise TypeError("ndindex types are not meant to be used as a shape - "
                        "did you mean to use the built-in tuple type?")

    if isinstance(shape, numbers.Number):
        if allow_int:
            shape = (operator_index(shape),)
        else:
            raise TypeError(f"expected sequence of integers, not {type(shape).__name__}")

    if not isinstance(shape, Sequence) or isinstance(shape, str):
        raise TypeError("expected sequence of integers" + allow_int*" or a single integer" + ", not " + type(shape).__name__)
    l = len(shape)

    newshape = []
    # numpy uses __getitem__ rather than __iter__ to index into shape, so we
    # match that
    for i in range(l):
        # Raise TypeError if invalid
        val = shape[i]
        if val is None:
            raise ValueError("unknonwn (None) dimensions are not supported")

        newshape.append(operator_index(shape[i]))

        if not allow_negative and val < 0:
            raise ValueError("unknown (negative) dimensions are not supported")

    if axis is not None:
        if len(newshape) <= axis:
            raise IndexError(f"too many indices for array: array is {len(shape)}-dimensional, but {axis + 1} were indexed")

    return tuple(newshape)

def associated_axis(broadcasted_shape, i, skip_axes):
    """
    Return the associated element of `broadcasted_shape` corresponding to
    `shape[i]` given `skip_axes`. If there is not such element (i.e., it's out
    of bounds), returns None.

    This function makes implicit assumptions about its input and is only
    designed for internal use.

    """
    skip_axes = sorted(skip_axes, reverse=True)
    if i >= 0:
        raise NotImplementedError
    if i in skip_axes:
        return None
    # We assume skip_axes are all negative and sorted
    j = i
    for sk in skip_axes:
        if sk >= i:
            j += 1
        else:
            break
    if ndindex(j).isvalid(len(broadcasted_shape)):
        return broadcasted_shape[j]
    return None

def remove_indices(x, idxes):
    """
    Return `x` with the indices `idxes` removed.

    This function is only intended for internal usage.
    """
    if isinstance(idxes, int):
        idxes = (idxes,)
    dim = len(x)
    _idxes = sorted({i if i >= 0 else i + dim for i in idxes})
    _idxes = [i - a for i, a in zip(_idxes, range(len(_idxes)))]
    _x = list(x)
    for i in _idxes:
        _x.pop(i)
    return tuple(_x)

def unremove_indices(x, idxes, *, val=None):
    """
    Insert `val` in `x` so that it appears at `idxes`.

    Note that idxes must be either all negative or all nonnegative.

    This function is only intended for internal usage.
    """
    if any(i >= 0 for i in idxes) and any(i < 0 for i in idxes):
        # A mix of positive and negative indices presents a fundamental
        # problem: sometimes the result is not unique. For example, x = [0];
        # idxes = [1, -1] could be satisfied by both [0, None] or [0, None,
        # None], depending on whether each index refers to a separate None or
        # not (note that both cases are supported by remove_indices(), because
        # there it is unambiguous). But even worse, in some cases, there may
        # be no way to satisfy the given requirement. For example, given x =
        # [0, 1, 2, 3]; idxes = [3, -3], there is no way to insert None into x
        # so that remove_indices(res, idxes) == x. To see this, simply observe
        # that there is no size list x such that remove_indices(x, [3, -3])
        # returns a tuple of size 4:
        #
        # >>> [len(remove_indices(list(range(n)), [3, -3])) for n in range(4, 10)]
        # [2, 3, 5, 5, 6, 7]
        raise NotImplementedError("Mixing both negative and nonnegative idxes is not yet supported")
    x = list(x)
    n = len(idxes) + len(x)
    _idxes = sorted({i if i >= 0 else i + n for i in idxes})
    for i in _idxes:
        x.insert(i, val)
    return tuple(x)

class ncycles:
    """
    Iterate `iterable` repeated `n` times.

    This is based on a recipe from the `Python itertools docs
    <https://docs.python.org/3/library/itertools.html#itertools-recipes>`_,
    but improved to give a repr, and to denest when it can. This makes
    debugging :func:`~.iter_indices` easier.

    This is only intended for internal usage.

    >>> from ndindex.shapetools import ncycles
    >>> ncycles(range(3), 2)
    ncycles(range(0, 3), 2)
    >>> list(_)
    [0, 1, 2, 0, 1, 2]
    >>> ncycles(ncycles(range(3), 3), 2)
    ncycles(range(0, 3), 6)

    """
    def __new__(cls, iterable, n):
        if n == 1:
            return iterable
        return object.__new__(cls)

    def __init__(self, iterable, n):
        if isinstance(iterable, ncycles):
            self.iterable = iterable.iterable
            self.n = iterable.n*n
        else:
            self.iterable = iterable
            self.n = n

    def __repr__(self):
        return f"ncycles({self.iterable!r}, {self.n!r})"

    def __iter__(self):
        return itertools.chain.from_iterable(itertools.repeat(tuple(self.iterable), self.n))

def normalize_skip_axes(shapes, skip_axes):
    """
    Return a canonical form of `skip_axes` corresponding to `shapes`.

    A canonical form of `skip_axes` is a list of tuples of integers, one for
    each shape in `shapes`, which are a unique set of axes for each
    corresponding shape.

    If `skip_axes` is an integer, this is basically `[(skip_axes,) for s
    in shapes]`. If `skip_axes` is a tuple, it is like `[skip_axes for s in
    shapes]`.

    The `skip_axes` must always refer to unique axes in each shape.

    The returned `skip_axes` will always be negative integers and will be
    sorted.

    This function is only intended for internal usage.

    """
    # Note: we assume asshape has already been called on the shapes in shapes
    if isinstance(skip_axes, Sequence):
        if skip_axes and all(isinstance(i, Sequence) for i in skip_axes):
            if len(skip_axes) != len(shapes):
                raise ValueError(f"Expected {len(shapes)} skip_axes")
            return [normalize_skip_axes([shape], skip_axis)[0] for shape, skip_axis in zip(shapes, skip_axes)]
        else:
            try:
                [operator_index(i) for i in skip_axes]
            except TypeError:
                raise TypeError("skip_axes must be an integer, a tuple of integers, or a list of tuples of integers")

    skip_axes = asshape(skip_axes, allow_negative=True)

    # From here, skip_axes is a single tuple of integers

    if not shapes and skip_axes:
        raise ValueError("skip_axes must be empty if there are no shapes")

    new_skip_axes = []
    err = None
    for shape in shapes:
        s = tuple(sorted(ndindex(i).reduce(len(shape), negative_int=True, axiserror=True).raw for i in skip_axes))
        if len(s) != len(set(s)):
            err = ValueError(f"skip_axes {skip_axes} are not unique for shape {shape}")
            # For testing
            err.skip_axes = skip_axes
            err.shape = shape
        new_skip_axes.append(s)
    if err:
        raise err
    return new_skip_axes
