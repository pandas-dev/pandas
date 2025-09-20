from .array import ArrayIndex
from .shapetools import asshape
from .subindex_helpers import subindex_slice

class IntegerArray(ArrayIndex):
    """
    Represents an integer array index.

    If `idx` is an n-dimensional integer array with shape `s = (s1, ..., sn)`
    and `a` is any array, `a[idx]` replaces the first dimension of `a` with
    dimensions of size `s1, ..., sn`, where each entry is indexed according to
    the entry in `idx` as an integer index.

    Integer arrays can also appear as part of tuple indices. In that case,
    they replace the axis being indexed. If more than one integer array
    appears inside of a tuple index, they are broadcast together and iterated
    as one. Furthermore, if an integer array appears in a tuple index, all
    integer indices in the tuple are treated as scalar integer arrays and are
    also broadcast. In general, an :any:`Integer` index semantically behaves
    the same as a scalar (`shape=()`) `IntegerArray`.

    A list (or list of lists) of integers may also be used in place of an
    array.

    See :doc:`../indexing-guide/multidimensional-indices/integer-arrays` for a
    description of the semantics of integer array indices.

    >>> from ndindex import IntegerArray
    >>> import numpy as np
    >>> idx = IntegerArray([[0, 1], [1, 2]])
    >>> a = np.arange(10)
    >>> a[idx.raw]
    array([[0, 1],
           [1, 2]])

    .. note::

       `IntegerArray` does *not* represent an array, but rather an *array
       index*. It does not have most methods that `numpy.ndarray` has, and
       should not be used in array contexts. See the document on
       :any:`type-confusion` for more details.

    """
    __slots__ = ()

    @property
    def dtype(self):
        """
        The dtype of `IntegerArray` is `np.intp`, which is typically either
        `np.int32` or `np.int64` depending on the platform.
        """
        from numpy import intp
        return intp

    def _raise_indexerror(self, shape, axis=0):
        size = shape[axis]
        out_of_bounds = (self.array >= size) | ((-size > self.array) & (self.array < 0))
        if out_of_bounds.any():
            raise IndexError(f"index {self.array[out_of_bounds].flat[0]} is out of bounds for axis {axis} with size {size}")

    def reduce(self, shape=None, *, axis=0, negative_int=False):
        """
        Reduce an `IntegerArray` index on an array of shape `shape`.

        The result will either be `IndexError` if the index is invalid for the
        given shape, an `IntegerArray` index where the values are all
        nonnegative, or, if `self` is a scalar array index (`self.shape ==
        ()`), an `Integer` whose value is nonnegative.

        If `negative_int` is `True` and a `shape` is provided, the result will
        be an `IntegerArray` with negative entries instead of positive
        entries.

        >>> from ndindex import IntegerArray
        >>> idx = IntegerArray([-5, 2])
        >>> idx.reduce((3,))
        Traceback (most recent call last):
        ...
        IndexError: index -5 is out of bounds for axis 0 with size 3
        >>> idx.reduce((9,))
        IntegerArray([4, 2])
        >>> idx.reduce((9,), negative_int=True)
        IntegerArray([-5, -7])

        See Also
        ========

        .NDIndex.reduce
        .Tuple.reduce
        .Slice.reduce
        .ellipsis.reduce
        .Newaxis.reduce
        .Integer.reduce
        .BooleanArray.reduce

        """
        if self.shape == ():
            return Integer(self.array).reduce(shape, axis=axis, negative_int=negative_int)

        if shape is None:
            return self

        shape = asshape(shape, axis=axis)

        self._raise_indexerror(shape, axis)

        size = shape[axis]
        new_array = self.array.copy()
        if negative_int:
            new_array[new_array >= 0] -= size
        else:
            new_array[new_array < 0] += size
        return IntegerArray(new_array)

    def newshape(self, shape):
        # The docstring for this method is on the NDIndex base class
        shape = asshape(shape)

        self._raise_indexerror(shape)
        return self.shape + shape[1:]

    def isempty(self, shape=None):
        if shape is not None:
            return 0 in self.newshape(shape)

        return 0 in self.shape

    def as_subindex(self, index):
        index = ndindex(index)

        if isinstance(index, Tuple):
            return Tuple(self).as_subindex(index)

        if isinstance(index, Slice):
            index = index.reduce()
            if (self.array < 0).any():
                raise NotImplementedError("IntegerArray.as_subindex() is only implemented for arrays with all nonnegative entries. Try calling reduce() with a shape first.")
            if index.step < 0:
                raise NotImplementedError("IntegerArray.as_subindex(Slice) is only implemented for slices with positive steps")

            # After reducing, start is not None when step > 0
            if index.stop is None or index.start < 0 or index.stop < 0:
                raise NotImplementedError("IntegerArray.as_subindex(Slice) is only implemented for slices with nonnegative start and stop. Try calling reduce() with a shape first.")

            # Equivalent to

            # res = []
            # for i in self.array.flat:
            #     s = Slice(i, i+1).as_subindex(index)
            #     if s == Slice(0, 0, 1):
            #         continue
            #     res.append(s.start)
            # return IntegerArray(res)

            # See also Integer.as_subindex().

            s = self.array
            start, stop, step = subindex_slice(s, s+1, 1,
                                               index.start, index.stop, index.step)
            if (stop <= 0).all():
                raise ValueError("Indices do not intersect")

            if start.shape == ():
                if start < stop:
                    return IntegerArray(start)
                else:
                    raise ValueError("Indices do not intersect")

            mask = start < stop
            if not mask.all():
                start = start[start < stop]

            if 0 in start.shape:
                raise ValueError("Indices do not intersect")
            return IntegerArray(start)

        raise NotImplementedError("IntegerArray.as_subindex is only implemented for slices")

    def selected_indices(self, shape, axis=None):
        if axis is None:
            yield from self.expand(shape).selected_indices(shape)
        else:
            shape = asshape(shape, axis=axis)
            for i in self.array.flat:
                yield Integer(i)

    def __eq__(self, other):
        from numpy import ndarray

        if isinstance(other, IntegerArray):
            b = other.array
        elif isinstance(other, ndarray):
            b = other
        elif isinstance(other, list):
            try:
                b = IntegerArray(other)
            except TypeError:
                return False
        else:
            return False
        a = self.array
        return a.shape == b.shape and (a == b).all()

    def __hash__(self):
        return super().__hash__()

# Imports at the bottom to avoid circular import issues
from .ndindex import ndindex
from .slice import Slice
from .tuple import Tuple
from .integer import Integer
