from .ndindex import NDIndex, operator_index
from .shapetools import AxisError, asshape

class Integer(NDIndex):
    """
    Represents an integer index on an axis of an nd-array.

    Any object that implements `__index__` can be used as an integer index.

    >>> from ndindex import Integer
    >>> idx = Integer(1)
    >>> [0, 1, 2][idx.raw]
    1
    >>> idx = Integer(-3)
    >>> [0, 1, 2][idx.raw]
    0

    Note that `Integer` itself implements `__index__`, so it can be used as an
    index directly. However, it is still recommended to use `raw` for
    consistency, as this only works for `Integer`.

    See :doc:`../indexing-guide/integer-indices` for a description of the
    semantics of integers as indices.

    .. note::

       `Integer` does *not* represent an integer, but rather an
       *integer index*. It does not have most methods that `int` has, and
       should not be used in non-indexing contexts. See the document on
       :any:`type-confusion` for more details.

    """
    __slots__ = ()

    def _typecheck(self, idx):
        idx = operator_index(idx)
        return (idx,)

    def __index__(self):
        return self.raw

    @property
    def raw(self):
        return self.args[0]

    def __len__(self):
        """
        Returns the number of elements indexed by `self`

        Since `self` is an integer index, this always returns 1. Note that
        integer indices always remove an axis.
        """
        return 1

    def isvalid(self, shape, _axis=0):
        # The docstring for this method is on the NDIndex base class
        shape = asshape(shape)
        if not shape:
            return False
        size = shape[_axis]
        return -size <= self.raw < size

    def _raise_indexerror(self, shape, axis=0):
        if not self.isvalid(shape, axis):
            size = shape[axis]
            raise IndexError(f"index {self.raw} is out of bounds for axis {axis} with size {size}")

    def reduce(self, shape=None, *, axis=0, negative_int=False, axiserror=False):
        """
        Reduce an Integer index on an array of shape `shape`.

        The result will either be `IndexError` if the index is invalid for the
        given shape, or an Integer index where the value is nonnegative.

        If `negative_int` is `True` and a `shape` is provided, then the result
        will be an Integer index where the value is negative.

        >>> from ndindex import Integer
        >>> idx = Integer(-5)
        >>> idx.reduce((3,))
        Traceback (most recent call last):
        ...
        IndexError: index -5 is out of bounds for axis 0 with size 3
        >>> idx.reduce((9,))
        Integer(4)

        See Also
        ========

        .NDIndex.reduce
        .Tuple.reduce
        .Slice.reduce
        .ellipsis.reduce
        .Newaxis.reduce
        .IntegerArray.reduce
        .BooleanArray.reduce

        """
        if shape is None:
            return self

        if axiserror:
            if not isinstance(shape, int): # pragma: no cover
                raise TypeError("axiserror=True requires shape to be an integer")
            if not self.isvalid(shape):
                raise AxisError(self.raw, shape)

        shape = asshape(shape, axis=axis)

        self._raise_indexerror(shape, axis)

        if self.raw < 0 and not negative_int:
            size = shape[axis]
            return self.__class__(size + self.raw)
        elif self.raw >= 0 and negative_int:
            size = shape[axis]
            return self.__class__(self.raw - size)

        return self

    def newshape(self, shape):
        # The docstring for this method is on the NDIndex base class
        shape = asshape(shape)

        self._raise_indexerror(shape)
        return shape[1:]

    def as_subindex(self, index):
        index = ndindex(index)

        if isinstance(index, Tuple):
            return Tuple(self).as_subindex(index)

        if not isinstance(index, Slice):
            raise NotImplementedError("Integer.as_subindex is only implemented for slices")

        if self == -1:
            s = Slice(self.args[0], None).as_subindex(index)
        else:
            s = Slice(self.args[0], self.args[0] + 1).as_subindex(index)
        if s == Slice(0, 0, 1):
            # The intersection is empty. There is no valid index we can return
            # here. We want an index that produces an empty array, but the
            # shape should be one less, to match a[self]. Since a[index] has
            # as many dimensions as a, there is no way to index a[index] so
            # that it gives one fewer dimension but is also empty. The best we
            # could do is to return a boolean array index array([False]),
            # which would replace the first dimension with a length 0
            # dimension. But
            #
            # 1. this isn't implemented yet,
            # 2. there are complications if this happens in multiple
            #    dimensions (it might not be possible to represent, I'm not
            #    sure), and
            # 3. Slice.as_subindex(Integer) also raises this exception in the
            #    case of an empty intersection (see the comment in that code).
            raise ValueError(f"{self} and {index} do not intersect")
        assert len(s) == 1
        return Integer(s.args[0])

    def isempty(self, shape=None):
        if shape is not None:
            return 0 in self.newshape(shape)

        return False

    def selected_indices(self, shape, axis=None):
        if axis is None:
            yield from self.expand(shape).selected_indices(shape)
        else:
            shape = asshape(shape, axis=axis)
            yield self

    def __eq__(self, other):
        if isinstance(other, Integer):
            return self.args == other.args
        try:
            other = operator_index(other)
        except TypeError:
            return False
        return self.args[0] == other

    def __hash__(self):
        return super().__hash__()


# Imports at the bottom to avoid circular import issues
from .ndindex import ndindex
from .slice import Slice
from .tuple import Tuple
