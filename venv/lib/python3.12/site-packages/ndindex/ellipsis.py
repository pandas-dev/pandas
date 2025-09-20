from .ndindex import NDIndex
from .tuple import Tuple
from .shapetools import asshape

class ellipsis(NDIndex):
    """
    Represents an ellipsis index, i.e., `...` (or `Ellipsis`).

    Ellipsis indices by themselves return the full array. Inside of a tuple
    index, an ellipsis skips 0 or more axes of the array so that everything
    after the ellipsis indexes the last axes of the array. A tuple index can
    have at most one ellipsis.

    See :doc:`../indexing-guide/multidimensional-indices/ellipses` for more
    details on the semantics of ellipsis indices.

    For example `a[(0, ..., -2)]` would index the first element on the first
    axis, the second-to-last element in the last axis, and include all the
    axes in between.

    >>> from numpy import arange
    >>> a = arange(2*3*4).reshape((2, 3, 4))
    >>> a
    array([[[ 0,  1,  2,  3],
            [ 4,  5,  6,  7],
            [ 8,  9, 10, 11]],
           [[12, 13, 14, 15],
            [16, 17, 18, 19],
            [20, 21, 22, 23]]])
    >>> a[0, ..., -2]
    array([ 2,  6, 10])

    An ellipsis can go at the beginning of end of a tuple index, and is
    allowed to match 0 axes.

    .. note::

       Unlike the standard Python `Ellipsis`, `ellipsis` is the type, not the
       object (the name is lowercase to avoid conflicting with the built-in).
       Use `ellipsis()` or `ndindex(...)` to create the object. In most
       ndindex contexts, `...` can be used instead of `ellipsis()`, for
       instance, when creating a `Tuple` object. Also unlike `Ellipsis`,
       `ellipsis()` is not singletonized, so you should not use `is` to
       compare it. See the document on :any:`type confusion
       <type-confusion-ellipsis>` for more details.

    """
    __slots__ = ()

    def _typecheck(self):
        return ()

    def reduce(self, shape=None, *, negative_int=False):
        """
        Reduce an ellipsis index

        Since an ellipsis by itself always returns the full array unchanged,
        `ellipsis().reduce()` returns `Tuple()` as a canonical form (the index
        `()` also always returns an array unchanged).

        >>> from ndindex import ellipsis
        >>> ellipsis().reduce()
        Tuple()

        See Also
        ========

        .NDIndex.reduce
        .Tuple.reduce
        .Slice.reduce
        .Newaxis.reduce
        .Integer.reduce
        .IntegerArray.reduce
        .BooleanArray.reduce

        """
        if shape is not None:
            shape = asshape(shape)
        return Tuple()

    @property
    def raw(self):
        return ...

    def isvalid(self, shape):
        shape = asshape(shape)
        return True

    def newshape(self, shape):
        # The docstring for this method is on the NDIndex base class
        shape = asshape(shape)

        return shape

    def as_subindex(self, index):
        return Tuple().as_subindex(index)

    def isempty(self, shape=None):
        return Tuple().isempty(shape=shape)

    def __eq__(self, other):
        return other is ... or isinstance(other, ellipsis)

    def __hash__(self):
        return super().__hash__()
