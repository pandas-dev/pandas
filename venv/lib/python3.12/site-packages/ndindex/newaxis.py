from .ndindex import NDIndex
from .shapetools import asshape

class Newaxis(NDIndex):
    """
    Represents a `np.newaxis` (i.e., `None`) index.

    `Newaxis` adds a shape 1 dimension to the array. If a `Newaxis` is inside
    of a tuple index, it adds a shape 1 dimension at that location in the
    index.

    For example, if `a` has shape `(2, 3)`, then `a[newaxis]` has shape `(1,
    2, 3)`, `a[:, newaxis]` has shape `(2, 1, 3)`, and so on.

    >>> from ndindex import Newaxis
    >>> from numpy import arange
    >>> a = arange(0,6).reshape(2,3)
    >>> a[Newaxis().raw].shape
    (1, 2, 3)
    >>> a[:, Newaxis().raw, :].shape
    (2, 1, 3)

    Using `Newaxis().raw` as an index is equivalent to using `numpy.newaxis`.

    See :doc:`../indexing-guide/multidimensional-indices/newaxis` for a
    description of the semantics of newaxis.

    .. note::

       Unlike the NumPy `newaxis`, `Newaxis` is the type, not the object (the
       name is uppercase to avoid conflicting with the NumPy type). Use
       `Newaxis()`, `ndindex(np.newaxis)`, or `ndindex(None)` to create the
       object. In most ndindex contexts, `np.newaxis` or `None` can be used
       instead of `Newaxis()`, for instance, when creating a `Tuple` object.
       Also unlike `None`, `Newaxis()` is not singletonized, so you should not
       use `is` to compare it. See the document on :any:`type-confusion` for
       more details.

    """
    __slots__ = ()

    def _typecheck(self):
        return ()

    @property
    def raw(self):
        return None

    def reduce(self, shape=None, *, axis=0, negative_int=False):
        """
        Reduce a `Newaxis` index

        There is no other index that is equivalent to a newaxis index by
        itself, so `Newaxis().reduce()` always returns `Newaxis()` unchanged.

        >>> from ndindex import Newaxis
        >>> Newaxis().reduce()
        Newaxis()

        See Also
        ========

        .NDIndex.reduce
        .Tuple.reduce
        .Slice.reduce
        .Integer.reduce
        .ellipsis.reduce
        .IntegerArray.reduce
        .BooleanArray.reduce

        """
        if shape is not None:
            shape = asshape(shape)
        return self

    def isvalid(self, shape):
        shape = asshape(shape)
        return True

    def newshape(self, shape):
        # The docstring for this method is on the NDIndex base class
        shape = asshape(shape)

        # reduce will raise IndexError if it should be raised
        self.reduce(shape)

        return (1,) + shape

    def isempty(self, shape=None):
        if shape is not None:
            return 0 in self.newshape(shape)

        return False

    def __eq__(self, other):
        return other is None or isinstance(other, Newaxis)

    def __hash__(self):
        return super().__hash__()
