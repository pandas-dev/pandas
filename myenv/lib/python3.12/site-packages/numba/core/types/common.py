"""
Helper classes / mixins for defining types.
"""

from .abstract import ArrayCompatible, Dummy, IterableType, IteratorType
from numba.core.errors import NumbaTypeError, NumbaValueError


class Opaque(Dummy):
    """
    A type that is a opaque pointer.
    """


class SimpleIterableType(IterableType):

    def __init__(self, name, iterator_type):
        self._iterator_type = iterator_type
        super(SimpleIterableType, self).__init__(name)

    @property
    def iterator_type(self):
        return self._iterator_type


class SimpleIteratorType(IteratorType):

    def __init__(self, name, yield_type):
        self._yield_type = yield_type
        super(SimpleIteratorType, self).__init__(name)

    @property
    def yield_type(self):
        return self._yield_type


class Buffer(IterableType, ArrayCompatible):
    """
    Type class for objects providing the buffer protocol.
    Derived classes exist for more specific cases.
    """
    mutable = True
    slice_is_copy = False
    aligned = True

    # CS and FS are not reserved for inner contig but strided
    LAYOUTS = frozenset(['C', 'F', 'CS', 'FS', 'A'])

    def __init__(self, dtype, ndim, layout, readonly=False, name=None):
        from .misc import unliteral

        if isinstance(dtype, Buffer):
            msg = ("The dtype of a Buffer type cannot itself be a Buffer type, "
                   "this is unsupported behaviour."
                   "\nThe dtype requested for the unsupported Buffer was: {}.")
            raise NumbaTypeError(msg.format(dtype))
        if layout not in self.LAYOUTS:
            raise NumbaValueError("Invalid layout '%s'" % layout)
        self.dtype = unliteral(dtype)
        self.ndim = ndim
        self.layout = layout
        if readonly:
            self.mutable = False
        if name is None:
            type_name = self.__class__.__name__.lower()
            if readonly:
                type_name = "readonly %s" % type_name
            name = "%s(%s, %sd, %s)" % (type_name, dtype, ndim, layout)
        super(Buffer, self).__init__(name)

    @property
    def iterator_type(self):
        from .iterators import ArrayIterator
        return ArrayIterator(self)

    @property
    def as_array(self):
        return self

    def copy(self, dtype=None, ndim=None, layout=None):
        if dtype is None:
            dtype = self.dtype
        if ndim is None:
            ndim = self.ndim
        if layout is None:
            layout = self.layout
        return self.__class__(dtype=dtype, ndim=ndim, layout=layout,
                              readonly=not self.mutable)

    @property
    def key(self):
        return self.dtype, self.ndim, self.layout, self.mutable

    @property
    def is_c_contig(self):
        return self.layout == 'C' or (self.ndim <= 1 and self.layout in 'CF')

    @property
    def is_f_contig(self):
        return self.layout == 'F' or (self.ndim <= 1 and self.layout in 'CF')

    @property
    def is_contig(self):
        return self.layout in 'CF'
