from .common import SimpleIterableType, SimpleIteratorType
from ..errors import TypingError


class RangeType(SimpleIterableType):

    def __init__(self, dtype):
        self.dtype = dtype
        name = "range_state_%s" % (dtype,)
        super(SimpleIterableType, self).__init__(name)
        self._iterator_type = RangeIteratorType(self.dtype)

    def unify(self, typingctx, other):
        if isinstance(other, RangeType):
            dtype = typingctx.unify_pairs(self.dtype, other.dtype)
            if dtype is not None:
                return RangeType(dtype)


class RangeIteratorType(SimpleIteratorType):

    def __init__(self, dtype):
        name = "range_iter_%s" % (dtype,)
        super(SimpleIteratorType, self).__init__(name)
        self._yield_type = dtype

    def unify(self, typingctx, other):
        if isinstance(other, RangeIteratorType):
            dtype = typingctx.unify_pairs(self.yield_type, other.yield_type)
            if dtype is not None:
                return RangeIteratorType(dtype)


class Generator(SimpleIteratorType):
    """
    Type class for Numba-compiled generator objects.
    """

    def __init__(self, gen_func, yield_type, arg_types, state_types,
                 has_finalizer):
        self.gen_func = gen_func
        self.arg_types = tuple(arg_types)
        self.state_types = tuple(state_types)
        self.has_finalizer = has_finalizer
        name = "%s generator(func=%s, args=%s, has_finalizer=%s)" % (
            yield_type, self.gen_func, self.arg_types,
            self.has_finalizer)
        super(Generator, self).__init__(name, yield_type)

    @property
    def key(self):
        return (self.gen_func, self.arg_types, self.yield_type,
                self.has_finalizer, self.state_types)


class EnumerateType(SimpleIteratorType):
    """
    Type class for `enumerate` objects.
    Type instances are parametered with the underlying source type.
    """

    def __init__(self, iterable_type):
        from numba.core.types import Tuple, intp
        self.source_type = iterable_type.iterator_type
        yield_type = Tuple([intp, self.source_type.yield_type])
        name = 'enumerate(%s)' % (self.source_type)
        super(EnumerateType, self).__init__(name, yield_type)

    @property
    def key(self):
        return self.source_type


class ZipType(SimpleIteratorType):
    """
    Type class for `zip` objects.
    Type instances are parametered with the underlying source types.
    """

    def __init__(self, iterable_types):
        from numba.core.types import Tuple
        self.source_types = tuple(tp.iterator_type for tp in iterable_types)
        yield_type = Tuple([tp.yield_type for tp in self.source_types])
        name = 'zip(%s)' % ', '.join(str(tp) for tp in self.source_types)
        super(ZipType, self).__init__(name, yield_type)

    @property
    def key(self):
        return self.source_types


class ArrayIterator(SimpleIteratorType):
    """
    Type class for iterators of array and buffer objects.
    """

    def __init__(self, array_type):
        self.array_type = array_type
        name = "iter(%s)" % (self.array_type,)
        nd = array_type.ndim
        if nd == 0:
            raise TypingError("iteration over a 0-d array")
        elif nd == 1:
            yield_type = array_type.dtype
        else:
            # iteration semantics leads to A order layout
            yield_type = array_type.copy(ndim=array_type.ndim - 1, layout='A')
        super(ArrayIterator, self).__init__(name, yield_type)
