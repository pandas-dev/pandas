from abc import ABCMeta, abstractmethod
from typing import Dict as ptDict, Type as ptType
import itertools
import weakref
from functools import cached_property

import numpy as np

from numba.core.utils import get_hashable_key

# Types are added to a global registry (_typecache) in order to assign
# them unique integer codes for fast matching in _dispatcher.c.
# However, we also want types to be disposable, therefore we ensure
# each type is interned as a weak reference, so that it lives only as
# long as necessary to keep a stable type code.
# NOTE: some types can still be made immortal elsewhere (for example
# in _dispatcher.c's internal caches).
_typecodes = itertools.count()

def _autoincr():
    n = next(_typecodes)
    # 4 billion types should be enough, right?
    assert n < 2 ** 32, "Limited to 4 billion types"
    return n

_typecache: ptDict[weakref.ref, weakref.ref] = {}

def _on_type_disposal(wr, _pop=_typecache.pop):
    _pop(wr, None)


class _TypeMetaclass(ABCMeta):
    """
    A metaclass that will intern instances after they are created.
    This is done by first creating a new instance (including calling
    __init__, which sets up the required attributes for equality
    and hashing), then looking it up in the _typecache registry.
    """

    def __init__(cls, name, bases, orig_vars):
        # __init__ is hooked to mark whether a Type class being defined is a
        # Numba internal type (one which is defined somewhere under the `numba`
        # module) or an external type (one which is defined elsewhere, for
        # example a user defined type).
        super(_TypeMetaclass, cls).__init__(name, bases, orig_vars)
        root = (cls.__module__.split('.'))[0]
        cls._is_internal = root == "numba"

    def _intern(cls, inst):
        # Try to intern the created instance
        wr = weakref.ref(inst, _on_type_disposal)
        orig = _typecache.get(wr)
        orig = orig and orig()
        if orig is not None:
            return orig
        else:
            inst._code = _autoincr()
            _typecache[wr] = wr
            return inst

    def __call__(cls, *args, **kwargs):
        """
        Instantiate *cls* (a Type subclass, presumably) and intern it.
        If an interned instance already exists, it is returned, otherwise
        the new instance is returned.
        """
        inst = type.__call__(cls, *args, **kwargs)
        return cls._intern(inst)


def _type_reconstructor(reconstructor, reconstructor_args, state):
    """
    Rebuild function for unpickling types.
    """
    obj = reconstructor(*reconstructor_args)
    if state:
        obj.__dict__.update(state)
    return type(obj)._intern(obj)


class Type(metaclass=_TypeMetaclass):
    """
    The base class for all Numba types.
    It is essential that proper equality comparison is implemented.  The
    default implementation uses the "key" property (overridable in subclasses)
    for both comparison and hashing, to ensure sane behaviour.
    """

    mutable = False
    # Rather the type is reflected at the python<->nopython boundary
    reflected = False

    def __init__(self, name):
        self.name = name

    @property
    def key(self):
        """
        A property used for __eq__, __ne__ and __hash__.  Can be overridden
        in subclasses.
        """
        return self.name

    @property
    def mangling_args(self):
        """
        Returns `(basename, args)` where `basename` is the name of the type
        and `args` is a sequence of parameters of the type.

        Subclass should override to specialize the behavior.
        By default, this returns `(self.name, ())`.
        """
        return self.name, ()

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

    def __hash__(self):
        return hash(self.key)

    def __eq__(self, other):
        return self.__class__ is other.__class__ and self.key == other.key

    def __ne__(self, other):
        return not (self == other)

    def __reduce__(self):
        reconstructor, args, state = super(Type, self).__reduce__()
        return (_type_reconstructor, (reconstructor, args, state))

    def unify(self, typingctx, other):
        """
        Try to unify this type with the *other*.  A third type must
        be returned, or None if unification is not possible.
        Only override this if the coercion logic cannot be expressed
        as simple casting rules.
        """
        return None

    def can_convert_to(self, typingctx, other):
        """
        Check whether this type can be converted to the *other*.
        If successful, must return a string describing the conversion, e.g.
        "exact", "promote", "unsafe", "safe"; otherwise None is returned.
        """
        return None

    def can_convert_from(self, typingctx, other):
        """
        Similar to *can_convert_to*, but in reverse.  Only needed if
        the type provides conversion from other types.
        """
        return None

    def is_precise(self):
        """
        Whether this type is precise, i.e. can be part of a successful
        type inference.  Default implementation returns True.
        """
        return True

    def augment(self, other):
        """
        Augment this type with the *other*.  Return the augmented type,
        or None if not supported.
        """
        return None

    # User-facing helpers.  These are not part of the core Type API but
    # are provided so that users can write e.g. `numba.boolean(1.5)`
    # (returns True) or `types.int32(types.int32[:])` (returns something
    # usable as a function signature).

    def __call__(self, *args):
        from numba.core.typing import signature
        if len(args) == 1 and not isinstance(args[0], Type):
            return self.cast_python_value(args[0])
        return signature(self, # return_type
                         *args)

    def __getitem__(self, args):
        """
        Return an array of this type.
        """
        from numba.core.types import Array
        ndim, layout = self._determine_array_spec(args)
        return Array(dtype=self, ndim=ndim, layout=layout)

    def _determine_array_spec(self, args):
        # XXX non-contiguous by default, even for 1d arrays,
        # doesn't sound very intuitive
        def validate_slice(s):
            return isinstance(s, slice) and s.start is None and s.stop is None

        if isinstance(args, (tuple, list)) and all(map(validate_slice, args)):
            ndim = len(args)
            if args[0].step == 1:
                layout = 'F'
            elif args[-1].step == 1:
                layout = 'C'
            else:
                layout = 'A'
        elif validate_slice(args):
            ndim = 1
            if args.step == 1:
                layout = 'C'
            else:
                layout = 'A'
        else:
            # Raise a KeyError to not be handled by collection constructors (e.g. list).
            raise KeyError(f"Can only index numba types with slices with no start or stop, got {args}.")

        return ndim, layout

    def cast_python_value(self, args):
        raise NotImplementedError


    @property
    def is_internal(self):
        """ Returns True if this class is an internally defined Numba type by
        virtue of the module in which it is instantiated, False else."""
        return self._is_internal

    def dump(self, tab=''):
        print(f'{tab}DUMP {type(self).__name__}[code={self._code}, name={self.name}]')

# XXX we should distinguish between Dummy (no meaningful
# representation, e.g. None or a builtin function) and Opaque (has a
# meaningful representation, e.g. ExternalFunctionPointer)

class Dummy(Type):
    """
    Base class for types that do not really have a representation and are
    compatible with a void*.
    """


class Hashable(Type):
    """
    Base class for hashable types.
    """


class Number(Hashable):
    """
    Base class for number types.
    """

    def unify(self, typingctx, other):
        """
        Unify the two number types using Numpy's rules.
        """
        from numba.np import numpy_support
        if isinstance(other, Number):
            # XXX: this can produce unsafe conversions,
            # e.g. would unify {int64, uint64} to float64
            a = numpy_support.as_dtype(self)
            b = numpy_support.as_dtype(other)
            sel = np.promote_types(a, b)
            return numpy_support.from_dtype(sel)


class Callable(Type):
    """
    Base class for callables.
    """

    @abstractmethod
    def get_call_type(self, context, args, kws):
        """
        Using the typing *context*, resolve the callable's signature for
        the given arguments.  A signature object is returned, or None.
        """

    @abstractmethod
    def get_call_signatures(self):
        """
        Returns a tuple of (list of signatures, parameterized)
        """

    @abstractmethod
    def get_impl_key(self, sig):
        """
        Returns the impl key for the given signature
        """


class DTypeSpec(Type):
    """
    Base class for types usable as "dtype" arguments to various Numpy APIs
    (e.g. np.empty()).
    """

    @property
    @abstractmethod
    def dtype(self):
        """
        The actual dtype denoted by this dtype spec (a Type instance).
        """


class IterableType(Type):
    """
    Base class for iterable types.
    """

    @property
    @abstractmethod
    def iterator_type(self):
        """
        The iterator type obtained when calling iter() (explicitly or implicitly).
        """


class Sized(Type):
    """
    Base class for objects that support len()
    """


class ConstSized(Sized):
    """
    For types that have a constant size
    """
    @abstractmethod
    def __len__(self):
        pass


class IteratorType(IterableType):
    """
    Base class for all iterator types.
    Derived classes should implement the *yield_type* attribute.
    """

    def __init__(self, name, **kwargs):
        super(IteratorType, self).__init__(name, **kwargs)

    @property
    @abstractmethod
    def yield_type(self):
        """
        The type of values yielded by the iterator.
        """

    # This is a property to avoid recursivity (for pickling)

    @property
    def iterator_type(self):
        return self


class Container(Sized, IterableType):
    """
    Base class for container types.
    """


class Sequence(Container):
    """
    Base class for 1d sequence types.  Instances should have the *dtype*
    attribute.
    """


class MutableSequence(Sequence):
    """
    Base class for 1d mutable sequence types.  Instances should have the
    *dtype* attribute.
    """

    mutable = True

class ArrayCompatible(Type):
    """
    Type class for Numpy array-compatible objects (typically, objects
    exposing an __array__ method).
    Derived classes should implement the *as_array* attribute.
    """
    # If overridden by a subclass, it should also implement typing
    # for '__array_wrap__' with arguments (input, formal result).
    array_priority = 0.0

    @property
    @abstractmethod
    def as_array(self):
        """
        The equivalent array type, for operations supporting array-compatible
        objects (such as ufuncs).
        """

    # For compatibility with types.Array

    @cached_property
    def ndim(self):
        return self.as_array.ndim

    @cached_property
    def layout(self):
        return self.as_array.layout

    @cached_property
    def dtype(self):
        return self.as_array.dtype


class Literal(Type):
    """Base class for Literal types.
    Literal types contain the original Python value in the type.

    A literal type should always be constructed from the `literal(val)`
    function.
    """

    # *ctor_map* is a dictionary mapping Python types to Literal subclasses
    # for constructing a numba type for a given Python type.
    # It is used in `literal(val)` function.
    # To add new Literal subclass, register a new mapping to this dict.
    ctor_map: ptDict[type, ptType['Literal']] = {}

    # *_literal_type_cache* is used to cache the numba type of the given value.
    _literal_type_cache = None

    def __init__(self, value):
        if type(self) is Literal:
            raise TypeError(
                "Cannot be constructed directly. "
                "Use `numba.types.literal(value)` instead",
            )
        self._literal_init(value)
        fmt = "Literal[{}]({})"
        super(Literal, self).__init__(fmt.format(type(value).__name__, value))

    def _literal_init(self, value):
        self._literal_value = value
        # We want to support constants of non-hashable values, therefore
        # fall back on the value's id() if necessary.
        self._key = get_hashable_key(value)

    @property
    def literal_value(self):
        return self._literal_value

    @property
    def literal_type(self):
        if self._literal_type_cache is None:
            from numba.core import typing
            ctx = typing.Context()
            try:
                res = ctx.resolve_value_type(self.literal_value)
            except ValueError as e:

                if "Int value is too large" in str(e):
                    # If a string literal cannot create an IntegerLiteral
                    # because of overflow we generate this message.
                    msg = f"Cannot create literal type. {str(e)}"
                    raise TypeError(msg)
                # Not all literal types have a literal_value that can be
                # resolved to a type, for example, LiteralStrKeyDict has a
                # literal_value that is a python dict for which there's no
                # `typeof` support.
                msg = "{} has no attribute 'literal_type'".format(self)
                raise AttributeError(msg)
            self._literal_type_cache = res

        return self._literal_type_cache



class TypeRef(Dummy):
    """Reference to a type.

    Used when a type is passed as a value.
    """
    def __init__(self, instance_type):
        self.instance_type = instance_type
        super(TypeRef, self).__init__('typeref[{}]'.format(self.instance_type))

    @property
    def key(self):
        return self.instance_type


class InitialValue(object):
    """
    Used as a mixin for a type will potentially have an initial value that will
    be carried in the .initial_value attribute.
    """
    def __init__(self, initial_value):
        self._initial_value = initial_value

    @property
    def initial_value(self):
        return self._initial_value


class Poison(Type):
    """
    This is the "bottom" type in the type system. It won't unify and it's
    unliteral version is Poison of itself. It's advisable for debugging purposes
    to call the constructor with the type that's being poisoned (for whatever
    reason) but this isn't strictly required.
    """
    def __init__(self, ty):
        self.ty = ty
        super(Poison, self).__init__(name="Poison<%s>" % ty)

    def __unliteral__(self):
        return Poison(self)

    def unify(self, typingctx, other):
        return None
