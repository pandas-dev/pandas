from numba.core.types.abstract import Callable, Literal, Type, Hashable
from numba.core.types.common import (Dummy, IterableType, Opaque,
                                     SimpleIteratorType)
from numba.core.typeconv import Conversion
from numba.core.errors import TypingError, LiteralTypingError
from numba.core.ir import UndefinedType
from numba.core.utils import get_hashable_key


class PyObject(Dummy):
    """
    A generic CPython object.
    """

    def is_precise(self):
        return False


class Phantom(Dummy):
    """
    A type that cannot be materialized.  A Phantom cannot be used as
    argument or return type.
    """


class Undefined(Dummy):
    """
    A type that is left imprecise.  This is used as a temporaray placeholder
    during type inference in the hope that the type can be later refined.
    """

    def is_precise(self):
        return False


class UndefVar(Dummy):
    """
    A type that is created by Expr.undef to represent an undefined variable.
    This type can be promoted to any other type.
    This is introduced to handle Python 3.12 LOAD_FAST_AND_CLEAR.
    """

    def can_convert_to(self, typingctx, other):
        return Conversion.promote


class RawPointer(Opaque):
    """
    A raw pointer without any specific meaning.
    """


class StringLiteral(Literal, Dummy):

    def can_convert_to(self, typingctx, other):
        if isinstance(other, UnicodeType):
            return Conversion.safe


Literal.ctor_map[str] = StringLiteral


def unliteral(lit_type):
    """
    Get base type from Literal type.
    """
    if hasattr(lit_type, '__unliteral__'):
        return lit_type.__unliteral__()
    return getattr(lit_type, 'literal_type', lit_type)


def literal(value):
    """Returns a Literal instance or raise LiteralTypingError
    """
    ty = type(value)
    if isinstance(value, Literal):
        msg = "the function does not accept a Literal type; got {} ({})"
        raise ValueError(msg.format(value, ty))
    try:
        ctor = Literal.ctor_map[ty]
    except KeyError:
        raise LiteralTypingError("{} cannot be used as a literal".format(ty))
    else:
        return ctor(value)


def maybe_literal(value):
    """Get a Literal type for the value or None.
    """
    try:
        return literal(value)
    except LiteralTypingError:
        return


class Omitted(Opaque):
    """
    An omitted function argument with a default value.
    """

    def __init__(self, value):
        self._value = value
        # Use helper function to support both hashable and non-hashable
        # values. See discussion in gh #6957.
        self._value_key = get_hashable_key(value)
        super(Omitted, self).__init__("omitted(default=%r)" % (value,))

    @property
    def key(self):
        return type(self._value), self._value_key

    @property
    def value(self):
        return self._value


class VarArg(Type):
    """
    Special type representing a variable number of arguments at the
    end of a function's signature.  Only used for signature matching,
    not for actual values.
    """

    def __init__(self, dtype):
        self.dtype = dtype
        super(VarArg, self).__init__("*%s" % dtype)

    @property
    def key(self):
        return self.dtype


class Module(Dummy):
    def __init__(self, pymod):
        self.pymod = pymod
        super(Module, self).__init__("Module(%s)" % pymod)

    @property
    def key(self):
        return self.pymod


class MemInfoPointer(Type):
    """
    Pointer to a Numba "meminfo" (i.e. the information for a managed
    piece of memory).
    """
    mutable = True

    def __init__(self, dtype):
        self.dtype = dtype
        name = "memory-managed *%s" % dtype
        super(MemInfoPointer, self).__init__(name)

    @property
    def key(self):
        return self.dtype


class CPointer(Type):
    """
    Type class for pointers to other types.

    Attributes
    ----------
        dtype : The pointee type
        addrspace : int
            The address space pointee belongs to.
    """
    mutable = True

    def __init__(self, dtype, addrspace=None):
        self.dtype = dtype
        self.addrspace = addrspace
        if addrspace is not None:
            name = "%s_%s*" % (dtype, addrspace)
        else:
            name = "%s*" % dtype
        super(CPointer, self).__init__(name)

    @property
    def key(self):
        return self.dtype, self.addrspace


class EphemeralPointer(CPointer):
    """
    Type class for pointers which aren't guaranteed to last long - e.g.
    stack-allocated slots.  The data model serializes such pointers
    by copying the data pointed to.
    """


class EphemeralArray(Type):
    """
    Similar to EphemeralPointer, but pointing to an array of elements,
    rather than a single one.  The array size must be known at compile-time.
    """

    def __init__(self, dtype, count):
        self.dtype = dtype
        self.count = count
        name = "*%s[%d]" % (dtype, count)
        super(EphemeralArray, self).__init__(name)

    @property
    def key(self):
        return self.dtype, self.count


class Object(Type):
    # XXX unused?
    mutable = True

    def __init__(self, clsobj):
        self.cls = clsobj
        name = "Object(%s)" % clsobj.__name__
        super(Object, self).__init__(name)

    @property
    def key(self):
        return self.cls


class Optional(Type):
    """
    Type class for optional types, i.e. union { some type, None }
    """

    def __init__(self, typ):
        assert not isinstance(typ, (Optional, NoneType))
        typ = unliteral(typ)
        self.type = typ
        name = "OptionalType(%s)" % self.type
        super(Optional, self).__init__(name)

    @property
    def key(self):
        return self.type

    def can_convert_to(self, typingctx, other):
        if isinstance(other, Optional):
            return typingctx.can_convert(self.type, other.type)
        else:
            conv = typingctx.can_convert(self.type, other)
            if conv is not None:
                return max(conv, Conversion.safe)

    def can_convert_from(self, typingctx, other):
        if isinstance(other, NoneType):
            return Conversion.promote
        elif isinstance(other, Optional):
            return typingctx.can_convert(other.type, self.type)
        else:
            conv = typingctx.can_convert(other, self.type)
            if conv is not None:
                return max(conv, Conversion.promote)

    def unify(self, typingctx, other):
        if isinstance(other, Optional):
            unified = typingctx.unify_pairs(self.type, other.type)
        else:
            unified = typingctx.unify_pairs(self.type, other)

        if unified is not None:
            if isinstance(unified, Optional):
                return unified
            else:
                return Optional(unified)


class NoneType(Opaque):
    """
    The type for None.
    """

    def unify(self, typingctx, other):
        """
        Turn anything to a Optional type;
        """
        if isinstance(other, (Optional, NoneType)):
            return other
        return Optional(other)


class EllipsisType(Opaque):
    """
    The type for the Ellipsis singleton.
    """


class ExceptionClass(Callable, Phantom):
    """
    The type of exception classes (not instances).
    """

    def __init__(self, exc_class):
        assert issubclass(exc_class, BaseException)
        name = "%s" % (exc_class.__name__)
        self.exc_class = exc_class
        super(ExceptionClass, self).__init__(name)

    def get_call_type(self, context, args, kws):
        return self.get_call_signatures()[0][0]

    def get_call_signatures(self):
        from numba.core import typing
        return_type = ExceptionInstance(self.exc_class)
        return [typing.signature(return_type)], False

    def get_impl_key(self, sig):
        return type(self)

    @property
    def key(self):
        return self.exc_class


class ExceptionInstance(Phantom):
    """
    The type of exception instances.  *exc_class* should be the
    exception class.
    """

    def __init__(self, exc_class):
        assert issubclass(exc_class, BaseException)
        name = "%s(...)" % (exc_class.__name__,)
        self.exc_class = exc_class
        super(ExceptionInstance, self).__init__(name)

    @property
    def key(self):
        return self.exc_class


class SliceType(Type):

    def __init__(self, name, members):
        assert members in (2, 3)
        self.members = members
        self.has_step = members >= 3
        super(SliceType, self).__init__(name)

    @property
    def key(self):
        return self.members


class SliceLiteral(Literal, SliceType):
    def __init__(self, value):
        self._literal_init(value)
        name = 'Literal[slice]({})'.format(value)
        members = 2 if value.step is None else 3
        SliceType.__init__(self, name=name, members=members)

    @property
    def key(self):
        sl = self.literal_value
        return sl.start, sl.stop, sl.step


Literal.ctor_map[slice] = SliceLiteral


class ClassInstanceType(Type):
    """
    The type of a jitted class *instance*.  It will be the return-type
    of the constructor of the class.
    """
    mutable = True
    name_prefix = "instance"

    def __init__(self, class_type):
        self.class_type = class_type
        name = "{0}.{1}".format(self.name_prefix, self.class_type.name)
        super(ClassInstanceType, self).__init__(name)

    def get_data_type(self):
        return ClassDataType(self)

    def get_reference_type(self):
        return self

    @property
    def key(self):
        return self.class_type.key

    @property
    def classname(self):
        return self.class_type.class_name

    @property
    def jit_props(self):
        return self.class_type.jit_props

    @property
    def jit_static_methods(self):
        return self.class_type.jit_static_methods

    @property
    def jit_methods(self):
        return self.class_type.jit_methods

    @property
    def struct(self):
        return self.class_type.struct

    @property
    def methods(self):
        return self.class_type.methods

    @property
    def static_methods(self):
        return self.class_type.static_methods


class ClassType(Callable, Opaque):
    """
    The type of the jitted class (not instance).  When the type of a class
    is called, its constructor is invoked.
    """
    mutable = True
    name_prefix = "jitclass"
    instance_type_class = ClassInstanceType

    def __init__(self, class_def, ctor_template_cls, struct, jit_methods,
                 jit_props, jit_static_methods):
        self.class_name = class_def.__name__
        self.class_doc = class_def.__doc__
        self._ctor_template_class = ctor_template_cls
        self.jit_methods = jit_methods
        self.jit_props = jit_props
        self.jit_static_methods = jit_static_methods
        self.struct = struct
        fielddesc = ','.join("{0}:{1}".format(k, v) for k, v in struct.items())
        name = "{0}.{1}#{2:x}<{3}>".format(self.name_prefix, self.class_name,
                                           id(self), fielddesc)
        super(ClassType, self).__init__(name)

    def get_call_type(self, context, args, kws):
        return self.ctor_template(context).apply(args, kws)

    def get_call_signatures(self):
        return (), True

    def get_impl_key(self, sig):
        return type(self)

    @property
    def methods(self):
        return {k: v.py_func for k, v in self.jit_methods.items()}

    @property
    def static_methods(self):
        return {k: v.py_func for k, v in self.jit_static_methods.items()}

    @property
    def instance_type(self):
        return ClassInstanceType(self)

    @property
    def ctor_template(self):
        return self._specialize_template(self._ctor_template_class)

    def _specialize_template(self, basecls):
        return type(basecls.__name__, (basecls,), dict(key=self))


class DeferredType(Type):
    """
    Represents a type that will be defined later.  It must be defined
    before it is materialized (used in the compiler).  Once defined, it
    behaves exactly as the type it is defining.
    """

    def __init__(self):
        self._define = None
        name = "{0}#{1}".format(type(self).__name__, id(self))
        super(DeferredType, self).__init__(name)

    def get(self):
        if self._define is None:
            raise RuntimeError("deferred type not defined")
        return self._define

    def define(self, typ):
        if self._define is not None:
            raise TypeError("deferred type already defined")
        if not isinstance(typ, Type):
            raise TypeError("arg is not a Type; got: {0}".format(type(typ)))
        self._define = typ

    def unify(self, typingctx, other):
        return typingctx.unify_pairs(self.get(), other)


class ClassDataType(Type):
    """
    Internal only.
    Represents the data of the instance.  The representation of
    ClassInstanceType contains a pointer to a ClassDataType which represents
    a C structure that contains all the data fields of the class instance.
    """

    def __init__(self, classtyp):
        self.class_type = classtyp
        name = "data.{0}".format(self.class_type.name)
        super(ClassDataType, self).__init__(name)


class ContextManager(Callable, Phantom):
    """
    An overly-simple ContextManager type that cannot be materialized.
    """

    def __init__(self, cm):
        self.cm = cm
        super(ContextManager, self).__init__("ContextManager({})".format(cm))

    def get_call_signatures(self):
        if not self.cm.is_callable:
            msg = "contextmanager {} is not callable".format(self.cm)
            raise TypingError(msg)

        return (), False

    def get_call_type(self, context, args, kws):
        from numba.core import typing

        if not self.cm.is_callable:
            msg = "contextmanager {} is not callable".format(self.cm)
            raise TypingError(msg)

        posargs = list(args) + [v for k, v in sorted(kws.items())]
        return typing.signature(self, *posargs)

    def get_impl_key(self, sig):
        return type(self)


class UnicodeType(IterableType, Hashable):

    def __init__(self, name):
        super(UnicodeType, self).__init__(name)

    @property
    def iterator_type(self):
        return UnicodeIteratorType(self)


class UnicodeIteratorType(SimpleIteratorType):

    def __init__(self, dtype):
        name = "iter_unicode"
        self.data = dtype
        super(UnicodeIteratorType, self).__init__(name, dtype)
