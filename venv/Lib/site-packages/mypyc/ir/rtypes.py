"""Types used in the intermediate representation.

These are runtime types (RTypes), as opposed to mypy Type objects.
The latter are only used during type checking and not directly used at
runtime.  Runtime types are derived from mypy types, but there's no
simple one-to-one correspondence. (Here 'runtime' means 'runtime
checked'.)

The generated IR ensures some runtime type safety properties based on
RTypes. Compiled code can assume that the runtime value matches the
static RType of a value. If the RType of a register is 'builtins.str'
(str_rprimitive), for example, the generated IR will ensure that the
register will have a 'str' object.

RTypes are simpler and less expressive than mypy (or PEP 484)
types. For example, all mypy types of form 'list[T]' (for arbitrary T)
are erased to the single RType 'builtins.list' (list_rprimitive).

mypyc.irbuild.mapper.Mapper.type_to_rtype converts mypy Types to mypyc
RTypes.

NOTE: As a convention, we don't create subclasses of concrete RType
      subclasses (e.g. you shouldn't define a subclass of RTuple, which
      is a concrete class). We prefer a flat class hierarchy.

      If you want to introduce a variant of an existing class, you'd
      typically add an attribute (e.g. a flag) to an existing concrete
      class to enable the new behavior. In rare cases, adding a new
      abstract base class could also be an option. Adding a completely
      separate class and sharing some functionality using module-level
      helper functions may also be reasonable.

      This makes it possible to use isinstance(x, <concrete RType
      subclass>) checks without worrying about potential subclasses
      and avoids most trouble caused by implementation inheritance.
"""

from __future__ import annotations

from abc import abstractmethod
from typing import TYPE_CHECKING, ClassVar, Final, Generic, TypeGuard, TypeVar, Union, final

from mypyc.common import HAVE_IMMORTAL, IS_32_BIT_PLATFORM, PLATFORM_SIZE, JsonDict, short_name
from mypyc.ir.deps import LIBRT_RANDOM, LIBRT_STRINGS, LIBRT_VECS, Dependency
from mypyc.namegen import NameGenerator

if TYPE_CHECKING:
    from mypyc.ir.class_ir import ClassIR
    from mypyc.ir.ops import DeserMaps

T = TypeVar("T")


class RType:
    """Abstract base class for runtime types (erased, only concrete; no generics)."""

    name: str
    # If True, the type has a special unboxed representation. If False, the
    # type is represented as PyObject *. Even if True, the representation
    # may contain pointers.
    is_unboxed = False
    # This is the C undefined value for this type. It's used for initialization
    # if there's no value yet, and for function return value on error/exception.
    #
    # TODO: This shouldn't be specific to C or a string
    c_undefined: str
    # If unboxed: does the unboxed version use reference counting?
    is_refcounted = True
    # C type; use Emitter.ctype() to access
    _ctype: str
    # If True, error/undefined value overlaps with a valid value. To
    # detect an exception, PyErr_Occurred() must be used in addition
    # to checking for error value as the return value of a function.
    #
    # For example, no i64 value can be reserved for error value, so we
    # pick an arbitrary value (-113) to signal error, but this is
    # also a valid non-error value. The chosen value is rare as a
    # normal, non-error value, so most of the time we can avoid calling
    # PyErr_Occurred() when checking for errors raised by called
    # functions.
    #
    # This also means that if an attribute with this type might be
    # undefined, we can't just rely on the error value to signal this.
    # Instead, we add a bitfield to keep track whether attributes with
    # "error overlap" have a value. If there is no value, AttributeError
    # is raised on attribute read. Parameters with default values also
    # use the bitfield trick to indicate whether the caller passed a
    # value. (If we can determine that an attribute is "always defined",
    # we never raise an AttributeError and don't need the bitfield
    # entry.)
    error_overlap = False
    dependencies: tuple[Dependency, ...] | None = None

    @abstractmethod
    def accept(self, visitor: RTypeVisitor[T]) -> T:
        raise NotImplementedError()

    def short_name(self) -> str:
        return short_name(self.name)

    @property
    @abstractmethod
    def may_be_immortal(self) -> bool:
        raise NotImplementedError

    def __str__(self) -> str:
        return short_name(self.name)

    def __repr__(self) -> str:
        return "<%s>" % self.__class__.__name__

    def serialize(self) -> JsonDict | str:
        raise NotImplementedError(f"Cannot serialize {self.__class__.__name__} instance")


def deserialize_type(data: JsonDict | str, ctx: DeserMaps) -> RType:
    """Deserialize a JSON-serialized RType.

    Arguments:
        data: The decoded JSON of the serialized type
        ctx: The deserialization maps to use
    """
    # Since there are so few types, we just case on them directly.  If
    # more get added we should switch to a system like mypy.types
    # uses.
    if isinstance(data, str):
        if data in ctx.classes:
            return RInstance(ctx.classes[data])
        elif data in RPrimitive.primitive_map:
            return RPrimitive.primitive_map[data]
        elif data == "void":
            return RVoid()
        else:
            assert False, f"Can't find class {data}"
    elif data[".class"] == "RTuple":
        return RTuple.deserialize(data, ctx)
    elif data[".class"] == "RUnion":
        return RUnion.deserialize(data, ctx)
    elif data[".class"] == "RVec":
        return RVec.deserialize(data, ctx)
    raise NotImplementedError("unexpected .class {}".format(data[".class"]))


class RTypeVisitor(Generic[T]):
    """Generic visitor over RTypes (uses the visitor design pattern)."""

    @abstractmethod
    def visit_rprimitive(self, typ: RPrimitive, /) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_rinstance(self, typ: RInstance, /) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_rvec(self, typ: RVec, /) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_runion(self, typ: RUnion, /) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_rtuple(self, typ: RTuple, /) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_rstruct(self, typ: RStruct, /) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_rarray(self, typ: RArray, /) -> T:
        raise NotImplementedError

    @abstractmethod
    def visit_rvoid(self, typ: RVoid, /) -> T:
        raise NotImplementedError


# These types are implemented as C structs.
RStructLike = Union["RStruct", "RVec"]


@final
class RVoid(RType):
    """The void type (no value).

    This is a singleton -- use void_rtype (below) to refer to this instead of
    constructing a new instance.
    """

    is_unboxed = False
    name = "void"
    ctype = "void"

    def accept(self, visitor: RTypeVisitor[T]) -> T:
        return visitor.visit_rvoid(self)

    @property
    def may_be_immortal(self) -> bool:
        return False

    def serialize(self) -> str:
        return "void"

    def __eq__(self, other: object) -> TypeGuard[RVoid]:
        return isinstance(other, RVoid)

    def __hash__(self) -> int:
        return hash(RVoid)


# Singleton instance of RVoid
void_rtype: Final = RVoid()


@final
class RPrimitive(RType):
    """Primitive type such as 'object' or 'int'.

    These often have custom ops associated with them. The 'object'
    primitive type can be used to hold arbitrary Python objects.

    Different primitive types have different representations, and
    primitives may be unboxed or boxed. Primitive types don't need to
    directly correspond to Python types, but most do.

    NOTE: All supported primitive types are defined below
    (e.g. object_rprimitive).
    """

    # Map from primitive names to primitive types and is used by deserialization
    primitive_map: ClassVar[dict[str, RPrimitive]] = {}

    def __init__(
        self,
        name: str,
        *,
        is_unboxed: bool,
        is_refcounted: bool,
        is_native_int: bool = False,
        is_signed: bool = False,
        ctype: str = "PyObject *",
        size: int = PLATFORM_SIZE,
        error_overlap: bool = False,
        may_be_immortal: bool = True,
        dependencies: tuple[Dependency, ...] | None = None,
    ) -> None:
        RPrimitive.primitive_map[name] = self

        self.name = name
        self.is_unboxed = is_unboxed
        self.is_refcounted = is_refcounted
        self.is_native_int = is_native_int
        self.is_signed = is_signed
        self._ctype = ctype
        self.size = size
        self.error_overlap = error_overlap
        self._may_be_immortal = may_be_immortal and HAVE_IMMORTAL
        self.dependencies = dependencies
        if ctype == "CPyTagged":
            self.c_undefined = "CPY_INT_TAG"
        elif ctype in ("int16_t", "int32_t", "int64_t"):
            # This is basically an arbitrary value that is pretty
            # unlikely to overlap with a real value.
            self.c_undefined = "-113"
        elif ctype == "CPyPtr":
            # TODO: Invent an overlapping error value?
            self.c_undefined = "0"
        elif ctype.endswith("*"):
            # Boxed and pointer types use the null pointer as the error value.
            self.c_undefined = "NULL"
        elif ctype == "char":
            self.c_undefined = "2"
        elif ctype == "double":
            self.c_undefined = "-113.0"
        elif ctype in ("uint8_t", "uint16_t", "uint32_t", "uint64_t"):
            self.c_undefined = "239"  # An arbitrary number
        else:
            assert False, "Unrecognized ctype: %r" % ctype

    def accept(self, visitor: RTypeVisitor[T]) -> T:
        return visitor.visit_rprimitive(self)

    @property
    def may_be_immortal(self) -> bool:
        return self._may_be_immortal

    def serialize(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return "<RPrimitive %s>" % self.name

    def __eq__(self, other: object) -> TypeGuard[RPrimitive]:
        return isinstance(other, RPrimitive) and other.name == self.name

    def __hash__(self) -> int:
        return hash(self.name)


# NOTE: All the supported instances of RPrimitive are defined
# below. Use these instead of creating new instances.

# Used to represent arbitrary objects and dynamically typed (Any)
# values. There are various ops that let you perform generic, runtime
# checked operations on these (that match Python semantics). See the
# ops in mypyc.primitives.misc_ops, including py_getattr_op,
# py_call_op, and many others.
#
# If there is no more specific RType available for some value, we fall
# back to using this type.
#
# NOTE: Even though this is very flexible, this type should be used as
# little as possible, as generic ops are typically slow. Other types,
# including other primitive types and RInstance, are usually much
# faster.
object_rprimitive: Final = RPrimitive("builtins.object", is_unboxed=False, is_refcounted=True)

# represents a low level pointer of an object
object_pointer_rprimitive: Final = RPrimitive(
    "object_ptr", is_unboxed=False, is_refcounted=False, ctype="PyObject **"
)

# Similar to object_primitive, but does not use automatic reference
# counting. Useful for temporaries.
object_non_refcounted_rprimitive: Final = RPrimitive(
    "builtins.object_nrc", is_unboxed=False, is_refcounted=False
)

# Arbitrary-precision integer (corresponds to Python 'int'). Small
# enough values are stored unboxed, while large integers are
# represented as a tagged pointer to a Python 'int' PyObject. The
# lowest bit is used as the tag to decide whether it is a signed
# unboxed value (shifted left by one) or a PyObject * pointing to an
# 'int' object. Pointers have the least significant bit set.
#
# The undefined/error value is the null pointer (1 -- only the least
# significant bit is set)).
#
# This cannot represent a subclass of int. An instance of a subclass
# of int is coerced to the corresponding 'int' value.
int_rprimitive: Final = RPrimitive(
    "builtins.int", is_unboxed=True, is_refcounted=True, ctype="CPyTagged"
)

# An unboxed integer. The representation is the same as for unboxed
# int_rprimitive (shifted left by one). These can be used when an
# integer is known to be small enough to fit size_t (CPyTagged).
short_int_rprimitive: Final = RPrimitive(
    "short_int", is_unboxed=True, is_refcounted=False, ctype="CPyTagged"
)

# Low level integer types (correspond to C integer types)

int16_rprimitive: Final = RPrimitive(
    "i16",
    is_unboxed=True,
    is_refcounted=False,
    is_native_int=True,
    is_signed=True,
    ctype="int16_t",
    size=2,
    error_overlap=True,
)
int32_rprimitive: Final = RPrimitive(
    "i32",
    is_unboxed=True,
    is_refcounted=False,
    is_native_int=True,
    is_signed=True,
    ctype="int32_t",
    size=4,
    error_overlap=True,
)
int64_rprimitive: Final = RPrimitive(
    "i64",
    is_unboxed=True,
    is_refcounted=False,
    is_native_int=True,
    is_signed=True,
    ctype="int64_t",
    size=8,
    error_overlap=True,
)
uint8_rprimitive: Final = RPrimitive(
    "u8",
    is_unboxed=True,
    is_refcounted=False,
    is_native_int=True,
    is_signed=False,
    ctype="uint8_t",
    size=1,
    error_overlap=True,
)

# The following unsigned native int types (u16, u32, u64) are not
# exposed to the user. They are for internal use within mypyc only.

u16_rprimitive: Final = RPrimitive(
    "u16",
    is_unboxed=True,
    is_refcounted=False,
    is_native_int=True,
    is_signed=False,
    ctype="uint16_t",
    size=2,
    error_overlap=True,
)
uint32_rprimitive: Final = RPrimitive(
    "u32",
    is_unboxed=True,
    is_refcounted=False,
    is_native_int=True,
    is_signed=False,
    ctype="uint32_t",
    size=4,
    error_overlap=True,
)
uint64_rprimitive: Final = RPrimitive(
    "u64",
    is_unboxed=True,
    is_refcounted=False,
    is_native_int=True,
    is_signed=False,
    ctype="uint64_t",
    size=8,
    error_overlap=True,
)

# The C 'int' type
c_int_rprimitive = int32_rprimitive

if IS_32_BIT_PLATFORM:
    c_size_t_rprimitive = uint32_rprimitive
    c_pyssize_t_rprimitive = RPrimitive(
        "native_int",
        is_unboxed=True,
        is_refcounted=False,
        is_native_int=True,
        is_signed=True,
        ctype="int32_t",
        size=4,
    )
else:
    c_size_t_rprimitive = uint64_rprimitive
    c_pyssize_t_rprimitive = RPrimitive(
        "native_int",
        is_unboxed=True,
        is_refcounted=False,
        is_native_int=True,
        is_signed=True,
        ctype="int64_t",
        size=8,
    )

# Untyped pointer, represented as integer in the C backend
pointer_rprimitive: Final = RPrimitive("ptr", is_unboxed=True, is_refcounted=False, ctype="CPyPtr")

# Untyped pointer, represented as void * in the C backend
c_pointer_rprimitive: Final = RPrimitive(
    "c_ptr", is_unboxed=False, is_refcounted=False, ctype="void *"
)

cstring_rprimitive: Final = RPrimitive(
    "cstring", is_unboxed=True, is_refcounted=False, ctype="const char *"
)

# The type corresponding to mypyc.common.BITMAP_TYPE
bitmap_rprimitive: Final = uint32_rprimitive

# Floats are represent as 'float' PyObject * values. (In the future
# we'll likely switch to a more efficient, unboxed representation.)
float_rprimitive: Final = RPrimitive(
    "builtins.float",
    is_unboxed=True,
    is_refcounted=False,
    ctype="double",
    size=8,
    error_overlap=True,
)

# An unboxed Python bool value. This actually has three possible values
# (0 -> False, 1 -> True, 2 -> error). If you only need True/False, use
# bit_rprimitive instead.
bool_rprimitive: Final = RPrimitive(
    "builtins.bool", is_unboxed=True, is_refcounted=False, ctype="char", size=1
)

# A low-level boolean value with two possible values: 0 and 1. Any
# other value results in undefined behavior. Undefined or error values
# are not supported.
bit_rprimitive: Final = RPrimitive(
    "bit", is_unboxed=True, is_refcounted=False, ctype="char", size=1
)

# The 'None' value. The possible values are 0 -> None and 2 -> error.
none_rprimitive: Final = RPrimitive(
    "builtins.None", is_unboxed=True, is_refcounted=False, ctype="char", size=1
)

# Python list object (or an instance of a subclass of list). These could be
# immortal, but since this is expected to be very rare, and the immortality checks
# can be pretty expensive for lists, we treat lists as non-immortal.
list_rprimitive: Final = RPrimitive(
    "builtins.list", is_unboxed=False, is_refcounted=True, may_be_immortal=False
)

# Python dict object (or an instance of a subclass of dict).
dict_rprimitive: Final = RPrimitive("builtins.dict", is_unboxed=False, is_refcounted=True)

# Python set object (or an instance of a subclass of set).
set_rprimitive: Final = RPrimitive("builtins.set", is_unboxed=False, is_refcounted=True)

# Python frozenset object (or an instance of a subclass of frozenset).
frozenset_rprimitive: Final = RPrimitive(
    "builtins.frozenset", is_unboxed=False, is_refcounted=True
)

# Python str object. At the C layer, str is referred to as unicode
# (PyUnicode).
str_rprimitive: Final = RPrimitive("builtins.str", is_unboxed=False, is_refcounted=True)

# Python bytes object.
bytes_rprimitive: Final = RPrimitive("builtins.bytes", is_unboxed=False, is_refcounted=True)

# Python bytearray object.
bytearray_rprimitive: Final = RPrimitive(
    "builtins.bytearray", is_unboxed=False, is_refcounted=True
)

# Tuple of an arbitrary length (corresponds to Tuple[t, ...], with
# explicit '...').
tuple_rprimitive: Final = RPrimitive("builtins.tuple", is_unboxed=False, is_refcounted=True)

# Python range object.
range_rprimitive: Final = RPrimitive("builtins.range", is_unboxed=False, is_refcounted=True)

KNOWN_NATIVE_TYPES: Final = {
    name: RPrimitive(name, is_unboxed=False, is_refcounted=True, dependencies=deps)
    for name, deps in [
        ("librt.internal.WriteBuffer", ()),
        ("librt.internal.ReadBuffer", ()),
        ("librt.strings.BytesWriter", (LIBRT_STRINGS,)),
        ("librt.strings.StringWriter", (LIBRT_STRINGS,)),
    ]
} | {
    "librt.random.Random": RPrimitive(
        "librt.random.Random", is_unboxed=False, is_refcounted=True, dependencies=(LIBRT_RANDOM,)
    )
}

bytes_writer_rprimitive: Final = KNOWN_NATIVE_TYPES["librt.strings.BytesWriter"]
string_writer_rprimitive: Final = KNOWN_NATIVE_TYPES["librt.strings.StringWriter"]
random_rprimitive: Final = KNOWN_NATIVE_TYPES["librt.random.Random"]


def is_native_rprimitive(rtype: RType) -> bool:
    return isinstance(rtype, RPrimitive) and rtype.name in KNOWN_NATIVE_TYPES


def is_tagged(rtype: RType) -> TypeGuard[RPrimitive]:
    return rtype is int_rprimitive or rtype is short_int_rprimitive


def is_any_int(rtype: RType) -> bool:
    return is_tagged(rtype) or is_int32_rprimitive(rtype) or is_int64_rprimitive(rtype)


def is_int_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return rtype is int_rprimitive


def is_short_int_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return rtype is short_int_rprimitive


def is_int16_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return rtype is int16_rprimitive


def is_int32_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return rtype is int32_rprimitive or (
        rtype is c_pyssize_t_rprimitive and rtype._ctype == "int32_t"
    )


def is_int64_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return rtype is int64_rprimitive or (
        rtype is c_pyssize_t_rprimitive and rtype._ctype == "int64_t"
    )


def is_fixed_width_rtype(rtype: RType) -> TypeGuard[RPrimitive]:
    return (
        is_int64_rprimitive(rtype)
        or is_int32_rprimitive(rtype)
        or is_int16_rprimitive(rtype)
        or is_uint8_rprimitive(rtype)
    )


def is_uint8_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return rtype is uint8_rprimitive


def is_uint32_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return rtype is uint32_rprimitive


def is_uint64_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return rtype is uint64_rprimitive


def is_c_py_ssize_t_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return rtype is c_pyssize_t_rprimitive


def is_pointer_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return rtype is pointer_rprimitive


def is_float_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return isinstance(rtype, RPrimitive) and rtype.name == "builtins.float"


def is_bool_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return isinstance(rtype, RPrimitive) and rtype.name == "builtins.bool"


def is_bit_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return isinstance(rtype, RPrimitive) and rtype.name == "bit"


def is_bool_or_bit_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return is_bool_rprimitive(rtype) or is_bit_rprimitive(rtype)


def is_object_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return isinstance(rtype, RPrimitive) and rtype.name == "builtins.object"


def is_none_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return isinstance(rtype, RPrimitive) and rtype.name == "builtins.None"


def is_list_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return isinstance(rtype, RPrimitive) and rtype.name == "builtins.list"


def is_dict_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return isinstance(rtype, RPrimitive) and rtype.name == "builtins.dict"


def is_set_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return isinstance(rtype, RPrimitive) and rtype.name == "builtins.set"


def is_frozenset_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return isinstance(rtype, RPrimitive) and rtype.name == "builtins.frozenset"


def is_str_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return isinstance(rtype, RPrimitive) and rtype.name == "builtins.str"


def is_bytes_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return isinstance(rtype, RPrimitive) and rtype.name == "builtins.bytes"


def is_bytearray_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return isinstance(rtype, RPrimitive) and rtype.name == "builtins.bytearray"


def is_tuple_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return isinstance(rtype, RPrimitive) and rtype.name == "builtins.tuple"


def is_range_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return isinstance(rtype, RPrimitive) and rtype.name == "builtins.range"


def is_sequence_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return isinstance(rtype, RPrimitive) and (
        is_list_rprimitive(rtype)
        or is_tuple_rprimitive(rtype)
        or is_str_rprimitive(rtype)
        or is_bytes_rprimitive(rtype)
        or is_bytearray_rprimitive(rtype)
    )


def is_immutable_rprimitive(rtype: RType) -> TypeGuard[RPrimitive]:
    return (
        is_str_rprimitive(rtype)
        or is_bytes_rprimitive(rtype)
        or is_tuple_rprimitive(rtype)
        or is_frozenset_rprimitive(rtype)
    )


class TupleNameVisitor(RTypeVisitor[str]):
    """Produce a tuple name based on the concrete representations of types."""

    def visit_rinstance(self, t: RInstance) -> str:
        return "O"

    def visit_rvec(self, t: RVec) -> str:
        if isinstance(t.item_type, RVec):
            # All nested vecs share a representation
            return "Vv"
        return "V" + t.item_type.accept(self)

    def visit_runion(self, t: RUnion) -> str:
        return "O"

    def visit_rprimitive(self, t: RPrimitive) -> str:
        if t._ctype == "CPyTagged":
            return "I"
        elif t._ctype == "char":
            return "C"
        elif t._ctype == "int64_t":
            return "8"  # "8 byte integer"
        elif t._ctype == "int32_t":
            return "4"  # "4 byte integer"
        elif t._ctype == "int16_t":
            return "2"  # "2 byte integer"
        elif t._ctype == "uint8_t":
            return "U1"  # "1 byte unsigned integer"
        elif t._ctype == "double":
            return "F"
        assert not t.is_unboxed, f"{t} unexpected unboxed type"
        return "O"

    def visit_rtuple(self, t: RTuple) -> str:
        parts = [elem.accept(self) for elem in t.types]
        return "T{}{}".format(len(parts), "".join(parts))

    def visit_rstruct(self, t: RStruct) -> str:
        if t.name == "VecNestedBufItem":
            return "Vi"
        assert False, "RStruct not supported in tuple"

    def visit_rarray(self, t: RArray) -> str:
        assert False, "RArray not supported in tuple"

    def visit_rvoid(self, t: RVoid) -> str:
        assert False, "rvoid in tuple?"


@final
class RTuple(RType):
    """Fixed-length unboxed tuple (represented as a C struct).

    These are used to represent mypy TupleType values (fixed-length
    Python tuples). Since this is unboxed, the identity of a tuple
    object is not preserved within compiled code. If the identity of a
    tuple is important, or there is a need to have multiple references
    to a single tuple object, a variable-length tuple should be used
    (tuple_rprimitive or Tuple[T, ...]  with explicit '...'), as they
    are boxed.

    These aren't immutable. However, user code won't be able to mutate
    individual tuple items.
    """

    is_unboxed = True

    def __init__(self, types: list[RType]) -> None:
        self.name = "tuple"
        self.types = tuple(types)
        self.is_refcounted = any(t.is_refcounted for t in self.types)
        # Generate a unique id which is used in naming corresponding C identifiers.
        # This is necessary since C does not have anonymous structural type equivalence
        # in the same way python can just assign a Tuple[int, bool] to a Tuple[int, bool].
        self.unique_id = self.accept(TupleNameVisitor())
        # Nominally the max c length is 31 chars, but I'm not honestly worried about this.
        self.struct_name = f"tuple_{self.unique_id}"
        self._ctype = f"{self.struct_name}"
        self.error_overlap = all(t.error_overlap for t in self.types) and bool(self.types)

    def accept(self, visitor: RTypeVisitor[T]) -> T:
        return visitor.visit_rtuple(self)

    @property
    def may_be_immortal(self) -> bool:
        return False

    def __str__(self) -> str:
        return "tuple[%s]" % ", ".join(str(typ) for typ in self.types)

    def __repr__(self) -> str:
        return "<RTuple %s>" % ", ".join(repr(typ) for typ in self.types)

    def __eq__(self, other: object) -> TypeGuard[RTuple]:
        return isinstance(other, RTuple) and self.types == other.types

    def __hash__(self) -> int:
        return hash((self.name, self.types))

    def serialize(self) -> JsonDict:
        types = [x.serialize() for x in self.types]
        return {".class": "RTuple", "types": types}

    @classmethod
    def deserialize(cls, data: JsonDict, ctx: DeserMaps) -> RTuple:
        types = [deserialize_type(t, ctx) for t in data["types"]]
        return RTuple(types)


# Exception tuple: (exception class, exception instance, traceback object)
exc_rtuple = RTuple([object_rprimitive, object_rprimitive, object_rprimitive])

# Dictionary iterator tuple: (should continue, internal offset, key, value)
# See mypyc.irbuild.for_helpers.ForDictionaryCommon for more details.
dict_next_rtuple_pair = RTuple(
    [bool_rprimitive, short_int_rprimitive, object_rprimitive, object_rprimitive]
)
# Same as above but just for key or value.
dict_next_rtuple_single = RTuple([bool_rprimitive, short_int_rprimitive, object_rprimitive])


def compute_rtype_alignment(typ: RType) -> int:
    """Compute alignment of a given type based on platform alignment rule"""
    platform_alignment = PLATFORM_SIZE
    if isinstance(typ, RPrimitive):
        return typ.size
    elif isinstance(typ, RInstance):
        return platform_alignment
    elif isinstance(typ, RUnion):
        return platform_alignment
    elif isinstance(typ, RArray):
        return compute_rtype_alignment(typ.item_type)
    else:
        if isinstance(typ, RTuple):
            items = list(typ.types)
        elif isinstance(typ, RStruct):
            items = typ.types
        else:
            assert False, "invalid rtype for computing alignment"
        max_alignment = max(compute_rtype_alignment(item) for item in items)
        return max_alignment


def compute_rtype_size(typ: RType) -> int:
    """Compute unaligned size of rtype"""
    if isinstance(typ, RPrimitive):
        return typ.size
    elif isinstance(typ, RTuple):
        return compute_aligned_offsets_and_size(list(typ.types))[1]
    elif isinstance(typ, RUnion):
        return PLATFORM_SIZE
    elif isinstance(typ, RStruct):
        return compute_aligned_offsets_and_size(typ.types)[1]
    elif isinstance(typ, RInstance):
        return PLATFORM_SIZE
    elif isinstance(typ, RArray):
        alignment = compute_rtype_alignment(typ)
        aligned_size = (compute_rtype_size(typ.item_type) + (alignment - 1)) & ~(alignment - 1)
        return aligned_size * typ.length
    else:
        assert False, "invalid rtype for computing size"


def compute_aligned_offsets_and_size(types: list[RType]) -> tuple[list[int], int]:
    """Compute offsets and total size of a list of types after alignment

    Note that the types argument are types of values that are stored
    sequentially with platform default alignment.
    """
    unaligned_sizes = [compute_rtype_size(typ) for typ in types]
    alignments = [compute_rtype_alignment(typ) for typ in types]

    current_offset = 0
    offsets = []
    final_size = 0
    for i in range(len(unaligned_sizes)):
        offsets.append(current_offset)
        if i + 1 < len(unaligned_sizes):
            cur_size = unaligned_sizes[i]
            current_offset += cur_size
            next_alignment = alignments[i + 1]
            # compute aligned offset,
            # check https://en.wikipedia.org/wiki/Data_structure_alignment for more information
            current_offset = (current_offset + (next_alignment - 1)) & -next_alignment
        else:
            struct_alignment = max(alignments)
            final_size = current_offset + unaligned_sizes[i]
            final_size = (final_size + (struct_alignment - 1)) & -struct_alignment
    return offsets, final_size


@final
class RStruct(RType):
    """C struct type"""

    def __init__(self, name: str, names: list[str], types: list[RType]) -> None:
        self.name = name
        self.names = names
        self.types = types
        self.is_refcounted = any(t.is_refcounted for t in self.types)

        # generate dummy names
        if len(self.names) < len(self.types):
            for i in range(len(self.types) - len(self.names)):
                self.names.append("_item" + str(i))
        self.offsets, self.size = compute_aligned_offsets_and_size(types)
        self._ctype = name

    def field_type(self, name: str) -> RType:
        for n, t in zip(self.names, self.types):
            if n == name:
                return t
        assert False, f"{self.name} has no field '{name}'"

    def accept(self, visitor: RTypeVisitor[T]) -> T:
        return visitor.visit_rstruct(self)

    @property
    def may_be_immortal(self) -> bool:
        return False

    def __str__(self) -> str:
        # if not tuple(unnamed structs)
        return "{}{{{}}}".format(
            self.name,
            ", ".join(name + ":" + str(typ) for name, typ in zip(self.names, self.types)),
        )

    def __repr__(self) -> str:
        return "<RStruct {}{{{}}}>".format(
            self.name,
            ", ".join(name + ":" + repr(typ) for name, typ in zip(self.names, self.types)),
        )

    def __eq__(self, other: object) -> TypeGuard[RStruct]:
        return (
            isinstance(other, RStruct)
            and self.name == other.name
            and self.names == other.names
            and self.types == other.types
        )

    def __hash__(self) -> int:
        return hash((self.name, tuple(self.names), tuple(self.types)))

    def serialize(self) -> JsonDict:
        assert False

    @classmethod
    def deserialize(cls, data: JsonDict, ctx: DeserMaps) -> RStruct:
        assert False


@final
class RInstance(RType):
    """Instance of user-defined class (compiled to C extension class).

    The runtime representation is 'PyObject *', and these are always
    boxed and thus reference-counted.

    These support fast method calls and fast attribute access using
    vtables, and they usually use a dict-free, struct-based
    representation of attributes. Method calls and attribute access
    can skip the vtable if we know that there is no overriding.

    These are also sometimes called 'native' types, since these have
    the most efficient representation and ops (along with certain
    RPrimitive types and RTuple).
    """

    is_unboxed = False

    def __init__(self, class_ir: ClassIR) -> None:
        # name is used for formatting the name in messages and debug output
        # so we want the fullname for precision.
        self.name = class_ir.fullname
        self.class_ir = class_ir
        self._ctype = "PyObject *"

    def accept(self, visitor: RTypeVisitor[T]) -> T:
        return visitor.visit_rinstance(self)

    @property
    def may_be_immortal(self) -> bool:
        return False

    def struct_name(self, names: NameGenerator) -> str:
        return self.class_ir.struct_name(names)

    def getter_index(self, name: str) -> int:
        return self.class_ir.vtable_entry(name)

    def setter_index(self, name: str) -> int:
        return self.getter_index(name) + 1

    def method_index(self, name: str) -> int:
        return self.class_ir.vtable_entry(name)

    def attr_type(self, name: str) -> RType:
        return self.class_ir.attr_type(name)

    def __repr__(self) -> str:
        return "<RInstance %s>" % self.name

    def __eq__(self, other: object) -> TypeGuard[RInstance]:
        return isinstance(other, RInstance) and other.name == self.name

    def __hash__(self) -> int:
        return hash(self.name)

    def serialize(self) -> str:
        return self.name


@final
class RVec(RType):
    """librt.vecs.vec[T]"""

    is_unboxed = True

    def __init__(self, item_type: RType) -> None:
        self.name = "vec[%s]" % item_type
        self.item_type = item_type
        self.names = ["len", "items"]
        self.dependencies = (LIBRT_VECS,)
        if isinstance(item_type, RUnion):
            non_opt = optional_value_type(item_type)
        else:
            non_opt = item_type
        if item_type in vec_buf_types:
            self._ctype = vec_c_types[item_type]
            self.buf_type = vec_buf_types[item_type]
            self.types = [c_pyssize_t_rprimitive, pointer_rprimitive]
        elif isinstance(non_opt, RVec):
            self._ctype = "VecNested"
            self.types = [c_pyssize_t_rprimitive, pointer_rprimitive]
            self.buf_type = VecNestedBufObject
        else:
            self._ctype = "VecT"
            self.types = [c_pyssize_t_rprimitive, pointer_rprimitive]
            self.buf_type = VecTBufObject

    @property
    def may_be_immortal(self) -> bool:
        return False

    def unwrap_item_type(self) -> RPrimitive | RInstance:
        """Return the non-optional value (non-vec) item type in a potentially nested vec."""
        item_type = self.item_type
        while True:
            if isinstance(item_type, RUnion):
                value_type = optional_value_type(item_type)
                assert value_type is not None
                item_type = value_type
            elif isinstance(item_type, RVec):
                item_type = item_type.item_type
            elif isinstance(item_type, (RPrimitive, RInstance)):
                return item_type
            else:
                assert False, f"unexpected item type: {self.item_type}"

    def is_optional(self) -> bool:
        item_type = self.item_type
        if isinstance(item_type, RUnion):
            return True
        elif isinstance(item_type, RVec):
            return item_type.is_optional()
        return False

    def depth(self) -> int:
        item_type = self.item_type
        if isinstance(item_type, RUnion):
            value_type = optional_value_type(item_type)
            assert value_type is not None
            item_type = value_type
        if isinstance(item_type, RVec):
            return 1 + item_type.depth()
        return 0

    def field_type(self, name: str) -> RType:
        if name == "len":
            return c_pyssize_t_rprimitive
        elif name == "items":
            return pointer_rprimitive
        assert False, f"RVec has no field '{name}'"

    def accept(self, visitor: RTypeVisitor[T]) -> T:
        return visitor.visit_rvec(self)

    def __str__(self) -> str:
        if self.is_optional() and self.depth() == 0:
            type_str = f"{self.unwrap_item_type()} | None"
        else:
            type_str = str(self.item_type)
        return f"vec[{type_str}]"

    def __repr__(self) -> str:
        return "<RVec[%s]>" % self.item_type

    def __eq__(self, other: object) -> bool:
        return isinstance(other, RVec) and other.item_type == self.item_type

    def __hash__(self) -> int:
        return hash(self.item_type) ^ 1

    def serialize(self) -> JsonDict:
        return {".class": "RVec", "item_type": self.item_type.serialize()}

    @classmethod
    def deserialize(cls, data: JsonDict, ctx: DeserMaps) -> RVec:
        return RVec(deserialize_type(data["item_type"], ctx))


@final
class RUnion(RType):
    """union[x, ..., y]"""

    is_unboxed = False

    def __init__(self, items: list[RType]) -> None:
        self.name = "union"
        self.items = items
        self.items_set = frozenset(items)
        self._ctype = "PyObject *"

    @staticmethod
    def make_simplified_union(items: list[RType]) -> RType:
        """Return a normalized union that covers the given items.

        Flatten nested unions and remove duplicate items.

        Overlapping items are *not* simplified. For example,
        [object, str] will not be simplified.
        """
        items = flatten_nested_unions(items)
        assert items

        unique_items = dict.fromkeys(items)
        if len(unique_items) > 1:
            return RUnion(list(unique_items))
        else:
            return next(iter(unique_items))

    def accept(self, visitor: RTypeVisitor[T]) -> T:
        return visitor.visit_runion(self)

    @property
    def may_be_immortal(self) -> bool:
        return any(item.may_be_immortal for item in self.items)

    def __repr__(self) -> str:
        return "<RUnion %s>" % ", ".join(str(item) for item in self.items)

    def __str__(self) -> str:
        return "union[%s]" % ", ".join(str(item) for item in self.items)

    # We compare based on the set because order in a union doesn't matter
    def __eq__(self, other: object) -> TypeGuard[RUnion]:
        return isinstance(other, RUnion) and self.items_set == other.items_set

    def __hash__(self) -> int:
        return hash(("union", self.items_set))

    def serialize(self) -> JsonDict:
        types = [x.serialize() for x in self.items]
        return {".class": "RUnion", "types": types}

    @classmethod
    def deserialize(cls, data: JsonDict, ctx: DeserMaps) -> RUnion:
        types = [deserialize_type(t, ctx) for t in data["types"]]
        return RUnion(types)


def flatten_nested_unions(types: list[RType]) -> list[RType]:
    if not any(isinstance(t, RUnion) for t in types):
        return types  # Fast path

    flat_items: list[RType] = []
    for t in types:
        if isinstance(t, RUnion):
            flat_items.extend(flatten_nested_unions(t.items))
        else:
            flat_items.append(t)
    return flat_items


def optional_value_type(rtype: RType) -> RType | None:
    """If rtype is the union of none_rprimitive and another type X, return X.

    Otherwise, return None.
    """
    if isinstance(rtype, RUnion) and len(rtype.items) == 2:
        if rtype.items[0] == none_rprimitive:
            return rtype.items[1]
        elif rtype.items[1] == none_rprimitive:
            return rtype.items[0]
    return None


def is_optional_type(rtype: RType) -> TypeGuard[RUnion]:
    """Is rtype an optional type with exactly two union items?"""
    return optional_value_type(rtype) is not None


@final
class RArray(RType):
    """Fixed-length C array type (for example, int[5]).

    Note that the implementation is a bit limited, and these can basically
    be only used for local variables that are initialized in one location.
    """

    def __init__(self, item_type: RType, length: int) -> None:
        self.item_type = item_type
        # Number of items
        self.length = length
        self.is_refcounted = False

    def accept(self, visitor: RTypeVisitor[T]) -> T:
        return visitor.visit_rarray(self)

    @property
    def may_be_immortal(self) -> bool:
        return False

    def __str__(self) -> str:
        return f"{self.item_type}[{self.length}]"

    def __repr__(self) -> str:
        return f"<RArray {self.item_type!r}[{self.length}]>"

    def __eq__(self, other: object) -> TypeGuard[RArray]:
        return (
            isinstance(other, RArray)
            and self.item_type == other.item_type
            and self.length == other.length
        )

    def __hash__(self) -> int:
        return hash((self.item_type, self.length))

    def serialize(self) -> JsonDict:
        assert False

    @classmethod
    def deserialize(cls, data: JsonDict, ctx: DeserMaps) -> RArray:
        assert False


PyObject = RStruct(
    name="PyObject",
    names=["ob_refcnt", "ob_type"],
    types=[c_pyssize_t_rprimitive, pointer_rprimitive],
)

PyVarObject = RStruct(
    name="PyVarObject", names=["ob_base", "ob_size"], types=[PyObject, c_pyssize_t_rprimitive]
)

setentry = RStruct(
    name="setentry", names=["key", "hash"], types=[pointer_rprimitive, c_pyssize_t_rprimitive]
)

smalltable = RStruct(name="smalltable", names=[], types=[setentry] * 8)

PySetObject = RStruct(
    name="PySetObject",
    names=[
        "ob_base",
        "fill",
        "used",
        "mask",
        "table",
        "hash",
        "finger",
        "smalltable",
        "weakreflist",
    ],
    types=[
        PyObject,
        c_pyssize_t_rprimitive,
        c_pyssize_t_rprimitive,
        c_pyssize_t_rprimitive,
        pointer_rprimitive,
        c_pyssize_t_rprimitive,
        c_pyssize_t_rprimitive,
        smalltable,
        pointer_rprimitive,
    ],
)

PyListObject = RStruct(
    name="PyListObject",
    names=["ob_base", "ob_item", "allocated"],
    types=[PyVarObject, pointer_rprimitive, c_pyssize_t_rprimitive],
)


def check_native_int_range(rtype: RPrimitive, n: int) -> bool:
    """Is n within the range of a native, fixed-width int type?

    Assume the type is a fixed-width int type.
    """
    if not rtype.is_signed:
        return 0 <= n < (1 << (8 * rtype.size))
    else:
        limit = 1 << (rtype.size * 8 - 1)
        return -limit <= n < limit


# Buffers for vec item types that have a packed representation
#
# Note that the 'items' fields are variable-length arrays, and mypyc IR isn't
# able to represent these, so the field type is omitted for now.

VecI64BufObject = RStruct(
    name="VecI64BufObject",
    names=["ob_base", "len", "items"],
    types=[PyVarObject, int64_rprimitive],
)

VecI32BufObject = RStruct(
    name="VecI32BufObject",
    names=["ob_base", "len", "items"],
    types=[PyVarObject, int64_rprimitive],
)

VecI16BufObject = RStruct(
    name="VecI16BufObject",
    names=["ob_base", "len", "items"],
    types=[PyVarObject, int64_rprimitive],
)

VecU8BufObject = RStruct(
    name="VecU8BufObject", names=["ob_base", "len", "items"], types=[PyVarObject, int64_rprimitive]
)

VecFloatBufObject = RStruct(
    name="VecFloatBufObject",
    names=["ob_base", "len", "items"],
    types=[PyVarObject, int64_rprimitive],
)

VecBoolBufObject = RStruct(
    name="VecBoolBufObject",
    names=["ob_base", "len", "items"],
    types=[PyVarObject, int64_rprimitive],
)


# Struct type for vec[i64] (in most cases use RVec instead).
VecI64 = RStruct(
    name="VecI64", names=["len", "items"], types=[c_pyssize_t_rprimitive, pointer_rprimitive]
)


# Buffer for vec[t]
VecTBufObject = RStruct(
    name="VecTBufObject",
    names=["ob_base", "item_type", "items"],
    types=[PyVarObject, c_pyssize_t_rprimitive, object_rprimitive],
)

# Struct type for vec[t] (in most cases use RVec instead).
VecT = RStruct(
    name="VecT", names=["len", "items"], types=[c_pyssize_t_rprimitive, pointer_rprimitive]
)

VecNestedBufItem = RStruct(
    name="VecNestedBufItem",
    names=["len", "items"],
    types=[c_pyssize_t_rprimitive, pointer_rprimitive],
)

# Buffer for vec[vec[t]]
VecNestedBufObject = RStruct(
    name="VecNestedBufObject",
    names=["ob_base", "item_type", "depth", "optionals", "items"],
    types=[
        PyVarObject,
        c_pyssize_t_rprimitive,
        int32_rprimitive,
        int32_rprimitive,
        VecNestedBufItem,
    ],
)

# Struct type for vec[vec[...]] (in most cases use RVec instead).
VecNested = RStruct(
    name="VecNested", names=["len", "items"], types=[c_pyssize_t_rprimitive, pointer_rprimitive]
)

VecNestedBufObject_rprimitive = RPrimitive(
    "VecNestedBufObject_ptr", is_unboxed=False, is_refcounted=True, ctype="VecNestedBufObject *"
)

VecNestedPopResult = RStruct(
    name="VecNestedPopResult", names=["vec", "item"], types=[VecNested, VecNestedBufItem]
)


vec_buf_types: Final[dict[RType, RStruct]] = {
    int64_rprimitive: VecI64BufObject,
    int32_rprimitive: VecI32BufObject,
    int16_rprimitive: VecI16BufObject,
    uint8_rprimitive: VecU8BufObject,
    float_rprimitive: VecFloatBufObject,
    bool_rprimitive: VecBoolBufObject,
}

vec_c_types: Final[dict[RType, str]] = {
    int64_rprimitive: "VecI64",
    int32_rprimitive: "VecI32",
    int16_rprimitive: "VecI16",
    uint8_rprimitive: "VecU8",
    float_rprimitive: "VecFloat",
    bool_rprimitive: "VecBool",
}

vec_api_fields: Final[dict[RType, str]] = {
    int64_rprimitive: "i64",
    int32_rprimitive: "i32",
    int16_rprimitive: "i16",
    uint8_rprimitive: "u8",
    float_rprimitive: "float_",
    bool_rprimitive: "bool_",
}

vec_api_by_item_type: Final[dict[RType, str]] = {
    int64_rprimitive: "VecI64Api",
    int32_rprimitive: "VecI32Api",
    int16_rprimitive: "VecI16Api",
    uint8_rprimitive: "VecU8Api",
    float_rprimitive: "VecFloatApi",
    bool_rprimitive: "VecBoolApi",
}

# These are special item type constants used in nested vecs to represent
# item types with specialized representations. These must match definitions
# in the vecs module (see VEC_ITEM_TYPE_I64 etc.).
vec_item_type_tags: Final[dict[RType, int]] = {
    int64_rprimitive: 2,
    int32_rprimitive: 6,
    int16_rprimitive: 10,
    uint8_rprimitive: 14,
    float_rprimitive: 18,
    bool_rprimitive: 22,
}
