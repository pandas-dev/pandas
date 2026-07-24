"""
Classes that are LLVM types
"""

import struct

from llvmlite import ir_layer_typed_pointers_enabled
from llvmlite.ir._utils import _StrCaching


def _wrapname(x):
    return '"{0}"'.format(x.replace('\\', '\\5c').replace('"', '\\22'))


class Type(_StrCaching):
    """
    The base class for all LLVM types.
    """
    is_pointer = False
    null = 'zeroinitializer'

    def __repr__(self):
        return "<%s %s>" % (type(self), str(self))

    def _to_string(self):
        raise NotImplementedError

    def as_pointer(self, addrspace=0):
        return PointerType(self, addrspace)

    def __ne__(self, other):
        return not (self == other)

    def _get_ll_global_value_type(self, target_data, context=None):
        """
        Convert this type object to an LLVM type.
        """
        from llvmlite.ir import Module, GlobalVariable
        from llvmlite.binding import parse_assembly

        if context is None:
            m = Module()
        else:
            m = Module(context=context)
        foo = GlobalVariable(m, self, name="foo")
        with parse_assembly(str(m)) as llmod:
            return llmod.get_global_variable(foo.name).global_value_type

    def get_abi_size(self, target_data, context=None):
        """
        Get the ABI size of this type according to data layout *target_data*.
        """
        llty = self._get_ll_global_value_type(target_data, context)
        return target_data.get_abi_size(llty)

    def get_element_offset(self, target_data, ndx, context=None):
        llty = self._get_ll_global_value_type(target_data, context)
        return target_data.get_element_offset(llty, ndx)

    def get_abi_alignment(self, target_data, context=None):
        """
        Get the minimum ABI alignment of this type according to data layout
        *target_data*.
        """
        llty = self._get_ll_global_value_type(target_data, context)
        return target_data.get_abi_alignment(llty)

    def format_constant(self, value):
        """
        Format constant *value* of this type.  This method may be overriden
        by subclasses.
        """
        return str(value)

    def wrap_constant_value(self, value):
        """
        Wrap constant *value* if necessary.  This method may be overriden
        by subclasses (especially aggregate types).
        """
        return value

    def __call__(self, value):
        """
        Create a LLVM constant of this type with the given Python value.
        """
        from llvmlite.ir import Constant
        return Constant(self, value)


class MetaDataType(Type):

    def _to_string(self):
        return "metadata"

    def as_pointer(self):
        raise TypeError

    def __eq__(self, other):
        return isinstance(other, MetaDataType)

    def __hash__(self):
        return hash(MetaDataType)


class LabelType(Type):
    """
    The label type is the type of e.g. basic blocks.
    """

    def _to_string(self):
        return "label"


class PointerType(Type):
    """
    The type of all pointer values.
    By default (without specialisation) represents an opaque pointer.
    """
    is_opaque = True
    is_pointer = True
    null = 'null'

    # Factory to create typed or opaque pointers based on `pointee'.
    def __new__(cls, pointee=None, addrspace=0):
        if cls is PointerType and pointee is not None and \
           type(pointee) is not PointerType:
            return super().__new__(_TypedPointerType)
        return super(PointerType, cls).__new__(cls)

    def __init__(self, pointee=None, addrspace=0):
        assert pointee is None or type(pointee) is PointerType
        self.addrspace = addrspace

    def _to_string(self):
        if self.addrspace != 0:
            return "ptr addrspace({0})".format(self.addrspace)
        else:
            return "ptr"

    def __eq__(self, other):
        return (isinstance(other, PointerType) and
                self.addrspace == other.addrspace)

    def __hash__(self):
        return hash(PointerType)

    @property
    def intrinsic_name(self):
        return 'p%d' % self.addrspace

    @classmethod
    def from_llvm(cls, typeref, ir_ctx):
        """
        Create from a llvmlite.binding.TypeRef
        """
        return cls()


class _TypedPointerType(PointerType):
    """
    The type of typed pointer values. To be removed eventually.
    """

    def __init__(self, pointee, addrspace=0):
        super(_TypedPointerType, self).__init__(None, addrspace)
        assert pointee is not None and type(pointee) is not PointerType
        assert not isinstance(pointee, VoidType)
        self.pointee = pointee
        self.is_opaque = False

    def _to_string(self):
        if ir_layer_typed_pointers_enabled:
            return "{0}*".format(self.pointee) if self.addrspace == 0 else \
                   "{0} addrspace({1})*".format(self.pointee, self.addrspace)
        return super(_TypedPointerType, self)._to_string()

    # This implements ``isOpaqueOrPointeeTypeEquals''.
    def __eq__(self, other):
        if isinstance(other, _TypedPointerType):
            return (self.pointee, self.addrspace) == (other.pointee,
                                                      other.addrspace)
        return (isinstance(other, PointerType) and
                self.addrspace == other.addrspace)

    def __hash__(self):
        return hash(_TypedPointerType)

    def gep(self, i):
        """
        Resolve the type of the i-th element (for getelementptr lookups).
        """
        if not isinstance(i.type, IntType):
            raise TypeError(i.type)
        return self.pointee

    @property
    def intrinsic_name(self):
        if ir_layer_typed_pointers_enabled:
            return 'p%d%s' % (self.addrspace, self.pointee.intrinsic_name)
        return super(_TypedPointerType, self).intrinsic_name


class VoidType(Type):
    """
    The type for empty values (e.g. a function returning no value).
    """

    def _to_string(self):
        return 'void'

    def __eq__(self, other):
        return isinstance(other, VoidType)

    def __hash__(self):
        return hash(VoidType)

    @classmethod
    def from_llvm(cls, typeref, ir_ctx):
        """
        Create from a llvmlite.binding.TypeRef
        """
        return cls()


class FunctionType(Type):
    """
    The type for functions.
    """

    def __init__(self, return_type, args, var_arg=False):
        self.return_type = return_type
        self.args = tuple(args)
        self.var_arg = var_arg

    def _to_string(self):
        if self.args:
            strargs = ', '.join([str(a) for a in self.args])
            if self.var_arg:
                return '{0} ({1}, ...)'.format(self.return_type, strargs)
            else:
                return '{0} ({1})'.format(self.return_type, strargs)
        elif self.var_arg:
            return '{0} (...)'.format(self.return_type)
        else:
            return '{0} ()'.format(self.return_type)

    def __eq__(self, other):
        if isinstance(other, FunctionType):
            return (self.return_type == other.return_type and
                    self.args == other.args and self.var_arg == other.var_arg)
        else:
            return False

    def __hash__(self):
        return hash(FunctionType)

    @classmethod
    def from_llvm(cls, typeref, ir_ctx):
        """
        Create from a llvmlite.binding.TypeRef
        """
        params = tuple(x.as_ir(ir_ctx=ir_ctx)
                       for x in typeref.get_function_parameters())
        ret = typeref.get_function_return().as_ir(ir_ctx=ir_ctx)
        is_vararg = typeref.is_function_vararg
        return cls(ret, params, is_vararg)


class IntType(Type):
    """
    The type for integers.
    """
    null = '0'
    _instance_cache = {}
    width: int

    def __new__(cls, bits):
        # Cache all common integer types
        if 0 <= bits <= 128:
            try:
                return cls._instance_cache[bits]
            except KeyError:
                inst = cls._instance_cache[bits] = cls.__new(bits)
                return inst
        return cls.__new(bits)

    @classmethod
    def __new(cls, bits):
        assert isinstance(bits, int) and bits >= 0
        self = super(IntType, cls).__new__(cls)
        self.width = bits
        return self

    def __getnewargs__(self):
        return self.width,

    def __copy__(self):
        return self

    def _to_string(self):
        return 'i%u' % (self.width,)

    def __eq__(self, other):
        if isinstance(other, IntType):
            return self.width == other.width
        else:
            return False

    def __hash__(self):
        return hash(IntType)

    def format_constant(self, val):
        if isinstance(val, bool):
            return str(val).lower()
        else:
            return str(val)

    def wrap_constant_value(self, val):
        if val is None:
            return 0
        return val

    @property
    def intrinsic_name(self):
        return str(self)

    @classmethod
    def from_llvm(cls, typeref, ir_ctx):
        """
        Create from a llvmlite.binding.TypeRef
        """
        return IntType(typeref.type_width)


def _as_float(value):
    """
    Truncate to single-precision float.
    """
    return struct.unpack('f', struct.pack('f', value))[0]


def _as_half(value):
    """
    Truncate to half-precision float.
    """
    try:
        return struct.unpack('e', struct.pack('e', value))[0]
    except struct.error:
        # 'e' only added in Python 3.6+
        return _as_float(value)


def _format_float_as_hex(value, packfmt, unpackfmt, numdigits):
    raw = struct.pack(packfmt, float(value))
    intrep = struct.unpack(unpackfmt, raw)[0]
    out = '{{0:#{0}x}}'.format(numdigits).format(intrep)
    return out


def _format_double(value):
    """
    Format *value* as a hexadecimal string of its IEEE double precision
    representation.
    """
    return _format_float_as_hex(value, 'd', 'Q', 16)


class _BaseFloatType(Type):

    def __new__(cls):
        return cls._instance_cache

    def __eq__(self, other):
        return isinstance(other, type(self))

    def __hash__(self):
        return hash(type(self))

    @classmethod
    def _create_instance(cls):
        cls._instance_cache = super(_BaseFloatType, cls).__new__(cls)

    @classmethod
    def from_llvm(cls, typeref, ir_ctx):
        """
        Create from a llvmlite.binding.TypeRef
        """
        return cls()


class HalfType(_BaseFloatType):
    """
    The type for single-precision floats.
    """
    null = '0.0'
    intrinsic_name = 'f16'

    def __str__(self):
        return 'half'

    def format_constant(self, value):
        return _format_double(_as_half(value))


class FloatType(_BaseFloatType):
    """
    The type for single-precision floats.
    """
    null = '0.0'
    intrinsic_name = 'f32'

    def __str__(self):
        return 'float'

    def format_constant(self, value):
        return _format_double(_as_float(value))


class DoubleType(_BaseFloatType):
    """
    The type for double-precision floats.
    """
    null = '0.0'
    intrinsic_name = 'f64'

    def __str__(self):
        return 'double'

    def format_constant(self, value):
        return _format_double(value)


for _cls in (HalfType, FloatType, DoubleType):
    _cls._create_instance()


class _Repeat(object):
    def __init__(self, value, size):
        self.value = value
        self.size = size

    def __len__(self):
        return self.size

    def __getitem__(self, item):
        if 0 <= item < self.size:
            return self.value
        else:
            raise IndexError(item)


class VectorType(Type):
    """
    The type for vectors of primitive data items (e.g. "<f32 x 4>").
    """

    def __init__(self, element, count):
        self.element = element
        self.count = count

    @property
    def elements(self):
        return _Repeat(self.element, self.count)

    def __len__(self):
        return self.count

    def _to_string(self):
        return "<%d x %s>" % (self.count, self.element)

    def __eq__(self, other):
        if isinstance(other, VectorType):
            return self.element == other.element and self.count == other.count

    def __hash__(self):
        # TODO: why does this not take self.element/self.count into account?
        return hash(VectorType)

    def __copy__(self):
        return self

    def format_constant(self, value):
        itemstring = ", " .join(["{0} {1}".format(x.type, x.get_reference())
                                 for x in value])
        return "<{0}>".format(itemstring)

    def wrap_constant_value(self, values):
        from . import Value, Constant
        if not isinstance(values, (list, tuple)):
            if isinstance(values, Constant):
                if values.type != self.element:
                    raise TypeError("expected {} for {}".format(
                        self.element, values.type))
                return (values, ) * self.count
            return (Constant(self.element, values), ) * self.count
        if len(values) != len(self):
            raise ValueError("wrong constant size for %s: got %d elements"
                             % (self, len(values)))
        return [Constant(ty, val) if not isinstance(val, Value) else val
                for ty, val in zip(self.elements, values)]

    @classmethod
    def from_llvm(cls, typeref, ir_ctx):
        """
        Create from a llvmlite.binding.TypeRef
        """
        [elemtyperef] = typeref.elements
        elemty = elemtyperef.as_ir(ir_ctx=ir_ctx)
        count = typeref.element_count
        return cls(elemty, count)


class Aggregate(Type):
    """
    Base class for aggregate types.
    See http://llvm.org/docs/LangRef.html#t-aggregate
    """

    def wrap_constant_value(self, values):
        from . import Value, Constant

        if not isinstance(values, (list, tuple)):
            return values
        if len(values) != len(self):
            raise ValueError("wrong constant size for %s: got %d elements"
                             % (self, len(values)))
        return [Constant(ty, val) if not isinstance(val, Value) else val
                for ty, val in zip(self.elements, values)]


class ArrayType(Aggregate):
    """
    The type for fixed-size homogenous arrays (e.g. "[f32 x 3]").
    """

    def __init__(self, element, count):
        self.element = element
        self.count = count

    @property
    def elements(self):
        return _Repeat(self.element, self.count)

    def __len__(self):
        return self.count

    def _to_string(self):
        return "[%d x %s]" % (self.count, self.element)

    def __eq__(self, other):
        if isinstance(other, ArrayType):
            return self.element == other.element and self.count == other.count

    def __hash__(self):
        return hash(ArrayType)

    def gep(self, i):
        """
        Resolve the type of the i-th element (for getelementptr lookups).
        """
        if not isinstance(i.type, IntType):
            raise TypeError(i.type)
        return self.element

    def format_constant(self, value):
        itemstring = ", " .join(["{0} {1}".format(x.type, x.get_reference())
                                 for x in value])
        return "[{0}]".format(itemstring)

    @classmethod
    def from_llvm(cls, typeref, ir_ctx):
        """
        Create from a llvmlite.binding.TypeRef
        """
        [elemtyperef] = typeref.elements
        elemty = elemtyperef.as_ir(ir_ctx=ir_ctx)
        count = typeref.element_count
        return cls(elemty, count)


class BaseStructType(Aggregate):
    """
    The base type for heterogenous struct types.
    """
    _packed = False

    @property
    def packed(self):
        """
        A boolean attribute that indicates whether the structure uses
        packed layout.
        """
        return self._packed

    @packed.setter
    def packed(self, val):
        self._packed = bool(val)

    def __len__(self):
        assert self.elements is not None
        return len(self.elements)

    def __iter__(self):
        assert self.elements is not None
        return iter(self.elements)

    @property
    def is_opaque(self):
        return self.elements is None

    def structure_repr(self):
        """
        Return the LLVM IR for the structure representation
        """
        ret = '{%s}' % ', '.join([str(x) for x in self.elements])
        return self._wrap_packed(ret)

    def format_constant(self, value):
        itemstring = ", " .join(["{0} {1}".format(x.type, x.get_reference())
                                 for x in value])
        ret = "{{{0}}}".format(itemstring)
        return self._wrap_packed(ret)

    def gep(self, i):
        """
        Resolve the type of the i-th element (for getelementptr lookups).

        *i* needs to be a LLVM constant, so that the type can be determined
        at compile-time.
        """
        if not isinstance(i.type, IntType):
            raise TypeError(i.type)
        return self.elements[i.constant]

    def _wrap_packed(self, textrepr):
        """
        Internal helper to wrap textual repr of struct type into packed struct
        """
        if self.packed:
            return '<{}>'.format(textrepr)
        else:
            return textrepr

    @classmethod
    def from_llvm(cls, typeref, ir_ctx):
        """
        Create from a llvmlite.binding.TypeRef
        """
        if typeref.is_literal_struct:
            elems = [el.as_ir(ir_ctx=ir_ctx) for el in typeref.elements]
            return cls(elems, typeref.is_packed_struct)
        else:
            return ir_ctx.get_identified_type(typeref.name)


class LiteralStructType(BaseStructType):
    """
    The type of "literal" structs, i.e. structs with a literally-defined
    type (by contrast with IdentifiedStructType).
    """

    null = 'zeroinitializer'

    def __init__(self, elems, packed=False):
        """
        *elems* is a sequence of types to be used as members.
        *packed* controls the use of packed layout.
        """
        self.elements = tuple(elems)
        self.packed = packed

    def _to_string(self):
        return self.structure_repr()

    def __eq__(self, other):
        if isinstance(other, LiteralStructType):
            return (self.elements == other.elements
                    and self.packed == other.packed)

    def __hash__(self):
        return hash(LiteralStructType)


class IdentifiedStructType(BaseStructType):
    """
    A type which is a named alias for another struct type, akin to a typedef.
    While literal struct types can be structurally equal (see
    LiteralStructType), identified struct types are compared by name.

    Do not use this directly.
    """
    null = 'zeroinitializer'

    def __init__(self, context, name, packed=False):
        """
        *context* is a llvmlite.ir.Context.
        *name* is the identifier for the new struct type.
        *packed* controls the use of packed layout.
        """
        assert name
        self.context = context
        self.name = name
        self.elements = None
        self.packed = packed

    def _to_string(self):
        return "%{name}".format(name=_wrapname(self.name))

    def get_declaration(self):
        """
        Returns the string for the declaration of the type
        """
        if self.is_opaque:
            out = "{strrep} = type opaque".format(strrep=str(self))
        else:
            out = "{strrep} = type {struct}".format(
                strrep=str(self), struct=self.structure_repr())
        return out

    def __eq__(self, other):
        if isinstance(other, IdentifiedStructType):
            return (self.name == other.name
                    and self.packed == other.packed)

    def __hash__(self):
        return hash(IdentifiedStructType)

    def set_body(self, *elems):
        if not self.is_opaque:
            raise RuntimeError("{name} is already defined".format(
                name=self.name))
        self.elements = tuple(elems)
