"""
Classes that are LLVM values: Value, Constant...
Instructions are in the instructions module.
"""

import functools
import string
import re
from types import MappingProxyType

from llvmlite.ir import values, types, _utils
from llvmlite.ir._utils import (_StrCaching, _StringReferenceCaching,
                                _HasMetadata)

_VALID_CHARS = (frozenset(map(ord, string.ascii_letters)) |
                frozenset(map(ord, string.digits)) |
                frozenset(map(ord, ' !#$%&\'()*+,-./:;<=>?@[]^_`{|}~')))

_SIMPLE_IDENTIFIER_RE = re.compile(r"[-a-zA-Z$._][-a-zA-Z$._0-9]*$")

_CMP_MAP = {
    '>': 'gt',
    '<': 'lt',
    '==': 'eq',
    '!=': 'ne',
    '>=': 'ge',
    '<=': 'le',
}


def _escape_string(text, _map={}):
    """
    Escape the given bytestring for safe use as a LLVM array constant.
    Any unicode string input is first encoded with utf8 into bytes.
    """
    if isinstance(text, str):
        text = text.encode()
    assert isinstance(text, (bytes, bytearray))

    if not _map:
        for ch in range(256):
            if ch in _VALID_CHARS:
                _map[ch] = chr(ch)
            else:
                _map[ch] = '\\%02x' % ch

    buf = [_map[ch] for ch in text]
    return ''.join(buf)


def _binop(opname):
    def wrap(fn):
        @functools.wraps(fn)
        def wrapped(lhs, rhs):
            if lhs.type != rhs.type:
                raise ValueError("Operands must be the same type, got (%s, %s)"
                                 % (lhs.type, rhs.type))

            fmt = "{0} ({1} {2}, {3} {4})".format(opname,
                                                  lhs.type, lhs.get_reference(),
                                                  rhs.type, rhs.get_reference())
            return FormattedConstant(lhs.type, fmt)

        return wrapped
    return wrap


def _castop(opname):
    def wrap(fn):
        @functools.wraps(fn)
        def wrapped(self, typ):
            fn(self, typ)
            if typ == self.type:
                return self

            op = "{0} ({1} {2} to {3})".format(opname, self.type,
                                               self.get_reference(), typ)
            return FormattedConstant(typ, op)

        return wrapped
    return wrap


class _ConstOpMixin(object):
    """
    A mixin defining constant operations, for use in constant-like classes.
    """

    #
    # Arithmetic APIs
    #

    @_binop('shl')
    def shl(self, other):
        """
        Left integer shift:
            lhs << rhs
        """

    @_binop('lshr')
    def lshr(self, other):
        """
        Logical (unsigned) right integer shift:
            lhs >> rhs
        """

    @_binop('ashr')
    def ashr(self, other):
        """
        Arithmetic (signed) right integer shift:
            lhs >> rhs
        """

    @_binop('add')
    def add(self, other):
        """
        Integer addition:
            lhs + rhs
        """

    @_binop('fadd')
    def fadd(self, other):
        """
        Floating-point addition:
            lhs + rhs
        """

    @_binop('sub')
    def sub(self, other):
        """
        Integer subtraction:
            lhs - rhs
        """

    @_binop('fsub')
    def fsub(self, other):
        """
        Floating-point subtraction:
            lhs - rhs
        """

    @_binop('mul')
    def mul(self, other):
        """
        Integer multiplication:
            lhs * rhs
        """

    @_binop('fmul')
    def fmul(self, other):
        """
        Floating-point multiplication:
            lhs * rhs
        """

    @_binop('udiv')
    def udiv(self, other):
        """
        Unsigned integer division:
            lhs / rhs
        """

    @_binop('sdiv')
    def sdiv(self, other):
        """
        Signed integer division:
            lhs / rhs
        """

    @_binop('fdiv')
    def fdiv(self, other):
        """
        Floating-point division:
            lhs / rhs
        """

    @_binop('urem')
    def urem(self, other):
        """
        Unsigned integer remainder:
            lhs % rhs
        """

    @_binop('srem')
    def srem(self, other):
        """
        Signed integer remainder:
            lhs % rhs
        """

    @_binop('frem')
    def frem(self, other):
        """
        Floating-point remainder:
            lhs % rhs
        """

    @_binop('or')
    def or_(self, other):
        """
        Bitwise integer OR:
            lhs | rhs
        """

    @_binop('and')
    def and_(self, other):
        """
        Bitwise integer AND:
            lhs & rhs
        """

    @_binop('xor')
    def xor(self, other):
        """
        Bitwise integer XOR:
            lhs ^ rhs
        """

    def _cmp(self, prefix, sign, cmpop, other):
        ins = prefix + 'cmp'
        try:
            op = _CMP_MAP[cmpop]
        except KeyError:
            raise ValueError("invalid comparison %r for %s" % (cmpop, ins))

        if not (prefix == 'i' and cmpop in ('==', '!=')):
            op = sign + op

        if self.type != other.type:
            raise ValueError("Operands must be the same type, got (%s, %s)"
                             % (self.type, other.type))

        fmt = "{0} {1} ({2} {3}, {4} {5})".format(
            ins, op,
            self.type, self.get_reference(),
            other.type, other.get_reference())

        return FormattedConstant(types.IntType(1), fmt)

    def icmp_signed(self, cmpop, other):
        """
        Signed integer comparison:
            lhs <cmpop> rhs

        where cmpop can be '==', '!=', '<', '<=', '>', '>='
        """
        return self._cmp('i', 's', cmpop, other)

    def icmp_unsigned(self, cmpop, other):
        """
        Unsigned integer (or pointer) comparison:
            lhs <cmpop> rhs

        where cmpop can be '==', '!=', '<', '<=', '>', '>='
        """
        return self._cmp('i', 'u', cmpop, other)

    def fcmp_ordered(self, cmpop, other):
        """
        Floating-point ordered comparison:
            lhs <cmpop> rhs

        where cmpop can be '==', '!=', '<', '<=', '>', '>=', 'ord', 'uno'
        """
        return self._cmp('f', 'o', cmpop, other)

    def fcmp_unordered(self, cmpop, other):
        """
        Floating-point unordered comparison:
            lhs <cmpop> rhs

        where cmpop can be '==', '!=', '<', '<=', '>', '>=', 'ord', 'uno'
        """
        return self._cmp('f', 'u', cmpop, other)

    #
    # Unary APIs
    #

    def not_(self):
        """
        Bitwise integer complement:
            ~value
        """
        if isinstance(self.type, types.VectorType):
            rhs = values.Constant(self.type, (-1,) * self.type.count)
        else:
            rhs = values.Constant(self.type, -1)

        return self.xor(rhs)

    def neg(self):
        """
        Integer negative:
            -value
        """
        zero = values.Constant(self.type, 0)
        return zero.sub(self)

    def fneg(self):
        """
        Floating-point negative:
            -value
        """
        fmt = "fneg ({0} {1})".format(self.type, self.get_reference())
        return FormattedConstant(self.type, fmt)

    #
    # Cast APIs
    #

    @_castop('trunc')
    def trunc(self, typ):
        """
        Truncating integer downcast to a smaller type.
        """

    @_castop('zext')
    def zext(self, typ):
        """
        Zero-extending integer upcast to a larger type
        """

    @_castop('sext')
    def sext(self, typ):
        """
        Sign-extending integer upcast to a larger type.
        """

    @_castop('fptrunc')
    def fptrunc(self, typ):
        """
        Floating-point downcast to a less precise type.
        """

    @_castop('fpext')
    def fpext(self, typ):
        """
        Floating-point upcast to a more precise type.
        """

    @_castop('bitcast')
    def bitcast(self, typ):
        """
        Pointer cast to a different pointer type.
        """

    @_castop('fptoui')
    def fptoui(self, typ):
        """
        Convert floating-point to unsigned integer.
        """

    @_castop('uitofp')
    def uitofp(self, typ):
        """
        Convert unsigned integer to floating-point.
        """

    @_castop('fptosi')
    def fptosi(self, typ):
        """
        Convert floating-point to signed integer.
        """

    @_castop('sitofp')
    def sitofp(self, typ):
        """
        Convert signed integer to floating-point.
        """

    @_castop('ptrtoint')
    def ptrtoint(self, typ):
        """
        Cast pointer to integer.
        """
        if not isinstance(self.type, types.PointerType):
            msg = "can only call ptrtoint() on pointer type, not '%s'"
            raise TypeError(msg % (self.type,))
        if not isinstance(typ, types.IntType):
            raise TypeError("can only ptrtoint() to integer type, not '%s'"
                            % (typ,))

    @_castop('inttoptr')
    def inttoptr(self, typ):
        """
        Cast integer to pointer.
        """
        if not isinstance(self.type, types.IntType):
            msg = "can only call inttoptr() on integer constants, not '%s'"
            raise TypeError(msg % (self.type,))
        if not isinstance(typ, types.PointerType):
            raise TypeError("can only inttoptr() to pointer type, not '%s'"
                            % (typ,))

    def gep(self, indices):
        """
        Call getelementptr on this pointer constant.
        """
        if not isinstance(self.type, types.PointerType):
            raise TypeError("can only call gep() on pointer constants, not '%s'"
                            % (self.type,))

        outtype = self.type
        for i in indices:
            outtype = outtype.gep(i)

        strindices = ["{0} {1}".format(idx.type, idx.get_reference())
                      for idx in indices]

        op = "getelementptr ({0}, {1} {2}, {3})".format(
            self.type.pointee, self.type,
            self.get_reference(), ', '.join(strindices))
        return FormattedConstant(outtype.as_pointer(self.addrspace), op)


class Value(object):
    """
    The base class for all values.
    """

    def __repr__(self):
        return "<ir.%s type='%s' ...>" % (self.__class__.__name__, self.type,)


class _Undefined(object):
    """
    'undef': a value for undefined values.
    """
    def __new__(cls):
        try:
            return Undefined
        except NameError:
            return object.__new__(_Undefined)


Undefined = _Undefined()


class Constant(_StrCaching, _StringReferenceCaching, _ConstOpMixin, Value):
    """
    A constant LLVM value.
    """

    def __init__(self, typ, constant):
        assert isinstance(typ, types.Type)
        assert not isinstance(typ, types.VoidType)
        self.type = typ
        constant = typ.wrap_constant_value(constant)
        self.constant = constant

    def _to_string(self):
        return '{0} {1}'.format(self.type, self.get_reference())

    def _get_reference(self):
        if self.constant is None:
            val = self.type.null

        elif self.constant is Undefined:
            val = "undef"

        elif isinstance(self.constant, bytearray):
            val = 'c"{0}"'.format(_escape_string(self.constant))

        else:
            val = self.type.format_constant(self.constant)

        return val

    @classmethod
    def literal_array(cls, elems):
        """
        Construct a literal array constant made of the given members.
        """
        tys = [el.type for el in elems]
        if len(tys) == 0:
            raise ValueError("need at least one element")
        ty = tys[0]
        for other in tys:
            if ty != other:
                raise TypeError("all elements must have the same type")
        return cls(types.ArrayType(ty, len(elems)), elems)

    @classmethod
    def literal_struct(cls, elems, packed=False):
        """
        Construct a literal structure constant made of the given members.
        """
        tys = [el.type for el in elems]
        return cls(types.LiteralStructType(tys, packed), elems)

    @property
    def addrspace(self):
        if not isinstance(self.type, types.PointerType):
            raise TypeError("Only pointer constant have address spaces")
        return self.type.addrspace

    def __eq__(self, other):
        if isinstance(other, Constant):
            return str(self) == str(other)
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(str(self))

    def __repr__(self):
        return "<ir.Constant type='%s' value=%r>" % (self.type, self.constant)


class FormattedConstant(Constant):
    """
    A constant with an already formatted IR representation.
    """

    def __init__(self, typ, constant):
        assert isinstance(constant, str)
        Constant.__init__(self, typ, constant)

    def _to_string(self):
        return self.constant

    def _get_reference(self):
        return self.constant


class NamedValue(_StrCaching, _StringReferenceCaching, Value):
    """
    The base class for named values.
    """
    name_prefix = '%'
    deduplicate_name = True

    def __init__(self, parent, type, name):
        assert parent is not None
        assert isinstance(type, types.Type)
        self.parent = parent
        self.type = type
        self._set_name(name)

    def _to_string(self):
        buf = []
        if not isinstance(self.type, types.VoidType):
            buf.append("{0} = ".format(self.get_reference()))
        self.descr(buf)
        return "".join(buf).rstrip()

    def descr(self, buf):
        raise NotImplementedError

    def _get_name(self):
        return self._name

    def _set_name(self, name):
        name = self.parent.scope.register(name,
                                          deduplicate=self.deduplicate_name)
        self._name = name

    name = property(_get_name, _set_name)

    def _get_reference(self):
        name = self.name
        # Quote and escape value name
        if '\\' in name or '"' in name:
            name = name.replace('\\', '\\5c').replace('"', '\\22')
        return '{0}"{1}"'.format(self.name_prefix, name)

    def __repr__(self):
        return "<ir.%s %r of type '%s'>" % (
            self.__class__.__name__, self.name, self.type)

    @property
    def function_type(self):
        ty = self.type
        if isinstance(ty, types.PointerType):
            ty = self.type.pointee
        if isinstance(ty, types.FunctionType):
            return ty
        else:
            raise TypeError("Not a function: {0}".format(self.type))


class MetaDataString(NamedValue):
    """
    A metadata string, i.e. a constant string used as a value in a metadata
    node.
    """

    def __init__(self, parent, string):
        super(MetaDataString, self).__init__(parent,
                                             types.MetaDataType(),
                                             name="")
        self.string = string

    def descr(self, buf):
        buf += (self.get_reference(), "\n")

    def _get_reference(self):
        return '!"{0}"'.format(_escape_string(self.string))

    _to_string = _get_reference

    def __eq__(self, other):
        if isinstance(other, MetaDataString):
            return self.string == other.string
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.string)


class MetaDataArgument(_StrCaching, _StringReferenceCaching, Value):
    """
    An argument value to a function taking metadata arguments.
    This can wrap any other kind of LLVM value.

    Do not instantiate directly, Builder.call() will create these
    automatically.
    """

    def __init__(self, value):
        assert isinstance(value, Value)
        assert not isinstance(value.type, types.MetaDataType)
        self.type = types.MetaDataType()
        self.wrapped_value = value

    def _get_reference(self):
        # e.g. "i32* %2"
        return "{0} {1}".format(self.wrapped_value.type,
                                self.wrapped_value.get_reference())

    _to_string = _get_reference


class NamedMetaData(object):
    """
    A named metadata node.

    Do not instantiate directly, use Module.add_named_metadata() instead.
    """

    def __init__(self, parent):
        self.parent = parent
        self.operands = []

    def add(self, md):
        self.operands.append(md)


class MDValue(NamedValue):
    """
    A metadata node's value, consisting of a sequence of elements ("operands").

    Do not instantiate directly, use Module.add_metadata() instead.
    """
    name_prefix = '!'

    def __init__(self, parent, values, name):
        super(MDValue, self).__init__(parent,
                                      types.MetaDataType(),
                                      name=name)
        self.operands = tuple(values)
        parent.metadata.append(self)

    def descr(self, buf):
        operands = []
        for op in self.operands:
            if isinstance(op.type, types.MetaDataType):
                if isinstance(op, Constant) and op.constant is None:
                    operands.append("null")
                else:
                    operands.append(op.get_reference())
            else:
                operands.append("{0} {1}".format(op.type, op.get_reference()))
        operands = ', '.join(operands)
        buf += ("!{{ {0} }}".format(operands), "\n")

    def _get_reference(self):
        return self.name_prefix + str(self.name)

    def __eq__(self, other):
        if isinstance(other, MDValue):
            return self.operands == other.operands
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.operands)


class DIToken:
    """
    A debug information enumeration value that should appear bare in
    the emitted metadata.

    Use this to wrap known constants, e.g. the DW_* enumerations.
    """

    def __init__(self, value):
        self.value = value


class DIValue(NamedValue):
    """
    A debug information descriptor, containing key-value pairs.

    Do not instantiate directly, use Module.add_debug_info() instead.
    """
    name_prefix = '!'

    def __init__(self, parent, is_distinct, kind, operands, name):
        super(DIValue, self).__init__(parent,
                                      types.MetaDataType(),
                                      name=name)
        self.is_distinct = is_distinct
        self.kind = kind
        self.operands = tuple(operands)
        parent.metadata.append(self)

    def descr(self, buf):
        if self.is_distinct:
            buf += ("distinct ",)
        operands = []
        for key, value in self.operands:
            if value is None:
                strvalue = "null"
            elif value is True:
                strvalue = "true"
            elif value is False:
                strvalue = "false"
            elif isinstance(value, DIToken):
                strvalue = value.value
            elif isinstance(value, str):
                strvalue = '"{}"'.format(_escape_string(value))
            elif isinstance(value, int):
                strvalue = str(value)
            elif isinstance(value, Constant):
                # Support for typed constants (e.g., i32 1, i8 2)
                # This calls Constant._to_string() which formats as "type value"
                strvalue = str(value)
            elif isinstance(value, NamedValue):
                strvalue = value.get_reference()
            else:
                raise TypeError("invalid operand type for debug info: %r"
                                % (value,))
            operands.append("{0}: {1}".format(key, strvalue))
        operands = ', '.join(operands)
        buf += ("!", self.kind, "(", operands, ")\n")

    def _get_reference(self):
        return self.name_prefix + str(self.name)

    def __eq__(self, other):
        if isinstance(other, DIValue):
            return self.is_distinct == other.is_distinct and \
                self.kind == other.kind and \
                self.operands == other.operands
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.is_distinct, self.kind, self.operands))


class GlobalValue(NamedValue, _ConstOpMixin, _HasMetadata):
    """
    A global value.
    """
    name_prefix = '@'
    deduplicate_name = False

    def __init__(self, *args, **kwargs):
        super(GlobalValue, self).__init__(*args, **kwargs)
        self.linkage = ''
        self.storage_class = ''
        self.section = ''
        self.metadata = {}


class GlobalVariable(GlobalValue):
    """
    A global variable.
    """

    def __init__(self, module, typ, name, addrspace=0):
        assert isinstance(typ, types.Type)
        super(GlobalVariable, self).__init__(module, typ.as_pointer(addrspace),
                                             name=name)
        self.value_type = typ
        self.initializer = None
        self.unnamed_addr = False
        self.global_constant = False
        self.addrspace = addrspace
        self.align = None
        self.parent.add_global(self)

    def descr(self, buf):
        if self.global_constant:
            kind = 'constant'
        else:
            kind = 'global'

        if not self.linkage:
            # Default to external linkage
            linkage = 'external' if self.initializer is None else ''
        else:
            linkage = self.linkage

        if linkage:
            buf.append(linkage + " ")
        if self.storage_class:
            buf.append(self.storage_class + " ")
        if self.unnamed_addr:
            buf.append("unnamed_addr ")
        if self.addrspace != 0:
            buf.append('addrspace({0:d}) '.format(self.addrspace))

        buf.append("{kind} {type}" .format(kind=kind, type=self.value_type))

        if self.initializer is not None:
            if self.initializer.type != self.value_type:
                raise TypeError("got initializer of type %s "
                                "for global value type %s"
                                % (self.initializer.type, self.value_type))
            buf.append(" " + self.initializer.get_reference())
        elif linkage not in ('external', 'extern_weak'):
            # emit 'undef' for non-external linkage GV
            buf.append(" " + self.value_type(Undefined).get_reference())

        if self.section:
            buf.append(", section \"%s\"" % (self.section,))

        if self.align is not None:
            buf.append(", align %d" % (self.align,))

        if self.metadata:
            buf.append(self._stringify_metadata(leading_comma=True))

        buf.append("\n")


class AttributeSet(set):
    """A set of string attribute.
    Only accept items listed in *_known*.

    Properties:
    * Iterate in sorted order
    """
    _known = ()

    def __init__(self, args=()):
        super().__init__()
        if isinstance(args, str):
            args = [args]
        for name in args:
            self.add(name)

    def _expand(self, name, typ):
        return name

    def add(self, name):
        if name not in self._known:
            raise ValueError('unknown attr {!r} for {}'.format(name, self))
        return super(AttributeSet, self).add(name)

    def _to_list(self, typ):
        return [self._expand(i, typ) for i in sorted(self)]


class FunctionAttributes(AttributeSet):
    _known = frozenset([
        'argmemonly', 'alwaysinline', 'builtin', 'cold', 'convergent',
        'inaccessiblememonly', 'inaccessiblemem_or_argmemonly', 'inlinehint',
        'jumptable', 'minsize', 'naked', 'nobuiltin', 'noduplicate',
        'noimplicitfloat', 'noinline', 'nonlazybind', 'norecurse',
        'noredzone', 'noreturn', 'nounwind', 'optnone', 'optsize',
        'readnone', 'readonly', 'returns_twice', 'sanitize_address',
        'sanitize_memory', 'sanitize_thread', 'ssp',
        'sspreg', 'sspstrong', 'uwtable'])

    def __init__(self, args=()):
        self._alignstack = 0
        self._personality = None
        super(FunctionAttributes, self).__init__(args)

    def add(self, name):
        if ((name == 'alwaysinline' and 'noinline' in self) or
                (name == 'noinline' and 'alwaysinline' in self)):
            raise ValueError("Can't have alwaysinline and noinline")

        super().add(name)

    @property
    def alignstack(self):
        return self._alignstack

    @alignstack.setter
    def alignstack(self, val):
        assert val >= 0
        self._alignstack = val

    @property
    def personality(self):
        return self._personality

    @personality.setter
    def personality(self, val):
        assert val is None or isinstance(val, GlobalValue)
        self._personality = val

    def _to_list(self, ret_type):
        attrs = super()._to_list(ret_type)
        if self.alignstack:
            attrs.append('alignstack({0:d})'.format(self.alignstack))
        if self.personality:
            attrs.append('personality {persty} {persfn}'.format(
                persty=self.personality.type,
                persfn=self.personality.get_reference()))
        return attrs


class Function(GlobalValue):
    """Represent a LLVM Function but does uses a Module as parent.
    Global Values are stored as a set of dependencies (attribute `depends`).
    """

    def __init__(self, module, ftype, name):
        assert isinstance(ftype, types.Type)
        super(Function, self).__init__(module, ftype.as_pointer(), name=name)
        self.ftype = ftype
        self.scope = _utils.NameScope()
        self.blocks = []
        self.attributes = FunctionAttributes()
        self.args = tuple([Argument(self, t)
                           for t in ftype.args])
        self.return_value = ReturnValue(self, ftype.return_type)
        self.parent.add_global(self)
        self.calling_convention = ''

    @property
    def module(self):
        return self.parent

    @property
    def entry_basic_block(self):
        return self.blocks[0]

    @property
    def basic_blocks(self):
        return self.blocks

    def append_basic_block(self, name=''):
        blk = Block(parent=self, name=name)
        self.blocks.append(blk)
        return blk

    def insert_basic_block(self, before, name=''):
        """Insert block before
        """
        blk = Block(parent=self, name=name)
        self.blocks.insert(before, blk)
        return blk

    def descr_prototype(self, buf):
        """
        Describe the prototype ("head") of the function.
        """
        state = "define" if self.blocks else "declare"
        ret = self.return_value
        args = ", ".join(str(a) for a in self.args)
        name = self.get_reference()
        attrs = ' ' + ' '.join(self.attributes._to_list(
            self.ftype.return_type)) if self.attributes else ''
        if any(self.args):
            vararg = ', ...' if self.ftype.var_arg else ''
        else:
            vararg = '...' if self.ftype.var_arg else ''
        linkage = self.linkage
        cconv = self.calling_convention
        prefix = " ".join(str(x) for x in [state, linkage, cconv, ret] if x)
        metadata = self._stringify_metadata()
        metadata = ' {}'.format(metadata) if metadata else ''
        section = ' section "{}"'.format(self.section) if self.section else ''
        pt_str = "{prefix} {name}({args}{vararg}){attrs}{section}{metadata}\n"
        prototype = pt_str.format(prefix=prefix, name=name, args=args,
                                  vararg=vararg, attrs=attrs, section=section,
                                  metadata=metadata)
        buf.append(prototype)

    def descr_body(self, buf):
        """
        Describe of the body of the function.
        """
        for blk in self.blocks:
            blk.descr(buf)

    def descr(self, buf):
        self.descr_prototype(buf)
        if self.blocks:
            buf.append("{\n")
            self.descr_body(buf)
            buf.append("}\n")

    def __str__(self):
        buf = []
        self.descr(buf)
        return "".join(buf)

    @property
    def is_declaration(self):
        return len(self.blocks) == 0


class ArgumentAttributes(AttributeSet):
    # List from
    # https://releases.llvm.org/14.0.0/docs/LangRef.html#parameter-attributes
    _known = MappingProxyType({
        # True (emit type),
        # False (emit name only)
        'byref': True,
        'byval': True,
        'elementtype': True,
        'immarg': False,
        'inalloca': True,
        'inreg': False,
        'nest': False,
        'noalias': False,
        'nocapture': False,
        'nofree': False,
        'nonnull': False,
        'noundef': False,
        'preallocated': True,
        'returned': False,
        'signext': False,
        'sret': True,
        'swiftasync': False,
        'swifterror': False,
        'swiftself': False,
        'zeroext': False,
    })

    def __init__(self, args=()):
        self._align = 0
        self._dereferenceable = 0
        self._dereferenceable_or_null = 0
        super(ArgumentAttributes, self).__init__(args)

    def _expand(self, name, typ):
        requires_type = self._known.get(name)
        if requires_type:
            return f"{name}({typ.pointee})"
        else:
            return name

    @property
    def align(self):
        return self._align

    @align.setter
    def align(self, val):
        assert isinstance(val, int) and val >= 0
        self._align = val

    @property
    def dereferenceable(self):
        return self._dereferenceable

    @dereferenceable.setter
    def dereferenceable(self, val):
        assert isinstance(val, int) and val >= 0
        self._dereferenceable = val

    @property
    def dereferenceable_or_null(self):
        return self._dereferenceable_or_null

    @dereferenceable_or_null.setter
    def dereferenceable_or_null(self, val):
        assert isinstance(val, int) and val >= 0
        self._dereferenceable_or_null = val

    def _to_list(self, typ):
        attrs = super()._to_list(typ)
        if self.align:
            attrs.append('align {0:d}'.format(self.align))
        if self.dereferenceable:
            attrs.append('dereferenceable({0:d})'.format(self.dereferenceable))
        if self.dereferenceable_or_null:
            dref = 'dereferenceable_or_null({0:d})'
            attrs.append(dref.format(self.dereferenceable_or_null))
        return attrs


class _BaseArgument(NamedValue):
    def __init__(self, parent, typ, name=''):
        assert isinstance(typ, types.Type)
        super(_BaseArgument, self).__init__(parent, typ, name=name)
        self.parent = parent
        self.attributes = ArgumentAttributes()

    def __repr__(self):
        return "<ir.%s %r of type %s>" % (self.__class__.__name__, self.name,
                                          self.type)

    def add_attribute(self, attr):
        self.attributes.add(attr)


class Argument(_BaseArgument):
    """
    The specification of a function argument.
    """

    def __str__(self):
        attrs = self.attributes._to_list(self.type)
        if attrs:
            return "{0} {1} {2}".format(self.type, ' '.join(attrs),
                                        self.get_reference())
        else:
            return "{0} {1}".format(self.type, self.get_reference())


class ReturnValue(_BaseArgument):
    """
    The specification of a function's return value.
    """

    def __str__(self):
        attrs = self.attributes._to_list(self.type)
        if attrs:
            return "{0} {1}".format(' '.join(attrs), self.type)
        else:
            return str(self.type)


class Block(NamedValue):
    """
    A LLVM IR basic block. A basic block is a sequence of
    instructions whose execution always goes from start to end.  That
    is, a control flow instruction (branch) can only appear as the
    last instruction, and incoming branches can only jump to the first
    instruction.
    """

    def __init__(self, parent, name=''):
        super(Block, self).__init__(parent, types.LabelType(), name=name)
        self.scope = parent.scope
        self.instructions = []
        self.terminator = None

    @property
    def is_terminated(self):
        return self.terminator is not None

    @property
    def function(self):
        return self.parent

    @property
    def module(self):
        return self.parent.module

    def descr(self, buf):
        buf.append("{0}:\n".format(self._format_name()))
        buf += ["  {0}\n".format(instr) for instr in self.instructions]

    def replace(self, old, new):
        """Replace an instruction"""
        if old.type != new.type:
            raise TypeError("new instruction has a different type")
        pos = self.instructions.index(old)
        self.instructions.remove(old)
        self.instructions.insert(pos, new)

        for bb in self.parent.basic_blocks:
            for instr in bb.instructions:
                instr.replace_usage(old, new)

    def _format_name(self):
        # Per the LLVM Language Ref on identifiers, names matching the following
        # regex do not need to be quoted: [%@][-a-zA-Z$._][-a-zA-Z$._0-9]*
        # Otherwise, the identifier must be quoted and escaped.
        name = self.name
        if not _SIMPLE_IDENTIFIER_RE.match(name):
            name = name.replace('\\', '\\5c').replace('"', '\\22')
            name = '"{0}"'.format(name)
        return name


class BlockAddress(Value):
    """
    The address of a basic block.
    """

    def __init__(self, function, basic_block):
        assert isinstance(function, Function)
        assert isinstance(basic_block, Block)
        self.type = types.IntType(8).as_pointer()
        self.function = function
        self.basic_block = basic_block

    def __str__(self):
        return '{0} {1}'.format(self.type, self.get_reference())

    def get_reference(self):
        return "blockaddress({0}, {1})".format(
            self.function.get_reference(),
            self.basic_block.get_reference())
