from functools import partial
from collections import deque

from llvmlite import ir

from numba.core.datamodel.registry import register_default
from numba.core import types, cgutils
from numba.np import numpy_support


class DataModel(object):
    """
    DataModel describe how a FE type is represented in the LLVM IR at
    different contexts.

    Contexts are:

    - value: representation inside function body.  Maybe stored in stack.
    The representation here are flexible.

    - data: representation used when storing into containers (e.g. arrays).

    - argument: representation used for function argument.  All composite
    types are unflattened into multiple primitive types.

    - return: representation used for return argument.

    Throughput the compiler pipeline, a LLVM value is usually passed around
    in the "value" representation.  All "as_" prefix function converts from
    "value" representation.  All "from_" prefix function converts to the
    "value"  representation.

    """
    def __init__(self, dmm, fe_type):
        self._dmm = dmm
        self._fe_type = fe_type

    @property
    def fe_type(self):
        return self._fe_type

    def get_value_type(self):
        raise NotImplementedError(self)

    def get_data_type(self):
        return self.get_value_type()

    def get_argument_type(self):
        """Return a LLVM type or nested tuple of LLVM type
        """
        return self.get_value_type()

    def get_return_type(self):
        return self.get_value_type()

    def as_data(self, builder, value):
        raise NotImplementedError(self)

    def as_argument(self, builder, value):
        """
        Takes one LLVM value
        Return a LLVM value or nested tuple of LLVM value
        """
        raise NotImplementedError(self)

    def as_return(self, builder, value):
        raise NotImplementedError(self)

    def from_data(self, builder, value):
        raise NotImplementedError(self)

    def from_argument(self, builder, value):
        """
        Takes a LLVM value or nested tuple of LLVM value
        Returns one LLVM value
        """
        raise NotImplementedError(self)

    def from_return(self, builder, value):
        raise NotImplementedError(self)

    def load_from_data_pointer(self, builder, ptr, align=None):
        """
        Load value from a pointer to data.
        This is the default implementation, sufficient for most purposes.
        """
        return self.from_data(builder, builder.load(ptr, align=align))

    def traverse(self, builder):
        """
        Traverse contained members.
        Returns a iterable of contained (types, getters).
        Each getter is a one-argument function accepting a LLVM value.
        """
        return []

    def traverse_models(self):
        """
        Recursively list all models involved in this model.
        """
        return [self._dmm[t] for t in self.traverse_types()]

    def traverse_types(self):
        """
        Recursively list all frontend types involved in this model.
        """
        types = [self._fe_type]
        queue = deque([self])
        while len(queue) > 0:
            dm = queue.popleft()

            for i_dm in dm.inner_models():
                if i_dm._fe_type not in types:
                    queue.append(i_dm)
                    types.append(i_dm._fe_type)

        return types

    def inner_models(self):
        """
        List all *inner* models.
        """
        return []

    def get_nrt_meminfo(self, builder, value):
        """
        Returns the MemInfo object or None if it is not tracked.
        It is only defined for types.meminfo_pointer
        """
        return None

    def has_nrt_meminfo(self):
        return False

    def contains_nrt_meminfo(self):
        """
        Recursively check all contained types for need for NRT meminfo.
        """
        return any(model.has_nrt_meminfo() for model in self.traverse_models())

    def _compared_fields(self):
        return (type(self), self._fe_type)

    def __hash__(self):
        return hash(tuple(self._compared_fields()))

    def __eq__(self, other):
        if type(self) is type(other):
            return self._compared_fields() == other._compared_fields()
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)


@register_default(types.Omitted)
class OmittedArgDataModel(DataModel):
    """
    A data model for omitted arguments.  Only the "argument" representation
    is defined, other representations raise a NotImplementedError.
    """
    # Omitted arguments are using a dummy value type
    def get_value_type(self):
        return ir.LiteralStructType([])

    # Omitted arguments don't produce any LLVM function argument.
    def get_argument_type(self):
        return ()

    def as_argument(self, builder, val):
        return ()

    def from_argument(self, builder, val):
        assert val == (), val
        return None


@register_default(types.Boolean)
@register_default(types.BooleanLiteral)
class BooleanModel(DataModel):
    _bit_type = ir.IntType(1)
    _byte_type = ir.IntType(8)

    def get_value_type(self):
        return self._bit_type

    def get_data_type(self):
        return self._byte_type

    def get_return_type(self):
        return self.get_data_type()

    def get_argument_type(self):
        return self.get_data_type()

    def as_data(self, builder, value):
        return builder.zext(value, self.get_data_type())

    def as_argument(self, builder, value):
        return self.as_data(builder, value)

    def as_return(self, builder, value):
        return self.as_data(builder, value)

    def from_data(self, builder, value):
        ty = self.get_value_type()
        resalloca = cgutils.alloca_once(builder, ty)
        cond = builder.icmp_unsigned('==', value, value.type(0))
        with builder.if_else(cond) as (then, otherwise):
            with then:
                builder.store(ty(0), resalloca)
            with otherwise:
                builder.store(ty(1), resalloca)
        return builder.load(resalloca)

    def from_argument(self, builder, value):
        return self.from_data(builder, value)

    def from_return(self, builder, value):
        return self.from_data(builder, value)


class PrimitiveModel(DataModel):
    """A primitive type can be represented natively in the target in all
    usage contexts.
    """

    def __init__(self, dmm, fe_type, be_type):
        super(PrimitiveModel, self).__init__(dmm, fe_type)
        self.be_type = be_type

    def get_value_type(self):
        return self.be_type

    def as_data(self, builder, value):
        return value

    def as_argument(self, builder, value):
        return value

    def as_return(self, builder, value):
        return value

    def from_data(self, builder, value):
        return value

    def from_argument(self, builder, value):
        return value

    def from_return(self, builder, value):
        return value


class ProxyModel(DataModel):
    """
    Helper class for models which delegate to another model.
    """

    def get_value_type(self):
        return self._proxied_model.get_value_type()

    def get_data_type(self):
        return self._proxied_model.get_data_type()

    def get_return_type(self):
        return self._proxied_model.get_return_type()

    def get_argument_type(self):
        return self._proxied_model.get_argument_type()

    def as_data(self, builder, value):
        return self._proxied_model.as_data(builder, value)

    def as_argument(self, builder, value):
        return self._proxied_model.as_argument(builder, value)

    def as_return(self, builder, value):
        return self._proxied_model.as_return(builder, value)

    def from_data(self, builder, value):
        return self._proxied_model.from_data(builder, value)

    def from_argument(self, builder, value):
        return self._proxied_model.from_argument(builder, value)

    def from_return(self, builder, value):
        return self._proxied_model.from_return(builder, value)


@register_default(types.EnumMember)
@register_default(types.IntEnumMember)
class EnumModel(ProxyModel):
    """
    Enum members are represented exactly like their values.
    """
    def __init__(self, dmm, fe_type):
        super(EnumModel, self).__init__(dmm, fe_type)
        self._proxied_model = dmm.lookup(fe_type.dtype)


@register_default(types.Opaque)
@register_default(types.PyObject)
@register_default(types.RawPointer)
@register_default(types.NoneType)
@register_default(types.StringLiteral)
@register_default(types.EllipsisType)
@register_default(types.Function)
@register_default(types.Type)
@register_default(types.Object)
@register_default(types.Module)
@register_default(types.Phantom)
@register_default(types.UndefVar)
@register_default(types.ContextManager)
@register_default(types.Dispatcher)
@register_default(types.ObjModeDispatcher)
@register_default(types.ExceptionClass)
@register_default(types.Dummy)
@register_default(types.ExceptionInstance)
@register_default(types.ExternalFunction)
@register_default(types.EnumClass)
@register_default(types.IntEnumClass)
@register_default(types.NumberClass)
@register_default(types.TypeRef)
@register_default(types.NamedTupleClass)
@register_default(types.DType)
@register_default(types.RecursiveCall)
@register_default(types.MakeFunctionLiteral)
@register_default(types.Poison)
class OpaqueModel(PrimitiveModel):
    """
    Passed as opaque pointers
    """
    _ptr_type = ir.IntType(8).as_pointer()

    def __init__(self, dmm, fe_type):
        be_type = self._ptr_type
        super(OpaqueModel, self).__init__(dmm, fe_type, be_type)


@register_default(types.MemInfoPointer)
class MemInfoModel(OpaqueModel):

    def inner_models(self):
        return [self._dmm.lookup(self._fe_type.dtype)]

    def has_nrt_meminfo(self):
        return True

    def get_nrt_meminfo(self, builder, value):
        return value


@register_default(types.Integer)
@register_default(types.IntegerLiteral)
class IntegerModel(PrimitiveModel):
    def __init__(self, dmm, fe_type):
        be_type = ir.IntType(fe_type.bitwidth)
        super(IntegerModel, self).__init__(dmm, fe_type, be_type)


@register_default(types.Float)
class FloatModel(PrimitiveModel):
    def __init__(self, dmm, fe_type):
        if fe_type == types.float32:
            be_type = ir.FloatType()
        elif fe_type == types.float64:
            be_type = ir.DoubleType()
        else:
            raise NotImplementedError(fe_type)
        super(FloatModel, self).__init__(dmm, fe_type, be_type)


@register_default(types.CPointer)
class PointerModel(PrimitiveModel):
    def __init__(self, dmm, fe_type):
        self._pointee_model = dmm.lookup(fe_type.dtype)
        self._pointee_be_type = self._pointee_model.get_data_type()
        be_type = self._pointee_be_type.as_pointer()
        super(PointerModel, self).__init__(dmm, fe_type, be_type)


@register_default(types.EphemeralPointer)
class EphemeralPointerModel(PointerModel):

    def get_data_type(self):
        return self._pointee_be_type

    def as_data(self, builder, value):
        value = builder.load(value)
        return self._pointee_model.as_data(builder, value)

    def from_data(self, builder, value):
        raise NotImplementedError("use load_from_data_pointer() instead")

    def load_from_data_pointer(self, builder, ptr, align=None):
        return builder.bitcast(ptr, self.get_value_type())


@register_default(types.EphemeralArray)
class EphemeralArrayModel(PointerModel):

    def __init__(self, dmm, fe_type):
        super(EphemeralArrayModel, self).__init__(dmm, fe_type)
        self._data_type = ir.ArrayType(self._pointee_be_type,
                                       self._fe_type.count)

    def get_data_type(self):
        return self._data_type

    def as_data(self, builder, value):
        values = [builder.load(cgutils.gep_inbounds(builder, value, i))
                  for i in range(self._fe_type.count)]
        return cgutils.pack_array(builder, values)

    def from_data(self, builder, value):
        raise NotImplementedError("use load_from_data_pointer() instead")

    def load_from_data_pointer(self, builder, ptr, align=None):
        return builder.bitcast(ptr, self.get_value_type())


@register_default(types.ExternalFunctionPointer)
class ExternalFuncPointerModel(PrimitiveModel):
    def __init__(self, dmm, fe_type):
        sig = fe_type.sig
        # Since the function is non-Numba, there is no adaptation
        # of arguments and return value, hence get_value_type().
        retty = dmm.lookup(sig.return_type).get_value_type()
        args = [dmm.lookup(t).get_value_type() for t in sig.args]
        be_type = ir.PointerType(ir.FunctionType(retty, args))
        super(ExternalFuncPointerModel, self).__init__(dmm, fe_type, be_type)


@register_default(types.UniTuple)
@register_default(types.NamedUniTuple)
@register_default(types.StarArgUniTuple)
class UniTupleModel(DataModel):
    def __init__(self, dmm, fe_type):
        super(UniTupleModel, self).__init__(dmm, fe_type)
        self._elem_model = dmm.lookup(fe_type.dtype)
        self._count = len(fe_type)
        self._value_type = ir.ArrayType(self._elem_model.get_value_type(),
                                        self._count)
        self._data_type = ir.ArrayType(self._elem_model.get_data_type(),
                                       self._count)

    def get_value_type(self):
        return self._value_type

    def get_data_type(self):
        return self._data_type

    def get_return_type(self):
        return self.get_value_type()

    def get_argument_type(self):
        return (self._elem_model.get_argument_type(),) * self._count

    def as_argument(self, builder, value):
        out = []
        for i in range(self._count):
            v = builder.extract_value(value, [i])
            v = self._elem_model.as_argument(builder, v)
            out.append(v)
        return out

    def from_argument(self, builder, value):
        out = ir.Constant(self.get_value_type(), ir.Undefined)
        for i, v in enumerate(value):
            v = self._elem_model.from_argument(builder, v)
            out = builder.insert_value(out, v, [i])
        return out

    def as_data(self, builder, value):
        out = ir.Constant(self.get_data_type(), ir.Undefined)
        for i in range(self._count):
            val = builder.extract_value(value, [i])
            dval = self._elem_model.as_data(builder, val)
            out = builder.insert_value(out, dval, [i])
        return out

    def from_data(self, builder, value):
        out = ir.Constant(self.get_value_type(), ir.Undefined)
        for i in range(self._count):
            val = builder.extract_value(value, [i])
            dval = self._elem_model.from_data(builder, val)
            out = builder.insert_value(out, dval, [i])
        return out

    def as_return(self, builder, value):
        return value

    def from_return(self, builder, value):
        return value

    def traverse(self, builder):
        def getter(i, value):
            return builder.extract_value(value, i)
        return [(self._fe_type.dtype, partial(getter, i))
                for i in range(self._count)]

    def inner_models(self):
        return [self._elem_model]


class CompositeModel(DataModel):
    """Any model that is composed of multiple other models should subclass from
    this.
    """
    pass


class StructModel(CompositeModel):
    _value_type = None
    _data_type = None

    def __init__(self, dmm, fe_type, members):
        super(StructModel, self).__init__(dmm, fe_type)
        if members:
            self._fields, self._members = zip(*members)
        else:
            self._fields = self._members = ()
        self._models = tuple([self._dmm.lookup(t) for t in self._members])

    def get_member_fe_type(self, name):
        """
        StructModel-specific: get the Numba type of the field named *name*.
        """
        pos = self.get_field_position(name)
        return self._members[pos]

    def get_value_type(self):
        if self._value_type is None:
            self._value_type = ir.LiteralStructType([t.get_value_type()
                                                    for t in self._models])
        return self._value_type

    def get_data_type(self):
        if self._data_type is None:
            self._data_type = ir.LiteralStructType([t.get_data_type()
                                                    for t in self._models])
        return self._data_type

    def get_argument_type(self):
        return tuple([t.get_argument_type() for t in self._models])

    def get_return_type(self):
        return self.get_data_type()

    def _as(self, methname, builder, value):
        extracted = []
        for i, dm in enumerate(self._models):
            extracted.append(getattr(dm, methname)(builder,
                                                   self.get(builder, value, i)))
        return tuple(extracted)

    def _from(self, methname, builder, value):
        struct = ir.Constant(self.get_value_type(), ir.Undefined)

        for i, (dm, val) in enumerate(zip(self._models, value)):
            v = getattr(dm, methname)(builder, val)
            struct = self.set(builder, struct, v, i)

        return struct

    def as_data(self, builder, value):
        """
        Converts the LLVM struct in `value` into a representation suited for
        storing into arrays.

        Note
        ----
        Current implementation rarely changes how types are represented for
        "value" and "data".  This is usually a pointless rebuild of the
        immutable LLVM struct value.  Luckily, LLVM optimization removes all
        redundancy.

        Sample usecase: Structures nested with pointers to other structures
        that can be serialized into  a flat representation when storing into
        array.
        """
        elems = self._as("as_data", builder, value)
        struct = ir.Constant(self.get_data_type(), ir.Undefined)
        for i, el in enumerate(elems):
            struct = builder.insert_value(struct, el, [i])
        return struct

    def from_data(self, builder, value):
        """
        Convert from "data" representation back into "value" representation.
        Usually invoked when loading from array.

        See notes in `as_data()`
        """
        vals = [builder.extract_value(value, [i])
                for i in range(len(self._members))]
        return self._from("from_data", builder, vals)

    def load_from_data_pointer(self, builder, ptr, align=None):
        values = []
        for i, model in enumerate(self._models):
            elem_ptr = cgutils.gep_inbounds(builder, ptr, 0, i)
            val = model.load_from_data_pointer(builder, elem_ptr, align)
            values.append(val)

        struct = ir.Constant(self.get_value_type(), ir.Undefined)
        for i, val in enumerate(values):
            struct = self.set(builder, struct, val, i)
        return struct

    def as_argument(self, builder, value):
        return self._as("as_argument", builder, value)

    def from_argument(self, builder, value):
        return self._from("from_argument", builder, value)

    def as_return(self, builder, value):
        elems = self._as("as_data", builder, value)
        struct = ir.Constant(self.get_data_type(), ir.Undefined)
        for i, el in enumerate(elems):
            struct = builder.insert_value(struct, el, [i])
        return struct

    def from_return(self, builder, value):
        vals = [builder.extract_value(value, [i])
                for i in range(len(self._members))]
        return self._from("from_data", builder, vals)

    def get(self, builder, val, pos):
        """Get a field at the given position or the fieldname

        Args
        ----
        builder:
            LLVM IRBuilder
        val:
            value to be inserted
        pos: int or str
            field index or field name

        Returns
        -------
        Extracted value
        """
        if isinstance(pos, str):
            pos = self.get_field_position(pos)
        return builder.extract_value(val, [pos],
                                     name="extracted." + self._fields[pos])

    def set(self, builder, stval, val, pos):
        """Set a field at the given position or the fieldname

        Args
        ----
        builder:
            LLVM IRBuilder
        stval:
            LLVM struct value
        val:
            value to be inserted
        pos: int or str
            field index or field name

        Returns
        -------
        A new LLVM struct with the value inserted
        """
        if isinstance(pos, str):
            pos = self.get_field_position(pos)
        return builder.insert_value(stval, val, [pos],
                                    name="inserted." + self._fields[pos])

    def get_field_position(self, field):
        try:
            return self._fields.index(field)
        except ValueError:
            raise KeyError("%s does not have a field named %r"
                           % (self.__class__.__name__, field))

    @property
    def field_count(self):
        return len(self._fields)

    def get_type(self, pos):
        """Get the frontend type (numba type) of a field given the position
         or the fieldname

        Args
        ----
        pos: int or str
            field index or field name
        """
        if isinstance(pos, str):
            pos = self.get_field_position(pos)
        return self._members[pos]

    def get_model(self, pos):
        """
        Get the datamodel of a field given the position or the fieldname.

        Args
        ----
        pos: int or str
            field index or field name
        """
        return self._models[pos]

    def traverse(self, builder):
        def getter(k, value):
            if value.type != self.get_value_type():
                args = self.get_value_type(), value.type
                raise TypeError("expecting {0} but got {1}".format(*args))
            return self.get(builder, value, k)

        return [(self.get_type(k), partial(getter, k)) for k in self._fields]

    def inner_models(self):
        return self._models


@register_default(types.Complex)
class ComplexModel(StructModel):
    _element_type = NotImplemented

    def __init__(self, dmm, fe_type):
        members = [
            ('real', fe_type.underlying_float),
            ('imag', fe_type.underlying_float),
        ]
        super(ComplexModel, self).__init__(dmm, fe_type, members)


@register_default(types.LiteralList)
@register_default(types.LiteralStrKeyDict)
@register_default(types.Tuple)
@register_default(types.NamedTuple)
@register_default(types.StarArgTuple)
class TupleModel(StructModel):
    def __init__(self, dmm, fe_type):
        members = [('f' + str(i), t) for i, t in enumerate(fe_type)]
        super(TupleModel, self).__init__(dmm, fe_type, members)


@register_default(types.UnionType)
class UnionModel(StructModel):
    def __init__(self, dmm, fe_type):
        members = [
            ('tag', types.uintp),
            # XXX: it should really be a MemInfoPointer(types.voidptr)
            ('payload', types.Tuple.from_types(fe_type.types)),
        ]
        super(UnionModel, self).__init__(dmm, fe_type, members)



@register_default(types.Pair)
class PairModel(StructModel):
    def __init__(self, dmm, fe_type):
        members = [('first', fe_type.first_type),
                   ('second', fe_type.second_type)]
        super(PairModel, self).__init__(dmm, fe_type, members)


@register_default(types.ListPayload)
class ListPayloadModel(StructModel):
    def __init__(self, dmm, fe_type):
        # The fields are mutable but the payload is always manipulated
        # by reference.  This scheme allows mutations of an array to
        # be seen by its iterators.
        members = [
            ('size', types.intp),
            ('allocated', types.intp),
            # This member is only used only for reflected lists
            ('dirty', types.boolean),
            # Actually an inlined var-sized array
            ('data', fe_type.container.dtype),
        ]
        super(ListPayloadModel, self).__init__(dmm, fe_type, members)


@register_default(types.List)
class ListModel(StructModel):
    def __init__(self, dmm, fe_type):
        payload_type = types.ListPayload(fe_type)
        members = [
            # The meminfo data points to a ListPayload
            ('meminfo', types.MemInfoPointer(payload_type)),
            # This member is only used only for reflected lists
            ('parent', types.pyobject),
        ]
        super(ListModel, self).__init__(dmm, fe_type, members)


@register_default(types.ListIter)
class ListIterModel(StructModel):
    def __init__(self, dmm, fe_type):
        payload_type = types.ListPayload(fe_type.container)
        members = [
            # The meminfo data points to a ListPayload (shared with the
            # original list object)
            ('meminfo', types.MemInfoPointer(payload_type)),
            ('index', types.EphemeralPointer(types.intp)),
            ]
        super(ListIterModel, self).__init__(dmm, fe_type, members)


@register_default(types.SetEntry)
class SetEntryModel(StructModel):
    def __init__(self, dmm, fe_type):
        dtype = fe_type.set_type.dtype
        members = [
            # -1 = empty, -2 = deleted
            ('hash', types.intp),
            ('key', dtype),
        ]
        super(SetEntryModel, self).__init__(dmm, fe_type, members)


@register_default(types.SetPayload)
class SetPayloadModel(StructModel):
    def __init__(self, dmm, fe_type):
        entry_type = types.SetEntry(fe_type.container)
        members = [
            # Number of active + deleted entries
            ('fill', types.intp),
            # Number of active entries
            ('used', types.intp),
            # Allocated size - 1 (size being a power of 2)
            ('mask', types.intp),
            # Search finger
            ('finger', types.intp),
            # This member is only used only for reflected sets
            ('dirty', types.boolean),
            # Actually an inlined var-sized array
            ('entries', entry_type),
        ]
        super(SetPayloadModel, self).__init__(dmm, fe_type, members)

@register_default(types.Set)
class SetModel(StructModel):
    def __init__(self, dmm, fe_type):
        payload_type = types.SetPayload(fe_type)
        members = [
            # The meminfo data points to a SetPayload
            ('meminfo', types.MemInfoPointer(payload_type)),
            # This member is only used only for reflected sets
            ('parent', types.pyobject),
        ]
        super(SetModel, self).__init__(dmm, fe_type, members)

@register_default(types.SetIter)
class SetIterModel(StructModel):
    def __init__(self, dmm, fe_type):
        payload_type = types.SetPayload(fe_type.container)
        members = [
            # The meminfo data points to a SetPayload (shared with the
            # original set object)
            ('meminfo', types.MemInfoPointer(payload_type)),
            # The index into the entries table
            ('index', types.EphemeralPointer(types.intp)),
            ]
        super(SetIterModel, self).__init__(dmm, fe_type, members)


@register_default(types.Array)
@register_default(types.Buffer)
@register_default(types.ByteArray)
@register_default(types.Bytes)
@register_default(types.MemoryView)
@register_default(types.PyArray)
class ArrayModel(StructModel):
    def __init__(self, dmm, fe_type):
        ndim = fe_type.ndim
        members = [
            ('meminfo', types.MemInfoPointer(fe_type.dtype)),
            ('parent', types.pyobject),
            ('nitems', types.intp),
            ('itemsize', types.intp),
            ('data', types.CPointer(fe_type.dtype)),
            ('shape', types.UniTuple(types.intp, ndim)),
            ('strides', types.UniTuple(types.intp, ndim)),

        ]
        super(ArrayModel, self).__init__(dmm, fe_type, members)


@register_default(types.ArrayFlags)
class ArrayFlagsModel(StructModel):
    def __init__(self, dmm, fe_type):
        members = [
            ('parent', fe_type.array_type),
        ]
        super(ArrayFlagsModel, self).__init__(dmm, fe_type, members)


@register_default(types.NestedArray)
class NestedArrayModel(ArrayModel):
    def __init__(self, dmm, fe_type):
        self._be_type = dmm.lookup(fe_type.dtype).get_data_type()
        super(NestedArrayModel, self).__init__(dmm, fe_type)

    def as_storage_type(self):
        """Return the LLVM type representation for the storage of
        the nestedarray.
        """
        ret = ir.ArrayType(self._be_type, self._fe_type.nitems)
        return ret


@register_default(types.Optional)
class OptionalModel(StructModel):
    def __init__(self, dmm, fe_type):
        members = [
            ('data', fe_type.type),
            ('valid', types.boolean),
        ]
        self._value_model = dmm.lookup(fe_type.type)
        super(OptionalModel, self).__init__(dmm, fe_type, members)

    def get_return_type(self):
        return self._value_model.get_return_type()

    def as_return(self, builder, value):
        raise NotImplementedError

    def from_return(self, builder, value):
        return self._value_model.from_return(builder, value)

    def traverse(self, builder):
        def get_data(value):
            valid = get_valid(value)
            data = self.get(builder, value, "data")
            return builder.select(valid, data, ir.Constant(data.type, None))
        def get_valid(value):
            return self.get(builder, value, "valid")

        return [(self.get_type("data"), get_data),
                (self.get_type("valid"), get_valid)]


@register_default(types.Record)
class RecordModel(CompositeModel):
    def __init__(self, dmm, fe_type):
        super(RecordModel, self).__init__(dmm, fe_type)
        self._models = [self._dmm.lookup(t) for _, t in fe_type.members]
        self._be_type = ir.ArrayType(ir.IntType(8), fe_type.size)
        self._be_ptr_type = self._be_type.as_pointer()

    def get_value_type(self):
        """Passed around as reference to underlying data
        """
        return self._be_ptr_type

    def get_argument_type(self):
        return self._be_ptr_type

    def get_return_type(self):
        return self._be_ptr_type

    def get_data_type(self):
        return self._be_type

    def as_data(self, builder, value):
        return builder.load(value)

    def from_data(self, builder, value):
        raise NotImplementedError("use load_from_data_pointer() instead")

    def as_argument(self, builder, value):
        return value

    def from_argument(self, builder, value):
        return value

    def as_return(self, builder, value):
        return value

    def from_return(self, builder, value):
        return value

    def load_from_data_pointer(self, builder, ptr, align=None):
        return builder.bitcast(ptr, self.get_value_type())


@register_default(types.UnicodeCharSeq)
class UnicodeCharSeq(DataModel):
    def __init__(self, dmm, fe_type):
        super(UnicodeCharSeq, self).__init__(dmm, fe_type)
        charty = ir.IntType(numpy_support.sizeof_unicode_char * 8)
        self._be_type = ir.ArrayType(charty, fe_type.count)

    def get_value_type(self):
        return self._be_type

    def get_data_type(self):
        return self._be_type

    def as_data(self, builder, value):
        return value

    def from_data(self, builder, value):
        return value

    def as_return(self, builder, value):
        return value

    def from_return(self, builder, value):
        return value

    def as_argument(self, builder, value):
        return value

    def from_argument(self, builder, value):
        return value


@register_default(types.CharSeq)
class CharSeq(DataModel):
    def __init__(self, dmm, fe_type):
        super(CharSeq, self).__init__(dmm, fe_type)
        charty = ir.IntType(8)
        self._be_type = ir.ArrayType(charty, fe_type.count)

    def get_value_type(self):
        return self._be_type

    def get_data_type(self):
        return self._be_type

    def as_data(self, builder, value):
        return value

    def from_data(self, builder, value):
        return value

    def as_return(self, builder, value):
        return value

    def from_return(self, builder, value):
        return value

    def as_argument(self, builder, value):
        return value

    def from_argument(self, builder, value):
        return value


class CContiguousFlatIter(StructModel):
    def __init__(self, dmm, fe_type, need_indices):
        assert fe_type.array_type.layout == 'C'
        array_type = fe_type.array_type
        dtype = array_type.dtype
        ndim = array_type.ndim
        members = [('array', array_type),
                   ('stride', types.intp),
                   ('index', types.EphemeralPointer(types.intp)),
                   ]
        if need_indices:
            # For ndenumerate()
            members.append(('indices', types.EphemeralArray(types.intp, ndim)))
        super(CContiguousFlatIter, self).__init__(dmm, fe_type, members)


class FlatIter(StructModel):
    def __init__(self, dmm, fe_type):
        array_type = fe_type.array_type
        dtype = array_type.dtype
        ndim = array_type.ndim
        members = [('array', array_type),
                   ('pointers', types.EphemeralArray(types.CPointer(dtype), ndim)),
                   ('indices', types.EphemeralArray(types.intp, ndim)),
                   ('exhausted', types.EphemeralPointer(types.boolean)),
        ]
        super(FlatIter, self).__init__(dmm, fe_type, members)


@register_default(types.UniTupleIter)
class UniTupleIter(StructModel):
    def __init__(self, dmm, fe_type):
        members = [('index', types.EphemeralPointer(types.intp)),
                   ('tuple', fe_type.container,)]
        super(UniTupleIter, self).__init__(dmm, fe_type, members)


@register_default(types.misc.SliceLiteral)
@register_default(types.SliceType)
class SliceModel(StructModel):
    def __init__(self, dmm, fe_type):
        members = [('start', types.intp),
                   ('stop', types.intp),
                   ('step', types.intp),
                   ]
        super(SliceModel, self).__init__(dmm, fe_type, members)


@register_default(types.NPDatetime)
@register_default(types.NPTimedelta)
class NPDatetimeModel(PrimitiveModel):
    def __init__(self, dmm, fe_type):
        be_type = ir.IntType(64)
        super(NPDatetimeModel, self).__init__(dmm, fe_type, be_type)


@register_default(types.ArrayIterator)
class ArrayIterator(StructModel):
    def __init__(self, dmm, fe_type):
        # We use an unsigned index to avoid the cost of negative index tests.
        members = [('index', types.EphemeralPointer(types.uintp)),
                   ('array', fe_type.array_type)]
        super(ArrayIterator, self).__init__(dmm, fe_type, members)


@register_default(types.EnumerateType)
class EnumerateType(StructModel):
    def __init__(self, dmm, fe_type):
        members = [('count', types.EphemeralPointer(types.intp)),
                   ('iter', fe_type.source_type)]

        super(EnumerateType, self).__init__(dmm, fe_type, members)


@register_default(types.ZipType)
class ZipType(StructModel):
    def __init__(self, dmm, fe_type):
        members = [('iter%d' % i, source_type.iterator_type)
                   for i, source_type in enumerate(fe_type.source_types)]
        super(ZipType, self).__init__(dmm, fe_type, members)


@register_default(types.RangeIteratorType)
class RangeIteratorType(StructModel):
    def __init__(self, dmm, fe_type):
        int_type = fe_type.yield_type
        members = [('iter', types.EphemeralPointer(int_type)),
                   ('stop', int_type),
                   ('step', int_type),
                   ('count', types.EphemeralPointer(int_type))]
        super(RangeIteratorType, self).__init__(dmm, fe_type, members)


@register_default(types.Generator)
class GeneratorModel(CompositeModel):
    def __init__(self, dmm, fe_type):
        super(GeneratorModel, self).__init__(dmm, fe_type)
        # XXX Fold this in DataPacker?
        self._arg_models = [self._dmm.lookup(t) for t in fe_type.arg_types
                            if not isinstance(t, types.Omitted)]
        self._state_models = [self._dmm.lookup(t) for t in fe_type.state_types]

        self._args_be_type = ir.LiteralStructType(
            [t.get_data_type() for t in self._arg_models])
        self._state_be_type = ir.LiteralStructType(
            [t.get_data_type() for t in self._state_models])
        # The whole generator closure
        self._be_type = ir.LiteralStructType(
            [self._dmm.lookup(types.int32).get_value_type(),
             self._args_be_type, self._state_be_type])
        self._be_ptr_type = self._be_type.as_pointer()

    def get_value_type(self):
        """
        The generator closure is passed around as a reference.
        """
        return self._be_ptr_type

    def get_argument_type(self):
        return self._be_ptr_type

    def get_return_type(self):
        return self._be_type

    def get_data_type(self):
        return self._be_type

    def as_argument(self, builder, value):
        return value

    def from_argument(self, builder, value):
        return value

    def as_return(self, builder, value):
        return self.as_data(builder, value)

    def from_return(self, builder, value):
        return self.from_data(builder, value)

    def as_data(self, builder, value):
        return builder.load(value)

    def from_data(self, builder, value):
        stack = cgutils.alloca_once(builder, value.type)
        builder.store(value, stack)
        return stack


@register_default(types.ArrayCTypes)
class ArrayCTypesModel(StructModel):
    def __init__(self, dmm, fe_type):
        # ndim = fe_type.ndim
        members = [('data', types.CPointer(fe_type.dtype)),
                   ('meminfo', types.MemInfoPointer(fe_type.dtype))]
        super(ArrayCTypesModel, self).__init__(dmm, fe_type, members)


@register_default(types.RangeType)
class RangeModel(StructModel):
    def __init__(self, dmm, fe_type):
        int_type = fe_type.iterator_type.yield_type
        members = [('start', int_type),
                   ('stop', int_type),
                   ('step', int_type)]
        super(RangeModel, self).__init__(dmm, fe_type, members)


# =============================================================================

@register_default(types.NumpyNdIndexType)
class NdIndexModel(StructModel):
    def __init__(self, dmm, fe_type):
        ndim = fe_type.ndim
        members = [('shape', types.UniTuple(types.intp, ndim)),
                   ('indices', types.EphemeralArray(types.intp, ndim)),
                   ('exhausted', types.EphemeralPointer(types.boolean)),
                   ]
        super(NdIndexModel, self).__init__(dmm, fe_type, members)


@register_default(types.NumpyFlatType)
def handle_numpy_flat_type(dmm, ty):
    if ty.array_type.layout == 'C':
        return CContiguousFlatIter(dmm, ty, need_indices=False)
    else:
        return FlatIter(dmm, ty)

@register_default(types.NumpyNdEnumerateType)
def handle_numpy_ndenumerate_type(dmm, ty):
    if ty.array_type.layout == 'C':
        return CContiguousFlatIter(dmm, ty, need_indices=True)
    else:
        return FlatIter(dmm, ty)

@register_default(types.BoundFunction)
def handle_bound_function(dmm, ty):
    # The same as the underlying type
    return dmm[ty.this]


@register_default(types.NumpyNdIterType)
class NdIter(StructModel):
    def __init__(self, dmm, fe_type):
        array_types = fe_type.arrays
        ndim = fe_type.ndim
        shape_len = ndim if fe_type.need_shaped_indexing else 1
        members = [('exhausted', types.EphemeralPointer(types.boolean)),
                   ('arrays', types.Tuple(array_types)),
                   # The iterator's main shape and indices
                   ('shape', types.UniTuple(types.intp, shape_len)),
                   ('indices', types.EphemeralArray(types.intp, shape_len)),
                   ]
        # Indexing state for the various sub-iterators
        # XXX use a tuple instead?
        for i, sub in enumerate(fe_type.indexers):
            kind, start_dim, end_dim, _ = sub
            member_name = 'index%d' % i
            if kind == 'flat':
                # A single index into the flattened array
                members.append((member_name, types.EphemeralPointer(types.intp)))
            elif kind in ('scalar', 'indexed', '0d'):
                # Nothing required
                pass
            else:
                assert 0
        # Slots holding values of the scalar args
        # XXX use a tuple instead?
        for i, ty in enumerate(fe_type.arrays):
            if not isinstance(ty, types.Array):
                member_name = 'scalar%d' % i
                members.append((member_name, types.EphemeralPointer(ty)))

        super(NdIter, self).__init__(dmm, fe_type, members)


@register_default(types.DeferredType)
class DeferredStructModel(CompositeModel):
    def __init__(self, dmm, fe_type):
        super(DeferredStructModel, self).__init__(dmm, fe_type)
        self.typename = "deferred.{0}".format(id(fe_type))
        self.actual_fe_type = fe_type.get()

    def get_value_type(self):
        return ir.global_context.get_identified_type(self.typename + '.value')

    def get_data_type(self):
        return ir.global_context.get_identified_type(self.typename + '.data')

    def get_argument_type(self):
        return self._actual_model.get_argument_type()

    def as_argument(self, builder, value):
        inner = self.get(builder, value)
        return self._actual_model.as_argument(builder, inner)

    def from_argument(self, builder, value):
        res = self._actual_model.from_argument(builder, value)
        return self.set(builder, self.make_uninitialized(), res)

    def from_data(self, builder, value):
        self._define()
        elem = self.get(builder, value)
        value = self._actual_model.from_data(builder, elem)
        out = self.make_uninitialized()
        return self.set(builder, out, value)

    def as_data(self, builder, value):
        self._define()
        elem = self.get(builder, value)
        value = self._actual_model.as_data(builder, elem)
        out = self.make_uninitialized(kind='data')
        return self.set(builder, out, value)

    def from_return(self, builder, value):
        return value

    def as_return(self, builder, value):
        return value

    def get(self, builder, value):
        return builder.extract_value(value, [0])

    def set(self, builder, value, content):
        return builder.insert_value(value, content, [0])

    def make_uninitialized(self, kind='value'):
        self._define()
        if kind == 'value':
            ty = self.get_value_type()
        else:
            ty = self.get_data_type()
        return ir.Constant(ty, ir.Undefined)

    def _define(self):
        valty = self.get_value_type()
        self._define_value_type(valty)
        datty = self.get_data_type()
        self._define_data_type(datty)

    def _define_value_type(self, value_type):
        if value_type.is_opaque:
            value_type.set_body(self._actual_model.get_value_type())

    def _define_data_type(self, data_type):
        if data_type.is_opaque:
            data_type.set_body(self._actual_model.get_data_type())

    @property
    def _actual_model(self):
        return self._dmm.lookup(self.actual_fe_type)

    def traverse(self, builder):
        return [(self.actual_fe_type,
                 lambda value: builder.extract_value(value, [0]))]


@register_default(types.StructRefPayload)
class StructPayloadModel(StructModel):
    """Model for the payload of a mutable struct
    """
    def __init__(self, dmm, fe_typ):
        members = tuple(fe_typ.field_dict.items())
        super().__init__(dmm, fe_typ, members)


class StructRefModel(StructModel):
    """Model for a mutable struct.
    A reference to the payload
    """
    def __init__(self, dmm, fe_typ):
        dtype = fe_typ.get_data_type()
        members = [
            ("meminfo", types.MemInfoPointer(dtype)),
        ]
        super().__init__(dmm, fe_typ, members)

