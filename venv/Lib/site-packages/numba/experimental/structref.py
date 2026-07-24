"""Utilities for defining a mutable struct.

A mutable struct is passed by reference;
hence, structref (a reference to a struct).

"""
import operator
from numba.core.cgutils import create_struct_proxy
from numba import njit
from numba.core import types, imputils, cgutils
from numba.core.datamodel import default_manager, models
from numba.core.extending import (
    infer_getattr,
    lower_getattr_generic,
    lower_setattr_generic,
    lower_builtin,
    box,
    unbox,
    NativeValue,
    intrinsic,
    overload,
)
from numba.core.typing.templates import AttributeTemplate


class _Utils:
    """Internal builder-code utils for structref definitions.
    """
    def __init__(self, context, builder, struct_type):
        """
        Parameters
        ----------
        context :
            a numba target context
        builder :
            a llvmlite IRBuilder
        struct_type : numba.core.types.StructRef
        """
        self.context = context
        self.builder = builder
        self.struct_type = struct_type

    def new_struct_ref(self, mi):
        """Encapsulate the MemInfo from a `StructRefPayload` in a `StructRef`
        """
        context = self.context
        builder = self.builder
        struct_type = self.struct_type

        st = cgutils.create_struct_proxy(struct_type)(context, builder)
        st.meminfo = mi
        return st

    def get_struct_ref(self, val):
        """Return a helper for accessing a StructRefType
        """
        context = self.context
        builder = self.builder
        struct_type = self.struct_type

        return cgutils.create_struct_proxy(struct_type)(
            context, builder, value=val
        )

    def get_data_pointer(self, val):
        """Get the data pointer to the payload from a `StructRefType`.
        """
        context = self.context
        builder = self.builder
        struct_type = self.struct_type

        structval = self.get_struct_ref(val)
        meminfo = structval.meminfo
        data_ptr = context.nrt.meminfo_data(builder, meminfo)

        valtype = struct_type.get_data_type()
        model = context.data_model_manager[valtype]
        alloc_type = model.get_value_type()
        data_ptr = builder.bitcast(data_ptr, alloc_type.as_pointer())
        return data_ptr

    def get_data_struct(self, val):
        """Get a getter/setter helper for accessing a `StructRefPayload`
        """
        context = self.context
        builder = self.builder
        struct_type = self.struct_type

        data_ptr = self.get_data_pointer(val)
        valtype = struct_type.get_data_type()
        dataval = cgutils.create_struct_proxy(valtype)(
            context, builder, ref=data_ptr
        )
        return dataval


def define_attributes(struct_typeclass):
    """Define attributes on `struct_typeclass`.

    Defines both setters and getters in jit-code.

    This is called directly in `register()`.
    """
    @infer_getattr
    class StructAttribute(AttributeTemplate):
        key = struct_typeclass

        def generic_resolve(self, typ, attr):
            if attr in typ.field_dict:
                attrty = typ.field_dict[attr]
                return attrty

    @lower_getattr_generic(struct_typeclass)
    def struct_getattr_impl(context, builder, typ, val, attr):
        utils = _Utils(context, builder, typ)
        dataval = utils.get_data_struct(val)
        ret = getattr(dataval, attr)
        fieldtype = typ.field_dict[attr]
        return imputils.impl_ret_borrowed(context, builder, fieldtype, ret)

    @lower_setattr_generic(struct_typeclass)
    def struct_setattr_impl(context, builder, sig, args, attr):
        [inst_type, val_type] = sig.args
        [instance, val] = args
        utils = _Utils(context, builder, inst_type)
        dataval = utils.get_data_struct(instance)
        # cast val to the correct type
        field_type = inst_type.field_dict[attr]
        casted = context.cast(builder, val, val_type, field_type)
        # read old
        old_value = getattr(dataval, attr)
        # incref new value
        context.nrt.incref(builder, val_type, casted)
        # decref old value (must be last in case new value is old value)
        context.nrt.decref(builder, val_type, old_value)
        # write new
        setattr(dataval, attr, casted)


def define_boxing(struct_type, obj_class):
    """Define the boxing & unboxing logic for `struct_type` to `obj_class`.

    Defines both boxing and unboxing.

    - boxing turns an instance of `struct_type` into a PyObject of `obj_class`
    - unboxing turns an instance of `obj_class` into an instance of
      `struct_type` in jit-code.


    Use this directly instead of `define_proxy()` when the user does not
    want any constructor to be defined.
    """
    if struct_type is types.StructRef:
        raise ValueError(f"cannot register {types.StructRef}")

    obj_ctor = obj_class._numba_box_

    @box(struct_type)
    def box_struct_ref(typ, val, c):
        """
        Convert a raw pointer to a Python int.
        """
        utils = _Utils(c.context, c.builder, typ)
        struct_ref = utils.get_struct_ref(val)
        meminfo = struct_ref.meminfo

        mip_type = types.MemInfoPointer(types.voidptr)
        boxed_meminfo = c.box(mip_type, meminfo)

        ctor_pyfunc = c.pyapi.unserialize(c.pyapi.serialize_object(obj_ctor))
        ty_pyobj = c.pyapi.unserialize(c.pyapi.serialize_object(typ))

        res = c.pyapi.call_function_objargs(
            ctor_pyfunc, [ty_pyobj, boxed_meminfo],
        )
        c.pyapi.decref(ctor_pyfunc)
        c.pyapi.decref(ty_pyobj)
        c.pyapi.decref(boxed_meminfo)
        return res

    @unbox(struct_type)
    def unbox_struct_ref(typ, obj, c):
        mi_obj = c.pyapi.object_getattr_string(obj, "_meminfo")

        mip_type = types.MemInfoPointer(types.voidptr)

        mi = c.unbox(mip_type, mi_obj).value

        utils = _Utils(c.context, c.builder, typ)
        struct_ref = utils.new_struct_ref(mi)
        out = struct_ref._getvalue()

        c.pyapi.decref(mi_obj)
        return NativeValue(out)


def define_constructor(py_class, struct_typeclass, fields):
    """Define the jit-code constructor for `struct_typeclass` using the
    Python type `py_class` and the required `fields`.

    Use this instead of `define_proxy()` if the user does not want boxing
    logic defined.
    """
    # Build source code for the constructor
    params = ', '.join(fields)
    indent = ' ' * 8
    init_fields_buf = []
    for k in fields:
        init_fields_buf.append(f"st.{k} = {k}")
    init_fields = f'\n{indent}'.join(init_fields_buf)

    source = f"""
def ctor({params}):
    struct_type = struct_typeclass(list(zip({list(fields)}, [{params}])))
    def impl({params}):
        st = new(struct_type)
        {init_fields}
        return st
    return impl
"""

    glbs = dict(struct_typeclass=struct_typeclass, new=new)
    exec(source, glbs)
    ctor = glbs['ctor']
    # Make it an overload
    overload(py_class)(ctor)


def define_proxy(py_class, struct_typeclass, fields):
    """Defines a PyObject proxy for a structref.

    This makes `py_class` a valid constructor for creating a instance of
    `struct_typeclass` that contains the members as defined by `fields`.

    Parameters
    ----------
    py_class : type
        The Python class for constructing an instance of `struct_typeclass`.
    struct_typeclass : numba.core.types.Type
        The structref type class to bind to.
    fields : Sequence[str]
        A sequence of field names.

    Returns
    -------
    None
    """
    define_constructor(py_class, struct_typeclass, fields)
    define_boxing(struct_typeclass, py_class)


def register(struct_type):
    """Register a `numba.core.types.StructRef` for use in jit-code.

    This defines the data-model for lowering an instance of `struct_type`.
    This defines attributes accessor and mutator for an instance of
    `struct_type`.

    Parameters
    ----------
    struct_type : type
        A subclass of `numba.core.types.StructRef`.

    Returns
    -------
    struct_type : type
        Returns the input argument so this can act like a decorator.

    Examples
    --------

    .. code-block::

        class MyStruct(numba.core.types.StructRef):
            ...  # the simplest subclass can be empty

        numba.experimental.structref.register(MyStruct)

    """
    if struct_type is types.StructRef:
        raise ValueError(f"cannot register {types.StructRef}")
    default_manager.register(struct_type, models.StructRefModel)
    define_attributes(struct_type)
    return struct_type


@intrinsic
def new(typingctx, struct_type):
    """new(struct_type)

    A jit-code only intrinsic. Used to allocate an **empty** mutable struct.
    The fields are zero-initialized and must be set manually after calling
    the function.

    Example:

        instance = new(MyStruct)
        instance.field = field_value
    """
    from numba.experimental.jitclass.base import imp_dtor

    inst_type = struct_type.instance_type

    def codegen(context, builder, signature, args):
        # FIXME: mostly the same as jitclass ctor_impl()
        model = context.data_model_manager[inst_type.get_data_type()]
        alloc_type = model.get_value_type()
        alloc_size = context.get_abi_sizeof(alloc_type)

        meminfo = context.nrt.meminfo_alloc_dtor(
            builder,
            context.get_constant(types.uintp, alloc_size),
            imp_dtor(context, builder.module, inst_type),
        )
        data_pointer = context.nrt.meminfo_data(builder, meminfo)
        data_pointer = builder.bitcast(data_pointer, alloc_type.as_pointer())

        # Nullify all data
        builder.store(cgutils.get_null_value(alloc_type), data_pointer)

        inst_struct = context.make_helper(builder, inst_type)
        inst_struct.meminfo = meminfo

        return inst_struct._getvalue()

    sig = inst_type(struct_type)
    return sig, codegen


class StructRefProxy:
    """A PyObject proxy to the Numba allocated structref data structure.

    Notes
    -----

    * Subclasses should not define ``__init__``.
    * Subclasses can override ``__new__``.
    """
    __slots__ = ('_type', '_meminfo')

    @classmethod
    def _numba_box_(cls, ty, mi):
        """Called by boxing logic, the conversion of Numba internal
        representation into a PyObject.

        Parameters
        ----------
        ty :
            a Numba type instance.
        mi :
            a wrapped MemInfoPointer.

        Returns
        -------
        instance :
             a StructRefProxy instance.
        """
        instance = super().__new__(cls)
        instance._type = ty
        instance._meminfo = mi
        return instance

    def __new__(cls, *args):
        """Construct a new instance of the structref.

        This takes positional-arguments only due to limitation of the compiler.
        The arguments are mapped to ``cls(*args)`` in jit-code.
        """
        try:
            # use cached ctor if available
            ctor = cls.__numba_ctor
        except AttributeError:
            # lazily define the ctor
            @njit
            def ctor(*args):
                return cls(*args)
            # cache it to attribute to avoid recompilation
            cls.__numba_ctor = ctor
        return ctor(*args)

    @property
    def _numba_type_(self):
        """Returns the Numba type instance for this structref instance.

        Subclasses should NOT override.
        """
        return self._type


@lower_builtin(operator.is_, types.StructRef, types.StructRef)
def structref_is(context, builder, sig, args):
    """
    Define the 'is' operator for structrefs by comparing the memory addresses.
    This is the identity check for structref objects.
    """
    a, b = args
    aty, bty = sig.args
    a_ptr = create_struct_proxy(aty)(context, builder, value=a).meminfo
    b_ptr = create_struct_proxy(bty)(context, builder, value=b).meminfo
    return builder.icmp_unsigned("==", a_ptr, b_ptr)
