"""Provides Numba type, FunctionType, that makes functions as
instances of a first-class function type.
"""

from numba.extending import typeof_impl
from numba.extending import models, register_model
from numba.extending import unbox, NativeValue, box
from numba.core.imputils import lower_constant, lower_cast
from numba.core.ccallback import CFunc
from numba.core import cgutils
from llvmlite import ir
from numba.core import types
from numba.core.types import (FunctionType, UndefinedFunctionType,
                              FunctionPrototype, WrapperAddressProtocol)
from numba.core.dispatcher import Dispatcher


@typeof_impl.register(WrapperAddressProtocol)
@typeof_impl.register(CFunc)
def typeof_function_type(val, c):
    if isinstance(val, CFunc):
        sig = val._sig
    elif isinstance(val, WrapperAddressProtocol):
        sig = val.signature()
    else:
        raise NotImplementedError(
            f'function type from {type(val).__name__}')
    return FunctionType(sig)


@register_model(FunctionPrototype)
class FunctionProtoModel(models.PrimitiveModel):
    """FunctionProtoModel describes the signatures of first-class functions
    """
    def __init__(self, dmm, fe_type):
        if isinstance(fe_type, FunctionType):
            ftype = fe_type.ftype
        elif isinstance(fe_type, FunctionPrototype):
            ftype = fe_type
        else:
            raise NotImplementedError((type(fe_type)))
        retty = dmm.lookup(ftype.rtype).get_value_type()
        args = [dmm.lookup(t).get_value_type() for t in ftype.atypes]
        be_type = ir.PointerType(ir.FunctionType(retty, args))
        super(FunctionProtoModel, self).__init__(dmm, fe_type, be_type)


@register_model(FunctionType)
@register_model(UndefinedFunctionType)
class FunctionModel(models.StructModel):
    """FunctionModel holds addresses of function implementations
    """
    def __init__(self, dmm, fe_type):
        members = [
            # address of cfunc wrapper function:
            ('addr', types.voidptr),
            # address of PyObject* referencing the Python function
            # object:
            ('pyaddr', types.voidptr),
        ]
        super(FunctionModel, self).__init__(dmm, fe_type, members)


@lower_constant(types.Dispatcher)
def lower_constant_dispatcher(context, builder, typ, pyval):
    return context.add_dynamic_addr(builder, id(pyval),
                                    info=type(pyval).__name__)


@lower_constant(FunctionType)
def lower_constant_function_type(context, builder, typ, pyval):
    typ = typ.get_precise()

    if isinstance(pyval, CFunc):
        addr = pyval._wrapper_address
        sfunc = cgutils.create_struct_proxy(typ)(context, builder)
        sfunc.addr = context.add_dynamic_addr(builder, addr,
                                              info=str(typ))
        sfunc.pyaddr = context.add_dynamic_addr(builder, id(pyval),
                                                info=type(pyval).__name__)
        return sfunc._getvalue()

    if isinstance(pyval, Dispatcher):
        sfunc = cgutils.create_struct_proxy(typ)(context, builder)
        sfunc.pyaddr = context.add_dynamic_addr(builder, id(pyval),
                                                info=type(pyval).__name__)
        return sfunc._getvalue()

    if isinstance(pyval, WrapperAddressProtocol):
        addr = pyval.__wrapper_address__()
        assert typ.check_signature(pyval.signature())
        sfunc = cgutils.create_struct_proxy(typ)(context, builder)
        sfunc.addr = context.add_dynamic_addr(builder, addr,
                                              info=str(typ))
        sfunc.pyaddr = context.add_dynamic_addr(builder, id(pyval),
                                                info=type(pyval).__name__)
        return sfunc._getvalue()

    # TODO: implement support for pytypes.FunctionType, ctypes.CFUNCTYPE
    raise NotImplementedError(
        'lower_constant_struct_function_type({}, {}, {}, {})'
        .format(context, builder, typ, pyval))


def _get_wrapper_address(func, sig):
    """Return the address of a compiled cfunc wrapper function of `func`.

    Warning: The compiled function must be compatible with the given
    signature `sig`. If it is not, then result of calling the compiled
    function is undefined. The compatibility is ensured when passing
    in a first-class function to a Numba njit compiled function either
    as an argument or via namespace scoping.

    Parameters
    ----------
    func : object
      A Numba cfunc or jit decoreated function or an object that
      implements the wrapper address protocol (see note below).
    sig : Signature
      The expected function signature.

    Returns
    -------
    addr : int
      An address in memory (pointer value) of the compiled function
      corresponding to the specified signature.

    Note: wrapper address protocol
    ------------------------------

    An object implements the wrapper address protocol iff the object
    provides a callable attribute named __wrapper_address__ that takes
    a Signature instance as the argument, and returns an integer
    representing the address or pointer value of a compiled function
    for the given signature.

    """
    if not sig.is_precise():
        # addr==-1 will indicate that no implementation is available
        # for cases where type-inference did not identified the
        # function type. For example, the type of an unused
        # jit-decorated function argument will be undefined but also
        # irrelevant.
        addr = -1
    elif hasattr(func, '__wrapper_address__'):
        # func can be any object that implements the
        # __wrapper_address__ protocol.
        addr = func.__wrapper_address__()
    elif isinstance(func, CFunc):
        assert sig == func._sig
        addr = func.address
    elif isinstance(func, Dispatcher):
        cres = func.get_compile_result(sig)
        wrapper_name = cres.fndesc.llvm_cfunc_wrapper_name
        addr = cres.library.get_pointer_to_function(wrapper_name)
    else:
        raise NotImplementedError(
            f'get wrapper address of {type(func)} instance with {sig!r}')
    if not isinstance(addr, int):
        raise TypeError(
            f'wrapper address must be integer, got {type(addr)} instance')
    if addr <= 0 and addr != -1:
        raise ValueError(f'wrapper address of {type(func)} instance must be'
                         f' a positive integer but got {addr} [sig={sig}]')
    # print(f'_get_wrapper_address[{func}]({sig=}) -> {addr}')
    return addr


def lower_get_wrapper_address(context, builder, func, sig,
                              failure_mode='return_exc'):
    """Low-level call to _get_wrapper_address(func, sig).

    When calling this function, GIL must be acquired.
    """
    pyapi = context.get_python_api(builder)

    # Get the cfunc wrapper address. The code below trusts that the
    # function numba.function._get_wrapper_address exists and can be
    # called with two arguments. However, if an exception is raised in
    # the function, then it will be caught and propagated to the
    # caller.

    modname = context.insert_const_string(builder.module, __name__)
    numba_mod = pyapi.import_module_noblock(modname)
    numba_func = pyapi.object_getattr_string(
        numba_mod, '_get_wrapper_address')
    pyapi.decref(numba_mod)
    sig_obj = pyapi.unserialize(pyapi.serialize_object(sig))

    addr = pyapi.call_function_objargs(numba_func, (func, sig_obj))

    if failure_mode != 'ignore':
        with builder.if_then(cgutils.is_null(builder, addr), likely=False):
            # _get_wrapper_address has raised an exception, propagate it
            # to the caller.
            if failure_mode == 'return_exc':
                context.call_conv.return_exc(builder)
            elif failure_mode == 'return_null':
                builder.ret(pyapi.get_null_object())
            else:
                raise NotImplementedError(failure_mode)
    # else the caller will handle addr == NULL

    return addr  # new reference or NULL


@unbox(FunctionType)
def unbox_function_type(typ, obj, c):
    typ = typ.get_precise()

    sfunc = cgutils.create_struct_proxy(typ)(c.context, c.builder)

    addr = lower_get_wrapper_address(
        c.context, c.builder, obj, typ.signature, failure_mode='return_null')
    sfunc.addr = c.pyapi.long_as_voidptr(addr)
    c.pyapi.decref(addr)

    llty = c.context.get_value_type(types.voidptr)
    sfunc.pyaddr = c.builder.ptrtoint(obj, llty)

    return NativeValue(sfunc._getvalue())


@box(FunctionType)
def box_function_type(typ, val, c):
    typ = typ.get_precise()

    sfunc = cgutils.create_struct_proxy(typ)(c.context, c.builder, value=val)
    pyaddr_ptr = cgutils.alloca_once(c.builder, c.pyapi.pyobj)
    raw_ptr = c.builder.inttoptr(sfunc.pyaddr, c.pyapi.pyobj)
    with c.builder.if_then(cgutils.is_null(c.builder, raw_ptr),
                           likely=False):
        cstr = f"first-class function {typ} parent object not set"
        c.pyapi.err_set_string("PyExc_MemoryError", cstr)
        c.builder.ret(c.pyapi.get_null_object())
    c.builder.store(raw_ptr, pyaddr_ptr)
    cfunc = c.builder.load(pyaddr_ptr)
    c.pyapi.incref(cfunc)
    return cfunc


@lower_cast(UndefinedFunctionType, FunctionType)
def lower_cast_function_type_to_function_type(
        context, builder, fromty, toty, val):
    return val


@lower_cast(types.Dispatcher, FunctionType)
def lower_cast_dispatcher_to_function_type(context, builder, fromty, toty, val):
    toty = toty.get_precise()

    pyapi = context.get_python_api(builder)
    sfunc = cgutils.create_struct_proxy(toty)(context, builder)

    gil_state = pyapi.gil_ensure()
    addr = lower_get_wrapper_address(
        context, builder, val, toty.signature,
        failure_mode='return_exc')
    sfunc.addr = pyapi.long_as_voidptr(addr)
    pyapi.decref(addr)
    pyapi.gil_release(gil_state)

    llty = context.get_value_type(types.voidptr)
    sfunc.pyaddr = builder.ptrtoint(val, llty)
    return sfunc._getvalue()
