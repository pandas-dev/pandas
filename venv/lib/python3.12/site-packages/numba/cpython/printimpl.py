"""
This file implements print functionality for the CPU.
"""
from numba.core import types, typing, cgutils
from numba.core.imputils import Registry, impl_ret_untracked

registry = Registry('printimpl')
lower = registry.lower


# NOTE: the current implementation relies on CPython API even in
#       nopython mode.


@lower("print_item", types.Literal)
def print_item_impl(context, builder, sig, args):
    """
    Print a single constant value.
    """
    ty, = sig.args
    val = ty.literal_value

    pyapi = context.get_python_api(builder)

    strobj = pyapi.unserialize(pyapi.serialize_object(val))
    pyapi.print_object(strobj)
    pyapi.decref(strobj)

    res = context.get_dummy_value()
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower("print_item", types.Any)
def print_item_impl(context, builder, sig, args):
    """
    Print a single native value by boxing it in a Python object and
    invoking the Python interpreter's print routine.
    """
    ty, = sig.args
    val, = args

    pyapi = context.get_python_api(builder)
    env_manager = context.get_env_manager(builder)

    if context.enable_nrt:
        context.nrt.incref(builder, ty, val)

    obj = pyapi.from_native_value(ty, val, env_manager)
    with builder.if_else(cgutils.is_not_null(builder, obj), likely=True) as (if_ok, if_error):
        with if_ok:
            pyapi.print_object(obj)
            pyapi.decref(obj)
        with if_error:
            cstr = context.insert_const_string(builder.module,
                                               "the print() function")
            strobj = pyapi.string_from_string(cstr)
            pyapi.err_write_unraisable(strobj)
            pyapi.decref(strobj)

    res = context.get_dummy_value()
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower(print, types.VarArg(types.Any))
def print_varargs_impl(context, builder, sig, args):
    """
    A entire print() call.
    """
    pyapi = context.get_python_api(builder)
    gil = pyapi.gil_ensure()

    for i, (argtype, argval) in enumerate(zip(sig.args, args)):
        signature = typing.signature(types.none, argtype)
        imp = context.get_function("print_item", signature)
        imp(builder, [argval])
        if i < len(args) - 1:
            pyapi.print_string(' ')
    pyapi.print_string('\n')

    pyapi.gil_release(gil)
    res = context.get_dummy_value()
    return impl_ret_untracked(context, builder, sig.return_type, res)
