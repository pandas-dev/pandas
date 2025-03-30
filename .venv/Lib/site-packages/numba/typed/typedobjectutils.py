""" Common compiler level utilities for typed dict and list. """

import operator
import warnings

from llvmlite import ir

from numba.core import types, cgutils
from numba.core import typing
from numba.core.registry import cpu_target
from numba.core.typeconv import Conversion
from numba.core.extending import intrinsic
from numba.core.errors import (TypingError, NumbaTypeSafetyWarning,
                               NumbaTypeError)


def _as_bytes(builder, ptr):
    """Helper to do (void*)ptr
    """
    return builder.bitcast(ptr, cgutils.voidptr_t)


@intrinsic
def _cast(typingctx, val, typ):
    """Cast *val* to *typ*
    """
    def codegen(context, builder, signature, args):
        [val, typ] = args
        context.nrt.incref(builder, signature.return_type, val)
        return val
    # Using implicit casting in argument types
    casted = typ.instance_type
    _sentry_safe_cast(val, casted)
    sig = casted(casted, typ)
    return sig, codegen


def _sentry_safe_cast(fromty, toty):
    """Check and raise TypingError if *fromty* cannot be safely cast to *toty*
    """
    tyctxt = cpu_target.typing_context
    fromty, toty = map(types.unliteral, (fromty, toty))
    by = tyctxt.can_convert(fromty, toty)

    def warn():
        m = 'unsafe cast from {} to {}. Precision may be lost.'
        warnings.warn(m.format(fromty, toty),
                      category=NumbaTypeSafetyWarning)

    isint = lambda x: isinstance(x, types.Integer)
    isflt = lambda x: isinstance(x, types.Float)
    iscmplx = lambda x: isinstance(x, types.Complex)
    isdict = lambda x: isinstance(x, types.DictType)
    # Only check against numeric types.
    if by is None or by > Conversion.safe:
        if isint(fromty) and isint(toty):
            # Accept if both types are ints
            warn()
        elif isint(fromty) and isflt(toty):
            # Accept if ints to floats
            warn()
        elif isflt(fromty) and isflt(toty):
            # Accept if floats to floats
            warn()
        elif iscmplx(fromty) and iscmplx(toty):
            # Accept if complex to complex
            warn()
        elif isdict(fromty) and isdict(toty):
            pass # it's complaining about initial values being different
        elif not isinstance(toty, types.Number):
            # Non-numbers
            warn()
        else:
            # Make it a hard error for numeric type that changes domain.
            m = 'cannot safely cast {} to {}. Please cast explicitly.'
            raise TypingError(m.format(fromty, toty))


def _sentry_safe_cast_default(default, valty):
    """Similar to _sentry_safe_cast but handle default value.
    """
    # Handle default values
    # TODO: simplify default values; too many possible way to spell None
    if default is None:
        return
    if isinstance(default, (types.Omitted, types.NoneType)):
        return
    return _sentry_safe_cast(default, valty)


@intrinsic
def _nonoptional(typingctx, val):
    """Typing trick to cast Optional[T] to T
    """
    if not isinstance(val, types.Optional):
        raise NumbaTypeError('expected an optional')

    def codegen(context, builder, sig, args):
        context.nrt.incref(builder, sig.return_type, args[0])
        return args[0]

    casted = val.type
    sig = casted(casted)
    return sig, codegen


def _container_get_data(context, builder, container_ty, c):
    """Helper to get the C list pointer in a numba containers.
    """
    ctor = cgutils.create_struct_proxy(container_ty)
    conatainer_struct = ctor(context, builder, value=c)
    return conatainer_struct.data


def _container_get_meminfo(context, builder, container_ty, c):
    """Helper to get the meminfo for a container
    """
    ctor = cgutils.create_struct_proxy(container_ty)
    conatainer_struct = ctor(context, builder, value=c)
    return conatainer_struct.meminfo


def _get_incref_decref(context, module, datamodel, container_element_type):
    assert datamodel.contains_nrt_meminfo()

    fe_type = datamodel.fe_type
    data_ptr_ty = datamodel.get_data_type().as_pointer()
    refct_fnty = ir.FunctionType(ir.VoidType(), [data_ptr_ty])
    incref_fn = cgutils.get_or_insert_function(
        module, refct_fnty, '.numba_{}.{}_incref'.format(
            context.fndesc.mangled_name, container_element_type),)

    builder = ir.IRBuilder(incref_fn.append_basic_block())
    context.nrt.incref(
        builder, fe_type,
        datamodel.load_from_data_pointer(builder, incref_fn.args[0]),
    )
    builder.ret_void()

    decref_fn = cgutils.get_or_insert_function(
        module, refct_fnty, name='.numba_{}.{}_decref'.format(
            context.fndesc.mangled_name, container_element_type),)
    builder = ir.IRBuilder(decref_fn.append_basic_block())
    context.nrt.decref(
        builder, fe_type,
        datamodel.load_from_data_pointer(builder, decref_fn.args[0]),
    )
    builder.ret_void()

    return incref_fn, decref_fn


def _get_equal(context, module, datamodel, container_element_type):
    assert datamodel.contains_nrt_meminfo()

    fe_type = datamodel.fe_type
    data_ptr_ty = datamodel.get_data_type().as_pointer()

    wrapfnty = context.call_conv.get_function_type(types.int32,
                                                   [fe_type, fe_type])
    argtypes = [fe_type, fe_type]

    def build_wrapper(fn):
        builder = ir.IRBuilder(fn.append_basic_block())
        args = context.call_conv.decode_arguments(builder, argtypes, fn)

        sig = typing.signature(types.boolean, fe_type, fe_type)
        op = operator.eq
        fnop = context.typing_context.resolve_value_type(op)
        fnop.get_call_type(context.typing_context, sig.args, {})
        eqfn = context.get_function(fnop, sig)
        res = eqfn(builder, args)
        intres = context.cast(builder, res, types.boolean, types.int32)
        context.call_conv.return_value(builder, intres)

    wrapfn = cgutils.get_or_insert_function(
        module, wrapfnty, name='.numba_{}.{}_equal.wrap'.format(
            context.fndesc.mangled_name, container_element_type))
    build_wrapper(wrapfn)

    equal_fnty = ir.FunctionType(ir.IntType(32), [data_ptr_ty, data_ptr_ty])
    equal_fn = cgutils.get_or_insert_function(
        module, equal_fnty, name='.numba_{}.{}_equal'.format(
            context.fndesc.mangled_name, container_element_type),)
    builder = ir.IRBuilder(equal_fn.append_basic_block())
    lhs = datamodel.load_from_data_pointer(builder, equal_fn.args[0])
    rhs = datamodel.load_from_data_pointer(builder, equal_fn.args[1])

    status, retval = context.call_conv.call_function(
        builder, wrapfn, types.int32, argtypes, [lhs, rhs],
    )
    with builder.if_then(status.is_ok, likely=True):
        with builder.if_then(status.is_none):
            builder.ret(context.get_constant(types.int32, 0))
        retval = context.cast(builder, retval, types.boolean, types.int32)
        builder.ret(retval)
    # Error out
    builder.ret(context.get_constant(types.int32, -1))

    return equal_fn
