from collections import namedtuple
import math
from functools import reduce

import numpy as np
import operator
import warnings

from llvmlite import ir

from numba.core.imputils import (lower_builtin, lower_getattr,
                                 lower_getattr_generic, lower_cast,
                                 lower_constant, iternext_impl,
                                 call_getiter, call_iternext, impl_ret_borrowed,
                                 impl_ret_untracked, numba_typeref_ctor)
from numba.core import typing, types, utils, cgutils
from numba.core.extending import overload, intrinsic
from numba.core.typeconv import Conversion
from numba.core.errors import (TypingError, LoweringError,
                               NumbaExperimentalFeatureWarning,
                               NumbaTypeError, RequireLiteralValue,
                               NumbaPerformanceWarning)
from numba.core.typing.templates import (AbstractTemplate, infer_global,
                                         signature)
from numba.misc.special import literal_unroll
from numba.core.typing.asnumbatype import as_numba_type


@overload(operator.truth)
def ol_truth(val):
    if isinstance(val, types.Boolean):
        def impl(val):
            return val
        return impl


@lower_builtin(operator.is_not, types.Any, types.Any)
def generic_is_not(context, builder, sig, args):
    """
    Implement `x is not y` as `not (x is y)`.
    """
    is_impl = context.get_function(operator.is_, sig)
    return builder.not_(is_impl(builder, args))


@lower_builtin(operator.is_, types.Any, types.Any)
def generic_is(context, builder, sig, args):
    """
    Default implementation for `x is y`
    """
    lhs_type, rhs_type = sig.args
    # the lhs and rhs have the same type
    if lhs_type == rhs_type:
            # mutable types
            if lhs_type.mutable:
                msg = 'no default `is` implementation'
                raise LoweringError(msg)
            # immutable types
            else:
                # fallbacks to `==`
                try:
                    eq_impl = context.get_function(operator.eq, sig)
                except NotImplementedError:
                    # no `==` implemented for this type
                    return cgutils.false_bit
                else:
                    return eq_impl(builder, args)
    else:
        return cgutils.false_bit


@lower_builtin(operator.is_, types.Opaque, types.Opaque)
def opaque_is(context, builder, sig, args):
    """
    Implementation for `x is y` for Opaque types.
    """
    lhs_type, rhs_type = sig.args
    # the lhs and rhs have the same type
    if lhs_type == rhs_type:
        lhs_ptr = builder.ptrtoint(args[0], cgutils.intp_t)
        rhs_ptr = builder.ptrtoint(args[1], cgutils.intp_t)

        return builder.icmp_unsigned('==', lhs_ptr, rhs_ptr)
    else:
        return cgutils.false_bit


@lower_builtin(operator.is_, types.Boolean, types.Boolean)
def bool_is_impl(context, builder, sig, args):
    """
    Implementation for `x is y` for types derived from types.Boolean
    (e.g. BooleanLiteral), and cross-checks between literal and non-literal
    booleans, to satisfy Python's behavior preserving identity for bools.
    """
    arg1, arg2 = args
    arg1_type, arg2_type = sig.args
    _arg1 = context.cast(builder, arg1, arg1_type, types.boolean)
    _arg2 = context.cast(builder, arg2, arg2_type, types.boolean)
    eq_impl = context.get_function(
        operator.eq,
        typing.signature(types.boolean, types.boolean, types.boolean)
    )
    return eq_impl(builder, (_arg1, _arg2))


# keep types.IntegerLiteral, as otherwise there's ambiguity between this and int_eq_impl
@lower_builtin(operator.eq, types.Literal, types.Literal)
@lower_builtin(operator.eq, types.IntegerLiteral, types.IntegerLiteral)
def const_eq_impl(context, builder, sig, args):
    arg1, arg2 = sig.args
    val = 0
    if arg1.literal_value == arg2.literal_value:
        val = 1
    res = ir.Constant(ir.IntType(1), val)
    return impl_ret_untracked(context, builder, sig.return_type, res)


# keep types.IntegerLiteral, as otherwise there's ambiguity between this and int_ne_impl
@lower_builtin(operator.ne, types.Literal, types.Literal)
@lower_builtin(operator.ne, types.IntegerLiteral, types.IntegerLiteral)
def const_ne_impl(context, builder, sig, args):
    arg1, arg2 = sig.args
    val = 0
    if arg1.literal_value != arg2.literal_value:
        val = 1
    res = ir.Constant(ir.IntType(1), val)
    return impl_ret_untracked(context, builder, sig.return_type, res)


def gen_non_eq(val):
    def none_equality(a, b):
        a_none = isinstance(a, types.NoneType)
        b_none = isinstance(b, types.NoneType)
        if a_none and b_none:
            def impl(a, b):
                return val
            return impl
        elif a_none ^ b_none:
            def impl(a, b):
                return not val
            return impl
    return none_equality

overload(operator.eq)(gen_non_eq(True))
overload(operator.ne)(gen_non_eq(False))

#-------------------------------------------------------------------------------

@lower_getattr_generic(types.DeferredType)
def deferred_getattr(context, builder, typ, value, attr):
    """
    Deferred.__getattr__ => redirect to the actual type.
    """
    inner_type = typ.get()
    val = context.cast(builder, value, typ, inner_type)
    imp = context.get_getattr(inner_type, attr)
    return imp(context, builder, inner_type, val, attr)

@lower_cast(types.Any, types.DeferredType)
@lower_cast(types.Optional, types.DeferredType)
@lower_cast(types.Boolean, types.DeferredType)
def any_to_deferred(context, builder, fromty, toty, val):
    actual = context.cast(builder, val, fromty, toty.get())
    model = context.data_model_manager[toty]
    return model.set(builder, model.make_uninitialized(), actual)

@lower_cast(types.DeferredType, types.Any)
@lower_cast(types.DeferredType, types.Boolean)
@lower_cast(types.DeferredType, types.Optional)
def deferred_to_any(context, builder, fromty, toty, val):
    model = context.data_model_manager[fromty]
    val = model.get(builder, val)
    return context.cast(builder, val, fromty.get(), toty)


#------------------------------------------------------------------------------

@lower_builtin(operator.getitem, types.CPointer, types.Integer)
def getitem_cpointer(context, builder, sig, args):
    base_ptr, idx = args
    elem_ptr = builder.gep(base_ptr, [idx])
    res = builder.load(elem_ptr)
    return impl_ret_borrowed(context, builder, sig.return_type, res)


@lower_builtin(operator.setitem, types.CPointer, types.Integer, types.Any)
def setitem_cpointer(context, builder, sig, args):
    base_ptr, idx, val = args
    elem_ptr = builder.gep(base_ptr, [idx])
    builder.store(val, elem_ptr)


#-------------------------------------------------------------------------------

def do_minmax(context, builder, argtys, args, cmpop):
    assert len(argtys) == len(args), (argtys, args)
    assert len(args) > 0

    def binary_minmax(accumulator, value):
        # This is careful to reproduce Python's algorithm, e.g.
        # max(1.5, nan, 2.5) should return 2.5 (not nan or 1.5)
        accty, acc = accumulator
        vty, v = value
        ty = context.typing_context.unify_types(accty, vty)
        assert ty is not None
        acc = context.cast(builder, acc, accty, ty)
        v = context.cast(builder, v, vty, ty)
        cmpsig = typing.signature(types.boolean, ty, ty)
        ge = context.get_function(cmpop, cmpsig)
        pred = ge(builder, (v, acc))
        res = builder.select(pred, v, acc)
        return ty, res

    typvals = zip(argtys, args)
    resty, resval = reduce(binary_minmax, typvals)
    return resval


@lower_builtin(max, types.BaseTuple)
def max_iterable(context, builder, sig, args):
    argtys = list(sig.args[0])
    args = cgutils.unpack_tuple(builder, args[0])
    return do_minmax(context, builder, argtys, args, operator.gt)

@lower_builtin(max, types.VarArg(types.Any))
def max_vararg(context, builder, sig, args):
    return do_minmax(context, builder, sig.args, args, operator.gt)

@lower_builtin(min, types.BaseTuple)
def min_iterable(context, builder, sig, args):
    argtys = list(sig.args[0])
    args = cgutils.unpack_tuple(builder, args[0])
    return do_minmax(context, builder, argtys, args, operator.lt)

@lower_builtin(min, types.VarArg(types.Any))
def min_vararg(context, builder, sig, args):
    return do_minmax(context, builder, sig.args, args, operator.lt)


def _round_intrinsic(tp):
    # round() rounds half to even
    return "llvm.rint.f%d" % (tp.bitwidth,)

@lower_builtin(round, types.Float)
def round_impl_unary(context, builder, sig, args):
    fltty = sig.args[0]
    llty = context.get_value_type(fltty)
    module = builder.module
    fnty = ir.FunctionType(llty, [llty])
    fn = cgutils.get_or_insert_function(module, fnty, _round_intrinsic(fltty))
    res = builder.call(fn, args)
    # unary round() returns an int
    res = builder.fptosi(res, context.get_value_type(sig.return_type))
    return impl_ret_untracked(context, builder, sig.return_type, res)

@lower_builtin(round, types.Float, types.Integer)
def round_impl_binary(context, builder, sig, args):
    fltty = sig.args[0]
    # Allow calling the intrinsic from the Python implementation below.
    # This avoids the conversion to an int in Python 3's unary round().
    _round = types.ExternalFunction(
        _round_intrinsic(fltty), typing.signature(fltty, fltty))

    def round_ndigits(x, ndigits):
        if math.isinf(x) or math.isnan(x):
            return x

        if ndigits >= 0:
            if ndigits > 22:
                # pow1 and pow2 are each safe from overflow, but
                # pow1*pow2 ~= pow(10.0, ndigits) might overflow.
                pow1 = 10.0 ** (ndigits - 22)
                pow2 = 1e22
            else:
                pow1 = 10.0 ** ndigits
                pow2 = 1.0
            y = (x * pow1) * pow2
            if math.isinf(y):
                return x
            return (_round(y) / pow2) / pow1

        else:
            pow1 = 10.0 ** (-ndigits)
            y = x / pow1
            return _round(y) * pow1

    res = context.compile_internal(builder, round_ndigits, sig, args)
    return impl_ret_untracked(context, builder, sig.return_type, res)


#-------------------------------------------------------------------------------
# Numeric constructors

@lower_builtin(int, types.Any)
@lower_builtin(float, types.Any)
def int_impl(context, builder, sig, args):
    [ty] = sig.args
    [val] = args
    res = context.cast(builder, val, ty, sig.return_type)
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower_builtin(float, types.StringLiteral)
def float_literal_impl(context, builder, sig, args):
    [ty] = sig.args
    res = context.get_constant(sig.return_type, float(ty.literal_value))
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower_builtin(complex, types.VarArg(types.Any))
def complex_impl(context, builder, sig, args):
    complex_type = sig.return_type
    float_type = complex_type.underlying_float
    if len(sig.args) == 1:
        [argty] = sig.args
        [arg] = args
        if isinstance(argty, types.Complex):
            # Cast Complex* to Complex*
            res = context.cast(builder, arg, argty, complex_type)
            return impl_ret_untracked(context, builder, sig.return_type, res)
        else:
            real = context.cast(builder, arg, argty, float_type)
            imag = context.get_constant(float_type, 0)

    elif len(sig.args) == 2:
        [realty, imagty] = sig.args
        [real, imag] = args
        real = context.cast(builder, real, realty, float_type)
        imag = context.cast(builder, imag, imagty, float_type)

    cmplx = context.make_complex(builder, complex_type)
    cmplx.real = real
    cmplx.imag = imag
    res = cmplx._getvalue()
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower_builtin(types.NumberClass, types.Any)
def number_constructor(context, builder, sig, args):
    """
    Call a number class, e.g. np.int32(...)
    """
    if isinstance(sig.return_type, types.Array):
        # Array constructor
        dt = sig.return_type.dtype
        def foo(*arg_hack):
            return np.array(arg_hack, dtype=dt)
        res = context.compile_internal(builder, foo, sig, args)
        return impl_ret_untracked(context, builder, sig.return_type, res)
    else:
        # Scalar constructor
        [val] = args
        [valty] = sig.args
        return context.cast(builder, val, valty, sig.return_type)


#-------------------------------------------------------------------------------
# Constants

@lower_constant(types.Dummy)
def constant_dummy(context, builder, ty, pyval):
    # This handles None, etc.
    return context.get_dummy_value()

@lower_constant(types.ExternalFunctionPointer)
def constant_function_pointer(context, builder, ty, pyval):
    ptrty = context.get_function_pointer_type(ty)
    ptrval = context.add_dynamic_addr(builder, ty.get_pointer(pyval),
                                      info=str(pyval))
    return builder.bitcast(ptrval, ptrty)


@lower_constant(types.Optional)
def constant_optional(context, builder, ty, pyval):
    if pyval is None:
        return context.make_optional_none(builder, ty.type)
    else:
        return context.make_optional_value(builder, ty.type, pyval)


# -----------------------------------------------------------------------------

@lower_builtin(type, types.Any)
def type_impl(context, builder, sig, args):
    """
    One-argument type() builtin.
    """
    return context.get_dummy_value()


@lower_builtin(iter, types.IterableType)
def iter_impl(context, builder, sig, args):
    ty, = sig.args
    val, = args
    iterval = call_getiter(context, builder, ty, val)
    return iterval


@lower_builtin(next, types.IteratorType)
def next_impl(context, builder, sig, args):
    iterty, = sig.args
    iterval, = args

    res = call_iternext(context, builder, iterty, iterval)

    with builder.if_then(builder.not_(res.is_valid()), likely=False):
        context.call_conv.return_user_exc(builder, StopIteration, ())

    return res.yielded_value()


# -----------------------------------------------------------------------------

@lower_builtin("not in", types.Any, types.Any)
def not_in(context, builder, sig, args):
    def in_impl(a, b):
        return operator.contains(b, a)

    res = context.compile_internal(builder, in_impl, sig, args)
    return builder.not_(res)


# -----------------------------------------------------------------------------

@lower_builtin(len, types.ConstSized)
def constsized_len(context, builder, sig, args):
    [ty] = sig.args
    retty = sig.return_type
    res = context.get_constant(retty, len(ty.types))
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower_builtin(bool, types.Sized)
def sized_bool(context, builder, sig, args):
    [ty] = sig.args
    if len(ty):
        return cgutils.true_bit
    else:
        return cgutils.false_bit

@lower_builtin(tuple)
def lower_empty_tuple(context, builder, sig, args):
    retty = sig.return_type
    res = context.get_constant_undef(retty)
    return impl_ret_untracked(context, builder, sig.return_type, res)

@lower_builtin(tuple, types.BaseTuple)
def lower_tuple(context, builder, sig, args):
    val, = args
    return impl_ret_borrowed(context, builder, sig.return_type, val)

@overload(bool)
def bool_sequence(x):
    valid_types = (
        types.CharSeq,
        types.UnicodeCharSeq,
        types.DictType,
        types.ListType,
        types.UnicodeType,
        types.Set,
    )

    if isinstance(x, valid_types):
        def bool_impl(x):
            return len(x) > 0
        return bool_impl

@overload(bool, inline='always')
def bool_none(x):
    if isinstance(x, types.NoneType) or x is None:
        return lambda x: False

# -----------------------------------------------------------------------------

def get_type_max_value(typ):
    if isinstance(typ, types.Float):
        return np.inf
    if isinstance(typ, types.Integer):
        return typ.maxval
    raise NotImplementedError("Unsupported type")

def get_type_min_value(typ):
    if isinstance(typ, types.Float):
        return -np.inf
    if isinstance(typ, types.Integer):
        return typ.minval
    raise NotImplementedError("Unsupported type")

@infer_global(get_type_min_value)
@infer_global(get_type_max_value)
class MinValInfer(AbstractTemplate):
    def generic(self, args, kws):
        assert not kws
        assert len(args) == 1
        if isinstance(args[0], (types.DType, types.NumberClass)):
            return signature(args[0].dtype, *args)

@lower_builtin(get_type_min_value, types.NumberClass)
@lower_builtin(get_type_min_value, types.DType)
def lower_get_type_min_value(context, builder, sig, args):
    typ = sig.args[0].dtype

    if isinstance(typ, types.Integer):
        bw = typ.bitwidth
        lty = ir.IntType(bw)
        val = typ.minval
        res = ir.Constant(lty, val)
    elif isinstance(typ, types.Float):
        bw = typ.bitwidth
        if bw == 32:
            lty = ir.FloatType()
        elif bw == 64:
            lty = ir.DoubleType()
        else:
            raise NotImplementedError("llvmlite only supports 32 and 64 bit floats")
        npty = getattr(np, 'float{}'.format(bw))
        res = ir.Constant(lty, -np.inf)
    elif isinstance(typ, (types.NPDatetime, types.NPTimedelta)):
        bw = 64
        lty = ir.IntType(bw)
        val = types.int64.minval + 1 # minval is NaT, so minval + 1 is the smallest value
        res = ir.Constant(lty, val)
    return impl_ret_untracked(context, builder, lty, res)

@lower_builtin(get_type_max_value, types.NumberClass)
@lower_builtin(get_type_max_value, types.DType)
def lower_get_type_max_value(context, builder, sig, args):
    typ = sig.args[0].dtype

    if isinstance(typ, types.Integer):
        bw = typ.bitwidth
        lty = ir.IntType(bw)
        val = typ.maxval
        res = ir.Constant(lty, val)
    elif isinstance(typ, types.Float):
        bw = typ.bitwidth
        if bw == 32:
            lty = ir.FloatType()
        elif bw == 64:
            lty = ir.DoubleType()
        else:
            raise NotImplementedError("llvmlite only supports 32 and 64 bit floats")
        npty = getattr(np, 'float{}'.format(bw))
        res = ir.Constant(lty, np.inf)
    elif isinstance(typ, (types.NPDatetime, types.NPTimedelta)):
        bw = 64
        lty = ir.IntType(bw)
        val = types.int64.maxval
        res = ir.Constant(lty, val)
    return impl_ret_untracked(context, builder, lty, res)

# -----------------------------------------------------------------------------

from numba.core.typing.builtins import IndexValue, IndexValueType
from numba.extending import overload, register_jitable

@lower_builtin(IndexValue, types.intp, types.Type)
@lower_builtin(IndexValue, types.uintp, types.Type)
def impl_index_value(context, builder, sig, args):
    typ = sig.return_type
    index, value = args
    index_value = cgutils.create_struct_proxy(typ)(context, builder)
    index_value.index = index
    index_value.value = value
    return index_value._getvalue()


@overload(min)
def indval_min(indval1, indval2):
    if isinstance(indval1, IndexValueType) and \
       isinstance(indval2, IndexValueType):
        def min_impl(indval1, indval2):
            if np.isnan(indval1.value):
                if np.isnan(indval2.value):
                    # both indval1 and indval2 are nans so order by index
                    if indval1.index < indval2.index:
                        return indval1
                    else:
                        return indval2
                else:
                    # comparing against one nan always considered less
                    return indval1
            elif np.isnan(indval2.value):
                # indval1 not a nan but indval2 is so consider indval2 less
                return indval2
            elif indval1.value > indval2.value:
                return indval2
            elif indval1.value == indval2.value:
                if indval1.index < indval2.index:
                    return indval1
                else:
                    return indval2
            return indval1
        return min_impl


@overload(min)
def boolval_min(val1, val2):
    if isinstance(val1, types.Boolean) and \
       isinstance(val2, types.Boolean):
        def bool_min_impl(val1, val2):
            return val1 and val2
        return bool_min_impl


@overload(max)
def indval_max(indval1, indval2):
    if isinstance(indval1, IndexValueType) and \
       isinstance(indval2, IndexValueType):
        def max_impl(indval1, indval2):
            if np.isnan(indval1.value):
                if np.isnan(indval2.value):
                    # both indval1 and indval2 are nans so order by index
                    if indval1.index < indval2.index:
                        return indval1
                    else:
                        return indval2
                else:
                    # comparing against one nan always considered larger
                    return indval1
            elif np.isnan(indval2.value):
                # indval1 not a nan but indval2 is so consider indval2 larger
                return indval2
            elif indval2.value > indval1.value:
                return indval2
            elif indval1.value == indval2.value:
                if indval1.index < indval2.index:
                    return indval1
                else:
                    return indval2
            return indval1
        return max_impl


@overload(max)
def boolval_max(val1, val2):
    if isinstance(val1, types.Boolean) and \
       isinstance(val2, types.Boolean):
        def bool_max_impl(val1, val2):
            return val1 or val2
        return bool_max_impl


greater_than = register_jitable(lambda a, b: a > b)
less_than = register_jitable(lambda a, b: a < b)


@register_jitable
def min_max_impl(iterable, op):
    if isinstance(iterable, types.IterableType):
        def impl(iterable):
            it = iter(iterable)
            return_val = next(it)
            for val in it:
                if op(val, return_val):
                    return_val = val
            return return_val
        return impl


@overload(min)
def iterable_min(iterable):
    return min_max_impl(iterable, less_than)


@overload(max)
def iterable_max(iterable):
    return min_max_impl(iterable, greater_than)


@lower_builtin(types.TypeRef, types.VarArg(types.Any))
def redirect_type_ctor(context, builder, sig, args):
    """Redirect constructor implementation to `numba_typeref_ctor(cls, *args)`,
    which should be overloaded by the type's implementation.

    For example:

        d = Dict()

    `d` will be typed as `TypeRef[DictType]()`.  Thus, it will call into this
    implementation.  We need to redirect the lowering to a function
    named ``numba_typeref_ctor``.
    """
    cls = sig.return_type

    def call_ctor(cls, *args):
        return numba_typeref_ctor(cls, *args)

    # Pack arguments into a tuple for `*args`
    ctor_args = types.Tuple.from_types(sig.args)
    # Make signature T(TypeRef[T], *args) where T is cls
    sig = typing.signature(cls, types.TypeRef(cls), ctor_args)
    if len(ctor_args) > 0:
        args = (context.get_dummy_value(),   # Type object has no runtime repr.
                context.make_tuple(builder, ctor_args, args))
    else:
        args = (context.get_dummy_value(),   # Type object has no runtime repr.
                context.make_tuple(builder, ctor_args, ()))

    return context.compile_internal(builder, call_ctor, sig, args)


@overload(sum)
def ol_sum(iterable, start=0):
    # Cpython explicitly rejects strings, bytes and bytearrays
    # https://github.com/python/cpython/blob/3.9/Python/bltinmodule.c#L2310-L2329 # noqa: E501
    error = None
    if isinstance(start, types.UnicodeType):
        error = ('strings', '')
    elif isinstance(start, types.Bytes):
        error = ('bytes', 'b')
    elif isinstance(start, types.ByteArray):
        error = ('bytearray', 'b')

    if error is not None:
        msg = "sum() can't sum {} [use {}''.join(seq) instead]".format(*error)
        raise TypingError(msg)

    # if the container is homogeneous then it's relatively easy to handle.
    if isinstance(iterable, (types.containers._HomogeneousTuple, types.List,
                             types.ListType, types.Array, types.RangeType)):
        iterator = iter
    elif isinstance(iterable, (types.containers._HeterogeneousTuple)):
        # if container is heterogeneous then literal unroll and hope for the
        # best.
        iterator = literal_unroll
    else:
        return None

    def impl(iterable, start=0):
        acc = start
        for x in iterator(iterable):
            # This most likely widens the type, this is expected Numba behaviour
            acc = acc + x
        return acc
    return impl


# ------------------------------------------------------------------------------
# map, filter, reduce


@overload(map)
def ol_map(func, iterable, *args):
    def impl(func, iterable, *args):
        for x in zip(iterable, *args):
            yield func(*x)
    return impl


@overload(filter)
def ol_filter(func, iterable):
    if (func is None) or isinstance(func, types.NoneType):
        def impl(func, iterable):
            for x in iterable:
                if x:
                    yield x
    else:
        def impl(func, iterable):
            for x in iterable:
                if func(x):
                    yield x
    return impl


@overload(isinstance)
def ol_isinstance(var, typs):

    def true_impl(var, typs):
        return True

    def false_impl(var, typs):
        return False

    var_ty = as_numba_type(var)

    if isinstance(var_ty, types.Optional):
        msg = f'isinstance cannot handle optional types. Found: "{var_ty}"'
        raise NumbaTypeError(msg)

    # NOTE: The current implementation of `isinstance` restricts the type of the
    # instance variable to types that are well known and in common use. The
    # danger of unrestricted type comparison is that a "default" of `False` is
    # required and this means that if there is a bug in the logic of the
    # comparison tree `isinstance` returns False! It's therefore safer to just
    # reject the compilation as untypable!
    supported_var_ty = (types.Number, types.Bytes, types.RangeType,
                        types.DictType, types.LiteralStrKeyDict, types.List,
                        types.ListType, types.Tuple, types.UniTuple, types.Set,
                        types.Function, types.ClassType, types.UnicodeType,
                        types.ClassInstanceType, types.NoneType, types.Array,
                        types.Boolean, types.Float, types.UnicodeCharSeq,
                        types.Complex, types.NPDatetime, types.NPTimedelta,)
    if not isinstance(var_ty, supported_var_ty):
        msg = f'isinstance() does not support variables of type "{var_ty}".'
        raise NumbaTypeError(msg)

    t_typs = typs

    # Check the types that the var can be an instance of, it'll be a scalar,
    # a unituple or a tuple.
    if isinstance(t_typs, types.UniTuple):
        # corner case - all types in isinstance are the same
        t_typs = (t_typs.key[0])

    if not isinstance(t_typs, types.Tuple):
        t_typs = (t_typs, )

    for typ in t_typs:

        if isinstance(typ, types.Function):
            key = typ.key[0]  # functions like int(..), float(..), str(..)
        elif isinstance(typ, types.ClassType):
            key = typ  # jitclasses
        else:
            key = typ.key

        # corner cases for bytes, range, ...
        # avoid registering those types on `as_numba_type`
        types_not_registered = {
            bytes: types.Bytes,
            range: types.RangeType,
            dict: (types.DictType, types.LiteralStrKeyDict),
            list: types.List,
            tuple: types.BaseTuple,
            set: types.Set,
        }
        if key in types_not_registered:
            if isinstance(var_ty, types_not_registered[key]):
                return true_impl
            continue

        if isinstance(typ, types.TypeRef):
            # Use of Numba type classes is in general not supported as they do
            # not work when the jit is disabled.
            if key not in (types.ListType, types.DictType):
                msg = ("Numba type classes (except numba.typed.* container "
                       "types) are not supported.")
                raise NumbaTypeError(msg)
            # Case for TypeRef (i.e. isinstance(var, typed.List))
            #      var_ty == ListType[int64] (instance)
            #         typ == types.ListType  (class)
            return true_impl if type(var_ty) is key else false_impl
        else:
            numba_typ = as_numba_type(key)
            if var_ty == numba_typ:
                return true_impl
            elif isinstance(numba_typ, (types.NPDatetime, types.NPTimedelta)):
                if isinstance(var_ty, type(numba_typ)):
                    return true_impl
            elif isinstance(numba_typ, types.ClassType) and \
                    isinstance(var_ty, types.ClassInstanceType) and \
                    var_ty.key == numba_typ.instance_type.key:
                # check for jitclasses
                return true_impl
            elif isinstance(numba_typ, types.Container) and \
                    numba_typ.key[0] == types.undefined:
                # check for containers (list, tuple, set, ...)
                if isinstance(var_ty, numba_typ.__class__) or \
                    (isinstance(var_ty, types.BaseTuple) and \
                        isinstance(numba_typ, types.BaseTuple)):
                    return true_impl

    return false_impl


# -- getattr implementation

def _getattr_raise_attr_exc(obj, name):
    # Dummy function for the purpose of creating an overloadable stub from
    # which to raise an AttributeError as needed
    pass


@overload(_getattr_raise_attr_exc)
def ol__getattr_raise_attr_exc(obj, name):
    if not isinstance(name, types.StringLiteral):
        raise RequireLiteralValue("argument 'name' must be a literal string")
    lname = name.literal_value
    message = f"'{obj}' has no attribute '{lname}'"

    def impl(obj, name):
        raise AttributeError(message)
    return impl


@intrinsic
def resolve_getattr(tyctx, obj, name, default):
    if not isinstance(name, types.StringLiteral):
        raise RequireLiteralValue("argument 'name' must be a literal string")
    lname = name.literal_value
    fn = tyctx.resolve_getattr(obj, lname)
    # Cannot handle things like `getattr(np, 'cos')` as the return type is
    # types.Function.
    if isinstance(fn, types.Function):
        msg = ("Returning function objects is not implemented. "
               f"getattr() was requested to return {fn} from attribute "
               f"'{lname}' of {obj}.")
        raise TypingError(msg)

    if fn is None:  # No attribute

        # if default is not _getattr_default then return the default
        if not (isinstance(default, types.NamedTuple) and
                default.instance_class == _getattr_default_type):
            # it's not the marker default value, so return it
            sig = default(obj, name, default)
            def impl(cgctx, builder, sig, llargs):
                tmp = llargs[-1]
                cgctx.nrt.incref(builder, default, tmp)
                return tmp
        else:
            # else wire in raising an AttributeError
            fnty = tyctx.resolve_value_type(_getattr_raise_attr_exc)
            raise_sig = fnty.get_call_type(tyctx, (obj, name), {})
            sig = types.none(obj, name, default)
            def impl(cgctx, builder, sig, llargs):
                native_impl = cgctx.get_function(fnty, raise_sig)
                return native_impl(builder, llargs[:-1])

    else:  # Attribute present, wire in handing it back to the overload(getattr)
        sig = fn(obj, name, default)
        if isinstance(fn, types.BoundFunction):
            # It's a method on an object
            def impl(cgctx, builder, sig, ll_args):
                cast_type = fn.this
                casted = cgctx.cast(builder, ll_args[0], obj, cast_type)
                res = cgctx.get_bound_function(builder, casted, cast_type)
                cgctx.nrt.incref(builder, fn, res)
                return res
        else:
            # Else it's some other type of attribute.
            # Ensure typing calls occur at typing time, not at lowering
            attrty = tyctx.resolve_getattr(obj, lname)
            def impl(cgctx, builder, sig, ll_args):
                attr_impl = cgctx.get_getattr(obj, lname)
                res = attr_impl(cgctx, builder, obj, ll_args[0], lname)
                casted = cgctx.cast(builder, res, attrty, fn)
                cgctx.nrt.incref(builder, fn, casted)
                return casted
    return sig, impl


# These are marker objects to indicate "no default has been provided" in a call
_getattr_default_type = namedtuple('_getattr_default_type', '')
_getattr_default = _getattr_default_type()


# getattr with no default arg, obj is an open type and name is forced as a
# literal string. The _getattr_default marker is used to indicate "no default
# was provided".
@overload(getattr, prefer_literal=True)
def ol_getattr_2(obj, name):
    def impl(obj, name):
        return resolve_getattr(obj, name, _getattr_default)
    return impl


# getattr with default arg present, obj is an open type, name is forced as a
# literal string, the "default" is again an open type. Note that the CPython
# definition is: `getattr(object, name[, default]) -> value`, the `default`
# is not a kwarg.
@overload(getattr)
def ol_getattr_3(obj, name, default):
    def impl(obj, name, default):
        return resolve_getattr(obj, name, default)
    return impl


@intrinsic
def resolve_hasattr(tyctx, obj, name):
    if not isinstance(name, types.StringLiteral):
        raise RequireLiteralValue("argument 'name' must be a literal string")
    lname = name.literal_value
    fn = tyctx.resolve_getattr(obj, lname)
    # Whilst technically the return type could be a types.bool_, the literal
    # value is resolvable at typing time. Propagating this literal information
    # into the type system allows the compiler to prune branches based on a
    # hasattr predicate. As a result the signature is based on literals. This is
    # "safe" because the overload requires a literal string so each will be a
    # different variant of (obj, literal(name)) -> literal(bool).
    if fn is None:
        retty = types.literal(False)
    else:
        retty = types.literal(True)
    sig = retty(obj, name)
    def impl(cgctx, builder, sig, ll_args):
       return cgutils.false_bit if fn is None else cgutils.true_bit
    return sig, impl

# hasattr cannot be implemented as a getattr call and then catching
# AttributeError because Numba doesn't support catching anything other than
# "Exception", so lacks the specificity required. Instead this implementation
# tries to resolve the attribute via typing information and returns True/False
# based on that.
@overload(hasattr)
def ol_hasattr(obj, name):
    def impl(obj, name):
        return resolve_hasattr(obj, name)
    return impl


@overload(repr)
def ol_repr_generic(obj):
    missing_repr_format = f"<object type:{obj}>"
    def impl(obj):
        attr = '__repr__'
        if hasattr(obj, attr) == True:
            return getattr(obj, attr)()
        else:
            # There's no __str__ or __repr__ defined for this object, return
            # something generic
            return missing_repr_format
    return impl


@overload(str)
def ol_str_generic(object=''):
    def impl(object=""):
        attr = '__str__'
        if hasattr(object, attr) == True:
            return getattr(object, attr)()
        else:
            return repr(object)
    return impl
