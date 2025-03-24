import itertools

import numpy as np
import operator

from numba.core import types, errors, config
from numba import prange
from numba.parfors.parfor import internal_prange

from numba.core.typing.templates import (AttributeTemplate, ConcreteTemplate,
                                         AbstractTemplate, infer_global, infer,
                                         infer_getattr, signature,
                                         bound_function, make_callable_template)


from numba.core.extending import (
    typeof_impl, type_callable, models, register_model, make_attribute_wrapper,
    )


@infer_global(print)
class Print(AbstractTemplate):
    def generic(self, args, kws):
        for a in args:
            sig = self.context.resolve_function_type("print_item", (a,), {})
            if sig is None:
                raise errors.TypingError("Type %s is not printable." % a)
            assert sig.return_type is types.none
        return signature(types.none, *args)

@infer
class PrintItem(AbstractTemplate):
    key = "print_item"

    def generic(self, args, kws):
        arg, = args
        return signature(types.none, *args)


@infer_global(abs)
class Abs(ConcreteTemplate):
    int_cases = [signature(ty, ty) for ty in sorted(types.py_signed_domain)]
    real_cases = [signature(ty, ty) for ty in sorted(types.py_real_domain)]
    complex_cases = [signature(ty.underlying_float, ty)
                     for ty in sorted(types.py_complex_domain)]
    cases = int_cases + real_cases + complex_cases


@infer_global(slice)
class Slice(ConcreteTemplate):
    cases = [
        signature(types.slice2_type, types.py_int),
        signature(types.slice2_type, types.none),
        signature(types.slice2_type, types.none, types.none),
        signature(types.slice2_type, types.none, types.py_int),
        signature(types.slice2_type, types.py_int, types.none),
        signature(types.slice2_type, types.py_int, types.py_int),
        signature(types.slice3_type, types.py_int, types.py_int, types.py_int),
        signature(types.slice3_type, types.none, types.py_int, types.py_int),
        signature(types.slice3_type, types.py_int, types.none, types.py_int),
        signature(types.slice3_type, types.py_int, types.py_int, types.none),
        signature(types.slice3_type, types.py_int, types.none, types.none),
        signature(types.slice3_type, types.none, types.py_int, types.none),
        signature(types.slice3_type, types.none, types.none, types.py_int),
        signature(types.slice3_type, types.none, types.none, types.none),
    ]


@infer_global(range, typing_key=range)
@infer_global(prange, typing_key=prange)
@infer_global(internal_prange, typing_key=internal_prange)
class Range(ConcreteTemplate):
    cases = [
        signature(types.range_state_type, types.py_int),
        signature(types.range_state_type, types.py_int, types.py_int),
        signature(types.range_state_type, types.py_int, types.py_int,
                types.py_int),
    ]


@infer
class GetIter(AbstractTemplate):
    key = "getiter"

    def generic(self, args, kws):
        assert not kws
        [obj] = args
        if isinstance(obj, types.IterableType):
            return signature(obj.iterator_type, obj)


@infer
class IterNext(AbstractTemplate):
    key = "iternext"

    def generic(self, args, kws):
        assert not kws
        [it] = args
        if isinstance(it, types.IteratorType):
            return signature(types.Pair(it.yield_type, types.py_bool), it)


@infer
class PairFirst(AbstractTemplate):
    """
    Given a heterogeneous pair, return the first element.
    """
    key = "pair_first"

    def generic(self, args, kws):
        assert not kws
        [pair] = args
        if isinstance(pair, types.Pair):
            return signature(pair.first_type, pair)


@infer
class PairSecond(AbstractTemplate):
    """
    Given a heterogeneous pair, return the second element.
    """
    key = "pair_second"

    def generic(self, args, kws):
        assert not kws
        [pair] = args
        if isinstance(pair, types.Pair):
            return signature(pair.second_type, pair)


def choose_result_bitwidth(*inputs):
    return max(tp.bitwidth for tp in inputs)

def choose_result_int(*inputs):
    """
    Choose the integer result type for an operation on integer inputs,
    according to the integer typing NBEP. In accordance with the new
    type system.
    """
    bitwidth = choose_result_bitwidth(*inputs)
    signed = any(tp.signed for tp in inputs)

    # If any integer is a NumPy integer, promotion should be to the
    # respective NumPy type.
    if any('np' in tp.name for tp in inputs):
        return types.NumPyInteger.from_bitwidth(bitwidth, signed)

    return types.py_int

all_ints = (
    sorted(set((types.py_int, types.py_int))) +
    sorted(set((types.np_int32, types.np_int64))) +
    sorted(set((types.np_uint32, types.np_uint64)))
    )
integer_binop_cases = tuple(
    signature(choose_result_int(op1, op2), op1, op2)
    for op1, op2 in itertools.product(all_ints, all_ints)
    )

class BinOp(ConcreteTemplate):
    cases = list(integer_binop_cases)

    cases += [signature(op, op, op) for op in sorted(types.py_real_domain)]
    cases += [signature(op, op, op) for op in sorted(types.py_complex_domain)]
    cases += [signature(op, op, op) for op in sorted(types.np_real_domain)]
    cases += [signature(op, op, op) for op in sorted(types.np_complex_domain)]


@infer_global(operator.add)
class BinOpAdd(BinOp):
    pass


@infer_global(operator.iadd)
class BinOpAdd(BinOp):
    pass


@infer_global(operator.sub)
class BinOpSub(BinOp):
    pass


@infer_global(operator.isub)
class BinOpSub(BinOp):
    pass


@infer_global(operator.mul)
class BinOpMul(BinOp):
    pass


@infer_global(operator.imul)
class BinOpMul(BinOp):
    pass


@infer_global(operator.mod)
class BinOpMod(ConcreteTemplate):
    cases = list(integer_binop_cases)
    cases += [signature(op, op, op) for op in sorted(types.np_real_domain)]


@infer_global(operator.imod)
class BinOpMod(ConcreteTemplate):
    cases = list(integer_binop_cases)
    cases += [signature(op, op, op) for op in sorted(types.np_real_domain)]


@infer_global(operator.truediv)
class BinOpTrueDiv(ConcreteTemplate):
    cases = [signature(types.np_float64, op1, op2)
             for op1, op2 in itertools.product(all_ints, all_ints)]
    cases += [signature(op, op, op) for op in sorted(types.np_real_domain)]
    cases += [signature(op, op, op) for op in sorted(types.np_complex_domain)]


@infer_global(operator.itruediv)
class BinOpTrueDiv(ConcreteTemplate):
    cases = [signature(types.np_float64, op1, op2)
             for op1, op2 in itertools.product(all_ints, all_ints)]
    cases += [signature(op, op, op) for op in sorted(types.np_real_domain)]
    cases += [signature(op, op, op) for op in sorted(types.np_complex_domain)]


@infer_global(operator.floordiv)
class BinOpFloorDiv(ConcreteTemplate):
    cases = list(integer_binop_cases)
    cases += [signature(op, op, op) for op in sorted(types.np_real_domain)]


@infer_global(operator.ifloordiv)
class BinOpFloorDiv(ConcreteTemplate):
    cases = list(integer_binop_cases)
    cases += [signature(op, op, op) for op in sorted(types.np_real_domain)]


@infer_global(divmod)
class DivMod(ConcreteTemplate):
    _tys = all_ints + sorted(types.np_real_domain)
    cases = [signature(types.UniTuple(ty, 2), ty, ty) for ty in _tys]


@infer_global(operator.pow)
class BinOpPower(ConcreteTemplate):
    cases = list(integer_binop_cases)
    # Ensure that float32 ** int doesn't go through DP computations
    cases += [signature(types.np_float32, types.np_float32, op)
              for op in (types.np_int32, types.np_int64, types.np_uint64)]
    cases += [signature(types.np_float64, types.np_float64, op)
              for op in (types.np_int32, types.np_int64, types.np_uint64)]
    cases += [signature(op, op, op)
              for op in sorted(types.np_real_domain)]
    cases += [signature(op, op, op)
              for op in sorted(types.np_complex_domain)]


@infer_global(operator.ipow)
class BinOpPower(ConcreteTemplate):
    cases = list(integer_binop_cases)
    # Ensure that float32 ** int doesn't go through DP computations
    cases += [signature(types.np_float32, types.np_float32, op)
              for op in (types.np_int32, types.np_int64, types.np_uint64)]
    cases += [signature(types.np_float64, types.np_float64, op)
              for op in (types.np_int32, types.np_int64, types.np_uint64)]
    cases += [signature(op, op, op)
              for op in sorted(types.np_real_domain)]
    cases += [signature(op, op, op)
              for op in sorted(types.np_complex_domain)]


@infer_global(pow)
class PowerBuiltin(BinOpPower):
    # TODO add 3 operand version
    pass


class BitwiseShiftOperation(ConcreteTemplate):
    # For bitshifts, only the first operand's signedness matters
    # to choose the operation's signedness (the second operand
    # should always be positive but will generally be considered
    # signed anyway, since it's often a constant integer).
    # (also, see issue #1995 for right-shifts)

    # The RHS type is fixed to 64-bit signed/unsigned ints.
    # The implementation will always cast the operands to the width of the
    # result type, which is the widest between the LHS type and (u)intp.
    cases = [signature(max(op, types.py_int), op, op2)
             for op in types.py_signed_domain
             for op2 in types.py_signed_domain]
    unsafe_casting = False


@infer_global(operator.lshift)
class BitwiseLeftShift(BitwiseShiftOperation):
    pass

@infer_global(operator.ilshift)
class BitwiseLeftShift(BitwiseShiftOperation):
    pass


@infer_global(operator.rshift)
class BitwiseRightShift(BitwiseShiftOperation):
    pass


@infer_global(operator.irshift)
class BitwiseRightShift(BitwiseShiftOperation):
    pass


class BitwiseLogicOperation(BinOp):
    cases = [signature(types.py_bool, types.py_bool, types.py_bool)]
    cases += [signature(types.np_bool_, types.np_bool_, types.np_bool_)]
    cases += list(integer_binop_cases)
    unsafe_casting = False


@infer_global(operator.and_)
class BitwiseAnd(BitwiseLogicOperation):
    pass


@infer_global(operator.iand)
class BitwiseAnd(BitwiseLogicOperation):
    pass


@infer_global(operator.or_)
class BitwiseOr(BitwiseLogicOperation):
    pass


@infer_global(operator.ior)
class BitwiseOr(BitwiseLogicOperation):
    pass


@infer_global(operator.xor)
class BitwiseXor(BitwiseLogicOperation):
    pass


@infer_global(operator.ixor)
class BitwiseXor(BitwiseLogicOperation):
    pass


# Bitwise invert and negate are special: we must not upcast the operand
# for unsigned numbers, as that would change the result.
# (i.e. ~np.int8(0) == 255 but ~np.int32(0) == 4294967295).

@infer_global(operator.invert)
class BitwiseInvert(ConcreteTemplate):
    # Note Numba follows the Numpy semantics of returning a bool,
    # while Python returns an int.  This makes it consistent with
    # np.invert() and makes array expressions correct.
    cases = [signature(types.py_bool, types.py_bool)]
    cases = [signature(types.np_bool_, types.np_bool_)]

    cases += [signature(choose_result_int(op), op) for op in sorted(types.np_unsigned_domain)]

    cases += [signature(choose_result_int(op), op) for op in sorted(types.py_signed_domain)]
    cases += [signature(choose_result_int(op), op) for op in sorted(types.np_signed_domain)]


    unsafe_casting = False


class UnaryOp(ConcreteTemplate):
    cases = [signature(choose_result_int(op), op) for op in sorted(types.np_unsigned_domain)]
    cases += [signature(choose_result_int(op), op) for op in sorted(types.py_signed_domain)]
    cases += [signature(choose_result_int(op), op) for op in sorted(types.np_signed_domain)]

    cases += [signature(op, op) for op in sorted(types.py_real_domain)]
    cases += [signature(op, op) for op in sorted(types.np_real_domain)]

    cases += [signature(op, op) for op in sorted(types.py_complex_domain)]
    cases += [signature(op, op) for op in sorted(types.np_complex_domain)]

    cases += [signature(types.py_int, types.py_bool)]
    cases += [signature(types.np_intp, types.np_bool_)]


@infer_global(operator.neg)
class UnaryNegate(UnaryOp):
    pass


@infer_global(operator.pos)
class UnaryPositive(UnaryOp):
   pass


@infer_global(operator.not_)
class UnaryNot(ConcreteTemplate):
    cases = [signature(types.np_bool_, types.np_bool_)]
    cases += [signature(types.np_bool_, op) for op in sorted(types.np_signed_domain)]
    cases += [signature(types.np_bool_, op) for op in sorted(types.np_unsigned_domain)]
    cases += [signature(types.np_bool_, op) for op in sorted(types.np_real_domain)]
    cases += [signature(types.np_bool_, op) for op in sorted(types.np_complex_domain)]


class OrderedCmpOp(ConcreteTemplate):
    cases = [signature(types.py_bool, types.py_bool, types.py_bool)]
    cases += [signature(types.py_bool, op, op) for op in sorted(types.py_signed_domain)]
    cases += [signature(types.py_bool, op, op) for op in sorted(types.py_real_domain)]
    cases = [signature(types.np_bool_, types.np_bool_, types.np_bool_)]
    cases += [signature(types.np_bool_, op, op) for op in sorted(types.np_signed_domain)]
    cases += [signature(types.np_bool_, op, op) for op in sorted(types.np_unsigned_domain)]
    cases += [signature(types.np_bool_, op, op) for op in sorted(types.np_real_domain)]


class UnorderedCmpOp(ConcreteTemplate):
    cases = OrderedCmpOp.cases + [
        signature(types.py_bool, op, op) for op in sorted(types.py_complex_domain)] + [
        signature(types.np_bool_, op, op) for op in sorted(types.np_complex_domain)]


@infer_global(operator.lt)
class CmpOpLt(OrderedCmpOp):
    pass


@infer_global(operator.le)
class CmpOpLe(OrderedCmpOp):
    pass


@infer_global(operator.gt)
class CmpOpGt(OrderedCmpOp):
    pass


@infer_global(operator.ge)
class CmpOpGe(OrderedCmpOp):
    pass


# more specific overloads should be registered first
@infer_global(operator.eq)
class ConstOpEq(AbstractTemplate):
    def generic(self, args, kws):
        assert not kws
        (arg1, arg2) = args
        if isinstance(arg1, types.Literal) and isinstance(arg2, types.Literal):
            return signature(types.np_bool_, arg1, arg2)


@infer_global(operator.ne)
class ConstOpNotEq(ConstOpEq):
    pass


@infer_global(operator.eq)
class CmpOpEq(UnorderedCmpOp):
    pass


@infer_global(operator.ne)
class CmpOpNe(UnorderedCmpOp):
    pass


class TupleCompare(AbstractTemplate):
    def generic(self, args, kws):
        [lhs, rhs] = args
        if isinstance(lhs, types.BaseTuple) and isinstance(rhs, types.BaseTuple):
            for u, v in zip(lhs, rhs):
                # Check element-wise comparability
                res = self.context.resolve_function_type(self.key, (u, v), {})
                if res is None:
                    break
            else:
                return signature(types.py_bool, lhs, rhs)


@infer_global(operator.eq)
class TupleEq(TupleCompare):
    pass


@infer_global(operator.ne)
class TupleNe(TupleCompare):
    pass


@infer_global(operator.ge)
class TupleGe(TupleCompare):
    pass


@infer_global(operator.gt)
class TupleGt(TupleCompare):
    pass


@infer_global(operator.le)
class TupleLe(TupleCompare):
    pass


@infer_global(operator.lt)
class TupleLt(TupleCompare):
    pass


@infer_global(operator.add)
class TupleAdd(AbstractTemplate):
    def generic(self, args, kws):
        if len(args) == 2:
            a, b = args
            if (isinstance(a, types.BaseTuple) and isinstance(b, types.BaseTuple)
                and not isinstance(a, types.BaseNamedTuple)
                and not isinstance(b, types.BaseNamedTuple)):
                res = types.BaseTuple.from_types(tuple(a) + tuple(b))
                return signature(res, a, b)


class CmpOpIdentity(AbstractTemplate):
    def generic(self, args, kws):
        [lhs, rhs] = args
        return signature(types.py_bool, lhs, rhs)


@infer_global(operator.is_)
class CmpOpIs(CmpOpIdentity):
    pass


@infer_global(operator.is_not)
class CmpOpIsNot(CmpOpIdentity):
    pass


def normalize_1d_index(index):
    """
    Normalize the *index* type (an integer or slice) for indexing a 1D
    sequence.
    """
    if isinstance(index, types.SliceType):
        return index

    elif isinstance(index, types.Integer):
        return types.np_intp if index.signed else types.uintp


@infer_global(operator.getitem)
class GetItemCPointer(AbstractTemplate):
    def generic(self, args, kws):
        assert not kws
        ptr, idx = args
        if isinstance(ptr, types.CPointer) and isinstance(idx, types.Integer):
            return signature(ptr.dtype, ptr, normalize_1d_index(idx))


@infer_global(operator.setitem)
class SetItemCPointer(AbstractTemplate):
    def generic(self, args, kws):
        assert not kws
        ptr, idx, val = args
        if isinstance(ptr, types.CPointer) and isinstance(idx, types.Integer):
            return signature(types.none, ptr, normalize_1d_index(idx), ptr.dtype)


@infer_global(len)
class Len(AbstractTemplate):
    def generic(self, args, kws):
        assert not kws
        (val,) = args
        if isinstance(val, (types.Buffer, types.BaseTuple)):
            return signature(types.py_int, val)
        elif isinstance(val, (types.RangeType)):
            return signature(val.dtype, val)

@infer_global(tuple)
class TupleConstructor(AbstractTemplate):
    def generic(self, args, kws):
        assert not kws
        # empty tuple case
        if len(args) == 0:
            return signature(types.Tuple(()))
        (val,) = args
        # tuple as input
        if isinstance(val, types.BaseTuple):
            return signature(val, val)


@infer_global(operator.contains)
class Contains(AbstractTemplate):
    def generic(self, args, kws):
        assert not kws
        (seq, val) = args

        if isinstance(seq, (types.Sequence)):
            return signature(types.py_bool, seq, val)

@infer_global(operator.truth)
class TupleBool(AbstractTemplate):
    def generic(self, args, kws):
        assert not kws
        (val,) = args
        if isinstance(val, (types.BaseTuple)):
            return signature(types.py_bool, val)


@infer
class StaticGetItemTuple(AbstractTemplate):
    key = "static_getitem"

    def generic(self, args, kws):
        tup, idx = args
        ret = None
        if not isinstance(tup, types.BaseTuple):
            return
        if isinstance(idx, int):
            try:
                ret = tup.types[idx]
            except IndexError:
                raise errors.NumbaIndexError("tuple index out of range")
        elif isinstance(idx, slice):
            ret = types.BaseTuple.from_types(tup.types[idx])
        if ret is not None:
            sig = signature(ret, *args)
            return sig


@infer
class StaticGetItemLiteralList(AbstractTemplate):
    key = "static_getitem"

    def generic(self, args, kws):
        tup, idx = args
        ret = None
        if not isinstance(tup, types.LiteralList):
            return
        if isinstance(idx, int):
            ret = tup.types[idx]
        if ret is not None:
            sig = signature(ret, *args)
            return sig


@infer
class StaticGetItemLiteralStrKeyDict(AbstractTemplate):
    key = "static_getitem"

    def generic(self, args, kws):
        tup, idx = args
        ret = None
        if not isinstance(tup, types.LiteralStrKeyDict):
            return
        if isinstance(idx, str):
            if idx in tup.fields:
                lookup = tup.fields.index(idx)
            else:
                raise errors.NumbaKeyError(f"Key '{idx}' is not in dict.")
            ret = tup.types[lookup]
        if ret is not None:
            sig = signature(ret, *args)
            return sig

@infer
class StaticGetItemClass(AbstractTemplate):
    """This handles the "static_getitem" when a Numba type is subscripted e.g:
    var = typed.List.empty_list(float64[::1, :])
    It only allows this on simple numerical types. Compound types, like
    records, are not supported.
    """
    key = "static_getitem"

    def generic(self, args, kws):
        clazz, idx = args
        if not isinstance(clazz, types.NumberClass):
            return
        ret = clazz.dtype[idx]
        sig = signature(ret, *args)
        return sig


# Generic implementation for "not in"

@infer
class GenericNotIn(AbstractTemplate):
    key = "not in"

    def generic(self, args, kws):
        args = args[::-1]
        sig = self.context.resolve_function_type(operator.contains, args, kws)
        return signature(sig.return_type, *sig.args[::-1])


#-------------------------------------------------------------------------------

@infer_getattr
class MemoryViewAttribute(AttributeTemplate):
    key = types.MemoryView

    def resolve_contiguous(self, buf):
        return types.py_bool

    def resolve_c_contiguous(self, buf):
        return types.py_bool

    def resolve_f_contiguous(self, buf):
        return types.py_bool

    def resolve_itemsize(self, buf):
        return types.py_int

    def resolve_nbytes(self, buf):
        return types.py_int

    def resolve_readonly(self, buf):
        return types.py_bool

    def resolve_shape(self, buf):
        return types.UniTuple(types.py_int, buf.ndim)

    def resolve_strides(self, buf):
        return types.UniTuple(types.py_int, buf.ndim)

    def resolve_ndim(self, buf):
        return types.py_int


#-------------------------------------------------------------------------------


@infer_getattr
class BooleanAttribute(AttributeTemplate):
    key = types.Boolean

    def resolve___class__(self, ty):
        return types.NumberClass(ty)

    @bound_function("number.item")
    def resolve_item(self, ty, args, kws):
        assert not kws
        if not args:
            return signature(ty)


@infer_getattr
class NumberAttribute(AttributeTemplate):
    key = types.Number

    def resolve___class__(self, ty):
        return types.NumberClass(ty)

    def resolve_real(self, ty):
        return getattr(ty, "underlying_float", ty)

    def resolve_imag(self, ty):
        return getattr(ty, "underlying_float", ty)

    @bound_function("complex.conjugate")
    def resolve_conjugate(self, ty, args, kws):
        assert not args
        assert not kws
        return signature(ty)

    @bound_function("number.item")
    def resolve_item(self, ty, args, kws):
        assert not kws
        if not args:
            return signature(ty)


@infer_getattr
class NPTimedeltaAttribute(AttributeTemplate):
    key = types.NPTimedelta

    def resolve___class__(self, ty):
        return types.NumberClass(ty)


@infer_getattr
class NPDatetimeAttribute(AttributeTemplate):
    key = types.NPDatetime

    def resolve___class__(self, ty):
        return types.NumberClass(ty)


@infer_getattr
class SliceAttribute(AttributeTemplate):
    key = types.SliceType

    def resolve_start(self, ty):
        return types.py_int

    def resolve_stop(self, ty):
        return types.py_int

    def resolve_step(self, ty):
        return types.py_int

    @bound_function("slice.indices")
    def resolve_indices(self, ty, args, kws):
        assert not kws
        if len(args) != 1:
            raise errors.NumbaTypeError(
                "indices() takes exactly one argument (%d given)" % len(args)
            )
        typ, = args
        if not isinstance(typ, types.Integer):
            raise errors.NumbaTypeError(
                "'%s' object cannot be interpreted as an integer" % typ
            )
        return signature(types.UniTuple(types.py_int, 3), types.py_int)


#-------------------------------------------------------------------------------


@infer_getattr
class NumberClassAttribute(AttributeTemplate):
    key = types.NumberClass

    def resolve___call__(self, classty):
        """
        Resolve a NumPy number class's constructor (e.g. calling numpy.int32(...))
        """
        ty = classty.instance_type

        def typer(val):
            if isinstance(val, (types.BaseTuple, types.Sequence)):
                # Array constructor, e.g. np.int32([1, 2])
                fnty = self.context.resolve_value_type(np.array)
                sig = fnty.get_call_type(self.context, (val, types.DType(ty)),
                                         {})
                return sig.return_type
            elif isinstance(val, (types.Number, types.Boolean, types.IntEnumMember)):
                 # Scalar constructor, e.g. np.int32(42)
                 return ty
            elif isinstance(val, (types.NPDatetime, types.NPTimedelta)):
                # Constructor cast from datetime-like, e.g.
                # > np.int64(np.datetime64("2000-01-01"))
                if ty.bitwidth == 64:
                    return ty
                else:
                    msg = (f"Cannot cast {val} to {ty} as {ty} is not 64 bits "
                           "wide.")
                    raise errors.TypingError(msg)
            else:
                if (isinstance(val, types.Array) and val.ndim == 0 and
                    val.dtype == ty):
                    # This is 0d array -> scalar degrading
                    return ty
                else:
                    # unsupported
                    msg = f"Casting {val} to {ty} directly is unsupported."
                    if isinstance(val, types.Array):
                        # array casts are supported a different way.
                        msg += f" Try doing '<array>.astype(np.{ty})' instead"
                    raise errors.TypingError(msg)

        return types.Function(make_callable_template(key=ty, typer=typer))


@infer_getattr
class TypeRefAttribute(AttributeTemplate):
    key = types.TypeRef

    def resolve___call__(self, classty):
        """
        Resolve a core number's constructor (e.g. calling int(...))

        Note:

        This is needed because of the limitation of the current type-system
        implementation.  Specifically, the lack of a higher-order type
        (i.e. passing the ``DictType`` vs ``DictType(key_type, value_type)``)
        """
        ty = classty.instance_type

        if isinstance(ty, type) and issubclass(ty, types.Type):
            # Redirect the typing to a:
            #   @type_callable(ty)
            #   def typeddict_call(context):
            #        ...
            # For example, see numba/typed/typeddict.py
            #   @type_callable(DictType)
            #   def typeddict_call(context):
            class Redirect(object):

                def __init__(self, context):
                    self.context =  context

                def __call__(self, *args, **kwargs):
                    result = self.context.resolve_function_type(ty, args, kwargs)
                    if hasattr(result, "pysig"):
                        self.pysig = result.pysig
                    return result

            return types.Function(make_callable_template(key=ty,
                                                         typer=Redirect(self.context)))


#------------------------------------------------------------------------------


class MinMaxBase(AbstractTemplate):

    def _unify_minmax(self, tys):
        for ty in tys:
            if not isinstance(ty, (types.Number, types.NPDatetime, types.NPTimedelta)):
                return
        return self.context.unify_types(*tys)

    def generic(self, args, kws):
        """
        Resolve a min() or max() call.
        """
        assert not kws

        if not args:
            return
        if len(args) == 1:
            # max(arg) only supported if arg is an iterable
            if isinstance(args[0], types.BaseTuple):
                tys = list(args[0])
                if not tys:
                    raise errors.TypingError("%s() argument is an empty tuple"
                                             % (self.key.__name__,))
            else:
                return
        else:
            # max(*args)
            tys = args
        retty = self._unify_minmax(tys)
        if retty is not None:
            return signature(retty, *args)


@infer_global(max)
class Max(MinMaxBase):
    pass


@infer_global(min)
class Min(MinMaxBase):
    pass


@infer_global(round)
class Round(ConcreteTemplate):
    cases = [
        signature(types.py_int, types.py_float),
        signature(types.py_float, types.py_float, types.py_int)
    ]


#------------------------------------------------------------------------------


@infer_global(bool)
class Bool(AbstractTemplate):

    def generic(self, args, kws):
        assert not kws
        [arg] = args
        if isinstance(arg, (types.Boolean, types.Number)):
            return signature(types.py_bool, arg)
        # XXX typing for bool cannot be polymorphic because of the
        # types.Function thing, so we redirect to the operator.truth
        # intrinsic.
        return self.context.resolve_function_type(operator.truth, args, kws)


@infer_global(int)
class Int(AbstractTemplate):

    def generic(self, args, kws):
        if kws:
            raise errors.NumbaAssertionError('kws not supported')

        [arg] = args

        if isinstance(arg, (types.Integer, types.Float,
                            types.Boolean, types.NPDatetime,
                            types.NPTimedelta)):
            return signature(types.py_int, arg)


@infer_global(float)
class Float(AbstractTemplate):

    def generic(self, args, kws):
        assert not kws

        [arg] = args

        if arg not in types.py_number_domain or arg not in types.np_number_domain:
            raise errors.NumbaTypeError("float() only support for numbers")

        if arg in types.py_complex_domain or arg in types.np_complex_domain:
            raise errors.NumbaTypeError("float() does not support complex")

        return signature(types.py_float, arg)



@infer_global(complex)
class Complex(AbstractTemplate):

    def generic(self, args, kws):
        assert not kws
        number_domain = types.py_number_domain | types.np_number_domain

        if len(args) == 1:
            [arg] = args
            if arg not in number_domain:
                raise errors.NumbaTypeError("complex() only support for numbers")

            return signature(types.py_complex, arg)

        elif len(args) == 2:
            [real, imag] = args
            if (real not in number_domain or
                imag not in number_domain):
                raise errors.NumbaTypeError("complex() only support for numbers")

            return signature(types.py_complex, real, imag)


#------------------------------------------------------------------------------

@infer_global(enumerate)
class Enumerate(AbstractTemplate):

    def generic(self, args, kws):
        assert not kws
        it = args[0]
        if len(args) > 1 and not isinstance(args[1], types.Integer):
            raise errors.NumbaTypeError("Only integers supported as start "
                                        "value in enumerate")
        elif len(args) > 2:
            #let python raise its own error
            enumerate(*args)

        if isinstance(it, types.IterableType):
            enumerate_type = types.EnumerateType(it)
            return signature(enumerate_type, *args)


@infer_global(zip)
class Zip(AbstractTemplate):

    def generic(self, args, kws):
        assert not kws
        if all(isinstance(it, types.IterableType) for it in args):
            zip_type = types.ZipType(args)
            return signature(zip_type, *args)


@infer_global(iter)
class Iter(AbstractTemplate):

    def generic(self, args, kws):
        assert not kws
        if len(args) == 1:
            it = args[0]
            if isinstance(it, types.IterableType):
                return signature(it.iterator_type, *args)


@infer_global(next)
class Next(AbstractTemplate):

    def generic(self, args, kws):
        assert not kws
        if len(args) == 1:
            it = args[0]
            if isinstance(it, types.IteratorType):
                return signature(it.yield_type, *args)


#------------------------------------------------------------------------------

@infer_global(type)
class TypeBuiltin(AbstractTemplate):

    def generic(self, args, kws):
        assert not kws
        if len(args) == 1:
            # One-argument type() -> return the __class__
            # Avoid literal types
            arg = types.unliteral(args[0])
            classty = self.context.resolve_getattr(arg, "__class__")
            if classty is not None:
                return signature(classty, *args)


#------------------------------------------------------------------------------

@infer_getattr
class OptionalAttribute(AttributeTemplate):
    key = types.Optional

    def generic_resolve(self, optional, attr):
        return self.context.resolve_getattr(optional.type, attr)

#------------------------------------------------------------------------------

@infer_getattr
class DeferredAttribute(AttributeTemplate):
    key = types.DeferredType

    def generic_resolve(self, deferred, attr):
        return self.context.resolve_getattr(deferred.get(), attr)

#------------------------------------------------------------------------------


class IndexValue(object):
    """
    Index and value
    """
    def __init__(self, ind, val):
        self.index = ind
        self.value = val

    def __repr__(self):
        return 'IndexValue(%f, %f)' % (self.index, self.value)


class IndexValueType(types.Type):
    def __init__(self, val_typ):
        self.val_typ = val_typ
        super(IndexValueType, self).__init__(
                                    name='IndexValueType({})'.format(val_typ))


@typeof_impl.register(IndexValue)
def typeof_index(val, c):
    val_typ = typeof_impl(val.value, c)
    return IndexValueType(val_typ)


@type_callable(IndexValue)
def type_index_value(context):
    def typer(ind, mval):
        if ind == types.np_intp or ind == types.uintp:
            return IndexValueType(mval)
    return typer


@register_model(IndexValueType)
class IndexValueModel(models.StructModel):
    def __init__(self, dmm, fe_type):
        members = [
            ('index', types.np_intp),
            ('value', fe_type.val_typ),
            ]
        models.StructModel.__init__(self, dmm, fe_type, members)


make_attribute_wrapper(IndexValueType, 'index', 'index')
make_attribute_wrapper(IndexValueType, 'value', 'value')
