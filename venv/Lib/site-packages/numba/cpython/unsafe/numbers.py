""" This module provides the unsafe things for targets/numbers.py
"""
from numba.core import types, errors
from numba.core.extending import intrinsic

from llvmlite import ir


@intrinsic
def viewer(tyctx, val, viewty):
    """ Bitcast a scalar 'val' to the given type 'viewty'. """
    bits = val.bitwidth
    if isinstance(viewty.dtype, types.Integer):
        bitcastty = ir.IntType(bits)
    elif isinstance(viewty.dtype, types.Float):
        bitcastty = ir.FloatType() if bits == 32 else ir.DoubleType()
    else:
        assert 0, "unreachable"

    def codegen(cgctx, builder, typ, args):
        flt = args[0]
        return builder.bitcast(flt, bitcastty)
    retty = viewty.dtype
    sig = retty(val, viewty)
    return sig, codegen


@intrinsic
def trailing_zeros(typeingctx, src):
    """Counts trailing zeros in the binary representation of an integer."""
    if not isinstance(src, types.Integer):
        msg = ("trailing_zeros is only defined for integers, but value passed "
               f"was '{src}'.")
        raise errors.NumbaTypeError(msg)

    def codegen(context, builder, signature, args):
        [src] = args
        return builder.cttz(src, ir.Constant(ir.IntType(1), 0))
    return src(src), codegen


@intrinsic
def leading_zeros(typeingctx, src):
    """Counts leading zeros in the binary representation of an integer."""
    if not isinstance(src, types.Integer):
        msg = ("leading_zeros is only defined for integers, but value passed "
               f"was '{src}'.")
        raise errors.NumbaTypeError(msg)

    def codegen(context, builder, signature, args):
        [src] = args
        return builder.ctlz(src, ir.Constant(ir.IntType(1), 0))
    return src(src), codegen
