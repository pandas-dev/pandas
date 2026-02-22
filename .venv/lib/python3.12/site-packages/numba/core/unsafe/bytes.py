"""
This file provides internal compiler utilities that support certain special
operations with bytes and workarounds for limitations enforced in userland.
"""

from numba.core.extending import intrinsic
from llvmlite import ir
from numba.core import types, cgutils


@intrinsic
def grab_byte(typingctx, data, offset):
    # returns a byte at a given offset in data
    def impl(context, builder, signature, args):
        data, idx = args
        ptr = builder.bitcast(data, ir.IntType(8).as_pointer())
        ch = builder.load(builder.gep(ptr, [idx]))
        return ch

    sig = types.uint8(types.voidptr, types.intp)
    return sig, impl


@intrinsic
def grab_uint64_t(typingctx, data, offset):
    # returns a uint64_t at a given offset in data
    def impl(context, builder, signature, args):
        data, idx = args
        ptr = builder.bitcast(data, ir.IntType(64).as_pointer())
        ch = builder.load(builder.gep(ptr, [idx]))
        return ch
    sig = types.uint64(types.voidptr, types.intp)
    return sig, impl


@intrinsic
def memcpy_region(typingctx, dst, dst_offset, src, src_offset, nbytes, align):
    '''Copy nbytes from *(src + src_offset) to *(dst + dst_offset)'''
    def codegen(context, builder, signature, args):
        [dst_val, dst_offset_val, src_val, src_offset_val, nbytes_val,
         align_val] = args
        src_ptr = builder.gep(src_val, [src_offset_val])
        dst_ptr = builder.gep(dst_val, [dst_offset_val])
        cgutils.raw_memcpy(builder, dst_ptr, src_ptr, nbytes_val, align_val)
        return context.get_dummy_value()

    sig = types.void(types.voidptr, types.intp, types.voidptr, types.intp,
                     types.intp, types.intp)
    return sig, codegen
