from functools import reduce
import operator
import math

from llvmlite import ir
import llvmlite.binding as ll

from numba.core.imputils import Registry, lower_cast
from numba.core.typing.npydecl import parse_dtype
from numba.core.datamodel import models
from numba.core import types, cgutils
from numba.np import ufunc_db
from numba.np.npyimpl import register_ufuncs
from .cudadrv import nvvm
from numba import cuda
from numba.cuda import nvvmutils, stubs, errors
from numba.cuda.types import dim3, CUDADispatcher

registry = Registry()
lower = registry.lower
lower_attr = registry.lower_getattr
lower_constant = registry.lower_constant


def initialize_dim3(builder, prefix):
    x = nvvmutils.call_sreg(builder, "%s.x" % prefix)
    y = nvvmutils.call_sreg(builder, "%s.y" % prefix)
    z = nvvmutils.call_sreg(builder, "%s.z" % prefix)
    return cgutils.pack_struct(builder, (x, y, z))


@lower_attr(types.Module(cuda), 'threadIdx')
def cuda_threadIdx(context, builder, sig, args):
    return initialize_dim3(builder, 'tid')


@lower_attr(types.Module(cuda), 'blockDim')
def cuda_blockDim(context, builder, sig, args):
    return initialize_dim3(builder, 'ntid')


@lower_attr(types.Module(cuda), 'blockIdx')
def cuda_blockIdx(context, builder, sig, args):
    return initialize_dim3(builder, 'ctaid')


@lower_attr(types.Module(cuda), 'gridDim')
def cuda_gridDim(context, builder, sig, args):
    return initialize_dim3(builder, 'nctaid')


@lower_attr(types.Module(cuda), 'laneid')
def cuda_laneid(context, builder, sig, args):
    return nvvmutils.call_sreg(builder, 'laneid')


@lower_attr(dim3, 'x')
def dim3_x(context, builder, sig, args):
    return builder.extract_value(args, 0)


@lower_attr(dim3, 'y')
def dim3_y(context, builder, sig, args):
    return builder.extract_value(args, 1)


@lower_attr(dim3, 'z')
def dim3_z(context, builder, sig, args):
    return builder.extract_value(args, 2)


# -----------------------------------------------------------------------------

@lower(cuda.const.array_like, types.Array)
def cuda_const_array_like(context, builder, sig, args):
    # This is a no-op because CUDATargetContext.make_constant_array already
    # created the constant array.
    return args[0]


_unique_smem_id = 0


def _get_unique_smem_id(name):
    """Due to bug with NVVM invalid internalizing of shared memory in the
    PTX output.  We can't mark shared memory to be internal. We have to
    ensure unique name is generated for shared memory symbol.
    """
    global _unique_smem_id
    _unique_smem_id += 1
    return "{0}_{1}".format(name, _unique_smem_id)


@lower(cuda.shared.array, types.IntegerLiteral, types.Any)
def cuda_shared_array_integer(context, builder, sig, args):
    length = sig.args[0].literal_value
    dtype = parse_dtype(sig.args[1])
    return _generic_array(context, builder, shape=(length,), dtype=dtype,
                          symbol_name=_get_unique_smem_id('_cudapy_smem'),
                          addrspace=nvvm.ADDRSPACE_SHARED,
                          can_dynsized=True)


@lower(cuda.shared.array, types.Tuple, types.Any)
@lower(cuda.shared.array, types.UniTuple, types.Any)
def cuda_shared_array_tuple(context, builder, sig, args):
    shape = [ s.literal_value for s in sig.args[0] ]
    dtype = parse_dtype(sig.args[1])
    return _generic_array(context, builder, shape=shape, dtype=dtype,
                          symbol_name=_get_unique_smem_id('_cudapy_smem'),
                          addrspace=nvvm.ADDRSPACE_SHARED,
                          can_dynsized=True)


@lower(cuda.local.array, types.IntegerLiteral, types.Any)
def cuda_local_array_integer(context, builder, sig, args):
    length = sig.args[0].literal_value
    dtype = parse_dtype(sig.args[1])
    return _generic_array(context, builder, shape=(length,), dtype=dtype,
                          symbol_name='_cudapy_lmem',
                          addrspace=nvvm.ADDRSPACE_LOCAL,
                          can_dynsized=False)


@lower(cuda.local.array, types.Tuple, types.Any)
@lower(cuda.local.array, types.UniTuple, types.Any)
def ptx_lmem_alloc_array(context, builder, sig, args):
    shape = [ s.literal_value for s in sig.args[0] ]
    dtype = parse_dtype(sig.args[1])
    return _generic_array(context, builder, shape=shape, dtype=dtype,
                          symbol_name='_cudapy_lmem',
                          addrspace=nvvm.ADDRSPACE_LOCAL,
                          can_dynsized=False)


@lower(stubs.threadfence_block)
def ptx_threadfence_block(context, builder, sig, args):
    assert not args
    fname = 'llvm.nvvm.membar.cta'
    lmod = builder.module
    fnty = ir.FunctionType(ir.VoidType(), ())
    sync = cgutils.get_or_insert_function(lmod, fnty, fname)
    builder.call(sync, ())
    return context.get_dummy_value()


@lower(stubs.threadfence_system)
def ptx_threadfence_system(context, builder, sig, args):
    assert not args
    fname = 'llvm.nvvm.membar.sys'
    lmod = builder.module
    fnty = ir.FunctionType(ir.VoidType(), ())
    sync = cgutils.get_or_insert_function(lmod, fnty, fname)
    builder.call(sync, ())
    return context.get_dummy_value()


@lower(stubs.threadfence)
def ptx_threadfence_device(context, builder, sig, args):
    assert not args
    fname = 'llvm.nvvm.membar.gl'
    lmod = builder.module
    fnty = ir.FunctionType(ir.VoidType(), ())
    sync = cgutils.get_or_insert_function(lmod, fnty, fname)
    builder.call(sync, ())
    return context.get_dummy_value()


@lower(stubs.syncwarp)
def ptx_syncwarp(context, builder, sig, args):
    mask = context.get_constant(types.int32, 0xFFFFFFFF)
    mask_sig = types.none(types.int32)
    return ptx_syncwarp_mask(context, builder, mask_sig, [mask])


@lower(stubs.syncwarp, types.i4)
def ptx_syncwarp_mask(context, builder, sig, args):
    fname = 'llvm.nvvm.bar.warp.sync'
    lmod = builder.module
    fnty = ir.FunctionType(ir.VoidType(), (ir.IntType(32),))
    sync = cgutils.get_or_insert_function(lmod, fnty, fname)
    builder.call(sync, args)
    return context.get_dummy_value()


@lower(stubs.shfl_sync_intrinsic, types.i4, types.i4, types.i4, types.i4,
       types.i4)
@lower(stubs.shfl_sync_intrinsic, types.i4, types.i4, types.i8, types.i4,
       types.i4)
@lower(stubs.shfl_sync_intrinsic, types.i4, types.i4, types.f4, types.i4,
       types.i4)
@lower(stubs.shfl_sync_intrinsic, types.i4, types.i4, types.f8, types.i4,
       types.i4)
def ptx_shfl_sync_i32(context, builder, sig, args):
    """
    The NVVM intrinsic for shfl only supports i32, but the cuda intrinsic
    function supports both 32 and 64 bit ints and floats, so for feature parity,
    i64, f32, and f64 are implemented. Floats by way of bitcasting the float to
    an int, then shuffling, then bitcasting back. And 64-bit values by packing
    them into 2 32bit values, shuffling thoose, and then packing back together.
    """
    mask, mode, value, index, clamp = args
    value_type = sig.args[2]
    if value_type in types.real_domain:
        value = builder.bitcast(value, ir.IntType(value_type.bitwidth))
    fname = 'llvm.nvvm.shfl.sync.i32'
    lmod = builder.module
    fnty = ir.FunctionType(
        ir.LiteralStructType((ir.IntType(32), ir.IntType(1))),
                            (ir.IntType(32), ir.IntType(32), ir.IntType(32),
                             ir.IntType(32), ir.IntType(32))
    )
    func = cgutils.get_or_insert_function(lmod, fnty, fname)
    if value_type.bitwidth == 32:
        ret = builder.call(func, (mask, mode, value, index, clamp))
        if value_type == types.float32:
            rv = builder.extract_value(ret, 0)
            pred = builder.extract_value(ret, 1)
            fv = builder.bitcast(rv, ir.FloatType())
            ret = cgutils.make_anonymous_struct(builder, (fv, pred))
    else:
        value1 = builder.trunc(value, ir.IntType(32))
        value_lshr = builder.lshr(value, context.get_constant(types.i8, 32))
        value2 = builder.trunc(value_lshr, ir.IntType(32))
        ret1 = builder.call(func, (mask, mode, value1, index, clamp))
        ret2 = builder.call(func, (mask, mode, value2, index, clamp))
        rv1 = builder.extract_value(ret1, 0)
        rv2 = builder.extract_value(ret2, 0)
        pred = builder.extract_value(ret1, 1)
        rv1_64 = builder.zext(rv1, ir.IntType(64))
        rv2_64 = builder.zext(rv2, ir.IntType(64))
        rv_shl = builder.shl(rv2_64, context.get_constant(types.i8, 32))
        rv = builder.or_(rv_shl, rv1_64)
        if value_type == types.float64:
            rv = builder.bitcast(rv, ir.DoubleType())
        ret = cgutils.make_anonymous_struct(builder, (rv, pred))
    return ret


@lower(stubs.vote_sync_intrinsic, types.i4, types.i4, types.boolean)
def ptx_vote_sync(context, builder, sig, args):
    fname = 'llvm.nvvm.vote.sync'
    lmod = builder.module
    fnty = ir.FunctionType(ir.LiteralStructType((ir.IntType(32),
                                                 ir.IntType(1))),
                           (ir.IntType(32), ir.IntType(32), ir.IntType(1)))
    func = cgutils.get_or_insert_function(lmod, fnty, fname)
    return builder.call(func, args)


@lower(stubs.match_any_sync, types.i4, types.i4)
@lower(stubs.match_any_sync, types.i4, types.i8)
@lower(stubs.match_any_sync, types.i4, types.f4)
@lower(stubs.match_any_sync, types.i4, types.f8)
def ptx_match_any_sync(context, builder, sig, args):
    mask, value = args
    width = sig.args[1].bitwidth
    if sig.args[1] in types.real_domain:
        value = builder.bitcast(value, ir.IntType(width))
    fname = 'llvm.nvvm.match.any.sync.i{}'.format(width)
    lmod = builder.module
    fnty = ir.FunctionType(ir.IntType(32), (ir.IntType(32), ir.IntType(width)))
    func = cgutils.get_or_insert_function(lmod, fnty, fname)
    return builder.call(func, (mask, value))


@lower(stubs.match_all_sync, types.i4, types.i4)
@lower(stubs.match_all_sync, types.i4, types.i8)
@lower(stubs.match_all_sync, types.i4, types.f4)
@lower(stubs.match_all_sync, types.i4, types.f8)
def ptx_match_all_sync(context, builder, sig, args):
    mask, value = args
    width = sig.args[1].bitwidth
    if sig.args[1] in types.real_domain:
        value = builder.bitcast(value, ir.IntType(width))
    fname = 'llvm.nvvm.match.all.sync.i{}'.format(width)
    lmod = builder.module
    fnty = ir.FunctionType(ir.LiteralStructType((ir.IntType(32),
                                                 ir.IntType(1))),
                           (ir.IntType(32), ir.IntType(width)))
    func = cgutils.get_or_insert_function(lmod, fnty, fname)
    return builder.call(func, (mask, value))


@lower(stubs.activemask)
def ptx_activemask(context, builder, sig, args):
    activemask = ir.InlineAsm(ir.FunctionType(ir.IntType(32), []),
                              "activemask.b32 $0;", '=r', side_effect=True)
    return builder.call(activemask, [])


@lower(stubs.lanemask_lt)
def ptx_lanemask_lt(context, builder, sig, args):
    activemask = ir.InlineAsm(ir.FunctionType(ir.IntType(32), []),
                              "mov.u32 $0, %lanemask_lt;", '=r',
                              side_effect=True)
    return builder.call(activemask, [])


@lower(stubs.popc, types.Any)
def ptx_popc(context, builder, sig, args):
    return builder.ctpop(args[0])


@lower(stubs.fma, types.Any, types.Any, types.Any)
def ptx_fma(context, builder, sig, args):
    return builder.fma(*args)


def float16_float_ty_constraint(bitwidth):
    typemap = {32: ('f32', 'f'), 64: ('f64', 'd')}

    try:
        return typemap[bitwidth]
    except KeyError:
        msg = f"Conversion between float16 and float{bitwidth} unsupported"
        raise errors.CudaLoweringError(msg)


@lower_cast(types.float16, types.Float)
def float16_to_float_cast(context, builder, fromty, toty, val):
    if fromty.bitwidth == toty.bitwidth:
        return val

    ty, constraint = float16_float_ty_constraint(toty.bitwidth)

    fnty = ir.FunctionType(context.get_value_type(toty), [ir.IntType(16)])
    asm = ir.InlineAsm(fnty, f"cvt.{ty}.f16 $0, $1;", f"={constraint},h")
    return builder.call(asm, [val])


@lower_cast(types.Float, types.float16)
def float_to_float16_cast(context, builder, fromty, toty, val):
    if fromty.bitwidth == toty.bitwidth:
        return val

    ty, constraint = float16_float_ty_constraint(fromty.bitwidth)

    fnty = ir.FunctionType(ir.IntType(16), [context.get_value_type(fromty)])
    asm = ir.InlineAsm(fnty, f"cvt.rn.f16.{ty} $0, $1;", f"=h,{constraint}")
    return builder.call(asm, [val])


def float16_int_constraint(bitwidth):
    typemap = { 8: 'c', 16: 'h', 32: 'r', 64: 'l' }

    try:
        return typemap[bitwidth]
    except KeyError:
        msg = f"Conversion between float16 and int{bitwidth} unsupported"
        raise errors.CudaLoweringError(msg)


@lower_cast(types.float16, types.Integer)
def float16_to_integer_cast(context, builder, fromty, toty, val):
    bitwidth = toty.bitwidth
    constraint = float16_int_constraint(bitwidth)
    signedness = 's' if toty.signed else 'u'

    fnty = ir.FunctionType(context.get_value_type(toty), [ir.IntType(16)])
    asm = ir.InlineAsm(fnty,
                       f"cvt.rni.{signedness}{bitwidth}.f16 $0, $1;",
                       f"={constraint},h")
    return builder.call(asm, [val])


@lower_cast(types.Integer, types.float16)
@lower_cast(types.IntegerLiteral, types.float16)
def integer_to_float16_cast(context, builder, fromty, toty, val):
    bitwidth = fromty.bitwidth
    constraint = float16_int_constraint(bitwidth)
    signedness = 's' if fromty.signed else 'u'

    fnty = ir.FunctionType(ir.IntType(16),
                           [context.get_value_type(fromty)])
    asm = ir.InlineAsm(fnty,
                       f"cvt.rn.f16.{signedness}{bitwidth} $0, $1;",
                       f"=h,{constraint}")
    return builder.call(asm, [val])


def lower_fp16_binary(fn, op):
    @lower(fn, types.float16, types.float16)
    def ptx_fp16_binary(context, builder, sig, args):
        fnty = ir.FunctionType(ir.IntType(16),
                               [ir.IntType(16), ir.IntType(16)])
        asm = ir.InlineAsm(fnty, f'{op}.f16 $0,$1,$2;', '=h,h,h')
        return builder.call(asm, args)


lower_fp16_binary(stubs.fp16.hadd, 'add')
lower_fp16_binary(operator.add, 'add')
lower_fp16_binary(operator.iadd, 'add')
lower_fp16_binary(stubs.fp16.hsub, 'sub')
lower_fp16_binary(operator.sub, 'sub')
lower_fp16_binary(operator.isub, 'sub')
lower_fp16_binary(stubs.fp16.hmul, 'mul')
lower_fp16_binary(operator.mul, 'mul')
lower_fp16_binary(operator.imul, 'mul')


@lower(stubs.fp16.hneg, types.float16)
def ptx_fp16_hneg(context, builder, sig, args):
    fnty = ir.FunctionType(ir.IntType(16), [ir.IntType(16)])
    asm = ir.InlineAsm(fnty, 'neg.f16 $0, $1;', '=h,h')
    return builder.call(asm, args)


@lower(operator.neg, types.float16)
def operator_hneg(context, builder, sig, args):
    return ptx_fp16_hneg(context, builder, sig, args)


@lower(stubs.fp16.habs, types.float16)
def ptx_fp16_habs(context, builder, sig, args):
    fnty = ir.FunctionType(ir.IntType(16), [ir.IntType(16)])
    asm = ir.InlineAsm(fnty, 'abs.f16 $0, $1;', '=h,h')
    return builder.call(asm, args)


@lower(abs, types.float16)
def operator_habs(context, builder, sig, args):
    return ptx_fp16_habs(context, builder, sig, args)


@lower(stubs.fp16.hfma, types.float16, types.float16, types.float16)
def ptx_hfma(context, builder, sig, args):
    argtys = [ir.IntType(16), ir.IntType(16), ir.IntType(16)]
    fnty = ir.FunctionType(ir.IntType(16), argtys)
    asm = ir.InlineAsm(fnty, "fma.rn.f16 $0,$1,$2,$3;", "=h,h,h,h")
    return builder.call(asm, args)


@lower(operator.truediv, types.float16, types.float16)
@lower(operator.itruediv, types.float16, types.float16)
def fp16_div_impl(context, builder, sig, args):
    def fp16_div(x, y):
        return cuda.fp16.hdiv(x, y)

    return context.compile_internal(builder, fp16_div, sig, args)


_fp16_cmp = """{{
          .reg .pred __$$f16_cmp_tmp;
          setp.{op}.f16 __$$f16_cmp_tmp, $1, $2;
          selp.u16 $0, 1, 0, __$$f16_cmp_tmp;
        }}"""


def _gen_fp16_cmp(op):
    def ptx_fp16_comparison(context, builder, sig, args):
        fnty = ir.FunctionType(ir.IntType(16), [ir.IntType(16), ir.IntType(16)])
        asm = ir.InlineAsm(fnty, _fp16_cmp.format(op=op), '=h,h,h')
        result = builder.call(asm, args)

        zero = context.get_constant(types.int16, 0)
        int_result = builder.bitcast(result, ir.IntType(16))
        return builder.icmp_unsigned("!=", int_result, zero)
    return ptx_fp16_comparison


lower(stubs.fp16.heq, types.float16, types.float16)(_gen_fp16_cmp('eq'))
lower(operator.eq, types.float16, types.float16)(_gen_fp16_cmp('eq'))
lower(stubs.fp16.hne, types.float16, types.float16)(_gen_fp16_cmp('ne'))
lower(operator.ne, types.float16, types.float16)(_gen_fp16_cmp('ne'))
lower(stubs.fp16.hge, types.float16, types.float16)(_gen_fp16_cmp('ge'))
lower(operator.ge, types.float16, types.float16)(_gen_fp16_cmp('ge'))
lower(stubs.fp16.hgt, types.float16, types.float16)(_gen_fp16_cmp('gt'))
lower(operator.gt, types.float16, types.float16)(_gen_fp16_cmp('gt'))
lower(stubs.fp16.hle, types.float16, types.float16)(_gen_fp16_cmp('le'))
lower(operator.le, types.float16, types.float16)(_gen_fp16_cmp('le'))
lower(stubs.fp16.hlt, types.float16, types.float16)(_gen_fp16_cmp('lt'))
lower(operator.lt, types.float16, types.float16)(_gen_fp16_cmp('lt'))


def lower_fp16_minmax(fn, fname, op):
    @lower(fn, types.float16, types.float16)
    def ptx_fp16_minmax(context, builder, sig, args):
        choice = _gen_fp16_cmp(op)(context, builder, sig, args)
        return builder.select(choice, args[0], args[1])


lower_fp16_minmax(stubs.fp16.hmax, 'max', 'gt')
lower_fp16_minmax(stubs.fp16.hmin, 'min', 'lt')

# See:
# https://docs.nvidia.com/cuda/libdevice-users-guide/__nv_cbrt.html#__nv_cbrt
# https://docs.nvidia.com/cuda/libdevice-users-guide/__nv_cbrtf.html#__nv_cbrtf


cbrt_funcs = {
    types.float32: '__nv_cbrtf',
    types.float64: '__nv_cbrt',
}


@lower(stubs.cbrt, types.float32)
@lower(stubs.cbrt, types.float64)
def ptx_cbrt(context, builder, sig, args):
    ty = sig.return_type
    fname = cbrt_funcs[ty]
    fty = context.get_value_type(ty)
    lmod = builder.module
    fnty = ir.FunctionType(fty, [fty])
    fn = cgutils.get_or_insert_function(lmod, fnty, fname)
    return builder.call(fn, args)


@lower(stubs.brev, types.u4)
def ptx_brev_u4(context, builder, sig, args):
    # FIXME the llvm.bitreverse.i32 intrinsic isn't supported by nvcc
    # return builder.bitreverse(args[0])

    fn = cgutils.get_or_insert_function(
        builder.module,
        ir.FunctionType(ir.IntType(32), (ir.IntType(32),)),
        '__nv_brev')
    return builder.call(fn, args)


@lower(stubs.brev, types.u8)
def ptx_brev_u8(context, builder, sig, args):
    # FIXME the llvm.bitreverse.i64 intrinsic isn't supported by nvcc
    # return builder.bitreverse(args[0])

    fn = cgutils.get_or_insert_function(
        builder.module,
        ir.FunctionType(ir.IntType(64), (ir.IntType(64),)),
        '__nv_brevll')
    return builder.call(fn, args)


@lower(stubs.clz, types.Any)
def ptx_clz(context, builder, sig, args):
    return builder.ctlz(
        args[0],
        context.get_constant(types.boolean, 0))


@lower(stubs.ffs, types.i4)
@lower(stubs.ffs, types.u4)
def ptx_ffs_32(context, builder, sig, args):
    fn = cgutils.get_or_insert_function(
        builder.module,
        ir.FunctionType(ir.IntType(32), (ir.IntType(32),)),
        '__nv_ffs')
    return builder.call(fn, args)


@lower(stubs.ffs, types.i8)
@lower(stubs.ffs, types.u8)
def ptx_ffs_64(context, builder, sig, args):
    fn = cgutils.get_or_insert_function(
        builder.module,
        ir.FunctionType(ir.IntType(32), (ir.IntType(64),)),
        '__nv_ffsll')
    return builder.call(fn, args)


@lower(stubs.selp, types.Any, types.Any, types.Any)
def ptx_selp(context, builder, sig, args):
    test, a, b = args
    return builder.select(test, a, b)


@lower(max, types.f4, types.f4)
def ptx_max_f4(context, builder, sig, args):
    fn = cgutils.get_or_insert_function(
        builder.module,
        ir.FunctionType(
            ir.FloatType(),
            (ir.FloatType(), ir.FloatType())),
        '__nv_fmaxf')
    return builder.call(fn, args)


@lower(max, types.f8, types.f4)
@lower(max, types.f4, types.f8)
@lower(max, types.f8, types.f8)
def ptx_max_f8(context, builder, sig, args):
    fn = cgutils.get_or_insert_function(
        builder.module,
        ir.FunctionType(
            ir.DoubleType(),
            (ir.DoubleType(), ir.DoubleType())),
        '__nv_fmax')

    return builder.call(fn, [
        context.cast(builder, args[0], sig.args[0], types.double),
        context.cast(builder, args[1], sig.args[1], types.double),
    ])


@lower(min, types.f4, types.f4)
def ptx_min_f4(context, builder, sig, args):
    fn = cgutils.get_or_insert_function(
        builder.module,
        ir.FunctionType(
            ir.FloatType(),
            (ir.FloatType(), ir.FloatType())),
        '__nv_fminf')
    return builder.call(fn, args)


@lower(min, types.f8, types.f4)
@lower(min, types.f4, types.f8)
@lower(min, types.f8, types.f8)
def ptx_min_f8(context, builder, sig, args):
    fn = cgutils.get_or_insert_function(
        builder.module,
        ir.FunctionType(
            ir.DoubleType(),
            (ir.DoubleType(), ir.DoubleType())),
        '__nv_fmin')

    return builder.call(fn, [
        context.cast(builder, args[0], sig.args[0], types.double),
        context.cast(builder, args[1], sig.args[1], types.double),
    ])


@lower(round, types.f4)
@lower(round, types.f8)
def ptx_round(context, builder, sig, args):
    fn = cgutils.get_or_insert_function(
        builder.module,
        ir.FunctionType(
            ir.IntType(64),
            (ir.DoubleType(),)),
        '__nv_llrint')
    return builder.call(fn, [
        context.cast(builder, args[0], sig.args[0], types.double),
    ])


# This rounding implementation follows the algorithm used in the "fallback
# version" of double_round in CPython.
# https://github.com/python/cpython/blob/a755410e054e1e2390de5830befc08fe80706c66/Objects/floatobject.c#L964-L1007

@lower(round, types.f4, types.Integer)
@lower(round, types.f8, types.Integer)
def round_to_impl(context, builder, sig, args):
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

        else:
            pow1 = 10.0 ** (-ndigits)
            y = x / pow1

        z = round(y)
        if (math.fabs(y - z) == 0.5):
            # halfway between two integers; use round-half-even
            z = 2.0 * round(y / 2.0)

        if ndigits >= 0:
            z = (z / pow2) / pow1
        else:
            z *= pow1

        return z

    return context.compile_internal(builder, round_ndigits, sig, args, )


def gen_deg_rad(const):
    def impl(context, builder, sig, args):
        argty, = sig.args
        factor = context.get_constant(argty, const)
        return builder.fmul(factor, args[0])
    return impl


_deg2rad = math.pi / 180.
_rad2deg = 180. / math.pi
lower(math.radians, types.f4)(gen_deg_rad(_deg2rad))
lower(math.radians, types.f8)(gen_deg_rad(_deg2rad))
lower(math.degrees, types.f4)(gen_deg_rad(_rad2deg))
lower(math.degrees, types.f8)(gen_deg_rad(_rad2deg))


def _normalize_indices(context, builder, indty, inds, aryty, valty):
    """
    Convert integer indices into tuple of intp
    """
    if indty in types.integer_domain:
        indty = types.UniTuple(dtype=indty, count=1)
        indices = [inds]
    else:
        indices = cgutils.unpack_tuple(builder, inds, count=len(indty))
    indices = [context.cast(builder, i, t, types.intp)
               for t, i in zip(indty, indices)]

    dtype = aryty.dtype
    if dtype != valty:
        raise TypeError("expect %s but got %s" % (dtype, valty))

    if aryty.ndim != len(indty):
        raise TypeError("indexing %d-D array with %d-D index" %
                        (aryty.ndim, len(indty)))

    return indty, indices


def _atomic_dispatcher(dispatch_fn):
    def imp(context, builder, sig, args):
        # The common argument handling code
        aryty, indty, valty = sig.args
        ary, inds, val = args
        dtype = aryty.dtype

        indty, indices = _normalize_indices(context, builder, indty, inds,
                                            aryty, valty)

        lary = context.make_array(aryty)(context, builder, ary)
        ptr = cgutils.get_item_pointer(context, builder, aryty, lary, indices,
                                       wraparound=True)
        # dispatcher to implementation base on dtype
        return dispatch_fn(context, builder, dtype, ptr, val)
    return imp


@lower(stubs.atomic.add, types.Array, types.intp, types.Any)
@lower(stubs.atomic.add, types.Array, types.UniTuple, types.Any)
@lower(stubs.atomic.add, types.Array, types.Tuple, types.Any)
@_atomic_dispatcher
def ptx_atomic_add_tuple(context, builder, dtype, ptr, val):
    if dtype == types.float32:
        lmod = builder.module
        return builder.call(nvvmutils.declare_atomic_add_float32(lmod),
                            (ptr, val))
    elif dtype == types.float64:
        lmod = builder.module
        return builder.call(nvvmutils.declare_atomic_add_float64(lmod),
                            (ptr, val))
    else:
        return builder.atomic_rmw('add', ptr, val, 'monotonic')


@lower(stubs.atomic.sub, types.Array, types.intp, types.Any)
@lower(stubs.atomic.sub, types.Array, types.UniTuple, types.Any)
@lower(stubs.atomic.sub, types.Array, types.Tuple, types.Any)
@_atomic_dispatcher
def ptx_atomic_sub(context, builder, dtype, ptr, val):
    if dtype == types.float32:
        lmod = builder.module
        return builder.call(nvvmutils.declare_atomic_sub_float32(lmod),
                            (ptr, val))
    elif dtype == types.float64:
        lmod = builder.module
        return builder.call(nvvmutils.declare_atomic_sub_float64(lmod),
                            (ptr, val))
    else:
        return builder.atomic_rmw('sub', ptr, val, 'monotonic')


@lower(stubs.atomic.inc, types.Array, types.intp, types.Any)
@lower(stubs.atomic.inc, types.Array, types.UniTuple, types.Any)
@lower(stubs.atomic.inc, types.Array, types.Tuple, types.Any)
@_atomic_dispatcher
def ptx_atomic_inc(context, builder, dtype, ptr, val):
    if dtype in cuda.cudadecl.unsigned_int_numba_types:
        bw = dtype.bitwidth
        lmod = builder.module
        fn = getattr(nvvmutils, f'declare_atomic_inc_int{bw}')
        return builder.call(fn(lmod), (ptr, val))
    else:
        raise TypeError(f'Unimplemented atomic inc with {dtype} array')


@lower(stubs.atomic.dec, types.Array, types.intp, types.Any)
@lower(stubs.atomic.dec, types.Array, types.UniTuple, types.Any)
@lower(stubs.atomic.dec, types.Array, types.Tuple, types.Any)
@_atomic_dispatcher
def ptx_atomic_dec(context, builder, dtype, ptr, val):
    if dtype in cuda.cudadecl.unsigned_int_numba_types:
        bw = dtype.bitwidth
        lmod = builder.module
        fn = getattr(nvvmutils, f'declare_atomic_dec_int{bw}')
        return builder.call(fn(lmod), (ptr, val))
    else:
        raise TypeError(f'Unimplemented atomic dec with {dtype} array')


def ptx_atomic_bitwise(stub, op):
    @_atomic_dispatcher
    def impl_ptx_atomic(context, builder, dtype, ptr, val):
        if dtype in (cuda.cudadecl.integer_numba_types):
            return builder.atomic_rmw(op, ptr, val, 'monotonic')
        else:
            raise TypeError(f'Unimplemented atomic {op} with {dtype} array')

    for ty in (types.intp, types.UniTuple, types.Tuple):
        lower(stub, types.Array, ty, types.Any)(impl_ptx_atomic)


ptx_atomic_bitwise(stubs.atomic.and_, 'and')
ptx_atomic_bitwise(stubs.atomic.or_, 'or')
ptx_atomic_bitwise(stubs.atomic.xor, 'xor')


@lower(stubs.atomic.exch, types.Array, types.intp, types.Any)
@lower(stubs.atomic.exch, types.Array, types.UniTuple, types.Any)
@lower(stubs.atomic.exch, types.Array, types.Tuple, types.Any)
@_atomic_dispatcher
def ptx_atomic_exch(context, builder, dtype, ptr, val):
    if dtype in (cuda.cudadecl.integer_numba_types):
        return builder.atomic_rmw('xchg', ptr, val, 'monotonic')
    else:
        raise TypeError(f'Unimplemented atomic exch with {dtype} array')


@lower(stubs.atomic.max, types.Array, types.intp, types.Any)
@lower(stubs.atomic.max, types.Array, types.Tuple, types.Any)
@lower(stubs.atomic.max, types.Array, types.UniTuple, types.Any)
@_atomic_dispatcher
def ptx_atomic_max(context, builder, dtype, ptr, val):
    lmod = builder.module
    if dtype == types.float64:
        return builder.call(nvvmutils.declare_atomic_max_float64(lmod),
                            (ptr, val))
    elif dtype == types.float32:
        return builder.call(nvvmutils.declare_atomic_max_float32(lmod),
                            (ptr, val))
    elif dtype in (types.int32, types.int64):
        return builder.atomic_rmw('max', ptr, val, ordering='monotonic')
    elif dtype in (types.uint32, types.uint64):
        return builder.atomic_rmw('umax', ptr, val, ordering='monotonic')
    else:
        raise TypeError('Unimplemented atomic max with %s array' % dtype)


@lower(stubs.atomic.min, types.Array, types.intp, types.Any)
@lower(stubs.atomic.min, types.Array, types.Tuple, types.Any)
@lower(stubs.atomic.min, types.Array, types.UniTuple, types.Any)
@_atomic_dispatcher
def ptx_atomic_min(context, builder, dtype, ptr, val):
    lmod = builder.module
    if dtype == types.float64:
        return builder.call(nvvmutils.declare_atomic_min_float64(lmod),
                            (ptr, val))
    elif dtype == types.float32:
        return builder.call(nvvmutils.declare_atomic_min_float32(lmod),
                            (ptr, val))
    elif dtype in (types.int32, types.int64):
        return builder.atomic_rmw('min', ptr, val, ordering='monotonic')
    elif dtype in (types.uint32, types.uint64):
        return builder.atomic_rmw('umin', ptr, val, ordering='monotonic')
    else:
        raise TypeError('Unimplemented atomic min with %s array' % dtype)


@lower(stubs.atomic.nanmax, types.Array, types.intp, types.Any)
@lower(stubs.atomic.nanmax, types.Array, types.Tuple, types.Any)
@lower(stubs.atomic.nanmax, types.Array, types.UniTuple, types.Any)
@_atomic_dispatcher
def ptx_atomic_nanmax(context, builder, dtype, ptr, val):
    lmod = builder.module
    if dtype == types.float64:
        return builder.call(nvvmutils.declare_atomic_nanmax_float64(lmod),
                            (ptr, val))
    elif dtype == types.float32:
        return builder.call(nvvmutils.declare_atomic_nanmax_float32(lmod),
                            (ptr, val))
    elif dtype in (types.int32, types.int64):
        return builder.atomic_rmw('max', ptr, val, ordering='monotonic')
    elif dtype in (types.uint32, types.uint64):
        return builder.atomic_rmw('umax', ptr, val, ordering='monotonic')
    else:
        raise TypeError('Unimplemented atomic max with %s array' % dtype)


@lower(stubs.atomic.nanmin, types.Array, types.intp, types.Any)
@lower(stubs.atomic.nanmin, types.Array, types.Tuple, types.Any)
@lower(stubs.atomic.nanmin, types.Array, types.UniTuple, types.Any)
@_atomic_dispatcher
def ptx_atomic_nanmin(context, builder, dtype, ptr, val):
    lmod = builder.module
    if dtype == types.float64:
        return builder.call(nvvmutils.declare_atomic_nanmin_float64(lmod),
                            (ptr, val))
    elif dtype == types.float32:
        return builder.call(nvvmutils.declare_atomic_nanmin_float32(lmod),
                            (ptr, val))
    elif dtype in (types.int32, types.int64):
        return builder.atomic_rmw('min', ptr, val, ordering='monotonic')
    elif dtype in (types.uint32, types.uint64):
        return builder.atomic_rmw('umin', ptr, val, ordering='monotonic')
    else:
        raise TypeError('Unimplemented atomic min with %s array' % dtype)


@lower(stubs.atomic.compare_and_swap, types.Array, types.Any, types.Any)
def ptx_atomic_compare_and_swap(context, builder, sig, args):
    sig = sig.return_type(sig.args[0], types.intp, sig.args[1], sig.args[2])
    args = (args[0], context.get_constant(types.intp, 0), args[1], args[2])
    return ptx_atomic_cas(context, builder, sig, args)


@lower(stubs.atomic.cas, types.Array, types.intp, types.Any, types.Any)
@lower(stubs.atomic.cas, types.Array, types.Tuple, types.Any, types.Any)
@lower(stubs.atomic.cas, types.Array, types.UniTuple, types.Any, types.Any)
def ptx_atomic_cas(context, builder, sig, args):
    aryty, indty, oldty, valty = sig.args
    ary, inds, old, val = args

    indty, indices = _normalize_indices(context, builder, indty, inds, aryty,
                                        valty)

    lary = context.make_array(aryty)(context, builder, ary)
    ptr = cgutils.get_item_pointer(context, builder, aryty, lary, indices,
                                   wraparound=True)

    if aryty.dtype in (cuda.cudadecl.integer_numba_types):
        lmod = builder.module
        bitwidth = aryty.dtype.bitwidth
        return nvvmutils.atomic_cmpxchg(builder, lmod, bitwidth, ptr, old, val)
    else:
        raise TypeError('Unimplemented atomic cas with %s array' % aryty.dtype)


# -----------------------------------------------------------------------------

@lower(stubs.nanosleep, types.uint32)
def ptx_nanosleep(context, builder, sig, args):
    nanosleep = ir.InlineAsm(ir.FunctionType(ir.VoidType(), [ir.IntType(32)]),
                             "nanosleep.u32 $0;", 'r', side_effect=True)
    ns = args[0]
    builder.call(nanosleep, [ns])


# -----------------------------------------------------------------------------


def _generic_array(context, builder, shape, dtype, symbol_name, addrspace,
                   can_dynsized=False):
    elemcount = reduce(operator.mul, shape, 1)

    # Check for valid shape for this type of allocation.
    # Only 1d arrays can be dynamic.
    dynamic_smem = elemcount <= 0 and can_dynsized and len(shape) == 1
    if elemcount <= 0 and not dynamic_smem:
        raise ValueError("array length <= 0")

    # Check that we support the requested dtype
    data_model = context.data_model_manager[dtype]
    other_supported_type = (
        isinstance(dtype, (types.Record, types.Boolean))
        or isinstance(data_model, models.StructModel)
        or dtype == types.float16
    )
    if dtype not in types.number_domain and not other_supported_type:
        raise TypeError("unsupported type: %s" % dtype)

    lldtype = context.get_data_type(dtype)
    laryty = ir.ArrayType(lldtype, elemcount)

    if addrspace == nvvm.ADDRSPACE_LOCAL:
        # Special case local address space allocation to use alloca
        # NVVM is smart enough to only use local memory if no register is
        # available
        dataptr = cgutils.alloca_once(builder, laryty, name=symbol_name)
    else:
        lmod = builder.module

        # Create global variable in the requested address space
        gvmem = cgutils.add_global_variable(lmod, laryty, symbol_name,
                                            addrspace)
        # Specify alignment to avoid misalignment bug
        align = context.get_abi_sizeof(lldtype)
        # Alignment is required to be a power of 2 for shared memory. If it is
        # not a power of 2 (e.g. for a Record array) then round up accordingly.
        gvmem.align = 1 << (align - 1 ).bit_length()

        if dynamic_smem:
            gvmem.linkage = 'external'
        else:
            ## Comment out the following line to workaround a NVVM bug
            ## which generates a invalid symbol name when the linkage
            ## is internal and in some situation.
            ## See _get_unique_smem_id()
            # gvmem.linkage = lc.LINKAGE_INTERNAL

            gvmem.initializer = ir.Constant(laryty, ir.Undefined)

        # Convert to generic address-space
        dataptr = builder.addrspacecast(gvmem, ir.PointerType(ir.IntType(8)),
                                        'generic')

    targetdata = ll.create_target_data(nvvm.NVVM().data_layout)
    lldtype = context.get_data_type(dtype)
    itemsize = lldtype.get_abi_size(targetdata)

    # Compute strides
    laststride = itemsize
    rstrides = []
    for i, lastsize in enumerate(reversed(shape)):
        rstrides.append(laststride)
        laststride *= lastsize
    strides = [s for s in reversed(rstrides)]
    kstrides = [context.get_constant(types.intp, s) for s in strides]

    # Compute shape
    if dynamic_smem:
        # Compute the shape based on the dynamic shared memory configuration.
        # Unfortunately NVVM does not provide an intrinsic for the
        # %dynamic_smem_size register, so we must read it using inline
        # assembly.
        get_dynshared_size = ir.InlineAsm(ir.FunctionType(ir.IntType(32), []),
                                          "mov.u32 $0, %dynamic_smem_size;",
                                          '=r', side_effect=True)
        dynsmem_size = builder.zext(builder.call(get_dynshared_size, []),
                                    ir.IntType(64))
        # Only 1-D dynamic shared memory is supported so the following is a
        # sufficient construction of the shape
        kitemsize = context.get_constant(types.intp, itemsize)
        kshape = [builder.udiv(dynsmem_size, kitemsize)]
    else:
        kshape = [context.get_constant(types.intp, s) for s in shape]

    # Create array object
    ndim = len(shape)
    aryty = types.Array(dtype=dtype, ndim=ndim, layout='C')
    ary = context.make_array(aryty)(context, builder)

    context.populate_array(ary,
                           data=builder.bitcast(dataptr, ary.data.type),
                           shape=kshape,
                           strides=kstrides,
                           itemsize=context.get_constant(types.intp, itemsize),
                           meminfo=None)
    return ary._getvalue()


@lower_constant(CUDADispatcher)
def cuda_dispatcher_const(context, builder, ty, pyval):
    return context.get_dummy_value()


# NumPy

register_ufuncs(ufunc_db.get_ufuncs(), lower)
