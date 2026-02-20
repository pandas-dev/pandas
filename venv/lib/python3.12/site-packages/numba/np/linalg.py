"""
Implementation of linear algebra operations.
"""


import contextlib
import warnings

from llvmlite import ir

import numpy as np
import operator

from numba.core.imputils import (lower_builtin, impl_ret_borrowed,
                                    impl_ret_new_ref, impl_ret_untracked)
from numba.core.typing import signature
from numba.core.extending import intrinsic, overload, register_jitable
from numba.core import types, cgutils, config
from numba.core.errors import TypingError, NumbaTypeError, \
    NumbaPerformanceWarning
from .arrayobj import make_array, _empty_nd_impl, array_copy
from numba.np import numpy_support as np_support

ll_char = ir.IntType(8)
ll_char_p = ll_char.as_pointer()
ll_void_p = ll_char_p
ll_intc = ir.IntType(32)
ll_intc_p = ll_intc.as_pointer()
intp_t = cgutils.intp_t
ll_intp_p = intp_t.as_pointer()


# fortran int type, this needs to match the F_INT C declaration in
# _lapack.c and is present to accommodate potential future 64bit int
# based LAPACK use.
F_INT_nptype = np.int32
if config.USE_LEGACY_TYPE_SYSTEM:
    F_INT_nbtype = types.int32

    # BLAS kinds as letters
    _blas_kinds = {
        types.float32: 's',
        types.float64: 'd',
        types.complex64: 'c',
        types.complex128: 'z',
    }
else:
    F_INT_nbtype = types.np_int32

    # BLAS kinds as letters
    _blas_kinds = {
        types.np_float32: 's',
        types.np_float64: 'd',
        types.np_complex64: 'c',
        types.np_complex128: 'z',
    }


def get_blas_kind(dtype, func_name="<BLAS function>"):
    kind = _blas_kinds.get(dtype)
    if kind is None:
        raise NumbaTypeError("unsupported dtype for %s()" % (func_name,))
    return kind


def ensure_blas():
    try:
        import scipy.linalg.cython_blas
    except ImportError:
        raise ImportError("scipy 0.16+ is required for linear algebra")


def ensure_lapack():
    try:
        import scipy.linalg.cython_lapack
    except ImportError:
        raise ImportError("scipy 0.16+ is required for linear algebra")


def make_constant_slot(context, builder, ty, val):
    const = context.get_constant_generic(builder, ty, val)
    return cgutils.alloca_once_value(builder, const)


class _BLAS:
    """
    Functions to return type signatures for wrapped
    BLAS functions.
    """

    def __init__(self):
        ensure_blas()

    @classmethod
    def numba_xxnrm2(cls, dtype):
        rtype = getattr(dtype, "underlying_float", dtype)
        sig = types.intc(types.char,             # kind
                         types.intp,             # n
                         types.CPointer(dtype),  # x
                         types.intp,             # incx
                         types.CPointer(rtype))  # returned

        return types.ExternalFunction("numba_xxnrm2", sig)

    @classmethod
    def numba_xxgemm(cls, dtype):
        sig = types.intc(
            types.char,             # kind
            types.char,             # transa
            types.char,             # transb
            types.intp,             # m
            types.intp,             # n
            types.intp,             # k
            types.CPointer(dtype),  # alpha
            types.CPointer(dtype),  # a
            types.intp,             # lda
            types.CPointer(dtype),  # b
            types.intp,             # ldb
            types.CPointer(dtype),  # beta
            types.CPointer(dtype),  # c
            types.intp              # ldc
        )
        return types.ExternalFunction("numba_xxgemm", sig)


class _LAPACK:
    """
    Functions to return type signatures for wrapped
    LAPACK functions.
    """

    def __init__(self):
        ensure_lapack()

    @classmethod
    def numba_xxgetrf(cls, dtype):
        sig = types.intc(types.char,                   # kind
                         types.intp,                   # m
                         types.intp,                   # n
                         types.CPointer(dtype),        # a
                         types.intp,                   # lda
                         types.CPointer(F_INT_nbtype)  # ipiv
                         )
        return types.ExternalFunction("numba_xxgetrf", sig)

    @classmethod
    def numba_ez_xxgetri(cls, dtype):
        sig = types.intc(types.char,                   # kind
                         types.intp,                   # n
                         types.CPointer(dtype),        # a
                         types.intp,                   # lda
                         types.CPointer(F_INT_nbtype)  # ipiv
                         )
        return types.ExternalFunction("numba_ez_xxgetri", sig)

    @classmethod
    def numba_ez_rgeev(cls, dtype):
        sig = types.intc(types.char,             # kind
                         types.char,             # jobvl
                         types.char,             # jobvr
                         types.intp,             # n
                         types.CPointer(dtype),  # a
                         types.intp,             # lda
                         types.CPointer(dtype),  # wr
                         types.CPointer(dtype),  # wi
                         types.CPointer(dtype),  # vl
                         types.intp,             # ldvl
                         types.CPointer(dtype),  # vr
                         types.intp              # ldvr
                         )
        return types.ExternalFunction("numba_ez_rgeev", sig)

    @classmethod
    def numba_ez_cgeev(cls, dtype):
        sig = types.intc(types.char,             # kind
                         types.char,             # jobvl
                         types.char,             # jobvr
                         types.intp,             # n
                         types.CPointer(dtype),  # a
                         types.intp,             # lda
                         types.CPointer(dtype),  # w
                         types.CPointer(dtype),  # vl
                         types.intp,             # ldvl
                         types.CPointer(dtype),  # vr
                         types.intp              # ldvr
                         )
        return types.ExternalFunction("numba_ez_cgeev", sig)

    @classmethod
    def numba_ez_xxxevd(cls, dtype):
        wtype = getattr(dtype, "underlying_float", dtype)
        sig = types.intc(types.char,             # kind
                         types.char,             # jobz
                         types.char,             # uplo
                         types.intp,             # n
                         types.CPointer(dtype),  # a
                         types.intp,             # lda
                         types.CPointer(wtype),  # w
                         )
        return types.ExternalFunction("numba_ez_xxxevd", sig)

    @classmethod
    def numba_xxpotrf(cls, dtype):
        sig = types.intc(types.char,             # kind
                         types.char,             # uplo
                         types.intp,             # n
                         types.CPointer(dtype),  # a
                         types.intp              # lda
                         )
        return types.ExternalFunction("numba_xxpotrf", sig)

    @classmethod
    def numba_ez_gesdd(cls, dtype):
        stype = getattr(dtype, "underlying_float", dtype)
        sig = types.intc(
            types.char,             # kind
            types.char,             # jobz
            types.intp,             # m
            types.intp,             # n
            types.CPointer(dtype),  # a
            types.intp,             # lda
            types.CPointer(stype),  # s
            types.CPointer(dtype),  # u
            types.intp,             # ldu
            types.CPointer(dtype),  # vt
            types.intp              # ldvt
        )

        return types.ExternalFunction("numba_ez_gesdd", sig)

    @classmethod
    def numba_ez_geqrf(cls, dtype):
        sig = types.intc(
            types.char,             # kind
            types.intp,             # m
            types.intp,             # n
            types.CPointer(dtype),  # a
            types.intp,             # lda
            types.CPointer(dtype),  # tau
        )
        return types.ExternalFunction("numba_ez_geqrf", sig)

    @classmethod
    def numba_ez_xxgqr(cls, dtype):
        sig = types.intc(
            types.char,             # kind
            types.intp,             # m
            types.intp,             # n
            types.intp,             # k
            types.CPointer(dtype),  # a
            types.intp,             # lda
            types.CPointer(dtype),  # tau
        )
        return types.ExternalFunction("numba_ez_xxgqr", sig)

    @classmethod
    def numba_ez_gelsd(cls, dtype):
        rtype = getattr(dtype, "underlying_float", dtype)
        sig = types.intc(
            types.char,                 # kind
            types.intp,                 # m
            types.intp,                 # n
            types.intp,                 # nrhs
            types.CPointer(dtype),      # a
            types.intp,                 # lda
            types.CPointer(dtype),      # b
            types.intp,                 # ldb
            types.CPointer(rtype),      # S
            types.float64,              # rcond
            types.CPointer(types.intc)  # rank
        )
        return types.ExternalFunction("numba_ez_gelsd", sig)

    @classmethod
    def numba_xgesv(cls, dtype):
        sig = types.intc(
            types.char,                    # kind
            types.intp,                    # n
            types.intp,                    # nhrs
            types.CPointer(dtype),         # a
            types.intp,                    # lda
            types.CPointer(F_INT_nbtype),  # ipiv
            types.CPointer(dtype),         # b
            types.intp                     # ldb
        )
        return types.ExternalFunction("numba_xgesv", sig)


@contextlib.contextmanager
def make_contiguous(context, builder, sig, args):
    """
    Ensure that all array arguments are contiguous, if necessary by
    copying them.
    A new (sig, args) tuple is yielded.
    """
    newtys = []
    newargs = []
    copies = []
    for ty, val in zip(sig.args, args):
        if not isinstance(ty, types.Array) or ty.layout in 'CF':
            newty, newval = ty, val
        else:
            newty = ty.copy(layout='C')
            copysig = signature(newty, ty)
            newval = array_copy(context, builder, copysig, (val,))
            copies.append((newty, newval))
        newtys.append(newty)
        newargs.append(newval)
    yield signature(sig.return_type, *newtys), tuple(newargs)
    for ty, val in copies:
        context.nrt.decref(builder, ty, val)


def check_c_int(context, builder, n):
    """
    Check whether *n* fits in a C `int`.
    """
    _maxint = 2**31 - 1

    def impl(n):
        if n > _maxint:
            raise OverflowError("array size too large to fit in C int")

    context.compile_internal(builder, impl,
                             signature(types.none, types.intp), (n,))


def check_blas_return(context, builder, res):
    """
    Check the integer error return from one of the BLAS wrappers in
    _helperlib.c.
    """
    with builder.if_then(cgutils.is_not_null(builder, res), likely=False):
        # Those errors shouldn't happen, it's easier to just abort the process
        pyapi = context.get_python_api(builder)
        pyapi.gil_ensure()
        pyapi.fatal_error("BLAS wrapper returned with an error")


def check_lapack_return(context, builder, res):
    """
    Check the integer error return from one of the LAPACK wrappers in
    _helperlib.c.
    """
    with builder.if_then(cgutils.is_not_null(builder, res), likely=False):
        # Those errors shouldn't happen, it's easier to just abort the process
        pyapi = context.get_python_api(builder)
        pyapi.gil_ensure()
        pyapi.fatal_error("LAPACK wrapper returned with an error")


def call_xxdot(context, builder, conjugate, dtype,
               n, a_data, b_data, out_data):
    """
    Call the BLAS vector * vector product function for the given arguments.
    """
    fnty = ir.FunctionType(ir.IntType(32),
                           [ll_char, ll_char, intp_t,    # kind, conjugate, n
                            ll_void_p, ll_void_p, ll_void_p,  # a, b, out
                            ])
    fn = cgutils.get_or_insert_function(builder.module, fnty, "numba_xxdot")

    kind = get_blas_kind(dtype)
    kind_val = ir.Constant(ll_char, ord(kind))
    conjugate = ir.Constant(ll_char, int(conjugate))

    res = builder.call(fn, (kind_val, conjugate, n,
                            builder.bitcast(a_data, ll_void_p),
                            builder.bitcast(b_data, ll_void_p),
                            builder.bitcast(out_data, ll_void_p)))
    check_blas_return(context, builder, res)


def call_xxgemv(context, builder, do_trans,
                m_type, m_shapes, m_data, v_data, out_data):
    """
    Call the BLAS matrix * vector product function for the given arguments.
    """
    fnty = ir.FunctionType(ir.IntType(32),
                           [ll_char, ll_char,                 # kind, trans
                            intp_t, intp_t,                   # m, n
                            ll_void_p, ll_void_p, intp_t,     # alpha, a, lda
                            ll_void_p, ll_void_p, ll_void_p,  # x, beta, y
                            ])
    fn = cgutils.get_or_insert_function(builder.module, fnty, "numba_xxgemv")

    dtype = m_type.dtype
    alpha = make_constant_slot(context, builder, dtype, 1.0)
    beta = make_constant_slot(context, builder, dtype, 0.0)

    if m_type.layout == 'F':
        m, n = m_shapes
        lda = m_shapes[0]
    else:
        n, m = m_shapes
        lda = m_shapes[1]

    kind = get_blas_kind(dtype)
    kind_val = ir.Constant(ll_char, ord(kind))
    trans = ir.Constant(ll_char, ord('t') if do_trans else ord('n'))

    res = builder.call(fn, (kind_val, trans, m, n,
                            builder.bitcast(alpha, ll_void_p),
                            builder.bitcast(m_data, ll_void_p), lda,
                            builder.bitcast(v_data, ll_void_p),
                            builder.bitcast(beta, ll_void_p),
                            builder.bitcast(out_data, ll_void_p)))
    check_blas_return(context, builder, res)


def call_xxgemm(context, builder,
                x_type, x_shapes, x_data,
                y_type, y_shapes, y_data,
                out_type, out_shapes, out_data):
    """
    Call the BLAS matrix * matrix product function for the given arguments.
    """
    fnty = ir.FunctionType(ir.IntType(32),
                           [ll_char,                       # kind
                            ll_char, ll_char,              # transa, transb
                            intp_t, intp_t, intp_t,        # m, n, k
                            ll_void_p, ll_void_p, intp_t,  # alpha, a, lda
                            ll_void_p, intp_t, ll_void_p,  # b, ldb, beta
                            ll_void_p, intp_t,             # c, ldc
                            ])
    fn = cgutils.get_or_insert_function(builder.module, fnty, "numba_xxgemm")

    m, k = x_shapes
    _k, n = y_shapes
    dtype = x_type.dtype
    alpha = make_constant_slot(context, builder, dtype, 1.0)
    beta = make_constant_slot(context, builder, dtype, 0.0)

    trans = ir.Constant(ll_char, ord('t'))
    notrans = ir.Constant(ll_char, ord('n'))

    def get_array_param(ty, shapes, data):
        return (
            # Transpose if layout different from result's
            notrans if ty.layout == out_type.layout else trans,
            # Size of the inner dimension in physical array order
            shapes[1] if ty.layout == 'C' else shapes[0],
            # The data pointer, unit-less
            builder.bitcast(data, ll_void_p),
        )

    transa, lda, data_a = get_array_param(y_type, y_shapes, y_data)
    transb, ldb, data_b = get_array_param(x_type, x_shapes, x_data)
    _, ldc, data_c = get_array_param(out_type, out_shapes, out_data)

    kind = get_blas_kind(dtype)
    kind_val = ir.Constant(ll_char, ord(kind))

    res = builder.call(fn, (kind_val, transa, transb, n, m, k,
                            builder.bitcast(alpha, ll_void_p), data_a, lda,
                            data_b, ldb, builder.bitcast(beta, ll_void_p),
                            data_c, ldc))
    check_blas_return(context, builder, res)


def dot_2_mm(context, builder, sig, args):
    """
    np.dot(matrix, matrix)
    """
    def dot_impl(a, b):
        m, k = a.shape
        _k, n = b.shape
        if k == 0:
            return np.zeros((m, n), a.dtype)
        out = np.empty((m, n), a.dtype)
        return np.dot(a, b, out)

    res = context.compile_internal(builder, dot_impl, sig, args)
    return impl_ret_new_ref(context, builder, sig.return_type, res)


def dot_2_vm(context, builder, sig, args):
    """
    np.dot(vector, matrix)
    """
    def dot_impl(a, b):
        m, = a.shape
        _m, n = b.shape
        if m == 0:
            return np.zeros((n, ), a.dtype)
        out = np.empty((n, ), a.dtype)
        return np.dot(a, b, out)

    res = context.compile_internal(builder, dot_impl, sig, args)
    return impl_ret_new_ref(context, builder, sig.return_type, res)


def dot_2_mv(context, builder, sig, args):
    """
    np.dot(matrix, vector)
    """
    def dot_impl(a, b):
        m, n = a.shape
        _n, = b.shape
        if n == 0:
            return np.zeros((m, ), a.dtype)
        out = np.empty((m, ), a.dtype)
        return np.dot(a, b, out)

    res = context.compile_internal(builder, dot_impl, sig, args)
    return impl_ret_new_ref(context, builder, sig.return_type, res)


def dot_2_vv(context, builder, sig, args, conjugate=False):
    """
    np.dot(vector, vector)
    np.vdot(vector, vector)
    """
    aty, bty = sig.args
    dtype = sig.return_type
    a = make_array(aty)(context, builder, args[0])
    b = make_array(bty)(context, builder, args[1])
    n, = cgutils.unpack_tuple(builder, a.shape)

    def check_args(a, b):
        m, = a.shape
        n, = b.shape
        if m != n:
            raise ValueError("incompatible array sizes for np.dot(a, b) "
                             "(vector * vector)")

    context.compile_internal(builder, check_args,
                             signature(types.none, *sig.args), args)
    check_c_int(context, builder, n)

    out = cgutils.alloca_once(builder, context.get_value_type(dtype))
    call_xxdot(context, builder, conjugate, dtype, n, a.data, b.data, out)
    return builder.load(out)


@overload(np.dot)
def dot_2(left, right):
    """
    np.dot(a, b)
    """
    return dot_2_impl('np.dot()', left, right)


@overload(operator.matmul)
def matmul_2(left, right):
    """
    a @ b
    """
    return dot_2_impl("'@'", left, right)


def dot_2_impl(name, left, right):
    if isinstance(left, types.Array) and isinstance(right, types.Array):
        @intrinsic
        def _impl(typingcontext, left, right):
            ndims = (left.ndim, right.ndim)

            def _dot2_codegen(context, builder, sig, args):
                ensure_blas()

                with make_contiguous(context, builder, sig, args) as (sig, args):
                    if ndims == (2, 2):
                        return dot_2_mm(context, builder, sig, args)
                    elif ndims == (2, 1):
                        return dot_2_mv(context, builder, sig, args)
                    elif ndims == (1, 2):
                        return dot_2_vm(context, builder, sig, args)
                    elif ndims == (1, 1):
                        return dot_2_vv(context, builder, sig, args)
                    else:
                        raise AssertionError('unreachable')

            if left.dtype != right.dtype:
                raise TypingError(
                    "%s arguments must all have the same dtype" % name)

            if ndims == (2, 2):
                return_type = types.Array(left.dtype, 2, 'C')
            elif ndims == (2, 1) or ndims == (1, 2):
                return_type = types.Array(left.dtype, 1, 'C')
            elif ndims == (1, 1):
                return_type = left.dtype
            else:
                raise TypingError(("%s: inputs must have compatible "
                                   "dimensions") % name)
            return signature(return_type, left, right), _dot2_codegen

        if left.layout not in 'CF' or right.layout not in 'CF':
            warnings.warn(
                "%s is faster on contiguous arrays, called on %s" % (
                    name, (left, right),), NumbaPerformanceWarning)

        return lambda left, right: _impl(left, right)


@overload(np.vdot)
def vdot(left, right):
    """
    np.vdot(a, b)
    """
    if isinstance(left, types.Array) and isinstance(right, types.Array):
        @intrinsic
        def _impl(typingcontext, left, right):
            def codegen(context, builder, sig, args):
                ensure_blas()

                with make_contiguous(context, builder, sig, args) as\
                        (sig, args):
                    return dot_2_vv(context, builder, sig, args, conjugate=True)

            if left.ndim != 1 or right.ndim != 1:
                raise TypingError("np.vdot() only supported on 1-D arrays")

            if left.dtype != right.dtype:
                raise TypingError(
                    "np.vdot() arguments must all have the same dtype")
            return signature(left.dtype, left, right), codegen

        if left.layout not in 'CF' or right.layout not in 'CF':
            warnings.warn(
                "np.vdot() is faster on contiguous arrays, called on %s"
                % ((left, right),), NumbaPerformanceWarning)

        return lambda left, right: _impl(left, right)


def dot_3_vm_check_args(a, b, out):
    m, = a.shape
    _m, n = b.shape
    if m != _m:
        raise ValueError("incompatible array sizes for "
                         "np.dot(a, b) (vector * matrix)")
    if out.shape != (n,):
        raise ValueError("incompatible output array size for "
                         "np.dot(a, b, out) (vector * matrix)")


def dot_3_mv_check_args(a, b, out):
    m, _n = a.shape
    n, = b.shape
    if n != _n:
        raise ValueError("incompatible array sizes for np.dot(a, b) "
                         "(matrix * vector)")
    if out.shape != (m,):
        raise ValueError("incompatible output array size for "
                         "np.dot(a, b, out) (matrix * vector)")


def dot_3_vm(context, builder, sig, args):
    """
    np.dot(vector, matrix, out)
    np.dot(matrix, vector, out)
    """
    xty, yty, outty = sig.args
    assert outty == sig.return_type
    dtype = xty.dtype

    x = make_array(xty)(context, builder, args[0])
    y = make_array(yty)(context, builder, args[1])
    out = make_array(outty)(context, builder, args[2])
    x_shapes = cgutils.unpack_tuple(builder, x.shape)
    y_shapes = cgutils.unpack_tuple(builder, y.shape)
    out_shapes = cgutils.unpack_tuple(builder, out.shape)
    if xty.ndim < yty.ndim:
        # Vector * matrix
        # Asked for x * y, we will compute y.T * x
        mty = yty
        m_shapes = y_shapes
        v_shape = x_shapes[0]
        lda = m_shapes[1]
        do_trans = yty.layout == 'F'
        m_data, v_data = y.data, x.data
        check_args = dot_3_vm_check_args
    else:
        # Matrix * vector
        # We will compute x * y
        mty = xty
        m_shapes = x_shapes
        v_shape = y_shapes[0]
        lda = m_shapes[0]
        do_trans = xty.layout == 'C'
        m_data, v_data = x.data, y.data
        check_args = dot_3_mv_check_args

    context.compile_internal(builder, check_args,
                             signature(types.none, *sig.args), args)
    for val in m_shapes:
        check_c_int(context, builder, val)

    zero = context.get_constant(types.intp, 0)
    both_empty = builder.icmp_signed('==', v_shape, zero)
    matrix_empty = builder.icmp_signed('==', lda, zero)
    is_empty = builder.or_(both_empty, matrix_empty)
    with builder.if_else(is_empty, likely=False) as (empty, nonempty):
        with empty:
            cgutils.memset(builder, out.data,
                           builder.mul(out.itemsize, out.nitems), 0)
        with nonempty:
            call_xxgemv(context, builder, do_trans, mty, m_shapes, m_data,
                        v_data, out.data)

    return impl_ret_borrowed(context, builder, sig.return_type,
                             out._getvalue())


def dot_3_mm(context, builder, sig, args):
    """
    np.dot(matrix, matrix, out)
    """
    xty, yty, outty = sig.args
    assert outty == sig.return_type
    dtype = xty.dtype

    x = make_array(xty)(context, builder, args[0])
    y = make_array(yty)(context, builder, args[1])
    out = make_array(outty)(context, builder, args[2])
    x_shapes = cgutils.unpack_tuple(builder, x.shape)
    y_shapes = cgutils.unpack_tuple(builder, y.shape)
    out_shapes = cgutils.unpack_tuple(builder, out.shape)
    m, k = x_shapes
    _k, n = y_shapes

    # The only case Numpy supports
    assert outty.layout == 'C'

    def check_args(a, b, out):
        m, k = a.shape
        _k, n = b.shape
        if k != _k:
            raise ValueError("incompatible array sizes for np.dot(a, b) "
                             "(matrix * matrix)")
        if out.shape != (m, n):
            raise ValueError("incompatible output array size for "
                             "np.dot(a, b, out) (matrix * matrix)")

    context.compile_internal(builder, check_args,
                             signature(types.none, *sig.args), args)

    check_c_int(context, builder, m)
    check_c_int(context, builder, k)
    check_c_int(context, builder, n)

    x_data = x.data
    y_data = y.data
    out_data = out.data

    # If eliminated dimension is zero, set all entries to zero and return
    zero = context.get_constant(types.intp, 0)
    both_empty = builder.icmp_signed('==', k, zero)
    x_empty = builder.icmp_signed('==', m, zero)
    y_empty = builder.icmp_signed('==', n, zero)
    is_empty = builder.or_(both_empty, builder.or_(x_empty, y_empty))
    with builder.if_else(is_empty, likely=False) as (empty, nonempty):
        with empty:
            cgutils.memset(builder, out.data,
                           builder.mul(out.itemsize, out.nitems), 0)
        with nonempty:
            # Check if any of the operands is really a 1-d vector represented
            # as a (1, k) or (k, 1) 2-d array.  In those cases, it is pessimal
            # to call the generic matrix * matrix product BLAS function.
            one = context.get_constant(types.intp, 1)
            is_left_vec = builder.icmp_signed('==', m, one)
            is_right_vec = builder.icmp_signed('==', n, one)

            with builder.if_else(is_right_vec) as (r_vec, r_mat):
                with r_vec:
                    with builder.if_else(is_left_vec) as (v_v, m_v):
                        with v_v:
                            # V * V
                            call_xxdot(context, builder, False, dtype,
                                       k, x_data, y_data, out_data)
                        with m_v:
                            # M * V
                            do_trans = xty.layout == outty.layout
                            call_xxgemv(context, builder, do_trans,
                                        xty, x_shapes, x_data, y_data, out_data)
                with r_mat:
                    with builder.if_else(is_left_vec) as (v_m, m_m):
                        with v_m:
                            # V * M
                            do_trans = yty.layout != outty.layout
                            call_xxgemv(context, builder, do_trans,
                                        yty, y_shapes, y_data, x_data, out_data)
                        with m_m:
                            # M * M
                            call_xxgemm(context, builder,
                                        xty, x_shapes, x_data,
                                        yty, y_shapes, y_data,
                                        outty, out_shapes, out_data)

    return impl_ret_borrowed(context, builder, sig.return_type,
                             out._getvalue())


@overload(np.dot)
def dot_3(left, right, out):
    """
    np.dot(a, b, out)
    """
    if (isinstance(left, types.Array) and isinstance(right, types.Array) and
            isinstance(out, types.Array)):
        @intrinsic
        def _impl(typingcontext, left, right, out):
            def codegen(context, builder, sig, args):
                ensure_blas()

                with make_contiguous(context, builder, sig, args) as (sig,
                                                                      args):
                    ndims = set(x.ndim for x in sig.args[:2])
                    if ndims == {2}:
                        return dot_3_mm(context, builder, sig, args)
                    elif ndims == {1, 2}:
                        return dot_3_vm(context, builder, sig, args)
                    else:
                        raise AssertionError('unreachable')
            if left.dtype != right.dtype or left.dtype != out.dtype:
                raise TypingError(
                    "np.dot() arguments must all have the same dtype")

            return signature(out, left, right, out), codegen

        if left.layout not in 'CF' or right.layout not in 'CF' or out.layout\
            not in 'CF':
            warnings.warn(
                "np.vdot() is faster on contiguous arrays, called on %s"
                % ((left, right),), NumbaPerformanceWarning)

        return lambda left, right, out: _impl(left, right, out)


if config.USE_LEGACY_TYPE_SYSTEM:
    fatal_error_func = types.ExternalFunction("numba_fatal_error", types.intc())
else:
    fatal_error_func = types.ExternalFunction("numba_fatal_error", types.c_intp())


@register_jitable
def _check_finite_matrix(a):
    for v in np.nditer(a):
        if not np.isfinite(v.item()):
            raise np.linalg.LinAlgError(
                "Array must not contain infs or NaNs.")


def _check_linalg_matrix(a, func_name, la_prefix=True):
    # la_prefix is present as some functions, e.g. np.trace()
    # are documented under "linear algebra" but aren't in the
    # module
    prefix = "np.linalg" if la_prefix else "np"
    interp = (prefix, func_name)
    # Unpack optional type
    if isinstance(a, types.Optional):
        a = a.type
    if not isinstance(a, types.Array):
        msg = "%s.%s() only supported for array types" % interp
        raise TypingError(msg, highlighting=False)
    if not a.ndim == 2:
        msg = "%s.%s() only supported on 2-D arrays." % interp
        raise TypingError(msg, highlighting=False)
    if not isinstance(a.dtype, (types.Float, types.Complex)):
        msg = "%s.%s() only supported on "\
            "float and complex arrays." % interp
        raise TypingError(msg, highlighting=False)


def _check_homogeneous_types(func_name, *types):
    t0 = types[0].dtype
    for t in types[1:]:
        if t.dtype != t0:
            msg = "np.linalg.%s() only supports inputs that have homogeneous dtypes." % func_name
            raise TypingError(msg, highlighting=False)


def _copy_to_fortran_order():
    pass


@overload(_copy_to_fortran_order)
def ol_copy_to_fortran_order(a):
    # This function copies the array 'a' into a new array with fortran order.
    # This exists because the copy routines don't take order flags yet.
    F_layout = a.layout == 'F'
    A_layout = a.layout == 'A'
    def impl(a):
        if F_layout:
            # it's F ordered at compile time, just copy
            acpy = np.copy(a)
        elif A_layout:
            # decide based on runtime value
            flag_f = a.flags.f_contiguous
            if flag_f:
                # it's already F ordered, so copy but in a round about way to
                # ensure that the copy is also F ordered
                acpy = np.copy(a.T).T
            else:
                # it's something else ordered, so let asfortranarray deal with
                # copying and making it fortran ordered
                acpy = np.asfortranarray(a)
        else:
            # it's C ordered at compile time, asfortranarray it.
            acpy = np.asfortranarray(a)
        return acpy
    return impl


@register_jitable
def _inv_err_handler(r):
    if r != 0:
        if r < 0:
            fatal_error_func()
            assert 0   # unreachable
        if r > 0:
            raise np.linalg.LinAlgError(
                "Matrix is singular to machine precision.")

@register_jitable
def _dummy_liveness_func(a):
    """pass a list of variables to be preserved through dead code elimination"""
    return a[0]


@overload(np.linalg.inv)
def inv_impl(a):
    ensure_lapack()

    _check_linalg_matrix(a, "inv")

    numba_xxgetrf = _LAPACK().numba_xxgetrf(a.dtype)

    numba_xxgetri = _LAPACK().numba_ez_xxgetri(a.dtype)

    kind = ord(get_blas_kind(a.dtype, "inv"))

    def inv_impl(a):
        n = a.shape[-1]
        if a.shape[-2] != n:
            msg = "Last 2 dimensions of the array must be square."
            raise np.linalg.LinAlgError(msg)

        _check_finite_matrix(a)

        acpy = _copy_to_fortran_order(a)

        if n == 0:
            return acpy

        ipiv = np.empty(n, dtype=F_INT_nptype)

        r = numba_xxgetrf(kind, n, n, acpy.ctypes, n, ipiv.ctypes)
        _inv_err_handler(r)

        r = numba_xxgetri(kind, n, acpy.ctypes, n, ipiv.ctypes)
        _inv_err_handler(r)

        # help liveness analysis
        _dummy_liveness_func([acpy.size, ipiv.size])
        return acpy

    return inv_impl


@register_jitable
def _handle_err_maybe_convergence_problem(r):
    if r != 0:
        if r < 0:
            fatal_error_func()
            assert 0   # unreachable
        if r > 0:
            raise ValueError("Internal algorithm failed to converge.")


def _check_linalg_1_or_2d_matrix(a, func_name, la_prefix=True):
    # la_prefix is present as some functions, e.g. np.trace()
    # are documented under "linear algebra" but aren't in the
    # module
    prefix = "np.linalg" if la_prefix else "np"
    interp = (prefix, func_name)
    # checks that a matrix is 1 or 2D
    if not isinstance(a, types.Array):
        raise TypingError("%s.%s() only supported for array types "
                          % interp)
    if not a.ndim <= 2:
        raise TypingError("%s.%s() only supported on 1 and 2-D arrays "
                          % interp)
    if not isinstance(a.dtype, (types.Float, types.Complex)):
        raise TypingError("%s.%s() only supported on "
                          "float and complex arrays." % interp)


@overload(np.linalg.cholesky)
def cho_impl(a):
    ensure_lapack()

    _check_linalg_matrix(a, "cholesky")

    numba_xxpotrf = _LAPACK().numba_xxpotrf(a.dtype)

    kind = ord(get_blas_kind(a.dtype, "cholesky"))
    UP = ord('U')
    LO = ord('L')

    def cho_impl(a):
        n = a.shape[-1]
        if a.shape[-2] != n:
            msg = "Last 2 dimensions of the array must be square."
            raise np.linalg.LinAlgError(msg)

        # The output is allocated in C order
        out = a.copy()

        if n == 0:
            return out

        # Pass UP since xxpotrf() operates in F order
        # The semantics ensure this works fine
        # (out is really its Hermitian in F order, but UP instructs
        #  xxpotrf to compute the Hermitian of the upper triangle
        #  => they cancel each other)
        r = numba_xxpotrf(kind, UP, n, out.ctypes, n)
        if r != 0:
            if r < 0:
                fatal_error_func()
                assert 0   # unreachable
            if r > 0:
                raise np.linalg.LinAlgError(
                    "Matrix is not positive definite.")
        # Zero out upper triangle, in F order
        for col in range(n):
            out[:col, col] = 0
        return out

    return cho_impl

@overload(np.linalg.eig)
def eig_impl(a):
    ensure_lapack()

    _check_linalg_matrix(a, "eig")

    numba_ez_rgeev = _LAPACK().numba_ez_rgeev(a.dtype)
    numba_ez_cgeev = _LAPACK().numba_ez_cgeev(a.dtype)

    kind = ord(get_blas_kind(a.dtype, "eig"))

    JOBVL = ord('N')
    JOBVR = ord('V')

    def real_eig_impl(a):
        """
        eig() implementation for real arrays.
        """
        n = a.shape[-1]
        if a.shape[-2] != n:
            msg = "Last 2 dimensions of the array must be square."
            raise np.linalg.LinAlgError(msg)

        _check_finite_matrix(a)

        acpy = _copy_to_fortran_order(a)

        ldvl = 1
        ldvr = n
        wr = np.empty(n, dtype=a.dtype)
        wi = np.empty(n, dtype=a.dtype)
        vl = np.empty((n, ldvl), dtype=a.dtype)
        vr = np.empty((n, ldvr), dtype=a.dtype)

        if n == 0:
            return (wr, vr.T)

        r = numba_ez_rgeev(kind,
                            JOBVL,
                            JOBVR,
                            n,
                            acpy.ctypes,
                            n,
                            wr.ctypes,
                            wi.ctypes,
                            vl.ctypes,
                            ldvl,
                            vr.ctypes,
                            ldvr)
        _handle_err_maybe_convergence_problem(r)

        # By design numba does not support dynamic return types, however,
        # Numpy does. Numpy uses this ability in the case of returning
        # eigenvalues/vectors of a real matrix. The return type of
        # np.linalg.eig(), when operating on a matrix in real space
        # depends on the values present in the matrix itself (recalling
        # that eigenvalues are the roots of the characteristic polynomial
        # of the system matrix, which will by construction depend on the
        # values present in the system matrix). As numba cannot handle
        # the case of a runtime decision based domain change relative to
        # the input type, if it is required numba raises as below.
        if np.any(wi):
            raise ValueError(
                "eig() argument must not cause a domain change.")

        # put these in to help with liveness analysis,
        # `.ctypes` doesn't keep the vars alive
        _dummy_liveness_func([acpy.size, vl.size, vr.size, wr.size, wi.size])
        return (wr, vr.T)

    def cmplx_eig_impl(a):
        """
        eig() implementation for complex arrays.
        """
        n = a.shape[-1]
        if a.shape[-2] != n:
            msg = "Last 2 dimensions of the array must be square."
            raise np.linalg.LinAlgError(msg)

        _check_finite_matrix(a)

        acpy = _copy_to_fortran_order(a)

        ldvl = 1
        ldvr = n
        w = np.empty(n, dtype=a.dtype)
        vl = np.empty((n, ldvl), dtype=a.dtype)
        vr = np.empty((n, ldvr), dtype=a.dtype)

        if n == 0:
            return (w, vr.T)

        r = numba_ez_cgeev(kind,
                            JOBVL,
                            JOBVR,
                            n,
                            acpy.ctypes,
                            n,
                            w.ctypes,
                            vl.ctypes,
                            ldvl,
                            vr.ctypes,
                            ldvr)
        _handle_err_maybe_convergence_problem(r)

        # put these in to help with liveness analysis,
        # `.ctypes` doesn't keep the vars alive
        _dummy_liveness_func([acpy.size, vl.size, vr.size, w.size])
        return (w, vr.T)

    if isinstance(a.dtype, types.scalars.Complex):
        return cmplx_eig_impl
    else:
        return real_eig_impl

@overload(np.linalg.eigvals)
def eigvals_impl(a):
    ensure_lapack()

    _check_linalg_matrix(a, "eigvals")

    numba_ez_rgeev = _LAPACK().numba_ez_rgeev(a.dtype)
    numba_ez_cgeev = _LAPACK().numba_ez_cgeev(a.dtype)

    kind = ord(get_blas_kind(a.dtype, "eigvals"))

    JOBVL = ord('N')
    JOBVR = ord('N')

    def real_eigvals_impl(a):
        """
        eigvals() implementation for real arrays.
        """
        n = a.shape[-1]
        if a.shape[-2] != n:
            msg = "Last 2 dimensions of the array must be square."
            raise np.linalg.LinAlgError(msg)

        _check_finite_matrix(a)

        acpy = _copy_to_fortran_order(a)

        ldvl = 1
        ldvr = 1
        wr = np.empty(n, dtype=a.dtype)

        if n == 0:
            return wr

        wi = np.empty(n, dtype=a.dtype)

        # not referenced but need setting for MKL null check
        vl = np.empty((1), dtype=a.dtype)
        vr = np.empty((1), dtype=a.dtype)

        r = numba_ez_rgeev(kind,
                            JOBVL,
                            JOBVR,
                            n,
                            acpy.ctypes,
                            n,
                            wr.ctypes,
                            wi.ctypes,
                            vl.ctypes,
                            ldvl,
                            vr.ctypes,
                            ldvr)
        _handle_err_maybe_convergence_problem(r)

        # By design numba does not support dynamic return types, however,
        # Numpy does. Numpy uses this ability in the case of returning
        # eigenvalues/vectors of a real matrix. The return type of
        # np.linalg.eigvals(), when operating on a matrix in real space
        # depends on the values present in the matrix itself (recalling
        # that eigenvalues are the roots of the characteristic polynomial
        # of the system matrix, which will by construction depend on the
        # values present in the system matrix). As numba cannot handle
        # the case of a runtime decision based domain change relative to
        # the input type, if it is required numba raises as below.
        if np.any(wi):
            raise ValueError(
                "eigvals() argument must not cause a domain change.")

        # put these in to help with liveness analysis,
        # `.ctypes` doesn't keep the vars alive
        _dummy_liveness_func([acpy.size, vl.size, vr.size, wr.size, wi.size])
        return wr

    def cmplx_eigvals_impl(a):
        """
        eigvals() implementation for complex arrays.
        """
        n = a.shape[-1]
        if a.shape[-2] != n:
            msg = "Last 2 dimensions of the array must be square."
            raise np.linalg.LinAlgError(msg)

        _check_finite_matrix(a)

        acpy = _copy_to_fortran_order(a)

        ldvl = 1
        ldvr = 1
        w = np.empty(n, dtype=a.dtype)

        if n == 0:
            return w

        vl = np.empty((1), dtype=a.dtype)
        vr = np.empty((1), dtype=a.dtype)

        r = numba_ez_cgeev(kind,
                            JOBVL,
                            JOBVR,
                            n,
                            acpy.ctypes,
                            n,
                            w.ctypes,
                            vl.ctypes,
                            ldvl,
                            vr.ctypes,
                            ldvr)
        _handle_err_maybe_convergence_problem(r)

        # put these in to help with liveness analysis,
        # `.ctypes` doesn't keep the vars alive
        _dummy_liveness_func([acpy.size, vl.size, vr.size, w.size])
        return w

    if isinstance(a.dtype, types.scalars.Complex):
        return cmplx_eigvals_impl
    else:
        return real_eigvals_impl

@overload(np.linalg.eigh)
def eigh_impl(a):
    ensure_lapack()

    _check_linalg_matrix(a, "eigh")

    # convert typing floats to numpy floats for use in the impl
    w_type = getattr(a.dtype, "underlying_float", a.dtype)
    w_dtype = np_support.as_dtype(w_type)

    numba_ez_xxxevd = _LAPACK().numba_ez_xxxevd(a.dtype)

    kind = ord(get_blas_kind(a.dtype, "eigh"))

    JOBZ = ord('V')
    UPLO = ord('L')

    def eigh_impl(a):
        n = a.shape[-1]

        if a.shape[-2] != n:
            msg = "Last 2 dimensions of the array must be square."
            raise np.linalg.LinAlgError(msg)

        _check_finite_matrix(a)

        acpy = _copy_to_fortran_order(a)

        w = np.empty(n, dtype=w_dtype)

        if n == 0:
            return (w, acpy)

        r = numba_ez_xxxevd(kind,  # kind
                            JOBZ,  # jobz
                            UPLO,  # uplo
                            n,  # n
                            acpy.ctypes,  # a
                            n,  # lda
                            w.ctypes  # w
                            )
        _handle_err_maybe_convergence_problem(r)

        # help liveness analysis
        _dummy_liveness_func([acpy.size, w.size])
        return (w, acpy)

    return eigh_impl

@overload(np.linalg.eigvalsh)
def eigvalsh_impl(a):
    ensure_lapack()

    _check_linalg_matrix(a, "eigvalsh")

    # convert typing floats to numpy floats for use in the impl
    w_type = getattr(a.dtype, "underlying_float", a.dtype)
    w_dtype = np_support.as_dtype(w_type)

    numba_ez_xxxevd = _LAPACK().numba_ez_xxxevd(a.dtype)

    kind = ord(get_blas_kind(a.dtype, "eigvalsh"))

    JOBZ = ord('N')
    UPLO = ord('L')

    def eigvalsh_impl(a):
        n = a.shape[-1]

        if a.shape[-2] != n:
            msg = "Last 2 dimensions of the array must be square."
            raise np.linalg.LinAlgError(msg)

        _check_finite_matrix(a)

        acpy = _copy_to_fortran_order(a)

        w = np.empty(n, dtype=w_dtype)

        if n == 0:
            return w

        r = numba_ez_xxxevd(kind,  # kind
                            JOBZ,  # jobz
                            UPLO,  # uplo
                            n,  # n
                            acpy.ctypes,  # a
                            n,  # lda
                            w.ctypes  # w
                            )
        _handle_err_maybe_convergence_problem(r)

        # help liveness analysis
        _dummy_liveness_func([acpy.size, w.size])
        return w

    return eigvalsh_impl

@overload(np.linalg.svd)
def svd_impl(a, full_matrices=1):
    ensure_lapack()

    _check_linalg_matrix(a, "svd")

    # convert typing floats to numpy floats for use in the impl
    s_type = getattr(a.dtype, "underlying_float", a.dtype)
    s_dtype = np_support.as_dtype(s_type)

    numba_ez_gesdd = _LAPACK().numba_ez_gesdd(a.dtype)

    kind = ord(get_blas_kind(a.dtype, "svd"))

    JOBZ_A = ord('A')
    JOBZ_S = ord('S')

    def svd_impl(a, full_matrices=1):
        n = a.shape[-1]
        m = a.shape[-2]

        if n == 0 or m == 0:
            raise np.linalg.LinAlgError("Arrays cannot be empty")

        _check_finite_matrix(a)

        acpy = _copy_to_fortran_order(a)

        ldu = m
        minmn = min(m, n)

        if full_matrices:
            JOBZ = JOBZ_A
            ucol = m
            ldvt = n
        else:
            JOBZ = JOBZ_S
            ucol = minmn
            ldvt = minmn

        u = np.empty((ucol, ldu), dtype=a.dtype)
        s = np.empty(minmn, dtype=s_dtype)
        vt = np.empty((n, ldvt), dtype=a.dtype)

        r = numba_ez_gesdd(
            kind,  # kind
            JOBZ,  # jobz
            m,  # m
            n,  # n
            acpy.ctypes,  # a
            m,  # lda
            s.ctypes,  # s
            u.ctypes,  # u
            ldu,  # ldu
            vt.ctypes,  # vt
            ldvt          # ldvt
        )
        _handle_err_maybe_convergence_problem(r)

        # help liveness analysis
        _dummy_liveness_func([acpy.size, vt.size, u.size, s.size])
        return (u.T, s, vt.T)

    return svd_impl


@overload(np.linalg.qr)
def qr_impl(a):
    ensure_lapack()

    _check_linalg_matrix(a, "qr")

    # Need two functions, the first computes R, storing it in the upper
    # triangle of A with the below diagonal part of A containing elementary
    # reflectors needed to construct Q. The second turns the below diagonal
    # entries of A into Q, storing Q in A (creates orthonormal columns from
    # the elementary reflectors).

    numba_ez_geqrf = _LAPACK().numba_ez_geqrf(a.dtype)
    numba_ez_xxgqr = _LAPACK().numba_ez_xxgqr(a.dtype)

    kind = ord(get_blas_kind(a.dtype, "qr"))

    def qr_impl(a):
        n = a.shape[-1]
        m = a.shape[-2]

        if n == 0 or m == 0:
            raise np.linalg.LinAlgError("Arrays cannot be empty")

        _check_finite_matrix(a)

        # copy A as it will be destroyed
        q = _copy_to_fortran_order(a)

        lda = m

        minmn = min(m, n)
        tau = np.empty((minmn), dtype=a.dtype)

        ret = numba_ez_geqrf(
            kind,  # kind
            m,  # m
            n,  # n
            q.ctypes,  # a
            m,  # lda
            tau.ctypes  # tau
        )
        if ret < 0:
            fatal_error_func()
            assert 0   # unreachable

        # pull out R, this is transposed because of Fortran
        r = np.zeros((n, minmn), dtype=a.dtype).T

        # the triangle in R
        for i in range(minmn):
            for j in range(i + 1):
                r[j, i] = q[j, i]

        # and the possible square in R
        for i in range(minmn, n):
            for j in range(minmn):
                r[j, i] = q[j, i]

        ret = numba_ez_xxgqr(
            kind,  # kind
            m,  # m
            minmn,  # n
            minmn,  # k
            q.ctypes,  # a
            m,  # lda
            tau.ctypes  # tau
        )
        _handle_err_maybe_convergence_problem(ret)

        # help liveness analysis
        _dummy_liveness_func([tau.size, q.size])
        return (q[:, :minmn], r)

    return qr_impl


# helpers and jitted specialisations required for np.linalg.lstsq
# and np.linalg.solve. These functions have "system" in their name
# as a differentiator.

def _system_copy_in_b(bcpy, b, nrhs):
    """
    Correctly copy 'b' into the 'bcpy' scratch space.
    """
    raise NotImplementedError


@overload(_system_copy_in_b)
def _system_copy_in_b_impl(bcpy, b, nrhs):
    if b.ndim == 1:
        def oneD_impl(bcpy, b, nrhs):
            bcpy[:b.shape[-1], 0] = b
        return oneD_impl
    else:
        def twoD_impl(bcpy, b, nrhs):
            bcpy[:b.shape[-2], :nrhs] = b
        return twoD_impl


def _system_compute_nrhs(b):
    """
    Compute the number of right hand sides in the system of equations
    """
    raise NotImplementedError


@overload(_system_compute_nrhs)
def _system_compute_nrhs_impl(b):
    if b.ndim == 1:
        def oneD_impl(b):
            return 1
        return oneD_impl
    else:
        def twoD_impl(b):
            return b.shape[-1]
        return twoD_impl


def _system_check_dimensionally_valid(a, b):
    """
    Check that AX=B style system input is dimensionally valid.
    """
    raise NotImplementedError


@overload(_system_check_dimensionally_valid)
def _system_check_dimensionally_valid_impl(a, b):
    ndim = b.ndim
    if ndim == 1:
        def oneD_impl(a, b):
            am = a.shape[-2]
            bm = b.shape[-1]
            if am != bm:
                raise np.linalg.LinAlgError(
                    "Incompatible array sizes, system is not dimensionally valid.")
        return oneD_impl
    else:
        def twoD_impl(a, b):
            am = a.shape[-2]
            bm = b.shape[-2]
            if am != bm:
                raise np.linalg.LinAlgError(
                    "Incompatible array sizes, system is not dimensionally valid.")
        return twoD_impl


def _system_check_non_empty(a, b):
    """
    Check that AX=B style system input is not empty.
    """
    raise NotImplementedError


@overload(_system_check_non_empty)
def _system_check_non_empty_impl(a, b):
    ndim = b.ndim
    if ndim == 1:
        def oneD_impl(a, b):
            am = a.shape[-2]
            an = a.shape[-1]
            bm = b.shape[-1]
            if am == 0 or bm == 0 or an == 0:
                raise np.linalg.LinAlgError('Arrays cannot be empty')
        return oneD_impl
    else:
        def twoD_impl(a, b):
            am = a.shape[-2]
            an = a.shape[-1]
            bm = b.shape[-2]
            bn = b.shape[-1]
            if am == 0 or bm == 0 or an == 0 or bn == 0:
                raise np.linalg.LinAlgError('Arrays cannot be empty')
        return twoD_impl


def _lstsq_residual(b, n, nrhs):
    """
    Compute the residual from the 'b' scratch space.
    """
    raise NotImplementedError


@overload(_lstsq_residual)
def _lstsq_residual_impl(b, n, nrhs):
    ndim = b.ndim
    dtype = b.dtype
    real_dtype = np_support.as_dtype(getattr(dtype, "underlying_float", dtype))

    if ndim == 1:
        if isinstance(dtype, (types.Complex)):
            def cmplx_impl(b, n, nrhs):
                res = np.empty((1,), dtype=real_dtype)
                res[0] = np.sum(np.abs(b[n:, 0])**2)
                return res
            return cmplx_impl
        else:
            def real_impl(b, n, nrhs):
                res = np.empty((1,), dtype=real_dtype)
                res[0] = np.sum(b[n:, 0]**2)
                return res
            return real_impl
    else:
        assert ndim == 2
        if isinstance(dtype, (types.Complex)):
            def cmplx_impl(b, n, nrhs):
                res = np.empty((nrhs), dtype=real_dtype)
                for k in range(nrhs):
                    res[k] = np.sum(np.abs(b[n:, k])**2)
                return res
            return cmplx_impl
        else:
            def real_impl(b, n, nrhs):
                res = np.empty((nrhs), dtype=real_dtype)
                for k in range(nrhs):
                    res[k] = np.sum(b[n:, k]**2)
                return res
            return real_impl


def _lstsq_solution(b, bcpy, n):
    """
    Extract 'x' (the lstsq solution) from the 'bcpy' scratch space.
    Note 'b' is only used to check the system input dimension...
    """
    raise NotImplementedError


@overload(_lstsq_solution)
def _lstsq_solution_impl(b, bcpy, n):
    if b.ndim == 1:
        def oneD_impl(b, bcpy, n):
            return bcpy.T.ravel()[:n]
        return oneD_impl
    else:
        def twoD_impl(b, bcpy, n):
            return bcpy[:n, :].copy()
        return twoD_impl


@overload(np.linalg.lstsq)
def lstsq_impl(a, b, rcond=-1.0):
    ensure_lapack()

    _check_linalg_matrix(a, "lstsq")

    # B can be 1D or 2D.
    _check_linalg_1_or_2d_matrix(b, "lstsq")

    _check_homogeneous_types("lstsq", a, b)

    np_dt = np_support.as_dtype(a.dtype)
    nb_dt = a.dtype

    # convert typing floats to np floats for use in the impl
    r_type = getattr(nb_dt, "underlying_float", nb_dt)
    real_dtype = np_support.as_dtype(r_type)

    # lapack solver
    numba_ez_gelsd = _LAPACK().numba_ez_gelsd(a.dtype)

    kind = ord(get_blas_kind(nb_dt, "lstsq"))

    # The following functions select specialisations based on
    # information around 'b', a lot of this effort is required
    # as 'b' can be either 1D or 2D, and then there are
    # some optimisations available depending on real or complex
    # space.

    def lstsq_impl(a, b, rcond=-1.0):
        n = a.shape[-1]
        m = a.shape[-2]
        nrhs = _system_compute_nrhs(b)

        # check the systems have no inf or NaN
        _check_finite_matrix(a)
        _check_finite_matrix(b)

        # check the system is not empty
        _system_check_non_empty(a, b)

        # check the systems are dimensionally valid
        _system_check_dimensionally_valid(a, b)

        minmn = min(m, n)
        maxmn = max(m, n)

        # a is destroyed on exit, copy it
        acpy = _copy_to_fortran_order(a)

        # b is overwritten on exit with the solution, copy allocate
        bcpy = np.empty((nrhs, maxmn), dtype=np_dt).T
        # specialised copy in due to b being 1 or 2D
        _system_copy_in_b(bcpy, b, nrhs)

        # Allocate returns
        s = np.empty(minmn, dtype=real_dtype)
        rank_ptr = np.empty(1, dtype=np.int32)

        r = numba_ez_gelsd(
            kind,  # kind
            m,  # m
            n,  # n
            nrhs,  # nrhs
            acpy.ctypes,  # a
            m,  # lda
            bcpy.ctypes,  # a
            maxmn,  # ldb
            s.ctypes,  # s
            rcond,  # rcond
            rank_ptr.ctypes  # rank
        )
        _handle_err_maybe_convergence_problem(r)

        # set rank to that which was computed
        rank = rank_ptr[0]

        # compute residuals
        if rank < n or m <= n:
            res = np.empty((0), dtype=real_dtype)
        else:
            # this requires additional dispatch as there's a faster
            # impl if the result is in the real domain (no abs() required)
            res = _lstsq_residual(bcpy, n, nrhs)

        # extract 'x', the solution
        x = _lstsq_solution(b, bcpy, n)

        # help liveness analysis
        _dummy_liveness_func([acpy.size, bcpy.size, s.size, rank_ptr.size])
        return (x, res, rank, s[:minmn])

    return lstsq_impl


def _solve_compute_return(b, bcpy):
    """
    Extract 'x' (the solution) from the 'bcpy' scratch space.
    Note 'b' is only used to check the system input dimension...
    """
    raise NotImplementedError


@overload(_solve_compute_return)
def _solve_compute_return_impl(b, bcpy):
    if b.ndim == 1:
        def oneD_impl(b, bcpy):
            return bcpy.T.ravel()
        return oneD_impl
    else:
        def twoD_impl(b, bcpy):
            return bcpy
        return twoD_impl


@overload(np.linalg.solve)
def solve_impl(a, b):
    ensure_lapack()

    _check_linalg_matrix(a, "solve")
    _check_linalg_1_or_2d_matrix(b, "solve")

    _check_homogeneous_types("solve", a, b)

    np_dt = np_support.as_dtype(a.dtype)
    nb_dt = a.dtype

    # the lapack solver
    numba_xgesv = _LAPACK().numba_xgesv(a.dtype)

    kind = ord(get_blas_kind(nb_dt, "solve"))

    def solve_impl(a, b):
        n = a.shape[-1]
        nrhs = _system_compute_nrhs(b)

        # check the systems have no inf or NaN
        _check_finite_matrix(a)
        _check_finite_matrix(b)

        # check the systems are dimensionally valid
        _system_check_dimensionally_valid(a, b)

        # a is destroyed on exit, copy it
        acpy = _copy_to_fortran_order(a)

        # b is overwritten on exit with the solution, copy allocate
        bcpy = np.empty((nrhs, n), dtype=np_dt).T
        if n == 0:
            return _solve_compute_return(b, bcpy)

        # specialised copy in due to b being 1 or 2D
        _system_copy_in_b(bcpy, b, nrhs)

        # allocate pivot array (needs to be fortran int size)
        ipiv = np.empty(n, dtype=F_INT_nptype)

        r = numba_xgesv(
            kind,        # kind
            n,           # n
            nrhs,        # nhrs
            acpy.ctypes,  # a
            n,           # lda
            ipiv.ctypes,  # ipiv
            bcpy.ctypes,  # b
            n            # ldb
        )
        _inv_err_handler(r)

        # help liveness analysis
        _dummy_liveness_func([acpy.size, bcpy.size, ipiv.size])
        return _solve_compute_return(b, bcpy)

    return solve_impl


@overload(np.linalg.pinv)
def pinv_impl(a, rcond=1.e-15):
    ensure_lapack()

    _check_linalg_matrix(a, "pinv")

    # convert typing floats to numpy floats for use in the impl
    s_type = getattr(a.dtype, "underlying_float", a.dtype)
    s_dtype = np_support.as_dtype(s_type)

    numba_ez_gesdd = _LAPACK().numba_ez_gesdd(a.dtype)

    numba_xxgemm = _BLAS().numba_xxgemm(a.dtype)

    kind = ord(get_blas_kind(a.dtype, "pinv"))
    JOB = ord('S')

    # need conjugate transposes
    TRANSA = ord('C')
    TRANSB = ord('C')

    # scalar constants
    dt = np_support.as_dtype(a.dtype)
    zero = np.array([0.], dtype=dt)
    one = np.array([1.], dtype=dt)

    def pinv_impl(a, rcond=1.e-15):

        # The idea is to build the pseudo-inverse via inverting the singular
        # value decomposition of a matrix `A`. Mathematically, this is roughly
        # A = U*S*V^H        [The SV decomposition of A]
        # A^+ = V*(S^+)*U^H  [The inverted SV decomposition of A]
        # where ^+ is pseudo inversion and ^H is Hermitian transpose.
        # As V and U are unitary, their inverses are simply their Hermitian
        # transpose. S has singular values on its diagonal and zero elsewhere,
        # it is inverted trivially by reciprocal of the diagonal values with
        # the exception that zero singular values remain as zero.
        #
        # The practical implementation can take advantage of a few things to
        # gain a few % performance increase:
        # * A is destroyed by the SVD algorithm from LAPACK so a copy is
        #   required, this memory is exactly the right size in which to return
        #   the pseudo-inverse and so can be reused for this purpose.
        # * The pseudo-inverse of S can be applied to either V or U^H, this
        #   then leaves a GEMM operation to compute the inverse via either:
        #   A^+ = (V*(S^+))*U^H
        #   or
        #   A^+ = V*((S^+)*U^H)
        #   however application of S^+ to V^H or U is more convenient as they
        #   are the result of the SVD algorithm. The application of the
        #   diagonal system is just a matrix multiplication which results in a
        #   row/column scaling (direction depending). To save effort, this
        #   "matrix multiplication" is applied to the smallest of U or V^H and
        #   only up to the point of "cut-off" (see next note) just as a direct
        #   scaling.
        # * The cut-off level for application of S^+ can be used to reduce
        #   total effort, this cut-off can come via rcond or may just naturally
        #   be present as a result of zeros in the singular values. Regardless
        #   there's no need to multiply by zeros in the application of S^+ to
        #   V^H or U as above. Further, the GEMM operation can be shrunk in
        #   effort by noting that the possible zero block generated by the
        #   presence of zeros in S^+ has no effect apart from wasting cycles as
        #   it is all fmadd()s where one operand is zero. The inner dimension
        #   of the GEMM operation can therefore be set as shrunk accordingly!

        n = a.shape[-1]
        m = a.shape[-2]

        _check_finite_matrix(a)

        acpy = _copy_to_fortran_order(a)

        if m == 0 or n == 0:
            return acpy.T.ravel().reshape(a.shape).T

        minmn = min(m, n)

        u = np.empty((minmn, m), dtype=a.dtype)
        s = np.empty(minmn, dtype=s_dtype)
        vt = np.empty((n, minmn), dtype=a.dtype)

        r = numba_ez_gesdd(
            kind,         # kind
            JOB,          # job
            m,            # m
            n,            # n
            acpy.ctypes,  # a
            m,            # lda
            s.ctypes,     # s
            u.ctypes,     # u
            m,            # ldu
            vt.ctypes,    # vt
            minmn         # ldvt
        )
        _handle_err_maybe_convergence_problem(r)

        # Invert singular values under threshold. Also find the index of
        # the threshold value as this is the upper limit for the application
        # of the inverted singular values. Finding this value saves
        # multiplication by a block of zeros that would be created by the
        # application of these values to either U or V^H ahead of multiplying
        # them together. This is done by simply in BLAS parlance via
        # restricting the `k` dimension to `cut_idx` in `xgemm` whilst keeping
        # the leading dimensions correct.

        cut_at = s[0] * rcond
        cut_idx = 0
        for k in range(minmn):
            if s[k] > cut_at:
                s[k] = 1. / s[k]
                cut_idx = k
        cut_idx += 1

        # Use cut_idx so there's no scaling by 0.
        if m >= n:
            # U is largest so apply S^+ to V^H.
            for i in range(n):
                for j in range(cut_idx):
                    vt[i, j] = vt[i, j] * s[j]
        else:
            # V^H is largest so apply S^+ to U.
            for i in range(cut_idx):
                s_local = s[i]
                for j in range(minmn):
                    u[i, j] = u[i, j] * s_local

        # Do (v^H)^H*U^H (obviously one of the matrices includes the S^+
        # scaling) and write back to acpy. Note the innner dimension of cut_idx
        # taking account of the possible zero block.
        # We can store the result in acpy, given we had to create it
        # for use in the SVD, and it is now redundant and the right size
        # but wrong shape.

        r = numba_xxgemm(
            kind,
            TRANSA,       # TRANSA
            TRANSB,       # TRANSB
            n,            # M
            m,            # N
            cut_idx,      # K
            one.ctypes,   # ALPHA
            vt.ctypes,    # A
            minmn,        # LDA
            u.ctypes,     # B
            m,            # LDB
            zero.ctypes,  # BETA
            acpy.ctypes,  # C
            n             # LDC
        )

        # help liveness analysis
        #acpy.size
        #vt.size
        #u.size
        #s.size
        #one.size
        #zero.size
        _dummy_liveness_func([acpy.size, vt.size, u.size, s.size, one.size,
            zero.size])
        return acpy.T.ravel().reshape(a.shape).T

    return pinv_impl


def _get_slogdet_diag_walker(a):
    """
    Walks the diag of a LUP decomposed matrix
    uses that det(A) = prod(diag(lup(A)))
    and also that log(a)+log(b) = log(a*b)
    The return sign is adjusted based on the values found
    such that the log(value) stays in the real domain.
    """
    if isinstance(a.dtype, types.Complex):
        @register_jitable
        def cmplx_diag_walker(n, a, sgn):
            # walk diagonal
            csgn = sgn + 0.j
            acc = 0.
            for k in range(n):
                absel = np.abs(a[k, k])
                csgn = csgn * (a[k, k] / absel)
                acc = acc + np.log(absel)
            return (csgn, acc)
        return cmplx_diag_walker
    else:
        @register_jitable
        def real_diag_walker(n, a, sgn):
            # walk diagonal
            acc = 0.
            for k in range(n):
                v = a[k, k]
                if v < 0.:
                    sgn = -sgn
                    v = -v
                acc = acc + np.log(v)
            # sgn is a float dtype
            return (sgn + 0., acc)
        return real_diag_walker


@overload(np.linalg.slogdet)
def slogdet_impl(a):
    ensure_lapack()

    _check_linalg_matrix(a, "slogdet")

    numba_xxgetrf = _LAPACK().numba_xxgetrf(a.dtype)

    kind = ord(get_blas_kind(a.dtype, "slogdet"))

    diag_walker = _get_slogdet_diag_walker(a)

    ONE = a.dtype(1)
    ZERO = getattr(a.dtype, "underlying_float", a.dtype)(0)

    def slogdet_impl(a):
        n = a.shape[-1]
        if a.shape[-2] != n:
            msg = "Last 2 dimensions of the array must be square."
            raise np.linalg.LinAlgError(msg)

        if n == 0:
            return (ONE, ZERO)

        _check_finite_matrix(a)

        acpy = _copy_to_fortran_order(a)

        ipiv = np.empty(n, dtype=F_INT_nptype)

        r = numba_xxgetrf(kind, n, n, acpy.ctypes, n, ipiv.ctypes)

        if r > 0:
            # factorisation failed, return same defaults as np
            return (0., -np.inf)
        _inv_err_handler(r)  # catch input-to-lapack problem

        # The following, prior to the call to diag_walker, is present
        # to account for the effect of possible permutations to the
        # sign of the determinant.
        # This is the same idea as in numpy:
        # File name `umath_linalg.c.src` e.g.
        # https://github.com/numpy/numpy/blob/master/numpy/linalg/umath_linalg.c.src
        # in function `@TYPE@_slogdet_single_element`.
        sgn = 1
        for k in range(n):
            sgn = sgn + (ipiv[k] != (k + 1))

        sgn = sgn & 1
        if sgn == 0:
            sgn = -1

        # help liveness analysis
        _dummy_liveness_func([ipiv.size])
        return diag_walker(n, acpy, sgn)

    return slogdet_impl


@overload(np.linalg.det)
def det_impl(a):

    ensure_lapack()

    _check_linalg_matrix(a, "det")

    def det_impl(a):
        (sgn, slogdet) = np.linalg.slogdet(a)
        return sgn * np.exp(slogdet)

    return det_impl


def _compute_singular_values(a):
    """
    Compute singular values of *a*.
    """
    raise NotImplementedError


@overload(_compute_singular_values)
def _compute_singular_values_impl(a):
    """
    Returns a function to compute singular values of `a`
    """
    numba_ez_gesdd = _LAPACK().numba_ez_gesdd(a.dtype)

    kind = ord(get_blas_kind(a.dtype, "svd"))

    # Flag for "only compute `S`" to give to xgesdd
    JOBZ_N = ord('N')

    nb_ret_type = getattr(a.dtype, "underlying_float", a.dtype)
    np_ret_type = np_support.as_dtype(nb_ret_type)
    np_dtype = np_support.as_dtype(a.dtype)

    # These are not referenced in the computation but must be set
    # for MKL.
    u = np.empty((1, 1), dtype=np_dtype)
    vt = np.empty((1, 1), dtype=np_dtype)

    def sv_function(a):
        """
        Computes singular values.
        """
        # Don't use the np.linalg.svd impl instead
        # call LAPACK to shortcut doing the "reconstruct
        # singular vectors from reflectors" step and just
        # get back the singular values.
        n = a.shape[-1]
        m = a.shape[-2]
        if m == 0 or n == 0:
            raise np.linalg.LinAlgError('Arrays cannot be empty')
        _check_finite_matrix(a)

        ldu = m
        minmn = min(m, n)

        # need to be >=1 but aren't referenced
        ucol = 1
        ldvt = 1

        acpy = _copy_to_fortran_order(a)

        # u and vt are not referenced however need to be
        # allocated (as done above) for MKL as it
        # checks for ref is nullptr.
        s = np.empty(minmn, dtype=np_ret_type)

        r = numba_ez_gesdd(
            kind,        # kind
            JOBZ_N,      # jobz
            m,           # m
            n,           # n
            acpy.ctypes,  # a
            m,           # lda
            s.ctypes,    # s
            u.ctypes,    # u
            ldu,         # ldu
            vt.ctypes,   # vt
            ldvt         # ldvt
        )
        _handle_err_maybe_convergence_problem(r)

        # help liveness analysis
        _dummy_liveness_func([acpy.size, vt.size, u.size, s.size])
        return s

    return sv_function


def _oneD_norm_2(a):
    """
    Compute the L2-norm of 1D-array *a*.
    """
    raise NotImplementedError


@overload(_oneD_norm_2)
def _oneD_norm_2_impl(a):
    nb_ret_type = getattr(a.dtype, "underlying_float", a.dtype)
    np_ret_type = np_support.as_dtype(nb_ret_type)

    xxnrm2 = _BLAS().numba_xxnrm2(a.dtype)

    kind = ord(get_blas_kind(a.dtype, "norm"))

    def impl(a):
        # Just ignore order, calls are guarded to only come
        # from cases where order=None or order=2.
        n = len(a)
        # Call L2-norm routine from BLAS
        ret = np.empty((1,), dtype=np_ret_type)
        jmp = int(a.strides[0] / a.itemsize)
        r = xxnrm2(
            kind,      # kind
            n,         # n
            a.ctypes,  # x
            jmp,       # incx
            ret.ctypes  # result
        )
        if r < 0:
            fatal_error_func()
            assert 0   # unreachable

        # help liveness analysis
        #ret.size
        #a.size
        _dummy_liveness_func([ret.size, a.size])
        return ret[0]

    return impl


def _get_norm_impl(x, ord_flag):
    # This function is quite involved as norm supports a large
    # range of values to select different norm types via kwarg `ord`.
    # The implementation below branches on dimension of the input
    # (1D or 2D). The default for `ord` is `None` which requires
    # special handling in numba, this is dealt with first in each of
    # the dimension branches. Following this the various norms are
    # computed via code that is in most cases simply a loop version
    # of a ufunc based version as found in numpy.

    # The following is common to both 1D and 2D cases.
    # Convert typing floats to numpy floats for use in the impl.
    # The return type is always a float, numba differs from numpy in
    # that it returns an input precision specific value whereas numpy
    # always returns np.float64.
    nb_ret_type = getattr(x.dtype, "underlying_float", x.dtype)
    np_ret_type = np_support.as_dtype(nb_ret_type)

    np_dtype = np_support.as_dtype(x.dtype)

    xxnrm2 = _BLAS().numba_xxnrm2(x.dtype)

    kind = ord(get_blas_kind(x.dtype, "norm"))

    if x.ndim == 1:
        # 1D cases

        # handle "ord" being "None", must be done separately
        if ord_flag in (None, types.none):
            def oneD_impl(x, ord=None):
                return _oneD_norm_2(x)
        else:
            def oneD_impl(x, ord=None):
                n = len(x)

                # Shortcut to handle zero length arrays
                # this differs slightly to numpy in that
                # numpy raises a ValueError for kwarg ord=
                # +/-np.inf as the reduction operations like
                # max() and min() don't accept zero length
                # arrays
                if n == 0:
                    return 0.0

                # Note: on order == 2
                # This is the same as for ord=="None" but because
                # we have to handle "None" specially this condition
                # is separated
                if ord == 2:
                    return _oneD_norm_2(x)
                elif ord == np.inf:
                    # max(abs(x))
                    ret = abs(x[0])
                    for k in range(1, n):
                        val = abs(x[k])
                        if val > ret:
                            ret = val
                    return ret

                elif ord == -np.inf:
                    # min(abs(x))
                    ret = abs(x[0])
                    for k in range(1, n):
                        val = abs(x[k])
                        if val < ret:
                            ret = val
                    return ret

                elif ord == 0:
                    # sum(x != 0)
                    ret = 0.0
                    for k in range(n):
                        if x[k] != 0.:
                            ret += 1.
                    return ret

                elif ord == 1:
                    # sum(abs(x))
                    ret = 0.0
                    for k in range(n):
                        ret += abs(x[k])
                    return ret

                else:
                    # sum(abs(x)**ord)**(1./ord)
                    ret = 0.0
                    for k in range(n):
                        ret += abs(x[k])**ord
                    return ret**(1. / ord)
        return oneD_impl

    elif x.ndim == 2:
        # 2D cases

        # handle "ord" being "None"
        if ord_flag in (None, types.none):
            # Force `x` to be C-order, so that we can take a contiguous
            # 1D view.
            if x.layout == 'C':
                @register_jitable
                def array_prepare(x):
                    return x
            elif x.layout == 'F':
                @register_jitable
                def array_prepare(x):
                    # Legal since L2(x) == L2(x.T)
                    return x.T
            else:
                @register_jitable
                def array_prepare(x):
                    return x.copy()

            # Compute the Frobenius norm, this is the L2,2 induced norm of `x`
            # which is the L2-norm of x.ravel() and so can be computed via BLAS
            def twoD_impl(x, ord=None):
                n = x.size
                if n == 0:
                    # reshape() currently doesn't support zero-sized arrays
                    return 0.0
                x_c = array_prepare(x)
                return _oneD_norm_2(x_c.reshape(n))
        else:
            # max value for this dtype
            max_val = np.finfo(np_ret_type.type).max

            def twoD_impl(x, ord=None):
                n = x.shape[-1]
                m = x.shape[-2]

                # Shortcut to handle zero size arrays
                # this differs slightly to numpy in that
                # numpy raises errors for some ord values
                # and in other cases returns zero.
                if x.size == 0:
                    return 0.0

                if ord == np.inf:
                    # max of sum of abs across rows
                    # max(sum(abs(x)), axis=1)
                    global_max = 0.
                    for ii in range(m):
                        tmp = 0.
                        for jj in range(n):
                            tmp += abs(x[ii, jj])
                        if tmp > global_max:
                            global_max = tmp
                    return global_max

                elif ord == -np.inf:
                    # min of sum of abs across rows
                    # min(sum(abs(x)), axis=1)
                    global_min = max_val
                    for ii in range(m):
                        tmp = 0.
                        for jj in range(n):
                            tmp += abs(x[ii, jj])
                        if tmp < global_min:
                            global_min = tmp
                    return global_min
                elif ord == 1:
                    # max of sum of abs across cols
                    # max(sum(abs(x)), axis=0)
                    global_max = 0.
                    for ii in range(n):
                        tmp = 0.
                        for jj in range(m):
                            tmp += abs(x[jj, ii])
                        if tmp > global_max:
                            global_max = tmp
                    return global_max

                elif ord == -1:
                    # min of sum of abs across cols
                    # min(sum(abs(x)), axis=0)
                    global_min = max_val
                    for ii in range(n):
                        tmp = 0.
                        for jj in range(m):
                            tmp += abs(x[jj, ii])
                        if tmp < global_min:
                            global_min = tmp
                    return global_min

                # Results via SVD, singular values are sorted on return
                # by definition.
                elif ord == 2:
                    # max SV
                    return _compute_singular_values(x)[0]
                elif ord == -2:
                    # min SV
                    return _compute_singular_values(x)[-1]
                else:
                    # replicate numpy error
                    raise ValueError("Invalid norm order for matrices.")
        return twoD_impl
    else:
        assert 0  # unreachable


@overload(np.linalg.norm)
def norm_impl(x, ord=None):
    ensure_lapack()

    _check_linalg_1_or_2d_matrix(x, "norm")

    return _get_norm_impl(x, ord)


@overload(np.linalg.cond)
def cond_impl(x, p=None):
    ensure_lapack()

    _check_linalg_matrix(x, "cond")

    def impl(x, p=None):
        # This is extracted for performance, numpy does approximately:
        # `condition = norm(x) * norm(inv(x))`
        # in the cases of `p == 2` or `p ==-2` singular values are used
        # for computing norms. This costs numpy an svd of `x` then an
        # inversion of `x` and another svd of `x`.
        # Below is a different approach, which also gives a more
        # accurate answer as there is no inversion involved.
        # Recall that the singular values of an inverted matrix are the
        # reciprocal of singular values of the original matrix.
        # Therefore calling `svd(x)` once yields all the information
        # needed about both `x` and `inv(x)` without the cost or
        # potential loss of accuracy incurred through inversion.
        # For the case of `p == 2`, the result is just the ratio of
        # `largest singular value/smallest singular value`, and for the
        # case of `p==-2` the result is simply the
        # `smallest singular value/largest singular value`.
        # As a result of this, numba accepts non-square matrices as
        # input when p==+/-2 as well as when p==None.
        if p == 2 or p == -2 or p is None:
            s = _compute_singular_values(x)
            if p == 2 or p is None:
                r = np.divide(s[0], s[-1])
            else:
                r = np.divide(s[-1], s[0])
        else:  # cases np.inf, -np.inf, 1, -1
            norm_x = np.linalg.norm(x, p)
            norm_inv_x = np.linalg.norm(np.linalg.inv(x), p)
            r = norm_x * norm_inv_x
        # NumPy uses a NaN mask, if the input has a NaN, it will return NaN,
        # Numba calls ban NaN through the use of _check_finite_matrix but this
        # catches cases where NaN occurs through floating point use
        if np.isnan(r):
            return np.inf
        else:
            return r
    return impl


@register_jitable
def _get_rank_from_singular_values(sv, t):
    """
    Gets rank from singular values with cut-off at a given tolerance
    """
    rank = 0
    for k in range(len(sv)):
        if sv[k] > t:
            rank = rank + 1
        else:  # sv is ordered big->small so break on condition not met
            break
    return rank


@overload(np.linalg.matrix_rank)
def matrix_rank_impl(A, tol=None):
    """
    Computes rank for matrices and vectors.
    The only issue that may arise is that because numpy uses double
    precision lapack calls whereas numba uses type specific lapack
    calls, some singular values may differ and therefore counting the
    number of them above a tolerance may lead to different counts,
    and therefore rank, in some cases.
    """
    ensure_lapack()

    _check_linalg_1_or_2d_matrix(A, "matrix_rank")

    def _2d_matrix_rank_impl(A, tol):

        # handle the tol==None case separately for type inference to work
        if tol in (None, types.none):
            nb_type = getattr(A.dtype, "underlying_float", A.dtype)
            np_type = np_support.as_dtype(nb_type)
            eps_val = np.finfo(np_type).eps

            def _2d_tol_none_impl(A, tol=None):
                s = _compute_singular_values(A)
                # replicate numpy default tolerance calculation
                r = A.shape[0]
                c = A.shape[1]
                l = max(r, c)
                t = s[0] * l * eps_val
                return _get_rank_from_singular_values(s, t)
            return _2d_tol_none_impl
        else:
            def _2d_tol_not_none_impl(A, tol=None):
                s = _compute_singular_values(A)
                return _get_rank_from_singular_values(s, tol)
            return _2d_tol_not_none_impl

    def _get_matrix_rank_impl(A, tol):
        ndim = A.ndim
        if ndim == 1:
            # NOTE: Technically, the numpy implementation could be argued as
            # incorrect for the case of a vector (1D matrix). If a tolerance
            # is provided and a vector with a singular value below tolerance is
            # encountered this should report a rank of zero, the numpy
            # implementation does not do this and instead elects to report that
            # if any value in the vector is nonzero then the rank is 1.
            # An example would be [0, 1e-15, 0, 2e-15] which numpy reports as
            # rank 1 invariant of `tol`. The singular value for this vector is
            # obviously sqrt(5)*1e-15 and so a tol of e.g. sqrt(6)*1e-15 should
            # lead to a reported rank of 0 whereas a tol of 1e-15 should lead
            # to a reported rank of 1, numpy reports 1 regardless.
            # The code below replicates the numpy behaviour.
            def _1d_matrix_rank_impl(A, tol=None):
                for k in range(len(A)):
                    if A[k] != 0.:
                        return 1
                return 0
            return _1d_matrix_rank_impl
        elif ndim == 2:
            return _2d_matrix_rank_impl(A, tol)
        else:
            assert 0  # unreachable

    return _get_matrix_rank_impl(A, tol)


@overload(np.linalg.matrix_power)
def matrix_power_impl(a, n):
    """
    Computes matrix power. Only integer powers are supported in numpy.
    """

    _check_linalg_matrix(a, "matrix_power")
    np_dtype = np_support.as_dtype(a.dtype)

    nt = getattr(n, 'dtype', n)
    if not isinstance(nt, types.Integer):
        raise NumbaTypeError("Exponent must be an integer.")

    def matrix_power_impl(a, n):

        if n == 0:
            # this should be eye() but it doesn't support
            # the dtype kwarg yet so do it manually to save
            # the copy required by eye(a.shape[0]).asdtype()
            A = np.zeros(a.shape, dtype=np_dtype)
            for k in range(a.shape[0]):
                A[k, k] = 1.
            return A

        am, an = a.shape[-1], a.shape[-2]
        if am != an:
            raise ValueError('input must be a square array')

        # empty, return a copy
        if am == 0:
            return a.copy()

        # note: to be consistent over contiguousness, C order is
        # returned as that is what dot() produces and the most common
        # paths through matrix_power will involve that. Therefore
        # copies are made here to ensure the data ordering is
        # correct for paths not going via dot().

        if n < 0:
            A = np.linalg.inv(a).copy()
            if n == -1:  # return now
                return A
            n = -n
        else:
            if n == 1:  # return a copy now
                return a.copy()
            A = a  # this is safe, `a` is only read

        if n < 4:
            if n == 2:
                return np.dot(A, A)
            if n == 3:
                return np.dot(np.dot(A, A), A)
        else:

            acc = A
            exp = n

            # Initialise ret, SSA cannot see the loop will execute, without this
            # it appears as uninitialised.
            ret = acc
            # tried a loop split and branchless using identity matrix as
            # input but it seems like having a "first entry" flag is quicker
            flag = True
            while exp != 0:
                if exp & 1:
                    if flag:
                        ret = acc
                        flag = False
                    else:
                        ret = np.dot(ret, acc)
                acc = np.dot(acc, acc)
                exp = exp >> 1

            return ret

    return matrix_power_impl

# This is documented under linalg despite not being in the module


@overload(np.trace)
def matrix_trace_impl(a, offset=0):
    """
    Computes the trace of an array.
    """

    _check_linalg_matrix(a, "trace", la_prefix=False)

    if not isinstance(offset, (int, types.Integer)):
        raise NumbaTypeError("integer argument expected, got %s" % offset)

    def matrix_trace_impl(a, offset=0):
        rows, cols = a.shape
        k = offset
        if k < 0:
            rows = rows + k
        if k > 0:
            cols = cols - k
        n = max(min(rows, cols), 0)
        ret = 0
        if k >= 0:
            for i in range(n):
                ret += a[i, k + i]
        else:
            for i in range(n):
                ret += a[i - k, i]
        return ret

    return matrix_trace_impl


def _check_scalar_or_lt_2d_mat(a, func_name, la_prefix=True):
    prefix = "np.linalg" if la_prefix else "np"
    interp = (prefix, func_name)
    # checks that a matrix is 1 or 2D
    if isinstance(a, types.Array):
        if not a.ndim <= 2:
            raise TypingError("%s.%s() only supported on 1 and 2-D arrays "
                              % interp, highlighting=False)


@register_jitable
def outer_impl_none(a, b, out):
    aa = np.asarray(a)
    bb = np.asarray(b)
    return np.multiply(aa.ravel().reshape((aa.size, 1)),
                        bb.ravel().reshape((1, bb.size)))


@register_jitable
def outer_impl_arr(a, b, out):
    aa = np.asarray(a)
    bb = np.asarray(b)
    np.multiply(aa.ravel().reshape((aa.size, 1)),
                bb.ravel().reshape((1, bb.size)),
                out)
    return out


def _get_outer_impl(a, b, out):
    if out in (None, types.none):
        return outer_impl_none
    else:
        return outer_impl_arr


@overload(np.outer)
def outer_impl(a, b, out=None):

    _check_scalar_or_lt_2d_mat(a, "outer", la_prefix=False)
    _check_scalar_or_lt_2d_mat(b, "outer", la_prefix=False)

    impl = _get_outer_impl(a, b, out)

    def outer_impl(a, b, out=None):
        return impl(a, b, out)

    return outer_impl


def _kron_normaliser_impl(x):
    # makes x into a 2d array
    if isinstance(x, types.Array):
        if x.layout not in ('C', 'F'):
            raise TypingError("np.linalg.kron only supports 'C' or 'F' layout "
                              "input arrays. Received an input of "
                              "layout '{}'.".format(x.layout))
        elif x.ndim == 2:
            @register_jitable
            def nrm_shape(x):
                xn = x.shape[-1]
                xm = x.shape[-2]
                return x.reshape(xm, xn)
            return nrm_shape
        else:
            @register_jitable
            def nrm_shape(x):
                xn = x.shape[-1]
                return x.reshape(1, xn)
            return nrm_shape
    else:  # assume its a scalar
        @register_jitable
        def nrm_shape(x):
            a = np.empty((1, 1), type(x))
            a[0] = x
            return a
        return nrm_shape


def _kron_return(a, b):
    # transforms c into something that kron would return
    # based on the shapes of a and b
    a_is_arr = isinstance(a, types.Array)
    b_is_arr = isinstance(b, types.Array)
    if a_is_arr and b_is_arr:
        if a.ndim == 2 or b.ndim == 2:
            @register_jitable
            def ret(a, b, c):
                return c
            return ret
        else:
            @register_jitable
            def ret(a, b, c):
                return c.reshape(c.size)
            return ret
    else:  # at least one of (a, b) is a scalar
        if a_is_arr:
            @register_jitable
            def ret(a, b, c):
                return c.reshape(a.shape)
            return ret
        elif b_is_arr:
            @register_jitable
            def ret(a, b, c):
                return c.reshape(b.shape)
            return ret
        else:  # both scalars
            @register_jitable
            def ret(a, b, c):
                return c[0]
            return ret


@overload(np.kron)
def kron_impl(a, b):

    _check_scalar_or_lt_2d_mat(a, "kron", la_prefix=False)
    _check_scalar_or_lt_2d_mat(b, "kron", la_prefix=False)

    fix_a = _kron_normaliser_impl(a)
    fix_b = _kron_normaliser_impl(b)
    ret_c = _kron_return(a, b)

    # this is fine because the ufunc for the Hadamard product
    # will reject differing dtypes in a and b.
    dt = getattr(a, 'dtype', a)

    def kron_impl(a, b):

        aa = fix_a(a)
        bb = fix_b(b)

        am = aa.shape[-2]
        an = aa.shape[-1]
        bm = bb.shape[-2]
        bn = bb.shape[-1]

        cm = am * bm
        cn = an * bn

        # allocate c
        C = np.empty((cm, cn), dtype=dt)

        # In practice this is runs quicker than the more obvious
        # `each element of A multiplied by B and assigned to
        # a block in C` like alg.

        # loop over rows of A
        for i in range(am):
            # compute the column offset into C
            rjmp = i * bm
            # loop over rows of B
            for k in range(bm):
                # compute row the offset into C
                irjmp = rjmp + k
                # slice a given row of B
                slc = bb[k, :]
                # loop over columns of A
                for j in range(an):
                    # vectorized assignment of an element of A
                    # multiplied by the current row of B into
                    # a slice of a row of C
                    cjmp = j * bn
                    C[irjmp, cjmp:cjmp + bn] = aa[i, j] * slc

        return ret_c(a, b, C)

    return kron_impl
