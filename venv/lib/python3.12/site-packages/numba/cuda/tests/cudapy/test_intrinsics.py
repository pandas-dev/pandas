import itertools
import numpy as np
import operator
import re
from numba import cuda, int64
from numba.cuda import compile_ptx
from numba.core.errors import TypingError
from numba.core.types import f2
from numba.cuda.testing import (unittest, CUDATestCase, skip_on_cudasim,
                                skip_unless_cc_53)


def simple_threadidx(ary):
    i = cuda.threadIdx.x
    ary[0] = i


def fill_threadidx(ary):
    i = cuda.threadIdx.x
    ary[i] = i


def fill3d_threadidx(ary):
    i = cuda.threadIdx.x
    j = cuda.threadIdx.y
    k = cuda.threadIdx.z

    ary[i, j, k] = (i + 1) * (j + 1) * (k + 1)


def simple_grid1d(ary):
    i = cuda.grid(1)
    ary[i] = i


def simple_grid2d(ary):
    i, j = cuda.grid(2)
    ary[i, j] = i + j


def simple_gridsize1d(ary):
    i = cuda.grid(1)
    x = cuda.gridsize(1)
    if i == 0:
        ary[0] = x


def simple_gridsize2d(ary):
    i, j = cuda.grid(2)
    x, y = cuda.gridsize(2)
    if i == 0 and j == 0:
        ary[0] = x
        ary[1] = y


def intrinsic_forloop_step(c):
    startX, startY = cuda.grid(2)
    gridX = cuda.gridDim.x * cuda.blockDim.x
    gridY = cuda.gridDim.y * cuda.blockDim.y
    height, width = c.shape

    for x in range(startX, width, gridX):
        for y in range(startY, height, gridY):
            c[y, x] = x + y


def simple_popc(ary, c):
    ary[0] = cuda.popc(c)


def simple_fma(ary, a, b, c):
    ary[0] = cuda.fma(a, b, c)


def simple_hadd(ary, a, b):
    ary[0] = cuda.fp16.hadd(a[0], b[0])


def simple_hadd_scalar(ary, a, b):
    ary[0] = cuda.fp16.hadd(a, b)


def simple_hfma(ary, a, b, c):
    ary[0] = cuda.fp16.hfma(a[0], b[0], c[0])


def simple_hfma_scalar(ary, a, b, c):
    ary[0] = cuda.fp16.hfma(a, b, c)


def simple_hsub(ary, a, b):
    ary[0] = cuda.fp16.hsub(a[0], b[0])


def simple_hsub_scalar(ary, a, b):
    ary[0] = cuda.fp16.hsub(a, b)


def simple_hmul(ary, a, b):
    ary[0] = cuda.fp16.hmul(a[0], b[0])


def simple_hmul_scalar(ary, a, b):
    ary[0] = cuda.fp16.hmul(a, b)


def simple_hdiv_scalar(ary, a, b):
    ary[0] = cuda.fp16.hdiv(a, b)


def simple_hdiv_kernel(ary, array_a, array_b):
    i = cuda.grid(1)
    if i < ary.size:
        a = array_a[i]
        b = array_b[i]
        ary[i] = cuda.fp16.hdiv(a, b)


def simple_hneg(ary, a):
    ary[0] = cuda.fp16.hneg(a[0])


def simple_hneg_scalar(ary, a):
    ary[0] = cuda.fp16.hneg(a)


def simple_habs(ary, a):
    ary[0] = cuda.fp16.habs(a[0])


def simple_habs_scalar(ary, a):
    ary[0] = cuda.fp16.habs(a)


def simple_heq_scalar(ary, a, b):
    ary[0] = cuda.fp16.heq(a, b)


def simple_hne_scalar(ary, a, b):
    ary[0] = cuda.fp16.hne(a, b)


def simple_hge_scalar(ary, a, b):
    ary[0] = cuda.fp16.hge(a, b)


def simple_hgt_scalar(ary, a, b):
    ary[0] = cuda.fp16.hgt(a, b)


def simple_hle_scalar(ary, a, b):
    ary[0] = cuda.fp16.hle(a, b)


def simple_hlt_scalar(ary, a, b):
    ary[0] = cuda.fp16.hlt(a, b)


@cuda.jit(device=True)
def hlt_func_1(x, y):
    return cuda.fp16.hlt(x, y)


@cuda.jit(device=True)
def hlt_func_2(x, y):
    return cuda.fp16.hlt(x, y)


def test_multiple_hcmp_1(r, a, b, c):
    # float16 predicates used in two separate functions
    r[0] = hlt_func_1(a, b) and hlt_func_2(b, c)


def test_multiple_hcmp_2(r, a, b, c):
    # The same float16 predicate used in the caller and callee
    r[0] = hlt_func_1(a, b) and cuda.fp16.hlt(b, c)


def test_multiple_hcmp_3(r, a, b, c):
    # Different float16 predicates used in the caller and callee
    r[0] = hlt_func_1(a, b) and cuda.fp16.hge(c, b)


def test_multiple_hcmp_4(r, a, b, c):
    # The same float16 predicates used twice in a function
    r[0] = cuda.fp16.hlt(a, b) and cuda.fp16.hlt(b, c)


def test_multiple_hcmp_5(r, a, b, c):
    # Different float16 predicates used in a function
    r[0] = cuda.fp16.hlt(a, b) and cuda.fp16.hge(c, b)


def simple_hmax_scalar(ary, a, b):
    ary[0] = cuda.fp16.hmax(a, b)


def simple_hmin_scalar(ary, a, b):
    ary[0] = cuda.fp16.hmin(a, b)


def simple_hsin(r, x):
    i = cuda.grid(1)

    if i < len(r):
        r[i] = cuda.fp16.hsin(x[i])


def simple_hcos(r, x):
    i = cuda.grid(1)

    if i < len(r):
        r[i] = cuda.fp16.hcos(x[i])


def simple_hlog(r, x):
    i = cuda.grid(1)

    if i < len(r):
        r[i] = cuda.fp16.hlog(x[i])


def simple_hlog2(r, x):
    i = cuda.grid(1)

    if i < len(r):
        r[i] = cuda.fp16.hlog2(x[i])


def simple_hlog10(r, x):
    i = cuda.grid(1)

    if i < len(r):
        r[i] = cuda.fp16.hlog10(x[i])


def simple_hexp(r, x):
    i = cuda.grid(1)

    if i < len(r):
        r[i] = cuda.fp16.hexp(x[i])


def simple_hexp2(r, x):
    i = cuda.grid(1)

    if i < len(r):
        r[i] = cuda.fp16.hexp2(x[i])


def simple_hsqrt(r, x):
    i = cuda.grid(1)

    if i < len(r):
        r[i] = cuda.fp16.hsqrt(x[i])


def simple_hrsqrt(r, x):

    i = cuda.grid(1)

    if i < len(r):
        r[i] = cuda.fp16.hrsqrt(x[i])


def numpy_hrsqrt(x, dtype):
    return x ** -0.5


def simple_hceil(r, x):
    i = cuda.grid(1)

    if i < len(r):
        r[i] = cuda.fp16.hceil(x[i])


def simple_hfloor(r, x):
    i = cuda.grid(1)

    if i < len(r):
        r[i] = cuda.fp16.hfloor(x[i])


def simple_hrcp(r, x):
    i = cuda.grid(1)

    if i < len(r):
        r[i] = cuda.fp16.hrcp(x[i])


def simple_htrunc(r, x):
    i = cuda.grid(1)

    if i < len(r):
        r[i] = cuda.fp16.htrunc(x[i])


def simple_hrint(r, x):
    i = cuda.grid(1)

    if i < len(r):
        r[i] = cuda.fp16.hrint(x[i])


def simple_cbrt(ary, a):
    ary[0] = cuda.cbrt(a)


def simple_brev(ary, c):
    ary[0] = cuda.brev(c)


def simple_clz(ary, c):
    ary[0] = cuda.clz(c)


def simple_ffs(ary, c):
    ary[0] = cuda.ffs(c)


def simple_round(ary, c):
    ary[0] = round(c)


def simple_round_to(ary, c, ndigits):
    ary[0] = round(c, ndigits)


def branching_with_ifs(a, b, c):
    i = cuda.grid(1)

    if a[i] > 4:
        if b % 2 == 0:
            a[i] = c[i]
        else:
            a[i] = 13
    else:
        a[i] = 3


def branching_with_selps(a, b, c):
    i = cuda.grid(1)

    inner = cuda.selp(b % 2 == 0, c[i], 13)
    a[i] = cuda.selp(a[i] > 4, inner, 3)


def simple_laneid(ary):
    i = cuda.grid(1)
    ary[i] = cuda.laneid


def simple_warpsize(ary):
    ary[0] = cuda.warpsize


def nonliteral_grid(x):
    cuda.grid(x)


def nonliteral_gridsize(x):
    cuda.gridsize(x)


class TestCudaIntrinsic(CUDATestCase):
    def setUp(self):
        super().setUp()
        np.random.seed(0)

    def test_simple_threadidx(self):
        compiled = cuda.jit("void(int32[:])")(simple_threadidx)
        ary = np.ones(1, dtype=np.int32)
        compiled[1, 1](ary)
        self.assertTrue(ary[0] == 0)

    def test_fill_threadidx(self):
        compiled = cuda.jit("void(int32[:])")(fill_threadidx)
        N = 10
        ary = np.ones(N, dtype=np.int32)
        exp = np.arange(N, dtype=np.int32)
        compiled[1, N](ary)
        self.assertTrue(np.all(ary == exp))

    def test_fill3d_threadidx(self):
        X, Y, Z = 4, 5, 6

        def c_contigous():
            compiled = cuda.jit("void(int32[:,:,::1])")(fill3d_threadidx)
            ary = np.zeros((X, Y, Z), dtype=np.int32)
            compiled[1, (X, Y, Z)](ary)
            return ary

        def f_contigous():
            compiled = cuda.jit("void(int32[::1,:,:])")(fill3d_threadidx)
            ary = np.asfortranarray(np.zeros((X, Y, Z), dtype=np.int32))
            compiled[1, (X, Y, Z)](ary)
            return ary

        c_res = c_contigous()
        f_res = f_contigous()
        self.assertTrue(np.all(c_res == f_res))

    @skip_on_cudasim('Cudasim does not check types')
    def test_nonliteral_grid_error(self):
        with self.assertRaisesRegex(TypingError, 'RequireLiteralValue'):
            cuda.jit('void(int32)')(nonliteral_grid)

    @skip_on_cudasim('Cudasim does not check types')
    def test_nonliteral_gridsize_error(self):
        with self.assertRaisesRegex(TypingError, 'RequireLiteralValue'):
            cuda.jit('void(int32)')(nonliteral_gridsize)

    def test_simple_grid1d(self):
        compiled = cuda.jit("void(int32[::1])")(simple_grid1d)
        ntid, nctaid = 3, 7
        nelem = ntid * nctaid
        ary = np.empty(nelem, dtype=np.int32)
        compiled[nctaid, ntid](ary)
        self.assertTrue(np.all(ary == np.arange(nelem)))

    def test_simple_grid2d(self):
        compiled = cuda.jit("void(int32[:,::1])")(simple_grid2d)
        ntid = (4, 3)
        nctaid = (5, 6)
        shape = (ntid[0] * nctaid[0], ntid[1] * nctaid[1])
        ary = np.empty(shape, dtype=np.int32)
        exp = ary.copy()
        compiled[nctaid, ntid](ary)

        for i in range(ary.shape[0]):
            for j in range(ary.shape[1]):
                exp[i, j] = i + j

        self.assertTrue(np.all(ary == exp))

    def test_simple_gridsize1d(self):
        compiled = cuda.jit("void(int32[::1])")(simple_gridsize1d)
        ntid, nctaid = 3, 7
        ary = np.zeros(1, dtype=np.int32)
        compiled[nctaid, ntid](ary)
        self.assertEqual(ary[0], nctaid * ntid)

    @skip_on_cudasim('Requires too many threads')
    def test_issue_9229(self):
        # Ensure that grid and grid size are correct - #9229 showed that they
        # overflowed an int32.
        @cuda.jit
        def f(grid_error, gridsize_error):
            i1 = cuda.grid(1)
            i2 = cuda.blockIdx.x * cuda.blockDim.x + cuda.threadIdx.x
            gs1 = cuda.gridsize(1)
            gs2 = cuda.blockDim.x * cuda.gridDim.x
            if i1 != i2:
                grid_error[0] = 1
            if gs1 != gs2:
                gridsize_error[0] = 1

        grid_error = np.zeros(1, dtype=np.uint64)
        gridsize_error = np.zeros(1, dtype=np.uint64)

        # A large enough grid for thread IDs to overflow an int32
        # (22121216 * 256 = 5663031296, which is greater than 2 ** 32)
        f[22121216, 256](grid_error, gridsize_error)

        self.assertEqual(grid_error[0], 0)
        self.assertEqual(gridsize_error[0], 0)

    @skip_on_cudasim('Tests PTX emission')
    def test_selp(self):
        sig = (int64[:], int64, int64[:])
        cu_branching_with_ifs = cuda.jit(sig)(branching_with_ifs)
        cu_branching_with_selps = cuda.jit(sig)(branching_with_selps)

        n = 32
        b = 6
        c = np.full(shape=32, fill_value=17, dtype=np.int64)

        expected = c.copy()
        expected[:5] = 3

        a = np.arange(n, dtype=np.int64)
        cu_branching_with_ifs[n, 1](a, b, c)
        ptx = cu_branching_with_ifs.inspect_asm(sig)
        self.assertEqual(2, len(re.findall(r'\s+bra\s+', ptx)))
        np.testing.assert_array_equal(a, expected, err_msg='branching')

        a = np.arange(n, dtype=np.int64)
        cu_branching_with_selps[n, 1](a, b, c)
        ptx = cu_branching_with_selps.inspect_asm(sig)
        self.assertEqual(0, len(re.findall(r'\s+bra\s+', ptx)))
        np.testing.assert_array_equal(a, expected, err_msg='selp')

    def test_simple_gridsize2d(self):
        compiled = cuda.jit("void(int32[::1])")(simple_gridsize2d)
        ntid = (4, 3)
        nctaid = (5, 6)
        ary = np.zeros(2, dtype=np.int32)
        compiled[nctaid, ntid](ary)

        self.assertEqual(ary[0], nctaid[0] * ntid[0])
        self.assertEqual(ary[1], nctaid[1] * ntid[1])

    def test_intrinsic_forloop_step(self):
        compiled = cuda.jit("void(int32[:,::1])")(intrinsic_forloop_step)
        ntid = (4, 3)
        nctaid = (5, 6)
        shape = (ntid[0] * nctaid[0], ntid[1] * nctaid[1])
        ary = np.empty(shape, dtype=np.int32)

        compiled[nctaid, ntid](ary)

        gridX, gridY = shape
        height, width = ary.shape
        for i, j in zip(range(ntid[0]), range(ntid[1])):
            startX, startY = gridX + i, gridY + j
            for x in range(startX, width, gridX):
                for y in range(startY, height, gridY):
                    self.assertTrue(ary[y, x] == x + y, (ary[y, x], x + y))

    def test_3dgrid(self):
        @cuda.jit
        def foo(out):
            x, y, z = cuda.grid(3)
            a, b, c = cuda.gridsize(3)
            out[x, y, z] = a * b * c

        arr = np.zeros(9 ** 3, dtype=np.int32).reshape(9, 9, 9)
        foo[(3, 3, 3), (3, 3, 3)](arr)

        np.testing.assert_equal(arr, 9 ** 3)

    def test_3dgrid_2(self):
        @cuda.jit
        def foo(out):
            x, y, z = cuda.grid(3)
            a, b, c = cuda.gridsize(3)
            grid_is_right = (
                x == cuda.threadIdx.x + cuda.blockIdx.x * cuda.blockDim.x and
                y == cuda.threadIdx.y + cuda.blockIdx.y * cuda.blockDim.y and
                z == cuda.threadIdx.z + cuda.blockIdx.z * cuda.blockDim.z
            )
            gridsize_is_right = (a == cuda.blockDim.x * cuda.gridDim.x and
                                 b == cuda.blockDim.y * cuda.gridDim.y and
                                 c == cuda.blockDim.z * cuda.gridDim.z)
            out[x, y, z] = grid_is_right and gridsize_is_right

        x, y, z = (4 * 3, 3 * 2, 2 * 4)
        arr = np.zeros((x * y * z), dtype=np.bool_).reshape(x, y, z)
        foo[(4, 3, 2), (3, 2, 4)](arr)

        self.assertTrue(np.all(arr))

    def test_popc_u4(self):
        compiled = cuda.jit("void(int32[:], uint32)")(simple_popc)
        ary = np.zeros(1, dtype=np.int32)
        compiled[1, 1](ary, 0xF0)
        self.assertEqual(ary[0], 4)

    def test_popc_u8(self):
        compiled = cuda.jit("void(int32[:], uint64)")(simple_popc)
        ary = np.zeros(1, dtype=np.int32)
        compiled[1, 1](ary, 0xF00000000000)
        self.assertEqual(ary[0], 4)

    def test_fma_f4(self):
        compiled = cuda.jit("void(f4[:], f4, f4, f4)")(simple_fma)
        ary = np.zeros(1, dtype=np.float32)
        compiled[1, 1](ary, 2., 3., 4.)
        np.testing.assert_allclose(ary[0], 2 * 3 + 4)

    def test_fma_f8(self):
        compiled = cuda.jit("void(f8[:], f8, f8, f8)")(simple_fma)
        ary = np.zeros(1, dtype=np.float64)
        compiled[1, 1](ary, 2., 3., 4.)
        np.testing.assert_allclose(ary[0], 2 * 3 + 4)

    @skip_unless_cc_53
    def test_hadd(self):
        compiled = cuda.jit("void(f2[:], f2[:], f2[:])")(simple_hadd)
        ary = np.zeros(1, dtype=np.float16)
        arg1 = np.array([3.], dtype=np.float16)
        arg2 = np.array([4.], dtype=np.float16)
        compiled[1, 1](ary, arg1, arg2)
        np.testing.assert_allclose(ary[0], arg1 + arg2)

    @skip_unless_cc_53
    def test_hadd_scalar(self):
        compiled = cuda.jit("void(f2[:], f2, f2)")(simple_hadd_scalar)
        ary = np.zeros(1, dtype=np.float16)
        arg1 = np.float16(3.1415926)
        arg2 = np.float16(3.)
        compiled[1, 1](ary, arg1, arg2)
        ref = arg1 + arg2
        np.testing.assert_allclose(ary[0], ref)

    @skip_on_cudasim('Compilation unsupported in the simulator')
    def test_hadd_ptx(self):
        args = (f2[:], f2, f2)
        ptx, _ = compile_ptx(simple_hadd_scalar, args, cc=(5, 3))
        self.assertIn('add.f16', ptx)

    @skip_unless_cc_53
    def test_hfma(self):
        compiled = cuda.jit("void(f2[:], f2[:], f2[:], f2[:])")(simple_hfma)
        ary = np.zeros(1, dtype=np.float16)
        arg1 = np.array([2.], dtype=np.float16)
        arg2 = np.array([3.], dtype=np.float16)
        arg3 = np.array([4.], dtype=np.float16)
        compiled[1, 1](ary, arg1, arg2, arg3)
        np.testing.assert_allclose(ary[0], arg1 * arg2 + arg3)

    @skip_unless_cc_53
    def test_hfma_scalar(self):
        compiled = cuda.jit("void(f2[:], f2, f2, f2)")(simple_hfma_scalar)
        ary = np.zeros(1, dtype=np.float16)
        arg1 = np.float16(2.)
        arg2 = np.float16(3.)
        arg3 = np.float16(4.)
        compiled[1, 1](ary, arg1, arg2, arg3)
        ref = arg1 * arg2 + arg3
        np.testing.assert_allclose(ary[0], ref)

    @skip_on_cudasim('Compilation unsupported in the simulator')
    def test_hfma_ptx(self):
        args = (f2[:], f2, f2, f2)
        ptx, _ = compile_ptx(simple_hfma_scalar, args, cc=(5, 3))
        self.assertIn('fma.rn.f16', ptx)

    @skip_unless_cc_53
    def test_hsub(self):
        compiled = cuda.jit("void(f2[:], f2[:], f2[:])")(simple_hsub)
        ary = np.zeros(1, dtype=np.float16)
        arg1 = np.array([3.], dtype=np.float16)
        arg2 = np.array([4.], dtype=np.float16)
        compiled[1, 1](ary, arg1, arg2)
        np.testing.assert_allclose(ary[0], arg1 - arg2)

    @skip_unless_cc_53
    def test_hsub_scalar(self):
        compiled = cuda.jit("void(f2[:], f2, f2)")(simple_hsub_scalar)
        ary = np.zeros(1, dtype=np.float16)
        arg1 = np.float16(3.1415926)
        arg2 = np.float16(1.57)
        compiled[1, 1](ary, arg1, arg2)
        ref = arg1 - arg2
        np.testing.assert_allclose(ary[0], ref)

    @skip_on_cudasim('Compilation unsupported in the simulator')
    def test_hsub_ptx(self):
        args = (f2[:], f2, f2)
        ptx, _ = compile_ptx(simple_hsub_scalar, args, cc=(5, 3))
        self.assertIn('sub.f16', ptx)

    @skip_unless_cc_53
    def test_hmul(self):
        compiled = cuda.jit()(simple_hmul)
        ary = np.zeros(1, dtype=np.float16)
        arg1 = np.array([3.], dtype=np.float16)
        arg2 = np.array([4.], dtype=np.float16)
        compiled[1, 1](ary, arg1, arg2)
        np.testing.assert_allclose(ary[0], arg1 * arg2)

    @skip_unless_cc_53
    def test_hmul_scalar(self):
        compiled = cuda.jit("void(f2[:], f2, f2)")(simple_hmul_scalar)
        ary = np.zeros(1, dtype=np.float16)
        arg1 = np.float16(3.1415926)
        arg2 = np.float16(1.57)
        compiled[1, 1](ary, arg1, arg2)
        ref = arg1 * arg2
        np.testing.assert_allclose(ary[0], ref)

    @skip_on_cudasim('Compilation unsupported in the simulator')
    def test_hmul_ptx(self):
        args = (f2[:], f2, f2)
        ptx, _ = compile_ptx(simple_hmul_scalar, args, cc=(5, 3))
        self.assertIn('mul.f16', ptx)

    @skip_unless_cc_53
    def test_hdiv_scalar(self):
        compiled = cuda.jit("void(f2[:], f2, f2)")(simple_hdiv_scalar)
        ary = np.zeros(1, dtype=np.float16)
        arg1 = np.float16(3.1415926)
        arg2 = np.float16(1.57)

        compiled[1, 1](ary, arg1, arg2)
        ref = arg1 / arg2
        np.testing.assert_allclose(ary[0], ref)

    @skip_unless_cc_53
    def test_hdiv(self):
        compiled = cuda.jit("void(f2[:], f2[:], f2[:])")(simple_hdiv_kernel)
        arry1 = np.random.randint(-65504, 65505, size=500).astype(np.float16)
        arry2 = np.random.randint(-65504, 65505, size=500).astype(np.float16)
        ary = np.zeros_like(arry1, dtype=np.float16)

        compiled.forall(ary.size)(ary, arry1, arry2)
        ref = arry1 / arry2
        np.testing.assert_allclose(ary, ref)

    @skip_unless_cc_53
    def test_hneg(self):
        compiled = cuda.jit("void(f2[:], f2[:])")(simple_hneg)
        ary = np.zeros(1, dtype=np.float16)
        arg1 = np.array([3.], dtype=np.float16)
        compiled[1, 1](ary, arg1)
        np.testing.assert_allclose(ary[0], -arg1)

    @skip_unless_cc_53
    def test_hneg_scalar(self):
        compiled = cuda.jit("void(f2[:], f2)")(simple_hneg_scalar)
        ary = np.zeros(1, dtype=np.float16)
        arg1 = np.float16(3.1415926)
        compiled[1, 1](ary, arg1)
        ref = -arg1
        np.testing.assert_allclose(ary[0], ref)

    @skip_on_cudasim('Compilation unsupported in the simulator')
    def test_hneg_ptx(self):
        args = (f2[:], f2)
        ptx, _ = compile_ptx(simple_hneg_scalar, args, cc=(5, 3))
        self.assertIn('neg.f16', ptx)

    @skip_unless_cc_53
    def test_habs(self):
        compiled = cuda.jit()(simple_habs)
        ary = np.zeros(1, dtype=np.float16)
        arg1 = np.array([-3.], dtype=np.float16)
        compiled[1, 1](ary, arg1)
        np.testing.assert_allclose(ary[0], abs(arg1))

    @skip_unless_cc_53
    def test_habs_scalar(self):
        compiled = cuda.jit("void(f2[:], f2)")(simple_habs_scalar)
        ary = np.zeros(1, dtype=np.float16)
        arg1 = np.float16(-3.1415926)
        compiled[1, 1](ary, arg1)
        ref = abs(arg1)
        np.testing.assert_allclose(ary[0], ref)

    @skip_on_cudasim('Compilation unsupported in the simulator')
    def test_habs_ptx(self):
        args = (f2[:], f2)
        ptx, _ = compile_ptx(simple_habs_scalar, args, cc=(5, 3))
        self.assertIn('abs.f16', ptx)

    @skip_unless_cc_53
    def test_fp16_intrinsics_common(self):
        kernels = (simple_hsin, simple_hcos,
                   simple_hlog, simple_hlog2, simple_hlog10,
                   simple_hsqrt, simple_hceil, simple_hfloor,
                   simple_hrcp, simple_htrunc, simple_hrint,
                   simple_hrsqrt)
        exp_kernels = (simple_hexp, simple_hexp2)
        expected_functions = (np.sin, np.cos,
                              np.log, np.log2, np.log10,
                              np.sqrt, np.ceil, np.floor,
                              np.reciprocal, np.trunc, np.rint,
                              numpy_hrsqrt)
        expected_exp_functions = (np.exp, np.exp2)

        # Generate random data
        N = 32
        np.random.seed(1)
        x = np.random.randint(1, 65505, size=N).astype(np.float16)
        r = np.zeros_like(x)
        for kernel, fn in zip(kernels, expected_functions):
            with self.subTest(fn=fn):
                kernel = cuda.jit("void(f2[:], f2[:])")(kernel)
                kernel[1,N](r, x)
                expected = fn(x, dtype=np.float16)
                np.testing.assert_allclose(r, expected)

        x2 = np.random.randint(1, 10, size=N).astype(np.float16)
        for kernel, fn in zip(exp_kernels, expected_exp_functions):
            with self.subTest(fn=fn):
                kernel = cuda.jit("void(f2[:], f2[:])")(kernel)
                kernel[1,N](r, x2)
                expected = fn(x2, dtype=np.float16)
                np.testing.assert_allclose(r, expected)

    @skip_unless_cc_53
    def test_hexp10(self):
        @cuda.jit()
        def hexp10_vectors(r, x):
            i = cuda.grid(1)

            if i < len(r):
                r[i] = cuda.fp16.hexp10(x[i])

        # Generate random data
        N = 32
        np.random.seed(1)
        x = np.random.rand(N).astype(np.float16)
        r = np.zeros_like(x)

        # Run the kernel
        hexp10_vectors[1, N](r, x)
        np.testing.assert_allclose(r, 10 ** x)

    @skip_unless_cc_53
    def test_fp16_comparison(self):
        fns = (simple_heq_scalar, simple_hne_scalar, simple_hge_scalar,
               simple_hgt_scalar, simple_hle_scalar, simple_hlt_scalar)
        ops = (operator.eq, operator.ne, operator.ge,
               operator.gt, operator.le, operator.lt)

        for fn, op in zip(fns, ops):
            with self.subTest(op=op):
                kernel = cuda.jit("void(b1[:], f2, f2)")(fn)

                expected = np.zeros(1, dtype=np.bool_)
                got = np.zeros(1, dtype=np.bool_)
                arg2 = np.float16(2)
                arg3 = np.float16(3)
                arg4 = np.float16(4)

                # Check with equal arguments
                kernel[1, 1](got, arg3, arg3)
                expected = op(arg3, arg3)
                self.assertEqual(expected, got[0])

                # Check with LHS < RHS
                kernel[1, 1](got, arg3, arg4)
                expected = op(arg3, arg4)
                self.assertEqual(expected, got[0])

                # Check with LHS > RHS
                kernel[1, 1](got, arg3, arg2)
                expected = op(arg3, arg2)
                self.assertEqual(expected, got[0])

    @skip_unless_cc_53
    def test_multiple_float16_comparisons(self):
        functions = (test_multiple_hcmp_1,
                     test_multiple_hcmp_2,
                     test_multiple_hcmp_3,
                     test_multiple_hcmp_4,
                     test_multiple_hcmp_5)
        for fn in functions:
            with self.subTest(fn=fn):
                compiled = cuda.jit("void(b1[:], f2, f2, f2)")(fn)
                ary = np.zeros(1, dtype=np.bool_)
                arg1 = np.float16(2.)
                arg2 = np.float16(3.)
                arg3 = np.float16(4.)
                compiled[1, 1](ary, arg1, arg2, arg3)
                self.assertTrue(ary[0])

    @skip_unless_cc_53
    def test_hmax(self):
        compiled = cuda.jit("void(f2[:], f2, f2)")(simple_hmax_scalar)
        ary = np.zeros(1, dtype=np.float16)
        arg1 = np.float16(3.)
        arg2 = np.float16(4.)
        compiled[1, 1](ary, arg1, arg2)
        np.testing.assert_allclose(ary[0], arg2)
        arg1 = np.float16(5.)
        compiled[1, 1](ary, arg1, arg2)
        np.testing.assert_allclose(ary[0], arg1)

    @skip_unless_cc_53
    def test_hmin(self):
        compiled = cuda.jit("void(f2[:], f2, f2)")(simple_hmin_scalar)
        ary = np.zeros(1, dtype=np.float16)
        arg1 = np.float16(3.)
        arg2 = np.float16(4.)
        compiled[1, 1](ary, arg1, arg2)
        np.testing.assert_allclose(ary[0], arg1)
        arg1 = np.float16(5.)
        compiled[1, 1](ary, arg1, arg2)
        np.testing.assert_allclose(ary[0], arg2)

    def test_cbrt_f32(self):
        compiled = cuda.jit("void(float32[:], float32)")(simple_cbrt)
        ary = np.zeros(1, dtype=np.float32)
        cbrt_arg = 2.
        compiled[1, 1](ary, cbrt_arg)
        np.testing.assert_allclose(ary[0], cbrt_arg ** (1 / 3))

    def test_cbrt_f64(self):
        compiled = cuda.jit("void(float64[:], float64)")(simple_cbrt)
        ary = np.zeros(1, dtype=np.float64)
        cbrt_arg = 6.
        compiled[1, 1](ary, cbrt_arg)
        np.testing.assert_allclose(ary[0], cbrt_arg ** (1 / 3))

    def test_brev_u4(self):
        compiled = cuda.jit("void(uint32[:], uint32)")(simple_brev)
        ary = np.zeros(1, dtype=np.uint32)
        compiled[1, 1](ary, 0x000030F0)
        self.assertEqual(ary[0], 0x0F0C0000)

    @skip_on_cudasim('only get given a Python "int", assumes 32 bits')
    def test_brev_u8(self):
        compiled = cuda.jit("void(uint64[:], uint64)")(simple_brev)
        ary = np.zeros(1, dtype=np.uint64)
        compiled[1, 1](ary, 0x000030F0000030F0)
        self.assertEqual(ary[0], 0x0F0C00000F0C0000)

    def test_clz_i4(self):
        compiled = cuda.jit("void(int32[:], int32)")(simple_clz)
        ary = np.zeros(1, dtype=np.int32)
        compiled[1, 1](ary, 0x00100000)
        self.assertEqual(ary[0], 11)

    def test_clz_u4(self):
        """
        Although the CUDA Math API
        (http://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__INTRINSIC__INT.html)
        only says int32 & int64 arguments are supported in C code, the LLVM
        IR input supports i8, i16, i32 & i64 (LLVM doesn't have a concept of
        unsigned integers, just unsigned operations on integers).
        http://docs.nvidia.com/cuda/nvvm-ir-spec/index.html#bit-manipulations-intrinics
        """
        compiled = cuda.jit("void(int32[:], uint32)")(simple_clz)
        ary = np.zeros(1, dtype=np.int32)
        compiled[1, 1](ary, 0x00100000)
        self.assertEqual(ary[0], 11)

    def test_clz_i4_1s(self):
        compiled = cuda.jit("void(int32[:], int32)")(simple_clz)
        ary = np.zeros(1, dtype=np.int32)
        compiled[1, 1](ary, 0xFFFFFFFF)
        self.assertEqual(ary[0], 0)

    def test_clz_i4_0s(self):
        compiled = cuda.jit("void(int32[:], int32)")(simple_clz)
        ary = np.zeros(1, dtype=np.int32)
        compiled[1, 1](ary, 0x0)
        self.assertEqual(ary[0], 32, "CUDA semantics")

    @skip_on_cudasim('only get given a Python "int", assumes 32 bits')
    def test_clz_i8(self):
        compiled = cuda.jit("void(int32[:], int64)")(simple_clz)
        ary = np.zeros(1, dtype=np.int32)
        compiled[1, 1](ary, 0x000000000010000)
        self.assertEqual(ary[0], 47)

    def test_ffs_i4(self):
        compiled = cuda.jit("void(int32[:], int32)")(simple_ffs)
        ary = np.zeros(1, dtype=np.int32)
        compiled[1, 1](ary, 0x00100000)
        self.assertEqual(ary[0], 21)
        compiled[1, 1](ary, 0x80000000)
        self.assertEqual(ary[0], 32)

    def test_ffs_u4(self):
        compiled = cuda.jit("void(int32[:], uint32)")(simple_ffs)
        ary = np.zeros(1, dtype=np.int32)
        compiled[1, 1](ary, 0x00100000)
        self.assertEqual(ary[0], 21)
        compiled[1, 1](ary, 0x80000000)
        self.assertEqual(ary[0], 32)

    def test_ffs_i4_1s(self):
        compiled = cuda.jit("void(int32[:], int32)")(simple_ffs)
        ary = np.zeros(1, dtype=np.int32)
        compiled[1, 1](ary, 0xFFFFFFFF)
        self.assertEqual(ary[0], 1)

    def test_ffs_i4_0s(self):
        compiled = cuda.jit("void(int32[:], int32)")(simple_ffs)
        ary = np.zeros(1, dtype=np.int32)
        compiled[1, 1](ary, 0x0)
        self.assertEqual(ary[0], 0)

    @skip_on_cudasim('only get given a Python "int", assumes 32 bits')
    def test_ffs_i8(self):
        compiled = cuda.jit("void(int32[:], int64)")(simple_ffs)
        ary = np.zeros(1, dtype=np.int32)
        compiled[1, 1](ary, 0x000000000010000)
        self.assertEqual(ary[0], 17)
        compiled[1, 1](ary, 0x100000000)
        self.assertEqual(ary[0], 33)

    def test_simple_laneid(self):
        compiled = cuda.jit("void(int32[:])")(simple_laneid)
        count = 2
        ary = np.zeros(count * 32, dtype=np.int32)
        exp = np.tile(np.arange(32, dtype=np.int32), count)
        compiled[1, count * 32](ary)
        self.assertTrue(np.all(ary == exp))

    def test_simple_warpsize(self):
        compiled = cuda.jit("void(int32[:])")(simple_warpsize)
        ary = np.zeros(1, dtype=np.int32)
        compiled[1, 1](ary)
        self.assertEqual(ary[0], 32, "CUDA semantics")

    def test_round_f4(self):
        compiled = cuda.jit("void(int64[:], float32)")(simple_round)
        ary = np.zeros(1, dtype=np.int64)

        for i in [-3.0, -2.5, -2.25, -1.5, 1.5, 2.25, 2.5, 2.75]:
            compiled[1, 1](ary, i)
            self.assertEqual(ary[0], round(i))

    def test_round_f8(self):
        compiled = cuda.jit("void(int64[:], float64)")(simple_round)
        ary = np.zeros(1, dtype=np.int64)

        for i in [-3.0, -2.5, -2.25, -1.5, 1.5, 2.25, 2.5, 2.75]:
            compiled[1, 1](ary, i)
            self.assertEqual(ary[0], round(i))

    def test_round_to_f4(self):
        compiled = cuda.jit("void(float32[:], float32, int32)")(simple_round_to)
        ary = np.zeros(1, dtype=np.float32)
        np.random.seed(123)
        vals = np.random.random(32).astype(np.float32)
        np.concatenate((vals, np.array([np.inf, -np.inf, np.nan])))
        digits = (
            # Common case branch of round_to_impl
            -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5,
            # The algorithm currently implemented can only round to 13 digits
            # with single precision. Note that this doesn't trigger the
            # "overflow safe" branch of the implementation, which can only be
            # hit when using double precision.
            13
        )
        for val, ndigits in itertools.product(vals, digits):
            with self.subTest(val=val, ndigits=ndigits):
                compiled[1, 1](ary, val, ndigits)
                self.assertPreciseEqual(ary[0], round(val, ndigits),
                                        prec='single')

    # CPython on most platforms uses rounding based on dtoa.c, whereas the CUDA
    # round-to implementation uses CPython's fallback implementation, which has
    # slightly different behavior at the edges of the domain. Since the CUDA
    # simulator executes using CPython, we need to skip this test when the
    # simulator is active.
    @skip_on_cudasim('Overflow behavior differs on CPython')
    def test_round_to_f4_overflow(self):
        # Test that the input value is returned when y in round_ndigits
        # overflows.
        compiled = cuda.jit("void(float32[:], float32, int32)")(simple_round_to)
        ary = np.zeros(1, dtype=np.float32)
        val = np.finfo(np.float32).max
        # An unusually large number of digits is required to hit the "y
        # overflows" branch of the implementation because the typing results in
        # the computation of y as float64.
        ndigits = 300
        compiled[1, 1](ary, val, ndigits)
        self.assertEqual(ary[0], val)

    def test_round_to_f4_halfway(self):
        compiled = cuda.jit("void(float32[:], float32, int32)")(simple_round_to)
        ary = np.zeros(1, dtype=np.float32)
        # Value chosen to trigger the "round to even" branch of the
        # implementation
        val = 0.3425
        ndigits = 3
        compiled[1, 1](ary, val, ndigits)
        self.assertPreciseEqual(ary[0], round(val, ndigits), prec='single')

    def test_round_to_f8(self):
        compiled = cuda.jit("void(float64[:], float64, int32)")(simple_round_to)
        ary = np.zeros(1, dtype=np.float64)
        np.random.seed(123)
        vals = np.random.random(32)
        np.concatenate((vals, np.array([np.inf, -np.inf, np.nan])))
        digits = (-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)

        for val, ndigits in itertools.product(vals, digits):
            with self.subTest(val=val, ndigits=ndigits):
                compiled[1, 1](ary, val, ndigits)
                self.assertPreciseEqual(ary[0], round(val, ndigits),
                                        prec='exact')

        # Trigger the "overflow safe" branch of the implementation
        val = 0.12345678987654321 * 10e-15
        ndigits = 23
        with self.subTest(val=val, ndigits=ndigits):
            compiled[1, 1](ary, val, ndigits)
            self.assertPreciseEqual(ary[0], round(val, ndigits),
                                    prec='double')

    # Skipped on cudasim for the same reasons as test_round_to_f4 above.
    @skip_on_cudasim('Overflow behavior differs on CPython')
    def test_round_to_f8_overflow(self):
        # Test that the input value is returned when y in round_ndigits
        # overflows.
        compiled = cuda.jit("void(float64[:], float64, int32)")(simple_round_to)
        ary = np.zeros(1, dtype=np.float64)
        val = np.finfo(np.float64).max
        # Unlike test_round_to_f4_overflow, a reasonable number of digits can
        # be used for this test to overflow y in round_ndigits.
        ndigits = 12
        compiled[1, 1](ary, val, ndigits)
        self.assertEqual(ary[0], val)

    def test_round_to_f8_halfway(self):
        compiled = cuda.jit("void(float64[:], float64, int32)")(simple_round_to)
        ary = np.zeros(1, dtype=np.float64)
        # Value chosen to trigger the "round to even" branch of the
        # implementation, with a value that is not exactly representable with a
        # float32, but only a float64.
        val = 0.5425
        ndigits = 3
        compiled[1, 1](ary, val, ndigits)
        self.assertPreciseEqual(ary[0], round(val, ndigits), prec='double')


if __name__ == '__main__':
    unittest.main()
