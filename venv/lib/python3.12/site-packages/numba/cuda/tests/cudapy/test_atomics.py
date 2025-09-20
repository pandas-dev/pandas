import numpy as np
from textwrap import dedent

from numba import cuda, uint32, uint64, float32, float64
from numba.cuda.testing import unittest, CUDATestCase, cc_X_or_above
from numba.core import config


@cuda.jit(device=True)
def atomic_cast_to_uint64(num):
    return uint64(num)


@cuda.jit(device=True)
def atomic_cast_to_int(num):
    return int(num)


@cuda.jit(device=True)
def atomic_cast_none(num):
    return num


@cuda.jit(device=True)
def atomic_binary_1dim_shared(ary, idx, op2, ary_dtype, ary_nelements,
                              binop_func, cast_func, initializer,
                              neg_idx):
    tid = cuda.threadIdx.x
    sm = cuda.shared.array(ary_nelements, ary_dtype)
    sm[tid] = initializer
    cuda.syncthreads()
    bin = cast_func(idx[tid] % ary_nelements)
    if neg_idx:
        bin = bin % ary_nelements
    binop_func(sm, bin, op2)
    cuda.syncthreads()
    ary[tid] = sm[tid]


@cuda.jit(device=True)
def atomic_binary_1dim_shared2(ary, idx, op2, ary_dtype, ary_nelements,
                               binop_func, cast_func):
    tid = cuda.threadIdx.x
    sm = cuda.shared.array(ary_nelements, ary_dtype)
    sm[tid] = ary[tid]
    cuda.syncthreads()
    bin = cast_func(idx[tid] % ary_nelements)
    binop_func(sm, bin, op2)
    cuda.syncthreads()
    ary[tid] = sm[tid]


@cuda.jit(device=True)
def atomic_binary_2dim_shared(ary, op2, ary_dtype, ary_shape,
                              binop_func, y_cast_func, neg_idx):
    tx = cuda.threadIdx.x
    ty = cuda.threadIdx.y
    sm = cuda.shared.array(ary_shape, ary_dtype)
    sm[tx, ty] = ary[tx, ty]
    cuda.syncthreads()
    bin = (tx, y_cast_func(ty))
    if neg_idx:
        bin = (bin[0] % ary_shape[0], bin[1] % ary_shape[1])
    binop_func(sm, bin, op2)
    cuda.syncthreads()
    ary[tx, ty] = sm[tx, ty]


@cuda.jit(device=True)
def atomic_binary_2dim_global(ary, op2, binop_func, y_cast_func, neg_idx):
    tx = cuda.threadIdx.x
    ty = cuda.threadIdx.y
    bin = (tx, y_cast_func(ty))
    if neg_idx:
        bin = (bin[0] % ary.shape[0], bin[1] % ary.shape[1])
    binop_func(ary, bin, op2)


@cuda.jit(device=True)
def atomic_binary_1dim_global(ary, idx, ary_nelements, op2,
                              binop_func, neg_idx):
    tid = cuda.threadIdx.x
    bin = int(idx[tid] % ary_nelements)
    if neg_idx:
        bin = bin % ary_nelements
    binop_func(ary, bin, op2)


def atomic_add(ary):
    atomic_binary_1dim_shared(ary, ary, 1, uint32, 32,
                              cuda.atomic.add, atomic_cast_none, 0, False)


def atomic_add_wrap(ary):
    atomic_binary_1dim_shared(ary, ary, 1, uint32, 32,
                              cuda.atomic.add, atomic_cast_none, 0, True)


def atomic_add2(ary):
    atomic_binary_2dim_shared(ary, 1, uint32, (4, 8),
                              cuda.atomic.add, atomic_cast_none, False)


def atomic_add2_wrap(ary):
    atomic_binary_2dim_shared(ary, 1, uint32, (4, 8),
                              cuda.atomic.add, atomic_cast_none, True)


def atomic_add3(ary):
    atomic_binary_2dim_shared(ary, 1, uint32, (4, 8),
                              cuda.atomic.add, atomic_cast_to_uint64, False)


def atomic_add_float(ary):
    atomic_binary_1dim_shared(ary, ary, 1.0, float32, 32,
                              cuda.atomic.add, atomic_cast_to_int, 0.0, False)


def atomic_add_float_wrap(ary):
    atomic_binary_1dim_shared(ary, ary, 1.0, float32, 32,
                              cuda.atomic.add, atomic_cast_to_int, 0.0, True)


def atomic_add_float_2(ary):
    atomic_binary_2dim_shared(ary, 1.0, float32, (4, 8),
                              cuda.atomic.add, atomic_cast_none, False)


def atomic_add_float_2_wrap(ary):
    atomic_binary_2dim_shared(ary, 1.0, float32, (4, 8),
                              cuda.atomic.add, atomic_cast_none, True)


def atomic_add_float_3(ary):
    atomic_binary_2dim_shared(ary, 1.0, float32, (4, 8),
                              cuda.atomic.add, atomic_cast_to_uint64, False)


def atomic_add_double_global(idx, ary):
    atomic_binary_1dim_global(ary, idx, 32, 1.0, cuda.atomic.add, False)


def atomic_add_double_global_wrap(idx, ary):
    atomic_binary_1dim_global(ary, idx, 32, 1.0, cuda.atomic.add, True)


def atomic_add_double_global_2(ary):
    atomic_binary_2dim_global(ary, 1, cuda.atomic.add, atomic_cast_none, False)


def atomic_add_double_global_2_wrap(ary):
    atomic_binary_2dim_global(ary, 1, cuda.atomic.add, atomic_cast_none, True)


def atomic_add_double_global_3(ary):
    atomic_binary_2dim_global(ary, 1, cuda.atomic.add, atomic_cast_to_uint64,
                              False)


def atomic_add_double(idx, ary):
    atomic_binary_1dim_shared(ary, idx, 1.0, float64, 32,
                              cuda.atomic.add, atomic_cast_none, 0.0, False)


def atomic_add_double_wrap(idx, ary):
    atomic_binary_1dim_shared(ary, idx, 1.0, float64, 32,
                              cuda.atomic.add, atomic_cast_none, 0.0, True)


def atomic_add_double_2(ary):
    atomic_binary_2dim_shared(ary, 1.0, float64, (4, 8),
                              cuda.atomic.add, atomic_cast_none, False)


def atomic_add_double_2_wrap(ary):
    atomic_binary_2dim_shared(ary, 1.0, float64, (4, 8),
                              cuda.atomic.add, atomic_cast_none, True)


def atomic_add_double_3(ary):
    atomic_binary_2dim_shared(ary, 1.0, float64, (4, 8),
                              cuda.atomic.add, atomic_cast_to_uint64, False)


def atomic_sub(ary):
    atomic_binary_1dim_shared(ary, ary, 1, uint32, 32,
                              cuda.atomic.sub, atomic_cast_none, 0, False)


def atomic_sub2(ary):
    atomic_binary_2dim_shared(ary, 1, uint32, (4, 8),
                              cuda.atomic.sub, atomic_cast_none, False)


def atomic_sub3(ary):
    atomic_binary_2dim_shared(ary, 1, uint32, (4, 8),
                              cuda.atomic.sub, atomic_cast_to_uint64, False)


def atomic_sub_float(ary):
    atomic_binary_1dim_shared(ary, ary, 1.0, float32, 32,
                              cuda.atomic.sub, atomic_cast_to_int, 0.0, False)


def atomic_sub_float_2(ary):
    atomic_binary_2dim_shared(ary, 1.0, float32, (4, 8),
                              cuda.atomic.sub, atomic_cast_none, False)


def atomic_sub_float_3(ary):
    atomic_binary_2dim_shared(ary, 1.0, float32, (4, 8),
                              cuda.atomic.sub, atomic_cast_to_uint64, False)


def atomic_sub_double(idx, ary):
    atomic_binary_1dim_shared(ary, idx, 1.0, float64, 32,
                              cuda.atomic.sub, atomic_cast_none, 0.0, False)


def atomic_sub_double_2(ary):
    atomic_binary_2dim_shared(ary, 1.0, float64, (4, 8),
                              cuda.atomic.sub, atomic_cast_none, False)


def atomic_sub_double_3(ary):
    atomic_binary_2dim_shared(ary, 1.0, float64, (4, 8),
                              cuda.atomic.sub, atomic_cast_to_uint64, False)


def atomic_sub_double_global(idx, ary):
    atomic_binary_1dim_global(ary, idx, 32, 1.0, cuda.atomic.sub, False)


def atomic_sub_double_global_2(ary):
    atomic_binary_2dim_global(ary, 1.0, cuda.atomic.sub, atomic_cast_none,
                              False)


def atomic_sub_double_global_3(ary):
    atomic_binary_2dim_shared(ary, 1.0, float64, (4, 8),
                              cuda.atomic.sub, atomic_cast_to_uint64, False)


def atomic_and(ary, op2):
    atomic_binary_1dim_shared(ary, ary, op2, uint32, 32,
                              cuda.atomic.and_, atomic_cast_none, 1, False)


def atomic_and2(ary, op2):
    atomic_binary_2dim_shared(ary, op2, uint32, (4, 8),
                              cuda.atomic.and_, atomic_cast_none, False)


def atomic_and3(ary, op2):
    atomic_binary_2dim_shared(ary, op2, uint32, (4, 8),
                              cuda.atomic.and_, atomic_cast_to_uint64, False)


def atomic_and_global(idx, ary, op2):
    atomic_binary_1dim_global(ary, idx, 32, op2, cuda.atomic.and_, False)


def atomic_and_global_2(ary, op2):
    atomic_binary_2dim_global(ary, op2, cuda.atomic.and_,
                              atomic_cast_none, False)


def atomic_or(ary, op2):
    atomic_binary_1dim_shared(ary, ary, op2, uint32, 32,
                              cuda.atomic.or_, atomic_cast_none, 0, False)


def atomic_or2(ary, op2):
    atomic_binary_2dim_shared(ary, op2, uint32, (4, 8),
                              cuda.atomic.or_, atomic_cast_none, False)


def atomic_or3(ary, op2):
    atomic_binary_2dim_shared(ary, op2, uint32, (4, 8),
                              cuda.atomic.or_, atomic_cast_to_uint64, False)


def atomic_or_global(idx, ary, op2):
    atomic_binary_1dim_global(ary, idx, 32, op2, cuda.atomic.or_, False)


def atomic_or_global_2(ary, op2):
    atomic_binary_2dim_global(ary, op2, cuda.atomic.or_,
                              atomic_cast_none, False)


def atomic_xor(ary, op2):
    atomic_binary_1dim_shared(ary, ary, op2, uint32, 32,
                              cuda.atomic.xor, atomic_cast_none, 0, False)


def atomic_xor2(ary, op2):
    atomic_binary_2dim_shared(ary, op2, uint32, (4, 8),
                              cuda.atomic.xor, atomic_cast_none, False)


def atomic_xor3(ary, op2):
    atomic_binary_2dim_shared(ary, op2, uint32, (4, 8),
                              cuda.atomic.xor, atomic_cast_to_uint64, False)


def atomic_xor_global(idx, ary, op2):
    atomic_binary_1dim_global(ary, idx, 32, op2, cuda.atomic.xor, False)


def atomic_xor_global_2(ary, op2):
    atomic_binary_2dim_global(ary, op2, cuda.atomic.xor,
                              atomic_cast_none, False)


def atomic_inc32(ary, idx, op2):
    atomic_binary_1dim_shared2(ary, idx, op2, uint32, 32,
                               cuda.atomic.inc, atomic_cast_none)


def atomic_inc64(ary, idx, op2):
    atomic_binary_1dim_shared2(ary, idx, op2, uint64, 32,
                               cuda.atomic.inc, atomic_cast_to_int)


def atomic_inc2_32(ary, op2):
    atomic_binary_2dim_shared(ary, op2, uint32, (4, 8),
                              cuda.atomic.inc, atomic_cast_none, False)


def atomic_inc2_64(ary, op2):
    atomic_binary_2dim_shared(ary, op2, uint64, (4, 8),
                              cuda.atomic.inc, atomic_cast_none, False)


def atomic_inc3(ary, op2):
    atomic_binary_2dim_shared(ary, op2, uint32, (4, 8),
                              cuda.atomic.inc, atomic_cast_to_uint64, False)


def atomic_inc_global(idx, ary, op2):
    atomic_binary_1dim_global(ary, idx, 32, op2, cuda.atomic.inc, False)


def atomic_inc_global_2(ary, op2):
    atomic_binary_2dim_global(ary, op2, cuda.atomic.inc,
                              atomic_cast_none, False)


def atomic_dec32(ary, idx, op2):
    atomic_binary_1dim_shared2(ary, idx, op2, uint32, 32,
                               cuda.atomic.dec, atomic_cast_none)


def atomic_dec64(ary, idx, op2):
    atomic_binary_1dim_shared2(ary, idx, op2, uint64, 32,
                               cuda.atomic.dec, atomic_cast_to_int)


def atomic_dec2_32(ary, op2):
    atomic_binary_2dim_shared(ary, op2, uint32, (4, 8),
                              cuda.atomic.dec, atomic_cast_none, False)


def atomic_dec2_64(ary, op2):
    atomic_binary_2dim_shared(ary, op2, uint64, (4, 8),
                              cuda.atomic.dec, atomic_cast_none, False)


def atomic_dec3(ary, op2):
    atomic_binary_2dim_shared(ary, op2, uint32, (4, 8),
                              cuda.atomic.dec, atomic_cast_to_uint64, False)


def atomic_dec_global(idx, ary, op2):
    atomic_binary_1dim_global(ary, idx, 32, op2, cuda.atomic.dec, False)


def atomic_dec_global_2(ary, op2):
    atomic_binary_2dim_global(ary, op2, cuda.atomic.dec,
                              atomic_cast_none, False)


def atomic_exch(ary, idx, op2):
    atomic_binary_1dim_shared2(ary, idx, op2, uint32, 32,
                               cuda.atomic.exch, atomic_cast_none)


def atomic_exch2(ary, op2):
    atomic_binary_2dim_shared(ary, op2, uint32, (4, 8),
                              cuda.atomic.exch, atomic_cast_none, False)


def atomic_exch3(ary, op2):
    atomic_binary_2dim_shared(ary, op2, uint64, (4, 8),
                              cuda.atomic.exch, atomic_cast_none, False)


def atomic_exch_global(idx, ary, op2):
    atomic_binary_1dim_global(ary, idx, 32, op2, cuda.atomic.exch, False)


def gen_atomic_extreme_funcs(func):

    fns = dedent("""
    def atomic(res, ary):
        tx = cuda.threadIdx.x
        bx = cuda.blockIdx.x
        {func}(res, 0, ary[tx, bx])

    def atomic_double_normalizedindex(res, ary):
        tx = cuda.threadIdx.x
        bx = cuda.blockIdx.x
        {func}(res, 0, ary[tx, uint64(bx)])

    def atomic_double_oneindex(res, ary):
        tx = cuda.threadIdx.x
        {func}(res, 0, ary[tx])

    def atomic_double_shared(res, ary):
        tid = cuda.threadIdx.x
        smary = cuda.shared.array(32, float64)
        smary[tid] = ary[tid]
        smres = cuda.shared.array(1, float64)
        if tid == 0:
            smres[0] = res[0]
        cuda.syncthreads()
        {func}(smres, 0, smary[tid])
        cuda.syncthreads()
        if tid == 0:
            res[0] = smres[0]
    """).format(func=func)
    ld = {}
    exec(fns, {'cuda': cuda, 'float64': float64, 'uint64': uint64}, ld)
    return (ld['atomic'], ld['atomic_double_normalizedindex'],
            ld['atomic_double_oneindex'], ld['atomic_double_shared'])


(atomic_max, atomic_max_double_normalizedindex, atomic_max_double_oneindex,
 atomic_max_double_shared) = gen_atomic_extreme_funcs('cuda.atomic.max')
(atomic_min, atomic_min_double_normalizedindex, atomic_min_double_oneindex,
 atomic_min_double_shared) = gen_atomic_extreme_funcs('cuda.atomic.min')
(atomic_nanmax, atomic_nanmax_double_normalizedindex,
 atomic_nanmax_double_oneindex, atomic_nanmax_double_shared) = \
    gen_atomic_extreme_funcs('cuda.atomic.nanmax')
(atomic_nanmin, atomic_nanmin_double_normalizedindex,
 atomic_nanmin_double_oneindex, atomic_nanmin_double_shared) = \
    gen_atomic_extreme_funcs('cuda.atomic.nanmin')


def atomic_compare_and_swap(res, old, ary, fill_val):
    gid = cuda.grid(1)
    if gid < res.size:
        old[gid] = cuda.atomic.compare_and_swap(res[gid:], fill_val, ary[gid])


def atomic_cas_1dim(res, old, ary, fill_val):
    gid = cuda.grid(1)
    if gid < res.size:
        old[gid] = cuda.atomic.cas(res, gid, fill_val, ary[gid])


def atomic_cas_2dim(res, old, ary, fill_val):
    gid = cuda.grid(2)
    if gid[0] < res.shape[0] and gid[1] < res.shape[1]:
        old[gid] = cuda.atomic.cas(res, gid, fill_val, ary[gid])


class TestCudaAtomics(CUDATestCase):
    def setUp(self):
        super().setUp()
        np.random.seed(0)

    def test_atomic_add(self):
        ary = np.random.randint(0, 32, size=32).astype(np.uint32)
        ary_wrap = ary.copy()
        orig = ary.copy()

        cuda_atomic_add = cuda.jit('void(uint32[:])')(atomic_add)
        cuda_atomic_add[1, 32](ary)

        cuda_atomic_add_wrap = cuda.jit('void(uint32[:])')(atomic_add_wrap)
        cuda_atomic_add_wrap[1, 32](ary_wrap)

        gold = np.zeros(32, dtype=np.uint32)
        for i in range(orig.size):
            gold[orig[i]] += 1

        self.assertTrue(np.all(ary == gold))
        self.assertTrue(np.all(ary_wrap == gold))

    def test_atomic_add2(self):
        ary = np.random.randint(0, 32, size=32).astype(np.uint32).reshape(4, 8)
        ary_wrap = ary.copy()
        orig = ary.copy()

        cuda_atomic_add2 = cuda.jit('void(uint32[:,:])')(atomic_add2)
        cuda_atomic_add2[1, (4, 8)](ary)

        cuda_atomic_add2_wrap = cuda.jit('void(uint32[:,:])')(atomic_add2_wrap)
        cuda_atomic_add2_wrap[1, (4, 8)](ary_wrap)

        self.assertTrue(np.all(ary == orig + 1))
        self.assertTrue(np.all(ary_wrap == orig + 1))

    def test_atomic_add3(self):
        ary = np.random.randint(0, 32, size=32).astype(np.uint32).reshape(4, 8)
        orig = ary.copy()
        cuda_atomic_add3 = cuda.jit('void(uint32[:,:])')(atomic_add3)
        cuda_atomic_add3[1, (4, 8)](ary)

        self.assertTrue(np.all(ary == orig + 1))

    def test_atomic_add_float(self):
        ary = np.random.randint(0, 32, size=32).astype(np.float32)
        ary_wrap = ary.copy()
        orig = ary.copy().astype(np.intp)

        cuda_atomic_add_float = cuda.jit('void(float32[:])')(atomic_add_float)
        cuda_atomic_add_float[1, 32](ary)

        add_float_wrap = cuda.jit('void(float32[:])')(atomic_add_float_wrap)
        add_float_wrap[1, 32](ary_wrap)

        gold = np.zeros(32, dtype=np.uint32)
        for i in range(orig.size):
            gold[orig[i]] += 1.0

        self.assertTrue(np.all(ary == gold))
        self.assertTrue(np.all(ary_wrap == gold))

    def test_atomic_add_float_2(self):
        ary = np.random.randint(0, 32, size=32).astype(np.float32).reshape(4, 8)
        ary_wrap = ary.copy()
        orig = ary.copy()

        cuda_atomic_add2 = cuda.jit('void(float32[:,:])')(atomic_add_float_2)
        cuda_atomic_add2[1, (4, 8)](ary)

        cuda_func_wrap = cuda.jit('void(float32[:,:])')(atomic_add_float_2_wrap)
        cuda_func_wrap[1, (4, 8)](ary_wrap)

        self.assertTrue(np.all(ary == orig + 1))
        self.assertTrue(np.all(ary_wrap == orig + 1))

    def test_atomic_add_float_3(self):
        ary = np.random.randint(0, 32, size=32).astype(np.float32).reshape(4, 8)
        orig = ary.copy()
        cuda_atomic_add3 = cuda.jit('void(float32[:,:])')(atomic_add_float_3)
        cuda_atomic_add3[1, (4, 8)](ary)

        self.assertTrue(np.all(ary == orig + 1))

    def assertCorrectFloat64Atomics(self, kernel, shared=True):
        if config.ENABLE_CUDASIM:
            return

        # Use the first (and only) definition
        asm = next(iter(kernel.inspect_asm().values()))
        if cc_X_or_above(6, 0):
            if cuda.runtime.get_version() > (12, 1):
                # CUDA 12.2 and above generate a more optimized reduction
                # instruction, because the result does not need to be
                # placed in a register.
                inst = 'red'
            else:
                inst = 'atom'

            if shared:
                inst = f'{inst}.shared'

            self.assertIn(f'{inst}.add.f64', asm)
        else:
            if shared:
                self.assertIn('atom.shared.cas.b64', asm)
            else:
                self.assertIn('atom.cas.b64', asm)

    def test_atomic_add_double(self):
        idx = np.random.randint(0, 32, size=32, dtype=np.int64)
        ary = np.zeros(32, np.float64)
        ary_wrap = ary.copy()

        cuda_fn = cuda.jit('void(int64[:], float64[:])')(atomic_add_double)
        cuda_fn[1, 32](idx, ary)

        wrap_fn = cuda.jit('void(int64[:], float64[:])')(atomic_add_double_wrap)
        wrap_fn[1, 32](idx, ary_wrap)

        gold = np.zeros(32, dtype=np.uint32)
        for i in range(idx.size):
            gold[idx[i]] += 1.0

        np.testing.assert_equal(ary, gold)
        np.testing.assert_equal(ary_wrap, gold)
        self.assertCorrectFloat64Atomics(cuda_fn)
        self.assertCorrectFloat64Atomics(wrap_fn)

    def test_atomic_add_double_2(self):
        ary = np.random.randint(0, 32, size=32).astype(np.float64).reshape(4, 8)
        ary_wrap = ary.copy()
        orig = ary.copy()

        cuda_fn = cuda.jit('void(float64[:,:])')(atomic_add_double_2)
        cuda_fn[1, (4, 8)](ary)

        cuda_fn_wrap = cuda.jit('void(float64[:,:])')(atomic_add_double_2_wrap)
        cuda_fn_wrap[1, (4, 8)](ary_wrap)

        np.testing.assert_equal(ary, orig + 1)
        np.testing.assert_equal(ary_wrap, orig + 1)
        self.assertCorrectFloat64Atomics(cuda_fn)
        self.assertCorrectFloat64Atomics(cuda_fn_wrap)

    def test_atomic_add_double_3(self):
        ary = np.random.randint(0, 32, size=32).astype(np.float64).reshape(4, 8)
        orig = ary.copy()
        cuda_func = cuda.jit('void(float64[:,:])')(atomic_add_double_3)
        cuda_func[1, (4, 8)](ary)

        np.testing.assert_equal(ary, orig + 1)
        self.assertCorrectFloat64Atomics(cuda_func)

    def test_atomic_add_double_global(self):
        idx = np.random.randint(0, 32, size=32, dtype=np.int64)
        ary = np.zeros(32, np.float64)
        ary_wrap = ary.copy()

        sig = 'void(int64[:], float64[:])'
        cuda_func = cuda.jit(sig)(atomic_add_double_global)
        wrap_cuda_func = cuda.jit(sig)(atomic_add_double_global_wrap)

        cuda_func[1, 32](idx, ary)
        wrap_cuda_func[1, 32](idx, ary_wrap)

        gold = np.zeros(32, dtype=np.uint32)
        for i in range(idx.size):
            gold[idx[i]] += 1.0

        np.testing.assert_equal(ary, gold)
        np.testing.assert_equal(ary_wrap, gold)
        self.assertCorrectFloat64Atomics(cuda_func, shared=False)
        self.assertCorrectFloat64Atomics(wrap_cuda_func, shared=False)

    def test_atomic_add_double_global_2(self):
        ary = np.random.randint(0, 32, size=32).astype(np.float64).reshape(4, 8)
        ary_wrap = ary.copy()
        orig = ary.copy()

        sig = 'void(float64[:,:])'
        cuda_func = cuda.jit(sig)(atomic_add_double_global_2)
        wrap_cuda_func = cuda.jit(sig)(atomic_add_double_global_2_wrap)

        cuda_func[1, (4, 8)](ary)
        wrap_cuda_func[1, (4, 8)](ary_wrap)

        np.testing.assert_equal(ary, orig + 1)
        np.testing.assert_equal(ary_wrap, orig + 1)
        self.assertCorrectFloat64Atomics(cuda_func, shared=False)
        self.assertCorrectFloat64Atomics(wrap_cuda_func, shared=False)

    def test_atomic_add_double_global_3(self):
        ary = np.random.randint(0, 32, size=32).astype(np.float64).reshape(4, 8)
        orig = ary.copy()
        cuda_func = cuda.jit('void(float64[:,:])')(atomic_add_double_global_3)
        cuda_func[1, (4, 8)](ary)

        np.testing.assert_equal(ary, orig + 1)
        self.assertCorrectFloat64Atomics(cuda_func, shared=False)

    def test_atomic_sub(self):
        ary = np.random.randint(0, 32, size=32).astype(np.uint32)
        orig = ary.copy()
        cuda_atomic_sub = cuda.jit('void(uint32[:])')(atomic_sub)
        cuda_atomic_sub[1, 32](ary)

        gold = np.zeros(32, dtype=np.uint32)
        for i in range(orig.size):
            gold[orig[i]] -= 1

        self.assertTrue(np.all(ary == gold))

    def test_atomic_sub2(self):
        ary = np.random.randint(0, 32, size=32).astype(np.uint32).reshape(4, 8)
        orig = ary.copy()
        cuda_atomic_sub2 = cuda.jit('void(uint32[:,:])')(atomic_sub2)
        cuda_atomic_sub2[1, (4, 8)](ary)
        self.assertTrue(np.all(ary == orig - 1))

    def test_atomic_sub3(self):
        ary = np.random.randint(0, 32, size=32).astype(np.uint32).reshape(4, 8)
        orig = ary.copy()
        cuda_atomic_sub3 = cuda.jit('void(uint32[:,:])')(atomic_sub3)
        cuda_atomic_sub3[1, (4, 8)](ary)
        self.assertTrue(np.all(ary == orig - 1))

    def test_atomic_sub_float(self):
        ary = np.random.randint(0, 32, size=32).astype(np.float32)
        orig = ary.copy().astype(np.intp)
        cuda_atomic_sub_float = cuda.jit('void(float32[:])')(atomic_sub_float)
        cuda_atomic_sub_float[1, 32](ary)

        gold = np.zeros(32, dtype=np.float32)
        for i in range(orig.size):
            gold[orig[i]] -= 1.0

        self.assertTrue(np.all(ary == gold))

    def test_atomic_sub_float_2(self):
        ary = np.random.randint(0, 32, size=32).astype(np.float32).reshape(4, 8)
        orig = ary.copy()
        cuda_atomic_sub2 = cuda.jit('void(float32[:,:])')(atomic_sub_float_2)
        cuda_atomic_sub2[1, (4, 8)](ary)
        self.assertTrue(np.all(ary == orig - 1))

    def test_atomic_sub_float_3(self):
        ary = np.random.randint(0, 32, size=32).astype(np.float32).reshape(4, 8)
        orig = ary.copy()
        cuda_atomic_sub3 = cuda.jit('void(float32[:,:])')(atomic_sub_float_3)
        cuda_atomic_sub3[1, (4, 8)](ary)
        self.assertTrue(np.all(ary == orig - 1))

    def test_atomic_sub_double(self):
        idx = np.random.randint(0, 32, size=32, dtype=np.int64)
        ary = np.zeros(32, np.float64)
        cuda_func = cuda.jit('void(int64[:], float64[:])')(atomic_sub_double)
        cuda_func[1, 32](idx, ary)

        gold = np.zeros(32, dtype=np.float64)
        for i in range(idx.size):
            gold[idx[i]] -= 1.0

        np.testing.assert_equal(ary, gold)

    def test_atomic_sub_double_2(self):
        ary = np.random.randint(0, 32, size=32).astype(np.float64).reshape(4, 8)
        orig = ary.copy()
        cuda_func = cuda.jit('void(float64[:,:])')(atomic_sub_double_2)
        cuda_func[1, (4, 8)](ary)
        np.testing.assert_equal(ary, orig - 1)

    def test_atomic_sub_double_3(self):
        ary = np.random.randint(0, 32, size=32).astype(np.float64).reshape(4, 8)
        orig = ary.copy()
        cuda_func = cuda.jit('void(float64[:,:])')(atomic_sub_double_3)
        cuda_func[1, (4, 8)](ary)
        np.testing.assert_equal(ary, orig - 1)

    def test_atomic_sub_double_global(self):
        idx = np.random.randint(0, 32, size=32, dtype=np.int64)
        ary = np.zeros(32, np.float64)
        sig = 'void(int64[:], float64[:])'
        cuda_func = cuda.jit(sig)(atomic_sub_double_global)
        cuda_func[1, 32](idx, ary)

        gold = np.zeros(32, dtype=np.float64)
        for i in range(idx.size):
            gold[idx[i]] -= 1.0

        np.testing.assert_equal(ary, gold)

    def test_atomic_sub_double_global_2(self):
        ary = np.random.randint(0, 32, size=32).astype(np.float64).reshape(4, 8)
        orig = ary.copy()
        cuda_func = cuda.jit('void(float64[:,:])')(atomic_sub_double_global_2)
        cuda_func[1, (4, 8)](ary)
        np.testing.assert_equal(ary, orig - 1)

    def test_atomic_sub_double_global_3(self):
        ary = np.random.randint(0, 32, size=32).astype(np.float64).reshape(4, 8)
        orig = ary.copy()
        cuda_func = cuda.jit('void(float64[:,:])')(atomic_sub_double_global_3)
        cuda_func[1, (4, 8)](ary)
        np.testing.assert_equal(ary, orig - 1)

    def test_atomic_and(self):
        rand_const = np.random.randint(500)
        ary = np.random.randint(0, 32, size=32).astype(np.uint32)
        orig = ary.copy()
        cuda_func = cuda.jit('void(uint32[:], uint32)')(atomic_and)
        cuda_func[1, 32](ary, rand_const)

        gold = ary.copy()
        for i in range(orig.size):
            gold[orig[i]] &= rand_const

        self.assertTrue(np.all(ary == gold))

    def test_atomic_and2(self):
        rand_const = np.random.randint(500)
        ary = np.random.randint(0, 32, size=32).astype(np.uint32).reshape(4, 8)
        orig = ary.copy()
        cuda_atomic_and2 = cuda.jit('void(uint32[:,:], uint32)')(atomic_and2)
        cuda_atomic_and2[1, (4, 8)](ary, rand_const)
        self.assertTrue(np.all(ary == orig & rand_const))

    def test_atomic_and3(self):
        rand_const = np.random.randint(500)
        ary = np.random.randint(0, 32, size=32).astype(np.uint32).reshape(4, 8)
        orig = ary.copy()
        cuda_atomic_and3 = cuda.jit('void(uint32[:,:], uint32)')(atomic_and3)
        cuda_atomic_and3[1, (4, 8)](ary, rand_const)
        self.assertTrue(np.all(ary == orig & rand_const))

    def test_atomic_and_global(self):
        rand_const = np.random.randint(500)
        idx = np.random.randint(0, 32, size=32, dtype=np.int32)
        ary = np.random.randint(0, 32, size=32, dtype=np.int32)
        sig = 'void(int32[:], int32[:], int32)'
        cuda_func = cuda.jit(sig)(atomic_and_global)
        cuda_func[1, 32](idx, ary, rand_const)

        gold = ary.copy()
        for i in range(idx.size):
            gold[idx[i]] &= rand_const

        np.testing.assert_equal(ary, gold)

    def test_atomic_and_global_2(self):
        rand_const = np.random.randint(500)
        ary = np.random.randint(0, 32, size=32).astype(np.uint32).reshape(4, 8)
        orig = ary.copy()
        cuda_func = cuda.jit('void(uint32[:,:], uint32)')(atomic_and_global_2)
        cuda_func[1, (4, 8)](ary, rand_const)
        np.testing.assert_equal(ary, orig & rand_const)

    def test_atomic_or(self):
        rand_const = np.random.randint(500)
        ary = np.random.randint(0, 32, size=32).astype(np.uint32)
        orig = ary.copy()
        cuda_func = cuda.jit('void(uint32[:], uint32)')(atomic_or)
        cuda_func[1, 32](ary, rand_const)

        gold = np.zeros(32, dtype=np.uint32)
        for i in range(orig.size):
            gold[orig[i]] |= rand_const

        self.assertTrue(np.all(ary == gold))

    def test_atomic_or2(self):
        rand_const = np.random.randint(500)
        ary = np.random.randint(0, 32, size=32).astype(np.uint32).reshape(4, 8)
        orig = ary.copy()
        cuda_atomic_and2 = cuda.jit('void(uint32[:,:], uint32)')(atomic_or2)
        cuda_atomic_and2[1, (4, 8)](ary, rand_const)
        self.assertTrue(np.all(ary == orig | rand_const))

    def test_atomic_or3(self):
        rand_const = np.random.randint(500)
        ary = np.random.randint(0, 32, size=32).astype(np.uint32).reshape(4, 8)
        orig = ary.copy()
        cuda_atomic_and3 = cuda.jit('void(uint32[:,:], uint32)')(atomic_or3)
        cuda_atomic_and3[1, (4, 8)](ary, rand_const)
        self.assertTrue(np.all(ary == orig | rand_const))

    def test_atomic_or_global(self):
        rand_const = np.random.randint(500)
        idx = np.random.randint(0, 32, size=32, dtype=np.int32)
        ary = np.random.randint(0, 32, size=32, dtype=np.int32)
        sig = 'void(int32[:], int32[:], int32)'
        cuda_func = cuda.jit(sig)(atomic_or_global)
        cuda_func[1, 32](idx, ary, rand_const)

        gold = ary.copy()
        for i in range(idx.size):
            gold[idx[i]] |= rand_const

        np.testing.assert_equal(ary, gold)

    def test_atomic_or_global_2(self):
        rand_const = np.random.randint(500)
        ary = np.random.randint(0, 32, size=32).astype(np.uint32).reshape(4, 8)
        orig = ary.copy()
        cuda_func = cuda.jit('void(uint32[:,:], uint32)')(atomic_or_global_2)
        cuda_func[1, (4, 8)](ary, rand_const)
        np.testing.assert_equal(ary, orig | rand_const)

    def test_atomic_xor(self):
        rand_const = np.random.randint(500)
        ary = np.random.randint(0, 32, size=32).astype(np.uint32)
        orig = ary.copy()
        cuda_func = cuda.jit('void(uint32[:], uint32)')(atomic_xor)
        cuda_func[1, 32](ary, rand_const)

        gold = np.zeros(32, dtype=np.uint32)
        for i in range(orig.size):
            gold[orig[i]] ^= rand_const

        self.assertTrue(np.all(ary == gold))

    def test_atomic_xor2(self):
        rand_const = np.random.randint(500)
        ary = np.random.randint(0, 32, size=32).astype(np.uint32).reshape(4, 8)
        orig = ary.copy()
        cuda_atomic_xor2 = cuda.jit('void(uint32[:,:], uint32)')(atomic_xor2)
        cuda_atomic_xor2[1, (4, 8)](ary, rand_const)
        self.assertTrue(np.all(ary == orig ^ rand_const))

    def test_atomic_xor3(self):
        rand_const = np.random.randint(500)
        ary = np.random.randint(0, 32, size=32).astype(np.uint32).reshape(4, 8)
        orig = ary.copy()
        cuda_atomic_xor3 = cuda.jit('void(uint32[:,:], uint32)')(atomic_xor3)
        cuda_atomic_xor3[1, (4, 8)](ary, rand_const)
        self.assertTrue(np.all(ary == orig ^ rand_const))

    def test_atomic_xor_global(self):
        rand_const = np.random.randint(500)
        idx = np.random.randint(0, 32, size=32, dtype=np.int32)
        ary = np.random.randint(0, 32, size=32, dtype=np.int32)
        gold = ary.copy()
        sig = 'void(int32[:], int32[:], int32)'
        cuda_func = cuda.jit(sig)(atomic_xor_global)
        cuda_func[1, 32](idx, ary, rand_const)

        for i in range(idx.size):
            gold[idx[i]] ^= rand_const

        np.testing.assert_equal(ary, gold)

    def test_atomic_xor_global_2(self):
        rand_const = np.random.randint(500)
        ary = np.random.randint(0, 32, size=32).astype(np.uint32).reshape(4, 8)
        orig = ary.copy()
        cuda_func = cuda.jit('void(uint32[:,:], uint32)')(atomic_xor_global_2)
        cuda_func[1, (4, 8)](ary, rand_const)
        np.testing.assert_equal(ary, orig ^ rand_const)

    def inc_dec_1dim_setup(self, dtype):
        rconst = np.random.randint(32,  dtype=dtype)
        rary = np.random.randint(0, 32, size=32).astype(dtype)
        ary_idx = np.arange(32, dtype=dtype)
        return rconst, rary, ary_idx

    def inc_dec_2dim_setup(self, dtype):
        rconst = np.random.randint(32, dtype=dtype)
        rary = np.random.randint(0, 32, size=32).astype(dtype).reshape(4, 8)
        return rconst, rary

    def check_inc_index(self, ary, idx, rconst, sig, nblocks, blksize, func):
        orig = ary.copy()
        cuda_func = cuda.jit(sig)(func)
        cuda_func[nblocks, blksize](ary, idx, rconst)
        np.testing.assert_equal(ary, np.where(orig >= rconst, 0, orig + 1))

    def check_inc_index2(self, ary, idx, rconst, sig, nblocks, blksize, func):
        orig = ary.copy()
        cuda_func = cuda.jit(sig)(func)
        cuda_func[nblocks, blksize](idx, ary, rconst)
        np.testing.assert_equal(ary, np.where(orig >= rconst, 0, orig + 1))

    def check_inc(self, ary, rconst, sig, nblocks, blksize, func):
        orig = ary.copy()
        cuda_func = cuda.jit(sig)(func)
        cuda_func[nblocks, blksize](ary, rconst)
        np.testing.assert_equal(ary, np.where(orig >= rconst, 0, orig + 1))

    def test_atomic_inc_32(self):
        rand_const, ary, idx = self.inc_dec_1dim_setup(dtype=np.uint32)
        sig = 'void(uint32[:], uint32[:], uint32)'
        self.check_inc_index(ary, idx, rand_const, sig, 1, 32, atomic_inc32)

    def test_atomic_inc_64(self):
        rand_const, ary, idx = self.inc_dec_1dim_setup(dtype=np.uint64)
        sig = 'void(uint64[:], uint64[:], uint64)'
        self.check_inc_index(ary, idx, rand_const, sig, 1, 32, atomic_inc64)

    def test_atomic_inc2_32(self):
        rand_const, ary = self.inc_dec_2dim_setup(np.uint32)
        sig = 'void(uint32[:,:], uint32)'
        self.check_inc(ary, rand_const, sig, 1, (4,8), atomic_inc2_32)

    def test_atomic_inc2_64(self):
        rand_const, ary = self.inc_dec_2dim_setup(np.uint64)
        sig = 'void(uint64[:,:], uint64)'
        self.check_inc(ary, rand_const, sig, 1, (4,8), atomic_inc2_64)

    def test_atomic_inc3(self):
        rand_const, ary = self.inc_dec_2dim_setup(np.uint32)
        sig = 'void(uint32[:,:], uint32)'
        self.check_inc(ary, rand_const, sig, 1, (4,8), atomic_inc3)

    def test_atomic_inc_global_32(self):
        rand_const, ary, idx = self.inc_dec_1dim_setup(dtype=np.uint32)
        sig = 'void(uint32[:], uint32[:], uint32)'
        self.check_inc_index2(ary, idx, rand_const, sig, 1, 32,
                              atomic_inc_global)

    def test_atomic_inc_global_64(self):
        rand_const, ary, idx = self.inc_dec_1dim_setup(dtype=np.uint64)
        sig = 'void(uint64[:], uint64[:], uint64)'
        self.check_inc_index2(ary, idx, rand_const, sig, 1, 32,
                              atomic_inc_global)

    def test_atomic_inc_global_2_32(self):
        rand_const, ary = self.inc_dec_2dim_setup(np.uint32)
        sig = 'void(uint32[:,:], uint32)'
        self.check_inc(ary, rand_const, sig, 1, (4,8), atomic_inc_global_2)

    def test_atomic_inc_global_2_64(self):
        rand_const, ary = self.inc_dec_2dim_setup(np.uint64)
        sig = 'void(uint64[:,:], uint64)'
        self.check_inc(ary, rand_const, sig, 1, (4,8), atomic_inc_global_2)

    def check_dec_index(self, ary, idx, rconst, sig, nblocks, blksize, func):
        orig = ary.copy()
        cuda_func = cuda.jit(sig)(func)
        cuda_func[nblocks, blksize](ary, idx, rconst)
        np.testing.assert_equal(ary, np.where(orig == 0, rconst,
                                              np.where(orig > rconst,
                                                       rconst,
                                                       orig - 1)))

    def check_dec_index2(self, ary, idx, rconst, sig, nblocks, blksize, func):
        orig = ary.copy()
        cuda_func = cuda.jit(sig)(func)
        cuda_func[nblocks, blksize](idx, ary, rconst)
        np.testing.assert_equal(ary, np.where(orig == 0, rconst,
                                              np.where(orig > rconst,
                                                       rconst,
                                                       orig - 1)))

    def check_dec(self, ary, rconst, sig, nblocks, blksize, func):
        orig = ary.copy()
        cuda_func = cuda.jit(sig)(func)
        cuda_func[nblocks, blksize](ary, rconst)
        np.testing.assert_equal(ary, np.where(orig == 0, rconst,
                                              np.where(orig > rconst,
                                                       rconst,
                                                       orig - 1)))

    def test_atomic_dec_32(self):
        rand_const, ary, idx = self.inc_dec_1dim_setup(dtype=np.uint32)
        sig = 'void(uint32[:], uint32[:], uint32)'
        self.check_dec_index(ary, idx, rand_const, sig, 1, 32, atomic_dec32)

    def test_atomic_dec_64(self):
        rand_const, ary, idx = self.inc_dec_1dim_setup(dtype=np.uint64)
        sig = 'void(uint64[:], uint64[:], uint64)'
        self.check_dec_index(ary, idx, rand_const, sig, 1, 32, atomic_dec64)

    def test_atomic_dec2_32(self):
        rand_const, ary = self.inc_dec_2dim_setup(np.uint32)
        sig = 'void(uint32[:,:], uint32)'
        self.check_dec(ary, rand_const, sig, 1, (4,8), atomic_dec2_32)

    def test_atomic_dec2_64(self):
        rand_const, ary = self.inc_dec_2dim_setup(np.uint64)
        sig = 'void(uint64[:,:], uint64)'
        self.check_dec(ary, rand_const, sig, 1, (4,8), atomic_dec2_64)

    def test_atomic_dec3_new(self):
        rand_const, ary = self.inc_dec_2dim_setup(np.uint32)
        sig = 'void(uint32[:,:], uint32)'
        self.check_dec(ary, rand_const, sig, 1, (4,8), atomic_dec3)

    def test_atomic_dec_global_32(self):
        rand_const, ary, idx = self.inc_dec_1dim_setup(dtype=np.uint32)
        sig = 'void(uint32[:], uint32[:], uint32)'
        self.check_dec_index2(ary, idx, rand_const, sig, 1, 32,
                              atomic_dec_global)

    def test_atomic_dec_global_64(self):
        rand_const, ary, idx = self.inc_dec_1dim_setup(dtype=np.uint64)
        sig = 'void(uint64[:], uint64[:], uint64)'
        self.check_dec_index2(ary, idx, rand_const, sig, 1, 32,
                              atomic_dec_global)

    def test_atomic_dec_global2_32(self):
        rand_const, ary = self.inc_dec_2dim_setup(np.uint32)
        sig = 'void(uint32[:,:], uint32)'
        self.check_dec(ary, rand_const, sig, 1, (4,8), atomic_dec_global_2)

    def test_atomic_dec_global2_64(self):
        rand_const, ary = self.inc_dec_2dim_setup(np.uint64)
        sig = 'void(uint64[:,:], uint64)'
        self.check_dec(ary, rand_const, sig, 1, (4,8), atomic_dec_global_2)

    def test_atomic_exch(self):
        rand_const = np.random.randint(50, 100, dtype=np.uint32)
        ary = np.random.randint(0, 32, size=32).astype(np.uint32)
        idx = np.arange(32, dtype=np.uint32)

        cuda_func = cuda.jit('void(uint32[:], uint32[:], uint32)')(atomic_exch)
        cuda_func[1, 32](ary, idx, rand_const)

        np.testing.assert_equal(ary, rand_const)

    def test_atomic_exch2(self):
        rand_const = np.random.randint(50, 100, dtype=np.uint32)
        ary = np.random.randint(0, 32, size=32).astype(np.uint32).reshape(4, 8)

        cuda_func = cuda.jit('void(uint32[:,:], uint32)')(atomic_exch2)
        cuda_func[1, (4, 8)](ary, rand_const)
        np.testing.assert_equal(ary, rand_const)

    def test_atomic_exch3(self):
        rand_const = np.random.randint(50, 100, dtype=np.uint64)
        ary = np.random.randint(0, 32, size=32).astype(np.uint64).reshape(4, 8)

        cuda_func = cuda.jit('void(uint64[:,:], uint64)')(atomic_exch3)
        cuda_func[1, (4, 8)](ary, rand_const)
        np.testing.assert_equal(ary, rand_const)

    def test_atomic_exch_global(self):
        rand_const = np.random.randint(50, 100, dtype=np.uint32)
        idx = np.arange(32, dtype=np.uint32)
        ary = np.random.randint(0, 32, size=32, dtype=np.uint32)

        sig = 'void(uint32[:], uint32[:], uint32)'
        cuda_func = cuda.jit(sig)(atomic_exch_global)
        cuda_func[1, 32](idx, ary, rand_const)
        np.testing.assert_equal(ary, rand_const)

    def check_atomic_max(self, dtype, lo, hi):
        vals = np.random.randint(lo, hi, size=(32, 32)).astype(dtype)
        res = np.zeros(1, dtype=vals.dtype)
        cuda_func = cuda.jit(atomic_max)
        cuda_func[32, 32](res, vals)
        gold = np.max(vals)
        np.testing.assert_equal(res, gold)

    def test_atomic_max_int32(self):
        self.check_atomic_max(dtype=np.int32, lo=-65535, hi=65535)

    def test_atomic_max_uint32(self):
        self.check_atomic_max(dtype=np.uint32, lo=0, hi=65535)

    def test_atomic_max_int64(self):
        self.check_atomic_max(dtype=np.int64, lo=-65535, hi=65535)

    def test_atomic_max_uint64(self):
        self.check_atomic_max(dtype=np.uint64, lo=0, hi=65535)

    def test_atomic_max_float32(self):
        self.check_atomic_max(dtype=np.float32, lo=-65535, hi=65535)

    def test_atomic_max_double(self):
        self.check_atomic_max(dtype=np.float64, lo=-65535, hi=65535)

    def test_atomic_max_double_normalizedindex(self):
        vals = np.random.randint(0, 65535, size=(32, 32)).astype(np.float64)
        res = np.zeros(1, np.float64)
        cuda_func = cuda.jit('void(float64[:], float64[:,:])')(
            atomic_max_double_normalizedindex)
        cuda_func[32, 32](res, vals)

        gold = np.max(vals)
        np.testing.assert_equal(res, gold)

    def test_atomic_max_double_oneindex(self):
        vals = np.random.randint(0, 128, size=32).astype(np.float64)
        res = np.zeros(1, np.float64)
        cuda_func = cuda.jit('void(float64[:], float64[:])')(
            atomic_max_double_oneindex)
        cuda_func[1, 32](res, vals)

        gold = np.max(vals)
        np.testing.assert_equal(res, gold)

    def check_atomic_min(self, dtype, lo, hi):
        vals = np.random.randint(lo, hi, size=(32, 32)).astype(dtype)
        res = np.array([65535], dtype=vals.dtype)
        cuda_func = cuda.jit(atomic_min)
        cuda_func[32, 32](res, vals)

        gold = np.min(vals)
        np.testing.assert_equal(res, gold)

    def test_atomic_min_int32(self):
        self.check_atomic_min(dtype=np.int32, lo=-65535, hi=65535)

    def test_atomic_min_uint32(self):
        self.check_atomic_min(dtype=np.uint32, lo=0, hi=65535)

    def test_atomic_min_int64(self):
        self.check_atomic_min(dtype=np.int64, lo=-65535, hi=65535)

    def test_atomic_min_uint64(self):
        self.check_atomic_min(dtype=np.uint64, lo=0, hi=65535)

    def test_atomic_min_float(self):
        self.check_atomic_min(dtype=np.float32, lo=-65535, hi=65535)

    def test_atomic_min_double(self):
        self.check_atomic_min(dtype=np.float64, lo=-65535, hi=65535)

    def test_atomic_min_double_normalizedindex(self):
        vals = np.random.randint(0, 65535, size=(32, 32)).astype(np.float64)
        res = np.ones(1, np.float64) * 65535
        cuda_func = cuda.jit('void(float64[:], float64[:,:])')(
            atomic_min_double_normalizedindex)
        cuda_func[32, 32](res, vals)

        gold = np.min(vals)
        np.testing.assert_equal(res, gold)

    def test_atomic_min_double_oneindex(self):
        vals = np.random.randint(0, 128, size=32).astype(np.float64)
        res = np.ones(1, np.float64) * 128
        cuda_func = cuda.jit('void(float64[:], float64[:])')(
            atomic_min_double_oneindex)
        cuda_func[1, 32](res, vals)

        gold = np.min(vals)
        np.testing.assert_equal(res, gold)

    # Taken together, _test_atomic_minmax_nan_location and
    # _test_atomic_minmax_nan_val check that NaNs are treated similarly to the
    # way they are in Python / NumPy - that is, {min,max}(a, b) == a if either
    # a or b is a NaN. For the atomics, this means that the max is taken as the
    # value stored in the memory location rather than the value supplied - i.e.
    # for:
    #
    #    cuda.atomic.{min,max}(ary, idx, val)
    #
    # the result will be ary[idx] for either of ary[idx] or val being NaN.

    def _test_atomic_minmax_nan_location(self, func):

        cuda_func = cuda.jit('void(float64[:], float64[:,:])')(func)

        vals = np.random.randint(0, 128, size=(1,1)).astype(np.float64)
        res = np.zeros(1, np.float64) + np.nan
        cuda_func[1, 1](res, vals)
        np.testing.assert_equal(res, [np.nan])

    def _test_atomic_minmax_nan_val(self, func):
        cuda_func = cuda.jit('void(float64[:], float64[:,:])')(func)

        res = np.random.randint(0, 128, size=1).astype(np.float64)
        gold = res.copy()
        vals = np.zeros((1, 1), np.float64) + np.nan
        cuda_func[1, 1](res, vals)

        np.testing.assert_equal(res, gold)

    def test_atomic_min_nan_location(self):
        self._test_atomic_minmax_nan_location(atomic_min)

    def test_atomic_max_nan_location(self):
        self._test_atomic_minmax_nan_location(atomic_max)

    def test_atomic_min_nan_val(self):
        self._test_atomic_minmax_nan_val(atomic_min)

    def test_atomic_max_nan_val(self):
        self._test_atomic_minmax_nan_val(atomic_max)

    def test_atomic_max_double_shared(self):
        vals = np.random.randint(0, 32, size=32).astype(np.float64)
        res = np.zeros(1, np.float64)
        sig = 'void(float64[:], float64[:])'
        cuda_func = cuda.jit(sig)(atomic_max_double_shared)
        cuda_func[1, 32](res, vals)

        gold = np.max(vals)
        np.testing.assert_equal(res, gold)

    def test_atomic_min_double_shared(self):
        vals = np.random.randint(0, 32, size=32).astype(np.float64)
        res = np.ones(1, np.float64) * 32
        sig = 'void(float64[:], float64[:])'
        cuda_func = cuda.jit(sig)(atomic_min_double_shared)
        cuda_func[1, 32](res, vals)

        gold = np.min(vals)
        np.testing.assert_equal(res, gold)

    def check_cas(self, n, fill, unfill, dtype, cas_func, ndim=1):
        res = [fill] * (n // 2) + [unfill] * (n // 2)
        np.random.shuffle(res)
        res = np.asarray(res, dtype=dtype)
        if ndim == 2:
            res.shape = (10, -1)
        out = np.zeros_like(res)
        ary = np.random.randint(1, 10, size=res.shape).astype(res.dtype)

        fill_mask = res == fill
        unfill_mask = res == unfill

        expect_res = np.zeros_like(res)
        expect_res[fill_mask] = ary[fill_mask]
        expect_res[unfill_mask] = unfill

        expect_out = res.copy()

        cuda_func = cuda.jit(cas_func)
        if ndim == 1:
            cuda_func[10, 10](res, out, ary, fill)
        else:
            cuda_func[(10, 10), (10, 10)](res, out, ary, fill)

        np.testing.assert_array_equal(expect_res, res)
        np.testing.assert_array_equal(expect_out, out)

    def test_atomic_compare_and_swap(self):
        self.check_cas(n=100, fill=-99, unfill=-1, dtype=np.int32,
                       cas_func=atomic_compare_and_swap)

    def test_atomic_compare_and_swap2(self):
        self.check_cas(n=100, fill=-45, unfill=-1, dtype=np.int64,
                       cas_func=atomic_compare_and_swap)

    def test_atomic_compare_and_swap3(self):
        rfill = np.random.randint(50, 500, dtype=np.uint32)
        runfill = np.random.randint(1, 25, dtype=np.uint32)
        self.check_cas(n=100, fill=rfill, unfill=runfill, dtype=np.uint32,
                       cas_func=atomic_compare_and_swap)

    def test_atomic_compare_and_swap4(self):
        rfill = np.random.randint(50, 500, dtype=np.uint64)
        runfill = np.random.randint(1, 25, dtype=np.uint64)
        self.check_cas(n=100, fill=rfill, unfill=runfill, dtype=np.uint64,
                       cas_func=atomic_compare_and_swap)

    def test_atomic_cas_1dim(self):
        self.check_cas(n=100, fill=-99, unfill=-1, dtype=np.int32,
                       cas_func=atomic_cas_1dim)

    def test_atomic_cas_2dim(self):
        self.check_cas(n=100, fill=-99, unfill=-1, dtype=np.int32,
                       cas_func=atomic_cas_2dim, ndim=2)

    def test_atomic_cas2_1dim(self):
        self.check_cas(n=100, fill=-45, unfill=-1, dtype=np.int64,
                       cas_func=atomic_cas_1dim)

    def test_atomic_cas2_2dim(self):
        self.check_cas(n=100, fill=-45, unfill=-1, dtype=np.int64,
                       cas_func=atomic_cas_2dim, ndim=2)

    def test_atomic_cas3_1dim(self):
        rfill = np.random.randint(50, 500, dtype=np.uint32)
        runfill = np.random.randint(1, 25, dtype=np.uint32)
        self.check_cas(n=100, fill=rfill, unfill=runfill, dtype=np.uint32,
                       cas_func=atomic_cas_1dim)

    def test_atomic_cas3_2dim(self):
        rfill = np.random.randint(50, 500, dtype=np.uint32)
        runfill = np.random.randint(1, 25, dtype=np.uint32)
        self.check_cas(n=100, fill=rfill, unfill=runfill, dtype=np.uint32,
                       cas_func=atomic_cas_2dim, ndim=2)

    def test_atomic_cas4_1dim(self):
        rfill = np.random.randint(50, 500, dtype=np.uint64)
        runfill = np.random.randint(1, 25, dtype=np.uint64)
        self.check_cas(n=100, fill=rfill, unfill=runfill, dtype=np.uint64,
                       cas_func=atomic_cas_1dim)

    def test_atomic_cas4_2dim(self):
        rfill = np.random.randint(50, 500, dtype=np.uint64)
        runfill = np.random.randint(1, 25, dtype=np.uint64)
        self.check_cas(n=100, fill=rfill, unfill=runfill, dtype=np.uint64,
                       cas_func=atomic_cas_2dim, ndim=2)

    # Tests that the atomic add, min, and max operations return the old value -
    # in the simulator, they did not (see Issue #5458). The max and min have
    # special handling for NaN values, so we explicitly test with a NaN in the
    # array being modified and the value provided.

    def _test_atomic_returns_old(self, kernel, initial):
        x = np.zeros(2, dtype=np.float32)
        x[0] = initial
        kernel[1, 1](x)
        if np.isnan(initial):
            self.assertTrue(np.isnan(x[1]))
        else:
            self.assertEqual(x[1], initial)

    def test_atomic_add_returns_old(self):
        @cuda.jit
        def kernel(x):
            x[1] = cuda.atomic.add(x, 0, 1)

        self._test_atomic_returns_old(kernel, 10)

    def test_atomic_max_returns_no_replace(self):
        @cuda.jit
        def kernel(x):
            x[1] = cuda.atomic.max(x, 0, 1)

        self._test_atomic_returns_old(kernel, 10)

    def test_atomic_max_returns_old_replace(self):
        @cuda.jit
        def kernel(x):
            x[1] = cuda.atomic.max(x, 0, 10)

        self._test_atomic_returns_old(kernel, 1)

    def test_atomic_max_returns_old_nan_in_array(self):
        @cuda.jit
        def kernel(x):
            x[1] = cuda.atomic.max(x, 0, 1)

        self._test_atomic_returns_old(kernel, np.nan)

    def test_atomic_max_returns_old_nan_val(self):
        @cuda.jit
        def kernel(x):
            x[1] = cuda.atomic.max(x, 0, np.nan)

        self._test_atomic_returns_old(kernel, 10)

    def test_atomic_min_returns_old_no_replace(self):
        @cuda.jit
        def kernel(x):
            x[1] = cuda.atomic.min(x, 0, 11)

        self._test_atomic_returns_old(kernel, 10)

    def test_atomic_min_returns_old_replace(self):
        @cuda.jit
        def kernel(x):
            x[1] = cuda.atomic.min(x, 0, 10)

        self._test_atomic_returns_old(kernel, 11)

    def test_atomic_min_returns_old_nan_in_array(self):
        @cuda.jit
        def kernel(x):
            x[1] = cuda.atomic.min(x, 0, 11)

        self._test_atomic_returns_old(kernel, np.nan)

    def test_atomic_min_returns_old_nan_val(self):
        @cuda.jit
        def kernel(x):
            x[1] = cuda.atomic.min(x, 0, np.nan)

        self._test_atomic_returns_old(kernel, 11)

    # Tests for atomic nanmin/nanmax

    # nanmax tests
    def check_atomic_nanmax(self, dtype, lo, hi, init_val):
        vals = np.random.randint(lo, hi, size=(32, 32)).astype(dtype)
        vals[1::2] = init_val
        res = np.zeros(1, dtype=vals.dtype)
        cuda_func = cuda.jit(atomic_nanmax)
        cuda_func[32, 32](res, vals)
        gold = np.nanmax(vals)
        np.testing.assert_equal(res, gold)

    def test_atomic_nanmax_int32(self):
        self.check_atomic_nanmax(dtype=np.int32, lo=-65535, hi=65535,
                                 init_val=0)

    def test_atomic_nanmax_uint32(self):
        self.check_atomic_nanmax(dtype=np.uint32, lo=0, hi=65535,
                                 init_val=0)

    def test_atomic_nanmax_int64(self):
        self.check_atomic_nanmax(dtype=np.int64, lo=-65535, hi=65535,
                                 init_val=0)

    def test_atomic_nanmax_uint64(self):
        self.check_atomic_nanmax(dtype=np.uint64, lo=0, hi=65535,
                                 init_val=0)

    def test_atomic_nanmax_float32(self):
        self.check_atomic_nanmax(dtype=np.float32, lo=-65535, hi=65535,
                                 init_val=np.nan)

    def test_atomic_nanmax_double(self):
        self.check_atomic_nanmax(dtype=np.float64, lo=-65535, hi=65535,
                                 init_val=np.nan)

    def test_atomic_nanmax_double_shared(self):
        vals = np.random.randint(0, 32, size=32).astype(np.float64)
        vals[1::2] = np.nan
        res = np.array([0], dtype=vals.dtype)
        sig = 'void(float64[:], float64[:])'
        cuda_func = cuda.jit(sig)(atomic_nanmax_double_shared)
        cuda_func[1, 32](res, vals)

        gold = np.nanmax(vals)
        np.testing.assert_equal(res, gold)

    def test_atomic_nanmax_double_oneindex(self):
        vals = np.random.randint(0, 128, size=32).astype(np.float64)
        vals[1::2] = np.nan
        res = np.zeros(1, np.float64)
        cuda_func = cuda.jit('void(float64[:], float64[:])')(
            atomic_max_double_oneindex)
        cuda_func[1, 32](res, vals)

        gold = np.nanmax(vals)
        np.testing.assert_equal(res, gold)

    # nanmin tests
    def check_atomic_nanmin(self, dtype, lo, hi, init_val):
        vals = np.random.randint(lo, hi, size=(32, 32)).astype(dtype)
        vals[1::2] = init_val
        res = np.array([65535], dtype=vals.dtype)
        cuda_func = cuda.jit(atomic_nanmin)
        cuda_func[32, 32](res, vals)

        gold = np.nanmin(vals)
        np.testing.assert_equal(res, gold)

    def test_atomic_nanmin_int32(self):
        self.check_atomic_nanmin(dtype=np.int32, lo=-65535, hi=65535,
                                 init_val=0)

    def test_atomic_nanmin_uint32(self):
        self.check_atomic_nanmin(dtype=np.uint32, lo=0, hi=65535,
                                 init_val=0)

    def test_atomic_nanmin_int64(self):
        self.check_atomic_nanmin(dtype=np.int64, lo=-65535, hi=65535,
                                 init_val=0)

    def test_atomic_nanmin_uint64(self):
        self.check_atomic_nanmin(dtype=np.uint64, lo=0, hi=65535,
                                 init_val=0)

    def test_atomic_nanmin_float(self):
        self.check_atomic_nanmin(dtype=np.float32, lo=-65535, hi=65535,
                                 init_val=np.nan)

    def test_atomic_nanmin_double(self):
        self.check_atomic_nanmin(dtype=np.float64, lo=-65535, hi=65535,
                                 init_val=np.nan)

    def test_atomic_nanmin_double_shared(self):
        vals = np.random.randint(0, 32, size=32).astype(np.float64)
        vals[1::2] = np.nan
        res = np.array([32], dtype=vals.dtype)
        sig = 'void(float64[:], float64[:])'
        cuda_func = cuda.jit(sig)(atomic_nanmin_double_shared)
        cuda_func[1, 32](res, vals)

        gold = np.nanmin(vals)
        np.testing.assert_equal(res, gold)

    def test_atomic_nanmin_double_oneindex(self):
        vals = np.random.randint(0, 128, size=32).astype(np.float64)
        vals[1::2] = np.nan
        res = np.array([128], np.float64)
        cuda_func = cuda.jit('void(float64[:], float64[:])')(
            atomic_min_double_oneindex)
        cuda_func[1, 32](res, vals)

        gold = np.nanmin(vals)
        np.testing.assert_equal(res, gold)

    # Returning old value tests

    def _test_atomic_nan_returns_old(self, kernel, initial):
        x = np.zeros(2, dtype=np.float32)
        x[0] = initial
        x[1] = np.nan
        kernel[1, 1](x)
        if np.isnan(initial):
            self.assertFalse(np.isnan(x[0]))
            self.assertTrue(np.isnan(x[1]))
        else:
            self.assertEqual(x[1], initial)

    def test_atomic_nanmax_returns_old_no_replace(self):
        @cuda.jit
        def kernel(x):
            x[1] = cuda.atomic.nanmax(x, 0, 1)

        self._test_atomic_nan_returns_old(kernel, 10)

    def test_atomic_nanmax_returns_old_replace(self):
        @cuda.jit
        def kernel(x):
            x[1] = cuda.atomic.nanmax(x, 0, 10)

        self._test_atomic_nan_returns_old(kernel, 1)

    def test_atomic_nanmax_returns_old_nan_in_array(self):
        @cuda.jit
        def kernel(x):
            x[1] = cuda.atomic.nanmax(x, 0, 1)

        self._test_atomic_nan_returns_old(kernel, np.nan)

    def test_atomic_nanmax_returns_old_nan_val(self):
        @cuda.jit
        def kernel(x):
            x[1] = cuda.atomic.nanmax(x, 0, np.nan)

        self._test_atomic_nan_returns_old(kernel, 10)

    def test_atomic_nanmin_returns_old_no_replace(self):
        @cuda.jit
        def kernel(x):
            x[1] = cuda.atomic.nanmin(x, 0, 11)

        self._test_atomic_nan_returns_old(kernel, 10)

    def test_atomic_nanmin_returns_old_replace(self):
        @cuda.jit
        def kernel(x):
            x[1] = cuda.atomic.nanmin(x, 0, 10)

        self._test_atomic_nan_returns_old(kernel, 11)

    def test_atomic_nanmin_returns_old_nan_in_array(self):
        @cuda.jit
        def kernel(x):
            x[1] = cuda.atomic.nanmin(x, 0, 11)

        self._test_atomic_nan_returns_old(kernel, np.nan)

    def test_atomic_nanmin_returns_old_nan_val(self):
        @cuda.jit
        def kernel(x):
            x[1] = cuda.atomic.nanmin(x, 0, np.nan)

        self._test_atomic_nan_returns_old(kernel, 11)


if __name__ == '__main__':
    unittest.main()
