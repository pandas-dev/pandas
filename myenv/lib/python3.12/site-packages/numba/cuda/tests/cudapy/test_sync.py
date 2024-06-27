import numpy as np
from numba import cuda, int32, float32
from numba.cuda.testing import skip_on_cudasim, unittest, CUDATestCase
from numba.core.config import ENABLE_CUDASIM


def useless_syncthreads(ary):
    i = cuda.grid(1)
    cuda.syncthreads()
    ary[i] = i


def useless_syncwarp(ary):
    i = cuda.grid(1)
    cuda.syncwarp()
    ary[i] = i


def useless_syncwarp_with_mask(ary):
    i = cuda.grid(1)
    cuda.syncwarp(0xFFFF)
    ary[i] = i


def coop_syncwarp(res):
    sm = cuda.shared.array(32, int32)
    i = cuda.grid(1)

    sm[i] = i
    cuda.syncwarp()

    if i < 16:
        sm[i] = sm[i] + sm[i + 16]
        cuda.syncwarp(0xFFFF)

    if i < 8:
        sm[i] = sm[i] + sm[i + 8]
        cuda.syncwarp(0xFF)

    if i < 4:
        sm[i] = sm[i] + sm[i + 4]
        cuda.syncwarp(0xF)

    if i < 2:
        sm[i] = sm[i] + sm[i + 2]
        cuda.syncwarp(0x3)

    if i == 0:
        res[0] = sm[0] + sm[1]


def simple_smem(ary):
    N = 100
    sm = cuda.shared.array(N, int32)
    i = cuda.grid(1)
    if i == 0:
        for j in range(N):
            sm[j] = j
    cuda.syncthreads()
    ary[i] = sm[i]


def coop_smem2d(ary):
    i, j = cuda.grid(2)
    sm = cuda.shared.array((10, 20), float32)
    sm[i, j] = (i + 1) / (j + 1)
    cuda.syncthreads()
    ary[i, j] = sm[i, j]


def dyn_shared_memory(ary):
    i = cuda.grid(1)
    sm = cuda.shared.array(0, float32)
    sm[i] = i * 2
    cuda.syncthreads()
    ary[i] = sm[i]


def use_threadfence(ary):
    ary[0] += 123
    cuda.threadfence()
    ary[0] += 321


def use_threadfence_block(ary):
    ary[0] += 123
    cuda.threadfence_block()
    ary[0] += 321


def use_threadfence_system(ary):
    ary[0] += 123
    cuda.threadfence_system()
    ary[0] += 321


def use_syncthreads_count(ary_in, ary_out):
    i = cuda.grid(1)
    ary_out[i] = cuda.syncthreads_count(ary_in[i])


def use_syncthreads_and(ary_in, ary_out):
    i = cuda.grid(1)
    ary_out[i] = cuda.syncthreads_and(ary_in[i])


def use_syncthreads_or(ary_in, ary_out):
    i = cuda.grid(1)
    ary_out[i] = cuda.syncthreads_or(ary_in[i])


def _safe_cc_check(cc):
    if ENABLE_CUDASIM:
        return True
    else:
        return cuda.get_current_device().compute_capability >= cc


class TestCudaSync(CUDATestCase):
    def _test_useless(self, kernel):
        compiled = cuda.jit("void(int32[::1])")(kernel)
        nelem = 10
        ary = np.empty(nelem, dtype=np.int32)
        exp = np.arange(nelem, dtype=np.int32)
        compiled[1, nelem](ary)
        np.testing.assert_equal(ary, exp)

    def test_useless_syncthreads(self):
        self._test_useless(useless_syncthreads)

    @skip_on_cudasim("syncwarp not implemented on cudasim")
    def test_useless_syncwarp(self):
        self._test_useless(useless_syncwarp)

    @skip_on_cudasim("syncwarp not implemented on cudasim")
    @unittest.skipUnless(_safe_cc_check((7, 0)),
                         "Partial masks require CC 7.0 or greater")
    def test_useless_syncwarp_with_mask(self):
        self._test_useless(useless_syncwarp_with_mask)

    @skip_on_cudasim("syncwarp not implemented on cudasim")
    @unittest.skipUnless(_safe_cc_check((7, 0)),
                         "Partial masks require CC 7.0 or greater")
    def test_coop_syncwarp(self):
        # coop_syncwarp computes the sum of all integers from 0 to 31 (496)
        # using a single warp
        expected = 496
        nthreads = 32
        nblocks = 1

        compiled = cuda.jit("void(int32[::1])")(coop_syncwarp)
        res = np.zeros(1, dtype=np.int32)
        compiled[nblocks, nthreads](res)
        np.testing.assert_equal(expected, res[0])

    def test_simple_smem(self):
        compiled = cuda.jit("void(int32[::1])")(simple_smem)
        nelem = 100
        ary = np.empty(nelem, dtype=np.int32)
        compiled[1, nelem](ary)
        self.assertTrue(np.all(ary == np.arange(nelem, dtype=np.int32)))

    def test_coop_smem2d(self):
        compiled = cuda.jit("void(float32[:,::1])")(coop_smem2d)
        shape = 10, 20
        ary = np.empty(shape, dtype=np.float32)
        compiled[1, shape](ary)
        exp = np.empty_like(ary)
        for i in range(ary.shape[0]):
            for j in range(ary.shape[1]):
                exp[i, j] = (i + 1) / (j + 1)
        self.assertTrue(np.allclose(ary, exp))

    def test_dyn_shared_memory(self):
        compiled = cuda.jit("void(float32[::1])")(dyn_shared_memory)
        shape = 50
        ary = np.empty(shape, dtype=np.float32)
        compiled[1, shape, 0, ary.size * 4](ary)
        self.assertTrue(np.all(ary == 2 * np.arange(ary.size, dtype=np.int32)))

    def test_threadfence_codegen(self):
        # Does not test runtime behavior, just the code generation.
        sig = (int32[:],)
        compiled = cuda.jit(sig)(use_threadfence)
        ary = np.zeros(10, dtype=np.int32)
        compiled[1, 1](ary)
        self.assertEqual(123 + 321, ary[0])
        if not ENABLE_CUDASIM:
            self.assertIn("membar.gl;", compiled.inspect_asm(sig))

    def test_threadfence_block_codegen(self):
        # Does not test runtime behavior, just the code generation.
        sig = (int32[:],)
        compiled = cuda.jit(sig)(use_threadfence_block)
        ary = np.zeros(10, dtype=np.int32)
        compiled[1, 1](ary)
        self.assertEqual(123 + 321, ary[0])
        if not ENABLE_CUDASIM:
            self.assertIn("membar.cta;", compiled.inspect_asm(sig))

    def test_threadfence_system_codegen(self):
        # Does not test runtime behavior, just the code generation.
        sig = (int32[:],)
        compiled = cuda.jit(sig)(use_threadfence_system)
        ary = np.zeros(10, dtype=np.int32)
        compiled[1, 1](ary)
        self.assertEqual(123 + 321, ary[0])
        if not ENABLE_CUDASIM:
            self.assertIn("membar.sys;", compiled.inspect_asm(sig))

    def _test_syncthreads_count(self, in_dtype):
        compiled = cuda.jit(use_syncthreads_count)
        ary_in = np.ones(72, dtype=in_dtype)
        ary_out = np.zeros(72, dtype=np.int32)
        ary_in[31] = 0
        ary_in[42] = 0
        compiled[1, 72](ary_in, ary_out)
        self.assertTrue(np.all(ary_out == 70))

    def test_syncthreads_count(self):
        self._test_syncthreads_count(np.int32)

    def test_syncthreads_count_upcast(self):
        self._test_syncthreads_count(np.int16)

    def test_syncthreads_count_downcast(self):
        self._test_syncthreads_count(np.int64)

    def _test_syncthreads_and(self, in_dtype):
        compiled = cuda.jit(use_syncthreads_and)
        nelem = 100
        ary_in = np.ones(nelem, dtype=in_dtype)
        ary_out = np.zeros(nelem, dtype=np.int32)
        compiled[1, nelem](ary_in, ary_out)
        self.assertTrue(np.all(ary_out == 1))
        ary_in[31] = 0
        compiled[1, nelem](ary_in, ary_out)
        self.assertTrue(np.all(ary_out == 0))

    def test_syncthreads_and(self):
        self._test_syncthreads_and(np.int32)

    def test_syncthreads_and_upcast(self):
        self._test_syncthreads_and(np.int16)

    def test_syncthreads_and_downcast(self):
        self._test_syncthreads_and(np.int64)

    def _test_syncthreads_or(self, in_dtype):
        compiled = cuda.jit(use_syncthreads_or)
        nelem = 100
        ary_in = np.zeros(nelem, dtype=in_dtype)
        ary_out = np.zeros(nelem, dtype=np.int32)
        compiled[1, nelem](ary_in, ary_out)
        self.assertTrue(np.all(ary_out == 0))
        ary_in[31] = 1
        compiled[1, nelem](ary_in, ary_out)
        self.assertTrue(np.all(ary_out == 1))

    def test_syncthreads_or(self):
        self._test_syncthreads_or(np.int32)

    def test_syncthreads_or_upcast(self):
        self._test_syncthreads_or(np.int16)

    def test_syncthreads_or_downcast(self):
        self._test_syncthreads_or(np.int64)


if __name__ == '__main__':
    unittest.main()
