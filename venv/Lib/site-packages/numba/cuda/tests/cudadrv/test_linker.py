import numpy as np
import warnings
from numba.cuda.testing import unittest
from numba.cuda.testing import (skip_on_cudasim, skip_if_cuda_includes_missing)
from numba.cuda.testing import CUDATestCase, test_data_dir
from numba.cuda.cudadrv.driver import (CudaAPIError, Linker,
                                       LinkerError)
from numba.cuda.cudadrv.error import NvrtcError
from numba.cuda import require_context
from numba.tests.support import ignore_internal_warnings
from numba import cuda, void, float64, int64, int32, typeof, float32


CONST1D = np.arange(10, dtype=np.float64)


def simple_const_mem(A):
    C = cuda.const.array_like(CONST1D)
    i = cuda.grid(1)

    A[i] = C[i] + 1.0


def func_with_lots_of_registers(x, a, b, c, d, e, f):
    a1 = 1.0
    a2 = 1.0
    a3 = 1.0
    a4 = 1.0
    a5 = 1.0
    b1 = 1.0
    b2 = 1.0
    b3 = 1.0
    b4 = 1.0
    b5 = 1.0
    c1 = 1.0
    c2 = 1.0
    c3 = 1.0
    c4 = 1.0
    c5 = 1.0
    d1 = 10
    d2 = 10
    d3 = 10
    d4 = 10
    d5 = 10
    for i in range(a):
        a1 += b
        a2 += c
        a3 += d
        a4 += e
        a5 += f
        b1 *= b
        b2 *= c
        b3 *= d
        b4 *= e
        b5 *= f
        c1 /= b
        c2 /= c
        c3 /= d
        c4 /= e
        c5 /= f
        d1 <<= b
        d2 <<= c
        d3 <<= d
        d4 <<= e
        d5 <<= f
    x[cuda.grid(1)] = a1 + a2 + a3 + a4 + a5
    x[cuda.grid(1)] += b1 + b2 + b3 + b4 + b5
    x[cuda.grid(1)] += c1 + c2 + c3 + c4 + c5
    x[cuda.grid(1)] += d1 + d2 + d3 + d4 + d5


def simple_smem(ary, dty):
    sm = cuda.shared.array(100, dty)
    i = cuda.grid(1)
    if i == 0:
        for j in range(100):
            sm[j] = j
    cuda.syncthreads()
    ary[i] = sm[i]


def coop_smem2d(ary):
    i, j = cuda.grid(2)
    sm = cuda.shared.array((10, 20), float32)
    sm[i, j] = (i + 1) / (j + 1)
    cuda.syncthreads()
    ary[i, j] = sm[i, j]


def simple_maxthreads(ary):
    i = cuda.grid(1)
    ary[i] = i


LMEM_SIZE = 1000


def simple_lmem(A, B, dty):
    C = cuda.local.array(LMEM_SIZE, dty)
    for i in range(C.shape[0]):
        C[i] = A[i]
    for i in range(C.shape[0]):
        B[i] = C[i]


@skip_on_cudasim('Linking unsupported in the simulator')
class TestLinker(CUDATestCase):
    _NUMBA_NVIDIA_BINDING_0_ENV = {'NUMBA_CUDA_USE_NVIDIA_BINDING': '0'}

    @require_context
    def test_linker_basic(self):
        '''Simply go through the constructor and destructor
        '''
        linker = Linker.new(cc=(5, 3))
        del linker

    def _test_linking(self, eager):
        global bar  # must be a global; other it is recognized as a freevar
        bar = cuda.declare_device('bar', 'int32(int32)')

        link = str(test_data_dir / 'jitlink.ptx')

        if eager:
            args = ['void(int32[:], int32[:])']
        else:
            args = []

        @cuda.jit(*args, link=[link])
        def foo(x, y):
            i = cuda.grid(1)
            x[i] += bar(y[i])

        A = np.array([123], dtype=np.int32)
        B = np.array([321], dtype=np.int32)

        foo[1, 1](A, B)

        self.assertTrue(A[0] == 123 + 2 * 321)

    def test_linking_lazy_compile(self):
        self._test_linking(eager=False)

    def test_linking_eager_compile(self):
        self._test_linking(eager=True)

    def test_linking_cu(self):
        bar = cuda.declare_device('bar', 'int32(int32)')

        link = str(test_data_dir / 'jitlink.cu')

        @cuda.jit(link=[link])
        def kernel(r, x):
            i = cuda.grid(1)

            if i < len(r):
                r[i] = bar(x[i])

        x = np.arange(10, dtype=np.int32)
        r = np.zeros_like(x)

        kernel[1, 32](r, x)

        # Matches the operation of bar() in jitlink.cu
        expected = x * 2
        np.testing.assert_array_equal(r, expected)

    def test_linking_cu_log_warning(self):
        bar = cuda.declare_device('bar', 'int32(int32)')

        link = str(test_data_dir / 'warn.cu')

        with warnings.catch_warnings(record=True) as w:
            ignore_internal_warnings()

            @cuda.jit('void(int32)', link=[link])
            def kernel(x):
                bar(x)

        self.assertEqual(len(w), 1, 'Expected warnings from NVRTC')
        # Check the warning refers to the log messages
        self.assertIn('NVRTC log messages', str(w[0].message))
        # Check the message pertaining to the unused variable is provided
        self.assertIn('declared but never referenced', str(w[0].message))

    def test_linking_cu_error(self):
        bar = cuda.declare_device('bar', 'int32(int32)')

        link = str(test_data_dir / 'error.cu')

        with self.assertRaises(NvrtcError) as e:
            @cuda.jit('void(int32)', link=[link])
            def kernel(x):
                bar(x)

        msg = e.exception.args[0]
        # Check the error message refers to the NVRTC compile
        self.assertIn('NVRTC Compilation failure', msg)
        # Check the expected error in the CUDA source is reported
        self.assertIn('identifier "SYNTAX" is undefined', msg)
        # Check the filename is reported correctly
        self.assertIn('in the compilation of "error.cu"', msg)

    def test_linking_unknown_filetype_error(self):
        expected_err = "Don't know how to link file with extension .cuh"
        with self.assertRaisesRegex(RuntimeError, expected_err):
            @cuda.jit('void()', link=['header.cuh'])
            def kernel():
                pass

    def test_linking_file_with_no_extension_error(self):
        expected_err = "Don't know how to link file with no extension"
        with self.assertRaisesRegex(RuntimeError, expected_err):
            @cuda.jit('void()', link=['data'])
            def kernel():
                pass

    @skip_if_cuda_includes_missing
    def test_linking_cu_cuda_include(self):
        link = str(test_data_dir / 'cuda_include.cu')

        # An exception will be raised when linking this kernel due to the
        # compile failure if CUDA includes cannot be found by Nvrtc.
        @cuda.jit('void()', link=[link])
        def kernel():
            pass

    def test_try_to_link_nonexistent(self):
        with self.assertRaises(LinkerError) as e:
            @cuda.jit('void(int32[::1])', link=['nonexistent.a'])
            def f(x):
                x[0] = 0
        self.assertIn('nonexistent.a not found', e.exception.args)

    def test_set_registers_no_max(self):
        """Ensure that the jitted kernel used in the test_set_registers_* tests
        uses more than 57 registers - this ensures that test_set_registers_*
        are really checking that they reduced the number of registers used from
        something greater than the maximum."""
        compiled = cuda.jit(func_with_lots_of_registers)
        compiled = compiled.specialize(np.empty(32), *range(6))
        self.assertGreater(compiled.get_regs_per_thread(), 57)

    def test_set_registers_57(self):
        compiled = cuda.jit(max_registers=57)(func_with_lots_of_registers)
        compiled = compiled.specialize(np.empty(32), *range(6))
        self.assertLessEqual(compiled.get_regs_per_thread(), 57)

    def test_set_registers_38(self):
        compiled = cuda.jit(max_registers=38)(func_with_lots_of_registers)
        compiled = compiled.specialize(np.empty(32), *range(6))
        self.assertLessEqual(compiled.get_regs_per_thread(), 38)

    def test_set_registers_eager(self):
        sig = void(float64[::1], int64, int64, int64, int64, int64, int64)
        compiled = cuda.jit(sig, max_registers=38)(func_with_lots_of_registers)
        self.assertLessEqual(compiled.get_regs_per_thread(), 38)

    def test_get_const_mem_size(self):
        sig = void(float64[::1])
        compiled = cuda.jit(sig)(simple_const_mem)
        const_mem_size = compiled.get_const_mem_size()
        self.assertGreaterEqual(const_mem_size, CONST1D.nbytes)

    def test_get_no_shared_memory(self):
        compiled = cuda.jit(func_with_lots_of_registers)
        compiled = compiled.specialize(np.empty(32), *range(6))
        shared_mem_size = compiled.get_shared_mem_per_block()
        self.assertEqual(shared_mem_size, 0)

    def test_get_shared_mem_per_block(self):
        sig = void(int32[::1], typeof(np.int32))
        compiled = cuda.jit(sig)(simple_smem)
        shared_mem_size = compiled.get_shared_mem_per_block()
        self.assertEqual(shared_mem_size, 400)

    def test_get_shared_mem_per_specialized(self):
        compiled = cuda.jit(simple_smem)
        compiled_specialized = compiled.specialize(
            np.zeros(100, dtype=np.int32), np.float64)
        shared_mem_size = compiled_specialized.get_shared_mem_per_block()
        self.assertEqual(shared_mem_size, 800)

    def test_get_max_threads_per_block(self):
        compiled = cuda.jit("void(float32[:,::1])")(coop_smem2d)
        max_threads = compiled.get_max_threads_per_block()
        self.assertGreater(max_threads, 0)

    def test_max_threads_exceeded(self):
        compiled = cuda.jit("void(int32[::1])")(simple_maxthreads)
        max_threads = compiled.get_max_threads_per_block()
        nelem = max_threads + 1
        ary = np.empty(nelem, dtype=np.int32)
        try:
            compiled[1, nelem](ary)
        except CudaAPIError as e:
            self.assertIn("cuLaunchKernel", e.msg)

    def test_get_local_mem_per_thread(self):
        sig = void(int32[::1], int32[::1], typeof(np.int32))
        compiled = cuda.jit(sig)(simple_lmem)
        local_mem_size = compiled.get_local_mem_per_thread()
        calc_size = np.dtype(np.int32).itemsize * LMEM_SIZE
        self.assertGreaterEqual(local_mem_size, calc_size)

    def test_get_local_mem_per_specialized(self):
        compiled = cuda.jit(simple_lmem)
        compiled_specialized = compiled.specialize(
            np.zeros(LMEM_SIZE, dtype=np.int32),
            np.zeros(LMEM_SIZE, dtype=np.int32),
            np.float64)
        local_mem_size = compiled_specialized.get_local_mem_per_thread()
        calc_size = np.dtype(np.float64).itemsize * LMEM_SIZE
        self.assertGreaterEqual(local_mem_size, calc_size)


if __name__ == '__main__':
    unittest.main()
