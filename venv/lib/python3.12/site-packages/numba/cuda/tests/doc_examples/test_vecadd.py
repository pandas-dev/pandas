import unittest

from numba.cuda.testing import CUDATestCase, skip_on_cudasim
from numba.tests.support import captured_stdout


@skip_on_cudasim("cudasim doesn't support cuda import at non-top-level")
class TestVecAdd(CUDATestCase):
    """
    Test simple vector addition
    """

    def setUp(self):
        # Prevent output from this test showing
        # up when running the test suite
        self._captured_stdout = captured_stdout()
        self._captured_stdout.__enter__()
        super().setUp()

    def tearDown(self):
        # No exception type, value, or traceback
        self._captured_stdout.__exit__(None, None, None)
        super().tearDown()

    def test_ex_vecadd(self):
        # ex_vecadd.import.begin
        import numpy as np
        from numba import cuda
        # ex_vecadd.import.end

        # ex_vecadd.kernel.begin
        @cuda.jit
        def f(a, b, c):
            # like threadIdx.x + (blockIdx.x * blockDim.x)
            tid = cuda.grid(1)
            size = len(c)

            if tid < size:
                c[tid] = a[tid] + b[tid]
        # ex_vecadd.kernel.end

        # Seed RNG for test repeatability
        np.random.seed(1)

        # ex_vecadd.allocate.begin
        N = 100000
        a = cuda.to_device(np.random.random(N))
        b = cuda.to_device(np.random.random(N))
        c = cuda.device_array_like(a)
        # ex_vecadd.allocate.end

        # ex_vecadd.forall.begin
        f.forall(len(a))(a, b, c)
        print(c.copy_to_host())
        # ex_vecadd.forall.end

        # ex_vecadd.launch.begin
        # Enough threads per block for several warps per block
        nthreads = 256
        # Enough blocks to cover the entire vector depending on its length
        nblocks = (len(a) // nthreads) + 1
        f[nblocks, nthreads](a, b, c)
        print(c.copy_to_host())
        # ex_vecadd.launch.end

        np.testing.assert_equal(
            c.copy_to_host(),
            a.copy_to_host() + b.copy_to_host()
        )


if __name__ == "__main__":
    unittest.main()
