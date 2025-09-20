import unittest

from numba.cuda.testing import CUDATestCase, skip_on_cudasim
from numba.tests.support import captured_stdout
import numpy as np


@skip_on_cudasim("cudasim doesn't support cuda import at non-top-level")
class TestCpuGpuCompat(CUDATestCase):
    """
    Test compatibility of CPU and GPU functions
    """

    def setUp(self):
        # Prevent output from this test showing up when running the test suite
        self._captured_stdout = captured_stdout()
        self._captured_stdout.__enter__()
        super().setUp()

    def tearDown(self):
        # No exception type, value, or traceback
        self._captured_stdout.__exit__(None, None, None)
        super().tearDown()

    def test_ex_cpu_gpu_compat(self):
        # ex_cpu_gpu_compat.import.begin
        from math import pi

        import numba
        from numba import cuda
        # ex_cpu_gpu_compat.import.end

        # ex_cpu_gpu_compat.allocate.begin
        X = cuda.to_device([1, 10, 234])
        Y = cuda.to_device([2, 2, 4014])
        Z = cuda.to_device([3, 14, 2211])
        results = cuda.to_device([0.0, 0.0, 0.0])
        # ex_cpu_gpu_compat.allocate.end

        # ex_cpu_gpu_compat.define.begin
        @numba.jit
        def business_logic(x, y, z):
            return 4 * z * (2 * x - (4 * y) / 2 * pi)
        # ex_cpu_gpu_compat.define.end

        # ex_cpu_gpu_compat.cpurun.begin
        print(business_logic(1, 2, 3))  # -126.79644737231007
        # ex_cpu_gpu_compat.cpurun.end

        # ex_cpu_gpu_compat.usegpu.begin
        @cuda.jit
        def f(res, xarr, yarr, zarr):
            tid = cuda.grid(1)
            if tid < len(xarr):
                # The function decorated with numba.jit may be directly reused
                res[tid] = business_logic(xarr[tid], yarr[tid], zarr[tid])
        # ex_cpu_gpu_compat.usegpu.end

        # ex_cpu_gpu_compat.launch.begin
        f.forall(len(X))(results, X, Y, Z)
        print(results)
        # [-126.79644737231007, 416.28324559588634, -218912930.2987788]
        # ex_cpu_gpu_compat.launch.end

        expect = [
            business_logic(x, y, z) for x, y, z in zip(X, Y, Z)
        ]

        np.testing.assert_equal(
            expect,
            results.copy_to_host()
        )


if __name__ == "__main__":
    unittest.main()
