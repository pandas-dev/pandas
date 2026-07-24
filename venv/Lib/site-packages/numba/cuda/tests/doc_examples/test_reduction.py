import unittest

from numba.cuda.testing import CUDATestCase, skip_on_cudasim
from numba.tests.support import captured_stdout


@skip_on_cudasim("cudasim doesn't support cuda import at non-top-level")
class TestReduction(CUDATestCase):
    """
    Test shared memory reduction
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

    def test_ex_reduction(self):
        # ex_reduction.import.begin
        import numpy as np
        from numba import cuda
        from numba.types import int32
        # ex_reduction.import.end

        # ex_reduction.allocate.begin
        # generate data
        a = cuda.to_device(np.arange(1024))
        nelem = len(a)
        # ex_reduction.allocate.end

        # ex_reduction.kernel.begin
        @cuda.jit
        def array_sum(data):
            tid = cuda.threadIdx.x
            size = len(data)
            if tid < size:
                i = cuda.grid(1)

                # Declare an array in shared memory
                shr = cuda.shared.array(nelem, int32)
                shr[tid] = data[i]

                # Ensure writes to shared memory are visible
                # to all threads before reducing
                cuda.syncthreads()

                s = 1
                while s < cuda.blockDim.x:
                    if tid % (2 * s) == 0:
                        # Stride by `s` and add
                        shr[tid] += shr[tid + s]
                    s *= 2
                    cuda.syncthreads()

                # After the loop, the zeroth  element contains the sum
                if tid == 0:
                    data[tid] = shr[tid]
        # ex_reduction.kernel.end

        # ex_reduction.launch.begin
        array_sum[1, nelem](a)
        print(a[0])                  # 523776
        print(sum(np.arange(1024)))  # 523776
        # ex_reduction.launch.end

        np.testing.assert_equal(a[0], sum(np.arange(1024)))


if __name__ == "__main__":
    unittest.main()
