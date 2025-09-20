import unittest

from numba.cuda.testing import CUDATestCase, skip_on_cudasim
from numba.tests.support import captured_stdout


@skip_on_cudasim("cudasim doesn't support cuda import at non-top-level")
class TestUFunc(CUDATestCase):
    """
    Test calling a UFunc
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

    def test_ex_cuda_ufunc_call(self):
        # ex_cuda_ufunc.begin
        import numpy as np
        from numba import cuda

        # A kernel calling a ufunc (sin, in this case)
        @cuda.jit
        def f(r, x):
            # Compute sin(x) with result written to r
            np.sin(x, r)

        # Declare input and output arrays
        x = np.arange(10, dtype=np.float32) - 5
        r = np.zeros_like(x)

        # Launch kernel that calls the ufunc
        f[1, 1](r, x)

        # A quick sanity check demonstrating equality of the sine computed by
        # the sin ufunc inside the kernel, and NumPy's sin ufunc
        np.testing.assert_allclose(r, np.sin(x))
        # ex_cuda_ufunc.end


if __name__ == "__main__":
    unittest.main()
