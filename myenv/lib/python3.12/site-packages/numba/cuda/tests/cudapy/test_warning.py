import numpy as np
from numba import cuda
from numba.cuda.testing import unittest, CUDATestCase, skip_on_cudasim
from numba.tests.support import linux_only, override_config
from numba.core.errors import NumbaPerformanceWarning
import warnings


@skip_on_cudasim('cudasim does not raise performance warnings')
class TestWarnings(CUDATestCase):
    def test_inefficient_launch_configuration(self):
        @cuda.jit
        def kernel():
            pass

        with override_config('CUDA_LOW_OCCUPANCY_WARNINGS', 1):
            with warnings.catch_warnings(record=True) as w:
                kernel[1, 1]()

        self.assertEqual(w[0].category, NumbaPerformanceWarning)
        self.assertIn('Grid size', str(w[0].message))
        self.assertIn('low occupancy', str(w[0].message))

    def test_efficient_launch_configuration(self):
        @cuda.jit
        def kernel():
            pass

        with override_config('CUDA_LOW_OCCUPANCY_WARNINGS', 1):
            with warnings.catch_warnings(record=True) as w:
                kernel[256, 256]()

        self.assertEqual(len(w), 0)

    def test_warn_on_host_array(self):
        @cuda.jit
        def foo(r, x):
            r[0] = x + 1

        N = 10
        arr_f32 = np.zeros(N, dtype=np.float32)
        with override_config('CUDA_WARN_ON_IMPLICIT_COPY', 1):
            with warnings.catch_warnings(record=True) as w:
                foo[1, N](arr_f32, N)

        self.assertEqual(w[0].category, NumbaPerformanceWarning)
        self.assertIn('Host array used in CUDA kernel will incur',
                      str(w[0].message))
        self.assertIn('copy overhead', str(w[0].message))

    def test_pinned_warn_on_host_array(self):
        @cuda.jit
        def foo(r, x):
            r[0] = x + 1

        N = 10
        ary = cuda.pinned_array(N, dtype=np.float32)

        with override_config('CUDA_WARN_ON_IMPLICIT_COPY', 1):
            with warnings.catch_warnings(record=True) as w:
                foo[1, N](ary, N)

        self.assertEqual(w[0].category, NumbaPerformanceWarning)
        self.assertIn('Host array used in CUDA kernel will incur',
                      str(w[0].message))
        self.assertIn('copy overhead', str(w[0].message))

    def test_nowarn_on_mapped_array(self):
        @cuda.jit
        def foo(r, x):
            r[0] = x + 1

        N = 10
        ary = cuda.mapped_array(N, dtype=np.float32)

        with override_config('CUDA_WARN_ON_IMPLICIT_COPY', 1):
            with warnings.catch_warnings(record=True) as w:
                foo[1, N](ary, N)

        self.assertEqual(len(w), 0)

    @linux_only
    def test_nowarn_on_managed_array(self):
        @cuda.jit
        def foo(r, x):
            r[0] = x + 1

        N = 10
        ary = cuda.managed_array(N, dtype=np.float32)

        with override_config('CUDA_WARN_ON_IMPLICIT_COPY', 1):
            with warnings.catch_warnings(record=True) as w:
                foo[1, N](ary, N)

        self.assertEqual(len(w), 0)

    def test_nowarn_on_device_array(self):
        @cuda.jit
        def foo(r, x):
            r[0] = x + 1

        N = 10
        ary = cuda.device_array(N, dtype=np.float32)

        with override_config('CUDA_WARN_ON_IMPLICIT_COPY', 1):
            with warnings.catch_warnings(record=True) as w:
                foo[1, N](ary, N)

        self.assertEqual(len(w), 0)

    def test_warn_on_debug_and_opt(self):
        with warnings.catch_warnings(record=True) as w:
            cuda.jit(debug=True, opt=True)

        self.assertEqual(len(w), 1)
        self.assertIn('not supported by CUDA', str(w[0].message))

    def test_warn_on_debug_and_opt_default(self):
        with warnings.catch_warnings(record=True) as w:
            cuda.jit(debug=True)

        self.assertEqual(len(w), 1)
        self.assertIn('not supported by CUDA', str(w[0].message))

    def test_no_warn_on_debug_and_no_opt(self):
        with warnings.catch_warnings(record=True) as w:
            cuda.jit(debug=True, opt=False)

        self.assertEqual(len(w), 0)

    def test_no_warn_with_no_debug_and_opt_kwargs(self):
        with warnings.catch_warnings(record=True) as w:
            cuda.jit()

        self.assertEqual(len(w), 0)


if __name__ == '__main__':
    unittest.main()
