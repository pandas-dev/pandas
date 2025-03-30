import numpy as np

from numba.cuda.testing import skip_on_cudasim, CUDATestCase
from numba import cuda, float64
import unittest


def kernel_func(x):
    x[0] = 1


def device_func(x, y, z):
    return x * y + z


# Fragments of code that are removed from kernel_func's PTX when optimization
# is on
removed_by_opt = ( '__local_depot0', 'call.uni', 'st.param.b64')


@skip_on_cudasim('Simulator does not optimize code')
class TestOptimization(CUDATestCase):
    def test_eager_opt(self):
        # Optimization should occur by default
        sig = (float64[::1],)
        kernel = cuda.jit(sig)(kernel_func)
        ptx = kernel.inspect_asm()

        for fragment in removed_by_opt:
            with self.subTest(fragment=fragment):
                self.assertNotIn(fragment, ptx[sig])

    def test_eager_noopt(self):
        # Optimization disabled
        sig = (float64[::1],)
        kernel = cuda.jit(sig, opt=False)(kernel_func)
        ptx = kernel.inspect_asm()

        for fragment in removed_by_opt:
            with self.subTest(fragment=fragment):
                self.assertIn(fragment, ptx[sig])

    def test_lazy_opt(self):
        # Optimization should occur by default
        kernel = cuda.jit(kernel_func)
        x = np.zeros(1, dtype=np.float64)
        kernel[1, 1](x)

        # Grab the PTX for the one definition that has just been jitted
        ptx = next(iter(kernel.inspect_asm().items()))[1]

        for fragment in removed_by_opt:
            with self.subTest(fragment=fragment):
                self.assertNotIn(fragment, ptx)

    def test_lazy_noopt(self):
        # Optimization disabled
        kernel = cuda.jit(opt=False)(kernel_func)
        x = np.zeros(1, dtype=np.float64)
        kernel[1, 1](x)

        # Grab the PTX for the one definition that has just been jitted
        ptx = next(iter(kernel.inspect_asm().items()))[1]

        for fragment in removed_by_opt:
            with self.subTest(fragment=fragment):
                self.assertIn(fragment, ptx)

    def test_device_opt(self):
        # Optimization should occur by default
        sig = (float64, float64, float64)
        device = cuda.jit(sig, device=True)(device_func)
        ptx = device.inspect_asm(sig)
        self.assertIn('fma.rn.f64', ptx)

    def test_device_noopt(self):
        # Optimization disabled
        sig = (float64, float64, float64)
        device = cuda.jit(sig, device=True, opt=False)(device_func)
        ptx = device.inspect_asm(sig)
        # Fused-multiply adds should be disabled when not optimizing
        self.assertNotIn('fma.rn.f64', ptx)


if __name__ == '__main__':
    unittest.main()
