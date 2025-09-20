from numba import cuda
from numba.core.errors import TypingError
from numba.cuda.testing import unittest, CUDATestCase, skip_on_cudasim


def noop(x):
    pass


class TestJitErrors(CUDATestCase):
    """
    Test compile-time errors with @jit.
    """

    def test_too_many_dims(self):
        kernfunc = cuda.jit(noop)

        with self.assertRaises(ValueError) as raises:
            kernfunc[(1, 2, 3, 4), (5, 6)]
        self.assertIn("griddim must be a sequence of 1, 2 or 3 integers, "
                      "got [1, 2, 3, 4]",
                      str(raises.exception))

        with self.assertRaises(ValueError) as raises:
            kernfunc[(1, 2,), (3, 4, 5, 6)]
        self.assertIn("blockdim must be a sequence of 1, 2 or 3 integers, "
                      "got [3, 4, 5, 6]",
                      str(raises.exception))

    def test_non_integral_dims(self):
        kernfunc = cuda.jit(noop)

        with self.assertRaises(TypeError) as raises:
            kernfunc[2.0, 3]
        self.assertIn("griddim must be a sequence of integers, got [2.0]",
                      str(raises.exception))

        with self.assertRaises(TypeError) as raises:
            kernfunc[2, 3.0]
        self.assertIn("blockdim must be a sequence of integers, got [3.0]",
                      str(raises.exception))

    def _test_unconfigured(self, kernfunc):
        with self.assertRaises(ValueError) as raises:
            kernfunc(0)
        self.assertIn("launch configuration was not specified",
                      str(raises.exception))

    def test_unconfigured_typed_cudakernel(self):
        kernfunc = cuda.jit("void(int32)")(noop)
        self._test_unconfigured(kernfunc)

    def test_unconfigured_untyped_cudakernel(self):
        kernfunc = cuda.jit(noop)
        self._test_unconfigured(kernfunc)

    @skip_on_cudasim('TypingError does not occur on simulator')
    def test_typing_error(self):
        # see #5860, this is present to catch changes to error reporting
        # accidentally breaking the CUDA target

        @cuda.jit(device=True)
        def dev_func(x):
            # floor is deliberately not imported for the purpose of this test.
            return floor(x)  # noqa: F821

        @cuda.jit
        def kernel_func():
            dev_func(1.5)

        with self.assertRaises(TypingError) as raises:
            kernel_func[1, 1]()
        excstr = str(raises.exception)
        self.assertIn("resolving callee type: type(CUDADispatcher", excstr)
        self.assertIn("NameError: name 'floor' is not defined", excstr)


if __name__ == '__main__':
    unittest.main()
