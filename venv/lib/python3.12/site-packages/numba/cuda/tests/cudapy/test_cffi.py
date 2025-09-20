import numpy as np

from numba import cuda, types
from numba.cuda.testing import (skip_on_cudasim, test_data_dir, unittest,
                                CUDATestCase)
from numba.tests.support import skip_unless_cffi


@skip_unless_cffi
@skip_on_cudasim('Simulator does not support linking')
class TestCFFI(CUDATestCase):
    def test_from_buffer(self):
        import cffi
        ffi = cffi.FFI()

        link = str(test_data_dir / 'jitlink.ptx')
        sig = types.void(types.CPointer(types.int32))
        array_mutator = cuda.declare_device('array_mutator', sig)

        @cuda.jit(link=[link])
        def mutate_array(x):
            x_ptr = ffi.from_buffer(x)
            array_mutator(x_ptr)

        x = np.arange(2).astype(np.int32)
        mutate_array[1, 1](x)

        # The foreign function should have copied element 1 to element 0
        self.assertEqual(x[0], x[1])


if __name__ == '__main__':
    unittest.main()
