import numpy as np
from numba import cuda
from numba.cuda.testing import CUDATestCase
import unittest


def reinterpret_array_type(byte_arr, start, stop, output):
    # Tested with just one thread
    val = byte_arr[start:stop].view(np.int32)[0]
    output[0] = val


class TestCudaArrayMethods(CUDATestCase):
    def test_reinterpret_array_type(self):
        """
        Reinterpret byte array as int32 in the GPU.
        """
        pyfunc = reinterpret_array_type
        kernel = cuda.jit(pyfunc)

        byte_arr = np.arange(256, dtype=np.uint8)
        itemsize = np.dtype(np.int32).itemsize
        for start in range(0, 256, itemsize):
            stop = start + itemsize
            expect = byte_arr[start:stop].view(np.int32)[0]

            output = np.zeros(1, dtype=np.int32)
            kernel[1, 1](byte_arr, start, stop, output)

            got = output[0]
            self.assertEqual(expect, got)


if __name__ == '__main__':
    unittest.main()
