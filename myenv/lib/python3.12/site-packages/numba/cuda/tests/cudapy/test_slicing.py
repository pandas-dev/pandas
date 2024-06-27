import numpy as np
from numba import cuda
from numba.cuda.testing import unittest, CUDATestCase


def foo(inp, out):
    for i in range(out.shape[0]):
        out[i] = inp[i]


def copy(inp, out):
    i = cuda.grid(1)
    cufoo(inp[i, :], out[i, :])


class TestCudaSlicing(CUDATestCase):
    def test_slice_as_arg(self):
        global cufoo
        cufoo = cuda.jit("void(int32[:], int32[:])", device=True)(foo)
        cucopy = cuda.jit("void(int32[:,:], int32[:,:])")(copy)

        inp = np.arange(100, dtype=np.int32).reshape(10, 10)
        out = np.zeros_like(inp)

        cucopy[1, 10](inp, out)

    def test_assign_empty_slice(self):
        # Issue #5017. Assigning to an empty slice should not result in a
        # CudaAPIError.
        N = 0
        a = range(N)
        arr = cuda.device_array(len(a))
        arr[:] = cuda.to_device(a)


if __name__ == '__main__':
    unittest.main()
