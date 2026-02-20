import pickle
import numpy as np
from numba import cuda, vectorize
from numba.core import types
from numba.cuda.testing import skip_on_cudasim, CUDATestCase
import unittest
from numba.np import numpy_support


@skip_on_cudasim('pickling not supported in CUDASIM')
class TestPickle(CUDATestCase):

    def check_call(self, callee):
        arr = np.array([100])
        expected = callee[1, 1](arr)

        # serialize and rebuild
        foo1 = pickle.loads(pickle.dumps(callee))
        del callee
        # call rebuild function
        got1 = foo1[1, 1](arr)
        np.testing.assert_equal(got1, expected)
        del got1

        # test serialization of previously serialized object
        foo2 = pickle.loads(pickle.dumps(foo1))
        del foo1
        # call rebuild function
        got2 = foo2[1, 1](arr)
        np.testing.assert_equal(got2, expected)
        del got2

        # test propagation of thread, block config
        foo3 = pickle.loads(pickle.dumps(foo2[5, 8]))
        del foo2
        self.assertEqual(foo3.griddim, (5, 1, 1))
        self.assertEqual(foo3.blockdim, (8, 1, 1))

    def test_pickling_jit_typing(self):
        @cuda.jit(device=True)
        def inner(a):
            return a + 1

        @cuda.jit('void(intp[:])')
        def foo(arr):
            arr[0] = inner(arr[0])

        self.check_call(foo)

    def test_pickling_jit(self):

        @cuda.jit(device=True)
        def inner(a):
            return a + 1

        @cuda.jit
        def foo(arr):
            arr[0] = inner(arr[0])

        self.check_call(foo)

    def test_pickling_vectorize(self):
        @vectorize(['intp(intp)', 'float64(float64)'], target='cuda')
        def cuda_vect(x):
            return x * 2

        # accommodate int representations in np.arange
        npty = numpy_support.as_dtype(types.intp)
        # get expected result
        ary = np.arange(10, dtype=npty)
        expected = cuda_vect(ary)
        # first pickle
        foo1 = pickle.loads(pickle.dumps(cuda_vect))
        del cuda_vect
        got1 = foo1(ary)
        np.testing.assert_equal(expected, got1)
        # second pickle
        foo2 = pickle.loads(pickle.dumps(foo1))
        del foo1
        got2 = foo2(ary)
        np.testing.assert_equal(expected, got2)


if __name__ == '__main__':
    unittest.main()
