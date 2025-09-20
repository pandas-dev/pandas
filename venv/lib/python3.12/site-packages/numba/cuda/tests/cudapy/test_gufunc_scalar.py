"""Example: sum each row using guvectorize

See Numpy documentation for detail about gufunc:
    http://docs.scipy.org/doc/numpy/reference/c-api.generalized-ufuncs.html
"""
import numpy as np
from numba import guvectorize, cuda
from numba.cuda.testing import skip_on_cudasim, CUDATestCase
import unittest


@skip_on_cudasim('ufunc API unsupported in the simulator')
class TestGUFuncScalar(CUDATestCase):
    def test_gufunc_scalar_output(self):
        #    function type:
        #        - has no void return type
        #        - array argument is one dimension fewer than the source array
        #        - scalar output is passed as a 1-element array.
        #
        #    signature: (n)->()
        #        - the function takes an array of n-element and output a scalar.

        @guvectorize(['void(int32[:], int32[:])'], '(n)->()', target='cuda')
        def sum_row(inp, out):
            tmp = 0.
            for i in range(inp.shape[0]):
                tmp += inp[i]
            out[0] = tmp

        # inp is (10000, 3)
        # out is (10000)
        # The outer (leftmost) dimension must match or numpy broadcasting
        # is performed. But, broadcasting on CUDA arrays is not supported.

        inp = np.arange(300, dtype=np.int32).reshape(100, 3)

        # invoke on CUDA with manually managed memory
        out1 = np.empty(100, dtype=inp.dtype)
        out2 = np.empty(100, dtype=inp.dtype)

        dev_inp = cuda.to_device(
            inp)                 # alloc and copy input data
        dev_out1 = cuda.to_device(out1, copy=False)   # alloc only

        sum_row(dev_inp, out=dev_out1)                # invoke the gufunc
        dev_out2 = sum_row(dev_inp)                   # invoke the gufunc

        dev_out1.copy_to_host(out1)                 # retrieve the result
        dev_out2.copy_to_host(out2)                 # retrieve the result

        # verify result
        for i in range(inp.shape[0]):
            self.assertTrue(out1[i] == inp[i].sum())
            self.assertTrue(out2[i] == inp[i].sum())

    def test_gufunc_scalar_output_bug(self):
        # Issue 2812: Error due to using input argument types as output argument
        @guvectorize(['void(int32, int32[:])'], '()->()', target='cuda')
        def twice(inp, out):
            out[0] = inp * 2

        self.assertEqual(twice(10), 20)
        arg = np.arange(10).astype(np.int32)
        self.assertPreciseEqual(twice(arg), arg * 2)

    def test_gufunc_scalar_input_saxpy(self):
        @guvectorize(['void(float32, float32[:], float32[:], float32[:])'],
                     '(),(t),(t)->(t)', target='cuda')
        def saxpy(a, x, y, out):
            for i in range(out.shape[0]):
                out[i] = a * x[i] + y[i]

        A = np.float32(2)
        X = np.arange(10, dtype=np.float32).reshape(5, 2)
        Y = np.arange(10, dtype=np.float32).reshape(5, 2)
        out = saxpy(A, X, Y)

        for j in range(5):
            for i in range(2):
                exp = A * X[j, i] + Y[j, i]
                self.assertTrue(exp == out[j, i])

        X = np.arange(10, dtype=np.float32)
        Y = np.arange(10, dtype=np.float32)
        out = saxpy(A, X, Y)

        for j in range(10):
            exp = A * X[j] + Y[j]
            self.assertTrue(exp == out[j], (exp, out[j]))

        A = np.arange(5, dtype=np.float32)
        X = np.arange(10, dtype=np.float32).reshape(5, 2)
        Y = np.arange(10, dtype=np.float32).reshape(5, 2)
        out = saxpy(A, X, Y)

        for j in range(5):
            for i in range(2):
                exp = A[j] * X[j, i] + Y[j, i]
                self.assertTrue(exp == out[j, i], (exp, out[j, i]))

    def test_gufunc_scalar_cast(self):
        @guvectorize(['void(int32, int32[:], int32[:])'], '(),(t)->(t)',
                     target='cuda')
        def foo(a, b, out):
            for i in range(b.size):
                out[i] = a * b[i]

        a = np.int64(2)  # type does not match signature (int32)
        b = np.arange(10).astype(np.int32)
        out = foo(a, b)
        np.testing.assert_equal(out, a * b)

        # test error
        a = np.array(a)
        da = cuda.to_device(a)
        self.assertEqual(da.dtype, np.int64)
        with self.assertRaises(TypeError) as raises:
            foo(da, b)

        self.assertIn("does not support .astype()", str(raises.exception))

    def test_gufunc_old_style_scalar_as_array(self):
        # Example from issue #2579
        @guvectorize(['void(int32[:],int32[:],int32[:])'], '(n),()->(n)',
                     target='cuda')
        def gufunc(x, y, res):
            for i in range(x.shape[0]):
                res[i] = x[i] + y[0]

        # Case 1
        a = np.array([1, 2, 3, 4], dtype=np.int32)
        b = np.array([2], dtype=np.int32)

        res = np.zeros(4, dtype=np.int32)

        expected = res.copy()
        expected = a + b

        gufunc(a, b, out=res)

        np.testing.assert_almost_equal(expected, res)

        # Case 2
        a = np.array([1, 2, 3, 4] * 2, dtype=np.int32).reshape(2, 4)
        b = np.array([2, 10], dtype=np.int32)

        res = np.zeros((2, 4), dtype=np.int32)

        expected = res.copy()
        expected[0] = a[0] + b[0]
        expected[1] = a[1] + b[1]

        gufunc(a, b, res)

        np.testing.assert_almost_equal(expected, res)


if __name__ == '__main__':
    unittest.main()
