import numpy as np

from collections import namedtuple
from numba import void, int32, float32, float64
from numba import guvectorize
from numba import cuda
from numba.cuda.testing import skip_on_cudasim, CUDATestCase
import unittest
import warnings
from numba.core.errors import NumbaPerformanceWarning, TypingError
from numba.tests.support import override_config


def _get_matmulcore_gufunc(dtype=float32):
    @guvectorize([void(dtype[:, :], dtype[:, :], dtype[:, :])],
                 '(m,n),(n,p)->(m,p)',
                 target='cuda')
    def matmulcore(A, B, C):
        m, n = A.shape
        n, p = B.shape
        for i in range(m):
            for j in range(p):
                C[i, j] = 0
                for k in range(n):
                    C[i, j] += A[i, k] * B[k, j]

    return matmulcore


@skip_on_cudasim('ufunc API unsupported in the simulator')
class TestCUDAGufunc(CUDATestCase):

    def test_gufunc_small(self):

        gufunc = _get_matmulcore_gufunc()

        matrix_ct = 2
        A = np.arange(matrix_ct * 2 * 4, dtype=np.float32).reshape(matrix_ct, 2,
                                                                   4)
        B = np.arange(matrix_ct * 4 * 5, dtype=np.float32).reshape(matrix_ct, 4,
                                                                   5)

        C = gufunc(A, B)
        Gold = np.matmul(A, B)
        self.assertTrue(np.allclose(C, Gold))

    def test_gufunc_auto_transfer(self):

        gufunc = _get_matmulcore_gufunc()

        matrix_ct = 2
        A = np.arange(matrix_ct * 2 * 4, dtype=np.float32).reshape(matrix_ct, 2,
                                                                   4)
        B = np.arange(matrix_ct * 4 * 5, dtype=np.float32).reshape(matrix_ct, 4,
                                                                   5)

        dB = cuda.to_device(B)

        C = gufunc(A, dB).copy_to_host()
        Gold = np.matmul(A, B)
        self.assertTrue(np.allclose(C, Gold))

    def test_gufunc(self):

        gufunc = _get_matmulcore_gufunc()

        matrix_ct = 1001 # an odd number to test thread/block division in CUDA
        A = np.arange(matrix_ct * 2 * 4, dtype=np.float32).reshape(matrix_ct, 2,
                                                                   4)
        B = np.arange(matrix_ct * 4 * 5, dtype=np.float32).reshape(matrix_ct, 4,
                                                                   5)

        C = gufunc(A, B)
        Gold = np.matmul(A, B)
        self.assertTrue(np.allclose(C, Gold))

    def test_gufunc_hidim(self):

        gufunc = _get_matmulcore_gufunc()

        matrix_ct = 100 # an odd number to test thread/block division in CUDA
        A = np.arange(matrix_ct * 2 * 4, dtype=np.float32).reshape(4, 25, 2, 4)
        B = np.arange(matrix_ct * 4 * 5, dtype=np.float32).reshape(4, 25, 4, 5)

        C = gufunc(A, B)
        Gold = np.matmul(A, B)
        self.assertTrue(np.allclose(C, Gold))

    def test_gufunc_new_axis(self):

        gufunc = _get_matmulcore_gufunc(dtype=float64)

        X = np.random.randn(10, 3, 3)
        Y = np.random.randn(3, 3)

        gold = np.matmul(X, Y)

        res1 = gufunc(X, Y)
        np.testing.assert_allclose(gold, res1)

        res2 = gufunc(X, np.tile(Y, (10, 1, 1)))
        np.testing.assert_allclose(gold, res2)

    def test_gufunc_stream(self):

        gufunc = _get_matmulcore_gufunc()

        #cuda.driver.flush_pending_free()
        matrix_ct = 1001 # an odd number to test thread/block division in CUDA
        A = np.arange(matrix_ct * 2 * 4, dtype=np.float32).reshape(matrix_ct, 2,
                                                                   4)
        B = np.arange(matrix_ct * 4 * 5, dtype=np.float32).reshape(matrix_ct, 4,
                                                                   5)

        stream = cuda.stream()
        dA = cuda.to_device(A, stream)
        dB = cuda.to_device(B, stream)

        dC = cuda.device_array(shape=(1001, 2, 5), dtype=A.dtype, stream=stream)
        dC = gufunc(dA, dB, out=dC, stream=stream)
        C = dC.copy_to_host(stream=stream)
        stream.synchronize()

        Gold = np.matmul(A, B)

        self.assertTrue(np.allclose(C, Gold))

    def test_copy(self):

        @guvectorize([void(float32[:], float32[:])],
                     '(x)->(x)',
                     target='cuda')
        def copy(A, B):
            for i in range(B.size):
                B[i] = A[i]

        A = np.arange(10, dtype=np.float32) + 1
        B = np.zeros_like(A)
        copy(A, out=B)
        np.testing.assert_allclose(A, B)

    def test_copy_unspecified_return(self):
        # Ensure that behaviour is correct when the return type is not
        # specified in the signature.
        @guvectorize([(float32[:], float32[:])],
                     '(x)->(x)',
                     target='cuda')
        def copy(A, B):
            for i in range(B.size):
                B[i] = A[i]

        A = np.arange(10, dtype=np.float32) + 1
        B = np.zeros_like(A)
        copy(A, out=B)
        self.assertTrue(np.allclose(A, B))

    def test_copy_odd(self):

        @guvectorize([void(float32[:], float32[:])],
                     '(x)->(x)',
                     target='cuda')
        def copy(A, B):
            for i in range(B.size):
                B[i] = A[i]

        A = np.arange(11, dtype=np.float32) + 1
        B = np.zeros_like(A)
        copy(A, out=B)
        self.assertTrue(np.allclose(A, B))

    def test_copy2d(self):

        @guvectorize([void(float32[:, :], float32[:, :])],
                     '(x, y)->(x, y)',
                     target='cuda')
        def copy2d(A, B):
            for x in range(B.shape[0]):
                for y in range(B.shape[1]):
                    B[x, y] = A[x, y]

        A = np.arange(30, dtype=np.float32).reshape(5, 6) + 1
        B = np.zeros_like(A)
        copy2d(A, out=B)
        self.assertTrue(np.allclose(A, B))

    def test_not_supported_call_from_jit(self):
        # not supported
        @guvectorize([void(int32[:], int32[:])],
                     '(n)->(n)', target='cuda')
        def gufunc_copy(A, b):
            for i in range(A.shape[0]):
                b[i] = A[i]

        @cuda.jit
        def cuda_jit(A, b):
            return gufunc_copy(A, b)

        A = np.arange(1024 * 32).astype('int32')
        b = np.zeros_like(A)
        msg = "Untyped global name 'gufunc_copy'.*"
        with self.assertRaisesRegex(TypingError, msg):
            cuda_jit[1, 1](A, b)

    # Test inefficient use of the GPU where the inputs are all mapped onto a
    # single thread in a single block.
    def test_inefficient_launch_configuration(self):
        @guvectorize(['void(float32[:], float32[:], float32[:])'],
                     '(n),(n)->(n)', target='cuda')
        def numba_dist_cuda(a, b, dist):
            len = a.shape[0]
            for i in range(len):
                dist[i] = a[i] * b[i]

        a = np.random.rand(1024 * 32).astype('float32')
        b = np.random.rand(1024 * 32).astype('float32')
        dist = np.zeros(a.shape[0]).astype('float32')

        with override_config('CUDA_LOW_OCCUPANCY_WARNINGS', 1):
            with warnings.catch_warnings(record=True) as w:
                numba_dist_cuda(a, b, dist)
                self.assertEqual(w[0].category, NumbaPerformanceWarning)
                self.assertIn('Grid size', str(w[0].message))
                self.assertIn('low occupancy', str(w[0].message))

    def test_efficient_launch_configuration(self):
        @guvectorize(['void(float32[:], float32[:], float32[:])'],
                     '(n),(n)->(n)', nopython=True, target='cuda')
        def numba_dist_cuda2(a, b, dist):
            len = a.shape[0]
            for i in range(len):
                dist[i] = a[i] * b[i]

        a = np.random.rand(524288 * 2).astype('float32').\
            reshape((524288, 2))
        b = np.random.rand(524288 * 2).astype('float32').\
            reshape((524288, 2))
        dist = np.zeros_like(a)

        with override_config('CUDA_LOW_OCCUPANCY_WARNINGS', 1):
            with warnings.catch_warnings(record=True) as w:
                numba_dist_cuda2(a, b, dist)
                self.assertEqual(len(w), 0)

    def test_nopython_flag(self):

        def foo(A, B):
            pass

        # nopython = True is fine
        guvectorize([void(float32[:], float32[:])], '(x)->(x)', target='cuda',
                    nopython=True)(foo)

        # nopython = False is bad
        with self.assertRaises(TypeError) as raises:
            guvectorize([void(float32[:], float32[:])], '(x)->(x)',
                        target='cuda', nopython=False)(foo)
        self.assertEqual("nopython flag must be True", str(raises.exception))

    def test_invalid_flags(self):
        # Check invalid flags
        def foo(A, B):
            pass

        with self.assertRaises(TypeError) as raises:
            guvectorize([void(float32[:], float32[:])], '(x)->(x)',
                        target='cuda', what1=True, ever2=False)(foo)
        head = "The following target options are not supported:"
        msg = str(raises.exception)
        self.assertEqual(msg[:len(head)], head)
        items = msg[len(head):].strip().split(',')
        items = [i.strip("'\" ") for i in items]
        self.assertEqual(set(['what1', 'ever2']), set(items))

    def test_duplicated_output(self):
        @guvectorize([void(float32[:], float32[:])], '(x)->(x)', target='cuda')
        def foo(inp, out):
            pass  # intentionally empty; never executed

        inp = out = np.zeros(10, dtype=np.float32)
        with self.assertRaises(ValueError) as raises:
            foo(inp, out, out=out)

        msg = "cannot specify argument 'out' as both positional and keyword"
        self.assertEqual(str(raises.exception), msg)

    def check_tuple_arg(self, a, b):
        @guvectorize([(float64[:], float64[:], float64[:])], '(n),(n)->()',
                     target='cuda')
        def gu_reduce(x, y, r):
            s = 0
            for i in range(len(x)):
                s += x[i] * y[i]
            r[0] = s

        r = gu_reduce(a, b)
        expected = np.sum(np.asarray(a) * np.asarray(b), axis=1)
        np.testing.assert_equal(expected, r)

    def test_tuple_of_tuple_arg(self):
        a = ((1.0, 2.0, 3.0),
             (4.0, 5.0, 6.0))
        b = ((1.5, 2.5, 3.5),
             (4.5, 5.5, 6.5))
        self.check_tuple_arg(a, b)

    def test_tuple_of_namedtuple_arg(self):
        Point = namedtuple('Point', ('x', 'y', 'z'))
        a = (Point(x=1.0, y=2.0, z=3.0),
             Point(x=4.0, y=5.0, z=6.0))
        b = (Point(x=1.5, y=2.5, z=3.5),
             Point(x=4.5, y=5.5, z=6.5))
        self.check_tuple_arg(a, b)

    def test_tuple_of_array_arg(self):
        a = (np.asarray((1.0, 2.0, 3.0)),
             np.asarray((4.0, 5.0, 6.0)))
        b = (np.asarray((1.5, 2.5, 3.5)),
             np.asarray((4.5, 5.5, 6.5)))
        self.check_tuple_arg(a, b)

    def test_gufunc_name(self):
        gufunc = _get_matmulcore_gufunc()
        self.assertEqual(gufunc.__name__, 'matmulcore')

    def test_bad_return_type(self):
        with self.assertRaises(TypeError) as te:
            @guvectorize([int32(int32[:], int32[:])], '(m)->(m)', target='cuda')
            def f(x, y):
                pass

        msg = str(te.exception)
        self.assertIn('guvectorized functions cannot return values', msg)
        self.assertIn('specifies int32 return type', msg)

    def test_incorrect_number_of_pos_args(self):
        @guvectorize([(int32[:], int32[:], int32[:])],
                     '(m),(m)->(m)', target='cuda')
        def f(x, y, z):
            pass

        arr = np.arange(5)

        # Inputs only, too few
        with self.assertRaises(TypeError) as te:
            f(arr)

        msg = str(te.exception)
        self.assertIn('gufunc accepts 2 positional arguments', msg)
        self.assertIn('or 3 positional arguments', msg)
        self.assertIn('Got 1 positional argument.', msg)

        # Inputs and outputs, too many
        with self.assertRaises(TypeError) as te:
            f(arr, arr, arr, arr)

        msg = str(te.exception)
        self.assertIn('gufunc accepts 2 positional arguments', msg)
        self.assertIn('or 3 positional arguments', msg)
        self.assertIn('Got 4 positional arguments.', msg)


@skip_on_cudasim('ufunc API unsupported in the simulator')
class TestMultipleOutputs(CUDATestCase):
    def test_multiple_outputs_same_type_passed_in(self):
        @guvectorize([void(float32[:], float32[:], float32[:])],
                     '(x)->(x),(x)',
                     target='cuda')
        def copy(A, B, C):
            for i in range(B.size):
                B[i] = A[i]
                C[i] = A[i]

        A = np.arange(10, dtype=np.float32) + 1
        B = np.zeros_like(A)
        C = np.zeros_like(A)
        copy(A, B, C)
        np.testing.assert_allclose(A, B)
        np.testing.assert_allclose(A, C)

    def test_multiple_outputs_distinct_values(self):

        @guvectorize([void(float32[:], float32[:], float32[:])],
                     '(x)->(x),(x)',
                     target='cuda')
        def copy_and_double(A, B, C):
            for i in range(B.size):
                B[i] = A[i]
                C[i] = A[i] * 2

        A = np.arange(10, dtype=np.float32) + 1
        B = np.zeros_like(A)
        C = np.zeros_like(A)
        copy_and_double(A, B, C)
        np.testing.assert_allclose(A, B)
        np.testing.assert_allclose(A * 2, C)

    def test_multiple_output_allocation(self):
        @guvectorize([void(float32[:], float32[:], float32[:])],
                     '(x)->(x),(x)',
                     target='cuda')
        def copy_and_double(A, B, C):
            for i in range(B.size):
                B[i] = A[i]
                C[i] = A[i] * 2

        A = np.arange(10, dtype=np.float32) + 1
        B, C = copy_and_double(A)
        np.testing.assert_allclose(A, B)
        np.testing.assert_allclose(A * 2, C)

    def test_multiple_output_dtypes(self):

        @guvectorize([void(int32[:], int32[:], float64[:])],
                     '(x)->(x),(x)',
                     target='cuda')
        def copy_and_multiply(A, B, C):
            for i in range(B.size):
                B[i] = A[i]
                C[i] = A[i] * 1.5

        A = np.arange(10, dtype=np.int32) + 1
        B = np.zeros_like(A)
        C = np.zeros_like(A, dtype=np.float64)
        copy_and_multiply(A, B, C)
        np.testing.assert_allclose(A, B)
        np.testing.assert_allclose(A * np.float64(1.5), C)

    def test_incorrect_number_of_pos_args(self):
        @guvectorize([(int32[:], int32[:], int32[:], int32[:])],
                     '(m),(m)->(m),(m)', target='cuda')
        def f(x, y, z, w):
            pass

        arr = np.arange(5)

        # Inputs only, too few
        with self.assertRaises(TypeError) as te:
            f(arr)

        msg = str(te.exception)
        self.assertIn('gufunc accepts 2 positional arguments', msg)
        self.assertIn('or 4 positional arguments', msg)
        self.assertIn('Got 1 positional argument.', msg)

        # Inputs and outputs, too many
        with self.assertRaises(TypeError) as te:
            f(arr, arr, arr, arr, arr)

        msg = str(te.exception)
        self.assertIn('gufunc accepts 2 positional arguments', msg)
        self.assertIn('or 4 positional arguments', msg)
        self.assertIn('Got 5 positional arguments.', msg)


if __name__ == '__main__':
    unittest.main()
