import numpy as np

from numba.cuda.testing import unittest, CUDATestCase
from numba.cuda.testing import skip_on_cudasim, skip_unless_cudasim
from numba import config, cuda


if config.ENABLE_CUDASIM:
    ARRAY_LIKE_FUNCTIONS = (cuda.device_array_like, cuda.pinned_array_like)
else:
    ARRAY_LIKE_FUNCTIONS = (cuda.device_array_like, cuda.mapped_array_like,
                            cuda.pinned_array_like)


class TestCudaArray(CUDATestCase):
    def test_gpu_array_zero_length(self):
        x = np.arange(0)
        dx = cuda.to_device(x)
        hx = dx.copy_to_host()
        self.assertEqual(x.shape, dx.shape)
        self.assertEqual(x.size, dx.size)
        self.assertEqual(x.shape, hx.shape)
        self.assertEqual(x.size, hx.size)

    def test_null_shape(self):
        null_shape = ()
        shape1 = cuda.device_array(()).shape
        shape2 = cuda.device_array_like(np.ndarray(())).shape
        self.assertEqual(shape1, null_shape)
        self.assertEqual(shape2, null_shape)

    def test_gpu_array_strided(self):

        @cuda.jit('void(double[:])')
        def kernel(x):
            i = cuda.grid(1)
            if i < x.shape[0]:
                x[i] = i

        x = np.arange(10, dtype=np.double)
        y = np.ndarray(shape=10 * 8, buffer=x, dtype=np.byte)
        z = np.ndarray(9, buffer=y[4:-4], dtype=np.double)
        kernel[10, 10](z)
        self.assertTrue(np.allclose(z, list(range(9))))

    def test_gpu_array_interleaved(self):

        @cuda.jit('void(double[:], double[:])')
        def copykernel(x, y):
            i = cuda.grid(1)
            if i < x.shape[0]:
                x[i] = i
                y[i] = i

        x = np.arange(10, dtype=np.double)
        y = x[:-1:2]
        # z = x[1::2]
        # n = y.size
        try:
            cuda.devicearray.auto_device(y)
        except ValueError:
            pass
        else:
            raise AssertionError("Should raise exception complaining the "
                                 "contiguous-ness of the array.")
            # Should we handle this use case?
            # assert z.size == y.size
            # copykernel[1, n](y, x)
            # print(y, z)
            # assert np.all(y == z)
            # assert np.all(y == list(range(n)))

    def test_auto_device_const(self):
        d, _ = cuda.devicearray.auto_device(2)
        self.assertTrue(np.all(d.copy_to_host() == np.array(2)))

    def _test_array_like_same(self, like_func, array):
        """
        Tests of *_array_like where shape, strides, dtype, and flags should
        all be equal.
        """
        array_like = like_func(array)
        self.assertEqual(array.shape, array_like.shape)
        self.assertEqual(array.strides, array_like.strides)
        self.assertEqual(array.dtype, array_like.dtype)
        self.assertEqual(array.flags['C_CONTIGUOUS'],
                         array_like.flags['C_CONTIGUOUS'])
        self.assertEqual(array.flags['F_CONTIGUOUS'],
                         array_like.flags['F_CONTIGUOUS'])

    def test_array_like_1d(self):
        d_a = cuda.device_array(10, order='C')
        for like_func in ARRAY_LIKE_FUNCTIONS:
            with self.subTest(like_func=like_func):
                self._test_array_like_same(like_func, d_a)

    def test_array_like_2d(self):
        d_a = cuda.device_array((10, 12), order='C')
        for like_func in ARRAY_LIKE_FUNCTIONS:
            with self.subTest(like_func=like_func):
                self._test_array_like_same(like_func, d_a)

    def test_array_like_2d_transpose(self):
        d_a = cuda.device_array((10, 12), order='C')
        for like_func in ARRAY_LIKE_FUNCTIONS:
            with self.subTest(like_func=like_func):
                self._test_array_like_same(like_func, d_a)

    def test_array_like_3d(self):
        d_a = cuda.device_array((10, 12, 14), order='C')
        for like_func in ARRAY_LIKE_FUNCTIONS:
            with self.subTest(like_func=like_func):
                self._test_array_like_same(like_func, d_a)

    def test_array_like_1d_f(self):
        d_a = cuda.device_array(10, order='F')
        for like_func in ARRAY_LIKE_FUNCTIONS:
            with self.subTest(like_func=like_func):
                self._test_array_like_same(like_func, d_a)

    def test_array_like_2d_f(self):
        d_a = cuda.device_array((10, 12), order='F')
        for like_func in ARRAY_LIKE_FUNCTIONS:
            with self.subTest(like_func=like_func):
                self._test_array_like_same(like_func, d_a)

    def test_array_like_2d_f_transpose(self):
        d_a = cuda.device_array((10, 12), order='F')
        for like_func in ARRAY_LIKE_FUNCTIONS:
            with self.subTest(like_func=like_func):
                self._test_array_like_same(like_func, d_a)

    def test_array_like_3d_f(self):
        d_a = cuda.device_array((10, 12, 14), order='F')
        for like_func in ARRAY_LIKE_FUNCTIONS:
            with self.subTest(like_func=like_func):
                self._test_array_like_same(like_func, d_a)

    def _test_array_like_view(self, like_func, view, d_view):
        """
        Tests of device_array_like where the original array is a view - the
        strides should not be equal because a contiguous array is expected.
        """
        nb_like = like_func(d_view)
        self.assertEqual(d_view.shape, nb_like.shape)
        self.assertEqual(d_view.dtype, nb_like.dtype)

        # Use NumPy as a reference for the expected strides
        np_like = np.zeros_like(view)
        self.assertEqual(nb_like.strides, np_like.strides)
        self.assertEqual(nb_like.flags['C_CONTIGUOUS'],
                         np_like.flags['C_CONTIGUOUS'])
        self.assertEqual(nb_like.flags['F_CONTIGUOUS'],
                         np_like.flags['F_CONTIGUOUS'])

    def test_array_like_1d_view(self):
        shape = 10
        view = np.zeros(shape)[::2]
        d_view = cuda.device_array(shape)[::2]
        for like_func in ARRAY_LIKE_FUNCTIONS:
            with self.subTest(like_func=like_func):
                self._test_array_like_view(like_func, view, d_view)

    def test_array_like_1d_view_f(self):
        shape = 10
        view = np.zeros(shape, order='F')[::2]
        d_view = cuda.device_array(shape, order='F')[::2]
        for like_func in ARRAY_LIKE_FUNCTIONS:
            with self.subTest(like_func=like_func):
                self._test_array_like_view(like_func, view, d_view)

    def test_array_like_2d_view(self):
        shape = (10, 12)
        view = np.zeros(shape)[::2, ::2]
        d_view = cuda.device_array(shape)[::2, ::2]
        for like_func in ARRAY_LIKE_FUNCTIONS:
            with self.subTest(like_func=like_func):
                self._test_array_like_view(like_func, view, d_view)

    def test_array_like_2d_view_f(self):
        shape = (10, 12)
        view = np.zeros(shape, order='F')[::2, ::2]
        d_view = cuda.device_array(shape, order='F')[::2, ::2]
        for like_func in ARRAY_LIKE_FUNCTIONS:
            with self.subTest(like_func=like_func):
                self._test_array_like_view(like_func, view, d_view)

    @skip_on_cudasim('Numba and NumPy stride semantics differ for transpose')
    def test_array_like_2d_view_transpose_device(self):
        shape = (10, 12)
        d_view = cuda.device_array(shape)[::2, ::2].T
        for like_func in ARRAY_LIKE_FUNCTIONS:
            with self.subTest(like_func=like_func):
                # This is a special case (see issue #4974) because creating the
                # transpose creates a new contiguous allocation with different
                # strides.  In this case, rather than comparing against NumPy,
                # we can only compare against expected values.
                like = like_func(d_view)
                self.assertEqual(d_view.shape, like.shape)
                self.assertEqual(d_view.dtype, like.dtype)
                self.assertEqual((40, 8), like.strides)
                self.assertTrue(like.flags['C_CONTIGUOUS'])
                self.assertFalse(like.flags['F_CONTIGUOUS'])

    @skip_unless_cudasim('Numba and NumPy stride semantics differ for '
                         'transpose')
    def test_array_like_2d_view_transpose_simulator(self):
        shape = (10, 12)
        view = np.zeros(shape)[::2, ::2].T
        d_view = cuda.device_array(shape)[::2, ::2].T
        for like_func in ARRAY_LIKE_FUNCTIONS:
            with self.subTest(like_func=like_func):
                # On the simulator, the transpose has different strides to on a
                # CUDA device (See issue #4974). Here we can compare strides
                # against NumPy as a reference.
                np_like = np.zeros_like(view)
                nb_like = like_func(d_view)
                self.assertEqual(d_view.shape, nb_like.shape)
                self.assertEqual(d_view.dtype, nb_like.dtype)
                self.assertEqual(np_like.strides, nb_like.strides)
                self.assertEqual(np_like.flags['C_CONTIGUOUS'],
                                 nb_like.flags['C_CONTIGUOUS'])
                self.assertEqual(np_like.flags['F_CONTIGUOUS'],
                                 nb_like.flags['F_CONTIGUOUS'])

    def test_array_like_2d_view_f_transpose(self):
        shape = (10, 12)
        view = np.zeros(shape, order='F')[::2, ::2].T
        d_view = cuda.device_array(shape, order='F')[::2, ::2].T
        for like_func in ARRAY_LIKE_FUNCTIONS:
            with self.subTest(like_func=like_func):
                self._test_array_like_view(like_func, view, d_view)

    @skip_on_cudasim('Kernel overloads not created in the simulator')
    def test_issue_4628(self):
        # CUDA Device arrays were reported as always being typed with 'A' order
        # so launching the kernel with a host array and then a device array
        # resulted in two overloads being compiled - one for 'C' order from
        # the host array, and one for 'A' order from the device array. With the
        # resolution of this issue, the order of the device array is also 'C',
        # so after the kernel launches there should only be one overload of
        # the function.
        @cuda.jit
        def func(A, out):
            i = cuda.grid(1)
            out[i] = A[i] * 2

        n = 128
        a = np.ones((n,))
        d_a = cuda.to_device(a)
        result = np.zeros((n,))

        func[1, 128](a, result)
        func[1, 128](d_a, result)

        self.assertEqual(1, len(func.overloads))


if __name__ == '__main__':
    unittest.main()
