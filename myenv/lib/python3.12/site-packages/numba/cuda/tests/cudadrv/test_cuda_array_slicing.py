from itertools import product

import numpy as np

from numba import cuda
from numba.cuda.testing import unittest, CUDATestCase, skip_on_cudasim
from unittest.mock import patch


class CudaArrayIndexing(CUDATestCase):
    def test_index_1d(self):
        arr = np.arange(10)
        darr = cuda.to_device(arr)
        x, = arr.shape
        for i in range(-x, x):
            self.assertEqual(arr[i], darr[i])
        with self.assertRaises(IndexError):
            darr[-x - 1]
        with self.assertRaises(IndexError):
            darr[x]

    def test_index_2d(self):
        arr = np.arange(3 * 4).reshape(3, 4)
        darr = cuda.to_device(arr)
        x, y = arr.shape
        for i in range(-x, x):
            for j in range(-y, y):
                self.assertEqual(arr[i, j], darr[i, j])
        with self.assertRaises(IndexError):
            darr[-x - 1, 0]
        with self.assertRaises(IndexError):
            darr[x, 0]
        with self.assertRaises(IndexError):
            darr[0, -y - 1]
        with self.assertRaises(IndexError):
            darr[0, y]

    def test_index_3d(self):
        arr = np.arange(3 * 4 * 5).reshape(3, 4, 5)
        darr = cuda.to_device(arr)
        x, y, z = arr.shape
        for i in range(-x, x):
            for j in range(-y, y):
                for k in range(-z, z):
                    self.assertEqual(arr[i, j, k], darr[i, j, k])
        with self.assertRaises(IndexError):
            darr[-x - 1, 0, 0]
        with self.assertRaises(IndexError):
            darr[x, 0, 0]
        with self.assertRaises(IndexError):
            darr[0, -y - 1, 0]
        with self.assertRaises(IndexError):
            darr[0, y, 0]
        with self.assertRaises(IndexError):
            darr[0, 0, -z - 1]
        with self.assertRaises(IndexError):
            darr[0, 0, z]


class CudaArrayStridedSlice(CUDATestCase):

    def test_strided_index_1d(self):
        arr = np.arange(10)
        darr = cuda.to_device(arr)
        for i in range(arr.size):
            np.testing.assert_equal(arr[i::2], darr[i::2].copy_to_host())

    def test_strided_index_2d(self):
        arr = np.arange(6 * 7).reshape(6, 7)
        darr = cuda.to_device(arr)

        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                np.testing.assert_equal(arr[i::2, j::2],
                                        darr[i::2, j::2].copy_to_host())

    def test_strided_index_3d(self):
        arr = np.arange(6 * 7 * 8).reshape(6, 7, 8)
        darr = cuda.to_device(arr)

        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                for k in range(arr.shape[2]):
                    np.testing.assert_equal(
                        arr[i::2, j::2, k::2],
                        darr[i::2, j::2, k::2].copy_to_host())


class CudaArraySlicing(CUDATestCase):
    def test_prefix_1d(self):
        arr = np.arange(5)
        darr = cuda.to_device(arr)
        for i in range(arr.size):
            expect = arr[i:]
            got = darr[i:].copy_to_host()
            self.assertTrue(np.all(expect == got))

    def test_prefix_2d(self):
        arr = np.arange(3 ** 2).reshape(3, 3)
        darr = cuda.to_device(arr)
        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                expect = arr[i:, j:]
                sliced = darr[i:, j:]
                self.assertEqual(expect.shape, sliced.shape)
                self.assertEqual(expect.strides, sliced.strides)
                got = sliced.copy_to_host()
                self.assertTrue(np.all(expect == got))

    def test_select_3d_first_two_dim(self):
        arr = np.arange(3 * 4 * 5).reshape(3, 4, 5)
        darr = cuda.to_device(arr)
        # Select first dimension
        for i in range(arr.shape[0]):
            expect = arr[i]
            sliced = darr[i]
            self.assertEqual(expect.shape, sliced.shape)
            self.assertEqual(expect.strides, sliced.strides)
            got = sliced.copy_to_host()
            self.assertTrue(np.all(expect == got))
        # Select second dimension
        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                expect = arr[i, j]
                sliced = darr[i, j]
                self.assertEqual(expect.shape, sliced.shape)
                self.assertEqual(expect.strides, sliced.strides)
                got = sliced.copy_to_host()
                self.assertTrue(np.all(expect == got))

    def test_select_f(self):
        a = np.arange(5 * 6 * 7).reshape(5, 6, 7, order='F')
        da = cuda.to_device(a)

        for i in range(a.shape[0]):
            for j in range(a.shape[1]):
                self.assertTrue(np.array_equal(da[i, j, :].copy_to_host(),
                                               a[i, j, :]))
            for j in range(a.shape[2]):
                self.assertTrue(np.array_equal(da[i, :, j].copy_to_host(),
                                               a[i, :, j]))
        for i in range(a.shape[1]):
            for j in range(a.shape[2]):
                self.assertTrue(np.array_equal(da[:, i, j].copy_to_host(),
                                               a[:, i, j]))

    def test_select_c(self):
        a = np.arange(5 * 6 * 7).reshape(5, 6, 7, order='C')
        da = cuda.to_device(a)

        for i in range(a.shape[0]):
            for j in range(a.shape[1]):
                self.assertTrue(np.array_equal(da[i, j, :].copy_to_host(),
                                               a[i, j, :]))
            for j in range(a.shape[2]):
                self.assertTrue(np.array_equal(da[i, :, j].copy_to_host(),
                                               a[i, :, j]))
        for i in range(a.shape[1]):
            for j in range(a.shape[2]):
                self.assertTrue(np.array_equal(da[:, i, j].copy_to_host(),
                                               a[:, i, j]))

    def test_prefix_select(self):
        arr = np.arange(5 * 7).reshape(5, 7, order='F')

        darr = cuda.to_device(arr)
        self.assertTrue(np.all(darr[:1, 1].copy_to_host() == arr[:1, 1]))

    def test_negative_slicing_1d(self):
        arr = np.arange(10)
        darr = cuda.to_device(arr)
        for i, j in product(range(-10, 10), repeat=2):
            np.testing.assert_array_equal(arr[i:j],
                                          darr[i:j].copy_to_host())

    def test_negative_slicing_2d(self):
        arr = np.arange(12).reshape(3, 4)
        darr = cuda.to_device(arr)
        for x, y, w, s in product(range(-4, 4), repeat=4):
            np.testing.assert_array_equal(arr[x:y, w:s],
                                          darr[x:y, w:s].copy_to_host())

    def test_empty_slice_1d(self):
        arr = np.arange(5)
        darr = cuda.to_device(arr)
        for i in range(darr.shape[0]):
            np.testing.assert_array_equal(darr[i:i].copy_to_host(), arr[i:i])
        # empty slice of empty slice
        self.assertFalse(darr[:0][:0].copy_to_host())
        # out-of-bound slice just produces empty slices
        np.testing.assert_array_equal(darr[:0][:1].copy_to_host(),
                                      arr[:0][:1])
        np.testing.assert_array_equal(darr[:0][-1:].copy_to_host(),
                                      arr[:0][-1:])

    def test_empty_slice_2d(self):
        arr = np.arange(5 * 7).reshape(5, 7)
        darr = cuda.to_device(arr)
        np.testing.assert_array_equal(darr[:0].copy_to_host(), arr[:0])
        np.testing.assert_array_equal(darr[3, :0].copy_to_host(), arr[3, :0])
        # empty slice of empty slice
        self.assertFalse(darr[:0][:0].copy_to_host())
        # out-of-bound slice just produces empty slices
        np.testing.assert_array_equal(darr[:0][:1].copy_to_host(), arr[:0][:1])
        np.testing.assert_array_equal(darr[:0][-1:].copy_to_host(),
                                      arr[:0][-1:])


class CudaArraySetting(CUDATestCase):
    """
    Most of the slicing logic is tested in the cases above, so these
    tests focus on the setting logic.
    """

    def test_scalar(self):
        arr = np.arange(5 * 7).reshape(5, 7)
        darr = cuda.to_device(arr)
        arr[2, 2] = 500
        darr[2, 2] = 500
        np.testing.assert_array_equal(darr.copy_to_host(), arr)

    def test_rank(self):
        arr = np.arange(5 * 7).reshape(5, 7)
        darr = cuda.to_device(arr)
        arr[2] = 500
        darr[2] = 500
        np.testing.assert_array_equal(darr.copy_to_host(), arr)

    def test_broadcast(self):
        arr = np.arange(5 * 7).reshape(5, 7)
        darr = cuda.to_device(arr)
        arr[:, 2] = 500
        darr[:, 2] = 500
        np.testing.assert_array_equal(darr.copy_to_host(), arr)

    def test_array_assign_column(self):
        arr = np.arange(5 * 7).reshape(5, 7)
        darr = cuda.to_device(arr)
        _400 = np.full(shape=7, fill_value=400)
        arr[2] = _400
        darr[2] = _400
        np.testing.assert_array_equal(darr.copy_to_host(), arr)

    def test_array_assign_row(self):
        arr = np.arange(5 * 7).reshape(5, 7)
        darr = cuda.to_device(arr)
        _400 = np.full(shape=5, fill_value=400)
        arr[:, 2] = _400
        darr[:, 2] = _400
        np.testing.assert_array_equal(darr.copy_to_host(), arr)

    def test_array_assign_subarray(self):
        arr = np.arange(5 * 6 * 7).reshape(5, 6, 7)
        darr = cuda.to_device(arr)
        _400 = np.full(shape=(6, 7), fill_value=400)
        arr[2] = _400
        darr[2] = _400
        np.testing.assert_array_equal(darr.copy_to_host(), arr)

    def test_array_assign_deep_subarray(self):
        arr = np.arange(5 * 6 * 7 * 8).reshape(5, 6, 7, 8)
        darr = cuda.to_device(arr)
        _400 = np.full(shape=(5, 6, 8), fill_value=400)
        arr[:, :, 2] = _400
        darr[:, :, 2] = _400
        np.testing.assert_array_equal(darr.copy_to_host(), arr)

    def test_array_assign_all(self):
        arr = np.arange(5 * 7).reshape(5, 7)
        darr = cuda.to_device(arr)
        _400 = np.full(shape=(5, 7), fill_value=400)
        arr[:] = _400
        darr[:] = _400
        np.testing.assert_array_equal(darr.copy_to_host(), arr)

    def test_strides(self):
        arr = np.ones(20)
        darr = cuda.to_device(arr)
        arr[::2] = 500
        darr[::2] = 500
        np.testing.assert_array_equal(darr.copy_to_host(), arr)

    def test_incompatible_highdim(self):
        darr = cuda.to_device(np.arange(5 * 7))

        with self.assertRaises(ValueError) as e:
            darr[:] = np.ones(shape=(1, 2, 3))

        self.assertIn(
            member=str(e.exception),
            container=[
                "Can't assign 3-D array to 1-D self",  # device
                "could not broadcast input array from shape (2,3) "
                "into shape (35,)",  # simulator, NP >= 1.20
            ])

    def test_incompatible_shape(self):
        darr = cuda.to_device(np.arange(5))

        with self.assertRaises(ValueError) as e:
            darr[:] = [1, 3]

        self.assertIn(
            member=str(e.exception),
            container=[
                "Can't copy sequence with size 2 to array axis 0 with "
                "dimension 5",  # device
                "could not broadcast input array from shape (2,) into "
                "shape (5,)",   # simulator, NP >= 1.20
            ])

    @skip_on_cudasim('cudasim does not use streams and operates synchronously')
    def test_sync(self):
        # There should be a synchronization when no stream is supplied
        darr = cuda.to_device(np.arange(5))

        with patch.object(cuda.cudadrv.driver.Stream, 'synchronize',
                          return_value=None) as mock_sync:
            darr[0] = 10

        mock_sync.assert_called_once()

    @skip_on_cudasim('cudasim does not use streams and operates synchronously')
    def test_no_sync_default_stream(self):
        # There should not be a synchronization when the array has a default
        # stream, whether it is the default stream, the legacy default stream,
        # the per-thread default stream, or another stream.
        streams = (cuda.stream(), cuda.default_stream(),
                   cuda.legacy_default_stream(),
                   cuda.per_thread_default_stream())

        for stream in streams:
            darr = cuda.to_device(np.arange(5), stream=stream)

            with patch.object(cuda.cudadrv.driver.Stream, 'synchronize',
                              return_value=None) as mock_sync:
                darr[0] = 10

            mock_sync.assert_not_called()

    @skip_on_cudasim('cudasim does not use streams and operates synchronously')
    def test_no_sync_supplied_stream(self):
        # There should not be a synchronization when a stream is supplied for
        # the setitem call, whether it is the default stream, the legacy default
        # stream, the per-thread default stream, or another stream.
        streams = (cuda.stream(), cuda.default_stream(),
                   cuda.legacy_default_stream(),
                   cuda.per_thread_default_stream())

        for stream in streams:
            darr = cuda.to_device(np.arange(5))

            with patch.object(cuda.cudadrv.driver.Stream, 'synchronize',
                              return_value=None) as mock_sync:
                darr.setitem(0, 10, stream=stream)

            mock_sync.assert_not_called()

    @unittest.skip('Requires PR #6367')
    def test_issue_6505(self):
        # On Windows, the writes to ary_v would not be visible prior to the
        # assertion, due to the assignment being done with a kernel launch that
        # returns asynchronously - there should now be a sync after the kernel
        # launch to ensure that the writes are always visible.
        ary = cuda.mapped_array(2, dtype=np.int32)
        ary[:] = 0

        ary_v = ary.view('u1')
        ary_v[1] = 1
        ary_v[5] = 1
        self.assertEqual(sum(ary), 512)


if __name__ == '__main__':
    unittest.main()
