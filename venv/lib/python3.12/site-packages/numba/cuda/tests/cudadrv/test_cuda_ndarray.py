import itertools
import numpy as np
from numba.cuda.cudadrv import devicearray
from numba import cuda
from numba.cuda.testing import unittest, CUDATestCase
from numba.cuda.testing import skip_on_cudasim


class TestCudaNDArray(CUDATestCase):
    def test_device_array_interface(self):
        dary = cuda.device_array(shape=100)
        devicearray.verify_cuda_ndarray_interface(dary)

        ary = np.empty(100)
        dary = cuda.to_device(ary)
        devicearray.verify_cuda_ndarray_interface(dary)

        ary = np.asarray(1.234)
        dary = cuda.to_device(ary)
        self.assertEqual(dary.ndim, 0)
        devicearray.verify_cuda_ndarray_interface(dary)

    def test_device_array_from_readonly(self):
        ary = np.arange(100, dtype=np.float32)
        # Make the array readonly
        ary.flags.writeable = False
        self.assertFalse(ary.flags.writeable)
        # Ensure that we can copy the readonly array
        dary = cuda.to_device(ary)
        retr = dary.copy_to_host()
        np.testing.assert_array_equal(retr, ary)

    def test_devicearray_dtype(self):
        dary = cuda.device_array(shape=(100,), dtype="f4")
        self.assertEqual(dary.dtype, np.dtype("f4"))

    def test_devicearray_no_copy(self):
        array = np.arange(100, dtype=np.float32)
        cuda.to_device(array, copy=False)

    def test_devicearray_shape(self):
        ary = np.arange(2 * 3 * 4).reshape(2, 3, 4)
        dary = cuda.to_device(ary)
        self.assertEqual(ary.shape, dary.shape)
        self.assertEqual(ary.shape[1:], dary.shape[1:])

    def test_devicearray(self):
        array = np.arange(100, dtype=np.int32)
        original = array.copy()
        gpumem = cuda.to_device(array)
        array[:] = 0
        gpumem.copy_to_host(array)

        np.testing.assert_array_equal(array, original)

    def test_stream_bind(self):
        stream = cuda.stream()
        with stream.auto_synchronize():
            arr = cuda.device_array(
                (3, 3),
                dtype=np.float64,
                stream=stream)
            self.assertEqual(arr.bind(stream).stream, stream)
            self.assertEqual(arr.stream, stream)

    def test_len_1d(self):
        ary = np.empty((3,))
        dary = cuda.device_array(3)
        self.assertEqual(len(ary), len(dary))

    def test_len_2d(self):
        ary = np.empty((3, 5))
        dary = cuda.device_array((3, 5))
        self.assertEqual(len(ary), len(dary))

    def test_len_3d(self):
        ary = np.empty((3, 5, 7))
        dary = cuda.device_array((3, 5, 7))
        self.assertEqual(len(ary), len(dary))

    def test_devicearray_partition(self):
        N = 100
        array = np.arange(N, dtype=np.int32)
        original = array.copy()
        gpumem = cuda.to_device(array)
        left, right = gpumem.split(N // 2)

        array[:] = 0

        self.assertTrue(np.all(array == 0))

        right.copy_to_host(array[N // 2:])
        left.copy_to_host(array[:N // 2])

        self.assertTrue(np.all(array == original))

    def test_devicearray_replace(self):
        N = 100
        array = np.arange(N, dtype=np.int32)
        original = array.copy()
        gpumem = cuda.to_device(array)
        cuda.to_device(array * 2, to=gpumem)
        gpumem.copy_to_host(array)
        np.testing.assert_array_equal(array, original * 2)

    @skip_on_cudasim('This works in the simulator')
    def test_devicearray_transpose_wrongdim(self):
        gpumem = cuda.to_device(np.array(np.arange(12)).reshape(3, 4, 1))

        with self.assertRaises(NotImplementedError) as e:
            np.transpose(gpumem)

        self.assertEqual(
            "transposing a non-2D DeviceNDArray isn't supported",
            str(e.exception))

    def test_devicearray_transpose_identity(self):
        # any-shape identities should work
        original = np.array(np.arange(24)).reshape(3, 4, 2)
        array = np.transpose(cuda.to_device(original),
                             axes=(0, 1, 2)).copy_to_host()
        self.assertTrue(np.all(array == original))

    def test_devicearray_transpose_duplicatedaxis(self):
        gpumem = cuda.to_device(np.array(np.arange(12)).reshape(3, 4))

        with self.assertRaises(ValueError) as e:
            np.transpose(gpumem, axes=(0, 0))

        self.assertIn(
            str(e.exception),
            container=[
                'invalid axes list (0, 0)',  # GPU
                'repeated axis in transpose',  # sim
            ])

    def test_devicearray_transpose_wrongaxis(self):
        gpumem = cuda.to_device(np.array(np.arange(12)).reshape(3, 4))

        with self.assertRaises(ValueError) as e:
            np.transpose(gpumem, axes=(0, 2))

        self.assertIn(
            str(e.exception),
            container=[
                'invalid axes list (0, 2)',  # GPU
                'invalid axis for this array',
                'axis 2 is out of bounds for array of dimension 2',  # sim
            ])

    def test_devicearray_view_ok(self):
        original = np.array(np.arange(12), dtype="i2").reshape(3, 4)
        array = cuda.to_device(original)
        for dtype in ("i4", "u4", "i8", "f8"):
            with self.subTest(dtype=dtype):
                np.testing.assert_array_equal(
                    array.view(dtype).copy_to_host(),
                    original.view(dtype)
                )

    def test_devicearray_view_ok_not_c_contig(self):
        original = np.array(np.arange(32), dtype="i2").reshape(4, 8)
        array = cuda.to_device(original)[:, ::2]
        original = original[:, ::2]
        np.testing.assert_array_equal(
            array.view("u2").copy_to_host(),
            original.view("u2")
        )

    def test_devicearray_view_bad_not_c_contig(self):
        original = np.array(np.arange(32), dtype="i2").reshape(4, 8)
        array = cuda.to_device(original)[:, ::2]
        with self.assertRaises(ValueError) as e:
            array.view("i4")

        msg = str(e.exception)
        self.assertIn('To change to a dtype of a different size,', msg)

        contiguous_pre_np123 = 'the array must be C-contiguous' in msg
        contiguous_post_np123 = 'the last axis must be contiguous' in msg
        self.assertTrue(contiguous_pre_np123 or contiguous_post_np123,
                        'Expected message to mention contiguity')

    def test_devicearray_view_bad_itemsize(self):
        original = np.array(np.arange(12), dtype="i2").reshape(4, 3)
        array = cuda.to_device(original)
        with self.assertRaises(ValueError) as e:
            array.view("i4")
        self.assertEqual(
            "When changing to a larger dtype,"
            " its size must be a divisor of the total size in bytes"
            " of the last axis of the array.",
            str(e.exception))

    def test_devicearray_transpose_ok(self):
        original = np.array(np.arange(12)).reshape(3, 4)
        array = np.transpose(cuda.to_device(original)).copy_to_host()
        self.assertTrue(np.all(array == original.T))

    def test_devicearray_transpose_T(self):
        original = np.array(np.arange(12)).reshape(3, 4)
        array = cuda.to_device(original).T.copy_to_host()
        self.assertTrue(np.all(array == original.T))

    def test_devicearray_contiguous_slice(self):
        # memcpys are dumb ranges of bytes, so trying to
        # copy to a non-contiguous range shouldn't work!
        a = np.arange(25).reshape(5, 5, order='F')
        s = np.full(fill_value=5, shape=(5,))

        d = cuda.to_device(a)
        a[2] = s

        # d is in F-order (not C-order), so d[2] is not contiguous
        # (40-byte strides). This means we can't memcpy to it!
        with self.assertRaises(ValueError) as e:
            d[2].copy_to_device(s)
        self.assertEqual(
            devicearray.errmsg_contiguous_buffer,
            str(e.exception))

        # if d[2].copy_to_device(s), then this would pass:
        # self.assertTrue((a == d.copy_to_host()).all())

    def _test_devicearray_contiguous_host_copy(self, a_c, a_f):
        """
        Checks host->device memcpys
        """
        self.assertTrue(a_c.flags.c_contiguous)
        self.assertTrue(a_f.flags.f_contiguous)

        for original, copy in [
            (a_f, a_f),
            (a_f, a_c),
            (a_c, a_f),
            (a_c, a_c),
        ]:
            msg = '%s => %s' % (
                'C' if original.flags.c_contiguous else 'F',
                'C' if copy.flags.c_contiguous else 'F',
            )

            d = cuda.to_device(original)
            d.copy_to_device(copy)
            self.assertTrue(np.all(d.copy_to_host() == a_c), msg=msg)
            self.assertTrue(np.all(d.copy_to_host() == a_f), msg=msg)

    def test_devicearray_contiguous_copy_host_3d(self):
        a_c = np.arange(5 * 5 * 5).reshape(5, 5, 5)
        a_f = np.array(a_c, order='F')
        self._test_devicearray_contiguous_host_copy(a_c, a_f)

    def test_devicearray_contiguous_copy_host_1d(self):
        a_c = np.arange(5)
        a_f = np.array(a_c, order='F')
        self._test_devicearray_contiguous_host_copy(a_c, a_f)

    def test_devicearray_contiguous_copy_device(self):
        a_c = np.arange(5 * 5 * 5).reshape(5, 5, 5)
        a_f = np.array(a_c, order='F')
        self.assertTrue(a_c.flags.c_contiguous)
        self.assertTrue(a_f.flags.f_contiguous)

        d = cuda.to_device(a_c)

        with self.assertRaises(ValueError) as e:
            d.copy_to_device(cuda.to_device(a_f))
        self.assertEqual(
            "incompatible strides: {} vs. {}".format(a_c.strides, a_f.strides),
            str(e.exception))

        d.copy_to_device(cuda.to_device(a_c))
        self.assertTrue(np.all(d.copy_to_host() == a_c))

        d = cuda.to_device(a_f)

        with self.assertRaises(ValueError) as e:
            d.copy_to_device(cuda.to_device(a_c))
        self.assertEqual(
            "incompatible strides: {} vs. {}".format(a_f.strides, a_c.strides),
            str(e.exception))

        d.copy_to_device(cuda.to_device(a_f))
        self.assertTrue(np.all(d.copy_to_host() == a_f))

    def test_devicearray_broadcast_host_copy(self):
        broadsize = 4
        coreshape = (2, 3)
        coresize = np.prod(coreshape)
        core_c = np.arange(coresize).reshape(coreshape, order='C')
        core_f = np.arange(coresize).reshape(coreshape, order='F')
        for dim in range(len(coreshape)):
            newindex = (slice(None),) * dim + (np.newaxis,)
            broadshape = coreshape[:dim] + (broadsize,) + coreshape[dim:]
            broad_c = np.broadcast_to(core_c[newindex], broadshape)
            broad_f = np.broadcast_to(core_f[newindex], broadshape)
            dbroad_c = cuda.to_device(broad_c)
            dbroad_f = cuda.to_device(broad_f)
            np.testing.assert_array_equal(dbroad_c.copy_to_host(), broad_c)
            np.testing.assert_array_equal(dbroad_f.copy_to_host(), broad_f)
            # Also test copying across different core orderings
            dbroad_c.copy_to_device(broad_f)
            dbroad_f.copy_to_device(broad_c)
            np.testing.assert_array_equal(dbroad_c.copy_to_host(), broad_f)
            np.testing.assert_array_equal(dbroad_f.copy_to_host(), broad_c)

    def test_devicearray_contiguous_host_strided(self):
        a_c = np.arange(10)
        d = cuda.to_device(a_c)
        arr = np.arange(20)[::2]
        d.copy_to_device(arr)
        np.testing.assert_array_equal(d.copy_to_host(), arr)

    def test_devicearray_contiguous_device_strided(self):
        d = cuda.to_device(np.arange(20))
        arr = np.arange(20)

        with self.assertRaises(ValueError) as e:
            d.copy_to_device(cuda.to_device(arr)[::2])
        self.assertEqual(
            devicearray.errmsg_contiguous_buffer,
            str(e.exception))

    @skip_on_cudasim('DeviceNDArray class not present in simulator')
    def test_devicearray_relaxed_strides(self):
        # From the reproducer in Issue #6824.

        # Construct a device array that is contiguous even though
        # the strides for the first axis (800) are not equal to
        # the strides * size (10 * 8 = 80) for the previous axis,
        # because the first axis size is 1.
        arr = devicearray.DeviceNDArray((1, 10), (800, 8), np.float64)

        # Ensure we still believe the array to be contiguous because
        # strides checking is relaxed.
        self.assertTrue(arr.flags['C_CONTIGUOUS'])
        self.assertTrue(arr.flags['F_CONTIGUOUS'])

    def test_c_f_contiguity_matches_numpy(self):
        # From the reproducer in Issue #4943.

        shapes = ((1, 4), (4, 1))
        orders = ('C', 'F')

        for shape, order in itertools.product(shapes, orders):
            arr = np.ndarray(shape, order=order)
            d_arr = cuda.to_device(arr)
            self.assertEqual(arr.flags['C_CONTIGUOUS'],
                             d_arr.flags['C_CONTIGUOUS'])
            self.assertEqual(arr.flags['F_CONTIGUOUS'],
                             d_arr.flags['F_CONTIGUOUS'])

    @skip_on_cudasim('Typing not done in the simulator')
    def test_devicearray_typing_order_simple_c(self):
        # C-order 1D array
        a = np.zeros(10, order='C')
        d = cuda.to_device(a)
        self.assertEqual(d._numba_type_.layout, 'C')

    @skip_on_cudasim('Typing not done in the simulator')
    def test_devicearray_typing_order_simple_f(self):
        # F-order array that is also C layout.
        a = np.zeros(10, order='F')
        d = cuda.to_device(a)
        self.assertEqual(d._numba_type_.layout, 'C')

    @skip_on_cudasim('Typing not done in the simulator')
    def test_devicearray_typing_order_2d_c(self):
        # C-order 2D array
        a = np.zeros((2, 10), order='C')
        d = cuda.to_device(a)
        self.assertEqual(d._numba_type_.layout, 'C')

    @skip_on_cudasim('Typing not done in the simulator')
    def test_devicearray_typing_order_2d_f(self):
        # F-order array that can only be F layout
        a = np.zeros((2, 10), order='F')
        d = cuda.to_device(a)
        self.assertEqual(d._numba_type_.layout, 'F')

    @skip_on_cudasim('Typing not done in the simulator')
    def test_devicearray_typing_order_noncontig_slice_c(self):
        # Non-contiguous slice of C-order array
        a = np.zeros((5, 5), order='C')
        d = cuda.to_device(a)[:,2]
        self.assertEqual(d._numba_type_.layout, 'A')

    @skip_on_cudasim('Typing not done in the simulator')
    def test_devicearray_typing_order_noncontig_slice_f(self):
        # Non-contiguous slice of F-order array
        a = np.zeros((5, 5), order='F')
        d = cuda.to_device(a)[2,:]
        self.assertEqual(d._numba_type_.layout, 'A')

    @skip_on_cudasim('Typing not done in the simulator')
    def test_devicearray_typing_order_contig_slice_c(self):
        # Contiguous slice of C-order array
        a = np.zeros((5, 5), order='C')
        d = cuda.to_device(a)[2,:]
        self.assertEqual(d._numba_type_.layout, 'C')

    @skip_on_cudasim('Typing not done in the simulator')
    def test_devicearray_typing_order_contig_slice_f(self):
        # Contiguous slice of F-order array - is both C- and F-contiguous, so
        # types as 'C' layout
        a = np.zeros((5, 5), order='F')
        d = cuda.to_device(a)[:,2]
        self.assertEqual(d._numba_type_.layout, 'C')

    @skip_on_cudasim('Typing not done in the simulator')
    def test_devicearray_typing_order_broadcasted(self):
        # Broadcasted array, similar to that used for passing scalars to ufuncs
        a = np.broadcast_to(np.array([1]), (10,))
        d = cuda.to_device(a)
        self.assertEqual(d._numba_type_.layout, 'A')

    def test_bug6697(self):
        ary = np.arange(10, dtype=np.int16)
        dary = cuda.to_device(ary)
        got = np.asarray(dary)
        self.assertEqual(got.dtype, dary.dtype)

    @skip_on_cudasim('DeviceNDArray class not present in simulator')
    def test_issue_8477(self):
        # Ensure that we can copy a zero-length device array to a zero-length
        # host array when the strides of the device and host arrays differ -
        # this should be possible because the strides are irrelevant when the
        # length is zero. For more info see
        # https://github.com/numba/numba/issues/8477.

        # Create a device array with shape (0,) and strides (8,)
        dev_array = devicearray.DeviceNDArray(shape=(0,), strides=(8,),
                                              dtype=np.int8)

        # Create a host array with shape (0,) and strides (0,)
        host_array = np.ndarray(shape=(0,), strides=(0,), dtype=np.int8)

        # Sanity check for this test - ensure our destination has the strides
        # we expect, because strides can be ignored in some cases by the
        # ndarray constructor - checking here ensures that we haven't failed to
        # account for unexpected behaviour across different versions of NumPy
        self.assertEqual(host_array.strides, (0,))

        # Ensure that the copy succeeds in both directions
        dev_array.copy_to_host(host_array)
        dev_array.copy_to_device(host_array)

        # Ensure that a device-to-device copy also succeeds when the strides
        # differ - one way of doing this is to copy the host array across and
        # use that for copies in both directions.
        dev_array_from_host = cuda.to_device(host_array)
        self.assertEqual(dev_array_from_host.shape, (0,))
        self.assertEqual(dev_array_from_host.strides, (0,))

        dev_array.copy_to_device(dev_array_from_host)
        dev_array_from_host.copy_to_device(dev_array)


class TestRecarray(CUDATestCase):
    def test_recarray(self):
        # From issue #4111
        a = np.recarray((16,), dtype=[
            ("value1", np.int64),
            ("value2", np.float64),
        ])
        a.value1 = np.arange(a.size, dtype=np.int64)
        a.value2 = np.arange(a.size, dtype=np.float64) / 100

        expect1 = a.value1
        expect2 = a.value2

        def test(x, out1, out2):
            i = cuda.grid(1)
            if i < x.size:
                out1[i] = x.value1[i]
                out2[i] = x.value2[i]

        got1 = np.zeros_like(expect1)
        got2 = np.zeros_like(expect2)
        cuda.jit(test)[1, a.size](a, got1, got2)

        np.testing.assert_array_equal(expect1, got1)
        np.testing.assert_array_equal(expect2, got2)


class TestCoreContiguous(CUDATestCase):
    def _test_against_array_core(self, view):
        self.assertEqual(
            devicearray.is_contiguous(view),
            devicearray.array_core(view).flags['C_CONTIGUOUS']
        )

    def test_device_array_like_1d(self):
        d_a = cuda.device_array(10, order='C')
        self._test_against_array_core(d_a)

    def test_device_array_like_2d(self):
        d_a = cuda.device_array((10, 12), order='C')
        self._test_against_array_core(d_a)

    def test_device_array_like_2d_transpose(self):
        d_a = cuda.device_array((10, 12), order='C')
        self._test_against_array_core(d_a.T)

    def test_device_array_like_3d(self):
        d_a = cuda.device_array((10, 12, 14), order='C')
        self._test_against_array_core(d_a)

    def test_device_array_like_1d_f(self):
        d_a = cuda.device_array(10, order='F')
        self._test_against_array_core(d_a)

    def test_device_array_like_2d_f(self):
        d_a = cuda.device_array((10, 12), order='F')
        self._test_against_array_core(d_a)

    def test_device_array_like_2d_f_transpose(self):
        d_a = cuda.device_array((10, 12), order='F')
        self._test_against_array_core(d_a.T)

    def test_device_array_like_3d_f(self):
        d_a = cuda.device_array((10, 12, 14), order='F')
        self._test_against_array_core(d_a)

    def test_1d_view(self):
        shape = 10
        view = np.zeros(shape)[::2]
        self._test_against_array_core(view)

    def test_1d_view_f(self):
        shape = 10
        view = np.zeros(shape, order='F')[::2]
        self._test_against_array_core(view)

    def test_2d_view(self):
        shape = (10, 12)
        view = np.zeros(shape)[::2, ::2]
        self._test_against_array_core(view)

    def test_2d_view_f(self):
        shape = (10, 12)
        view = np.zeros(shape, order='F')[::2, ::2]
        self._test_against_array_core(view)


if __name__ == '__main__':
    unittest.main()
