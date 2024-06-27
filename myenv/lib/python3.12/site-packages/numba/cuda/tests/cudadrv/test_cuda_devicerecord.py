import numpy as np
import ctypes
from numba.cuda.cudadrv.devicearray import (DeviceRecord, from_record_like,
                                            auto_device)
from numba.cuda.testing import unittest, CUDATestCase
from numba.cuda.testing import skip_on_cudasim
from numba.np import numpy_support
from numba import cuda

N_CHARS = 5

recordtype = np.dtype(
    [
        ('a', np.float64),
        ('b', np.int32),
        ('c', np.complex64),
        ('d', (np.str_, N_CHARS))
    ],
    align=True
)

recordwitharray = np.dtype(
    [
        ('g', np.int32),
        ('h', np.float32, 2)
    ],
    align=True
)

recwithmat = np.dtype([('i', np.int32),
                       ('j', np.float32, (3, 3))])

recwithrecwithmat = np.dtype([('x', np.int32), ('y', recwithmat)])


@skip_on_cudasim('Device Record API unsupported in the simulator')
class TestCudaDeviceRecord(CUDATestCase):
    """
    Tests the DeviceRecord class with np.void host types.
    """
    def setUp(self):
        super().setUp()
        self._create_data(np.zeros)

    def _create_data(self, array_ctor):
        self.dtype = np.dtype([('a', np.int32), ('b', np.float32)], align=True)
        self.hostz = array_ctor(1, self.dtype)[0]
        self.hostnz = array_ctor(1, self.dtype)[0]
        self.hostnz['a'] = 10
        self.hostnz['b'] = 11.0

    def _check_device_record(self, reference, rec):
        self.assertEqual(rec.shape, tuple())
        self.assertEqual(rec.strides, tuple())
        self.assertEqual(rec.dtype, reference.dtype)
        self.assertEqual(rec.alloc_size, reference.dtype.itemsize)
        self.assertIsNotNone(rec.gpu_data)
        self.assertNotEqual(rec.device_ctypes_pointer, ctypes.c_void_p(0))

        numba_type = numpy_support.from_dtype(reference.dtype)
        self.assertEqual(rec._numba_type_, numba_type)

    def test_device_record_interface(self):
        hostrec = self.hostz.copy()
        devrec = DeviceRecord(self.dtype)
        self._check_device_record(hostrec, devrec)

    def test_device_record_copy(self):
        hostrec = self.hostz.copy()
        devrec = DeviceRecord(self.dtype)
        devrec.copy_to_device(hostrec)

        # Copy back and check values are all zeros
        hostrec2 = self.hostnz.copy()
        devrec.copy_to_host(hostrec2)
        np.testing.assert_equal(self.hostz, hostrec2)

        # Copy non-zero values to GPU and back and check values
        hostrec3 = self.hostnz.copy()
        devrec.copy_to_device(hostrec3)

        hostrec4 = self.hostz.copy()
        devrec.copy_to_host(hostrec4)
        np.testing.assert_equal(hostrec4, self.hostnz)

    def test_from_record_like(self):
        # Create record from host record
        hostrec = self.hostz.copy()
        devrec = from_record_like(hostrec)
        self._check_device_record(hostrec, devrec)

        # Create record from device record and check for distinct data
        devrec2 = from_record_like(devrec)
        self._check_device_record(devrec, devrec2)
        self.assertNotEqual(devrec.gpu_data, devrec2.gpu_data)

    def test_auto_device(self):
        # Create record from host record
        hostrec = self.hostnz.copy()
        devrec, new_gpu_obj = auto_device(hostrec)
        self._check_device_record(hostrec, devrec)
        self.assertTrue(new_gpu_obj)

        # Copy data back and check it is equal to auto_device arg
        hostrec2 = self.hostz.copy()
        devrec.copy_to_host(hostrec2)
        np.testing.assert_equal(hostrec2, hostrec)


class TestCudaDeviceRecordWithRecord(TestCudaDeviceRecord):
    """
    Tests the DeviceRecord class with np.record host types
    """
    def setUp(self):
        CUDATestCase.setUp(self)
        self._create_data(np.recarray)


@skip_on_cudasim('Structured array attr access not supported in simulator')
class TestRecordDtypeWithStructArrays(CUDATestCase):
    '''
    Test operation of device arrays on structured arrays.
    '''

    def _createSampleArrays(self):
        self.sample1d = cuda.device_array(3, dtype=recordtype)
        self.samplerec1darr = cuda.device_array(1, dtype=recordwitharray)[0]
        self.samplerecmat = cuda.device_array(1,dtype=recwithmat)[0]

    def setUp(self):
        super().setUp()
        self._createSampleArrays()

        ary = self.sample1d
        for i in range(ary.size):
            x = i + 1
            ary[i]['a'] = x / 2
            ary[i]['b'] = x
            ary[i]['c'] = x * 1j
            ary[i]['d'] = str(x) * N_CHARS

    def test_structured_array1(self):
        ary = self.sample1d
        for i in range(self.sample1d.size):
            x = i + 1
            self.assertEqual(ary[i]['a'], x / 2)
            self.assertEqual(ary[i]['b'], x)
            self.assertEqual(ary[i]['c'], x * 1j)
            self.assertEqual(ary[i]['d'], str(x) * N_CHARS)

    def test_structured_array2(self):
        ary = self.samplerec1darr
        ary['g'] = 2
        ary['h'][0] = 3.0
        ary['h'][1] = 4.0
        self.assertEqual(ary['g'], 2)
        self.assertEqual(ary['h'][0], 3.0)
        self.assertEqual(ary['h'][1], 4.0)

    def test_structured_array3(self):
        ary = self.samplerecmat
        mat = np.array([[5.0, 10.0, 15.0],
                       [20.0, 25.0, 30.0],
                       [35.0, 40.0, 45.0]],
                       dtype=np.float32).reshape(3,3)
        ary['j'][:] = mat
        np.testing.assert_equal(ary['j'], mat)

    def test_structured_array4(self):
        arr = np.zeros(1, dtype=recwithrecwithmat)
        d_arr = cuda.to_device(arr)
        d_arr[0]['y']['i'] = 1
        self.assertEqual(d_arr[0]['y']['i'], 1)
        d_arr[0]['y']['j'][0, 0] = 2.0
        self.assertEqual(d_arr[0]['y']['j'][0, 0], 2.0)


if __name__ == '__main__':
    unittest.main()
