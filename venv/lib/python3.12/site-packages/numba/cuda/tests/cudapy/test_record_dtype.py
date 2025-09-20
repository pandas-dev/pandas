import numpy as np
from numba import cuda
from numba.core import types
from numba.cuda.testing import skip_on_cudasim, CUDATestCase
import unittest
from numba.np import numpy_support


def set_a(ary, i, v):
    ary[i].a = v


def set_b(ary, i, v):
    ary[i].b = v


def set_c(ary, i, v):
    ary[i].c = v


def set_record(ary, i, j):
    ary[i] = ary[j]


def record_set_a(r, v):
    r.a = v


def record_set_b(r, v):
    r.b = v


def record_set_c(r, v):
    r.c = v


def record_read_a(r, arr):
    arr[0] = r.a


def record_read_b(r, arr):
    arr[0] = r.b


def record_read_c(r, arr):
    arr[0] = r.c


def record_write_array(r):
    r.g = 2
    r.h[0] = 3.0
    r.h[1] = 4.0


def record_write_2d_array(r):
    r.i = 3
    r.j[0, 0] = 5.0
    r.j[0, 1] = 6.0
    r.j[1, 0] = 7.0
    r.j[1, 1] = 8.0
    r.j[2, 0] = 9.0
    r.j[2, 1] = 10.0


def record_read_array(r, a):
    a[0] = r.h[0]
    a[1] = r.h[1]


def record_read_2d_array(r, a):
    a[0, 0] = r.j[0, 0]
    a[0, 1] = r.j[0, 1]
    a[1, 0] = r.j[1, 0]
    a[1, 1] = r.j[1, 1]
    a[2, 0] = r.j[2, 0]
    a[2, 1] = r.j[2, 1]


recordtype = np.dtype(
    [
        ('a', np.float64),
        ('b', np.int32),
        ('c', np.complex64),
        ('d', (np.uint8, 5))
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

recordwith2darray = np.dtype([('i', np.int32),
                              ('j', np.float32, (3, 2))])

nested_array1_dtype = np.dtype([("array1", np.int16, (3,))], align=True)

nested_array2_dtype = np.dtype([("array2", np.int16, (3, 2))], align=True)


# Functions used for "full array" tests

def record_write_full_array(rec):
    rec.j[:, :] = np.ones((3, 2))


def record_write_full_array_alt(rec):
    rec['j'][:, :] = np.ones((3, 2))


def recarray_set_record(ary, rec):
    ary[0] = rec


def recarray_write_array_of_nestedarray_broadcast(ary):
    ary.j[:, :, :] = 1
    return ary


def record_setitem_array(rec_source, rec_dest):
    rec_dest['j'] = rec_source['j']


def recarray_write_array_of_nestedarray(ary):
    ary.j[:, :, :] = np.ones((2, 3, 2))
    return ary


def recarray_getitem_return(ary):
    return ary[0]


def recarray_getitem_field_return(ary):
    return ary['h']


def recarray_getitem_field_return2(ary):
    return ary.h


def recarray_getitem_field_return2_2d(ary):
    return ary.j


def record_read_array0(ary):
    return ary.h[0]


def record_read_array1(ary):
    return ary.h[1]


def record_read_whole_array(ary):
    return ary.h


def record_read_2d_array00(ary):
    return ary.j[0, 0]


def record_read_2d_array10(ary):
    return ary.j[1, 0]


def record_read_2d_array01(ary):
    return ary.j[0, 1]


def assign_array_to_nested(dest, src):
    dest['array1'] = src


def assign_array_to_nested_2d(dest, src):
    dest['array2'] = src


class TestRecordDtype(CUDATestCase):

    def _createSampleArrays(self):
        self.sample1d = np.recarray(3, dtype=recordtype)
        self.samplerec1darr = np.recarray(1, dtype=recordwitharray)[0]
        self.samplerec2darr = np.recarray(1, dtype=recordwith2darray)[0]

    def setUp(self):
        super().setUp()
        self._createSampleArrays()

        ary = self.sample1d
        for i in range(ary.size):
            x = i + 1
            ary[i]['a'] = x / 2
            ary[i]['b'] = x
            ary[i]['c'] = x * 1j
            ary[i]['d'] = "%d" % x

    def get_cfunc(self, pyfunc, argspec):
        return cuda.jit()(pyfunc)

    def _test_set_equal(self, pyfunc, value, valuetype):
        rec = numpy_support.from_dtype(recordtype)
        cfunc = self.get_cfunc(pyfunc, (rec[:], types.intp, valuetype))

        for i in range(self.sample1d.size):
            got = self.sample1d.copy()

            # Force the argument to the pure Python function to be
            # a recarray, as attribute access isn't supported on
            # structured arrays.
            expect = got.copy().view(np.recarray)

            cfunc[1, 1](got, i, value)
            pyfunc(expect, i, value)

            # Match the entire array to ensure no memory corruption
            self.assertTrue(np.all(expect == got))

    def test_set_a(self):
        self._test_set_equal(set_a, 3.1415, types.float64)
        # Test again to check if coercion works
        self._test_set_equal(set_a, 3., types.float32)

    def test_set_b(self):
        self._test_set_equal(set_b, 123, types.int32)
        # Test again to check if coercion works
        self._test_set_equal(set_b, 123, types.float64)

    def test_set_c(self):
        self._test_set_equal(set_c, 43j, types.complex64)
        # Test again to check if coercion works
        self._test_set_equal(set_c, 43j, types.complex128)

    def test_set_record(self):
        pyfunc = set_record
        rec = numpy_support.from_dtype(recordtype)
        cfunc = self.get_cfunc(pyfunc, (rec[:], types.intp, types.intp))

        test_indices = [(0, 1), (1, 2), (0, 2)]
        for i, j in test_indices:
            expect = self.sample1d.copy()
            pyfunc(expect, i, j)

            got = self.sample1d.copy()
            cfunc[1, 1](got, i, j)

            # Match the entire array to ensure no memory corruption
            self.assertEqual(expect[i], expect[j])
            self.assertEqual(got[i], got[j])
            self.assertTrue(np.all(expect == got))

    def _test_rec_set(self, v, pyfunc, f):
        rec = self.sample1d.copy()[0]
        nbrecord = numpy_support.from_dtype(recordtype)
        cfunc = self.get_cfunc(pyfunc, (nbrecord,))
        cfunc[1, 1](rec, v)
        np.testing.assert_equal(rec[f], v)

    def test_rec_set_a(self):
        self._test_rec_set(np.float64(1.5), record_set_a, 'a')

    def test_rec_set_b(self):
        self._test_rec_set(np.int32(2), record_set_b, 'b')

    def test_rec_set_c(self):
        self._test_rec_set(np.complex64(4.0 + 5.0j), record_set_c, 'c')

    def _test_rec_read(self, v, pyfunc, f):
        rec = self.sample1d.copy()[0]
        rec[f] = v
        arr = np.zeros(1, v.dtype)
        nbrecord = numpy_support.from_dtype(recordtype)
        cfunc = self.get_cfunc(pyfunc, (nbrecord,))
        cfunc[1, 1](rec, arr)
        np.testing.assert_equal(arr[0], v)

    def test_rec_read_a(self):
        self._test_rec_read(np.float64(1.5), record_read_a, 'a')

    def test_rec_read_b(self):
        self._test_rec_read(np.int32(2), record_read_b, 'b')

    def test_rec_read_c(self):
        self._test_rec_read(np.complex64(4.0 + 5.0j), record_read_c, 'c')

    def test_record_write_1d_array(self):
        '''
        Test writing to a 1D array within a structured type
        '''
        rec = self.samplerec1darr.copy()
        nbrecord = numpy_support.from_dtype(recordwitharray)
        cfunc = self.get_cfunc(record_write_array, (nbrecord,))

        cfunc[1, 1](rec)
        expected = self.samplerec1darr.copy()
        expected['g'] = 2
        expected['h'][0] = 3.0
        expected['h'][1] = 4.0

        np.testing.assert_equal(expected, rec)

    def test_record_write_2d_array(self):
        '''
        Test writing to a 2D array within a structured type
        '''
        rec = self.samplerec2darr.copy()
        nbrecord = numpy_support.from_dtype(recordwith2darray)
        cfunc = self.get_cfunc(record_write_2d_array, (nbrecord,))
        cfunc[1, 1](rec)

        expected = self.samplerec2darr.copy()
        expected['i'] = 3
        expected['j'][:] = np.asarray([5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
                                      np.float32).reshape(3, 2)
        np.testing.assert_equal(expected, rec)

    def test_record_read_1d_array(self):
        '''
        Test reading from a 1D array within a structured type
        '''
        rec = self.samplerec1darr.copy()
        rec['h'][0] = 4.0
        rec['h'][1] = 5.0

        nbrecord = numpy_support.from_dtype(recordwitharray)
        cfunc = self.get_cfunc(record_read_array, (nbrecord,))
        arr = np.zeros(2, dtype=rec['h'].dtype)
        cfunc[1, 1](rec, arr)

        np.testing.assert_equal(rec['h'], arr)

    def test_record_read_2d_array(self):
        '''
        Test reading from a 2D array within a structured type
        '''
        rec = self.samplerec2darr.copy()
        rec['j'][:] = np.asarray([5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
                                 np.float32).reshape(3, 2)

        nbrecord = numpy_support.from_dtype(recordwith2darray)
        cfunc = self.get_cfunc(record_read_2d_array, (nbrecord,))
        arr = np.zeros((3,2), dtype=rec['j'].dtype)
        cfunc[1, 1](rec, arr)

        np.testing.assert_equal(rec['j'], arr)


@skip_on_cudasim('Structured array attr access not supported in simulator')
class TestRecordDtypeWithStructArrays(TestRecordDtype):
    '''
    Same as TestRecordDtype, but using structured arrays instead of recarrays.
    '''

    def _createSampleArrays(self):
        self.sample1d = np.zeros(3, dtype=recordtype)
        self.samplerec1darr = np.zeros(1, dtype=recordwitharray)[0]
        self.samplerec2darr = np.zeros(1, dtype=recordwith2darray)[0]


class TestNestedArrays(CUDATestCase):

    # These tests mirror those from
    # numba.tests.test_record_dtype.TestNestedArrays added in PR
    # #7359: https://github.com/numba/numba/pull/7359

    # The code cannot be shared between the two classes without modification,
    # as the CUDA test implementations need to be launched (and in some cases
    # wrapped in an outer function to handle the return value). Otherwise, the
    # code here is kept as similar to that in the equivalent CPU tests as
    # possible.

    # Reading records / recarrays

    def get_cfunc(self, pyfunc, retty):
        # Create a host-callable function for testing CUDA device functions
        # that get a value from a record array
        inner = cuda.jit(device=True)(pyfunc)

        @cuda.jit
        def outer(arg0, res):
            res[0] = inner(arg0)

        def host(arg0):
            res = np.zeros(1, dtype=retty)
            outer[1, 1](arg0, res)
            return res[0]

        return host

    def test_record_read_array(self):
        # Test reading from a 1D array within a structured type
        nbval = np.recarray(1, dtype=recordwitharray)
        nbval[0].h[0] = 15.0
        nbval[0].h[1] = 25.0
        cfunc = self.get_cfunc(record_read_array0, np.float32)
        res = cfunc(nbval[0])
        np.testing.assert_equal(res, nbval[0].h[0])

        cfunc = self.get_cfunc(record_read_array1, np.float32)
        res = cfunc(nbval[0])
        np.testing.assert_equal(res, nbval[0].h[1])

    def test_record_read_2d_array(self):
        # Test reading from a 2D array within a structured type
        nbval = np.recarray(1, dtype=recordwith2darray)
        nbval[0].j = np.asarray([1.5, 2.5, 3.5, 4.5, 5.5, 6.5],
                                np.float32).reshape(3, 2)
        cfunc = self.get_cfunc(record_read_2d_array00, np.float32)
        res = cfunc(nbval[0])
        np.testing.assert_equal(res, nbval[0].j[0, 0])

        cfunc = self.get_cfunc(record_read_2d_array01, np.float32)
        res = cfunc(nbval[0])
        np.testing.assert_equal(res, nbval[0].j[0, 1])

        cfunc = self.get_cfunc(record_read_2d_array10, np.float32)
        res = cfunc(nbval[0])
        np.testing.assert_equal(res, nbval[0].j[1, 0])

    def test_setitem(self):
        def gen():
            nbarr1 = np.recarray(1, dtype=recordwith2darray)
            nbarr1[0] = np.array([(1, ((1, 2), (4, 5), (2, 3)))],
                                 dtype=recordwith2darray)[0]
            nbarr2 = np.recarray(1, dtype=recordwith2darray)
            nbarr2[0] = np.array([(10, ((10, 20), (40, 50), (20, 30)))],
                                 dtype=recordwith2darray)[0]
            return nbarr1[0], nbarr2[0]
        pyfunc = record_setitem_array
        pyargs = gen()
        pyfunc(*pyargs)

        cfunc = cuda.jit(pyfunc)
        cuargs = gen()
        cfunc[1, 1](*cuargs)
        np.testing.assert_equal(pyargs, cuargs)

    def test_getitem_idx(self):
        # Test __getitem__ with numerical index

        # This tests returning a record when passing an array and
        # returning the first item when passing a record
        nbarr = np.recarray(2, dtype=recordwitharray)
        nbarr[0] = np.array([(1, (2, 3))], dtype=recordwitharray)[0]
        for arg, retty in [(nbarr, recordwitharray), (nbarr[0], np.int32)]:
            pyfunc = recarray_getitem_return
            arr_expected = pyfunc(arg)
            cfunc = self.get_cfunc(pyfunc, retty)
            arr_res = cfunc(arg)
            np.testing.assert_equal(arr_res, arr_expected)

    # Writing to records / recarrays

    @skip_on_cudasim('Structured array attr access not supported in simulator')
    def test_set_record(self):
        # Test setting an entire record
        rec = np.ones(2, dtype=recordwith2darray).view(np.recarray)[0]
        nbarr = np.zeros(2, dtype=recordwith2darray).view(np.recarray)
        arr = np.zeros(2, dtype=recordwith2darray).view(np.recarray)
        pyfunc = recarray_set_record
        pyfunc(arr, rec)
        kernel = cuda.jit(pyfunc)
        kernel[1, 1](nbarr, rec)
        np.testing.assert_equal(nbarr, arr)

    def test_assign_array_to_nested(self):
        src = (np.arange(3) + 1).astype(np.int16)
        got = np.zeros(2, dtype=nested_array1_dtype)
        expected = np.zeros(2, dtype=nested_array1_dtype)

        pyfunc = assign_array_to_nested
        kernel = cuda.jit(pyfunc)

        kernel[1, 1](got[0], src)
        pyfunc(expected[0], src)

        np.testing.assert_array_equal(expected, got)

    def test_assign_array_to_nested_2d(self):
        src = (np.arange(6) + 1).astype(np.int16).reshape((3, 2))
        got = np.zeros(2, dtype=nested_array2_dtype)
        expected = np.zeros(2, dtype=nested_array2_dtype)

        pyfunc = assign_array_to_nested_2d
        kernel = cuda.jit(pyfunc)

        kernel[1, 1](got[0], src)
        pyfunc(expected[0], src)

        np.testing.assert_array_equal(expected, got)

    def test_issue_7693(self):
        src_dtype = np.dtype([
            ("user", np.float64),
            ("array", np.int16, (3,))],
            align=True)

        dest_dtype = np.dtype([
            ("user1", np.float64),
            ("array1", np.int16, (3,))],
            align=True)

        @cuda.jit
        def copy(index, src, dest):
            dest['user1'] = src[index]['user']
            dest['array1'] = src[index]['array']

        source = np.zeros(2, dtype=src_dtype)
        got = np.zeros(2, dtype=dest_dtype)
        expected = np.zeros(2, dtype=dest_dtype)

        source[0] = (1.2, [1, 2, 3])
        copy[1, 1](0, source, got[0])
        copy.py_func(0, source, expected[0])

        np.testing.assert_array_equal(expected, got)

    # Reading and returning arrays from recarrays - the following functions are
    # all xfailed because CUDA cannot handle returning arrays from device
    # functions (or creating arrays in general).

    @unittest.expectedFailure
    def test_getitem_idx_2darray(self):
        # Test __getitem__ with numerical index
        #
        # This test returning a record when passing an array and
        # return the first item when passing a record
        nbarr = np.recarray(2, dtype=recordwith2darray)
        nbarr[0] = np.array([(1, ((1,2),(4,5),(2,3)))],
                            dtype=recordwith2darray)[0]
        for arg, retty in [(nbarr, recordwith2darray),
                           (nbarr[0], (np.float32, (3, 2)))]:
            pyfunc = recarray_getitem_field_return2_2d
            arr_expected = pyfunc(arg)
            cfunc = self.get_cfunc(pyfunc, retty)
            arr_res = cfunc(arg)
            np.testing.assert_equal(arr_res, arr_expected)

    @unittest.expectedFailure
    def test_return_getattr_getitem_fieldname(self):
        # Test __getitem__ with field name and getattr .field_name
        #
        # This tests returning a array of nestedarrays when passing an array and
        # returning a nestedarray when passing a record
        nbarr = np.recarray(2, dtype=recordwitharray)
        nbarr[0] = np.array([(1, (2,3))], dtype=recordwitharray)[0]
        for arg, retty in [(nbarr, recordwitharray), (nbarr[0], np.float32)]:
            for pyfunc in [recarray_getitem_field_return,
                           recarray_getitem_field_return2]:
                arr_expected = pyfunc(arg)
                cfunc = self.get_cfunc(pyfunc, retty)
                arr_res = cfunc(arg)
                np.testing.assert_equal(arr_res, arr_expected)

    @unittest.expectedFailure
    def test_record_read_arrays(self):
        # Test reading from a 1D array within a structured type
        nbval = np.recarray(2, dtype=recordwitharray)
        nbval[0].h[0] = 15.0
        nbval[0].h[1] = 25.0
        nbval[1].h[0] = 35.0
        nbval[1].h[1] = 45.4
        cfunc = self.get_cfunc(record_read_whole_array, np.float32)
        res = cfunc(nbval)
        np.testing.assert_equal(res, nbval.h)

    @unittest.expectedFailure
    def test_return_array(self):
        # Test getitem record AND array within record and returning it
        nbval = np.recarray(2, dtype=recordwitharray)
        nbval[0] = np.array([(1, (2,3))], dtype=recordwitharray)[0]
        pyfunc = record_read_array0
        arr_expected = pyfunc(nbval)
        cfunc = self.get_cfunc(pyfunc, np.float32)
        arr_res = cfunc(nbval)
        np.testing.assert_equal(arr_expected, arr_res)

    @skip_on_cudasim('Will unexpectedly pass on cudasim')
    @unittest.expectedFailure
    def test_set_array(self):
        #Test setting an entire array within one record
        arr = np.zeros(2, dtype=recordwith2darray).view(np.recarray)
        rec = arr[0]
        nbarr = np.zeros(2, dtype=recordwith2darray).view(np.recarray)
        nbrec = nbarr[0]
        for pyfunc in (record_write_full_array, record_write_full_array_alt):
            pyfunc(rec)
            kernel = cuda.jit(pyfunc)
            kernel[1, 1](nbrec)
            np.testing.assert_equal(nbarr, arr)

    @unittest.expectedFailure
    def test_set_arrays(self):
        # Test setting an entire array of arrays (multiple records)
        arr = np.zeros(2, dtype=recordwith2darray).view(np.recarray)
        nbarr = np.zeros(2, dtype=recordwith2darray).view(np.recarray)
        for pyfunc in (
                recarray_write_array_of_nestedarray_broadcast,
                recarray_write_array_of_nestedarray,
        ):
            arr_expected = pyfunc(arr)
            cfunc = self.get_cfunc(pyfunc, nbarr.dtype)
            arr_res = cfunc(nbarr)
            np.testing.assert_equal(arr_res, arr_expected)


if __name__ == '__main__':
    unittest.main()
