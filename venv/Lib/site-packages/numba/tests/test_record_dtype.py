import numpy as np
import ctypes
from numba import jit, literal_unroll, njit, typeof
from numba.core import types
from numba.core.itanium_mangler import mangle_type
from numba.core.errors import TypingError
import unittest
from numba.np import numpy_support
from numba.tests.support import TestCase, skip_ppc64le_issue6465

_FS = ('e', 'f')


def get_a(ary, i):
    return ary[i].a


def get_b(ary, i):
    return ary[i].b


def get_c(ary, i):
    return ary[i].c


def make_getitem(item):
    # This also exercises constant lookup from a closure variable
    def get_xx(ary, i):
        return ary[i][item]
    return get_xx


# Issue #1664: constant index lookup should fall back to regular getitem
def get_zero_a(ary, _unused):
    return ary[0].a


getitem_a = make_getitem('a')
getitem_b = make_getitem('b')
getitem_c = make_getitem('c')
getitem_0 = make_getitem(0)
getitem_1 = make_getitem(1)
getitem_2 = make_getitem(2)
getitem_10 = make_getitem(10)  # OOB


def get_a_subarray(ary, i):
    return ary.a[i]


def get_b_subarray(ary, i):
    return ary.b[i]


def get_c_subarray(ary, i):
    return ary.c[i]


def get_a_zero(ary, _unused):
    return ary.a[0]


def make_getitem_subarray(item):
    # This also exercises constant lookup from a closure variable
    def get_xx_subarray(ary, i):
        return ary[item][i]
    return get_xx_subarray


getitem_a_subarray = make_getitem_subarray('a')
getitem_b_subarray = make_getitem_subarray('b')
getitem_c_subarray = make_getitem_subarray('c')


def get_two_arrays_a(ary1, ary2, i):
    return ary1[i].a + ary2[i].a


def get_two_arrays_b(ary1, ary2, i):
    return ary1[i].b + ary2[i].b


def get_two_arrays_c(ary1, ary2, i):
    return ary1[i].c + ary2[i].c


def get_two_arrays_distinct(ary1, ary2, i):
    return ary1[i].a + ary2[i].f


def set_a(ary, i, v):
    ary[i].a = v


def set_b(ary, i, v):
    ary[i].b = v


def set_c(ary, i, v):
    ary[i].c = v


def make_setitem(item):
    def set_xx(ary, i, v):
        ary[i][item] = v
    return set_xx


setitem_a = make_setitem('a')
setitem_b = make_setitem('b')
setitem_c = make_setitem('c')
setitem_0 = make_setitem(0)
setitem_1 = make_setitem(1)
setitem_2 = make_setitem(2)
setitem_10 = make_setitem(10)  # OOB


def set_a_subarray(ary, i, v):
    ary.a[i] = v


def set_b_subarray(ary, i, v):
    ary.b[i] = v


def set_c_subarray(ary, i, v):
    ary.c[i] = v


def make_setitem_subarray(item):
    def set_xx_subarray(ary, i, v):
        ary[item][i] = v
    return set_xx_subarray


setitem_a_subarray = make_setitem('a')
setitem_b_subarray = make_setitem('b')
setitem_c_subarray = make_setitem('c')


def set_record(ary, i, j):
    ary[i] = ary[j]


def get_record_a(rec, val):
    x = rec.a
    rec.a = val
    return x


def get_record_b(rec, val):
    x = rec.b
    rec.b = val
    return x


def get_record_c(rec, val):
    x = rec.c
    rec.c = val
    return x


def get_record_rev_a(val, rec):
    x = rec.a
    rec.a = val
    return x


def get_record_rev_b(val, rec):
    x = rec.b
    rec.b = val
    return x


def get_record_rev_c(val, rec):
    x = rec.c
    rec.c = val
    return x


def get_two_records_a(rec1, rec2):
    x = rec1.a + rec2.a
    return x


def get_two_records_b(rec1, rec2):
    x = rec1.b + rec2.b
    return x


def get_two_records_c(rec1, rec2):
    x = rec1.c + rec2.c
    return x


def get_two_records_distinct(rec1, rec2):
    x = rec1.a + rec2.f
    return x


def record_return(ary, i):
    return ary[i]


def record_write_array(ary):
    ary.g = 2
    ary.h[0] = 3.0
    ary.h[1] = 4.0


def record_write_2d_array(ary):
    ary.i = 3
    ary.j[0, 0] = 5.0
    ary.j[0, 1] = 6.0
    ary.j[1, 0] = 7.0
    ary.j[1, 1] = 8.0
    ary.j[2, 0] = 9.0
    ary.j[2, 1] = 10.0


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
    ary.j[:, :, :] = np.ones((2, 3, 2), dtype=np.float64)
    return ary


def recarray_getitem_return(ary):
    return ary[0]


def recarray_getitem_field_return(ary):
    return ary['h']


def recarray_getitem_field_return2(ary):
    return ary.h


def recarray_getitem_field_return2_2d(ary):
    return ary.j


def rec_getitem_field_slice_2d(rec):
    return rec.j[0]


def recarray_getitem_field_slice_2d(ary):
    return ary.j[0][0]


def array_rec_getitem_field_slice_2d_0(rec):
    return rec['j'][0]


def array_getitem_field_slice_2d_0(ary):
    return ary['j'][0][0]


def array_rec_getitem_field_slice_2d_1(rec):
    return rec['j'][1]


def array_getitem_field_slice_2d_1(ary):
    return ary['j'][1][0]


def rec_getitem_range_slice_4d(rec):
    return rec.p[0:2, 0, 1:4, 3:6]


def recarray_getitem_range_slice_4d(ary):
    return ary[0].p[0:2, 0, 1:4, 3:6]


def record_read_array0(ary):
    return ary.h[0]


def record_read_array0_alt(ary):
    return ary[0].h


def record_read_array1(ary):
    return ary.h[1]


def record_read_whole_array(ary):
    return ary.h


def record_read_2d_array00(ary):
    return ary.j[0,0]


def record_read_2d_array10(ary):
    return ary.j[1,0]


def record_read_2d_array01(ary):
    return ary.j[0,1]


def record_read_first_arr(ary):
    return ary.k[2, 2]


def record_read_second_arr(ary):
    return ary.l[2, 2]


def get_shape(rec):
    return np.shape(rec.j)


def get_size(rec):
    return np.size(rec.j)


def get_charseq(ary, i):
    return ary[i].n


def set_charseq(ary, i, cs):
    ary[i].n = cs


def get_charseq_tuple(ary, i):
    return ary[i].m, ary[i].n


def get_field1(rec):
    fs = ('e', 'f')
    f = fs[1]
    return rec[f]


def get_field2(rec):
    fs = ('e', 'f')
    out = 0
    for f in literal_unroll(fs):
        out += rec[f]
    return out


def get_field3(rec):
    f = _FS[1]
    return rec[f]


def get_field4(rec):
    out = 0
    for f in literal_unroll(_FS):
        out += rec[f]
    return out


def set_field1(rec):
    fs = ('e', 'f')
    f = fs[1]
    rec[f] = 10
    return rec


def set_field2(rec):
    fs = ('e', 'f')
    for f in literal_unroll(fs):
        rec[f] = 10
    return rec


def set_field3(rec):
    f = _FS[1]
    rec[f] = 10
    return rec


def set_field4(rec):
    for f in literal_unroll(_FS):
        rec[f] = 10
    return rec


def set_field_slice(arr):
    arr['k'][:] = 0.0
    return arr


def assign_array_to_nested(dest):
    tmp = (np.arange(3) + 1).astype(np.int16)
    dest['array1'] = tmp


def assign_array_to_nested_2d(dest):
    tmp = (np.arange(6) + 1).astype(np.int16).reshape((3, 2))
    dest['array2'] = tmp


recordtype = np.dtype([('a', np.float64),
                       ('b', np.int16),
                       ('c', np.complex64),
                       ('d', (np.str_, 5))])

recordtype2 = np.dtype([('e', np.int32),
                        ('f', np.float64)], align=True)

recordtype3 = np.dtype([('first', np.float32),
                        ('second', np.float64)])

recordwitharray = np.dtype([('g', np.int32),
                            ('h', np.float32, 2)])

recordwith2darray = np.dtype([('i', np.int32),
                              ('j', np.float32, (3, 2))])

recordwith2arrays = np.dtype([('k', np.int32, (10, 20)),
                              ('l', np.int32, (6, 12))])

recordwithcharseq = np.dtype([('m', np.int32),
                              ('n', 'S5')])

recordwith4darray = np.dtype([('o', np.int64),
                              ('p', np.float32, (3, 2, 5, 7)),
                              ('q', 'U10'),])

nested_array1_dtype = np.dtype([("array1", np.int16, (3,))], align=True)

nested_array2_dtype = np.dtype([("array2", np.int16, (3, 2))], align=True)


class TestRecordDtypeMakeCStruct(TestCase):
    def test_two_scalars(self):

        class Ref(ctypes.Structure):
            _fields_ = [
                ('apple', ctypes.c_int32),
                ('orange', ctypes.c_float),
            ]

        ty = types.Record.make_c_struct([
            ('apple', types.int32),
            ('orange', types.float32),
        ])
        # Correct offsets
        self.assertEqual(len(ty), 2)
        self.assertEqual(ty.offset('apple'), Ref.apple.offset)
        self.assertEqual(ty.offset('orange'), Ref.orange.offset)
        # Correct size
        self.assertEqual(ty.size, ctypes.sizeof(Ref))
        # Is aligned
        dtype = ty.dtype
        self.assertTrue(dtype.isalignedstruct)

    def test_three_scalars(self):

        class Ref(ctypes.Structure):
            _fields_ = [
                ('apple', ctypes.c_int32),
                ('mango', ctypes.c_int8),
                ('orange', ctypes.c_float),
            ]

        ty = types.Record.make_c_struct([
            ('apple', types.int32),
            ('mango', types.int8),
            ('orange', types.float32),
        ])
        # Correct offsets
        self.assertEqual(len(ty), 3)
        self.assertEqual(ty.offset('apple'), Ref.apple.offset)
        self.assertEqual(ty.offset('mango'), Ref.mango.offset)
        self.assertEqual(ty.offset('orange'), Ref.orange.offset)
        # Correct size
        self.assertEqual(ty.size, ctypes.sizeof(Ref))
        # Is aligned
        dtype = ty.dtype
        self.assertTrue(dtype.isalignedstruct)

    def test_complex_struct(self):
        class Complex(ctypes.Structure):
            _fields_ = [
                ('real', ctypes.c_double),
                ('imag', ctypes.c_double),
            ]

        class Ref(ctypes.Structure):
            _fields_ = [
                ('apple', ctypes.c_int32),
                ('mango', Complex),
            ]

        ty = types.Record.make_c_struct([
            ('apple', types.intc),
            ('mango', types.complex128),
        ])
        # Correct offsets
        self.assertEqual(len(ty), 2)
        self.assertEqual(ty.offset('apple'), Ref.apple.offset)
        self.assertEqual(ty.offset('mango'), Ref.mango.offset)
        # Correct size
        self.assertEqual(ty.size, ctypes.sizeof(Ref))
        # Is aligned?
        self.assertTrue(ty.dtype.isalignedstruct)

    def test_nestedarray_issue_8132(self):
        # issue#8132 is caused by misrepresenting the NestedArray. Instead of
        # using the storage layout for the contained data, it is represented
        # as an array object structure.

        # Make an array that is longer than the array object structure.
        data = np.arange(27 * 2, dtype=np.float64).reshape(27, 2)
        recty = types.Record.make_c_struct([
            ('data', types.NestedArray(dtype=types.float64, shape=data.shape)),
        ])
        arr = np.array((data,), dtype=recty.dtype)
        # unpack the nestedarray as a normal array
        [extracted_array] = arr.tolist()
        # check that is matches the original values
        np.testing.assert_array_equal(extracted_array, data)


class TestRecordDtype(TestCase):

    def _createSampleArrays(self):
        '''
        Set up the data structures to be used with the Numpy and Numba
        versions of functions.

        In this case, both accept recarrays.
        '''
        self.refsample1d = np.recarray(3, dtype=recordtype)
        self.refsample1d2 = np.recarray(3, dtype=recordtype2)
        self.refsample1d3 = np.recarray(3, dtype=recordtype)

        self.nbsample1d = np.recarray(3, dtype=recordtype)
        self.nbsample1d2 = np.recarray(3, dtype=recordtype2)
        self.nbsample1d3 = np.recarray(3, dtype=recordtype)

    def setUp(self):

        self._createSampleArrays()

        for ary in (self.refsample1d, self.nbsample1d):
            for i in range(ary.size):
                x = i + 1
                ary[i]['a'] = x / 2
                ary[i]['b'] = x
                ary[i]['c'] = x * 1j
                ary[i]['d'] = "%d" % x

        for ary2 in (self.refsample1d2, self.nbsample1d2):
            for i in range(ary2.size):
                x = i + 5
                ary2[i]['e'] = x
                ary2[i]['f'] = x / 2

        for ary3 in (self.refsample1d3, self.nbsample1d3):
            for i in range(ary3.size):
                x = i + 10
                ary3[i]['a'] = x / 2
                ary3[i]['b'] = x
                ary3[i]['c'] = x * 1j
                ary3[i]['d'] = "%d" % x

    def get_cfunc(self, pyfunc, argspec):
        return njit(argspec)(pyfunc)

    def test_from_dtype(self):
        rec = numpy_support.from_dtype(recordtype)
        self.assertEqual(rec.typeof('a'), types.float64)
        self.assertEqual(rec.typeof('b'), types.int16)
        self.assertEqual(rec.typeof('c'), types.complex64)
        self.assertEqual(rec.typeof('d'), types.UnicodeCharSeq(5))
        self.assertEqual(rec.offset('a'), recordtype.fields['a'][1])
        self.assertEqual(rec.offset('b'), recordtype.fields['b'][1])
        self.assertEqual(rec.offset('c'), recordtype.fields['c'][1])
        self.assertEqual(rec.offset('d'), recordtype.fields['d'][1])
        self.assertEqual(recordtype.itemsize, rec.size)

    def _test_get_equal(self, pyfunc):
        rec = numpy_support.from_dtype(recordtype)
        cfunc = self.get_cfunc(pyfunc, (rec[:], types.intp))
        for i in range(self.refsample1d.size):
            self.assertEqual(pyfunc(self.refsample1d, i),
                             cfunc(self.nbsample1d, i))

    def test_get_a(self):
        self._test_get_equal(get_a)
        self._test_get_equal(get_a_subarray)
        self._test_get_equal(getitem_a)
        self._test_get_equal(getitem_a_subarray)
        self._test_get_equal(get_a_zero)
        self._test_get_equal(get_zero_a)

    def test_get_b(self):
        self._test_get_equal(get_b)
        self._test_get_equal(get_b_subarray)
        self._test_get_equal(getitem_b)
        self._test_get_equal(getitem_b_subarray)

    def test_get_c(self):
        self._test_get_equal(get_c)
        self._test_get_equal(get_c_subarray)
        self._test_get_equal(getitem_c)
        self._test_get_equal(getitem_c_subarray)

    def test_getitem_static_int_index(self):
        self._test_get_equal(getitem_0)
        self._test_get_equal(getitem_1)
        self._test_get_equal(getitem_2)

        # this exception on accessing OOB integer index
        rec = numpy_support.from_dtype(recordtype)
        with self.assertRaises(TypingError) as raises:
            self.get_cfunc(getitem_10, (rec[:], types.intp))

        msg = "Requested index 10 is out of range"
        self.assertIn(msg, str(raises.exception))

    def _test_get_two_equal(self, pyfunc):
        '''
        Test with two arrays of the same type
        '''
        rec = numpy_support.from_dtype(recordtype)
        cfunc = self.get_cfunc(pyfunc, (rec[:], rec[:], types.intp))
        for i in range(self.refsample1d.size):
            self.assertEqual(pyfunc(self.refsample1d, self.refsample1d3, i),
                             cfunc(self.nbsample1d, self.nbsample1d3, i))

    def test_two_distinct_arrays(self):
        '''
        Test with two arrays of distinct record types
        '''
        pyfunc = get_two_arrays_distinct
        rec1 = numpy_support.from_dtype(recordtype)
        rec2 = numpy_support.from_dtype(recordtype2)
        cfunc = self.get_cfunc(pyfunc, (rec1[:], rec2[:], types.intp))
        for i in range(self.refsample1d.size):
            pres = pyfunc(self.refsample1d, self.refsample1d2, i)
            cres = cfunc(self.nbsample1d, self.nbsample1d2, i)
            self.assertEqual(pres,cres)

    def test_get_two_a(self):
        self._test_get_two_equal(get_two_arrays_a)

    def test_get_two_b(self):
        self._test_get_two_equal(get_two_arrays_b)

    def test_get_two_c(self):
        self._test_get_two_equal(get_two_arrays_c)

    def _test_set_equal(self, pyfunc, value, valuetype):
        rec = numpy_support.from_dtype(recordtype)
        cfunc = self.get_cfunc(pyfunc, (rec[:], types.intp, valuetype))

        for i in range(self.refsample1d.size):
            expect = self.refsample1d.copy()
            pyfunc(expect, i, value)

            got = self.nbsample1d.copy()
            cfunc(got, i, value)

            # Match the entire array to ensure no memory corruption
            np.testing.assert_equal(expect, got)

    def test_set_a(self):
        def check(pyfunc):
            self._test_set_equal(pyfunc, 3.1415, types.float64)
            # Test again to check if coercion works
            self._test_set_equal(pyfunc, 3., types.float32)
        check(set_a)
        check(set_a_subarray)
        check(setitem_a)
        check(setitem_a_subarray)

    def test_set_b(self):
        def check(pyfunc):
            self._test_set_equal(pyfunc, 123, types.int32)
            # Test again to check if coercion works
            self._test_set_equal(pyfunc, 123, types.float64)
        check(set_b)
        check(set_b_subarray)
        check(setitem_b)
        check(setitem_b_subarray)

    def test_set_c(self):
        def check(pyfunc):
            self._test_set_equal(pyfunc, 43j, types.complex64)
            # Test again to check if coercion works
            self._test_set_equal(pyfunc, 43j, types.complex128)
        check(set_c)
        check(set_c_subarray)
        check(setitem_c)
        check(setitem_c_subarray)

    def test_setitem_static_int_index(self):
        def check(pyfunc):
            self._test_set_equal(pyfunc, 3.1415, types.float64)
            # Test again to check if coercion works
            self._test_set_equal(pyfunc, 3., types.float32)

        check(setitem_0)
        check(setitem_1)
        check(setitem_2)

        # this exception on accessing OOB integer index
        rec = numpy_support.from_dtype(recordtype)
        with self.assertRaises(TypingError) as raises:
            self.get_cfunc(setitem_10, (rec[:], types.intp, types.float64))

        msg = "Requested index 10 is out of range"
        self.assertIn(msg, str(raises.exception))

    def test_set_record(self):
        pyfunc = set_record
        rec = numpy_support.from_dtype(recordtype)
        cfunc = self.get_cfunc(pyfunc, (rec[:], types.intp, types.intp))

        test_indices = [(0, 1), (1, 2), (0, 2)]
        for i, j in test_indices:
            expect = self.refsample1d.copy()
            pyfunc(expect, i, j)

            got = self.nbsample1d.copy()
            cfunc(got, i, j)

            # Match the entire array to ensure no memory corruption
            self.assertEqual(expect[i], expect[j])
            self.assertEqual(got[i], got[j])
            np.testing.assert_equal(expect, got)

    def _test_record_args(self, revargs):
        # Testing scalar record value as argument

        npval = self.refsample1d.copy()[0]
        nbval = self.nbsample1d.copy()[0]
        attrs = 'abc'
        valtypes = types.float64, types.int16, types.complex64
        values = 1.23, 12345, 123 + 456j
        # Check for potential leaks (issue #441)
        with self.assertRefCount(nbval):
            for attr, valtyp, val in zip(attrs, valtypes, values):
                expected = getattr(npval, attr)
                nbrecord = numpy_support.from_dtype(recordtype)

                # Test with a record as either the first argument or the second
                # argument (issue #870)
                if revargs:
                    prefix = 'get_record_rev_'
                    argtypes = (valtyp, nbrecord)
                    args = (val, nbval)
                else:
                    prefix = 'get_record_'
                    argtypes = (nbrecord, valtyp)
                    args = (nbval, val)

                pyfunc = globals()[prefix + attr]
                cfunc = self.get_cfunc(pyfunc, argtypes)

                got = cfunc(*args)
                try:
                    self.assertEqual(expected, got)
                except AssertionError:
                    # On ARM, a LLVM misoptimization can produce buggy code,
                    # see https://llvm.org/bugs/show_bug.cgi?id=24669
                    import llvmlite.binding as ll
                    if attr != 'c':
                        raise
                    triple = 'armv7l-unknown-linux-gnueabihf'
                    if ll.get_default_triple() != triple:
                        raise
                    self.assertEqual(val, got)
                else:
                    self.assertEqual(nbval[attr], val)
                del got, expected, args

    def test_record_args(self):
        self._test_record_args(False)

    def test_record_args_reverse(self):
        self._test_record_args(True)

    def test_two_records(self):
        '''
        Testing the use of two scalar records of the same type
        '''
        npval1 = self.refsample1d.copy()[0]
        npval2 = self.refsample1d.copy()[1]
        nbval1 = self.nbsample1d.copy()[0]
        nbval2 = self.nbsample1d.copy()[1]
        attrs = 'abc'
        valtypes = types.float64, types.int32, types.complex64

        for attr, valtyp in zip(attrs, valtypes):
            expected = getattr(npval1, attr) + getattr(npval2, attr)

            nbrecord = numpy_support.from_dtype(recordtype)
            pyfunc = globals()['get_two_records_' + attr]
            cfunc = self.get_cfunc(pyfunc, (nbrecord, nbrecord))

            got = cfunc(nbval1, nbval2)
            self.assertEqual(expected, got)

    def test_two_distinct_records(self):
        '''
        Testing the use of two scalar records of differing type
        '''
        nbval1 = self.nbsample1d.copy()[0]
        nbval2 = self.refsample1d2.copy()[0]
        expected = nbval1['a'] + nbval2['f']

        nbrecord1 = numpy_support.from_dtype(recordtype)
        nbrecord2 = numpy_support.from_dtype(recordtype2)
        cfunc = self.get_cfunc(get_two_records_distinct, (nbrecord1, nbrecord2))

        got = cfunc(nbval1, nbval2)
        self.assertEqual(expected, got)

    def test_record_return(self):
        # Testing scalar record value as return value.
        # We can only return a copy of the record.
        pyfunc = record_return
        recty = numpy_support.from_dtype(recordtype)
        cfunc = self.get_cfunc(pyfunc, (recty[:], types.intp))

        attrs = 'abc'
        indices = [0, 1, 2]
        for index, attr in zip(indices, attrs):
            nbary = self.nbsample1d.copy()
            # Check for potential leaks
            with self.assertRefCount(nbary):
                res = cfunc(nbary, index)
                self.assertEqual(nbary[index], res)
                # Prove that this is a by-value copy
                setattr(res, attr, 0)
                self.assertNotEqual(nbary[index], res)
                del res

    def test_record_arg_transform(self):
        # Testing that transforming the name of a record type argument to a
        # function does not result in the fields of the record being used to
        # uniquely identify them, and that no other condition results in the
        # transformed name being excessively long.

        rec = numpy_support.from_dtype(recordtype3)
        transformed = mangle_type(rec)
        self.assertNotIn('first', transformed)
        self.assertNotIn('second', transformed)
        # len(transformed) is generally 10, but could be longer if a large
        # number of typecodes are in use. Checking <20 should provide enough
        # tolerance.
        self.assertLess(len(transformed), 20)

        struct_arr = types.Array(rec, 1, 'C')
        transformed = mangle_type(struct_arr)
        self.assertIn('Array', transformed)
        self.assertNotIn('first', transformed)
        self.assertNotIn('second', transformed)
        # Length is usually 50 - 5 chars tolerance as above.
        self.assertLess(len(transformed), 50)

    def test_record_two_arrays(self):
        # Tests that comparison of NestedArrays by key is working correctly. If
        # the two NestedArrays in recordwith2arrays compare equal (same length
        # and ndim but different shape) incorrect code will be generated for one
        # of the functions.

        nbrecord = numpy_support.from_dtype(recordwith2arrays)
        rec = np.recarray(1, dtype=recordwith2arrays)[0]
        rec.k[:] = np.arange(200).reshape(10,20)
        rec.l[:] = np.arange(72).reshape(6,12)

        pyfunc = record_read_first_arr
        cfunc = self.get_cfunc(pyfunc, (nbrecord,))
        self.assertEqual(cfunc(rec), pyfunc(rec))

        pyfunc = record_read_second_arr
        cfunc = self.get_cfunc(pyfunc, (nbrecord,))
        self.assertEqual(cfunc(rec), pyfunc(rec))

        pyfunc = set_field_slice
        cfunc = self.get_cfunc(pyfunc, (nbrecord,))
        np.testing.assert_array_equal(cfunc(rec), pyfunc(rec))

    def test_structure_dtype_with_titles(self):
        # the following is the definition of int4 vector type from pyopencl
        vecint4 = np.dtype([(('x', 's0'), 'i4'), (('y', 's1'), 'i4'),
                            (('z', 's2'), 'i4'), (('w', 's3'), 'i4')])
        nbtype = numpy_support.from_dtype(vecint4)
        self.assertEqual(len(nbtype.fields), len(vecint4.fields))

        arr = np.zeros(10, dtype=vecint4)

        def pyfunc(a):
            for i in range(a.size):
                j = i + 1
                a[i]['s0'] = j * 2
                a[i]['x'] += -1

                a[i]['s1'] = j * 3
                a[i]['y'] += -2

                a[i]['s2'] = j * 4
                a[i]['z'] += -3

                a[i]['s3'] = j * 5
                a[i]['w'] += -4

            return a

        expect = pyfunc(arr.copy())
        cfunc = self.get_cfunc(pyfunc, (nbtype[:],))
        got = cfunc(arr.copy())
        np.testing.assert_equal(expect, got)

    def test_record_dtype_with_titles_roundtrip(self):
        recdtype = np.dtype([(("title a", 'a'), np.float64), ('b', np.float64)])
        nbtype = numpy_support.from_dtype(recdtype)
        self.assertTrue(nbtype.is_title('title a'))
        self.assertFalse(nbtype.is_title('a'))
        self.assertFalse(nbtype.is_title('b'))
        got = numpy_support.as_dtype(nbtype)
        self.assertTrue(got, recdtype)


def _get_cfunc_nopython(pyfunc, argspec):
    return jit(argspec, nopython=True)(pyfunc)


class TestRecordDtypeWithDispatcher(TestRecordDtype):
    '''
    Same as TestRecordDtype, but stressing the Dispatcher's type dispatch
    mechanism (issue #384). Note that this does not stress caching of ndarray
    typecodes as the path that uses the cache is not taken with recarrays.
    '''

    def get_cfunc(self, pyfunc, argspec):
        return _get_cfunc_nopython(pyfunc, argspec)


class TestRecordDtypeWithStructArrays(TestRecordDtype):
    '''
    Same as TestRecordDtype, but using structured arrays instead of recarrays.
    '''

    def _createSampleArrays(self):
        '''
        Two different versions of the data structures are required because Numba
        supports attribute access on structured arrays, whereas Numpy does not.

        However, the semantics of recarrays and structured arrays are equivalent
        for these tests so Numpy with recarrays can be used for comparison with
        Numba using structured arrays.
        '''

        self.refsample1d = np.recarray(3, dtype=recordtype)
        self.refsample1d2 = np.recarray(3, dtype=recordtype2)
        self.refsample1d3 = np.recarray(3, dtype=recordtype)

        self.nbsample1d = np.zeros(3, dtype=recordtype)
        self.nbsample1d2 = np.zeros(3, dtype=recordtype2)
        self.nbsample1d3 = np.zeros(3, dtype=recordtype)


class TestRecordDtypeWithStructArraysAndDispatcher(TestRecordDtypeWithStructArrays):    # noqa: E501
    '''
    Same as TestRecordDtypeWithStructArrays, stressing the Dispatcher's type
    dispatch mechanism (issue #384) and caching of ndarray typecodes for void
    types (which occur in structured arrays).
    '''

    def get_cfunc(self, pyfunc, argspec):
        return _get_cfunc_nopython(pyfunc, argspec)


@skip_ppc64le_issue6465
class TestRecordDtypeWithCharSeq(TestCase):

    def _createSampleaArray(self):
        self.refsample1d = np.recarray(3, dtype=recordwithcharseq)
        self.nbsample1d = np.zeros(3, dtype=recordwithcharseq)

    def _fillData(self, arr):
        for i in range(arr.size):
            arr[i]['m'] = i

        arr[0]['n'] = 'abcde'  # no null-byte
        arr[1]['n'] = 'xyz'  # null-byte
        arr[2]['n'] = 'u\x00v\x00\x00'  # null-byte at the middle and at the end

    def setUp(self):
        self._createSampleaArray()
        self._fillData(self.refsample1d)
        self._fillData(self.nbsample1d)

    def get_cfunc(self, pyfunc):
        rectype = numpy_support.from_dtype(recordwithcharseq)
        return njit((rectype[:], types.intp))(pyfunc)

    def test_return_charseq(self):
        pyfunc = get_charseq
        cfunc = self.get_cfunc(pyfunc)
        for i in range(self.refsample1d.size):
            expected = pyfunc(self.refsample1d, i)
            got = cfunc(self.nbsample1d, i)
            self.assertEqual(expected, got)

    def test_npm_argument_charseq(self):
        # Test CharSeq as NPM argument

        def pyfunc(arr, i):
            return arr[i].n

        identity = njit(lambda x: x)   # an identity function

        @jit(nopython=True)
        def cfunc(arr, i):
            return identity(arr[i].n)

        for i in range(self.refsample1d.size):
            expected = pyfunc(self.refsample1d, i)
            got = cfunc(self.nbsample1d, i)
            self.assertEqual(expected, got)

    def test_py_argument_charseq(self):
        # Test CharSeq as python wrapper argument

        pyfunc = set_charseq

        # compile
        rectype = numpy_support.from_dtype(recordwithcharseq)
        sig = (rectype[::1], types.intp, rectype.typeof('n'))
        cfunc = njit(sig)(pyfunc).overloads[sig].entry_point

        for i in range(self.refsample1d.size):
            chars = "{0}".format(hex(i + 10))
            pyfunc(self.refsample1d, i, chars)
            cfunc(self.nbsample1d, i, chars)
            np.testing.assert_equal(self.refsample1d, self.nbsample1d)

    def test_py_argument_char_seq_near_overflow(self):
        # Test strings that are as long as the charseq capacity

        pyfunc = set_charseq
        # compile
        rectype = numpy_support.from_dtype(recordwithcharseq)
        sig = (rectype[::1], types.intp, rectype.typeof('n'))
        cfunc = njit(sig)(pyfunc).overloads[sig].entry_point

        cs_near_overflow = "abcde"

        self.assertEqual(len(cs_near_overflow),
                         recordwithcharseq['n'].itemsize)

        cfunc(self.nbsample1d, 0, cs_near_overflow)
        self.assertEqual(self.nbsample1d[0]['n'].decode('ascii'),
                         cs_near_overflow)
        # Check that we didn't overwrite
        np.testing.assert_equal(self.refsample1d[1:], self.nbsample1d[1:])

    def test_py_argument_char_seq_truncate(self):
        # NumPy silently truncates strings to fix inside charseq

        pyfunc = set_charseq
        # compile
        rectype = numpy_support.from_dtype(recordwithcharseq)
        sig = (rectype[::1], types.intp, rectype.typeof('n'))
        cfunc = njit(sig)(pyfunc).overloads[sig].entry_point

        cs_overflowed = "abcdef"

        pyfunc(self.refsample1d, 1, cs_overflowed)
        cfunc(self.nbsample1d, 1, cs_overflowed)
        np.testing.assert_equal(self.refsample1d, self.nbsample1d)
        self.assertEqual(self.refsample1d[1].n,
                         cs_overflowed[:-1].encode("ascii"))

    def test_return_charseq_tuple(self):
        pyfunc = get_charseq_tuple
        cfunc = self.get_cfunc(pyfunc)
        for i in range(self.refsample1d.size):
            expected = pyfunc(self.refsample1d, i)
            got = cfunc(self.nbsample1d, i)
            self.assertEqual(expected, got)


class TestRecordArrayGetItem(TestCase):

    # Test getitem when index is Literal[str]

    def test_literal_variable(self):
        arr = np.array([1, 2], dtype=recordtype2)
        pyfunc = get_field1
        jitfunc = njit(pyfunc)
        self.assertEqual(pyfunc(arr[0]), jitfunc(arr[0]))

    def test_literal_unroll(self):
        arr = np.array([1, 2], dtype=recordtype2)
        pyfunc = get_field2
        jitfunc = njit(pyfunc)
        self.assertEqual(pyfunc(arr[0]), jitfunc(arr[0]))

    def test_literal_variable_global_tuple(self):
        # This tests the getitem of record array when the indexes come from a
        # global tuple. It tests getitem behaviour but also tests that a global
        # tuple is being typed as a tuple of constants.

        arr = np.array([1, 2], dtype=recordtype2)
        pyfunc = get_field3
        jitfunc = njit(pyfunc)
        self.assertEqual(pyfunc(arr[0]), jitfunc(arr[0]))

    def test_literal_unroll_global_tuple(self):
        # This tests the getitem of record array when the indexes come from a
        # global tuple and are being unrolled.
        # It tests getitem behaviour but also tests that literal_unroll accepts
        # a global tuple as argument

        arr = np.array([1, 2], dtype=recordtype2)
        pyfunc = get_field4
        jitfunc = njit(pyfunc)
        self.assertEqual(pyfunc(arr[0]), jitfunc(arr[0]))

    def test_literal_unroll_free_var_tuple(self):
        # This tests the getitem of record array when the indexes come from a
        # free variable tuple (not local, not global) and are being unrolled.
        # It tests getitem behaviour but also tests that literal_unroll accepts
        # a free variable tuple as argument

        fs = ('e', 'f')
        arr = np.array([1, 2], dtype=recordtype2)

        def get_field(rec):
            out = 0
            for f in literal_unroll(fs):
                out += rec[f]
            return out

        jitfunc = njit(get_field)
        self.assertEqual(get_field(arr[0]), jitfunc(arr[0]))

    def test_error_w_invalid_field(self):
        arr = np.array([1, 2], dtype=recordtype3)
        jitfunc = njit(get_field1)
        with self.assertRaises(TypingError) as raises:
            jitfunc(arr[0])
        self.assertIn("Field 'f' was not found in record with fields "
                      "('first', 'second')", str(raises.exception))

    def test_literal_unroll_dynamic_to_static_getitem_transform(self):
        # See issue #6634
        keys = ('a', 'b', 'c')
        n = 5

        def pyfunc(rec):
            x = np.zeros((n,))
            for o in literal_unroll(keys):
                x += rec[o]
            return x

        dt = np.float64
        ldd = [np.arange(dt(n)) for x in keys]
        ldk = [(x, np.float64,) for x in keys]
        rec = np.rec.fromarrays(ldd, dtype=ldk)

        expected = pyfunc(rec)
        got = njit(pyfunc)(rec)
        np.testing.assert_allclose(expected, got)


class TestRecordArraySetItem(TestCase):
    """
    Test setitem when index is Literal[str]
    """
    def test_literal_variable(self):
        arr = np.array([1, 2], dtype=recordtype2)
        pyfunc = set_field1
        jitfunc = njit(pyfunc)
        self.assertEqual(pyfunc(arr[0].copy()), jitfunc(arr[0].copy()))

    def test_literal_unroll(self):
        arr = np.array([1, 2], dtype=recordtype2)
        pyfunc = set_field2
        jitfunc = njit(pyfunc)
        self.assertEqual(pyfunc(arr[0].copy()), jitfunc(arr[0].copy()))

    def test_literal_variable_global_tuple(self):
        # This tests the setitem of record array when the indexes come from a
        # global tuple. It tests getitem behaviour but also tests that a global
        # tuple is being typed as a tuple of constants.

        arr = np.array([1, 2], dtype=recordtype2)
        pyfunc = set_field3
        jitfunc = njit(pyfunc)
        self.assertEqual(pyfunc(arr[0].copy()), jitfunc(arr[0].copy()))

    def test_literal_unroll_global_tuple(self):
        # This tests the setitem of record array when the indexes come from a
        # global tuple and are being unrolled.
        # It tests setitem behaviour but also tests that literal_unroll accepts
        # a global tuple as argument

        arr = np.array([1, 2], dtype=recordtype2)
        pyfunc = set_field4
        jitfunc = njit(pyfunc)
        self.assertEqual(pyfunc(arr[0].copy()), jitfunc(arr[0].copy()))

    def test_literal_unroll_free_var_tuple(self):

        # This tests the setitem of record array when the indexes come from a
        # free variable tuple (not local, not global) and are being unrolled.
        # It tests setitem behaviour but also tests that literal_unroll accepts
        # a free variable tuple as argument

        arr = np.array([1, 2], dtype=recordtype2)
        fs = arr.dtype.names

        def set_field(rec):
            for f in literal_unroll(fs):
                rec[f] = 10
            return rec

        jitfunc = njit(set_field)
        self.assertEqual(set_field(arr[0].copy()), jitfunc(arr[0].copy()))

    def test_error_w_invalid_field(self):
        arr = np.array([1, 2], dtype=recordtype3)
        jitfunc = njit(set_field1)
        with self.assertRaises(TypingError) as raises:
            jitfunc(arr[0])
        self.assertIn("Field 'f' was not found in record with fields "
                      "('first', 'second')", str(raises.exception))


class TestSubtyping(TestCase):
    def setUp(self):
        self.value = 2
        a_dtype = np.dtype([('a', 'f8')])
        ab_dtype = np.dtype([('a', 'f8'), ('b', 'f8')])
        self.a_rec1 = np.array([1], dtype=a_dtype)[0]
        self.a_rec2 = np.array([2], dtype=a_dtype)[0]
        self.ab_rec1 = np.array([(self.value, 3)], dtype=ab_dtype)[0]
        self.ab_rec2 = np.array([(self.value + 1, 3)], dtype=ab_dtype)[0]
        self.func = lambda rec: rec['a']

    def test_common_field(self):
        # Test that subtypes do not require new compilations

        njit_sig = njit(types.float64(typeof(self.a_rec1)))
        functions = [
            njit(self.func),  # jitted function with open njit
            njit_sig(self.func)  # jitted fc with closed signature
        ]

        for fc in functions:
            fc(self.a_rec1)
            fc.disable_compile()
            y = fc(self.ab_rec1)
            self.assertEqual(self.value, y)

    def test_tuple_of_records(self):

        @njit
        def foo(rec_tup):
            x = 0
            for i in range(len(rec_tup)):
                x += rec_tup[i]['a']
            return x

        foo((self.a_rec1, self.a_rec2))
        foo.disable_compile()
        y = foo((self.ab_rec1, self.ab_rec2))
        self.assertEqual(2 * self.value + 1, y)

    def test_array_field(self):
        # Tests subtyping with array fields

        rec1 = np.empty(1, dtype=[('a', 'f8', (4,))])[0]
        rec1['a'][0] = 1
        rec2 = np.empty(1, dtype=[('a', 'f8', (4,)), ('b', 'f8')])[0]
        rec2['a'][0] = self.value

        @njit
        def foo(rec):
            return rec['a'][0]

        foo(rec1)
        foo.disable_compile()
        y = foo(rec2)
        self.assertEqual(self.value, y)

    def test_no_subtyping1(self):
        # test that conversion rules don't allow subtypes with different field
        # names

        c_dtype = np.dtype([('c', 'f8')])
        c_rec1 = np.array([1], dtype=c_dtype)[0]

        @njit
        def foo(rec):
            return rec['c']

        foo(c_rec1)
        foo.disable_compile()
        with self.assertRaises(TypeError) as err:
            foo(self.a_rec1)
            self.assertIn("No matching definition for argument type(s) Record",
                          str(err.exception))

    def test_no_subtyping2(self):
        # test that conversion rules don't allow smaller records as subtypes

        jit_fc = njit(self.func)
        jit_fc(self.ab_rec1)
        jit_fc.disable_compile()
        with self.assertRaises(TypeError) as err:
            jit_fc(self.a_rec1)
            self.assertIn("No matching definition for argument type(s) Record",
                          str(err.exception))

    def test_no_subtyping3(self):
        # test that conversion rules don't allow records with fields with same
        # name but incompatible type

        other_a_rec = np.array(['a'], dtype=np.dtype([('a', 'U25')]))[0]
        jit_fc = njit(self.func)
        jit_fc(self.a_rec1)
        jit_fc.disable_compile()
        with self.assertRaises(TypeError) as err:
            jit_fc(other_a_rec)
            self.assertIn("No matching definition for argument type(s) Record",
                          str(err.exception))

    def test_branch_pruning(self):
        # test subtyping behaviour in a case with a dead branch

        @njit
        def foo(rec, flag=None):
            n = 0
            n += rec['a']
            if flag is not None:
                # Dead branch pruning will hide this branch
                n += rec['b']
                rec['b'] += 20
            return n

        self.assertEqual(foo(self.a_rec1), self.a_rec1[0])

        # storing value because it will be mutated
        k = self.ab_rec1[1]
        self.assertEqual(foo(self.ab_rec1, flag=1), self.ab_rec1[0] + k)
        self.assertEqual(self.ab_rec1[1], k + 20)

        foo.disable_compile()
        self.assertEqual(len(foo.nopython_signatures), 2)
        self.assertEqual(foo(self.a_rec1) + 1, foo(self.ab_rec1))
        self.assertEqual(foo(self.ab_rec1, flag=1), self.ab_rec1[0] + k + 20)


class TestNestedArrays(TestCase):

    def get_cfunc(self, pyfunc, argspec):
        return njit(argspec)(pyfunc)

    def test_record_write_array(self):
        # Testing writing to a 1D array within a structured type

        nbval = np.recarray(1, dtype=recordwitharray)
        nbrecord = numpy_support.from_dtype(recordwitharray)
        cfunc = self.get_cfunc(record_write_array, (nbrecord,))
        cfunc(nbval[0])

        expected = np.recarray(1, dtype=recordwitharray)
        expected[0].g = 2
        expected[0].h[0] = 3.0
        expected[0].h[1] = 4.0
        np.testing.assert_equal(expected, nbval)

    def test_record_write_2d_array(self):
        # Test writing to a 2D array within a structured type

        nbval = np.recarray(1, dtype=recordwith2darray)
        nbrecord = numpy_support.from_dtype(recordwith2darray)
        cfunc = self.get_cfunc(record_write_2d_array, (nbrecord,))
        cfunc(nbval[0])

        expected = np.recarray(1, dtype=recordwith2darray)
        expected[0].i = 3
        expected[0].j[:] = np.asarray([5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
                                      np.float32).reshape(3, 2)
        np.testing.assert_equal(expected, nbval)

    def test_record_read_array(self):
        # Test reading from a 1D array within a structured type

        nbval = np.recarray(1, dtype=recordwitharray)
        nbval[0].h[0] = 15.0
        nbval[0].h[1] = 25.0
        nbrecord = numpy_support.from_dtype(recordwitharray)
        cfunc = self.get_cfunc(record_read_array0, (nbrecord,))
        res = cfunc(nbval[0])
        np.testing.assert_equal(res, nbval[0].h[0])

        cfunc = self.get_cfunc(record_read_array1, (nbrecord,))
        res = cfunc(nbval[0])
        np.testing.assert_equal(res, nbval[0].h[1])

    def test_record_read_arrays(self):
        # Test reading the array of 1D array
        #
        # This array or arrays is created by getitem on the
        # record array (as opposed to getitem on a record)

        nbval = np.recarray(2, dtype=recordwitharray)
        nbval[0].h[0] = 15.0
        nbval[0].h[1] = 25.0
        nbval[1].h[0] = 35.0
        nbval[1].h[1] = 45.4
        ty_arr = typeof(nbval)
        cfunc = self.get_cfunc(record_read_whole_array, (ty_arr,))
        res = cfunc(nbval)
        np.testing.assert_equal(res, nbval.h)

    def test_record_read_2d_array(self):
        # Test reading from a 2D array within a structured type

        nbval = np.recarray(1, dtype=recordwith2darray)
        nbval[0].j = np.asarray([1.5, 2.5, 3.5, 4.5, 5.5, 6.5],
                                np.float32).reshape(3, 2)
        nbrecord = numpy_support.from_dtype(recordwith2darray)
        cfunc = self.get_cfunc(record_read_2d_array00, (nbrecord,))
        res = cfunc(nbval[0])
        np.testing.assert_equal(res, nbval[0].j[0, 0])

        cfunc = self.get_cfunc(record_read_2d_array01, (nbrecord,))
        res = cfunc(nbval[0])
        np.testing.assert_equal(res, nbval[0].j[0, 1])

        cfunc = self.get_cfunc(record_read_2d_array10, (nbrecord,))
        res = cfunc(nbval[0])
        np.testing.assert_equal(res, nbval[0].j[1, 0])

    def test_set_record(self):
        # Test setting an entire record

        rec = np.ones(2, dtype=recordwith2darray).view(np.recarray)[0]
        nbarr = np.zeros(2, dtype=recordwith2darray).view(np.recarray)
        arr = np.zeros(2, dtype=recordwith2darray).view(np.recarray)
        ty_arr = typeof(nbarr)
        ty_rec = typeof(rec)
        pyfunc = recarray_set_record
        pyfunc(arr, rec)
        cfunc = self.get_cfunc(pyfunc, (ty_arr, ty_rec))
        cfunc(nbarr, rec)
        np.testing.assert_equal(nbarr, arr)

    def test_set_array(self):
        # Test setting an entire array within one record

        arr = np.zeros(2, dtype=recordwith2darray).view(np.recarray)
        rec = arr[0]
        nbarr = np.zeros(2, dtype=recordwith2darray).view(np.recarray)
        nbrec = nbarr[0]
        ty = typeof(nbrec)
        for pyfunc in (record_write_full_array, record_write_full_array_alt):
            pyfunc(rec)
            cfunc = self.get_cfunc(pyfunc, (ty,))
            cfunc(nbrec)
            np.testing.assert_equal(nbarr, arr)

    def test_set_arrays(self):
        # Test setting an entire array of arrays (multiple records)

        arr = np.zeros(2, dtype=recordwith2darray).view(np.recarray)
        nbarr = np.zeros(2, dtype=recordwith2darray).view(np.recarray)
        ty = typeof(nbarr)
        for pyfunc in (
                recarray_write_array_of_nestedarray_broadcast,
                recarray_write_array_of_nestedarray,
        ):
            arr_expected = pyfunc(arr)
            cfunc = self.get_cfunc(pyfunc, (ty,))
            arr_res = cfunc(nbarr)
            np.testing.assert_equal(arr_res, arr_expected)

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

        nbargs = gen()
        cfunc = self.get_cfunc(pyfunc, tuple((typeof(arg) for arg in nbargs)))
        cfunc(*nbargs)
        np.testing.assert_equal(pyargs, nbargs)

    def test_setitem_whole_array_error(self):
        # Ensure we raise a suitable error when attempting to assign an
        # array to a whole array's worth of nested arrays.
        nbarr1 = np.recarray(1, dtype=recordwith2darray)
        nbarr2 = np.recarray(1, dtype=recordwith2darray)
        args = (nbarr1, nbarr2)
        pyfunc = record_setitem_array
        errmsg = "Unsupported array index type"
        with self.assertRaisesRegex(TypingError, errmsg):
            self.get_cfunc(pyfunc, tuple((typeof(arg) for arg in args)))

    def test_getitem_idx(self):
        # Test __getitem__ with numerical index

        # This test returning a record when passing an array and
        # return the first item when passing a record

        nbarr = np.recarray(2, dtype=recordwitharray)
        nbarr[0] = np.array([(1, (2, 3))], dtype=recordwitharray)[0]
        for arg in [nbarr, nbarr[0]]:
            ty = typeof(arg)
            pyfunc = recarray_getitem_return
            arr_expected = pyfunc(arg)
            cfunc = self.get_cfunc(pyfunc, (ty,))
            arr_res = cfunc(arg)
            np.testing.assert_equal(arr_res, arr_expected)

    def test_getitem_idx_2darray(self):
        # Test __getitem__ with numerical index
        #
        # This test returning a record when passing an array and
        # return the first item when passing a record

        nbarr = np.recarray(2, dtype=recordwith2darray)
        nbarr[0] = np.array([(1, ((1, 2), (4, 5), (2, 3)))],
                            dtype=recordwith2darray)[0]
        for arg in [nbarr, nbarr[0]]:
            ty = typeof(arg)
            pyfunc = recarray_getitem_field_return2_2d
            arr_expected = pyfunc(arg)
            cfunc = self.get_cfunc(pyfunc, (ty,))
            arr_res = cfunc(arg)
            np.testing.assert_equal(arr_res, arr_expected)

    def test_return_getattr_getitem_fieldname(self):
        # Test __getitem__ with field name and getattr .field_name
        #
        # This tests returning a array of nestedarrays when passing an array and
        # returning a nestedarray when passing a record

        nbarr = np.recarray(2, dtype=recordwitharray)
        nbarr[0] = np.array([(1, (2, 3))], dtype=recordwitharray)[0]
        for arg in [nbarr, nbarr[0]]:
            for pyfunc in [recarray_getitem_field_return,
                           recarray_getitem_field_return2]:
                ty = typeof(arg)
                arr_expected = pyfunc(arg)
                cfunc = self.get_cfunc(pyfunc, (ty,))
                arr_res = cfunc(arg)
                np.testing.assert_equal(arr_res, arr_expected)

    def test_return_array(self):
        # Test getitem record AND array within record and returning it

        nbval = np.recarray(2, dtype=recordwitharray)
        nbval[0] = np.array([(1, (2, 3))], dtype=recordwitharray)[0]
        for pyfunc in [record_read_array0, record_read_array0_alt]:
            arr_expected = pyfunc(nbval)
            cfunc = self.get_cfunc(pyfunc, (typeof(nbval),))
            arr_res = cfunc(nbval)
            np.testing.assert_equal(arr_expected, arr_res)

    def test_slice_2d_array(self):
        # test slicing the nestedarray inside a record
        nbarr = np.recarray(2, dtype=recordwith2darray)
        nbarr[0] = np.array([(1, ((1, 2), (4, 5), (2, 3)))],
                            dtype=recordwith2darray)[0]

        funcs = rec_getitem_field_slice_2d, recarray_getitem_field_slice_2d
        for arg, pyfunc in zip([nbarr[0], nbarr], funcs):
            ty = typeof(arg)
            arr_expected = pyfunc(arg)
            cfunc = self.get_cfunc(pyfunc, (ty,))
            arr_res = cfunc(arg)
            np.testing.assert_equal(arr_res, arr_expected)

    def test_shape(self):
        # test getting the shape of a nestedarray inside a record
        nbarr = np.recarray(2, dtype=recordwith2darray)
        nbarr[0] = np.array([(1, ((1, 2), (4, 5), (2, 3)))],
                            dtype=recordwith2darray)[0]

        arg = nbarr[0]
        pyfunc = get_shape
        ty = typeof(arg)
        arr_expected = pyfunc(arg)
        cfunc = self.get_cfunc(pyfunc, (ty,))
        arr_res = cfunc(arg)
        np.testing.assert_equal(arr_res, arr_expected)

    def test_size(self):
        # test getting the size of a nestedarray inside a record
        nbarr = np.recarray(2, dtype=recordwith2darray)
        nbarr[0] = np.array([(1, ((1, 2), (4, 5), (2, 3)))],
                            dtype=recordwith2darray)[0]

        arg = nbarr[0]
        pyfunc = get_size
        ty = typeof(arg)
        arr_expected = pyfunc(arg)
        cfunc = self.get_cfunc(pyfunc, (ty,))
        arr_res = cfunc(arg)
        np.testing.assert_equal(arr_res, arr_expected)

    def test_corner_slice(self):
        # testing corner cases while slicing nested arrays
        nbarr = np.recarray((1, 2, 3, 5, 7, 13, 17), dtype=recordwith4darray,
                            order='F')
        np.random.seed(1)
        for index, _ in np.ndenumerate(nbarr):
            nbarr[index].p = np.random.randint(0, 1000, (3, 2, 5, 7),
                                               np.int64).astype(np.float32)

        funcs = rec_getitem_range_slice_4d, recarray_getitem_range_slice_4d
        for arg, pyfunc in zip([nbarr[0], nbarr], funcs):
            ty = typeof(arg)
            arr_expected = pyfunc(arg)
            cfunc = self.get_cfunc(pyfunc, (ty,))
            arr_res = cfunc(arg)
            np.testing.assert_equal(arr_res, arr_expected)

    def test_broadcast_slice(self):
        nbarr = np.recarray(2, dtype=recordwith2darray)
        nbarr[0] = np.array([(1, ((1, 2), (4, 5), (2, 3)))],
                            dtype=recordwith2darray)[0]
        nbarr[1] = np.array([(10, ((10, 20), (40, 50), (20, 30)))],
                            dtype=recordwith2darray)[0]
        nbarr = np.broadcast_to(nbarr, (3, 2))

        funcs = (
            array_rec_getitem_field_slice_2d_0,
            array_getitem_field_slice_2d_0,
            array_rec_getitem_field_slice_2d_1,
            array_getitem_field_slice_2d_1
        )

        for arg, pyfunc in zip([nbarr[0], nbarr, nbarr[1], nbarr], funcs):
            ty = typeof(arg)
            arr_expected = pyfunc(arg)
            cfunc = self.get_cfunc(pyfunc, (ty,))
            arr_res = cfunc(arg)
            np.testing.assert_equal(arr_res, arr_expected)

    def test_assign_array_to_nested(self):
        got = np.zeros(2, dtype=nested_array1_dtype)
        expected = np.zeros(2, dtype=nested_array1_dtype)

        cfunc = njit(assign_array_to_nested)
        cfunc(got[0])
        assign_array_to_nested(expected[0])

        np.testing.assert_array_equal(expected, got)

    def test_assign_array_to_nested_2d(self):
        got = np.zeros(2, dtype=nested_array2_dtype)
        expected = np.zeros(2, dtype=nested_array2_dtype)

        cfunc = njit(assign_array_to_nested_2d)
        cfunc(got[0])
        assign_array_to_nested_2d(expected[0])

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

        @njit
        def copy(index, src, dest):
            dest['user1'] = src[index]['user']
            dest['array1'] = src[index]['array']

        source = np.zeros(2, dtype=src_dtype)
        got = np.zeros(2, dtype=dest_dtype)
        expected = np.zeros(2, dtype=dest_dtype)

        source[0] = (1.2, [1, 2, 3])
        copy(0, source, got[0])
        copy.py_func(0, source, expected[0])

        np.testing.assert_array_equal(expected, got)

    def test_issue_1469_1(self):
        # Dimensions of nested arrays as a dtype are concatenated with the
        # dimensions of the array.
        nptype = np.dtype((np.float64, (4,)))
        nbtype = types.Array(numpy_support.from_dtype(nptype), 2, 'C')
        expected = types.Array(types.float64, 3, 'C')
        self.assertEqual(nbtype, expected)

    def test_issue_1469_2(self):
        # Nesting nested arrays results in concatentated dimensions -
        # In this example a 3D array of a 2D nested array of 2D nested arrays
        # results in a 7D type.
        nptype = np.dtype((np.dtype((np.float64, (5, 2))), (3, 6)))
        nbtype = types.Array(numpy_support.from_dtype(nptype), 3, 'C')
        expected = types.Array(types.float64, 7, 'C')
        self.assertEqual(nbtype, expected)

    def test_issue_1469_3(self):
        # Nested arrays in record dtypes are accepted, but alignment is not
        # guaranteed.
        nptype = np.dtype([('a', np.float64,(4,))])
        nbtype = types.Array(numpy_support.from_dtype(nptype), 2, 'C')

        # Manual construction of expected array of record type
        natype = types.NestedArray(types.float64, (4,))
        fields = [('a', {'type': natype, 'offset': 0})]
        rectype = types.Record(fields=fields, size=32, aligned=False)
        expected = types.Array(rectype, 2, 'C', aligned=False)

        self.assertEqual(nbtype, expected)

    def test_issue_3158_1(self):
        # A nested array dtype.
        item = np.dtype([('some_field', np.int32)])
        items = np.dtype([('items', item, 3)])

        @njit
        def fn(x):
            return x[0]

        arr = np.asarray([([(0,), (1,), (2,)],),
                          ([(3,), (4,), (5,)],)],
                         dtype=items)

        expected = fn.py_func(arr)
        actual = fn(arr)

        self.assertEqual(expected, actual)

    def test_issue_3158_2(self):
        # A slightly more complex nested array dtype example.
        dtype1 = np.dtype([('a', 'i8'), ('b', 'i4')])
        dtype2 = np.dtype((dtype1, (2, 2)))
        dtype3 = np.dtype([('x', '?'), ('y', dtype2)])

        @njit
        def fn(arr):
            return arr[0]

        arr = np.asarray([(False, [[(0, 1), (2, 3)], [(4, 5), (6, 7)]]),
                          (True, [[(8, 9), (10, 11)], [(12, 13), (14, 15)]])],
                         dtype=dtype3)
        expected = fn.py_func(arr)
        actual = fn(arr)

        self.assertEqual(expected, actual)


if __name__ == '__main__':
    unittest.main()
