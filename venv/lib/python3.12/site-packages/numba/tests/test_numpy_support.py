"""
Test helper functions from numba.numpy_support.
"""


import sys
from itertools import product

import numpy as np

import unittest
from numba.core import types
from numba.core.errors import NumbaNotImplementedError
from numba.tests.support import TestCase
from numba.tests.enum_usecases import Shake, RequestError
from numba.np import numpy_support


class TestFromDtype(TestCase):

    def test_number_types(self):
        """
        Test from_dtype() and as_dtype() with the various scalar number types.
        """
        f = numpy_support.from_dtype

        def check(typechar, numba_type):
            # Only native ordering and alignment is supported
            dtype = np.dtype(typechar)
            self.assertIs(f(dtype), numba_type)
            self.assertIs(f(np.dtype('=' + typechar)), numba_type)
            self.assertEqual(dtype, numpy_support.as_dtype(numba_type))

        check('?', types.bool_)
        check('f', types.float32)
        check('f4', types.float32)
        check('d', types.float64)
        check('f8', types.float64)

        check('F', types.complex64)
        check('c8', types.complex64)
        check('D', types.complex128)
        check('c16', types.complex128)

        check('O', types.pyobject)

        check('b', types.int8)
        check('i1', types.int8)
        check('B', types.uint8)
        check('u1', types.uint8)

        check('h', types.int16)
        check('i2', types.int16)
        check('H', types.uint16)
        check('u2', types.uint16)

        check('i', types.int32)
        check('i4', types.int32)
        check('I', types.uint32)
        check('u4', types.uint32)

        check('q', types.int64)
        check('Q', types.uint64)
        for name in ('int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32',
                     'int64', 'uint64', 'intp', 'uintp'):
            self.assertIs(f(np.dtype(name)), getattr(types, name))

        # Non-native alignments are unsupported (except for 1-byte types)
        foreign_align = '>' if sys.byteorder == 'little' else '<'
        for letter in 'hHiIlLqQfdFD':
            self.assertRaises(NumbaNotImplementedError, f,
                              np.dtype(foreign_align + letter))

    def test_string_types(self):
        """
        Test from_dtype() and as_dtype() with the character string types.
        """
        def check(typestring, numba_type):
            # Only native ordering and alignment is supported
            dtype = np.dtype(typestring)
            self.assertEqual(numpy_support.from_dtype(dtype), numba_type)
            self.assertEqual(dtype, numpy_support.as_dtype(numba_type))

        check('S10', types.CharSeq(10))
        check('a11', types.CharSeq(11))
        check('U12', types.UnicodeCharSeq(12))

    def check_datetime_types(self, letter, nb_class):
        def check(dtype, numba_type, code):
            tp = numpy_support.from_dtype(dtype)
            self.assertEqual(tp, numba_type)
            self.assertEqual(tp.unit_code, code)
            self.assertEqual(numpy_support.as_dtype(numba_type), dtype)
            self.assertEqual(numpy_support.as_dtype(tp), dtype)

        # Unit-less ("generic") type
        check(np.dtype(letter), nb_class(''), 14)

    def test_datetime_types(self):
        """
        Test from_dtype() and as_dtype() with the datetime types.
        """
        self.check_datetime_types('M', types.NPDatetime)

    def test_timedelta_types(self):
        """
        Test from_dtype() and as_dtype() with the timedelta types.
        """
        self.check_datetime_types('m', types.NPTimedelta)

    def test_struct_types(self):
        def check(dtype, fields, size, aligned):
            tp = numpy_support.from_dtype(dtype)
            self.assertIsInstance(tp, types.Record)
            # Only check for dtype equality, as the Numba type may be interned
            self.assertEqual(tp.dtype, dtype)
            self.assertEqual(tp.fields, fields)
            self.assertEqual(tp.size, size)
            self.assertEqual(tp.aligned, aligned)

        dtype = np.dtype([('a', np.int16), ('b', np.int32)])
        check(dtype,
              fields={'a': (types.int16, 0, None, None),
                      'b': (types.int32, 2, None, None)},
              size=6, aligned=False)

        dtype = np.dtype([('a', np.int16), ('b', np.int32)], align=True)
        check(dtype,
              fields={'a': (types.int16, 0, None, None),
                      'b': (types.int32, 4, None, None)},
              size=8, aligned=True)

        dtype = np.dtype([('m', np.int32), ('n', 'S5')])
        check(dtype,
              fields={'m': (types.int32, 0, None, None),
                      'n': (types.CharSeq(5), 4, None, None)},
              size=9, aligned=False)

    def test_enum_type(self):

        def check(base_inst, enum_def, type_class):
            np_dt = np.dtype(base_inst)
            nb_ty = numpy_support.from_dtype(np_dt)
            inst = type_class(enum_def, nb_ty)
            recovered = numpy_support.as_dtype(inst)
            self.assertEqual(np_dt, recovered)

        dts = [np.float64, np.int32, np.complex128, np.bool_]
        enums = [Shake, RequestError]

        for dt, enum in product(dts, enums):
            check(dt, enum, types.EnumMember)

        for dt, enum in product(dts, enums):
            check(dt, enum, types.IntEnumMember)


class ValueTypingTestBase(object):
    """
    Common tests for the typing of values.  Also used by test_special.
    """

    def check_number_values(self, func):
        """
        Test *func*() with scalar numeric values.
        """
        f = func
        # Standard Python types get inferred by numpy
        self.assertIn(f(1), (types.int32, types.int64))
        self.assertIn(f(2**31 - 1), (types.int32, types.int64))
        self.assertIn(f(-2**31), (types.int32, types.int64))
        self.assertIs(f(1.0), types.float64)
        self.assertIs(f(1.0j), types.complex128)
        self.assertIs(f(True), types.bool_)
        self.assertIs(f(False), types.bool_)
        # Numpy scalar types get converted by from_dtype()
        for name in ('int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32',
                     'int64', 'uint64', 'intc', 'uintc', 'intp', 'uintp',
                     'float32', 'float64', 'complex64', 'complex128',
                     'bool_'):
            val = getattr(np, name)()
            self.assertIs(f(val), getattr(types, name))

    def _base_check_datetime_values(self, func, np_type, nb_type):
        f = func
        for unit in [
            '', 'Y', 'M', 'D', 'h', 'm', 's',
            'ms', 'us', 'ns', 'ps', 'fs', 'as',
        ]:
            if unit:
                t = np_type(3, unit)
            else:
                # "generic" datetime / timedelta
                t = np_type('Nat')
            tp = f(t)
            # This ensures the unit hasn't been lost
            self.assertEqual(tp, nb_type(unit))

    def check_datetime_values(self, func):
        """
        Test *func*() with np.datetime64 values.
        """
        self._base_check_datetime_values(func, np.datetime64, types.NPDatetime)

    def check_timedelta_values(self, func):
        """
        Test *func*() with np.timedelta64 values.
        """
        self._base_check_datetime_values(func, np.timedelta64,
                                         types.NPTimedelta)


class TestArrayScalars(ValueTypingTestBase, TestCase):

    def test_number_values(self):
        """
        Test map_arrayscalar_type() with scalar number values.
        """
        self.check_number_values(numpy_support.map_arrayscalar_type)

    def test_datetime_values(self):
        """
        Test map_arrayscalar_type() with np.datetime64 values.
        """
        f = numpy_support.map_arrayscalar_type
        self.check_datetime_values(f)
        # datetime64s with a non-one factor shouldn't be supported
        t = np.datetime64('2014', '10Y')
        with self.assertRaises(NumbaNotImplementedError):
            f(t)

    def test_timedelta_values(self):
        """
        Test map_arrayscalar_type() with np.timedelta64 values.
        """
        f = numpy_support.map_arrayscalar_type
        self.check_timedelta_values(f)
        # timedelta64s with a non-one factor shouldn't be supported
        t = np.timedelta64(10, '10Y')
        with self.assertRaises(NumbaNotImplementedError):
            f(t)


class FakeUFunc(object):
    __slots__ = ('nin', 'nout', 'types', 'ntypes')
    __name__ = "fake ufunc"

    def __init__(self, types):
        self.types = types
        in_, out = self.types[0].split('->')
        self.nin = len(in_)
        self.nout = len(out)
        self.ntypes = len(types)
        for tp in types:
            in_, out = self.types[0].split('->')
            assert len(in_) == self.nin
            assert len(out) == self.nout


# Typical types for np.add, np.multiply, np.isnan
_add_types = ['??->?', 'bb->b', 'BB->B', 'hh->h', 'HH->H', 'ii->i', 'II->I',
              'll->l', 'LL->L', 'qq->q', 'QQ->Q', 'ee->e', 'ff->f', 'dd->d',
              'gg->g', 'FF->F', 'DD->D', 'GG->G', 'Mm->M', 'mm->m', 'mM->M',
              'OO->O']

_mul_types = ['??->?', 'bb->b', 'BB->B', 'hh->h', 'HH->H', 'ii->i', 'II->I',
              'll->l', 'LL->L', 'qq->q', 'QQ->Q', 'ee->e', 'ff->f', 'dd->d',
              'gg->g', 'FF->F', 'DD->D', 'GG->G', 'mq->m', 'qm->m', 'md->m',
              'dm->m', 'OO->O']

# Those ones only have floating-point loops
_isnan_types = ['e->?', 'f->?', 'd->?', 'g->?', 'F->?', 'D->?', 'G->?']
_sqrt_types = ['e->e', 'f->f', 'd->d', 'g->g', 'F->F', 'D->D', 'G->G', 'O->O']


class TestUFuncs(TestCase):
    """
    Test ufunc helpers.
    """

    def test_ufunc_find_matching_loop(self):
        f = numpy_support.ufunc_find_matching_loop
        np_add = FakeUFunc(_add_types)
        np_mul = FakeUFunc(_mul_types)
        np_isnan = FakeUFunc(_isnan_types)
        np_sqrt = FakeUFunc(_sqrt_types)

        def check(ufunc, input_types, sigs, output_types=()):
            """
            Check that ufunc_find_matching_loop() finds one of the given
            *sigs* for *ufunc*, *input_types* and optional *output_types*.
            """
            loop = f(ufunc, input_types + output_types)
            self.assertTrue(loop)
            if isinstance(sigs, str):
                sigs = (sigs,)
            self.assertIn(loop.ufunc_sig, sigs,
                          "inputs=%s and outputs=%s should have selected "
                          "one of %s, got %s"
                          % (input_types, output_types, sigs, loop.ufunc_sig))
            self.assertEqual(len(loop.numpy_inputs), len(loop.inputs))
            self.assertEqual(len(loop.numpy_outputs), len(loop.outputs))
            if not output_types:
                # Add explicit outputs and check the result is the same
                loop_explicit = f(ufunc, list(input_types) + loop.outputs)
                self.assertEqual(loop_explicit, loop)
            else:
                self.assertEqual(loop.outputs, list(output_types))
            # Round-tripping inputs and outputs
            loop_rt = f(ufunc, loop.inputs + loop.outputs)
            self.assertEqual(loop_rt, loop)
            return loop

        def check_exact(ufunc, input_types, sigs, output_types=()):
            """
            Like check(), but also ensure no casting of inputs occurred.
            """
            loop = check(ufunc, input_types, sigs, output_types)
            self.assertEqual(loop.inputs, list(input_types))

        def check_no_match(ufunc, input_types):
            loop = f(ufunc, input_types)
            self.assertIs(loop, None)

        # Exact matching for number types
        check_exact(np_add, (types.bool_, types.bool_), '??->?')
        check_exact(np_add, (types.int8, types.int8), 'bb->b')
        check_exact(np_add, (types.uint8, types.uint8), 'BB->B')
        check_exact(np_add, (types.int64, types.int64), ('ll->l', 'qq->q'))
        check_exact(np_add, (types.uint64, types.uint64), ('LL->L', 'QQ->Q'))
        check_exact(np_add, (types.float32, types.float32), 'ff->f')
        check_exact(np_add, (types.float64, types.float64), 'dd->d')
        check_exact(np_add, (types.complex64, types.complex64), 'FF->F')
        check_exact(np_add, (types.complex128, types.complex128), 'DD->D')

        # Exact matching for datetime64 and timedelta64 types
        check_exact(np_add, (types.NPTimedelta('s'), types.NPTimedelta('s')),
                    'mm->m', output_types=(types.NPTimedelta('s'),))
        check_exact(np_add, (types.NPTimedelta('ms'), types.NPDatetime('s')),
                    'mM->M', output_types=(types.NPDatetime('ms'),))
        check_exact(np_add, (types.NPDatetime('s'), types.NPTimedelta('s')),
                    'Mm->M', output_types=(types.NPDatetime('s'),))
        check_exact(np_add, (types.NPDatetime('s'), types.NPTimedelta('')),
                    'Mm->M', output_types=(types.NPDatetime('s'),))
        check_exact(np_add, (types.NPDatetime('ns'), types.NPTimedelta('')),
                    'Mm->M', output_types=(types.NPDatetime('ns'),))
        check_exact(np_add, (types.NPTimedelta(''), types.NPDatetime('s')),
                    'mM->M', output_types=(types.NPDatetime('s'),))
        check_exact(np_add, (types.NPTimedelta(''), types.NPDatetime('ns')),
                    'mM->M', output_types=(types.NPDatetime('ns'),))
        check_exact(np_mul, (types.NPTimedelta('s'), types.int64),
                    'mq->m', output_types=(types.NPTimedelta('s'),))
        check_exact(np_mul, (types.float64, types.NPTimedelta('s')),
                    'dm->m', output_types=(types.NPTimedelta('s'),))

        # Mix and match number types, with casting
        check(np_add, (types.bool_, types.int8), 'bb->b')
        check(np_add, (types.uint8, types.bool_), 'BB->B')
        check(np_add, (types.int16, types.uint16), 'ii->i')
        check(np_add, (types.complex64, types.float64), 'DD->D')
        check(np_add, (types.float64, types.complex64), 'DD->D')
        # Integers, when used together with floating-point numbers,
        # should cast to any real or complex (see #2006)
        int_types = [types.int32, types.uint32, types.int64, types.uint64]
        for intty in int_types:
            check(np_add, (types.float32, intty), 'ff->f')
            check(np_add, (types.float64, intty), 'dd->d')
            check(np_add, (types.complex64, intty), 'FF->F')
            check(np_add, (types.complex128, intty), 'DD->D')
        # However, when used alone, they should cast only to
        # floating-point types of sufficient precision
        # (typical use case: np.sqrt(2) should give an accurate enough value)
        for intty in int_types:
            check(np_sqrt, (intty,), 'd->d')
            check(np_isnan, (intty,), 'd->?')

        # With some timedelta64 arguments as well
        check(np_mul, (types.NPTimedelta('s'), types.int32),
              'mq->m', output_types=(types.NPTimedelta('s'),))
        check(np_mul, (types.NPTimedelta('s'), types.uint32),
              'mq->m', output_types=(types.NPTimedelta('s'),))
        check(np_mul, (types.NPTimedelta('s'), types.float32),
              'md->m', output_types=(types.NPTimedelta('s'),))
        check(np_mul, (types.float32, types.NPTimedelta('s')),
              'dm->m', output_types=(types.NPTimedelta('s'),))

        # No match
        check_no_match(np_add, (types.NPDatetime('s'), types.NPDatetime('s')))
        # No implicit casting from int64 to timedelta64 (Numpy would allow
        # this).
        check_no_match(np_add, (types.NPTimedelta('s'), types.int64))

    def test_layout_checker(self):
        def check_arr(arr):
            dims = arr.shape
            strides = arr.strides
            itemsize = arr.dtype.itemsize
            is_c = numpy_support.is_contiguous(dims, strides, itemsize)
            is_f = numpy_support.is_fortran(dims, strides, itemsize)
            expect_c = arr.flags['C_CONTIGUOUS']
            expect_f = arr.flags['F_CONTIGUOUS']
            self.assertEqual(is_c, expect_c)
            self.assertEqual(is_f, expect_f)

        arr = np.arange(24)
        # 1D
        check_arr(arr)
        # 2D
        check_arr(arr.reshape((3, 8)))
        check_arr(arr.reshape((3, 8)).T)
        check_arr(arr.reshape((3, 8))[::2])
        # 3D
        check_arr(arr.reshape((2, 3, 4)))
        check_arr(arr.reshape((2, 3, 4)).T)
        # middle axis is shape 1
        check_arr(arr.reshape((2, 3, 4))[:, ::3])
        check_arr(arr.reshape((2, 3, 4)).T[:, ::3])

        # leading axis is shape 1
        check_arr(arr.reshape((2, 3, 4))[::2])
        check_arr(arr.reshape((2, 3, 4)).T[:, :, ::2])
        # 2 leading axis are shape 1
        check_arr(arr.reshape((2, 3, 4))[::2, ::3])
        check_arr(arr.reshape((2, 3, 4)).T[:, ::3, ::2])
        # single item slices for all axis
        check_arr(arr.reshape((2, 3, 4))[::2, ::3, ::4])
        check_arr(arr.reshape((2, 3, 4)).T[::4, ::3, ::2])
        # 4D
        check_arr(arr.reshape((2, 2, 3, 2))[::2, ::2, ::3])
        check_arr(arr.reshape((2, 2, 3, 2)).T[:, ::3, ::2, ::2])
        # outer zero dims
        check_arr(arr.reshape((2, 2, 3, 2))[::5, ::2, ::3])
        check_arr(arr.reshape((2, 2, 3, 2)).T[:, ::3, ::2, ::5])


if __name__ == '__main__':
    unittest.main()
