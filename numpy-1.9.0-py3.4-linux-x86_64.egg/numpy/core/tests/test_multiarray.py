from __future__ import division, absolute_import, print_function

import tempfile
import sys
import os
import shutil
import warnings
import operator
import io
if sys.version_info[0] >= 3:
    import builtins
else:
    import __builtin__ as builtins
from decimal import Decimal


import numpy as np
from nose import SkipTest
from numpy.core import *
from numpy.compat import asbytes, getexception, strchar, sixu
from test_print import in_foreign_locale
from numpy.core.multiarray_tests import (
        test_neighborhood_iterator, test_neighborhood_iterator_oob,
        test_pydatamem_seteventhook_start, test_pydatamem_seteventhook_end,
        test_inplace_increment, get_buffer_info
        )
from numpy.testing import (
        TestCase, run_module_suite, assert_, assert_raises,
        assert_equal, assert_almost_equal, assert_array_equal,
        assert_array_almost_equal, assert_allclose,
        assert_array_less, runstring, dec
        )

# Need to test an object that does not fully implement math interface
from datetime import timedelta


if sys.version_info[:2] > (3, 2):
    # In Python 3.3 the representation of empty shape, strides and suboffsets
    # is an empty tuple instead of None.
    # http://docs.python.org/dev/whatsnew/3.3.html#api-changes
    EMPTY = ()
else:
    EMPTY = None


class TestFlags(TestCase):
    def setUp(self):
        self.a = arange(10)

    def test_writeable(self):
        mydict = locals()
        self.a.flags.writeable = False
        self.assertRaises(ValueError, runstring, 'self.a[0] = 3', mydict)
        self.assertRaises(ValueError, runstring, 'self.a[0:1].itemset(3)', mydict)
        self.a.flags.writeable = True
        self.a[0] = 5
        self.a[0] = 0

    def test_otherflags(self):
        assert_equal(self.a.flags.carray, True)
        assert_equal(self.a.flags.farray, False)
        assert_equal(self.a.flags.behaved, True)
        assert_equal(self.a.flags.fnc, False)
        assert_equal(self.a.flags.forc, True)
        assert_equal(self.a.flags.owndata, True)
        assert_equal(self.a.flags.writeable, True)
        assert_equal(self.a.flags.aligned, True)
        assert_equal(self.a.flags.updateifcopy, False)

class TestHash(TestCase):
    # see #3793
    def test_int(self):
        for st, ut, s in [(np.int8, np.uint8, 8),
                          (np.int16, np.uint16, 16),
                          (np.int32, np.uint32, 32),
                          (np.int64, np.uint64, 64)]:
            for i in range(1, s):
                assert_equal(hash(st(-2**i)), hash(-2**i),
                             err_msg="%r: -2**%d" % (st, i))
                assert_equal(hash(st(2**(i - 1))), hash(2**(i - 1)),
                             err_msg="%r: 2**%d" % (st, i - 1))
                assert_equal(hash(st(2**i - 1)), hash(2**i - 1),
                             err_msg="%r: 2**%d - 1" % (st, i))

                i = max(i - 1, 1)
                assert_equal(hash(ut(2**(i - 1))), hash(2**(i - 1)),
                             err_msg="%r: 2**%d" % (ut, i - 1))
                assert_equal(hash(ut(2**i - 1)), hash(2**i - 1),
                             err_msg="%r: 2**%d - 1" % (ut, i))

class TestAttributes(TestCase):
    def setUp(self):
        self.one = arange(10)
        self.two = arange(20).reshape(4, 5)
        self.three = arange(60, dtype=float64).reshape(2, 5, 6)

    def test_attributes(self):
        assert_equal(self.one.shape, (10,))
        assert_equal(self.two.shape, (4, 5))
        assert_equal(self.three.shape, (2, 5, 6))
        self.three.shape = (10, 3, 2)
        assert_equal(self.three.shape, (10, 3, 2))
        self.three.shape = (2, 5, 6)
        assert_equal(self.one.strides, (self.one.itemsize,))
        num = self.two.itemsize
        assert_equal(self.two.strides, (5*num, num))
        num = self.three.itemsize
        assert_equal(self.three.strides, (30*num, 6*num, num))
        assert_equal(self.one.ndim, 1)
        assert_equal(self.two.ndim, 2)
        assert_equal(self.three.ndim, 3)
        num = self.two.itemsize
        assert_equal(self.two.size, 20)
        assert_equal(self.two.nbytes, 20*num)
        assert_equal(self.two.itemsize, self.two.dtype.itemsize)
        assert_equal(self.two.base, arange(20))

    def test_dtypeattr(self):
        assert_equal(self.one.dtype, dtype(int_))
        assert_equal(self.three.dtype, dtype(float_))
        assert_equal(self.one.dtype.char, 'l')
        assert_equal(self.three.dtype.char, 'd')
        self.assertTrue(self.three.dtype.str[0] in '<>')
        assert_equal(self.one.dtype.str[1], 'i')
        assert_equal(self.three.dtype.str[1], 'f')

    def test_int_subclassing(self):
        # Regression test for https://github.com/numpy/numpy/pull/3526

        numpy_int = np.int_(0)

        if sys.version_info[0] >= 3:
            # On Py3k int_ should not inherit from int, because it's not fixed-width anymore
            assert_equal(isinstance(numpy_int, int), False)
        else:
            # Otherwise, it should inherit from int...
            assert_equal(isinstance(numpy_int, int), True)

            # ... and fast-path checks on C-API level should also work
            from numpy.core.multiarray_tests import test_int_subclass
            assert_equal(test_int_subclass(numpy_int), True)

    def test_stridesattr(self):
        x = self.one
        def make_array(size, offset, strides):
            return ndarray(size, buffer=x, dtype=int,
                           offset=offset*x.itemsize,
                           strides=strides*x.itemsize)
        assert_equal(make_array(4, 4, -1), array([4, 3, 2, 1]))
        self.assertRaises(ValueError, make_array, 4, 4, -2)
        self.assertRaises(ValueError, make_array, 4, 2, -1)
        self.assertRaises(ValueError, make_array, 8, 3, 1)
        assert_equal(make_array(8, 3, 0), np.array([3]*8))
        # Check behavior reported in gh-2503:
        self.assertRaises(ValueError, make_array, (2, 3), 5, array([-2, -3]))
        make_array(0, 0, 10)

    def test_set_stridesattr(self):
        x = self.one
        def make_array(size, offset, strides):
            try:
                r = ndarray([size], dtype=int, buffer=x, offset=offset*x.itemsize)
            except:
                raise RuntimeError(getexception())
            r.strides = strides=strides*x.itemsize
            return r
        assert_equal(make_array(4, 4, -1), array([4, 3, 2, 1]))
        assert_equal(make_array(7, 3, 1), array([3, 4, 5, 6, 7, 8, 9]))
        self.assertRaises(ValueError, make_array, 4, 4, -2)
        self.assertRaises(ValueError, make_array, 4, 2, -1)
        self.assertRaises(RuntimeError, make_array, 8, 3, 1)
        # Check that the true extent of the array is used.
        # Test relies on as_strided base not exposing a buffer.
        x = np.lib.stride_tricks.as_strided(arange(1), (10, 10), (0, 0))
        def set_strides(arr, strides):
            arr.strides = strides
        self.assertRaises(ValueError, set_strides, x, (10*x.itemsize, x.itemsize))

        # Test for offset calculations:
        x = np.lib.stride_tricks.as_strided(np.arange(10, dtype=np.int8)[-1],
                                                    shape=(10,), strides=(-1,))
        self.assertRaises(ValueError, set_strides, x[::-1], -1)
        a = x[::-1]
        a.strides = 1
        a[::2].strides = 2

    def test_fill(self):
        for t in "?bhilqpBHILQPfdgFDGO":
            x = empty((3, 2, 1), t)
            y = empty((3, 2, 1), t)
            x.fill(1)
            y[...] = 1
            assert_equal(x, y)

    def test_fill_struct_array(self):
        # Filling from a scalar
        x = array([(0, 0.0), (1, 1.0)], dtype='i4,f8')
        x.fill(x[0])
        assert_equal(x['f1'][1], x['f1'][0])
        # Filling from a tuple that can be converted
        # to a scalar
        x = np.zeros(2, dtype=[('a', 'f8'), ('b', 'i4')])
        x.fill((3.5, -2))
        assert_array_equal(x['a'], [3.5, 3.5])
        assert_array_equal(x['b'], [-2, -2])


class TestArrayConstruction(TestCase):
    def test_array(self):
        d = np.ones(6)
        r = np.array([d, d])
        assert_equal(r, np.ones((2, 6)))

        d = np.ones(6)
        tgt = np.ones((2, 6))
        r = np.array([d, d])
        assert_equal(r, tgt)
        tgt[1] = 2
        r = np.array([d, d + 1])
        assert_equal(r, tgt)

        d = np.ones(6)
        r = np.array([[d, d]])
        assert_equal(r, np.ones((1, 2, 6)))

        d = np.ones(6)
        r = np.array([[d, d], [d, d]])
        assert_equal(r, np.ones((2, 2, 6)))

        d = np.ones((6, 6))
        r = np.array([d, d])
        assert_equal(r, np.ones((2, 6, 6)))

        d = np.ones((6, ))
        r = np.array([[d, d + 1], d + 2])
        assert_equal(len(r), 2)
        assert_equal(r[0], [d, d + 1])
        assert_equal(r[1], d + 2)

        tgt = np.ones((2, 3), dtype=np.bool)
        tgt[0, 2] = False
        tgt[1, 0:2] = False
        r = np.array([[True, True, False], [False, False, True]])
        assert_equal(r, tgt)
        r = np.array([[True, False], [True, False], [False, True]])
        assert_equal(r, tgt.T)


class TestAssignment(TestCase):
    def test_assignment_broadcasting(self):
        a = np.arange(6).reshape(2, 3)

        # Broadcasting the input to the output
        a[...] = np.arange(3)
        assert_equal(a, [[0, 1, 2], [0, 1, 2]])
        a[...] = np.arange(2).reshape(2, 1)
        assert_equal(a, [[0, 0, 0], [1, 1, 1]])

        # For compatibility with <= 1.5, a limited version of broadcasting
        # the output to the input.
        #
        # This behavior is inconsistent with NumPy broadcasting
        # in general, because it only uses one of the two broadcasting
        # rules (adding a new "1" dimension to the left of the shape),
        # applied to the output instead of an input. In NumPy 2.0, this kind
        # of broadcasting assignment will likely be disallowed.
        a[...] = np.arange(6)[::-1].reshape(1, 2, 3)
        assert_equal(a, [[5, 4, 3], [2, 1, 0]])
        # The other type of broadcasting would require a reduction operation.
        def assign(a, b):
            a[...] = b
        assert_raises(ValueError, assign, a, np.arange(12).reshape(2, 2, 3))

    def test_assignment_errors(self):
        # Address issue #2276
        class C:
            pass
        a = np.zeros(1)
        def assign(v):
            a[0] = v
        assert_raises((AttributeError, TypeError), assign, C())
        assert_raises(ValueError, assign, [1])

class TestDtypedescr(TestCase):
    def test_construction(self):
        d1 = dtype('i4')
        assert_equal(d1, dtype(int32))
        d2 = dtype('f8')
        assert_equal(d2, dtype(float64))

class TestZeroRank(TestCase):
    def setUp(self):
        self.d = array(0), array('x', object)

    def test_ellipsis_subscript(self):
        a, b = self.d
        self.assertEqual(a[...], 0)
        self.assertEqual(b[...], 'x')
        self.assertTrue(a[...].base is a) # `a[...] is a` in numpy <1.9.
        self.assertTrue(b[...].base is b) # `b[...] is b` in numpy <1.9.

    def test_empty_subscript(self):
        a, b = self.d
        self.assertEqual(a[()], 0)
        self.assertEqual(b[()], 'x')
        self.assertTrue(type(a[()]) is a.dtype.type)
        self.assertTrue(type(b[()]) is str)

    def test_invalid_subscript(self):
        a, b = self.d
        self.assertRaises(IndexError, lambda x: x[0], a)
        self.assertRaises(IndexError, lambda x: x[0], b)
        self.assertRaises(IndexError, lambda x: x[array([], int)], a)
        self.assertRaises(IndexError, lambda x: x[array([], int)], b)

    def test_ellipsis_subscript_assignment(self):
        a, b = self.d
        a[...] = 42
        self.assertEqual(a, 42)
        b[...] = ''
        self.assertEqual(b.item(), '')

    def test_empty_subscript_assignment(self):
        a, b = self.d
        a[()] = 42
        self.assertEqual(a, 42)
        b[()] = ''
        self.assertEqual(b.item(), '')

    def test_invalid_subscript_assignment(self):
        a, b = self.d
        def assign(x, i, v):
            x[i] = v
        self.assertRaises(IndexError, assign, a, 0, 42)
        self.assertRaises(IndexError, assign, b, 0, '')
        self.assertRaises(ValueError, assign, a, (), '')

    def test_newaxis(self):
        a, b = self.d
        self.assertEqual(a[newaxis].shape, (1,))
        self.assertEqual(a[..., newaxis].shape, (1,))
        self.assertEqual(a[newaxis, ...].shape, (1,))
        self.assertEqual(a[..., newaxis].shape, (1,))
        self.assertEqual(a[newaxis, ..., newaxis].shape, (1, 1))
        self.assertEqual(a[..., newaxis, newaxis].shape, (1, 1))
        self.assertEqual(a[newaxis, newaxis, ...].shape, (1, 1))
        self.assertEqual(a[(newaxis,)*10].shape, (1,)*10)

    def test_invalid_newaxis(self):
        a, b = self.d
        def subscript(x, i): x[i]
        self.assertRaises(IndexError, subscript, a, (newaxis, 0))
        self.assertRaises(IndexError, subscript, a, (newaxis,)*50)

    def test_constructor(self):
        x = ndarray(())
        x[()] = 5
        self.assertEqual(x[()], 5)
        y = ndarray((), buffer=x)
        y[()] = 6
        self.assertEqual(x[()], 6)

    def test_output(self):
        x = array(2)
        self.assertRaises(ValueError, add, x, [1], x)


class TestScalarIndexing(TestCase):
    def setUp(self):
        self.d = array([0, 1])[0]

    def test_ellipsis_subscript(self):
        a = self.d
        self.assertEqual(a[...], 0)
        self.assertEqual(a[...].shape, ())

    def test_empty_subscript(self):
        a = self.d
        self.assertEqual(a[()], 0)
        self.assertEqual(a[()].shape, ())

    def test_invalid_subscript(self):
        a = self.d
        self.assertRaises(IndexError, lambda x: x[0], a)
        self.assertRaises(IndexError, lambda x: x[array([], int)], a)

    def test_invalid_subscript_assignment(self):
        a = self.d
        def assign(x, i, v):
            x[i] = v
        self.assertRaises(TypeError, assign, a, 0, 42)

    def test_newaxis(self):
        a = self.d
        self.assertEqual(a[newaxis].shape, (1,))
        self.assertEqual(a[..., newaxis].shape, (1,))
        self.assertEqual(a[newaxis, ...].shape, (1,))
        self.assertEqual(a[..., newaxis].shape, (1,))
        self.assertEqual(a[newaxis, ..., newaxis].shape, (1, 1))
        self.assertEqual(a[..., newaxis, newaxis].shape, (1, 1))
        self.assertEqual(a[newaxis, newaxis, ...].shape, (1, 1))
        self.assertEqual(a[(newaxis,)*10].shape, (1,)*10)

    def test_invalid_newaxis(self):
        a = self.d
        def subscript(x, i): x[i]
        self.assertRaises(IndexError, subscript, a, (newaxis, 0))
        self.assertRaises(IndexError, subscript, a, (newaxis,)*50)

    def test_overlapping_assignment(self):
        # With positive strides
        a = np.arange(4)
        a[:-1] = a[1:]
        assert_equal(a, [1, 2, 3, 3])

        a = np.arange(4)
        a[1:] = a[:-1]
        assert_equal(a, [0, 0, 1, 2])

        # With positive and negative strides
        a = np.arange(4)
        a[:] = a[::-1]
        assert_equal(a, [3, 2, 1, 0])

        a = np.arange(6).reshape(2, 3)
        a[::-1,:] = a[:, ::-1]
        assert_equal(a, [[5, 4, 3], [2, 1, 0]])

        a = np.arange(6).reshape(2, 3)
        a[::-1, ::-1] = a[:, ::-1]
        assert_equal(a, [[3, 4, 5], [0, 1, 2]])

        # With just one element overlapping
        a = np.arange(5)
        a[:3] = a[2:]
        assert_equal(a, [2, 3, 4, 3, 4])

        a = np.arange(5)
        a[2:] = a[:3]
        assert_equal(a, [0, 1, 0, 1, 2])

        a = np.arange(5)
        a[2::-1] = a[2:]
        assert_equal(a, [4, 3, 2, 3, 4])

        a = np.arange(5)
        a[2:] = a[2::-1]
        assert_equal(a, [0, 1, 2, 1, 0])

        a = np.arange(5)
        a[2::-1] = a[:1:-1]
        assert_equal(a, [2, 3, 4, 3, 4])

        a = np.arange(5)
        a[:1:-1] = a[2::-1]
        assert_equal(a, [0, 1, 0, 1, 2])

class TestCreation(TestCase):
    def test_from_attribute(self):
        class x(object):
            def __array__(self, dtype=None):
                pass
        self.assertRaises(ValueError, array, x())

    def test_from_string(self) :
        types = np.typecodes['AllInteger'] + np.typecodes['Float']
        nstr = ['123', '123']
        result = array([123, 123], dtype=int)
        for type in types :
            msg = 'String conversion for %s' % type
            assert_equal(array(nstr, dtype=type), result, err_msg=msg)

    def test_void(self):
        arr = np.array([], dtype='V')
        assert_equal(arr.dtype.kind, 'V')

    def test_zeros(self):
        types = np.typecodes['AllInteger'] + np.typecodes['AllFloat']
        for dt in types:
            d = np.zeros((13,), dtype=dt)
            assert_equal(np.count_nonzero(d), 0)
            # true for ieee floats
            assert_equal(d.sum(), 0)
            assert_(not d.any())

            d = np.zeros(2, dtype='(2,4)i4')
            assert_equal(np.count_nonzero(d), 0)
            assert_equal(d.sum(), 0)
            assert_(not d.any())

            d = np.zeros(2, dtype='4i4')
            assert_equal(np.count_nonzero(d), 0)
            assert_equal(d.sum(), 0)
            assert_(not d.any())

            d = np.zeros(2, dtype='(2,4)i4, (2,4)i4')
            assert_equal(np.count_nonzero(d), 0)

    @dec.slow
    def test_zeros_big(self):
        # test big array as they might be allocated different by the sytem
        types = np.typecodes['AllInteger'] + np.typecodes['AllFloat']
        for dt in types:
            d = np.zeros((30 * 1024**2,), dtype=dt)
            assert_(not d.any())

    def test_zeros_obj(self):
        # test initialization from PyLong(0)
        d = np.zeros((13,), dtype=object)
        assert_array_equal(d, [0] * 13)
        assert_equal(np.count_nonzero(d), 0)

    def test_zeros_obj_obj(self):
        d = zeros(10, dtype=[('k', object, 2)])
        assert_array_equal(d['k'], 0)

    def test_zeros_like_like_zeros(self):
        # test zeros_like returns the same as zeros
        for c in np.typecodes['All']:
            if c == 'V':
                continue
            d = zeros((3,3), dtype=c)
            assert_array_equal(zeros_like(d), d)
            assert_equal(zeros_like(d).dtype, d.dtype)
        # explicitly check some special cases
        d = zeros((3,3), dtype='S5')
        assert_array_equal(zeros_like(d), d)
        assert_equal(zeros_like(d).dtype, d.dtype)
        d = zeros((3,3), dtype='U5')
        assert_array_equal(zeros_like(d), d)
        assert_equal(zeros_like(d).dtype, d.dtype)

        d = zeros((3,3), dtype='<i4')
        assert_array_equal(zeros_like(d), d)
        assert_equal(zeros_like(d).dtype, d.dtype)
        d = zeros((3,3), dtype='>i4')
        assert_array_equal(zeros_like(d), d)
        assert_equal(zeros_like(d).dtype, d.dtype)

        d = zeros((3,3), dtype='<M8[s]')
        assert_array_equal(zeros_like(d), d)
        assert_equal(zeros_like(d).dtype, d.dtype)
        d = zeros((3,3), dtype='>M8[s]')
        assert_array_equal(zeros_like(d), d)
        assert_equal(zeros_like(d).dtype, d.dtype)

        d = zeros((3,3), dtype='f4,f4')
        assert_array_equal(zeros_like(d), d)
        assert_equal(zeros_like(d).dtype, d.dtype)

    def test_empty_unicode(self):
        # don't throw decode errors on garbage memory
        for i in range(5, 100, 5):
            d = np.empty(i, dtype='U')
            str(d)

    def test_sequence_non_homogenous(self):
        assert_equal(np.array([4, 2**80]).dtype, np.object)
        assert_equal(np.array([4, 2**80, 4]).dtype, np.object)
        assert_equal(np.array([2**80, 4]).dtype, np.object)
        assert_equal(np.array([2**80] * 3).dtype, np.object)
        assert_equal(np.array([[1, 1],[1j, 1j]]).dtype, np.complex)
        assert_equal(np.array([[1j, 1j],[1, 1]]).dtype, np.complex)
        assert_equal(np.array([[1, 1, 1],[1, 1j, 1.], [1, 1, 1]]).dtype, np.complex)

    @dec.skipif(sys.version_info[0] >= 3)
    def test_sequence_long(self):
        assert_equal(np.array([long(4), long(4)]).dtype, np.long)
        assert_equal(np.array([long(4), 2**80]).dtype, np.object)
        assert_equal(np.array([long(4), 2**80, long(4)]).dtype, np.object)
        assert_equal(np.array([2**80, long(4)]).dtype, np.object)

    def test_non_sequence_sequence(self):
        """Should not segfault.

        Class Fail breaks the sequence protocol for new style classes, i.e.,
        those derived from object. Class Map is a mapping type indicated by
        raising a ValueError. At some point we may raise a warning instead
        of an error in the Fail case.

        """
        class Fail(object):
            def __len__(self):
                return 1

            def __getitem__(self, index):
                raise ValueError()

        class Map(object):
            def __len__(self):
                return 1

            def __getitem__(self, index):
                raise KeyError()

        a = np.array([Map()])
        assert_(a.shape == (1,))
        assert_(a.dtype == np.dtype(object))
        assert_raises(ValueError, np.array, [Fail()])


class TestStructured(TestCase):
    def test_subarray_field_access(self):
        a = np.zeros((3, 5), dtype=[('a', ('i4', (2, 2)))])
        a['a'] = np.arange(60).reshape(3, 5, 2, 2)

        # Since the subarray is always in C-order, a transpose
        # does not swap the subarray:
        assert_array_equal(a.T['a'], a['a'].transpose(1, 0, 2, 3))

        # In Fortran order, the subarray gets appended
        # like in all other cases, not prepended as a special case
        b = a.copy(order='F')
        assert_equal(a['a'].shape, b['a'].shape)
        assert_equal(a.T['a'].shape, a.T.copy()['a'].shape)

    def test_subarray_comparison(self):
        # Check that comparisons between record arrays with
        # multi-dimensional field types work properly
        a = np.rec.fromrecords(
            [([1, 2, 3], 'a', [[1, 2], [3, 4]]), ([3, 3, 3], 'b', [[0, 0], [0, 0]])],
            dtype=[('a', ('f4', 3)), ('b', np.object), ('c', ('i4', (2, 2)))])
        b = a.copy()
        assert_equal(a==b, [True, True])
        assert_equal(a!=b, [False, False])
        b[1].b = 'c'
        assert_equal(a==b, [True, False])
        assert_equal(a!=b, [False, True])
        for i in range(3):
            b[0].a = a[0].a
            b[0].a[i] = 5
            assert_equal(a==b, [False, False])
            assert_equal(a!=b, [True, True])
        for i in range(2):
            for j in range(2):
                b = a.copy()
                b[0].c[i, j] = 10
                assert_equal(a==b, [False, True])
                assert_equal(a!=b, [True, False])

        # Check that broadcasting with a subarray works
        a = np.array([[(0,)], [(1,)]], dtype=[('a', 'f8')])
        b = np.array([(0,), (0,), (1,)], dtype=[('a', 'f8')])
        assert_equal(a==b, [[True, True, False], [False, False, True]])
        assert_equal(b==a, [[True, True, False], [False, False, True]])
        a = np.array([[(0,)], [(1,)]], dtype=[('a', 'f8', (1,))])
        b = np.array([(0,), (0,), (1,)], dtype=[('a', 'f8', (1,))])
        assert_equal(a==b, [[True, True, False], [False, False, True]])
        assert_equal(b==a, [[True, True, False], [False, False, True]])
        a = np.array([[([0, 0],)], [([1, 1],)]], dtype=[('a', 'f8', (2,))])
        b = np.array([([0, 0],), ([0, 1],), ([1, 1],)], dtype=[('a', 'f8', (2,))])
        assert_equal(a==b, [[True, False, False], [False, False, True]])
        assert_equal(b==a, [[True, False, False], [False, False, True]])

        # Check that broadcasting Fortran-style arrays with a subarray work
        a = np.array([[([0, 0],)], [([1, 1],)]], dtype=[('a', 'f8', (2,))], order='F')
        b = np.array([([0, 0],), ([0, 1],), ([1, 1],)], dtype=[('a', 'f8', (2,))])
        assert_equal(a==b, [[True, False, False], [False, False, True]])
        assert_equal(b==a, [[True, False, False], [False, False, True]])

        # Check that incompatible sub-array shapes don't result to broadcasting
        x = np.zeros((1,), dtype=[('a', ('f4', (1, 2))), ('b', 'i1')])
        y = np.zeros((1,), dtype=[('a', ('f4', (2,))), ('b', 'i1')])
        assert_equal(x == y, False)

        x = np.zeros((1,), dtype=[('a', ('f4', (2, 1))), ('b', 'i1')])
        y = np.zeros((1,), dtype=[('a', ('f4', (2,))), ('b', 'i1')])
        assert_equal(x == y, False)

        # Check that structured arrays that are different only in
        # byte-order work
        a = np.array([(5, 42), (10, 1)], dtype=[('a', '>i8'), ('b', '<f8')])
        b = np.array([(5, 43), (10, 1)], dtype=[('a', '<i8'), ('b', '>f8')])
        assert_equal(a == b, [False, True])


class TestBool(TestCase):
    def test_test_interning(self):
        a0 = bool_(0)
        b0 = bool_(False)
        self.assertTrue(a0 is b0)
        a1 = bool_(1)
        b1 = bool_(True)
        self.assertTrue(a1 is b1)
        self.assertTrue(array([True])[0] is a1)
        self.assertTrue(array(True)[()] is a1)

    def test_sum(self):
        d = np.ones(101, dtype=np.bool);
        assert_equal(d.sum(), d.size)
        assert_equal(d[::2].sum(), d[::2].size)
        assert_equal(d[::-2].sum(), d[::-2].size)

        d = np.frombuffer(b'\xff\xff' * 100, dtype=bool)
        assert_equal(d.sum(), d.size)
        assert_equal(d[::2].sum(), d[::2].size)
        assert_equal(d[::-2].sum(), d[::-2].size)

    def check_count_nonzero(self, power, length):
        powers = [2 ** i for i in range(length)]
        for i in range(2**power):
            l = [(i & x) != 0 for x in powers]
            a = np.array(l, dtype=np.bool)
            c = builtins.sum(l)
            self.assertEqual(np.count_nonzero(a), c)
            av = a.view(np.uint8)
            av *= 3
            self.assertEqual(np.count_nonzero(a), c)
            av *= 4
            self.assertEqual(np.count_nonzero(a), c)
            av[av != 0] = 0xFF
            self.assertEqual(np.count_nonzero(a), c)

    def test_count_nonzero(self):
        # check all 12 bit combinations in a length 17 array
        # covers most cases of the 16 byte unrolled code
        self.check_count_nonzero(12, 17)

    @dec.slow
    def test_count_nonzero_all(self):
        # check all combinations in a length 17 array
        # covers all cases of the 16 byte unrolled code
        self.check_count_nonzero(17, 17)

    def test_count_nonzero_unaligned(self):
        # prevent mistakes as e.g. gh-4060
        for o in range(7):
            a = np.zeros((18,), dtype=np.bool)[o+1:]
            a[:o] = True
            self.assertEqual(np.count_nonzero(a), builtins.sum(a.tolist()))
            a = np.ones((18,), dtype=np.bool)[o+1:]
            a[:o] = False
            self.assertEqual(np.count_nonzero(a), builtins.sum(a.tolist()))

class TestMethods(TestCase):
    def test_test_round(self):
        assert_equal(array([1.2, 1.5]).round(), [1, 2])
        assert_equal(array(1.5).round(), 2)
        assert_equal(array([12.2, 15.5]).round(-1), [10, 20])
        assert_equal(array([12.15, 15.51]).round(1), [12.2, 15.5])

    def test_transpose(self):
        a = array([[1, 2], [3, 4]])
        assert_equal(a.transpose(), [[1, 3], [2, 4]])
        self.assertRaises(ValueError, lambda: a.transpose(0))
        self.assertRaises(ValueError, lambda: a.transpose(0, 0))
        self.assertRaises(ValueError, lambda: a.transpose(0, 1, 2))

    def test_sort(self):
        # test ordering for floats and complex containing nans. It is only
        # necessary to check the lessthan comparison, so sorts that
        # only follow the insertion sort path are sufficient. We only
        # test doubles and complex doubles as the logic is the same.

        # check doubles
        msg = "Test real sort order with nans"
        a = np.array([np.nan, 1, 0])
        b = sort(a)
        assert_equal(b, a[::-1], msg)
        # check complex
        msg = "Test complex sort order with nans"
        a = np.zeros(9, dtype=np.complex128)
        a.real += [np.nan, np.nan, np.nan, 1, 0, 1, 1, 0, 0]
        a.imag += [np.nan, 1, 0, np.nan, np.nan, 1, 0, 1, 0]
        b = sort(a)
        assert_equal(b, a[::-1], msg)

        # all c scalar sorts use the same code with different types
        # so it suffices to run a quick check with one type. The number
        # of sorted items must be greater than ~50 to check the actual
        # algorithm because quick and merge sort fall over to insertion
        # sort for small arrays.
        a = np.arange(101)
        b = a[::-1].copy()
        for kind in ['q', 'm', 'h'] :
            msg = "scalar sort, kind=%s" % kind
            c = a.copy();
            c.sort(kind=kind)
            assert_equal(c, a, msg)
            c = b.copy();
            c.sort(kind=kind)
            assert_equal(c, a, msg)

        # test complex sorts. These use the same code as the scalars
        # but the compare fuction differs.
        ai = a*1j + 1
        bi = b*1j + 1
        for kind in ['q', 'm', 'h'] :
            msg = "complex sort, real part == 1, kind=%s" % kind
            c = ai.copy();
            c.sort(kind=kind)
            assert_equal(c, ai, msg)
            c = bi.copy();
            c.sort(kind=kind)
            assert_equal(c, ai, msg)
        ai = a + 1j
        bi = b + 1j
        for kind in ['q', 'm', 'h'] :
            msg = "complex sort, imag part == 1, kind=%s" % kind
            c = ai.copy();
            c.sort(kind=kind)
            assert_equal(c, ai, msg)
            c = bi.copy();
            c.sort(kind=kind)
            assert_equal(c, ai, msg)

        # test string sorts.
        s = 'aaaaaaaa'
        a = np.array([s + chr(i) for i in range(101)])
        b = a[::-1].copy()
        for kind in ['q', 'm', 'h'] :
            msg = "string sort, kind=%s" % kind
            c = a.copy();
            c.sort(kind=kind)
            assert_equal(c, a, msg)
            c = b.copy();
            c.sort(kind=kind)
            assert_equal(c, a, msg)

        # test unicode sorts.
        s = 'aaaaaaaa'
        a = np.array([s + chr(i) for i in range(101)], dtype=np.unicode)
        b = a[::-1].copy()
        for kind in ['q', 'm', 'h'] :
            msg = "unicode sort, kind=%s" % kind
            c = a.copy();
            c.sort(kind=kind)
            assert_equal(c, a, msg)
            c = b.copy();
            c.sort(kind=kind)
            assert_equal(c, a, msg)

        # test object array sorts.
        a = np.empty((101,), dtype=np.object)
        a[:] = list(range(101))
        b = a[::-1]
        for kind in ['q', 'h', 'm'] :
            msg = "object sort, kind=%s" % kind
            c = a.copy();
            c.sort(kind=kind)
            assert_equal(c, a, msg)
            c = b.copy();
            c.sort(kind=kind)
            assert_equal(c, a, msg)

        # test record array sorts.
        dt = np.dtype([('f', float), ('i', int)])
        a = array([(i, i) for i in range(101)], dtype = dt)
        b = a[::-1]
        for kind in ['q', 'h', 'm'] :
            msg = "object sort, kind=%s" % kind
            c = a.copy();
            c.sort(kind=kind)
            assert_equal(c, a, msg)
            c = b.copy();
            c.sort(kind=kind)
            assert_equal(c, a, msg)

        # test datetime64 sorts.
        a = np.arange(0, 101, dtype='datetime64[D]')
        b = a[::-1]
        for kind in ['q', 'h', 'm'] :
            msg = "datetime64 sort, kind=%s" % kind
            c = a.copy();
            c.sort(kind=kind)
            assert_equal(c, a, msg)
            c = b.copy();
            c.sort(kind=kind)
            assert_equal(c, a, msg)

        # test timedelta64 sorts.
        a = np.arange(0, 101, dtype='timedelta64[D]')
        b = a[::-1]
        for kind in ['q', 'h', 'm'] :
            msg = "timedelta64 sort, kind=%s" % kind
            c = a.copy();
            c.sort(kind=kind)
            assert_equal(c, a, msg)
            c = b.copy();
            c.sort(kind=kind)
            assert_equal(c, a, msg)

        # check axis handling. This should be the same for all type
        # specific sorts, so we only check it for one type and one kind
        a = np.array([[3, 2], [1, 0]])
        b = np.array([[1, 0], [3, 2]])
        c = np.array([[2, 3], [0, 1]])
        d = a.copy()
        d.sort(axis=0)
        assert_equal(d, b, "test sort with axis=0")
        d = a.copy()
        d.sort(axis=1)
        assert_equal(d, c, "test sort with axis=1")
        d = a.copy()
        d.sort()
        assert_equal(d, c, "test sort with default axis")

    def test_copy(self):
        def assert_fortran(arr):
            assert_(arr.flags.fortran)
            assert_(arr.flags.f_contiguous)
            assert_(not arr.flags.c_contiguous)

        def assert_c(arr):
            assert_(not arr.flags.fortran)
            assert_(not arr.flags.f_contiguous)
            assert_(arr.flags.c_contiguous)

        a = np.empty((2, 2), order='F')
        # Test copying a Fortran array
        assert_c(a.copy())
        assert_c(a.copy('C'))
        assert_fortran(a.copy('F'))
        assert_fortran(a.copy('A'))

        # Now test starting with a C array.
        a = np.empty((2, 2), order='C')
        assert_c(a.copy())
        assert_c(a.copy('C'))
        assert_fortran(a.copy('F'))
        assert_c(a.copy('A'))

    def test_sort_order(self):
        # Test sorting an array with fields
        x1=np.array([21, 32, 14])
        x2=np.array(['my', 'first', 'name'])
        x3=np.array([3.1, 4.5, 6.2])
        r=np.rec.fromarrays([x1, x2, x3], names='id,word,number')

        r.sort(order=['id'])
        assert_equal(r.id, array([14, 21, 32]))
        assert_equal(r.word, array(['name', 'my', 'first']))
        assert_equal(r.number, array([6.2, 3.1, 4.5]))

        r.sort(order=['word'])
        assert_equal(r.id, array([32, 21, 14]))
        assert_equal(r.word, array(['first', 'my', 'name']))
        assert_equal(r.number, array([4.5, 3.1, 6.2]))

        r.sort(order=['number'])
        assert_equal(r.id, array([21, 32, 14]))
        assert_equal(r.word, array(['my', 'first', 'name']))
        assert_equal(r.number, array([3.1, 4.5, 6.2]))

        if sys.byteorder == 'little':
            strtype = '>i2'
        else:
            strtype = '<i2'
        mydtype = [('name', strchar + '5'), ('col2', strtype)]
        r = np.array([('a', 1), ('b', 255), ('c', 3), ('d', 258)],
                     dtype= mydtype)
        r.sort(order='col2')
        assert_equal(r['col2'], [1, 3, 255, 258])
        assert_equal(r, np.array([('a', 1), ('c', 3), ('b', 255), ('d', 258)],
                                 dtype=mydtype))

    def test_argsort(self):
        # all c scalar argsorts use the same code with different types
        # so it suffices to run a quick check with one type. The number
        # of sorted items must be greater than ~50 to check the actual
        # algorithm because quick and merge sort fall over to insertion
        # sort for small arrays.
        a = np.arange(101)
        b = a[::-1].copy()
        for kind in ['q', 'm', 'h'] :
            msg = "scalar argsort, kind=%s" % kind
            assert_equal(a.copy().argsort(kind=kind), a, msg)
            assert_equal(b.copy().argsort(kind=kind), b, msg)

        # test complex argsorts. These use the same code as the scalars
        # but the compare fuction differs.
        ai = a*1j + 1
        bi = b*1j + 1
        for kind in ['q', 'm', 'h'] :
            msg = "complex argsort, kind=%s" % kind
            assert_equal(ai.copy().argsort(kind=kind), a, msg)
            assert_equal(bi.copy().argsort(kind=kind), b, msg)
        ai = a + 1j
        bi = b + 1j
        for kind in ['q', 'm', 'h'] :
            msg = "complex argsort, kind=%s" % kind
            assert_equal(ai.copy().argsort(kind=kind), a, msg)
            assert_equal(bi.copy().argsort(kind=kind), b, msg)

        # test string argsorts.
        s = 'aaaaaaaa'
        a = np.array([s + chr(i) for i in range(101)])
        b = a[::-1].copy()
        r = np.arange(101)
        rr = r[::-1]
        for kind in ['q', 'm', 'h'] :
            msg = "string argsort, kind=%s" % kind
            assert_equal(a.copy().argsort(kind=kind), r, msg)
            assert_equal(b.copy().argsort(kind=kind), rr, msg)

        # test unicode argsorts.
        s = 'aaaaaaaa'
        a = np.array([s + chr(i) for i in range(101)], dtype=np.unicode)
        b = a[::-1]
        r = np.arange(101)
        rr = r[::-1]
        for kind in ['q', 'm', 'h'] :
            msg = "unicode argsort, kind=%s" % kind
            assert_equal(a.copy().argsort(kind=kind), r, msg)
            assert_equal(b.copy().argsort(kind=kind), rr, msg)

        # test object array argsorts.
        a = np.empty((101,), dtype=np.object)
        a[:] = list(range(101))
        b = a[::-1]
        r = np.arange(101)
        rr = r[::-1]
        for kind in ['q', 'm', 'h'] :
            msg = "object argsort, kind=%s" % kind
            assert_equal(a.copy().argsort(kind=kind), r, msg)
            assert_equal(b.copy().argsort(kind=kind), rr, msg)

        # test structured array argsorts.
        dt = np.dtype([('f', float), ('i', int)])
        a = array([(i, i) for i in range(101)], dtype = dt)
        b = a[::-1]
        r = np.arange(101)
        rr = r[::-1]
        for kind in ['q', 'm', 'h'] :
            msg = "structured array argsort, kind=%s" % kind
            assert_equal(a.copy().argsort(kind=kind), r, msg)
            assert_equal(b.copy().argsort(kind=kind), rr, msg)

        # test datetime64 argsorts.
        a = np.arange(0, 101, dtype='datetime64[D]')
        b = a[::-1]
        r = np.arange(101)
        rr = r[::-1]
        for kind in ['q', 'h', 'm'] :
            msg = "datetime64 argsort, kind=%s" % kind
            assert_equal(a.copy().argsort(kind=kind), r, msg)
            assert_equal(b.copy().argsort(kind=kind), rr, msg)

        # test timedelta64 argsorts.
        a = np.arange(0, 101, dtype='timedelta64[D]')
        b = a[::-1]
        r = np.arange(101)
        rr = r[::-1]
        for kind in ['q', 'h', 'm'] :
            msg = "timedelta64 argsort, kind=%s" % kind
            assert_equal(a.copy().argsort(kind=kind), r, msg)
            assert_equal(b.copy().argsort(kind=kind), rr, msg)


        # check axis handling. This should be the same for all type
        # specific argsorts, so we only check it for one type and one kind
        a = np.array([[3, 2], [1, 0]])
        b = np.array([[1, 1], [0, 0]])
        c = np.array([[1, 0], [1, 0]])
        assert_equal(a.copy().argsort(axis=0), b)
        assert_equal(a.copy().argsort(axis=1), c)
        assert_equal(a.copy().argsort(), c)
        # using None is known fail at this point
        #assert_equal(a.copy().argsort(axis=None, c)

        # check that stable argsorts are stable
        r = np.arange(100)
        # scalars
        a = np.zeros(100)
        assert_equal(a.argsort(kind='m'), r)
        # complex
        a = np.zeros(100, dtype=np.complex)
        assert_equal(a.argsort(kind='m'), r)
        # string
        a = np.array(['aaaaaaaaa' for i in range(100)])
        assert_equal(a.argsort(kind='m'), r)
        # unicode
        a = np.array(['aaaaaaaaa' for i in range(100)], dtype=np.unicode)
        assert_equal(a.argsort(kind='m'), r)

    def test_sort_unicode_kind(self):
        d = np.arange(10)
        k = b'\xc3\xa4'.decode("UTF8")
        assert_raises(ValueError, d.sort, kind=k)
        assert_raises(ValueError, d.argsort, kind=k)

    def test_searchsorted(self):
        # test for floats and complex containing nans. The logic is the
        # same for all float types so only test double types for now.
        # The search sorted routines use the compare functions for the
        # array type, so this checks if that is consistent with the sort
        # order.

        # check double
        a = np.array([0, 1, np.nan])
        msg = "Test real searchsorted with nans, side='l'"
        b = a.searchsorted(a, side='l')
        assert_equal(b, np.arange(3), msg)
        msg = "Test real searchsorted with nans, side='r'"
        b = a.searchsorted(a, side='r')
        assert_equal(b, np.arange(1, 4), msg)
        # check double complex
        a = np.zeros(9, dtype=np.complex128)
        a.real += [0, 0, 1, 1, 0, 1, np.nan, np.nan, np.nan]
        a.imag += [0, 1, 0, 1, np.nan, np.nan, 0, 1, np.nan]
        msg = "Test complex searchsorted with nans, side='l'"
        b = a.searchsorted(a, side='l')
        assert_equal(b, np.arange(9), msg)
        msg = "Test complex searchsorted with nans, side='r'"
        b = a.searchsorted(a, side='r')
        assert_equal(b, np.arange(1, 10), msg)
        msg = "Test searchsorted with little endian, side='l'"
        a = np.array([0, 128], dtype='<i4')
        b = a.searchsorted(np.array(128, dtype='<i4'))
        assert_equal(b, 1, msg)
        msg = "Test searchsorted with big endian, side='l'"
        a = np.array([0, 128], dtype='>i4')
        b = a.searchsorted(np.array(128, dtype='>i4'))
        assert_equal(b, 1, msg)

        # Check 0 elements
        a = np.ones(0)
        b = a.searchsorted([0, 1, 2], 'l')
        assert_equal(b, [0, 0, 0])
        b = a.searchsorted([0, 1, 2], 'r')
        assert_equal(b, [0, 0, 0])
        a = np.ones(1)
        # Check 1 element
        b = a.searchsorted([0, 1, 2], 'l')
        assert_equal(b, [0, 0, 1])
        b = a.searchsorted([0, 1, 2], 'r')
        assert_equal(b, [0, 1, 1])
        # Check all elements equal
        a = np.ones(2)
        b = a.searchsorted([0, 1, 2], 'l')
        assert_equal(b, [0, 0, 2])
        b = a.searchsorted([0, 1, 2], 'r')
        assert_equal(b, [0, 2, 2])

        # Test searching unaligned array
        a = np.arange(10)
        aligned = np.empty(a.itemsize * a.size + 1, 'uint8')
        unaligned = aligned[1:].view(a.dtype)
        unaligned[:] = a
        # Test searching unaligned array
        b = unaligned.searchsorted(a, 'l')
        assert_equal(b, a)
        b = unaligned.searchsorted(a, 'r')
        assert_equal(b, a + 1)
        # Test searching for unaligned keys
        b = a.searchsorted(unaligned, 'l')
        assert_equal(b, a)
        b = a.searchsorted(unaligned, 'r')
        assert_equal(b, a + 1)

        # Test smart resetting of binsearch indices
        a = np.arange(5)
        b = a.searchsorted([6, 5, 4], 'l')
        assert_equal(b, [5, 5, 4])
        b = a.searchsorted([6, 5, 4], 'r')
        assert_equal(b, [5, 5, 5])

        # Test all type specific binary search functions
        types = ''.join((np.typecodes['AllInteger'], np.typecodes['AllFloat'],
                         np.typecodes['Datetime'], '?O'))
        for dt in types:
            if dt == 'M':
                dt = 'M8[D]'
            if dt == '?':
                a = np.arange(2, dtype=dt)
                out = np.arange(2)
            else:
                a = np.arange(0, 5, dtype=dt)
                out = np.arange(5)
            b = a.searchsorted(a, 'l')
            assert_equal(b, out)
            b = a.searchsorted(a, 'r')
            assert_equal(b, out + 1)


    def test_searchsorted_unicode(self):
        # Test searchsorted on unicode strings.

        # 1.6.1 contained a string length miscalculation in
        # arraytypes.c.src:UNICODE_compare() which manifested as
        # incorrect/inconsistent results from searchsorted.
        a = np.array(['P:\\20x_dapi_cy3\\20x_dapi_cy3_20100185_1',
                      'P:\\20x_dapi_cy3\\20x_dapi_cy3_20100186_1',
                      'P:\\20x_dapi_cy3\\20x_dapi_cy3_20100187_1',
                      'P:\\20x_dapi_cy3\\20x_dapi_cy3_20100189_1',
                      'P:\\20x_dapi_cy3\\20x_dapi_cy3_20100190_1',
                      'P:\\20x_dapi_cy3\\20x_dapi_cy3_20100191_1',
                      'P:\\20x_dapi_cy3\\20x_dapi_cy3_20100192_1',
                      'P:\\20x_dapi_cy3\\20x_dapi_cy3_20100193_1',
                      'P:\\20x_dapi_cy3\\20x_dapi_cy3_20100194_1',
                      'P:\\20x_dapi_cy3\\20x_dapi_cy3_20100195_1',
                      'P:\\20x_dapi_cy3\\20x_dapi_cy3_20100196_1',
                      'P:\\20x_dapi_cy3\\20x_dapi_cy3_20100197_1',
                      'P:\\20x_dapi_cy3\\20x_dapi_cy3_20100198_1',
                      'P:\\20x_dapi_cy3\\20x_dapi_cy3_20100199_1'],
                      dtype=np.unicode)
        ind = np.arange(len(a))
        assert_equal([a.searchsorted(v, 'left') for v in a], ind)
        assert_equal([a.searchsorted(v, 'right') for v in a], ind + 1)
        assert_equal([a.searchsorted(a[i], 'left') for i in ind], ind)
        assert_equal([a.searchsorted(a[i], 'right') for i in ind], ind + 1)

    def test_searchsorted_with_sorter(self):
        a = np.array([5, 2, 1, 3, 4])
        s = np.argsort(a)
        assert_raises(TypeError, np.searchsorted, a, 0, sorter=(1, (2, 3)))
        assert_raises(TypeError, np.searchsorted, a, 0, sorter=[1.1])
        assert_raises(ValueError, np.searchsorted, a, 0, sorter=[1, 2, 3, 4])
        assert_raises(ValueError, np.searchsorted, a, 0, sorter=[1, 2, 3, 4, 5, 6])

        # bounds check
        assert_raises(ValueError, np.searchsorted, a, 4, sorter=[0, 1, 2, 3, 5])
        assert_raises(ValueError, np.searchsorted, a, 0, sorter=[-1, 0, 1, 2, 3])
        assert_raises(ValueError, np.searchsorted, a, 0, sorter=[4, 0, -1, 2, 3])

        a = np.random.rand(300)
        s = a.argsort()
        b = np.sort(a)
        k = np.linspace(0, 1, 20)
        assert_equal(b.searchsorted(k), a.searchsorted(k, sorter=s))

        a = np.array([0, 1, 2, 3, 5]*20)
        s = a.argsort()
        k = [0, 1, 2, 3, 5]
        expected = [0, 20, 40, 60, 80]
        assert_equal(a.searchsorted(k, side='l', sorter=s), expected)
        expected = [20, 40, 60, 80, 100]
        assert_equal(a.searchsorted(k, side='r', sorter=s), expected)

        # Test searching unaligned array
        keys = np.arange(10)
        a = keys.copy()
        np.random.shuffle(s)
        s = a.argsort()
        aligned = np.empty(a.itemsize * a.size + 1, 'uint8')
        unaligned = aligned[1:].view(a.dtype)
        # Test searching unaligned array
        unaligned[:] = a
        b = unaligned.searchsorted(keys, 'l', s)
        assert_equal(b, keys)
        b = unaligned.searchsorted(keys, 'r', s)
        assert_equal(b, keys + 1)
        # Test searching for unaligned keys
        unaligned[:] = keys
        b = a.searchsorted(unaligned, 'l', s)
        assert_equal(b, keys)
        b = a.searchsorted(unaligned, 'r', s)
        assert_equal(b, keys + 1)

        # Test all type specific indirect binary search functions
        types = ''.join((np.typecodes['AllInteger'], np.typecodes['AllFloat'],
                         np.typecodes['Datetime'], '?O'))
        for dt in types:
            if dt == 'M':
                dt = 'M8[D]'
            if dt == '?':
                a = np.array([1, 0], dtype=dt)
                # We want the sorter array to be of a type that is different
                # from np.intp in all platforms, to check for #4698
                s = np.array([1, 0], dtype=np.int16)
                out = np.array([1, 0])
            else:
                a = np.array([3, 4, 1, 2, 0], dtype=dt)
                # We want the sorter array to be of a type that is different
                # from np.intp in all platforms, to check for #4698
                s = np.array([4, 2, 3, 0, 1], dtype=np.int16)
                out = np.array([3, 4, 1, 2, 0], dtype=np.intp)
            b = a.searchsorted(a, 'l', s)
            assert_equal(b, out)
            b = a.searchsorted(a, 'r', s)
            assert_equal(b, out + 1)

        # Test non-contiguous sorter array
        a = np.array([3, 4, 1, 2, 0])
        srt = np.empty((10,), dtype=np.intp)
        srt[1::2] = -1
        srt[::2] = [4, 2, 3, 0, 1]
        s = srt[::2]
        out = np.array([3, 4, 1, 2, 0], dtype=np.intp)
        b = a.searchsorted(a, 'l', s)
        assert_equal(b, out)
        b = a.searchsorted(a, 'r', s)
        assert_equal(b, out + 1)

    def test_partition(self):
        d = np.arange(10)
        assert_raises(TypeError, np.partition, d, 2, kind=1)
        assert_raises(ValueError, np.partition, d, 2, kind="nonsense")
        assert_raises(ValueError, np.argpartition, d, 2, kind="nonsense")
        assert_raises(ValueError, d.partition, 2, axis=0, kind="nonsense")
        assert_raises(ValueError, d.argpartition, 2, axis=0, kind="nonsense")
        for k in ("introselect",):
            d = np.array([])
            assert_array_equal(np.partition(d, 0, kind=k), d)
            assert_array_equal(np.argpartition(d, 0, kind=k), d)
            d = np.ones((1))
            assert_array_equal(np.partition(d, 0, kind=k)[0], d)
            assert_array_equal(d[np.argpartition(d, 0, kind=k)],
                               np.partition(d, 0, kind=k))

            # kth not modified
            kth = np.array([30, 15, 5])
            okth = kth.copy()
            np.partition(np.arange(40), kth)
            assert_array_equal(kth, okth)

            for r in ([2, 1], [1, 2], [1, 1]):
                d = np.array(r)
                tgt = np.sort(d)
                assert_array_equal(np.partition(d, 0, kind=k)[0], tgt[0])
                assert_array_equal(np.partition(d, 1, kind=k)[1], tgt[1])
                assert_array_equal(d[np.argpartition(d, 0, kind=k)],
                                   np.partition(d, 0, kind=k))
                assert_array_equal(d[np.argpartition(d, 1, kind=k)],
                                   np.partition(d, 1, kind=k))
                for i in range(d.size):
                    d[i:].partition(0, kind=k)
                assert_array_equal(d, tgt)

            for r in ([3, 2, 1], [1, 2, 3], [2, 1, 3], [2, 3, 1],
                      [1, 1, 1], [1, 2, 2], [2, 2, 1], [1, 2, 1]):
                d = np.array(r)
                tgt = np.sort(d)
                assert_array_equal(np.partition(d, 0, kind=k)[0], tgt[0])
                assert_array_equal(np.partition(d, 1, kind=k)[1], tgt[1])
                assert_array_equal(np.partition(d, 2, kind=k)[2], tgt[2])
                assert_array_equal(d[np.argpartition(d, 0, kind=k)],
                                   np.partition(d, 0, kind=k))
                assert_array_equal(d[np.argpartition(d, 1, kind=k)],
                                   np.partition(d, 1, kind=k))
                assert_array_equal(d[np.argpartition(d, 2, kind=k)],
                                   np.partition(d, 2, kind=k))
                for i in range(d.size):
                    d[i:].partition(0, kind=k)
                assert_array_equal(d, tgt)

            d = np.ones((50))
            assert_array_equal(np.partition(d, 0, kind=k), d)
            assert_array_equal(d[np.argpartition(d, 0, kind=k)],
                               np.partition(d, 0, kind=k))

            # sorted
            d = np.arange((49))
            self.assertEqual(np.partition(d, 5, kind=k)[5], 5)
            self.assertEqual(np.partition(d, 15, kind=k)[15], 15)
            assert_array_equal(d[np.argpartition(d, 5, kind=k)],
                               np.partition(d, 5, kind=k))
            assert_array_equal(d[np.argpartition(d, 15, kind=k)],
                               np.partition(d, 15, kind=k))

            # rsorted
            d = np.arange((47))[::-1]
            self.assertEqual(np.partition(d, 6, kind=k)[6], 6)
            self.assertEqual(np.partition(d, 16, kind=k)[16], 16)
            assert_array_equal(d[np.argpartition(d, 6, kind=k)],
                               np.partition(d, 6, kind=k))
            assert_array_equal(d[np.argpartition(d, 16, kind=k)],
                               np.partition(d, 16, kind=k))

            assert_array_equal(np.partition(d, -6, kind=k),
                               np.partition(d, 41, kind=k))
            assert_array_equal(np.partition(d, -16, kind=k),
                               np.partition(d, 31, kind=k))
            assert_array_equal(d[np.argpartition(d, -6, kind=k)],
                               np.partition(d, 41, kind=k))

            # median of 3 killer, O(n^2) on pure median 3 pivot quickselect
            # exercises the median of median of 5 code used to keep O(n)
            d = np.arange(1000000)
            x = np.roll(d, d.size // 2)
            mid = x.size // 2 + 1
            assert_equal(np.partition(x, mid)[mid], mid)
            d = np.arange(1000001)
            x = np.roll(d, d.size // 2 + 1)
            mid = x.size // 2 + 1
            assert_equal(np.partition(x, mid)[mid], mid)

            # max
            d = np.ones(10); d[1] = 4;
            assert_equal(np.partition(d, (2, -1))[-1], 4)
            assert_equal(np.partition(d, (2, -1))[2], 1)
            assert_equal(d[np.argpartition(d, (2, -1))][-1], 4)
            assert_equal(d[np.argpartition(d, (2, -1))][2], 1)
            d[1] = np.nan
            assert_(np.isnan(d[np.argpartition(d, (2, -1))][-1]))
            assert_(np.isnan(np.partition(d, (2, -1))[-1]))

            # equal elements
            d = np.arange((47)) % 7
            tgt = np.sort(np.arange((47)) % 7)
            np.random.shuffle(d)
            for i in range(d.size):
                self.assertEqual(np.partition(d, i, kind=k)[i], tgt[i])
            assert_array_equal(d[np.argpartition(d, 6, kind=k)],
                               np.partition(d, 6, kind=k))
            assert_array_equal(d[np.argpartition(d, 16, kind=k)],
                               np.partition(d, 16, kind=k))
            for i in range(d.size):
                d[i:].partition(0, kind=k)
            assert_array_equal(d, tgt)

            d = np.array([0, 1, 2, 3, 4, 5, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
                          7, 7, 7, 7, 7, 9])
            kth = [0, 3, 19, 20]
            assert_equal(np.partition(d, kth, kind=k)[kth], (0, 3, 7, 7))
            assert_equal(d[np.argpartition(d, kth, kind=k)][kth], (0, 3, 7, 7))

            d = np.array([2, 1])
            d.partition(0, kind=k)
            assert_raises(ValueError, d.partition, 2)
            assert_raises(ValueError, d.partition, 3, axis=1)
            assert_raises(ValueError, np.partition, d, 2)
            assert_raises(ValueError, np.partition, d, 2, axis=1)
            assert_raises(ValueError, d.argpartition, 2)
            assert_raises(ValueError, d.argpartition, 3, axis=1)
            assert_raises(ValueError, np.argpartition, d, 2)
            assert_raises(ValueError, np.argpartition, d, 2, axis=1)
            d = np.arange(10).reshape((2, 5))
            d.partition(1, axis=0, kind=k)
            d.partition(4, axis=1, kind=k)
            np.partition(d, 1, axis=0, kind=k)
            np.partition(d, 4, axis=1, kind=k)
            np.partition(d, 1, axis=None, kind=k)
            np.partition(d, 9, axis=None, kind=k)
            d.argpartition(1, axis=0, kind=k)
            d.argpartition(4, axis=1, kind=k)
            np.argpartition(d, 1, axis=0, kind=k)
            np.argpartition(d, 4, axis=1, kind=k)
            np.argpartition(d, 1, axis=None, kind=k)
            np.argpartition(d, 9, axis=None, kind=k)
            assert_raises(ValueError, d.partition, 2, axis=0)
            assert_raises(ValueError, d.partition, 11, axis=1)
            assert_raises(TypeError, d.partition, 2, axis=None)
            assert_raises(ValueError, np.partition, d, 9, axis=1)
            assert_raises(ValueError, np.partition, d, 11, axis=None)
            assert_raises(ValueError, d.argpartition, 2, axis=0)
            assert_raises(ValueError, d.argpartition, 11, axis=1)
            assert_raises(ValueError, np.argpartition, d, 9, axis=1)
            assert_raises(ValueError, np.argpartition, d, 11, axis=None)

            td = [(dt, s) for dt in [np.int32, np.float32, np.complex64]
                          for s in (9, 16)]
            for dt, s in td:
                aae = assert_array_equal
                at = self.assertTrue

                d = np.arange(s, dtype=dt)
                np.random.shuffle(d)
                d1 = np.tile(np.arange(s, dtype=dt), (4, 1))
                map(np.random.shuffle, d1)
                d0 = np.transpose(d1)
                for i in range(d.size):
                    p = np.partition(d, i, kind=k)
                    self.assertEqual(p[i], i)
                    # all before are smaller
                    assert_array_less(p[:i], p[i])
                    # all after are larger
                    assert_array_less(p[i], p[i + 1:])
                    aae(p, d[np.argpartition(d, i, kind=k)])

                    p = np.partition(d1, i, axis=1, kind=k)
                    aae(p[:, i], np.array([i] * d1.shape[0], dtype=dt))
                    # array_less does not seem to work right
                    at((p[:, :i].T <= p[:, i]).all(),
                       msg="%d: %r <= %r" % (i, p[:, i], p[:, :i].T))
                    at((p[:, i + 1:].T > p[:, i]).all(),
                       msg="%d: %r < %r" % (i, p[:, i], p[:, i + 1:].T))
                    aae(p, d1[np.arange(d1.shape[0])[:, None],
                        np.argpartition(d1, i, axis=1, kind=k)])

                    p = np.partition(d0, i, axis=0, kind=k)
                    aae(p[i,:], np.array([i] * d1.shape[0],
                                                         dtype=dt))
                    # array_less does not seem to work right
                    at((p[:i,:] <= p[i,:]).all(),
                       msg="%d: %r <= %r" % (i, p[i,:], p[:i,:]))
                    at((p[i + 1:,:] > p[i,:]).all(),
                       msg="%d: %r < %r" % (i, p[i,:], p[:, i + 1:]))
                    aae(p, d0[np.argpartition(d0, i, axis=0, kind=k),
                        np.arange(d0.shape[1])[None,:]])

                    # check inplace
                    dc = d.copy()
                    dc.partition(i, kind=k)
                    assert_equal(dc, np.partition(d, i, kind=k))
                    dc = d0.copy()
                    dc.partition(i, axis=0, kind=k)
                    assert_equal(dc, np.partition(d0, i, axis=0, kind=k))
                    dc = d1.copy()
                    dc.partition(i, axis=1, kind=k)
                    assert_equal(dc, np.partition(d1, i, axis=1, kind=k))


    def assert_partitioned(self, d, kth):
        prev = 0
        for k in np.sort(kth):
            assert_array_less(d[prev:k], d[k], err_msg='kth %d' % k)
            assert_((d[k:] >= d[k]).all(),
                    msg="kth %d, %r not greater equal %d" % (k, d[k:], d[k]))
            prev = k + 1


    def test_partition_iterative(self):
            d = np.arange(17)
            kth = (0, 1, 2, 429, 231)
            assert_raises(ValueError, d.partition, kth)
            assert_raises(ValueError, d.argpartition, kth)
            d = np.arange(10).reshape((2, 5))
            assert_raises(ValueError, d.partition, kth, axis=0)
            assert_raises(ValueError, d.partition, kth, axis=1)
            assert_raises(ValueError, np.partition, d, kth, axis=1)
            assert_raises(ValueError, np.partition, d, kth, axis=None)

            d = np.array([3, 4, 2, 1])
            p = np.partition(d, (0, 3))
            self.assert_partitioned(p, (0, 3))
            self.assert_partitioned(d[np.argpartition(d, (0, 3))], (0, 3))

            assert_array_equal(p, np.partition(d, (-3, -1)))
            assert_array_equal(p, d[np.argpartition(d, (-3, -1))])

            d = np.arange(17)
            np.random.shuffle(d)
            d.partition(range(d.size))
            assert_array_equal(np.arange(17), d)
            np.random.shuffle(d)
            assert_array_equal(np.arange(17), d[d.argpartition(range(d.size))])

            # test unsorted kth
            d = np.arange(17)
            np.random.shuffle(d)
            keys = np.array([1, 3, 8, -2])
            np.random.shuffle(d)
            p = np.partition(d, keys)
            self.assert_partitioned(p, keys)
            p = d[np.argpartition(d, keys)]
            self.assert_partitioned(p, keys)
            np.random.shuffle(keys)
            assert_array_equal(np.partition(d, keys), p)
            assert_array_equal(d[np.argpartition(d, keys)], p)

            # equal kth
            d = np.arange(20)[::-1]
            self.assert_partitioned(np.partition(d, [5]*4), [5])
            self.assert_partitioned(np.partition(d, [5]*4 + [6, 13]),
                                    [5]*4 + [6, 13])
            self.assert_partitioned(d[np.argpartition(d, [5]*4)], [5])
            self.assert_partitioned(d[np.argpartition(d, [5]*4 + [6, 13])],
                                    [5]*4 + [6, 13])

            d = np.arange(12)
            np.random.shuffle(d)
            d1 = np.tile(np.arange(12), (4, 1))
            map(np.random.shuffle, d1)
            d0 = np.transpose(d1)

            kth = (1, 6, 7, -1)
            p = np.partition(d1, kth, axis=1)
            pa = d1[np.arange(d1.shape[0])[:, None],
                    d1.argpartition(kth, axis=1)]
            assert_array_equal(p, pa)
            for i in range(d1.shape[0]):
                self.assert_partitioned(p[i,:], kth)
            p = np.partition(d0, kth, axis=0)
            pa = d0[np.argpartition(d0, kth, axis=0),
                    np.arange(d0.shape[1])[None,:]]
            assert_array_equal(p, pa)
            for i in range(d0.shape[1]):
                self.assert_partitioned(p[:, i], kth)


    def test_partition_cdtype(self):
        d = array([('Galahad', 1.7, 38), ('Arthur', 1.8, 41),
                   ('Lancelot', 1.9, 38)],
                  dtype=[('name', '|S10'), ('height', '<f8'), ('age', '<i4')])

        tgt = np.sort(d, order=['age', 'height'])
        assert_array_equal(np.partition(d, range(d.size),
                                        order=['age', 'height']),
                           tgt)
        assert_array_equal(d[np.argpartition(d, range(d.size),
                                             order=['age', 'height'])],
                           tgt)
        for k in range(d.size):
            assert_equal(np.partition(d, k, order=['age', 'height'])[k],
                        tgt[k])
            assert_equal(d[np.argpartition(d, k, order=['age', 'height'])][k],
                         tgt[k])

        d = array(['Galahad', 'Arthur', 'zebra', 'Lancelot'])
        tgt = np.sort(d)
        assert_array_equal(np.partition(d, range(d.size)), tgt)
        for k in range(d.size):
            assert_equal(np.partition(d, k)[k], tgt[k])
            assert_equal(d[np.argpartition(d, k)][k], tgt[k])

    def test_partition_unicode_kind(self):
        d = np.arange(10)
        k = b'\xc3\xa4'.decode("UTF8")
        assert_raises(ValueError, d.partition, 2, kind=k)
        assert_raises(ValueError, d.argpartition, 2, kind=k)

    def test_partition_fuzz(self):
        # a few rounds of random data testing
       for j in range(10, 30):
           for i in range(1, j - 2):
               d = np.arange(j)
               np.random.shuffle(d)
               d = d % np.random.randint(2, 30)
               idx = np.random.randint(d.size)
               kth = [0, idx, i, i + 1]
               tgt = np.sort(d)[kth]
               assert_array_equal(np.partition(d, kth)[kth], tgt,
                                  err_msg="data: %r\n kth: %r" % (d, kth))

    def test_flatten(self):
        x0 = np.array([[1, 2, 3], [4, 5, 6]], np.int32)
        x1 = np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]], np.int32)
        y0 = np.array([1, 2, 3, 4, 5, 6], np.int32)
        y0f = np.array([1, 4, 2, 5, 3, 6], np.int32)
        y1 = np.array([1, 2, 3, 4, 5, 6, 7, 8], np.int32)
        y1f = np.array([1, 5, 3, 7, 2, 6, 4, 8], np.int32)
        assert_equal(x0.flatten(), y0)
        assert_equal(x0.flatten('F'), y0f)
        assert_equal(x0.flatten('F'), x0.T.flatten())
        assert_equal(x1.flatten(), y1)
        assert_equal(x1.flatten('F'), y1f)
        assert_equal(x1.flatten('F'), x1.T.flatten())

    def test_dot(self):
        a = np.array([[1, 0], [0, 1]])
        b = np.array([[0, 1], [1, 0]])
        c = np.array([[9, 1], [1, -9]])

        assert_equal(np.dot(a, b), a.dot(b))
        assert_equal(np.dot(np.dot(a, b), c), a.dot(b).dot(c))

        # test passing in an output array
        c = np.zeros_like(a)
        a.dot(b, c)
        assert_equal(c, np.dot(a, b))

        # test keyword args
        c = np.zeros_like(a)
        a.dot(b=b, out=c)
        assert_equal(c, np.dot(a, b))

    @dec.skipif(True) # ufunc override disabled for 1.9
    def test_dot_override(self):
        class A(object):
            def __numpy_ufunc__(self, ufunc, method, pos, inputs, **kwargs):
                return "A"

        class B(object):
            def __numpy_ufunc__(self, ufunc, method, pos, inputs, **kwargs):
                return NotImplemented

        a = A()
        b = B()
        c = np.array([[1]])

        assert_equal(np.dot(a, b), "A")
        assert_equal(c.dot(a), "A")
        assert_raises(TypeError, np.dot, b, c)
        assert_raises(TypeError, c.dot, b)

    def test_diagonal(self):
        a = np.arange(12).reshape((3, 4))
        assert_equal(a.diagonal(), [0, 5, 10])
        assert_equal(a.diagonal(0), [0, 5, 10])
        assert_equal(a.diagonal(1), [1, 6, 11])
        assert_equal(a.diagonal(-1), [4, 9])

        b = np.arange(8).reshape((2, 2, 2))
        assert_equal(b.diagonal(), [[0, 6], [1, 7]])
        assert_equal(b.diagonal(0), [[0, 6], [1, 7]])
        assert_equal(b.diagonal(1), [[2], [3]])
        assert_equal(b.diagonal(-1), [[4], [5]])
        assert_raises(ValueError, b.diagonal, axis1=0, axis2=0)
        assert_equal(b.diagonal(0, 1, 2), [[0, 3], [4, 7]])
        assert_equal(b.diagonal(0, 0, 1), [[0, 6], [1, 7]])
        assert_equal(b.diagonal(offset=1, axis1=0, axis2=2), [[1], [3]])
        # Order of axis argument doesn't matter:
        assert_equal(b.diagonal(0, 2, 1), [[0, 3], [4, 7]])

    def test_diagonal_view_notwriteable(self):
        # this test is only for 1.9, the diagonal view will be
        # writeable in 1.10.
        a = np.eye(3).diagonal()
        assert_(not a.flags.writeable)
        assert_(not a.flags.owndata)

        a = np.diagonal(np.eye(3))
        assert_(not a.flags.writeable)
        assert_(not a.flags.owndata)

        a = np.diag(np.eye(3))
        assert_(not a.flags.writeable)
        assert_(not a.flags.owndata)

    def test_diagonal_memleak(self):
        # Regression test for a bug that crept in at one point
        a = np.zeros((100, 100))
        assert_(sys.getrefcount(a) < 50)
        for i in range(100):
            a.diagonal()
        assert_(sys.getrefcount(a) < 50)

    def test_ravel(self):
        a = np.array([[0, 1], [2, 3]])
        assert_equal(a.ravel(), [0, 1, 2, 3])
        assert_(not a.ravel().flags.owndata)
        assert_equal(a.ravel('F'), [0, 2, 1, 3])
        assert_equal(a.ravel(order='C'), [0, 1, 2, 3])
        assert_equal(a.ravel(order='F'), [0, 2, 1, 3])
        assert_equal(a.ravel(order='A'), [0, 1, 2, 3])
        assert_(not a.ravel(order='A').flags.owndata)
        assert_equal(a.ravel(order='K'), [0, 1, 2, 3])
        assert_(not a.ravel(order='K').flags.owndata)
        assert_equal(a.ravel(), a.reshape(-1))

        a = np.array([[0, 1], [2, 3]], order='F')
        assert_equal(a.ravel(), [0, 1, 2, 3])
        assert_equal(a.ravel(order='A'), [0, 2, 1, 3])
        assert_equal(a.ravel(order='K'), [0, 2, 1, 3])
        assert_(not a.ravel(order='A').flags.owndata)
        assert_(not a.ravel(order='K').flags.owndata)
        assert_equal(a.ravel(), a.reshape(-1))
        assert_equal(a.ravel(order='A'), a.reshape(-1, order='A'))

        a = np.array([[0, 1], [2, 3]])[::-1,:]
        assert_equal(a.ravel(), [2, 3, 0, 1])
        assert_equal(a.ravel(order='C'), [2, 3, 0, 1])
        assert_equal(a.ravel(order='F'), [2, 0, 3, 1])
        assert_equal(a.ravel(order='A'), [2, 3, 0, 1])
        # 'K' doesn't reverse the axes of negative strides
        assert_equal(a.ravel(order='K'), [2, 3, 0, 1])
        assert_(a.ravel(order='K').flags.owndata)

    def test_conjugate(self):
        a = np.array([1-1j, 1+1j, 23+23.0j])
        ac = a.conj()
        assert_equal(a.real, ac.real)
        assert_equal(a.imag, -ac.imag)
        assert_equal(ac, a.conjugate())
        assert_equal(ac, np.conjugate(a))

        a = np.array([1-1j, 1+1j, 23+23.0j], 'F')
        ac = a.conj()
        assert_equal(a.real, ac.real)
        assert_equal(a.imag, -ac.imag)
        assert_equal(ac, a.conjugate())
        assert_equal(ac, np.conjugate(a))

        a = np.array([1, 2, 3])
        ac = a.conj()
        assert_equal(a, ac)
        assert_equal(ac, a.conjugate())
        assert_equal(ac, np.conjugate(a))

        a = np.array([1.0, 2.0, 3.0])
        ac = a.conj()
        assert_equal(a, ac)
        assert_equal(ac, a.conjugate())
        assert_equal(ac, np.conjugate(a))

        a = np.array([1-1j, 1+1j, 1, 2.0], object)
        ac = a.conj()
        assert_equal(ac, [k.conjugate() for k in a])
        assert_equal(ac, a.conjugate())
        assert_equal(ac, np.conjugate(a))

        a = np.array([1-1j, 1, 2.0, 'f'], object)
        assert_raises(AttributeError, lambda: a.conj())
        assert_raises(AttributeError, lambda: a.conjugate()) 


class TestBinop(object):
    @dec.skipif(True) # ufunc override disabled for 1.9
    def test_ufunc_override_rop_precedence(self):
        # Check that __rmul__ and other right-hand operations have
        # precedence over __numpy_ufunc__

        ops = {
            '__add__':      ('__radd__', np.add, True),
            '__sub__':      ('__rsub__', np.subtract, True),
            '__mul__':      ('__rmul__', np.multiply, True),
            '__truediv__':  ('__rtruediv__', np.true_divide, True),
            '__floordiv__': ('__rfloordiv__', np.floor_divide, True),
            '__mod__':      ('__rmod__', np.remainder, True),
            '__divmod__':   ('__rdivmod__', None, False),
            '__pow__':      ('__rpow__', np.power, True),
            '__lshift__':   ('__rlshift__', np.left_shift, True),
            '__rshift__':   ('__rrshift__', np.right_shift, True),
            '__and__':      ('__rand__', np.bitwise_and, True),
            '__xor__':      ('__rxor__', np.bitwise_xor, True),
            '__or__':       ('__ror__', np.bitwise_or, True),
            '__ge__':       ('__le__', np.less_equal, False),
            '__gt__':       ('__lt__', np.less, False),
            '__le__':       ('__ge__', np.greater_equal, False),
            '__lt__':       ('__gt__', np.greater, False),
            '__eq__':       ('__eq__', np.equal, False),
            '__ne__':       ('__ne__', np.not_equal, False),
        }

        class OtherNdarraySubclass(ndarray):
            pass

        class OtherNdarraySubclassWithOverride(ndarray):
            def __numpy_ufunc__(self, *a, **kw):
                raise AssertionError(("__numpy_ufunc__ %r %r shouldn't have "
                                      "been called!") % (a, kw))

        def check(op_name, ndsubclass):
            rop_name, np_op, has_iop = ops[op_name]

            if has_iop:
                iop_name = '__i' + op_name[2:]
                iop = getattr(operator, iop_name)

            if op_name == "__divmod__":
                op = divmod
            else:
                op = getattr(operator, op_name)

            # Dummy class
            def __init__(self, *a, **kw):
                pass

            def __numpy_ufunc__(self, *a, **kw):
                raise AssertionError(("__numpy_ufunc__ %r %r shouldn't have "
                                      "been called!") % (a, kw))

            def __op__(self, *other):
                return "op"

            def __rop__(self, *other):
                return "rop"

            if ndsubclass:
                bases = (ndarray,)
            else:
                bases = (object,)

            dct = {'__init__': __init__,
                   '__numpy_ufunc__': __numpy_ufunc__,
                   op_name: __op__}
            if op_name != rop_name:
                dct[rop_name] = __rop__

            cls = type("Rop" + rop_name, bases, dct)

            # Check behavior against both bare ndarray objects and a
            # ndarray subclasses with and without their own override
            obj = cls((1,), buffer=np.ones(1,))

            arr_objs = [np.array([1]),
                        np.array([2]).view(OtherNdarraySubclass),
                        np.array([3]).view(OtherNdarraySubclassWithOverride),
                        ]

            for arr in arr_objs:
                err_msg = "%r %r" % (op_name, arr,)

                # Check that ndarray op gives up if it sees a non-subclass
                if not isinstance(obj, arr.__class__):
                    assert_equal(getattr(arr, op_name)(obj),
                                 NotImplemented, err_msg=err_msg)

                # Check that the Python binops have priority
                assert_equal(op(obj, arr), "op", err_msg=err_msg)
                if op_name == rop_name:
                    assert_equal(op(arr, obj), "op", err_msg=err_msg)
                else:
                    assert_equal(op(arr, obj), "rop", err_msg=err_msg)

                # Check that Python binops have priority also for in-place ops
                if has_iop:
                    assert_equal(getattr(arr, iop_name)(obj),
                                 NotImplemented, err_msg=err_msg)
                    if op_name != "__pow__":
                        # inplace pow requires the other object to be
                        # integer-like?
                        assert_equal(iop(arr, obj), "rop", err_msg=err_msg)

                # Check that ufunc call __numpy_ufunc__ normally
                if np_op is not None:
                    assert_raises(AssertionError, np_op, arr, obj,
                                  err_msg=err_msg)
                    assert_raises(AssertionError, np_op, obj, arr,
                                  err_msg=err_msg)

        # Check all binary operations
        for op_name in sorted(ops.keys()):
            yield check, op_name, True
            yield check, op_name, False

    @dec.skipif(True) # ufunc override disabled for 1.9
    def test_ufunc_override_rop_simple(self):
        # Check parts of the binary op overriding behavior in an
        # explicit test case that is easier to understand.

        class SomeClass(object):
            def __numpy_ufunc__(self, *a, **kw):
                return "ufunc"
            def __mul__(self, other):
                return 123
            def __rmul__(self, other):
                return 321
            def __gt__(self, other):
                return "yep"
            def __lt__(self, other):
                return "nope"

        class SomeClass2(SomeClass, ndarray):
            def __numpy_ufunc__(self, ufunc, method, i, inputs, **kw):
                if ufunc is np.multiply:
                    return "ufunc"
                else:
                    inputs = list(inputs)
                    inputs[i] = np.asarray(self)
                    func = getattr(ufunc, method)
                    r = func(*inputs, **kw)
                    if 'out' in kw:
                        return r
                    else:
                        x = SomeClass2(r.shape, dtype=r.dtype)
                        x[...] = r
                        return x

        arr = np.array([0])
        obj = SomeClass()
        obj2 = SomeClass2((1,), dtype=np.int_)
        obj2[0] = 9

        assert_equal(obj * arr, 123)
        assert_equal(arr * obj, 321)
        assert_equal(arr > obj, "nope")
        assert_equal(arr < obj, "yep")
        assert_equal(np.multiply(arr, obj), "ufunc")
        arr *= obj
        assert_equal(arr, 321)

        assert_equal(obj2 * arr, 123)
        assert_equal(arr * obj2, 321)
        assert_equal(arr > obj2, "nope")
        assert_equal(arr < obj2, "yep")
        assert_equal(np.multiply(arr, obj2), "ufunc")
        arr *= obj2
        assert_equal(arr, 321)

        obj2 += 33
        assert_equal(obj2[0], 42)
        assert_equal(obj2.sum(), 42)
        assert_(isinstance(obj2, SomeClass2))


class TestSubscripting(TestCase):
    def test_test_zero_rank(self):
        x = array([1, 2, 3])
        self.assertTrue(isinstance(x[0], np.int_))
        if sys.version_info[0] < 3:
            self.assertTrue(isinstance(x[0], int))
        self.assertTrue(type(x[0, ...]) is ndarray)


class TestPickling(TestCase):
    def test_roundtrip(self):
        import pickle
        carray = array([[2, 9], [7, 0], [3, 8]])
        DATA = [
            carray,
            transpose(carray),
            array([('xxx', 1, 2.0)], dtype=[('a', (str, 3)), ('b', int),
                                            ('c', float)])
        ]

        for a in DATA:
            assert_equal(a, pickle.loads(a.dumps()), err_msg="%r" % a)

    def _loads(self, obj):
        if sys.version_info[0] >= 3:
            return loads(obj, encoding='latin1')
        else:
            return loads(obj)

    # version 0 pickles, using protocol=2 to pickle
    # version 0 doesn't have a version field
    def test_version0_int8(self):
        s = '\x80\x02cnumpy.core._internal\n_reconstruct\nq\x01cnumpy\nndarray\nq\x02K\x00\x85U\x01b\x87Rq\x03(K\x04\x85cnumpy\ndtype\nq\x04U\x02i1K\x00K\x01\x87Rq\x05(U\x01|NNJ\xff\xff\xff\xffJ\xff\xff\xff\xfftb\x89U\x04\x01\x02\x03\x04tb.'
        a = array([1, 2, 3, 4], dtype=int8)
        p = self._loads(asbytes(s))
        assert_equal(a, p)

    def test_version0_float32(self):
        s = '\x80\x02cnumpy.core._internal\n_reconstruct\nq\x01cnumpy\nndarray\nq\x02K\x00\x85U\x01b\x87Rq\x03(K\x04\x85cnumpy\ndtype\nq\x04U\x02f4K\x00K\x01\x87Rq\x05(U\x01<NNJ\xff\xff\xff\xffJ\xff\xff\xff\xfftb\x89U\x10\x00\x00\x80?\x00\x00\x00@\x00\x00@@\x00\x00\x80@tb.'
        a = array([1.0, 2.0, 3.0, 4.0], dtype=float32)
        p = self._loads(asbytes(s))
        assert_equal(a, p)

    def test_version0_object(self):
        s = '\x80\x02cnumpy.core._internal\n_reconstruct\nq\x01cnumpy\nndarray\nq\x02K\x00\x85U\x01b\x87Rq\x03(K\x02\x85cnumpy\ndtype\nq\x04U\x02O8K\x00K\x01\x87Rq\x05(U\x01|NNJ\xff\xff\xff\xffJ\xff\xff\xff\xfftb\x89]q\x06(}q\x07U\x01aK\x01s}q\x08U\x01bK\x02setb.'
        a = np.array([{'a':1}, {'b':2}])
        p = self._loads(asbytes(s))
        assert_equal(a, p)

    # version 1 pickles, using protocol=2 to pickle
    def test_version1_int8(self):
        s = '\x80\x02cnumpy.core._internal\n_reconstruct\nq\x01cnumpy\nndarray\nq\x02K\x00\x85U\x01b\x87Rq\x03(K\x01K\x04\x85cnumpy\ndtype\nq\x04U\x02i1K\x00K\x01\x87Rq\x05(K\x01U\x01|NNJ\xff\xff\xff\xffJ\xff\xff\xff\xfftb\x89U\x04\x01\x02\x03\x04tb.'
        a = array([1, 2, 3, 4], dtype=int8)
        p = self._loads(asbytes(s))
        assert_equal(a, p)

    def test_version1_float32(self):
        s = '\x80\x02cnumpy.core._internal\n_reconstruct\nq\x01cnumpy\nndarray\nq\x02K\x00\x85U\x01b\x87Rq\x03(K\x01K\x04\x85cnumpy\ndtype\nq\x04U\x02f4K\x00K\x01\x87Rq\x05(K\x01U\x01<NNJ\xff\xff\xff\xffJ\xff\xff\xff\xfftb\x89U\x10\x00\x00\x80?\x00\x00\x00@\x00\x00@@\x00\x00\x80@tb.'
        a = array([1.0, 2.0, 3.0, 4.0], dtype=float32)
        p = self._loads(asbytes(s))
        assert_equal(a, p)

    def test_version1_object(self):
        s = '\x80\x02cnumpy.core._internal\n_reconstruct\nq\x01cnumpy\nndarray\nq\x02K\x00\x85U\x01b\x87Rq\x03(K\x01K\x02\x85cnumpy\ndtype\nq\x04U\x02O8K\x00K\x01\x87Rq\x05(K\x01U\x01|NNJ\xff\xff\xff\xffJ\xff\xff\xff\xfftb\x89]q\x06(}q\x07U\x01aK\x01s}q\x08U\x01bK\x02setb.'
        a = array([{'a':1}, {'b':2}])
        p = self._loads(asbytes(s))
        assert_equal(a, p)

    def test_subarray_int_shape(self):
        s = "cnumpy.core.multiarray\n_reconstruct\np0\n(cnumpy\nndarray\np1\n(I0\ntp2\nS'b'\np3\ntp4\nRp5\n(I1\n(I1\ntp6\ncnumpy\ndtype\np7\n(S'V6'\np8\nI0\nI1\ntp9\nRp10\n(I3\nS'|'\np11\nN(S'a'\np12\ng3\ntp13\n(dp14\ng12\n(g7\n(S'V4'\np15\nI0\nI1\ntp16\nRp17\n(I3\nS'|'\np18\n(g7\n(S'i1'\np19\nI0\nI1\ntp20\nRp21\n(I3\nS'|'\np22\nNNNI-1\nI-1\nI0\ntp23\nb(I2\nI2\ntp24\ntp25\nNNI4\nI1\nI0\ntp26\nbI0\ntp27\nsg3\n(g7\n(S'V2'\np28\nI0\nI1\ntp29\nRp30\n(I3\nS'|'\np31\n(g21\nI2\ntp32\nNNI2\nI1\nI0\ntp33\nbI4\ntp34\nsI6\nI1\nI0\ntp35\nbI00\nS'\\x01\\x01\\x01\\x01\\x01\\x02'\np36\ntp37\nb."
        a = np.array([(1, (1, 2))], dtype=[('a', 'i1', (2, 2)), ('b', 'i1', 2)])
        p = self._loads(asbytes(s))
        assert_equal(a, p)


class TestFancyIndexing(TestCase):
    def test_list(self):
        x = ones((1, 1))
        x[:, [0]] = 2.0
        assert_array_equal(x, array([[2.0]]))

        x = ones((1, 1, 1))
        x[:,:, [0]] = 2.0
        assert_array_equal(x, array([[[2.0]]]))

    def test_tuple(self):
        x = ones((1, 1))
        x[:, (0,)] = 2.0
        assert_array_equal(x, array([[2.0]]))
        x = ones((1, 1, 1))
        x[:,:, (0,)] = 2.0
        assert_array_equal(x, array([[[2.0]]]))

    def test_mask(self):
        x = array([1, 2, 3, 4])
        m = array([0, 1], bool)
        assert_array_equal(x[m], array([2]))

    def test_mask2(self):
        x = array([[1, 2, 3, 4], [5, 6, 7, 8]])
        m = array([0, 1], bool)
        m2 = array([[0, 1], [1, 0]], bool)
        m3 = array([[0, 1]], bool)
        assert_array_equal(x[m], array([[5, 6, 7, 8]]))
        assert_array_equal(x[m2], array([2, 5]))
        assert_array_equal(x[m3], array([2]))

    def test_assign_mask(self):
        x = array([1, 2, 3, 4])
        m = array([0, 1], bool)
        x[m] = 5
        assert_array_equal(x, array([1, 5, 3, 4]))

    def test_assign_mask2(self):
        xorig = array([[1, 2, 3, 4], [5, 6, 7, 8]])
        m = array([0, 1], bool)
        m2 = array([[0, 1], [1, 0]], bool)
        m3 = array([[0, 1]], bool)
        x = xorig.copy()
        x[m] = 10
        assert_array_equal(x, array([[1, 2, 3, 4], [10, 10, 10, 10]]))
        x = xorig.copy()
        x[m2] = 10
        assert_array_equal(x, array([[1, 10, 3, 4], [10, 6, 7, 8]]))
        x = xorig.copy()
        x[m3] = 10
        assert_array_equal(x, array([[1, 10, 3, 4], [5, 6, 7, 8]]))


class TestStringCompare(TestCase):
    def test_string(self):
        g1 = array(["This", "is", "example"])
        g2 = array(["This", "was", "example"])
        assert_array_equal(g1 == g2, [g1[i] == g2[i] for i in [0, 1, 2]])
        assert_array_equal(g1 != g2, [g1[i] != g2[i] for i in [0, 1, 2]])
        assert_array_equal(g1 <= g2, [g1[i] <= g2[i] for i in [0, 1, 2]])
        assert_array_equal(g1 >= g2, [g1[i] >= g2[i] for i in [0, 1, 2]])
        assert_array_equal(g1 < g2, [g1[i] < g2[i] for i in [0, 1, 2]])
        assert_array_equal(g1 > g2, [g1[i] > g2[i] for i in [0, 1, 2]])

    def test_mixed(self):
        g1 = array(["spam", "spa", "spammer", "and eggs"])
        g2 = "spam"
        assert_array_equal(g1 == g2, [x == g2 for x in g1])
        assert_array_equal(g1 != g2, [x != g2 for x in g1])
        assert_array_equal(g1 < g2, [x < g2 for x in g1])
        assert_array_equal(g1 > g2, [x > g2 for x in g1])
        assert_array_equal(g1 <= g2, [x <= g2 for x in g1])
        assert_array_equal(g1 >= g2, [x >= g2 for x in g1])


    def test_unicode(self):
        g1 = array([sixu("This"), sixu("is"), sixu("example")])
        g2 = array([sixu("This"), sixu("was"), sixu("example")])
        assert_array_equal(g1 == g2, [g1[i] == g2[i] for i in [0, 1, 2]])
        assert_array_equal(g1 != g2, [g1[i] != g2[i] for i in [0, 1, 2]])
        assert_array_equal(g1 <= g2, [g1[i] <= g2[i] for i in [0, 1, 2]])
        assert_array_equal(g1 >= g2, [g1[i] >= g2[i] for i in [0, 1, 2]])
        assert_array_equal(g1 < g2,  [g1[i] < g2[i] for i in [0, 1, 2]])
        assert_array_equal(g1 > g2,  [g1[i] > g2[i] for i in [0, 1, 2]])


class TestArgmax(TestCase):

    nan_arr = [
        ([0, 1, 2, 3, np.nan], 4),
        ([0, 1, 2, np.nan, 3], 3),
        ([np.nan, 0, 1, 2, 3], 0),
        ([np.nan, 0, np.nan, 2, 3], 0),
        ([0, 1, 2, 3, complex(0, np.nan)], 4),
        ([0, 1, 2, 3, complex(np.nan, 0)], 4),
        ([0, 1, 2, complex(np.nan, 0), 3], 3),
        ([0, 1, 2, complex(0, np.nan), 3], 3),
        ([complex(0, np.nan), 0, 1, 2, 3], 0),
        ([complex(np.nan, np.nan), 0, 1, 2, 3], 0),
        ([complex(np.nan, 0), complex(np.nan, 2), complex(np.nan, 1)], 0),
        ([complex(np.nan, np.nan), complex(np.nan, 2), complex(np.nan, 1)], 0),
        ([complex(np.nan, 0), complex(np.nan, 2), complex(np.nan, np.nan)], 0),

        ([complex(0, 0), complex(0, 2), complex(0, 1)], 1),
        ([complex(1, 0), complex(0, 2), complex(0, 1)], 0),
        ([complex(1, 0), complex(0, 2), complex(1, 1)], 2),

        ([np.datetime64('1923-04-14T12:43:12'),
          np.datetime64('1994-06-21T14:43:15'),
          np.datetime64('2001-10-15T04:10:32'),
          np.datetime64('1995-11-25T16:02:16'),
          np.datetime64('2005-01-04T03:14:12'),
          np.datetime64('2041-12-03T14:05:03')], 5),
        ([np.datetime64('1935-09-14T04:40:11'),
          np.datetime64('1949-10-12T12:32:11'),
          np.datetime64('2010-01-03T05:14:12'),
          np.datetime64('2015-11-20T12:20:59'),
          np.datetime64('1932-09-23T10:10:13'),
          np.datetime64('2014-10-10T03:50:30')], 3),
        ([np.datetime64('2059-03-14T12:43:12'),
          np.datetime64('1996-09-21T14:43:15'),
          np.datetime64('2001-10-15T04:10:32'),
          np.datetime64('2022-12-25T16:02:16'),
          np.datetime64('1963-10-04T03:14:12'),
          np.datetime64('2013-05-08T18:15:23')], 0),

        ([timedelta(days=5, seconds=14), timedelta(days=2, seconds=35),
          timedelta(days=-1, seconds=23)], 0),
        ([timedelta(days=1, seconds=43), timedelta(days=10, seconds=5),
          timedelta(days=5, seconds=14)], 1),
        ([timedelta(days=10, seconds=24), timedelta(days=10, seconds=5),
          timedelta(days=10, seconds=43)], 2),

        ([False, False, False, False, True], 4),
        ([False, False, False, True, False], 3),
        ([True, False, False, False, False], 0),
        ([True, False, True, False, False], 0),

        # Can't reduce a "flexible type"
        #(['a', 'z', 'aa', 'zz'], 3),
        #(['zz', 'a', 'aa', 'a'], 0),
        #(['aa', 'z', 'zz', 'a'], 2),
    ]

    def test_all(self):
        a = np.random.normal(0, 1, (4, 5, 6, 7, 8))
        for i in range(a.ndim):
            amax = a.max(i)
            aargmax = a.argmax(i)
            axes = list(range(a.ndim))
            axes.remove(i)
            assert_(all(amax == aargmax.choose(*a.transpose(i,*axes))))

    def test_combinations(self):
        for arr, pos in self.nan_arr:
            assert_equal(np.argmax(arr), pos, err_msg="%r"%arr)
            assert_equal(arr[np.argmax(arr)], np.max(arr), err_msg="%r"%arr)

    def test_output_shape(self):
        # see also gh-616
        a = np.ones((10, 5))
        # Check some simple shape mismatches
        out = np.ones(11, dtype=np.int_)
        assert_raises(ValueError, a.argmax, -1, out)

        out = np.ones((2, 5), dtype=np.int_)
        assert_raises(ValueError, a.argmax, -1, out)

        # these could be relaxed possibly (used to allow even the previous)
        out = np.ones((1, 10), dtype=np.int_)
        assert_raises(ValueError, a.argmax, -1, np.ones((1, 10)))

        out = np.ones(10, dtype=np.int_)
        a.argmax(-1, out=out)
        assert_equal(out, a.argmax(-1))


class TestArgmin(TestCase):

    nan_arr = [
        ([0, 1, 2, 3, np.nan], 4),
        ([0, 1, 2, np.nan, 3], 3),
        ([np.nan, 0, 1, 2, 3], 0),
        ([np.nan, 0, np.nan, 2, 3], 0),
        ([0, 1, 2, 3, complex(0, np.nan)], 4),
        ([0, 1, 2, 3, complex(np.nan, 0)], 4),
        ([0, 1, 2, complex(np.nan, 0), 3], 3),
        ([0, 1, 2, complex(0, np.nan), 3], 3),
        ([complex(0, np.nan), 0, 1, 2, 3], 0),
        ([complex(np.nan, np.nan), 0, 1, 2, 3], 0),
        ([complex(np.nan, 0), complex(np.nan, 2), complex(np.nan, 1)], 0),
        ([complex(np.nan, np.nan), complex(np.nan, 2), complex(np.nan, 1)], 0),
        ([complex(np.nan, 0), complex(np.nan, 2), complex(np.nan, np.nan)], 0),

        ([complex(0, 0), complex(0, 2), complex(0, 1)], 0),
        ([complex(1, 0), complex(0, 2), complex(0, 1)], 2),
        ([complex(1, 0), complex(0, 2), complex(1, 1)], 1),

        ([np.datetime64('1923-04-14T12:43:12'),
          np.datetime64('1994-06-21T14:43:15'),
          np.datetime64('2001-10-15T04:10:32'),
          np.datetime64('1995-11-25T16:02:16'),
          np.datetime64('2005-01-04T03:14:12'),
          np.datetime64('2041-12-03T14:05:03')], 0),
        ([np.datetime64('1935-09-14T04:40:11'),
          np.datetime64('1949-10-12T12:32:11'),
          np.datetime64('2010-01-03T05:14:12'),
          np.datetime64('2014-11-20T12:20:59'),
          np.datetime64('2015-09-23T10:10:13'),
          np.datetime64('1932-10-10T03:50:30')], 5),
        ([np.datetime64('2059-03-14T12:43:12'),
          np.datetime64('1996-09-21T14:43:15'),
          np.datetime64('2001-10-15T04:10:32'),
          np.datetime64('2022-12-25T16:02:16'),
          np.datetime64('1963-10-04T03:14:12'),
          np.datetime64('2013-05-08T18:15:23')], 4),

        ([timedelta(days=5, seconds=14), timedelta(days=2, seconds=35),
          timedelta(days=-1, seconds=23)], 2),
        ([timedelta(days=1, seconds=43), timedelta(days=10, seconds=5),
          timedelta(days=5, seconds=14)], 0),
        ([timedelta(days=10, seconds=24), timedelta(days=10, seconds=5),
          timedelta(days=10, seconds=43)], 1),

        ([True, True, True, True, False], 4),
        ([True, True, True, False, True], 3),
        ([False, True, True, True, True], 0),
        ([False, True, False, True, True], 0),

        # Can't reduce a "flexible type"
        #(['a', 'z', 'aa', 'zz'], 0),
        #(['zz', 'a', 'aa', 'a'], 1),
        #(['aa', 'z', 'zz', 'a'], 3),
    ]

    def test_all(self):
        a = np.random.normal(0, 1, (4, 5, 6, 7, 8))
        for i in range(a.ndim):
            amin = a.min(i)
            aargmin = a.argmin(i)
            axes = list(range(a.ndim))
            axes.remove(i)
            assert_(all(amin == aargmin.choose(*a.transpose(i,*axes))))

    def test_combinations(self):
        for arr, pos in self.nan_arr:
            assert_equal(np.argmin(arr), pos, err_msg="%r"%arr)
            assert_equal(arr[np.argmin(arr)], np.min(arr), err_msg="%r"%arr)

    def test_minimum_signed_integers(self):

        a = np.array([1, -2**7, -2**7 + 1], dtype=np.int8)
        assert_equal(np.argmin(a), 1)

        a = np.array([1, -2**15, -2**15 + 1], dtype=np.int16)
        assert_equal(np.argmin(a), 1)

        a = np.array([1, -2**31, -2**31 + 1], dtype=np.int32)
        assert_equal(np.argmin(a), 1)

        a = np.array([1, -2**63, -2**63 + 1], dtype=np.int64)
        assert_equal(np.argmin(a), 1)

    def test_output_shape(self):
        # see also gh-616
        a = np.ones((10, 5))
        # Check some simple shape mismatches
        out = np.ones(11, dtype=np.int_)
        assert_raises(ValueError, a.argmin, -1, out)

        out = np.ones((2, 5), dtype=np.int_)
        assert_raises(ValueError, a.argmin, -1, out)

        # these could be relaxed possibly (used to allow even the previous)
        out = np.ones((1, 10), dtype=np.int_)
        assert_raises(ValueError, a.argmin, -1, np.ones((1, 10)))

        out = np.ones(10, dtype=np.int_)
        a.argmin(-1, out=out)
        assert_equal(out, a.argmin(-1))


class TestMinMax(TestCase):
    def test_scalar(self):
        assert_raises(ValueError, np.amax, 1, 1)
        assert_raises(ValueError, np.amin, 1, 1)

        assert_equal(np.amax(1, axis=0), 1)
        assert_equal(np.amin(1, axis=0), 1)
        assert_equal(np.amax(1, axis=None), 1)
        assert_equal(np.amin(1, axis=None), 1)

    def test_axis(self):
        assert_raises(ValueError, np.amax, [1, 2, 3], 1000)
        assert_equal(np.amax([[1, 2, 3]], axis=1), 3)

class TestNewaxis(TestCase):
    def test_basic(self):
        sk = array([0, -0.1, 0.1])
        res = 250*sk[:, newaxis]
        assert_almost_equal(res.ravel(), 250*sk)


class TestClip(TestCase):
    def _check_range(self, x, cmin, cmax):
        assert_(np.all(x >= cmin))
        assert_(np.all(x <= cmax))

    def _clip_type(self,type_group,array_max,
                   clip_min,clip_max,inplace=False,
                   expected_min=None,expected_max=None):
        if expected_min is None:
            expected_min = clip_min
        if expected_max is None:
            expected_max = clip_max

        for T in np.sctypes[type_group]:
            if sys.byteorder == 'little':
                byte_orders = ['=', '>']
            else:
                byte_orders = ['<', '=']

            for byteorder in byte_orders:
                dtype = np.dtype(T).newbyteorder(byteorder)

                x = (np.random.random(1000) * array_max).astype(dtype)
                if inplace:
                    x.clip(clip_min, clip_max, x)
                else:
                    x = x.clip(clip_min, clip_max)
                    byteorder = '='

                if x.dtype.byteorder == '|': byteorder = '|'
                assert_equal(x.dtype.byteorder, byteorder)
                self._check_range(x, expected_min, expected_max)
        return x

    def test_basic(self):
        for inplace in [False, True]:
            self._clip_type('float', 1024, -12.8, 100.2, inplace=inplace)
            self._clip_type('float', 1024, 0, 0, inplace=inplace)

            self._clip_type('int', 1024, -120, 100.5, inplace=inplace)
            self._clip_type('int', 1024, 0, 0, inplace=inplace)

            x = self._clip_type('uint', 1024, -120, 100, expected_min=0,
                    inplace=inplace)
            x = self._clip_type('uint', 1024, 0, 0, inplace=inplace)

    def test_record_array(self):
        rec = np.array([(-5, 2.0, 3.0), (5.0, 4.0, 3.0)],
                      dtype=[('x', '<f8'), ('y', '<f8'), ('z', '<f8')])
        y = rec['x'].clip(-0.3, 0.5)
        self._check_range(y, -0.3, 0.5)

    def test_max_or_min(self):
        val = np.array([0, 1, 2, 3, 4, 5, 6, 7])
        x = val.clip(3)
        assert_(np.all(x >= 3))
        x = val.clip(min=3)
        assert_(np.all(x >= 3))
        x = val.clip(max=4)
        assert_(np.all(x <= 4))


class TestPutmask(object):
    def tst_basic(self, x, T, mask, val):
        np.putmask(x, mask, val)
        assert_(np.all(x[mask] == T(val)))
        assert_(x.dtype == T)

    def test_ip_types(self):
        unchecked_types = [str, unicode, np.void, object]

        x = np.random.random(1000)*100
        mask = x < 40

        for val in [-100, 0, 15]:
            for types in np.sctypes.values():
                for T in types:
                    if T not in unchecked_types:
                        yield self.tst_basic, x.copy().astype(T), T, mask, val

    def test_mask_size(self):
        assert_raises(ValueError, np.putmask, np.array([1, 2, 3]), [True], 5)

    def tst_byteorder(self, dtype):
        x = np.array([1, 2, 3], dtype)
        np.putmask(x, [True, False, True], -1)
        assert_array_equal(x, [-1, 2, -1])

    def test_ip_byteorder(self):
        for dtype in ('>i4', '<i4'):
            yield self.tst_byteorder, dtype

    def test_record_array(self):
        # Note mixed byteorder.
        rec = np.array([(-5, 2.0, 3.0), (5.0, 4.0, 3.0)],
                      dtype=[('x', '<f8'), ('y', '>f8'), ('z', '<f8')])
        np.putmask(rec['x'], [True, False], 10)
        assert_array_equal(rec['x'], [10, 5])
        assert_array_equal(rec['y'], [2, 4])
        assert_array_equal(rec['z'], [3, 3])
        np.putmask(rec['y'], [True, False], 11)
        assert_array_equal(rec['x'], [10, 5])
        assert_array_equal(rec['y'], [11, 4])
        assert_array_equal(rec['z'], [3, 3])

    def test_masked_array(self):
        ## x = np.array([1,2,3])
        ## z = np.ma.array(x,mask=[True,False,False])
        ## np.putmask(z,[True,True,True],3)
        pass


class TestTake(object):
    def tst_basic(self, x):
        ind = list(range(x.shape[0]))
        assert_array_equal(x.take(ind, axis=0), x)

    def test_ip_types(self):
        unchecked_types = [str, unicode, np.void, object]

        x = np.random.random(24)*100
        x.shape = 2, 3, 4
        for types in np.sctypes.values():
            for T in types:
                if T not in unchecked_types:
                    yield self.tst_basic, x.copy().astype(T)

    def test_raise(self):
        x = np.random.random(24)*100
        x.shape = 2, 3, 4
        assert_raises(IndexError, x.take, [0, 1, 2], axis=0)
        assert_raises(IndexError, x.take, [-3], axis=0)
        assert_array_equal(x.take([-1], axis=0)[0], x[1])

    def test_clip(self):
        x = np.random.random(24)*100
        x.shape = 2, 3, 4
        assert_array_equal(x.take([-1], axis=0, mode='clip')[0], x[0])
        assert_array_equal(x.take([2], axis=0, mode='clip')[0], x[1])

    def test_wrap(self):
        x = np.random.random(24)*100
        x.shape = 2, 3, 4
        assert_array_equal(x.take([-1], axis=0, mode='wrap')[0], x[1])
        assert_array_equal(x.take([2], axis=0, mode='wrap')[0], x[0])
        assert_array_equal(x.take([3], axis=0, mode='wrap')[0], x[1])

    def tst_byteorder(self, dtype):
        x = np.array([1, 2, 3], dtype)
        assert_array_equal(x.take([0, 2, 1]), [1, 3, 2])

    def test_ip_byteorder(self):
        for dtype in ('>i4', '<i4'):
            yield self.tst_byteorder, dtype

    def test_record_array(self):
        # Note mixed byteorder.
        rec = np.array([(-5, 2.0, 3.0), (5.0, 4.0, 3.0)],
                      dtype=[('x', '<f8'), ('y', '>f8'), ('z', '<f8')])
        rec1 = rec.take([1])
        assert_(rec1['x'] == 5.0 and rec1['y'] == 4.0)


class TestLexsort(TestCase):
    def test_basic(self):
        a = [1, 2, 1, 3, 1, 5]
        b = [0, 4, 5, 6, 2, 3]
        idx = np.lexsort((b, a))
        expected_idx = np.array([0, 4, 2, 1, 3, 5])
        assert_array_equal(idx, expected_idx)

        x = np.vstack((b, a))
        idx = np.lexsort(x)
        assert_array_equal(idx, expected_idx)

        assert_array_equal(x[1][idx], np.sort(x[1]))

    def test_datetime(self):
        a = np.array([0,0,0], dtype='datetime64[D]')
        b = np.array([2,1,0], dtype='datetime64[D]')
        idx = np.lexsort((b, a))
        expected_idx = np.array([2, 1, 0])
        assert_array_equal(idx, expected_idx)

        a = np.array([0,0,0], dtype='timedelta64[D]')
        b = np.array([2,1,0], dtype='timedelta64[D]')
        idx = np.lexsort((b, a))
        expected_idx = np.array([2, 1, 0])
        assert_array_equal(idx, expected_idx)


class TestIO(object):
    """Test tofile, fromfile, tobytes, and fromstring"""

    def setUp(self):
        shape = (2, 4, 3)
        rand = np.random.random
        self.x = rand(shape) + rand(shape).astype(np.complex)*1j
        self.x[0,:, 1] = [nan, inf, -inf, nan]
        self.dtype = self.x.dtype
        self.tempdir = tempfile.mkdtemp()
        self.filename = tempfile.mktemp(dir=self.tempdir)

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_bool_fromstring(self):
        v = np.array([True, False, True, False], dtype=np.bool_)
        y = np.fromstring('1 0 -2.3 0.0', sep=' ', dtype=np.bool_)
        assert_array_equal(v, y)

    def test_uint64_fromstring(self):
        d = np.fromstring("9923372036854775807 104783749223640",
                          dtype=np.uint64, sep=' ');
        e = np.array([9923372036854775807, 104783749223640], dtype=np.uint64)
        assert_array_equal(d, e)

    def test_int64_fromstring(self):
        d = np.fromstring("-25041670086757 104783749223640",
                          dtype=np.int64, sep=' ');
        e = np.array([-25041670086757, 104783749223640], dtype=np.int64)
        assert_array_equal(d, e)

    def test_empty_files_binary(self):
        f = open(self.filename, 'w')
        f.close()
        y = fromfile(self.filename)
        assert_(y.size == 0, "Array not empty")

    def test_empty_files_text(self):
        f = open(self.filename, 'w')
        f.close()
        y = fromfile(self.filename, sep=" ")
        assert_(y.size == 0, "Array not empty")

    def test_roundtrip_file(self):
        f = open(self.filename, 'wb')
        self.x.tofile(f)
        f.close()
        # NB. doesn't work with flush+seek, due to use of C stdio
        f = open(self.filename, 'rb')
        y = np.fromfile(f, dtype=self.dtype)
        f.close()
        assert_array_equal(y, self.x.flat)

    def test_roundtrip_filename(self):
        self.x.tofile(self.filename)
        y = np.fromfile(self.filename, dtype=self.dtype)
        assert_array_equal(y, self.x.flat)

    def test_roundtrip_binary_str(self):
        s = self.x.tobytes()
        y = np.fromstring(s, dtype=self.dtype)
        assert_array_equal(y, self.x.flat)

        s = self.x.tobytes('F')
        y = np.fromstring(s, dtype=self.dtype)
        assert_array_equal(y, self.x.flatten('F'))

    def test_roundtrip_str(self):
        x = self.x.real.ravel()
        s = "@".join(map(str, x))
        y = np.fromstring(s, sep="@")
        # NB. str imbues less precision
        nan_mask = ~np.isfinite(x)
        assert_array_equal(x[nan_mask], y[nan_mask])
        assert_array_almost_equal(x[~nan_mask], y[~nan_mask], decimal=5)

    def test_roundtrip_repr(self):
        x = self.x.real.ravel()
        s = "@".join(map(repr, x))
        y = np.fromstring(s, sep="@")
        assert_array_equal(x, y)

    def test_file_position_after_fromfile(self):
        # gh-4118
        sizes = [io.DEFAULT_BUFFER_SIZE//8,
                 io.DEFAULT_BUFFER_SIZE,
                 io.DEFAULT_BUFFER_SIZE*8]

        for size in sizes:
            f = open(self.filename, 'wb')
            f.seek(size-1)
            f.write(b'\0')
            f.close()

            for mode in ['rb', 'r+b']:
                err_msg = "%d %s" % (size, mode)

                f = open(self.filename, mode)
                f.read(2)
                np.fromfile(f, dtype=np.float64, count=1)
                pos = f.tell()
                f.close()
                assert_equal(pos, 10, err_msg=err_msg)

    def test_file_position_after_tofile(self):
        # gh-4118
        sizes = [io.DEFAULT_BUFFER_SIZE//8,
                 io.DEFAULT_BUFFER_SIZE,
                 io.DEFAULT_BUFFER_SIZE*8]

        for size in sizes:
            err_msg = "%d" % (size,)

            f = open(self.filename, 'wb')
            f.seek(size-1)
            f.write(b'\0')
            f.seek(10)
            f.write(b'12')
            np.array([0], dtype=np.float64).tofile(f)
            pos = f.tell()
            f.close()
            assert_equal(pos, 10 + 2 + 8, err_msg=err_msg)

            f = open(self.filename, 'r+b')
            f.read(2)
            f.seek(0, 1) # seek between read&write required by ANSI C
            np.array([0], dtype=np.float64).tofile(f)
            pos = f.tell()
            f.close()
            assert_equal(pos, 10, err_msg=err_msg)

    def _check_from(self, s, value, **kw):
        y = np.fromstring(asbytes(s), **kw)
        assert_array_equal(y, value)

        f = open(self.filename, 'wb')
        f.write(asbytes(s))
        f.close()
        y = np.fromfile(self.filename, **kw)
        assert_array_equal(y, value)

    def test_nan(self):
        self._check_from("nan +nan -nan NaN nan(foo) +NaN(BAR) -NAN(q_u_u_x_)",
                         [nan, nan, nan, nan, nan, nan, nan],
                         sep=' ')

    def test_inf(self):
        self._check_from("inf +inf -inf infinity -Infinity iNfInItY -inF",
                         [inf, inf, -inf, inf, -inf, inf, -inf], sep=' ')

    def test_numbers(self):
        self._check_from("1.234 -1.234 .3 .3e55 -123133.1231e+133",
                         [1.234, -1.234, .3, .3e55, -123133.1231e+133], sep=' ')

    def test_binary(self):
        self._check_from('\x00\x00\x80?\x00\x00\x00@\x00\x00@@\x00\x00\x80@',
                         array([1, 2, 3, 4]),
                         dtype='<f4')

    @dec.slow # takes > 1 minute on mechanical hard drive
    def test_big_binary(self):
        """Test workarounds for 32-bit limited fwrite, fseek, and ftell
        calls in windows. These normally would hang doing something like this.
        See http://projects.scipy.org/numpy/ticket/1660"""
        if sys.platform != 'win32':
            return
        try:
            # before workarounds, only up to 2**32-1 worked
            fourgbplus = 2**32 + 2**16
            testbytes = np.arange(8, dtype=np.int8)
            n = len(testbytes)
            flike = tempfile.NamedTemporaryFile()
            f = flike.file
            np.tile(testbytes, fourgbplus // testbytes.nbytes).tofile(f)
            flike.seek(0)
            a = np.fromfile(f, dtype=np.int8)
            flike.close()
            assert_(len(a) == fourgbplus)
            # check only start and end for speed:
            assert_((a[:n] == testbytes).all())
            assert_((a[-n:] == testbytes).all())
        except (MemoryError, ValueError):
            pass

    def test_string(self):
        self._check_from('1,2,3,4', [1., 2., 3., 4.], sep=',')

    def test_counted_string(self):
        self._check_from('1,2,3,4', [1., 2., 3., 4.], count=4, sep=',')
        self._check_from('1,2,3,4', [1., 2., 3.], count=3, sep=',')
        self._check_from('1,2,3,4', [1., 2., 3., 4.], count=-1, sep=',')

    def test_string_with_ws(self):
        self._check_from('1 2  3     4   ', [1, 2, 3, 4], dtype=int, sep=' ')

    def test_counted_string_with_ws(self):
        self._check_from('1 2  3     4   ', [1, 2, 3], count=3, dtype=int,
                         sep=' ')

    def test_ascii(self):
        self._check_from('1 , 2 , 3 , 4', [1., 2., 3., 4.], sep=',')
        self._check_from('1,2,3,4', [1., 2., 3., 4.], dtype=float, sep=',')

    def test_malformed(self):
        self._check_from('1.234 1,234', [1.234, 1.], sep=' ')

    def test_long_sep(self):
        self._check_from('1_x_3_x_4_x_5', [1, 3, 4, 5], sep='_x_')

    def test_dtype(self):
        v = np.array([1, 2, 3, 4], dtype=np.int_)
        self._check_from('1,2,3,4', v, sep=',', dtype=np.int_)

    def test_dtype_bool(self):
        # can't use _check_from because fromstring can't handle True/False
        v = np.array([True, False, True, False], dtype=np.bool_)
        s = '1,0,-2.3,0'
        f = open(self.filename, 'wb')
        f.write(asbytes(s))
        f.close()
        y = np.fromfile(self.filename, sep=',', dtype=np.bool_)
        assert_(y.dtype == '?')
        assert_array_equal(y, v)

    def test_tofile_sep(self):
        x = np.array([1.51, 2, 3.51, 4], dtype=float)
        f = open(self.filename, 'w')
        x.tofile(f, sep=',')
        f.close()
        f = open(self.filename, 'r')
        s = f.read()
        f.close()
        assert_equal(s, '1.51,2.0,3.51,4.0')

    def test_tofile_format(self):
        x = np.array([1.51, 2, 3.51, 4], dtype=float)
        f = open(self.filename, 'w')
        x.tofile(f, sep=',', format='%.2f')
        f.close()
        f = open(self.filename, 'r')
        s = f.read()
        f.close()
        assert_equal(s, '1.51,2.00,3.51,4.00')

    def test_locale(self):
        in_foreign_locale(self.test_numbers)()
        in_foreign_locale(self.test_nan)()
        in_foreign_locale(self.test_inf)()
        in_foreign_locale(self.test_counted_string)()
        in_foreign_locale(self.test_ascii)()
        in_foreign_locale(self.test_malformed)()
        in_foreign_locale(self.test_tofile_sep)()
        in_foreign_locale(self.test_tofile_format)()


class TestFromBuffer(object):
    def tst_basic(self, buffer, expected, kwargs):
        assert_array_equal(np.frombuffer(buffer,**kwargs), expected)

    def test_ip_basic(self):
        for byteorder in ['<', '>']:
            for dtype in [float, int, np.complex]:
                dt = np.dtype(dtype).newbyteorder(byteorder)
                x = (np.random.random((4, 7))*5).astype(dt)
                buf = x.tobytes()
                yield self.tst_basic, buf, x.flat, {'dtype':dt}

    def test_empty(self):
        yield self.tst_basic, asbytes(''), np.array([]), {}


class TestFlat(TestCase):
    def setUp(self):
        a0 = arange(20.0)
        a = a0.reshape(4, 5)
        a0.shape = (4, 5)
        a.flags.writeable = False
        self.a = a
        self.b = a[::2, ::2]
        self.a0 = a0
        self.b0 = a0[::2, ::2]

    def test_contiguous(self):
        testpassed = False
        try:
            self.a.flat[12] = 100.0
        except ValueError:
            testpassed = True
        assert testpassed
        assert self.a.flat[12] == 12.0

    def test_discontiguous(self):
        testpassed = False
        try:
            self.b.flat[4] = 100.0
        except ValueError:
            testpassed = True
        assert testpassed
        assert self.b.flat[4] == 12.0

    def test___array__(self):
        c = self.a.flat.__array__()
        d = self.b.flat.__array__()
        e = self.a0.flat.__array__()
        f = self.b0.flat.__array__()

        assert c.flags.writeable is False
        assert d.flags.writeable is False
        assert e.flags.writeable is True
        assert f.flags.writeable is True

        assert c.flags.updateifcopy is False
        assert d.flags.updateifcopy is False
        assert e.flags.updateifcopy is False
        assert f.flags.updateifcopy is True
        assert f.base is self.b0

class TestResize(TestCase):
    def test_basic(self):
        x = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        x.resize((5, 5))
        assert_array_equal(x.flat[:9],
                np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]).flat)
        assert_array_equal(x[9:].flat, 0)

    def test_check_reference(self):
        x = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        y = x
        self.assertRaises(ValueError, x.resize, (5, 1))

    def test_int_shape(self):
        x = np.eye(3)
        x.resize(3)
        assert_array_equal(x, np.eye(3)[0,:])

    def test_none_shape(self):
        x = np.eye(3)
        x.resize(None)
        assert_array_equal(x, np.eye(3))
        x.resize()
        assert_array_equal(x, np.eye(3))

    def test_invalid_arguements(self):
        self.assertRaises(TypeError, np.eye(3).resize, 'hi')
        self.assertRaises(ValueError, np.eye(3).resize, -1)
        self.assertRaises(TypeError, np.eye(3).resize, order=1)
        self.assertRaises(TypeError, np.eye(3).resize, refcheck='hi')

    def test_freeform_shape(self):
        x = np.eye(3)
        x.resize(3, 2, 1)
        assert_(x.shape == (3, 2, 1))

    def test_zeros_appended(self):
        x = np.eye(3)
        x.resize(2, 3, 3)
        assert_array_equal(x[0], np.eye(3))
        assert_array_equal(x[1], np.zeros((3, 3)))

    def test_obj_obj(self):
        # check memory is initialized on resize, gh-4857
        a = ones(10, dtype=[('k', object, 2)])
        a.resize(15,)
        assert_equal(a.shape, (15,))
        assert_array_equal(a['k'][-5:], 0)
        assert_array_equal(a['k'][:-5], 1)


class TestRecord(TestCase):
    def test_field_rename(self):
        dt = np.dtype([('f', float), ('i', int)])
        dt.names = ['p', 'q']
        assert_equal(dt.names, ['p', 'q'])

    if sys.version_info[0] >= 3:
        def test_bytes_fields(self):
            # Bytes are not allowed in field names and not recognized in titles
            # on Py3
            assert_raises(TypeError, np.dtype, [(asbytes('a'), int)])
            assert_raises(TypeError, np.dtype, [(('b', asbytes('a')), int)])

            dt = np.dtype([((asbytes('a'), 'b'), int)])
            assert_raises(ValueError, dt.__getitem__, asbytes('a'))

            x = np.array([(1,), (2,), (3,)], dtype=dt)
            assert_raises(ValueError, x.__getitem__, asbytes('a'))

            y = x[0]
            assert_raises(IndexError, y.__getitem__, asbytes('a'))
    else:
        def test_unicode_field_titles(self):
            # Unicode field titles are added to field dict on Py2
            title = unicode('b')
            dt = np.dtype([((title, 'a'), int)])
            dt[title]
            dt['a']
            x = np.array([(1,), (2,), (3,)], dtype=dt)
            x[title]
            x['a']
            y = x[0]
            y[title]
            y['a']

        def test_unicode_field_names(self):
            # Unicode field names are not allowed on Py2
            title = unicode('b')
            assert_raises(TypeError, np.dtype, [(title, int)])
            assert_raises(TypeError, np.dtype, [(('a', title), int)])

    def test_field_names(self):
        # Test unicode and 8-bit / byte strings can be used
        a = np.zeros((1,), dtype=[('f1', 'i4'),
                                  ('f2', 'i4'),
                                  ('f3', [('sf1', 'i4')])])
        is_py3 = sys.version_info[0] >= 3
        if is_py3:
            funcs = (str,)
            # byte string indexing fails gracefully
            assert_raises(ValueError, a.__setitem__, asbytes('f1'), 1)
            assert_raises(ValueError, a.__getitem__, asbytes('f1'))
            assert_raises(IndexError, a['f1'].__setitem__, asbytes('sf1'), 1)
            assert_raises(IndexError, a['f1'].__getitem__, asbytes('sf1'))
        else:
            funcs = (str, unicode)
        for func in funcs:
            b = a.copy()
            fn1 = func('f1')
            b[fn1] = 1
            assert_equal(b[fn1], 1)
            fnn = func('not at all')
            assert_raises(ValueError, b.__setitem__, fnn, 1)
            assert_raises(ValueError, b.__getitem__, fnn)
            b[0][fn1] = 2
            assert_equal(b[fn1], 2)
            # Subfield
            assert_raises(IndexError, b[0].__setitem__, fnn, 1)
            assert_raises(IndexError, b[0].__getitem__, fnn)
            # Subfield
            fn3 = func('f3')
            sfn1 = func('sf1')
            b[fn3][sfn1] = 1
            assert_equal(b[fn3][sfn1], 1)
            assert_raises(ValueError, b[fn3].__setitem__, fnn, 1)
            assert_raises(ValueError, b[fn3].__getitem__, fnn)
            # multiple Subfields
            fn2 = func('f2')
            b[fn2] = 3
            assert_equal(b[['f1', 'f2']][0].tolist(), (2, 3))
            assert_equal(b[['f2', 'f1']][0].tolist(), (3, 2))
            assert_equal(b[['f1', 'f3']][0].tolist(), (2, (1,)))
            # view of subfield view/copy
            assert_equal(b[['f1', 'f2']][0].view(('i4', 2)).tolist(), (2, 3))
            assert_equal(b[['f2', 'f1']][0].view(('i4', 2)).tolist(), (3, 2))
            view_dtype=[('f1', 'i4'), ('f3', [('', 'i4')])]
            assert_equal(b[['f1', 'f3']][0].view(view_dtype).tolist(), (2, (1,)))
        # non-ascii unicode field indexing is well behaved
        if not is_py3:
            raise SkipTest('non ascii unicode field indexing skipped; '
                           'raises segfault on python 2.x')
        else:
            assert_raises(ValueError, a.__setitem__, sixu('\u03e0'), 1)
            assert_raises(ValueError, a.__getitem__, sixu('\u03e0'))

    def test_field_names_deprecation(self):

        def collect_warning_types(f, *args, **kwargs):
            with warnings.catch_warnings(record=True) as log:
                warnings.simplefilter("always")
                f(*args, **kwargs)
            return [w.category for w in log]

        a = np.zeros((1,), dtype=[('f1', 'i4'),
                                  ('f2', 'i4'),
                                  ('f3', [('sf1', 'i4')])])
        a['f1'][0] = 1
        a['f2'][0] = 2
        a['f3'][0] = (3,)
        b = np.zeros((1,), dtype=[('f1', 'i4'),
                                  ('f2', 'i4'),
                                  ('f3', [('sf1', 'i4')])])
        b['f1'][0] = 1
        b['f2'][0] = 2
        b['f3'][0] = (3,)

        # All the different functions raise a warning, but not an error, and
        # 'a' is not modified:
        assert_equal(collect_warning_types(a[['f1', 'f2']].__setitem__, 0, (10, 20)),
                     [FutureWarning])
        assert_equal(a, b)
        # Views also warn
        subset = a[['f1', 'f2']]
        subset_view = subset.view()
        assert_equal(collect_warning_types(subset_view['f1'].__setitem__, 0, 10),
                     [FutureWarning])
        # But the write goes through:
        assert_equal(subset['f1'][0], 10)
        # Only one warning per multiple field indexing, though (even if there are
        # multiple views involved):
        assert_equal(collect_warning_types(subset['f1'].__setitem__, 0, 10),
                     [])

    def test_record_hash(self):
        a = np.array([(1, 2), (1, 2)], dtype='i1,i2')
        a.flags.writeable = False
        b = np.array([(1, 2), (3, 4)], dtype=[('num1', 'i1'), ('num2', 'i2')])
        b.flags.writeable = False
        c = np.array([(1, 2), (3, 4)], dtype='i1,i2')
        c.flags.writeable = False
        self.assertTrue(hash(a[0]) == hash(a[1]))
        self.assertTrue(hash(a[0]) == hash(b[0]))
        self.assertTrue(hash(a[0]) != hash(b[1]))
        self.assertTrue(hash(c[0]) == hash(a[0]) and c[0] == a[0])

    def test_record_no_hash(self):
        a = np.array([(1, 2), (1, 2)], dtype='i1,i2')
        self.assertRaises(TypeError, hash, a[0])

class TestView(TestCase):
    def test_basic(self):
        x = np.array([(1, 2, 3, 4), (5, 6, 7, 8)], dtype=[('r', np.int8), ('g', np.int8),
                                                  ('b', np.int8), ('a', np.int8)])
        # We must be specific about the endianness here:
        y = x.view(dtype='<i4')
        # ... and again without the keyword.
        z = x.view('<i4')
        assert_array_equal(y, z)
        assert_array_equal(y, [67305985, 134678021])

def _mean(a, **args):
    return a.mean(**args)

def _var(a, **args):
    return a.var(**args)

def _std(a, **args):
    return a.std(**args)

class TestStats(TestCase):

    funcs = [_mean, _var, _std]

    def setUp(self):
        np.random.seed(range(3))
        self.rmat = np.random.random((4, 5))
        self.cmat = self.rmat + 1j * self.rmat
        self.omat = np.array([Decimal(repr(r)) for r in self.rmat.flat])
        self.omat = self.omat.reshape(4, 5)

    def test_keepdims(self):
        mat = np.eye(3)
        for f in self.funcs:
            for axis in [0, 1]:
                res = f(mat, axis=axis, keepdims=True)
                assert_(res.ndim == mat.ndim)
                assert_(res.shape[axis] == 1)
            for axis in [None]:
                res = f(mat, axis=axis, keepdims=True)
                assert_(res.shape == (1, 1))

    def test_out(self):
        mat = np.eye(3)
        for f in self.funcs:
            out = np.zeros(3)
            tgt = f(mat, axis=1)
            res = f(mat, axis=1, out=out)
            assert_almost_equal(res, out)
            assert_almost_equal(res, tgt)
        out = np.empty(2)
        assert_raises(ValueError, f, mat, axis=1, out=out)
        out = np.empty((2, 2))
        assert_raises(ValueError, f, mat, axis=1, out=out)

    def test_dtype_from_input(self):

        icodes = np.typecodes['AllInteger']
        fcodes = np.typecodes['AllFloat']

        # object type
        for f in self.funcs:
            mat = np.array([[Decimal(1)]*3]*3)
            tgt = mat.dtype.type
            res = f(mat, axis=1).dtype.type
            assert_(res is tgt)
            # scalar case
            res = type(f(mat, axis=None))
            assert_(res is Decimal)

        # integer types
        for f in self.funcs:
            for c in icodes:
                mat = np.eye(3, dtype=c)
                tgt = np.float64
                res = f(mat, axis=1).dtype.type
                assert_(res is tgt)
                # scalar case
                res = f(mat, axis=None).dtype.type
                assert_(res is tgt)

        # mean for float types
        for f in [_mean]:
            for c in fcodes:
                mat = np.eye(3, dtype=c)
                tgt = mat.dtype.type
                res = f(mat, axis=1).dtype.type
                assert_(res is tgt)
                # scalar case
                res = f(mat, axis=None).dtype.type
                assert_(res is tgt)

        # var, std for float types
        for f in [_var, _std]:
            for c in fcodes:
                mat = np.eye(3, dtype=c)
                # deal with complex types
                tgt = mat.real.dtype.type
                res = f(mat, axis=1).dtype.type
                assert_(res is tgt)
                # scalar case
                res = f(mat, axis=None).dtype.type
                assert_(res is tgt)

    def test_dtype_from_dtype(self):
        icodes = np.typecodes['AllInteger']
        fcodes = np.typecodes['AllFloat']
        mat = np.eye(3)

        # stats for integer types
        # fixme:
        # this needs definition as there are lots places along the line
        # where type casting may take place.
        #for f in self.funcs:
            #for c in icodes:
                #tgt = np.dtype(c).type
                #res = f(mat, axis=1, dtype=c).dtype.type
                #assert_(res is tgt)
                ## scalar case
                #res = f(mat, axis=None, dtype=c).dtype.type
                #assert_(res is tgt)

        # stats for float types
        for f in self.funcs:
            for c in fcodes:
                tgt = np.dtype(c).type
                res = f(mat, axis=1, dtype=c).dtype.type
                assert_(res is tgt)
                # scalar case
                res = f(mat, axis=None, dtype=c).dtype.type
                assert_(res is tgt)

    def test_ddof(self):
        for f in [_var]:
            for ddof in range(3):
                dim = self.rmat.shape[1]
                tgt = f(self.rmat, axis=1) * dim
                res = f(self.rmat, axis=1, ddof=ddof) * (dim - ddof)
        for f in [_std]:
            for ddof in range(3):
                dim = self.rmat.shape[1]
                tgt = f(self.rmat, axis=1) * np.sqrt(dim)
                res = f(self.rmat, axis=1, ddof=ddof) * np.sqrt(dim - ddof)
                assert_almost_equal(res, tgt)
                assert_almost_equal(res, tgt)

    def test_ddof_too_big(self):
        dim = self.rmat.shape[1]
        for f in [_var, _std]:
            for ddof in range(dim, dim + 2):
                with warnings.catch_warnings(record=True) as w:
                    warnings.simplefilter('always')
                    res = f(self.rmat, axis=1, ddof=ddof)
                    assert_(not (res < 0).any())
                    assert_(len(w) > 0)
                    assert_(issubclass(w[0].category, RuntimeWarning))

    def test_empty(self):
        A = np.zeros((0, 3))
        for f in self.funcs:
            for axis in [0, None]:
                with warnings.catch_warnings(record=True) as w:
                    warnings.simplefilter('always')
                    assert_(np.isnan(f(A, axis=axis)).all())
                    assert_(len(w) > 0)
                    assert_(issubclass(w[0].category, RuntimeWarning))
            for axis in [1]:
                with warnings.catch_warnings(record=True) as w:
                    warnings.simplefilter('always')
                    assert_equal(f(A, axis=axis), np.zeros([]))

    def test_mean_values(self):
        for mat in [self.rmat, self.cmat, self.omat]:
            for axis in [0, 1]:
                tgt = mat.sum(axis=axis)
                res = _mean(mat, axis=axis) * mat.shape[axis]
                assert_almost_equal(res, tgt)
            for axis in [None]:
                tgt = mat.sum(axis=axis)
                res = _mean(mat, axis=axis) * np.prod(mat.shape)
                assert_almost_equal(res, tgt)

    def test_var_values(self):
        for mat in [self.rmat, self.cmat, self.omat]:
            for axis in [0, 1, None]:
                msqr = _mean(mat * mat.conj(), axis=axis)
                mean = _mean(mat, axis=axis)
                tgt = msqr - mean * mean.conjugate()
                res = _var(mat, axis=axis)
                assert_almost_equal(res, tgt)

    def test_std_values(self):
        for mat in [self.rmat, self.cmat, self.omat]:
            for axis in [0, 1, None]:
                tgt = np.sqrt(_var(mat, axis=axis))
                res = _std(mat, axis=axis)
                assert_almost_equal(res, tgt)


    def test_subclass(self):
        class TestArray(np.ndarray):
            def __new__(cls, data, info):
                result = np.array(data)
                result = result.view(cls)
                result.info = info
                return result
            def __array_finalize__(self, obj):
                self.info = getattr(obj, "info", '')

        dat = TestArray([[1, 2, 3, 4], [5, 6, 7, 8]], 'jubba')
        res = dat.mean(1)
        assert_(res.info == dat.info)
        res = dat.std(1)
        assert_(res.info == dat.info)
        res = dat.var(1)
        assert_(res.info == dat.info)

class TestDot(TestCase):
    def test_dot_2args(self):
        from numpy.core.multiarray import dot

        a = np.array([[1, 2], [3, 4]], dtype=float)
        b = np.array([[1, 0], [1, 1]], dtype=float)
        c = np.array([[3, 2], [7, 4]], dtype=float)

        d = dot(a, b)
        assert_allclose(c, d)

    def test_dot_3args(self):
        from numpy.core.multiarray import dot

        np.random.seed(22)
        f = np.random.random_sample((1024, 16))
        v = np.random.random_sample((16, 32))

        r = np.empty((1024, 32))
        for i in range(12):
            dot(f, v, r)
        assert_equal(sys.getrefcount(r), 2)
        r2 = dot(f, v, out=None)
        assert_array_equal(r2, r)
        assert_(r is dot(f, v, out=r))

        v = v[:, 0].copy() # v.shape == (16,)
        r = r[:, 0].copy() # r.shape == (1024,)
        r2 = dot(f, v)
        assert_(r is dot(f, v, r))
        assert_array_equal(r2, r)

    def test_dot_3args_errors(self):
        from numpy.core.multiarray import dot

        np.random.seed(22)
        f = np.random.random_sample((1024, 16))
        v = np.random.random_sample((16, 32))

        r = np.empty((1024, 31))
        assert_raises(ValueError, dot, f, v, r)

        r = np.empty((1024,))
        assert_raises(ValueError, dot, f, v, r)

        r = np.empty((32,))
        assert_raises(ValueError, dot, f, v, r)

        r = np.empty((32, 1024))
        assert_raises(ValueError, dot, f, v, r)
        assert_raises(ValueError, dot, f, v, r.T)

        r = np.empty((1024, 64))
        assert_raises(ValueError, dot, f, v, r[:, ::2])
        assert_raises(ValueError, dot, f, v, r[:, :32])

        r = np.empty((1024, 32), dtype=np.float32)
        assert_raises(ValueError, dot, f, v, r)

        r = np.empty((1024, 32), dtype=int)
        assert_raises(ValueError, dot, f, v, r)

    def test_dot_scalar_and_matrix_of_objects(self):
        # Ticket #2469
        arr = np.matrix([1, 2], dtype=object)
        desired = np.matrix([[3, 6]], dtype=object)
        assert_equal(np.dot(arr, 3), desired)
        assert_equal(np.dot(3, arr), desired)


class TestInner(TestCase):

    def test_inner_scalar_and_matrix_of_objects(self):
        # Ticket #4482
        arr = np.matrix([1, 2], dtype=object)
        desired = np.matrix([[3, 6]], dtype=object)
        assert_equal(np.inner(arr, 3), desired)
        assert_equal(np.inner(3, arr), desired)


class TestSummarization(TestCase):
    def test_1d(self):
        A = np.arange(1001)
        strA = '[   0    1    2 ...,  998  999 1000]'
        assert_(str(A) == strA)

        reprA = 'array([   0,    1,    2, ...,  998,  999, 1000])'
        assert_(repr(A) == reprA)

    def test_2d(self):
        A = np.arange(1002).reshape(2, 501)
        strA = '[[   0    1    2 ...,  498  499  500]\n' \
               ' [ 501  502  503 ...,  999 1000 1001]]'
        assert_(str(A) == strA)

        reprA = 'array([[   0,    1,    2, ...,  498,  499,  500],\n' \
                '       [ 501,  502,  503, ...,  999, 1000, 1001]])'
        assert_(repr(A) == reprA)


class TestChoose(TestCase):
    def setUp(self):
        self.x = 2*ones((3,), dtype=int)
        self.y = 3*ones((3,), dtype=int)
        self.x2 = 2*ones((2, 3), dtype=int)
        self.y2 = 3*ones((2, 3), dtype=int)
        self.ind = [0, 0, 1]

    def test_basic(self):
        A = np.choose(self.ind, (self.x, self.y))
        assert_equal(A, [2, 2, 3])

    def test_broadcast1(self):
        A = np.choose(self.ind, (self.x2, self.y2))
        assert_equal(A, [[2, 2, 3], [2, 2, 3]])

    def test_broadcast2(self):
        A = np.choose(self.ind, (self.x, self.y2))
        assert_equal(A, [[2, 2, 3], [2, 2, 3]])


# TODO: test for multidimensional
NEIGH_MODE = {'zero': 0, 'one': 1, 'constant': 2, 'circular': 3, 'mirror': 4}
class TestNeighborhoodIter(TestCase):
    # Simple, 2d tests
    def _test_simple2d(self, dt):
        # Test zero and one padding for simple data type
        x = np.array([[0, 1], [2, 3]], dtype=dt)
        r = [np.array([[0, 0, 0], [0, 0, 1]], dtype=dt),
             np.array([[0, 0, 0], [0, 1, 0]], dtype=dt),
             np.array([[0, 0, 1], [0, 2, 3]], dtype=dt),
             np.array([[0, 1, 0], [2, 3, 0]], dtype=dt)]
        l = test_neighborhood_iterator(x, [-1, 0, -1, 1], x[0],
                NEIGH_MODE['zero'])
        assert_array_equal(l, r)

        r = [np.array([[1, 1, 1], [1, 0, 1]], dtype=dt),
             np.array([[1, 1, 1], [0, 1, 1]], dtype=dt),
             np.array([[1, 0, 1], [1, 2, 3]], dtype=dt),
             np.array([[0, 1, 1], [2, 3, 1]], dtype=dt)]
        l = test_neighborhood_iterator(x, [-1, 0, -1, 1], x[0],
                NEIGH_MODE['one'])
        assert_array_equal(l, r)

        r = [np.array([[4, 4, 4], [4, 0, 1]], dtype=dt),
             np.array([[4, 4, 4], [0, 1, 4]], dtype=dt),
             np.array([[4, 0, 1], [4, 2, 3]], dtype=dt),
             np.array([[0, 1, 4], [2, 3, 4]], dtype=dt)]
        l = test_neighborhood_iterator(x, [-1, 0, -1, 1], 4,
                NEIGH_MODE['constant'])
        assert_array_equal(l, r)

    def test_simple2d(self):
        self._test_simple2d(np.float)

    def test_simple2d_object(self):
        self._test_simple2d(Decimal)

    def _test_mirror2d(self, dt):
        x = np.array([[0, 1], [2, 3]], dtype=dt)
        r = [np.array([[0, 0, 1], [0, 0, 1]], dtype=dt),
             np.array([[0, 1, 1], [0, 1, 1]], dtype=dt),
             np.array([[0, 0, 1], [2, 2, 3]], dtype=dt),
             np.array([[0, 1, 1], [2, 3, 3]], dtype=dt)]
        l = test_neighborhood_iterator(x, [-1, 0, -1, 1], x[0],
                NEIGH_MODE['mirror'])
        assert_array_equal(l, r)

    def test_mirror2d(self):
        self._test_mirror2d(np.float)

    def test_mirror2d_object(self):
        self._test_mirror2d(Decimal)

    # Simple, 1d tests
    def _test_simple(self, dt):
        # Test padding with constant values
        x = np.linspace(1, 5, 5).astype(dt)
        r = [[0, 1, 2], [1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 0]]
        l = test_neighborhood_iterator(x, [-1, 1], x[0], NEIGH_MODE['zero'])
        assert_array_equal(l, r)

        r = [[1, 1, 2], [1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 1]]
        l = test_neighborhood_iterator(x, [-1, 1], x[0], NEIGH_MODE['one'])
        assert_array_equal(l, r)

        r = [[x[4], 1, 2], [1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, x[4]]]
        l = test_neighborhood_iterator(x, [-1, 1], x[4], NEIGH_MODE['constant'])
        assert_array_equal(l, r)

    def test_simple_float(self):
        self._test_simple(np.float)

    def test_simple_object(self):
        self._test_simple(Decimal)

    # Test mirror modes
    def _test_mirror(self, dt):
        x = np.linspace(1, 5, 5).astype(dt)
        r = np.array([[2, 1, 1, 2, 3], [1, 1, 2, 3, 4], [1, 2, 3, 4, 5],
                [2, 3, 4, 5, 5], [3, 4, 5, 5, 4]], dtype=dt)
        l = test_neighborhood_iterator(x, [-2, 2], x[1], NEIGH_MODE['mirror'])
        self.assertTrue([i.dtype == dt for i in l])
        assert_array_equal(l, r)

    def test_mirror(self):
        self._test_mirror(np.float)

    def test_mirror_object(self):
        self._test_mirror(Decimal)

    # Circular mode
    def _test_circular(self, dt):
        x = np.linspace(1, 5, 5).astype(dt)
        r = np.array([[4, 5, 1, 2, 3], [5, 1, 2, 3, 4], [1, 2, 3, 4, 5],
                [2, 3, 4, 5, 1], [3, 4, 5, 1, 2]], dtype=dt)
        l = test_neighborhood_iterator(x, [-2, 2], x[0], NEIGH_MODE['circular'])
        assert_array_equal(l, r)

    def test_circular(self):
        self._test_circular(np.float)

    def test_circular_object(self):
        self._test_circular(Decimal)

# Test stacking neighborhood iterators
class TestStackedNeighborhoodIter(TestCase):
    # Simple, 1d test: stacking 2 constant-padded neigh iterators
    def test_simple_const(self):
        dt = np.float64
        # Test zero and one padding for simple data type
        x = np.array([1, 2, 3], dtype=dt)
        r = [np.array([0], dtype=dt),
             np.array([0], dtype=dt),
             np.array([1], dtype=dt),
             np.array([2], dtype=dt),
             np.array([3], dtype=dt),
             np.array([0], dtype=dt),
             np.array([0], dtype=dt)]
        l = test_neighborhood_iterator_oob(x, [-2, 4], NEIGH_MODE['zero'],
                [0, 0], NEIGH_MODE['zero'])
        assert_array_equal(l, r)

        r = [np.array([1, 0, 1], dtype=dt),
             np.array([0, 1, 2], dtype=dt),
             np.array([1, 2, 3], dtype=dt),
             np.array([2, 3, 0], dtype=dt),
             np.array([3, 0, 1], dtype=dt)]
        l = test_neighborhood_iterator_oob(x, [-1, 3], NEIGH_MODE['zero'],
                [-1, 1], NEIGH_MODE['one'])
        assert_array_equal(l, r)

    # 2nd simple, 1d test: stacking 2 neigh iterators, mixing const padding and
    # mirror padding
    def test_simple_mirror(self):
        dt = np.float64
        # Stacking zero on top of mirror
        x = np.array([1, 2, 3], dtype=dt)
        r = [np.array([0, 1, 1], dtype=dt),
             np.array([1, 1, 2], dtype=dt),
             np.array([1, 2, 3], dtype=dt),
             np.array([2, 3, 3], dtype=dt),
             np.array([3, 3, 0], dtype=dt)]
        l = test_neighborhood_iterator_oob(x, [-1, 3], NEIGH_MODE['mirror'],
                [-1, 1], NEIGH_MODE['zero'])
        assert_array_equal(l, r)

        # Stacking mirror on top of zero
        x = np.array([1, 2, 3], dtype=dt)
        r = [np.array([1, 0, 0], dtype=dt),
             np.array([0, 0, 1], dtype=dt),
             np.array([0, 1, 2], dtype=dt),
             np.array([1, 2, 3], dtype=dt),
             np.array([2, 3, 0], dtype=dt)]
        l = test_neighborhood_iterator_oob(x, [-1, 3], NEIGH_MODE['zero'],
                [-2, 0], NEIGH_MODE['mirror'])
        assert_array_equal(l, r)

        # Stacking mirror on top of zero: 2nd
        x = np.array([1, 2, 3], dtype=dt)
        r = [np.array([0, 1, 2], dtype=dt),
             np.array([1, 2, 3], dtype=dt),
             np.array([2, 3, 0], dtype=dt),
             np.array([3, 0, 0], dtype=dt),
             np.array([0, 0, 3], dtype=dt)]
        l = test_neighborhood_iterator_oob(x, [-1, 3], NEIGH_MODE['zero'],
                [0, 2], NEIGH_MODE['mirror'])
        assert_array_equal(l, r)

        # Stacking mirror on top of zero: 3rd
        x = np.array([1, 2, 3], dtype=dt)
        r = [np.array([1, 0, 0, 1, 2], dtype=dt),
             np.array([0, 0, 1, 2, 3], dtype=dt),
             np.array([0, 1, 2, 3, 0], dtype=dt),
             np.array([1, 2, 3, 0, 0], dtype=dt),
             np.array([2, 3, 0, 0, 3], dtype=dt)]
        l = test_neighborhood_iterator_oob(x, [-1, 3], NEIGH_MODE['zero'],
                [-2, 2], NEIGH_MODE['mirror'])
        assert_array_equal(l, r)

    # 3rd simple, 1d test: stacking 2 neigh iterators, mixing const padding and
    # circular padding
    def test_simple_circular(self):
        dt = np.float64
        # Stacking zero on top of mirror
        x = np.array([1, 2, 3], dtype=dt)
        r = [np.array([0, 3, 1], dtype=dt),
             np.array([3, 1, 2], dtype=dt),
             np.array([1, 2, 3], dtype=dt),
             np.array([2, 3, 1], dtype=dt),
             np.array([3, 1, 0], dtype=dt)]
        l = test_neighborhood_iterator_oob(x, [-1, 3], NEIGH_MODE['circular'],
                [-1, 1], NEIGH_MODE['zero'])
        assert_array_equal(l, r)

        # Stacking mirror on top of zero
        x = np.array([1, 2, 3], dtype=dt)
        r = [np.array([3, 0, 0], dtype=dt),
             np.array([0, 0, 1], dtype=dt),
             np.array([0, 1, 2], dtype=dt),
             np.array([1, 2, 3], dtype=dt),
             np.array([2, 3, 0], dtype=dt)]
        l = test_neighborhood_iterator_oob(x, [-1, 3], NEIGH_MODE['zero'],
                [-2, 0], NEIGH_MODE['circular'])
        assert_array_equal(l, r)

        # Stacking mirror on top of zero: 2nd
        x = np.array([1, 2, 3], dtype=dt)
        r = [np.array([0, 1, 2], dtype=dt),
             np.array([1, 2, 3], dtype=dt),
             np.array([2, 3, 0], dtype=dt),
             np.array([3, 0, 0], dtype=dt),
             np.array([0, 0, 1], dtype=dt)]
        l = test_neighborhood_iterator_oob(x, [-1, 3], NEIGH_MODE['zero'],
                [0, 2], NEIGH_MODE['circular'])
        assert_array_equal(l, r)

        # Stacking mirror on top of zero: 3rd
        x = np.array([1, 2, 3], dtype=dt)
        r = [np.array([3, 0, 0, 1, 2], dtype=dt),
             np.array([0, 0, 1, 2, 3], dtype=dt),
             np.array([0, 1, 2, 3, 0], dtype=dt),
             np.array([1, 2, 3, 0, 0], dtype=dt),
             np.array([2, 3, 0, 0, 1], dtype=dt)]
        l = test_neighborhood_iterator_oob(x, [-1, 3], NEIGH_MODE['zero'],
                [-2, 2], NEIGH_MODE['circular'])
        assert_array_equal(l, r)

    # 4th simple, 1d test: stacking 2 neigh iterators, but with lower iterator
    # being strictly within the array
    def test_simple_strict_within(self):
        dt = np.float64
        # Stacking zero on top of zero, first neighborhood strictly inside the
        # array
        x = np.array([1, 2, 3], dtype=dt)
        r = [np.array([1, 2, 3, 0], dtype=dt)]
        l = test_neighborhood_iterator_oob(x, [1, 1], NEIGH_MODE['zero'],
                [-1, 2], NEIGH_MODE['zero'])
        assert_array_equal(l, r)

        # Stacking mirror on top of zero, first neighborhood strictly inside the
        # array
        x = np.array([1, 2, 3], dtype=dt)
        r = [np.array([1, 2, 3, 3], dtype=dt)]
        l = test_neighborhood_iterator_oob(x, [1, 1], NEIGH_MODE['zero'],
                [-1, 2], NEIGH_MODE['mirror'])
        assert_array_equal(l, r)

        # Stacking mirror on top of zero, first neighborhood strictly inside the
        # array
        x = np.array([1, 2, 3], dtype=dt)
        r = [np.array([1, 2, 3, 1], dtype=dt)]
        l = test_neighborhood_iterator_oob(x, [1, 1], NEIGH_MODE['zero'],
                [-1, 2], NEIGH_MODE['circular'])
        assert_array_equal(l, r)

class TestWarnings(object):

    def test_complex_warning(self):
        x = np.array([1, 2])
        y = np.array([1-2j, 1+2j])

        with warnings.catch_warnings():
            warnings.simplefilter("error", np.ComplexWarning)
            assert_raises(np.ComplexWarning, x.__setitem__, slice(None), y)
            assert_equal(x, [1, 2])

class TestMinScalarType(object):

    def test_usigned_shortshort(self):
        dt = np.min_scalar_type(2**8-1)
        wanted = np.dtype('uint8')
        assert_equal(wanted, dt)

    def test_usigned_short(self):
        dt = np.min_scalar_type(2**16-1)
        wanted = np.dtype('uint16')
        assert_equal(wanted, dt)

    def test_usigned_int(self):
        dt = np.min_scalar_type(2**32-1)
        wanted = np.dtype('uint32')
        assert_equal(wanted, dt)

    def test_usigned_longlong(self):
        dt = np.min_scalar_type(2**63-1)
        wanted = np.dtype('uint64')
        assert_equal(wanted, dt)

    def test_object(self):
        dt = np.min_scalar_type(2**64)
        wanted = np.dtype('O')
        assert_equal(wanted, dt)


if sys.version_info[:2] == (2, 6):
    from numpy.core.multiarray import memorysimpleview as memoryview

from numpy.core._internal import _dtype_from_pep3118

class TestPEP3118Dtype(object):
    def _check(self, spec, wanted):
        dt = np.dtype(wanted)
        if isinstance(wanted, list) and isinstance(wanted[-1], tuple):
            if wanted[-1][0] == '':
                names = list(dt.names)
                names[-1] = ''
                dt.names = tuple(names)
        assert_equal(_dtype_from_pep3118(spec), dt,
                     err_msg="spec %r != dtype %r" % (spec, wanted))

    def test_native_padding(self):
        align = np.dtype('i').alignment
        for j in range(8):
            if j == 0:
                s = 'bi'
            else:
                s = 'b%dxi' % j
            self._check('@'+s, {'f0': ('i1', 0),
                                'f1': ('i', align*(1 + j//align))})
            self._check('='+s, {'f0': ('i1', 0),
                                'f1': ('i', 1+j)})

    def test_native_padding_2(self):
        # Native padding should work also for structs and sub-arrays
        self._check('x3T{xi}', {'f0': (({'f0': ('i', 4)}, (3,)), 4)})
        self._check('^x3T{xi}', {'f0': (({'f0': ('i', 1)}, (3,)), 1)})

    def test_trailing_padding(self):
        # Trailing padding should be included, *and*, the item size
        # should match the alignment if in aligned mode
        align = np.dtype('i').alignment
        def VV(n):
            return 'V%d' % (align*(1 + (n-1)//align))

        self._check('ix', [('f0', 'i'), ('', VV(1))])
        self._check('ixx', [('f0', 'i'), ('', VV(2))])
        self._check('ixxx', [('f0', 'i'), ('', VV(3))])
        self._check('ixxxx', [('f0', 'i'), ('', VV(4))])
        self._check('i7x', [('f0', 'i'), ('', VV(7))])

        self._check('^ix', [('f0', 'i'), ('', 'V1')])
        self._check('^ixx', [('f0', 'i'), ('', 'V2')])
        self._check('^ixxx', [('f0', 'i'), ('', 'V3')])
        self._check('^ixxxx', [('f0', 'i'), ('', 'V4')])
        self._check('^i7x', [('f0', 'i'), ('', 'V7')])

    def test_native_padding_3(self):
        dt = np.dtype(
                [('a', 'b'), ('b', 'i'),
                    ('sub', np.dtype('b,i')), ('c', 'i')],
                align=True)
        self._check("T{b:a:xxxi:b:T{b:f0:=i:f1:}:sub:xxxi:c:}", dt)

        dt = np.dtype(
                [('a', 'b'), ('b', 'i'), ('c', 'b'), ('d', 'b'),
                    ('e', 'b'), ('sub', np.dtype('b,i', align=True))])
        self._check("T{b:a:=i:b:b:c:b:d:b:e:T{b:f0:xxxi:f1:}:sub:}", dt)

    def test_padding_with_array_inside_struct(self):
        dt = np.dtype(
                [('a', 'b'), ('b', 'i'), ('c', 'b', (3,)),
                    ('d', 'i')],
                align=True)
        self._check("T{b:a:xxxi:b:3b:c:xi:d:}", dt)

    def test_byteorder_inside_struct(self):
        # The byte order after @T{=i} should be '=', not '@'.
        # Check this by noting the absence of native alignment.
        self._check('@T{^i}xi', {'f0': ({'f0': ('i', 0)}, 0),
                                 'f1': ('i', 5)})

    def test_intra_padding(self):
        # Natively aligned sub-arrays may require some internal padding
        align = np.dtype('i').alignment
        def VV(n):
            return 'V%d' % (align*(1 + (n-1)//align))

        self._check('(3)T{ix}', ({'f0': ('i', 0), '': (VV(1), 4)}, (3,)))

class TestNewBufferProtocol(object):
    def _check_roundtrip(self, obj):
        obj = np.asarray(obj)
        x = memoryview(obj)
        y = np.asarray(x)
        y2 = np.array(x)
        assert_(not y.flags.owndata)
        assert_(y2.flags.owndata)

        assert_equal(y.dtype, obj.dtype)
        assert_equal(y.shape, obj.shape)
        assert_array_equal(obj, y)

        assert_equal(y2.dtype, obj.dtype)
        assert_equal(y2.shape, obj.shape)
        assert_array_equal(obj, y2)

    def test_roundtrip(self):
        x = np.array([1, 2, 3, 4, 5], dtype='i4')
        self._check_roundtrip(x)

        x = np.array([[1, 2], [3, 4]], dtype=np.float64)
        self._check_roundtrip(x)

        x = np.zeros((3, 3, 3), dtype=np.float32)[:, 0,:]
        self._check_roundtrip(x)

        dt = [('a', 'b'),
              ('b', 'h'),
              ('c', 'i'),
              ('d', 'l'),
              ('dx', 'q'),
              ('e', 'B'),
              ('f', 'H'),
              ('g', 'I'),
              ('h', 'L'),
              ('hx', 'Q'),
              ('i', np.single),
              ('j', np.double),
              ('k', np.longdouble),
              ('ix', np.csingle),
              ('jx', np.cdouble),
              ('kx', np.clongdouble),
              ('l', 'S4'),
              ('m', 'U4'),
              ('n', 'V3'),
              ('o', '?'),
              ('p', np.half),
             ]
        x = np.array(
                [(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    asbytes('aaaa'), 'bbbb', asbytes('xxx'), True, 1.0)],
                dtype=dt)
        self._check_roundtrip(x)

        x = np.array(([[1, 2], [3, 4]],), dtype=[('a', (int, (2, 2)))])
        self._check_roundtrip(x)

        x = np.array([1, 2, 3], dtype='>i2')
        self._check_roundtrip(x)

        x = np.array([1, 2, 3], dtype='<i2')
        self._check_roundtrip(x)

        x = np.array([1, 2, 3], dtype='>i4')
        self._check_roundtrip(x)

        x = np.array([1, 2, 3], dtype='<i4')
        self._check_roundtrip(x)

        # check long long can be represented as non-native
        x = np.array([1, 2, 3], dtype='>q')
        self._check_roundtrip(x)

        # Native-only data types can be passed through the buffer interface
        # only in native byte order
        if sys.byteorder == 'little':
            x = np.array([1, 2, 3], dtype='>g')
            assert_raises(ValueError, self._check_roundtrip, x)
            x = np.array([1, 2, 3], dtype='<g')
            self._check_roundtrip(x)
        else:
            x = np.array([1, 2, 3], dtype='>g')
            self._check_roundtrip(x)
            x = np.array([1, 2, 3], dtype='<g')
            assert_raises(ValueError, self._check_roundtrip, x)

    def test_roundtrip_half(self):
        half_list = [
            1.0,
            -2.0,
            6.5504 * 10**4, #  (max half precision)
            2**-14, # ~= 6.10352 * 10**-5 (minimum positive normal)
            2**-24, # ~= 5.96046 * 10**-8 (minimum strictly positive subnormal)
            0.0,
            -0.0,
            float('+inf'),
            float('-inf'),
            0.333251953125, # ~= 1/3
        ]

        x = np.array(half_list, dtype='>e')
        self._check_roundtrip(x)
        x = np.array(half_list, dtype='<e')
        self._check_roundtrip(x)

    def test_roundtrip_single_types(self):
        for typ in np.typeDict.values():
            dtype = np.dtype(typ)

            if dtype.char in 'Mm':
                # datetimes cannot be used in buffers
                continue
            if dtype.char == 'V':
                # skip void
                continue

            x = np.zeros(4, dtype=dtype)
            self._check_roundtrip(x)

            if dtype.char not in 'qQgG':
                dt = dtype.newbyteorder('<')
                x = np.zeros(4, dtype=dt)
                self._check_roundtrip(x)

                dt = dtype.newbyteorder('>')
                x = np.zeros(4, dtype=dt)
                self._check_roundtrip(x)

    def test_roundtrip_scalar(self):
        # Issue #4015.
        self._check_roundtrip(0)

    def test_export_simple_1d(self):
        x = np.array([1, 2, 3, 4, 5], dtype='i')
        y = memoryview(x)
        assert_equal(y.format, 'i')
        assert_equal(y.shape, (5,))
        assert_equal(y.ndim, 1)
        assert_equal(y.strides, (4,))
        assert_equal(y.suboffsets, EMPTY)
        assert_equal(y.itemsize, 4)

    def test_export_simple_nd(self):
        x = np.array([[1, 2], [3, 4]], dtype=np.float64)
        y = memoryview(x)
        assert_equal(y.format, 'd')
        assert_equal(y.shape, (2, 2))
        assert_equal(y.ndim, 2)
        assert_equal(y.strides, (16, 8))
        assert_equal(y.suboffsets, EMPTY)
        assert_equal(y.itemsize, 8)

    def test_export_discontiguous(self):
        x = np.zeros((3, 3, 3), dtype=np.float32)[:, 0,:]
        y = memoryview(x)
        assert_equal(y.format, 'f')
        assert_equal(y.shape, (3, 3))
        assert_equal(y.ndim, 2)
        assert_equal(y.strides, (36, 4))
        assert_equal(y.suboffsets, EMPTY)
        assert_equal(y.itemsize, 4)

    def test_export_record(self):
        dt = [('a', 'b'),
              ('b', 'h'),
              ('c', 'i'),
              ('d', 'l'),
              ('dx', 'q'),
              ('e', 'B'),
              ('f', 'H'),
              ('g', 'I'),
              ('h', 'L'),
              ('hx', 'Q'),
              ('i', np.single),
              ('j', np.double),
              ('k', np.longdouble),
              ('ix', np.csingle),
              ('jx', np.cdouble),
              ('kx', np.clongdouble),
              ('l', 'S4'),
              ('m', 'U4'),
              ('n', 'V3'),
              ('o', '?'),
              ('p', np.half),
             ]
        x = np.array(
                [(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    asbytes('aaaa'), 'bbbb', asbytes('   '), True, 1.0)],
                dtype=dt)
        y = memoryview(x)
        assert_equal(y.shape, (1,))
        assert_equal(y.ndim, 1)
        assert_equal(y.suboffsets, EMPTY)

        sz = sum([dtype(b).itemsize for a, b in dt])
        if dtype('l').itemsize == 4:
            assert_equal(y.format, 'T{b:a:=h:b:i:c:l:d:q:dx:B:e:@H:f:=I:g:L:h:Q:hx:f:i:d:j:^g:k:=Zf:ix:Zd:jx:^Zg:kx:4s:l:=4w:m:3x:n:?:o:@e:p:}')
        else:
            assert_equal(y.format, 'T{b:a:=h:b:i:c:q:d:q:dx:B:e:@H:f:=I:g:Q:h:Q:hx:f:i:d:j:^g:k:=Zf:ix:Zd:jx:^Zg:kx:4s:l:=4w:m:3x:n:?:o:@e:p:}')
        # Cannot test if NPY_RELAXED_STRIDES_CHECKING changes the strides
        if not (np.ones(1).strides[0] == np.iinfo(np.intp).max):
            assert_equal(y.strides, (sz,))
        assert_equal(y.itemsize, sz)

    def test_export_subarray(self):
        x = np.array(([[1, 2], [3, 4]],), dtype=[('a', ('i', (2, 2)))])
        y = memoryview(x)
        assert_equal(y.format, 'T{(2,2)i:a:}')
        assert_equal(y.shape, EMPTY)
        assert_equal(y.ndim, 0)
        assert_equal(y.strides, EMPTY)
        assert_equal(y.suboffsets, EMPTY)
        assert_equal(y.itemsize, 16)

    def test_export_endian(self):
        x = np.array([1, 2, 3], dtype='>i')
        y = memoryview(x)
        if sys.byteorder == 'little':
            assert_equal(y.format, '>i')
        else:
            assert_equal(y.format, 'i')

        x = np.array([1, 2, 3], dtype='<i')
        y = memoryview(x)
        if sys.byteorder == 'little':
            assert_equal(y.format, 'i')
        else:
            assert_equal(y.format, '<i')

    def test_export_flags(self):
        # Check SIMPLE flag, see also gh-3613 (exception should be BufferError)
        assert_raises(ValueError, get_buffer_info, np.arange(5)[::2], ('SIMPLE',))

    def test_padding(self):
        for j in range(8):
            x = np.array([(1,), (2,)], dtype={'f0': (int, j)})
            self._check_roundtrip(x)

    def test_reference_leak(self):
        count_1 = sys.getrefcount(np.core._internal)
        a = np.zeros(4)
        b = memoryview(a)
        c = np.asarray(b)
        count_2 = sys.getrefcount(np.core._internal)
        assert_equal(count_1, count_2)

    def test_padded_struct_array(self):
        dt1 = np.dtype(
                [('a', 'b'), ('b', 'i'), ('sub', np.dtype('b,i')), ('c', 'i')],
                align=True)
        x1 = np.arange(dt1.itemsize, dtype=np.int8).view(dt1)
        self._check_roundtrip(x1)

        dt2 = np.dtype(
                [('a', 'b'), ('b', 'i'), ('c', 'b', (3,)), ('d', 'i')],
                align=True)
        x2 = np.arange(dt2.itemsize, dtype=np.int8).view(dt2)
        self._check_roundtrip(x2)

        dt3 = np.dtype(
                [('a', 'b'), ('b', 'i'), ('c', 'b'), ('d', 'b'),
                    ('e', 'b'), ('sub', np.dtype('b,i', align=True))])
        x3 = np.arange(dt3.itemsize, dtype=np.int8).view(dt3)
        self._check_roundtrip(x3)


class TestArrayAttributeDeletion(object):

    def test_multiarray_writable_attributes_deletion(self):
        """ticket #2046, should not seqfault, raise AttributeError"""
        a = np.ones(2)
        attr =  ['shape', 'strides', 'data', 'dtype', 'real', 'imag', 'flat']
        for s in attr:
            assert_raises(AttributeError, delattr, a, s)


    def test_multiarray_not_writable_attributes_deletion(self):
        a = np.ones(2)
        attr = ["ndim", "flags", "itemsize", "size", "nbytes", "base",
                "ctypes", "T", "__array_interface__", "__array_struct__",
                "__array_priority__", "__array_finalize__"]
        for s in attr:
            assert_raises(AttributeError, delattr, a, s)


    def test_multiarray_flags_writable_attribute_deletion(self):
        a = np.ones(2).flags
        attr = ['updateifcopy', 'aligned', 'writeable']
        for s in attr:
            assert_raises(AttributeError, delattr, a, s)


    def test_multiarray_flags_not_writable_attribute_deletion(self):
        a = np.ones(2).flags
        attr = ["contiguous", "c_contiguous", "f_contiguous", "fortran",
                "owndata", "fnc", "forc", "behaved", "carray", "farray",
                "num"]
        for s in attr:
            assert_raises(AttributeError, delattr, a, s)


def test_array_interface():
    # Test scalar coercion within the array interface
    class Foo(object):
        def __init__(self, value):
            self.value = value
            self.iface = {'typestr' : '=f8'}
        def __float__(self):
            return float(self.value)
        @property
        def __array_interface__(self):
            return self.iface
    f = Foo(0.5)
    assert_equal(np.array(f), 0.5)
    assert_equal(np.array([f]), [0.5])
    assert_equal(np.array([f, f]), [0.5, 0.5])
    assert_equal(np.array(f).dtype, np.dtype('=f8'))
    # Test various shape definitions
    f.iface['shape'] = ()
    assert_equal(np.array(f), 0.5)
    f.iface['shape'] = None
    assert_raises(TypeError, np.array, f)
    f.iface['shape'] = (1, 1)
    assert_equal(np.array(f), [[0.5]])
    f.iface['shape'] = (2,)
    assert_raises(ValueError, np.array, f)

    # test scalar with no shape
    class ArrayLike(object):
        array = np.array(1)
        __array_interface__ = array.__array_interface__
    assert_equal(np.array(ArrayLike()), 1)


def test_flat_element_deletion():
    it = np.ones(3).flat
    try:
        del it[1]
        del it[1:2]
    except TypeError:
        pass
    except:
        raise AssertionError

def test_scalar_element_deletion():
    a = np.zeros(2, dtype=[('x', 'int'), ('y', 'int')])
    assert_raises(ValueError, a[0].__delitem__, 'x')

class TestMemEventHook(TestCase):
    def test_mem_seteventhook(self):
        # The actual tests are within the C code in
        # multiarray/multiarray_tests.c.src
        test_pydatamem_seteventhook_start()
        # force an allocation and free of a numpy array
        # needs to be larger then limit of small memory cacher in ctors.c
        a = np.zeros(1000)
        del a
        test_pydatamem_seteventhook_end()

class TestMapIter(TestCase):
    def test_mapiter(self):
        # The actual tests are within the C code in
        # multiarray/multiarray_tests.c.src

        a = arange(12).reshape((3, 4)).astype(float)
        index = ([1, 1, 2, 0],
                 [0, 0, 2, 3])
        vals = [50, 50, 30, 16]

        test_inplace_increment(a, index, vals)
        assert_equal(a, [[   0.,   1.,    2.,   19.,],
                         [ 104.,    5.,    6.,    7.,],
                         [   8.,    9.,   40.,   11.,]])

        b = arange(6).astype(float)
        index = (array([1, 2, 0]),)
        vals = [50, 4, 100.1]
        test_inplace_increment(b, index, vals)
        assert_equal(b, [ 100.1,  51.,   6.,   3.,   4.,   5. ])


class PriorityNdarray():
    __array_priority__ = 1000

    def __init__(self, array):
        self.array = array

    def __lt__(self, array):
        if isinstance(array, PriorityNdarray):
            array = array.array
        return PriorityNdarray(self.array < array)

    def __gt__(self, array):
        if isinstance(array, PriorityNdarray):
            array = array.array
        return PriorityNdarray(self.array > array)

    def __le__(self, array):
        if isinstance(array, PriorityNdarray):
            array = array.array
        return PriorityNdarray(self.array <= array)

    def __ge__(self, array):
        if isinstance(array, PriorityNdarray):
            array = array.array
        return PriorityNdarray(self.array >= array)

    def __eq__(self, array):
        if isinstance(array, PriorityNdarray):
            array = array.array
        return PriorityNdarray(self.array == array)

    def __ne__(self, array):
        if isinstance(array, PriorityNdarray):
            array = array.array
        return PriorityNdarray(self.array != array)


class TestArrayPriority(TestCase):
    def test_lt(self):
        l = np.asarray([0., -1., 1.], dtype=dtype)
        r = np.asarray([0., 1., -1.], dtype=dtype)
        lp = PriorityNdarray(l)
        rp = PriorityNdarray(r)
        res1 = l < r
        res2 = l < rp
        res3 = lp < r
        res4 = lp < rp

        assert_array_equal(res1, res2.array)
        assert_array_equal(res1, res3.array)
        assert_array_equal(res1, res4.array)
        assert_(isinstance(res1, np.ndarray))
        assert_(isinstance(res2, PriorityNdarray))
        assert_(isinstance(res3, PriorityNdarray))
        assert_(isinstance(res4, PriorityNdarray))

    def test_gt(self):
        l = np.asarray([0., -1., 1.], dtype=dtype)
        r = np.asarray([0., 1., -1.], dtype=dtype)
        lp = PriorityNdarray(l)
        rp = PriorityNdarray(r)
        res1 = l > r
        res2 = l > rp
        res3 = lp > r
        res4 = lp > rp

        assert_array_equal(res1, res2.array)
        assert_array_equal(res1, res3.array)
        assert_array_equal(res1, res4.array)
        assert_(isinstance(res1, np.ndarray))
        assert_(isinstance(res2, PriorityNdarray))
        assert_(isinstance(res3, PriorityNdarray))
        assert_(isinstance(res4, PriorityNdarray))

    def test_le(self):
        l = np.asarray([0., -1., 1.], dtype=dtype)
        r = np.asarray([0., 1., -1.], dtype=dtype)
        lp = PriorityNdarray(l)
        rp = PriorityNdarray(r)
        res1 = l <= r
        res2 = l <= rp
        res3 = lp <= r
        res4 = lp <= rp

        assert_array_equal(res1, res2.array)
        assert_array_equal(res1, res3.array)
        assert_array_equal(res1, res4.array)
        assert_(isinstance(res1, np.ndarray))
        assert_(isinstance(res2, PriorityNdarray))
        assert_(isinstance(res3, PriorityNdarray))
        assert_(isinstance(res4, PriorityNdarray))

    def test_ge(self):
        l = np.asarray([0., -1., 1.], dtype=dtype)
        r = np.asarray([0., 1., -1.], dtype=dtype)
        lp = PriorityNdarray(l)
        rp = PriorityNdarray(r)
        res1 = l >= r
        res2 = l >= rp
        res3 = lp >= r
        res4 = lp >= rp

        assert_array_equal(res1, res2.array)
        assert_array_equal(res1, res3.array)
        assert_array_equal(res1, res4.array)
        assert_(isinstance(res1, np.ndarray))
        assert_(isinstance(res2, PriorityNdarray))
        assert_(isinstance(res3, PriorityNdarray))
        assert_(isinstance(res4, PriorityNdarray))

    def test_eq(self):
        l = np.asarray([0., -1., 1.], dtype=dtype)
        r = np.asarray([0., 1., -1.], dtype=dtype)
        lp = PriorityNdarray(l)
        rp = PriorityNdarray(r)
        res1 = l == r
        res2 = l == rp
        res3 = lp == r
        res4 = lp == rp

        assert_array_equal(res1, res2.array)
        assert_array_equal(res1, res3.array)
        assert_array_equal(res1, res4.array)
        assert_(isinstance(res1, np.ndarray))
        assert_(isinstance(res2, PriorityNdarray))
        assert_(isinstance(res3, PriorityNdarray))
        assert_(isinstance(res4, PriorityNdarray))

    def test_ne(self):
        l = np.asarray([0., -1., 1.], dtype=dtype)
        r = np.asarray([0., 1., -1.], dtype=dtype)
        lp = PriorityNdarray(l)
        rp = PriorityNdarray(r)
        res1 = l != r
        res2 = l != rp
        res3 = lp != r
        res4 = lp != rp

        assert_array_equal(res1, res2.array)
        assert_array_equal(res1, res3.array)
        assert_array_equal(res1, res4.array)
        assert_(isinstance(res1, np.ndarray))
        assert_(isinstance(res2, PriorityNdarray))
        assert_(isinstance(res3, PriorityNdarray))
        assert_(isinstance(res4, PriorityNdarray))


class TestConversion(TestCase):
    def test_array_scalar_relational_operation(self):
        #All integer
        for dt1 in np.typecodes['AllInteger']:
            assert_(1 > np.array(0, dtype=dt1), "type %s failed" % (dt1,))
            assert_(not 1 < np.array(0, dtype=dt1), "type %s failed" % (dt1,))

            for dt2 in np.typecodes['AllInteger']:
                assert_(np.array(1, dtype=dt1) > np.array(0, dtype=dt2),
                        "type %s and %s failed" % (dt1, dt2))
                assert_(not np.array(1, dtype=dt1) < np.array(0, dtype=dt2),
                        "type %s and %s failed" % (dt1, dt2))

        #Unsigned integers
        for dt1 in 'BHILQP':
            assert_(-1 < np.array(1, dtype=dt1), "type %s failed" % (dt1,))
            assert_(not -1 > np.array(1, dtype=dt1), "type %s failed" % (dt1,))
            assert_(-1 != np.array(1, dtype=dt1), "type %s failed" % (dt1,))

            #unsigned vs signed
            for dt2 in 'bhilqp':
                assert_(np.array(1, dtype=dt1) > np.array(-1, dtype=dt2),
                        "type %s and %s failed" % (dt1, dt2))
                assert_(not np.array(1, dtype=dt1) < np.array(-1, dtype=dt2),
                        "type %s and %s failed" % (dt1, dt2))
                assert_(np.array(1, dtype=dt1) != np.array(-1, dtype=dt2),
                        "type %s and %s failed" % (dt1, dt2))

        #Signed integers and floats
        for dt1 in 'bhlqp' + np.typecodes['Float']:
            assert_(1 > np.array(-1, dtype=dt1), "type %s failed" % (dt1,))
            assert_(not 1 < np.array(-1, dtype=dt1), "type %s failed" % (dt1,))
            assert_(-1 == np.array(-1, dtype=dt1), "type %s failed" % (dt1,))

            for dt2 in 'bhlqp' + np.typecodes['Float']:
                assert_(np.array(1, dtype=dt1) > np.array(-1, dtype=dt2),
                        "type %s and %s failed" % (dt1, dt2))
                assert_(not np.array(1, dtype=dt1) < np.array(-1, dtype=dt2),
                        "type %s and %s failed" % (dt1, dt2))
                assert_(np.array(-1, dtype=dt1) == np.array(-1, dtype=dt2),
                        "type %s and %s failed" % (dt1, dt2))


class TestWhere(TestCase):
    def test_basic(self):
        dts = [np.bool, np.int16, np.int32, np.int64, np.double, np.complex128,
               np.longdouble, np.clongdouble]
        for dt in dts:
            c = np.ones(53, dtype=np.bool)
            assert_equal(np.where( c, dt(0), dt(1)), dt(0))
            assert_equal(np.where(~c, dt(0), dt(1)), dt(1))
            assert_equal(np.where(True, dt(0), dt(1)), dt(0))
            assert_equal(np.where(False, dt(0), dt(1)), dt(1))
            d = np.ones_like(c).astype(dt)
            e = np.zeros_like(d)
            r = d.astype(dt)
            c[7] = False
            r[7] = e[7]
            assert_equal(np.where(c, e, e), e)
            assert_equal(np.where(c, d, e), r)
            assert_equal(np.where(c, d, e[0]), r)
            assert_equal(np.where(c, d[0], e), r)
            assert_equal(np.where(c[::2], d[::2], e[::2]), r[::2])
            assert_equal(np.where(c[1::2], d[1::2], e[1::2]), r[1::2])
            assert_equal(np.where(c[::3], d[::3], e[::3]), r[::3])
            assert_equal(np.where(c[1::3], d[1::3], e[1::3]), r[1::3])
            assert_equal(np.where(c[::-2], d[::-2], e[::-2]), r[::-2])
            assert_equal(np.where(c[::-3], d[::-3], e[::-3]), r[::-3])
            assert_equal(np.where(c[1::-3], d[1::-3], e[1::-3]), r[1::-3])

    def test_exotic(self):
        # object
        assert_array_equal(np.where(True, None, None), np.array(None))
        # zero sized
        m = np.array([], dtype=bool).reshape(0, 3)
        b = np.array([], dtype=np.float64).reshape(0, 3)
        assert_array_equal(np.where(m, 0, b), np.array([]).reshape(0, 3))

        # object cast
        d = np.array([-1.34, -0.16, -0.54, -0.31, -0.08, -0.95, 0.000, 0.313,
                      0.547, -0.18, 0.876, 0.236, 1.969, 0.310, 0.699, 1.013,
                      1.267, 0.229, -1.39, 0.487])
        nan = float('NaN')
        e = np.array(['5z', '0l', nan, 'Wz', nan, nan, 'Xq', 'cs', nan, nan,
                     'QN', nan, nan, 'Fd', nan, nan, 'kp', nan, '36', 'i1'],
                     dtype=object);
        m = np.array([0,0,1,0,1,1,0,0,1,1,0,1,1,0,1,1,0,1,0,0], dtype=bool)

        r = e[:]
        r[np.where(m)] = d[np.where(m)]
        assert_array_equal(np.where(m, d, e), r)

        r = e[:]
        r[np.where(~m)] = d[np.where(~m)]
        assert_array_equal(np.where(m, e, d), r)

        assert_array_equal(np.where(m, e, e), e)

        # minimal dtype result with NaN scalar (e.g required by pandas)
        d =  np.array([1., 2.], dtype=np.float32)
        e =  float('NaN')
        assert_equal(np.where(True, d, e).dtype, np.float32)
        e =  float('Infinity')
        assert_equal(np.where(True, d, e).dtype, np.float32)
        e =  float('-Infinity')
        assert_equal(np.where(True, d, e).dtype, np.float32)
        # also check upcast
        e =  float(1e150)
        assert_equal(np.where(True, d, e).dtype, np.float64)

    def test_ndim(self):
        c = [True, False]
        a = np.zeros((2, 25))
        b = np.ones((2, 25))
        r = np.where(np.array(c)[:,np.newaxis], a, b)
        assert_array_equal(r[0], a[0])
        assert_array_equal(r[1], b[0])

        a = a.T
        b = b.T
        r = np.where(c, a, b)
        assert_array_equal(r[:,0], a[:,0])
        assert_array_equal(r[:,1], b[:,0])

    def test_dtype_mix(self):
        c = np.array([False, True, False, False, False, False, True, False,
                     False, False, True, False])
        a = np.uint32(1)
        b =  np.array([5., 0., 3., 2., -1., -4., 0., -10., 10., 1., 0., 3.],
                      dtype=np.float64)
        r = np.array([5., 1., 3., 2., -1., -4., 1., -10., 10., 1., 1., 3.],
                     dtype=np.float64)
        assert_equal(np.where(c, a, b), r)

        a = a.astype(np.float32)
        b = b.astype(np.int64)
        assert_equal(np.where(c, a, b), r)

        # non bool mask
        c = c.astype(np.int)
        c[c != 0] = 34242324
        assert_equal(np.where(c, a, b), r)
        # invert
        tmpmask = c != 0
        c[c == 0] = 41247212
        c[tmpmask] = 0
        assert_equal(np.where(c, b, a), r)

    def test_foreign(self):
        c = np.array([False, True, False, False, False, False, True, False,
                     False, False, True, False])
        r = np.array([5., 1., 3., 2., -1., -4., 1., -10., 10., 1., 1., 3.],
                     dtype=np.float64)
        a = np.ones(1, dtype='>i4')
        b =  np.array([5., 0., 3., 2., -1., -4., 0., -10., 10., 1., 0., 3.],
                      dtype=np.float64)
        assert_equal(np.where(c, a, b), r)

        b = b.astype('>f8')
        assert_equal(np.where(c, a, b), r)

        a = a.astype('<i4')
        assert_equal(np.where(c, a, b), r)

        c = c.astype('>i4')
        assert_equal(np.where(c, a, b), r)

    def test_error(self):
        c = [True, True]
        a = np.ones((4, 5))
        b = np.ones((5, 5))
        assert_raises(ValueError, np.where, c, a, a)
        assert_raises(ValueError, np.where, c[0], a, b)

    def test_string(self):
        # gh-4778 check strings are properly filled with nulls
        a = np.array("abc")
        b = np.array("x" * 753)
        assert_equal(np.where(True, a, b), "abc")
        assert_equal(np.where(False, b, a), "abc")

        # check native datatype sized strings
        a = np.array("abcd")
        b = np.array("x" * 8)
        assert_equal(np.where(True, a, b), "abcd")
        assert_equal(np.where(False, b, a), "abcd")


if __name__ == "__main__":
    run_module_suite()
