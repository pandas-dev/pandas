# pylint: disable-msg=W0401,W0511,W0611,W0612,W0614,R0201,E1102
"""Tests suite for MaskedArray & subclassing.

:author: Pierre Gerard-Marchant
:contact: pierregm_at_uga_dot_edu
"""
from __future__ import division, absolute_import, print_function

__author__ = "Pierre GF Gerard-Marchant"

import warnings
import sys
import pickle
from functools import reduce

from nose.tools import assert_raises

import numpy as np
import numpy.ma.core
import numpy.core.fromnumeric as fromnumeric
from numpy import ndarray
from numpy.ma.testutils import *
from numpy.ma.core import *
from numpy.compat import asbytes, asbytes_nested

pi = np.pi


#..............................................................................
class TestMaskedArray(TestCase):
    # Base test class for MaskedArrays.

    def setUp(self):
        # Base data definition.
        x = np.array([1., 1., 1., -2., pi/2.0, 4., 5., -10., 10., 1., 2., 3.])
        y = np.array([5., 0., 3., 2., -1., -4., 0., -10., 10., 1., 0., 3.])
        a10 = 10.
        m1 = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]
        m2 = [0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1]
        xm = masked_array(x, mask=m1)
        ym = masked_array(y, mask=m2)
        z = np.array([-.5, 0., .5, .8])
        zm = masked_array(z, mask=[0, 1, 0, 0])
        xf = np.where(m1, 1e+20, x)
        xm.set_fill_value(1e+20)
        self.d = (x, y, a10, m1, m2, xm, ym, z, zm, xf)

    def test_basicattributes(self):
        # Tests some basic array attributes.
        a = array([1, 3, 2])
        b = array([1, 3, 2], mask=[1, 0, 1])
        assert_equal(a.ndim, 1)
        assert_equal(b.ndim, 1)
        assert_equal(a.size, 3)
        assert_equal(b.size, 3)
        assert_equal(a.shape, (3,))
        assert_equal(b.shape, (3,))

    def test_basic0d(self):
        # Checks masking a scalar
        x = masked_array(0)
        assert_equal(str(x), '0')
        x = masked_array(0, mask=True)
        assert_equal(str(x), str(masked_print_option))
        x = masked_array(0, mask=False)
        assert_equal(str(x), '0')
        x = array(0, mask=1)
        self.assertTrue(x.filled().dtype is x._data.dtype)

    def test_basic1d(self):
        # Test of basic array creation and properties in 1 dimension.
        (x, y, a10, m1, m2, xm, ym, z, zm, xf) = self.d
        self.assertTrue(not isMaskedArray(x))
        self.assertTrue(isMaskedArray(xm))
        self.assertTrue((xm - ym).filled(0).any())
        fail_if_equal(xm.mask.astype(int), ym.mask.astype(int))
        s = x.shape
        assert_equal(np.shape(xm), s)
        assert_equal(xm.shape, s)
        assert_equal(xm.dtype, x.dtype)
        assert_equal(zm.dtype, z.dtype)
        assert_equal(xm.size, reduce(lambda x, y:x * y, s))
        assert_equal(count(xm), len(m1) - reduce(lambda x, y:x + y, m1))
        assert_array_equal(xm, xf)
        assert_array_equal(filled(xm, 1.e20), xf)
        assert_array_equal(x, xm)

    def test_basic2d(self):
        # Test of basic array creation and properties in 2 dimensions.
        (x, y, a10, m1, m2, xm, ym, z, zm, xf) = self.d
        for s in [(4, 3), (6, 2)]:
            x.shape = s
            y.shape = s
            xm.shape = s
            ym.shape = s
            xf.shape = s
            #
            self.assertTrue(not isMaskedArray(x))
            self.assertTrue(isMaskedArray(xm))
            assert_equal(shape(xm), s)
            assert_equal(xm.shape, s)
            assert_equal(xm.size, reduce(lambda x, y:x * y, s))
            assert_equal(count(xm), len(m1) - reduce(lambda x, y:x + y, m1))
            assert_equal(xm, xf)
            assert_equal(filled(xm, 1.e20), xf)
            assert_equal(x, xm)

    def test_concatenate_basic(self):
        # Tests concatenations.
        (x, y, a10, m1, m2, xm, ym, z, zm, xf) = self.d
        # basic concatenation
        assert_equal(np.concatenate((x, y)), concatenate((xm, ym)))
        assert_equal(np.concatenate((x, y)), concatenate((x, y)))
        assert_equal(np.concatenate((x, y)), concatenate((xm, y)))
        assert_equal(np.concatenate((x, y, x)), concatenate((x, ym, x)))

    def test_concatenate_alongaxis(self):
        # Tests concatenations.
        (x, y, a10, m1, m2, xm, ym, z, zm, xf) = self.d
        # Concatenation along an axis
        s = (3, 4)
        x.shape = y.shape = xm.shape = ym.shape = s
        assert_equal(xm.mask, np.reshape(m1, s))
        assert_equal(ym.mask, np.reshape(m2, s))
        xmym = concatenate((xm, ym), 1)
        assert_equal(np.concatenate((x, y), 1), xmym)
        assert_equal(np.concatenate((xm.mask, ym.mask), 1), xmym._mask)
        #
        x = zeros(2)
        y = array(ones(2), mask=[False, True])
        z = concatenate((x, y))
        assert_array_equal(z, [0, 0, 1, 1])
        assert_array_equal(z.mask, [False, False, False, True])
        z = concatenate((y, x))
        assert_array_equal(z, [1, 1, 0, 0])
        assert_array_equal(z.mask, [False, True, False, False])

    def test_concatenate_flexible(self):
        # Tests the concatenation on flexible arrays.
        data = masked_array(list(zip(np.random.rand(10),
                                     np.arange(10))),
                            dtype=[('a', float), ('b', int)])
        #
        test = concatenate([data[:5], data[5:]])
        assert_equal_records(test, data)

    def test_creation_ndmin(self):
        # Check the use of ndmin
        x = array([1, 2, 3], mask=[1, 0, 0], ndmin=2)
        assert_equal(x.shape, (1, 3))
        assert_equal(x._data, [[1, 2, 3]])
        assert_equal(x._mask, [[1, 0, 0]])

    def test_creation_ndmin_from_maskedarray(self):
        # Make sure we're not losing the original mask w/ ndmin
        x = array([1, 2, 3])
        x[-1] = masked
        xx = array(x, ndmin=2, dtype=float)
        assert_equal(x.shape, x._mask.shape)
        assert_equal(xx.shape, xx._mask.shape)

    def test_creation_maskcreation(self):
        # Tests how masks are initialized at the creation of Maskedarrays.
        data = arange(24, dtype=float)
        data[[3, 6, 15]] = masked
        dma_1 = MaskedArray(data)
        assert_equal(dma_1.mask, data.mask)
        dma_2 = MaskedArray(dma_1)
        assert_equal(dma_2.mask, dma_1.mask)
        dma_3 = MaskedArray(dma_1, mask=[1, 0, 0, 0] * 6)
        fail_if_equal(dma_3.mask, dma_1.mask)

    def test_creation_with_list_of_maskedarrays(self):
        # Tests creaating a masked array from alist of masked arrays.
        x = array(np.arange(5), mask=[1, 0, 0, 0, 0])
        data = array((x, x[::-1]))
        assert_equal(data, [[0, 1, 2, 3, 4], [4, 3, 2, 1, 0]])
        assert_equal(data._mask, [[1, 0, 0, 0, 0], [0, 0, 0, 0, 1]])
        #
        x.mask = nomask
        data = array((x, x[::-1]))
        assert_equal(data, [[0, 1, 2, 3, 4], [4, 3, 2, 1, 0]])
        self.assertTrue(data.mask is nomask)

    def test_asarray(self):
        (x, y, a10, m1, m2, xm, ym, z, zm, xf) = self.d
        xm.fill_value = -9999
        xm._hardmask = True
        xmm = asarray(xm)
        assert_equal(xmm._data, xm._data)
        assert_equal(xmm._mask, xm._mask)
        assert_equal(xmm.fill_value, xm.fill_value)
        assert_equal(xmm._hardmask, xm._hardmask)

    def test_fix_invalid(self):
        # Checks fix_invalid.
        with np.errstate(invalid='ignore'):
            data = masked_array([np.nan, 0., 1.], mask=[0, 0, 1])
            data_fixed = fix_invalid(data)
            assert_equal(data_fixed._data, [data.fill_value, 0., 1.])
            assert_equal(data_fixed._mask, [1., 0., 1.])

    def test_maskedelement(self):
        # Test of masked element
        x = arange(6)
        x[1] = masked
        self.assertTrue(str(masked) == '--')
        self.assertTrue(x[1] is masked)
        assert_equal(filled(x[1], 0), 0)
        # don't know why these should raise an exception...
        #self.assertRaises(Exception, lambda x,y: x+y, masked, masked)
        #self.assertRaises(Exception, lambda x,y: x+y, masked, 2)
        #self.assertRaises(Exception, lambda x,y: x+y, masked, xx)
        #self.assertRaises(Exception, lambda x,y: x+y, xx, masked)

    def test_set_element_as_object(self):
        # Tests setting elements with object
        a = empty(1, dtype=object)
        x = (1, 2, 3, 4, 5)
        a[0] = x
        assert_equal(a[0], x)
        self.assertTrue(a[0] is x)
        #
        import datetime
        dt = datetime.datetime.now()
        a[0] = dt
        self.assertTrue(a[0] is dt)

    def test_indexing(self):
        # Tests conversions and indexing
        x1 = np.array([1, 2, 4, 3])
        x2 = array(x1, mask=[1, 0, 0, 0])
        x3 = array(x1, mask=[0, 1, 0, 1])
        x4 = array(x1)
        # test conversion to strings
        junk, garbage = str(x2), repr(x2)
        assert_equal(np.sort(x1), sort(x2, endwith=False))
        # tests of indexing
        assert_(type(x2[1]) is type(x1[1]))
        assert_(x1[1] == x2[1])
        assert_(x2[0] is masked)
        assert_equal(x1[2], x2[2])
        assert_equal(x1[2:5], x2[2:5])
        assert_equal(x1[:], x2[:])
        assert_equal(x1[1:], x3[1:])
        x1[2] = 9
        x2[2] = 9
        assert_equal(x1, x2)
        x1[1:3] = 99
        x2[1:3] = 99
        assert_equal(x1, x2)
        x2[1] = masked
        assert_equal(x1, x2)
        x2[1:3] = masked
        assert_equal(x1, x2)
        x2[:] = x1
        x2[1] = masked
        assert_(allequal(getmask(x2), array([0, 1, 0, 0])))
        x3[:] = masked_array([1, 2, 3, 4], [0, 1, 1, 0])
        assert_(allequal(getmask(x3), array([0, 1, 1, 0])))
        x4[:] = masked_array([1, 2, 3, 4], [0, 1, 1, 0])
        assert_(allequal(getmask(x4), array([0, 1, 1, 0])))
        assert_(allequal(x4, array([1, 2, 3, 4])))
        x1 = np.arange(5) * 1.0
        x2 = masked_values(x1, 3.0)
        assert_equal(x1, x2)
        assert_(allequal(array([0, 0, 0, 1, 0], MaskType), x2.mask))
        assert_equal(3.0, x2.fill_value)
        x1 = array([1, 'hello', 2, 3], object)
        x2 = np.array([1, 'hello', 2, 3], object)
        s1 = x1[1]
        s2 = x2[1]
        assert_equal(type(s2), str)
        assert_equal(type(s1), str)
        assert_equal(s1, s2)
        assert_(x1[1:1].shape == (0,))

    def test_copy(self):
        # Tests of some subtle points of copying and sizing.
        n = [0, 0, 1, 0, 0]
        m = make_mask(n)
        m2 = make_mask(m)
        self.assertTrue(m is m2)
        m3 = make_mask(m, copy=1)
        self.assertTrue(m is not m3)

        x1 = np.arange(5)
        y1 = array(x1, mask=m)
        #self.assertTrue( y1._data is x1)
        assert_equal(y1._data.__array_interface__, x1.__array_interface__)
        self.assertTrue(allequal(x1, y1.data))
        #self.assertTrue( y1.mask is m)
        assert_equal(y1._mask.__array_interface__, m.__array_interface__)

        y1a = array(y1)
        self.assertTrue(y1a._data.__array_interface__ ==
                        y1._data.__array_interface__)
        self.assertTrue(y1a.mask is y1.mask)

        y2 = array(x1, mask=m)
        self.assertTrue(y2._data.__array_interface__ == x1.__array_interface__)
        #self.assertTrue( y2.mask is m)
        self.assertTrue(y2._mask.__array_interface__ == m.__array_interface__)
        self.assertTrue(y2[2] is masked)
        y2[2] = 9
        self.assertTrue(y2[2] is not masked)
        #self.assertTrue( y2.mask is not m)
        self.assertTrue(y2._mask.__array_interface__ != m.__array_interface__)
        self.assertTrue(allequal(y2.mask, 0))

        y3 = array(x1 * 1.0, mask=m)
        self.assertTrue(filled(y3).dtype is (x1 * 1.0).dtype)

        x4 = arange(4)
        x4[2] = masked
        y4 = resize(x4, (8,))
        assert_equal(concatenate([x4, x4]), y4)
        assert_equal(getmask(y4), [0, 0, 1, 0, 0, 0, 1, 0])
        y5 = repeat(x4, (2, 2, 2, 2), axis=0)
        assert_equal(y5, [0, 0, 1, 1, 2, 2, 3, 3])
        y6 = repeat(x4, 2, axis=0)
        assert_equal(y5, y6)
        y7 = x4.repeat((2, 2, 2, 2), axis=0)
        assert_equal(y5, y7)
        y8 = x4.repeat(2, 0)
        assert_equal(y5, y8)

        y9 = x4.copy()
        assert_equal(y9._data, x4._data)
        assert_equal(y9._mask, x4._mask)
        #
        x = masked_array([1, 2, 3], mask=[0, 1, 0])
        # Copy is False by default
        y = masked_array(x)
        assert_equal(y._data.ctypes.data, x._data.ctypes.data)
        assert_equal(y._mask.ctypes.data, x._mask.ctypes.data)
        y = masked_array(x, copy=True)
        assert_not_equal(y._data.ctypes.data, x._data.ctypes.data)
        assert_not_equal(y._mask.ctypes.data, x._mask.ctypes.data)

    def test_deepcopy(self):
        from copy import deepcopy
        a = array([0, 1, 2], mask=[False, True, False])
        copied = deepcopy(a)
        assert_equal(copied.mask, a.mask)
        assert_not_equal(id(a._mask), id(copied._mask))
        #
        copied[1] = 1
        assert_equal(copied.mask, [0, 0, 0])
        assert_equal(a.mask, [0, 1, 0])
        #
        copied = deepcopy(a)
        assert_equal(copied.mask, a.mask)
        copied.mask[1] = False
        assert_equal(copied.mask, [0, 0, 0])
        assert_equal(a.mask, [0, 1, 0])

    def test_str_repr(self):
        a = array([0, 1, 2], mask=[False, True, False])
        assert_equal(str(a), '[0 -- 2]')
        assert_equal(repr(a), 'masked_array(data = [0 -- 2],\n'
                              '             mask = [False  True False],\n'
                              '       fill_value = 999999)\n')

    def test_pickling(self):
        # Tests pickling
        a = arange(10)
        a[::3] = masked
        a.fill_value = 999
        a_pickled = pickle.loads(a.dumps())
        assert_equal(a_pickled._mask, a._mask)
        assert_equal(a_pickled._data, a._data)
        assert_equal(a_pickled.fill_value, 999)

    def test_pickling_subbaseclass(self):
        # Test pickling w/ a subclass of ndarray
        a = array(np.matrix(list(range(10))), mask=[1, 0, 1, 0, 0] * 2)
        a_pickled = pickle.loads(a.dumps())
        assert_equal(a_pickled._mask, a._mask)
        assert_equal(a_pickled, a)
        self.assertTrue(isinstance(a_pickled._data, np.matrix))

    def test_pickling_maskedconstant(self):
        # Test pickling MaskedConstant
        mc = np.ma.masked
        mc_pickled = pickle.loads(mc.dumps())
        assert_equal(mc_pickled._baseclass, mc._baseclass)
        assert_equal(mc_pickled._mask, mc._mask)
        assert_equal(mc_pickled._data, mc._data)

    def test_pickling_wstructured(self):
        # Tests pickling w/ structured array
        a = array([(1, 1.), (2, 2.)], mask=[(0, 0), (0, 1)],
                  dtype=[('a', int), ('b', float)])
        a_pickled = pickle.loads(a.dumps())
        assert_equal(a_pickled._mask, a._mask)
        assert_equal(a_pickled, a)

    def test_pickling_keepalignment(self):
        # Tests pickling w/ F_CONTIGUOUS arrays
        a = arange(10)
        a.shape = (-1, 2)
        b = a.T
        test = pickle.loads(pickle.dumps(b))
        assert_equal(test, b)

    def test_single_element_subscript(self):
        # Tests single element subscripts of Maskedarrays.
        a = array([1, 3, 2])
        b = array([1, 3, 2], mask=[1, 0, 1])
        assert_equal(a[0].shape, ())
        assert_equal(b[0].shape, ())
        assert_equal(b[1].shape, ())

    def test_topython(self):
        # Tests some communication issues with Python.
        assert_equal(1, int(array(1)))
        assert_equal(1.0, float(array(1)))
        assert_equal(1, int(array([[[1]]])))
        assert_equal(1.0, float(array([[1]])))
        self.assertRaises(TypeError, float, array([1, 1]))
        #
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            assert_(np.isnan(float(array([1], mask=[1]))))
        #
        a = array([1, 2, 3], mask=[1, 0, 0])
        self.assertRaises(TypeError, lambda:float(a))
        assert_equal(float(a[-1]), 3.)
        self.assertTrue(np.isnan(float(a[0])))
        self.assertRaises(TypeError, int, a)
        assert_equal(int(a[-1]), 3)
        self.assertRaises(MAError, lambda:int(a[0]))

    def test_oddfeatures_1(self):
        # Test of other odd features
        x = arange(20)
        x = x.reshape(4, 5)
        x.flat[5] = 12
        assert_(x[1, 0] == 12)
        z = x + 10j * x
        assert_equal(z.real, x)
        assert_equal(z.imag, 10 * x)
        assert_equal((z * conjugate(z)).real, 101 * x * x)
        z.imag[...] = 0.0
        #
        x = arange(10)
        x[3] = masked
        assert_(str(x[3]) == str(masked))
        c = x >= 8
        assert_(count(where(c, masked, masked)) == 0)
        assert_(shape(where(c, masked, masked)) == c.shape)
        #
        z = masked_where(c, x)
        assert_(z.dtype is x.dtype)
        assert_(z[3] is masked)
        assert_(z[4] is not masked)
        assert_(z[7] is not masked)
        assert_(z[8] is masked)
        assert_(z[9] is masked)
        assert_equal(x, z)

    def test_oddfeatures_2(self):
        # Tests some more features.
        x = array([1., 2., 3., 4., 5.])
        c = array([1, 1, 1, 0, 0])
        x[2] = masked
        z = where(c, x, -x)
        assert_equal(z, [1., 2., 0., -4., -5])
        c[0] = masked
        z = where(c, x, -x)
        assert_equal(z, [1., 2., 0., -4., -5])
        assert_(z[0] is masked)
        assert_(z[1] is not masked)
        assert_(z[2] is masked)

    def test_oddfeatures_3(self):
        # Tests some generic features
        atest = array([10], mask=True)
        btest = array([20])
        idx = atest.mask
        atest[idx] = btest[idx]
        assert_equal(atest, [20])

    def test_filled_w_object_dtype(self):
        a = np.ma.masked_all(1, dtype='O')
        assert_equal(a.filled('x')[0], 'x')

    def test_filled_w_flexible_dtype(self):
        # Test filled w/ flexible dtype
        flexi = array([(1, 1, 1)],
                      dtype=[('i', int), ('s', '|S8'), ('f', float)])
        flexi[0] = masked
        assert_equal(flexi.filled(),
                     np.array([(default_fill_value(0),
                                default_fill_value('0'),
                                default_fill_value(0.),)], dtype=flexi.dtype))
        flexi[0] = masked
        assert_equal(flexi.filled(1),
                     np.array([(1, '1', 1.)], dtype=flexi.dtype))

    def test_filled_w_mvoid(self):
        # Test filled w/ mvoid
        ndtype = [('a', int), ('b', float)]
        a = mvoid((1, 2.), mask=[(0, 1)], dtype=ndtype)
        # Filled using default
        test = a.filled()
        assert_equal(tuple(test), (1, default_fill_value(1.)))
        # Explicit fill_value
        test = a.filled((-1, -1))
        assert_equal(tuple(test), (1, -1))
        # Using predefined filling values
        a.fill_value = (-999, -999)
        assert_equal(tuple(a.filled()), (1, -999))

    def test_filled_w_nested_dtype(self):
        # Test filled w/ nested dtype
        ndtype = [('A', int), ('B', [('BA', int), ('BB', int)])]
        a = array([(1, (1, 1)), (2, (2, 2))],
                  mask=[(0, (1, 0)), (0, (0, 1))], dtype=ndtype)
        test = a.filled(0)
        control = np.array([(1, (0, 1)), (2, (2, 0))], dtype=ndtype)
        assert_equal(test, control)
        #
        test = a['B'].filled(0)
        control = np.array([(0, 1), (2, 0)], dtype=a['B'].dtype)
        assert_equal(test, control)

    def test_filled_w_f_order(self):
        # Test filled w/ F-contiguous array
        a = array(np.array([(0, 1, 2), (4, 5, 6)], order='F'),
                  mask=np.array([(0, 0, 1), (1, 0, 0)], order='F'),
                  order='F')  # this is currently ignored
        self.assertTrue(a.flags['F_CONTIGUOUS'])
        self.assertTrue(a.filled(0).flags['F_CONTIGUOUS'])

    def test_optinfo_propagation(self):
        # Checks that _optinfo dictionary isn't back-propagated
        x = array([1, 2, 3, ], dtype=float)
        x._optinfo['info'] = '???'
        y = x.copy()
        assert_equal(y._optinfo['info'], '???')
        y._optinfo['info'] = '!!!'
        assert_equal(x._optinfo['info'], '???')

    def test_fancy_printoptions(self):
        # Test printing a masked array w/ fancy dtype.
        fancydtype = np.dtype([('x', int), ('y', [('t', int), ('s', float)])])
        test = array([(1, (2, 3.0)), (4, (5, 6.0))],
                     mask=[(1, (0, 1)), (0, (1, 0))],
                     dtype=fancydtype)
        control = "[(--, (2, --)) (4, (--, 6.0))]"
        assert_equal(str(test), control)

    def test_flatten_structured_array(self):
        # Test flatten_structured_array on arrays
        # On ndarray
        ndtype = [('a', int), ('b', float)]
        a = np.array([(1, 1), (2, 2)], dtype=ndtype)
        test = flatten_structured_array(a)
        control = np.array([[1., 1.], [2., 2.]], dtype=np.float)
        assert_equal(test, control)
        assert_equal(test.dtype, control.dtype)
        # On masked_array
        a = array([(1, 1), (2, 2)], mask=[(0, 1), (1, 0)], dtype=ndtype)
        test = flatten_structured_array(a)
        control = array([[1., 1.], [2., 2.]],
                        mask=[[0, 1], [1, 0]], dtype=np.float)
        assert_equal(test, control)
        assert_equal(test.dtype, control.dtype)
        assert_equal(test.mask, control.mask)
        # On masked array with nested structure
        ndtype = [('a', int), ('b', [('ba', int), ('bb', float)])]
        a = array([(1, (1, 1.1)), (2, (2, 2.2))],
                  mask=[(0, (1, 0)), (1, (0, 1))], dtype=ndtype)
        test = flatten_structured_array(a)
        control = array([[1., 1., 1.1], [2., 2., 2.2]],
                        mask=[[0, 1, 0], [1, 0, 1]], dtype=np.float)
        assert_equal(test, control)
        assert_equal(test.dtype, control.dtype)
        assert_equal(test.mask, control.mask)
        # Keeping the initial shape
        ndtype = [('a', int), ('b', float)]
        a = np.array([[(1, 1), ], [(2, 2), ]], dtype=ndtype)
        test = flatten_structured_array(a)
        control = np.array([[[1., 1.], ], [[2., 2.], ]], dtype=np.float)
        assert_equal(test, control)
        assert_equal(test.dtype, control.dtype)

    def test_void0d(self):
        # Test creating a mvoid object
        ndtype = [('a', int), ('b', int)]
        a = np.array([(1, 2,)], dtype=ndtype)[0]
        f = mvoid(a)
        assert_(isinstance(f, mvoid))
        #
        a = masked_array([(1, 2)], mask=[(1, 0)], dtype=ndtype)[0]
        assert_(isinstance(a, mvoid))
        #
        a = masked_array([(1, 2), (1, 2)], mask=[(1, 0), (0, 0)], dtype=ndtype)
        f = mvoid(a._data[0], a._mask[0])
        assert_(isinstance(f, mvoid))

    def test_mvoid_getitem(self):
        # Test mvoid.__getitem__
        ndtype = [('a', int), ('b', int)]
        a = masked_array([(1, 2,), (3, 4)], mask=[(0, 0), (1, 0)],
                         dtype=ndtype)
        # w/o mask
        f = a[0]
        self.assertTrue(isinstance(f, mvoid))
        assert_equal((f[0], f['a']), (1, 1))
        assert_equal(f['b'], 2)
        # w/ mask
        f = a[1]
        self.assertTrue(isinstance(f, mvoid))
        self.assertTrue(f[0] is masked)
        self.assertTrue(f['a'] is masked)
        assert_equal(f[1], 4)

    def test_mvoid_iter(self):
        # Test iteration on __getitem__
        ndtype = [('a', int), ('b', int)]
        a = masked_array([(1, 2,), (3, 4)], mask=[(0, 0), (1, 0)],
                         dtype=ndtype)
        # w/o mask
        assert_equal(list(a[0]), [1, 2])
        # w/ mask
        assert_equal(list(a[1]), [masked, 4])

    def test_mvoid_print(self):
        # Test printing a mvoid
        mx = array([(1, 1), (2, 2)], dtype=[('a', int), ('b', int)])
        assert_equal(str(mx[0]), "(1, 1)")
        mx['b'][0] = masked
        ini_display = masked_print_option._display
        masked_print_option.set_display("-X-")
        try:
            assert_equal(str(mx[0]), "(1, -X-)")
            assert_equal(repr(mx[0]), "(1, -X-)")
        finally:
            masked_print_option.set_display(ini_display)


#------------------------------------------------------------------------------
class TestMaskedArrayArithmetic(TestCase):
    # Base test class for MaskedArrays.

    def setUp(self):
        # Base data definition.
        x = np.array([1., 1., 1., -2., pi/2.0, 4., 5., -10., 10., 1., 2., 3.])
        y = np.array([5., 0., 3., 2., -1., -4., 0., -10., 10., 1., 0., 3.])
        a10 = 10.
        m1 = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]
        m2 = [0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1]
        xm = masked_array(x, mask=m1)
        ym = masked_array(y, mask=m2)
        z = np.array([-.5, 0., .5, .8])
        zm = masked_array(z, mask=[0, 1, 0, 0])
        xf = np.where(m1, 1e+20, x)
        xm.set_fill_value(1e+20)
        self.d = (x, y, a10, m1, m2, xm, ym, z, zm, xf)
        self.err_status = np.geterr()
        np.seterr(divide='ignore', invalid='ignore')

    def tearDown(self):
        np.seterr(**self.err_status)

    def test_basic_arithmetic(self):
        # Test of basic arithmetic.
        (x, y, a10, m1, m2, xm, ym, z, zm, xf) = self.d
        a2d = array([[1, 2], [0, 4]])
        a2dm = masked_array(a2d, [[0, 0], [1, 0]])
        assert_equal(a2d * a2d, a2d * a2dm)
        assert_equal(a2d + a2d, a2d + a2dm)
        assert_equal(a2d - a2d, a2d - a2dm)
        for s in [(12,), (4, 3), (2, 6)]:
            x = x.reshape(s)
            y = y.reshape(s)
            xm = xm.reshape(s)
            ym = ym.reshape(s)
            xf = xf.reshape(s)
            assert_equal(-x, -xm)
            assert_equal(x + y, xm + ym)
            assert_equal(x - y, xm - ym)
            assert_equal(x * y, xm * ym)
            assert_equal(x / y, xm / ym)
            assert_equal(a10 + y, a10 + ym)
            assert_equal(a10 - y, a10 - ym)
            assert_equal(a10 * y, a10 * ym)
            assert_equal(a10 / y, a10 / ym)
            assert_equal(x + a10, xm + a10)
            assert_equal(x - a10, xm - a10)
            assert_equal(x * a10, xm * a10)
            assert_equal(x / a10, xm / a10)
            assert_equal(x ** 2, xm ** 2)
            assert_equal(abs(x) ** 2.5, abs(xm) ** 2.5)
            assert_equal(x ** y, xm ** ym)
            assert_equal(np.add(x, y), add(xm, ym))
            assert_equal(np.subtract(x, y), subtract(xm, ym))
            assert_equal(np.multiply(x, y), multiply(xm, ym))
            assert_equal(np.divide(x, y), divide(xm, ym))

    def test_divide_on_different_shapes(self):
        x = arange(6, dtype=float)
        x.shape = (2, 3)
        y = arange(3, dtype=float)
        #
        z = x / y
        assert_equal(z, [[-1., 1., 1.], [-1., 4., 2.5]])
        assert_equal(z.mask, [[1, 0, 0], [1, 0, 0]])
        #
        z = x / y[None,:]
        assert_equal(z, [[-1., 1., 1.], [-1., 4., 2.5]])
        assert_equal(z.mask, [[1, 0, 0], [1, 0, 0]])
        #
        y = arange(2, dtype=float)
        z = x / y[:, None]
        assert_equal(z, [[-1., -1., -1.], [3., 4., 5.]])
        assert_equal(z.mask, [[1, 1, 1], [0, 0, 0]])

    def test_mixed_arithmetic(self):
        # Tests mixed arithmetics.
        na = np.array([1])
        ma = array([1])
        self.assertTrue(isinstance(na + ma, MaskedArray))
        self.assertTrue(isinstance(ma + na, MaskedArray))

    def test_limits_arithmetic(self):
        tiny = np.finfo(float).tiny
        a = array([tiny, 1. / tiny, 0.])
        assert_equal(getmaskarray(a / 2), [0, 0, 0])
        assert_equal(getmaskarray(2 / a), [1, 0, 1])

    def test_masked_singleton_arithmetic(self):
        # Tests some scalar arithmetics on MaskedArrays.
        # Masked singleton should remain masked no matter what
        xm = array(0, mask=1)
        self.assertTrue((1 / array(0)).mask)
        self.assertTrue((1 + xm).mask)
        self.assertTrue((-xm).mask)
        self.assertTrue(maximum(xm, xm).mask)
        self.assertTrue(minimum(xm, xm).mask)

    def test_masked_singleton_equality(self):
        # Tests (in)equality on masked snigleton
        a = array([1, 2, 3], mask=[1, 1, 0])
        assert_((a[0] == 0) is masked)
        assert_((a[0] != 0) is masked)
        assert_equal((a[-1] == 0), False)
        assert_equal((a[-1] != 0), True)

    def test_arithmetic_with_masked_singleton(self):
        # Checks that there's no collapsing to masked
        x = masked_array([1, 2])
        y = x * masked
        assert_equal(y.shape, x.shape)
        assert_equal(y._mask, [True, True])
        y = x[0] * masked
        assert_(y is masked)
        y = x + masked
        assert_equal(y.shape, x.shape)
        assert_equal(y._mask, [True, True])

    def test_arithmetic_with_masked_singleton_on_1d_singleton(self):
        # Check that we're not losing the shape of a singleton
        x = masked_array([1, ])
        y = x + masked
        assert_equal(y.shape, x.shape)
        assert_equal(y.mask, [True, ])

    def test_scalar_arithmetic(self):
        x = array(0, mask=0)
        assert_equal(x.filled().ctypes.data, x.ctypes.data)
        # Make sure we don't lose the shape in some circumstances
        xm = array((0, 0)) / 0.
        assert_equal(xm.shape, (2,))
        assert_equal(xm.mask, [1, 1])

    def test_basic_ufuncs(self):
        # Test various functions such as sin, cos.
        (x, y, a10, m1, m2, xm, ym, z, zm, xf) = self.d
        assert_equal(np.cos(x), cos(xm))
        assert_equal(np.cosh(x), cosh(xm))
        assert_equal(np.sin(x), sin(xm))
        assert_equal(np.sinh(x), sinh(xm))
        assert_equal(np.tan(x), tan(xm))
        assert_equal(np.tanh(x), tanh(xm))
        assert_equal(np.sqrt(abs(x)), sqrt(xm))
        assert_equal(np.log(abs(x)), log(xm))
        assert_equal(np.log10(abs(x)), log10(xm))
        assert_equal(np.exp(x), exp(xm))
        assert_equal(np.arcsin(z), arcsin(zm))
        assert_equal(np.arccos(z), arccos(zm))
        assert_equal(np.arctan(z), arctan(zm))
        assert_equal(np.arctan2(x, y), arctan2(xm, ym))
        assert_equal(np.absolute(x), absolute(xm))
        assert_equal(np.angle(x + 1j*y), angle(xm + 1j*ym))
        assert_equal(np.angle(x + 1j*y, deg=True), angle(xm + 1j*ym, deg=True))
        assert_equal(np.equal(x, y), equal(xm, ym))
        assert_equal(np.not_equal(x, y), not_equal(xm, ym))
        assert_equal(np.less(x, y), less(xm, ym))
        assert_equal(np.greater(x, y), greater(xm, ym))
        assert_equal(np.less_equal(x, y), less_equal(xm, ym))
        assert_equal(np.greater_equal(x, y), greater_equal(xm, ym))
        assert_equal(np.conjugate(x), conjugate(xm))

    def test_count_func(self):
        # Tests count
        assert_equal(1, count(1))
        assert_equal(0, array(1, mask=[1]))

        ott = array([0., 1., 2., 3.], mask=[1, 0, 0, 0])
        res = count(ott)
        self.assertTrue(res.dtype.type is np.intp)
        assert_equal(3, res)

        ott = ott.reshape((2, 2))
        res = count(ott)
        assert_(res.dtype.type is np.intp)
        assert_equal(3, res)
        res = count(ott, 0)
        assert_(isinstance(res, ndarray))
        assert_equal([1, 2], res)
        assert_(getmask(res) is nomask)

        ott= array([0., 1., 2., 3.])
        res = count(ott, 0)
        assert_(isinstance(res, ndarray))
        assert_(res.dtype.type is np.intp)

        assert_raises(IndexError, ott.count, 1)

    def test_minmax_func(self):
        # Tests minimum and maximum.
        (x, y, a10, m1, m2, xm, ym, z, zm, xf) = self.d
        # max doesn't work if shaped
        xr = np.ravel(x)
        xmr = ravel(xm)
        # following are true because of careful selection of data
        assert_equal(max(xr), maximum(xmr))
        assert_equal(min(xr), minimum(xmr))
        #
        assert_equal(minimum([1, 2, 3], [4, 0, 9]), [1, 0, 3])
        assert_equal(maximum([1, 2, 3], [4, 0, 9]), [4, 2, 9])
        x = arange(5)
        y = arange(5) - 2
        x[3] = masked
        y[0] = masked
        assert_equal(minimum(x, y), where(less(x, y), x, y))
        assert_equal(maximum(x, y), where(greater(x, y), x, y))
        assert_(minimum(x) == 0)
        assert_(maximum(x) == 4)
        #
        x = arange(4).reshape(2, 2)
        x[-1, -1] = masked
        assert_equal(maximum(x), 2)

    def test_minimummaximum_func(self):
        a = np.ones((2, 2))
        aminimum = minimum(a, a)
        self.assertTrue(isinstance(aminimum, MaskedArray))
        assert_equal(aminimum, np.minimum(a, a))
        #
        aminimum = minimum.outer(a, a)
        self.assertTrue(isinstance(aminimum, MaskedArray))
        assert_equal(aminimum, np.minimum.outer(a, a))
        #
        amaximum = maximum(a, a)
        self.assertTrue(isinstance(amaximum, MaskedArray))
        assert_equal(amaximum, np.maximum(a, a))
        #
        amaximum = maximum.outer(a, a)
        self.assertTrue(isinstance(amaximum, MaskedArray))
        assert_equal(amaximum, np.maximum.outer(a, a))

    def test_minmax_reduce(self):
        # Test np.min/maximum.reduce on array w/ full False mask
        a = array([1, 2, 3], mask=[False, False, False])
        b = np.maximum.reduce(a)
        assert_equal(b, 3)

    def test_minmax_funcs_with_output(self):
        # Tests the min/max functions with explicit outputs
        mask = np.random.rand(12).round()
        xm = array(np.random.uniform(0, 10, 12), mask=mask)
        xm.shape = (3, 4)
        for funcname in ('min', 'max'):
            # Initialize
            npfunc = getattr(np, funcname)
            mafunc = getattr(numpy.ma.core, funcname)
            # Use the np version
            nout = np.empty((4,), dtype=int)
            try:
                result = npfunc(xm, axis=0, out=nout)
            except MaskError:
                pass
            nout = np.empty((4,), dtype=float)
            result = npfunc(xm, axis=0, out=nout)
            self.assertTrue(result is nout)
            # Use the ma version
            nout.fill(-999)
            result = mafunc(xm, axis=0, out=nout)
            self.assertTrue(result is nout)

    def test_minmax_methods(self):
        # Additional tests on max/min
        (_, _, _, _, _, xm, _, _, _, _) = self.d
        xm.shape = (xm.size,)
        assert_equal(xm.max(), 10)
        self.assertTrue(xm[0].max() is masked)
        self.assertTrue(xm[0].max(0) is masked)
        self.assertTrue(xm[0].max(-1) is masked)
        assert_equal(xm.min(), -10.)
        self.assertTrue(xm[0].min() is masked)
        self.assertTrue(xm[0].min(0) is masked)
        self.assertTrue(xm[0].min(-1) is masked)
        assert_equal(xm.ptp(), 20.)
        self.assertTrue(xm[0].ptp() is masked)
        self.assertTrue(xm[0].ptp(0) is masked)
        self.assertTrue(xm[0].ptp(-1) is masked)
        #
        x = array([1, 2, 3], mask=True)
        self.assertTrue(x.min() is masked)
        self.assertTrue(x.max() is masked)
        self.assertTrue(x.ptp() is masked)

    def test_addsumprod(self):
        # Tests add, sum, product.
        (x, y, a10, m1, m2, xm, ym, z, zm, xf) = self.d
        assert_equal(np.add.reduce(x), add.reduce(x))
        assert_equal(np.add.accumulate(x), add.accumulate(x))
        assert_equal(4, sum(array(4), axis=0))
        assert_equal(4, sum(array(4), axis=0))
        assert_equal(np.sum(x, axis=0), sum(x, axis=0))
        assert_equal(np.sum(filled(xm, 0), axis=0), sum(xm, axis=0))
        assert_equal(np.sum(x, 0), sum(x, 0))
        assert_equal(np.product(x, axis=0), product(x, axis=0))
        assert_equal(np.product(x, 0), product(x, 0))
        assert_equal(np.product(filled(xm, 1), axis=0), product(xm, axis=0))
        s = (3, 4)
        x.shape = y.shape = xm.shape = ym.shape = s
        if len(s) > 1:
            assert_equal(np.concatenate((x, y), 1), concatenate((xm, ym), 1))
            assert_equal(np.add.reduce(x, 1), add.reduce(x, 1))
            assert_equal(np.sum(x, 1), sum(x, 1))
            assert_equal(np.product(x, 1), product(x, 1))

    def test_binops_d2D(self):
        # Test binary operations on 2D data
        a = array([[1.], [2.], [3.]], mask=[[False], [True], [True]])
        b = array([[2., 3.], [4., 5.], [6., 7.]])
        #
        test = a * b
        control = array([[2., 3.], [2., 2.], [3., 3.]],
                        mask=[[0, 0], [1, 1], [1, 1]])
        assert_equal(test, control)
        assert_equal(test.data, control.data)
        assert_equal(test.mask, control.mask)
        #
        test = b * a
        control = array([[2., 3.], [4., 5.], [6., 7.]],
                        mask=[[0, 0], [1, 1], [1, 1]])
        assert_equal(test, control)
        assert_equal(test.data, control.data)
        assert_equal(test.mask, control.mask)
        #
        a = array([[1.], [2.], [3.]])
        b = array([[2., 3.], [4., 5.], [6., 7.]],
                  mask=[[0, 0], [0, 0], [0, 1]])
        test = a * b
        control = array([[2, 3], [8, 10], [18, 3]],
                        mask=[[0, 0], [0, 0], [0, 1]])
        assert_equal(test, control)
        assert_equal(test.data, control.data)
        assert_equal(test.mask, control.mask)
        #
        test = b * a
        control = array([[2, 3], [8, 10], [18, 7]],
                        mask=[[0, 0], [0, 0], [0, 1]])
        assert_equal(test, control)
        assert_equal(test.data, control.data)
        assert_equal(test.mask, control.mask)

    def test_domained_binops_d2D(self):
        # Test domained binary operations on 2D data
        a = array([[1.], [2.], [3.]], mask=[[False], [True], [True]])
        b = array([[2., 3.], [4., 5.], [6., 7.]])
        #
        test = a / b
        control = array([[1. / 2., 1. / 3.], [2., 2.], [3., 3.]],
                        mask=[[0, 0], [1, 1], [1, 1]])
        assert_equal(test, control)
        assert_equal(test.data, control.data)
        assert_equal(test.mask, control.mask)
        #
        test = b / a
        control = array([[2. / 1., 3. / 1.], [4., 5.], [6., 7.]],
                        mask=[[0, 0], [1, 1], [1, 1]])
        assert_equal(test, control)
        assert_equal(test.data, control.data)
        assert_equal(test.mask, control.mask)
        #
        a = array([[1.], [2.], [3.]])
        b = array([[2., 3.], [4., 5.], [6., 7.]],
                  mask=[[0, 0], [0, 0], [0, 1]])
        test = a / b
        control = array([[1. / 2, 1. / 3], [2. / 4, 2. / 5], [3. / 6, 3]],
                        mask=[[0, 0], [0, 0], [0, 1]])
        assert_equal(test, control)
        assert_equal(test.data, control.data)
        assert_equal(test.mask, control.mask)
        #
        test = b / a
        control = array([[2 / 1., 3 / 1.], [4 / 2., 5 / 2.], [6 / 3., 7]],
                        mask=[[0, 0], [0, 0], [0, 1]])
        assert_equal(test, control)
        assert_equal(test.data, control.data)
        assert_equal(test.mask, control.mask)

    def test_noshrinking(self):
        # Check that we don't shrink a mask when not wanted
        # Binary operations
        a = masked_array([1., 2., 3.], mask=[False, False, False],
                         shrink=False)
        b = a + 1
        assert_equal(b.mask, [0, 0, 0])
        # In place binary operation
        a += 1
        assert_equal(a.mask, [0, 0, 0])
        # Domained binary operation
        b = a / 1.
        assert_equal(b.mask, [0, 0, 0])
        # In place binary operation
        a /= 1.
        assert_equal(a.mask, [0, 0, 0])

    def test_mod(self):
        # Tests mod
        (x, y, a10, m1, m2, xm, ym, z, zm, xf) = self.d
        assert_equal(mod(x, y), mod(xm, ym))
        test = mod(ym, xm)
        assert_equal(test, np.mod(ym, xm))
        assert_equal(test.mask, mask_or(xm.mask, ym.mask))
        test = mod(xm, ym)
        assert_equal(test, np.mod(xm, ym))
        assert_equal(test.mask, mask_or(mask_or(xm.mask, ym.mask), (ym == 0)))

    def test_TakeTransposeInnerOuter(self):
        # Test of take, transpose, inner, outer products
        x = arange(24)
        y = np.arange(24)
        x[5:6] = masked
        x = x.reshape(2, 3, 4)
        y = y.reshape(2, 3, 4)
        assert_equal(np.transpose(y, (2, 0, 1)), transpose(x, (2, 0, 1)))
        assert_equal(np.take(y, (2, 0, 1), 1), take(x, (2, 0, 1), 1))
        assert_equal(np.inner(filled(x, 0), filled(y, 0)),
                     inner(x, y))
        assert_equal(np.outer(filled(x, 0), filled(y, 0)),
                     outer(x, y))
        y = array(['abc', 1, 'def', 2, 3], object)
        y[2] = masked
        t = take(y, [0, 3, 4])
        assert_(t[0] == 'abc')
        assert_(t[1] == 2)
        assert_(t[2] == 3)

    def test_imag_real(self):
        # Check complex
        xx = array([1 + 10j, 20 + 2j], mask=[1, 0])
        assert_equal(xx.imag, [10, 2])
        assert_equal(xx.imag.filled(), [1e+20, 2])
        assert_equal(xx.imag.dtype, xx._data.imag.dtype)
        assert_equal(xx.real, [1, 20])
        assert_equal(xx.real.filled(), [1e+20, 20])
        assert_equal(xx.real.dtype, xx._data.real.dtype)

    def test_methods_with_output(self):
        xm = array(np.random.uniform(0, 10, 12)).reshape(3, 4)
        xm[:, 0] = xm[0] = xm[-1, -1] = masked
        #
        funclist = ('sum', 'prod', 'var', 'std', 'max', 'min', 'ptp', 'mean',)
        #
        for funcname in funclist:
            npfunc = getattr(np, funcname)
            xmmeth = getattr(xm, funcname)
            # A ndarray as explicit input
            output = np.empty(4, dtype=float)
            output.fill(-9999)
            result = npfunc(xm, axis=0, out=output)
            # ... the result should be the given output
            assert_(result is output)
            assert_equal(result, xmmeth(axis=0, out=output))
            #
            output = empty(4, dtype=int)
            result = xmmeth(axis=0, out=output)
            assert_(result is output)
            assert_(output[0] is masked)

    def test_eq_on_structured(self):
        # Test the equality of structured arrays
        ndtype = [('A', int), ('B', int)]
        a = array([(1, 1), (2, 2)], mask=[(0, 1), (0, 0)], dtype=ndtype)
        test = (a == a)
        assert_equal(test, [True, True])
        assert_equal(test.mask, [False, False])
        b = array([(1, 1), (2, 2)], mask=[(1, 0), (0, 0)], dtype=ndtype)
        test = (a == b)
        assert_equal(test, [False, True])
        assert_equal(test.mask, [True, False])
        b = array([(1, 1), (2, 2)], mask=[(0, 1), (1, 0)], dtype=ndtype)
        test = (a == b)
        assert_equal(test, [True, False])
        assert_equal(test.mask, [False, False])

    def test_ne_on_structured(self):
        # Test the equality of structured arrays
        ndtype = [('A', int), ('B', int)]
        a = array([(1, 1), (2, 2)], mask=[(0, 1), (0, 0)], dtype=ndtype)
        test = (a != a)
        assert_equal(test, [False, False])
        assert_equal(test.mask, [False, False])
        b = array([(1, 1), (2, 2)], mask=[(1, 0), (0, 0)], dtype=ndtype)
        test = (a != b)
        assert_equal(test, [True, False])
        assert_equal(test.mask, [True, False])
        b = array([(1, 1), (2, 2)], mask=[(0, 1), (1, 0)], dtype=ndtype)
        test = (a != b)
        assert_equal(test, [False, True])
        assert_equal(test.mask, [False, False])

    def test_eq_w_None(self):
        # Really, comparisons with None should not be done, but
        # check them anyway
        # With partial mask
        a = array([1, 2], mask=[0, 1])
        assert_equal(a == None, False)
        assert_equal(a.data == None, False)
        assert_equal(a.mask == None, False)
        assert_equal(a != None, True)
        # With nomask
        a = array([1, 2], mask=False)
        assert_equal(a == None, False)
        assert_equal(a != None, True)
        # With complete mask
        a = array([1, 2], mask=True)
        assert_equal(a == None, False)
        assert_equal(a != None, True)
        # Fully masked, even comparison to None should return "masked"
        a = masked
        assert_equal(a == None, masked)

    def test_eq_w_scalar(self):
        a = array(1)
        assert_equal(a == 1, True)
        assert_equal(a == 0, False)
        assert_equal(a != 1, False)
        assert_equal(a != 0, True)

    def test_numpyarithmetics(self):
        # Check that the mask is not back-propagated when using numpy functions
        a = masked_array([-1, 0, 1, 2, 3], mask=[0, 0, 0, 0, 1])
        control = masked_array([np.nan, np.nan, 0, np.log(2), -1],
                               mask=[1, 1, 0, 0, 1])
        #
        test = log(a)
        assert_equal(test, control)
        assert_equal(test.mask, control.mask)
        assert_equal(a.mask, [0, 0, 0, 0, 1])
        #
        test = np.log(a)
        assert_equal(test, control)
        assert_equal(test.mask, control.mask)
        assert_equal(a.mask, [0, 0, 0, 0, 1])


#------------------------------------------------------------------------------
class TestMaskedArrayAttributes(TestCase):

    def test_keepmask(self):
        # Tests the keep mask flag
        x = masked_array([1, 2, 3], mask=[1, 0, 0])
        mx = masked_array(x)
        assert_equal(mx.mask, x.mask)
        mx = masked_array(x, mask=[0, 1, 0], keep_mask=False)
        assert_equal(mx.mask, [0, 1, 0])
        mx = masked_array(x, mask=[0, 1, 0], keep_mask=True)
        assert_equal(mx.mask, [1, 1, 0])
        # We default to true
        mx = masked_array(x, mask=[0, 1, 0])
        assert_equal(mx.mask, [1, 1, 0])

    def test_hardmask(self):
        # Test hard_mask
        d = arange(5)
        n = [0, 0, 0, 1, 1]
        m = make_mask(n)
        xh = array(d, mask=m, hard_mask=True)
        # We need to copy, to avoid updating d in xh !
        xs = array(d, mask=m, hard_mask=False, copy=True)
        xh[[1, 4]] = [10, 40]
        xs[[1, 4]] = [10, 40]
        assert_equal(xh._data, [0, 10, 2, 3, 4])
        assert_equal(xs._data, [0, 10, 2, 3, 40])
        #assert_equal(xh.mask.ctypes._data, m.ctypes._data)
        assert_equal(xs.mask, [0, 0, 0, 1, 0])
        self.assertTrue(xh._hardmask)
        self.assertTrue(not xs._hardmask)
        xh[1:4] = [10, 20, 30]
        xs[1:4] = [10, 20, 30]
        assert_equal(xh._data, [0, 10, 20, 3, 4])
        assert_equal(xs._data, [0, 10, 20, 30, 40])
        #assert_equal(xh.mask.ctypes._data, m.ctypes._data)
        assert_equal(xs.mask, nomask)
        xh[0] = masked
        xs[0] = masked
        assert_equal(xh.mask, [1, 0, 0, 1, 1])
        assert_equal(xs.mask, [1, 0, 0, 0, 0])
        xh[:] = 1
        xs[:] = 1
        assert_equal(xh._data, [0, 1, 1, 3, 4])
        assert_equal(xs._data, [1, 1, 1, 1, 1])
        assert_equal(xh.mask, [1, 0, 0, 1, 1])
        assert_equal(xs.mask, nomask)
        # Switch to soft mask
        xh.soften_mask()
        xh[:] = arange(5)
        assert_equal(xh._data, [0, 1, 2, 3, 4])
        assert_equal(xh.mask, nomask)
        # Switch back to hard mask
        xh.harden_mask()
        xh[xh < 3] = masked
        assert_equal(xh._data, [0, 1, 2, 3, 4])
        assert_equal(xh._mask, [1, 1, 1, 0, 0])
        xh[filled(xh > 1, False)] = 5
        assert_equal(xh._data, [0, 1, 2, 5, 5])
        assert_equal(xh._mask, [1, 1, 1, 0, 0])
        #
        xh = array([[1, 2], [3, 4]], mask=[[1, 0], [0, 0]], hard_mask=True)
        xh[0] = 0
        assert_equal(xh._data, [[1, 0], [3, 4]])
        assert_equal(xh._mask, [[1, 0], [0, 0]])
        xh[-1, -1] = 5
        assert_equal(xh._data, [[1, 0], [3, 5]])
        assert_equal(xh._mask, [[1, 0], [0, 0]])
        xh[filled(xh < 5, False)] = 2
        assert_equal(xh._data, [[1, 2], [2, 5]])
        assert_equal(xh._mask, [[1, 0], [0, 0]])

    def test_hardmask_again(self):
        # Another test of hardmask
        d = arange(5)
        n = [0, 0, 0, 1, 1]
        m = make_mask(n)
        xh = array(d, mask=m, hard_mask=True)
        xh[4:5] = 999
        #assert_equal(xh.mask.ctypes._data, m.ctypes._data)
        xh[0:1] = 999
        assert_equal(xh._data, [999, 1, 2, 3, 4])

    def test_hardmask_oncemore_yay(self):
        # OK, yet another test of hardmask
        # Make sure that harden_mask/soften_mask//unshare_mask returns self
        a = array([1, 2, 3], mask=[1, 0, 0])
        b = a.harden_mask()
        assert_equal(a, b)
        b[0] = 0
        assert_equal(a, b)
        assert_equal(b, array([1, 2, 3], mask=[1, 0, 0]))
        a = b.soften_mask()
        a[0] = 0
        assert_equal(a, b)
        assert_equal(b, array([0, 2, 3], mask=[0, 0, 0]))

    def test_smallmask(self):
        # Checks the behaviour of _smallmask
        a = arange(10)
        a[1] = masked
        a[1] = 1
        assert_equal(a._mask, nomask)
        a = arange(10)
        a._smallmask = False
        a[1] = masked
        a[1] = 1
        assert_equal(a._mask, zeros(10))

    def test_shrink_mask(self):
        # Tests .shrink_mask()
        a = array([1, 2, 3], mask=[0, 0, 0])
        b = a.shrink_mask()
        assert_equal(a, b)
        assert_equal(a.mask, nomask)

    def test_flat(self):
        # Test that flat can return all types of items [#4585, #4615]
        # test simple access
        test = masked_array(np.matrix([[1, 2, 3]]), mask=[0, 0, 1])
        assert_equal(test.flat[1], 2)
        assert_equal(test.flat[2], masked)
        self.assertTrue(np.all(test.flat[0:2] == test[0, 0:2]))
        # Test flat on masked_matrices
        test = masked_array(np.matrix([[1, 2, 3]]), mask=[0, 0, 1])
        test.flat = masked_array([3, 2, 1], mask=[1, 0, 0])
        control = masked_array(np.matrix([[3, 2, 1]]), mask=[1, 0, 0])
        assert_equal(test, control)
        # Test setting
        test = masked_array(np.matrix([[1, 2, 3]]), mask=[0, 0, 1])
        testflat = test.flat
        testflat[:] = testflat[[2, 1, 0]]
        assert_equal(test, control)
        testflat[0] = 9
        assert_equal(test[0, 0], 9)
        # test 2-D record array
        # ... on structured array w/ masked records
        x = array([[(1, 1.1, 'one'), (2, 2.2, 'two'), (3, 3.3, 'thr')],
                   [(4, 4.4, 'fou'), (5, 5.5, 'fiv'), (6, 6.6, 'six')]],
                  dtype=[('a', int), ('b', float), ('c', '|S8')])
        x['a'][0, 1] = masked
        x['b'][1, 0] = masked
        x['c'][0, 2] = masked
        x[-1, -1] = masked
        xflat = x.flat
        assert_equal(xflat[0], x[0, 0])
        assert_equal(xflat[1], x[0, 1])
        assert_equal(xflat[2], x[0, 2])
        assert_equal(xflat[:3], x[0])
        assert_equal(xflat[3], x[1, 0])
        assert_equal(xflat[4], x[1, 1])
        assert_equal(xflat[5], x[1, 2])
        assert_equal(xflat[3:], x[1])
        assert_equal(xflat[-1], x[-1, -1])
        i = 0
        j = 0
        for xf in xflat:
            assert_equal(xf, x[j, i])
            i += 1
            if i >= x.shape[-1]:
                i = 0
                j += 1
        # test that matrices keep the correct shape (#4615)
        a = masked_array(np.matrix(np.eye(2)), mask=0)
        b = a.flat
        b01 = b[:2]
        assert_equal(b01.data, array([[1., 0.]]))
        assert_equal(b01.mask, array([[False, False]]))


#------------------------------------------------------------------------------
class TestFillingValues(TestCase):

    def test_check_on_scalar(self):
        # Test _check_fill_value set to valid and invalid values
        _check_fill_value = np.ma.core._check_fill_value
        #
        fval = _check_fill_value(0, int)
        assert_equal(fval, 0)
        fval = _check_fill_value(None, int)
        assert_equal(fval, default_fill_value(0))
        #
        fval = _check_fill_value(0, "|S3")
        assert_equal(fval, asbytes("0"))
        fval = _check_fill_value(None, "|S3")
        assert_equal(fval, default_fill_value("|S3"))
        self.assertRaises(TypeError, _check_fill_value, 1e+20, int)
        self.assertRaises(TypeError, _check_fill_value, 'stuff', int)

    def test_check_on_fields(self):
        # Tests _check_fill_value with records
        _check_fill_value = np.ma.core._check_fill_value
        ndtype = [('a', int), ('b', float), ('c', "|S3")]
        # A check on a list should return a single record
        fval = _check_fill_value([-999, -12345678.9, "???"], ndtype)
        self.assertTrue(isinstance(fval, ndarray))
        assert_equal(fval.item(), [-999, -12345678.9, asbytes("???")])
        # A check on None should output the defaults
        fval = _check_fill_value(None, ndtype)
        self.assertTrue(isinstance(fval, ndarray))
        assert_equal(fval.item(), [default_fill_value(0),
                                   default_fill_value(0.),
                                   asbytes(default_fill_value("0"))])
        #.....Using a structured type as fill_value should work
        fill_val = np.array((-999, -12345678.9, "???"), dtype=ndtype)
        fval = _check_fill_value(fill_val, ndtype)
        self.assertTrue(isinstance(fval, ndarray))
        assert_equal(fval.item(), [-999, -12345678.9, asbytes("???")])

        #.....Using a flexible type w/ a different type shouldn't matter
        # BEHAVIOR in 1.5 and earlier: match structured types by position
        #fill_val = np.array((-999, -12345678.9, "???"),
        #                    dtype=[("A", int), ("B", float), ("C", "|S3")])
        # BEHAVIOR in 1.6 and later: match structured types by name
        fill_val = np.array(("???", -999, -12345678.9),
                            dtype=[("c", "|S3"), ("a", int), ("b", float), ])
        fval = _check_fill_value(fill_val, ndtype)
        self.assertTrue(isinstance(fval, ndarray))
        assert_equal(fval.item(), [-999, -12345678.9, asbytes("???")])

        #.....Using an object-array shouldn't matter either
        fill_val = np.ndarray(shape=(1,), dtype=object)
        fill_val[0] = (-999, -12345678.9, asbytes("???"))
        fval = _check_fill_value(fill_val, object)
        self.assertTrue(isinstance(fval, ndarray))
        assert_equal(fval.item(), [-999, -12345678.9, asbytes("???")])
        # NOTE: This test was never run properly as "fill_value" rather than
        # "fill_val" was assigned.  Written properly, it fails.
        #fill_val = np.array((-999, -12345678.9, "???"))
        #fval = _check_fill_value(fill_val, ndtype)
        #self.assertTrue(isinstance(fval, ndarray))
        #assert_equal(fval.item(), [-999, -12345678.9, asbytes("???")])
        #.....One-field-only flexible type should work as well
        ndtype = [("a", int)]
        fval = _check_fill_value(-999999999, ndtype)
        self.assertTrue(isinstance(fval, ndarray))
        assert_equal(fval.item(), (-999999999,))

    def test_fillvalue_conversion(self):
        # Tests the behavior of fill_value during conversion
        # We had a tailored comment to make sure special attributes are
        # properly dealt with
        a = array(asbytes_nested(['3', '4', '5']))
        a._optinfo.update({'comment':"updated!"})
        #
        b = array(a, dtype=int)
        assert_equal(b._data, [3, 4, 5])
        assert_equal(b.fill_value, default_fill_value(0))
        #
        b = array(a, dtype=float)
        assert_equal(b._data, [3, 4, 5])
        assert_equal(b.fill_value, default_fill_value(0.))
        #
        b = a.astype(int)
        assert_equal(b._data, [3, 4, 5])
        assert_equal(b.fill_value, default_fill_value(0))
        assert_equal(b._optinfo['comment'], "updated!")
        #
        b = a.astype([('a', '|S3')])
        assert_equal(b['a']._data, a._data)
        assert_equal(b['a'].fill_value, a.fill_value)

    def test_fillvalue(self):
        # Yet more fun with the fill_value
        data = masked_array([1, 2, 3], fill_value=-999)
        series = data[[0, 2, 1]]
        assert_equal(series._fill_value, data._fill_value)
        #
        mtype = [('f', float), ('s', '|S3')]
        x = array([(1, 'a'), (2, 'b'), (pi, 'pi')], dtype=mtype)
        x.fill_value = 999
        assert_equal(x.fill_value.item(), [999., asbytes('999')])
        assert_equal(x['f'].fill_value, 999)
        assert_equal(x['s'].fill_value, asbytes('999'))
        #
        x.fill_value = (9, '???')
        assert_equal(x.fill_value.item(), (9, asbytes('???')))
        assert_equal(x['f'].fill_value, 9)
        assert_equal(x['s'].fill_value, asbytes('???'))
        #
        x = array([1, 2, 3.1])
        x.fill_value = 999
        assert_equal(np.asarray(x.fill_value).dtype, float)
        assert_equal(x.fill_value, 999.)
        assert_equal(x._fill_value, np.array(999.))

    def test_fillvalue_exotic_dtype(self):
        # Tests yet more exotic flexible dtypes
        _check_fill_value = np.ma.core._check_fill_value
        ndtype = [('i', int), ('s', '|S8'), ('f', float)]
        control = np.array((default_fill_value(0),
                            default_fill_value('0'),
                            default_fill_value(0.),),
                           dtype=ndtype)
        assert_equal(_check_fill_value(None, ndtype), control)
        # The shape shouldn't matter
        ndtype = [('f0', float, (2, 2))]
        control = np.array((default_fill_value(0.),),
                           dtype=[('f0', float)]).astype(ndtype)
        assert_equal(_check_fill_value(None, ndtype), control)
        control = np.array((0,), dtype=[('f0', float)]).astype(ndtype)
        assert_equal(_check_fill_value(0, ndtype), control)
        #
        ndtype = np.dtype("int, (2,3)float, float")
        control = np.array((default_fill_value(0),
                            default_fill_value(0.),
                            default_fill_value(0.),),
                           dtype="int, float, float").astype(ndtype)
        test = _check_fill_value(None, ndtype)
        assert_equal(test, control)
        control = np.array((0, 0, 0), dtype="int, float, float").astype(ndtype)
        assert_equal(_check_fill_value(0, ndtype), control)

    def test_extremum_fill_value(self):
        # Tests extremum fill values for flexible type.
        a = array([(1, (2, 3)), (4, (5, 6))],
                  dtype=[('A', int), ('B', [('BA', int), ('BB', int)])])
        test = a.fill_value
        assert_equal(test['A'], default_fill_value(a['A']))
        assert_equal(test['B']['BA'], default_fill_value(a['B']['BA']))
        assert_equal(test['B']['BB'], default_fill_value(a['B']['BB']))
        #
        test = minimum_fill_value(a)
        assert_equal(test[0], minimum_fill_value(a['A']))
        assert_equal(test[1][0], minimum_fill_value(a['B']['BA']))
        assert_equal(test[1][1], minimum_fill_value(a['B']['BB']))
        assert_equal(test[1], minimum_fill_value(a['B']))
        #
        test = maximum_fill_value(a)
        assert_equal(test[0], maximum_fill_value(a['A']))
        assert_equal(test[1][0], maximum_fill_value(a['B']['BA']))
        assert_equal(test[1][1], maximum_fill_value(a['B']['BB']))
        assert_equal(test[1], maximum_fill_value(a['B']))

    def test_fillvalue_individual_fields(self):
        # Test setting fill_value on individual fields
        ndtype = [('a', int), ('b', int)]
        # Explicit fill_value
        a = array(list(zip([1, 2, 3], [4, 5, 6])),
                  fill_value=(-999, -999), dtype=ndtype)
        aa = a['a']
        aa.set_fill_value(10)
        assert_equal(aa._fill_value, np.array(10))
        assert_equal(tuple(a.fill_value), (10, -999))
        a.fill_value['b'] = -10
        assert_equal(tuple(a.fill_value), (10, -10))
        # Implicit fill_value
        t = array(list(zip([1, 2, 3], [4, 5, 6])), dtype=ndtype)
        tt = t['a']
        tt.set_fill_value(10)
        assert_equal(tt._fill_value, np.array(10))
        assert_equal(tuple(t.fill_value), (10, default_fill_value(0)))

    def test_fillvalue_implicit_structured_array(self):
        # Check that fill_value is always defined for structured arrays
        ndtype = ('b', float)
        adtype = ('a', float)
        a = array([(1.,), (2.,)], mask=[(False,), (False,)],
                  fill_value=(np.nan,), dtype=np.dtype([adtype]))
        b = empty(a.shape, dtype=[adtype, ndtype])
        b['a'] = a['a']
        b['a'].set_fill_value(a['a'].fill_value)
        f = b._fill_value[()]
        assert_(np.isnan(f[0]))
        assert_equal(f[-1], default_fill_value(1.))

    def test_fillvalue_as_arguments(self):
        # Test adding a fill_value parameter to empty/ones/zeros
        a = empty(3, fill_value=999.)
        assert_equal(a.fill_value, 999.)
        #
        a = ones(3, fill_value=999., dtype=float)
        assert_equal(a.fill_value, 999.)
        #
        a = zeros(3, fill_value=0., dtype=complex)
        assert_equal(a.fill_value, 0.)
        #
        a = identity(3, fill_value=0., dtype=complex)
        assert_equal(a.fill_value, 0.)

    def test_fillvalue_in_view(self):
        # Test the behavior of fill_value in view

        # Create initial masked array
        x = array([1, 2, 3], fill_value=1, dtype=np.int64)

        # Check that fill_value is preserved by default
        y = x.view()
        assert_(y.fill_value == 1)

        # Check that fill_value is preserved if dtype is specified and the
        # dtype is an ndarray sub-class and has a _fill_value attribute
        y = x.view(MaskedArray)
        assert_(y.fill_value == 1)

        # Check that fill_value is preserved if type is specified and the
        # dtype is an ndarray sub-class and has a _fill_value attribute (by
        # default, the first argument is dtype, not type)
        y = x.view(type=MaskedArray)
        assert_(y.fill_value == 1)

        # Check that code does not crash if passed an ndarray sub-class that
        # does not have a _fill_value attribute
        y = x.view(np.ndarray)
        y = x.view(type=np.ndarray)

        # Check that fill_value can be overriden with view
        y = x.view(MaskedArray, fill_value=2)
        assert_(y.fill_value == 2)

        # Check that fill_value can be overriden with view (using type=)
        y = x.view(type=MaskedArray, fill_value=2)
        assert_(y.fill_value == 2)

        # Check that fill_value gets reset if passed a dtype but not a
        # fill_value. This is because even though in some cases one can safely
        # cast the fill_value, e.g. if taking an int64 view of an int32 array,
        # in other cases, this cannot be done (e.g. int32 view of an int64
        # array with a large fill_value).
        y = x.view(dtype=np.int32)
        assert_(y.fill_value == 999999)


#------------------------------------------------------------------------------
class TestUfuncs(TestCase):
    # Test class for the application of ufuncs on MaskedArrays.

    def setUp(self):
        # Base data definition.
        self.d = (array([1.0, 0, -1, pi / 2] * 2, mask=[0, 1] + [0] * 6),
                  array([1.0, 0, -1, pi / 2] * 2, mask=[1, 0] + [0] * 6),)
        self.err_status = np.geterr()
        np.seterr(divide='ignore', invalid='ignore')

    def tearDown(self):
        np.seterr(**self.err_status)

    def test_testUfuncRegression(self):
        # Tests new ufuncs on MaskedArrays.
        for f in ['sqrt', 'log', 'log10', 'exp', 'conjugate',
                  'sin', 'cos', 'tan',
                  'arcsin', 'arccos', 'arctan',
                  'sinh', 'cosh', 'tanh',
                  'arcsinh',
                  'arccosh',
                  'arctanh',
                  'absolute', 'fabs', 'negative',
                  # 'nonzero', 'around',
                  'floor', 'ceil',
                  # 'sometrue', 'alltrue',
                  'logical_not',
                  'add', 'subtract', 'multiply',
                  'divide', 'true_divide', 'floor_divide',
                  'remainder', 'fmod', 'hypot', 'arctan2',
                  'equal', 'not_equal', 'less_equal', 'greater_equal',
                  'less', 'greater',
                  'logical_and', 'logical_or', 'logical_xor',
                  ]:
            try:
                uf = getattr(umath, f)
            except AttributeError:
                uf = getattr(fromnumeric, f)
            mf = getattr(numpy.ma.core, f)
            args = self.d[:uf.nin]
            ur = uf(*args)
            mr = mf(*args)
            assert_equal(ur.filled(0), mr.filled(0), f)
            assert_mask_equal(ur.mask, mr.mask, err_msg=f)

    def test_reduce(self):
        # Tests reduce on MaskedArrays.
        a = self.d[0]
        self.assertTrue(not alltrue(a, axis=0))
        self.assertTrue(sometrue(a, axis=0))
        assert_equal(sum(a[:3], axis=0), 0)
        assert_equal(product(a, axis=0), 0)
        assert_equal(add.reduce(a), pi)

    def test_minmax(self):
        # Tests extrema on MaskedArrays.
        a = arange(1, 13).reshape(3, 4)
        amask = masked_where(a < 5, a)
        assert_equal(amask.max(), a.max())
        assert_equal(amask.min(), 5)
        assert_equal(amask.max(0), a.max(0))
        assert_equal(amask.min(0), [5, 6, 7, 8])
        self.assertTrue(amask.max(1)[0].mask)
        self.assertTrue(amask.min(1)[0].mask)

    def test_ndarray_mask(self):
        # Check that the mask of the result is a ndarray (not a MaskedArray...)
        a = masked_array([-1, 0, 1, 2, 3], mask=[0, 0, 0, 0, 1])
        test = np.sqrt(a)
        control = masked_array([-1, 0, 1, np.sqrt(2), -1],
                               mask=[1, 0, 0, 0, 1])
        assert_equal(test, control)
        assert_equal(test.mask, control.mask)
        self.assertTrue(not isinstance(test.mask, MaskedArray))

    def test_treatment_of_NotImplemented(self):
        # Check any NotImplemented returned by umath.<ufunc> is passed on
        a = masked_array([1., 2.], mask=[1, 0])
        # basic tests for _MaskedBinaryOperation
        assert_(a.__mul__('abc') is NotImplemented)
        assert_(multiply.outer(a, 'abc') is NotImplemented)
        # and for _DomainedBinaryOperation
        assert_(a.__div__('abc') is NotImplemented)

        # also check explicitly that rmul of another class can be accessed
        class MyClass(str):
            def __mul__(self, other):
                return "My mul"

            def __rmul__(self, other):
                return "My rmul"

        me = MyClass()
        assert_(me * a == "My mul")
        assert_(a * me == "My rmul")


#------------------------------------------------------------------------------
class TestMaskedArrayInPlaceArithmetics(TestCase):
    # Test MaskedArray Arithmetics

    def setUp(self):
        x = arange(10)
        y = arange(10)
        xm = arange(10)
        xm[2] = masked
        self.intdata = (x, y, xm)
        self.floatdata = (x.astype(float), y.astype(float), xm.astype(float))

    def test_inplace_addition_scalar(self):
        # Test of inplace additions
        (x, y, xm) = self.intdata
        xm[2] = masked
        x += 1
        assert_equal(x, y + 1)
        xm += 1
        assert_equal(xm, y + 1)
        #
        (x, _, xm) = self.floatdata
        id1 = x.data.ctypes._data
        x += 1.
        assert_(id1 == x.data.ctypes._data)
        assert_equal(x, y + 1.)

    def test_inplace_addition_array(self):
        # Test of inplace additions
        (x, y, xm) = self.intdata
        m = xm.mask
        a = arange(10, dtype=np.int16)
        a[-1] = masked
        x += a
        xm += a
        assert_equal(x, y + a)
        assert_equal(xm, y + a)
        assert_equal(xm.mask, mask_or(m, a.mask))

    def test_inplace_subtraction_scalar(self):
        # Test of inplace subtractions
        (x, y, xm) = self.intdata
        x -= 1
        assert_equal(x, y - 1)
        xm -= 1
        assert_equal(xm, y - 1)

    def test_inplace_subtraction_array(self):
        # Test of inplace subtractions
        (x, y, xm) = self.floatdata
        m = xm.mask
        a = arange(10, dtype=float)
        a[-1] = masked
        x -= a
        xm -= a
        assert_equal(x, y - a)
        assert_equal(xm, y - a)
        assert_equal(xm.mask, mask_or(m, a.mask))

    def test_inplace_multiplication_scalar(self):
        # Test of inplace multiplication
        (x, y, xm) = self.floatdata
        x *= 2.0
        assert_equal(x, y * 2)
        xm *= 2.0
        assert_equal(xm, y * 2)

    def test_inplace_multiplication_array(self):
        # Test of inplace multiplication
        (x, y, xm) = self.floatdata
        m = xm.mask
        a = arange(10, dtype=float)
        a[-1] = masked
        x *= a
        xm *= a
        assert_equal(x, y * a)
        assert_equal(xm, y * a)
        assert_equal(xm.mask, mask_or(m, a.mask))

    def test_inplace_division_scalar_int(self):
        # Test of inplace division
        (x, y, xm) = self.intdata
        x = arange(10) * 2
        xm = arange(10) * 2
        xm[2] = masked
        x //= 2
        assert_equal(x, y)
        xm //= 2
        assert_equal(xm, y)

    def test_inplace_division_scalar_float(self):
        # Test of inplace division
        (x, y, xm) = self.floatdata
        x /= 2.0
        assert_equal(x, y / 2.0)
        xm /= arange(10)
        assert_equal(xm, ones((10,)))

    def test_inplace_division_array_float(self):
        # Test of inplace division
        (x, y, xm) = self.floatdata
        m = xm.mask
        a = arange(10, dtype=float)
        a[-1] = masked
        x /= a
        xm /= a
        assert_equal(x, y / a)
        assert_equal(xm, y / a)
        assert_equal(xm.mask, mask_or(mask_or(m, a.mask), (a == 0)))

    def test_inplace_division_misc(self):
        #
        x = [1., 1., 1., -2., pi / 2., 4., 5., -10., 10., 1., 2., 3.]
        y = [5., 0., 3., 2., -1., -4., 0., -10., 10., 1., 0., 3.]
        m1 = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]
        m2 = [0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1]
        xm = masked_array(x, mask=m1)
        ym = masked_array(y, mask=m2)
        #
        z = xm / ym
        assert_equal(z._mask, [1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1])
        assert_equal(z._data,
                     [1., 1., 1., -1., -pi / 2., 4., 5., 1., 1., 1., 2., 3.])
        #assert_equal(z._data, [0.2,1.,1./3.,-1.,-pi/2.,-1.,5.,1.,1.,1.,2.,1.])
        #
        xm = xm.copy()
        xm /= ym
        assert_equal(xm._mask, [1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1])
        assert_equal(z._data,
                     [1., 1., 1., -1., -pi / 2., 4., 5., 1., 1., 1., 2., 3.])
        #assert_equal(xm._data,
        #             [1/5.,1.,1./3.,-1.,-pi/2.,-1.,5.,1.,1.,1.,2.,1.])

    def test_datafriendly_add(self):
        # Test keeping data w/ (inplace) addition
        x = array([1, 2, 3], mask=[0, 0, 1])
        # Test add w/ scalar
        xx = x + 1
        assert_equal(xx.data, [2, 3, 3])
        assert_equal(xx.mask, [0, 0, 1])
        # Test iadd w/ scalar
        x += 1
        assert_equal(x.data, [2, 3, 3])
        assert_equal(x.mask, [0, 0, 1])
        # Test add w/ array
        x = array([1, 2, 3], mask=[0, 0, 1])
        xx = x + array([1, 2, 3], mask=[1, 0, 0])
        assert_equal(xx.data, [1, 4, 3])
        assert_equal(xx.mask, [1, 0, 1])
        # Test iadd w/ array
        x = array([1, 2, 3], mask=[0, 0, 1])
        x += array([1, 2, 3], mask=[1, 0, 0])
        assert_equal(x.data, [1, 4, 3])
        assert_equal(x.mask, [1, 0, 1])

    def test_datafriendly_sub(self):
        # Test keeping data w/ (inplace) subtraction
        # Test sub w/ scalar
        x = array([1, 2, 3], mask=[0, 0, 1])
        xx = x - 1
        assert_equal(xx.data, [0, 1, 3])
        assert_equal(xx.mask, [0, 0, 1])
        # Test isub w/ scalar
        x = array([1, 2, 3], mask=[0, 0, 1])
        x -= 1
        assert_equal(x.data, [0, 1, 3])
        assert_equal(x.mask, [0, 0, 1])
        # Test sub w/ array
        x = array([1, 2, 3], mask=[0, 0, 1])
        xx = x - array([1, 2, 3], mask=[1, 0, 0])
        assert_equal(xx.data, [1, 0, 3])
        assert_equal(xx.mask, [1, 0, 1])
        # Test isub w/ array
        x = array([1, 2, 3], mask=[0, 0, 1])
        x -= array([1, 2, 3], mask=[1, 0, 0])
        assert_equal(x.data, [1, 0, 3])
        assert_equal(x.mask, [1, 0, 1])

    def test_datafriendly_mul(self):
        # Test keeping data w/ (inplace) multiplication
        # Test mul w/ scalar
        x = array([1, 2, 3], mask=[0, 0, 1])
        xx = x * 2
        assert_equal(xx.data, [2, 4, 3])
        assert_equal(xx.mask, [0, 0, 1])
        # Test imul w/ scalar
        x = array([1, 2, 3], mask=[0, 0, 1])
        x *= 2
        assert_equal(x.data, [2, 4, 3])
        assert_equal(x.mask, [0, 0, 1])
        # Test mul w/ array
        x = array([1, 2, 3], mask=[0, 0, 1])
        xx = x * array([10, 20, 30], mask=[1, 0, 0])
        assert_equal(xx.data, [1, 40, 3])
        assert_equal(xx.mask, [1, 0, 1])
        # Test imul w/ array
        x = array([1, 2, 3], mask=[0, 0, 1])
        x *= array([10, 20, 30], mask=[1, 0, 0])
        assert_equal(x.data, [1, 40, 3])
        assert_equal(x.mask, [1, 0, 1])

    def test_datafriendly_div(self):
        # Test keeping data w/ (inplace) division
        # Test div on scalar
        x = array([1, 2, 3], mask=[0, 0, 1])
        xx = x / 2.
        assert_equal(xx.data, [1 / 2., 2 / 2., 3])
        assert_equal(xx.mask, [0, 0, 1])
        # Test idiv on scalar
        x = array([1., 2., 3.], mask=[0, 0, 1])
        x /= 2.
        assert_equal(x.data, [1 / 2., 2 / 2., 3])
        assert_equal(x.mask, [0, 0, 1])
        # Test div on array
        x = array([1., 2., 3.], mask=[0, 0, 1])
        xx = x / array([10., 20., 30.], mask=[1, 0, 0])
        assert_equal(xx.data, [1., 2. / 20., 3.])
        assert_equal(xx.mask, [1, 0, 1])
        # Test idiv on array
        x = array([1., 2., 3.], mask=[0, 0, 1])
        x /= array([10., 20., 30.], mask=[1, 0, 0])
        assert_equal(x.data, [1., 2 / 20., 3.])
        assert_equal(x.mask, [1, 0, 1])

    def test_datafriendly_pow(self):
        # Test keeping data w/ (inplace) power
        # Test pow on scalar
        x = array([1., 2., 3.], mask=[0, 0, 1])
        xx = x ** 2.5
        assert_equal(xx.data, [1., 2. ** 2.5, 3.])
        assert_equal(xx.mask, [0, 0, 1])
        # Test ipow on scalar
        x **= 2.5
        assert_equal(x.data, [1., 2. ** 2.5, 3])
        assert_equal(x.mask, [0, 0, 1])

    def test_datafriendly_add_arrays(self):
        a = array([[1, 1], [3, 3]])
        b = array([1, 1], mask=[0, 0])
        a += b
        assert_equal(a, [[2, 2], [4, 4]])
        if a.mask is not nomask:
            assert_equal(a.mask, [[0, 0], [0, 0]])
        #
        a = array([[1, 1], [3, 3]])
        b = array([1, 1], mask=[0, 1])
        a += b
        assert_equal(a, [[2, 2], [4, 4]])
        assert_equal(a.mask, [[0, 1], [0, 1]])

    def test_datafriendly_sub_arrays(self):
        a = array([[1, 1], [3, 3]])
        b = array([1, 1], mask=[0, 0])
        a -= b
        assert_equal(a, [[0, 0], [2, 2]])
        if a.mask is not nomask:
            assert_equal(a.mask, [[0, 0], [0, 0]])
        #
        a = array([[1, 1], [3, 3]])
        b = array([1, 1], mask=[0, 1])
        a -= b
        assert_equal(a, [[0, 0], [2, 2]])
        assert_equal(a.mask, [[0, 1], [0, 1]])

    def test_datafriendly_mul_arrays(self):
        a = array([[1, 1], [3, 3]])
        b = array([1, 1], mask=[0, 0])
        a *= b
        assert_equal(a, [[1, 1], [3, 3]])
        if a.mask is not nomask:
            assert_equal(a.mask, [[0, 0], [0, 0]])
        #
        a = array([[1, 1], [3, 3]])
        b = array([1, 1], mask=[0, 1])
        a *= b
        assert_equal(a, [[1, 1], [3, 3]])
        assert_equal(a.mask, [[0, 1], [0, 1]])


#------------------------------------------------------------------------------
class TestMaskedArrayMethods(TestCase):
    # Test class for miscellaneous MaskedArrays methods.
    def setUp(self):
        # Base data definition.
        x = np.array([8.375, 7.545, 8.828, 8.5, 1.757, 5.928,
                      8.43, 7.78, 9.865, 5.878, 8.979, 4.732,
                      3.012, 6.022, 5.095, 3.116, 5.238, 3.957,
                      6.04, 9.63, 7.712, 3.382, 4.489, 6.479,
                      7.189, 9.645, 5.395, 4.961, 9.894, 2.893,
                      7.357, 9.828, 6.272, 3.758, 6.693, 0.993])
        X = x.reshape(6, 6)
        XX = x.reshape(3, 2, 2, 3)

        m = np.array([0, 1, 0, 1, 0, 0,
                     1, 0, 1, 1, 0, 1,
                     0, 0, 0, 1, 0, 1,
                     0, 0, 0, 1, 1, 1,
                     1, 0, 0, 1, 0, 0,
                     0, 0, 1, 0, 1, 0])
        mx = array(data=x, mask=m)
        mX = array(data=X, mask=m.reshape(X.shape))
        mXX = array(data=XX, mask=m.reshape(XX.shape))

        m2 = np.array([1, 1, 0, 1, 0, 0,
                      1, 1, 1, 1, 0, 1,
                      0, 0, 1, 1, 0, 1,
                      0, 0, 0, 1, 1, 1,
                      1, 0, 0, 1, 1, 0,
                      0, 0, 1, 0, 1, 1])
        m2x = array(data=x, mask=m2)
        m2X = array(data=X, mask=m2.reshape(X.shape))
        m2XX = array(data=XX, mask=m2.reshape(XX.shape))
        self.d = (x, X, XX, m, mx, mX, mXX, m2x, m2X, m2XX)

    def test_generic_methods(self):
        # Tests some MaskedArray methods.
        a = array([1, 3, 2])
        assert_equal(a.any(), a._data.any())
        assert_equal(a.all(), a._data.all())
        assert_equal(a.argmax(), a._data.argmax())
        assert_equal(a.argmin(), a._data.argmin())
        assert_equal(a.choose(0, 1, 2, 3, 4), a._data.choose(0, 1, 2, 3, 4))
        assert_equal(a.compress([1, 0, 1]), a._data.compress([1, 0, 1]))
        assert_equal(a.conj(), a._data.conj())
        assert_equal(a.conjugate(), a._data.conjugate())
        #
        m = array([[1, 2], [3, 4]])
        assert_equal(m.diagonal(), m._data.diagonal())
        assert_equal(a.sum(), a._data.sum())
        assert_equal(a.take([1, 2]), a._data.take([1, 2]))
        assert_equal(m.transpose(), m._data.transpose())

    def test_allclose(self):
        # Tests allclose on arrays
        a = np.random.rand(10)
        b = a + np.random.rand(10) * 1e-8
        self.assertTrue(allclose(a, b))
        # Test allclose w/ infs
        a[0] = np.inf
        self.assertTrue(not allclose(a, b))
        b[0] = np.inf
        self.assertTrue(allclose(a, b))
        # Test all close w/ masked
        a = masked_array(a)
        a[-1] = masked
        self.assertTrue(allclose(a, b, masked_equal=True))
        self.assertTrue(not allclose(a, b, masked_equal=False))
        # Test comparison w/ scalar
        a *= 1e-8
        a[0] = 0
        self.assertTrue(allclose(a, 0, masked_equal=True))

        # Test that the function works for MIN_INT integer typed arrays
        a = masked_array([np.iinfo(np.int_).min], dtype=np.int_)
        self.assertTrue(allclose(a, a))

    def test_allany(self):
        # Checks the any/all methods/functions.
        x = np.array([[0.13, 0.26, 0.90],
                      [0.28, 0.33, 0.63],
                      [0.31, 0.87, 0.70]])
        m = np.array([[True, False, False],
                      [False, False, False],
                      [True, True, False]], dtype=np.bool_)
        mx = masked_array(x, mask=m)
        mxbig = (mx > 0.5)
        mxsmall = (mx < 0.5)
        #
        self.assertFalse(mxbig.all())
        self.assertTrue(mxbig.any())
        assert_equal(mxbig.all(0), [False, False, True])
        assert_equal(mxbig.all(1), [False, False, True])
        assert_equal(mxbig.any(0), [False, False, True])
        assert_equal(mxbig.any(1), [True, True, True])
        #
        self.assertFalse(mxsmall.all())
        self.assertTrue(mxsmall.any())
        assert_equal(mxsmall.all(0), [True, True, False])
        assert_equal(mxsmall.all(1), [False, False, False])
        assert_equal(mxsmall.any(0), [True, True, False])
        assert_equal(mxsmall.any(1), [True, True, False])

    def test_allany_onmatrices(self):
        x = np.array([[0.13, 0.26, 0.90],
                      [0.28, 0.33, 0.63],
                      [0.31, 0.87, 0.70]])
        X = np.matrix(x)
        m = np.array([[True, False, False],
                      [False, False, False],
                      [True, True, False]], dtype=np.bool_)
        mX = masked_array(X, mask=m)
        mXbig = (mX > 0.5)
        mXsmall = (mX < 0.5)
        #
        self.assertFalse(mXbig.all())
        self.assertTrue(mXbig.any())
        assert_equal(mXbig.all(0), np.matrix([False, False, True]))
        assert_equal(mXbig.all(1), np.matrix([False, False, True]).T)
        assert_equal(mXbig.any(0), np.matrix([False, False, True]))
        assert_equal(mXbig.any(1), np.matrix([True, True, True]).T)
        #
        self.assertFalse(mXsmall.all())
        self.assertTrue(mXsmall.any())
        assert_equal(mXsmall.all(0), np.matrix([True, True, False]))
        assert_equal(mXsmall.all(1), np.matrix([False, False, False]).T)
        assert_equal(mXsmall.any(0), np.matrix([True, True, False]))
        assert_equal(mXsmall.any(1), np.matrix([True, True, False]).T)

    def test_allany_oddities(self):
        # Some fun with all and any
        store = empty((), dtype=bool)
        full = array([1, 2, 3], mask=True)
        #
        self.assertTrue(full.all() is masked)
        full.all(out=store)
        self.assertTrue(store)
        self.assertTrue(store._mask, True)
        self.assertTrue(store is not masked)
        #
        store = empty((), dtype=bool)
        self.assertTrue(full.any() is masked)
        full.any(out=store)
        self.assertTrue(not store)
        self.assertTrue(store._mask, True)
        self.assertTrue(store is not masked)

    def test_argmax_argmin(self):
        # Tests argmin & argmax on MaskedArrays.
        (x, X, XX, m, mx, mX, mXX, m2x, m2X, m2XX) = self.d
        #
        assert_equal(mx.argmin(), 35)
        assert_equal(mX.argmin(), 35)
        assert_equal(m2x.argmin(), 4)
        assert_equal(m2X.argmin(), 4)
        assert_equal(mx.argmax(), 28)
        assert_equal(mX.argmax(), 28)
        assert_equal(m2x.argmax(), 31)
        assert_equal(m2X.argmax(), 31)
        #
        assert_equal(mX.argmin(0), [2, 2, 2, 5, 0, 5])
        assert_equal(m2X.argmin(0), [2, 2, 4, 5, 0, 4])
        assert_equal(mX.argmax(0), [0, 5, 0, 5, 4, 0])
        assert_equal(m2X.argmax(0), [5, 5, 0, 5, 1, 0])
        #
        assert_equal(mX.argmin(1), [4, 1, 0, 0, 5, 5, ])
        assert_equal(m2X.argmin(1), [4, 4, 0, 0, 5, 3])
        assert_equal(mX.argmax(1), [2, 4, 1, 1, 4, 1])
        assert_equal(m2X.argmax(1), [2, 4, 1, 1, 1, 1])

    def test_clip(self):
        # Tests clip on MaskedArrays.
        x = np.array([8.375, 7.545, 8.828, 8.5, 1.757, 5.928,
                      8.43, 7.78, 9.865, 5.878, 8.979, 4.732,
                      3.012, 6.022, 5.095, 3.116, 5.238, 3.957,
                      6.04, 9.63, 7.712, 3.382, 4.489, 6.479,
                      7.189, 9.645, 5.395, 4.961, 9.894, 2.893,
                      7.357, 9.828, 6.272, 3.758, 6.693, 0.993])
        m = np.array([0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1,
                      0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1,
                      1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0])
        mx = array(x, mask=m)
        clipped = mx.clip(2, 8)
        assert_equal(clipped.mask, mx.mask)
        assert_equal(clipped._data, x.clip(2, 8))
        assert_equal(clipped._data, mx._data.clip(2, 8))

    def test_compress(self):
        # test compress
        a = masked_array([1., 2., 3., 4., 5.], fill_value=9999)
        condition = (a > 1.5) & (a < 3.5)
        assert_equal(a.compress(condition), [2., 3.])
        #
        a[[2, 3]] = masked
        b = a.compress(condition)
        assert_equal(b._data, [2., 3.])
        assert_equal(b._mask, [0, 1])
        assert_equal(b.fill_value, 9999)
        assert_equal(b, a[condition])
        #
        condition = (a < 4.)
        b = a.compress(condition)
        assert_equal(b._data, [1., 2., 3.])
        assert_equal(b._mask, [0, 0, 1])
        assert_equal(b.fill_value, 9999)
        assert_equal(b, a[condition])
        #
        a = masked_array([[10, 20, 30], [40, 50, 60]],
                         mask=[[0, 0, 1], [1, 0, 0]])
        b = a.compress(a.ravel() >= 22)
        assert_equal(b._data, [30, 40, 50, 60])
        assert_equal(b._mask, [1, 1, 0, 0])
        #
        x = np.array([3, 1, 2])
        b = a.compress(x >= 2, axis=1)
        assert_equal(b._data, [[10, 30], [40, 60]])
        assert_equal(b._mask, [[0, 1], [1, 0]])

    def test_compressed(self):
        # Tests compressed
        a = array([1, 2, 3, 4], mask=[0, 0, 0, 0])
        b = a.compressed()
        assert_equal(b, a)
        a[0] = masked
        b = a.compressed()
        assert_equal(b, [2, 3, 4])
        #
        a = array(np.matrix([1, 2, 3, 4]), mask=[0, 0, 0, 0])
        b = a.compressed()
        assert_equal(b, a)
        self.assertTrue(isinstance(b, np.matrix))
        a[0, 0] = masked
        b = a.compressed()
        assert_equal(b, [[2, 3, 4]])

    def test_empty(self):
        # Tests empty/like
        datatype = [('a', int), ('b', float), ('c', '|S8')]
        a = masked_array([(1, 1.1, '1.1'), (2, 2.2, '2.2'), (3, 3.3, '3.3')],
                         dtype=datatype)
        assert_equal(len(a.fill_value.item()), len(datatype))
        #
        b = empty_like(a)
        assert_equal(b.shape, a.shape)
        assert_equal(b.fill_value, a.fill_value)
        #
        b = empty(len(a), dtype=datatype)
        assert_equal(b.shape, a.shape)
        assert_equal(b.fill_value, a.fill_value)

    def test_put(self):
        # Tests put.
        d = arange(5)
        n = [0, 0, 0, 1, 1]
        m = make_mask(n)
        x = array(d, mask=m)
        self.assertTrue(x[3] is masked)
        self.assertTrue(x[4] is masked)
        x[[1, 4]] = [10, 40]
        #self.assertTrue(x.mask is not m)
        self.assertTrue(x[3] is masked)
        self.assertTrue(x[4] is not masked)
        assert_equal(x, [0, 10, 2, -1, 40])
        #
        x = masked_array(arange(10), mask=[1, 0, 0, 0, 0] * 2)
        i = [0, 2, 4, 6]
        x.put(i, [6, 4, 2, 0])
        assert_equal(x, asarray([6, 1, 4, 3, 2, 5, 0, 7, 8, 9, ]))
        assert_equal(x.mask, [0, 0, 0, 0, 0, 1, 0, 0, 0, 0])
        x.put(i, masked_array([0, 2, 4, 6], [1, 0, 1, 0]))
        assert_array_equal(x, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, ])
        assert_equal(x.mask, [1, 0, 0, 0, 1, 1, 0, 0, 0, 0])
        #
        x = masked_array(arange(10), mask=[1, 0, 0, 0, 0] * 2)
        put(x, i, [6, 4, 2, 0])
        assert_equal(x, asarray([6, 1, 4, 3, 2, 5, 0, 7, 8, 9, ]))
        assert_equal(x.mask, [0, 0, 0, 0, 0, 1, 0, 0, 0, 0])
        put(x, i, masked_array([0, 2, 4, 6], [1, 0, 1, 0]))
        assert_array_equal(x, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, ])
        assert_equal(x.mask, [1, 0, 0, 0, 1, 1, 0, 0, 0, 0])

    def test_put_hardmask(self):
        # Tests put on hardmask
        d = arange(5)
        n = [0, 0, 0, 1, 1]
        m = make_mask(n)
        xh = array(d + 1, mask=m, hard_mask=True, copy=True)
        xh.put([4, 2, 0, 1, 3], [1, 2, 3, 4, 5])
        assert_equal(xh._data, [3, 4, 2, 4, 5])

    def test_putmask(self):
        x = arange(6) + 1
        mx = array(x, mask=[0, 0, 0, 1, 1, 1])
        mask = [0, 0, 1, 0, 0, 1]
        # w/o mask, w/o masked values
        xx = x.copy()
        putmask(xx, mask, 99)
        assert_equal(xx, [1, 2, 99, 4, 5, 99])
        # w/ mask, w/o masked values
        mxx = mx.copy()
        putmask(mxx, mask, 99)
        assert_equal(mxx._data, [1, 2, 99, 4, 5, 99])
        assert_equal(mxx._mask, [0, 0, 0, 1, 1, 0])
        # w/o mask, w/ masked values
        values = array([10, 20, 30, 40, 50, 60], mask=[1, 1, 1, 0, 0, 0])
        xx = x.copy()
        putmask(xx, mask, values)
        assert_equal(xx._data, [1, 2, 30, 4, 5, 60])
        assert_equal(xx._mask, [0, 0, 1, 0, 0, 0])
        # w/ mask, w/ masked values
        mxx = mx.copy()
        putmask(mxx, mask, values)
        assert_equal(mxx._data, [1, 2, 30, 4, 5, 60])
        assert_equal(mxx._mask, [0, 0, 1, 1, 1, 0])
        # w/ mask, w/ masked values + hardmask
        mxx = mx.copy()
        mxx.harden_mask()
        putmask(mxx, mask, values)
        assert_equal(mxx, [1, 2, 30, 4, 5, 60])

    def test_ravel(self):
        # Tests ravel
        a = array([[1, 2, 3, 4, 5]], mask=[[0, 1, 0, 0, 0]])
        aravel = a.ravel()
        assert_equal(aravel._mask.shape, aravel.shape)
        a = array([0, 0], mask=[1, 1])
        aravel = a.ravel()
        assert_equal(aravel._mask.shape, a.shape)
        a = array(np.matrix([1, 2, 3, 4, 5]), mask=[[0, 1, 0, 0, 0]])
        aravel = a.ravel()
        assert_equal(aravel.shape, (1, 5))
        assert_equal(aravel._mask.shape, a.shape)
        # Checks that small_mask is preserved
        a = array([1, 2, 3, 4], mask=[0, 0, 0, 0], shrink=False)
        assert_equal(a.ravel()._mask, [0, 0, 0, 0])
        # Test that the fill_value is preserved
        a.fill_value = -99
        a.shape = (2, 2)
        ar = a.ravel()
        assert_equal(ar._mask, [0, 0, 0, 0])
        assert_equal(ar._data, [1, 2, 3, 4])
        assert_equal(ar.fill_value, -99)

    def test_reshape(self):
        # Tests reshape
        x = arange(4)
        x[0] = masked
        y = x.reshape(2, 2)
        assert_equal(y.shape, (2, 2,))
        assert_equal(y._mask.shape, (2, 2,))
        assert_equal(x.shape, (4,))
        assert_equal(x._mask.shape, (4,))

    def test_sort(self):
        # Test sort
        x = array([1, 4, 2, 3], mask=[0, 1, 0, 0], dtype=np.uint8)
        #
        sortedx = sort(x)
        assert_equal(sortedx._data, [1, 2, 3, 4])
        assert_equal(sortedx._mask, [0, 0, 0, 1])
        #
        sortedx = sort(x, endwith=False)
        assert_equal(sortedx._data, [4, 1, 2, 3])
        assert_equal(sortedx._mask, [1, 0, 0, 0])
        #
        x.sort()
        assert_equal(x._data, [1, 2, 3, 4])
        assert_equal(x._mask, [0, 0, 0, 1])
        #
        x = array([1, 4, 2, 3], mask=[0, 1, 0, 0], dtype=np.uint8)
        x.sort(endwith=False)
        assert_equal(x._data, [4, 1, 2, 3])
        assert_equal(x._mask, [1, 0, 0, 0])
        #
        x = [1, 4, 2, 3]
        sortedx = sort(x)
        self.assertTrue(not isinstance(sorted, MaskedArray))
        #
        x = array([0, 1, -1, -2, 2], mask=nomask, dtype=np.int8)
        sortedx = sort(x, endwith=False)
        assert_equal(sortedx._data, [-2, -1, 0, 1, 2])
        x = array([0, 1, -1, -2, 2], mask=[0, 1, 0, 0, 1], dtype=np.int8)
        sortedx = sort(x, endwith=False)
        assert_equal(sortedx._data, [1, 2, -2, -1, 0])
        assert_equal(sortedx._mask, [1, 1, 0, 0, 0])

    def test_sort_2d(self):
        # Check sort of 2D array.
        # 2D array w/o mask
        a = masked_array([[8, 4, 1], [2, 0, 9]])
        a.sort(0)
        assert_equal(a, [[2, 0, 1], [8, 4, 9]])
        a = masked_array([[8, 4, 1], [2, 0, 9]])
        a.sort(1)
        assert_equal(a, [[1, 4, 8], [0, 2, 9]])
        # 2D array w/mask
        a = masked_array([[8, 4, 1], [2, 0, 9]], mask=[[1, 0, 0], [0, 0, 1]])
        a.sort(0)
        assert_equal(a, [[2, 0, 1], [8, 4, 9]])
        assert_equal(a._mask, [[0, 0, 0], [1, 0, 1]])
        a = masked_array([[8, 4, 1], [2, 0, 9]], mask=[[1, 0, 0], [0, 0, 1]])
        a.sort(1)
        assert_equal(a, [[1, 4, 8], [0, 2, 9]])
        assert_equal(a._mask, [[0, 0, 1], [0, 0, 1]])
        # 3D
        a = masked_array([[[7, 8, 9], [4, 5, 6], [1, 2, 3]],
                          [[1, 2, 3], [7, 8, 9], [4, 5, 6]],
                          [[7, 8, 9], [1, 2, 3], [4, 5, 6]],
                          [[4, 5, 6], [1, 2, 3], [7, 8, 9]]])
        a[a % 4 == 0] = masked
        am = a.copy()
        an = a.filled(99)
        am.sort(0)
        an.sort(0)
        assert_equal(am, an)
        am = a.copy()
        an = a.filled(99)
        am.sort(1)
        an.sort(1)
        assert_equal(am, an)
        am = a.copy()
        an = a.filled(99)
        am.sort(2)
        an.sort(2)
        assert_equal(am, an)

    def test_sort_flexible(self):
        # Test sort on flexible dtype.
        a = array(
            data=[(3, 3), (3, 2), (2, 2), (2, 1), (1, 0), (1, 1), (1, 2)],
            mask=[(0, 0), (0, 1), (0, 0), (0, 0), (1, 0), (0, 0), (0, 0)],
            dtype=[('A', int), ('B', int)])
        #
        test = sort(a)
        b = array(
            data=[(1, 1), (1, 2), (2, 1), (2, 2), (3, 3), (3, 2), (1, 0)],
            mask=[(0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 1), (1, 0)],
            dtype=[('A', int), ('B', int)])
        assert_equal(test, b)
        assert_equal(test.mask, b.mask)
        #
        test = sort(a, endwith=False)
        b = array(
            data=[(1, 0), (1, 1), (1, 2), (2, 1), (2, 2), (3, 2), (3, 3), ],
            mask=[(1, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 1), (0, 0), ],
            dtype=[('A', int), ('B', int)])
        assert_equal(test, b)
        assert_equal(test.mask, b.mask)

    def test_argsort(self):
        # Test argsort
        a = array([1, 5, 2, 4, 3], mask=[1, 0, 0, 1, 0])
        assert_equal(np.argsort(a), argsort(a))

    def test_squeeze(self):
        # Check squeeze
        data = masked_array([[1, 2, 3]])
        assert_equal(data.squeeze(), [1, 2, 3])
        data = masked_array([[1, 2, 3]], mask=[[1, 1, 1]])
        assert_equal(data.squeeze(), [1, 2, 3])
        assert_equal(data.squeeze()._mask, [1, 1, 1])
        data = masked_array([[1]], mask=True)
        self.assertTrue(data.squeeze() is masked)

    def test_swapaxes(self):
        # Tests swapaxes on MaskedArrays.
        x = np.array([8.375, 7.545, 8.828, 8.5, 1.757, 5.928,
                      8.43, 7.78, 9.865, 5.878, 8.979, 4.732,
                      3.012, 6.022, 5.095, 3.116, 5.238, 3.957,
                      6.04, 9.63, 7.712, 3.382, 4.489, 6.479,
                      7.189, 9.645, 5.395, 4.961, 9.894, 2.893,
                      7.357, 9.828, 6.272, 3.758, 6.693, 0.993])
        m = np.array([0, 1, 0, 1, 0, 0,
                      1, 0, 1, 1, 0, 1,
                      0, 0, 0, 1, 0, 1,
                      0, 0, 0, 1, 1, 1,
                      1, 0, 0, 1, 0, 0,
                      0, 0, 1, 0, 1, 0])
        mX = array(x, mask=m).reshape(6, 6)
        mXX = mX.reshape(3, 2, 2, 3)
        #
        mXswapped = mX.swapaxes(0, 1)
        assert_equal(mXswapped[-1], mX[:, -1])

        mXXswapped = mXX.swapaxes(0, 2)
        assert_equal(mXXswapped.shape, (2, 2, 3, 3))

    def test_take(self):
        # Tests take
        x = masked_array([10, 20, 30, 40], [0, 1, 0, 1])
        assert_equal(x.take([0, 0, 3]), masked_array([10, 10, 40], [0, 0, 1]))
        assert_equal(x.take([0, 0, 3]), x[[0, 0, 3]])
        assert_equal(x.take([[0, 1], [0, 1]]),
                     masked_array([[10, 20], [10, 20]], [[0, 1], [0, 1]]))
        #
        x = array([[10, 20, 30], [40, 50, 60]], mask=[[0, 0, 1], [1, 0, 0, ]])
        assert_equal(x.take([0, 2], axis=1),
                     array([[10, 30], [40, 60]], mask=[[0, 1], [1, 0]]))
        assert_equal(take(x, [0, 2], axis=1),
                     array([[10, 30], [40, 60]], mask=[[0, 1], [1, 0]]))

    def test_take_masked_indices(self):
        # Test take w/ masked indices
        a = np.array((40, 18, 37, 9, 22))
        indices = np.arange(3)[None,:] + np.arange(5)[:, None]
        mindices = array(indices, mask=(indices >= len(a)))
        # No mask
        test = take(a, mindices, mode='clip')
        ctrl = array([[40, 18, 37],
                      [18, 37, 9],
                      [37, 9, 22],
                      [9, 22, 22],
                      [22, 22, 22]])
        assert_equal(test, ctrl)
        # Masked indices
        test = take(a, mindices)
        ctrl = array([[40, 18, 37],
                      [18, 37, 9],
                      [37, 9, 22],
                      [9, 22, 40],
                      [22, 40, 40]])
        ctrl[3, 2] = ctrl[4, 1] = ctrl[4, 2] = masked
        assert_equal(test, ctrl)
        assert_equal(test.mask, ctrl.mask)
        # Masked input + masked indices
        a = array((40, 18, 37, 9, 22), mask=(0, 1, 0, 0, 0))
        test = take(a, mindices)
        ctrl[0, 1] = ctrl[1, 0] = masked
        assert_equal(test, ctrl)
        assert_equal(test.mask, ctrl.mask)

    def test_tolist(self):
        # Tests to list
        # ... on 1D
        x = array(np.arange(12))
        x[[1, -2]] = masked
        xlist = x.tolist()
        self.assertTrue(xlist[1] is None)
        self.assertTrue(xlist[-2] is None)
        # ... on 2D
        x.shape = (3, 4)
        xlist = x.tolist()
        ctrl = [[0, None, 2, 3], [4, 5, 6, 7], [8, 9, None, 11]]
        assert_equal(xlist[0], [0, None, 2, 3])
        assert_equal(xlist[1], [4, 5, 6, 7])
        assert_equal(xlist[2], [8, 9, None, 11])
        assert_equal(xlist, ctrl)
        # ... on structured array w/ masked records
        x = array(list(zip([1, 2, 3],
                           [1.1, 2.2, 3.3],
                           ['one', 'two', 'thr'])),
                  dtype=[('a', int), ('b', float), ('c', '|S8')])
        x[-1] = masked
        assert_equal(x.tolist(),
                     [(1, 1.1, asbytes('one')),
                      (2, 2.2, asbytes('two')),
                      (None, None, None)])
        # ... on structured array w/ masked fields
        a = array([(1, 2,), (3, 4)], mask=[(0, 1), (0, 0)],
                  dtype=[('a', int), ('b', int)])
        test = a.tolist()
        assert_equal(test, [[1, None], [3, 4]])
        # ... on mvoid
        a = a[0]
        test = a.tolist()
        assert_equal(test, [1, None])

    def test_tolist_specialcase(self):
        # Test mvoid.tolist: make sure we return a standard Python object
        a = array([(0, 1), (2, 3)], dtype=[('a', int), ('b', int)])
        # w/o mask: each entry is a np.void whose elements are standard Python
        for entry in a:
            for item in entry.tolist():
                assert_(not isinstance(item, np.generic))
        # w/ mask: each entry is a ma.void whose elements should be
        # standard Python
        a.mask[0] = (0, 1)
        for entry in a:
            for item in entry.tolist():
                assert_(not isinstance(item, np.generic))

    def test_toflex(self):
        # Test the conversion to records
        data = arange(10)
        record = data.toflex()
        assert_equal(record['_data'], data._data)
        assert_equal(record['_mask'], data._mask)
        #
        data[[0, 1, 2, -1]] = masked
        record = data.toflex()
        assert_equal(record['_data'], data._data)
        assert_equal(record['_mask'], data._mask)
        #
        ndtype = [('i', int), ('s', '|S3'), ('f', float)]
        data = array([(i, s, f) for (i, s, f) in zip(np.arange(10),
                                                     'ABCDEFGHIJKLM',
                                                     np.random.rand(10))],
                     dtype=ndtype)
        data[[0, 1, 2, -1]] = masked
        record = data.toflex()
        assert_equal(record['_data'], data._data)
        assert_equal(record['_mask'], data._mask)
        #
        ndtype = np.dtype("int, (2,3)float, float")
        data = array([(i, f, ff) for (i, f, ff) in zip(np.arange(10),
                                                       np.random.rand(10),
                                                       np.random.rand(10))],
                     dtype=ndtype)
        data[[0, 1, 2, -1]] = masked
        record = data.toflex()
        assert_equal_records(record['_data'], data._data)
        assert_equal_records(record['_mask'], data._mask)

    def test_fromflex(self):
        # Test the reconstruction of a masked_array from a record
        a = array([1, 2, 3])
        test = fromflex(a.toflex())
        assert_equal(test, a)
        assert_equal(test.mask, a.mask)
        #
        a = array([1, 2, 3], mask=[0, 0, 1])
        test = fromflex(a.toflex())
        assert_equal(test, a)
        assert_equal(test.mask, a.mask)
        #
        a = array([(1, 1.), (2, 2.), (3, 3.)], mask=[(1, 0), (0, 0), (0, 1)],
                  dtype=[('A', int), ('B', float)])
        test = fromflex(a.toflex())
        assert_equal(test, a)
        assert_equal(test.data, a.data)

    def test_arraymethod(self):
        # Test a _arraymethod w/ n argument
        marray = masked_array([[1, 2, 3, 4, 5]], mask=[0, 0, 1, 0, 0])
        control = masked_array([[1], [2], [3], [4], [5]],
                               mask=[0, 0, 1, 0, 0])
        assert_equal(marray.T, control)
        assert_equal(marray.transpose(), control)
        #
        assert_equal(MaskedArray.cumsum(marray.T, 0), control.cumsum(0))


#------------------------------------------------------------------------------
class TestMaskedArrayMathMethods(TestCase):

    def setUp(self):
        # Base data definition.
        x = np.array([8.375, 7.545, 8.828, 8.5, 1.757, 5.928,
                      8.43, 7.78, 9.865, 5.878, 8.979, 4.732,
                      3.012, 6.022, 5.095, 3.116, 5.238, 3.957,
                      6.04, 9.63, 7.712, 3.382, 4.489, 6.479,
                      7.189, 9.645, 5.395, 4.961, 9.894, 2.893,
                      7.357, 9.828, 6.272, 3.758, 6.693, 0.993])
        X = x.reshape(6, 6)
        XX = x.reshape(3, 2, 2, 3)

        m = np.array([0, 1, 0, 1, 0, 0,
                     1, 0, 1, 1, 0, 1,
                     0, 0, 0, 1, 0, 1,
                     0, 0, 0, 1, 1, 1,
                     1, 0, 0, 1, 0, 0,
                     0, 0, 1, 0, 1, 0])
        mx = array(data=x, mask=m)
        mX = array(data=X, mask=m.reshape(X.shape))
        mXX = array(data=XX, mask=m.reshape(XX.shape))

        m2 = np.array([1, 1, 0, 1, 0, 0,
                      1, 1, 1, 1, 0, 1,
                      0, 0, 1, 1, 0, 1,
                      0, 0, 0, 1, 1, 1,
                      1, 0, 0, 1, 1, 0,
                      0, 0, 1, 0, 1, 1])
        m2x = array(data=x, mask=m2)
        m2X = array(data=X, mask=m2.reshape(X.shape))
        m2XX = array(data=XX, mask=m2.reshape(XX.shape))
        self.d = (x, X, XX, m, mx, mX, mXX, m2x, m2X, m2XX)

    def test_cumsumprod(self):
        # Tests cumsum & cumprod on MaskedArrays.
        (x, X, XX, m, mx, mX, mXX, m2x, m2X, m2XX) = self.d
        mXcp = mX.cumsum(0)
        assert_equal(mXcp._data, mX.filled(0).cumsum(0))
        mXcp = mX.cumsum(1)
        assert_equal(mXcp._data, mX.filled(0).cumsum(1))
        #
        mXcp = mX.cumprod(0)
        assert_equal(mXcp._data, mX.filled(1).cumprod(0))
        mXcp = mX.cumprod(1)
        assert_equal(mXcp._data, mX.filled(1).cumprod(1))

    def test_cumsumprod_with_output(self):
        # Tests cumsum/cumprod w/ output
        xm = array(np.random.uniform(0, 10, 12)).reshape(3, 4)
        xm[:, 0] = xm[0] = xm[-1, -1] = masked
        #
        for funcname in ('cumsum', 'cumprod'):
            npfunc = getattr(np, funcname)
            xmmeth = getattr(xm, funcname)

            # A ndarray as explicit input
            output = np.empty((3, 4), dtype=float)
            output.fill(-9999)
            result = npfunc(xm, axis=0, out=output)
            # ... the result should be the given output
            self.assertTrue(result is output)
            assert_equal(result, xmmeth(axis=0, out=output))
            #
            output = empty((3, 4), dtype=int)
            result = xmmeth(axis=0, out=output)
            self.assertTrue(result is output)

    def test_ptp(self):
        # Tests ptp on MaskedArrays.
        (x, X, XX, m, mx, mX, mXX, m2x, m2X, m2XX) = self.d
        (n, m) = X.shape
        assert_equal(mx.ptp(), mx.compressed().ptp())
        rows = np.zeros(n, np.float)
        cols = np.zeros(m, np.float)
        for k in range(m):
            cols[k] = mX[:, k].compressed().ptp()
        for k in range(n):
            rows[k] = mX[k].compressed().ptp()
        assert_equal(mX.ptp(0), cols)
        assert_equal(mX.ptp(1), rows)

    def test_sum_object(self):
        # Test sum on object dtype
        a = masked_array([1, 2, 3], mask=[1, 0, 0], dtype=np.object)
        assert_equal(a.sum(), 5)
        a = masked_array([[1, 2, 3], [4, 5, 6]], dtype=object)
        assert_equal(a.sum(axis=0), [5, 7, 9])

    def test_prod_object(self):
        # Test prod on object dtype
        a = masked_array([1, 2, 3], mask=[1, 0, 0], dtype=np.object)
        assert_equal(a.prod(), 2 * 3)
        a = masked_array([[1, 2, 3], [4, 5, 6]], dtype=object)
        assert_equal(a.prod(axis=0), [4, 10, 18])

    def test_meananom_object(self):
        # Test mean/anom on object dtype
        a = masked_array([1, 2, 3], dtype=np.object)
        assert_equal(a.mean(), 2)
        assert_equal(a.anom(), [-1, 0, 1])

    def test_trace(self):
        # Tests trace on MaskedArrays.
        (x, X, XX, m, mx, mX, mXX, m2x, m2X, m2XX) = self.d
        mXdiag = mX.diagonal()
        assert_equal(mX.trace(), mX.diagonal().compressed().sum())
        assert_almost_equal(mX.trace(),
                            X.trace() - sum(mXdiag.mask * X.diagonal(),
                                            axis=0))

    def test_varstd(self):
        # Tests var & std on MaskedArrays.
        (x, X, XX, m, mx, mX, mXX, m2x, m2X, m2XX) = self.d
        assert_almost_equal(mX.var(axis=None), mX.compressed().var())
        assert_almost_equal(mX.std(axis=None), mX.compressed().std())
        assert_almost_equal(mX.std(axis=None, ddof=1),
                            mX.compressed().std(ddof=1))
        assert_almost_equal(mX.var(axis=None, ddof=1),
                            mX.compressed().var(ddof=1))
        assert_equal(mXX.var(axis=3).shape, XX.var(axis=3).shape)
        assert_equal(mX.var().shape, X.var().shape)
        (mXvar0, mXvar1) = (mX.var(axis=0), mX.var(axis=1))
        assert_almost_equal(mX.var(axis=None, ddof=2),
                            mX.compressed().var(ddof=2))
        assert_almost_equal(mX.std(axis=None, ddof=2),
                            mX.compressed().std(ddof=2))
        for k in range(6):
            assert_almost_equal(mXvar1[k], mX[k].compressed().var())
            assert_almost_equal(mXvar0[k], mX[:, k].compressed().var())
            assert_almost_equal(np.sqrt(mXvar0[k]),
                                mX[:, k].compressed().std())

    def test_varstd_specialcases(self):
        # Test a special case for var
        nout = np.array(-1, dtype=float)
        mout = array(-1, dtype=float)
        #
        x = array(arange(10), mask=True)
        for methodname in ('var', 'std'):
            method = getattr(x, methodname)
            self.assertTrue(method() is masked)
            self.assertTrue(method(0) is masked)
            self.assertTrue(method(-1) is masked)
            # Using a masked array as explicit output
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                _ = method(out=mout)
            self.assertTrue(mout is not masked)
            assert_equal(mout.mask, True)
            # Using a ndarray as explicit output
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                _ = method(out=nout)
            self.assertTrue(np.isnan(nout))
        #
        x = array(arange(10), mask=True)
        x[-1] = 9
        for methodname in ('var', 'std'):
            method = getattr(x, methodname)
            self.assertTrue(method(ddof=1) is masked)
            self.assertTrue(method(0, ddof=1) is masked)
            self.assertTrue(method(-1, ddof=1) is masked)
            # Using a masked array as explicit output
            method(out=mout, ddof=1)
            self.assertTrue(mout is not masked)
            assert_equal(mout.mask, True)
            # Using a ndarray as explicit output
            method(out=nout, ddof=1)
            self.assertTrue(np.isnan(nout))

    def test_varstd_ddof(self):
        a = array([[1, 1, 0], [1, 1, 0]], mask=[[0, 0, 1], [0, 0, 1]])
        test = a.std(axis=0, ddof=0)
        assert_equal(test.filled(0), [0, 0, 0])
        assert_equal(test.mask, [0, 0, 1])
        test = a.std(axis=0, ddof=1)
        assert_equal(test.filled(0), [0, 0, 0])
        assert_equal(test.mask, [0, 0, 1])
        test = a.std(axis=0, ddof=2)
        assert_equal(test.filled(0), [0, 0, 0])
        assert_equal(test.mask, [1, 1, 1])

    def test_diag(self):
        # Test diag
        x = arange(9).reshape((3, 3))
        x[1, 1] = masked
        out = np.diag(x)
        assert_equal(out, [0, 4, 8])
        out = diag(x)
        assert_equal(out, [0, 4, 8])
        assert_equal(out.mask, [0, 1, 0])
        out = diag(out)
        control = array([[0, 0, 0], [0, 4, 0], [0, 0, 8]],
                        mask=[[0, 0, 0], [0, 1, 0], [0, 0, 0]])
        assert_equal(out, control)

    def test_axis_methods_nomask(self):
        # Test the combination nomask & methods w/ axis
        a = array([[1, 2, 3], [4, 5, 6]])
        #
        assert_equal(a.sum(0), [5, 7, 9])
        assert_equal(a.sum(-1), [6, 15])
        assert_equal(a.sum(1), [6, 15])
        #
        assert_equal(a.prod(0), [4, 10, 18])
        assert_equal(a.prod(-1), [6, 120])
        assert_equal(a.prod(1), [6, 120])
        #
        assert_equal(a.min(0), [1, 2, 3])
        assert_equal(a.min(-1), [1, 4])
        assert_equal(a.min(1), [1, 4])
        #
        assert_equal(a.max(0), [4, 5, 6])
        assert_equal(a.max(-1), [3, 6])
        assert_equal(a.max(1), [3, 6])


#------------------------------------------------------------------------------
class TestMaskedArrayMathMethodsComplex(TestCase):
    # Test class for miscellaneous MaskedArrays methods.
    def setUp(self):
        # Base data definition.
        x = np.array([8.375j, 7.545j, 8.828j, 8.5j, 1.757j, 5.928,
                      8.43, 7.78, 9.865, 5.878, 8.979, 4.732,
                      3.012, 6.022, 5.095, 3.116, 5.238, 3.957,
                      6.04, 9.63, 7.712, 3.382, 4.489, 6.479j,
                      7.189j, 9.645, 5.395, 4.961, 9.894, 2.893,
                      7.357, 9.828, 6.272, 3.758, 6.693, 0.993j])
        X = x.reshape(6, 6)
        XX = x.reshape(3, 2, 2, 3)

        m = np.array([0, 1, 0, 1, 0, 0,
                     1, 0, 1, 1, 0, 1,
                     0, 0, 0, 1, 0, 1,
                     0, 0, 0, 1, 1, 1,
                     1, 0, 0, 1, 0, 0,
                     0, 0, 1, 0, 1, 0])
        mx = array(data=x, mask=m)
        mX = array(data=X, mask=m.reshape(X.shape))
        mXX = array(data=XX, mask=m.reshape(XX.shape))

        m2 = np.array([1, 1, 0, 1, 0, 0,
                      1, 1, 1, 1, 0, 1,
                      0, 0, 1, 1, 0, 1,
                      0, 0, 0, 1, 1, 1,
                      1, 0, 0, 1, 1, 0,
                      0, 0, 1, 0, 1, 1])
        m2x = array(data=x, mask=m2)
        m2X = array(data=X, mask=m2.reshape(X.shape))
        m2XX = array(data=XX, mask=m2.reshape(XX.shape))
        self.d = (x, X, XX, m, mx, mX, mXX, m2x, m2X, m2XX)

    def test_varstd(self):
        # Tests var & std on MaskedArrays.
        (x, X, XX, m, mx, mX, mXX, m2x, m2X, m2XX) = self.d
        assert_almost_equal(mX.var(axis=None), mX.compressed().var())
        assert_almost_equal(mX.std(axis=None), mX.compressed().std())
        assert_equal(mXX.var(axis=3).shape, XX.var(axis=3).shape)
        assert_equal(mX.var().shape, X.var().shape)
        (mXvar0, mXvar1) = (mX.var(axis=0), mX.var(axis=1))
        assert_almost_equal(mX.var(axis=None, ddof=2),
                            mX.compressed().var(ddof=2))
        assert_almost_equal(mX.std(axis=None, ddof=2),
                            mX.compressed().std(ddof=2))
        for k in range(6):
            assert_almost_equal(mXvar1[k], mX[k].compressed().var())
            assert_almost_equal(mXvar0[k], mX[:, k].compressed().var())
            assert_almost_equal(np.sqrt(mXvar0[k]),
                                mX[:, k].compressed().std())


#------------------------------------------------------------------------------
class TestMaskedArrayFunctions(TestCase):
    # Test class for miscellaneous functions.

    def setUp(self):
        x = np.array([1., 1., 1., -2., pi/2.0, 4., 5., -10., 10., 1., 2., 3.])
        y = np.array([5., 0., 3., 2., -1., -4., 0., -10., 10., 1., 0., 3.])
        m1 = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]
        m2 = [0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1]
        xm = masked_array(x, mask=m1)
        ym = masked_array(y, mask=m2)
        xm.set_fill_value(1e+20)
        self.info = (xm, ym)

    def test_masked_where_bool(self):
        x = [1, 2]
        y = masked_where(False, x)
        assert_equal(y, [1, 2])
        assert_equal(y[1], 2)

    def test_masked_equal_wlist(self):
        x = [1, 2, 3]
        mx = masked_equal(x, 3)
        assert_equal(mx, x)
        assert_equal(mx._mask, [0, 0, 1])
        mx = masked_not_equal(x, 3)
        assert_equal(mx, x)
        assert_equal(mx._mask, [1, 1, 0])

    def test_masked_equal_fill_value(self):
        x = [1, 2, 3]
        mx = masked_equal(x, 3)
        assert_equal(mx._mask, [0, 0, 1])
        assert_equal(mx.fill_value, 3)

    def test_masked_where_condition(self):
        # Tests masking functions.
        x = array([1., 2., 3., 4., 5.])
        x[2] = masked
        assert_equal(masked_where(greater(x, 2), x), masked_greater(x, 2))
        assert_equal(masked_where(greater_equal(x, 2), x),
                     masked_greater_equal(x, 2))
        assert_equal(masked_where(less(x, 2), x), masked_less(x, 2))
        assert_equal(masked_where(less_equal(x, 2), x),
                     masked_less_equal(x, 2))
        assert_equal(masked_where(not_equal(x, 2), x), masked_not_equal(x, 2))
        assert_equal(masked_where(equal(x, 2), x), masked_equal(x, 2))
        assert_equal(masked_where(not_equal(x, 2), x), masked_not_equal(x, 2))
        assert_equal(masked_where([1, 1, 0, 0, 0], [1, 2, 3, 4, 5]),
                     [99, 99, 3, 4, 5])

    def test_masked_where_oddities(self):
        # Tests some generic features.
        atest = ones((10, 10, 10), dtype=float)
        btest = zeros(atest.shape, MaskType)
        ctest = masked_where(btest, atest)
        assert_equal(atest, ctest)

    def test_masked_where_shape_constraint(self):
        a = arange(10)
        try:
            test = masked_equal(1, a)
        except IndexError:
            pass
        else:
            raise AssertionError("Should have failed...")
        test = masked_equal(a, 1)
        assert_equal(test.mask, [0, 1, 0, 0, 0, 0, 0, 0, 0, 0])

    def test_masked_otherfunctions(self):
        assert_equal(masked_inside(list(range(5)), 1, 3),
                     [0, 199, 199, 199, 4])
        assert_equal(masked_outside(list(range(5)), 1, 3), [199, 1, 2, 3, 199])
        assert_equal(masked_inside(array(list(range(5)),
                                         mask=[1, 0, 0, 0, 0]), 1, 3).mask,
                     [1, 1, 1, 1, 0])
        assert_equal(masked_outside(array(list(range(5)),
                                          mask=[0, 1, 0, 0, 0]), 1, 3).mask,
                     [1, 1, 0, 0, 1])
        assert_equal(masked_equal(array(list(range(5)),
                                        mask=[1, 0, 0, 0, 0]), 2).mask,
                     [1, 0, 1, 0, 0])
        assert_equal(masked_not_equal(array([2, 2, 1, 2, 1],
                                            mask=[1, 0, 0, 0, 0]), 2).mask,
                     [1, 0, 1, 0, 1])

    def test_round(self):
        a = array([1.23456, 2.34567, 3.45678, 4.56789, 5.67890],
                  mask=[0, 1, 0, 0, 0])
        assert_equal(a.round(), [1., 2., 3., 5., 6.])
        assert_equal(a.round(1), [1.2, 2.3, 3.5, 4.6, 5.7])
        assert_equal(a.round(3), [1.235, 2.346, 3.457, 4.568, 5.679])
        b = empty_like(a)
        a.round(out=b)
        assert_equal(b, [1., 2., 3., 5., 6.])

        x = array([1., 2., 3., 4., 5.])
        c = array([1, 1, 1, 0, 0])
        x[2] = masked
        z = where(c, x, -x)
        assert_equal(z, [1., 2., 0., -4., -5])
        c[0] = masked
        z = where(c, x, -x)
        assert_equal(z, [1., 2., 0., -4., -5])
        assert_(z[0] is masked)
        assert_(z[1] is not masked)
        assert_(z[2] is masked)

    def test_round_with_output(self):
        # Testing round with an explicit output

        xm = array(np.random.uniform(0, 10, 12)).reshape(3, 4)
        xm[:, 0] = xm[0] = xm[-1, -1] = masked

        # A ndarray as explicit input
        output = np.empty((3, 4), dtype=float)
        output.fill(-9999)
        result = np.round(xm, decimals=2, out=output)
        # ... the result should be the given output
        self.assertTrue(result is output)
        assert_equal(result, xm.round(decimals=2, out=output))
        #
        output = empty((3, 4), dtype=float)
        result = xm.round(decimals=2, out=output)
        self.assertTrue(result is output)

    def test_identity(self):
        a = identity(5)
        self.assertTrue(isinstance(a, MaskedArray))
        assert_equal(a, np.identity(5))

    def test_power(self):
        x = -1.1
        assert_almost_equal(power(x, 2.), 1.21)
        self.assertTrue(power(x, masked) is masked)
        x = array([-1.1, -1.1, 1.1, 1.1, 0.])
        b = array([0.5, 2., 0.5, 2., -1.], mask=[0, 0, 0, 0, 1])
        y = power(x, b)
        assert_almost_equal(y, [0, 1.21, 1.04880884817, 1.21, 0.])
        assert_equal(y._mask, [1, 0, 0, 0, 1])
        b.mask = nomask
        y = power(x, b)
        assert_equal(y._mask, [1, 0, 0, 0, 1])
        z = x ** b
        assert_equal(z._mask, y._mask)
        assert_almost_equal(z, y)
        assert_almost_equal(z._data, y._data)
        x **= b
        assert_equal(x._mask, y._mask)
        assert_almost_equal(x, y)
        assert_almost_equal(x._data, y._data)

    def test_power_w_broadcasting(self):
        # Test power w/ broadcasting
        a2 = np.array([[1., 2., 3.], [4., 5., 6.]])
        a2m = array(a2, mask=[[1, 0, 0], [0, 0, 1]])
        b1 = np.array([2, 4, 3])
        b2 = np.array([b1, b1])
        b2m = array(b2, mask=[[0, 1, 0], [0, 1, 0]])
        #
        ctrl = array([[1 ** 2, 2 ** 4, 3 ** 3], [4 ** 2, 5 ** 4, 6 ** 3]],
                     mask=[[1, 1, 0], [0, 1, 1]])
        # No broadcasting, base & exp w/ mask
        test = a2m ** b2m
        assert_equal(test, ctrl)
        assert_equal(test.mask, ctrl.mask)
        # No broadcasting, base w/ mask, exp w/o mask
        test = a2m ** b2
        assert_equal(test, ctrl)
        assert_equal(test.mask, a2m.mask)
        # No broadcasting, base w/o mask, exp w/ mask
        test = a2 ** b2m
        assert_equal(test, ctrl)
        assert_equal(test.mask, b2m.mask)
        #
        ctrl = array([[2 ** 2, 4 ** 4, 3 ** 3], [2 ** 2, 4 ** 4, 3 ** 3]],
                     mask=[[0, 1, 0], [0, 1, 0]])
        test = b1 ** b2m
        assert_equal(test, ctrl)
        assert_equal(test.mask, ctrl.mask)
        test = b2m ** b1
        assert_equal(test, ctrl)
        assert_equal(test.mask, ctrl.mask)

    def test_where(self):
        # Test the where function
        x = np.array([1., 1., 1., -2., pi/2.0, 4., 5., -10., 10., 1., 2., 3.])
        y = np.array([5., 0., 3., 2., -1., -4., 0., -10., 10., 1., 0., 3.])
        m1 = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]
        m2 = [0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1]
        xm = masked_array(x, mask=m1)
        ym = masked_array(y, mask=m2)
        xm.set_fill_value(1e+20)
        #
        d = where(xm > 2, xm, -9)
        assert_equal(d, [-9., -9., -9., -9., -9., 4.,
                         -9., -9., 10., -9., -9., 3.])
        assert_equal(d._mask, xm._mask)
        d = where(xm > 2, -9, ym)
        assert_equal(d, [5., 0., 3., 2., -1., -9.,
                         -9., -10., -9., 1., 0., -9.])
        assert_equal(d._mask, [1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0])
        d = where(xm > 2, xm, masked)
        assert_equal(d, [-9., -9., -9., -9., -9., 4.,
                         -9., -9., 10., -9., -9., 3.])
        tmp = xm._mask.copy()
        tmp[(xm <= 2).filled(True)] = True
        assert_equal(d._mask, tmp)
        #
        ixm = xm.astype(int)
        d = where(ixm > 2, ixm, masked)
        assert_equal(d, [-9, -9, -9, -9, -9, 4, -9, -9, 10, -9, -9, 3])
        assert_equal(d.dtype, ixm.dtype)

    def test_where_with_masked_choice(self):
        x = arange(10)
        x[3] = masked
        c = x >= 8
        # Set False to masked
        z = where(c, x, masked)
        assert_(z.dtype is x.dtype)
        assert_(z[3] is masked)
        assert_(z[4] is masked)
        assert_(z[7] is masked)
        assert_(z[8] is not masked)
        assert_(z[9] is not masked)
        assert_equal(x, z)
        # Set True to masked
        z = where(c, masked, x)
        assert_(z.dtype is x.dtype)
        assert_(z[3] is masked)
        assert_(z[4] is not masked)
        assert_(z[7] is not masked)
        assert_(z[8] is masked)
        assert_(z[9] is masked)

    def test_where_with_masked_condition(self):
        x = array([1., 2., 3., 4., 5.])
        c = array([1, 1, 1, 0, 0])
        x[2] = masked
        z = where(c, x, -x)
        assert_equal(z, [1., 2., 0., -4., -5])
        c[0] = masked
        z = where(c, x, -x)
        assert_equal(z, [1., 2., 0., -4., -5])
        assert_(z[0] is masked)
        assert_(z[1] is not masked)
        assert_(z[2] is masked)
        #
        x = arange(1, 6)
        x[-1] = masked
        y = arange(1, 6) * 10
        y[2] = masked
        c = array([1, 1, 1, 0, 0], mask=[1, 0, 0, 0, 0])
        cm = c.filled(1)
        z = where(c, x, y)
        zm = where(cm, x, y)
        assert_equal(z, zm)
        assert_(getmask(zm) is nomask)
        assert_equal(zm, [1, 2, 3, 40, 50])
        z = where(c, masked, 1)
        assert_equal(z, [99, 99, 99, 1, 1])
        z = where(c, 1, masked)
        assert_equal(z, [99, 1, 1, 99, 99])

    def test_where_type(self):
        # Test the type conservation with where
        x = np.arange(4, dtype=np.int32)
        y = np.arange(4, dtype=np.float32) * 2.2
        test = where(x > 1.5, y, x).dtype
        control = np.find_common_type([np.int32, np.float32], [])
        assert_equal(test, control)

    def test_choose(self):
        # Test choose
        choices = [[0, 1, 2, 3], [10, 11, 12, 13],
                   [20, 21, 22, 23], [30, 31, 32, 33]]
        chosen = choose([2, 3, 1, 0], choices)
        assert_equal(chosen, array([20, 31, 12, 3]))
        chosen = choose([2, 4, 1, 0], choices, mode='clip')
        assert_equal(chosen, array([20, 31, 12, 3]))
        chosen = choose([2, 4, 1, 0], choices, mode='wrap')
        assert_equal(chosen, array([20, 1, 12, 3]))
        # Check with some masked indices
        indices_ = array([2, 4, 1, 0], mask=[1, 0, 0, 1])
        chosen = choose(indices_, choices, mode='wrap')
        assert_equal(chosen, array([99, 1, 12, 99]))
        assert_equal(chosen.mask, [1, 0, 0, 1])
        # Check with some masked choices
        choices = array(choices, mask=[[0, 0, 0, 1], [1, 1, 0, 1],
                                       [1, 0, 0, 0], [0, 0, 0, 0]])
        indices_ = [2, 3, 1, 0]
        chosen = choose(indices_, choices, mode='wrap')
        assert_equal(chosen, array([20, 31, 12, 3]))
        assert_equal(chosen.mask, [1, 0, 0, 1])

    def test_choose_with_out(self):
        # Test choose with an explicit out keyword
        choices = [[0, 1, 2, 3], [10, 11, 12, 13],
                   [20, 21, 22, 23], [30, 31, 32, 33]]
        store = empty(4, dtype=int)
        chosen = choose([2, 3, 1, 0], choices, out=store)
        assert_equal(store, array([20, 31, 12, 3]))
        self.assertTrue(store is chosen)
        # Check with some masked indices + out
        store = empty(4, dtype=int)
        indices_ = array([2, 3, 1, 0], mask=[1, 0, 0, 1])
        chosen = choose(indices_, choices, mode='wrap', out=store)
        assert_equal(store, array([99, 31, 12, 99]))
        assert_equal(store.mask, [1, 0, 0, 1])
        # Check with some masked choices + out ina ndarray !
        choices = array(choices, mask=[[0, 0, 0, 1], [1, 1, 0, 1],
                                       [1, 0, 0, 0], [0, 0, 0, 0]])
        indices_ = [2, 3, 1, 0]
        store = empty(4, dtype=int).view(ndarray)
        chosen = choose(indices_, choices, mode='wrap', out=store)
        assert_equal(store, array([999999, 31, 12, 999999]))

    def test_reshape(self):
        a = arange(10)
        a[0] = masked
        # Try the default
        b = a.reshape((5, 2))
        assert_equal(b.shape, (5, 2))
        self.assertTrue(b.flags['C'])
        # Try w/ arguments as list instead of tuple
        b = a.reshape(5, 2)
        assert_equal(b.shape, (5, 2))
        self.assertTrue(b.flags['C'])
        # Try w/ order
        b = a.reshape((5, 2), order='F')
        assert_equal(b.shape, (5, 2))
        self.assertTrue(b.flags['F'])
        # Try w/ order
        b = a.reshape(5, 2, order='F')
        assert_equal(b.shape, (5, 2))
        self.assertTrue(b.flags['F'])
        #
        c = np.reshape(a, (2, 5))
        self.assertTrue(isinstance(c, MaskedArray))
        assert_equal(c.shape, (2, 5))
        self.assertTrue(c[0, 0] is masked)
        self.assertTrue(c.flags['C'])

    def test_make_mask_descr(self):
        # Test make_mask_descr
        # Flexible
        ntype = [('a', np.float), ('b', np.float)]
        test = make_mask_descr(ntype)
        assert_equal(test, [('a', np.bool), ('b', np.bool)])
        # Standard w/ shape
        ntype = (np.float, 2)
        test = make_mask_descr(ntype)
        assert_equal(test, (np.bool, 2))
        # Standard standard
        ntype = np.float
        test = make_mask_descr(ntype)
        assert_equal(test, np.dtype(np.bool))
        # Nested
        ntype = [('a', np.float), ('b', [('ba', np.float), ('bb', np.float)])]
        test = make_mask_descr(ntype)
        control = np.dtype([('a', 'b1'), ('b', [('ba', 'b1'), ('bb', 'b1')])])
        assert_equal(test, control)
        # Named+ shape
        ntype = [('a', (np.float, 2))]
        test = make_mask_descr(ntype)
        assert_equal(test, np.dtype([('a', (np.bool, 2))]))
        # 2 names
        ntype = [(('A', 'a'), float)]
        test = make_mask_descr(ntype)
        assert_equal(test, np.dtype([(('A', 'a'), bool)]))

    def test_make_mask(self):
        # Test make_mask
        # w/ a list as an input
        mask = [0, 1]
        test = make_mask(mask)
        assert_equal(test.dtype, MaskType)
        assert_equal(test, [0, 1])
        # w/ a ndarray as an input
        mask = np.array([0, 1], dtype=np.bool)
        test = make_mask(mask)
        assert_equal(test.dtype, MaskType)
        assert_equal(test, [0, 1])
        # w/ a flexible-type ndarray as an input - use default
        mdtype = [('a', np.bool), ('b', np.bool)]
        mask = np.array([(0, 0), (0, 1)], dtype=mdtype)
        test = make_mask(mask)
        assert_equal(test.dtype, MaskType)
        assert_equal(test, [1, 1])
        # w/ a flexible-type ndarray as an input - use input dtype
        mdtype = [('a', np.bool), ('b', np.bool)]
        mask = np.array([(0, 0), (0, 1)], dtype=mdtype)
        test = make_mask(mask, dtype=mask.dtype)
        assert_equal(test.dtype, mdtype)
        assert_equal(test, mask)
        # w/ a flexible-type ndarray as an input - use input dtype
        mdtype = [('a', np.float), ('b', np.float)]
        bdtype = [('a', np.bool), ('b', np.bool)]
        mask = np.array([(0, 0), (0, 1)], dtype=mdtype)
        test = make_mask(mask, dtype=mask.dtype)
        assert_equal(test.dtype, bdtype)
        assert_equal(test, np.array([(0, 0), (0, 1)], dtype=bdtype))

    def test_mask_or(self):
        # Initialize
        mtype = [('a', np.bool), ('b', np.bool)]
        mask = np.array([(0, 0), (0, 1), (1, 0), (0, 0)], dtype=mtype)
        # Test using nomask as input
        test = mask_or(mask, nomask)
        assert_equal(test, mask)
        test = mask_or(nomask, mask)
        assert_equal(test, mask)
        # Using False as input
        test = mask_or(mask, False)
        assert_equal(test, mask)
        # Using True as input. Won't work, but keep it for the kicks
        # test = mask_or(mask, True)
        # control = np.array([(1, 1), (1, 1), (1, 1), (1, 1)], dtype=mtype)
        # assert_equal(test, control)
        # Using another array w / the same dtype
        other = np.array([(0, 1), (0, 1), (0, 1), (0, 1)], dtype=mtype)
        test = mask_or(mask, other)
        control = np.array([(0, 1), (0, 1), (1, 1), (0, 1)], dtype=mtype)
        assert_equal(test, control)
        # Using another array w / a different dtype
        othertype = [('A', np.bool), ('B', np.bool)]
        other = np.array([(0, 1), (0, 1), (0, 1), (0, 1)], dtype=othertype)
        try:
            test = mask_or(mask, other)
        except ValueError:
            pass
        # Using nested arrays
        dtype = [('a', np.bool), ('b', [('ba', np.bool), ('bb', np.bool)])]
        amask = np.array([(0, (1, 0)), (0, (1, 0))], dtype=dtype)
        bmask = np.array([(1, (0, 1)), (0, (0, 0))], dtype=dtype)
        cntrl = np.array([(1, (1, 1)), (0, (1, 0))], dtype=dtype)
        assert_equal(mask_or(amask, bmask), cntrl)

    def test_flatten_mask(self):
        # Tests flatten mask
        # Standarad dtype
        mask = np.array([0, 0, 1], dtype=np.bool)
        assert_equal(flatten_mask(mask), mask)
        # Flexible dtype
        mask = np.array([(0, 0), (0, 1)], dtype=[('a', bool), ('b', bool)])
        test = flatten_mask(mask)
        control = np.array([0, 0, 0, 1], dtype=bool)
        assert_equal(test, control)

        mdtype = [('a', bool), ('b', [('ba', bool), ('bb', bool)])]
        data = [(0, (0, 0)), (0, (0, 1))]
        mask = np.array(data, dtype=mdtype)
        test = flatten_mask(mask)
        control = np.array([0, 0, 0, 0, 0, 1], dtype=bool)
        assert_equal(test, control)

    def test_on_ndarray(self):
        # Test functions on ndarrays
        a = np.array([1, 2, 3, 4])
        m = array(a, mask=False)
        test = anom(a)
        assert_equal(test, m.anom())
        test = reshape(a, (2, 2))
        assert_equal(test, m.reshape(2, 2))

    def test_compress(self):
        # Test compress function on ndarray and masked array
        # Address Github #2495.
        arr = np.arange(8)
        arr.shape = 4, 2
        cond = np.array([True, False, True, True])
        control = arr[[0, 2, 3]]
        test = np.ma.compress(cond, arr, axis=0)
        assert_equal(test, control)
        marr = np.ma.array(arr)
        test = np.ma.compress(cond, marr, axis=0)
        assert_equal(test, control)

    def test_compressed(self):
        # Test ma.compressed function.
        # Address gh-4026
        a = np.ma.array([1, 2])
        test = np.ma.compressed(a)
        assert_(type(test) is np.ndarray)
        # Test case when input data is ndarray subclass
        class A(np.ndarray):
            pass
        a = np.ma.array(A(shape=0))
        test = np.ma.compressed(a)
        assert_(type(test) is A)
        # Test that compress flattens
        test = np.ma.compressed([[1],[2]])
        assert_equal(test.ndim, 1)
        test = np.ma.compressed([[[[[1]]]]])
        assert_equal(test.ndim, 1)
        # Test case when input is MaskedArray subclass
        class M(MaskedArray):
            pass
        test = np.ma.compressed(M(shape=(0,1,2)))
        assert_equal(test.ndim, 1)
        # with .compessed() overriden
        class M(MaskedArray):
            def compressed(self):
                return 42
        test = np.ma.compressed(M(shape=(0,1,2)))
        assert_equal(test, 42)

#------------------------------------------------------------------------------
class TestMaskedFields(TestCase):
    #
    def setUp(self):
        ilist = [1, 2, 3, 4, 5]
        flist = [1.1, 2.2, 3.3, 4.4, 5.5]
        slist = ['one', 'two', 'three', 'four', 'five']
        ddtype = [('a', int), ('b', float), ('c', '|S8')]
        mdtype = [('a', bool), ('b', bool), ('c', bool)]
        mask = [0, 1, 0, 0, 1]
        base = array(list(zip(ilist, flist, slist)), mask=mask, dtype=ddtype)
        self.data = dict(base=base, mask=mask, ddtype=ddtype, mdtype=mdtype)

    def test_set_records_masks(self):
        base = self.data['base']
        mdtype = self.data['mdtype']
        # Set w/ nomask or masked
        base.mask = nomask
        assert_equal_records(base._mask, np.zeros(base.shape, dtype=mdtype))
        base.mask = masked
        assert_equal_records(base._mask, np.ones(base.shape, dtype=mdtype))
        # Set w/ simple boolean
        base.mask = False
        assert_equal_records(base._mask, np.zeros(base.shape, dtype=mdtype))
        base.mask = True
        assert_equal_records(base._mask, np.ones(base.shape, dtype=mdtype))
        # Set w/ list
        base.mask = [0, 0, 0, 1, 1]
        assert_equal_records(base._mask,
                             np.array([(x, x, x) for x in [0, 0, 0, 1, 1]],
                                      dtype=mdtype))

    def test_set_record_element(self):
        # Check setting an element of a record)
        base = self.data['base']
        (base_a, base_b, base_c) = (base['a'], base['b'], base['c'])
        base[0] = (pi, pi, 'pi')

        assert_equal(base_a.dtype, int)
        assert_equal(base_a._data, [3, 2, 3, 4, 5])

        assert_equal(base_b.dtype, float)
        assert_equal(base_b._data, [pi, 2.2, 3.3, 4.4, 5.5])

        assert_equal(base_c.dtype, '|S8')
        assert_equal(base_c._data,
                     asbytes_nested(['pi', 'two', 'three', 'four', 'five']))

    def test_set_record_slice(self):
        base = self.data['base']
        (base_a, base_b, base_c) = (base['a'], base['b'], base['c'])
        base[:3] = (pi, pi, 'pi')

        assert_equal(base_a.dtype, int)
        assert_equal(base_a._data, [3, 3, 3, 4, 5])

        assert_equal(base_b.dtype, float)
        assert_equal(base_b._data, [pi, pi, pi, 4.4, 5.5])

        assert_equal(base_c.dtype, '|S8')
        assert_equal(base_c._data,
                     asbytes_nested(['pi', 'pi', 'pi', 'four', 'five']))

    def test_mask_element(self):
        "Check record access"
        base = self.data['base']
        (base_a, base_b, base_c) = (base['a'], base['b'], base['c'])
        base[0] = masked
        #
        for n in ('a', 'b', 'c'):
            assert_equal(base[n].mask, [1, 1, 0, 0, 1])
            assert_equal(base[n]._data, base._data[n])

    def test_getmaskarray(self):
        # Test getmaskarray on flexible dtype
        ndtype = [('a', int), ('b', float)]
        test = empty(3, dtype=ndtype)
        assert_equal(getmaskarray(test),
                     np.array([(0, 0), (0, 0), (0, 0)],
                              dtype=[('a', '|b1'), ('b', '|b1')]))
        test[:] = masked
        assert_equal(getmaskarray(test),
                     np.array([(1, 1), (1, 1), (1, 1)],
                              dtype=[('a', '|b1'), ('b', '|b1')]))

    def test_view(self):
        # Test view w/ flexible dtype
        iterator = list(zip(np.arange(10), np.random.rand(10)))
        data = np.array(iterator)
        a = array(iterator, dtype=[('a', float), ('b', float)])
        a.mask[0] = (1, 0)
        controlmask = np.array([1] + 19 * [0], dtype=bool)
        # Transform globally to simple dtype
        test = a.view(float)
        assert_equal(test, data.ravel())
        assert_equal(test.mask, controlmask)
        # Transform globally to dty
        test = a.view((float, 2))
        assert_equal(test, data)
        assert_equal(test.mask, controlmask.reshape(-1, 2))
        #
        test = a.view((float, 2), np.matrix)
        assert_equal(test, data)
        self.assertTrue(isinstance(test, np.matrix))

    def test_getitem(self):
        ndtype = [('a', float), ('b', float)]
        a = array(list(zip(np.random.rand(10), np.arange(10))), dtype=ndtype)
        a.mask = np.array(list(zip([0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
                                   [1, 0, 0, 0, 0, 0, 0, 0, 1, 0])),
                          dtype=[('a', bool), ('b', bool)])
        # No mask
        self.assertTrue(isinstance(a[1], MaskedArray))
        # One element masked
        self.assertTrue(isinstance(a[0], MaskedArray))
        assert_equal_records(a[0]._data, a._data[0])
        assert_equal_records(a[0]._mask, a._mask[0])
        # All element masked
        self.assertTrue(isinstance(a[-2], MaskedArray))
        assert_equal_records(a[-2]._data, a._data[-2])
        assert_equal_records(a[-2]._mask, a._mask[-2])

    def test_setitem(self):
        # Issue 4866: check that one can set individual items in [record][col]
        # and [col][record] order
        ndtype = np.dtype([('a', float), ('b', int)])
        ma = np.ma.MaskedArray([(1.0, 1), (2.0, 2)], dtype=ndtype)
        ma['a'][1] = 3.0
        assert_equal(ma['a'], np.array([1.0, 3.0]))
        ma[1]['a'] = 4.0
        assert_equal(ma['a'], np.array([1.0, 4.0]))
        # Issue 2403
        mdtype = np.dtype([('a', bool), ('b', bool)])
        # soft mask
        control = np.array([(False, True), (True, True)], dtype=mdtype)
        a = np.ma.masked_all((2,), dtype=ndtype)
        a['a'][0] = 2
        assert_equal(a.mask, control)
        a = np.ma.masked_all((2,), dtype=ndtype)
        a[0]['a'] = 2
        assert_equal(a.mask, control)
        # hard mask
        control = np.array([(True, True), (True, True)], dtype=mdtype)
        a = np.ma.masked_all((2,), dtype=ndtype)
        a.harden_mask()
        a['a'][0] = 2
        assert_equal(a.mask, control)
        a = np.ma.masked_all((2,), dtype=ndtype)
        a.harden_mask()
        a[0]['a'] = 2
        assert_equal(a.mask, control)

    def test_element_len(self):
        # check that len() works for mvoid (Github issue #576)
        for rec in self.data['base']:
            assert_equal(len(rec), len(self.data['ddtype']))


#------------------------------------------------------------------------------
class TestMaskedView(TestCase):
    #
    def setUp(self):
        iterator = list(zip(np.arange(10), np.random.rand(10)))
        data = np.array(iterator)
        a = array(iterator, dtype=[('a', float), ('b', float)])
        a.mask[0] = (1, 0)
        controlmask = np.array([1] + 19 * [0], dtype=bool)
        self.data = (data, a, controlmask)

    def test_view_to_nothing(self):
        (data, a, controlmask) = self.data
        test = a.view()
        self.assertTrue(isinstance(test, MaskedArray))
        assert_equal(test._data, a._data)
        assert_equal(test._mask, a._mask)

    def test_view_to_type(self):
        (data, a, controlmask) = self.data
        test = a.view(np.ndarray)
        self.assertTrue(not isinstance(test, MaskedArray))
        assert_equal(test, a._data)
        assert_equal_records(test, data.view(a.dtype).squeeze())

    def test_view_to_simple_dtype(self):
        (data, a, controlmask) = self.data
        # View globally
        test = a.view(float)
        self.assertTrue(isinstance(test, MaskedArray))
        assert_equal(test, data.ravel())
        assert_equal(test.mask, controlmask)

    def test_view_to_flexible_dtype(self):
        (data, a, controlmask) = self.data
        #
        test = a.view([('A', float), ('B', float)])
        assert_equal(test.mask.dtype.names, ('A', 'B'))
        assert_equal(test['A'], a['a'])
        assert_equal(test['B'], a['b'])
        #
        test = a[0].view([('A', float), ('B', float)])
        self.assertTrue(isinstance(test, MaskedArray))
        assert_equal(test.mask.dtype.names, ('A', 'B'))
        assert_equal(test['A'], a['a'][0])
        assert_equal(test['B'], a['b'][0])
        #
        test = a[-1].view([('A', float), ('B', float)])
        self.assertTrue(isinstance(test, MaskedArray))
        assert_equal(test.dtype.names, ('A', 'B'))
        assert_equal(test['A'], a['a'][-1])
        assert_equal(test['B'], a['b'][-1])

    def test_view_to_subdtype(self):
        (data, a, controlmask) = self.data
        # View globally
        test = a.view((float, 2))
        self.assertTrue(isinstance(test, MaskedArray))
        assert_equal(test, data)
        assert_equal(test.mask, controlmask.reshape(-1, 2))
        # View on 1 masked element
        test = a[0].view((float, 2))
        self.assertTrue(isinstance(test, MaskedArray))
        assert_equal(test, data[0])
        assert_equal(test.mask, (1, 0))
        # View on 1 unmasked element
        test = a[-1].view((float, 2))
        self.assertTrue(isinstance(test, MaskedArray))
        assert_equal(test, data[-1])

    def test_view_to_dtype_and_type(self):
        (data, a, controlmask) = self.data
        #
        test = a.view((float, 2), np.matrix)
        assert_equal(test, data)
        self.assertTrue(isinstance(test, np.matrix))
        self.assertTrue(not isinstance(test, MaskedArray))


def test_masked_array():
    a = np.ma.array([0, 1, 2, 3], mask=[0, 0, 1, 0])
    assert_equal(np.argwhere(a), [[1], [3]])

def test_append_masked_array():
    a = np.ma.masked_equal([1,2,3], value=2)
    b = np.ma.masked_equal([4,3,2], value=2)

    result = np.ma.append(a, b)
    expected_data = [1, 2, 3, 4, 3, 2]
    expected_mask = [False, True, False, False, False, True]
    assert_array_equal(result.data, expected_data)
    assert_array_equal(result.mask, expected_mask)

    a = np.ma.masked_all((2,2))
    b = np.ma.ones((3,1))

    result = np.ma.append(a, b)
    expected_data = [1] * 3
    expected_mask = [True] * 4 + [False] * 3
    assert_array_equal(result.data[-3], expected_data)
    assert_array_equal(result.mask, expected_mask)

    result = np.ma.append(a, b, axis=None)
    assert_array_equal(result.data[-3], expected_data)
    assert_array_equal(result.mask, expected_mask)


def test_append_masked_array_along_axis():
    a = np.ma.masked_equal([1,2,3], value=2)
    b = np.ma.masked_values([[4, 5, 6], [7, 8, 9]], 7)

    # When `axis` is specified, `values` must have the correct shape.
    assert_raises(ValueError, np.ma.append, a, b, axis=0)

    result = np.ma.append(a[np.newaxis,:], b, axis=0)
    expected = np.ma.arange(1, 10)
    expected[[1, 6]] = np.ma.masked
    expected = expected.reshape((3,3))
    assert_array_equal(result.data, expected.data)
    assert_array_equal(result.mask, expected.mask)


###############################################################################
if __name__ == "__main__":
    run_module_suite()
