from __future__ import division, absolute_import, print_function

import warnings

import numpy as np
from numpy.testing import (
    run_module_suite, TestCase, assert_, assert_equal, assert_array_equal,
    assert_almost_equal, assert_array_almost_equal, assert_raises,
    assert_allclose, assert_array_max_ulp, assert_warns, assert_raises_regex
    )
from numpy.random import rand
from numpy.lib import *
from numpy.compat import long


class TestAny(TestCase):
    def test_basic(self):
        y1 = [0, 0, 1, 0]
        y2 = [0, 0, 0, 0]
        y3 = [1, 0, 1, 0]
        assert_(np.any(y1))
        assert_(np.any(y3))
        assert_(not np.any(y2))

    def test_nd(self):
        y1 = [[0, 0, 0], [0, 1, 0], [1, 1, 0]]
        assert_(np.any(y1))
        assert_array_equal(np.sometrue(y1, axis=0), [1, 1, 0])
        assert_array_equal(np.sometrue(y1, axis=1), [0, 1, 1])


class TestAll(TestCase):
    def test_basic(self):
        y1 = [0, 1, 1, 0]
        y2 = [0, 0, 0, 0]
        y3 = [1, 1, 1, 1]
        assert_(not np.all(y1))
        assert_(np.all(y3))
        assert_(not np.all(y2))
        assert_(np.all(~np.array(y2)))

    def test_nd(self):
        y1 = [[0, 0, 1], [0, 1, 1], [1, 1, 1]]
        assert_(not np.all(y1))
        assert_array_equal(np.alltrue(y1, axis=0), [0, 0, 1])
        assert_array_equal(np.alltrue(y1, axis=1), [0, 0, 1])


class TestCopy(TestCase):
    def test_basic(self):
        a = np.array([[1, 2], [3, 4]])
        a_copy = np.copy(a)
        assert_array_equal(a, a_copy)
        a_copy[0, 0] = 10
        assert_equal(a[0, 0], 1)
        assert_equal(a_copy[0, 0], 10)

    def test_order(self):
        # It turns out that people rely on np.copy() preserving order by
        # default; changing this broke scikit-learn:
        #   https://github.com/scikit-learn/scikit-learn/commit/7842748cf777412c506a8c0ed28090711d3a3783
        a = np.array([[1, 2], [3, 4]])
        assert_(a.flags.c_contiguous)
        assert_(not a.flags.f_contiguous)
        a_fort = np.array([[1, 2], [3, 4]], order="F")
        assert_(not a_fort.flags.c_contiguous)
        assert_(a_fort.flags.f_contiguous)
        a_copy = np.copy(a)
        assert_(a_copy.flags.c_contiguous)
        assert_(not a_copy.flags.f_contiguous)
        a_fort_copy = np.copy(a_fort)
        assert_(not a_fort_copy.flags.c_contiguous)
        assert_(a_fort_copy.flags.f_contiguous)


class TestAverage(TestCase):
    def test_basic(self):
        y1 = np.array([1, 2, 3])
        assert_(average(y1, axis=0) == 2.)
        y2 = np.array([1., 2., 3.])
        assert_(average(y2, axis=0) == 2.)
        y3 = [0., 0., 0.]
        assert_(average(y3, axis=0) == 0.)

        y4 = np.ones((4, 4))
        y4[0, 1] = 0
        y4[1, 0] = 2
        assert_almost_equal(y4.mean(0), average(y4, 0))
        assert_almost_equal(y4.mean(1), average(y4, 1))

        y5 = rand(5, 5)
        assert_almost_equal(y5.mean(0), average(y5, 0))
        assert_almost_equal(y5.mean(1), average(y5, 1))

        y6 = np.matrix(rand(5, 5))
        assert_array_equal(y6.mean(0), average(y6, 0))

    def test_weights(self):
        y = np.arange(10)
        w = np.arange(10)
        actual = average(y, weights=w)
        desired = (np.arange(10) ** 2).sum()*1. / np.arange(10).sum()
        assert_almost_equal(actual, desired)

        y1 = np.array([[1, 2, 3], [4, 5, 6]])
        w0 = [1, 2]
        actual = average(y1, weights=w0, axis=0)
        desired = np.array([3., 4., 5.])
        assert_almost_equal(actual, desired)

        w1 = [0, 0, 1]
        actual = average(y1, weights=w1, axis=1)
        desired = np.array([3., 6.])
        assert_almost_equal(actual, desired)

        # This should raise an error. Can we test for that ?
        # assert_equal(average(y1, weights=w1), 9./2.)

        # 2D Case
        w2 = [[0, 0, 1], [0, 0, 2]]
        desired = np.array([3., 6.])
        assert_array_equal(average(y1, weights=w2, axis=1), desired)
        assert_equal(average(y1, weights=w2), 5.)

    def test_returned(self):
        y = np.array([[1, 2, 3], [4, 5, 6]])

        # No weights
        avg, scl = average(y, returned=True)
        assert_equal(scl, 6.)

        avg, scl = average(y, 0, returned=True)
        assert_array_equal(scl, np.array([2., 2., 2.]))

        avg, scl = average(y, 1, returned=True)
        assert_array_equal(scl, np.array([3., 3.]))

        # With weights
        w0 = [1, 2]
        avg, scl = average(y, weights=w0, axis=0, returned=True)
        assert_array_equal(scl, np.array([3., 3., 3.]))

        w1 = [1, 2, 3]
        avg, scl = average(y, weights=w1, axis=1, returned=True)
        assert_array_equal(scl, np.array([6., 6.]))

        w2 = [[0, 0, 1], [1, 2, 3]]
        avg, scl = average(y, weights=w2, axis=1, returned=True)
        assert_array_equal(scl, np.array([1., 6.]))


class TestSelect(TestCase):
    choices = [np.array([1, 2, 3]),
               np.array([4, 5, 6]),
               np.array([7, 8, 9])]
    conditions = [np.array([False, False, False]),
                  np.array([False, True, False]),
                  np.array([False, False, True])]

    def _select(self, cond, values, default=0):
        output = []
        for m in range(len(cond)):
            output += [V[m] for V, C in zip(values, cond) if C[m]] or [default]
        return output

    def test_basic(self):
        choices = self.choices
        conditions = self.conditions
        assert_array_equal(select(conditions, choices, default=15),
                           self._select(conditions, choices, default=15))

        assert_equal(len(choices), 3)
        assert_equal(len(conditions), 3)

    def test_broadcasting(self):
        conditions = [np.array(True), np.array([False, True, False])]
        choices = [1, np.arange(12).reshape(4, 3)]
        assert_array_equal(select(conditions, choices), np.ones((4, 3)))
        # default can broadcast too:
        assert_equal(select([True], [0], default=[0]).shape, (1,))

    def test_return_dtype(self):
        assert_equal(select(self.conditions, self.choices, 1j).dtype,
                     np.complex_)
        # But the conditions need to be stronger then the scalar default
        # if it is scalar.
        choices = [choice.astype(np.int8) for choice in self.choices]
        assert_equal(select(self.conditions, choices).dtype, np.int8)

        d = np.array([1, 2, 3, np.nan, 5, 7])
        m = np.isnan(d)
        assert_equal(select([m], [d]), [0, 0, 0, np.nan, 0, 0])

    def test_deprecated_empty(self):
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            assert_equal(select([], [], 3j), 3j)

        with warnings.catch_warnings():
            warnings.simplefilter("always")
            assert_warns(DeprecationWarning, select, [], [])
            warnings.simplefilter("error")
            assert_raises(DeprecationWarning, select, [], [])

    def test_non_bool_deprecation(self):
        choices = self.choices
        conditions = self.conditions[:]
        with warnings.catch_warnings():
            warnings.filterwarnings("always")
            conditions[0] = conditions[0].astype(np.int_)
            assert_warns(DeprecationWarning, select, conditions, choices)
            conditions[0] = conditions[0].astype(np.uint8)
            assert_warns(DeprecationWarning, select, conditions, choices)
            warnings.filterwarnings("error")
            assert_raises(DeprecationWarning, select, conditions, choices)

    def test_many_arguments(self):
        # This used to be limited by NPY_MAXARGS == 32
        conditions = [np.array([False])] * 100
        choices = [np.array([1])] * 100
        select(conditions, choices)


class TestInsert(TestCase):
    def test_basic(self):
        a = [1, 2, 3]
        assert_equal(insert(a, 0, 1), [1, 1, 2, 3])
        assert_equal(insert(a, 3, 1), [1, 2, 3, 1])
        assert_equal(insert(a, [1, 1, 1], [1, 2, 3]), [1, 1, 2, 3, 2, 3])
        assert_equal(insert(a, 1, [1, 2, 3]), [1, 1, 2, 3, 2, 3])
        assert_equal(insert(a, [1, -1, 3], 9), [1, 9, 2, 9, 3, 9])
        assert_equal(insert(a, slice(-1, None, -1), 9), [9, 1, 9, 2, 9, 3])
        assert_equal(insert(a, [-1, 1, 3], [7, 8, 9]), [1, 8, 2, 7, 3, 9])
        b = np.array([0, 1], dtype=np.float64)
        assert_equal(insert(b, 0, b[0]), [0., 0., 1.])
        assert_equal(insert(b, [], []), b)
        # Bools will be treated differently in the future:
        #assert_equal(insert(a, np.array([True]*4), 9), [9,1,9,2,9,3,9])
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', '', FutureWarning)
            assert_equal(
                insert(a, np.array([True]*4), 9), [1, 9, 9, 9, 9, 2, 3])
            assert_(w[0].category is FutureWarning)

    def test_multidim(self):
        a = [[1, 1, 1]]
        r = [[2, 2, 2],
             [1, 1, 1]]
        assert_equal(insert(a, 0, [1]), [1, 1, 1, 1])
        assert_equal(insert(a, 0, [2, 2, 2], axis=0), r)
        assert_equal(insert(a, 0, 2, axis=0), r)
        assert_equal(insert(a, 2, 2, axis=1), [[1, 1, 2, 1]])

        a = np.array([[1, 1], [2, 2], [3, 3]])
        b = np.arange(1, 4).repeat(3).reshape(3, 3)
        c = np.concatenate(
            (a[:, 0:1], np.arange(1, 4).repeat(3).reshape(3, 3).T,
             a[:, 1:2]), axis=1)
        assert_equal(insert(a, [1], [[1], [2], [3]], axis=1), b)
        assert_equal(insert(a, [1], [1, 2, 3], axis=1), c)
        # scalars behave differently, in this case exactly opposite:
        assert_equal(insert(a, 1, [1, 2, 3], axis=1), b)
        assert_equal(insert(a, 1, [[1], [2], [3]], axis=1), c)

        a = np.arange(4).reshape(2, 2)
        assert_equal(insert(a[:, :1], 1, a[:, 1], axis=1), a)
        assert_equal(insert(a[:1, :], 1, a[1, :], axis=0), a)

        # negative axis value
        a = np.arange(24).reshape((2, 3, 4))
        assert_equal(insert(a, 1, a[:, :, 3], axis=-1),
                     insert(a, 1, a[:, :, 3], axis=2))
        assert_equal(insert(a, 1, a[:, 2, :], axis=-2),
                     insert(a, 1, a[:, 2, :], axis=1))

        # invalid axis value
        assert_raises(IndexError, insert, a, 1, a[:, 2, :], axis=3)
        assert_raises(IndexError, insert, a, 1, a[:, 2, :], axis=-4)

        # negative axis value
        a = np.arange(24).reshape((2,3,4))
        assert_equal(insert(a, 1, a[:,:,3], axis=-1),
                     insert(a, 1, a[:,:,3], axis=2))
        assert_equal(insert(a, 1, a[:,2,:], axis=-2),
                     insert(a, 1, a[:,2,:], axis=1))

    def test_0d(self):
        # This is an error in the future
        a = np.array(1)
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', '', DeprecationWarning)
            assert_equal(insert(a, [], 2, axis=0), np.array(2))
            assert_(w[0].category is DeprecationWarning)

    def test_subclass(self):
        class SubClass(np.ndarray):
            pass
        a = np.arange(10).view(SubClass)
        assert_(isinstance(np.insert(a, 0, [0]), SubClass))
        assert_(isinstance(np.insert(a, [], []), SubClass))
        assert_(isinstance(np.insert(a, [0, 1], [1, 2]), SubClass))
        assert_(isinstance(np.insert(a, slice(1, 2), [1, 2]), SubClass))
        assert_(isinstance(np.insert(a, slice(1, -2, -1), []), SubClass))
        # This is an error in the future:
        a = np.array(1).view(SubClass)
        assert_(isinstance(np.insert(a, 0, [0]), SubClass))

    def test_index_array_copied(self):
        x = np.array([1, 1, 1])
        np.insert([0, 1, 2], x, [3, 4, 5])
        assert_equal(x, np.array([1, 1, 1]))

    def test_structured_array(self):
        a = np.array([(1, 'a'), (2, 'b'), (3, 'c')],
                     dtype=[('foo', 'i'), ('bar', 'a1')])
        val = (4, 'd')
        b = np.insert(a, 0, val)
        assert_array_equal(b[0], np.array(val, dtype=b.dtype))
        val = [(4, 'd')] * 2
        b = np.insert(a, [0, 2], val)
        assert_array_equal(b[[0, 3]], np.array(val, dtype=b.dtype))


class TestAmax(TestCase):
    def test_basic(self):
        a = [3, 4, 5, 10, -3, -5, 6.0]
        assert_equal(np.amax(a), 10.0)
        b = [[3, 6.0, 9.0],
             [4, 10.0, 5.0],
             [8, 3.0, 2.0]]
        assert_equal(np.amax(b, axis=0), [8.0, 10.0, 9.0])
        assert_equal(np.amax(b, axis=1), [9.0, 10.0, 8.0])


class TestAmin(TestCase):
    def test_basic(self):
        a = [3, 4, 5, 10, -3, -5, 6.0]
        assert_equal(np.amin(a), -5.0)
        b = [[3, 6.0, 9.0],
             [4, 10.0, 5.0],
             [8, 3.0, 2.0]]
        assert_equal(np.amin(b, axis=0), [3.0, 3.0, 2.0])
        assert_equal(np.amin(b, axis=1), [3.0, 4.0, 2.0])


class TestPtp(TestCase):
    def test_basic(self):
        a = [3, 4, 5, 10, -3, -5, 6.0]
        assert_equal(np.ptp(a, axis=0), 15.0)
        b = [[3, 6.0, 9.0],
             [4, 10.0, 5.0],
             [8, 3.0, 2.0]]
        assert_equal(np.ptp(b, axis=0), [5.0, 7.0, 7.0])
        assert_equal(np.ptp(b, axis=-1), [6.0, 6.0, 6.0])


class TestCumsum(TestCase):
    def test_basic(self):
        ba = [1, 2, 10, 11, 6, 5, 4]
        ba2 = [[1, 2, 3, 4], [5, 6, 7, 9], [10, 3, 4, 5]]
        for ctype in [np.int8, np.uint8, np.int16, np.uint16, np.int32,
                      np.uint32, np.float32, np.float64, np.complex64, np.complex128]:
            a = np.array(ba, ctype)
            a2 = np.array(ba2, ctype)

            tgt = np.array([1, 3, 13, 24, 30, 35, 39], ctype)
            assert_array_equal(np.cumsum(a, axis=0), tgt)

            tgt = np.array(
                [[1, 2, 3, 4], [6, 8, 10, 13], [16, 11, 14, 18]], ctype)
            assert_array_equal(np.cumsum(a2, axis=0), tgt)

            tgt = np.array(
                [[1, 3, 6, 10], [5, 11, 18, 27], [10, 13, 17, 22]], ctype)
            assert_array_equal(np.cumsum(a2, axis=1), tgt)


class TestProd(TestCase):
    def test_basic(self):
        ba = [1, 2, 10, 11, 6, 5, 4]
        ba2 = [[1, 2, 3, 4], [5, 6, 7, 9], [10, 3, 4, 5]]
        for ctype in [np.int16, np.uint16, np.int32, np.uint32,
                      np.float32, np.float64, np.complex64, np.complex128]:
            a = np.array(ba, ctype)
            a2 = np.array(ba2, ctype)
            if ctype in ['1', 'b']:
                self.assertRaises(ArithmeticError, prod, a)
                self.assertRaises(ArithmeticError, prod, a2, 1)
                self.assertRaises(ArithmeticError, prod, a)
            else:
                assert_equal(np.prod(a, axis=0), 26400)
                assert_array_equal(np.prod(a2, axis=0),
                                   np.array([50, 36, 84, 180], ctype))
                assert_array_equal(np.prod(a2, axis=-1),
                                   np.array([24, 1890, 600], ctype))


class TestCumprod(TestCase):
    def test_basic(self):
        ba = [1, 2, 10, 11, 6, 5, 4]
        ba2 = [[1, 2, 3, 4], [5, 6, 7, 9], [10, 3, 4, 5]]
        for ctype in [np.int16, np.uint16, np.int32, np.uint32,
                      np.float32, np.float64, np.complex64, np.complex128]:
            a = np.array(ba, ctype)
            a2 = np.array(ba2, ctype)
            if ctype in ['1', 'b']:
                self.assertRaises(ArithmeticError, cumprod, a)
                self.assertRaises(ArithmeticError, cumprod, a2, 1)
                self.assertRaises(ArithmeticError, cumprod, a)
            else:
                assert_array_equal(np.cumprod(a, axis=-1),
                                   np.array([1, 2, 20, 220,
                                             1320, 6600, 26400], ctype))
                assert_array_equal(np.cumprod(a2, axis=0),
                                   np.array([[1, 2, 3, 4],
                                             [5, 12, 21, 36],
                                             [50, 36, 84, 180]], ctype))
                assert_array_equal(np.cumprod(a2, axis=-1),
                                   np.array([[1, 2, 6, 24],
                                             [5, 30, 210, 1890],
                                             [10, 30, 120, 600]], ctype))


class TestDiff(TestCase):
    def test_basic(self):
        x = [1, 4, 6, 7, 12]
        out = np.array([3, 2, 1, 5])
        out2 = np.array([-1, -1, 4])
        out3 = np.array([0, 5])
        assert_array_equal(diff(x), out)
        assert_array_equal(diff(x, n=2), out2)
        assert_array_equal(diff(x, n=3), out3)

    def test_nd(self):
        x = 20 * rand(10, 20, 30)
        out1 = x[:, :, 1:] - x[:, :, :-1]
        out2 = out1[:, :, 1:] - out1[:, :, :-1]
        out3 = x[1:, :, :] - x[:-1, :, :]
        out4 = out3[1:, :, :] - out3[:-1, :, :]
        assert_array_equal(diff(x), out1)
        assert_array_equal(diff(x, n=2), out2)
        assert_array_equal(diff(x, axis=0), out3)
        assert_array_equal(diff(x, n=2, axis=0), out4)


class TestDelete(TestCase):
    def setUp(self):
        self.a = np.arange(5)
        self.nd_a = np.arange(5).repeat(2).reshape(1, 5, 2)

    def _check_inverse_of_slicing(self, indices):
        a_del = delete(self.a, indices)
        nd_a_del = delete(self.nd_a, indices, axis=1)
        msg = 'Delete failed for obj: %r' % indices
        # NOTE: The cast should be removed after warning phase for bools
        if not isinstance(indices, (slice, int, long, np.integer)):
            indices = np.asarray(indices, dtype=np.intp)
            indices = indices[(indices >= 0) & (indices < 5)]
        assert_array_equal(setxor1d(a_del, self.a[indices, ]), self.a,
                           err_msg=msg)
        xor = setxor1d(nd_a_del[0, :, 0], self.nd_a[0, indices, 0])
        assert_array_equal(xor, self.nd_a[0, :, 0], err_msg=msg)

    def test_slices(self):
        lims = [-6, -2, 0, 1, 2, 4, 5]
        steps = [-3, -1, 1, 3]
        for start in lims:
            for stop in lims:
                for step in steps:
                    s = slice(start, stop, step)
                    self._check_inverse_of_slicing(s)

    def test_fancy(self):
        # Deprecation/FutureWarning tests should be kept after change.
        self._check_inverse_of_slicing(np.array([[0, 1], [2, 1]]))
        with warnings.catch_warnings():
            warnings.filterwarnings('error', category=DeprecationWarning)
            assert_raises(DeprecationWarning, delete, self.a, [100])
            assert_raises(DeprecationWarning, delete, self.a, [-100])
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', category=FutureWarning)
            self._check_inverse_of_slicing([0, -1, 2, 2])
            obj = np.array([True, False, False], dtype=bool)
            self._check_inverse_of_slicing(obj)
            assert_(w[0].category is FutureWarning)
            assert_(w[1].category is FutureWarning)

    def test_single(self):
        self._check_inverse_of_slicing(0)
        self._check_inverse_of_slicing(-4)

    def test_0d(self):
        a = np.array(1)
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', '', DeprecationWarning)
            assert_equal(delete(a, [], axis=0), a)
            assert_(w[0].category is DeprecationWarning)

    def test_subclass(self):
        class SubClass(np.ndarray):
            pass
        a = self.a.view(SubClass)
        assert_(isinstance(delete(a, 0), SubClass))
        assert_(isinstance(delete(a, []), SubClass))
        assert_(isinstance(delete(a, [0, 1]), SubClass))
        assert_(isinstance(delete(a, slice(1, 2)), SubClass))
        assert_(isinstance(delete(a, slice(1, -2)), SubClass))


class TestGradient(TestCase):
    def test_basic(self):
        v = [[1, 1], [3, 4]]
        x = np.array(v)
        dx = [np.array([[2., 3.], [2., 3.]]),
              np.array([[0., 0.], [1., 1.]])]
        assert_array_equal(gradient(x), dx)
        assert_array_equal(gradient(v), dx)

    def test_badargs(self):
        # for 2D array, gradient can take 0, 1, or 2 extra args
        x = np.array([[1, 1], [3, 4]])
        assert_raises(SyntaxError, gradient, x, np.array([1., 1.]),
                      np.array([1., 1.]), np.array([1., 1.]))

    def test_masked(self):
        # Make sure that gradient supports subclasses like masked arrays
        x = np.ma.array([[1, 1], [3, 4]])
        assert_equal(type(gradient(x)[0]), type(x))

    def test_datetime64(self):
        # Make sure gradient() can handle special types like datetime64
        x = np.array(
            ['1910-08-16', '1910-08-11', '1910-08-10', '1910-08-12',
             '1910-10-12', '1910-12-12', '1912-12-12'],
            dtype='datetime64[D]')
        dx = np.array(
            [-7, -3, 0, 31, 61, 396, 1066],
            dtype='timedelta64[D]')
        assert_array_equal(gradient(x), dx)
        assert_(dx.dtype == np.dtype('timedelta64[D]'))

    def test_timedelta64(self):
        # Make sure gradient() can handle special types like timedelta64
        x = np.array(
            [-5, -3, 10, 12, 61, 321, 300],
            dtype='timedelta64[D]')
        dx = np.array(
            [-3, 7, 7, 25, 154, 119, -161],
            dtype='timedelta64[D]')
        assert_array_equal(gradient(x), dx)
        assert_(dx.dtype == np.dtype('timedelta64[D]'))

    def test_second_order_accurate(self):
        # Testing that the relative numerical error is less that 3% for
        # this example problem. This corresponds to second order
        # accurate finite differences for all interior and boundary
        # points.
        x = np.linspace(0, 1, 10)
        dx = x[1] - x[0]
        y = 2 * x ** 3 + 4 * x ** 2 + 2 * x
        analytical = 6 * x ** 2 + 8 * x + 2
        num_error = np.abs((np.gradient(y, dx) / analytical) - 1)
        assert_(np.all(num_error < 0.03) == True)


class TestAngle(TestCase):
    def test_basic(self):
        x = [1 + 3j, np.sqrt(2) / 2.0 + 1j * np.sqrt(2) / 2,
             1, 1j, -1, -1j, 1 - 3j, -1 + 3j]
        y = angle(x)
        yo = [
            np.arctan(3.0 / 1.0),
            np.arctan(1.0), 0, np.pi / 2, np.pi, -np.pi / 2.0,
            -np.arctan(3.0 / 1.0), np.pi - np.arctan(3.0 / 1.0)]
        z = angle(x, deg=1)
        zo = np.array(yo) * 180 / np.pi
        assert_array_almost_equal(y, yo, 11)
        assert_array_almost_equal(z, zo, 11)


class TestTrimZeros(TestCase):
    """ only testing for integer splits.
    """
    def test_basic(self):
        a = np.array([0, 0, 1, 2, 3, 4, 0])
        res = trim_zeros(a)
        assert_array_equal(res, np.array([1, 2, 3, 4]))

    def test_leading_skip(self):
        a = np.array([0, 0, 1, 0, 2, 3, 4, 0])
        res = trim_zeros(a)
        assert_array_equal(res, np.array([1, 0, 2, 3, 4]))

    def test_trailing_skip(self):
        a = np.array([0, 0, 1, 0, 2, 3, 0, 4, 0])
        res = trim_zeros(a)
        assert_array_equal(res, np.array([1, 0, 2, 3, 0, 4]))


class TestExtins(TestCase):
    def test_basic(self):
        a = np.array([1, 3, 2, 1, 2, 3, 3])
        b = extract(a > 1, a)
        assert_array_equal(b, [3, 2, 2, 3, 3])

    def test_place(self):
        a = np.array([1, 4, 3, 2, 5, 8, 7])
        place(a, [0, 1, 0, 1, 0, 1, 0], [2, 4, 6])
        assert_array_equal(a, [1, 2, 3, 4, 5, 6, 7])

    def test_both(self):
        a = rand(10)
        mask = a > 0.5
        ac = a.copy()
        c = extract(mask, a)
        place(a, mask, 0)
        place(a, mask, c)
        assert_array_equal(a, ac)


class TestVectorize(TestCase):
    def test_simple(self):
        def addsubtract(a, b):
            if a > b:
                return a - b
            else:
                return a + b
        f = vectorize(addsubtract)
        r = f([0, 3, 6, 9], [1, 3, 5, 7])
        assert_array_equal(r, [1, 6, 1, 2])

    def test_scalar(self):
        def addsubtract(a, b):
            if a > b:
                return a - b
            else:
                return a + b
        f = vectorize(addsubtract)
        r = f([0, 3, 6, 9], 5)
        assert_array_equal(r, [5, 8, 1, 4])

    def test_large(self):
        x = np.linspace(-3, 2, 10000)
        f = vectorize(lambda x: x)
        y = f(x)
        assert_array_equal(y, x)

    def test_ufunc(self):
        import math
        f = vectorize(math.cos)
        args = np.array([0, 0.5*np.pi, np.pi, 1.5*np.pi, 2*np.pi])
        r1 = f(args)
        r2 = np.cos(args)
        assert_array_equal(r1, r2)

    def test_keywords(self):
        import math

        def foo(a, b=1):
            return a + b
        f = vectorize(foo)
        args = np.array([1, 2, 3])
        r1 = f(args)
        r2 = np.array([2, 3, 4])
        assert_array_equal(r1, r2)
        r1 = f(args, 2)
        r2 = np.array([3, 4, 5])
        assert_array_equal(r1, r2)

    def test_keywords_no_func_code(self):
        # This needs to test a function that has keywords but
        # no func_code attribute, since otherwise vectorize will
        # inspect the func_code.
        import random
        try:
            f = vectorize(random.randrange)
        except:
            raise AssertionError()

    def test_keywords2_ticket_2100(self):
        r"""Test kwarg support: enhancement ticket 2100"""
        import math

        def foo(a, b=1):
            return a + b
        f = vectorize(foo)
        args = np.array([1, 2, 3])
        r1 = f(a=args)
        r2 = np.array([2, 3, 4])
        assert_array_equal(r1, r2)
        r1 = f(b=1, a=args)
        assert_array_equal(r1, r2)
        r1 = f(args, b=2)
        r2 = np.array([3, 4, 5])
        assert_array_equal(r1, r2)

    def test_keywords3_ticket_2100(self):
        """Test excluded with mixed positional and kwargs: ticket 2100"""
        def mypolyval(x, p):
            _p = list(p)
            res = _p.pop(0)
            while _p:
                res = res*x + _p.pop(0)
            return res
        vpolyval = np.vectorize(mypolyval, excluded=['p', 1])
        ans = [3, 6]
        assert_array_equal(ans, vpolyval(x=[0, 1], p=[1, 2, 3]))
        assert_array_equal(ans, vpolyval([0, 1], p=[1, 2, 3]))
        assert_array_equal(ans, vpolyval([0, 1], [1, 2, 3]))

    def test_keywords4_ticket_2100(self):
        """Test vectorizing function with no positional args."""
        @vectorize
        def f(**kw):
            res = 1.0
            for _k in kw:
                res *= kw[_k]
            return res
        assert_array_equal(f(a=[1, 2], b=[3, 4]), [3, 8])

    def test_keywords5_ticket_2100(self):
        """Test vectorizing function with no kwargs args."""
        @vectorize
        def f(*v):
            return np.prod(v)
        assert_array_equal(f([1, 2], [3, 4]), [3, 8])

    def test_coverage1_ticket_2100(self):
        def foo():
            return 1
        f = vectorize(foo)
        assert_array_equal(f(), 1)

    def test_assigning_docstring(self):
        def foo(x):
            return x
        doc = "Provided documentation"
        f = vectorize(foo, doc=doc)
        assert_equal(f.__doc__, doc)

    def test_UnboundMethod_ticket_1156(self):
        """Regression test for issue 1156"""
        class Foo:
            b = 2

            def bar(self, a):
                return a**self.b
        assert_array_equal(vectorize(Foo().bar)(np.arange(9)),
                           np.arange(9)**2)
        assert_array_equal(vectorize(Foo.bar)(Foo(), np.arange(9)),
                           np.arange(9)**2)

    def test_execution_order_ticket_1487(self):
        """Regression test for dependence on execution order: issue 1487"""
        f1 = vectorize(lambda x: x)
        res1a = f1(np.arange(3))
        res1b = f1(np.arange(0.1, 3))
        f2 = vectorize(lambda x: x)
        res2b = f2(np.arange(0.1, 3))
        res2a = f2(np.arange(3))
        assert_equal(res1a, res2a)
        assert_equal(res1b, res2b)

    def test_string_ticket_1892(self):
        """Test vectorization over strings: issue 1892."""
        f = np.vectorize(lambda x: x)
        s = '0123456789'*10
        assert_equal(s, f(s))
        #z = f(np.array([s,s]))
        #assert_array_equal([s,s], f(s))

    def test_cache(self):
        """Ensure that vectorized func called exactly once per argument."""
        _calls = [0]

        @vectorize
        def f(x):
            _calls[0] += 1
            return x**2
        f.cache = True
        x = np.arange(5)
        assert_array_equal(f(x), x*x)
        assert_equal(_calls[0], len(x))

    def test_otypes(self):
        f = np.vectorize(lambda x: x)
        f.otypes = 'i'
        x = np.arange(5)
        assert_array_equal(f(x), x)


class TestDigitize(TestCase):
    def test_forward(self):
        x = np.arange(-6, 5)
        bins = np.arange(-5, 5)
        assert_array_equal(digitize(x, bins), np.arange(11))

    def test_reverse(self):
        x = np.arange(5, -6, -1)
        bins = np.arange(5, -5, -1)
        assert_array_equal(digitize(x, bins), np.arange(11))

    def test_random(self):
        x = rand(10)
        bin = np.linspace(x.min(), x.max(), 10)
        assert_(np.all(digitize(x, bin) != 0))

    def test_right_basic(self):
        x = [1, 5, 4, 10, 8, 11, 0]
        bins = [1, 5, 10]
        default_answer = [1, 2, 1, 3, 2, 3, 0]
        assert_array_equal(digitize(x, bins), default_answer)
        right_answer = [0, 1, 1, 2, 2, 3, 0]
        assert_array_equal(digitize(x, bins, True), right_answer)

    def test_right_open(self):
        x = np.arange(-6, 5)
        bins = np.arange(-6, 4)
        assert_array_equal(digitize(x, bins, True), np.arange(11))

    def test_right_open_reverse(self):
        x = np.arange(5, -6, -1)
        bins = np.arange(4, -6, -1)
        assert_array_equal(digitize(x, bins, True), np.arange(11))

    def test_right_open_random(self):
        x = rand(10)
        bins = np.linspace(x.min(), x.max(), 10)
        assert_(np.all(digitize(x, bins, True) != 10))

    def test_monotonic(self):
        x = [-1, 0, 1, 2]
        bins = [0, 0, 1]
        assert_array_equal(digitize(x, bins, False), [0, 2, 3, 3])
        assert_array_equal(digitize(x, bins, True), [0, 0, 2, 3])
        bins = [1, 1, 0]
        assert_array_equal(digitize(x, bins, False), [3, 2, 0, 0])
        assert_array_equal(digitize(x, bins, True), [3, 3, 2, 0])
        bins = [1, 1, 1, 1]
        assert_array_equal(digitize(x, bins, False), [0, 0, 4, 4])
        assert_array_equal(digitize(x, bins, True), [0, 0, 0, 4])
        bins = [0, 0, 1, 0]
        assert_raises(ValueError, digitize, x, bins)
        bins = [1, 1, 0, 1]
        assert_raises(ValueError, digitize, x, bins)


class TestUnwrap(TestCase):
    def test_simple(self):
                #check that unwrap removes jumps greather that 2*pi
        assert_array_equal(unwrap([1, 1 + 2 * np.pi]), [1, 1])
        #check that unwrap maintans continuity
        assert_(np.all(diff(unwrap(rand(10) * 100)) < np.pi))


class TestFilterwindows(TestCase):
    def test_hanning(self):
        #check symmetry
        w = hanning(10)
        assert_array_almost_equal(w, flipud(w), 7)
        #check known value
        assert_almost_equal(np.sum(w, axis=0), 4.500, 4)

    def test_hamming(self):
        #check symmetry
        w = hamming(10)
        assert_array_almost_equal(w, flipud(w), 7)
        #check known value
        assert_almost_equal(np.sum(w, axis=0), 4.9400, 4)

    def test_bartlett(self):
        #check symmetry
        w = bartlett(10)
        assert_array_almost_equal(w, flipud(w), 7)
        #check known value
        assert_almost_equal(np.sum(w, axis=0), 4.4444, 4)

    def test_blackman(self):
        #check symmetry
        w = blackman(10)
        assert_array_almost_equal(w, flipud(w), 7)
        #check known value
        assert_almost_equal(np.sum(w, axis=0), 3.7800, 4)


class TestTrapz(TestCase):
    def test_simple(self):
        x = np.arange(-10, 10, .1)
        r = trapz(np.exp(-.5*x**2) / np.sqrt(2*np.pi), dx=0.1)
        #check integral of normal equals 1
        assert_almost_equal(r, 1, 7)

    def test_ndim(self):
        x = np.linspace(0, 1, 3)
        y = np.linspace(0, 2, 8)
        z = np.linspace(0, 3, 13)

        wx = np.ones_like(x) * (x[1] - x[0])
        wx[0] /= 2
        wx[-1] /= 2
        wy = np.ones_like(y) * (y[1] - y[0])
        wy[0] /= 2
        wy[-1] /= 2
        wz = np.ones_like(z) * (z[1] - z[0])
        wz[0] /= 2
        wz[-1] /= 2

        q = x[:, None, None] + y[None, :, None] + z[None, None, :]

        qx = (q * wx[:, None, None]).sum(axis=0)
        qy = (q * wy[None, :, None]).sum(axis=1)
        qz = (q * wz[None, None, :]).sum(axis=2)

        # n-d `x`
        r = trapz(q, x=x[:, None, None], axis=0)
        assert_almost_equal(r, qx)
        r = trapz(q, x=y[None, :, None], axis=1)
        assert_almost_equal(r, qy)
        r = trapz(q, x=z[None, None, :], axis=2)
        assert_almost_equal(r, qz)

        # 1-d `x`
        r = trapz(q, x=x, axis=0)
        assert_almost_equal(r, qx)
        r = trapz(q, x=y, axis=1)
        assert_almost_equal(r, qy)
        r = trapz(q, x=z, axis=2)
        assert_almost_equal(r, qz)

    def test_masked(self):
        #Testing that masked arrays behave as if the function is 0 where
        #masked
        x = np.arange(5)
        y = x * x
        mask = x == 2
        ym = np.ma.array(y, mask=mask)
        r = 13.0  # sum(0.5 * (0 + 1) * 1.0 + 0.5 * (9 + 16))
        assert_almost_equal(trapz(ym, x), r)

        xm = np.ma.array(x, mask=mask)
        assert_almost_equal(trapz(ym, xm), r)

        xm = np.ma.array(x, mask=mask)
        assert_almost_equal(trapz(y, xm), r)

    def test_matrix(self):
        #Test to make sure matrices give the same answer as ndarrays
        x = np.linspace(0, 5)
        y = x * x
        r = trapz(y, x)
        mx = np.matrix(x)
        my = np.matrix(y)
        mr = trapz(my, mx)
        assert_almost_equal(mr, r)


class TestSinc(TestCase):
    def test_simple(self):
        assert_(sinc(0) == 1)
        w = sinc(np.linspace(-1, 1, 100))
        #check symmetry
        assert_array_almost_equal(w, flipud(w), 7)

    def test_array_like(self):
        x = [0, 0.5]
        y1 = sinc(np.array(x))
        y2 = sinc(list(x))
        y3 = sinc(tuple(x))
        assert_array_equal(y1, y2)
        assert_array_equal(y1, y3)


class TestHistogram(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_simple(self):
        n = 100
        v = rand(n)
        (a, b) = histogram(v)
        #check if the sum of the bins equals the number of samples
        assert_equal(np.sum(a, axis=0), n)
        #check that the bin counts are evenly spaced when the data is from a
        # linear function
        (a, b) = histogram(np.linspace(0, 10, 100))
        assert_array_equal(a, 10)

    def test_one_bin(self):
        # Ticket 632
        hist, edges = histogram([1, 2, 3, 4], [1, 2])
        assert_array_equal(hist, [2, ])
        assert_array_equal(edges, [1, 2])
        assert_raises(ValueError, histogram, [1, 2], bins=0)
        h, e = histogram([1, 2], bins=1)
        assert_equal(h, np.array([2]))
        assert_allclose(e, np.array([1., 2.]))

    def test_normed(self):
        # Check that the integral of the density equals 1.
        n = 100
        v = rand(n)
        a, b = histogram(v, normed=True)
        area = np.sum(a * diff(b))
        assert_almost_equal(area, 1)

        # Check with non-constant bin widths (buggy but backwards compatible)
        v = np.arange(10)
        bins = [0, 1, 5, 9, 10]
        a, b = histogram(v, bins, normed=True)
        area = np.sum(a * diff(b))
        assert_almost_equal(area, 1)

    def test_density(self):
        # Check that the integral of the density equals 1.
        n = 100
        v = rand(n)
        a, b = histogram(v, density=True)
        area = np.sum(a * diff(b))
        assert_almost_equal(area, 1)

        # Check with non-constant bin widths
        v = np.arange(10)
        bins = [0, 1, 3, 6, 10]
        a, b = histogram(v, bins, density=True)
        assert_array_equal(a, .1)
        assert_equal(np.sum(a*diff(b)), 1)

        # Variale bin widths are especially useful to deal with
        # infinities.
        v = np.arange(10)
        bins = [0, 1, 3, 6, np.inf]
        a, b = histogram(v, bins, density=True)
        assert_array_equal(a, [.1, .1, .1, 0.])

        # Taken from a bug report from N. Becker on the numpy-discussion
        # mailing list Aug. 6, 2010.
        counts, dmy = np.histogram(
            [1, 2, 3, 4], [0.5, 1.5, np.inf], density=True)
        assert_equal(counts, [.25, 0])

    def test_outliers(self):
        # Check that outliers are not tallied
        a = np.arange(10) + .5

        # Lower outliers
        h, b = histogram(a, range=[0, 9])
        assert_equal(h.sum(), 9)

        # Upper outliers
        h, b = histogram(a, range=[1, 10])
        assert_equal(h.sum(), 9)

        # Normalization
        h, b = histogram(a, range=[1, 9], normed=True)
        assert_almost_equal((h * diff(b)).sum(), 1, decimal=15)

        # Weights
        w = np.arange(10) + .5
        h, b = histogram(a, range=[1, 9], weights=w, normed=True)
        assert_equal((h * diff(b)).sum(), 1)

        h, b = histogram(a, bins=8, range=[1, 9], weights=w)
        assert_equal(h, w[1:-1])

    def test_type(self):
        # Check the type of the returned histogram
        a = np.arange(10) + .5
        h, b = histogram(a)
        assert_(issubdtype(h.dtype, int))

        h, b = histogram(a, normed=True)
        assert_(issubdtype(h.dtype, float))

        h, b = histogram(a, weights=np.ones(10, int))
        assert_(issubdtype(h.dtype, int))

        h, b = histogram(a, weights=np.ones(10, float))
        assert_(issubdtype(h.dtype, float))

    def test_f32_rounding(self):
        # gh-4799, check that the rounding of the edges works with float32
        x = np.array([276.318359  , -69.593948  , 21.329449], dtype=np.float32)
        y = np.array([5005.689453, 4481.327637, 6010.369629], dtype=np.float32)
        counts_hist, xedges, yedges = np.histogram2d(x, y, bins=100)
        assert_equal(counts_hist.sum(), 3.)

    def test_weights(self):
        v = rand(100)
        w = np.ones(100) * 5
        a, b = histogram(v)
        na, nb = histogram(v, normed=True)
        wa, wb = histogram(v, weights=w)
        nwa, nwb = histogram(v, weights=w, normed=True)
        assert_array_almost_equal(a * 5, wa)
        assert_array_almost_equal(na, nwa)

        # Check weights are properly applied.
        v = np.linspace(0, 10, 10)
        w = np.concatenate((np.zeros(5), np.ones(5)))
        wa, wb = histogram(v, bins=np.arange(11), weights=w)
        assert_array_almost_equal(wa, w)

        # Check with integer weights
        wa, wb = histogram([1, 2, 2, 4], bins=4, weights=[4, 3, 2, 1])
        assert_array_equal(wa, [4, 5, 0, 1])
        wa, wb = histogram(
            [1, 2, 2, 4], bins=4, weights=[4, 3, 2, 1], normed=True)
        assert_array_almost_equal(wa, np.array([4, 5, 0, 1]) / 10. / 3. * 4)

        # Check weights with non-uniform bin widths
        a, b = histogram(
            np.arange(9), [0, 1, 3, 6, 10],
            weights=[2, 1, 1, 1, 1, 1, 1, 1, 1], density=True)
        assert_almost_equal(a, [.2, .1, .1, .075])

    def test_empty(self):
        a, b = histogram([], bins=([0, 1]))
        assert_array_equal(a, np.array([0]))
        assert_array_equal(b, np.array([0, 1]))


class TestHistogramdd(TestCase):
    def test_simple(self):
        x = np.array([[-.5, .5, 1.5], [-.5, 1.5, 2.5], [-.5, 2.5, .5],
                      [.5,  .5, 1.5], [.5,  1.5, 2.5], [.5,  2.5, 2.5]])
        H, edges = histogramdd(x, (2, 3, 3),
                               range=[[-1, 1], [0, 3], [0, 3]])
        answer = np.array([[[0, 1, 0], [0, 0, 1], [1, 0, 0]],
                           [[0, 1, 0], [0, 0, 1], [0, 0, 1]]])
        assert_array_equal(H, answer)

        # Check normalization
        ed = [[-2, 0, 2], [0, 1, 2, 3], [0, 1, 2, 3]]
        H, edges = histogramdd(x, bins=ed, normed=True)
        assert_(np.all(H == answer / 12.))

        # Check that H has the correct shape.
        H, edges = histogramdd(x, (2, 3, 4),
                               range=[[-1, 1], [0, 3], [0, 4]],
                               normed=True)
        answer = np.array([[[0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0]],
                           [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 1, 0]]])
        assert_array_almost_equal(H, answer / 6., 4)
        # Check that a sequence of arrays is accepted and H has the correct
        # shape.
        z = [np.squeeze(y) for y in split(x, 3, axis=1)]
        H, edges = histogramdd(
            z, bins=(4, 3, 2), range=[[-2, 2], [0, 3], [0, 2]])
        answer = np.array([[[0, 0], [0, 0], [0, 0]],
                           [[0, 1], [0, 0], [1, 0]],
                           [[0, 1], [0, 0], [0, 0]],
                           [[0, 0], [0, 0], [0, 0]]])
        assert_array_equal(H, answer)

        Z = np.zeros((5, 5, 5))
        Z[list(range(5)), list(range(5)), list(range(5))] = 1.
        H, edges = histogramdd([np.arange(5), np.arange(5), np.arange(5)], 5)
        assert_array_equal(H, Z)

    def test_shape_3d(self):
        # All possible permutations for bins of different lengths in 3D.
        bins = ((5, 4, 6), (6, 4, 5), (5, 6, 4), (4, 6, 5), (6, 5, 4),
                (4, 5, 6))
        r = rand(10, 3)
        for b in bins:
            H, edges = histogramdd(r, b)
            assert_(H.shape == b)

    def test_shape_4d(self):
        # All possible permutations for bins of different lengths in 4D.
        bins = ((7, 4, 5, 6), (4, 5, 7, 6), (5, 6, 4, 7), (7, 6, 5, 4),
                (5, 7, 6, 4), (4, 6, 7, 5), (6, 5, 7, 4), (7, 5, 4, 6),
                (7, 4, 6, 5), (6, 4, 7, 5), (6, 7, 5, 4), (4, 6, 5, 7),
                (4, 7, 5, 6), (5, 4, 6, 7), (5, 7, 4, 6), (6, 7, 4, 5),
                (6, 5, 4, 7), (4, 7, 6, 5), (4, 5, 6, 7), (7, 6, 4, 5),
                (5, 4, 7, 6), (5, 6, 7, 4), (6, 4, 5, 7), (7, 5, 6, 4))

        r = rand(10, 4)
        for b in bins:
            H, edges = histogramdd(r, b)
            assert_(H.shape == b)

    def test_weights(self):
        v = rand(100, 2)
        hist, edges = histogramdd(v)
        n_hist, edges = histogramdd(v, normed=True)
        w_hist, edges = histogramdd(v, weights=np.ones(100))
        assert_array_equal(w_hist, hist)
        w_hist, edges = histogramdd(v, weights=np.ones(100) * 2, normed=True)
        assert_array_equal(w_hist, n_hist)
        w_hist, edges = histogramdd(v, weights=np.ones(100, int) * 2)
        assert_array_equal(w_hist, 2 * hist)

    def test_identical_samples(self):
        x = np.zeros((10, 2), int)
        hist, edges = histogramdd(x, bins=2)
        assert_array_equal(edges[0], np.array([-0.5, 0., 0.5]))

    def test_empty(self):
        a, b = histogramdd([[], []], bins=([0, 1], [0, 1]))
        assert_array_max_ulp(a, np.array([[0.]]))
        a, b = np.histogramdd([[], [], []], bins=2)
        assert_array_max_ulp(a, np.zeros((2, 2, 2)))

    def test_bins_errors(self):
        """There are two ways to specify bins. Check for the right errors when
        mixing those."""
        x = np.arange(8).reshape(2, 4)
        assert_raises(ValueError, np.histogramdd, x, bins=[-1, 2, 4, 5])
        assert_raises(ValueError, np.histogramdd, x, bins=[1, 0.99, 1, 1])
        assert_raises(
            ValueError, np.histogramdd, x, bins=[1, 1, 1, [1, 2, 2, 3]])
        assert_raises(
            ValueError, np.histogramdd, x, bins=[1, 1, 1, [1, 2, 3, -3]])
        assert_(np.histogramdd(x, bins=[1, 1, 1, [1, 2, 3, 4]]))

    def test_inf_edges(self):
        """Test using +/-inf bin edges works. See #1788."""
        with np.errstate(invalid='ignore'):
            x = np.arange(6).reshape(3, 2)
            expected = np.array([[1, 0], [0, 1], [0, 1]])
            h, e = np.histogramdd(x, bins=[3, [-np.inf, 2, 10]])
            assert_allclose(h, expected)
            h, e = np.histogramdd(x, bins=[3, np.array([-1, 2, np.inf])])
            assert_allclose(h, expected)
            h, e = np.histogramdd(x, bins=[3, [-np.inf, 3, np.inf]])
            assert_allclose(h, expected)

    def test_rightmost_binedge(self):
        """Test event very close to rightmost binedge.
        See Github issue #4266"""
        x = [0.9999999995]
        bins = [[0.,0.5,1.0]]
        hist, _ = histogramdd(x, bins=bins)
        assert_(hist[0] == 0.0)
        assert_(hist[1] == 1.)
        x = [1.0]
        bins = [[0.,0.5,1.0]]
        hist, _ = histogramdd(x, bins=bins)
        assert_(hist[0] == 0.0)
        assert_(hist[1] == 1.)
        x = [1.0000000001]
        bins = [[0.,0.5,1.0]]
        hist, _ = histogramdd(x, bins=bins)
        assert_(hist[0] == 0.0)
        assert_(hist[1] == 1.)
        x = [1.0001]
        bins = [[0.,0.5,1.0]]
        hist, _ = histogramdd(x, bins=bins)
        assert_(hist[0] == 0.0)
        assert_(hist[1] == 0.0)


class TestUnique(TestCase):
    def test_simple(self):
        x = np.array([4, 3, 2, 1, 1, 2, 3, 4, 0])
        assert_(np.all(unique(x) == [0, 1, 2, 3, 4]))
        assert_(unique(np.array([1, 1, 1, 1, 1])) == np.array([1]))
        x = ['widget', 'ham', 'foo', 'bar', 'foo', 'ham']
        assert_(np.all(unique(x) == ['bar', 'foo', 'ham', 'widget']))
        x = np.array([5 + 6j, 1 + 1j, 1 + 10j, 10, 5 + 6j])
        assert_(np.all(unique(x) == [1 + 1j, 1 + 10j, 5 + 6j, 10]))


class TestCheckFinite(TestCase):
    def test_simple(self):
        a = [1, 2, 3]
        b = [1, 2, np.inf]
        c = [1, 2, np.nan]
        np.lib.asarray_chkfinite(a)
        assert_raises(ValueError, np.lib.asarray_chkfinite, b)
        assert_raises(ValueError, np.lib.asarray_chkfinite, c)

    def test_dtype_order(self):
        """Regression test for missing dtype and order arguments"""
        a = [1, 2, 3]
        a = np.lib.asarray_chkfinite(a, order='F', dtype=np.float64)
        assert_(a.dtype == np.float64)


class TestCorrCoef(TestCase):
    A = np.array(
        [[0.15391142, 0.18045767, 0.14197213],
         [0.70461506, 0.96474128, 0.27906989],
         [0.9297531, 0.32296769, 0.19267156]])
    B = np.array(
        [[0.10377691, 0.5417086, 0.49807457],
         [0.82872117, 0.77801674, 0.39226705],
         [0.9314666, 0.66800209, 0.03538394]])
    res1 = np.array(
        [[1., 0.9379533, -0.04931983],
         [0.9379533, 1., 0.30007991],
         [-0.04931983, 0.30007991, 1.]])
    res2 = np.array(
        [[1., 0.9379533, -0.04931983, 0.30151751, 0.66318558, 0.51532523],
         [0.9379533, 1., 0.30007991, -0.04781421, 0.88157256, 0.78052386],
         [-0.04931983, 0.30007991, 1., -0.96717111, 0.71483595, 0.83053601],
         [0.30151751, -0.04781421, -0.96717111, 1., -0.51366032, -0.66173113],
         [0.66318558, 0.88157256, 0.71483595, -0.51366032, 1., 0.98317823],
         [0.51532523, 0.78052386, 0.83053601, -0.66173113, 0.98317823, 1.]])

    def test_non_array(self):
        assert_almost_equal(np.corrcoef([0, 1, 0], [1, 0, 1]),
                            [[1., -1.], [-1.,  1.]])

    def test_simple(self):
        assert_almost_equal(corrcoef(self.A), self.res1)
        assert_almost_equal(corrcoef(self.A, self.B), self.res2)

    def test_ddof(self):
        assert_almost_equal(corrcoef(self.A, ddof=-1), self.res1)
        assert_almost_equal(corrcoef(self.A, self.B, ddof=-1), self.res2)

    def test_complex(self):
        x = np.array([[1, 2, 3], [1j, 2j, 3j]])
        assert_allclose(corrcoef(x), np.array([[1., -1.j], [1.j, 1.]]))

    def test_xy(self):
        x = np.array([[1, 2, 3]])
        y = np.array([[1j, 2j, 3j]])
        assert_allclose(np.corrcoef(x, y), np.array([[1., -1.j], [1.j, 1.]]))

    def test_empty(self):
        with warnings.catch_warnings(record=True):
            warnings.simplefilter('always', RuntimeWarning)
            assert_array_equal(corrcoef(np.array([])), np.nan)
            assert_array_equal(corrcoef(np.array([]).reshape(0, 2)),
                               np.array([]).reshape(0, 0))
            assert_array_equal(corrcoef(np.array([]).reshape(2, 0)),
                               np.array([[np.nan, np.nan], [np.nan, np.nan]]))

    def test_wrong_ddof(self):
        x = np.array([[0, 2], [1, 1], [2, 0]]).T
        with warnings.catch_warnings(record=True):
            warnings.simplefilter('always', RuntimeWarning)
            assert_array_equal(corrcoef(x, ddof=5),
                               np.array([[np.nan, np.nan], [np.nan, np.nan]]))


class TestCov(TestCase):
    def test_basic(self):
        x = np.array([[0, 2], [1, 1], [2, 0]]).T
        assert_allclose(cov(x), np.array([[1., -1.], [-1., 1.]]))

    def test_complex(self):
        x = np.array([[1, 2, 3], [1j, 2j, 3j]])
        assert_allclose(cov(x), np.array([[1., -1.j], [1.j, 1.]]))

    def test_xy(self):
        x = np.array([[1, 2, 3]])
        y = np.array([[1j, 2j, 3j]])
        assert_allclose(cov(x, y), np.array([[1., -1.j], [1.j, 1.]]))

    def test_empty(self):
        with warnings.catch_warnings(record=True):
            warnings.simplefilter('always', RuntimeWarning)
            assert_array_equal(cov(np.array([])), np.nan)
            assert_array_equal(cov(np.array([]).reshape(0, 2)),
                               np.array([]).reshape(0, 0))
            assert_array_equal(cov(np.array([]).reshape(2, 0)),
                               np.array([[np.nan, np.nan], [np.nan, np.nan]]))

    def test_wrong_ddof(self):
        x = np.array([[0, 2], [1, 1], [2, 0]]).T
        with warnings.catch_warnings(record=True):
            warnings.simplefilter('always', RuntimeWarning)
            assert_array_equal(cov(x, ddof=5),
                               np.array([[np.inf, -np.inf], [-np.inf, np.inf]]))


class Test_I0(TestCase):
    def test_simple(self):
        assert_almost_equal(
            i0(0.5),
            np.array(1.0634833707413234))

        A = np.array([0.49842636, 0.6969809, 0.22011976, 0.0155549])
        assert_almost_equal(
            i0(A),
            np.array([1.06307822, 1.12518299, 1.01214991, 1.00006049]))

        B = np.array([[0.827002, 0.99959078],
                      [0.89694769, 0.39298162],
                      [0.37954418, 0.05206293],
                      [0.36465447, 0.72446427],
                      [0.48164949, 0.50324519]])
        assert_almost_equal(
            i0(B),
            np.array([[1.17843223, 1.26583466],
                      [1.21147086, 1.03898290],
                      [1.03633899, 1.00067775],
                      [1.03352052, 1.13557954],
                      [1.05884290, 1.06432317]]))


class TestKaiser(TestCase):
    def test_simple(self):
        assert_(np.isfinite(kaiser(1, 1.0)))
        assert_almost_equal(kaiser(0, 1.0),
                            np.array([]))
        assert_almost_equal(kaiser(2, 1.0),
                            np.array([0.78984831, 0.78984831]))
        assert_almost_equal(kaiser(5, 1.0),
                            np.array([0.78984831, 0.94503323, 1.,
                                      0.94503323, 0.78984831]))
        assert_almost_equal(kaiser(5, 1.56789),
                            np.array([0.58285404, 0.88409679, 1.,
                                      0.88409679, 0.58285404]))

    def test_int_beta(self):
        kaiser(3, 4)


class TestMsort(TestCase):
    def test_simple(self):
        A = np.array([[0.44567325, 0.79115165, 0.54900530],
                      [0.36844147, 0.37325583, 0.96098397],
                      [0.64864341, 0.52929049, 0.39172155]])
        assert_almost_equal(
            msort(A),
            np.array([[0.36844147, 0.37325583, 0.39172155],
                      [0.44567325, 0.52929049, 0.54900530],
                      [0.64864341, 0.79115165, 0.96098397]]))


class TestMeshgrid(TestCase):
    def test_simple(self):
        [X, Y] = meshgrid([1, 2, 3], [4, 5, 6, 7])
        assert_array_equal(X, np.array([[1, 2, 3],
                                        [1, 2, 3],
                                        [1, 2, 3],
                                        [1, 2, 3]]))
        assert_array_equal(Y, np.array([[4, 4, 4],
                                        [5, 5, 5],
                                        [6, 6, 6],
                                        [7, 7, 7]]))

    def test_single_input(self):
        [X] = meshgrid([1, 2, 3, 4])
        assert_array_equal(X, np.array([1, 2, 3, 4]))

    def test_no_input(self):
        args = []
        assert_array_equal([], meshgrid(*args))

    def test_indexing(self):
        x = [1, 2, 3]
        y = [4, 5, 6, 7]
        [X, Y] = meshgrid(x, y, indexing='ij')
        assert_array_equal(X, np.array([[1, 1, 1, 1],
                                        [2, 2, 2, 2],
                                        [3, 3, 3, 3]]))
        assert_array_equal(Y, np.array([[4, 5, 6, 7],
                                        [4, 5, 6, 7],
                                        [4, 5, 6, 7]]))

        # Test expected shapes:
        z = [8, 9]
        assert_(meshgrid(x, y)[0].shape == (4, 3))
        assert_(meshgrid(x, y, indexing='ij')[0].shape == (3, 4))
        assert_(meshgrid(x, y, z)[0].shape == (4, 3, 2))
        assert_(meshgrid(x, y, z, indexing='ij')[0].shape == (3, 4, 2))

        assert_raises(ValueError, meshgrid, x, y, indexing='notvalid')

    def test_sparse(self):
        [X, Y] = meshgrid([1, 2, 3], [4, 5, 6, 7], sparse=True)
        assert_array_equal(X, np.array([[1, 2, 3]]))
        assert_array_equal(Y, np.array([[4], [5], [6], [7]]))

    def test_invalid_arguments(self):
        # Test that meshgrid complains about invalid arguments
        # Regression test for issue #4755:
        # https://github.com/numpy/numpy/issues/4755
        assert_raises(TypeError, meshgrid,
                      [1, 2, 3], [4, 5, 6, 7], indices='ij')


class TestPiecewise(TestCase):
    def test_simple(self):
        # Condition is single bool list
        x = piecewise([0, 0], [True, False], [1])
        assert_array_equal(x, [1, 0])

        # List of conditions: single bool list
        x = piecewise([0, 0], [[True, False]], [1])
        assert_array_equal(x, [1, 0])

        # Conditions is single bool array
        x = piecewise([0, 0], np.array([True, False]), [1])
        assert_array_equal(x, [1, 0])

        # Condition is single int array
        x = piecewise([0, 0], np.array([1, 0]), [1])
        assert_array_equal(x, [1, 0])

        # List of conditions: int array
        x = piecewise([0, 0], [np.array([1, 0])], [1])
        assert_array_equal(x, [1, 0])

        x = piecewise([0, 0], [[False, True]], [lambda x:-1])
        assert_array_equal(x, [0, -1])

    def test_two_conditions(self):
        x = piecewise([1, 2], [[True, False], [False, True]], [3, 4])
        assert_array_equal(x, [3, 4])

    def test_default(self):
        # No value specified for x[1], should be 0
        x = piecewise([1, 2], [True, False], [2])
        assert_array_equal(x, [2, 0])

        # Should set x[1] to 3
        x = piecewise([1, 2], [True, False], [2, 3])
        assert_array_equal(x, [2, 3])

    def test_0d(self):
        x = np.array(3)
        y = piecewise(x, x > 3, [4, 0])
        assert_(y.ndim == 0)
        assert_(y == 0)

        x = 5
        y = piecewise(x, [[True], [False]], [1, 0])
        assert_(y.ndim == 0)
        assert_(y == 1)

    def test_0d_comparison(self):
        x = 3
        y = piecewise(x, [x <= 3, x > 3], [4, 0])


class TestBincount(TestCase):
    def test_simple(self):
        y = np.bincount(np.arange(4))
        assert_array_equal(y, np.ones(4))

    def test_simple2(self):
        y = np.bincount(np.array([1, 5, 2, 4, 1]))
        assert_array_equal(y, np.array([0, 2, 1, 0, 1, 1]))

    def test_simple_weight(self):
        x = np.arange(4)
        w = np.array([0.2, 0.3, 0.5, 0.1])
        y = np.bincount(x, w)
        assert_array_equal(y, w)

    def test_simple_weight2(self):
        x = np.array([1, 2, 4, 5, 2])
        w = np.array([0.2, 0.3, 0.5, 0.1, 0.2])
        y = np.bincount(x, w)
        assert_array_equal(y, np.array([0, 0.2, 0.5, 0, 0.5, 0.1]))

    def test_with_minlength(self):
        x = np.array([0, 1, 0, 1, 1])
        y = np.bincount(x, minlength=3)
        assert_array_equal(y, np.array([2, 3, 0]))

    def test_with_minlength_smaller_than_maxvalue(self):
        x = np.array([0, 1, 1, 2, 2, 3, 3])
        y = np.bincount(x, minlength=2)
        assert_array_equal(y, np.array([1, 2, 2, 2]))

    def test_with_minlength_and_weights(self):
        x = np.array([1, 2, 4, 5, 2])
        w = np.array([0.2, 0.3, 0.5, 0.1, 0.2])
        y = np.bincount(x, w, 8)
        assert_array_equal(y, np.array([0, 0.2, 0.5, 0, 0.5, 0.1, 0, 0]))

    def test_empty(self):
        x = np.array([], dtype=int)
        y = np.bincount(x)
        assert_array_equal(x, y)

    def test_empty_with_minlength(self):
        x = np.array([], dtype=int)
        y = np.bincount(x, minlength=5)
        assert_array_equal(y, np.zeros(5, dtype=int))

    def test_with_incorrect_minlength(self):
        x = np.array([], dtype=int)
        assert_raises_regex(TypeError, "an integer is required",
                            lambda: np.bincount(x, minlength="foobar"))
        assert_raises_regex(ValueError, "must be positive",
                            lambda: np.bincount(x, minlength=-1))
        assert_raises_regex(ValueError, "must be positive",
                            lambda: np.bincount(x, minlength=0))

        x = np.arange(5)
        assert_raises_regex(TypeError, "an integer is required",
                            lambda: np.bincount(x, minlength="foobar"))
        assert_raises_regex(ValueError, "minlength must be positive",
                            lambda: np.bincount(x, minlength=-1))
        assert_raises_regex(ValueError, "minlength must be positive",
                            lambda: np.bincount(x, minlength=0))


class TestInterp(TestCase):
    def test_exceptions(self):
        assert_raises(ValueError, interp, 0, [], [])
        assert_raises(ValueError, interp, 0, [0], [1, 2])

    def test_basic(self):
        x = np.linspace(0, 1, 5)
        y = np.linspace(0, 1, 5)
        x0 = np.linspace(0, 1, 50)
        assert_almost_equal(np.interp(x0, x, y), x0)

    def test_right_left_behavior(self):
        assert_equal(interp([-1, 0, 1], [0], [1]), [1, 1, 1])
        assert_equal(interp([-1, 0, 1], [0], [1], left=0), [0, 1, 1])
        assert_equal(interp([-1, 0, 1], [0], [1], right=0), [1, 1, 0])
        assert_equal(interp([-1, 0, 1], [0], [1], left=0, right=0), [0, 1, 0])

    def test_scalar_interpolation_point(self):
        x = np.linspace(0, 1, 5)
        y = np.linspace(0, 1, 5)
        x0 = 0
        assert_almost_equal(np.interp(x0, x, y), x0)
        x0 = .3
        assert_almost_equal(np.interp(x0, x, y), x0)
        x0 = np.float32(.3)
        assert_almost_equal(np.interp(x0, x, y), x0)
        x0 = np.float64(.3)
        assert_almost_equal(np.interp(x0, x, y), x0)
        x0 = np.nan
        assert_almost_equal(np.interp(x0, x, y), x0)

    def test_zero_dimensional_interpolation_point(self):
        x = np.linspace(0, 1, 5)
        y = np.linspace(0, 1, 5)
        x0 = np.array(.3)
        assert_almost_equal(np.interp(x0, x, y), x0)
        x0 = np.array(.3, dtype=object)
        assert_almost_equal(np.interp(x0, x, y), .3)

    def test_if_len_x_is_small(self):
        xp = np.arange(0, 10, 0.0001)
        fp = np.sin(xp)
        assert_almost_equal(np.interp(np.pi, xp, fp), 0.0)


def compare_results(res, desired):
    for i in range(len(desired)):
        assert_array_equal(res[i], desired[i])


class TestScoreatpercentile(TestCase):

    def test_basic(self):
        x = np.arange(8) * 0.5
        assert_equal(np.percentile(x, 0), 0.)
        assert_equal(np.percentile(x, 100), 3.5)
        assert_equal(np.percentile(x, 50), 1.75)

    def test_api(self):
        d = np.ones(5)
        np.percentile(d, 5, None, None, False)
        np.percentile(d, 5, None, None, False, 'linear')
        o = np.ones((1,))
        np.percentile(d, 5, None, o, False, 'linear')

    def test_2D(self):
        x = np.array([[1, 1, 1],
                     [1, 1, 1],
                     [4, 4, 3],
                     [1, 1, 1],
                     [1, 1, 1]])
        assert_array_equal(np.percentile(x, 50, axis=0), [1, 1, 1])

    def test_linear(self):

        # Test defaults
        assert_equal(np.percentile(range(10), 50), 4.5)

        # explicitly specify interpolation_method 'fraction' (the default)
        assert_equal(np.percentile(range(10), 50,
                                   interpolation='linear'), 4.5)

    def test_lower_higher(self):

        # interpolation_method 'lower'/'higher'
        assert_equal(np.percentile(range(10), 50,
                                   interpolation='lower'), 4)
        assert_equal(np.percentile(range(10), 50,
                                   interpolation='higher'), 5)

    def test_midpoint(self):
        assert_equal(np.percentile(range(10), 51,
                                   interpolation='midpoint'), 4.5)

    def test_nearest(self):
        assert_equal(np.percentile(range(10), 51,
                                   interpolation='nearest'), 5)
        assert_equal(np.percentile(range(10), 49,
                                   interpolation='nearest'), 4)

    def test_sequence(self):
        x = np.arange(8) * 0.5
        assert_equal(np.percentile(x, [0, 100, 50]), [0, 3.5, 1.75])

    def test_axis(self):
        x = np.arange(12).reshape(3, 4)

        assert_equal(np.percentile(x, (25, 50, 100)), [2.75, 5.5, 11.0])

        r0 = [[2, 3, 4, 5], [4, 5, 6, 7], [8, 9, 10, 11]]
        assert_equal(np.percentile(x, (25, 50, 100), axis=0), r0)

        r1 = [[0.75, 1.5, 3], [4.75, 5.5, 7], [8.75, 9.5, 11]]
        assert_equal(np.percentile(x, (25, 50, 100), axis=1), np.array(r1).T)

        # ensure qth axis is always first as with np.array(old_percentile(..))
        x = np.arange(3 * 4 * 5 * 6).reshape(3, 4, 5, 6)
        assert_equal(np.percentile(x, (25, 50)).shape, (2,))
        assert_equal(np.percentile(x, (25, 50, 75)).shape, (3,))
        assert_equal(np.percentile(x, (25, 50), axis=0).shape, (2, 4, 5, 6))
        assert_equal(np.percentile(x, (25, 50), axis=1).shape, (2, 3, 5, 6))
        assert_equal(np.percentile(x, (25, 50), axis=2).shape, (2, 3, 4, 6))
        assert_equal(np.percentile(x, (25, 50), axis=3).shape, (2, 3, 4, 5))
        assert_equal(np.percentile(x, (25, 50, 75), axis=1).shape, (3, 3, 5, 6))
        assert_equal(np.percentile(x, (25, 50),
                                   interpolation="higher").shape, (2,))
        assert_equal(np.percentile(x, (25, 50, 75),
                                   interpolation="higher").shape, (3,))
        assert_equal(np.percentile(x, (25, 50), axis=0,
                                   interpolation="higher").shape, (2, 4, 5, 6))
        assert_equal(np.percentile(x, (25, 50), axis=1,
                                   interpolation="higher").shape, (2, 3, 5, 6))
        assert_equal(np.percentile(x, (25, 50), axis=2,
                                   interpolation="higher").shape, (2, 3, 4, 6))
        assert_equal(np.percentile(x, (25, 50), axis=3,
                                   interpolation="higher").shape, (2, 3, 4, 5))
        assert_equal(np.percentile(x, (25, 50, 75), axis=1,
                                   interpolation="higher").shape, (3, 3, 5, 6))

    def test_scalar_q(self):
        # test for no empty dimensions for compatiblity with old percentile
        x = np.arange(12).reshape(3, 4)
        assert_equal(np.percentile(x, 50), 5.5)
        self.assertTrue(np.isscalar(np.percentile(x, 50)))
        r0 = np.array([ 4.,  5.,  6.,  7.])
        assert_equal(np.percentile(x, 50, axis=0), r0)
        assert_equal(np.percentile(x, 50, axis=0).shape, r0.shape)
        r1 = np.array([ 1.5,  5.5,  9.5])
        assert_almost_equal(np.percentile(x, 50, axis=1), r1)
        assert_equal(np.percentile(x, 50, axis=1).shape, r1.shape)

        out = np.empty(1)
        assert_equal(np.percentile(x, 50, out=out), 5.5)
        assert_equal(out, 5.5)
        out = np.empty(4)
        assert_equal(np.percentile(x, 50, axis=0, out=out), r0)
        assert_equal(out, r0)
        out = np.empty(3)
        assert_equal(np.percentile(x, 50, axis=1, out=out), r1)
        assert_equal(out, r1)

        # test for no empty dimensions for compatiblity with old percentile
        x = np.arange(12).reshape(3, 4)
        assert_equal(np.percentile(x, 50, interpolation='lower'), 5.)
        self.assertTrue(np.isscalar(np.percentile(x, 50)))
        r0 = np.array([ 4.,  5.,  6.,  7.])
        c0 = np.percentile(x, 50, interpolation='lower', axis=0)
        assert_equal(c0, r0)
        assert_equal(c0.shape, r0.shape)
        r1 = np.array([ 1.,  5.,  9.])
        c1 = np.percentile(x, 50, interpolation='lower', axis=1)
        assert_almost_equal(c1, r1)
        assert_equal(c1.shape, r1.shape)

        out = np.empty((), dtype=x.dtype)
        c = np.percentile(x, 50, interpolation='lower', out=out)
        assert_equal(c, 5)
        assert_equal(out, 5)
        out = np.empty(4, dtype=x.dtype)
        c = np.percentile(x, 50, interpolation='lower', axis=0, out=out)
        assert_equal(c, r0)
        assert_equal(out, r0)
        out = np.empty(3, dtype=x.dtype)
        c = np.percentile(x, 50, interpolation='lower', axis=1, out=out)
        assert_equal(c, r1)
        assert_equal(out, r1)

    def test_exception(self):
        assert_raises(ValueError, np.percentile, [1, 2], 56,
                      interpolation='foobar')
        assert_raises(ValueError, np.percentile, [1], 101)
        assert_raises(ValueError, np.percentile, [1], -1)
        assert_raises(ValueError, np.percentile, [1], list(range(50)) + [101])
        assert_raises(ValueError, np.percentile, [1], list(range(50)) + [-0.1])

    def test_percentile_list(self):
        assert_equal(np.percentile([1, 2, 3], 0), 1)

    def test_percentile_out(self):
        x = np.array([1, 2, 3])
        y = np.zeros((3,))
        p = (1, 2, 3)
        np.percentile(x, p, out=y)
        assert_equal(y, np.percentile(x, p))

        x = np.array([[1, 2, 3],
                      [4, 5, 6]])

        y = np.zeros((3, 3))
        np.percentile(x, p, axis=0, out=y)
        assert_equal(y, np.percentile(x, p, axis=0))

        y = np.zeros((3, 2))
        np.percentile(x, p, axis=1, out=y)
        assert_equal(y, np.percentile(x, p, axis=1))

        x = np.arange(12).reshape(3, 4)
        # q.dim > 1, float
        r0 = np.array([[2.,  3.,  4., 5.], [4., 5., 6., 7.]])
        out = np.empty((2, 4))
        assert_equal(np.percentile(x, (25, 50), axis=0, out=out), r0)
        assert_equal(out, r0)
        r1 = np.array([[0.75,  4.75,  8.75], [1.5,  5.5,  9.5]])
        out = np.empty((2, 3))
        assert_equal(np.percentile(x, (25, 50), axis=1, out=out), r1)
        assert_equal(out, r1)

        # q.dim > 1, int
        r0 = np.array([[0,  1,  2, 3], [4, 5, 6, 7]])
        out = np.empty((2, 4), dtype=x.dtype)
        c = np.percentile(x, (25, 50), interpolation='lower', axis=0, out=out)
        assert_equal(c, r0)
        assert_equal(out, r0)
        r1 = np.array([[0,  4,  8], [1,  5,  9]])
        out = np.empty((2, 3), dtype=x.dtype)
        c = np.percentile(x, (25, 50), interpolation='lower', axis=1, out=out)
        assert_equal(c, r1)
        assert_equal(out, r1)

    def test_percentile_empty_dim(self):
        # empty dims are preserved
        d = np.arange(11*2).reshape(11, 1, 2, 1)
        assert_array_equal(np.percentile(d, 50, axis=0).shape, (1, 2, 1))
        assert_array_equal(np.percentile(d, 50, axis=1).shape, (11, 2, 1))
        assert_array_equal(np.percentile(d, 50, axis=2).shape, (11, 1, 1))
        assert_array_equal(np.percentile(d, 50, axis=3).shape, (11, 1, 2))
        assert_array_equal(np.percentile(d, 50, axis=-1).shape, (11, 1, 2))
        assert_array_equal(np.percentile(d, 50, axis=-2).shape, (11, 1, 1))
        assert_array_equal(np.percentile(d, 50, axis=-3).shape, (11, 2, 1))
        assert_array_equal(np.percentile(d, 50, axis=-4).shape, (1, 2, 1))

        assert_array_equal(np.percentile(d, 50, axis=2,
                                         interpolation='midpoint').shape,
                           (11, 1, 1))
        assert_array_equal(np.percentile(d, 50, axis=-2,
                                         interpolation='midpoint').shape,
                           (11, 1, 1))

        assert_array_equal(np.array(np.percentile(d, [10, 50], axis=0)).shape,
                           (2, 1, 2, 1))
        assert_array_equal(np.array(np.percentile(d, [10, 50], axis=1)).shape,
                           (2, 11, 2, 1))
        assert_array_equal(np.array(np.percentile(d, [10, 50], axis=2)).shape,
                           (2, 11, 1, 1))
        assert_array_equal(np.array(np.percentile(d, [10, 50], axis=3)).shape,
                           (2, 11, 1, 2))


    def test_percentile_no_overwrite(self):
        a = np.array([2, 3, 4, 1])
        np.percentile(a, [50], overwrite_input=False)
        assert_equal(a, np.array([2, 3, 4, 1]))

        a = np.array([2, 3, 4, 1])
        np.percentile(a, [50])
        assert_equal(a, np.array([2, 3, 4, 1]))

    def test_no_p_overwrite(self):
        p = np.linspace(0., 100., num=5)
        np.percentile(np.arange(100.), p, interpolation="midpoint")
        assert_array_equal(p, np.linspace(0., 100., num=5))
        p = np.linspace(0., 100., num=5).tolist()
        np.percentile(np.arange(100.), p, interpolation="midpoint")
        assert_array_equal(p, np.linspace(0., 100., num=5).tolist())

    def test_percentile_overwrite(self):
        a = np.array([2, 3, 4, 1])
        b = np.percentile(a, [50], overwrite_input=True)
        assert_equal(b, np.array([2.5]))

        b = np.percentile([2, 3, 4, 1], [50], overwrite_input=True)
        assert_equal(b, np.array([2.5]))

    def test_extended_axis(self):
        o = np.random.normal(size=(71, 23))
        x = np.dstack([o] * 10)
        assert_equal(np.percentile(x, 30, axis=(0, 1)), np.percentile(o, 30))
        x = np.rollaxis(x, -1, 0)
        assert_equal(np.percentile(x, 30, axis=(-2, -1)), np.percentile(o, 30))
        x = x.swapaxes(0, 1).copy()
        assert_equal(np.percentile(x, 30, axis=(0, -1)), np.percentile(o, 30))
        x = x.swapaxes(0, 1).copy()

        assert_equal(np.percentile(x, [25, 60], axis=(0, 1, 2)),
                     np.percentile(x, [25, 60], axis=None))
        assert_equal(np.percentile(x, [25, 60], axis=(0,)),
                     np.percentile(x, [25, 60], axis=0))

        d = np.arange(3 * 5 * 7 * 11).reshape(3, 5, 7, 11)
        np.random.shuffle(d)
        assert_equal(np.percentile(d, 25,  axis=(0, 1, 2))[0],
                     np.percentile(d[:, :, :, 0].flatten(), 25))
        assert_equal(np.percentile(d, [10, 90], axis=(0, 1, 3))[:, 1],
                     np.percentile(d[:, :, 1, :].flatten(), [10, 90]))
        assert_equal(np.percentile(d, 25, axis=(3, 1, -4))[2],
                     np.percentile(d[:, :, 2, :].flatten(), 25))
        assert_equal(np.percentile(d, 25, axis=(3, 1, 2))[2],
                     np.percentile(d[2, :, :, :].flatten(), 25))
        assert_equal(np.percentile(d, 25, axis=(3, 2))[2, 1],
                     np.percentile(d[2, 1, :, :].flatten(), 25))
        assert_equal(np.percentile(d, 25, axis=(1, -2))[2, 1],
                     np.percentile(d[2, :, :, 1].flatten(), 25))
        assert_equal(np.percentile(d, 25, axis=(1, 3))[2, 2],
                     np.percentile(d[2, :, 2, :].flatten(), 25))

    def test_extended_axis_invalid(self):
        d = np.ones((3, 5, 7, 11))
        assert_raises(IndexError, np.percentile, d, axis=-5, q=25)
        assert_raises(IndexError, np.percentile, d, axis=(0, -5), q=25)
        assert_raises(IndexError, np.percentile, d, axis=4, q=25)
        assert_raises(IndexError, np.percentile, d, axis=(0, 4), q=25)
        assert_raises(ValueError, np.percentile, d, axis=(1, 1), q=25)

    def test_keepdims(self):
        d = np.ones((3, 5, 7, 11))
        assert_equal(np.percentile(d, 7, axis=None, keepdims=True).shape,
                     (1, 1, 1, 1))
        assert_equal(np.percentile(d, 7, axis=(0, 1), keepdims=True).shape,
                     (1, 1, 7, 11))
        assert_equal(np.percentile(d, 7, axis=(0, 3), keepdims=True).shape,
                     (1, 5, 7, 1))
        assert_equal(np.percentile(d, 7, axis=(1,), keepdims=True).shape,
                     (3, 1, 7, 11))
        assert_equal(np.percentile(d, 7, (0, 1, 2, 3), keepdims=True).shape,
                     (1, 1, 1, 1))
        assert_equal(np.percentile(d, 7, axis=(0, 1, 3), keepdims=True).shape,
                     (1, 1, 7, 1))

        assert_equal(np.percentile(d, [1, 7], axis=(0, 1, 3),
                                   keepdims=True).shape, (2, 1, 1, 7, 1))
        assert_equal(np.percentile(d, [1, 7], axis=(0, 3),
                                   keepdims=True).shape, (2, 1, 5, 7, 1))


class TestMedian(TestCase):
    def test_basic(self):
        a0 = np.array(1)
        a1 = np.arange(2)
        a2 = np.arange(6).reshape(2, 3)
        assert_equal(np.median(a0), 1)
        assert_allclose(np.median(a1), 0.5)
        assert_allclose(np.median(a2), 2.5)
        assert_allclose(np.median(a2, axis=0), [1.5,  2.5,  3.5])
        assert_equal(np.median(a2, axis=1), [1, 4])
        assert_allclose(np.median(a2, axis=None), 2.5)

        a = np.array([0.0444502, 0.0463301, 0.141249, 0.0606775])
        assert_almost_equal((a[1] + a[3]) / 2., np.median(a))
        a = np.array([0.0463301, 0.0444502, 0.141249])
        assert_equal(a[0], np.median(a))
        a = np.array([0.0444502, 0.141249, 0.0463301])
        assert_equal(a[-1], np.median(a))
        # check array scalar result
        assert_equal(np.median(a).ndim, 0)
        a[1] = np.nan
        assert_equal(np.median(a).ndim, 0)

    def test_axis_keyword(self):
        a3 = np.array([[2, 3],
                       [0, 1],
                       [6, 7],
                       [4, 5]])
        for a in [a3, np.random.randint(0, 100, size=(2, 3, 4))]:
            orig = a.copy()
            np.median(a, axis=None)
            for ax in range(a.ndim):
                np.median(a, axis=ax)
            assert_array_equal(a, orig)

        assert_allclose(np.median(a3, axis=0), [3,  4])
        assert_allclose(np.median(a3.T, axis=1), [3,  4])
        assert_allclose(np.median(a3), 3.5)
        assert_allclose(np.median(a3, axis=None), 3.5)
        assert_allclose(np.median(a3.T), 3.5)

    def test_overwrite_keyword(self):
        a3 = np.array([[2, 3],
                       [0, 1],
                       [6, 7],
                       [4, 5]])
        a0 = np.array(1)
        a1 = np.arange(2)
        a2 = np.arange(6).reshape(2, 3)
        assert_allclose(np.median(a0.copy(), overwrite_input=True), 1)
        assert_allclose(np.median(a1.copy(), overwrite_input=True), 0.5)
        assert_allclose(np.median(a2.copy(), overwrite_input=True), 2.5)
        assert_allclose(np.median(a2.copy(), overwrite_input=True, axis=0),
                        [1.5,  2.5,  3.5])
        assert_allclose(
            np.median(a2.copy(), overwrite_input=True, axis=1), [1, 4])
        assert_allclose(
            np.median(a2.copy(), overwrite_input=True, axis=None), 2.5)
        assert_allclose(
            np.median(a3.copy(), overwrite_input=True, axis=0), [3,  4])
        assert_allclose(np.median(a3.T.copy(), overwrite_input=True, axis=1),
                        [3,  4])

        a4 = np.arange(3 * 4 * 5, dtype=np.float32).reshape((3, 4, 5))
        map(np.random.shuffle, a4)
        assert_allclose(np.median(a4, axis=None),
                        np.median(a4.copy(), axis=None, overwrite_input=True))
        assert_allclose(np.median(a4, axis=0),
                        np.median(a4.copy(), axis=0, overwrite_input=True))
        assert_allclose(np.median(a4, axis=1),
                        np.median(a4.copy(), axis=1, overwrite_input=True))
        assert_allclose(np.median(a4, axis=2),
                        np.median(a4.copy(), axis=2, overwrite_input=True))

    def test_array_like(self):
        x = [1, 2, 3]
        assert_almost_equal(np.median(x), 2)
        x2 = [x]
        assert_almost_equal(np.median(x2), 2)
        assert_allclose(np.median(x2, axis=0), x)

    def test_subclass(self):
        # gh-3846
        class MySubClass(np.ndarray):
            def __new__(cls, input_array, info=None):
                obj = np.asarray(input_array).view(cls)
                obj.info = info
                return obj

            def mean(self, axis=None, dtype=None, out=None):
                return -7

        a = MySubClass([1,2,3])
        assert_equal(np.median(a), -7)

    def test_object(self):
        o = np.arange(7.);
        assert_(type(np.median(o.astype(object))), float)
        o[2] = np.nan
        assert_(type(np.median(o.astype(object))), float)

    def test_extended_axis(self):
        o = np.random.normal(size=(71, 23))
        x = np.dstack([o] * 10)
        assert_equal(np.median(x, axis=(0, 1)), np.median(o))
        x = np.rollaxis(x, -1, 0)
        assert_equal(np.median(x, axis=(-2, -1)), np.median(o))
        x = x.swapaxes(0, 1).copy()
        assert_equal(np.median(x, axis=(0, -1)), np.median(o))

        assert_equal(np.median(x, axis=(0, 1, 2)), np.median(x, axis=None))
        assert_equal(np.median(x, axis=(0, )), np.median(x, axis=0))
        assert_equal(np.median(x, axis=(-1, )), np.median(x, axis=-1))

        d = np.arange(3 * 5 * 7 * 11).reshape(3, 5, 7, 11)
        np.random.shuffle(d)
        assert_equal(np.median(d, axis=(0, 1, 2))[0],
                     np.median(d[:, :, :, 0].flatten()))
        assert_equal(np.median(d, axis=(0, 1, 3))[1],
                     np.median(d[:, :, 1, :].flatten()))
        assert_equal(np.median(d, axis=(3, 1, -4))[2],
                     np.median(d[:, :, 2, :].flatten()))
        assert_equal(np.median(d, axis=(3, 1, 2))[2],
                     np.median(d[2, :, :, :].flatten()))
        assert_equal(np.median(d, axis=(3, 2))[2, 1],
                     np.median(d[2, 1, :, :].flatten()))
        assert_equal(np.median(d, axis=(1, -2))[2, 1],
                     np.median(d[2, :, :, 1].flatten()))
        assert_equal(np.median(d, axis=(1, 3))[2, 2],
                     np.median(d[2, :, 2, :].flatten()))

    def test_extended_axis_invalid(self):
        d = np.ones((3, 5, 7, 11))
        assert_raises(IndexError, np.median, d, axis=-5)
        assert_raises(IndexError, np.median, d, axis=(0, -5))
        assert_raises(IndexError, np.median, d, axis=4)
        assert_raises(IndexError, np.median, d, axis=(0, 4))
        assert_raises(ValueError, np.median, d, axis=(1, 1))

    def test_keepdims(self):
        d = np.ones((3, 5, 7, 11))
        assert_equal(np.median(d, axis=None, keepdims=True).shape,
                     (1, 1, 1, 1))
        assert_equal(np.median(d, axis=(0, 1), keepdims=True).shape,
                     (1, 1, 7, 11))
        assert_equal(np.median(d, axis=(0, 3), keepdims=True).shape,
                     (1, 5, 7, 1))
        assert_equal(np.median(d, axis=(1,), keepdims=True).shape,
                     (3, 1, 7, 11))
        assert_equal(np.median(d, axis=(0, 1, 2, 3), keepdims=True).shape,
                     (1, 1, 1, 1))
        assert_equal(np.median(d, axis=(0, 1, 3), keepdims=True).shape,
                     (1, 1, 7, 1))



class TestAdd_newdoc_ufunc(TestCase):

    def test_ufunc_arg(self):
        assert_raises(TypeError, add_newdoc_ufunc, 2, "blah")
        assert_raises(ValueError, add_newdoc_ufunc, np.add, "blah")

    def test_string_arg(self):
        assert_raises(TypeError, add_newdoc_ufunc, np.add, 3)


class TestAdd_newdoc(TestCase):
    def test_add_doc(self):
        # test np.add_newdoc
        tgt = "Current flat index into the array."
        self.assertEqual(np.core.flatiter.index.__doc__[:len(tgt)], tgt)
        self.assertTrue(len(np.core.ufunc.identity.__doc__) > 300)
        self.assertTrue(len(np.lib.index_tricks.mgrid.__doc__) > 300)


if __name__ == "__main__":
    run_module_suite()
