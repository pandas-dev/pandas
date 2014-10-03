from __future__ import division, absolute_import, print_function

import warnings
import numpy as np
from numpy.testing import (TestCase, assert_, assert_raises, assert_array_equal,
                           assert_equal, run_module_suite)
from numpy.core import (array, arange, atleast_1d, atleast_2d, atleast_3d,
                        vstack, hstack, newaxis, concatenate)
from numpy.compat import long

class TestAtleast1d(TestCase):
    def test_0D_array(self):
        a = array(1)
        b = array(2)
        res = [atleast_1d(a), atleast_1d(b)]
        desired = [array([1]), array([2])]
        assert_array_equal(res, desired)

    def test_1D_array(self):
        a = array([1, 2])
        b = array([2, 3])
        res = [atleast_1d(a), atleast_1d(b)]
        desired = [array([1, 2]), array([2, 3])]
        assert_array_equal(res, desired)

    def test_2D_array(self):
        a = array([[1, 2], [1, 2]])
        b = array([[2, 3], [2, 3]])
        res = [atleast_1d(a), atleast_1d(b)]
        desired = [a, b]
        assert_array_equal(res, desired)

    def test_3D_array(self):
        a = array([[1, 2], [1, 2]])
        b = array([[2, 3], [2, 3]])
        a = array([a, a])
        b = array([b, b])
        res = [atleast_1d(a), atleast_1d(b)]
        desired = [a, b]
        assert_array_equal(res, desired)

    def test_r1array(self):
        """ Test to make sure equivalent Travis O's r1array function
        """
        assert_(atleast_1d(3).shape == (1,))
        assert_(atleast_1d(3j).shape == (1,))
        assert_(atleast_1d(long(3)).shape == (1,))
        assert_(atleast_1d(3.0).shape == (1,))
        assert_(atleast_1d([[2, 3], [4, 5]]).shape == (2, 2))


class TestAtleast2d(TestCase):
    def test_0D_array(self):
        a = array(1)
        b = array(2)
        res = [atleast_2d(a), atleast_2d(b)]
        desired = [array([[1]]), array([[2]])]
        assert_array_equal(res, desired)

    def test_1D_array(self):
        a = array([1, 2])
        b = array([2, 3])
        res = [atleast_2d(a), atleast_2d(b)]
        desired = [array([[1, 2]]), array([[2, 3]])]
        assert_array_equal(res, desired)

    def test_2D_array(self):
        a = array([[1, 2], [1, 2]])
        b = array([[2, 3], [2, 3]])
        res = [atleast_2d(a), atleast_2d(b)]
        desired = [a, b]
        assert_array_equal(res, desired)

    def test_3D_array(self):
        a = array([[1, 2], [1, 2]])
        b = array([[2, 3], [2, 3]])
        a = array([a, a])
        b = array([b, b])
        res = [atleast_2d(a), atleast_2d(b)]
        desired = [a, b]
        assert_array_equal(res, desired)

    def test_r2array(self):
        """ Test to make sure equivalent Travis O's r2array function
        """
        assert_(atleast_2d(3).shape == (1, 1))
        assert_(atleast_2d([3j, 1]).shape == (1, 2))
        assert_(atleast_2d([[[3, 1], [4, 5]], [[3, 5], [1, 2]]]).shape == (2, 2, 2))


class TestAtleast3d(TestCase):
    def test_0D_array(self):
        a = array(1)
        b = array(2)
        res = [atleast_3d(a), atleast_3d(b)]
        desired = [array([[[1]]]), array([[[2]]])]
        assert_array_equal(res, desired)

    def test_1D_array(self):
        a = array([1, 2])
        b = array([2, 3])
        res = [atleast_3d(a), atleast_3d(b)]
        desired = [array([[[1], [2]]]), array([[[2], [3]]])]
        assert_array_equal(res, desired)

    def test_2D_array(self):
        a = array([[1, 2], [1, 2]])
        b = array([[2, 3], [2, 3]])
        res = [atleast_3d(a), atleast_3d(b)]
        desired = [a[:,:, newaxis], b[:,:, newaxis]]
        assert_array_equal(res, desired)

    def test_3D_array(self):
        a = array([[1, 2], [1, 2]])
        b = array([[2, 3], [2, 3]])
        a = array([a, a])
        b = array([b, b])
        res = [atleast_3d(a), atleast_3d(b)]
        desired = [a, b]
        assert_array_equal(res, desired)


class TestHstack(TestCase):
    def test_0D_array(self):
        a = array(1)
        b = array(2)
        res=hstack([a, b])
        desired = array([1, 2])
        assert_array_equal(res, desired)

    def test_1D_array(self):
        a = array([1])
        b = array([2])
        res=hstack([a, b])
        desired = array([1, 2])
        assert_array_equal(res, desired)

    def test_2D_array(self):
        a = array([[1], [2]])
        b = array([[1], [2]])
        res=hstack([a, b])
        desired = array([[1, 1], [2, 2]])
        assert_array_equal(res, desired)


class TestVstack(TestCase):
    def test_0D_array(self):
        a = array(1)
        b = array(2)
        res=vstack([a, b])
        desired = array([[1], [2]])
        assert_array_equal(res, desired)

    def test_1D_array(self):
        a = array([1])
        b = array([2])
        res=vstack([a, b])
        desired = array([[1], [2]])
        assert_array_equal(res, desired)

    def test_2D_array(self):
        a = array([[1], [2]])
        b = array([[1], [2]])
        res=vstack([a, b])
        desired = array([[1], [2], [1], [2]])
        assert_array_equal(res, desired)

    def test_2D_array2(self):
        a = array([1, 2])
        b = array([1, 2])
        res=vstack([a, b])
        desired = array([[1, 2], [1, 2]])
        assert_array_equal(res, desired)

def test_concatenate_axis_None():
    a = np.arange(4, dtype=np.float64).reshape((2, 2))
    b = list(range(3))
    c = ['x']
    r = np.concatenate((a, a), axis=None)
    assert_equal(r.dtype, a.dtype)
    assert_equal(r.ndim, 1)
    r = np.concatenate((a, b), axis=None)
    assert_equal(r.size, a.size + len(b))
    assert_equal(r.dtype, a.dtype)
    r = np.concatenate((a, b, c), axis=None)
    d = array(['0.0', '1.0', '2.0', '3.0',
               '0', '1', '2', 'x'])
    assert_array_equal(r, d)


def test_concatenate():
    # Test concatenate function
    # No arrays raise ValueError
    assert_raises(ValueError, concatenate, ())
    # Scalars cannot be concatenated
    assert_raises(ValueError, concatenate, (0,))
    assert_raises(ValueError, concatenate, (array(0),))
    # One sequence returns unmodified (but as array)
    r4 = list(range(4))
    assert_array_equal(concatenate((r4,)), r4)
    # Any sequence
    assert_array_equal(concatenate((tuple(r4),)), r4)
    assert_array_equal(concatenate((array(r4),)), r4)
    # 1D default concatenation
    r3 = list(range(3))
    assert_array_equal(concatenate((r4, r3)), r4 + r3)
    # Mixed sequence types
    assert_array_equal(concatenate((tuple(r4), r3)), r4 + r3)
    assert_array_equal(concatenate((array(r4), r3)), r4 + r3)
    # Explicit axis specification
    assert_array_equal(concatenate((r4, r3), 0), r4 + r3)
    # Including negative
    assert_array_equal(concatenate((r4, r3), -1), r4 + r3)
    # 2D
    a23 = array([[10, 11, 12], [13, 14, 15]])
    a13 = array([[0, 1, 2]])
    res = array([[10, 11, 12], [13, 14, 15], [0, 1, 2]])
    assert_array_equal(concatenate((a23, a13)), res)
    assert_array_equal(concatenate((a23, a13), 0), res)
    assert_array_equal(concatenate((a23.T, a13.T), 1), res.T)
    assert_array_equal(concatenate((a23.T, a13.T), -1), res.T)
    # Arrays much match shape
    assert_raises(ValueError, concatenate, (a23.T, a13.T), 0)
    # 3D
    res = arange(2 * 3 * 7).reshape((2, 3, 7))
    a0 = res[..., :4]
    a1 = res[..., 4:6]
    a2 = res[..., 6:]
    assert_array_equal(concatenate((a0, a1, a2), 2), res)
    assert_array_equal(concatenate((a0, a1, a2), -1), res)
    assert_array_equal(concatenate((a0.T, a1.T, a2.T), 0), res.T)


def test_concatenate_sloppy0():
    # Versions of numpy < 1.7.0 ignored axis argument value for 1D arrays.  We
    # allow this for now, but in due course we will raise an error
    r4 = list(range(4))
    r3 = list(range(3))
    assert_array_equal(concatenate((r4, r3), 0), r4 + r3)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', DeprecationWarning)
        assert_array_equal(concatenate((r4, r3), -10), r4 + r3)
        assert_array_equal(concatenate((r4, r3), 10), r4 + r3)
        # Confirm DeprecationWarning raised
        warnings.simplefilter('error', DeprecationWarning)
        assert_raises(DeprecationWarning, concatenate, (r4, r3), 10)


if __name__ == "__main__":
    run_module_suite()
