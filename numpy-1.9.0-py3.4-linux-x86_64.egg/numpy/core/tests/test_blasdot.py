from __future__ import division, absolute_import, print_function

import numpy as np
import sys
from numpy.core import zeros, float64
from numpy.testing import dec, TestCase, assert_almost_equal, assert_, \
     assert_raises, assert_array_equal, assert_allclose, assert_equal
from numpy.core.multiarray import inner as inner_

DECPREC = 14

class TestInner(TestCase):
    def test_vecself(self):
        """Ticket 844."""
        # Inner product of a vector with itself segfaults or give meaningless
        # result
        a = zeros(shape = (1, 80), dtype = float64)
        p = inner_(a, a)
        assert_almost_equal(p, 0, decimal = DECPREC)

try:
    import numpy.core._dotblas as _dotblas
except ImportError:
    _dotblas = None

@dec.skipif(_dotblas is None, "Numpy is not compiled with _dotblas")
def test_blasdot_used():
    from numpy.core import dot, vdot, inner, alterdot, restoredot
    assert_(dot is _dotblas.dot)
    assert_(vdot is _dotblas.vdot)
    assert_(inner is _dotblas.inner)
    assert_(alterdot is _dotblas.alterdot)
    assert_(restoredot is _dotblas.restoredot)


def test_dot_2args():
    from numpy.core import dot

    a = np.array([[1, 2], [3, 4]], dtype=float)
    b = np.array([[1, 0], [1, 1]], dtype=float)
    c = np.array([[3, 2], [7, 4]], dtype=float)

    d = dot(a, b)
    assert_allclose(c, d)

def test_dot_3args():
    np.random.seed(22)
    f = np.random.random_sample((1024, 16))
    v = np.random.random_sample((16, 32))

    r = np.empty((1024, 32))
    for i in range(12):
        np.dot(f, v, r)
    assert_equal(sys.getrefcount(r), 2)
    r2 = np.dot(f, v, out=None)
    assert_array_equal(r2, r)
    assert_(r is np.dot(f, v, out=r))

    v = v[:, 0].copy() # v.shape == (16,)
    r = r[:, 0].copy() # r.shape == (1024,)
    r2 = np.dot(f, v)
    assert_(r is np.dot(f, v, r))
    assert_array_equal(r2, r)

def test_dot_3args_errors():
    np.random.seed(22)
    f = np.random.random_sample((1024, 16))
    v = np.random.random_sample((16, 32))

    r = np.empty((1024, 31))
    assert_raises(ValueError, np.dot, f, v, r)

    r = np.empty((1024,))
    assert_raises(ValueError, np.dot, f, v, r)

    r = np.empty((32,))
    assert_raises(ValueError, np.dot, f, v, r)

    r = np.empty((32, 1024))
    assert_raises(ValueError, np.dot, f, v, r)
    assert_raises(ValueError, np.dot, f, v, r.T)

    r = np.empty((1024, 64))
    assert_raises(ValueError, np.dot, f, v, r[:, ::2])
    assert_raises(ValueError, np.dot, f, v, r[:, :32])

    r = np.empty((1024, 32), dtype=np.float32)
    assert_raises(ValueError, np.dot, f, v, r)

    r = np.empty((1024, 32), dtype=int)
    assert_raises(ValueError, np.dot, f, v, r)

def test_dot_array_order():
    """ Test numpy dot with different order C, F

    Comparing results with multiarray dot.
    Double and single precisions array are compared using relative
    precision of 7 and 5 decimals respectively.
    Use 30 decimal when comparing exact operations like:
        (a.b)' = b'.a'
    """
    _dot = np.core.multiarray.dot
    a_dim, b_dim, c_dim = 10, 4, 7
    orders = ["C", "F"]
    dtypes_prec = {np.float64: 7, np.float32: 5}
    np.random.seed(7)

    for arr_type, prec in dtypes_prec.items():
        for a_order in orders:
            a = np.asarray(np.random.randn(a_dim, a_dim),
                dtype=arr_type, order=a_order)
            assert_array_equal(np.dot(a, a), a.dot(a))
            # (a.a)' = a'.a', note that mse~=1e-31 needs almost_equal
            assert_almost_equal(a.dot(a), a.T.dot(a.T).T, decimal=prec)

            #
            # Check with making explicit copy
            #
            a_T = a.T.copy(order=a_order)
            assert_almost_equal(a_T.dot(a_T), a.T.dot(a.T), decimal=prec)
            assert_almost_equal(a.dot(a_T), a.dot(a.T), decimal=prec)
            assert_almost_equal(a_T.dot(a), a.T.dot(a), decimal=prec)

            #
            # Compare with multiarray dot
            #
            assert_almost_equal(a.dot(a), _dot(a, a), decimal=prec)
            assert_almost_equal(a.T.dot(a), _dot(a.T, a), decimal=prec)
            assert_almost_equal(a.dot(a.T), _dot(a, a.T), decimal=prec)
            assert_almost_equal(a.T.dot(a.T), _dot(a.T, a.T), decimal=prec)
            for res in a.dot(a), a.T.dot(a), a.dot(a.T), a.T.dot(a.T):
                assert res.flags.c_contiguous

            for b_order in orders:
                b = np.asarray(np.random.randn(a_dim, b_dim),
                    dtype=arr_type, order=b_order)
                b_T = b.T.copy(order=b_order)
                assert_almost_equal(a_T.dot(b), a.T.dot(b), decimal=prec)
                assert_almost_equal(b_T.dot(a), b.T.dot(a), decimal=prec)
                # (b'.a)' = a'.b
                assert_almost_equal(b.T.dot(a), a.T.dot(b).T, decimal=prec)
                assert_almost_equal(a.dot(b), _dot(a, b), decimal=prec)
                assert_almost_equal(b.T.dot(a), _dot(b.T, a), decimal=prec)


                for c_order in orders:
                    c = np.asarray(np.random.randn(b_dim, c_dim),
                        dtype=arr_type, order=c_order)
                    c_T = c.T.copy(order=c_order)
                    assert_almost_equal(c.T.dot(b.T), c_T.dot(b_T), decimal=prec)
                    assert_almost_equal(c.T.dot(b.T).T, b.dot(c), decimal=prec)
                    assert_almost_equal(b.dot(c), _dot(b, c), decimal=prec)
                    assert_almost_equal(c.T.dot(b.T), _dot(c.T, b.T), decimal=prec)

@dec.skipif(True) # ufunc override disabled for 1.9
def test_dot_override():
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
