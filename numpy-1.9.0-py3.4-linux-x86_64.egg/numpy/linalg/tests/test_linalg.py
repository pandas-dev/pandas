""" Test functions for linalg module

"""
from __future__ import division, absolute_import, print_function

import os
import sys
import itertools
import traceback

import numpy as np
from numpy import array, single, double, csingle, cdouble, dot, identity
from numpy import multiply, atleast_2d, inf, asarray, matrix
from numpy import linalg
from numpy.linalg import matrix_power, norm, matrix_rank
from numpy.testing import (
    assert_, assert_equal, assert_raises, assert_array_equal,
    assert_almost_equal, assert_allclose, run_module_suite,
    dec
)


def ifthen(a, b):
    return not a or b


def imply(a, b):
    return not a or b


old_assert_almost_equal = assert_almost_equal


def assert_almost_equal(a, b, **kw):
    if asarray(a).dtype.type in (single, csingle):
        decimal = 6
    else:
        decimal = 12
    old_assert_almost_equal(a, b, decimal=decimal, **kw)


def get_real_dtype(dtype):
    return {single: single, double: double,
            csingle: single, cdouble: double}[dtype]


def get_complex_dtype(dtype):
    return {single: csingle, double: cdouble,
            csingle: csingle, cdouble: cdouble}[dtype]

def get_rtol(dtype):
    # Choose a safe rtol
    if dtype in (single, csingle):
        return 1e-5
    else:
        return 1e-11

class LinalgCase(object):
    def __init__(self, name, a, b, exception_cls=None):
        assert isinstance(name, str)
        self.name = name
        self.a = a
        self.b = b
        self.exception_cls = exception_cls

    def check(self, do):
        if self.exception_cls is None:
            do(self.a, self.b)
        else:
            assert_raises(self.exception_cls, do, self.a, self.b)

    def __repr__(self):
        return "<LinalgCase: %s>" % (self.name,)


#
# Base test cases
#

np.random.seed(1234)

SQUARE_CASES = [
    LinalgCase("single",
               array([[1., 2.], [3., 4.]], dtype=single),
               array([2., 1.], dtype=single)),
    LinalgCase("double",
               array([[1., 2.], [3., 4.]], dtype=double),
               array([2., 1.], dtype=double)),
    LinalgCase("double_2",
               array([[1., 2.], [3., 4.]], dtype=double),
               array([[2., 1., 4.], [3., 4., 6.]], dtype=double)),
    LinalgCase("csingle",
               array([[1.+2j, 2+3j], [3+4j, 4+5j]], dtype=csingle),
               array([2.+1j, 1.+2j], dtype=csingle)),
    LinalgCase("cdouble",
               array([[1.+2j, 2+3j], [3+4j, 4+5j]], dtype=cdouble),
               array([2.+1j, 1.+2j], dtype=cdouble)),
    LinalgCase("cdouble_2",
               array([[1.+2j, 2+3j], [3+4j, 4+5j]], dtype=cdouble),
               array([[2.+1j, 1.+2j, 1+3j], [1-2j, 1-3j, 1-6j]], dtype=cdouble)),
    LinalgCase("empty",
               atleast_2d(array([], dtype = double)),
               atleast_2d(array([], dtype = double)),
               linalg.LinAlgError),
    LinalgCase("8x8",
               np.random.rand(8, 8),
               np.random.rand(8)),
    LinalgCase("1x1",
               np.random.rand(1, 1),
               np.random.rand(1)),
    LinalgCase("nonarray",
               [[1, 2], [3, 4]],
               [2, 1]),
    LinalgCase("matrix_b_only",
               array([[1., 2.], [3., 4.]]),
               matrix([2., 1.]).T),
    LinalgCase("matrix_a_and_b",
               matrix([[1., 2.], [3., 4.]]),
               matrix([2., 1.]).T),
]

NONSQUARE_CASES = [
    LinalgCase("single_nsq_1",
               array([[1., 2., 3.], [3., 4., 6.]], dtype=single),
               array([2., 1.], dtype=single)),
    LinalgCase("single_nsq_2",
               array([[1., 2.], [3., 4.], [5., 6.]], dtype=single),
               array([2., 1., 3.], dtype=single)),
    LinalgCase("double_nsq_1",
               array([[1., 2., 3.], [3., 4., 6.]], dtype=double),
               array([2., 1.], dtype=double)),
    LinalgCase("double_nsq_2",
               array([[1., 2.], [3., 4.], [5., 6.]], dtype=double),
               array([2., 1., 3.], dtype=double)),
    LinalgCase("csingle_nsq_1",
               array([[1.+1j, 2.+2j, 3.-3j], [3.-5j, 4.+9j, 6.+2j]], dtype=csingle),
               array([2.+1j, 1.+2j], dtype=csingle)),
    LinalgCase("csingle_nsq_2",
               array([[1.+1j, 2.+2j], [3.-3j, 4.-9j], [5.-4j, 6.+8j]], dtype=csingle),
               array([2.+1j, 1.+2j, 3.-3j], dtype=csingle)),
    LinalgCase("cdouble_nsq_1",
               array([[1.+1j, 2.+2j, 3.-3j], [3.-5j, 4.+9j, 6.+2j]], dtype=cdouble),
               array([2.+1j, 1.+2j], dtype=cdouble)),
    LinalgCase("cdouble_nsq_2",
               array([[1.+1j, 2.+2j], [3.-3j, 4.-9j], [5.-4j, 6.+8j]], dtype=cdouble),
               array([2.+1j, 1.+2j, 3.-3j], dtype=cdouble)),
    LinalgCase("cdouble_nsq_1_2",
               array([[1.+1j, 2.+2j, 3.-3j], [3.-5j, 4.+9j, 6.+2j]], dtype=cdouble),
               array([[2.+1j, 1.+2j], [1-1j, 2-2j]], dtype=cdouble)),
    LinalgCase("cdouble_nsq_2_2",
               array([[1.+1j, 2.+2j], [3.-3j, 4.-9j], [5.-4j, 6.+8j]], dtype=cdouble),
               array([[2.+1j, 1.+2j], [1-1j, 2-2j], [1-1j, 2-2j]], dtype=cdouble)),
    LinalgCase("8x11",
               np.random.rand(8, 11),
               np.random.rand(11)),
    LinalgCase("1x5",
               np.random.rand(1, 5),
               np.random.rand(5)),
    LinalgCase("5x1",
               np.random.rand(5, 1),
               np.random.rand(1)),
]

HERMITIAN_CASES = [
    LinalgCase("hsingle",
               array([[1., 2.], [2., 1.]], dtype=single),
               None),
    LinalgCase("hdouble",
               array([[1., 2.], [2., 1.]], dtype=double),
               None),
    LinalgCase("hcsingle",
               array([[1., 2+3j], [2-3j, 1]], dtype=csingle),
               None),
    LinalgCase("hcdouble",
               array([[1., 2+3j], [2-3j, 1]], dtype=cdouble),
               None),
    LinalgCase("hempty",
               atleast_2d(array([], dtype = double)),
               None,
               linalg.LinAlgError),
    LinalgCase("hnonarray",
               [[1, 2], [2, 1]],
               None),
    LinalgCase("matrix_b_only",
               array([[1., 2.], [2., 1.]]),
               None),
    LinalgCase("hmatrix_a_and_b",
               matrix([[1., 2.], [2., 1.]]),
               None),
    LinalgCase("hmatrix_1x1",
               np.random.rand(1, 1),
               None),
]


#
# Gufunc test cases
#

GENERALIZED_SQUARE_CASES = []
GENERALIZED_NONSQUARE_CASES = []
GENERALIZED_HERMITIAN_CASES = []

for tgt, src in ((GENERALIZED_SQUARE_CASES, SQUARE_CASES),
                 (GENERALIZED_NONSQUARE_CASES, NONSQUARE_CASES),
                 (GENERALIZED_HERMITIAN_CASES, HERMITIAN_CASES)):
    for case in src:
        if not isinstance(case.a, np.ndarray):
            continue
        
        a = np.array([case.a, 2*case.a, 3*case.a])
        if case.b is None:
            b = None
        else:
            b = np.array([case.b, 7*case.b, 6*case.b])
        new_case = LinalgCase(case.name + "_tile3", a, b,
                              case.exception_cls)
        tgt.append(new_case)

        a = np.array([case.a]*2*3).reshape((3, 2) + case.a.shape)
        if case.b is None:
            b = None
        else:
            b = np.array([case.b]*2*3).reshape((3, 2) + case.b.shape)
        new_case = LinalgCase(case.name + "_tile213", a, b,
                              case.exception_cls)
        tgt.append(new_case)

#
# Generate stride combination variations of the above
#

def _stride_comb_iter(x):
    """
    Generate cartesian product of strides for all axes
    """

    if not isinstance(x, np.ndarray):
        yield x, "nop"
        return

    stride_set = [(1,)]*x.ndim
    stride_set[-1] = (1, 3, -4)
    if x.ndim > 1:
        stride_set[-2] = (1, 3, -4)
    if x.ndim > 2:
        stride_set[-3] = (1, -4)

    for repeats in itertools.product(*tuple(stride_set)):
        new_shape = [abs(a*b) for a, b in zip(x.shape, repeats)]
        slices = tuple([slice(None, None, repeat) for repeat in repeats])

        # new array with different strides, but same data
        xi = np.empty(new_shape, dtype=x.dtype)
        xi.view(np.uint32).fill(0xdeadbeef)
        xi = xi[slices]
        xi[...] = x
        xi = xi.view(x.__class__)
        assert np.all(xi == x)
        yield xi, "stride_" + "_".join(["%+d" % j for j in repeats])

        # generate also zero strides if possible
        if x.ndim >= 1 and x.shape[-1] == 1:
            s = list(x.strides)
            s[-1] = 0
            xi = np.lib.stride_tricks.as_strided(x, strides=s)
            yield xi, "stride_xxx_0"
        if x.ndim >= 2 and x.shape[-2] == 1:
            s = list(x.strides)
            s[-2] = 0
            xi = np.lib.stride_tricks.as_strided(x, strides=s)
            yield xi, "stride_xxx_0_x"
        if x.ndim >= 2 and x.shape[:-2] == (1, 1):
            s = list(x.strides)
            s[-1] = 0
            s[-2] = 0
            xi = np.lib.stride_tricks.as_strided(x, strides=s)
            yield xi, "stride_xxx_0_0"

for src in (SQUARE_CASES,
            NONSQUARE_CASES,
            HERMITIAN_CASES,
            GENERALIZED_SQUARE_CASES,
            GENERALIZED_NONSQUARE_CASES,
            GENERALIZED_HERMITIAN_CASES):

    new_cases = []
    for case in src:
        for a, a_tag in _stride_comb_iter(case.a):
            for b, b_tag in _stride_comb_iter(case.b):
                new_case = LinalgCase(case.name + "_" + a_tag + "_" + b_tag, a, b,
                                      exception_cls=case.exception_cls)
                new_cases.append(new_case)
    src.extend(new_cases)


#
# Test different routines against the above cases
#

def _check_cases(func, cases):
    for case in cases:
        try:
            case.check(func)
        except Exception:
            msg = "In test case: %r\n\n" % case
            msg += traceback.format_exc()
            raise AssertionError(msg)

class LinalgTestCase(object):
    def test_sq_cases(self):
        _check_cases(self.do, SQUARE_CASES)


class LinalgNonsquareTestCase(object):
    def test_sq_cases(self):
        _check_cases(self.do, NONSQUARE_CASES)


class LinalgGeneralizedTestCase(object):
    @dec.slow
    def test_generalized_sq_cases(self):
        _check_cases(self.do, GENERALIZED_SQUARE_CASES)


class LinalgGeneralizedNonsquareTestCase(object):
    @dec.slow
    def test_generalized_nonsq_cases(self):
        _check_cases(self.do, GENERALIZED_NONSQUARE_CASES)


class HermitianTestCase(object):
    def test_herm_cases(self):
        _check_cases(self.do, HERMITIAN_CASES)


class HermitianGeneralizedTestCase(object):
    @dec.slow
    def test_generalized_herm_cases(self):
        _check_cases(self.do, GENERALIZED_HERMITIAN_CASES)


def dot_generalized(a, b):
    a = asarray(a)
    if a.ndim >= 3:
        if a.ndim == b.ndim:
            # matrix x matrix
            new_shape = a.shape[:-1] + b.shape[-1:]
        elif a.ndim == b.ndim + 1:
            # matrix x vector
            new_shape = a.shape[:-1]
        else:
            raise ValueError("Not implemented...")
        r = np.empty(new_shape, dtype=np.common_type(a, b))
        for c in itertools.product(*map(range, a.shape[:-2])):
            r[c] = dot(a[c], b[c])
        return r
    else:
        return dot(a, b)


def identity_like_generalized(a):
    a = asarray(a)
    if a.ndim >= 3:
        r = np.empty(a.shape, dtype=a.dtype)
        for c in itertools.product(*map(range, a.shape[:-2])):
            r[c] = identity(a.shape[-2])
        return r
    else:
        return identity(a.shape[0])


class TestSolve(LinalgTestCase, LinalgGeneralizedTestCase):
    def do(self, a, b):
        x = linalg.solve(a, b)
        assert_almost_equal(b, dot_generalized(a, x))
        assert_(imply(isinstance(b, matrix), isinstance(x, matrix)))

    def test_types(self):
        def check(dtype):
            x = np.array([[1, 0.5], [0.5, 1]], dtype=dtype)
            assert_equal(linalg.solve(x, x).dtype, dtype)
        for dtype in [single, double, csingle, cdouble]:
            yield check, dtype

    def test_0_size(self):
        class ArraySubclass(np.ndarray):
            pass
        # Test system of 0x0 matrices
        a = np.arange(8).reshape(2, 2, 2)
        b = np.arange(6).reshape(1, 2, 3).view(ArraySubclass)

        expected = linalg.solve(a, b)[:, 0:0,:]
        result = linalg.solve(a[:, 0:0, 0:0], b[:, 0:0,:])
        assert_array_equal(result, expected)
        assert_(isinstance(result, ArraySubclass))

        # Test errors for non-square and only b's dimension being 0
        assert_raises(linalg.LinAlgError, linalg.solve, a[:, 0:0, 0:1], b)
        assert_raises(ValueError, linalg.solve, a, b[:, 0:0,:])

        # Test broadcasting error
        b = np.arange(6).reshape(1, 3, 2) # broadcasting error
        assert_raises(ValueError, linalg.solve, a, b)
        assert_raises(ValueError, linalg.solve, a[0:0], b[0:0])

        # Test zero "single equations" with 0x0 matrices.
        b = np.arange(2).reshape(1, 2).view(ArraySubclass)
        expected = linalg.solve(a, b)[:, 0:0]
        result = linalg.solve(a[:, 0:0, 0:0], b[:, 0:0])
        assert_array_equal(result, expected)
        assert_(isinstance(result, ArraySubclass))

        b = np.arange(3).reshape(1, 3)
        assert_raises(ValueError, linalg.solve, a, b)
        assert_raises(ValueError, linalg.solve, a[0:0], b[0:0])
        assert_raises(ValueError, linalg.solve, a[:, 0:0, 0:0], b)

    def test_0_size_k(self):
        # test zero multiple equation (K=0) case.
        class ArraySubclass(np.ndarray):
            pass
        a = np.arange(4).reshape(1, 2, 2)
        b = np.arange(6).reshape(3, 2, 1).view(ArraySubclass)

        expected = linalg.solve(a, b)[:,:, 0:0]
        result = linalg.solve(a, b[:,:, 0:0])
        assert_array_equal(result, expected)
        assert_(isinstance(result, ArraySubclass))

        # test both zero.
        expected = linalg.solve(a, b)[:, 0:0, 0:0]
        result = linalg.solve(a[:, 0:0, 0:0], b[:,0:0, 0:0])
        assert_array_equal(result, expected)
        assert_(isinstance(result, ArraySubclass))


class TestInv(LinalgTestCase, LinalgGeneralizedTestCase):
    def do(self, a, b):
        a_inv = linalg.inv(a)
        assert_almost_equal(dot_generalized(a, a_inv),
                            identity_like_generalized(a))
        assert_(imply(isinstance(a, matrix), isinstance(a_inv, matrix)))

    def test_types(self):
        def check(dtype):
            x = np.array([[1, 0.5], [0.5, 1]], dtype=dtype)
            assert_equal(linalg.inv(x).dtype, dtype)
        for dtype in [single, double, csingle, cdouble]:
            yield check, dtype

    def test_0_size(self):
        # Check that all kinds of 0-sized arrays work
        class ArraySubclass(np.ndarray):
            pass
        a = np.zeros((0, 1, 1), dtype=np.int_).view(ArraySubclass)
        res = linalg.inv(a)
        assert_(res.dtype.type is np.float64)
        assert_equal(a.shape, res.shape)
        assert_(isinstance(a, ArraySubclass))

        a = np.zeros((0, 0), dtype=np.complex64).view(ArraySubclass)
        res = linalg.inv(a)
        assert_(res.dtype.type is np.complex64)
        assert_equal(a.shape, res.shape)


class TestEigvals(LinalgTestCase, LinalgGeneralizedTestCase):
    def do(self, a, b):
        ev = linalg.eigvals(a)
        evalues, evectors = linalg.eig(a)
        assert_almost_equal(ev, evalues)

    def test_types(self):
        def check(dtype):
            x = np.array([[1, 0.5], [0.5, 1]], dtype=dtype)
            assert_equal(linalg.eigvals(x).dtype, dtype)
            x = np.array([[1, 0.5], [-1, 1]], dtype=dtype)
            assert_equal(linalg.eigvals(x).dtype, get_complex_dtype(dtype))
        for dtype in [single, double, csingle, cdouble]:
            yield check, dtype


class TestEig(LinalgTestCase, LinalgGeneralizedTestCase):
    def do(self, a, b):
        evalues, evectors = linalg.eig(a)
        assert_allclose(dot_generalized(a, evectors),
                        np.asarray(evectors) * np.asarray(evalues)[...,None,:],
                        rtol=get_rtol(evalues.dtype))
        assert_(imply(isinstance(a, matrix), isinstance(evectors, matrix)))

    def test_types(self):
        def check(dtype):
            x = np.array([[1, 0.5], [0.5, 1]], dtype=dtype)
            w, v = np.linalg.eig(x)
            assert_equal(w.dtype, dtype)
            assert_equal(v.dtype, dtype)

            x = np.array([[1, 0.5], [-1, 1]], dtype=dtype)
            w, v = np.linalg.eig(x)
            assert_equal(w.dtype, get_complex_dtype(dtype))
            assert_equal(v.dtype, get_complex_dtype(dtype))

        for dtype in [single, double, csingle, cdouble]:
            yield check, dtype


class TestSVD(LinalgTestCase, LinalgGeneralizedTestCase):
    def do(self, a, b):
        u, s, vt = linalg.svd(a, 0)
        assert_allclose(a, dot_generalized(np.asarray(u) * np.asarray(s)[...,None,:],
                                           np.asarray(vt)),
                        rtol=get_rtol(u.dtype))
        assert_(imply(isinstance(a, matrix), isinstance(u, matrix)))
        assert_(imply(isinstance(a, matrix), isinstance(vt, matrix)))

    def test_types(self):
        def check(dtype):
            x = np.array([[1, 0.5], [0.5, 1]], dtype=dtype)
            u, s, vh = linalg.svd(x)
            assert_equal(u.dtype, dtype)
            assert_equal(s.dtype, get_real_dtype(dtype))
            assert_equal(vh.dtype, dtype)
            s = linalg.svd(x, compute_uv=False)
            assert_equal(s.dtype, get_real_dtype(dtype))

        for dtype in [single, double, csingle, cdouble]:
            yield check, dtype


class TestCondSVD(LinalgTestCase, LinalgGeneralizedTestCase):
    def do(self, a, b):
        c = asarray(a) # a might be a matrix
        s = linalg.svd(c, compute_uv=False)
        old_assert_almost_equal(s[0]/s[-1], linalg.cond(a), decimal=5)


class TestCond2(LinalgTestCase):
    def do(self, a, b):
        c = asarray(a) # a might be a matrix
        s = linalg.svd(c, compute_uv=False)
        old_assert_almost_equal(s[0]/s[-1], linalg.cond(a, 2), decimal=5)


class TestCondInf(object):
    def test(self):
        A = array([[1., 0, 0], [0, -2., 0], [0, 0, 3.]])
        assert_almost_equal(linalg.cond(A, inf), 3.)


class TestPinv(LinalgTestCase):
    def do(self, a, b):
        a_ginv = linalg.pinv(a)
        assert_almost_equal(dot(a, a_ginv), identity(asarray(a).shape[0]))
        assert_(imply(isinstance(a, matrix), isinstance(a_ginv, matrix)))


class TestDet(LinalgTestCase, LinalgGeneralizedTestCase):
    def do(self, a, b):
        d = linalg.det(a)
        (s, ld) = linalg.slogdet(a)
        if asarray(a).dtype.type in (single, double):
            ad = asarray(a).astype(double)
        else:
            ad = asarray(a).astype(cdouble)
        ev = linalg.eigvals(ad)
        assert_almost_equal(d, multiply.reduce(ev, axis=-1))
        assert_almost_equal(s * np.exp(ld), multiply.reduce(ev, axis=-1))

        s = np.atleast_1d(s)
        ld = np.atleast_1d(ld)
        m = (s != 0)
        assert_almost_equal(np.abs(s[m]), 1)
        assert_equal(ld[~m], -inf)

    def test_zero(self):
        assert_equal(linalg.det([[0.0]]), 0.0)
        assert_equal(type(linalg.det([[0.0]])), double)
        assert_equal(linalg.det([[0.0j]]), 0.0)
        assert_equal(type(linalg.det([[0.0j]])), cdouble)

        assert_equal(linalg.slogdet([[0.0]]), (0.0, -inf))
        assert_equal(type(linalg.slogdet([[0.0]])[0]), double)
        assert_equal(type(linalg.slogdet([[0.0]])[1]), double)
        assert_equal(linalg.slogdet([[0.0j]]), (0.0j, -inf))
        assert_equal(type(linalg.slogdet([[0.0j]])[0]), cdouble)
        assert_equal(type(linalg.slogdet([[0.0j]])[1]), double)

    def test_types(self):
        def check(dtype):
            x = np.array([[1, 0.5], [0.5, 1]], dtype=dtype)
            assert_equal(np.linalg.det(x).dtype, dtype)
            ph, s = np.linalg.slogdet(x)
            assert_equal(s.dtype, get_real_dtype(dtype))
            assert_equal(ph.dtype, dtype)
        for dtype in [single, double, csingle, cdouble]:
            yield check, dtype


class TestLstsq(LinalgTestCase, LinalgNonsquareTestCase):
    def do(self, a, b):
        arr = np.asarray(a)
        m, n = arr.shape
        u, s, vt = linalg.svd(a, 0)
        x, residuals, rank, sv = linalg.lstsq(a, b)
        if m <= n:
            assert_almost_equal(b, dot(a, x))
            assert_equal(rank, m)
        else:
            assert_equal(rank, n)
        assert_almost_equal(sv, sv.__array_wrap__(s))
        if rank == n and m > n:
            expect_resids = (np.asarray(abs(np.dot(a, x) - b))**2).sum(axis=0)
            expect_resids = np.asarray(expect_resids)
            if len(np.asarray(b).shape) == 1:
                expect_resids.shape = (1,)
                assert_equal(residuals.shape, expect_resids.shape)
        else:
            expect_resids = np.array([]).view(type(x))
        assert_almost_equal(residuals, expect_resids)
        assert_(np.issubdtype(residuals.dtype, np.floating))
        assert_(imply(isinstance(b, matrix), isinstance(x, matrix)))
        assert_(imply(isinstance(b, matrix), isinstance(residuals, matrix)))


class TestMatrixPower(object):
    R90 = array([[0, 1], [-1, 0]])
    Arb22 = array([[4, -7], [-2, 10]])
    noninv = array([[1, 0], [0, 0]])
    arbfloat = array([[0.1, 3.2], [1.2, 0.7]])

    large = identity(10)
    t = large[1,:].copy()
    large[1,:] = large[0,:]
    large[0,:] = t

    def test_large_power(self):
        assert_equal(matrix_power(self.R90, 2**100+2**10+2**5+1), self.R90)

    def test_large_power_trailing_zero(self):
        assert_equal(matrix_power(self.R90, 2**100+2**10+2**5), identity(2))

    def testip_zero(self):
        def tz(M):
            mz = matrix_power(M, 0)
            assert_equal(mz, identity(M.shape[0]))
            assert_equal(mz.dtype, M.dtype)
        for M in [self.Arb22, self.arbfloat, self.large]:
            yield tz, M

    def testip_one(self):
        def tz(M):
            mz = matrix_power(M, 1)
            assert_equal(mz, M)
            assert_equal(mz.dtype, M.dtype)
        for M in [self.Arb22, self.arbfloat, self.large]:
            yield tz, M

    def testip_two(self):
        def tz(M):
            mz = matrix_power(M, 2)
            assert_equal(mz, dot(M, M))
            assert_equal(mz.dtype, M.dtype)
        for M in [self.Arb22, self.arbfloat, self.large]:
            yield tz, M

    def testip_invert(self):
        def tz(M):
            mz = matrix_power(M, -1)
            assert_almost_equal(identity(M.shape[0]), dot(mz, M))
        for M in [self.R90, self.Arb22, self.arbfloat, self.large]:
            yield tz, M

    def test_invert_noninvertible(self):
        import numpy.linalg
        assert_raises(numpy.linalg.linalg.LinAlgError,
                      lambda: matrix_power(self.noninv, -1))


class TestBoolPower(object):
    def test_square(self):
        A = array([[True, False], [True, True]])
        assert_equal(matrix_power(A, 2), A)


class TestEigvalsh(HermitianTestCase, HermitianGeneralizedTestCase):
    def do(self, a, b):
        # note that eigenvalue arrays must be sorted since
        # their order isn't guaranteed.
        ev = linalg.eigvalsh(a, 'L')
        evalues, evectors = linalg.eig(a)
        ev.sort(axis=-1)
        evalues.sort(axis=-1)
        assert_allclose(ev, evalues,
                        rtol=get_rtol(ev.dtype))

        ev2 = linalg.eigvalsh(a, 'U')
        ev2.sort(axis=-1)
        assert_allclose(ev2, evalues,
                        rtol=get_rtol(ev.dtype))

    def test_types(self):
        def check(dtype):
            x = np.array([[1, 0.5], [0.5, 1]], dtype=dtype)
            w = np.linalg.eigvalsh(x)
            assert_equal(w.dtype, get_real_dtype(dtype))
        for dtype in [single, double, csingle, cdouble]:
            yield check, dtype

    def test_invalid(self):
        x = np.array([[1, 0.5], [0.5, 1]], dtype=np.float32)
        assert_raises(ValueError, np.linalg.eigvalsh, x, UPLO="lrong")
        assert_raises(ValueError, np.linalg.eigvalsh, x, "lower")
        assert_raises(ValueError, np.linalg.eigvalsh, x, "upper")

    def test_UPLO(self):
        Klo = np.array([[0, 0],[1, 0]], dtype=np.double)
        Kup = np.array([[0, 1],[0, 0]], dtype=np.double)
        tgt = np.array([-1, 1], dtype=np.double)
        rtol = get_rtol(np.double)

        # Check default is 'L'
        w = np.linalg.eigvalsh(Klo)
        assert_allclose(np.sort(w), tgt, rtol=rtol)
        # Check 'L'
        w = np.linalg.eigvalsh(Klo, UPLO='L')
        assert_allclose(np.sort(w), tgt, rtol=rtol)
        # Check 'l'
        w = np.linalg.eigvalsh(Klo, UPLO='l')
        assert_allclose(np.sort(w), tgt, rtol=rtol)
        # Check 'U'
        w = np.linalg.eigvalsh(Kup, UPLO='U')
        assert_allclose(np.sort(w), tgt, rtol=rtol)
        # Check 'u'
        w = np.linalg.eigvalsh(Kup, UPLO='u')
        assert_allclose(np.sort(w), tgt, rtol=rtol)


class TestEigh(HermitianTestCase, HermitianGeneralizedTestCase):
    def do(self, a, b):
        # note that eigenvalue arrays must be sorted since
        # their order isn't guaranteed.
        ev, evc = linalg.eigh(a)
        evalues, evectors = linalg.eig(a)
        ev.sort(axis=-1)
        evalues.sort(axis=-1)
        assert_almost_equal(ev, evalues)

        assert_allclose(dot_generalized(a, evc),
                        np.asarray(ev)[...,None,:] * np.asarray(evc),
                        rtol=get_rtol(ev.dtype))

        ev2, evc2 = linalg.eigh(a, 'U')
        ev2.sort(axis=-1)
        assert_almost_equal(ev2, evalues)

        assert_allclose(dot_generalized(a, evc2),
                        np.asarray(ev2)[...,None,:] * np.asarray(evc2),
                        rtol=get_rtol(ev.dtype), err_msg=repr(a))

    def test_types(self):
        def check(dtype):
            x = np.array([[1, 0.5], [0.5, 1]], dtype=dtype)
            w, v = np.linalg.eigh(x)
            assert_equal(w.dtype, get_real_dtype(dtype))
            assert_equal(v.dtype, dtype)
        for dtype in [single, double, csingle, cdouble]:
            yield check, dtype

    def test_invalid(self):
        x = np.array([[1, 0.5], [0.5, 1]], dtype=np.float32)
        assert_raises(ValueError, np.linalg.eigh, x, UPLO="lrong")
        assert_raises(ValueError, np.linalg.eigh, x, "lower")
        assert_raises(ValueError, np.linalg.eigh, x, "upper")

    def test_UPLO(self):
        Klo = np.array([[0, 0],[1, 0]], dtype=np.double)
        Kup = np.array([[0, 1],[0, 0]], dtype=np.double)
        tgt = np.array([-1, 1], dtype=np.double)
        rtol = get_rtol(np.double)

        # Check default is 'L'
        w, v = np.linalg.eigh(Klo)
        assert_allclose(np.sort(w), tgt, rtol=rtol)
        # Check 'L'
        w, v = np.linalg.eigh(Klo, UPLO='L')
        assert_allclose(np.sort(w), tgt, rtol=rtol)
        # Check 'l'
        w, v = np.linalg.eigh(Klo, UPLO='l')
        assert_allclose(np.sort(w), tgt, rtol=rtol)
        # Check 'U'
        w, v = np.linalg.eigh(Kup, UPLO='U')
        assert_allclose(np.sort(w), tgt, rtol=rtol)
        # Check 'u'
        w, v = np.linalg.eigh(Kup, UPLO='u')
        assert_allclose(np.sort(w), tgt, rtol=rtol)


class _TestNorm(object):

    dt = None
    dec = None

    def test_empty(self):
        assert_equal(norm([]), 0.0)
        assert_equal(norm(array([], dtype=self.dt)), 0.0)
        assert_equal(norm(atleast_2d(array([], dtype=self.dt))), 0.0)

    def test_vector(self):
        a = [1, 2, 3, 4]
        b = [-1, -2, -3, -4]
        c = [-1, 2, -3, 4]

        def _test(v):
            np.testing.assert_almost_equal(norm(v), 30**0.5,
                                           decimal=self.dec)
            np.testing.assert_almost_equal(norm(v, inf), 4.0,
                                           decimal=self.dec)
            np.testing.assert_almost_equal(norm(v, -inf), 1.0,
                                           decimal=self.dec)
            np.testing.assert_almost_equal(norm(v, 1), 10.0,
                                           decimal=self.dec)
            np.testing.assert_almost_equal(norm(v, -1), 12.0/25,
                                           decimal=self.dec)
            np.testing.assert_almost_equal(norm(v, 2), 30**0.5,
                                           decimal=self.dec)
            np.testing.assert_almost_equal(norm(v, -2), ((205./144)**-0.5),
                                           decimal=self.dec)
            np.testing.assert_almost_equal(norm(v, 0), 4,
                                           decimal=self.dec)

        for v in (a, b, c,):
            _test(v)

        for v in (array(a, dtype=self.dt), array(b, dtype=self.dt),
                  array(c, dtype=self.dt)):
            _test(v)

    def test_matrix(self):
        A = matrix([[1, 3], [5, 7]], dtype=self.dt)
        assert_almost_equal(norm(A), 84**0.5)
        assert_almost_equal(norm(A, 'fro'), 84**0.5)
        assert_almost_equal(norm(A, inf), 12.0)
        assert_almost_equal(norm(A, -inf), 4.0)
        assert_almost_equal(norm(A, 1), 10.0)
        assert_almost_equal(norm(A, -1), 6.0)
        assert_almost_equal(norm(A, 2), 9.1231056256176615)
        assert_almost_equal(norm(A, -2), 0.87689437438234041)

        assert_raises(ValueError, norm, A, 'nofro')
        assert_raises(ValueError, norm, A, -3)
        assert_raises(ValueError, norm, A, 0)

    def test_axis(self):
        # Vector norms.
        # Compare the use of `axis` with computing the norm of each row
        # or column separately.
        A = array([[1, 2, 3], [4, 5, 6]], dtype=self.dt)
        for order in [None, -1, 0, 1, 2, 3, np.Inf, -np.Inf]:
            expected0 = [norm(A[:, k], ord=order) for k in range(A.shape[1])]
            assert_almost_equal(norm(A, ord=order, axis=0), expected0)
            expected1 = [norm(A[k,:], ord=order) for k in range(A.shape[0])]
            assert_almost_equal(norm(A, ord=order, axis=1), expected1)

        # Matrix norms.
        B = np.arange(1, 25, dtype=self.dt).reshape(2, 3, 4)

        for order in [None, -2, 2, -1, 1, np.Inf, -np.Inf, 'fro']:
            assert_almost_equal(norm(A, ord=order), norm(A, ord=order,
                                                         axis=(0, 1)))

            n = norm(B, ord=order, axis=(1, 2))
            expected = [norm(B[k], ord=order) for k in range(B.shape[0])]
            assert_almost_equal(n, expected)

            n = norm(B, ord=order, axis=(2, 1))
            expected = [norm(B[k].T, ord=order) for k in range(B.shape[0])]
            assert_almost_equal(n, expected)

            n = norm(B, ord=order, axis=(0, 2))
            expected = [norm(B[:, k,:], ord=order) for k in range(B.shape[1])]
            assert_almost_equal(n, expected)

            n = norm(B, ord=order, axis=(0, 1))
            expected = [norm(B[:,:, k], ord=order) for k in range(B.shape[2])]
            assert_almost_equal(n, expected)

    def test_bad_args(self):
        # Check that bad arguments raise the appropriate exceptions.

        A = array([[1, 2, 3], [4, 5, 6]], dtype=self.dt)
        B = np.arange(1, 25, dtype=self.dt).reshape(2, 3, 4)

        # Using `axis=<integer>` or passing in a 1-D array implies vector
        # norms are being computed, so also using `ord='fro'` raises a
        # ValueError.
        assert_raises(ValueError, norm, A, 'fro', 0)
        assert_raises(ValueError, norm, [3, 4], 'fro', None)

        # Similarly, norm should raise an exception when ord is any finite
        # number other than 1, 2, -1 or -2 when computing matrix norms.
        for order in [0, 3]:
            assert_raises(ValueError, norm, A, order, None)
            assert_raises(ValueError, norm, A, order, (0, 1))
            assert_raises(ValueError, norm, B, order, (1, 2))

        # Invalid axis
        assert_raises(ValueError, norm, B, None, 3)
        assert_raises(ValueError, norm, B, None, (2, 3))
        assert_raises(ValueError, norm, B, None, (0, 1, 2))

    def test_longdouble_norm(self):
        # Non-regression test: p-norm of longdouble would previously raise
        # UnboundLocalError.
        x = np.arange(10, dtype=np.longdouble)
        old_assert_almost_equal(norm(x, ord=3), 12.65, decimal=2)

    def test_intmin(self):
        # Non-regression test: p-norm of signed integer would previously do
        # float cast and abs in the wrong order.
        x = np.array([-2 ** 31], dtype=np.int32)
        old_assert_almost_equal(norm(x, ord=3), 2 ** 31, decimal=5)

    def test_complex_high_ord(self):
        # gh-4156
        d = np.empty((2,), dtype=np.clongdouble)
        d[0] = 6+7j
        d[1] = -6+7j
        res = 11.615898132184
        old_assert_almost_equal(np.linalg.norm(d, ord=3), res, decimal=10)
        d = d.astype(np.complex128)
        old_assert_almost_equal(np.linalg.norm(d, ord=3), res, decimal=9)
        d = d.astype(np.complex64)
        old_assert_almost_equal(np.linalg.norm(d, ord=3), res, decimal=5)


class TestNormDouble(_TestNorm):
    dt = np.double
    dec = 12


class TestNormSingle(_TestNorm):
    dt = np.float32
    dec = 6


class TestNormInt64(_TestNorm):
    dt = np.int64
    dec = 12


class TestMatrixRank(object):
    def test_matrix_rank(self):
        # Full rank matrix
        yield assert_equal, 4, matrix_rank(np.eye(4))
        # rank deficient matrix
        I=np.eye(4); I[-1, -1] = 0.
        yield assert_equal, matrix_rank(I), 3
        # All zeros - zero rank
        yield assert_equal, matrix_rank(np.zeros((4, 4))), 0
        # 1 dimension - rank 1 unless all 0
        yield assert_equal, matrix_rank([1, 0, 0, 0]), 1
        yield assert_equal, matrix_rank(np.zeros((4,))), 0
        # accepts array-like
        yield assert_equal, matrix_rank([1]), 1
        # greater than 2 dimensions raises error
        yield assert_raises, TypeError, matrix_rank, np.zeros((2, 2, 2))
        # works on scalar
        yield assert_equal, matrix_rank(1), 1


def test_reduced_rank():
    # Test matrices with reduced rank
    rng = np.random.RandomState(20120714)
    for i in range(100):
        # Make a rank deficient matrix
        X = rng.normal(size=(40, 10))
        X[:, 0] = X[:, 1] + X[:, 2]
        # Assert that matrix_rank detected deficiency
        assert_equal(matrix_rank(X), 9)
        X[:, 3] = X[:, 4] + X[:, 5]
        assert_equal(matrix_rank(X), 8)


class TestQR(object):

    def check_qr(self, a):
        # This test expects the argument `a` to be an ndarray or
        # a subclass of an ndarray of inexact type.
        a_type = type(a)
        a_dtype = a.dtype
        m, n = a.shape
        k = min(m, n)

        # mode == 'complete'
        q, r = linalg.qr(a, mode='complete')
        assert_(q.dtype == a_dtype)
        assert_(r.dtype == a_dtype)
        assert_(isinstance(q, a_type))
        assert_(isinstance(r, a_type))
        assert_(q.shape == (m, m))
        assert_(r.shape == (m, n))
        assert_almost_equal(dot(q, r), a)
        assert_almost_equal(dot(q.T.conj(), q), np.eye(m))
        assert_almost_equal(np.triu(r), r)

        # mode == 'reduced'
        q1, r1 = linalg.qr(a, mode='reduced')
        assert_(q1.dtype == a_dtype)
        assert_(r1.dtype == a_dtype)
        assert_(isinstance(q1, a_type))
        assert_(isinstance(r1, a_type))
        assert_(q1.shape == (m, k))
        assert_(r1.shape == (k, n))
        assert_almost_equal(dot(q1, r1), a)
        assert_almost_equal(dot(q1.T.conj(), q1), np.eye(k))
        assert_almost_equal(np.triu(r1), r1)

        # mode == 'r'
        r2 = linalg.qr(a, mode='r')
        assert_(r2.dtype == a_dtype)
        assert_(isinstance(r2, a_type))
        assert_almost_equal(r2, r1)

    def test_qr_empty(self):
        a = np.zeros((0, 2))
        assert_raises(linalg.LinAlgError, linalg.qr, a)

    def test_mode_raw(self):
        # The factorization is not unique and varies between libraries,
        # so it is not possible to check against known values. Functional
        # testing is a possibility, but awaits the exposure of more
        # of the functions in lapack_lite. Consequently, this test is
        # very limited in scope. Note that the results are in FORTRAN
        # order, hence the h arrays are transposed.
        a = array([[1, 2], [3, 4], [5, 6]], dtype=np.double)
        b = a.astype(np.single)

        # Test double
        h, tau = linalg.qr(a, mode='raw')
        assert_(h.dtype == np.double)
        assert_(tau.dtype == np.double)
        assert_(h.shape == (2, 3))
        assert_(tau.shape == (2,))

        h, tau = linalg.qr(a.T, mode='raw')
        assert_(h.dtype == np.double)
        assert_(tau.dtype == np.double)
        assert_(h.shape == (3, 2))
        assert_(tau.shape == (2,))

    def test_mode_all_but_economic(self):
        a = array([[1, 2], [3, 4]])
        b = array([[1, 2], [3, 4], [5, 6]])
        for dt in "fd":
            m1 = a.astype(dt)
            m2 = b.astype(dt)
            self.check_qr(m1)
            self.check_qr(m2)
            self.check_qr(m2.T)
            self.check_qr(matrix(m1))
        for dt in "fd":
            m1 = 1 + 1j * a.astype(dt)
            m2 = 1 + 1j * b.astype(dt)
            self.check_qr(m1)
            self.check_qr(m2)
            self.check_qr(m2.T)
            self.check_qr(matrix(m1))


def test_byteorder_check():
    # Byte order check should pass for native order
    if sys.byteorder == 'little':
        native = '<'
    else:
        native = '>'

    for dtt in (np.float32, np.float64):
        arr = np.eye(4, dtype=dtt)
        n_arr = arr.newbyteorder(native)
        sw_arr = arr.newbyteorder('S').byteswap()
        assert_equal(arr.dtype.byteorder, '=')
        for routine in (linalg.inv, linalg.det, linalg.pinv):
            # Normal call
            res = routine(arr)
            # Native but not '='
            assert_array_equal(res, routine(n_arr))
            # Swapped
            assert_array_equal(res, routine(sw_arr))


def test_generalized_raise_multiloop():
    # It should raise an error even if the error doesn't occur in the
    # last iteration of the ufunc inner loop

    invertible = np.array([[1, 2], [3, 4]])
    non_invertible = np.array([[1, 1], [1, 1]])

    x = np.zeros([4, 4, 2, 2])[1::2]
    x[...] = invertible
    x[0, 0] = non_invertible

    assert_raises(np.linalg.LinAlgError, np.linalg.inv, x)

def test_xerbla_override():
    # Check that our xerbla has been successfully linked in. If it is not,
    # the default xerbla routine is called, which prints a message to stdout
    # and may, or may not, abort the process depending on the LAPACK package.
    from nose import SkipTest

    try:
        pid = os.fork()
    except (OSError, AttributeError):
        # fork failed, or not running on POSIX
        raise SkipTest("Not POSIX or fork failed.")

    if pid == 0:
        # child; close i/o file handles
        os.close(1)
        os.close(0)
        # Avoid producing core files.
        import resource
        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
        # These calls may abort.
        try:
            np.linalg.lapack_lite.xerbla()
        except ValueError:
            pass
        except:
            os._exit(os.EX_CONFIG)

        try:
            a = np.array([[1.]])
            np.linalg.lapack_lite.dorgqr(
                1, 1, 1, a,
                0, # <- invalid value
                a, a, 0, 0)
        except ValueError as e:
            if "DORGQR parameter number 5" in str(e):
                # success
                os._exit(os.EX_OK)

        # Did not abort, but our xerbla was not linked in.
        os._exit(os.EX_CONFIG)
    else:
        # parent
        pid, status = os.wait()
        if os.WEXITSTATUS(status) != os.EX_OK or os.WIFSIGNALED(status):
            raise SkipTest('Numpy xerbla not linked in.')


if __name__ == "__main__":
    run_module_suite()
