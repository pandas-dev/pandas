# Tests numpy methods of <class 'function'>

import itertools
import math
import platform
from functools import partial
from itertools import product
from textwrap import dedent

import numpy as np

from numba import jit, njit, typeof
from numba.core import types
from numba.typed import List, Dict
from numba.np.numpy_support import numpy_version
from numba.core.errors import TypingError
from numba.core.config import IS_32BITS
from numba.core.utils import pysignature
from numba.np.extensions import cross2d
from numba.tests.support import (TestCase, MemoryLeakMixin,
                                 needs_blas, run_in_subprocess,
                                 skip_if_numpy_2, IS_NUMPY_2,
                                 IS_MACOS_ARM64)
import unittest


def sinc(x):
    return np.sinc(x)


def angle1(x):
    return np.angle(x)


def angle2(x, deg):
    return np.angle(x, deg)


def array_equal(a, b):
    return np.array_equal(a, b)


def intersect1d_2(a, b):
    return np.intersect1d(a, b)


def intersect1d_3(a, b, assume_unique=False):
    return np.intersect1d(a, b, assume_unique)


def append(arr, values, axis):
    return np.append(arr, values, axis=axis)


def count_nonzero(arr, axis):
    return np.count_nonzero(arr, axis=axis)


def delete(arr, obj):
    return np.delete(arr, obj)


def diff1(a):
    return np.diff(a)


def diff2(a, n):
    return np.diff(a, n)


def bincount1(a):
    return np.bincount(a)


def bincount2(a, w):
    return np.bincount(a, weights=w)


def bincount3(a, w=None, minlength=0):
    return np.bincount(a, w, minlength)


def searchsorted(a, v):
    return np.searchsorted(a, v)


def searchsorted_left(a, v):
    return np.searchsorted(a, v, side='left')


def searchsorted_right(a, v):
    return np.searchsorted(a, v, side='right')


def digitize(*args):
    return np.digitize(*args)


def histogram(*args):
    return np.histogram(*args)


def machar(*args):
    return np.MachAr()


def iscomplex(x):
    return np.iscomplex(x)


def iscomplexobj(x):
    return np.iscomplexobj(x)


def isscalar(x):
    return np.isscalar(x)


def isreal(x):
    return np.isreal(x)


def isrealobj(x):
    return np.isrealobj(x)


def isneginf(x, out=None):
    return np.isneginf(x, out)


def isposinf(x, out=None):
    return np.isposinf(x, out)


def isclose(a, b, rtol=1e-05, atol=1e-08, equal_nan=False):
    return np.isclose(a, b, rtol, atol, equal_nan)


def isnat(x):
    return np.isnat(x)


def iinfo(*args):
    return np.iinfo(*args)


def finfo(*args):
    return np.finfo(*args)


def finfo_machar(*args):
    return np.finfo(*args).machar


def fliplr(a):
    return np.fliplr(a)


def flipud(a):
    return np.flipud(a)


def flip(a):
    return np.flip(a)


def logspace2(start, stop):
    return np.logspace(start, stop)


def logspace3(start, stop, num=50):
    return np.logspace(start, stop, num=num)


def geomspace2(start, stop):
    return np.geomspace(start, stop)


def geomspace3(start, stop, num=50):
    return np.geomspace(start, stop, num=num)


def rot90(a):
    return np.rot90(a)


def rot90_k(a, k=1):
    return np.rot90(a, k)


def array_split(a, indices, axis=0):
    return np.array_split(a, indices, axis=axis)


def split(a, indices, axis=0):
    return np.split(a, indices, axis=axis)


def vsplit(a, ind_or_sec):
    return np.vsplit(a, ind_or_sec)


def hsplit(a, ind_or_sec):
    return np.hsplit(a, ind_or_sec)


def dsplit(a, ind_or_sec):
    return np.dsplit(a, ind_or_sec)


def correlate(a, v, mode="valid"):
    return np.correlate(a, v, mode=mode)


def convolve(a, v, mode="full"):
    return np.convolve(a, v, mode=mode)


def tri_n(N):
    return np.tri(N)


def tri_n_m(N, M=None):
    return np.tri(N, M)


def tri_n_k(N, k=0):
    return np.tri(N, k)


def tri_n_m_k(N, M=None, k=0):
    return np.tri(N, M, k)


def tril_m(m):
    return np.tril(m)


def tril_m_k(m, k=0):
    return np.tril(m, k)


def tril_indices_n(n):
    return np.tril_indices(n)


def tril_indices_n_k(n, k=0):
    return np.tril_indices(n, k)


def tril_indices_n_m(n, m=None):
    return np.tril_indices(n, m=m)


def tril_indices_n_k_m(n, k=0, m=None):
    return np.tril_indices(n, k, m)


def tril_indices_from_arr(arr):
    return np.tril_indices_from(arr)


def tril_indices_from_arr_k(arr, k=0):
    return np.tril_indices_from(arr, k)


def triu_m(m):
    return np.triu(m)


def triu_m_k(m, k=0):
    return np.triu(m, k)


def triu_indices_n(n):
    return np.triu_indices(n)


def triu_indices_n_k(n, k=0):
    return np.triu_indices(n, k)


def triu_indices_n_m(n, m=None):
    return np.triu_indices(n, m=m)


def triu_indices_n_k_m(n, k=0, m=None):
    return np.triu_indices(n, k, m)


def triu_indices_from_arr(arr):
    return np.triu_indices_from(arr)


def triu_indices_from_arr_k(arr, k=0):
    return np.triu_indices_from(arr, k)


def vander(x, N=None, increasing=False):
    return np.vander(x, N, increasing)


def partition(a, kth):
    return np.partition(a, kth)


def argpartition(a, kth):
    return np.argpartition(a, kth)


def cov(m, y=None, rowvar=True, bias=False, ddof=None):
    return np.cov(m, y, rowvar, bias, ddof)


def corrcoef(x, y=None, rowvar=True):
    return np.corrcoef(x, y, rowvar)


def ediff1d(ary, to_end=None, to_begin=None):
    return np.ediff1d(ary, to_end, to_begin)


def roll(a, shift):
    return np.roll(a, shift)


def asarray(a):
    return np.asarray(a)


def asarray_kws(a, dtype):
    return np.asarray(a, dtype=dtype)


def asfarray(a, dtype=np.float64):
    return np.asfarray(a, dtype=dtype)


def asfarray_default_kwarg(a):
    return np.asfarray(a)


def extract(condition, arr):
    return np.extract(condition, arr)


def np_trapz(y):
    return np.trapz(y)


def np_trapz_x(y, x):
    return np.trapz(y, x)


def np_trapz_dx(y, dx):
    return np.trapz(y, dx=dx)


def np_trapz_x_dx(y, x, dx):
    return np.trapz(y, x, dx)


def np_trapezoid(y):
    return np.trapezoid(y)


def np_trapezoid_x(y, x):
    return np.trapezoid(y, x)


def np_trapezoid_dx(y, dx):
    return np.trapezoid(y, dx=dx)


def np_trapezoid_x_dx(y, x, dx):
    return np.trapezoid(y, x, dx)


def np_allclose(a, b, rtol=1e-05, atol=1e-08, equal_nan=False):
    return np.allclose(a, b, rtol, atol, equal_nan)


def np_average(a, axis=None, weights=None):
    return np.average(a, axis=axis, weights=weights)


def interp(x, xp, fp):
    return np.interp(x, xp, fp)


def np_repeat(a, repeats):
    return np.repeat(a, repeats)


def array_repeat(a, repeats):
    return np.asarray(a).repeat(repeats)


def np_select(condlist, choicelist, default=0):
    return np.select(condlist, choicelist, default=default)


def np_select_defaults(condlist, choicelist):
    return np.select(condlist, choicelist)


def np_bartlett(M):
    return np.bartlett(M)


def np_blackman(M):
    return np.blackman(M)


def np_hamming(M):
    return np.hamming(M)


def np_hanning(M):
    return np.hanning(M)


def np_kaiser(M, beta):
    return np.kaiser(M, beta)


def np_cross(a, b):
    return np.cross(a, b)


def np_trim_zeros(a, trim='fb'):
    return np.trim_zeros(a, trim)


def nb_cross2d(a, b):
    return cross2d(a, b)


def flip_lr(a):
    return np.fliplr(a)


def flip_ud(a):
    return np.flipud(a)


def np_union1d(a, b):
    return np.union1d(a,b)


def np_asarray_chkfinite(a, dtype=None):
    return np.asarray_chkfinite(a, dtype)


def unwrap(p, discont=None, axis=-1, period=6.283185307179586):
    return np.unwrap(p, discont, axis, period=period)


def unwrap1(p):
    return np.unwrap(p)


def unwrap13(p, period):
    return np.unwrap(p, period=period)


def unwrap123(p, period, discont):
    return np.unwrap(p, period=period, discont=discont)


def array_contains(a, key):
    return key in a


def swapaxes(a, a1, a2):
    return np.swapaxes(a, a1, a2)


def nan_to_num(X, copy=True, nan=0.0):
    return np.nan_to_num(X, copy=copy, nan=nan)


def np_indices(dimensions):
    return np.indices(dimensions)


def diagflat1(v):
    return np.diagflat(v)


def diagflat2(v, k=0):
    return np.diagflat(v, k)


def np_setxor1d_2(a, b):
    return np.setxor1d(a, b)


def np_setxor1d_3(a, b, assume_unique=False):
    return np.setxor1d(a, b, assume_unique)


def np_setdiff1d_2(a, b):
    return np.setdiff1d(a, b)


def np_setdiff1d_3(a, b, assume_unique=False):
    return np.setdiff1d(a, b, assume_unique)


def np_in1d_2(a, b):
    return np.in1d(a, b)


def np_in1d_3a(a, b, assume_unique=False):
    return np.in1d(a, b, assume_unique=assume_unique)


def np_in1d_3b(a, b, invert=False):
    return np.in1d(a, b, invert=invert)


def np_in1d_4(a, b, assume_unique=False, invert=False):
    return np.in1d(a, b, assume_unique, invert)


def np_isin_2(a, b):
    return np.isin(a, b)


def np_isin_3a(a, b, assume_unique=False):
    return np.isin(a, b, assume_unique=assume_unique)


def np_isin_3b(a, b, invert=False):
    return np.isin(a, b, invert=invert)


def np_isin_4(a, b, assume_unique=False, invert=False):
    return np.isin(a, b, assume_unique, invert)


class TestNPFunctions(MemoryLeakMixin, TestCase):
    """
    Tests for various Numpy functions.
    """

    def setUp(self):
        super(TestNPFunctions, self).setUp()
        self.rnd = np.random.RandomState(42)

    def run_unary(self, pyfunc, x_types, x_values, func_extra_types=None,
                  func_extra_args=None, ignore_sign_on_zero=False, abs_tol=None,
                  **kwargs):
        """
        Runs tests for a unary function operating in the numerical real space.

        Parameters
        ----------
        pyfunc : a python function definition holding that calls the numpy
                 functions to be tested.
        x_types: the types of the values being tested, see numba.types
        x_values: the numerical values of the values to be tested
        func_extra_types: the types of additional arguments to the numpy
                          function
        func_extra_args:  additional arguments to the numpy function
        ignore_sign_on_zero: boolean as to whether to allow zero values
        with incorrect signs to be considered equal
        prec: the required precision match, see assertPreciseEqual

        Notes:
        ------
        x_types and x_values must have the same length

        """
        for tx, vx in zip(x_types, x_values):
            if func_extra_args is None:
                func_extra_types = func_extra_args = [()]
            for xtypes, xargs in zip(func_extra_types, func_extra_args):
                cfunc = njit((tx,) + xtypes,)(pyfunc)
                got = cfunc(vx, *xargs)
                expected = pyfunc(vx, *xargs)
                try:
                    scalty = tx.dtype
                except AttributeError:
                    scalty = tx
                prec = ('single'
                        if scalty in (types.float32, types.complex64)
                        else 'double')
                msg = 'for input %r with prec %r' % (vx, prec)
                self.assertPreciseEqual(got, expected,
                                        prec=prec,
                                        msg=msg,
                                        ignore_sign_on_zero=ignore_sign_on_zero,
                                        abs_tol=abs_tol, **kwargs)

    def test_sinc(self):
        """
        Tests the sinc() function.
        This test is purely to assert numerical computations are correct.
        """

        # Ignore sign of zeros, this will need masking depending on numpy
        # version once the fix to numpy complex division is in upstream
        # See: https://github.com/numpy/numpy/pull/6699
        isoz = True

        # Testing sinc(1.) leads to sin(pi)/pi, which is below machine
        # precision in practice on most machines. Small floating point
        # differences in sin() etc. may lead to large differences in the result
        # that are at a range that is inaccessible using standard width
        # floating point representations.
        # e.g. Assume float64 type.
        # sin(pi) ~= 1e-16, but should be zero
        # sin(pi)/pi ~= 1e-17, should be zero, error carried from above
        # float64 has log10(2^53)~=15.9 digits of precision and the magnitude
        # change in the alg is > 16  digits (1.0...0 -> 0.0...0),
        # so comparison via ULP is invalid.
        # We therefore opt to assume that values under machine precision are
        # equal in this case.
        tol = "eps"

        pyfunc = sinc

        def check(x_types, x_values, **kwargs):
            self.run_unary(pyfunc, x_types, x_values,
                           ignore_sign_on_zero=isoz, abs_tol=tol,
                           **kwargs)

        # real domain scalar context
        x_values = [1., -1., 0.0, -0.0, 0.5, -0.5, 5, -5, 5e-21, -5e-21]
        x_types = [types.float32, types.float64] * (len(x_values) // 2)
        check(x_types, x_values)

        # real domain vector context
        x_values = [np.array(x_values, dtype=np.float64)]
        x_types = [typeof(v) for v in x_values]
        check(x_types, x_values)

        # complex domain scalar context
        x_values = [1.+0j, -1+0j, 0.0+0.0j, -0.0+0.0j, 0+1j, 0-1j, 0.5+0.0j, # noqa
                    -0.5+0.0j, 0.5+0.5j, -0.5-0.5j, 5+5j, -5-5j,             # noqa
                    # the following are to test sin(x)/x for small x
                    5e-21+0j, -5e-21+0j, 5e-21j, +(0-5e-21j)                 # noqa
                    ]
        x_types = [types.complex64, types.complex128] * (len(x_values) // 2)
        check(x_types, x_values, ulps=2)

        # complex domain vector context
        x_values = [np.array(x_values, dtype=np.complex128)]
        x_types = [typeof(v) for v in x_values]
        check(x_types, x_values, ulps=2)

    def test_sinc_exceptions(self):
        pyfunc = sinc
        cfunc = jit(nopython=True)(pyfunc)

        with self.assertRaises(TypingError) as raises:
            cfunc('str')
        self.assertIn('Argument "x" must be a Number or array-like',
                      str(raises.exception))
        # Exceptions leak references
        self.disable_leak_check()

    def test_contains(self):
        def arrs():
            a_0 = np.arange(10, 50)
            k_0 = 20

            yield a_0, k_0

            a_1 = np.arange(6)
            k_1 = 10

            yield a_1, k_1

            single_val_a = np.asarray([20])
            k_in = 20
            k_out = 13

            yield single_val_a, k_in
            yield single_val_a, k_out

            empty_arr = np.asarray([])
            yield empty_arr, k_out

            # np scalars

            bool_arr = np.array([True, False])
            yield bool_arr, True
            yield bool_arr, k_0

            np.random.seed(2)
            float_arr = np.random.rand(10)
            np.random.seed(2)
            rand_k = np.random.rand()
            present_k = float_arr[0]

            yield float_arr, rand_k
            yield float_arr, present_k

            complx_arr = float_arr.view(np.complex128)
            yield complx_arr, complx_arr[0]
            yield complx_arr, rand_k

            np.random.seed(2)
            uint_arr = np.random.randint(10, size=15, dtype=np.uint8)
            yield uint_arr, 5
            yield uint_arr, 25

        pyfunc = array_contains

        cfunc = jit(nopython=True)(pyfunc)

        for arr, key in arrs():
            expected = pyfunc(arr, key)
            received = cfunc(arr, key)

            self.assertPreciseEqual(expected, received)

    def test_angle(self):
        """
        Tests the angle() function.
        This test is purely to assert numerical computations are correct.
        """
        pyfunc1 = angle1
        pyfunc2 = angle2

        def check(x_types, x_values):
            # angle(x)
            self.run_unary(pyfunc1, x_types, x_values)
            # angle(x, deg)
            xtra_values = [(True,), (False,)]
            xtra_types = [(types.bool_,)] * len(xtra_values)
            self.run_unary(pyfunc2, x_types, x_values,
                           func_extra_types=xtra_types,
                           func_extra_args=xtra_values,)

        # real domain scalar context
        x_values = [1., -1., 0.0, -0.0, 0.5, -0.5, 5, -5]
        x_types = [types.float32, types.float64] * (len(x_values) // 2 + 1)
        check(x_types, x_values)

        # real domain vector context
        x_values = [np.array(x_values, dtype=np.float64)]
        x_types = [typeof(v) for v in x_values]
        check(x_types, x_values)

        # complex domain scalar context
        x_values = [1.+0j, -1+0j, 0.0+0.0j, -0.0+0.0j, 1j, -1j, 0.5+0.0j, # noqa
                    -0.5+0.0j, 0.5+0.5j, -0.5-0.5j, 5+5j, -5-5j]          # noqa
        x_types = [types.complex64, types.complex128] * (len(x_values) // 2 + 1)
        check(x_types, x_values)

        # complex domain vector context
        x_values = np.array(x_values)
        x_types = [types.complex64, types.complex128]
        check(x_types, x_values)

    def test_angle_return_type(self):
        # see issue #8949
        def numba_angle(x):
            r = np.angle(x)
            return r.dtype

        pyfunc = numba_angle
        x_values = [1., -1., 1. + 0j, -5 - 5j]
        x_types = ['f4', 'f8', 'c8', 'c16']
        for val, typ in zip(x_values, x_types):
            x = np.array([val], dtype=typ)
            cfunc = jit(nopython=True)(pyfunc)
            expected = pyfunc(x)
            got = cfunc(x)
            self.assertEqual(expected, got)

    def test_angle_exceptions(self):
        pyfunc = angle1
        cfunc = jit(nopython=True)(pyfunc)

        with self.assertRaises(TypingError) as raises:
            cfunc('hello')
        self.assertIn('Argument "z" must be a complex or Array[complex]',
                      str(raises.exception))
        # Exceptions leak references
        self.disable_leak_check()

    def test_array_equal(self):
        def arrays():
            yield np.array([]), np.array([])
            yield np.array([1, 2]), np.array([1, 2])
            yield np.array([]), np.array([1])
            x = np.arange(10).reshape(5, 2)
            x[1][1] = 30
            yield np.arange(10).reshape(5, 2), x
            yield x, x
            yield (1, 2, 3), (1, 2, 3)
            yield 2, 2
            yield 3, 2
            yield True, True
            yield True, False
            yield True, 2
            yield True, 1
            yield False, 0

        pyfunc = array_equal
        cfunc = jit(nopython=True)(pyfunc)

        for arr, obj in arrays():
            expected = pyfunc(arr, obj)
            got = cfunc(arr, obj)
            self.assertPreciseEqual(expected, got)

    def test_array_equal_exception(self):
        pyfunc = array_equal
        cfunc = jit(nopython=True)(pyfunc)

        with self.assertRaises(TypingError) as raises:
            cfunc(np.arange(3 * 4).reshape(3, 4), None)
        self.assertIn(
            'Both arguments to "array_equals" must be array-like',
            str(raises.exception)
        )

    def test_intersect1d_2(self):

        def arrays():
            yield (List.empty_list(types.float64),
                   List.empty_list(types.float64))  # two empty arrays
            yield [1], List.empty_list(types.float64)  # empty right
            yield List.empty_list(types.float64), [1]  # empty left
            yield [1], [2]  # singletons no intersection
            yield [1], [1]  # singletons one intersection
            yield [1, 2], [1]
            yield [1, 2, 2], [2, 2]
            yield [1, 2], [2, 1]
            yield [1, 2, 3], [1, 2, 3]
            # from numpy:
            # https://github.com/numpy/numpy/blob/b0371ef240560e78b651a5d7c9407ae3212a3d56/numpy/lib/tests/test_arraysetops.py#L17 # noqa: E501
            yield [5, 7, 1, 2], [2, 4, 3, 1, 5]
            yield [5, 5, 7, 1, 2], [2, 1, 4, 3, 3, 1, 5]

        pyfunc = intersect1d_2
        cfunc = jit(nopython=True)(pyfunc)

        for a, b in arrays():
            # a = np.array(a)
            # b = np.array(b)
            if isinstance(a, list):
                a = List(a)
            if isinstance(b, list):
                b = List(b)
            expected = pyfunc(a, b)
            got = cfunc(a, b)
            self.assertPreciseEqual(expected, got)

    def test_intersect1d_3(self):

        def arrays():
            yield (List.empty_list(types.float64),
                   List.empty_list(types.float64))  # two empty arrays
            yield [1], List.empty_list(types.float64)  # empty right
            yield List.empty_list(types.float64), [1]  # empty left
            yield [1], [2]  # singletons no intersection
            yield [1], [1]  # singletons one intersection
            yield [1, 2], [1]
            yield [1, 2, 2], [2, 2]
            yield [1, 2], [2, 1]
            yield [1, 2, 3], [1, 2, 3]
            # from numpy:
            # https://github.com/numpy/numpy/blob/b0371ef240560e78b651a5d7c9407ae3212a3d56/numpy/lib/tests/test_arraysetops.py#L17 # noqa: E501
            yield [5, 7, 1, 2], [2, 4, 3, 1, 5]
            yield [5, 5, 7, 1, 2], [2, 1, 4, 3, 3, 1, 5]

        pyfunc = intersect1d_3
        cfunc = jit(nopython=True)(pyfunc)

        for a, b in arrays():
            if isinstance(a, list):
                a = List(a)
            if isinstance(b, list):
                b = List(b)
            expected = pyfunc(a, b, assume_unique=False)
            got = cfunc(a, b, assume_unique=False)
            self.assertPreciseEqual(expected, got)
            if len(np.unique(a)) == len(a) and len(np.unique(b)) == len(b):
                expected = pyfunc(a, b, assume_unique=True)
                got = cfunc(a, b, assume_unique=True)
                self.assertPreciseEqual(expected, got)

    def test_intersect1d_errors(self):
        np_pyfunc = intersect1d_3
        np_nbfunc = njit(np_pyfunc)

        a = np.array([1])
        b = np.array([2])
        self.disable_leak_check()
        with self.assertRaises(TypingError):
            np_nbfunc(a, b, "foo")
        with self.assertRaises(TypingError):
            np_nbfunc("foo", b, True)
        with self.assertRaises(TypingError):
            np_nbfunc(a, "foo", True)

    def test_count_nonzero(self):

        def arrays():
            yield np.array([]), None
            yield np.zeros(10), None
            yield np.arange(10), None
            yield np.ones(10, dtype=np.bool_), 0
            yield np.arange(3 * 4 * 5).reshape(3, 4, 5), None
            yield np.arange(3 * 4).reshape(3, 4), 0
            yield np.arange(3 * 4).reshape(3, 4), 1

        pyfunc = count_nonzero
        cfunc = jit(nopython=True)(pyfunc)

        for arr, axis in arrays():
            expected = pyfunc(arr, axis)
            got = cfunc(arr, axis)
            self.assertPreciseEqual(expected, got)

    def test_np_append(self):
        def arrays():
            yield 2, 2, None
            yield np.arange(10), 3, None
            yield np.arange(10), np.arange(3), None
            yield np.arange(10).reshape(5, 2), np.arange(3), None
            yield np.array([[1, 2, 3], [4, 5, 6]]), np.array([[7, 8, 9]]), 0
            arr = np.array([[1, 2, 3], [4, 5, 6]])
            yield arr, arr, 1

        pyfunc = append
        cfunc = jit(nopython=True)(pyfunc)

        for arr, obj, axis in arrays():
            expected = pyfunc(arr, obj, axis)
            got = cfunc(arr, obj, axis)
            self.assertPreciseEqual(expected, got)

    def test_np_append_exceptions(self):
        pyfunc = append
        cfunc = jit(nopython=True)(pyfunc)
        arr = np.array([[1, 2, 3], [4, 5, 6]])
        values = np.array([[7, 8, 9]])
        axis = 0

        # first argument must be array-like
        with self.assertRaises(TypingError) as raises:
            cfunc(None, values, axis)
        self.assertIn(
            'The first argument "arr" must be array-like',
            str(raises.exception)
        )

        # second argument must also be array-like
        with self.assertRaises(TypingError) as raises:
            cfunc(arr, None, axis)
        self.assertIn(
            'The second argument "values" must be array-like',
            str(raises.exception)
        )

        # third argument must be either nonelike or an integer
        with self.assertRaises(TypingError) as raises:
            cfunc(arr, values, axis=0.0)
        self.assertIn(
            'The third argument "axis" must be an integer',
            str(raises.exception)
        )
        # Exceptions leak references
        self.disable_leak_check()

    def test_delete(self):

        def arrays():
            # array, obj
            #
            # an array-like type
            yield [1, 2, 3, 4, 5], 3
            yield [1, 2, 3, 4, 5], [2, 3]
            # 1d array, scalar
            yield np.arange(10), 3
            yield np.arange(10), -3 # Negative obj
            # 1d array, list
            yield np.arange(10), [3, 5, 6]
            yield np.arange(10), [2, 3, 4, 5]
            # 3d array, scalar
            yield np.arange(3 * 4 * 5).reshape(3, 4, 5), 2
            # 3d array, list
            yield np.arange(3 * 4 * 5).reshape(3, 4, 5), [5, 30, 27, 8]
            # slices
            yield [1, 2, 3, 4], slice(1, 3, 1)
            yield np.arange(10), slice(10)

        pyfunc = delete
        cfunc = jit(nopython=True)(pyfunc)

        for arr, obj in arrays():
            expected = pyfunc(arr, obj)
            got = cfunc(arr, obj)
            self.assertPreciseEqual(expected, got)

    def test_delete_exceptions(self):
        pyfunc = delete
        cfunc = jit(nopython=True)(pyfunc)
        self.disable_leak_check()

        with self.assertRaises(TypingError) as raises:
            cfunc([1, 2], 3.14)
        self.assertIn(
            'obj should be of Integer dtype',
            str(raises.exception)
        )

        with self.assertRaises(TypingError) as raises:
            cfunc(np.arange(10), [3.5, 5.6, 6.2])
        self.assertIn(
            'obj should be of Integer dtype',
            str(raises.exception)
        )

        with self.assertRaises(TypingError) as raises:
            cfunc(2, 3)
        self.assertIn(
            'arr must be either an Array or a Sequence',
            str(raises.exception)
        )

        with self.assertRaises(IndexError) as raises:
            cfunc([1, 2], 3)
        self.assertIn(
            'obj must be less than the len(arr)',
            str(raises.exception),
        )
        # Exceptions leak references
        self.disable_leak_check()

    def diff_arrays(self):
        """
        Some test arrays for np.diff()
        """
        a = np.arange(12) ** 3
        yield a
        b = a.reshape((3, 4))
        yield b
        c = np.arange(24).reshape((3, 2, 4)) ** 3
        yield c

    def test_diff1(self):
        pyfunc = diff1
        cfunc = jit(nopython=True)(pyfunc)
        for arr in self.diff_arrays():
            expected = pyfunc(arr)
            got = cfunc(arr)
            self.assertPreciseEqual(expected, got)

        # 0-dim array
        a = np.array(42)
        with self.assertTypingError():
            cfunc(a)

    def test_diff2(self):
        pyfunc = diff2
        cfunc = jit(nopython=True)(pyfunc)
        for arr in self.diff_arrays():
            size = arr.shape[-1]
            for n in (0, 1, 2, 3, size - 1, size, size + 1, 421):
                expected = pyfunc(arr, n)
                got = cfunc(arr, n)
                self.assertPreciseEqual(expected, got)

    def test_diff2_exceptions(self):
        pyfunc = diff2
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        # 0-dim array
        arr = np.array(42)
        with self.assertTypingError():
            cfunc(arr, 1)

        # Invalid `n`
        arr = np.arange(10)
        for n in (-1, -2, -42):
            with self.assertRaises(ValueError) as raises:
                cfunc(arr, n)
            self.assertIn("order must be non-negative", str(raises.exception))
        # Exceptions leak references
        self.disable_leak_check()

    def test_isscalar(self):
        def values():
            yield 3
            yield np.asarray([3])
            yield (3,)
            yield 3j
            yield 'numba'
            yield int(10)
            yield np.int16(12345)
            yield 4.234
            yield True
            yield None
            yield np.timedelta64(10, 'Y')
            yield np.datetime64('nat')
            yield np.datetime64(1, 'Y')

        pyfunc = isscalar
        cfunc = jit(nopython=True)(pyfunc)
        for x in values():
            expected = pyfunc(x)
            got = cfunc(x)
            self.assertEqual(expected, got, x)

    def test_isobj_functions(self):
        def values():
            yield 1
            yield 1 + 0j
            yield np.asarray([3, 1 + 0j, True])
            yield "hello world"

        @jit(nopython=True)
        def optional_fn(x, cond, cfunc):
            y = x if cond else None
            return cfunc(y)

        pyfuncs = [iscomplexobj, isrealobj]
        for pyfunc in pyfuncs:
            cfunc = jit(nopython=True)(pyfunc)
            for x in values():
                expected = pyfunc(x)
                got = cfunc(x)
                self.assertEqual(expected, got)

                # optional type
                expected_optional = optional_fn.py_func(x, True, pyfunc)
                got_optional = optional_fn(x, True, cfunc)
                self.assertEqual(expected_optional, got_optional)

                # none type
                expected_none = optional_fn.py_func(x, False, pyfunc)
                got_none = optional_fn(x, False, cfunc)
                self.assertEqual(expected_none, got_none)

            self.assertEqual(len(cfunc.signatures), 8)

    def test_is_real_or_complex(self):
        def values():
            yield np.array([1 + 1j, 1 + 0j, 4.5, 3, 2, 2j])
            yield np.array([1, 2, 3])
            yield 3
            yield 12j
            yield 1 + 4j
            yield 10 + 0j
            yield (1 + 4j, 2 + 0j)
            yield np.array([[1, 2], [3, 4], [5, 6], [7, 8]])

        pyfuncs = [iscomplex, isreal]
        for pyfunc in pyfuncs:
            cfunc = jit(nopython=True)(pyfunc)
            for x in values():
                expected = pyfunc(x)
                got = cfunc(x)
                self.assertPreciseEqual(expected, got)

    def test_isneg_or_ispos_inf(self):
        def values():
            yield -np.inf, None
            yield np.inf, None
            yield np.inf, None
            yield np.asarray([-np.inf, 0., np.inf]), None
            yield -np.inf, np.zeros(1, dtype=np.bool_)
            yield np.inf, np.zeros(1, dtype=np.bool_)
            yield np.inf, np.zeros(1, dtype=np.bool_)
            yield -np.inf, np.empty(12)
            yield np.asarray([-np.inf, 0., np.inf]), np.zeros(3, dtype=np.bool_)

        pyfuncs = [isneginf, isposinf]
        for pyfunc in pyfuncs:
            cfunc = jit(nopython=True)(pyfunc)
            for x, out in values():
                expected = pyfunc(x, out)
                got = cfunc(x, out)
                self.assertPreciseEqual(expected, got)

    def test_isclose(self):
        rtol = 1e-5
        atol = 1e-8
        arr = np.array([100, 1000])
        aran = np.arange(8).reshape((2, 2, 2))
        kw = {'rtol': rtol, 'atol': atol}

        def values():
            yield 1e10, 1.00001e10, {}
            yield 1e10, np.nan, {}
            yield np.array([1e-8, 1e-7]), np.array([0.0, 0.0]), {}
            yield np.array([1e10, 1e-7]), np.array([1.00001e10, 1e-8]), {}
            yield np.array([1e10, 1e-8]), np.array([1.00001e10, 1e-9]), {}
            yield np.array([1e10, 1e-8]), np.array([1.0001e10, 1e-9]), {}
            yield np.array([1.0, np.nan]), np.array([1.0, np.nan]), {}
            yield np.array([1.0, np.nan]), np.array([1.0, np.nan]), {'equal_nan': True}  # noqa
            yield np.array([np.nan, np.nan]), np.array([1.0, np.nan]), {'equal_nan': True}  # noqa
            yield np.array([1e-100, 1e-7]), np.array([0.0, 0.0]), {'atol': 0.0}
            yield np.array([1e-10, 1e-10]), np.array([1e-20, 0.0]), {}
            yield np.array([1e-10, 1e-10]), np.array([1e-20, 0.999999e-10]), {'atol': 0.0}  # noqa
            yield np.array([1, np.inf, 2]), np.array([3, np.inf, 4]), kw
            yield np.array([atol, np.inf, -np.inf, np.nan]), np.array([0]), kw
            yield np.array([atol, np.inf, -np.inf, np.nan]), 0, kw
            yield 0, np.array([atol, np.inf, -np.inf, np.nan]), kw

            # tests taken from
            # https://github.com/numpy/numpy/blob/aac965af6032b69d5cb515ad785cc9a331e816f4/numpy/core/tests/test_numeric.py#L2298-L2335  # noqa: E501

            # all close tests
            yield np.array([0, 1]), np.array([1, 0]), kw
            yield arr, arr, kw
            yield np.array([1]), np.array([1 + rtol + atol]), kw
            yield arr, arr + arr * rtol, kw
            yield arr, arr + arr * rtol + atol, kw
            yield aran, aran + aran * rtol, kw
            yield np.inf, np.inf, kw
            yield -np.inf, np.inf, kw
            yield np.inf, np.array([np.inf]), kw
            yield np.array([np.inf, -np.inf]), np.array([np.inf, -np.inf]), kw

            # none close tests
            yield np.array([np.inf, 0]), np.array([1, np.inf]), kw
            yield np.array([np.inf, -np.inf]), np.array([1, 0]), kw
            yield np.array([np.inf, np.inf]), np.array([1, -np.inf]), kw
            yield np.array([np.inf, np.inf]), np.array([1, 0]), kw
            yield np.array([np.nan, 0]), np.array([np.nan, -np.inf]), kw
            yield np.array([atol * 2]), np.array([0]), kw
            yield np.array([1]), np.array([1 + rtol + atol * 2]), kw
            yield aran, aran + rtol * 1.1 * aran + atol * 1.1, kw
            yield np.array(np.array([np.inf, 1])), np.array(np.array([0, np.inf])), kw  # noqa

            # some close tests
            yield np.array([np.inf, 0]), np.array([atol * 2, atol * 2]), kw
            yield np.array([np.inf, 0]), np.array([np.inf, atol * 2]), kw
            yield np.array([atol, 1, 1e6 * (1 + 2 * rtol) + atol]), np.array([0, np.nan, 1e6]), kw  # noqa
            yield np.arange(3), np.array([0, 1, 2.1]), kw
            yield np.nan, np.array([np.nan, np.nan, np.nan]), kw
            yield np.array([0]), np.array([atol, np.inf, -np.inf, np.nan]), kw
            yield 0, np.array([atol, np.inf, -np.inf, np.nan]), kw

        pyfunc = isclose
        cfunc = jit(nopython=True)(pyfunc)
        for a, b, kwargs in values():
            expected = pyfunc(a, b, **kwargs)
            got = cfunc(a, b, **kwargs)
            if isinstance(expected, np.bool_):
                self.assertEqual(expected, got)
            else:
                self.assertTrue(np.array_equal(expected, got))

    def isclose_exception(self):
        pyfunc = isclose
        cfunc = jit(nopython=True)(pyfunc)
        inps = [
            (np.asarray([1e10, 1e-9, np.nan]),
             np.asarray([1.0001e10, 1e-9]),
             1e-05, 1e-08, False,
             "shape mismatch: objects cannot be broadcast to a single shape",
             ValueError),
            ('hello', 3, False, 1e-08, False,
             'The first argument "a" must be array-like',
             TypingError),
            (3, 'hello', False, 1e-08, False,
             'The second argument "b" must be array-like',
             TypingError),
            (2, 3, False, 1e-08, False,
             'The third argument "rtol" must be a floating point',
             TypingError),
            (2, 3, 1e-05, False, False,
             'The fourth argument "atol" must be a floating point',
             TypingError),
            (2, 3, 1e-05, 1e-08, 1,
             'The fifth argument "equal_nan" must be a boolean',
             TypingError),
        ]

        for a, b, rtol, atol, equal_nan, exc_msg, exc in inps:
            with self.assertRaisesRegex(exc, exc_msg):
                cfunc(a, b, rtol, atol, equal_nan)

    def bincount_sequences(self):
        """
        Some test sequences for np.bincount()
        """
        a = [1, 2, 5, 2, 3, 20]
        b = np.array([5, 8, 42, 5])
        c = self.rnd.randint(0, 100, size=300).astype(np.int8)
        return (a, b, c)

    def test_bincount1(self):
        pyfunc = bincount1
        cfunc = jit(nopython=True)(pyfunc)
        for seq in self.bincount_sequences():
            expected = pyfunc(seq)
            got = cfunc(seq)
            self.assertPreciseEqual(expected, got)

    def test_bincount1_exceptions(self):
        pyfunc = bincount1
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        # Negative input
        with self.assertRaises(ValueError) as raises:
            cfunc([2, -1])
        self.assertIn("first argument must be non-negative",
                      str(raises.exception))
        # Exceptions leak references
        self.disable_leak_check()

    def test_bincount2(self):
        pyfunc = bincount2
        cfunc = jit(nopython=True)(pyfunc)
        for seq in self.bincount_sequences():
            w = [math.sqrt(x) - 2 for x in seq]
            # weights as list, then array, mixed types, check upcast is ok
            for weights in (w, np.array(w), seq, np.array(seq)):
                expected = pyfunc(seq, weights)
                got = cfunc(seq, weights)
                self.assertPreciseEqual(expected, got)

    def test_bincount2_exceptions(self):
        pyfunc = bincount2
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        # Negative input
        with self.assertRaises(ValueError) as raises:
            cfunc([2, -1], [0, 0])
        self.assertIn("first argument must be non-negative",
                      str(raises.exception))

        # Mismatching input sizes
        with self.assertRaises(ValueError) as raises:
            cfunc([2, -1], [0])
        self.assertIn("weights and list don't have the same length",
                      str(raises.exception))

    def test_bincount3(self):
        pyfunc = bincount3
        cfunc = jit(nopython=True)(pyfunc)
        for seq in self.bincount_sequences():
            a_max = max(seq)
            # Length should be a_max in the first case, minlength in the second
            for minlength in (a_max, a_max + 2):
                expected = pyfunc(seq, None, minlength)
                got = cfunc(seq, None, minlength)
                self.assertEqual(len(expected), len(got))
                self.assertPreciseEqual(expected, got)

    def test_bincount3_exceptions(self):
        pyfunc = bincount3
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        # Negative input
        with self.assertRaises(ValueError) as raises:
            cfunc([2, -1], [0, 0])
        self.assertIn("first argument must be non-negative",
                      str(raises.exception))

        # Negative minlength
        with self.assertRaises(ValueError) as raises:
            cfunc([17, 38], None, -1)
        self.assertIn("'minlength' must not be negative",
                      str(raises.exception))

    def test_searchsorted(self):
        pyfunc = searchsorted
        cfunc = jit(nopython=True)(pyfunc)

        pyfunc_left = searchsorted_left
        cfunc_left = jit(nopython=True)(pyfunc_left)

        pyfunc_right = searchsorted_right
        cfunc_right = jit(nopython=True)(pyfunc_right)

        def check(a, v):
            expected = pyfunc(a, v)
            got = cfunc(a, v)
            self.assertPreciseEqual(expected, got)

            expected = pyfunc_left(a, v)
            got = cfunc_left(a, v)
            self.assertPreciseEqual(expected, got)

            expected = pyfunc_right(a, v)
            got = cfunc_right(a, v)
            self.assertPreciseEqual(expected, got)

        # First with integer values (no NaNs)
        bins = np.arange(5) ** 2
        values = np.arange(20) - 1

        for a in (bins, list(bins)):
            # Scalar values
            for v in values:
                check(a, v)
            # Array values
            for v in (values, values.reshape((4, 5))):
                check(a, v)
            # Sequence values
            check(a, list(values))

        # Second with float values (including NaNs)
        bins = np.float64(list(bins) + [float('nan')] * 7) / 2.0
        values = np.arange(20) - 0.5

        for a in (bins, list(bins)):
            # Scalar values
            for v in values:
                check(a, v)
            # Array values
            for v in (values, values.reshape((4, 5))):
                check(a, v)
            # Sequence values
            check(a, list(values))

        # nonsense value for 'side' raises TypingError
        def bad_side(a, v):
            return np.searchsorted(a, v, side='nonsense')
        cfunc = jit(nopython=True)(bad_side)
        with self.assertTypingError():
            cfunc([1,2], 1)

        # non-constant value for 'side' raises TypingError
        def nonconst_side(a, v, side='left'):
            return np.searchsorted(a, v, side=side)
        cfunc = jit(nopython=True)(nonconst_side)
        with self.assertTypingError():
            cfunc([1,2], 1, side='right')

        # Test unordered values
        a = np.array([1, 2, 0])
        v = np.array(
            [
                [5, 4],
                [6, 7],
                [2, 1],
                [0, 3],
            ]
        )
        check(a, v)

        a = np.array([9, 1, 4, 2, 0, 3, 7, 6, 8])
        v = np.array(
            [
                [5, 10],
                [10, 5],
                [-1, 5],
            ]
        )
        check(a, v)

    def test_searchsorted_supplemental(self):
        pyfunc = searchsorted
        cfunc = jit(nopython=True)(pyfunc)

        pyfunc_left = searchsorted_left
        cfunc_left = jit(nopython=True)(pyfunc_left)

        pyfunc_right = searchsorted_right
        cfunc_right = jit(nopython=True)(pyfunc_right)

        def check(a, v):
            expected = pyfunc(a, v)
            got = cfunc(a, v)
            self.assertPreciseEqual(expected, got)

            expected = pyfunc_left(a, v)
            got = cfunc_left(a, v)
            self.assertPreciseEqual(expected, got)

            expected = pyfunc_right(a, v)
            got = cfunc_right(a, v)
            self.assertPreciseEqual(expected, got)

        element_pool = list(range(-5, 50))
        element_pool += [np.nan] * 5 + [np.inf] * 3 + [-np.inf] * 3

        for i in range(1000):
            sample_size = self.rnd.choice([5, 10, 25])

            # `a` and `v` not sorted; either may have repeating values
            a = self.rnd.choice(element_pool, sample_size)
            v = self.rnd.choice(element_pool, sample_size + (i % 3 - 1))

            # output should match numpy regardless of whether `a` is sorted
            check(a, v)
            check(np.sort(a), v)

        ones = np.ones(5)
        nans = np.full(len(ones), fill_value=np.nan)
        check(ones, ones)

        # `a` and / or `v` full of nans
        check(ones, nans)
        check(nans, ones)
        check(nans, nans)

        # `v` is zero size
        a = np.arange(1)
        v = np.arange(0)
        check(a, v)

        # `a` and `v` booleans
        a = np.array([False, False, True, True])
        v = np.array([False, True])
        check(a, v)

        # `v` is a (scalar) boolean
        a = [1, 2, 3]
        v = True
        check(a, v)

        # `a` and `v` arrays of strings
        a = np.array(['1', '2', '3'])
        v = np.array(['2', '4'])
        check(a, v)

    def test_searchsorted_complex(self):
        pyfunc = searchsorted
        cfunc = jit(nopython=True)(pyfunc)

        pyfunc_left = searchsorted_left
        cfunc_left = jit(nopython=True)(pyfunc_left)

        pyfunc_right = searchsorted_right
        cfunc_right = jit(nopython=True)(pyfunc_right)

        def check(a, v):
            expected = pyfunc(a, v)
            got = cfunc(a, v)
            self.assertPreciseEqual(expected, got)

            expected = pyfunc_left(a, v)
            got = cfunc_left(a, v)
            self.assertPreciseEqual(expected, got)

            expected = pyfunc_right(a, v)
            got = cfunc_right(a, v)
            self.assertPreciseEqual(expected, got)

        pool = [0, 1, np.nan]
        element_pool = [complex(*c) for c in itertools.product(pool, pool)]

        for i in range(100):
            sample_size = self.rnd.choice([3, 5, len(element_pool)])

            # `a` and `v` not sorted; either may have repeating values
            a = self.rnd.choice(element_pool, sample_size)
            v = self.rnd.choice(element_pool, sample_size + (i % 3 - 1))

            # output should match numpy regardless of whether `a` is sorted
            check(a, v)
            check(np.sort(a), v)

        # check type promotion (a complex; v not so much)
        check(a=np.array(element_pool), v=np.arange(2))

    def test_digitize(self):
        pyfunc = digitize
        cfunc = jit(nopython=True)(pyfunc)

        def check(*args):
            expected = pyfunc(*args)
            got = cfunc(*args)
            self.assertPreciseEqual(expected, got)

        values = np.float64((0, 0.99, 1, 4.4, 4.5, 7, 8, 9, 9.5,
                             float('inf'), float('-inf'), float('nan')))
        assert len(values) == 12
        self.rnd.shuffle(values)

        bins1 = np.float64([1, 3, 4.5, 8])
        bins2 = np.float64([1, 3, 4.5, 8, float('inf'), float('-inf')])
        bins3 = np.float64([1, 3, 4.5, 8, float('inf'), float('-inf')]
                           + [float('nan')] * 10)

        all_bins = [bins1, bins2, bins3]
        xs = [values, values.reshape((3, 4))]

        # 2-ary digitize()
        for bins in all_bins:
            bins.sort()
            for x in xs:
                check(x, bins)
                check(x, bins[::-1])

        # 3-ary digitize()
        for bins in all_bins:
            for right in (True, False):
                check(values, bins, right)
                check(values, bins[::-1], right)

        # Sequence input
        check(list(values), bins1)

        # per https://github.com/numba/numba/issues/8768
        check(np.array([np.nan, 1]), np.array([1.5, np.nan]))

    def test_digitize_non_monotonic_bins(self):
        # Exceptions leak references
        self.disable_leak_check()

        pyfunc = digitize
        cfunc = jit(nopython=True)(pyfunc)

        def check_error(*args):
            for fn in (pyfunc, cfunc):
                with self.assertRaises(ValueError) as raises:
                    fn(*args)

                msg = 'bins must be monotonically increasing or decreasing'
                self.assertIn(msg, str(raises.exception))

        x = np.array([np.nan, 1])
        bins = np.array([np.nan, 1.5, 2.3, np.nan])
        check_error(x, bins)

        x = [-1, 0, 1, 2]
        bins = [0, 0, 1, 0]
        check_error(x, bins)

        bins = [1, 1, 0, 1]
        check_error(x, bins)

    def test_digitize_supplemental(self):
        # inspired by the tests in
        # https://github.com/numpy/numpy/blob/a277f62/numpy/lib/tests/test_function_base.py

        pyfunc = digitize
        cfunc = jit(nopython=True)(pyfunc)

        def check(*args):
            expected = pyfunc(*args)
            got = cfunc(*args)
            self.assertPreciseEqual(expected, got)

        # forward
        x = np.arange(-6, 5)
        bins = np.arange(-5, 5)
        check(x, bins)

        # reverse
        x = np.arange(5, -6, -1)
        bins = np.arange(5, -5, -1)
        check(x, bins)

        # random
        x = self.rnd.rand(10)
        bins = np.linspace(x.min(), x.max(), 10)
        check(x, bins)

        # right_basic
        x = [1, 5, 4, 10, 8, 11, 0]
        bins = [1, 5, 10]
        check(x, bins)

        # right_open
        x = np.arange(-6, 5)
        bins = np.arange(-6, 4)
        check(x, bins, True)

        # right_open_reverse
        x = np.arange(5, -6, -1)
        bins = np.arange(4, -6, -1)
        check(x, bins, True)

        # right_open_random
        x = self.rnd.rand(10)
        bins = np.linspace(x.min(), x.max(), 10)
        check(x, bins, True)

        # monotonic
        x = [-1, 0, 1, 2]
        bins = [0, 0, 1]
        check(x, bins)

        bins = [1, 1, 0]
        check(x, bins)

        bins = [1, 1, 1, 1]
        check(x, bins)

        # large_integers_increasing
        x = 2 ** 54  # loses precision in a float
        check([x], [x - 1, x + 1])

    def test_digitize_raise_if_x_complex(self):
        # Exceptions leak references
        self.disable_leak_check()

        pyfunc = digitize
        cfunc = jit(nopython=True)(pyfunc)

        x = np.array([1 + 1j])
        y = np.array([1., 3., 4.5, 8.])
        msg = 'x may not be complex'

        for func in pyfunc, cfunc:
            with self.assertTypingError() as raises:
                func(x, y)
                self.assertIn(msg, str(raises.exception))

    def test_histogram(self):
        pyfunc = histogram
        cfunc = jit(nopython=True)(pyfunc)

        def check(*args):
            pyhist, pybins = pyfunc(*args)
            chist, cbins = cfunc(*args)
            self.assertPreciseEqual(pyhist, chist)
            # There can be a slight discrepancy in the linspace() result
            # when `bins` is an integer...
            self.assertPreciseEqual(pybins, cbins, prec='double', ulps=2)

        def check_values(values):
            # Explicit bins array
            # (note Numpy seems to not support NaN bins)
            bins = np.float64([1, 3, 4.5, 8])
            check(values, bins)
            check(values.reshape((3, 4)), bins)

            # Explicit number of bins
            check(values, 7)

            # Explicit number of bins and bins range
            check(values, 7, (1.0, 13.5))

            # Implicit bins=10
            check(values)

        values = np.float64((0, 0.99, 1, 4.4, 4.5, 7, 8,
                             9, 9.5, 42.5, -1.0, -0.0))
        assert len(values) == 12
        self.rnd.shuffle(values)

        check_values(values)

    def _test_correlate_convolve(self, pyfunc):
        cfunc = jit(nopython=True)(pyfunc)
        # only 1d arrays are accepted, test varying lengths
        # and varying dtype
        lengths = (1, 2, 3, 7)
        dts = [np.int8, np.int32, np.int64, np.float32, np.float64,
               np.complex64, np.complex128]
        modes = ["full", "valid", "same"]

        for dt1, dt2, n, m, mode in itertools.product(
            dts, dts, lengths, lengths, modes
        ):
            a = np.arange(n, dtype=dt1)
            v = np.arange(m, dtype=dt2)

            if np.issubdtype(dt1, np.complexfloating):
                a = (a + 1j * a).astype(dt1)
            if np.issubdtype(dt2, np.complexfloating):
                v = (v + 1j * v).astype(dt2)

            expected = pyfunc(a, v, mode=mode)
            got = cfunc(a, v, mode=mode)

            self.assertPreciseEqual(expected, got)

        _a = np.arange(12).reshape(4, 3)
        _b = np.arange(12)
        for x, y in [(_a, _b), (_b, _a)]:
            with self.assertRaises(TypingError) as raises:
                cfunc(x, y)
            msg = 'only supported on 1D arrays'
            self.assertIn(msg, str(raises.exception))

    def test_correlate(self):
        self._test_correlate_convolve(correlate)

    def _test_correlate_convolve_exceptions(self, fn):
        # Exceptions leak references
        self.disable_leak_check()

        # convolve raises if either array has a 0 dimension
        _a = np.ones(shape=(0,))
        _b = np.arange(5)
        cfunc = jit(nopython=True)(fn)
        for x, y in [(_a, _b), (_b, _a)]:
            with self.assertRaises(ValueError) as raises:
                cfunc(x, y)
            if len(x) == 0:
                self.assertIn("'a' cannot be empty", str(raises.exception))
            else:
                self.assertIn("'v' cannot be empty", str(raises.exception))

        with self.assertRaises(ValueError) as raises:
            cfunc(_b, _b, mode="invalid mode")

            self.assertIn("Invalid 'mode'", str(raises.exception))

    def test_correlate_exceptions(self):
        # correlate supported 0 dimension arrays until 1.18
        self._test_correlate_convolve_exceptions(correlate)

    def test_convolve(self):
        self._test_correlate_convolve(convolve)

    def test_convolve_exceptions(self):
        self._test_correlate_convolve_exceptions(convolve)

    def _check_output(self, pyfunc, cfunc, params, abs_tol=None):
        expected = pyfunc(**params)
        got = cfunc(**params)
        self.assertPreciseEqual(expected, got, abs_tol=abs_tol)

    def test_vander_basic(self):
        pyfunc = vander
        cfunc = jit(nopython=True)(pyfunc)
        _check_output = partial(self._check_output, pyfunc, cfunc)

        def _check(x):
            n_choices = [None, 0, 1, 2, 3, 4]
            increasing_choices = [True, False]

            # N and increasing defaulted
            params = {'x': x}
            _check_output(params)

            # N provided and increasing defaulted
            for n in n_choices:
                params = {'x': x, 'N': n}
                _check_output(params)

            # increasing provided and N defaulted:
            for increasing in increasing_choices:
                params = {'x': x, 'increasing': increasing}
                _check_output(params)

            # both n and increasing supplied
            for n in n_choices:
                for increasing in increasing_choices:
                    params = {'x': x, 'N': n, 'increasing': increasing}
                    _check_output(params)

        _check(np.array([1, 2, 3, 5]))
        _check(np.arange(7) - 10.5)
        _check(np.linspace(3, 10, 5))
        _check(np.array([1.2, np.nan, np.inf, -np.inf]))
        _check(np.array([]))
        _check(np.arange(-5, 5) - 0.3)

        # # boolean array
        _check(np.array([True] * 5 + [False] * 4))

        # cycle through dtypes to check type promotion a la numpy
        for dtype in np.int32, np.int64, np.float32, np.float64:
            _check(np.arange(10, dtype=dtype))

        # non array inputs
        _check([0, 1, 2, 3])
        _check((4, 5, 6, 7))
        _check((0.0, 1.0, 2.0))
        _check(())

        # edge cases
        _check((3, 4.444, 3.142))
        _check((True, False, 4))

    def test_vander_exceptions(self):
        pyfunc = vander
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        x = np.arange(5) - 0.5

        def _check_n(N):
            with self.assertTypingError() as raises:
                cfunc(x, N=N)
            self.assertIn("Second argument N must be None or an integer",
                          str(raises.exception))

        for N in 1.1, True, np.inf, [1, 2]:
            _check_n(N)

        with self.assertRaises(ValueError) as raises:
            cfunc(x, N=-1)
        self.assertIn("Negative dimensions are not allowed",
                      str(raises.exception))

        def _check_1d(x):
            with self.assertRaises(ValueError) as raises:
                cfunc(x)
            self.assertEqual("x must be a one-dimensional array or sequence.",
                             str(raises.exception))

        x = np.arange(27).reshape((3, 3, 3))
        _check_1d(x)

        x = ((2, 3), (4, 5))
        _check_1d(x)

    def test_tri_n_basic(self):
        pyfunc = tri_n
        cfunc = jit(nopython=True)(pyfunc)
        _check = partial(self._check_output, pyfunc, cfunc)

        def n_variations():
            return np.arange(-4, 8)  # number of rows

        # N supplied, M and k defaulted
        for n in n_variations():
            params = {'N': n}
            _check(params)

    def test_tri_n_m_basic(self):
        pyfunc = tri_n_m
        cfunc = jit(nopython=True)(pyfunc)
        _check = partial(self._check_output, pyfunc, cfunc)

        def n_variations():
            return np.arange(-4, 8)  # number of rows

        def m_variations():
            # number of columns
            return itertools.chain.from_iterable(([None], range(-5, 9)))

        # N supplied, M and k defaulted
        for n in n_variations():
            params = {'N': n}
            _check(params)

        # N and M supplied, k defaulted
        for n in n_variations():
            for m in m_variations():
                params = {'N': n, 'M': m}
                _check(params)

    def test_tri_n_k_basic(self):
        pyfunc = tri_n_k
        cfunc = jit(nopython=True)(pyfunc)
        _check = partial(self._check_output, pyfunc, cfunc)

        def n_variations():
            return np.arange(-4, 8)  # number of rows

        def k_variations():
            return np.arange(-10, 10)  # offset

        # N supplied, M and k defaulted
        for n in n_variations():
            params = {'N': n}
            _check(params)

        # N and k supplied, M defaulted
        for n in n_variations():
            for k in k_variations():
                params = {'N': n, 'k': k}
                _check(params)

    def test_tri_n_m_k_basic(self):
        pyfunc = tri_n_m_k
        cfunc = jit(nopython=True)(pyfunc)
        _check = partial(self._check_output, pyfunc, cfunc)

        def n_variations():
            return np.arange(-4, 8)  # number of rows

        def m_variations():
            # number of columns
            return itertools.chain.from_iterable(([None], range(-5, 9)))

        def k_variations():
            return np.arange(-10, 10)  # offset

        # N supplied, M and k defaulted
        for n in n_variations():
            params = {'N': n}
            _check(params)

        # N and M supplied, k defaulted
        for n in n_variations():
            for m in m_variations():
                params = {'N': n, 'M': m}
                _check(params)

        # N and k supplied, M defaulted
        for n in n_variations():
            for k in k_variations():
                params = {'N': n, 'k': k}
                _check(params)

        # N, M and k supplied
        for n in n_variations():
            for k in k_variations():
                for m in m_variations():
                    params = {'N': n, 'M': m, 'k': k}
                    _check(params)

    def test_tri_exceptions(self):
        pyfunc = tri_n_m_k
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        def _check(k):
            with self.assertTypingError() as raises:
                cfunc(5, 6, k=k)
            assert "k must be an integer" in str(raises.exception)

        for k in 1.5, True, np.inf, [1, 2]:
            _check(k)

    def _triangular_matrix_tests_m(self, pyfunc):
        cfunc = jit(nopython=True)(pyfunc)

        def _check(arr):
            expected = pyfunc(arr)
            got = cfunc(arr)
            # TODO: Contiguity of result not consistent with numpy
            self.assertEqual(got.dtype, expected.dtype)
            np.testing.assert_array_equal(got, expected)

        return self._triangular_matrix_tests_inner(self, pyfunc, _check)

    def _triangular_matrix_tests_m_k(self, pyfunc):
        cfunc = jit(nopython=True)(pyfunc)

        def _check(arr):
            for k in itertools.chain.from_iterable(([None], range(-10, 10))):
                if k is None:
                    params = {}
                else:
                    params = {'k': k}
                expected = pyfunc(arr, **params)
                got = cfunc(arr, **params)
                # TODO: Contiguity of result not consistent with numpy
                self.assertEqual(got.dtype, expected.dtype)
                np.testing.assert_array_equal(got, expected)

        return self._triangular_matrix_tests_inner(self, pyfunc, _check)

    @staticmethod
    def _triangular_matrix_tests_inner(self, pyfunc, _check):

        def check_odd(a):
            _check(a)
            a = a.reshape((9, 7))
            _check(a)
            a = a.reshape((7, 1, 3, 3))
            _check(a)
            _check(a.T)

        def check_even(a):
            _check(a)
            a = a.reshape((4, 16))
            _check(a)
            a = a.reshape((4, 2, 2, 4))
            _check(a)
            _check(a.T)

        check_odd(np.arange(63) + 10.5)
        check_even(np.arange(64) - 10.5)

        # edge cases
        _check(np.arange(360).reshape(3, 4, 5, 6))
        _check(np.array([]))
        _check(np.arange(9).reshape((3, 3))[::-1])
        _check(np.arange(9).reshape((3, 3), order='F'))

        arr = (np.arange(64) - 10.5).reshape((4, 2, 2, 4))
        _check(arr)
        _check(np.asfortranarray(arr))

    def _triangular_matrix_exceptions(self, pyfunc):
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        a = np.ones((5, 6))
        with self.assertTypingError() as raises:
            cfunc(a, k=1.5)
            self.assertIn("k must be an integer", str(raises.exception))

    def _triangular_indices_tests_base(self, pyfunc, args):
        cfunc = jit(nopython=True)(pyfunc)

        for x in args:
            expected = pyfunc(*x)
            got = cfunc(*x)
            self.assertEqual(type(expected), type(got))
            self.assertEqual(len(expected), len(got))
            for e, g in zip(expected, got):
                np.testing.assert_array_equal(e, g)

    def _triangular_indices_tests_n(self, pyfunc):
        self._triangular_indices_tests_base(
            pyfunc,
            [[n] for n in range(10)]
        )

    def _triangular_indices_tests_n_k(self, pyfunc):
        self._triangular_indices_tests_base(
            pyfunc,
            [[n, k] for n in range(10) for k in range(-n - 1, n + 2)]
        )

    def _triangular_indices_tests_n_m(self, pyfunc):
        self._triangular_indices_tests_base(
            pyfunc,
            [[n, m] for n in range(10) for m in range(2 * n)]
        )

    def _triangular_indices_tests_n_k_m(self, pyfunc):
        self._triangular_indices_tests_base(
            pyfunc,
            [[n, k, m] for n in range(10)
             for k in range(-n - 1, n + 2)
             for m in range(2 * n)]
        )

        # Check jitted version works with default values for kwargs
        cfunc = jit(nopython=True)(pyfunc)
        cfunc(1)

    def _triangular_indices_from_tests_arr(self, pyfunc):
        cfunc = jit(nopython=True)(pyfunc)

        for dtype in [int, float, bool]:
            for n,m in itertools.product(range(10), range(10)):
                arr = np.ones((n, m), dtype)
                expected = pyfunc(arr)
                got = cfunc(arr)
                self.assertEqual(type(expected), type(got))
                self.assertEqual(len(expected), len(got))
                for e, g in zip(expected, got):
                    np.testing.assert_array_equal(e, g)

    def _triangular_indices_from_tests_arr_k(self, pyfunc):
        cfunc = jit(nopython=True)(pyfunc)

        for dtype in [int, float, bool]:
            for n,m in itertools.product(range(10), range(10)):
                arr = np.ones((n, m), dtype)
                for k in range(-10, 10):
                    expected = pyfunc(arr)
                    got = cfunc(arr)
                    self.assertEqual(type(expected), type(got))
                    self.assertEqual(len(expected), len(got))
                    for e, g in zip(expected, got):
                        np.testing.assert_array_equal(e, g)

    def _triangular_indices_exceptions(self, pyfunc):
        cfunc = jit(nopython=True)(pyfunc)
        parameters = pysignature(pyfunc).parameters

        with self.assertTypingError() as raises:
            cfunc(1.0)
        self.assertIn("n must be an integer", str(raises.exception))

        if 'k' in parameters:
            with self.assertTypingError() as raises:
                cfunc(1, k=1.0)
            self.assertIn("k must be an integer", str(raises.exception))

        if 'm' in parameters:
            with self.assertTypingError() as raises:
                cfunc(1, m=1.0)
            self.assertIn("m must be an integer", str(raises.exception))

    def _triangular_indices_from_exceptions(self, pyfunc, test_k=True):
        cfunc = jit(nopython=True)(pyfunc)

        for ndims in [0, 1, 3]:
            a = np.ones([5] * ndims)
            with self.assertTypingError() as raises:
                cfunc(a)
            self.assertIn("input array must be 2-d", str(raises.exception))

        if test_k:
            a = np.ones([5, 5])
            with self.assertTypingError() as raises:
                cfunc(a, k=0.5)
            self.assertIn("k must be an integer", str(raises.exception))

    def test_tril_basic(self):
        self._triangular_matrix_tests_m(tril_m)
        self._triangular_matrix_tests_m_k(tril_m_k)

    def test_tril_exceptions(self):
        self._triangular_matrix_exceptions(tril_m_k)

    def test_tril_indices(self):
        self._triangular_indices_tests_n(tril_indices_n)
        self._triangular_indices_tests_n_k(tril_indices_n_k)
        self._triangular_indices_tests_n_m(tril_indices_n_m)
        self._triangular_indices_tests_n_k_m(tril_indices_n_k_m)
        self._triangular_indices_exceptions(tril_indices_n)
        self._triangular_indices_exceptions(tril_indices_n_k)
        self._triangular_indices_exceptions(tril_indices_n_m)
        self._triangular_indices_exceptions(tril_indices_n_k_m)

    def test_tril_indices_from(self):
        self._triangular_indices_from_tests_arr(tril_indices_from_arr)
        self._triangular_indices_from_tests_arr_k(tril_indices_from_arr_k)
        self._triangular_indices_from_exceptions(tril_indices_from_arr, False)
        self._triangular_indices_from_exceptions(tril_indices_from_arr_k, True)

    def test_triu_basic(self):
        self._triangular_matrix_tests_m(triu_m)
        self._triangular_matrix_tests_m_k(triu_m_k)

    def test_triu_exceptions(self):
        self._triangular_matrix_exceptions(triu_m_k)

    def test_triu_indices(self):
        self._triangular_indices_tests_n(triu_indices_n)
        self._triangular_indices_tests_n_k(triu_indices_n_k)
        self._triangular_indices_tests_n_m(triu_indices_n_m)
        self._triangular_indices_tests_n_k_m(triu_indices_n_k_m)
        self._triangular_indices_exceptions(triu_indices_n)
        self._triangular_indices_exceptions(triu_indices_n_k)
        self._triangular_indices_exceptions(triu_indices_n_m)
        self._triangular_indices_exceptions(triu_indices_n_k_m)

    def test_triu_indices_from(self):
        self._triangular_indices_from_tests_arr(triu_indices_from_arr)
        self._triangular_indices_from_tests_arr_k(triu_indices_from_arr_k)
        self._triangular_indices_from_exceptions(triu_indices_from_arr, False)
        self._triangular_indices_from_exceptions(triu_indices_from_arr_k, True)

    def test_indices_basic(self):
        pyfunc = np_indices
        cfunc = njit(np_indices)

        def inputs():
            # Taken from https://github.com/numpy/numpy/blob/db4f43983cb938f12c311e1f5b7165e270c393b4/numpy/core/tests/test_numeric.py#L3383-L3407 # noqa: E501
            yield (4, 3)
            yield (4,)
            yield (0,)
            yield (2, 2, 3, 5)

        for dims in inputs():
            self.assertPreciseEqual(pyfunc(dims), cfunc(dims))

    def test_indices_exception(self):
        cfunc = njit(np_indices)

        self.disable_leak_check()

        errmsg = 'The argument "dimensions" must be a tuple of integers'
        with self.assertRaises(TypingError) as raises:
            cfunc("abc")
        self.assertIn(errmsg,
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc((2.0, 3.0))
        self.assertIn(errmsg,
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc((2, 3.0))
        self.assertIn(errmsg,
                      str(raises.exception))

    def partition_sanity_check(self, pyfunc, cfunc, a, kth):
        # as NumPy uses a different algorithm, we do not expect to
        # match outputs exactly...
        expected = pyfunc(a, kth)
        got = cfunc(a, kth)

        # but we do expect the unordered collection of elements up to the
        # kth to tie out
        self.assertPreciseEqual(np.unique(expected[:kth]), np.unique(got[:kth]))

        # likewise the unordered collection of elements from the kth onwards
        self.assertPreciseEqual(np.unique(expected[kth:]), np.unique(got[kth:]))

    def argpartition_sanity_check(self, pyfunc, cfunc, a, kth):
        # as NumPy uses a different algorithm, we do not expect to
        # match outputs exactly...
        expected = pyfunc(a, kth)
        got = cfunc(a, kth)

        # but we do expect the unordered collection of elements up to the
        # kth to tie out
        self.assertPreciseEqual(np.unique(a[expected[:kth]]),
                                np.unique(a[got[:kth]]))

        # likewise the unordered collection of elements from the kth onwards
        self.assertPreciseEqual(np.unique(a[expected[kth:]]),
                                np.unique(a[got[kth:]]))

    def test_partition_fuzz(self):
        # inspired by the test of the same name in:
        # https://github.com/numpy/numpy/blob/043a840/numpy/core/tests/test_multiarray.py    # noqa: E501
        pyfunc = partition
        cfunc = jit(nopython=True)(pyfunc)

        for j in range(10, 30):
            for i in range(1, j - 2):
                d = np.arange(j)
                self.rnd.shuffle(d)
                d = d % self.rnd.randint(2, 30)
                idx = self.rnd.randint(d.size)
                kth = [0, idx, i, i + 1, -idx, -i]  # include negative kth's
                tgt = np.sort(d)[kth]
                self.assertPreciseEqual(cfunc(d, kth)[kth],
                                        tgt)  # a -> array
                self.assertPreciseEqual(cfunc(d.tolist(), kth)[kth],
                                        tgt)  # a -> list
                self.assertPreciseEqual(cfunc(tuple(d.tolist()), kth)[kth],
                                        tgt)  # a -> tuple

                for k in kth:
                    self.partition_sanity_check(pyfunc, cfunc, d, k)

    def test_argpartition_fuzz(self):
        # inspired by the test of the same name in:
        # https://github.com/numpy/numpy/blob/043a840/numpy/core/tests/test_multiarray.py    # noqa: E501
        pyfunc = argpartition
        cfunc = jit(nopython=True)(pyfunc)

        for j in range(10, 30):
            for i in range(1, j - 2):
                d = np.arange(j)
                self.rnd.shuffle(d)
                d = d % self.rnd.randint(2, 30)
                idx = self.rnd.randint(d.size)
                kth = [0, idx, i, i + 1, -idx, -i]  # include negative kth's
                tgt = np.argsort(d)[kth]
                self.assertPreciseEqual(d[cfunc(d, kth)[kth]],
                                        d[tgt])  # a -> array
                self.assertPreciseEqual(d[cfunc(d.tolist(), kth)[kth]],
                                        d[tgt])  # a -> list
                self.assertPreciseEqual(d[cfunc(tuple(d.tolist()), kth)[kth]],
                                        d[tgt])  # a -> tuple

                for k in kth:
                    self.argpartition_sanity_check(pyfunc, cfunc, d, k)

    def test_partition_exception_out_of_range(self):
        # inspired by the test of the same name in:
        # https://github.com/numpy/numpy/blob/043a840/numpy/core/tests/test_multiarray.py    # noqa: E501
        pyfunc = partition
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        # Test out of range values in kth raise an error
        a = np.arange(10)

        def _check(a, kth):
            with self.assertRaises(ValueError) as e:
                cfunc(a, kth)
            assert str(e.exception) == "kth out of bounds"

        _check(a, 10)
        _check(a, -11)
        _check(a, (3, 30))

    def test_argpartition_exception_out_of_range(self):
        # inspired by the test of the same name in:
        # https://github.com/numpy/numpy/blob/043a840/numpy/core/tests/test_multiarray.py    # noqa: E501
        pyfunc = argpartition
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        # Test out of range values in kth raise an error
        a = np.arange(10)

        def _check(a, kth):
            with self.assertRaises(ValueError) as e:
                cfunc(a, kth)
            assert str(e.exception) == "kth out of bounds"

        _check(a, 10)
        _check(a, -11)
        _check(a, (3, 30))

    def test_partition_exception_non_integer_kth(self):
        # inspired by the test of the same name in:
        # https://github.com/numpy/numpy/blob/043a840/numpy/core/tests/test_multiarray.py    # noqa: E501
        pyfunc = partition
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        def _check(a, kth):
            with self.assertTypingError() as raises:
                cfunc(a, kth)
            self.assertIn("Partition index must be integer",
                          str(raises.exception))

        a = np.arange(10)
        _check(a, 9.0)
        _check(a, (3.3, 4.4))
        _check(a, np.array((1, 2, np.nan)))

    def test_argpartition_exception_non_integer_kth(self):
        # inspired by the test of the same name in:
        # https://github.com/numpy/numpy/blob/043a840/numpy/core/tests/test_multiarray.py    # noqa: E501
        pyfunc = argpartition
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        def _check(a, kth):
            with self.assertTypingError() as raises:
                cfunc(a, kth)
            self.assertIn("Partition index must be integer",
                          str(raises.exception))

        a = np.arange(10)
        _check(a, 9.0)
        _check(a, (3.3, 4.4))
        _check(a, np.array((1, 2, np.nan)))

    def test_partition_exception_a_not_array_like(self):
        pyfunc = partition
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        def _check(a, kth):
            with self.assertTypingError() as raises:
                cfunc(a, kth)
            self.assertIn('The first argument must be an array-like',
                          str(raises.exception))

        _check(4, 0)
        _check('Sausages', 0)

    def test_argpartition_exception_a_not_array_like(self):
        pyfunc = argpartition
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        def _check(a, kth):
            with self.assertTypingError() as raises:
                cfunc(a, kth)
            self.assertIn('The first argument must be an array-like',
                          str(raises.exception))

        _check(4, 0)
        _check('Sausages', 0)

    def test_partition_exception_a_zero_dim(self):
        pyfunc = partition
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        def _check(a, kth):
            with self.assertTypingError() as raises:
                cfunc(a, kth)
            self.assertIn('The first argument must be at least 1-D (found 0-D)',
                          str(raises.exception))

        _check(np.array(1), 0)

    def test_argpartition_exception_a_zero_dim(self):
        pyfunc = argpartition
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        def _check(a, kth):
            with self.assertTypingError() as raises:
                cfunc(a, kth)
            self.assertIn('The first argument must be at least 1-D (found 0-D)',
                          str(raises.exception))

        _check(np.array(1), 0)

    def test_partition_exception_kth_multi_dimensional(self):
        pyfunc = partition
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        def _check(a, kth):
            with self.assertRaises(ValueError) as raises:
                cfunc(a, kth)
            self.assertIn('kth must be scalar or 1-D', str(raises.exception))

        _check(np.arange(10), kth=np.arange(6).reshape(3, 2))

    def test_argpartition_exception_kth_multi_dimensional(self):
        pyfunc = argpartition
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        def _check(a, kth):
            with self.assertRaises(ValueError) as raises:
                cfunc(a, kth)
            self.assertIn('kth must be scalar or 1-D', str(raises.exception))

        _check(np.arange(10), kth=np.arange(6).reshape(3, 2))

    def test_partition_empty_array(self):
        # inspired by the test of the same name in:
        # https://github.com/numpy/numpy/blob/043a840/numpy/core/tests/test_multiarray.py    # noqa: E501
        pyfunc = partition
        cfunc = jit(nopython=True)(pyfunc)

        def check(a, kth=0):
            expected = pyfunc(a, kth)
            got = cfunc(a, kth)
            self.assertPreciseEqual(expected, got)

        # check axis handling for multidimensional empty arrays
        a = np.array([])
        a.shape = (3, 2, 1, 0)

        # include this with some other empty data structures
        for arr in a, (), np.array([]):
            check(arr)

    def test_argpartition_empty_array(self):
        # inspired by the test of the same name in:
        # https://github.com/numpy/numpy/blob/043a840/numpy/core/tests/test_multiarray.py    # noqa: E501
        pyfunc = argpartition
        cfunc = jit(nopython=True)(pyfunc)

        def check(a, kth=0):
            expected = pyfunc(a, kth)
            got = cfunc(a, kth)
            self.assertPreciseEqual(expected, got)

        # check axis handling for multidimensional empty arrays
        a = np.array([])
        a.shape = (3, 2, 1, 0)

        # include this with some other empty data structures
        for arr in a, (), np.array([]):
            check(arr)

    def test_partition_basic(self):
        # inspired by the test of the same name in:
        # https://github.com/numpy/numpy/blob/043a840/numpy/core/tests/test_multiarray.py    # noqa: E501
        pyfunc = partition
        cfunc = jit(nopython=True)(pyfunc)

        d = np.array([])
        got = cfunc(d, 0)
        self.assertPreciseEqual(d, got)

        d = np.ones(1)
        got = cfunc(d, 0)
        self.assertPreciseEqual(d, got)

        # kth not modified
        kth = np.array([30, 15, 5])
        okth = kth.copy()
        cfunc(np.arange(40), kth)
        self.assertPreciseEqual(kth, okth)

        for r in ([2, 1], [1, 2], [1, 1]):
            d = np.array(r)
            tgt = np.sort(d)
            for k in 0, 1:
                self.assertPreciseEqual(cfunc(d, k)[k], tgt[k])
                self.partition_sanity_check(pyfunc, cfunc, d, k)

        for r in ([3, 2, 1], [1, 2, 3], [2, 1, 3], [2, 3, 1],
                  [1, 1, 1], [1, 2, 2], [2, 2, 1], [1, 2, 1]):
            d = np.array(r)
            tgt = np.sort(d)
            for k in 0, 1, 2:
                self.assertPreciseEqual(cfunc(d, k)[k], tgt[k])
                self.partition_sanity_check(pyfunc, cfunc, d, k)

        d = np.ones(50)
        self.assertPreciseEqual(cfunc(d, 0), d)

        # sorted
        d = np.arange(49)
        for k in 5, 15:
            self.assertEqual(cfunc(d, k)[k], k)
            self.partition_sanity_check(pyfunc, cfunc, d, k)

        # rsorted, with input flavours: array, list and tuple
        d = np.arange(47)[::-1]
        for a in d, d.tolist(), tuple(d.tolist()):
            self.assertEqual(cfunc(a, 6)[6], 6)
            self.assertEqual(cfunc(a, 16)[16], 16)
            self.assertPreciseEqual(cfunc(a, -6), cfunc(a, 41))
            self.assertPreciseEqual(cfunc(a, -16), cfunc(a, 31))
            self.partition_sanity_check(pyfunc, cfunc, d, -16)

        # median of 3 killer, O(n^2) on pure median 3 pivot quickselect
        # exercises the median of median of 5 code used to keep O(n)
        d = np.arange(1000000)
        x = np.roll(d, d.size // 2)
        mid = x.size // 2 + 1
        self.assertEqual(cfunc(x, mid)[mid], mid)
        d = np.arange(1000001)
        x = np.roll(d, d.size // 2 + 1)
        mid = x.size // 2 + 1
        self.assertEqual(cfunc(x, mid)[mid], mid)

        # max
        d = np.ones(10)
        d[1] = 4
        self.assertEqual(cfunc(d, (2, -1))[-1], 4)
        self.assertEqual(cfunc(d, (2, -1))[2], 1)
        d[1] = np.nan
        assert np.isnan(cfunc(d, (2, -1))[-1])

        # equal elements
        d = np.arange(47) % 7
        tgt = np.sort(np.arange(47) % 7)
        self.rnd.shuffle(d)
        for i in range(d.size):
            self.assertEqual(cfunc(d, i)[i], tgt[i])
            self.partition_sanity_check(pyfunc, cfunc, d, i)

        d = np.array([0, 1, 2, 3, 4, 5, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
                      7, 7, 7, 7, 7, 9])
        kth = [0, 3, 19, 20]
        self.assertEqual(tuple(cfunc(d, kth)[kth]), (0, 3, 7, 7))

        td = [(dt, s) for dt in [np.int32, np.float32] for s in (9, 16)]
        for dt, s in td:
            d = np.arange(s, dtype=dt)
            self.rnd.shuffle(d)
            d1 = np.tile(np.arange(s, dtype=dt), (4, 1))
            map(self.rnd.shuffle, d1)
            for i in range(d.size):
                p = cfunc(d, i)
                self.assertEqual(p[i], i)
                # all before are smaller
                np.testing.assert_array_less(p[:i], p[i])
                # all after are larger
                np.testing.assert_array_less(p[i], p[i + 1:])
                # sanity check
                self.partition_sanity_check(pyfunc, cfunc, d, i)

    def test_argpartition_basic(self):
        # inspired by the test of the same name in:
        # https://github.com/numpy/numpy/blob/043a840/numpy/core/tests/test_multiarray.py    # noqa: E501
        pyfunc = argpartition
        cfunc = jit(nopython=True)(pyfunc)

        d = np.array([], dtype=np.int64)
        expected = pyfunc(d, 0)
        got = cfunc(d, 0)
        self.assertPreciseEqual(expected, got)

        d = np.ones(1, dtype=np.int64)
        expected = pyfunc(d, 0)
        got = cfunc(d, 0)
        self.assertPreciseEqual(expected, got)

        # kth not modified
        kth = np.array([30, 15, 5])
        okth = kth.copy()
        cfunc(np.arange(40), kth)
        self.assertPreciseEqual(kth, okth)

        for r in ([2, 1], [1, 2], [1, 1]):
            d = np.array(r)
            tgt = np.argsort(d)
            for k in 0, 1:
                self.assertPreciseEqual(d[cfunc(d, k)[k]], d[tgt[k]])
                self.argpartition_sanity_check(pyfunc, cfunc, d, k)

        for r in ([3, 2, 1], [1, 2, 3], [2, 1, 3], [2, 3, 1],
                  [1, 1, 1], [1, 2, 2], [2, 2, 1], [1, 2, 1]):
            d = np.array(r)
            tgt = np.argsort(d)
            for k in 0, 1, 2:
                self.assertPreciseEqual(d[cfunc(d, k)[k]], d[tgt[k]])
                self.argpartition_sanity_check(pyfunc, cfunc, d, k)

        d = np.ones(50)
        self.assertPreciseEqual(d[cfunc(d, 0)], d)

        # sorted
        d = np.arange(49)
        for k in 5, 15:
            self.assertEqual(cfunc(d, k)[k], k)
            self.partition_sanity_check(pyfunc, cfunc, d, k)

        # rsorted, with input flavours: array, list and tuple
        d = np.arange(47)[::-1]
        for a in d, d.tolist(), tuple(d.tolist()):
            self.assertEqual(cfunc(a, 6)[6], 40)
            self.assertEqual(cfunc(a, 16)[16], 30)
            self.assertPreciseEqual(cfunc(a, -6), cfunc(a, 41))
            self.assertPreciseEqual(cfunc(a, -16), cfunc(a, 31))
            self.argpartition_sanity_check(pyfunc, cfunc, d, -16)

        # median of 3 killer, O(n^2) on pure median 3 pivot quickselect
        # exercises the median of median of 5 code used to keep O(n)
        d = np.arange(1000000)
        x = np.roll(d, d.size // 2)
        mid = x.size // 2 + 1
        self.assertEqual(x[cfunc(x, mid)[mid]], mid)
        d = np.arange(1000001)
        x = np.roll(d, d.size // 2 + 1)
        mid = x.size // 2 + 1
        self.assertEqual(x[cfunc(x, mid)[mid]], mid)

        # max
        d = np.ones(10)
        d[1] = 4
        self.assertEqual(d[cfunc(d, (2, -1))[-1]], 4)
        self.assertEqual(d[cfunc(d, (2, -1))[2]], 1)
        d[1] = np.nan
        assert np.isnan(d[cfunc(d, (2, -1))[-1]])

        # equal elements
        d = np.arange(47) % 7
        tgt = np.sort(np.arange(47) % 7)
        self.rnd.shuffle(d)
        for i in range(d.size):
            self.assertEqual(d[cfunc(d, i)[i]], tgt[i])
            self.argpartition_sanity_check(pyfunc, cfunc, d, i)

        d = np.array([0, 1, 2, 3, 4, 5, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
                      7, 7, 7, 7, 7, 9])
        kth = [0, 3, 19, 20]
        self.assertEqual(tuple(d[cfunc(d, kth)[kth]]), (0, 3, 7, 7))

        td = [(dt, s) for dt in [np.int32, np.float32] for s in (9, 16)]
        for dt, s in td:
            d = np.arange(s, dtype=dt)
            self.rnd.shuffle(d)
            d1 = np.tile(np.arange(s, dtype=dt), (4, 1))
            map(self.rnd.shuffle, d1)
            for i in range(d.size):
                p = d[cfunc(d, i)]
                self.assertEqual(p[i], i)
                # all before are smaller
                np.testing.assert_array_less(p[:i], p[i])
                # all after are larger
                np.testing.assert_array_less(p[i], p[i + 1:])
                # sanity check
                self.argpartition_sanity_check(pyfunc, cfunc, d, i)

    def assert_partitioned(self, pyfunc, cfunc, d, kth):
        prev = 0
        for k in np.sort(kth):
            np.testing.assert_array_less(d[prev:k], d[k],
                                         err_msg='kth %d' % k)
            self.assertTrue((d[k:] >= d[k]).all(),
                            msg=("kth %d, %r not greater equal "
                                 "%d" % (k, d[k:], d[k])))
            prev = k + 1
            self.partition_sanity_check(pyfunc, cfunc, d, k)

    def assert_argpartitioned(self, pyfunc, cfunc, d, kth):
        prev = 0
        for k in np.sort(kth):
            np.testing.assert_array_less(d[prev:k], d[k],
                                         err_msg='kth %d' % k)
            self.assertTrue((d[k:] >= d[k]).all(),
                            msg=("kth %d, %r not greater equal "
                                 "%d" % (k, d[k:], d[k])))
            prev = k + 1
            self.argpartition_sanity_check(pyfunc, cfunc, d, k)

    def test_partition_iterative(self):
        # inspired by the test of the same name in:
        # https://github.com/numpy/numpy/blob/043a840/numpy/core/tests/test_multiarray.py
        pyfunc = partition
        cfunc = jit(nopython=True)(pyfunc)

        assert_partitioned = partial(self.assert_partitioned, pyfunc, cfunc)

        d = np.array([3, 4, 2, 1])
        p = cfunc(d, (0, 3))
        assert_partitioned(p, (0, 3))
        assert_partitioned(d[np.argpartition(d, (0, 3))], (0, 3))

        self.assertPreciseEqual(p, cfunc(d, (-3, -1)))

        d = np.arange(17)
        self.rnd.shuffle(d)
        self.assertPreciseEqual(np.arange(17), cfunc(d, list(range(d.size))))

        # test unsorted kth
        d = np.arange(17)
        self.rnd.shuffle(d)
        keys = np.array([1, 3, 8, -2])
        self.rnd.shuffle(d)
        p = cfunc(d, keys)
        assert_partitioned(p, keys)
        self.rnd.shuffle(keys)
        self.assertPreciseEqual(cfunc(d, keys), p)

        # equal kth
        d = np.arange(20)[::-1]
        assert_partitioned(cfunc(d, [5] * 4), [5])
        assert_partitioned(cfunc(d, [5] * 4 + [6, 13]), [5] * 4 + [6, 13])

    def test_argpartition_iterative(self):
        # inspired by the test of the same name in:
        # https://github.com/numpy/numpy/blob/043a840/numpy/core/tests/test_multiarray.py
        pyfunc = argpartition
        cfunc = jit(nopython=True)(pyfunc)

        assert_argpartitioned = partial(self.assert_argpartitioned,
                                        pyfunc,
                                        cfunc)

        d = np.array([3, 4, 2, 1])
        p = d[cfunc(d, (0, 3))]
        assert_argpartitioned(p, (0, 3))
        assert_argpartitioned(d[np.argpartition(d, (0, 3))], (0, 3))

        self.assertPreciseEqual(p, d[cfunc(d, (-3, -1))])

        d = np.arange(17)
        self.rnd.shuffle(d)
        self.assertPreciseEqual(np.arange(17), d[cfunc(d, list(range(d.size)))])

        # test unsorted kth
        d = np.arange(17)
        self.rnd.shuffle(d)
        keys = np.array([1, 3, 8, -2])
        self.rnd.shuffle(d)
        p = d[cfunc(d, keys)]
        assert_argpartitioned(p, keys)
        self.rnd.shuffle(keys)
        self.assertPreciseEqual(d[cfunc(d, keys)], p)

        # equal kth
        d = np.arange(20)[::-1]
        assert_argpartitioned(d[cfunc(d, [5] * 4)], [5])
        assert_argpartitioned(d[cfunc(d, [5] * 4 + [6, 13])], [5] * 4 + [6, 13])

    def test_partition_multi_dim(self):
        pyfunc = partition
        cfunc = jit(nopython=True)(pyfunc)

        def check(a, kth):
            expected = pyfunc(a, kth)
            got = cfunc(a, kth)
            self.assertPreciseEqual(expected[:, :, kth], got[:, :, kth])

            for s in np.ndindex(expected.shape[:-1]):
                self.assertPreciseEqual(np.unique(expected[s][:kth]),
                                        np.unique(got[s][:kth]))
                self.assertPreciseEqual(np.unique(expected[s][kth:]),
                                        np.unique(got[s][kth:]))

        def a_variations(a):
            yield a
            yield a.T
            yield np.asfortranarray(a)
            yield np.full_like(a, fill_value=np.nan)
            yield np.full_like(a, fill_value=np.inf)
            # multi-dimensional tuple input
            yield (((1.0, 3.142, -np.inf, 3),),)

        a = np.linspace(1, 10, 48)
        a[4:7] = np.nan
        a[8] = -np.inf
        a[9] = np.inf
        a = a.reshape((4, 3, 4))

        for arr in a_variations(a):
            for k in range(-3, 3):
                check(arr, k)

    def test_argpartition_multi_dim(self):
        pyfunc = argpartition
        cfunc = jit(nopython=True)(pyfunc)

        def check(a, kth):
            expected = pyfunc(a, kth)
            got = cfunc(a, kth)
            a = np.asarray(a)
            idx = np.ndindex(a.shape[:-1])
            for s in idx:
                self.assertPreciseEqual(a[s][expected[s][kth]],
                                        a[s][got[s][kth]])

            for s in np.ndindex(expected.shape[:-1]):
                self.assertPreciseEqual(np.unique(a[s][expected[s][:kth]]),
                                        np.unique(a[s][got[s][:kth]]))
                self.assertPreciseEqual(np.unique(a[s][expected[s][kth:]]),
                                        np.unique(a[s][got[s][kth:]]))

        def a_variations(a):
            yield a
            yield a.T
            yield np.asfortranarray(a)
            yield np.full_like(a, fill_value=np.nan)
            yield np.full_like(a, fill_value=np.inf)
            # multi-dimensional tuple input
            yield (((1.0, 3.142, -np.inf, 3),),)

        a = np.linspace(1, 10, 48)
        a[4:7] = np.nan
        a[8] = -np.inf
        a[9] = np.inf
        a = a.reshape((4, 3, 4))

        for arr in a_variations(a):
            for k in range(-3, 3):
                check(arr, k)

    def test_partition_boolean_inputs(self):
        pyfunc = partition
        cfunc = jit(nopython=True)(pyfunc)

        for d in np.linspace(1, 10, 17), np.array((True, False, True)):
            for kth in True, False, -1, 0, 1:
                self.partition_sanity_check(pyfunc, cfunc, d, kth)

    def test_argpartition_boolean_inputs(self):
        pyfunc = argpartition
        cfunc = jit(nopython=True)(pyfunc)

        for d in np.linspace(1, 10, 17), np.array((True, False, True)):
            for kth in True, False, -1, 0, 1:
                self.argpartition_sanity_check(pyfunc, cfunc, d, kth)

    @needs_blas
    def test_cov_invalid_ddof(self):
        pyfunc = cov
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        m = np.array([[0, 2], [1, 1], [2, 0]]).T

        for ddof in np.arange(4), 4j:
            with self.assertTypingError() as raises:
                cfunc(m, ddof=ddof)
            self.assertIn('ddof must be a real numerical scalar type',
                          str(raises.exception))

        for ddof in np.nan, np.inf:
            with self.assertRaises(ValueError) as raises:
                cfunc(m, ddof=ddof)
            self.assertIn('Cannot convert non-finite ddof to integer',
                          str(raises.exception))

        for ddof in 1.1, -0.7:
            with self.assertRaises(ValueError) as raises:
                cfunc(m, ddof=ddof)
            self.assertIn('ddof must be integral value', str(raises.exception))

    def corr_corrcoef_basic(self, pyfunc, first_arg_name):
        cfunc = jit(nopython=True)(pyfunc)
        _check = partial(self._check_output, pyfunc, cfunc, abs_tol=1e-14)

        def input_variations():
            # array inputs
            yield np.array([[0, 2], [1, 1], [2, 0]]).T
            yield self.rnd.randn(100).reshape(5, 20)
            yield np.asfortranarray(np.array([[0, 2], [1, 1], [2, 0]]).T)
            yield self.rnd.randn(100).reshape(5, 20)[:, ::2]
            yield np.array([0.3942, 0.5969, 0.7730, 0.9918, 0.7964])
            yield np.full((4, 5), fill_value=True)
            yield np.array([np.nan, 0.5969, -np.inf, 0.9918, 0.7964])
            yield np.linspace(-3, 3, 33).reshape(33, 1)

            # non-array inputs
            yield ((0.1, 0.2), (0.11, 0.19), (0.09, 0.21))  # UniTuple
            yield ((0.1, 0.2), (0.11, 0.19), (0.09j, 0.21j))  # Tuple
            yield (-2.1, -1, 4.3)
            yield (1, 2, 3)
            yield [4, 5, 6]
            yield ((0.1, 0.2, 0.3), (0.1, 0.2, 0.3))
            yield [(1, 2, 3), (1, 3, 2)]
            yield 3.142
            yield ((1.1, 2.2, 1.5),)

            # empty data structures
            yield np.array([])
            yield np.array([]).reshape(0, 2)
            yield np.array([]).reshape(2, 0)
            yield ()

        # all inputs other than the first are defaulted
        for input_arr in input_variations():
            _check({first_arg_name: input_arr})

    @needs_blas
    def test_corrcoef_basic(self):
        pyfunc = corrcoef
        self.corr_corrcoef_basic(pyfunc, first_arg_name='x')

    @needs_blas
    def test_cov_basic(self):
        pyfunc = cov
        self.corr_corrcoef_basic(pyfunc, first_arg_name='m')

    @needs_blas
    def test_cov_explicit_arguments(self):
        pyfunc = cov
        cfunc = jit(nopython=True)(pyfunc)
        _check = partial(self._check_output, pyfunc, cfunc, abs_tol=1e-14)

        m = self.rnd.randn(105).reshape(15, 7)
        y_choices = None, m[::-1]
        rowvar_choices = False, True
        bias_choices = False, True
        ddof_choice = None, -1, 0, 1, 3.0, True

        products = itertools.product(y_choices, rowvar_choices,
                                     bias_choices, ddof_choice)
        for y, rowvar, bias, ddof in products:
            params = {'m': m, 'y': y, 'ddof': ddof,
                      'bias': bias, 'rowvar': rowvar}
            _check(params)

    @needs_blas
    def test_corrcoef_explicit_arguments(self):
        pyfunc = corrcoef
        cfunc = jit(nopython=True)(pyfunc)
        _check = partial(self._check_output, pyfunc, cfunc, abs_tol=1e-14)

        x = self.rnd.randn(105).reshape(15, 7)
        y_choices = None, x[::-1]
        rowvar_choices = False, True

        for y, rowvar in itertools.product(y_choices, rowvar_choices):
            params = {'x': x, 'y': y, 'rowvar': rowvar}
            _check(params)

    def cov_corrcoef_edge_cases(self, pyfunc, first_arg_name):
        cfunc = jit(nopython=True)(pyfunc)
        _check = partial(self._check_output, pyfunc, cfunc, abs_tol=1e-14)

        # some of these examples borrowed from numpy doc string examples:
        # https://github.com/numpy/numpy/blob/v1.15.0/numpy/lib/function_base.py#L2199-L2231    # noqa: E501
        # some borrowed from TestCov and TestCorrCoef:
        # https://github.com/numpy/numpy/blob/80d3a7a/numpy/lib/tests/test_function_base.py    # noqa: E501
        m = np.array([-2.1, -1, 4.3])
        y = np.array([3, 1.1, 0.12])
        params = {first_arg_name: m, 'y': y}
        _check(params)

        m = np.array([1, 2, 3])  # test case modified such that m is 1D
        y = np.array([[1j, 2j, 3j]])
        params = {first_arg_name: m, 'y': y}
        _check(params)

        m = np.array([1, 2, 3])
        y = (1j, 2j, 3j)
        params = {first_arg_name: m, 'y': y}
        _check(params)
        params = {first_arg_name: y, 'y': m}  # flip real and complex inputs
        _check(params)

        m = np.array([1, 2, 3])
        y = (1j, 2j, 3)  # note last item is not complex
        params = {first_arg_name: m, 'y': y}
        _check(params)
        params = {first_arg_name: y, 'y': m}  # flip real and complex inputs
        _check(params)

        m = np.array([])
        y = np.array([])
        params = {first_arg_name: m, 'y': y}
        _check(params)

        m = 1.1
        y = 2.2
        params = {first_arg_name: m, 'y': y}
        _check(params)

        m = self.rnd.randn(10, 3)
        y = np.array([-2.1, -1, 4.3]).reshape(1, 3) / 10
        params = {first_arg_name: m, 'y': y}
        _check(params)

        m = np.array([-2.1, -1, 4.3])
        y = np.array([[3, 1.1, 0.12], [3, 1.1, 0.12]])
        params = {first_arg_name: m, 'y': y}
        _check(params)

        for rowvar in False, True:
            m = np.array([-2.1, -1, 4.3])
            y = np.array([[3, 1.1, 0.12], [3, 1.1, 0.12], [4, 1.1, 0.12]])
            params = {first_arg_name: m, 'y': y, 'rowvar': rowvar}
            _check(params)
            # swap m and y
            params = {first_arg_name: y, 'y': m, 'rowvar': rowvar}
            _check(params)

    @needs_blas
    def test_corrcoef_edge_cases(self):
        pyfunc = corrcoef
        self.cov_corrcoef_edge_cases(pyfunc, first_arg_name='x')

        cfunc = jit(nopython=True)(pyfunc)
        _check = partial(self._check_output, pyfunc, cfunc, abs_tol=1e-14)

        for x in (np.nan, -np.inf, 3.142, 0):
            params = {'x': x}
            _check(params)

    @needs_blas
    def test_corrcoef_edge_case_extreme_values(self):
        pyfunc = corrcoef
        cfunc = jit(nopython=True)(pyfunc)
        _check = partial(self._check_output, pyfunc, cfunc, abs_tol=1e-14)

        # extreme values
        x = ((1e-100, 1e100), (1e100, 1e-100))
        params = {'x': x}
        _check(params)

    @needs_blas
    def test_cov_edge_cases(self):
        pyfunc = cov
        self.cov_corrcoef_edge_cases(pyfunc, first_arg_name='m')

        cfunc = jit(nopython=True)(pyfunc)
        _check = partial(self._check_output, pyfunc, cfunc, abs_tol=1e-14)

        # invalid ddof
        m = np.array([[0, 2], [1, 1], [2, 0]]).T
        params = {'m': m, 'ddof': 5}
        _check(params)

    @needs_blas
    def test_cov_exceptions(self):
        pyfunc = cov
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        def _check_m(m):
            with self.assertTypingError() as raises:
                cfunc(m)
            self.assertIn('m has more than 2 dimensions', str(raises.exception))

        m = np.ones((5, 6, 7))
        _check_m(m)

        m = ((((1, 2, 3), (2, 2, 2)),),)
        _check_m(m)

        m = [[[5, 6, 7]]]
        _check_m(m)

        def _check_y(m, y):
            with self.assertTypingError() as raises:
                cfunc(m, y=y)
            self.assertIn('y has more than 2 dimensions', str(raises.exception))

        m = np.ones((5, 6))
        y = np.ones((5, 6, 7))
        _check_y(m, y)

        m = np.array((1.1, 2.2, 1.1))
        y = (((1.2, 2.2, 2.3),),)
        _check_y(m, y)

        m = np.arange(3)
        y = np.arange(4)
        with self.assertRaises(ValueError) as raises:
            cfunc(m, y=y)
        self.assertIn('m and y have incompatible dimensions',
                      str(raises.exception))
        # Numpy raises ValueError: all the input array dimensions except for the
        # concatenation axis must match exactly.

        m = np.array([-2.1, -1, 4.3]).reshape(1, 3)
        with self.assertRaises(RuntimeError) as raises:
            cfunc(m)
        self.assertIn('2D array containing a single row is unsupported',
                      str(raises.exception))

    def test_ediff1d_basic(self):
        pyfunc = ediff1d
        cfunc = jit(nopython=True)(pyfunc)
        _check = partial(self._check_output, pyfunc, cfunc)

        def to_variations(a):
            yield None
            yield a
            yield a.astype(np.int16)

        def ary_variations(a):
            yield a
            yield a.reshape(3, 2, 2)
            yield a.astype(np.int32)

        for ary in ary_variations(np.linspace(-2, 7, 12)):
            params = {'ary': ary}
            _check(params)

            for a in to_variations(ary):
                params = {'ary': ary, 'to_begin': a}
                _check(params)

                params = {'ary': ary, 'to_end': a}
                _check(params)

                for b in to_variations(ary):
                    params = {'ary': ary, 'to_begin': a, 'to_end': b}
                    _check(params)

    def test_ediff1d_exceptions(self):
        pyfunc = ediff1d
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        with self.assertTypingError() as e:
            cfunc(np.array((True, True, False)))

        msg = "Boolean dtype is unsupported (as per NumPy)"
        assert msg in str(e.exception)

    def test_fliplr_basic(self):
        pyfunc = fliplr
        cfunc = jit(nopython=True)(pyfunc)

        def a_variations():
            yield np.arange(10).reshape(5, 2)
            yield np.arange(20).reshape(5, 2, 2)
            yield ((1, 2),)
            yield ([1, 2], [3, 4],)

        for a in a_variations():
            expected = pyfunc(a)
            got = cfunc(a)
            self.assertPreciseEqual(expected, got)

        with self.assertRaises(TypingError) as raises:
            cfunc("abc")

        self.assertIn("Cannot np.fliplr on %s type" % types.unicode_type,
                      str(raises.exception))

    def test_fliplr_exception(self):
        pyfunc = fliplr
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        with self.assertRaises(TypingError) as raises:
            cfunc(np.arange(3))

        self.assertIn("cannot index array", str(raises.exception))
        self.assertIn("with 2 indices", str(raises.exception))

    def test_flipud_basic(self):
        pyfunc = flipud
        cfunc = jit(nopython=True)(pyfunc)

        def a_variations():
            yield [1]
            yield np.arange(10)
            yield np.arange(10).reshape(5, 2)
            yield np.arange(20).reshape(5, 2, 2)
            yield ((1, 2),)
            yield ([1, 2], [3, 4],)

        for a in a_variations():
            expected = pyfunc(a)
            got = cfunc(a)
            self.assertPreciseEqual(expected, got)

        with self.assertRaises(TypingError) as raises:
            cfunc("abc")

        self.assertIn("Cannot np.flipud on %s type" % types.unicode_type,
                      str(raises.exception))

    def test_flipud_exception(self):
        pyfunc = flipud
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        with self.assertRaises(TypingError) as raises:
            cfunc(1)

        self.assertIn("cannot index array", str(raises.exception))
        self.assertIn("with 1 indices", str(raises.exception))

    def test_flip_basic(self):
        pyfunc = flip
        cfunc = jit(nopython=True)(pyfunc)

        def a_variations():
            yield np.array(1)
            yield np.arange(10)
            yield np.arange(10).reshape(5, 2)
            yield np.arange(20).reshape(5, 2, 2)

        for a in a_variations():
            expected = pyfunc(a)
            got = cfunc(a)
            self.assertPreciseEqual(expected, got)

        with self.assertRaises(TypingError) as raises:
            cfunc((1, 2, 3))

        self.assertIn("Cannot np.flip on UniTuple", str(raises.exception))

    def test_logspace2_basic(self):

        def inputs():
            #start, stop
            yield 1, 60
            yield -1, 60
            yield -60, -1
            yield -1, -60
            yield 60, -1
            yield 1.0, 60.0
            yield -60.0, -1.0
            yield -1.0, 60.0
            yield 0.0, np.e
            yield 0.0, np.pi
            if numpy_version < (2, 0):
                yield np.complex64(1), np.complex64(2)
                yield np.complex64(2j), np.complex64(4j)
                yield np.complex64(2), np.complex64(4j)
                yield np.complex64(1 + 2j), np.complex64(3 + 4j)
                yield np.complex64(1 - 2j), np.complex64(3 - 4j)
                yield np.complex64(-1 + 2j), np.complex64(3 + 4j)

        pyfunc = logspace2
        cfunc = jit(nopython=True)(pyfunc)

        for start, stop in inputs():
            np.testing.assert_allclose(pyfunc(start, stop), cfunc(start, stop))

    def test_logspace2_exception(self):
        cfunc = jit(nopython=True)(logspace2)

        self.disable_leak_check()

        with self.assertRaises(TypingError) as raises:
            cfunc("abc", 5)
        self.assertIn('The first argument "start" must be a number',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc(5, "abc")
        self.assertIn('The second argument "stop" must be a number',
                      str(raises.exception))

    def test_logspace3_basic(self):

        def inputs():
            #start, stop
            yield 1, 60
            yield -1, 60
            yield -60, -1
            yield -1, -60
            yield 60, -1
            yield 1.0, 60.0
            yield -60.0, -1.0
            yield -1.0, 60.0
            yield 0.0, np.e
            yield 0.0, np.pi
            if numpy_version < (2, 0):
                yield np.complex64(1), np.complex64(2)
                yield np.complex64(2j), np.complex64(4j)
                yield np.complex64(2), np.complex64(4j)
                yield np.complex64(1 + 2j), np.complex64(3 + 4j)
                yield np.complex64(1 - 2j), np.complex64(3 - 4j)
                yield np.complex64(-1 + 2j), np.complex64(3 + 4j)

        pyfunc = logspace3
        cfunc = jit(nopython=True)(pyfunc)

        for start, stop in inputs():
            np.testing.assert_allclose(pyfunc(start, stop), cfunc(start, stop))

    def test_logspace3_with_num_basic(self):

        def inputs():
            #start, stop, num
            yield 1, 60, 20
            yield -1, 60, 30
            yield -60, -1, 40
            yield -1, -60, 50
            yield 60, -1, 60
            yield 1.0, 60.0, 70
            yield -60.0, -1.0, 80
            yield -1.0, 60.0, 90
            yield 0.0, np.e, 20
            yield 0.0, np.pi, 30
            if numpy_version < (2, 0):
                yield np.complex64(1), np.complex64(2), 40
                yield np.complex64(2j), np.complex64(4j), 50
                yield np.complex64(2), np.complex64(4j), 60
                yield np.complex64(1 + 2j), np.complex64(3 + 4j), 70
                yield np.complex64(1 - 2j), np.complex64(3 - 4j), 80
                yield np.complex64(-1 + 2j), np.complex64(3 + 4j), 90

        pyfunc = logspace3
        cfunc = jit(nopython=True)(pyfunc)

        for start, stop, num in inputs():
            np.testing.assert_allclose(pyfunc(start, stop, num),
                                       cfunc(start, stop, num))

    def test_logspace3_exception(self):
        cfunc = jit(nopython=True)(logspace3)

        self.disable_leak_check()

        with self.assertRaises(TypingError) as raises:
            cfunc("abc", 5)
        self.assertIn('The first argument "start" must be a number',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc(5, "abc")
        self.assertIn('The second argument "stop" must be a number',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc(0, 5, "abc")
        self.assertIn('The third argument "num" must be an integer',
                      str(raises.exception))

    def test_geomspace2_basic(self):

        def inputs():
            #start, stop
            yield -1, -60
            yield 1.0, 60.0
            yield -60.0, -1.0
            yield 1, 1000
            yield 1000, 1
            yield 1, 256
            yield -1000, -1
            yield -1, np.complex64(2j)
            yield np.complex64(2j), -1
            yield -1.0, np.complex64(2j)
            yield np.complex64(1j), np.complex64(1000j)
            yield np.complex64(-1 + 0j), np.complex64(1 + 0j)
            yield np.complex64(1), np.complex64(2)
            yield np.complex64(2j), np.complex64(4j)
            yield np.complex64(2), np.complex64(4j)
            yield np.complex64(1 + 2j), np.complex64(3 + 4j)
            yield np.complex64(1 - 2j), np.complex64(3 - 4j)
            yield np.complex64(-1 + 2j), np.complex64(3 + 4j)

        pyfunc = geomspace2
        cfunc = jit(nopython=True)(pyfunc)

        for start, stop in inputs():
            self.assertPreciseEqual(pyfunc(start, stop),
                                    cfunc(start, stop),
                                    abs_tol=1e-12)

    def test_geomspace2_exception(self):
        cfunc = jit(nopython=True)(geomspace2)

        self.disable_leak_check()

        with self.assertRaises(TypingError) as raises:
            cfunc("abc", 5)
        self.assertIn('The argument "start" must be a number',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc(5, "abc")
        self.assertIn('The argument "stop" must be a number',
                      str(raises.exception))

        with self.assertRaises(ValueError) as raises:
            cfunc(0, 5)
        self.assertIn('Geometric sequence cannot include zero',
                      str(raises.exception))

        with self.assertRaises(ValueError) as raises:
            cfunc(5, 0)
        self.assertIn('Geometric sequence cannot include zero',
                      str(raises.exception))

    def test_geomspace3_basic(self):

        def inputs():
            #start, stop, num
            yield -1, -60, 50
            yield 1.0, 60.0, 70
            yield -60.0, -1.0, 80
            yield 1, 1000, 4
            yield 1, 1000, 3
            yield 1000, 1, 4
            yield 1, 256, 9
            yield -1000, -1, 4
            yield -1, np.complex64(2j), 10
            yield np.complex64(2j), -1, 20
            yield -1.0, np.complex64(2j), 30
            yield np.complex64(1j), np.complex64(1000j), 4
            yield np.complex64(-1 + 0j), np.complex64(1 + 0j), 5
            yield np.complex64(1), np.complex64(2), 40
            yield np.complex64(2j), np.complex64(4j), 50
            yield np.complex64(2), np.complex64(4j), 60
            yield np.complex64(1 + 2j), np.complex64(3 + 4j), 70
            yield np.complex64(1 - 2j), np.complex64(3 - 4j), 80
            yield np.complex64(-1 + 2j), np.complex64(3 + 4j), 90

        pyfunc = geomspace3
        cfunc = jit(nopython=True)(pyfunc)

        for start, stop, num in inputs():
            self.assertPreciseEqual(pyfunc(start, stop, num),
                                    cfunc(start, stop, num),
                                    abs_tol=1e-14)

    def test_geomspace3_exception(self):
        cfunc = jit(nopython=True)(geomspace3)

        self.disable_leak_check()

        with self.assertRaises(TypingError) as raises:
            cfunc("abc", 5, 10)
        self.assertIn('The argument "start" must be a number',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc(5, "abc", 10)
        self.assertIn('The argument "stop" must be a number',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc(5, 10, "abc")
        self.assertIn('The argument "num" must be an integer',
                      str(raises.exception))

        with self.assertRaises(ValueError) as raises:
            cfunc(0, 5, 5)
        self.assertIn('Geometric sequence cannot include zero',
                      str(raises.exception))

        with self.assertRaises(ValueError) as raises:
            cfunc(5, 0, 5)
        self.assertIn('Geometric sequence cannot include zero',
                      str(raises.exception))

    def test_geomspace_numpy(self):
        cfunc2 = jit(nopython=True)(geomspace2)
        cfunc3 = jit(nopython=True)(geomspace3)
        pfunc3 = geomspace3

        # https://github.com/numpy/numpy/blob/ab2178b47c0ee834180c318db196976623710691/numpy/core/tests/test_function_base.py#L122C3-L142
        # test_basic
        y = cfunc2(1, 1e6)
        self.assertEqual(len(y), 50)
        y = cfunc3(1, 1e6, num=100)
        self.assertEqual(y[-1], 10 ** 6)
        y = cfunc3(1, 1e6, num=7)
        self.assertPreciseEqual(y, pfunc3(1,1e6, num=7))

        y = cfunc3(8, 2, num=3)
        self.assertPreciseEqual(y, pfunc3(8, 2, num=3))
        self.assertTrue([x == 0 for x in y.imag])

        y = cfunc3(-1, -100, num=3)
        self.assertPreciseEqual(y, pfunc3(-1, -100, num=3))
        self.assertTrue([x == 0 for x in y.imag])

        y = cfunc3(-100, -1, num=3)
        self.assertPreciseEqual(y, pfunc3(-100, -1, num=3))
        self.assertTrue([x == 0 for x in y.imag])

        # test_boundaries_match_start_and_stop_exactly
        start = 0.3
        stop = 20.3

        y = cfunc3(start, stop, num=1)
        self.assertPreciseEqual(y[0], start)

        y = cfunc3(start, stop, num=3)
        self.assertPreciseEqual(y[0], start)
        self.assertPreciseEqual(y[-1], stop)

        # test_nan_interior
        with np.errstate(invalid='ignore'):
            y = cfunc3(-3, 3, num=4)

        self.assertPreciseEqual(y[0], -3.0)
        self.assertTrue(np.isnan(y[1:-1]).all())
        self.assertPreciseEqual(y[3], 3.0)

        # test_complex
        # Purely imaginary
        y = cfunc3(1j, 16j, num=5)
        self.assertPreciseEqual(y, pfunc3(1j, 16j, num=5), abs_tol=1e-14)
        self.assertTrue([x == 0 for x in y.real])

        y = cfunc3(-4j, -324j, num=5)
        self.assertPreciseEqual(y, pfunc3(-4j, -324j, num=5), abs_tol=1e-13)
        self.assertTrue([x == 0 for x in y.real])

        y = cfunc3(1 + 1j, 1000 + 1000j, num=4)
        self.assertPreciseEqual(y,
                                pfunc3(1 + 1j, 1000 + 1000j, num=4),
                                abs_tol=1e-13)

        y = cfunc3(-1 + 1j, -1000 + 1000j, num=4)
        self.assertPreciseEqual(y,
                                pfunc3(-1 + 1j, -1000 + 1000j, num=4),
                                abs_tol=1e-13)

        # Logarithmic spirals
        if numpy_version < (2, 0):
            y = cfunc3(-1 + 0j, 1 + 0j, num=3)
            self.assertPreciseEqual(y, pfunc3(-1 + 0j, 1 + 0j, num=3))

        y = cfunc3(0 + 3j, -3 + 0j, 3)
        self.assertPreciseEqual(y, pfunc3(0 + 3j, -3 + 0j, 3), abs_tol=1e-15)
        y = cfunc3(0 + 3j, 3 + 0j, 3)
        self.assertPreciseEqual(y, pfunc3(0 + 3j, 3 + 0j, 3), abs_tol=1e-15)
        y = cfunc3(-3 + 0j, 0 - 3j, 3)
        self.assertPreciseEqual(y, pfunc3(-3 + 0j, 0 - 3j, 3), abs_tol=1e-15)
        y = cfunc3(0 + 3j, -3 + 0j, 3)
        self.assertPreciseEqual(y, pfunc3(0 + 3j, -3 + 0j, 3), abs_tol=1e-15)
        y = cfunc3(-2 - 3j, 5 + 7j, 7)
        self.assertPreciseEqual(y, pfunc3(-2 - 3j, 5 + 7j, 7), abs_tol=1e-14)

        y = cfunc3(3j, -5, 2)
        self.assertPreciseEqual(y, pfunc3(3j, -5, 2))
        y = cfunc3(-5, 3j, 2)
        self.assertPreciseEqual(y, pfunc3(-5, 3j, 2))

    def test_rot90_basic(self):
        pyfunc = rot90
        cfunc = jit(nopython=True)(pyfunc)

        def a_variations():
            yield np.arange(10).reshape(5, 2)
            yield np.arange(20).reshape(5, 2, 2)
            yield np.arange(64).reshape(2, 2, 2, 2, 2, 2)

        for a in a_variations():
            expected = pyfunc(a)
            got = cfunc(a)
            self.assertPreciseEqual(expected, got)

    def test_rot90_with_k_basic(self):
        pyfunc = rot90_k
        cfunc = jit(nopython=True)(pyfunc)

        def a_variations():
            yield np.arange(10).reshape(5, 2)
            yield np.arange(20).reshape(5, 2, 2)
            yield np.arange(64).reshape(2, 2, 2, 2, 2, 2)

        for a in a_variations():
            for k in range(-5, 6):
                expected = pyfunc(a, k)
                got = cfunc(a, k)
                self.assertPreciseEqual(expected, got)

    def test_rot90_exception(self):
        pyfunc = rot90_k
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        with self.assertRaises(TypingError) as raises:
            cfunc("abc")

        self.assertIn('The first argument "m" must be an array',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc(np.arange(4).reshape(2, 2), k="abc")

        self.assertIn('The second argument "k" must be an integer',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc(np.arange(3))

        self.assertIn("Input must be >= 2-d.", str(raises.exception))

    def _check_split(self, func):
        # Since np.split and np.array_split are very similar
        pyfunc = func
        cfunc = jit(nopython=True)(pyfunc)

        def args_variations():
            a = np.arange(100)
            yield a, 2
            yield a, 2, 0
            yield a, [1, 4, 72]
            yield list(a), [1, 4, 72]
            yield tuple(a), [1, 4, 72]
            yield a, [1, 4, 72], 0
            yield list(a), [1, 4, 72], 0
            yield tuple(a), [1, 4, 72], 0

            a = np.arange(64).reshape(4, 4, 4)
            yield a, 2
            yield a, 2, 0
            yield a, 2, 1
            yield a, [2, 1, 5]
            yield a, [2, 1, 5], 1
            yield a, [2, 1, 5], 2
            yield a, [1, 3]
            yield a, [1, 3], 1
            yield a, [1, 3], 2
            yield a, [1], -1
            yield a, [1], -2
            yield a, [1], -3
            yield a, np.array([], dtype=np.int64), 0

            a = np.arange(100).reshape(2, -1)
            yield a, 1
            yield a, 1, 0
            yield a, [1], 0
            yield a, 50, 1
            yield a, np.arange(10, 50, 10), 1
            yield a, (1,)
            yield a, (np.int32(4), 10)

            a = np.array([])
            yield a, 1
            yield a, 2
            yield a, (2, 3), 0
            yield a, 1, 0

            a = np.array([[]])
            yield a, 1
            yield a, (2, 3), 1
            yield a, 1, 0
            yield a, 1, 1

        for args in args_variations():
            expected = pyfunc(*args)
            got = cfunc(*args)

            np.testing.assert_equal(expected, list(got))

    def _check_array_split(self, func):
        # array_split specific checks, mainly dealing with `int`s
        pyfunc = func
        cfunc = jit(nopython=True)(pyfunc)

        def args_variations():
            yield np.arange(8), 3
            yield list(np.arange(8)), 3
            yield tuple(np.arange(8)), 3
            yield np.arange(24).reshape(12, 2), 5

        for args in args_variations():
            expected = pyfunc(*args)
            got = cfunc(*args)

            np.testing.assert_equal(expected, list(got))

    def test_array_split_basic(self):
        self._check_split(array_split)
        self._check_array_split(array_split)

    def test_split_basic(self):
        self._check_split(split)

        self.disable_leak_check()  # The exception leaks
        with self.assertRaises(ValueError) as raises:
            njit(split)(np.ones(5), 2)
        self.assertIn(
            "array split does not result in an equal division",
            str(raises.exception)
        )

        with self.assertRaises(ValueError) as raises:
            njit(split)(np.ones(5), [3], axis=-3)
        self.assertIn("np.split: Argument axis out of bounds",
                      str(raises.exception))

    def test_vhdsplit_basic(self):
        # split and array_split have more comprehensive tests of splitting.
        # only do simple tests on vsplit, hsplit and dsplit
        # Based on tests from https://github.com/numpy/numpy/blob/f0befec40376fc46fdaceac2c49c7349ad671bde/numpy/lib/tests/test_shape_base.py#L538-L624 # noqa: E501
        def inputs1D():
            # test_1D_array
            yield np.array([1, 2, 3, 4]), 2
            yield np.array([1., 2., 3., 4.]), 2

        def inputs2D():
            # test_2D_array
            yield np.array([[1, 2, 3, 4], [1, 2, 3, 4]]), 2
            yield np.array([[1., 2., 3., 4.], [1., 2., 3., 4.]]), 2
            yield np.arange(16.0).reshape(4, 4), 2
            yield np.arange(16.0).reshape(4, 4), np.array([3, 6])
            yield np.arange(16.0).reshape(4, 4), [3, 6]
            yield np.arange(16.0).reshape(4, 4), (3, 6)
            yield np.arange(8.0).reshape(2, 2, 2), 2

        def inputs3D():
            # test_3D_array
            np.array([[[1, 2, 3, 4],
                       [1, 2, 3, 4]],
                      [[1, 2, 3, 4],
                       [1, 2, 3, 4]]]), 2
            yield np.arange(16.0).reshape(2, 2, 4), 2
            yield np.arange(16.0).reshape(2, 2, 4), np.array([3, 6])
            yield np.arange(16.0).reshape(2, 2, 4), [3, 6]
            yield np.arange(16.0).reshape(2, 2, 4), (3, 6)
            yield np.arange(8.0).reshape(2, 2, 2), 2

        inputs = [inputs1D(), inputs2D(), inputs3D()]
        for (f, mindim, name) in [(vsplit, 2, "vsplit"),
                                  (hsplit, 1, "hsplit"),
                                  (dsplit, 3, "dsplit")]:
            pyfunc = f
            cfunc = njit(pyfunc)
            for i in range(mindim, 4):
                for a, ind_or_sec in inputs[i - 1]:
                    self.assertPreciseEqual(pyfunc(a, ind_or_sec),
                                            cfunc(a, ind_or_sec))

    def test_vhdsplit_exception(self):
        # Single test method for vsplit, hsplit and dsplit exceptions
        for (f, mindim, name) in [(vsplit, 2, "vsplit"),
                                  (hsplit, 1, "hsplit"),
                                  (dsplit, 3, "dsplit")]:
            cfunc = jit(nopython=True)(f)
            self.disable_leak_check()
            with self.assertRaises(TypingError) as raises:
                cfunc(1, 2)
            self.assertIn('The argument "ary" must be an array',
                          str(raises.exception))
            with self.assertRaises(TypingError) as raises:
                cfunc("abc", 2)
            self.assertIn('The argument "ary" must be an array',
                          str(raises.exception))
            with self.assertRaises(TypingError) as raises:
                cfunc(np.array([[1, 2, 3, 4], [1, 2, 3, 4]]), "abc")
            self.assertIn(('The argument "indices_or_sections" must be int or '
                           '1d-array'),
                          str(raises.exception))
            with self.assertRaises(ValueError) as raises:
                cfunc(np.array(1), 2)
            self.assertIn(name + ' only works on arrays of ' + str(mindim) +
                          ' or more dimensions',
                          str(raises.exception))

    def test_roll_basic(self):
        pyfunc = roll
        cfunc = jit(nopython=True)(pyfunc)

        def a_variations():
            yield np.arange(7)
            yield np.arange(3 * 4 * 5).reshape(3, 4, 5)
            yield [1.1, 2.2, 3.3]
            yield (True, False, True)
            yield False
            yield 4
            yield (9,)
            yield np.asfortranarray(np.array([[1.1, np.nan], [np.inf, 7.8]]))
            yield np.array([])
            yield ()

        def shift_variations():
            return itertools.chain.from_iterable(((True, False),
                                                  range(-10, 10)))

        for a in a_variations():
            for shift in shift_variations():
                expected = pyfunc(a, shift)
                got = cfunc(a, shift)
                self.assertPreciseEqual(expected, got)

    def test_roll_exceptions(self):
        pyfunc = roll
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        for shift in 1.1, (1, 2):
            with self.assertTypingError() as e:
                cfunc(np.arange(10), shift)

            msg = "shift must be an integer"
            assert msg in str(e.exception)

    def test_extract_basic(self):
        pyfunc = extract
        cfunc = jit(nopython=True)(pyfunc)
        _check = partial(self._check_output, pyfunc, cfunc)

        a = np.arange(10)
        self.rnd.shuffle(a)
        for threshold in range(-3, 13):
            cond = a > threshold
            _check({'condition': cond, 'arr': a})

        a = np.arange(60).reshape(4, 5, 3)
        cond = a > 11.2
        _check({'condition': cond, 'arr': a})

        a = ((1, 2, 3), (3, 4, 5), (4, 5, 6))
        cond = np.eye(3).flatten()
        _check({'condition': cond, 'arr': a})

        a = [1.1, 2.2, 3.3, 4.4]
        cond = [1, 1, 0, 1]
        _check({'condition': cond, 'arr': a})

        a = np.linspace(-2, 10, 6)
        element_pool = (True, False, np.nan, -1, -1.0, -1.2, 1, 1.0, 1.5j)
        for cond in itertools.combinations_with_replacement(element_pool, 4):
            _check({'condition': cond, 'arr': a})
            _check({'condition': np.array(cond).reshape(2, 2), 'arr': a})

        a = np.array([1, 2, 3])
        cond = np.array([])
        _check({'condition': cond, 'arr': a})

        a = np.array([1, 2, 3])
        cond = np.array([1, 0, 1, 0])  # but [1, 0, 1, 0, 1] raises
        _check({'condition': cond, 'arr': a})

        a = np.array([[1, 2, 3], [4, 5, 6]])
        cond = [1, 0, 1, 0, 1, 0]  # but [1, 0, 1, 0, 1, 0, 1] raises
        _check({'condition': cond, 'arr': a})

        a = np.array([[1, 2, 3], [4, 5, 6]])
        cond = np.array([1, 0, 1, 0, 1, 0, 0, 0]).reshape(2, 2, 2)
        _check({'condition': cond, 'arr': a})

        a = np.asfortranarray(np.arange(60).reshape(3, 4, 5))
        cond = np.repeat((0, 1), 30)
        _check({'condition': cond, 'arr': a})
        _check({'condition': cond, 'arr': a[::-1]})

        a = np.array(4)
        for cond in 0, 1:
            _check({'condition': cond, 'arr': a})

        a = 1
        cond = 1
        _check({'condition': cond, 'arr': a})

        a = np.array(1)
        cond = np.array([True, False])
        _check({'condition': cond, 'arr': a})

        a = np.arange(4)
        cond = np.array([1, 0, 1, 0, 0, 0]).reshape(2, 3) * 1j
        _check({'condition': cond, 'arr': a})

    def test_extract_exceptions(self):
        pyfunc = extract
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        a = np.array([])
        cond = np.array([1, 2, 3])

        with self.assertRaises(ValueError) as e:
            cfunc(cond, a)
        self.assertIn('Cannot extract from an empty array', str(e.exception))

        def _check(cond, a):
            msg = 'condition shape inconsistent with arr shape'
            with self.assertRaises(ValueError) as e:
                cfunc(cond, a)
            self.assertIn(msg, str(e.exception))

        a = np.array([[1, 2, 3], [1, 2, 3]])
        cond = [1, 0, 1, 0, 1, 0, 1]
        _check(cond, a)

        a = np.array([1, 2, 3])
        cond = np.array([1, 0, 1, 0, 1])
        _check(cond, a)

        a = np.array(60)  # note, this is 0D
        cond = 0, 1
        _check(cond, a)

        a = np.arange(4)
        cond = np.array([True, False, False, False, True])
        _check(cond, a)

        a = np.arange(4)
        cond = np.array([True, False, True, False, False, True, False])
        _check(cond, a)

    @unittest.skipUnless(IS_NUMPY_2, "New in numpy 2.0+")
    def test_np_trapezoid_basic(self):
        self.test_np_trapz_basic(pyfunc=np_trapezoid)

    def test_np_trapz_basic(self, pyfunc=np_trapz):
        cfunc = jit(nopython=True)(pyfunc)
        _check = partial(self._check_output, pyfunc, cfunc)

        y = [1, 2, 3]
        _check({'y': y})

        y = (3, 1, 2, 2, 2)
        _check({'y': y})

        y = np.arange(15).reshape(3, 5)
        _check({'y': y})

        y = np.linspace(-10, 10, 60).reshape(4, 3, 5)
        _check({'y': y}, abs_tol=1e-13)

        self.rnd.shuffle(y)
        _check({'y': y}, abs_tol=1e-13)

        y = np.array([])
        _check({'y': y})

        y = np.array([3.142, np.nan, np.inf, -np.inf, 5])
        _check({'y': y})

        y = np.arange(20) + np.linspace(0, 10, 20) * 1j
        _check({'y': y})

        y = np.array([], dtype=np.complex128)
        _check({'y': y})

        y = (True, False, True)
        _check({'y': y})

    @unittest.skipUnless(IS_NUMPY_2, "New in numpy 2.0+")
    def test_np_trapezoid_x_basic(self):
        self.test_np_trapz_x_basic(pyfunc=np_trapezoid_x)

    def test_np_trapz_x_basic(self, pyfunc=np_trapz_x):
        cfunc = jit(nopython=True)(pyfunc)
        _check = partial(self._check_output, pyfunc, cfunc)

        y = [1, 2, 3]
        x = [4, 6, 8]
        _check({'y': y, 'x': x})

        y = [1, 2, 3, 4, 5]
        x = (4, 6)
        _check({'y': y, 'x': x})

        y = (1, 2, 3, 4, 5)
        x = [4, 5, 6, 7, 8]
        _check({'y': y, 'x': x})

        y = np.array([1, 2, 3, 4, 5])
        x = [4, 4]
        _check({'y': y, 'x': x})

        y = np.array([])
        x = np.array([2, 3])
        _check({'y': y, 'x': x})

        y = (1, 2, 3, 4, 5)
        x = None
        _check({'y': y, 'x': x})

        y = np.arange(20).reshape(5, 4)
        x = np.array([4, 5])
        _check({'y': y, 'x': x})

        y = np.arange(20).reshape(5, 4)
        x = np.array([4, 5, 6, 7])
        _check({'y': y, 'x': x})

        y = np.arange(60).reshape(5, 4, 3)
        x = np.array([4, 5])
        _check({'y': y, 'x': x})

        y = np.arange(60).reshape(5, 4, 3)
        x = np.array([4, 5, 7])
        _check({'y': y, 'x': x})

        y = np.arange(60).reshape(5, 4, 3)
        self.rnd.shuffle(y)
        x = y + 1.1
        self.rnd.shuffle(x)
        _check({'y': y, 'x': x})

        y = np.arange(20)
        x = y + np.linspace(0, 10, 20) * 1j
        _check({'y': y, 'x': x})

        y = np.array([1, 2, 3])
        x = np.array([1 + 1j, 1 + 2j])
        _check({'y': y, 'x': x})

    @unittest.skipUnless(IS_NUMPY_2, "New in numpy 2.0+")
    def test_trapezoid_numpy_questionable(self):
        self.test_trapz_numpy_questionable(pyfunc=np_trapezoid)

    @unittest.skip('NumPy behaviour questionable')
    def test_trapz_numpy_questionable(self, pyfunc=np_trapz):
        # https://github.com/numpy/numpy/issues/12858
        cfunc = jit(nopython=True)(pyfunc)
        _check = partial(self._check_output, pyfunc, cfunc)

        # passes (NumPy and Numba return 2.0)
        y = np.array([True, False, True, True]).astype(int)
        _check({'y': y})

        # fails (NumPy returns 1.5; Numba returns 2.0)
        y = np.array([True, False, True, True])
        _check({'y': y})

    @unittest.skipUnless(IS_NUMPY_2, "New in numpy 2.0+")
    def test_np_trapezoid_dx_basic(self):
        self.test_np_trapz_dx_basic(pyfunc=np_trapezoid_dx)

    def test_np_trapz_dx_basic(self, pyfunc=np_trapz_dx):
        cfunc = jit(nopython=True)(pyfunc)
        _check = partial(self._check_output, pyfunc, cfunc)

        y = [1, 2, 3]
        dx = 2
        _check({'y': y, 'dx': dx})

        y = [1, 2, 3, 4, 5]
        dx = [1, 4, 5, 6]
        _check({'y': y, 'dx': dx})

        y = [1, 2, 3, 4, 5]
        dx = [1, 4, 5, 6]
        _check({'y': y, 'dx': dx})

        y = np.linspace(-2, 5, 10)
        dx = np.nan
        _check({'y': y, 'dx': dx})

        y = np.linspace(-2, 5, 10)
        dx = np.inf
        _check({'y': y, 'dx': dx})

        y = np.linspace(-2, 5, 10)
        dx = np.linspace(-2, 5, 9)
        _check({'y': y, 'dx': dx}, abs_tol=1e-13)

        y = np.arange(60).reshape(4, 5, 3) * 1j
        dx = np.arange(40).reshape(4, 5, 2)
        _check({'y': y, 'dx': dx})

        x = np.arange(-10, 10, .1)
        r = cfunc(np.exp(-.5 * x ** 2) / np.sqrt(2 * np.pi), dx=0.1)
        # check integral of normal equals 1
        np.testing.assert_almost_equal(r, 1, 7)

        y = np.arange(20)
        dx = 1j
        _check({'y': y, 'dx': dx})

        y = np.arange(20)
        dx = np.array([5])
        _check({'y': y, 'dx': dx})

    @unittest.skipUnless(IS_NUMPY_2, "New in numpy 2.0+")
    def test_np_trapezoid_x_dx_basic(self):
        self.test_np_trapz_x_dx_basic(pyfunc=np_trapezoid_x_dx)

    def test_np_trapz_x_dx_basic(self, pyfunc=np_trapz_x_dx):
        cfunc = jit(nopython=True)(pyfunc)
        _check = partial(self._check_output, pyfunc, cfunc)

        # dx should be ignored
        for dx in (None, 2, np.array([1, 2, 3, 4, 5])):
            y = [1, 2, 3]
            x = [4, 6, 8]
            _check({'y': y, 'x': x, 'dx': dx})

            y = [1, 2, 3, 4, 5]
            x = [4, 6]
            _check({'y': y, 'x': x, 'dx': dx})

            y = [1, 2, 3, 4, 5]
            x = [4, 5, 6, 7, 8]
            _check({'y': y, 'x': x, 'dx': dx})

            y = np.arange(60).reshape(4, 5, 3)
            self.rnd.shuffle(y)
            x = y * 1.1
            x[2, 2, 2] = np.nan
            _check({'y': y, 'x': x, 'dx': dx})

    @unittest.skipUnless(IS_NUMPY_2, "New in numpy 2.0+")
    def test_np_trapezoid_x_dx_exceptions(self):
        self.test_np_trapz_x_dx_exceptions(pyfunc=np_trapezoid_x_dx)

    def test_np_trapz_x_dx_exceptions(self, pyfunc=np_trapz_x_dx):
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        def check_not_ok(params):
            with self.assertRaises(ValueError) as e:
                cfunc(*params)

            self.assertIn('unable to broadcast', str(e.exception))

        y = [1, 2, 3, 4, 5]
        for x in [4, 5, 6, 7, 8, 9], [4, 5, 6]:
            check_not_ok((y, x, 1.0))

        y = np.arange(60).reshape(3, 4, 5)
        x = np.arange(36).reshape(3, 4, 3)
        check_not_ok((y, x, 1.0))

        y = np.arange(60).reshape(3, 4, 5)
        x = np.array([4, 5, 6, 7])
        check_not_ok((y, x, 1.0))

        y = [1, 2, 3, 4, 5]
        dx = np.array([1.0, 2.0])
        check_not_ok((y, None, dx))

        y = np.arange(60).reshape(3, 4, 5)
        dx = np.arange(60).reshape(3, 4, 5)
        check_not_ok((y, None, dx))

        with self.assertTypingError() as e:
            y = np.array(4)
            check_not_ok((y, None, 1.0))

        self.assertIn('y cannot be 0D', str(e.exception))

        for y in 5, False, np.nan:
            with self.assertTypingError() as e:
                cfunc(y, None, 1.0)

            self.assertIn('y cannot be a scalar', str(e.exception))

    def test_average(self):

        #array of random numbers
        N = 100
        a = np.random.ranf(N) * 100
        w = np.random.ranf(N) * 100
        w0 = np.zeros(N)

        #boolean array and weights
        a_bool = np.random.ranf(N) > 0.5
        w_bool = np.random.ranf(N) > 0.5

        #array of random ints
        a_int = np.random.randint(101, size=N)
        w_int = np.random.randint(101, size=N)

        #3D array of random numbers
        d0 = 100
        d1 = 50
        d2 = 25
        a_3d = np.random.rand(d0,d1,d2) * 100
        w_3d = np.random.rand(d0,d1,d2) * 100

        pyfunc = np_average
        cfunc = jit(nopython=True)(pyfunc)

        #test case for average with weights
        #(number of elements in array and weight array are equal)
        self.assertAlmostEqual( pyfunc(a,weights=w),
                                cfunc(a,weights=w), places=10)
        self.assertAlmostEqual( pyfunc(a_3d,weights=w_3d),
                                cfunc(a_3d,weights=w_3d), places=10)

        #test case for average with array and weights with
        #int datatype (number of elements in array and weight array are equal)
        self.assertAlmostEqual( pyfunc(a_int,weights=w_int),
                                cfunc(a_int,weights=w_int), places=10)

        #test case for average with boolean weights
        self.assertAlmostEqual( pyfunc(a,weights=w_bool),
                                cfunc(a,weights=w_bool), places=10)
        self.assertAlmostEqual( pyfunc(a_bool,weights=w),
                                cfunc(a_bool,weights=w), places=10)
        self.assertAlmostEqual( pyfunc(a_bool, weights=w_bool),
                                cfunc(a_bool, weights=w_bool), places=10)

        #test case for average without weights
        self.assertAlmostEqual(pyfunc(a), cfunc(a), places=10)
        self.assertAlmostEqual(pyfunc(a_3d), cfunc(a_3d), places=10)

        def test_weights_zero_sum(data, weights):
            with self.assertRaises(ZeroDivisionError) as e:
                cfunc(data, weights=weights)
            err = e.exception
            self.assertEqual(str(err),
                             "Weights sum to zero, can't be normalized.")

        #test case when sum of weights is zero
        test_weights_zero_sum(a, weights=w0)

        def test_1D_weights(data, weights):
            with self.assertRaises(TypeError) as e:
                cfunc(data, weights=weights)
            err = e.exception
            self.assertEqual(str(err),
                             "Numba does not support average when shapes of "
                             "a and weights differ.")

        def test_1D_weights_axis(data, axis, weights):
            with self.assertRaises(TypeError) as e:
                cfunc(data,axis=axis, weights=weights)
            err = e.exception
            self.assertEqual(str(err),
                             "Numba does not support average with axis.")

        #small case to test exceptions for 2D array and 1D weights
        data = np.arange(6).reshape((3,2,1))
        w = np.asarray([1. / 4, 3. / 4])

        #test without axis argument
        test_1D_weights(data, weights=w)

        #test with axis argument
        test_1D_weights_axis(data, axis=1, weights=w)

    def test_allclose(self):

        pyfunc = np_allclose
        cfunc = jit(nopython=True)(pyfunc)

        min_int = np.iinfo(np.int_).min
        a = np.array([min_int], dtype=np.int_)

        simple_data = [
            (np.asarray([1e10, 1e-7]), np.asarray([1.00001e10, 1e-8])),
            (np.asarray([1e10, 1e-8]), np.asarray([1.00001e10, 1e-9])),
            (np.asarray([1e10, 1e-8]), np.asarray([1.0001e10, 1e-9])),
            (np.asarray([1e10]), np.asarray([1.0001e10, 1e-9])),
            (1.0, 1.0),
            (np.array([np.inf, 1]), np.array([0, np.inf])),
            (a, a)
        ]

        for a, b in simple_data:
            py_result = pyfunc(a, b)
            c_result = cfunc(a, b)
            self.assertEqual(py_result, c_result)

        a = np.asarray([1.0, np.nan])
        b = np.asarray([1.0, np.nan])
        self.assertFalse(cfunc(a, b))
        self.assertEqual(pyfunc(a, b, equal_nan=True),
                         cfunc(a, b, equal_nan=True))

        b = np.asarray([np.nan, 1.0])
        self.assertEqual(pyfunc(a, b), cfunc(a, b))

        noise_levels = [1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0.0]
        zero_array = np.zeros((25, 4))
        a = np.random.ranf((25, 4))
        for noise in noise_levels:
            for rtol in noise_levels:
                for atol in noise_levels:
                    py_result = pyfunc(zero_array, noise,
                                       atol=atol, rtol=rtol)
                    c_result = cfunc(zero_array, noise,
                                     atol=atol, rtol=rtol)
                    self.assertEqual(py_result, c_result)

                    py_result = pyfunc(noise, zero_array,
                                       atol=atol, rtol=rtol)
                    c_result = cfunc(noise, zero_array,
                                     atol=atol, rtol=rtol)
                    self.assertEqual(py_result, c_result)

                    py_result = pyfunc(np.asarray([noise]), zero_array,
                                       atol=atol, rtol=rtol)
                    c_result = cfunc(np.asarray([noise]), zero_array,
                                     atol=atol, rtol=rtol)
                    self.assertEqual(py_result, c_result)

                    py_result = pyfunc(a, a + noise, atol=atol, rtol=rtol)
                    c_result = cfunc(a, a + noise, atol=atol, rtol=rtol)
                    self.assertEqual(py_result, c_result)

                    py_result = pyfunc(a + noise, a, atol=atol, rtol=rtol)
                    c_result = cfunc(a + noise, a, atol=atol, rtol=rtol)
                    self.assertEqual(py_result, c_result)

    def test_ip_allclose_numpy(self):
        # https://github.com/numpy/numpy/blob/4adc87dff15a247e417d50f10cc4def8e1c17a03/numpy/core/tests/test_numeric.py#L2402-L2420    # noqa: E501
        pyfunc = np_allclose
        cfunc = jit(nopython=True)(pyfunc)

        arr = np.array([100.0, 1000.0])
        aran = np.arange(125).astype(dtype=np.float64).reshape((5, 5, 5))

        atol = 1e-8
        rtol = 1e-5

        numpy_data = [
            (np.asarray([1, 0]), np.asarray([1, 0])),
            (np.asarray([atol]), np.asarray([0.0])),
            (np.asarray([1.0]), np.asarray([1 + rtol + atol])),
            (arr, arr + arr * rtol),
            (arr, arr + arr * rtol + atol * 2),
            (aran, aran + aran * rtol),
            (np.inf, np.inf),
            (np.inf, np.asarray([np.inf]))
        ]

        for (x, y) in numpy_data:
            self.assertEqual(pyfunc(x, y), cfunc(x, y))

    def test_ip_not_allclose_numpy(self):
        # https://github.com/numpy/numpy/blob/4adc87dff15a247e417d50f10cc4def8e1c17a03/numpy/core/tests/test_numeric.py#L2422-L2441    # noqa: E501

        pyfunc = np_allclose
        cfunc = jit(nopython=True)(pyfunc)

        aran = np.arange(125).astype(dtype=np.float64).reshape((5, 5, 5))

        atol = 1e-8
        rtol = 1e-5

        numpy_data = [
            (np.asarray([np.inf, 0]), np.asarray([1.0, np.inf])),
            (np.asarray([np.inf, 0]), np.asarray([1.0, 0])),
            (np.asarray([np.inf, np.inf]), np.asarray([1.0, np.inf])),
            (np.asarray([np.inf, np.inf]), np.asarray([1.0, 0.0])),
            (np.asarray([-np.inf, 0.0]), np.asarray([np.inf, 0.0])),
            (np.asarray([np.nan, 0.0]), np.asarray([np.nan, 0.0])),
            (np.asarray([atol * 2]), np.asarray([0.0])),
            (np.asarray([1.0]), np.asarray([1 + rtol + atol * 2])),
            (aran, aran + aran * atol + atol * 2),
            (np.array([np.inf, 1.0]), np.array([0.0, np.inf]))
        ]

        for (x, y) in numpy_data:
            self.assertEqual(pyfunc(x, y), cfunc(x, y))

    def test_return_class_is_ndarray_numpy(self):
        # https://github.com/numpy/numpy/blob/4adc87dff15a247e417d50f10cc4def8e1c17a03/numpy/core/tests/test_numeric.py#L2460-L2468    # noqa: E501

        pyfunc = np_allclose
        cfunc = jit(nopython=True)(pyfunc)

        class Foo(np.ndarray):
            def __new__(cls, *args, **kwargs):
                return np.array(*args, **kwargs).view(cls)

        a = Foo([1])
        self.assertTrue(type(cfunc(a, a)) is bool)

    def test_equalnan_numpy(self):
        # https://github.com/numpy/numpy/blob/4adc87dff15a247e417d50f10cc4def8e1c17a03/numpy/core/tests/test_numeric.py#L2456-L2458    # noqa: E501
        pyfunc = np_allclose
        cfunc = jit(nopython=True)(pyfunc)

        x = np.array([1.0, np.nan])

        self.assertEqual(pyfunc(x, x, equal_nan=True),
                         cfunc(x, x, equal_nan=True))

    def test_no_parameter_modification_numpy(self):
        # https://github.com/numpy/numpy/blob/4adc87dff15a247e417d50f10cc4def8e1c17a03/numpy/core/tests/test_numeric.py#L2443-L2448    # noqa: E501

        pyfunc = np_allclose
        cfunc = jit(nopython=True)(pyfunc)

        x = np.array([np.inf, 1])
        y = np.array([0, np.inf])

        cfunc(x, y)
        np.testing.assert_array_equal(x, np.array([np.inf, 1]))
        np.testing.assert_array_equal(y, np.array([0, np.inf]))

    def test_min_int_numpy(self):
        # https://github.com/numpy/numpy/blob/4adc87dff15a247e417d50f10cc4def8e1c17a03/numpy/core/tests/test_numeric.py#L2450-L2454    # noqa: E501

        pyfunc = np_allclose
        cfunc = jit(nopython=True)(pyfunc)

        min_int = np.iinfo(np.int_).min
        a = np.array([min_int], dtype=np.int_)

        self.assertEqual(pyfunc(a, a), cfunc(a, a))

    def test_allclose_exception(self):
        self.disable_leak_check()

        pyfunc = np_allclose
        cfunc = jit(nopython=True)(pyfunc)

        inps = [
            (np.asarray([1e10, 1e-9, np.nan]),
             np.asarray([1.0001e10, 1e-9]),
             1e-05, 1e-08, False,
             "shape mismatch: objects cannot be broadcast to a single shape",
             ValueError),
            ('hello', 3, False, 1e-08, False,
             'The first argument "a" must be array-like',
             TypingError),
            (3, 'hello', False, 1e-08, False,
             'The second argument "b" must be array-like',
             TypingError),
            (2, 3, False, 1e-08, False,
             'The third argument "rtol" must be a floating point',
             TypingError),
            (2, 3, 1e-05, False, False,
             'The fourth argument "atol" must be a floating point',
             TypingError),
            (2, 3, 1e-05, 1e-08, 1,
             'The fifth argument "equal_nan" must be a boolean',
             TypingError),
        ]

        for a, b, rtol, atol, equal_nan, exc_msg, exc in inps:
            with self.assertRaisesRegex(exc, exc_msg):
                cfunc(a, b, rtol, atol, equal_nan)

    def test_interp_basic(self):
        pyfunc = interp
        cfunc = jit(nopython=True)(pyfunc)
        _check = partial(self._check_output, pyfunc, cfunc, abs_tol=1e-10)

        x = np.linspace(-5, 5, 25)
        xp = np.arange(-4, 8)
        fp = xp + 1.5
        _check(params={'x': x, 'xp': xp, 'fp': fp})
        self.rnd.shuffle(x)
        _check(params={'x': x, 'xp': xp, 'fp': fp})
        self.rnd.shuffle(fp)
        _check(params={'x': x, 'xp': xp, 'fp': fp})

        x[:5] = np.nan
        x[-5:] = np.inf
        self.rnd.shuffle(x)
        _check(params={'x': x, 'xp': xp, 'fp': fp})
        fp[:5] = np.nan
        fp[-5:] = -np.inf
        self.rnd.shuffle(fp)
        _check(params={'x': x, 'xp': xp, 'fp': fp})

        x = np.arange(-4, 8)
        xp = x + 1
        fp = x + 2
        _check(params={'x': x, 'xp': xp, 'fp': fp})

        x = (2.2, 3.3, -5.0)
        xp = (2, 3, 4)
        fp = (5, 6, 7)
        _check(params={'x': x, 'xp': xp, 'fp': fp})

        x = ((2.2, 3.3, -5.0), (1.2, 1.3, 4.0))
        xp = np.linspace(-4, 4, 10)
        fp = np.arange(-5, 5)
        _check(params={'x': x, 'xp': xp, 'fp': fp})

        x = np.array([1.4, np.nan, np.inf, -np.inf, 0.0, -9.1])
        x = x.reshape(3, 2, order='F')
        xp = np.linspace(-4, 4, 10)
        fp = np.arange(-5, 5)
        _check(params={'x': x, 'xp': xp, 'fp': fp})

        for x in range(-2, 4):
            xp = [0, 1, 2]
            fp = (3, 4, 5)
            _check(params={'x': x, 'xp': xp, 'fp': fp})

        x = np.array([])
        xp = [0, 1, 2]
        fp = (3, 4, 5)
        _check(params={'x': x, 'xp': xp, 'fp': fp})

        x = np.linspace(0, 25, 60).reshape(3, 4, 5)
        xp = np.arange(20)
        fp = xp - 10
        _check(params={'x': x, 'xp': xp, 'fp': fp})

        x = np.nan
        xp = np.arange(5)
        fp = np.full(5, np.nan)
        _check(params={'x': x, 'xp': xp, 'fp': fp})

        x = np.nan
        xp = [3]
        fp = [4]
        _check(params={'x': x, 'xp': xp, 'fp': fp})

        x = np.arange(-4, 8)
        xp = x
        fp = x
        _check(params={'x': x, 'xp': xp, 'fp': fp})

        x = [True, False]
        xp = np.arange(-4, 8)
        fp = xp
        _check(params={'x': x, 'xp': xp, 'fp': fp})

        x = [-np.inf, -1.0, 0.0, 1.0, np.inf]
        xp = np.arange(-4, 8)
        fp = xp * 2.2
        _check(params={'x': x, 'xp': xp, 'fp': fp})

        x = np.linspace(-10, 10, 10)
        xp = np.array([-np.inf, -1.0, 0.0, 1.0, np.inf])
        fp = xp * 2.2
        _check(params={'x': x, 'xp': xp, 'fp': fp})

        x = self.rnd.randn(100)
        xp = np.linspace(-3, 3, 100)
        fp = np.full(100, fill_value=3.142)
        _check(params={'x': x, 'xp': xp, 'fp': fp})

        for factor in 1, -1:
            x = np.array([5, 6, 7]) * factor
            xp = [1, 2]
            fp = [3, 4]
            _check(params={'x': x, 'xp': xp, 'fp': fp})

        x = 1
        xp = [1]
        fp = [True]
        _check(params={'x': x, 'xp': xp, 'fp': fp})

        x = np.linspace(0, 1, 5)
        y = np.linspace(0, 1, 5)
        x0 = np.linspace(0, 1, 50)
        out = cfunc(x0, x, y)
        np.testing.assert_almost_equal(out, x0)

        x = np.array([1, 2, 3, 4])
        xp = np.array([1, 2, 3, 4])
        fp = np.array([1, 2, 3.01, 4])
        _check(params={'x': x, 'xp': xp, 'fp': fp})

        xp = [1]
        fp = [np.inf]
        _check(params={'x': 1, 'xp': xp, 'fp': fp})

        x = np.array([1, 2, 2.5, 3, 4])
        xp = np.array([1, 2, 3, 4])
        fp = np.array([1, 2, np.nan, 4])
        _check({'x': x, 'xp': xp, 'fp': fp})

        x = np.array([1, 1.5, 2, 2.5, 3, 4, 4.5, 5, 5.5])
        xp = np.array([1, 2, 3, 4, 5])
        fp = np.array([np.nan, 2, np.nan, 4, np.nan])
        _check({'x': x, 'xp': xp, 'fp': fp})

        x = np.array([1, 2, 2.5, 3, 4])
        xp = np.array([1, 2, 3, 4])
        fp = np.array([1, 2, np.inf, 4])
        _check({'x': x, 'xp': xp, 'fp': fp})

        x = np.array([1, 1.5, np.nan, 2.5, -np.inf, 4, 4.5, 5, np.inf, 0, 7])
        xp = np.array([1, 2, 3, 4, 5, 6])
        fp = np.array([1, 2, np.nan, 4, 3, np.inf])
        _check({'x': x, 'xp': xp, 'fp': fp})

        x = np.array([3.10034867, 3.0999066, 3.10001529])
        xp = np.linspace(0, 10, 1 + 20000)
        fp = np.sin(xp / 2.0)
        _check({'x': x, 'xp': xp, 'fp': fp})

        x = self.rnd.uniform(0, 2 * np.pi, (100,))
        xp = np.linspace(0, 2 * np.pi, 1000)
        fp = np.cos(xp)
        exact = np.cos(x)
        got = cfunc(x, xp, fp)
        np.testing.assert_allclose(exact, got, atol=1e-5)

        # very dense calibration
        x = self.rnd.randn(10)
        xp = np.linspace(-10, 10, 1000)
        fp = np.ones_like(xp)
        _check({'x': x, 'xp': xp, 'fp': fp})

        # very sparse calibration
        x = self.rnd.randn(1000)
        xp = np.linspace(-10, 10, 10)
        fp = np.ones_like(xp)
        _check({'x': x, 'xp': xp, 'fp': fp})

    def _make_some_values_non_finite(self, a):
        p = a.size // 100
        np.put(a, self.rnd.choice(range(a.size), p, replace=False), np.nan)
        np.put(a, self.rnd.choice(range(a.size), p, replace=False), -np.inf)
        np.put(a, self.rnd.choice(range(a.size), p, replace=False), np.inf)

    def arrays(self, ndata):
        # much_finer_grid
        yield np.linspace(2.0, 7.0, 1 + ndata * 5)
        # finer_grid
        yield np.linspace(2.0, 7.0, 1 + ndata)
        # similar_grid
        yield np.linspace(2.1, 6.8, 1 + ndata // 2)
        # coarser_grid
        yield np.linspace(2.1, 7.5, 1 + ndata // 2)
        # much_coarser_grid
        yield np.linspace(1.1, 9.5, 1 + ndata // 5)
        # finer_stretched_grid
        yield np.linspace(3.1, 5.3, 1 + ndata) * 1.09
        # similar_stretched_grid
        yield np.linspace(3.1, 8.3, 1 + ndata // 2) * 1.09
        # finer_compressed_grid
        yield np.linspace(3.1, 5.3, 1 + ndata) * 0.91
        # similar_compressed_grid
        yield np.linspace(3.1, 8.3, 1 + ndata // 2) * 0.91
        # warped_grid
        yield np.linspace(3.1, 5.3, 1 + ndata // 2) + 0.3 * np.sin(
            np.arange(1 + ndata / 2) * np.pi / (1 + ndata / 2))
        # very_low_noise_grid
        yield np.linspace(3.1, 5.3, 1 + ndata) + self.rnd.normal(
            size=1 + ndata, scale=0.5 / ndata)
        # low_noise_grid
        yield np.linspace(3.1, 5.3, 1 + ndata) + self.rnd.normal(
            size=1 + ndata, scale=2.0 / ndata)
        # med_noise_grid
        yield np.linspace(3.1, 5.3, 1 + ndata) + self.rnd.normal(
            size=1 + ndata, scale=5.0 / ndata)
        # high_noise_grid
        yield np.linspace(3.1, 5.3, 1 + ndata) + self.rnd.normal(
            size=1 + ndata, scale=20.0 / ndata)
        # very_high_noise_grid
        yield np.linspace(3.1, 5.3, 1 + ndata) + self.rnd.normal(
            size=1 + ndata, scale=50.0 / ndata)
        # extreme_noise_grid
        yield np.linspace(3.1, 5.3, 1 + ndata) + self.rnd.normal(
            size=1 + ndata, scale=200.0 / ndata)
        # random_fine_grid
        yield self.rnd.rand(1 + ndata) * 9.0 + 0.6
        # random_grid
        yield self.rnd.rand(1 + ndata * 2) * 4.0 + 1.3

    def test_interp_stress_tests(self):
        pyfunc = interp
        cfunc = jit(nopython=True)(pyfunc)

        ndata = 20000
        xp = np.linspace(0, 10, 1 + ndata)
        fp = np.sin(xp / 2.0)

        for x in self.arrays(ndata):
            atol = 1e-14  # using abs_tol as otherwise fails on 32bit builds

            expected = pyfunc(x, xp, fp)
            got = cfunc(x, xp, fp)
            self.assertPreciseEqual(expected, got, abs_tol=atol)

            # no longer require xp to be monotonically increasing
            # (in keeping with numpy) even if the output might not
            # be meaningful; shuffle all inputs
            self.rnd.shuffle(x)
            expected = pyfunc(x, xp, fp)
            got = cfunc(x, xp, fp)
            self.assertPreciseEqual(expected, got, abs_tol=atol)

            self.rnd.shuffle(xp)
            expected = pyfunc(x, xp, fp)
            got = cfunc(x, xp, fp)
            self.assertPreciseEqual(expected, got, abs_tol=atol)

            self.rnd.shuffle(fp)
            expected = pyfunc(x, xp, fp)
            got = cfunc(x, xp, fp)
            self.assertPreciseEqual(expected, got, abs_tol=atol)

            # add some values non finite
            self._make_some_values_non_finite(x)
            expected = pyfunc(x, xp, fp)
            got = cfunc(x, xp, fp)
            self.assertPreciseEqual(expected, got, abs_tol=atol)

            self._make_some_values_non_finite(xp)
            expected = pyfunc(x, xp, fp)
            got = cfunc(x, xp, fp)
            self.assertPreciseEqual(expected, got, abs_tol=atol)

            self._make_some_values_non_finite(fp)
            expected = pyfunc(x, xp, fp)
            got = cfunc(x, xp, fp)
            self.assertPreciseEqual(expected, got, abs_tol=atol)

    @unittest.skipIf(IS_NUMPY_2 and IS_MACOS_ARM64, "NEP 50 interaction issue.")
    def test_interp_complex_stress_tests(self):
        pyfunc = interp
        cfunc = jit(nopython=True)(pyfunc)

        ndata = 2000
        xp = np.linspace(0, 10, 1 + ndata)

        real = np.sin(xp / 2.0)
        real[:200] = self.rnd.choice([np.inf, -np.inf, np.nan], 200)
        self.rnd.shuffle(real)

        imag = np.cos(xp / 2.0)
        imag[:200] = self.rnd.choice([np.inf, -np.inf, np.nan], 200)
        self.rnd.shuffle(imag)

        fp = real + 1j * imag

        for x in self.arrays(ndata):
            expected = pyfunc(x, xp, fp)
            got = cfunc(x, xp, fp)
            np.testing.assert_allclose(expected, got, equal_nan=True)

            self.rnd.shuffle(x)
            self.rnd.shuffle(xp)
            self.rnd.shuffle(fp)
            np.testing.assert_allclose(expected, got, equal_nan=True)

    def test_interp_exceptions(self):
        pyfunc = interp
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        x = np.array([1, 2, 3])
        xp = np.array([])
        fp = np.array([])

        with self.assertRaises(ValueError) as e:
            cfunc(x, xp, fp)

        msg = "array of sample points is empty"
        self.assertIn(msg, str(e.exception))

        x = 1
        xp = np.array([1, 2, 3])
        fp = np.array([1, 2])

        with self.assertRaises(ValueError) as e:
            cfunc(x, xp, fp)

        msg = "fp and xp are not of the same size."
        self.assertIn(msg, str(e.exception))

        x = 1
        xp = np.arange(6).reshape(3, 2)
        fp = np.arange(6)

        with self.assertTypingError() as e:
            cfunc(x, xp, fp)

        msg = "xp must be 1D"
        self.assertIn(msg, str(e.exception))

        x = 1
        xp = np.arange(6)
        fp = np.arange(6).reshape(3, 2)

        with self.assertTypingError() as e:
            cfunc(x, xp, fp)

        msg = "fp must be 1D"
        self.assertIn(msg, str(e.exception))

        x = 1 + 1j
        xp = np.arange(6)
        fp = np.arange(6)

        with self.assertTypingError() as e:
            cfunc(x, xp, fp)

        complex_dtype_msg = (
            "Cannot cast array data from complex dtype "
            "to float64 dtype"
        )
        self.assertIn(complex_dtype_msg, str(e.exception))

        x = 1
        xp = (np.arange(6) + 1j).astype(np.complex64)
        fp = np.arange(6)

        with self.assertTypingError() as e:
            cfunc(x, xp, fp)

        self.assertIn(complex_dtype_msg, str(e.exception))

    def test_interp_non_finite_calibration(self):
        # examples from
        # https://github.com/numpy/numpy/issues/12951
        pyfunc = interp
        cfunc = jit(nopython=True)(pyfunc)
        _check = partial(self._check_output, pyfunc, cfunc)

        xp = np.array([0, 1, 9, 10])
        fp = np.array([-np.inf, 0.1, 0.9, np.inf])
        x = np.array([0.2, 9.5])
        params = {'x': x, 'xp': xp, 'fp': fp}
        _check(params)

        xp = np.array([-np.inf, 1, 9, np.inf])
        fp = np.array([0, 0.1, 0.9, 1])
        x = np.array([0.2, 9.5])
        params = {'x': x, 'xp': xp, 'fp': fp}
        _check(params)

    def test_interp_supplemental_tests(self):
        # inspired by class TestInterp
        # https://github.com/numpy/numpy/blob/f5b6850f231/numpy/lib/tests/test_function_base.py    # noqa: E501
        pyfunc = interp
        cfunc = jit(nopython=True)(pyfunc)

        for size in range(1, 10):
            xp = np.arange(size, dtype=np.double)
            yp = np.ones(size, dtype=np.double)
            incpts = np.array([-1, 0, size - 1, size], dtype=np.double)
            decpts = incpts[::-1]

            incres = cfunc(incpts, xp, yp)
            decres = cfunc(decpts, xp, yp)
            inctgt = np.array([1, 1, 1, 1], dtype=float)
            dectgt = inctgt[::-1]
            np.testing.assert_almost_equal(incres, inctgt)
            np.testing.assert_almost_equal(decres, dectgt)

        x = np.linspace(0, 1, 5)
        y = np.linspace(0, 1, 5)
        x0 = 0
        np.testing.assert_almost_equal(cfunc(x0, x, y), x0)
        x0 = 0.3
        np.testing.assert_almost_equal(cfunc(x0, x, y), x0)
        x0 = np.float32(0.3)
        np.testing.assert_almost_equal(cfunc(x0, x, y), x0)
        x0 = np.float64(0.3)
        np.testing.assert_almost_equal(cfunc(x0, x, y), x0)
        x0 = np.nan
        np.testing.assert_almost_equal(cfunc(x0, x, y), x0)

        x = np.linspace(0, 1, 5)
        y = np.linspace(0, 1, 5)
        x0 = np.array(0.3)
        np.testing.assert_almost_equal(cfunc(x0, x, y), x0)

        xp = np.arange(0, 10, 0.0001)
        fp = np.sin(xp)
        np.testing.assert_almost_equal(cfunc(np.pi, xp, fp), 0.0)

    def test_interp_supplemental_complex_tests(self):
        # inspired by class TestInterp
        # https://github.com/numpy/numpy/blob/f5b6850f231/numpy/lib/tests/test_function_base.py    # noqa: E501
        pyfunc = interp
        cfunc = jit(nopython=True)(pyfunc)

        x = np.linspace(0, 1, 5)
        y = np.linspace(0, 1, 5) + (1 + np.linspace(0, 1, 5)) * 1.0j
        x0 = 0.3
        y0 = x0 + (1 + x0) * 1.0j
        np.testing.assert_almost_equal(cfunc(x0, x, y), y0)

    def test_interp_float_precision_handled_per_numpy(self):
        # test cases from https://github.com/numba/numba/issues/4890
        pyfunc = interp
        cfunc = jit(nopython=True)(pyfunc)
        dtypes = [np.float32, np.float64, np.int32, np.int64]

        for combo in itertools.combinations_with_replacement(dtypes, 3):
            xp_dtype, fp_dtype, x_dtype = combo
            xp = np.arange(10, dtype=xp_dtype)
            fp = (xp ** 2).astype(fp_dtype)
            x = np.linspace(2, 3, 10, dtype=x_dtype)

            expected = pyfunc(x, xp, fp)
            got = cfunc(x, xp, fp)
            self.assertPreciseEqual(expected, got)

    def test_isnat(self):
        def values():
            yield np.datetime64("2016-01-01")
            yield np.datetime64("NaT")
            yield np.datetime64('NaT', 'ms')
            yield np.datetime64('NaT', 'ns')
            yield np.datetime64('2038-01-19T03:14:07')

            yield np.timedelta64('NaT', "ms")
            yield np.timedelta64(34, "ms")

            for unit in ['Y', 'M', 'W', 'D',
                         'h', 'm', 's', 'ms', 'us',
                         'ns', 'ps', 'fs', 'as']:
                yield np.array([123, -321, "NaT"],
                               dtype='<datetime64[%s]' % unit)
                yield np.array([123, -321, "NaT"],
                               dtype='<timedelta64[%s]' % unit)

        pyfunc = isnat
        cfunc = jit(nopython=True)(pyfunc)

        for x in values():
            expected = pyfunc(x)
            got = cfunc(x)
            if isinstance(x, np.ndarray):
                self.assertPreciseEqual(expected, got, (x,))
            else:
                self.assertEqual(expected, got, x)

    def test_asarray(self):

        def input_variations():
            """
            To quote from: https://docs.scipy.org/doc/numpy/reference/generated/numpy.asarray.html    # noqa: E501
            Input data, in any form that can be converted to an array.
            This includes:
            * lists
            * lists of tuples
            * tuples
            * tuples of tuples
            * tuples of lists
            * ndarrays
            """
            yield 1j
            yield 1.2
            yield False
            yield 1
            yield [1, 2, 3]
            yield [(1, 2, 3), (1, 2, 3)]
            yield (1, 2, 3)
            yield ((1, 2, 3), (1, 2, 3))
            yield ([1, 2, 3], [1, 2, 3])
            yield np.array([])
            yield np.arange(4)
            yield np.arange(12).reshape(3, 4)
            yield np.arange(12).reshape(3, 4).T

            # Test cases for `numba.typed.List`
            def make_list(values):
                a = List()
                for i in values:
                    a.append(i)
                return a
            yield make_list((1, 2, 3))
            yield make_list((1.0, 2.0, 3.0))
            yield make_list((1j, 2j, 3j))
            yield make_list((True, False, True))

        # used to check that if the input is already an array and the dtype is
        # the same as that of the input/omitted then the array itself is
        # returned.
        def check_pass_through(jitted, expect_same, params):
            returned = jitted(**params)
            if expect_same:
                self.assertTrue(returned is params['a'])
            else:
                self.assertTrue(returned is not params['a'])
                # should be numerically the same, just different dtype
                np.testing.assert_allclose(returned, params['a'])
                self.assertTrue(returned.dtype == params['dtype'])

        for pyfunc in [asarray, asarray_kws]:
            cfunc = jit(nopython=True)(pyfunc)
            _check = partial(self._check_output, pyfunc, cfunc)

            for x in input_variations():
                params = {'a': x}
                if 'kws' in pyfunc.__name__:
                    for dt in [None, np.complex128]:
                        params['dtype'] = dt
                        _check(params)
                else:
                    _check(params)

                # check the behaviour over a dtype change (or not!)
                x = np.arange(10, dtype=np.float32)
                params = {'a': x}
                if 'kws' in pyfunc.__name__:
                    params['dtype'] = None
                    check_pass_through(cfunc, True, params)
                    params['dtype'] = np.complex128
                    check_pass_through(cfunc, False, params)
                    params['dtype'] = np.float32
                    check_pass_through(cfunc, True, params)
                else:
                    check_pass_through(cfunc, True, params)

    def test_asarray_literal(self):

        def case1():
            return np.asarray("hello world")

        def case2(): # kind1
            s = "hello world"
            return np.asarray(s)

        def case3(): # kind2
            s = '大处 着眼，小处着手。大大大处'
            return np.asarray(s)

        def case4():
            s = ''
            return np.asarray(s)

        funcs = [case1, case2, case3, case4]

        for pyfunc in funcs:
            cfunc = jit(nopython=True)(pyfunc)
            expected = pyfunc()
            got = cfunc()
            self.assertPreciseEqual(expected, got)

    def test_asarray_rejects_List_with_illegal_dtype(self):
        self.disable_leak_check()
        cfunc = jit(nopython=True)(asarray)

        def test_reject(alist):
            with self.assertRaises(TypingError) as e:
                cfunc(alist)
            self.assertIn(
                "asarray support for List is limited "
                "to Boolean and Number types",
                str(e.exception))

        def make_none_typed_list():
            l = List()
            l.append(None)
            return l

        def make_nested_list():
            l = List()
            m = List()
            m.append(1)
            l.append(m)
            return l

        def make_nested_list_with_dict():
            l = List()
            d = Dict()
            d[1] = "a"
            l.append(d)
            return l

        def make_unicode_list():
            l = List()
            for i in ("a", "bc", "def"):
                l.append(i)
            return l

        test_reject(make_none_typed_list())
        test_reject(make_nested_list())
        test_reject(make_nested_list_with_dict())
        test_reject(make_unicode_list())

    @skip_if_numpy_2
    def test_asfarray(self):
        def inputs():
            yield np.array([1, 2, 3]), None
            yield np.array([2, 3], dtype=np.float32), np.float32
            yield np.array([2, 3], dtype=np.int8), np.int8
            yield np.array([2, 3], dtype=np.int8), np.complex64
            yield np.array([2, 3], dtype=np.int8), np.complex128

        pyfunc = asfarray
        cfunc = jit(nopython=True)(pyfunc)

        for arr, dt in inputs():
            if dt is None:
                expected = pyfunc(arr)
                got = cfunc(arr)
            else:
                expected = pyfunc(arr, dtype=dt)
                got = cfunc(arr, dtype=dt)

            self.assertPreciseEqual(expected, got)
            self.assertTrue(np.issubdtype(got.dtype, np.inexact), got.dtype)

        # test default kwarg variant
        pyfunc = asfarray_default_kwarg
        cfunc = jit(nopython=True)(pyfunc)
        arr = np.array([1, 2, 3])
        expected = pyfunc(arr)
        got = cfunc(arr)
        self.assertPreciseEqual(expected, got)
        self.assertTrue(np.issubdtype(got.dtype, np.inexact), got.dtype)

    def test_repeat(self):
        # np.repeat(a, repeats)
        np_pyfunc = np_repeat
        np_nbfunc = njit(np_pyfunc)

        # a.repeat(repeats)
        array_pyfunc = array_repeat
        array_nbfunc = njit(array_pyfunc)

        for pyfunc, nbfunc in ((np_pyfunc, np_nbfunc),
                               (array_pyfunc, array_nbfunc)):

            def check(a, repeats):
                self.assertPreciseEqual(pyfunc(a, repeats), nbfunc(a, repeats))

            # test array arguments
            target_numpy_values = [
                np.ones(1),
                np.arange(1000),
                np.array([[0, 1], [2, 3]]),
                np.array([]),
                np.array([[], []]),
            ]

            target_numpy_types = [
                np.uint32,
                np.int32,
                np.uint64,
                np.int64,
                np.float32,
                np.float64,
                np.complex64,
                np.complex128,
            ]

            target_numpy_inputs = (np.array(a,dtype=t) for a,t in
                                   itertools.product(target_numpy_values,
                                                     target_numpy_types))

            target_non_numpy_inputs = [
                1,
                1.0,
                True,
                1j,
                [0, 1, 2],
                (0, 1, 2),
            ]
            for i in itertools.chain(target_numpy_inputs,
                                     target_non_numpy_inputs):
                check(i, repeats=0)
                check(i, repeats=1)
                check(i, repeats=2)
                check(i, repeats=3)
                check(i, repeats=100)

            # check broadcasting when repeats is an array/list
            one = np.arange(1)
            for i in ([0], [1], [2]):
                check(one, repeats=i)
                check(one, repeats=np.array(i))

            two = np.arange(2)
            for i in ([0, 0], [0, 1], [1, 0], [0, 1], [1, 2], [2, 1], [2, 2]):
                check(two, repeats=i)
                check(two, repeats=np.array(i))

            check(two, repeats=np.array([2, 2], dtype=np.int32))
            check(np.arange(10), repeats=np.arange(10))

    def test_repeat_exception(self):
        # np.repeat(a, repeats)
        np_pyfunc = np_repeat
        np_nbfunc = njit(np_pyfunc)

        # a.repeat(repeats)
        array_pyfunc = array_repeat
        array_nbfunc = njit(array_pyfunc)

        self.disable_leak_check()

        for pyfunc, nbfunc in ((np_pyfunc, np_nbfunc),
                               (array_pyfunc, array_nbfunc)):

            # negative repeat argument
            with self.assertRaises(ValueError) as e:
                nbfunc(np.ones(1), -1)
            self.assertIn("negative dimensions are not allowed",
                          str(e.exception))

            # float repeat argument has custom error message
            with self.assertRaises(TypingError) as e:
                nbfunc(np.ones(1), 1.0)
            self.assertIn(
                "The repeats argument must be an integer "
                "or an array-like of integer dtype",
                str(e.exception))

            # negative repeat argument as array
            with self.assertRaises(ValueError) as e:
                nbfunc(np.ones(2), np.array([1, -1]))
            self.assertIn("negative dimensions are not allowed",
                          str(e.exception))

            # broadcasting error, repeats too large
            with self.assertRaises(ValueError) as e:
                nbfunc(np.ones(2), np.array([1, 1, 1]))
            self.assertIn("operands could not be broadcast together",
                          str(e.exception))

            # broadcasting error, repeats too small
            with self.assertRaises(ValueError) as e:
                nbfunc(np.ones(5), np.array([1, 1, 1, 1]))
            self.assertIn("operands could not be broadcast together",
                          str(e.exception))

            # float repeat argument has custom error message
            with self.assertRaises(TypingError) as e:
                nbfunc(np.ones(2), [1.0, 1.0])
            self.assertIn(
                "The repeats argument must be an integer "
                "or an array-like of integer dtype",
                str(e.exception))

            for rep in [True, "a", "1"]:
                with self.assertRaises(TypingError):
                    nbfunc(np.ones(1), rep)

    def test_select(self):
        np_pyfunc = np_select
        np_nbfunc = njit(np_select)

        test_cases = [
            # Each test case below is one tuple.
            # Each tuple is separated by a description of what's being tested

            # test with arrays of length 3 instead of 2 and a different default
            ([np.array([False, False, False]),
              np.array([False, True, False]),
              np.array([False, False, True])],
             [np.array([1, 2, 3]),
              np.array([4, 5, 6]),
              np.array([7, 8, 9])], 15.3),
            # test with arrays of length 1 instead of 2
            ([np.array([True]),
              np.array([False])], [np.array([1]), np.array([2])], 0),
            # test with lists of length 100 of arrays of length 1
            ([np.array([False])] * 100, [np.array([1])] * 100, 0),
            # passing arrays with NaNs
            ([np.isnan(np.array([1, 2, 3, np.nan, 5, 7]))] * 2,
             [np.array([1, 2, 3, np.nan, 5, 7])] * 2, 0),
            # passing lists with 2d arrays
            ([np.isnan(np.array([[1, 2, 3, np.nan, 5, 7]]))] * 2,
             [np.array([[1, 2, 3, np.nan, 5, 7]])] * 2, 0),
            # passing arrays with complex numbers
            ([np.isnan(np.array([1, 2, 3 + 2j, np.nan, 5, 7]))] * 2,
             [np.array([1, 2, 3 + 2j, np.nan, 5, 7])] * 2, 0)
        ]

        for x in (np.arange(10), np.arange(10).reshape((5, 2))):
            # test with two lists
            test_cases.append(([x < 3, x > 5], [x, x ** 2], 0))
            # test with two tuples
            test_cases.append(((x < 3, x > 5), (x, x ** 2), 0))
            # test with one list and one tuple
            test_cases.append(([x < 3, x > 5], (x, x ** 2), 0))
            # test with one tuple and one list
            test_cases.append(((x < 3, x > 5), [x, x ** 2], 0))

        for condlist,  choicelist, default in test_cases:
            self.assertPreciseEqual(np_pyfunc(condlist, choicelist, default),
                                    np_nbfunc(condlist, choicelist, default))

        np_pyfunc_defaults = np_select_defaults
        np_nbfunc_defaults = njit(np_select_defaults)
        # check the defaults work, using whatever the last input was
        self.assertPreciseEqual(np_pyfunc_defaults(condlist, choicelist),
                                np_nbfunc_defaults(condlist, choicelist))

    def test_select_exception(self):
        np_nbfunc = njit(np_select)
        x = np.arange(10)
        self.disable_leak_check()
        for condlist, choicelist, default, expected_error, expected_text in [
            # Each test case below is one tuple.
            # Each tuple is separated by the description of the intended error

            # passing condlist of dim zero
            ([np.array(True), np.array([False, True, False])],
             [np.array(1), np.arange(12).reshape(4, 3)], 0,
             TypingError, "condlist arrays must be of at least dimension 1"),
            # condlist and choicelist with different dimensions
            ([np.array(True), np.array(False)], [np.array([1]), np.array([2])],
             0, TypingError, "condlist and choicelist elements must have the "
                             "same number of dimensions"),
            # condlist and choicelist with different dimensions
            ([np.array([True]), np.array([False])],
             [np.array([[1]]), np.array([[2]])], 0, TypingError,
             "condlist and choicelist elements must have the "
             "same number of dimensions"),
            # passing choicelist of dim zero
            ([np.array(True), np.array(False)], [np.array(1), np.array(2)], 0,
             TypingError, "condlist arrays must be of at least dimension 1"),
            # passing an array as condlist instead of a list or tuple
            (np.isnan(np.array([1, 2, 3, np.nan, 5, 7])),
             np.array([1, 2, 3, np.nan, 5, 7]), 0, TypingError,
             "condlist must be a List or a Tuple"),
            # default is a list
            ([True], [0], [0], TypingError,
             "default must be a scalar"),
            # condlist with ints instead of booleans
            ([(x < 3).astype(int), (x > 5).astype(int)], [x, x ** 2], 0,
             TypingError, "condlist arrays must contain booleans"),
            # condlist and choicelist of different length
            ([x > 9, x > 8, x > 7, x > 6], [x, x**2, x], 0, ValueError,
             "list of cases must be same length as list of conditions"),

            # condlist contains tuples instead of arrays
            # if in the future numba's np.where accepts tuples, the
            # implementation of np.select should also accept them and
            # the following two test cases should be normal tests
            # instead of negative tests

            # test with lists of length 100 of tuples of length 1 for condlist
            ([(False,)] * 100, [np.array([1])] * 100, 0, TypingError,
             'items of condlist must be arrays'),
            # test with lists of length 100 of tuples of length 1 for choicelist
            ([np.array([False])] * 100, [(1,)] * 100, 0, TypingError,
             'items of choicelist must be arrays'),
        ]:
            with self.assertRaises(expected_error) as e:
                np_nbfunc(condlist, choicelist, default)
            self.assertIn(expected_text, str(e.exception))

    def test_windowing(self):
        def check_window(func):
            np_pyfunc = func
            np_nbfunc = njit(func)

            for M in [0, 1, 5, 12]:
                expected = np_pyfunc(M)
                got = np_nbfunc(M)
                self.assertPreciseEqual(expected, got, prec='double')

            for M in ['a', 1.1, 1j]:
                with self.assertRaises(TypingError) as raises:
                    np_nbfunc(1.1)
                self.assertIn("M must be an integer", str(raises.exception))

        check_window(np_bartlett)
        check_window(np_blackman)
        check_window(np_hamming)
        check_window(np_hanning)

        # Test np.kaiser separately
        np_pyfunc = np_kaiser
        np_nbfunc = njit(np_kaiser)

        for M in [0, 1, 5, 12]:
            for beta in [0.0, 5.0, 14.0]:
                expected = np_pyfunc(M, beta)
                got = np_nbfunc(M, beta)

                if IS_32BITS or platform.machine() in ['ppc64le', 'aarch64']:
                    self.assertPreciseEqual(expected,
                                            got, prec='double', ulps=2)
                else:
                    self.assertPreciseEqual(expected, got, prec='double',
                                            ulps=2)

        for M in ['a', 1.1, 1j]:
            with self.assertRaises(TypingError) as raises:
                np_nbfunc(M, 1.0)
            self.assertIn("M must be an integer", str(raises.exception))

        for beta in ['a', 1j]:
            with self.assertRaises(TypingError) as raises:
                np_nbfunc(5, beta)
            self.assertIn("beta must be an integer or float",
                          str(raises.exception))

    def test_cross(self):
        pyfunc = np_cross
        cfunc = jit(nopython=True)(pyfunc)
        pairs = [
            # 3x3 (n-dims)
            (
                np.array([[1, 2, 3], [4, 5, 6]]),
                np.array([[4, 5, 6], [1, 2, 3]])
            ),
            # 2x3 array-like (n-dims)
            (
                np.array([[1, 2, 3], [4, 5, 6]]),
                ((4, 5), (1, 2))
            ),
            # 3x3 (1-dim) with type promotion
            (
                np.array([1, 2, 3], dtype=np.int64),
                np.array([4, 5, 6], dtype=np.float64)
            ),
            # 3x3 array-like (1-dim)
            (
                (1, 2, 3),
                (4, 5, 6)
            ),
            # 2x3 (1-dim)
            (
                np.array([1, 2]),
                np.array([4, 5, 6])
            ),
            # 3x3 (with broadcasting 1d x 2d)
            (
                np.array([1, 2, 3]),
                np.array([[4, 5, 6], [1, 2, 3]])
            ),
            # 3x3 (with broadcasting 2d x 1d)
            (
                np.array([[1, 2, 3], [4, 5, 6]]),
                np.array([1, 2, 3])
            ),
            # 3x2 (with higher order broadcasting)
            (
                np.arange(36).reshape(6, 2, 3),
                np.arange(4).reshape(2, 2)
            )
        ]

        for x, y in pairs:
            expected = pyfunc(x, y)
            got = cfunc(x, y)
            self.assertPreciseEqual(expected, got)

    def test_cross_exceptions(self):
        pyfunc = np_cross
        cfunc = jit(nopython=True)(pyfunc)
        self.disable_leak_check()

        # test incompatible dimensions for ndim == 1
        with self.assertRaises(ValueError) as raises:
            cfunc(
                np.arange(4),
                np.arange(3)
            )
        self.assertIn(
            'Incompatible dimensions for cross product',
            str(raises.exception)
        )

        # test 2d cross product error for ndim == 1
        with self.assertRaises(ValueError) as raises:
            cfunc(
                np.array((1, 2)),
                np.array((3, 4))
            )
        self.assertIn(
            'Dimensions for both inputs is 2.',
            str(raises.exception)
        )

        self.assertIn(
            '`cross2d(a, b)` from `numba.np.extensions`.',
            str(raises.exception)
        )

        # test incompatible dimensions for ndim > 1
        with self.assertRaises(ValueError) as raises:
            cfunc(
                np.arange(8).reshape((2, 4)),
                np.arange(6)[::-1].reshape((2, 3))
            )
        self.assertIn(
            'Incompatible dimensions for cross product',
            str(raises.exception)
        )

        # test 2d cross product error for ndim == 1
        with self.assertRaises(ValueError) as raises:
            cfunc(
                np.arange(8).reshape((4, 2)),
                np.arange(8)[::-1].reshape((4, 2))
            )
        self.assertIn(
            'Dimensions for both inputs is 2',
            str(raises.exception)
        )

        # test non-array-like input
        with self.assertRaises(TypingError) as raises:
            cfunc(
                set([1, 2, 3]),
                set([4, 5, 6])
            )
        self.assertIn(
            'Inputs must be array-like.',
            str(raises.exception)
        )

    def test_cross2d(self):
        pyfunc = np_cross
        cfunc = njit(nb_cross2d)
        pairs = [
            # 2x2 (n-dims)
            (
                np.array([[1, 2], [4, 5]]),
                np.array([[4, 5], [1, 2]])
            ),
            # 2x2 array-like (n-dims)
            (
                np.array([[1, 2], [4, 5]]),
                ((4, 5), (1, 2))
            ),
            # 2x2 (1-dim) with type promotion
            (
                np.array([1, 2], dtype=np.int64),
                np.array([4, 5], dtype=np.float64)
            ),
            # 2x2 array-like (1-dim)
            (
                (1, 2),
                (4, 5)
            ),
            # 2x2 (with broadcasting 1d x 2d)
            (
                np.array([1, 2]),
                np.array([[4, 5], [1, 2]])
            ),
            # 2x2 (with broadcasting 2d x 1d)
            (
                np.array([[1, 2], [4, 5]]),
                np.array([1, 2])
            ),
            # 2x2 (with higher order broadcasting)
            (
                np.arange(36).reshape(6, 3, 2),
                np.arange(6).reshape(3, 2)
            )
        ]

        for x, y in pairs:
            expected = pyfunc(x, y)
            got = cfunc(x, y)
            self.assertPreciseEqual(expected, got)

    def test_cross2d_exceptions(self):
        cfunc = njit(nb_cross2d)
        self.disable_leak_check()

        # test incompatible dimensions for ndim == 1
        with self.assertRaises(ValueError) as raises:
            cfunc(
                np.array((1, 2, 3)),
                np.array((4, 5, 6))
            )
        self.assertIn(
            'Incompatible dimensions for 2D cross product',
            str(raises.exception)
        )

        # test incompatible dimensions for ndim > 1
        with self.assertRaises(ValueError) as raises:
            cfunc(
                np.arange(6).reshape((2, 3)),
                np.arange(6)[::-1].reshape((2, 3))
            )
        self.assertIn(
            'Incompatible dimensions for 2D cross product',
            str(raises.exception)
        )

        # test non-array-like input
        with self.assertRaises(TypingError) as raises:
            cfunc(
                set([1, 2]),
                set([4, 5])
            )
        self.assertIn(
            'Inputs must be array-like.',
            str(raises.exception)
        )

    def test_trim_zeros(self):

        def arrays():
            yield np.array([])
            yield np.zeros(5)
            yield np.zeros(1)
            yield np.array([1, 2, 3])
            yield np.array([0, 1, 2, 3])
            yield np.array([0., 1., 2., np.nan, 0.])
            yield np.array(['0', 'Hello', 'world'])

        def explicit_trim():
            yield np.array([0, 1, 2, 0, 0]), 'FB'
            yield np.array([0, 1, 2]), 'B'
            yield np.array([np.nan, 0., 1.2, 2.3, 0.]), 'b'
            yield np.array([0, 0, 1, 2, 5]), 'f'
            if numpy_version < (2, 2):
                # abf and d are not supported in numpy >= 2.2
                yield np.array([0, 1, 2, 0]), 'abf'
                yield np.array([0, 4, 0]), 'd'

            yield np.array(['\0', '1', '2']), 'f'

        pyfunc = np_trim_zeros
        cfunc = jit(nopython=True)(pyfunc)

        for arr in arrays():
            expected = pyfunc(arr)
            got = cfunc(arr)
            self.assertPreciseEqual(expected, got)

        for arr, trim in explicit_trim():
            expected = pyfunc(arr, trim)
            got = cfunc(arr, trim)
            self.assertPreciseEqual(expected, got)

    def test_trim_zeros_numpy(self):
        # https://github.com/numpy/numpy/blob/9d8d46ad615a7e13256b930146ac369f651016c0/numpy/lib/tests/test_function_base.py#L1251-L1313
        a = np.array([0, 0, 1, 0, 2, 3, 4, 0])
        b = a.astype(float)
        c = a.astype(complex)
        # d = a.astype(object)
        values = [a, b, c]

        # test_basic
        slc = np.s_[2:-1]
        for arr in values:
            res = np_trim_zeros(arr)
            self.assertPreciseEqual(res, arr[slc])

        # test_leading_skip
        slc = np.s_[:-1]
        for arr in values:
            res = np_trim_zeros(arr, trim='b')
            self.assertPreciseEqual(res, arr[slc])

        # test_trailing_skip
        slc = np.s_[2:]
        for arr in values:
            res = np_trim_zeros(arr, trim='F')
            self.assertPreciseEqual(res, arr[slc])

        # test_all_zero
        for _arr in values:
            arr = np.zeros_like(_arr, dtype=_arr.dtype)

            res1 = np_trim_zeros(arr, trim='B')
            assert len(res1) == 0

            res2 = np_trim_zeros(arr, trim='f')
            assert len(res2) == 0

        # test_size_zero
        arr = np.zeros(0)
        res = np_trim_zeros(arr)
        self.assertPreciseEqual(arr, res)

        # test_overflow
        for arr in [np.array([0, 2**62, 0]), np.array([0, 2**63, 0]),
                    np.array([0, 2**64, 0])]:
            slc = np.s_[1:2]
            res = np_trim_zeros(arr)
            self.assertPreciseEqual(res, arr[slc])

        # test_no_trim
        arr = np.array([None, 1, None])
        res = np_trim_zeros(arr)
        self.assertPreciseEqual(arr, res)

        # test_list_to_list
        res = np_trim_zeros(a.tolist())
        assert isinstance(res, list)

    def test_trim_zeros_exceptions(self):
        self.disable_leak_check()
        cfunc = jit(nopython=True)(np_trim_zeros)

        with self.assertRaises(TypingError) as raises:
            cfunc(np.array([[1, 2, 3], [4, 5, 6]]))
        self.assertIn(
            'array must be 1D',
            str(raises.exception)
        )

        with self.assertRaises(TypingError) as raises:
            cfunc(3)
        self.assertIn(
            'The first argument must be an array',
            str(raises.exception)
        )

        with self.assertRaises(TypingError) as raises:
            cfunc({0, 1, 2})
        self.assertIn(
            'The first argument must be an array',
            str(raises.exception)
        )

        with self.assertRaises(TypingError) as raises:
            cfunc(np.array([0, 1, 2]), 1)
        self.assertIn(
            'The second argument must be a string',
            str(raises.exception)
        )

    def test_union1d(self):
        pyfunc = np_union1d
        cfunc = jit(nopython=True)(pyfunc)
        arrays = [
            # Test 1d arrays
            (
                np.array([1, 2, 3]),
                np.array([2, 3, 4])
            ),
            # Test 2d with 1d array
            (
                np.array([[1, 2, 3], [2, 3, 4]]),
                np.array([2, 5, 6])
            ),
            # Test 3d with 1d array
            (
                np.arange(0, 20).reshape(2,2,5),
                np.array([1, 20, 21])
            ),
            # Test 2d with 3d array
            (
                np.arange(0, 10).reshape(2,5),
                np.arange(0, 20).reshape(2,5,2)
            ),
            # Test other array-like
            (
                np.array([False, True, 7]),
                np.array([1, 2, 3])
            )
        ]

        for a, b in arrays:
            expected = pyfunc(a,b)
            got = cfunc(a,b)
            self.assertPreciseEqual(expected, got)

    def test_union1d_exceptions(self):
        cfunc = jit(nopython=True)(np_union1d)
        self.disable_leak_check()

        # Test inputs not array-like
        with self.assertRaises(TypingError) as raises:
            cfunc("Hello", np.array([1,2]))
        self.assertIn(
            "The arguments to np.union1d must be array-like",
            str(raises.exception)
        )
        with self.assertRaises(TypingError) as raises:
            cfunc(np.array([1,2]), "Hello")
        self.assertIn(
            "The arguments to np.union1d must be array-like",
            str(raises.exception)
        )
        with self.assertRaises(TypingError) as raises:
            cfunc("Hello", "World")
        self.assertIn(
            "The arguments to np.union1d must be array-like",
            str(raises.exception)
        )

        # Test Unicode array exceptions
        with self.assertRaises(TypingError) as raises:
            cfunc(np.array(['hello', 'world']), np.array(['a', 'b']))
        self.assertIn(
            "For Unicode arrays, arrays must have same dtype",
            str(raises.exception)
        )
        with self.assertRaises(TypingError) as raises:
            cfunc(np.array(['c', 'd']), np.array(['foo', 'bar']))
        self.assertIn(
            "For Unicode arrays, arrays must have same dtype",
            str(raises.exception)
        )
        with self.assertRaises(TypingError) as raises:
            cfunc(np.array(['c', 'd']), np.array([1, 2]))
        self.assertIn(
            "For Unicode arrays, arrays must have same dtype",
            str(raises.exception)
        )
        with self.assertRaises(TypingError) as raises:
            cfunc(np.array(['c', 'd']), np.array([1.1, 2.5]))
        self.assertIn(
            "For Unicode arrays, arrays must have same dtype",
            str(raises.exception)
        )

    def test_asarray_chkfinite(self):
        pyfunc = np_asarray_chkfinite
        cfunc = jit(nopython=True)(pyfunc)
        self.disable_leak_check()

        pairs = [
            #1D array with all args
            (
                np.array([1, 2, 3]),
                np.float32,
            ),
            #1D array
            (
                np.array([1, 2, 3]),
            ),
            #1D array-like
            (
                [1, 2, 3, 4],
            ),
            # 2x2 (n-dims)
            (
                np.array([[1, 2], [3, 4]]),
                np.float32,
            ),
            # 2x2 array-like (n-dims)
            (
                ((1, 2), (3, 4)),
                np.int64
            ),
            # 2x2 (1-dim) with type promotion
            (
                np.array([1, 2], dtype=np.int64),
            ),
            # 3x2 (with higher order broadcasting)
            (
                np.arange(36).reshape(6, 2, 3),
            ),
        ]

        for pair in pairs:
            expected = pyfunc(*pair)
            got = cfunc(*pair)
            self.assertPreciseEqual(expected, got)

    def test_asarray_chkfinite_exceptions(self):
        cfunc = jit(nopython=True)(np_asarray_chkfinite)
        self.disable_leak_check()

        #test for single value
        with self.assertRaises(TypingError) as e:
            cfunc(2)
        msg = "The argument to np.asarray_chkfinite must be array-like"
        self.assertIn(msg, str(e.exception))

        #test for NaNs
        with self.assertRaises(ValueError) as e:
            cfunc(np.array([2, 4, np.nan, 5]))
        self.assertIn("array must not contain infs or NaNs", str(e.exception))

        #test for infs
        with self.assertRaises(ValueError) as e:
            cfunc(np.array([1, 2, np.inf, 4]))
        self.assertIn("array must not contain infs or NaNs", str(e.exception))

        #test for dtype
        with self.assertRaises(TypingError) as e:
            cfunc(np.array([1, 2, 3, 4]), 'float32')
        self.assertIn("dtype must be a valid Numpy dtype", str(e.exception))

    def test_unwrap_basic(self):
        pyfunc = unwrap
        cfunc = njit(pyfunc)

        pyfunc1 = unwrap1
        cfunc1 = njit(pyfunc1)

        pyfunc13 = unwrap13
        cfunc13 = njit(pyfunc13)

        pyfunc123 = unwrap123
        cfunc123 = njit(pyfunc123)
        # Based on tests from https://github.com/numpy/numpy/blob/3032e84ff34f20def2ef4ebf9f8695947af3fd24/numpy/lib/tests/test_function_base.py#L1979-L2003 # noqa: E501
        # Additional tests are included to ensure proper support for
        # higher dimensional arrays

        # p only
        def inputs1():
            yield np.array([1, 1 + 2 * np.pi])
            phase = np.linspace(0, np.pi, num=5)
            phase[3:] += np.pi
            yield phase
            yield np.arange(16).reshape((4,4))
            yield np.arange(160, step=10).reshape((4,4))
            yield np.arange(240, step=10).reshape((2,3,4))

        for p in inputs1():
            self.assertPreciseEqual(pyfunc1(p), cfunc1(p))

        uneven_seq = np.array([0, 75, 150, 225, 300, 430])
        wrap_uneven = np.mod(uneven_seq, 250)

        # p and period only
        def inputs13():
            yield np.array([1, 1 + 256]), 255
            yield np.array([0, 75, 150, 225, 300]), 255
            yield np.array([0, 1, 2, -1, 0]), 4
            yield np.array([2, 3, 4, 5, 2, 3, 4, 5]), 4
            yield wrap_uneven, 250

        # check that you can set axis=-1 without errors
        self.assertPreciseEqual(pyfunc(wrap_uneven, axis=-1, period=250),
                                cfunc(wrap_uneven, axis=-1, period=250))

        for p, period in inputs13():
            self.assertPreciseEqual(pyfunc13(p, period=period),
                                    cfunc13(p, period=period))

        # p, period and discont
        def inputs123():
            yield wrap_uneven, 250, 140

        for p, period, discont in inputs123():
            self.assertPreciseEqual(pyfunc123(p, period=period,
                                              discont=discont),
                                    cfunc123(p, period=period,
                                             discont=discont))

    def test_unwrap_exception(self):
        cfunc = njit(unwrap)
        self.disable_leak_check()

        with self.assertRaises(TypingError) as e:
            cfunc('abc')
        self.assertIn('The argument "p" must be array-like',
                      str(e.exception))

        with self.assertRaises(TypingError) as e:
            cfunc(np.array([1, 2]), 'abc')
        self.assertIn('The argument "discont" must be a scalar',
                      str(e.exception))

        with self.assertRaises(TypingError) as e:
            cfunc(np.array([1, 2]), 3, period='abc')
        self.assertIn('The argument "period" must be a scalar',
                      str(e.exception))

        with self.assertRaises(TypingError) as e:
            cfunc(np.array([1, 2]), 3, axis='abc')
        self.assertIn('The argument "axis" must be an integer',
                      str(e.exception))

        with self.assertRaises(ValueError) as e:
            cfunc(np.array([1, 2]), 3, axis=2)
        self.assertIn('Value for argument "axis" is not supported',
                      str(e.exception))

    def test_swapaxes_basic(self):
        pyfunc = swapaxes
        cfunc = jit(nopython=True)(pyfunc)

        def a_variations():
            yield np.arange(10)
            yield np.arange(10).reshape(2, 5)
            yield np.arange(60).reshape(5, 4, 3)

        for a in a_variations():
            for a1 in range(-a.ndim, a.ndim):
                for a2 in range(-a.ndim, a.ndim):
                    expected = pyfunc(a, a1, a2)
                    got = cfunc(a, a1, a2)
                    self.assertPreciseEqual(expected, got)

    def test_swapaxes_exception(self):
        pyfunc = swapaxes
        cfunc = jit(nopython=True)(pyfunc)

        # Exceptions leak references
        self.disable_leak_check()

        with self.assertRaises(TypingError) as raises:
            cfunc('abc', 0, 0)

        self.assertIn('The first argument "a" must be an array',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc(np.arange(4), 'abc', 0)

        self.assertIn('The second argument "axis1" must be an integer',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc(np.arange(4), 0, 'abc')

        self.assertIn('The third argument "axis2" must be an integer',
                      str(raises.exception))

        with self.assertRaises(ValueError) as raises:
            cfunc(np.arange(4), 1, 0)

        self.assertIn('np.swapaxes: Argument axis1 out of bounds',
                      str(raises.exception))

        with self.assertRaises(ValueError) as raises:
            cfunc(np.arange(8).reshape(2, 4), 0, -3)

        self.assertIn('np.swapaxes: Argument axis2 out of bounds',
                      str(raises.exception))

    def test_take_along_axis(self):
        a = np.arange(24).reshape((3, 1, 4, 2))

        # For now axis must be literal, test explicitly defined implementations
        @njit
        def axis_none(a, i):
            return np.take_along_axis(a, i, axis=None)

        indices = np.array([1, 2], dtype=np.uint64)
        self.assertPreciseEqual(axis_none(a, indices),
                                axis_none.py_func(a, indices))

        def gen(axis):
            @njit
            def impl(a, i):
                return np.take_along_axis(a, i, axis)
            return impl

        for i in range(-1, a.ndim):
            jfunc = gen(i)
            ai = np.argsort(a, axis=i)
            self.assertPreciseEqual(jfunc(a, ai), jfunc.py_func(a, ai))

    def test_take_along_axis_broadcasting(self):
        # Based on
        # https://github.com/numpy/numpy/blob/v1.21.0/numpy/lib/tests/test_shape_base.py#L74-L79
        # This demonstrates that arrays are broadcast before the algorithm is
        # applied.
        arr = np.ones((3, 4, 1))
        ai = np.ones((1, 2, 5), dtype=np.intp)

        def gen(axis):
            @njit
            def impl(a, i):
                return np.take_along_axis(a, i, axis)
            return impl

        # Check same axis but expressed as positive/negative value
        for i in (1, -2):
            check = gen(i)
            expected = check.py_func(arr, ai)
            actual = check(arr, ai)
            self.assertPreciseEqual(expected, actual)
            self.assertEqual(actual.shape, (3, 2, 5))

    def test_take_along_axis_exceptions(self):
        arr2d = np.arange(8).reshape(2, 4)
        # Valid indices when axis=None is passed to take_along_axis:
        indices_none = np.array([0, 1], dtype=np.uint64)
        indices = np.ones((2, 4), dtype=np.uint64)

        # For now axis must be literal, so we need to construct functions with
        # explicit axis:
        def gen(axis):
            @njit
            def impl(a, i):
                return np.take_along_axis(a, i, axis)
            return impl

        with self.assertRaises(TypingError) as raises:
            gen("a")(arr2d, indices)
        self.assertIn("axis must be an integer", str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            gen(-3)(arr2d, indices)
        self.assertIn("axis is out of bounds", str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            gen(2)(arr2d, indices)
        self.assertIn("axis is out of bounds", str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            gen(None)(12, indices_none)
        self.assertIn('"arr" must be an array', str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            gen(None)(arr2d, 5)
        self.assertIn('"indices" must be an array', str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            gen(None)(arr2d, np.array([0.0, 1.0]))
        self.assertIn(
            'indices array must contain integers',
            str(raises.exception)
        )

        @njit
        def not_literal_axis(a, i, axis):
            return np.take_along_axis(a, i, axis)

        with self.assertRaises(TypingError) as raises:
            not_literal_axis(arr2d, indices, 0)
        self.assertIn("axis must be a literal value", str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            gen(0)(arr2d, np.array([0, 1], dtype=np.uint64))
        self.assertIn("must have the same number of dimensions",
                      str(raises.exception))

        # With axis None, array's ndim is implicitly 1.
        with self.assertRaises(TypingError) as raises:
            gen(None)(arr2d, arr2d)
        self.assertIn("must have the same number of dimensions",
                      str(raises.exception))

        with self.assertRaises(ValueError) as raises:
            gen(0)(arr2d, np.ones((2, 3), dtype=np.uint64))
        self.assertIn("dimensions don't match", str(raises.exception))
        # Exceptions leak references
        self.disable_leak_check()

    def test_nan_to_num(self):
        # Test cases are from
        # https://github.com/numpy/numpy/blob/8ff45c5bb520db04af8720bf1d34a392a8d2561a/numpy/lib/tests/test_type_check.py#L350-L452
        values = [
            np.nan,
            1,
            1.1,
            1 + 1j,
            complex(-np.inf, np.nan),
            complex(np.nan, np.nan),
            np.array([1], dtype=int),
            np.array([complex(-np.inf, np.inf), complex(1, np.nan),
                      complex(np.nan, 1), complex(np.inf, -np.inf)]),
            np.array([0.1, 1.0, 0.4]),
            np.array([1, 2, 3]),
            np.array([[0.1, 1.0, 0.4], [0.4, 1.2, 4.0]]),
            np.array([0.1, np.nan, 0.4]),
            np.array([[0.1, np.nan, 0.4], [np.nan, 1.2, 4.0]]),
            np.array([-np.inf, np.nan, np.inf]),
            np.array([-np.inf, np.nan, np.inf], dtype=np.float32)
        ]
        nans = [0.0, 10]

        pyfunc = nan_to_num
        cfunc = njit(nan_to_num)

        for value, nan in product(values, nans):
            expected = pyfunc(value, nan=nan)
            got = cfunc(value, nan=nan)
            self.assertPreciseEqual(expected, got)

    def test_nan_to_num_copy_false(self):
        # Check that copy=False operates in-place.
        cfunc = njit(nan_to_num)

        x = np.array([0.1, 0.4, np.nan])
        expected = 1.0
        cfunc(x, copy=False, nan=expected)
        self.assertPreciseEqual(x[-1], expected)

        x_complex = np.array([0.1, 0.4, complex(np.nan, np.nan)])
        cfunc(x_complex, copy=False, nan=expected)
        self.assertPreciseEqual(x_complex[-1], 1. + 1.j)

    def test_nan_to_num_invalid_argument(self):
        cfunc = njit(nan_to_num)

        with self.assertTypingError() as raises:
            cfunc("invalid_input")
        self.assertIn("The first argument must be a scalar or an array-like",
                      str(raises.exception))

    def test_diagflat_basic(self):
        pyfunc1 = diagflat1
        cfunc1 = njit(pyfunc1)
        pyfunc2 = diagflat2
        cfunc2 = njit(pyfunc2)

        def inputs():
            yield np.array([1,2]), 1
            yield np.array([[1,2],[3,4]]), -2
            yield np.arange(8).reshape((2,2,2)), 2
            yield [1, 2], 1
            yield np.array([]), 1

        for v, k in inputs():
            self.assertPreciseEqual(pyfunc1(v), cfunc1(v))
            self.assertPreciseEqual(pyfunc2(v, k), cfunc2(v, k))

    def test_diagflat1_exception(self):
        pyfunc = diagflat1
        cfunc = njit(pyfunc)

        self.disable_leak_check()

        with self.assertRaises(TypingError) as raises:
            cfunc("abc")
        self.assertIn('The argument "v" must be array-like',
                      str(raises.exception))

    def test_diagflat2_exception(self):
        pyfunc = diagflat2
        cfunc = njit(pyfunc)

        self.disable_leak_check()

        with self.assertRaises(TypingError) as raises:
            cfunc("abc", 2)
        self.assertIn('The argument "v" must be array-like',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc([1, 2], "abc")
        self.assertIn('The argument "k" must be an integer',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc([1, 2], 3.0)
        self.assertIn('The argument "k" must be an integer',
                      str(raises.exception))

    @staticmethod
    def _setxor_arrays():
        yield (List.empty_list(types.float64),
               List.empty_list(types.float64))  # two empty arrays
        yield [1], List.empty_list(types.float64)  # empty right
        yield List.empty_list(types.float64), [1]  # empty left
        yield [1], [2]  # singletons - xor == union
        yield [1], [1]  # singletons - xor == nothing
        yield [1, 2], [1]
        yield [1, 2, 2], [2, 2]
        yield [1, 2, 2], [2, 2, 3]
        yield [1, 2], [2, 1]
        yield [1, 2, 3], [1, 2, 3]
        yield [2, 3, 4, 0], [1, 3]
        # from numpy:
        # https://github.com/numpy/numpy/blob/b0371ef240560e78b651a5d7c9407ae3212a3d56/numpy/lib/tests/test_arraysetops.py#L86 # noqa: E501
        yield [5, 7, 1, 2], [2, 4, 3, 1, 5]
        yield [1, 2, 3], [6, 5, 4]
        yield [1, 8, 2, 3], [1, 2, 3, 4, 5, 6]

    def test_setxor1d_2(self):
        np_pyfunc = np_setxor1d_2
        np_nbfunc = njit(np_pyfunc)

        def check(ar1, ar2):
            if isinstance(ar1, list):
                ar1 = List(ar1)
            if isinstance(ar2, list):
                ar2 = List(ar2)
            expected = np_pyfunc(ar1, ar2)
            got = np_nbfunc(ar1, ar2)
            self.assertPreciseEqual(expected, got, msg=f"ar1={ar1}, ar2={ar2}")

        for a, b in self._setxor_arrays():
            check(a, b)

    def test_setxor1d_3(self):
        np_pyfunc = np_setxor1d_3
        np_nbfunc = njit(np_pyfunc)

        def check(ar1, ar2, assume_unique=False):
            if isinstance(ar1, list):
                ar1 = List(ar1)
            if isinstance(ar2, list):
                ar2 = List(ar2)
            expected = np_pyfunc(ar1, ar2, assume_unique)
            got = np_nbfunc(ar1, ar2, assume_unique)
            self.assertPreciseEqual(expected, got, msg=f"ar1={ar1}, ar2={ar2}")

        for a, b in self._setxor_arrays():
            check(a, b)
            if len(np.unique(a)) == len(a) and len(np.unique(b)) == len(b):
                check(a, b, True)

    def test_setxor1d_errors(self):
        np_pyfunc = np_setxor1d_3
        np_nbfunc = njit(np_pyfunc)

        a = np.array([1])
        b = np.array([2])
        self.disable_leak_check()
        with self.assertRaises(TypingError):
            np_nbfunc(a, b, "foo")
        with self.assertRaises(TypingError):
            np_nbfunc("foo", b, True)
        with self.assertRaises(TypingError):
            np_nbfunc(a, "foo", True)

    @staticmethod
    def _setdiff_arrays():
        yield (List.empty_list(types.float64),
               List.empty_list(types.float64))  # two empty arrays
        yield [1], List.empty_list(types.float64)  # empty right
        yield List.empty_list(types.float64), [1]  # empty left
        yield [1], [2]  # singletons - diff == [1]
        yield [1], [1]  # singletons - diff == nothing
        yield [1, 2], [1]
        yield [1, 2, 2], [2, 2]
        yield [1, 2, 2], [2, 2, 3]
        yield [1, 2], [2, 1]
        yield [1, 2, 3], [1, 2, 3]
        yield [2, 3, 4, 0], [1, 3]

        # https://github.com/numpy/numpy/blob/b0371ef240560e78b651a5d7c9407ae3212a3d56/numpy/lib/tests/test_arraysetops.py#L558 # noqa: E501
        yield (np.array([6, 5, 4, 7, 1, 2, 7, 4]),
               np.array([2, 4, 3, 3, 2, 1, 5]))
        yield np.arange(21), np.arange(19)
        yield np.array([3, 2, 1]), np.array([7, 5, 2])

    def test_setdiff1d_2(self):
        np_pyfunc = np_setdiff1d_2
        np_nbfunc = njit(np_pyfunc)

        def check(ar1, ar2):
            if isinstance(ar1, list):
                ar1 = List(ar1)
            if isinstance(ar2, list):
                ar2 = List(ar2)
            expected = np_pyfunc(ar1, ar2)
            got = np_nbfunc(ar1, ar2)
            self.assertPreciseEqual(expected, got, msg=f"ar1={ar1}, ar2={ar2}")

        for a, b in self._setdiff_arrays():
            check(a, b)

    def test_setdiff1d_3(self):
        np_pyfunc = np_setdiff1d_3
        np_nbfunc = njit(np_pyfunc)

        def check(ar1, ar2, assume_unique=False):
            if isinstance(ar1, list):
                ar1 = List(ar1)
            if isinstance(ar2, list):
                ar2 = List(ar2)
            expected = np_pyfunc(ar1, ar2, assume_unique)
            got = np_nbfunc(ar1, ar2, assume_unique)
            self.assertPreciseEqual(expected, got, msg=f"ar1={ar1}, ar2={ar2}")

        for a, b in self._setdiff_arrays():
            check(a, b)
            if len(np.unique(a)) == len(a) and len(np.unique(b)) == len(b):
                check(a, b, True)

    def test_setdiff1d_errors(self):
        np_pyfunc = np_setdiff1d_3
        np_nbfunc = njit(np_pyfunc)

        a = np.array([1])
        b = np.array([2])
        self.disable_leak_check()
        with self.assertRaises(TypingError):
            np_nbfunc(a, b, "foo")
        with self.assertRaises(TypingError):
            np_nbfunc("foo", b, True)
        with self.assertRaises(TypingError):
            np_nbfunc(a, "foo", True)

    @staticmethod
    def _in1d_arrays():
        yield (List.empty_list(types.float64),
               List.empty_list(types.float64))  # two empty arrays
        yield [1], List.empty_list(types.float64)  # empty right
        yield List.empty_list(types.float64), [1]  # empty left
        yield [1], [2]  # singletons - False
        yield [1], [1]  # singletons - True
        yield [1, 2], [1]
        yield [1, 2, 2], [2, 2]
        yield [1, 2, 2], [2, 2, 3]
        yield [1, 2], [2, 1]
        yield [1, 2, 3], [1, 2, 3]
        yield [2, 3, 4, 0], [3, 1]
        yield [2, 3], np.arange(20)  # Test the "sorting" method.
        yield [2, 3], np.tile(np.arange(5), 4)

    def test_in1d_2(self):
        np_pyfunc = np_in1d_2
        np_nbfunc = njit(np_pyfunc)

        def check(ar1, ar2):
            if isinstance(ar1, list):
                ar1 = List(ar1)
            if isinstance(ar2, list):
                ar2 = List(ar2)
            expected = np_pyfunc(ar1, ar2)
            got = np_nbfunc(ar1, ar2)
            self.assertPreciseEqual(expected, got, msg=f"ar1={ar1}, ar2={ar2}")

        for a, b in self._in1d_arrays():
            check(a, b)

    def test_in1d_3a(self):
        np_pyfunc = np_in1d_3a
        np_nbfunc = njit(np_pyfunc)

        def check(ar1, ar2, assume_unique=False):
            if isinstance(ar1, list):
                ar1 = List(ar1)
            if isinstance(ar2, list):
                ar2 = List(ar2)
            expected = np_pyfunc(ar1, ar2, assume_unique)
            got = np_nbfunc(ar1, ar2, assume_unique)
            self.assertPreciseEqual(expected, got, msg=f"ar1={ar1}, ar2={ar2}")

        for a, b in self._in1d_arrays():
            check(a, b)
            if len(np.unique(a)) == len(a) and len(np.unique(b)) == len(b):
                check(a, b, assume_unique=True)

    def test_in1d_3b(self):
        np_pyfunc = np_in1d_3b
        np_nbfunc = njit(np_pyfunc)

        def check(ar1, ar2, invert=False):
            if isinstance(ar1, list):
                ar1 = List(ar1)
            if isinstance(ar2, list):
                ar2 = List(ar2)
            expected = np_pyfunc(ar1, ar2, invert)
            got = np_nbfunc(ar1, ar2, invert)
            self.assertPreciseEqual(expected, got, msg=f"ar1={ar1}, ar2={ar2}")

        for a, b in self._in1d_arrays():
            check(a, b, invert=False)
            check(a, b, invert=True)

    def test_in1d_4(self):
        np_pyfunc = np_in1d_4
        np_nbfunc = njit(np_pyfunc)

        def check(ar1, ar2, assume_unique=False, invert=False):
            if isinstance(ar1, list):
                ar1 = List(ar1)
            if isinstance(ar2, list):
                ar2 = List(ar2)
            expected = np_pyfunc(ar1, ar2, assume_unique, invert)
            got = np_nbfunc(ar1, ar2, assume_unique, invert)
            self.assertPreciseEqual(expected, got, msg=f"ar1={ar1}, ar2={ar2}")

        for a, b in self._in1d_arrays():
            check(a, b, invert=False)
            check(a, b, invert=True)
            if len(np.unique(a)) == len(a) and len(np.unique(b)) == len(b):
                check(a, b, assume_unique=True, invert=False)
                check(a, b, assume_unique=True, invert=True)

    def test_in1d_errors(self):
        np_pyfunc = np_in1d_4
        np_nbfunc = njit(np_pyfunc)

        a = np.array([1])
        b = np.array([2])
        x = np_nbfunc(a, b)
        self.assertPreciseEqual(x, np.array([False]))

        self.disable_leak_check()
        with self.assertRaises(TypingError):
            np_nbfunc(a, b, "foo", False)
        with self.assertRaises(TypingError):
            np_nbfunc(a, b, False, "foo")
        with self.assertRaises(TypingError):
            np_nbfunc("foo", b, True, False)
        with self.assertRaises(TypingError):
            np_nbfunc(a, "foo", True, False)

        @njit()
        def np_in1d_kind(a, b, kind):
            return np.in1d(a, b, kind=kind)

        with self.assertRaises(TypingError):
            np_in1d_kind(a, b, kind=None)
        with self.assertRaises(TypingError):
            np_in1d_kind(a, b, kind="table")

    @staticmethod
    def _isin_arrays():
        yield (List.empty_list(types.float64),
               List.empty_list(types.float64))  # two empty arrays
        yield (np.zeros((1, 0), dtype=np.int64),
               List.empty_list(types.int64))  # two-dim array - shape (1, 0)
        yield (np.zeros((0, 0), dtype=np.int64),
               List.empty_list(types.int64))
        yield (np.zeros((0, 1), dtype=np.int64),
               List.empty_list(types.int64))
        yield [1], List.empty_list(types.float64)  # empty right
        yield List.empty_list(types.float64), [1]  # empty left
        yield [1], [2]  # singletons - False
        yield [1], [1]  # singletons - True
        yield [1, 2], [1]
        yield [1, 2, 2], [2, 2]
        yield [1, 2, 2], [2, 2, 3]
        yield [1, 2], [2, 1]
        yield [2, 3], np.arange(20)  # Test the "sorting" method.
        yield [2, 3], np.tile(np.arange(5), 4)
        yield np.arange(30).reshape(2, 3, 5), [5, 7, 10, 15]  # 3d

        # from numpy
        # https://github.com/numpy/numpy/blob/b0371ef240560e78b651a5d7c9407ae3212a3d56/numpy/lib/tests/test_arraysetops.py#L200 # noqa: E501
        a = np.arange(24).reshape([2, 3, 4])
        b = np.array([[10, 20, 30], [0, 1, 3], [11, 22, 33]])
        yield a, b
        yield np.array(3), b
        yield a, np.array(3)
        yield np.array(3), np.array(3)
        yield 5, b
        yield a, 6
        yield 5, 6
        yield List.empty_list(types.int64), b
        yield a, List.empty_list(types.int64)

        for dtype in [bool, np.int64, np.float64]:
            if dtype in {np.int64, np.float64}:
                ar = np.array([10, 20, 30], dtype=dtype)
            elif dtype in {bool}:
                ar = np.array([True, False, False])

            empty_array = np.array([], dtype=dtype)

            yield empty_array, ar
            yield ar, empty_array
            yield empty_array, empty_array

        for mult in (1, 10):
            yield [5, 7, 1, 2], [2, 4, 3, 1, 5] * mult
            yield [8, 7, 1, 2], [2, 4, 3, 1, 5] * mult
            yield [4, 7, 1, 8], [2, 4, 3, 1, 5] * mult
            a = [5, 4, 5, 3, 4, 4, 3, 4, 3, 5, 2, 1, 5, 5]
            yield a, [2, 3, 4] * mult
            yield a, [2, 3, 4] * mult + [5, 5, 4] * mult
            yield np.array([5, 7, 1, 2]), np.array([2, 4, 3, 1, 5] * mult)
            yield np.array([5, 7, 1, 1, 2]), np.array([2, 4, 3, 3, 1, 5] * mult)
            yield np.array([5, 5]), np.array([2, 2] * mult)

        yield np.array([5]), np.array([2])
        yield np.array([True, False]), np.array([False, False, False])

        for dtype1, dtype2 in [
            (np.int8, np.int16),
            (np.int16, np.int8),
            (np.uint8, np.uint16),
            (np.uint16, np.uint8),
            (np.uint8, np.int16),
            (np.int16, np.uint8),
        ]:
            is_dtype2_signed = np.issubdtype(dtype2, np.signedinteger)
            ar1 = np.array([0, 0, 1, 1], dtype=dtype1)

            if is_dtype2_signed:
                ar2 = np.array([-128, 0, 127], dtype=dtype2)
            else:
                ar2 = np.array([127, 0, 255], dtype=dtype2)

            yield ar1, ar2

        for dtype in np.typecodes["AllInteger"]:
            a = np.array([True, False, False], dtype=bool)
            b = np.array([0, 0, 0, 0], dtype=dtype)
            yield a, b
            yield b, a

    def test_isin_2(self):
        np_pyfunc = np_isin_2
        np_nbfunc = njit(np_pyfunc)

        def check(ar1, ar2):
            expected = np_pyfunc(ar1, ar2)
            if isinstance(ar1, list):
                ar1 = List(ar1)
            if isinstance(ar2, list):
                ar2 = List(ar2)
            got = np_nbfunc(ar1, ar2)
            self.assertPreciseEqual(expected, got, msg=f"ar1={ar1}, ar2={ar2}")

        for a, b in self._isin_arrays():
            check(a, b)

    def test_isin_3a(self):
        np_pyfunc = np_isin_3a
        np_nbfunc = njit(np_pyfunc)

        def check(ar1, ar2, assume_unique=False):
            expected = np_pyfunc(ar1, ar2, assume_unique)
            if isinstance(ar1, list):
                ar1 = List(ar1)
            if isinstance(ar2, list):
                ar2 = List(ar2)
            got = np_nbfunc(ar1, ar2, assume_unique)
            self.assertPreciseEqual(expected, got, msg=f"ar1={ar1}, ar2={ar2}")

        for a, b in self._isin_arrays():
            check(a, b)

            try:
                len_a = len(a)
            except TypeError:
                len_a = 1
            try:
                len_b = len(b)
            except TypeError:
                len_b = 1
            if len(np.unique(a)) == len_a and len(np.unique(b)) == len_b:
                check(a, b, assume_unique=True)

    def test_isin_3b(self):
        np_pyfunc = np_isin_3b
        np_nbfunc = njit(np_pyfunc)

        def check(ar1, ar2, invert=False):
            expected = np_pyfunc(ar1, ar2, invert)
            if isinstance(ar1, list):
                ar1 = List(ar1)
            if isinstance(ar2, list):
                ar2 = List(ar2)
            got = np_nbfunc(ar1, ar2, invert)
            self.assertPreciseEqual(expected, got, msg=f"ar1={ar1}, ar2={ar2}")

        for a, b in self._isin_arrays():
            check(a, b, invert=False)
            check(a, b, invert=True)

    def test_isin_4(self):
        np_pyfunc = np_isin_4
        np_nbfunc = njit(np_pyfunc)

        def check(ar1, ar2, assume_unique=False, invert=False):
            expected = np_pyfunc(ar1, ar2, assume_unique, invert)
            if isinstance(ar1, list):
                ar1 = List(ar1)
            if isinstance(ar2, list):
                ar2 = List(ar2)
            got = np_nbfunc(ar1, ar2, assume_unique, invert)
            self.assertPreciseEqual(expected, got, msg=f"ar1={ar1}, ar2={ar2}")

        for a, b in self._isin_arrays():
            check(a, b, invert=False)
            check(a, b, invert=True)

            try:
                len_a = len(a)
            except TypeError:
                len_a = 1
            try:
                len_b = len(b)
            except TypeError:
                len_b = 1
            if len(np.unique(a)) == len_a and len(np.unique(b)) == len_b:
                check(a, b, assume_unique=True, invert=False)
                check(a, b, assume_unique=True, invert=True)

    def test_isin_errors(self):
        np_pyfunc = np_isin_4
        np_nbfunc = njit(np_pyfunc)

        a = np.array([1])
        b = np.array([2])
        x = np_nbfunc(a, b)
        self.assertPreciseEqual(x, np.array([False]))

        self.disable_leak_check()
        with self.assertRaises(TypingError):
            np_nbfunc(a, b, "foo", False)
        with self.assertRaises(TypingError):
            np_nbfunc(a, b, False, "foo")
        with self.assertRaises(TypingError):
            np_nbfunc("foo", b, True, False)
        with self.assertRaises(TypingError):
            np_nbfunc(a, "foo", True, False)

        @njit()
        def np_isin_kind(a, b, kind):
            return np.isin(a, b, kind=kind)

        with self.assertRaises(TypingError):
            np_isin_kind(a, b,  kind=None)
        with self.assertRaises(TypingError):
            np_isin_kind(a, b, kind="table")

    def test_setops_manyways(self):
        # https://github.com/numpy/numpy/blob/b0371ef240560e78b651a5d7c9407ae3212a3d56/numpy/lib/tests/test_arraysetops.py#L588 # noqa: E501
        nb_setxor1d = njit(np_setxor1d_2)
        nb_intersect1d = njit(intersect1d_2)
        nb_union1d = njit(np_union1d)
        nb_setdiff1d = njit(np_setdiff1d_2)

        a = np.array([5, 7, 1, 2, 8])
        b = np.array([9, 8, 2, 4, 3, 1, 5])

        c1 = nb_setxor1d(a, b)
        aux1 = nb_intersect1d(a, b)
        aux2 = nb_union1d(a, b)
        c2 = nb_setdiff1d(aux2, aux1)
        self.assertPreciseEqual(c1, c2)


class TestNPMachineParameters(TestCase):
    # tests np.finfo, np.iinfo, np.MachAr

    template = '''
def foo():
    ty = np.%s
    return np.%s(ty)
'''

    def check(self, func, attrs, *args):
        pyfunc = func
        cfunc = jit(nopython=True)(pyfunc)

        expected = pyfunc(*args)
        got = cfunc(*args)

        # check result
        for attr in attrs:
            self.assertPreciseEqual(getattr(expected, attr),
                                    getattr(got, attr))

    def create_harcoded_variant(self, basefunc, ty):
        #create an instance of using the function with a hardcoded type
        #and eval it into existence, return the function for use
        tystr = ty.__name__
        basestr = basefunc.__name__
        funcstr = self.template % (tystr, basestr)
        dct = {}
        exec(compile(funcstr, '<string>', 'exec'), globals(), dct)
        return dct['foo']

    def test_finfo(self):
        types = [np.float32, np.float64, np.complex64, np.complex128]
        attrs = ('eps', 'epsneg', 'iexp', 'machep', 'max', 'maxexp', 'negep',
                 'nexp', 'nmant', 'precision', 'resolution', 'tiny', 'bits',)
        for ty in types:
            self.check(finfo, attrs, ty(1))
            hc_func = self.create_harcoded_variant(np.finfo, ty)
            self.check(hc_func, attrs)

        # check unsupported attr raises
        with self.assertRaises(TypingError) as raises:
            cfunc = jit(nopython=True)(finfo_machar)
            cfunc(7.)
        msg = "Unknown attribute 'machar' of type finfo"
        self.assertIn(msg, str(raises.exception))

        # check invalid type raises
        with self.assertTypingError():
            cfunc = jit(nopython=True)(finfo)
            cfunc(np.int32(7))

    def test_iinfo(self):
        # check types and instances of types
        types = [np.int8, np.int16, np.int32, np.int64, np.uint8, np.uint16,
                 np.uint32, np.uint64]
        attrs = ('min', 'max', 'bits',)
        for ty in types:
            self.check(iinfo, attrs, ty(1))
            hc_func = self.create_harcoded_variant(np.iinfo, ty)
            self.check(hc_func, attrs)

        # check invalid type raises
        with self.assertTypingError():
            cfunc = jit(nopython=True)(iinfo)
            cfunc(np.float64(7))


class TestRegistryImports(TestCase):

    def test_unsafe_import_in_registry(self):
        # See 8940
        # This should not fail
        code = dedent("""
            import numba
            import numpy as np
            @numba.njit
            def foo():
                np.array([1 for _ in range(1)])
            foo()
            print("OK")
        """)
        result, error = run_in_subprocess(code)
        # Assert that the bytestring "OK" was printed to stdout
        self.assertEqual(b"OK", result.strip())
        self.assertEqual(b"", error.strip(), msg=f"--ERROR--\n{error}\n")


if __name__ == '__main__':
    unittest.main()
