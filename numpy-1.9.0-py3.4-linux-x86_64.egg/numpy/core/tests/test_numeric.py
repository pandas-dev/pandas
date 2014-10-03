from __future__ import division, absolute_import, print_function

import sys
import platform
from decimal import Decimal
import warnings
import itertools
import platform

import numpy as np
from numpy.core import *
from numpy.core import umath
from numpy.random import rand, randint, randn
from numpy.testing import *
from numpy.core.multiarray import dot as dot_


class Vec(object):
    def __init__(self,sequence=None):
        if sequence is None:
            sequence=[]
        self.array=array(sequence)
    def __add__(self, other):
        out=Vec()
        out.array=self.array+other.array
        return out
    def __sub__(self, other):
        out=Vec()
        out.array=self.array-other.array
        return out
    def __mul__(self, other): # with scalar
        out=Vec(self.array.copy())
        out.array*=other
        return out
    def __rmul__(self, other):
        return self*other


class TestDot(TestCase):
    def setUp(self):
        self.A = rand(10, 8)
        self.b1 = rand(8, 1)
        self.b2 = rand(8)
        self.b3 = rand(1, 8)
        self.b4 = rand(10)
        self.N = 14

    def test_matmat(self):
        A = self.A
        c1 = dot(A.transpose(), A)
        c2 = dot_(A.transpose(), A)
        assert_almost_equal(c1, c2, decimal=self.N)

    def test_matvec(self):
        A, b1 = self.A, self.b1
        c1 = dot(A, b1)
        c2 = dot_(A, b1)
        assert_almost_equal(c1, c2, decimal=self.N)

    def test_matvec2(self):
        A, b2 = self.A, self.b2
        c1 = dot(A, b2)
        c2 = dot_(A, b2)
        assert_almost_equal(c1, c2, decimal=self.N)

    def test_vecmat(self):
        A, b4 = self.A, self.b4
        c1 = dot(b4, A)
        c2 = dot_(b4, A)
        assert_almost_equal(c1, c2, decimal=self.N)

    def test_vecmat2(self):
        b3, A = self.b3, self.A
        c1 = dot(b3, A.transpose())
        c2 = dot_(b3, A.transpose())
        assert_almost_equal(c1, c2, decimal=self.N)

    def test_vecmat3(self):
        A, b4 = self.A, self.b4
        c1 = dot(A.transpose(), b4)
        c2 = dot_(A.transpose(), b4)
        assert_almost_equal(c1, c2, decimal=self.N)

    def test_vecvecouter(self):
        b1, b3 = self.b1, self.b3
        c1 = dot(b1, b3)
        c2 = dot_(b1, b3)
        assert_almost_equal(c1, c2, decimal=self.N)

    def test_vecvecinner(self):
        b1, b3 = self.b1, self.b3
        c1 = dot(b3, b1)
        c2 = dot_(b3, b1)
        assert_almost_equal(c1, c2, decimal=self.N)

    def test_columnvect1(self):
        b1 = ones((3, 1))
        b2 = [5.3]
        c1 = dot(b1, b2)
        c2 = dot_(b1, b2)
        assert_almost_equal(c1, c2, decimal=self.N)

    def test_columnvect2(self):
        b1 = ones((3, 1)).transpose()
        b2 = [6.2]
        c1 = dot(b2, b1)
        c2 = dot_(b2, b1)
        assert_almost_equal(c1, c2, decimal=self.N)

    def test_vecscalar(self):
        b1 = rand(1, 1)
        b2 = rand(1, 8)
        c1 = dot(b1, b2)
        c2 = dot_(b1, b2)
        assert_almost_equal(c1, c2, decimal=self.N)

    def test_vecscalar2(self):
        b1 = rand(8, 1)
        b2 = rand(1, 1)
        c1 = dot(b1, b2)
        c2 = dot_(b1, b2)
        assert_almost_equal(c1, c2, decimal=self.N)

    def test_all(self):
        dims = [(), (1,), (1, 1)]
        for dim1 in dims:
            for dim2 in dims:
                arg1 = rand(*dim1)
                arg2 = rand(*dim2)
                c1 = dot(arg1, arg2)
                c2 = dot_(arg1, arg2)
                assert_(c1.shape == c2.shape)
                assert_almost_equal(c1, c2, decimal=self.N)

    def test_vecobject(self):
        U_non_cont = transpose([[1., 1.], [1., 2.]])
        U_cont = ascontiguousarray(U_non_cont)
        x = array([Vec([1., 0.]), Vec([0., 1.])])
        zeros = array([Vec([0., 0.]), Vec([0., 0.])])
        zeros_test = dot(U_cont, x) - dot(U_non_cont, x)
        assert_equal(zeros[0].array, zeros_test[0].array)
        assert_equal(zeros[1].array, zeros_test[1].array)


class TestResize(TestCase):
    def test_copies(self):
        A = array([[1, 2], [3, 4]])
        Ar1 = array([[1, 2, 3, 4], [1, 2, 3, 4]])
        assert_equal(resize(A, (2, 4)), Ar1)

        Ar2 = array([[1, 2], [3, 4], [1, 2], [3, 4]])
        assert_equal(resize(A, (4, 2)), Ar2)

        Ar3 = array([[1, 2, 3], [4, 1, 2], [3, 4, 1], [2, 3, 4]])
        assert_equal(resize(A, (4, 3)), Ar3)

    def test_zeroresize(self):
        A = array([[1, 2], [3, 4]])
        Ar = resize(A, (0,))
        assert_equal(Ar, array([]))

class TestNonarrayArgs(TestCase):
    # check that non-array arguments to functions wrap them in arrays
    def test_squeeze(self):
        A = [[[1, 1, 1], [2, 2, 2], [3, 3, 3]]]
        assert_(squeeze(A).shape == (3, 3))

    def test_cumproduct(self):
        A = [[1, 2, 3], [4, 5, 6]]
        assert_(all(cumproduct(A) == array([1, 2, 6, 24, 120, 720])))

    def test_size(self):
        A = [[1, 2, 3], [4, 5, 6]]
        assert_(size(A) == 6)
        assert_(size(A, 0) == 2)
        assert_(size(A, 1) == 3)

    def test_mean(self):
        A = [[1, 2, 3], [4, 5, 6]]
        assert_(mean(A) == 3.5)
        assert_(all(mean(A, 0) == array([2.5, 3.5, 4.5])))
        assert_(all(mean(A, 1) == array([2., 5.])))

        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', '', RuntimeWarning)
            assert_(isnan(mean([])))
            assert_(w[0].category is RuntimeWarning)

    def test_std(self):
        A = [[1, 2, 3], [4, 5, 6]]
        assert_almost_equal(std(A), 1.707825127659933)
        assert_almost_equal(std(A, 0), array([1.5, 1.5, 1.5]))
        assert_almost_equal(std(A, 1), array([0.81649658, 0.81649658]))

        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', '', RuntimeWarning)
            assert_(isnan(std([])))
            assert_(w[0].category is RuntimeWarning)

    def test_var(self):
        A = [[1, 2, 3], [4, 5, 6]]
        assert_almost_equal(var(A), 2.9166666666666665)
        assert_almost_equal(var(A, 0), array([2.25, 2.25, 2.25]))
        assert_almost_equal(var(A, 1), array([0.66666667, 0.66666667]))

        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', '', RuntimeWarning)
            assert_(isnan(var([])))
            assert_(w[0].category is RuntimeWarning)


class TestBoolScalar(TestCase):
    def test_logical(self):
        f = False_
        t = True_
        s = "xyz"
        self.assertTrue((t and s) is s)
        self.assertTrue((f and s) is f)

    def test_bitwise_or(self):
        f = False_
        t = True_
        self.assertTrue((t | t) is t)
        self.assertTrue((f | t) is t)
        self.assertTrue((t | f) is t)
        self.assertTrue((f | f) is f)

    def test_bitwise_and(self):
        f = False_
        t = True_
        self.assertTrue((t & t) is t)
        self.assertTrue((f & t) is f)
        self.assertTrue((t & f) is f)
        self.assertTrue((f & f) is f)

    def test_bitwise_xor(self):
        f = False_
        t = True_
        self.assertTrue((t ^ t) is f)
        self.assertTrue((f ^ t) is t)
        self.assertTrue((t ^ f) is t)
        self.assertTrue((f ^ f) is f)


class TestBoolArray(TestCase):
    def setUp(self):
        # offset for simd tests
        self.t = array([True] * 41, dtype=np.bool)[1::]
        self.f = array([False] * 41, dtype=np.bool)[1::]
        self.o = array([False] * 42, dtype=np.bool)[2::]
        self.nm = self.f.copy()
        self.im = self.t.copy()
        self.nm[3] = True
        self.nm[-2] = True
        self.im[3] = False
        self.im[-2] = False

    def test_all_any(self):
        self.assertTrue(self.t.all())
        self.assertTrue(self.t.any())
        self.assertFalse(self.f.all())
        self.assertFalse(self.f.any())
        self.assertTrue(self.nm.any())
        self.assertTrue(self.im.any())
        self.assertFalse(self.nm.all())
        self.assertFalse(self.im.all())
        # check bad element in all positions
        for i in range(256 - 7):
            d = array([False] * 256, dtype=np.bool)[7::]
            d[i] = True
            self.assertTrue(np.any(d))
            e = array([True] * 256, dtype=np.bool)[7::]
            e[i] = False
            self.assertFalse(np.all(e))
            assert_array_equal(e, ~d)
        # big array test for blocked libc loops
        for i in list(range(9, 6000, 507)) + [7764, 90021, -10]:
            d = array([False] * 100043, dtype=np.bool)
            d[i] = True
            self.assertTrue(np.any(d), msg="%r" % i)
            e = array([True] * 100043, dtype=np.bool)
            e[i] = False
            self.assertFalse(np.all(e), msg="%r" % i)

    def test_logical_not_abs(self):
        assert_array_equal(~self.t, self.f)
        assert_array_equal(np.abs(~self.t), self.f)
        assert_array_equal(np.abs(~self.f), self.t)
        assert_array_equal(np.abs(self.f), self.f)
        assert_array_equal(~np.abs(self.f), self.t)
        assert_array_equal(~np.abs(self.t), self.f)
        assert_array_equal(np.abs(~self.nm), self.im)
        np.logical_not(self.t, out=self.o)
        assert_array_equal(self.o, self.f)
        np.abs(self.t, out=self.o)
        assert_array_equal(self.o, self.t)

    def test_logical_and_or_xor(self):
        assert_array_equal(self.t | self.t, self.t)
        assert_array_equal(self.f | self.f, self.f)
        assert_array_equal(self.t | self.f, self.t)
        assert_array_equal(self.f | self.t, self.t)
        np.logical_or(self.t, self.t, out=self.o)
        assert_array_equal(self.o, self.t)
        assert_array_equal(self.t & self.t, self.t)
        assert_array_equal(self.f & self.f, self.f)
        assert_array_equal(self.t & self.f, self.f)
        assert_array_equal(self.f & self.t, self.f)
        np.logical_and(self.t, self.t, out=self.o)
        assert_array_equal(self.o, self.t)
        assert_array_equal(self.t ^ self.t, self.f)
        assert_array_equal(self.f ^ self.f, self.f)
        assert_array_equal(self.t ^ self.f, self.t)
        assert_array_equal(self.f ^ self.t, self.t)
        np.logical_xor(self.t, self.t, out=self.o)
        assert_array_equal(self.o, self.f)

        assert_array_equal(self.nm & self.t, self.nm)
        assert_array_equal(self.im & self.f, False)
        assert_array_equal(self.nm & True, self.nm)
        assert_array_equal(self.im & False, self.f)
        assert_array_equal(self.nm | self.t, self.t)
        assert_array_equal(self.im | self.f, self.im)
        assert_array_equal(self.nm | True, self.t)
        assert_array_equal(self.im | False, self.im)
        assert_array_equal(self.nm ^ self.t, self.im)
        assert_array_equal(self.im ^ self.f, self.im)
        assert_array_equal(self.nm ^ True, self.im)
        assert_array_equal(self.im ^ False, self.im)


class TestBoolCmp(TestCase):
    def setUp(self):
        self.f = ones(256, dtype=np.float32)
        self.ef = ones(self.f.size, dtype=np.bool)
        self.d = ones(128, dtype=np.float64)
        self.ed = ones(self.d.size, dtype=np.bool)
        # generate values for all permutation of 256bit simd vectors
        s = 0
        for i in range(32):
            self.f[s:s+8] = [i & 2**x for x in range(8)]
            self.ef[s:s+8] = [(i & 2**x) != 0 for x in range(8)]
            s += 8
        s = 0
        for i in range(16):
            self.d[s:s+4] = [i & 2**x for x in range(4)]
            self.ed[s:s+4] = [(i & 2**x) != 0 for x in range(4)]
            s += 4

        self.nf = self.f.copy()
        self.nd = self.d.copy()
        self.nf[self.ef] = np.nan
        self.nd[self.ed] = np.nan

    def test_float(self):
        # offset for alignment test
        for i in range(4):
            assert_array_equal(self.f[i:] > 0, self.ef[i:])
            assert_array_equal(self.f[i:] - 1 >= 0, self.ef[i:])
            assert_array_equal(self.f[i:] == 0, ~self.ef[i:])
            assert_array_equal(-self.f[i:] < 0, self.ef[i:])
            assert_array_equal(-self.f[i:] + 1 <= 0, self.ef[i:])
            r = self.f[i:] != 0
            assert_array_equal(r, self.ef[i:])
            r2 = self.f[i:] != np.zeros_like(self.f[i:])
            r3 = 0 != self.f[i:]
            assert_array_equal(r, r2)
            assert_array_equal(r, r3)
            # check bool == 0x1
            assert_array_equal(r.view(np.int8), r.astype(np.int8))
            assert_array_equal(r2.view(np.int8), r2.astype(np.int8))
            assert_array_equal(r3.view(np.int8), r3.astype(np.int8))

            # isnan on amd64 takes the same codepath
            assert_array_equal(np.isnan(self.nf[i:]), self.ef[i:])

    def test_double(self):
        # offset for alignment test
        for i in range(2):
            assert_array_equal(self.d[i:] > 0, self.ed[i:])
            assert_array_equal(self.d[i:] - 1 >= 0, self.ed[i:])
            assert_array_equal(self.d[i:] == 0, ~self.ed[i:])
            assert_array_equal(-self.d[i:] < 0, self.ed[i:])
            assert_array_equal(-self.d[i:] + 1 <= 0, self.ed[i:])
            r = self.d[i:] != 0
            assert_array_equal(r, self.ed[i:])
            r2 = self.d[i:] != np.zeros_like(self.d[i:])
            r3 = 0 != self.d[i:]
            assert_array_equal(r, r2)
            assert_array_equal(r, r3)
            # check bool == 0x1
            assert_array_equal(r.view(np.int8), r.astype(np.int8))
            assert_array_equal(r2.view(np.int8), r2.astype(np.int8))
            assert_array_equal(r3.view(np.int8), r3.astype(np.int8))

            # isnan on amd64 takes the same codepath
            assert_array_equal(np.isnan(self.nd[i:]), self.ed[i:])


class TestSeterr(TestCase):
    def test_default(self):
        err = geterr()
        self.assertEqual(err, dict(
            divide='warn',
            invalid='warn',
            over='warn',
            under='ignore',
        ))

    def test_set(self):
        with np.errstate():
            err = seterr()
            old = seterr(divide='print')
            self.assertTrue(err == old)
            new = seterr()
            self.assertTrue(new['divide'] == 'print')
            seterr(over='raise')
            self.assertTrue(geterr()['over'] == 'raise')
            self.assertTrue(new['divide'] == 'print')
            seterr(**old)
            self.assertTrue(geterr() == old)

    @dec.skipif(platform.machine() == "armv5tel", "See gh-413.")
    def test_divide_err(self):
        with errstate(divide='raise'):
            try:
                array([1.]) / array([0.])
            except FloatingPointError:
                pass
            else:
                self.fail()
            seterr(divide='ignore')
            array([1.]) / array([0.])

    def test_errobj(self):
        olderrobj = np.geterrobj()
        self.called = 0
        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                with errstate(divide='warn'):
                    np.seterrobj([20000, 1, None])
                    array([1.]) / array([0.])
                    self.assertEqual(len(w), 1)

            def log_err(*args):
                self.called += 1
                extobj_err = args
                assert (len(extobj_err) == 2)
                assert ("divide" in extobj_err[0])

            with errstate(divide='ignore'):
                np.seterrobj([20000, 3, log_err])
                array([1.]) / array([0.])
            self.assertEqual(self.called, 1)

            np.seterrobj(olderrobj)
            with errstate(divide='ignore'):
                np.divide(1., 0., extobj=[20000, 3, log_err])
            self.assertEqual(self.called, 2)
        finally:
            np.seterrobj(olderrobj)
            del self.called

    def test_errobj_noerrmask(self):
        # errmask = 0 has a special code path for the default
        olderrobj = np.geterrobj()
        try:
            # set errobj to something non default
            np.seterrobj([umath.UFUNC_BUFSIZE_DEFAULT,
                         umath.ERR_DEFAULT + 1, None])
            #call a ufunc
            np.isnan(np.array([6]))
            # same with the default, lots of times to get rid of possible
            # pre-existing stack in the code
            for i in range(10000):
                np.seterrobj([umath.UFUNC_BUFSIZE_DEFAULT, umath.ERR_DEFAULT,
                             None])
            np.isnan(np.array([6]))
        finally:
            np.seterrobj(olderrobj)


class TestFloatExceptions(TestCase):
    def assert_raises_fpe(self, fpeerr, flop, x, y):
        ftype = type(x)
        try:
            flop(x, y)
            assert_(False,
                    "Type %s did not raise fpe error '%s'." % (ftype, fpeerr))
        except FloatingPointError as exc:
            assert_(str(exc).find(fpeerr) >= 0,
                    "Type %s raised wrong fpe error '%s'." % (ftype, exc))

    def assert_op_raises_fpe(self, fpeerr, flop, sc1, sc2):
        # Check that fpe exception is raised.
        #
        # Given a floating operation `flop` and two scalar values, check that
        # the operation raises the floating point exception specified by
        #`fpeerr`. Tests all variants with 0-d array scalars as well.

        self.assert_raises_fpe(fpeerr, flop, sc1, sc2);
        self.assert_raises_fpe(fpeerr, flop, sc1[()], sc2);
        self.assert_raises_fpe(fpeerr, flop, sc1, sc2[()]);
        self.assert_raises_fpe(fpeerr, flop, sc1[()], sc2[()]);

    @dec.knownfailureif(True, "See ticket #2350")
    def test_floating_exceptions(self):
        # Test basic arithmetic function errors
        with np.errstate(all='raise'):
            # Test for all real and complex float types
            for typecode in np.typecodes['AllFloat']:
                ftype = np.obj2sctype(typecode)
                if np.dtype(ftype).kind == 'f':
                    # Get some extreme values for the type
                    fi = np.finfo(ftype)
                    ft_tiny = fi.tiny
                    ft_max = fi.max
                    ft_eps = fi.eps
                    underflow = 'underflow'
                    divbyzero = 'divide by zero'
                else:
                    # 'c', complex, corresponding real dtype
                    rtype = type(ftype(0).real)
                    fi = np.finfo(rtype)
                    ft_tiny = ftype(fi.tiny)
                    ft_max = ftype(fi.max)
                    ft_eps = ftype(fi.eps)
                    # The complex types raise different exceptions
                    underflow = ''
                    divbyzero = ''
                overflow = 'overflow'
                invalid = 'invalid'

                self.assert_raises_fpe(underflow,
                        lambda a, b:a/b, ft_tiny, ft_max)
                self.assert_raises_fpe(underflow,
                        lambda a, b:a*b, ft_tiny, ft_tiny)
                self.assert_raises_fpe(overflow,
                        lambda a, b:a*b, ft_max, ftype(2))
                self.assert_raises_fpe(overflow,
                        lambda a, b:a/b, ft_max, ftype(0.5))
                self.assert_raises_fpe(overflow,
                        lambda a, b:a+b, ft_max, ft_max*ft_eps)
                self.assert_raises_fpe(overflow,
                        lambda a, b:a-b, -ft_max, ft_max*ft_eps)
                self.assert_raises_fpe(overflow,
                        np.power, ftype(2), ftype(2**fi.nexp))
                self.assert_raises_fpe(divbyzero,
                        lambda a, b:a/b, ftype(1), ftype(0))
                self.assert_raises_fpe(invalid,
                        lambda a, b:a/b, ftype(np.inf), ftype(np.inf))
                self.assert_raises_fpe(invalid,
                        lambda a, b:a/b, ftype(0), ftype(0))
                self.assert_raises_fpe(invalid,
                        lambda a, b:a-b, ftype(np.inf), ftype(np.inf))
                self.assert_raises_fpe(invalid,
                        lambda a, b:a+b, ftype(np.inf), ftype(-np.inf))
                self.assert_raises_fpe(invalid,
                        lambda a, b:a*b, ftype(0), ftype(np.inf))

    def test_warnings(self):
        # test warning code path
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            with np.errstate(all="warn"):
                np.divide(1, 0.)
                self.assertEqual(len(w), 1)
                self.assertTrue("divide by zero" in str(w[0].message))
                np.array(1e300) * np.array(1e300)
                self.assertEqual(len(w), 2)
                self.assertTrue("overflow" in str(w[-1].message))
                np.array(np.inf) - np.array(np.inf)
                self.assertEqual(len(w), 3)
                self.assertTrue("invalid value" in str(w[-1].message))
                np.array(1e-300) * np.array(1e-300)
                self.assertEqual(len(w), 4)
                self.assertTrue("underflow" in str(w[-1].message))


class TestTypes(TestCase):
    def check_promotion_cases(self, promote_func):
        #Tests that the scalars get coerced correctly.
        b = np.bool_(0)
        i8, i16, i32, i64 = int8(0), int16(0), int32(0), int64(0)
        u8, u16, u32, u64 = uint8(0), uint16(0), uint32(0), uint64(0)
        f32, f64, fld = float32(0), float64(0), longdouble(0)
        c64, c128, cld = complex64(0), complex128(0), clongdouble(0)

        # coercion within the same kind
        assert_equal(promote_func(i8, i16), np.dtype(int16))
        assert_equal(promote_func(i32, i8), np.dtype(int32))
        assert_equal(promote_func(i16, i64), np.dtype(int64))
        assert_equal(promote_func(u8, u32), np.dtype(uint32))
        assert_equal(promote_func(f32, f64), np.dtype(float64))
        assert_equal(promote_func(fld, f32), np.dtype(longdouble))
        assert_equal(promote_func(f64, fld), np.dtype(longdouble))
        assert_equal(promote_func(c128, c64), np.dtype(complex128))
        assert_equal(promote_func(cld, c128), np.dtype(clongdouble))
        assert_equal(promote_func(c64, fld), np.dtype(clongdouble))

        # coercion between kinds
        assert_equal(promote_func(b, i32), np.dtype(int32))
        assert_equal(promote_func(b, u8), np.dtype(uint8))
        assert_equal(promote_func(i8, u8), np.dtype(int16))
        assert_equal(promote_func(u8, i32), np.dtype(int32))
        assert_equal(promote_func(i64, u32), np.dtype(int64))
        assert_equal(promote_func(u64, i32), np.dtype(float64))
        assert_equal(promote_func(i32, f32), np.dtype(float64))
        assert_equal(promote_func(i64, f32), np.dtype(float64))
        assert_equal(promote_func(f32, i16), np.dtype(float32))
        assert_equal(promote_func(f32, u32), np.dtype(float64))
        assert_equal(promote_func(f32, c64), np.dtype(complex64))
        assert_equal(promote_func(c128, f32), np.dtype(complex128))
        assert_equal(promote_func(cld, f64), np.dtype(clongdouble))

        # coercion between scalars and 1-D arrays
        assert_equal(promote_func(array([b]), i8), np.dtype(int8))
        assert_equal(promote_func(array([b]), u8), np.dtype(uint8))
        assert_equal(promote_func(array([b]), i32), np.dtype(int32))
        assert_equal(promote_func(array([b]), u32), np.dtype(uint32))
        assert_equal(promote_func(array([i8]), i64), np.dtype(int8))
        assert_equal(promote_func(u64, array([i32])), np.dtype(int32))
        assert_equal(promote_func(i64, array([u32])), np.dtype(uint32))
        assert_equal(promote_func(int32(-1), array([u64])), np.dtype(float64))
        assert_equal(promote_func(f64, array([f32])), np.dtype(float32))
        assert_equal(promote_func(fld, array([f32])), np.dtype(float32))
        assert_equal(promote_func(array([f64]), fld), np.dtype(float64))
        assert_equal(promote_func(fld, array([c64])), np.dtype(complex64))
        assert_equal(promote_func(c64, array([f64])), np.dtype(complex128))
        assert_equal(promote_func(complex64(3j), array([f64])),
                                                    np.dtype(complex128))

        # coercion between scalars and 1-D arrays, where
        # the scalar has greater kind than the array
        assert_equal(promote_func(array([b]), f64), np.dtype(float64))
        assert_equal(promote_func(array([b]), i64), np.dtype(int64))
        assert_equal(promote_func(array([b]), u64), np.dtype(uint64))
        assert_equal(promote_func(array([i8]), f64), np.dtype(float64))
        assert_equal(promote_func(array([u16]), f64), np.dtype(float64))

        # uint and int are treated as the same "kind" for
        # the purposes of array-scalar promotion.
        assert_equal(promote_func(array([u16]), i32), np.dtype(uint16))

        # float and complex are treated as the same "kind" for
        # the purposes of array-scalar promotion, so that you can do
        # (0j + float32array) to get a complex64 array instead of
        # a complex128 array.
        assert_equal(promote_func(array([f32]), c128), np.dtype(complex64))

    def test_coercion(self):
        def res_type(a, b):
            return np.add(a, b).dtype
        self.check_promotion_cases(res_type)

        # Use-case: float/complex scalar * bool/int8 array
        #           shouldn't narrow the float/complex type
        for a in [np.array([True, False]), np.array([-3, 12], dtype=np.int8)]:
            b = 1.234 * a
            assert_equal(b.dtype, np.dtype('f8'), "array type %s" % a.dtype)
            b = np.longdouble(1.234) * a
            assert_equal(b.dtype, np.dtype(np.longdouble),
                                                "array type %s" % a.dtype)
            b = np.float64(1.234) * a
            assert_equal(b.dtype, np.dtype('f8'), "array type %s" % a.dtype)
            b = np.float32(1.234) * a
            assert_equal(b.dtype, np.dtype('f4'), "array type %s" % a.dtype)
            b = np.float16(1.234) * a
            assert_equal(b.dtype, np.dtype('f2'), "array type %s" % a.dtype)

            b = 1.234j * a
            assert_equal(b.dtype, np.dtype('c16'), "array type %s" % a.dtype)
            b = np.clongdouble(1.234j) * a
            assert_equal(b.dtype, np.dtype(np.clongdouble),
                                                "array type %s" % a.dtype)
            b = np.complex128(1.234j) * a
            assert_equal(b.dtype, np.dtype('c16'), "array type %s" % a.dtype)
            b = np.complex64(1.234j) * a
            assert_equal(b.dtype, np.dtype('c8'), "array type %s" % a.dtype)

        # The following use-case is problematic, and to resolve its
        # tricky side-effects requires more changes.
        #
        ## Use-case: (1-t)*a, where 't' is a boolean array and 'a' is
        ##            a float32, shouldn't promote to float64
        #a = np.array([1.0, 1.5], dtype=np.float32)
        #t = np.array([True, False])
        #b = t*a
        #assert_equal(b, [1.0, 0.0])
        #assert_equal(b.dtype, np.dtype('f4'))
        #b = (1-t)*a
        #assert_equal(b, [0.0, 1.5])
        #assert_equal(b.dtype, np.dtype('f4'))
        ## Probably ~t (bitwise negation) is more proper to use here,
        ## but this is arguably less intuitive to understand at a glance, and
        ## would fail if 't' is actually an integer array instead of boolean:
        #b = (~t)*a
        #assert_equal(b, [0.0, 1.5])
        #assert_equal(b.dtype, np.dtype('f4'))

    def test_result_type(self):
        self.check_promotion_cases(np.result_type)
        assert_(np.result_type(None) == np.dtype(None))

    def test_promote_types_endian(self):
        # promote_types should always return native-endian types
        assert_equal(np.promote_types('<i8', '<i8'), np.dtype('i8'))
        assert_equal(np.promote_types('>i8', '>i8'), np.dtype('i8'))

        assert_equal(np.promote_types('>i8', '>U16'), np.dtype('U21'))
        assert_equal(np.promote_types('<i8', '<U16'), np.dtype('U21'))
        assert_equal(np.promote_types('>U16', '>i8'), np.dtype('U21'))
        assert_equal(np.promote_types('<U16', '<i8'), np.dtype('U21'))

        assert_equal(np.promote_types('<S5', '<U8'), np.dtype('U8'))
        assert_equal(np.promote_types('>S5', '>U8'), np.dtype('U8'))
        assert_equal(np.promote_types('<U8', '<S5'), np.dtype('U8'))
        assert_equal(np.promote_types('>U8', '>S5'), np.dtype('U8'))
        assert_equal(np.promote_types('<U5', '<U8'), np.dtype('U8'))
        assert_equal(np.promote_types('>U8', '>U5'), np.dtype('U8'))

        assert_equal(np.promote_types('<M8', '<M8'), np.dtype('M8'))
        assert_equal(np.promote_types('>M8', '>M8'), np.dtype('M8'))
        assert_equal(np.promote_types('<m8', '<m8'), np.dtype('m8'))
        assert_equal(np.promote_types('>m8', '>m8'), np.dtype('m8'))

    def test_promote_types_strings(self):
        assert_equal(np.promote_types('bool', 'S'), np.dtype('S5'))
        assert_equal(np.promote_types('b', 'S'), np.dtype('S4'))
        assert_equal(np.promote_types('u1', 'S'), np.dtype('S3'))
        assert_equal(np.promote_types('u2', 'S'), np.dtype('S5'))
        assert_equal(np.promote_types('u4', 'S'), np.dtype('S10'))
        assert_equal(np.promote_types('u8', 'S'), np.dtype('S20'))
        assert_equal(np.promote_types('i1', 'S'), np.dtype('S4'))
        assert_equal(np.promote_types('i2', 'S'), np.dtype('S6'))
        assert_equal(np.promote_types('i4', 'S'), np.dtype('S11'))
        assert_equal(np.promote_types('i8', 'S'), np.dtype('S21'))
        assert_equal(np.promote_types('bool', 'U'), np.dtype('U5'))
        assert_equal(np.promote_types('b', 'U'), np.dtype('U4'))
        assert_equal(np.promote_types('u1', 'U'), np.dtype('U3'))
        assert_equal(np.promote_types('u2', 'U'), np.dtype('U5'))
        assert_equal(np.promote_types('u4', 'U'), np.dtype('U10'))
        assert_equal(np.promote_types('u8', 'U'), np.dtype('U20'))
        assert_equal(np.promote_types('i1', 'U'), np.dtype('U4'))
        assert_equal(np.promote_types('i2', 'U'), np.dtype('U6'))
        assert_equal(np.promote_types('i4', 'U'), np.dtype('U11'))
        assert_equal(np.promote_types('i8', 'U'), np.dtype('U21'))
        assert_equal(np.promote_types('bool', 'S1'), np.dtype('S5'))
        assert_equal(np.promote_types('bool', 'S30'), np.dtype('S30'))
        assert_equal(np.promote_types('b', 'S1'), np.dtype('S4'))
        assert_equal(np.promote_types('b', 'S30'), np.dtype('S30'))
        assert_equal(np.promote_types('u1', 'S1'), np.dtype('S3'))
        assert_equal(np.promote_types('u1', 'S30'), np.dtype('S30'))
        assert_equal(np.promote_types('u2', 'S1'), np.dtype('S5'))
        assert_equal(np.promote_types('u2', 'S30'), np.dtype('S30'))
        assert_equal(np.promote_types('u4', 'S1'), np.dtype('S10'))
        assert_equal(np.promote_types('u4', 'S30'), np.dtype('S30'))
        assert_equal(np.promote_types('u8', 'S1'), np.dtype('S20'))
        assert_equal(np.promote_types('u8', 'S30'), np.dtype('S30'))

    def test_can_cast(self):
        assert_(np.can_cast(np.int32, np.int64))
        assert_(np.can_cast(np.float64, np.complex))
        assert_(not np.can_cast(np.complex, np.float))

        assert_(np.can_cast('i8', 'f8'))
        assert_(not np.can_cast('i8', 'f4'))
        assert_(np.can_cast('i4', 'S11'))

        assert_(np.can_cast('i8', 'i8', 'no'))
        assert_(not np.can_cast('<i8', '>i8', 'no'))

        assert_(np.can_cast('<i8', '>i8', 'equiv'))
        assert_(not np.can_cast('<i4', '>i8', 'equiv'))

        assert_(np.can_cast('<i4', '>i8', 'safe'))
        assert_(not np.can_cast('<i8', '>i4', 'safe'))

        assert_(np.can_cast('<i8', '>i4', 'same_kind'))
        assert_(not np.can_cast('<i8', '>u4', 'same_kind'))

        assert_(np.can_cast('<i8', '>u4', 'unsafe'))

        assert_(np.can_cast('bool', 'S5'))
        assert_(not np.can_cast('bool', 'S4'))

        assert_(np.can_cast('b', 'S4'))
        assert_(not np.can_cast('b', 'S3'))

        assert_(np.can_cast('u1', 'S3'))
        assert_(not np.can_cast('u1', 'S2'))
        assert_(np.can_cast('u2', 'S5'))
        assert_(not np.can_cast('u2', 'S4'))
        assert_(np.can_cast('u4', 'S10'))
        assert_(not np.can_cast('u4', 'S9'))
        assert_(np.can_cast('u8', 'S20'))
        assert_(not np.can_cast('u8', 'S19'))

        assert_(np.can_cast('i1', 'S4'))
        assert_(not np.can_cast('i1', 'S3'))
        assert_(np.can_cast('i2', 'S6'))
        assert_(not np.can_cast('i2', 'S5'))
        assert_(np.can_cast('i4', 'S11'))
        assert_(not np.can_cast('i4', 'S10'))
        assert_(np.can_cast('i8', 'S21'))
        assert_(not np.can_cast('i8', 'S20'))

        assert_(np.can_cast('bool', 'S5'))
        assert_(not np.can_cast('bool', 'S4'))

        assert_(np.can_cast('b', 'U4'))
        assert_(not np.can_cast('b', 'U3'))

        assert_(np.can_cast('u1', 'U3'))
        assert_(not np.can_cast('u1', 'U2'))
        assert_(np.can_cast('u2', 'U5'))
        assert_(not np.can_cast('u2', 'U4'))
        assert_(np.can_cast('u4', 'U10'))
        assert_(not np.can_cast('u4', 'U9'))
        assert_(np.can_cast('u8', 'U20'))
        assert_(not np.can_cast('u8', 'U19'))

        assert_(np.can_cast('i1', 'U4'))
        assert_(not np.can_cast('i1', 'U3'))
        assert_(np.can_cast('i2', 'U6'))
        assert_(not np.can_cast('i2', 'U5'))
        assert_(np.can_cast('i4', 'U11'))
        assert_(not np.can_cast('i4', 'U10'))
        assert_(np.can_cast('i8', 'U21'))
        assert_(not np.can_cast('i8', 'U20'))

        assert_raises(TypeError, np.can_cast, 'i4', None)
        assert_raises(TypeError, np.can_cast, None, 'i4')


# Custom exception class to test exception propagation in fromiter
class NIterError(Exception): pass


class TestFromiter(TestCase):
    def makegen(self):
        for x in range(24):
            yield x**2

    def test_types(self):
        ai32 = fromiter(self.makegen(), int32)
        ai64 = fromiter(self.makegen(), int64)
        af = fromiter(self.makegen(), float)
        self.assertTrue(ai32.dtype == dtype(int32))
        self.assertTrue(ai64.dtype == dtype(int64))
        self.assertTrue(af.dtype == dtype(float))

    def test_lengths(self):
        expected = array(list(self.makegen()))
        a = fromiter(self.makegen(), int)
        a20 = fromiter(self.makegen(), int, 20)
        self.assertTrue(len(a) == len(expected))
        self.assertTrue(len(a20) == 20)
        self.assertRaises(ValueError, fromiter,
                          self.makegen(), int, len(expected) + 10)

    def test_values(self):
        expected = array(list(self.makegen()))
        a = fromiter(self.makegen(), int)
        a20 = fromiter(self.makegen(), int, 20)
        self.assertTrue(alltrue(a == expected, axis=0))
        self.assertTrue(alltrue(a20 == expected[:20], axis=0))

    def load_data(self, n, eindex):
        # Utility method for the issue 2592 tests.
        # Raise an exception at the desired index in the iterator.
        for e in range(n):
            if e == eindex:
                raise NIterError('error at index %s' % eindex)
            yield e

    def test_2592(self):
        # Test iteration exceptions are correctly raised.
        count, eindex = 10, 5
        self.assertRaises(NIterError, np.fromiter,
                          self.load_data(count, eindex), dtype=int, count=count)

    def test_2592_edge(self):
        # Test iter. exceptions, edge case (exception at end of iterator).
        count = 10
        eindex = count-1
        self.assertRaises(NIterError, np.fromiter,
                          self.load_data(count, eindex), dtype=int, count=count)


class TestNonzero(TestCase):
    def test_nonzero_trivial(self):
        assert_equal(np.count_nonzero(array([])), 0)
        assert_equal(np.count_nonzero(array([], dtype='?')), 0)
        assert_equal(np.nonzero(array([])), ([],))

        assert_equal(np.count_nonzero(array(0)), 0)
        assert_equal(np.count_nonzero(array(0, dtype='?')), 0)
        assert_equal(np.nonzero(array(0)), ([],))
        assert_equal(np.count_nonzero(array(1)), 1)
        assert_equal(np.count_nonzero(array(1, dtype='?')), 1)
        assert_equal(np.nonzero(array(1)), ([0],))

    def test_nonzero_onedim(self):
        x = array([1, 0, 2, -1, 0, 0, 8])
        assert_equal(np.count_nonzero(x), 4)
        assert_equal(np.count_nonzero(x), 4)
        assert_equal(np.nonzero(x), ([0, 2, 3, 6],))

        x = array([(1, 2), (0, 0), (1, 1), (-1, 3), (0, 7)],
                            dtype=[('a', 'i4'), ('b', 'i2')])
        assert_equal(np.count_nonzero(x['a']), 3)
        assert_equal(np.count_nonzero(x['b']), 4)
        assert_equal(np.nonzero(x['a']), ([0, 2, 3],))
        assert_equal(np.nonzero(x['b']), ([0, 2, 3, 4],))

    def test_nonzero_twodim(self):
        x = array([[0, 1, 0], [2, 0, 3]])
        assert_equal(np.count_nonzero(x), 3)
        assert_equal(np.nonzero(x), ([0, 1, 1], [1, 0, 2]))

        x = np.eye(3)
        assert_equal(np.count_nonzero(x), 3)
        assert_equal(np.nonzero(x), ([0, 1, 2], [0, 1, 2]))

        x = array([[(0, 1), (0, 0), (1, 11)],
                   [(1, 1), (1, 0), (0, 0)],
                   [(0, 0), (1, 5), (0, 1)]], dtype=[('a', 'f4'), ('b', 'u1')])
        assert_equal(np.count_nonzero(x['a']), 4)
        assert_equal(np.count_nonzero(x['b']), 5)
        assert_equal(np.nonzero(x['a']), ([0, 1, 1, 2], [2, 0, 1, 1]))
        assert_equal(np.nonzero(x['b']), ([0, 0, 1, 2, 2], [0, 2, 0, 1, 2]))

        assert_(not x['a'].T.flags.aligned)
        assert_equal(np.count_nonzero(x['a'].T), 4)
        assert_equal(np.count_nonzero(x['b'].T), 5)
        assert_equal(np.nonzero(x['a'].T), ([0, 1, 1, 2], [1, 1, 2, 0]))
        assert_equal(np.nonzero(x['b'].T), ([0, 0, 1, 2, 2], [0, 1, 2, 0, 2]))

    def test_sparse(self):
        # test special sparse condition boolean code path
        for i in range(20):
            c = np.zeros(200, dtype=np.bool)
            c[i::20] = True
            assert_equal(np.nonzero(c)[0], np.arange(i, 200 + i, 20))

            c = np.zeros(400, dtype=np.bool)
            c[10 + i:20 + i] = True
            c[20 + i*2] = True
            assert_equal(np.nonzero(c)[0],
                         np.concatenate((np.arange(10 +i, 20 + i), [20 +i*2])))


class TestIndex(TestCase):
    def test_boolean(self):
        a = rand(3, 5, 8)
        V = rand(5, 8)
        g1 = randint(0, 5, size=15)
        g2 = randint(0, 8, size=15)
        V[g1, g2] = -V[g1, g2]
        assert_((array([a[0][V>0], a[1][V>0], a[2][V>0]]) == a[:, V>0]).all())

    def test_boolean_edgecase(self):
        a = np.array([], dtype='int32')
        b = np.array([], dtype='bool')
        c = a[b]
        assert_equal(c, [])
        assert_equal(c.dtype, np.dtype('int32'))


class TestBinaryRepr(TestCase):
    def test_zero(self):
        assert_equal(binary_repr(0), '0')

    def test_large(self):
        assert_equal(binary_repr(10736848), '101000111101010011010000')

    def test_negative(self):
        assert_equal(binary_repr(-1), '-1')
        assert_equal(binary_repr(-1, width=8), '11111111')

class TestBaseRepr(TestCase):
    def test_base3(self):
        assert_equal(base_repr(3**5, 3), '100000')

    def test_positive(self):
        assert_equal(base_repr(12, 10), '12')
        assert_equal(base_repr(12, 10, 4), '000012')
        assert_equal(base_repr(12, 4), '30')
        assert_equal(base_repr(3731624803700888, 36), '10QR0ROFCEW')

    def test_negative(self):
        assert_equal(base_repr(-12, 10), '-12')
        assert_equal(base_repr(-12, 10, 4), '-000012')
        assert_equal(base_repr(-12, 4), '-30')

class TestArrayComparisons(TestCase):
    def test_array_equal(self):
        res = array_equal(array([1, 2]), array([1, 2]))
        assert_(res)
        assert_(type(res) is bool)
        res = array_equal(array([1, 2]), array([1, 2, 3]))
        assert_(not res)
        assert_(type(res) is bool)
        res = array_equal(array([1, 2]), array([3, 4]))
        assert_(not res)
        assert_(type(res) is bool)
        res = array_equal(array([1, 2]), array([1, 3]))
        assert_(not res)
        assert_(type(res) is bool)
        res = array_equal(array(['a'], dtype='S1'), array(['a'], dtype='S1'))
        assert_(res)
        assert_(type(res) is bool)
        res = array_equal(array([('a', 1)], dtype='S1,u4'), array([('a', 1)], dtype='S1,u4'))
        assert_(res)
        assert_(type(res) is bool)

    def test_array_equiv(self):
        res = array_equiv(array([1, 2]), array([1, 2]))
        assert_(res)
        assert_(type(res) is bool)
        res = array_equiv(array([1, 2]), array([1, 2, 3]))
        assert_(not res)
        assert_(type(res) is bool)
        res = array_equiv(array([1, 2]), array([3, 4]))
        assert_(not res)
        assert_(type(res) is bool)
        res = array_equiv(array([1, 2]), array([1, 3]))
        assert_(not res)
        assert_(type(res) is bool)

        res = array_equiv(array([1, 1]), array([1]))
        assert_(res)
        assert_(type(res) is bool)
        res = array_equiv(array([1, 1]), array([[1], [1]]))
        assert_(res)
        assert_(type(res) is bool)
        res = array_equiv(array([1, 2]), array([2]))
        assert_(not res)
        assert_(type(res) is bool)
        res = array_equiv(array([1, 2]), array([[1], [2]]))
        assert_(not res)
        assert_(type(res) is bool)
        res = array_equiv(array([1, 2]), array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]))
        assert_(not res)
        assert_(type(res) is bool)


def assert_array_strict_equal(x, y):
    assert_array_equal(x, y)
    # Check flags, 32 bit arches typically don't provide 16 byte alignment
    if ((x.dtype.alignment <= 8 or
            np.intp().dtype.itemsize != 4) and
            sys.platform != 'win32'):
        assert_(x.flags == y.flags)
    else:
        assert_(x.flags.owndata == y.flags.owndata)
        assert_(x.flags.writeable == y.flags.writeable)
        assert_(x.flags.c_contiguous == y.flags.c_contiguous)
        assert_(x.flags.f_contiguous == y.flags.f_contiguous)
        assert_(x.flags.updateifcopy == y.flags.updateifcopy)
    # check endianness
    assert_(x.dtype.isnative == y.dtype.isnative)


class TestClip(TestCase):
    def setUp(self):
        self.nr = 5
        self.nc = 3

    def fastclip(self, a, m, M, out=None):
        if out is None:
            return a.clip(m, M)
        else:
            return a.clip(m, M, out)

    def clip(self, a, m, M, out=None):
        # use slow-clip
        selector = less(a, m)+2*greater(a, M)
        return selector.choose((a, m, M), out=out)

    # Handy functions
    def _generate_data(self, n, m):
        return randn(n, m)

    def _generate_data_complex(self, n, m):
        return randn(n, m) + 1.j *rand(n, m)

    def _generate_flt_data(self, n, m):
        return (randn(n, m)).astype(float32)

    def _neg_byteorder(self, a):
        a = asarray(a)
        if sys.byteorder == 'little':
            a = a.astype(a.dtype.newbyteorder('>'))
        else:
            a = a.astype(a.dtype.newbyteorder('<'))
        return a

    def _generate_non_native_data(self, n, m):
        data = randn(n, m)
        data = self._neg_byteorder(data)
        assert_(not data.dtype.isnative)
        return data

    def _generate_int_data(self, n, m):
        return (10 * rand(n, m)).astype(int64)

    def _generate_int32_data(self, n, m):
        return (10 * rand(n, m)).astype(int32)

    # Now the real test cases
    def test_simple_double(self):
        #Test native double input with scalar min/max.
        a   = self._generate_data(self.nr, self.nc)
        m   = 0.1
        M   = 0.6
        ac  = self.fastclip(a, m, M)
        act = self.clip(a, m, M)
        assert_array_strict_equal(ac, act)

    def test_simple_int(self):
        #Test native int input with scalar min/max.
        a   = self._generate_int_data(self.nr, self.nc)
        a   = a.astype(int)
        m   = -2
        M   = 4
        ac  = self.fastclip(a, m, M)
        act = self.clip(a, m, M)
        assert_array_strict_equal(ac, act)

    def test_array_double(self):
        #Test native double input with array min/max.
        a   = self._generate_data(self.nr, self.nc)
        m   = zeros(a.shape)
        M   = m + 0.5
        ac  = self.fastclip(a, m, M)
        act = self.clip(a, m, M)
        assert_array_strict_equal(ac, act)

    def test_simple_nonnative(self):
        #Test non native double input with scalar min/max.
        #Test native double input with non native double scalar min/max.
        a   = self._generate_non_native_data(self.nr, self.nc)
        m   = -0.5
        M   = 0.6
        ac  = self.fastclip(a, m, M)
        act = self.clip(a, m, M)
        assert_array_equal(ac, act)

        #Test native double input with non native double scalar min/max.
        a   = self._generate_data(self.nr, self.nc)
        m   = -0.5
        M   = self._neg_byteorder(0.6)
        assert_(not M.dtype.isnative)
        ac  = self.fastclip(a, m, M)
        act = self.clip(a, m, M)
        assert_array_equal(ac, act)

    def test_simple_complex(self):
        #Test native complex input with native double scalar min/max.
        #Test native input with complex double scalar min/max.
        a   = 3 * self._generate_data_complex(self.nr, self.nc)
        m   = -0.5
        M   = 1.
        ac  = self.fastclip(a, m, M)
        act = self.clip(a, m, M)
        assert_array_strict_equal(ac, act)

        #Test native input with complex double scalar min/max.
        a   = 3 * self._generate_data(self.nr, self.nc)
        m   = -0.5 + 1.j
        M   = 1. + 2.j
        ac  = self.fastclip(a, m, M)
        act = self.clip(a, m, M)
        assert_array_strict_equal(ac, act)

    def test_clip_non_contig(self):
        #Test clip for non contiguous native input and native scalar min/max.
        a   = self._generate_data(self.nr * 2, self.nc * 3)
        a   = a[::2, ::3]
        assert_(not a.flags['F_CONTIGUOUS'])
        assert_(not a.flags['C_CONTIGUOUS'])
        ac  = self.fastclip(a, -1.6, 1.7)
        act = self.clip(a, -1.6, 1.7)
        assert_array_strict_equal(ac, act)

    def test_simple_out(self):
        #Test native double input with scalar min/max.
        a   = self._generate_data(self.nr, self.nc)
        m   = -0.5
        M   = 0.6
        ac  = zeros(a.shape)
        act = zeros(a.shape)
        self.fastclip(a, m, M, ac)
        self.clip(a, m, M, act)
        assert_array_strict_equal(ac, act)

    def test_simple_int32_inout(self):
        #Test native int32 input with double min/max and int32 out.
        a   = self._generate_int32_data(self.nr, self.nc)
        m   = float64(0)
        M   = float64(2)
        ac  = zeros(a.shape, dtype = int32)
        act = ac.copy()
        self.fastclip(a, m, M, ac)
        self.clip(a, m, M, act)
        assert_array_strict_equal(ac, act)

    def test_simple_int64_out(self):
        #Test native int32 input with int32 scalar min/max and int64 out.
        a   = self._generate_int32_data(self.nr, self.nc)
        m   = int32(-1)
        M   = int32(1)
        ac  = zeros(a.shape, dtype = int64)
        act = ac.copy()
        self.fastclip(a, m, M, ac)
        self.clip(a, m, M, act)
        assert_array_strict_equal(ac, act)

    def test_simple_int64_inout(self):
        #Test native int32 input with double array min/max and int32 out.
        a   = self._generate_int32_data(self.nr, self.nc)
        m   = zeros(a.shape, float64)
        M   = float64(1)
        ac  = zeros(a.shape, dtype = int32)
        act = ac.copy()
        self.fastclip(a, m, M, ac)
        self.clip(a, m, M, act)
        assert_array_strict_equal(ac, act)

    def test_simple_int32_out(self):
        #Test native double input with scalar min/max and int out.
        a   = self._generate_data(self.nr, self.nc)
        m   = -1.0
        M   = 2.0
        ac  = zeros(a.shape, dtype = int32)
        act = ac.copy()
        self.fastclip(a, m, M, ac)
        self.clip(a, m, M, act)
        assert_array_strict_equal(ac, act)

    def test_simple_inplace_01(self):
        #Test native double input with array min/max in-place.
        a   = self._generate_data(self.nr, self.nc)
        ac  = a.copy()
        m   = zeros(a.shape)
        M   = 1.0
        self.fastclip(a, m, M, a)
        self.clip(a, m, M, ac)
        assert_array_strict_equal(a, ac)

    def test_simple_inplace_02(self):
        #Test native double input with scalar min/max in-place.
        a   = self._generate_data(self.nr, self.nc)
        ac  = a.copy()
        m   = -0.5
        M   = 0.6
        self.fastclip(a, m, M, a)
        self.clip(a, m, M, ac)
        assert_array_strict_equal(a, ac)

    def test_noncontig_inplace(self):
        #Test non contiguous double input with double scalar min/max in-place.
        a   = self._generate_data(self.nr * 2, self.nc * 3)
        a   = a[::2, ::3]
        assert_(not a.flags['F_CONTIGUOUS'])
        assert_(not a.flags['C_CONTIGUOUS'])
        ac  = a.copy()
        m   = -0.5
        M   = 0.6
        self.fastclip(a, m, M, a)
        self.clip(a, m, M, ac)
        assert_array_equal(a, ac)

    def test_type_cast_01(self):
        #Test native double input with scalar min/max.
        a   = self._generate_data(self.nr, self.nc)
        m   = -0.5
        M   = 0.6
        ac  = self.fastclip(a, m, M)
        act = self.clip(a, m, M)
        assert_array_strict_equal(ac, act)

    def test_type_cast_02(self):
        #Test native int32 input with int32 scalar min/max.
        a   = self._generate_int_data(self.nr, self.nc)
        a   = a.astype(int32)
        m   = -2
        M   = 4
        ac  = self.fastclip(a, m, M)
        act = self.clip(a, m, M)
        assert_array_strict_equal(ac, act)

    def test_type_cast_03(self):
        #Test native int32 input with float64 scalar min/max.
        a   = self._generate_int32_data(self.nr, self.nc)
        m   = -2
        M   = 4
        ac  = self.fastclip(a, float64(m), float64(M))
        act = self.clip(a, float64(m), float64(M))
        assert_array_strict_equal(ac, act)

    def test_type_cast_04(self):
        #Test native int32 input with float32 scalar min/max.
        a   = self._generate_int32_data(self.nr, self.nc)
        m   = float32(-2)
        M   = float32(4)
        act = self.fastclip(a, m, M)
        ac  = self.clip(a, m, M)
        assert_array_strict_equal(ac, act)

    def test_type_cast_05(self):
        #Test native int32 with double arrays min/max.
        a   = self._generate_int_data(self.nr, self.nc)
        m   = -0.5
        M   = 1.
        ac  = self.fastclip(a, m * zeros(a.shape), M)
        act = self.clip(a, m * zeros(a.shape), M)
        assert_array_strict_equal(ac, act)

    def test_type_cast_06(self):
        #Test native with NON native scalar min/max.
        a   = self._generate_data(self.nr, self.nc)
        m   = 0.5
        m_s = self._neg_byteorder(m)
        M   = 1.
        act = self.clip(a, m_s, M)
        ac  = self.fastclip(a, m_s, M)
        assert_array_strict_equal(ac, act)

    def test_type_cast_07(self):
        #Test NON native with native array min/max.
        a   = self._generate_data(self.nr, self.nc)
        m   = -0.5 * ones(a.shape)
        M   = 1.
        a_s = self._neg_byteorder(a)
        assert_(not a_s.dtype.isnative)
        act = a_s.clip(m, M)
        ac  = self.fastclip(a_s, m, M)
        assert_array_strict_equal(ac, act)

    def test_type_cast_08(self):
        #Test NON native with native scalar min/max.
        a   = self._generate_data(self.nr, self.nc)
        m   = -0.5
        M   = 1.
        a_s = self._neg_byteorder(a)
        assert_(not a_s.dtype.isnative)
        ac  = self.fastclip(a_s, m, M)
        act = a_s.clip(m, M)
        assert_array_strict_equal(ac, act)

    def test_type_cast_09(self):
        #Test native with NON native array min/max.
        a   = self._generate_data(self.nr, self.nc)
        m   = -0.5 * ones(a.shape)
        M   = 1.
        m_s = self._neg_byteorder(m)
        assert_(not m_s.dtype.isnative)
        ac  = self.fastclip(a, m_s, M)
        act = self.clip(a, m_s, M)
        assert_array_strict_equal(ac, act)

    def test_type_cast_10(self):
        #Test native int32 with float min/max and float out for output argument.
        a   = self._generate_int_data(self.nr, self.nc)
        b   = zeros(a.shape, dtype = float32)
        m   = float32(-0.5)
        M   = float32(1)
        act = self.clip(a, m, M, out = b)
        ac  = self.fastclip(a, m, M, out = b)
        assert_array_strict_equal(ac, act)

    def test_type_cast_11(self):
        #Test non native with native scalar, min/max, out non native
        a   = self._generate_non_native_data(self.nr, self.nc)
        b   = a.copy()
        b   = b.astype(b.dtype.newbyteorder('>'))
        bt  = b.copy()
        m   = -0.5
        M   = 1.
        self.fastclip(a, m, M, out = b)
        self.clip(a, m, M, out = bt)
        assert_array_strict_equal(b, bt)

    def test_type_cast_12(self):
        #Test native int32 input and min/max and float out
        a   = self._generate_int_data(self.nr, self.nc)
        b   = zeros(a.shape, dtype = float32)
        m   = int32(0)
        M   = int32(1)
        act = self.clip(a, m, M, out = b)
        ac  = self.fastclip(a, m, M, out = b)
        assert_array_strict_equal(ac, act)

    def test_clip_with_out_simple(self):
        #Test native double input with scalar min/max
        a   = self._generate_data(self.nr, self.nc)
        m   = -0.5
        M   = 0.6
        ac  = zeros(a.shape)
        act = zeros(a.shape)
        self.fastclip(a, m, M, ac)
        self.clip(a, m, M, act)
        assert_array_strict_equal(ac, act)

    def test_clip_with_out_simple2(self):
        #Test native int32 input with double min/max and int32 out
        a   = self._generate_int32_data(self.nr, self.nc)
        m   = float64(0)
        M   = float64(2)
        ac  = zeros(a.shape, dtype = int32)
        act = ac.copy()
        self.fastclip(a, m, M, ac)
        self.clip(a, m, M, act)
        assert_array_strict_equal(ac, act)

    def test_clip_with_out_simple_int32(self):
        #Test native int32 input with int32 scalar min/max and int64 out
        a   = self._generate_int32_data(self.nr, self.nc)
        m   = int32(-1)
        M   = int32(1)
        ac  = zeros(a.shape, dtype = int64)
        act = ac.copy()
        self.fastclip(a, m, M, ac)
        self.clip(a, m, M, act)
        assert_array_strict_equal(ac, act)

    def test_clip_with_out_array_int32(self):
        #Test native int32 input with double array min/max and int32 out
        a   = self._generate_int32_data(self.nr, self.nc)
        m   = zeros(a.shape, float64)
        M   = float64(1)
        ac  = zeros(a.shape, dtype = int32)
        act = ac.copy()
        self.fastclip(a, m, M, ac)
        self.clip(a, m, M, act)
        assert_array_strict_equal(ac, act)

    def test_clip_with_out_array_outint32(self):
        #Test native double input with scalar min/max and int out
        a   = self._generate_data(self.nr, self.nc)
        m   = -1.0
        M   = 2.0
        ac  = zeros(a.shape, dtype = int32)
        act = ac.copy()
        self.fastclip(a, m, M, ac)
        self.clip(a, m, M, act)
        assert_array_strict_equal(ac, act)

    def test_clip_inplace_array(self):
        #Test native double input with array min/max
        a   = self._generate_data(self.nr, self.nc)
        ac  = a.copy()
        m   = zeros(a.shape)
        M   = 1.0
        self.fastclip(a, m, M, a)
        self.clip(a, m, M, ac)
        assert_array_strict_equal(a, ac)

    def test_clip_inplace_simple(self):
        #Test native double input with scalar min/max
        a   = self._generate_data(self.nr, self.nc)
        ac  = a.copy()
        m   = -0.5
        M   = 0.6
        self.fastclip(a, m, M, a)
        self.clip(a, m, M, ac)
        assert_array_strict_equal(a, ac)

    def test_clip_func_takes_out(self):
        # Ensure that the clip() function takes an out= argument.
        a = self._generate_data(self.nr, self.nc)
        ac = a.copy()
        m = -0.5
        M = 0.6
        a2 = clip(a, m, M, out=a)
        self.clip(a, m, M, ac)
        assert_array_strict_equal(a2, ac)
        self.assertTrue(a2 is a)


class TestAllclose(object):
    rtol = 1e-5
    atol = 1e-8

    def setUp(self):
        self.olderr = np.seterr(invalid='ignore')

    def tearDown(self):
        np.seterr(**self.olderr)

    def tst_allclose(self, x, y):
        assert_(allclose(x, y), "%s and %s not close" % (x, y))

    def tst_not_allclose(self, x, y):
        assert_(not allclose(x, y), "%s and %s shouldn't be close" % (x, y))

    def test_ip_allclose(self):
        #Parametric test factory.
        arr = array([100, 1000])
        aran = arange(125).reshape((5, 5, 5))

        atol = self.atol
        rtol = self.rtol

        data = [([1, 0], [1, 0]),
                ([atol], [0]),
                ([1], [1+rtol+atol]),
                (arr, arr + arr*rtol),
                (arr, arr + arr*rtol + atol*2),
                (aran, aran + aran*rtol),
                (inf, inf),
                (inf, [inf])]

        for (x, y) in data:
            yield (self.tst_allclose, x, y)

    def test_ip_not_allclose(self):
        #Parametric test factory.
        aran = arange(125).reshape((5, 5, 5))

        atol = self.atol
        rtol = self.rtol

        data = [([inf, 0], [1, inf]),
                ([inf, 0], [1, 0]),
                ([inf, inf], [1, inf]),
                ([inf, inf], [1, 0]),
                ([-inf, 0], [inf, 0]),
                ([nan, 0], [nan, 0]),
                ([atol*2], [0]),
                ([1], [1+rtol+atol*2]),
                (aran, aran + aran*atol + atol*2),
                (array([inf, 1]), array([0, inf]))]

        for (x, y) in data:
            yield (self.tst_not_allclose, x, y)

    def test_no_parameter_modification(self):
        x = array([inf, 1])
        y = array([0, inf])
        allclose(x, y)
        assert_array_equal(x, array([inf, 1]))
        assert_array_equal(y, array([0, inf]))


    def test_min_int(self):
        # Could make problems because of abs(min_int) == min_int
        min_int = np.iinfo(np.int_).min
        a = np.array([min_int], dtype=np.int_)
        assert_(allclose(a, a))


class TestIsclose(object):
    rtol = 1e-5
    atol = 1e-8

    def setup(self):
        atol = self.atol
        rtol = self.rtol
        arr = array([100, 1000])
        aran = arange(125).reshape((5, 5, 5))

        self.all_close_tests = [
                ([1, 0], [1, 0]),
                ([atol], [0]),
                ([1], [1 + rtol + atol]),
                (arr, arr + arr*rtol),
                (arr, arr + arr*rtol + atol),
                (aran, aran + aran*rtol),
                (inf, inf),
                (inf, [inf]),
                ([inf, -inf], [inf, -inf]),
                ]
        self.none_close_tests = [
                ([inf, 0], [1, inf]),
                ([inf, -inf], [1, 0]),
                ([inf, inf], [1, -inf]),
                ([inf, inf], [1, 0]),
                ([nan, 0], [nan, -inf]),
                ([atol*2], [0]),
                ([1], [1 + rtol + atol*2]),
                (aran, aran + rtol*1.1*aran + atol*1.1),
                (array([inf, 1]), array([0, inf])),
                ]
        self.some_close_tests = [
                ([inf, 0], [inf, atol*2]),
                ([atol, 1, 1e6*(1 + 2*rtol) + atol], [0, nan, 1e6]),
                (arange(3), [0, 1, 2.1]),
                (nan, [nan, nan, nan]),
                ([0], [atol, inf, -inf, nan]),
                (0, [atol, inf, -inf, nan]),
                ]
        self.some_close_results = [
                [True, False],
                [True, False, False],
                [True, True, False],
                [False, False, False],
                [True, False, False, False],
                [True, False, False, False],
                ]

    def test_ip_isclose(self):
        self.setup()
        tests = self.some_close_tests
        results = self.some_close_results
        for (x, y), result in zip(tests, results):
            yield (assert_array_equal, isclose(x, y), result)

    def tst_all_isclose(self, x, y):
        assert_(all(isclose(x, y)), "%s and %s not close" % (x, y))

    def tst_none_isclose(self, x, y):
        msg = "%s and %s shouldn't be close"
        assert_(not any(isclose(x, y)), msg % (x, y))

    def tst_isclose_allclose(self, x, y):
        msg = "isclose.all() and allclose aren't same for %s and %s"
        assert_array_equal(isclose(x, y).all(), allclose(x, y), msg % (x, y))

    def test_ip_all_isclose(self):
        self.setup()
        for (x, y) in self.all_close_tests:
            yield (self.tst_all_isclose, x, y)

    def test_ip_none_isclose(self):
        self.setup()
        for (x, y) in self.none_close_tests:
            yield (self.tst_none_isclose, x, y)

    def test_ip_isclose_allclose(self):
        self.setup()
        tests = (self.all_close_tests + self.none_close_tests +
                 self.some_close_tests)
        for (x, y) in tests:
            yield (self.tst_isclose_allclose, x, y)

    def test_equal_nan(self):
        assert_array_equal(isclose(nan, nan, equal_nan=True), [True])
        arr = array([1.0, nan])
        assert_array_equal(isclose(arr, arr, equal_nan=True), [True, True])

    def test_masked_arrays(self):
        x = np.ma.masked_where([True, True, False], np.arange(3))
        assert_(type(x) is type(isclose(2, x)))

        x = np.ma.masked_where([True, True, False], [nan, inf, nan])
        assert_(type(x) is type(isclose(inf, x)))

        x = np.ma.masked_where([True, True, False], [nan, nan, nan])
        y = isclose(nan, x, equal_nan=True)
        assert_(type(x) is type(y))
        # Ensure that the mask isn't modified...
        assert_array_equal([True, True, False], y.mask)

        x = np.ma.masked_where([True, True, False], [nan, nan, nan])
        y = isclose(x, x, equal_nan=True)
        assert_(type(x) is type(y))
        # Ensure that the mask isn't modified...
        assert_array_equal([True, True, False], y.mask)

    def test_scalar_return(self):
        assert_(isscalar(isclose(1, 1)))

    def test_no_parameter_modification(self):
        x = array([inf, 1])
        y = array([0, inf])
        isclose(x, y)
        assert_array_equal(x, array([inf, 1]))
        assert_array_equal(y, array([0, inf]))

class TestStdVar(TestCase):
    def setUp(self):
        self.A = array([1, -1, 1, -1])
        self.real_var = 1

    def test_basic(self):
        assert_almost_equal(var(self.A), self.real_var)
        assert_almost_equal(std(self.A)**2, self.real_var)

    def test_scalars(self):
        assert_equal(var(1), 0)
        assert_equal(std(1), 0)

    def test_ddof1(self):
        assert_almost_equal(var(self.A, ddof=1),
                            self.real_var*len(self.A)/float(len(self.A)-1))
        assert_almost_equal(std(self.A, ddof=1)**2,
                            self.real_var*len(self.A)/float(len(self.A)-1))

    def test_ddof2(self):
        assert_almost_equal(var(self.A, ddof=2),
                            self.real_var*len(self.A)/float(len(self.A)-2))
        assert_almost_equal(std(self.A, ddof=2)**2,
                            self.real_var*len(self.A)/float(len(self.A)-2))

class TestStdVarComplex(TestCase):
    def test_basic(self):
        A = array([1, 1.j, -1, -1.j])
        real_var = 1
        assert_almost_equal(var(A), real_var)
        assert_almost_equal(std(A)**2, real_var)

    def test_scalars(self):
        assert_equal(var(1j), 0)
        assert_equal(std(1j), 0)


class TestCreationFuncs(TestCase):
    #Test ones, zeros, empty and filled

    def setUp(self):
        self.dtypes = ('b', 'i', 'u', 'f', 'c', 'S', 'a', 'U', 'V')
        self.orders = {'C': 'c_contiguous', 'F': 'f_contiguous'}
        self.ndims = 10

    def check_function(self, func, fill_value=None):
        par = (
            (0, 1, 2),
            range(self.ndims),
            self.orders,
            self.dtypes,
            2**np.arange(9)
        )
        fill_kwarg = {}
        if fill_value is not None:
            fill_kwarg = {'fill_value': fill_value}
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', DeprecationWarning)
            for size, ndims, order, type, bytes in itertools.product(*par):
                shape = ndims * [size]
                try:
                    dtype = np.dtype('{0}{1}'.format(type, bytes))
                except TypeError: # dtype combination does not exist
                    continue
                else:
                    # do not fill void type
                    if fill_value is not None and type in 'V':
                        continue

                    arr = func(shape, order=order, dtype=dtype,
                               **fill_kwarg)

                    assert_(arr.dtype == dtype)
                    assert_(getattr(arr.flags, self.orders[order]))

                    if fill_value is not None:
                        if dtype.str.startswith('|S'):
                            val = str(fill_value)
                        else:
                            val = fill_value
                        assert_equal(arr, dtype.type(val))

    def test_zeros(self):
        self.check_function(np.zeros)

    def test_ones(self):
        self.check_function(np.zeros)

    def test_empty(self):
        self.check_function(np.empty)

    def test_filled(self):
        self.check_function(np.full, 0)
        self.check_function(np.full, 1)

    def test_for_reference_leak(self):
        # Make sure we have an object for reference
        dim = 1
        beg = sys.getrefcount(dim)
        np.zeros([dim]*10)
        assert_(sys.getrefcount(dim) == beg)
        np.ones([dim]*10)
        assert_(sys.getrefcount(dim) == beg)
        np.empty([dim]*10)
        assert_(sys.getrefcount(dim) == beg)
        np.full([dim]*10, 0)
        assert_(sys.getrefcount(dim) == beg)



class TestLikeFuncs(TestCase):
    '''Test ones_like, zeros_like, empty_like and full_like'''

    def setUp(self):
        self.data = [
                # Array scalars
                (array(3.), None),
                (array(3), 'f8'),
                # 1D arrays
                (arange(6, dtype='f4'), None),
                (arange(6), 'c16'),
                # 2D C-layout arrays
                (arange(6).reshape(2, 3), None),
                (arange(6).reshape(3, 2), 'i1'),
                # 2D F-layout arrays
                (arange(6).reshape((2, 3), order='F'), None),
                (arange(6).reshape((3, 2), order='F'), 'i1'),
                # 3D C-layout arrays
                (arange(24).reshape(2, 3, 4), None),
                (arange(24).reshape(4, 3, 2), 'f4'),
                # 3D F-layout arrays
                (arange(24).reshape((2, 3, 4), order='F'), None),
                (arange(24).reshape((4, 3, 2), order='F'), 'f4'),
                # 3D non-C/F-layout arrays
                (arange(24).reshape(2, 3, 4).swapaxes(0, 1), None),
                (arange(24).reshape(4, 3, 2).swapaxes(0, 1), '?'),
                     ]

    def compare_array_value(self, dz, value, fill_value):
        if value is not None:
            if fill_value:
                try:
                    z = dz.dtype.type(value)
                except OverflowError:
                    pass
                else:
                    assert_(all(dz == z))
            else:
                assert_(all(dz == value))

    def check_like_function(self, like_function, value, fill_value=False):
        if fill_value:
            fill_kwarg = {'fill_value': value}
        else:
            fill_kwarg = {}
        for d, dtype in self.data:
            # default (K) order, dtype
            dz = like_function(d, dtype=dtype, **fill_kwarg)
            assert_equal(dz.shape, d.shape)
            assert_equal(array(dz.strides)*d.dtype.itemsize,
                         array(d.strides)*dz.dtype.itemsize)
            assert_equal(d.flags.c_contiguous, dz.flags.c_contiguous)
            assert_equal(d.flags.f_contiguous, dz.flags.f_contiguous)
            if dtype is None:
                assert_equal(dz.dtype, d.dtype)
            else:
                assert_equal(dz.dtype, np.dtype(dtype))
            self.compare_array_value(dz, value, fill_value)

            # C order, default dtype
            dz = like_function(d, order='C', dtype=dtype, **fill_kwarg)
            assert_equal(dz.shape, d.shape)
            assert_(dz.flags.c_contiguous)
            if dtype is None:
                assert_equal(dz.dtype, d.dtype)
            else:
                assert_equal(dz.dtype, np.dtype(dtype))
            self.compare_array_value(dz, value, fill_value)

            # F order, default dtype
            dz = like_function(d, order='F', dtype=dtype, **fill_kwarg)
            assert_equal(dz.shape, d.shape)
            assert_(dz.flags.f_contiguous)
            if dtype is None:
                assert_equal(dz.dtype, d.dtype)
            else:
                assert_equal(dz.dtype, np.dtype(dtype))
            self.compare_array_value(dz, value, fill_value)

            # A order
            dz = like_function(d, order='A', dtype=dtype, **fill_kwarg)
            assert_equal(dz.shape, d.shape)
            if d.flags.f_contiguous:
                assert_(dz.flags.f_contiguous)
            else:
                assert_(dz.flags.c_contiguous)
            if dtype is None:
                assert_equal(dz.dtype, d.dtype)
            else:
                assert_equal(dz.dtype, np.dtype(dtype))
            self.compare_array_value(dz, value, fill_value)

        # Test the 'subok' parameter
        a = np.matrix([[1, 2], [3, 4]])

        b = like_function(a, **fill_kwarg)
        assert_(type(b) is np.matrix)

        b = like_function(a, subok=False, **fill_kwarg)
        assert_(type(b) is not np.matrix)

    def test_ones_like(self):
        self.check_like_function(np.ones_like, 1)

    def test_zeros_like(self):
        self.check_like_function(np.zeros_like, 0)

    def test_empty_like(self):
        self.check_like_function(np.empty_like, None)

    def test_filled_like(self):
        self.check_like_function(np.full_like, 0, True)
        self.check_like_function(np.full_like, 1, True)
        self.check_like_function(np.full_like, 1000, True)
        self.check_like_function(np.full_like, 123.456, True)
        self.check_like_function(np.full_like, np.inf, True)

class _TestCorrelate(TestCase):
    def _setup(self, dt):
        self.x = np.array([1, 2, 3, 4, 5], dtype=dt)
        self.y = np.array([-1, -2, -3], dtype=dt)
        self.z1 = np.array([ -3.,  -8., -14., -20., -26., -14.,  -5.], dtype=dt)
        self.z2 = np.array([ -5.,  -14., -26., -20., -14., -8.,  -3.], dtype=dt)

    def test_float(self):
        self._setup(np.float)
        z = np.correlate(self.x, self.y, 'full', old_behavior=self.old_behavior)
        assert_array_almost_equal(z, self.z1)
        z = np.correlate(self.y, self.x, 'full', old_behavior=self.old_behavior)
        assert_array_almost_equal(z, self.z2)

    def test_object(self):
        self._setup(Decimal)
        z = np.correlate(self.x, self.y, 'full', old_behavior=self.old_behavior)
        assert_array_almost_equal(z, self.z1)
        z = np.correlate(self.y, self.x, 'full', old_behavior=self.old_behavior)
        assert_array_almost_equal(z, self.z2)

class TestCorrelate(_TestCorrelate):
    old_behavior = True
    def _setup(self, dt):
        # correlate uses an unconventional definition so that correlate(a, b)
        # == correlate(b, a), so force the corresponding outputs to be the same
        # as well
        _TestCorrelate._setup(self, dt)
        self.z2 = self.z1

    @dec.deprecated()
    def test_complex(self):
        x = np.array([1, 2, 3, 4+1j], dtype=np.complex)
        y = np.array([-1, -2j, 3+1j], dtype=np.complex)
        r_z = np.array([3+1j, 6, 8-1j, 9+1j, -1-8j, -4-1j], dtype=np.complex)
        z = np.correlate(x, y, 'full', old_behavior=self.old_behavior)
        assert_array_almost_equal(z, r_z)

    @dec.deprecated()
    def test_float(self):
        _TestCorrelate.test_float(self)

    @dec.deprecated()
    def test_object(self):
        _TestCorrelate.test_object(self)

class TestCorrelateNew(_TestCorrelate):
    old_behavior = False
    def test_complex(self):
        x = np.array([1, 2, 3, 4+1j], dtype=np.complex)
        y = np.array([-1, -2j, 3+1j], dtype=np.complex)
        r_z = np.array([3-1j, 6, 8+1j, 11+5j, -5+8j, -4-1j], dtype=np.complex)
        #z = np.acorrelate(x, y, 'full')
        #assert_array_almost_equal(z, r_z)

        r_z = r_z[::-1].conjugate()
        z = np.correlate(y, x, 'full', old_behavior=self.old_behavior)
        assert_array_almost_equal(z, r_z)

class TestArgwhere(object):
    def test_2D(self):
        x = np.arange(6).reshape((2, 3))
        assert_array_equal(np.argwhere(x > 1),
                           [[0, 2],
                            [1, 0],
                            [1, 1],
                            [1, 2]])

    def test_list(self):
        assert_equal(np.argwhere([4, 0, 2, 1, 3]), [[0], [2], [3], [4]])

class TestStringFunction(object):
    def test_set_string_function(self):
        a = np.array([1])
        np.set_string_function(lambda x: "FOO", repr=True)
        assert_equal(repr(a), "FOO")
        np.set_string_function(None, repr=True)
        assert_equal(repr(a), "array([1])")

        np.set_string_function(lambda x: "FOO", repr=False)
        assert_equal(str(a), "FOO")
        np.set_string_function(None, repr=False)
        assert_equal(str(a), "[1]")

class TestRoll(TestCase):
    def test_roll1d(self):
        x = np.arange(10)
        xr = np.roll(x, 2)
        assert_equal(xr, np.array([8, 9, 0, 1, 2, 3, 4, 5, 6, 7]))

    def test_roll2d(self):
        x2 = np.reshape(np.arange(10), (2, 5))
        x2r = np.roll(x2, 1)
        assert_equal(x2r, np.array([[9, 0, 1, 2, 3], [4, 5, 6, 7, 8]]))

        x2r = np.roll(x2, 1, axis=0)
        assert_equal(x2r, np.array([[5, 6, 7, 8, 9], [0, 1, 2, 3, 4]]))

        x2r = np.roll(x2, 1, axis=1)
        assert_equal(x2r, np.array([[4, 0, 1, 2, 3], [9, 5, 6, 7, 8]]))

    def test_roll_empty(self):
        x = np.array([])
        assert_equal(np.roll(x, 1), np.array([]))

class TestCross(TestCase):
    def test_2x2(self):
        u = [1, 2]
        v = [3, 4]
        z = -2
        cp = np.cross(u, v)
        assert_equal(cp, z)
        cp = np.cross(v, u)
        assert_equal(cp, -z)

    def test_2x3(self):
        u = [1, 2]
        v = [3, 4, 5]
        z = np.array([10, -5, -2])
        cp = np.cross(u, v)
        assert_equal(cp, z)
        cp = np.cross(v, u)
        assert_equal(cp, -z)

    def test_3x3(self):
        u = [1, 2, 3]
        v = [4, 5, 6]
        z = np.array([-3, 6, -3])
        cp = cross(u, v)
        assert_equal(cp, z)
        cp = np.cross(v, u)
        assert_equal(cp, -z)

    def test_broadcasting(self):
        # Ticket #2624 (Trac #2032)
        u = np.tile([1, 2], (11, 1))
        v = np.tile([3, 4], (11, 1))
        z = -2
        assert_equal(np.cross(u, v), z)
        assert_equal(np.cross(v, u), -z)
        assert_equal(np.cross(u, u), 0)

        u = np.tile([1, 2], (11, 1)).T
        v = np.tile([3, 4, 5], (11, 1))
        z = np.tile([10, -5, -2], (11, 1))
        assert_equal(np.cross(u, v, axisa=0), z)
        assert_equal(np.cross(v, u.T), -z)
        assert_equal(np.cross(v, v), 0)

        u = np.tile([1, 2, 3], (11, 1)).T
        v = np.tile([3, 4], (11, 1)).T
        z = np.tile([-12, 9, -2], (11, 1))
        assert_equal(np.cross(u, v, axisa=0, axisb=0), z)
        assert_equal(np.cross(v.T, u.T), -z)
        assert_equal(np.cross(u.T, u.T), 0)

        u = np.tile([1, 2, 3], (5, 1))
        v = np.tile([4, 5, 6], (5, 1)).T
        z = np.tile([-3, 6, -3], (5, 1))
        assert_equal(np.cross(u, v, axisb=0), z)
        assert_equal(np.cross(v.T, u), -z)
        assert_equal(np.cross(u, u), 0)

    def test_broadcasting_shapes(self):
        u = np.ones((2, 1, 3))
        v = np.ones((5, 3))
        assert_equal(np.cross(u, v).shape, (2, 5, 3))
        u = np.ones((10, 3, 5))
        v = np.ones((2, 5))
        assert_equal(np.cross(u, v, axisa=1, axisb=0).shape, (10, 5, 3))
        assert_raises(ValueError, np.cross, u, v, axisa=1, axisb=2)
        assert_raises(ValueError, np.cross, u, v, axisa=3, axisb=0)
        u = np.ones((10, 3, 5, 7))
        v = np.ones((5, 7, 2))
        assert_equal(np.cross(u, v, axisa=1, axisc=2).shape, (10, 5, 3, 7))
        assert_raises(ValueError, np.cross, u, v, axisa=-5, axisb=2)
        assert_raises(ValueError, np.cross, u, v, axisa=1, axisb=-4)

def test_outer_out_param():
    arr1 = np.ones((5,))
    arr2 = np.ones((2,))
    arr3 = np.linspace(-2, 2, 5)
    out1 = np.ndarray(shape=(5,5))
    out2 = np.ndarray(shape=(2, 5))
    res1 = np.outer(arr1, arr3, out1)
    assert_equal(res1, out1)
    assert_equal(np.outer(arr2, arr3, out2), out2)

if __name__ == "__main__":
    run_module_suite()
