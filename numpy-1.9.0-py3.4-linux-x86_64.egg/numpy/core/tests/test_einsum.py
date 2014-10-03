from __future__ import division, absolute_import, print_function

import sys
import warnings
from decimal import Decimal

import numpy as np
from numpy.testing import *

class TestEinSum(TestCase):
    def test_einsum_errors(self):
        # Need enough arguments
        assert_raises(ValueError, np.einsum)
        assert_raises(ValueError, np.einsum, "")

        # subscripts must be a string
        assert_raises(TypeError, np.einsum, 0, 0)

        # out parameter must be an array
        assert_raises(TypeError, np.einsum, "", 0, out='test')

        # order parameter must be a valid order
        assert_raises(TypeError, np.einsum, "", 0, order='W')

        # casting parameter must be a valid casting
        assert_raises(ValueError, np.einsum, "", 0, casting='blah')

        # dtype parameter must be a valid dtype
        assert_raises(TypeError, np.einsum, "", 0, dtype='bad_data_type')

        # other keyword arguments are rejected
        assert_raises(TypeError, np.einsum, "", 0, bad_arg=0)

        # issue 4528 revealed a segfault with this call
        assert_raises(TypeError, np.einsum, *(None,)*63)

        # number of operands must match count in subscripts string
        assert_raises(ValueError, np.einsum, "", 0, 0)
        assert_raises(ValueError, np.einsum, ",", 0, [0], [0])
        assert_raises(ValueError, np.einsum, ",", [0])

        # can't have more subscripts than dimensions in the operand
        assert_raises(ValueError, np.einsum, "i", 0)
        assert_raises(ValueError, np.einsum, "ij", [0, 0])
        assert_raises(ValueError, np.einsum, "...i", 0)
        assert_raises(ValueError, np.einsum, "i...j", [0, 0])
        assert_raises(ValueError, np.einsum, "i...", 0)
        assert_raises(ValueError, np.einsum, "ij...", [0, 0])

        # invalid ellipsis
        assert_raises(ValueError, np.einsum, "i..", [0, 0])
        assert_raises(ValueError, np.einsum, ".i...", [0, 0])
        assert_raises(ValueError, np.einsum, "j->..j", [0, 0])
        assert_raises(ValueError, np.einsum, "j->.j...", [0, 0])

        # invalid subscript character
        assert_raises(ValueError, np.einsum, "i%...", [0, 0])
        assert_raises(ValueError, np.einsum, "...j$", [0, 0])
        assert_raises(ValueError, np.einsum, "i->&", [0, 0])

        # output subscripts must appear in input
        assert_raises(ValueError, np.einsum, "i->ij", [0, 0])

        # output subscripts may only be specified once
        assert_raises(ValueError, np.einsum, "ij->jij", [[0, 0], [0, 0]])

        # dimensions much match when being collapsed
        assert_raises(ValueError, np.einsum, "ii", np.arange(6).reshape(2, 3))
        assert_raises(ValueError, np.einsum, "ii->i", np.arange(6).reshape(2, 3))

        # broadcasting to new dimensions must be enabled explicitly
        assert_raises(ValueError, np.einsum, "i", np.arange(6).reshape(2, 3))
        assert_raises(ValueError, np.einsum, "i->i", [[0, 1], [0, 1]],
                                            out=np.arange(4).reshape(2, 2))

    def test_einsum_views(self):
        # pass-through
        a = np.arange(6)
        a.shape = (2, 3)

        b = np.einsum("...", a)
        assert_(b.base is a)

        b = np.einsum(a, [Ellipsis])
        assert_(b.base is a)

        b = np.einsum("ij", a)
        assert_(b.base is a)
        assert_equal(b, a)

        b = np.einsum(a, [0, 1])
        assert_(b.base is a)
        assert_equal(b, a)

        # transpose
        a = np.arange(6)
        a.shape = (2, 3)

        b = np.einsum("ji", a)
        assert_(b.base is a)
        assert_equal(b, a.T)

        b = np.einsum(a, [1, 0])
        assert_(b.base is a)
        assert_equal(b, a.T)

        # diagonal
        a = np.arange(9)
        a.shape = (3, 3)

        b = np.einsum("ii->i", a)
        assert_(b.base is a)
        assert_equal(b, [a[i, i] for i in range(3)])

        b = np.einsum(a, [0, 0], [0])
        assert_(b.base is a)
        assert_equal(b, [a[i, i] for i in range(3)])

        # diagonal with various ways of broadcasting an additional dimension
        a = np.arange(27)
        a.shape = (3, 3, 3)

        b = np.einsum("...ii->...i", a)
        assert_(b.base is a)
        assert_equal(b, [[x[i, i] for i in range(3)] for x in a])

        b = np.einsum(a, [Ellipsis, 0, 0], [Ellipsis, 0])
        assert_(b.base is a)
        assert_equal(b, [[x[i, i] for i in range(3)] for x in a])

        b = np.einsum("ii...->...i", a)
        assert_(b.base is a)
        assert_equal(b, [[x[i, i] for i in range(3)]
                         for x in a.transpose(2, 0, 1)])

        b = np.einsum(a, [0, 0, Ellipsis], [Ellipsis, 0])
        assert_(b.base is a)
        assert_equal(b, [[x[i, i] for i in range(3)]
                         for x in a.transpose(2, 0, 1)])

        b = np.einsum("...ii->i...", a)
        assert_(b.base is a)
        assert_equal(b, [a[:, i, i] for i in range(3)])

        b = np.einsum(a, [Ellipsis, 0, 0], [0, Ellipsis])
        assert_(b.base is a)
        assert_equal(b, [a[:, i, i] for i in range(3)])

        b = np.einsum("jii->ij", a)
        assert_(b.base is a)
        assert_equal(b, [a[:, i, i] for i in range(3)])

        b = np.einsum(a, [1, 0, 0], [0, 1])
        assert_(b.base is a)
        assert_equal(b, [a[:, i, i] for i in range(3)])

        b = np.einsum("ii...->i...", a)
        assert_(b.base is a)
        assert_equal(b, [a.transpose(2, 0, 1)[:, i, i] for i in range(3)])

        b = np.einsum(a, [0, 0, Ellipsis], [0, Ellipsis])
        assert_(b.base is a)
        assert_equal(b, [a.transpose(2, 0, 1)[:, i, i] for i in range(3)])

        b = np.einsum("i...i->i...", a)
        assert_(b.base is a)
        assert_equal(b, [a.transpose(1, 0, 2)[:, i, i] for i in range(3)])

        b = np.einsum(a, [0, Ellipsis, 0], [0, Ellipsis])
        assert_(b.base is a)
        assert_equal(b, [a.transpose(1, 0, 2)[:, i, i] for i in range(3)])

        b = np.einsum("i...i->...i", a)
        assert_(b.base is a)
        assert_equal(b, [[x[i, i] for i in range(3)]
                         for x in a.transpose(1, 0, 2)])

        b = np.einsum(a, [0, Ellipsis, 0], [Ellipsis, 0])
        assert_(b.base is a)
        assert_equal(b, [[x[i, i] for i in range(3)]
                         for x in a.transpose(1, 0, 2)])

        # triple diagonal
        a = np.arange(27)
        a.shape = (3, 3, 3)

        b = np.einsum("iii->i", a)
        assert_(b.base is a)
        assert_equal(b, [a[i, i, i] for i in range(3)])

        b = np.einsum(a, [0, 0, 0], [0])
        assert_(b.base is a)
        assert_equal(b, [a[i, i, i] for i in range(3)])

        # swap axes
        a = np.arange(24)
        a.shape = (2, 3, 4)

        b = np.einsum("ijk->jik", a)
        assert_(b.base is a)
        assert_equal(b, a.swapaxes(0, 1))

        b = np.einsum(a, [0, 1, 2], [1, 0, 2])
        assert_(b.base is a)
        assert_equal(b, a.swapaxes(0, 1))

    def check_einsum_sums(self, dtype):
        # Check various sums.  Does many sizes to exercise unrolled loops.

        # sum(a, axis=-1)
        for n in range(1, 17):
            a = np.arange(n, dtype=dtype)
            assert_equal(np.einsum("i->", a), np.sum(a, axis=-1).astype(dtype))
            assert_equal(np.einsum(a, [0], []),
                         np.sum(a, axis=-1).astype(dtype))

        for n in range(1, 17):
            a = np.arange(2*3*n, dtype=dtype).reshape(2, 3, n)
            assert_equal(np.einsum("...i->...", a),
                         np.sum(a, axis=-1).astype(dtype))
            assert_equal(np.einsum(a, [Ellipsis, 0], [Ellipsis]),
                         np.sum(a, axis=-1).astype(dtype))

        # sum(a, axis=0)
        for n in range(1, 17):
            a = np.arange(2*n, dtype=dtype).reshape(2, n)
            assert_equal(np.einsum("i...->...", a),
                         np.sum(a, axis=0).astype(dtype))
            assert_equal(np.einsum(a, [0, Ellipsis], [Ellipsis]),
                         np.sum(a, axis=0).astype(dtype))

        for n in range(1, 17):
            a = np.arange(2*3*n, dtype=dtype).reshape(2, 3, n)
            assert_equal(np.einsum("i...->...", a),
                         np.sum(a, axis=0).astype(dtype))
            assert_equal(np.einsum(a, [0, Ellipsis], [Ellipsis]),
                         np.sum(a, axis=0).astype(dtype))

        # trace(a)
        for n in range(1, 17):
            a = np.arange(n*n, dtype=dtype).reshape(n, n)
            assert_equal(np.einsum("ii", a), np.trace(a).astype(dtype))
            assert_equal(np.einsum(a, [0, 0]), np.trace(a).astype(dtype))

        # multiply(a, b)
        assert_equal(np.einsum("..., ...", 3, 4), 12) # scalar case
        for n in range(1, 17):
            a = np.arange(3*n, dtype=dtype).reshape(3, n)
            b = np.arange(2*3*n, dtype=dtype).reshape(2, 3, n)
            assert_equal(np.einsum("..., ...", a, b), np.multiply(a, b))
            assert_equal(np.einsum(a, [Ellipsis], b, [Ellipsis]),
                         np.multiply(a, b))

        # inner(a,b)
        for n in range(1, 17):
            a = np.arange(2*3*n, dtype=dtype).reshape(2, 3, n)
            b = np.arange(n, dtype=dtype)
            assert_equal(np.einsum("...i, ...i", a, b), np.inner(a, b))
            assert_equal(np.einsum(a, [Ellipsis, 0], b, [Ellipsis, 0]),
                         np.inner(a, b))

        for n in range(1, 11):
            a = np.arange(n*3*2, dtype=dtype).reshape(n, 3, 2)
            b = np.arange(n, dtype=dtype)
            assert_equal(np.einsum("i..., i...", a, b), np.inner(a.T, b.T).T)
            assert_equal(np.einsum(a, [0, Ellipsis], b, [0, Ellipsis]),
                         np.inner(a.T, b.T).T)

        # outer(a,b)
        for n in range(1, 17):
            a = np.arange(3, dtype=dtype)+1
            b = np.arange(n, dtype=dtype)+1
            assert_equal(np.einsum("i,j", a, b), np.outer(a, b))
            assert_equal(np.einsum(a, [0], b, [1]), np.outer(a, b))

        # Suppress the complex warnings for the 'as f8' tests
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', np.ComplexWarning)

            # matvec(a,b) / a.dot(b) where a is matrix, b is vector
            for n in range(1, 17):
                a = np.arange(4*n, dtype=dtype).reshape(4, n)
                b = np.arange(n, dtype=dtype)
                assert_equal(np.einsum("ij, j", a, b), np.dot(a, b))
                assert_equal(np.einsum(a, [0, 1], b, [1]), np.dot(a, b))

                c = np.arange(4, dtype=dtype)
                np.einsum("ij,j", a, b, out=c,
                            dtype='f8', casting='unsafe')
                assert_equal(c,
                            np.dot(a.astype('f8'),
                                   b.astype('f8')).astype(dtype))
                c[...] = 0
                np.einsum(a, [0, 1], b, [1], out=c,
                            dtype='f8', casting='unsafe')
                assert_equal(c,
                            np.dot(a.astype('f8'),
                                   b.astype('f8')).astype(dtype))

            for n in range(1, 17):
                a = np.arange(4*n, dtype=dtype).reshape(4, n)
                b = np.arange(n, dtype=dtype)
                assert_equal(np.einsum("ji,j", a.T, b.T), np.dot(b.T, a.T))
                assert_equal(np.einsum(a.T, [1, 0], b.T, [1]), np.dot(b.T, a.T))

                c = np.arange(4, dtype=dtype)
                np.einsum("ji,j", a.T, b.T, out=c, dtype='f8', casting='unsafe')
                assert_equal(c,
                        np.dot(b.T.astype('f8'),
                               a.T.astype('f8')).astype(dtype))
                c[...] = 0
                np.einsum(a.T, [1, 0], b.T, [1], out=c,
                            dtype='f8', casting='unsafe')
                assert_equal(c,
                        np.dot(b.T.astype('f8'),
                               a.T.astype('f8')).astype(dtype))

            # matmat(a,b) / a.dot(b) where a is matrix, b is matrix
            for n in range(1, 17):
                if n < 8 or dtype != 'f2':
                    a = np.arange(4*n, dtype=dtype).reshape(4, n)
                    b = np.arange(n*6, dtype=dtype).reshape(n, 6)
                    assert_equal(np.einsum("ij,jk", a, b), np.dot(a, b))
                    assert_equal(np.einsum(a, [0, 1], b, [1, 2]), np.dot(a, b))

            for n in range(1, 17):
                a = np.arange(4*n, dtype=dtype).reshape(4, n)
                b = np.arange(n*6, dtype=dtype).reshape(n, 6)
                c = np.arange(24, dtype=dtype).reshape(4, 6)
                np.einsum("ij,jk", a, b, out=c, dtype='f8', casting='unsafe')
                assert_equal(c,
                            np.dot(a.astype('f8'),
                                   b.astype('f8')).astype(dtype))
                c[...] = 0
                np.einsum(a, [0, 1], b, [1, 2], out=c,
                                dtype='f8', casting='unsafe')
                assert_equal(c,
                            np.dot(a.astype('f8'),
                                   b.astype('f8')).astype(dtype))

            # matrix triple product (note this is not currently an efficient
            # way to multiply 3 matrices)
            a = np.arange(12, dtype=dtype).reshape(3, 4)
            b = np.arange(20, dtype=dtype).reshape(4, 5)
            c = np.arange(30, dtype=dtype).reshape(5, 6)
            if dtype != 'f2':
                assert_equal(np.einsum("ij,jk,kl", a, b, c),
                                    a.dot(b).dot(c))
                assert_equal(np.einsum(a, [0, 1], b, [1, 2], c, [2, 3]),
                                    a.dot(b).dot(c))

            d = np.arange(18, dtype=dtype).reshape(3, 6)
            np.einsum("ij,jk,kl", a, b, c, out=d,
                                dtype='f8', casting='unsafe')
            assert_equal(d, a.astype('f8').dot(b.astype('f8')
                        ).dot(c.astype('f8')).astype(dtype))
            d[...] = 0
            np.einsum(a, [0, 1], b, [1, 2], c, [2, 3], out=d,
                                dtype='f8', casting='unsafe')
            assert_equal(d, a.astype('f8').dot(b.astype('f8')
                        ).dot(c.astype('f8')).astype(dtype))

            # tensordot(a, b)
            if np.dtype(dtype) != np.dtype('f2'):
                a = np.arange(60, dtype=dtype).reshape(3, 4, 5)
                b = np.arange(24, dtype=dtype).reshape(4, 3, 2)
                assert_equal(np.einsum("ijk, jil -> kl", a, b),
                                np.tensordot(a, b, axes=([1, 0], [0, 1])))
                assert_equal(np.einsum(a, [0, 1, 2], b, [1, 0, 3], [2, 3]),
                                np.tensordot(a, b, axes=([1, 0], [0, 1])))

                c = np.arange(10, dtype=dtype).reshape(5, 2)
                np.einsum("ijk,jil->kl", a, b, out=c,
                                        dtype='f8', casting='unsafe')
                assert_equal(c, np.tensordot(a.astype('f8'), b.astype('f8'),
                                        axes=([1, 0], [0, 1])).astype(dtype))
                c[...] = 0
                np.einsum(a, [0, 1, 2], b, [1, 0, 3], [2, 3], out=c,
                                        dtype='f8', casting='unsafe')
                assert_equal(c, np.tensordot(a.astype('f8'), b.astype('f8'),
                                        axes=([1, 0], [0, 1])).astype(dtype))

        # logical_and(logical_and(a!=0, b!=0), c!=0)
        a = np.array([1,   3,   -2,   0,   12,  13,   0,   1], dtype=dtype)
        b = np.array([0,   3.5, 0.,   -2,  0,   1,    3,   12], dtype=dtype)
        c = np.array([True, True, False, True, True, False, True, True])
        assert_equal(np.einsum("i,i,i->i", a, b, c,
                                dtype='?', casting='unsafe'),
                            np.logical_and(np.logical_and(a!=0, b!=0), c!=0))
        assert_equal(np.einsum(a, [0], b, [0], c, [0], [0],
                                dtype='?', casting='unsafe'),
                            np.logical_and(np.logical_and(a!=0, b!=0), c!=0))

        a = np.arange(9, dtype=dtype)
        assert_equal(np.einsum(",i->", 3, a), 3*np.sum(a))
        assert_equal(np.einsum(3, [], a, [0], []), 3*np.sum(a))
        assert_equal(np.einsum("i,->", a, 3), 3*np.sum(a))
        assert_equal(np.einsum(a, [0], 3, [], []), 3*np.sum(a))

        # Various stride0, contiguous, and SSE aligned variants
        for n in range(1, 25):
            a = np.arange(n, dtype=dtype)
            if np.dtype(dtype).itemsize > 1:
                assert_equal(np.einsum("...,...", a, a), np.multiply(a, a))
                assert_equal(np.einsum("i,i", a, a), np.dot(a, a))
                assert_equal(np.einsum("i,->i", a, 2), 2*a)
                assert_equal(np.einsum(",i->i", 2, a), 2*a)
                assert_equal(np.einsum("i,->", a, 2), 2*np.sum(a))
                assert_equal(np.einsum(",i->", 2, a), 2*np.sum(a))

                assert_equal(np.einsum("...,...", a[1:], a[:-1]),
                             np.multiply(a[1:], a[:-1]))
                assert_equal(np.einsum("i,i", a[1:], a[:-1]),
                             np.dot(a[1:], a[:-1]))
                assert_equal(np.einsum("i,->i", a[1:], 2), 2*a[1:])
                assert_equal(np.einsum(",i->i", 2, a[1:]), 2*a[1:])
                assert_equal(np.einsum("i,->", a[1:], 2), 2*np.sum(a[1:]))
                assert_equal(np.einsum(",i->", 2, a[1:]), 2*np.sum(a[1:]))

        # An object array, summed as the data type
        a = np.arange(9, dtype=object)

        b = np.einsum("i->", a, dtype=dtype, casting='unsafe')
        assert_equal(b, np.sum(a))
        assert_equal(b.dtype, np.dtype(dtype))

        b = np.einsum(a, [0], [], dtype=dtype, casting='unsafe')
        assert_equal(b, np.sum(a))
        assert_equal(b.dtype, np.dtype(dtype))

        # A case which was failing (ticket #1885)
        p = np.arange(2) + 1
        q = np.arange(4).reshape(2, 2) + 3
        r = np.arange(4).reshape(2, 2) + 7
        assert_equal(np.einsum('z,mz,zm->', p, q, r), 253)

    def test_einsum_sums_int8(self):
        self.check_einsum_sums('i1');

    def test_einsum_sums_uint8(self):
        self.check_einsum_sums('u1');

    def test_einsum_sums_int16(self):
        self.check_einsum_sums('i2');

    def test_einsum_sums_uint16(self):
        self.check_einsum_sums('u2');

    def test_einsum_sums_int32(self):
        self.check_einsum_sums('i4');

    def test_einsum_sums_uint32(self):
        self.check_einsum_sums('u4');

    def test_einsum_sums_int64(self):
        self.check_einsum_sums('i8');

    def test_einsum_sums_uint64(self):
        self.check_einsum_sums('u8');

    def test_einsum_sums_float16(self):
        self.check_einsum_sums('f2');

    def test_einsum_sums_float32(self):
        self.check_einsum_sums('f4');

    def test_einsum_sums_float64(self):
        self.check_einsum_sums('f8');

    def test_einsum_sums_longdouble(self):
        self.check_einsum_sums(np.longdouble);

    def test_einsum_sums_cfloat64(self):
        self.check_einsum_sums('c8');

    def test_einsum_sums_cfloat128(self):
        self.check_einsum_sums('c16');

    def test_einsum_sums_clongdouble(self):
        self.check_einsum_sums(np.clongdouble);

    def test_einsum_misc(self):
        # This call used to crash because of a bug in
        # PyArray_AssignZero
        a = np.ones((1, 2))
        b = np.ones((2, 2, 1))
        assert_equal(np.einsum('ij...,j...->i...', a, b), [[[2], [2]]])

        # The iterator had an issue with buffering this reduction
        a = np.ones((5, 12, 4, 2, 3), np.int64)
        b = np.ones((5, 12, 11), np.int64)
        assert_equal(np.einsum('ijklm,ijn,ijn->', a, b, b),
                        np.einsum('ijklm,ijn->', a, b))

        # Issue #2027, was a problem in the contiguous 3-argument
        # inner loop implementation
        a = np.arange(1, 3)
        b = np.arange(1, 5).reshape(2, 2)
        c = np.arange(1, 9).reshape(4, 2)
        assert_equal(np.einsum('x,yx,zx->xzy', a, b, c),
                    [[[1,  3], [3,  9], [5, 15], [7, 21]],
                    [[8, 16], [16, 32], [24, 48], [32, 64]]])

    def test_einsum_broadcast(self):
        # Issue #2455 change in handling ellipsis
        # remove the 'middle broadcast' error
        # only use the 'RIGHT' iteration in prepare_op_axes
        # adds auto broadcast on left where it belongs
        # broadcast on right has to be explicit

        A = np.arange(2*3*4).reshape(2,3,4)
        B = np.arange(3)
        ref = np.einsum('ijk,j->ijk',A, B)
        assert_equal(np.einsum('ij...,j...->ij...',A, B), ref)
        assert_equal(np.einsum('ij...,...j->ij...',A, B), ref)
        assert_equal(np.einsum('ij...,j->ij...',A, B), ref) # used to raise error

        A = np.arange(12).reshape((4,3))
        B = np.arange(6).reshape((3,2))
        ref = np.einsum('ik,kj->ij', A, B)
        assert_equal(np.einsum('ik...,k...->i...', A, B), ref)
        assert_equal(np.einsum('ik...,...kj->i...j', A, B), ref)
        assert_equal(np.einsum('...k,kj', A, B), ref) # used to raise error
        assert_equal(np.einsum('ik,k...->i...', A, B), ref) # used to raise error

        dims=[2,3,4,5];
        a = np.arange(np.prod(dims)).reshape(dims)
        v = np.arange(dims[2])
        ref = np.einsum('ijkl,k->ijl', a, v)
        assert_equal(np.einsum('ijkl,k', a, v), ref)
        assert_equal(np.einsum('...kl,k', a, v), ref)  # used to raise error
        assert_equal(np.einsum('...kl,k...', a, v), ref)
        # no real diff from 1st

        J,K,M=160,160,120;
        A=np.arange(J*K*M).reshape(1,1,1,J,K,M)
        B=np.arange(J*K*M*3).reshape(J,K,M,3)
        ref = np.einsum('...lmn,...lmno->...o', A, B)
        assert_equal(np.einsum('...lmn,lmno->...o', A, B), ref)  # used to raise error

    def test_einsum_fixedstridebug(self):
        # Issue #4485 obscure einsum bug
        # This case revealed a bug in nditer where it reported a stride
        # as 'fixed' (0) when it was in fact not fixed during processing
        # (0 or 4). The reason for the bug was that the check for a fixed
        # stride was using the information from the 2D inner loop reuse
        # to restrict the iteration dimensions it had to validate to be
        # the same, but that 2D inner loop reuse logic is only triggered
        # during the buffer copying step, and hence it was invalid to
        # rely on those values. The fix is to check all the dimensions
        # of the stride in question, which in the test case reveals that
        # the stride is not fixed.
        #
        # NOTE: This test is triggered by the fact that the default buffersize,
        #       used by einsum, is 8192, and 3*2731 = 8193, is larger than that
        #       and results in a mismatch between the buffering and the
        #       striding for operand A.
        A = np.arange(2*3).reshape(2,3).astype(np.float32)
        B = np.arange(2*3*2731).reshape(2,3,2731).astype(np.int16)
        es = np.einsum('cl,cpx->lpx', A, B)
        tp = np.tensordot(A, B, axes=(0, 0))
        assert_equal(es, tp)
        # The following is the original test case from the bug report,
        # made repeatable by changing random arrays to aranges.
        A = np.arange(3*3).reshape(3,3).astype(np.float64)
        B = np.arange(3*3*64*64).reshape(3,3,64,64).astype(np.float32)
        es = np.einsum ('cl,cpxy->lpxy', A,B)
        tp = np.tensordot(A,B, axes=(0,0))
        assert_equal(es, tp)


if __name__ == "__main__":
    run_module_suite()
