import numpy as np

from numba import jit
from numba.core import types
from numba.tests.support import TestCase, tag
import unittest


# Array overlaps involving a displacement

def array_overlap1(src, dest, k=1):
    assert src.shape == dest.shape
    dest[k:] = src[:-k]

def array_overlap2(src, dest, k=1):
    assert src.shape == dest.shape
    dest[:-k] = src[k:]

def array_overlap3(src, dest, k=1):
    assert src.shape == dest.shape
    dest[:,:-k] = src[:,k:]

def array_overlap4(src, dest, k=1):
    assert src.shape == dest.shape
    dest[:,k:] = src[:,:-k]

def array_overlap5(src, dest, k=1):
    assert src.shape == dest.shape
    dest[...,:-k] = src[...,k:]

def array_overlap6(src, dest, k=1):
    assert src.shape == dest.shape
    dest[...,k:] = src[...,:-k]

# Array overlaps involving an in-place reversal

def array_overlap11(src, dest):
    assert src.shape == dest.shape
    dest[::-1] = src

def array_overlap12(src, dest):
    assert src.shape == dest.shape
    dest[:] = src[::-1]

def array_overlap13(src, dest):
    assert src.shape == dest.shape
    dest[:,::-1] = src

def array_overlap14(src, dest):
    assert src.shape == dest.shape
    dest[:] = src[:,::-1]

def array_overlap15(src, dest):
    assert src.shape == dest.shape
    dest[...,::-1] = src

def array_overlap16(src, dest):
    assert src.shape == dest.shape
    dest[:] = src[...,::-1]


class TestArrayOverlap(TestCase):

    def check_overlap(self, pyfunc, min_ndim, have_k_argument=False):
        N = 4

        def vary_layouts(orig):
            yield orig.copy(order='C')
            yield orig.copy(order='F')
            a = orig[::-1].copy()[::-1]
            assert not a.flags.c_contiguous and not a.flags.f_contiguous
            yield a

        def check(pyfunc, cfunc, pydest, cdest, kwargs):
            pyfunc(pydest, pydest, **kwargs)
            cfunc(cdest, cdest, **kwargs)
            self.assertPreciseEqual(pydest, cdest)

        cfunc = jit(nopython=True)(pyfunc)
        # Check for up to 3d arrays
        for ndim in range(min_ndim, 4):
            shape = (N,) * ndim
            orig = np.arange(0, N**ndim).reshape(shape)
            # Note we cannot copy a 'A' layout array exactly (bitwise),
            # so instead we call vary_layouts() twice
            for pydest, cdest in zip(vary_layouts(orig), vary_layouts(orig)):
                if have_k_argument:
                    for k in range(1, N):
                        check(pyfunc, cfunc, pydest, cdest, dict(k=k))
                else:
                    check(pyfunc, cfunc, pydest, cdest, {})

    def check_overlap_with_k(self, pyfunc, min_ndim):
        self.check_overlap(pyfunc, min_ndim=min_ndim, have_k_argument=True)

    def test_overlap1(self):
        self.check_overlap_with_k(array_overlap1, min_ndim=1)

    def test_overlap2(self):
        self.check_overlap_with_k(array_overlap2, min_ndim=1)

    def test_overlap3(self):
        self.check_overlap_with_k(array_overlap3, min_ndim=2)

    def test_overlap4(self):
        self.check_overlap_with_k(array_overlap4, min_ndim=2)

    def test_overlap5(self):
        self.check_overlap_with_k(array_overlap5, min_ndim=1)

    def test_overlap6(self):
        self.check_overlap_with_k(array_overlap6, min_ndim=1)

    def test_overlap11(self):
        self.check_overlap(array_overlap11, min_ndim=1)

    def test_overlap12(self):
        self.check_overlap(array_overlap12, min_ndim=1)

    def test_overlap13(self):
        self.check_overlap(array_overlap13, min_ndim=2)

    def test_overlap14(self):
        self.check_overlap(array_overlap14, min_ndim=2)

    def test_overlap15(self):
        self.check_overlap(array_overlap15, min_ndim=1)

    def test_overlap16(self):
        self.check_overlap(array_overlap16, min_ndim=1)


if __name__ == '__main__':
    unittest.main()
