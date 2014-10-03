#!/usr/bin/env python
"""Test functions for fftpack.helper module

Copied from fftpack.helper by Pearu Peterson, October 2005

"""
from __future__ import division, absolute_import, print_function

import numpy as np
from numpy.testing import TestCase, run_module_suite, assert_array_almost_equal
from numpy import fft
from numpy import pi


class TestFFTShift(TestCase):

    def test_definition(self):
        x = [0, 1, 2, 3, 4, -4, -3, -2, -1]
        y = [-4, -3, -2, -1, 0, 1, 2, 3, 4]
        assert_array_almost_equal(fft.fftshift(x), y)
        assert_array_almost_equal(fft.ifftshift(y), x)
        x = [0, 1, 2, 3, 4, -5, -4, -3, -2, -1]
        y = [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4]
        assert_array_almost_equal(fft.fftshift(x), y)
        assert_array_almost_equal(fft.ifftshift(y), x)

    def test_inverse(self):
        for n in [1, 4, 9, 100, 211]:
            x = np.random.random((n,))
            assert_array_almost_equal(fft.ifftshift(fft.fftshift(x)), x)

    def test_axes_keyword(self):
        freqs = [[ 0,  1,  2], [ 3,  4, -4], [-3, -2, -1]]
        shifted = [[-1, -3, -2], [ 2,  0,  1], [-4,  3,  4]]
        assert_array_almost_equal(fft.fftshift(freqs, axes=(0, 1)), shifted)
        assert_array_almost_equal(fft.fftshift(freqs, axes=0),
                fft.fftshift(freqs, axes=(0,)))
        assert_array_almost_equal(fft.ifftshift(shifted, axes=(0, 1)), freqs)
        assert_array_almost_equal(fft.ifftshift(shifted, axes=0),
                fft.ifftshift(shifted, axes=(0,)))


class TestFFTFreq(TestCase):

    def test_definition(self):
        x = [0, 1, 2, 3, 4, -4, -3, -2, -1]
        assert_array_almost_equal(9*fft.fftfreq(9), x)
        assert_array_almost_equal(9*pi*fft.fftfreq(9, pi), x)
        x = [0, 1, 2, 3, 4, -5, -4, -3, -2, -1]
        assert_array_almost_equal(10*fft.fftfreq(10), x)
        assert_array_almost_equal(10*pi*fft.fftfreq(10, pi), x)


class TestRFFTFreq(TestCase):

    def test_definition(self):
        x = [0, 1, 2, 3, 4]
        assert_array_almost_equal(9*fft.rfftfreq(9), x)
        assert_array_almost_equal(9*pi*fft.rfftfreq(9, pi), x)
        x = [0, 1, 2, 3, 4, 5]
        assert_array_almost_equal(10*fft.rfftfreq(10), x)
        assert_array_almost_equal(10*pi*fft.rfftfreq(10, pi), x)


class TestIRFFTN(TestCase):

    def test_not_last_axis_success(self):
        ar, ai = np.random.random((2, 16, 8, 32))
        a = ar + 1j*ai

        axes = (-2,)

        # Should not raise error
        fft.irfftn(a, axes=axes)


if __name__ == "__main__":
    run_module_suite()
