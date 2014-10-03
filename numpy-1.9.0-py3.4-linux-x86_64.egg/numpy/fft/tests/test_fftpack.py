from __future__ import division, absolute_import, print_function

import numpy as np
from numpy.testing import TestCase, run_module_suite, assert_array_almost_equal
from numpy.testing import assert_array_equal
import threading
import sys
if sys.version_info[0] >= 3:
    import queue
else:
    import Queue as queue


def fft1(x):
    L = len(x)
    phase = -2j*np.pi*(np.arange(L)/float(L))
    phase = np.arange(L).reshape(-1, 1) * phase
    return np.sum(x*np.exp(phase), axis=1)


class TestFFTShift(TestCase):

    def test_fft_n(self):
        self.assertRaises(ValueError, np.fft.fft, [1, 2, 3], 0)


class TestFFT1D(TestCase):

    def test_basic(self):
        rand = np.random.random
        x = rand(30) + 1j*rand(30)
        assert_array_almost_equal(fft1(x), np.fft.fft(x))


class TestFFTThreadSafe(TestCase):
    threads = 16
    input_shape = (800, 200)

    def _test_mtsame(self, func, *args):
        def worker(args, q):
            q.put(func(*args))

        q = queue.Queue()
        expected = func(*args)

        # Spin off a bunch of threads to call the same function simultaneously
        t = [threading.Thread(target=worker, args=(args, q))
             for i in range(self.threads)]
        [x.start() for x in t]

        [x.join() for x in t]
        # Make sure all threads returned the correct value
        for i in range(self.threads):
            assert_array_equal(q.get(timeout=5), expected,
                'Function returned wrong value in multithreaded context')

    def test_fft(self):
        a = np.ones(self.input_shape) * 1+0j
        self._test_mtsame(np.fft.fft, a)

    def test_ifft(self):
        a = np.ones(self.input_shape) * 1+0j
        self._test_mtsame(np.fft.ifft, a)

    def test_rfft(self):
        a = np.ones(self.input_shape)
        self._test_mtsame(np.fft.rfft, a)

    def test_irfft(self):
        a = np.ones(self.input_shape) * 1+0j
        self._test_mtsame(np.fft.irfft, a)


if __name__ == "__main__":
    run_module_suite()
