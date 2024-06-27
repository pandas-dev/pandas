import numpy as np
import platform

from numba import cuda
from numba.cuda.testing import unittest, ContextResettingTestCase


class TestPinned(ContextResettingTestCase):

    def _run_copies(self, A):
        A0 = np.copy(A)

        stream = cuda.stream()
        ptr = cuda.to_device(A, copy=False, stream=stream)
        ptr.copy_to_device(A, stream=stream)
        ptr.copy_to_host(A, stream=stream)
        stream.synchronize()

        self.assertTrue(np.allclose(A, A0))

    def test_pinned(self):
        machine = platform.machine()
        if machine.startswith('arm') or machine.startswith('aarch64'):
            count = 262144   # 2MB
        else:
            count = 2097152  # 16MB
        A = np.arange(count)
        with cuda.pinned(A):
            self._run_copies(A)

    def test_unpinned(self):
        A = np.arange(2 * 1024 * 1024) # 16 MB
        self._run_copies(A)


if __name__ == '__main__':
    unittest.main()
