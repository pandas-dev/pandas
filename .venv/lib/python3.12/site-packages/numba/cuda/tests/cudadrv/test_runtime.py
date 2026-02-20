import multiprocessing
import os
from numba.core import config
from numba.cuda.cudadrv.runtime import runtime
from numba.cuda.testing import unittest, SerialMixin, skip_on_cudasim
from unittest.mock import patch


def set_visible_devices_and_check(q):
    try:
        from numba import cuda
        import os

        os.environ['CUDA_VISIBLE_DEVICES'] = '0'
        q.put(len(cuda.gpus.lst))
    except: # noqa: E722
        # Sentinel value for error executing test code
        q.put(-1)


if config.ENABLE_CUDASIM:
    SUPPORTED_VERSIONS = (-1, -1),
else:
    SUPPORTED_VERSIONS = ((11, 0), (11, 1), (11, 2), (11, 3), (11, 4), (11, 5),
                          (11, 6), (11, 7))


class TestRuntime(unittest.TestCase):
    def test_is_supported_version_true(self):
        for v in SUPPORTED_VERSIONS:
            with patch.object(runtime, 'get_version', return_value=v):
                self.assertTrue(runtime.is_supported_version())

    @skip_on_cudasim('The simulator always simulates a supported runtime')
    def test_is_supported_version_false(self):
        # Check with an old unsupported version and some potential future
        # versions
        for v in ((10, 2), (11, 8), (12, 0)):
            with patch.object(runtime, 'get_version', return_value=v):
                self.assertFalse(runtime.is_supported_version())

    def test_supported_versions(self):
        self.assertEqual(SUPPORTED_VERSIONS, runtime.supported_versions)


class TestVisibleDevices(unittest.TestCase, SerialMixin):
    def test_visible_devices_set_after_import(self):
        # See Issue #6149. This test checks that we can set
        # CUDA_VISIBLE_DEVICES after importing Numba and have the value
        # reflected in the available list of GPUs. Prior to the fix for this
        # issue, Numba made a call to runtime.get_version() on import that
        # initialized the driver and froze the list of available devices before
        # CUDA_VISIBLE_DEVICES could be set by the user.

        # Avoid importing cuda at the top level so that
        # set_visible_devices_and_check gets to import it first in its process
        from numba import cuda

        if len(cuda.gpus.lst) in (0, 1):
            self.skipTest('This test requires multiple GPUs')

        if os.environ.get('CUDA_VISIBLE_DEVICES'):
            msg = 'Cannot test when CUDA_VISIBLE_DEVICES already set'
            self.skipTest(msg)

        ctx = multiprocessing.get_context('spawn')
        q = ctx.Queue()
        p = ctx.Process(target=set_visible_devices_and_check, args=(q,))
        p.start()
        try:
            visible_gpu_count = q.get()
        finally:
            p.join()

        # Make an obvious distinction between an error running the test code
        # and an incorrect number of GPUs in the list
        msg = 'Error running set_visible_devices_and_check'
        self.assertNotEqual(visible_gpu_count, -1, msg=msg)

        # The actual check that we see only one GPU
        self.assertEqual(visible_gpu_count, 1)


if __name__ == '__main__':
    unittest.main()
