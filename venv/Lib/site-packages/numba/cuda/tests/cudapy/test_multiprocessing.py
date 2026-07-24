import os
import multiprocessing as mp

import numpy as np

from numba import cuda
from numba.cuda.testing import skip_on_cudasim, CUDATestCase
import unittest

has_mp_get_context = hasattr(mp, 'get_context')
is_unix = os.name == 'posix'


def fork_test(q):
    from numba.cuda.cudadrv.error import CudaDriverError
    try:
        cuda.to_device(np.arange(1))
    except CudaDriverError as e:
        q.put(e)
    else:
        q.put(None)


@skip_on_cudasim('disabled for cudasim')
class TestMultiprocessing(CUDATestCase):
    @unittest.skipUnless(has_mp_get_context, 'requires mp.get_context')
    @unittest.skipUnless(is_unix, 'requires Unix')
    def test_fork(self):
        """
        Test fork detection.
        """
        cuda.current_context()  # force cuda initialize
        # fork in process that also uses CUDA
        ctx = mp.get_context('fork')
        q = ctx.Queue()
        proc = ctx.Process(target=fork_test, args=[q])
        proc.start()
        exc = q.get()
        proc.join()
        # there should be an exception raised in the child process
        self.assertIsNotNone(exc)
        self.assertIn('CUDA initialized before forking', str(exc))


if __name__ == '__main__':
    unittest.main()
