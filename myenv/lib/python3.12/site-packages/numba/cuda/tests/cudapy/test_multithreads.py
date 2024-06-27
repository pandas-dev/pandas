import traceback
import threading
import multiprocessing
import numpy as np
from numba import cuda
from numba.cuda.testing import (skip_on_cudasim, skip_under_cuda_memcheck,
                                CUDATestCase)
import unittest

try:
    from concurrent.futures import ThreadPoolExecutor
except ImportError:
    has_concurrent_futures = False
else:
    has_concurrent_futures = True


has_mp_get_context = hasattr(multiprocessing, 'get_context')


def check_concurrent_compiling():
    @cuda.jit
    def foo(x):
        x[0] += 1

    def use_foo(x):
        foo[1, 1](x)
        return x

    arrays = [cuda.to_device(np.arange(10)) for i in range(10)]
    expected = np.arange(10)
    expected[0] += 1
    with ThreadPoolExecutor(max_workers=4) as e:
        for ary in e.map(use_foo, arrays):
            np.testing.assert_equal(ary, expected)


def spawn_process_entry(q):
    try:
        check_concurrent_compiling()
    # Catch anything that goes wrong in the threads
    except:  # noqa: E722
        msg = traceback.format_exc()
        q.put('\n'.join(['', '=' * 80, msg]))
    else:
        q.put(None)


@skip_under_cuda_memcheck('Hangs cuda-memcheck')
@skip_on_cudasim('disabled for cudasim')
class TestMultiThreadCompiling(CUDATestCase):

    @unittest.skipIf(not has_concurrent_futures, "no concurrent.futures")
    def test_concurrent_compiling(self):
        check_concurrent_compiling()

    @unittest.skipIf(not has_mp_get_context, "no multiprocessing.get_context")
    def test_spawn_concurrent_compilation(self):
        # force CUDA context init
        cuda.get_current_device()
        # use "spawn" to avoid inheriting the CUDA context
        ctx = multiprocessing.get_context('spawn')

        q = ctx.Queue()
        p = ctx.Process(target=spawn_process_entry, args=(q,))
        p.start()
        try:
            err = q.get()
        finally:
            p.join()
        if err is not None:
            raise AssertionError(err)
        self.assertEqual(p.exitcode, 0, 'test failed in child process')

    def test_invalid_context_error_with_d2h(self):
        def d2h(arr, out):
            out[:] = arr.copy_to_host()

        arr = np.arange(1, 4)
        out = np.zeros_like(arr)
        darr = cuda.to_device(arr)
        th = threading.Thread(target=d2h, args=[darr, out])
        th.start()
        th.join()
        np.testing.assert_equal(arr, out)

    def test_invalid_context_error_with_d2d(self):
        def d2d(dst, src):
            dst.copy_to_device(src)

        arr = np.arange(100)
        common = cuda.to_device(arr)
        darr = cuda.to_device(np.zeros(common.shape, dtype=common.dtype))
        th = threading.Thread(target=d2d, args=[darr, common])
        th.start()
        th.join()
        np.testing.assert_equal(darr.copy_to_host(), arr)


if __name__ == '__main__':
    unittest.main()
