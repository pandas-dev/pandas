import multiprocessing as mp
import logging
import traceback
from numba.cuda.testing import unittest, CUDATestCase
from numba.cuda.testing import (skip_on_cudasim, skip_with_cuda_python,
                                skip_under_cuda_memcheck)
from numba.tests.support import linux_only


def child_test():
    from numba import cuda, int32, void
    from numba.core import config
    import io
    import numpy as np
    import threading

    # Enable PTDS before we make any CUDA driver calls.  Enabling it first
    # ensures that PTDS APIs are used because the CUDA driver looks up API
    # functions on first use and memoizes them.
    config.CUDA_PER_THREAD_DEFAULT_STREAM = 1

    # Set up log capture for the Driver API so we can see what API calls were
    # used.
    logbuf = io.StringIO()
    handler = logging.StreamHandler(logbuf)
    cudadrv_logger = logging.getLogger('numba.cuda.cudadrv.driver')
    cudadrv_logger.addHandler(handler)
    cudadrv_logger.setLevel(logging.DEBUG)

    # Set up data for our test, and copy over to the device
    N = 2 ** 16
    N_THREADS = 10
    N_ADDITIONS = 4096

    # Seed the RNG for repeatability
    np.random.seed(1)
    x = np.random.randint(low=0, high=1000, size=N, dtype=np.int32)
    r = np.zeros_like(x)

    # One input and output array for each thread
    xs = [cuda.to_device(x) for _ in range(N_THREADS)]
    rs = [cuda.to_device(r) for _ in range(N_THREADS)]

    # Compute the grid size and get the [per-thread] default stream
    n_threads = 256
    n_blocks = N // n_threads
    stream = cuda.default_stream()

    # A simple multiplication-by-addition kernel. What it does exactly is not
    # too important; only that we have a kernel that does something.
    @cuda.jit(void(int32[::1], int32[::1]))
    def f(r, x):
        i = cuda.grid(1)

        if i > len(r):
            return

        # Accumulate x into r
        for j in range(N_ADDITIONS):
            r[i] += x[i]

    # This function will be used to launch the kernel from each thread on its
    # own unique data.
    def kernel_thread(n):
        f[n_blocks, n_threads, stream](rs[n], xs[n])

    # Create threads
    threads = [threading.Thread(target=kernel_thread, args=(i,))
               for i in range(N_THREADS)]

    # Start all threads
    for thread in threads:
        thread.start()

    # Wait for all threads to finish, to ensure that we don't synchronize with
    # the device until all kernels are scheduled.
    for thread in threads:
        thread.join()

    # Synchronize with the device
    cuda.synchronize()

    # Check output is as expected
    expected = x * N_ADDITIONS
    for i in range(N_THREADS):
        np.testing.assert_equal(rs[i].copy_to_host(), expected)

    # Return the driver log output to the calling process for checking
    handler.flush()
    return logbuf.getvalue()


def child_test_wrapper(result_queue):
    try:
        output = child_test()
        success = True
    # Catch anything raised so it can be propagated
    except: # noqa: E722
        output = traceback.format_exc()
        success = False

    result_queue.put((success, output))


# Run on Linux only until the reason for test hangs on Windows (Issue #8635,
# https://github.com/numba/numba/issues/8635) is diagnosed
@linux_only
@skip_under_cuda_memcheck('Hangs cuda-memcheck')
@skip_on_cudasim('Streams not supported on the simulator')
class TestPTDS(CUDATestCase):
    @skip_with_cuda_python('Function names unchanged for PTDS with NV Binding')
    def test_ptds(self):
        # Run a test with PTDS enabled in a child process
        ctx = mp.get_context('spawn')
        result_queue = ctx.Queue()
        proc = ctx.Process(target=child_test_wrapper, args=(result_queue,))
        proc.start()
        proc.join()
        success, output = result_queue.get()

        # Ensure the child process ran to completion before checking its output
        if not success:
            self.fail(output)

        # Functions with a per-thread default stream variant that we expect to
        # see in the output
        ptds_functions = ('cuMemcpyHtoD_v2_ptds', 'cuLaunchKernel_ptsz',
                          'cuMemcpyDtoH_v2_ptds')

        for fn in ptds_functions:
            with self.subTest(fn=fn, expected=True):
                self.assertIn(fn, output)

        # Non-PTDS versions of the functions that we should not see in the
        # output:
        legacy_functions = ('cuMemcpyHtoD_v2', 'cuLaunchKernel',
                            'cuMemcpyDtoH_v2')

        for fn in legacy_functions:
            with self.subTest(fn=fn, expected=False):
                # Ensure we only spot these function names appearing without a
                # _ptds or _ptsz suffix by checking including the end of the
                # line in the log
                fn_at_end = f'{fn}\n'
                self.assertNotIn(fn_at_end, output)


if __name__ == '__main__':
    unittest.main()
