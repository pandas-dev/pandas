import time
import ctypes

import numpy as np

from numba.tests.support import captured_stdout
from numba import vectorize, guvectorize
import unittest


class TestParUfuncIssues(unittest.TestCase):

    _numba_parallel_test_ = False

    def test_thread_response(self):
        """
        Related to #89.
        This does not test #89 but tests the fix for it.
        We want to make sure the worker threads can be used multiple times
        and with different time gap between each execution.
        """

        @vectorize('float64(float64, float64)', target='parallel')
        def fnv(a, b):
            return a + b

        sleep_time = 1   # 1 second
        while sleep_time > 0.00001:    # 10us
            time.sleep(sleep_time)
            a = b = np.arange(10**5)
            np.testing.assert_equal(a + b, fnv(a, b))
            # Reduce sleep time
            sleep_time /= 2

    def test_gil_reacquire_deadlock(self):
        """
        Testing issue #1998 due to GIL reacquiring
        """
        # make a ctypes callback that requires the GIL
        proto = ctypes.CFUNCTYPE(None, ctypes.c_int32)
        characters = 'abcdefghij'

        def bar(x):
            print(characters[x])

        cbar = proto(bar)

        # our unit under test
        @vectorize(['int32(int32)'], target='parallel', nopython=True)
        def foo(x):
            print(x % 10)  # this reacquires the GIL
            cbar(x % 10)   # this reacquires the GIL
            return x * 2

        # Numpy ufunc has a heuristic to determine whether to release the GIL
        # during execution.  Small input size (10) seems to not release the GIL.
        # Large input size (1000) seems to release the GIL.
        for nelem in [1, 10, 100, 1000]:
            # inputs
            a = np.arange(nelem, dtype=np.int32)
            acopy = a.copy()
            # run and capture stdout
            with captured_stdout() as buf:
                got = foo(a)
            stdout = buf.getvalue()
            buf.close()
            # process outputs from print
            got_output = sorted(map(lambda x: x.strip(), stdout.splitlines()))
            # build expected output
            expected_output = [str(x % 10) for x in range(nelem)]
            expected_output += [characters[x % 10] for x in range(nelem)]
            expected_output = sorted(expected_output)
            # verify
            self.assertEqual(got_output, expected_output)
            np.testing.assert_equal(got, 2 * acopy)



class TestParGUfuncIssues(unittest.TestCase):

    _numba_parallel_test_ = False

    def test_gil_reacquire_deadlock(self):
        """
        Testing similar issue to #1998 due to GIL reacquiring for Gufunc
        """
        # make a ctypes callback that requires the GIL
        proto = ctypes.CFUNCTYPE(None, ctypes.c_int32)
        characters = 'abcdefghij'

        def bar(x):
            print(characters[x])

        cbar = proto(bar)

        # our unit under test
        @guvectorize(['(int32, int32[:])'], "()->()",
                     target='parallel', nopython=True)
        def foo(x, out):
            print(x % 10)  # this reacquires the GIL
            cbar(x % 10)   # this reacquires the GIL
            out[0] = x * 2

        # Numpy ufunc has a heuristic to determine whether to release the GIL
        # during execution.  Small input size (10) seems to not release the GIL.
        # Large input size (1000) seems to release the GIL.
        for nelem in [1, 10, 100, 1000]:
            # inputs
            a = np.arange(nelem, dtype=np.int32)
            acopy = a.copy()
            # run and capture stdout
            with captured_stdout() as buf:
                got = foo(a)
            stdout = buf.getvalue()
            buf.close()
            # process outputs from print
            got_output = sorted(map(lambda x: x.strip(), stdout.splitlines()))
            # build expected output
            expected_output = [str(x % 10) for x in range(nelem)]
            expected_output += [characters[x % 10] for x in range(nelem)]
            expected_output = sorted(expected_output)
            # verify
            self.assertEqual(got_output, expected_output)
            np.testing.assert_equal(got, 2 * acopy)


if __name__ == '__main__':
    unittest.main()
