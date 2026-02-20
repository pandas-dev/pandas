import cProfile as profiler
import os
import pstats
import subprocess
import sys

import numpy as np

from numba import jit, types
from numba.tests.support import needs_blas
import unittest


def generate_standard_dot_case():

    @jit((types.float32[::1], types.float32[::1],))
    def dot(a, b):
        sum = 0
        for i in range(len(a)):
            sum += a[i]*b[i]
        return sum

    dot._enable_sysmon = True
    return dot, dot


def generate_raising_dot_case():

    @jit((types.float32[::1], types.float32[::1],))
    def raising_dot(a, b):
        # this is like dot above, it does all the work, but the raises
        sum = 0
        for i in range(len(a)):
            sum += a[i]*b[i]
        raise ValueError("problem with dot")

    raising_dot._enable_sysmon = True

    def call_raising_dot(a, b):
        try:
            raising_dot(a, b)
        except ValueError:
            pass

    return raising_dot, call_raising_dot


def np_dot(a, b):
    return np.dot(a, b)


class TestProfiler(unittest.TestCase):

    def check_profiler_dot(self, caller, cfunc):
        """
        Make sure the jit-compiled function shows up in the profile stats
        as a regular Python function.
        """
        a = np.arange(16, dtype=np.float32)
        b = np.arange(16, dtype=np.float32)
        n_calls = 123
        p = profiler.Profile()
        p.enable()
        try:
            for _ in range(n_calls):
                caller(a, b)
        finally:
            p.disable()
        stats = pstats.Stats(p).strip_dirs()

        def check_stats_for_key(stats, code, n_calls):
            expected_key = (os.path.basename(code.co_filename),
                            code.co_firstlineno,
                            code.co_name,
                            )
            # check the key is in the stats
            self.assertIn(expected_key, stats.stats)
            # check that call for the key has been made `n_calls` times.
            func_stats = stats.stats[expected_key]
            self.assertEqual(func_stats[:2], (n_calls, n_calls))

        # check the JIT compiled function
        check_stats_for_key(stats, cfunc.py_func.__code__, n_calls)

        # check the caller if it's not the same as the cfunc
        if caller is not cfunc:
            check_stats_for_key(stats, caller.__code__, n_calls)

    def test_profiler(self):
        dot, _ = generate_standard_dot_case()
        self.check_profiler_dot(dot, dot)

    def test_profiler_for_raising_function(self):
        raising_dot, call_raising_dot = generate_raising_dot_case()
        self.check_profiler_dot(call_raising_dot, raising_dot)

    @needs_blas
    def test_profiler_np_dot(self):
        # Issue #1786: initializing BLAS would crash when profiling
        code = """if 1:
            import cProfile as profiler

            import numpy as np

            from numba import jit
            from numba.tests.test_profiler import np_dot

            cfunc = jit(nopython=True)(np_dot)

            a = np.arange(16, dtype=np.float32)
            b = np.arange(16, dtype=np.float32)

            p = profiler.Profile()
            p.enable()
            cfunc(a, b)
            cfunc(a, b)
            p.disable()
            """
        subprocess.check_call([sys.executable, "-c", code])

    def test_issue_3229(self):
        # Issue #3229: Seemingly random segfaults when profiling due to
        # frame injection.
        # numba.tests.npyufunc.test_dufunc.TestDUFunc.test_npm_call is the
        # first test case crashing when profiling. Fingers crossed fixing
        # this is sufficient proof for the general case.

        code = """if 1:
            import cProfile as profiler
            p = profiler.Profile()
            p.enable()

            from numba.tests.npyufunc.test_dufunc import TestDUFunc
            t = TestDUFunc('test_npm_call')
            t.test_npm_call()

            p.disable()
            """
        subprocess.check_call([sys.executable, "-c", code])

if __name__ == '__main__':
    unittest.main()
