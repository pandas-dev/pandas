import os
import json
import unittest
from textwrap import dedent
from tempfile import TemporaryDirectory

from numba.tests.support import TestCase, run_in_subprocess


class TestChromeTraceModule(TestCase):
    """
    Test chrome tracing generated file(s).
    """

    def test_trace_output(self):
        code = """
            from numba import njit
            import numpy as np

            x = np.arange(100).reshape(10, 10)

            @njit
            def go_fast(a):
                trace = 0.0
                for i in range(a.shape[0]):
                    trace += np.tanh(a[i, i])
                return a + trace

            go_fast(x)
        """

        src = dedent(code)
        with TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "test_trace.json")
            env = os.environ.copy()
            env['NUMBA_CHROME_TRACE'] = path
            run_in_subprocess(src, env=env)
            with open(path) as file:
                events = json.load(file)
                self.assertIsInstance(events, list)
                for ev in events:
                    self.assertIsInstance(ev, dict)
                    # check that each record has the right fields.
                    self.assertEqual(
                        set(ev.keys()),
                        {"cat", "pid", "tid", "ph", "name", "args", "ts"},
                    )


if __name__ == "__main__":
    unittest.main()
