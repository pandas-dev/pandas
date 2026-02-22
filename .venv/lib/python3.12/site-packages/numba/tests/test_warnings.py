import os
import subprocess
import sys
import warnings
import numpy as np

import unittest
from numba import jit
from numba.core.errors import (
    NumbaWarning,
    deprecated,
    NumbaDeprecationWarning,
    NumbaPendingDeprecationWarning,
)
from numba.core import errors
from numba.tests.support import ignore_internal_warnings


class TestBuiltins(unittest.TestCase):

    def check_objmode_deprecation_warning(self, w):
        # Object mode fall-back is slated for deprecation, check the warning
        msg = ("Fall-back from the nopython compilation path to the object "
               "mode compilation path has been detected")
        self.assertEqual(w.category, NumbaDeprecationWarning)
        self.assertIn(msg, str(w.message))

    def check_nopython_kwarg_missing_warning(self, w):
        # nopython default is scheduled to change when objmode fall-back is
        # removed, check warning.
        msg = ("The \'nopython\' keyword argument was not supplied")
        self.assertEqual(w.category, NumbaDeprecationWarning)
        self.assertIn(msg, str(w.message))

    def test_return_type_warning_with_nrt(self):
        """
        Rerun test_return_type_warning with nrt
        """
        y = np.ones(4, dtype=np.float32)

        def return_external_array():
            return y

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', NumbaWarning)
            ignore_internal_warnings()

            cfunc = jit(nopython=True)(return_external_array)
            cfunc()
            # No more warning
            self.assertEqual(len(w), 0)

    def test_no_warning_with_forceobj(self):
        def add(x, y):
            a = [] # noqa dead
            return x + y

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', NumbaWarning)
            ignore_internal_warnings()

            cfunc = jit(add, forceobj=True)
            cfunc(1, 2)

            self.assertEqual(len(w), 0)

    def test_deprecated(self):
        @deprecated('foo')
        def bar():
            pass

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            ignore_internal_warnings()
            bar()

            self.assertEqual(len(w), 1)
            self.assertEqual(w[0].category, DeprecationWarning)
            self.assertIn('bar', str(w[0].message))
            self.assertIn('foo', str(w[0].message))

    def test_warnings_fixer(self):
        # For some context, see #4083

        wfix = errors.WarningsFixer(errors.NumbaWarning)
        with wfix.catch_warnings('foo', 10):
            warnings.warn(errors.NumbaWarning('same'))
            warnings.warn(errors.NumbaDeprecationWarning('same'))
            ignore_internal_warnings()

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            ignore_internal_warnings()
            wfix.flush()

            self.assertEqual(len(w), 2)
            # the order of these will be backwards to the above, the
            # WarningsFixer flush method sorts with a key based on str
            # comparison
            self.assertEqual(w[0].category, NumbaDeprecationWarning)
            self.assertEqual(w[1].category, NumbaWarning)
            self.assertIn('same', str(w[0].message))
            self.assertIn('same', str(w[1].message))

    def test_disable_performance_warnings(self):

        not_found_ret_code = 55
        found_ret_code = 99
        expected = "'parallel=True' was specified but no transformation"

        # NOTE: the error_usecases is needed as the NumbaPerformanceWarning's
        # for parallel=True failing to parallelise do not appear for functions
        # defined by string eval/exec etc.
        parallel_code = """if 1:
            import warnings
            from numba.tests.error_usecases import foo
            import numba
            from numba.tests.support import ignore_internal_warnings
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter('always')
                ignore_internal_warnings()
                foo()
            for x in w:
                if x.category == numba.errors.NumbaPerformanceWarning:
                    if "%s" in str(x.message):
                        exit(%s)
            exit(%s)
        """ % (expected, found_ret_code, not_found_ret_code)

        # run in the standard env, warning should raise
        popen = subprocess.Popen([sys.executable, "-c", parallel_code])
        out, err = popen.communicate()
        self.assertEqual(popen.returncode, found_ret_code)

        # run in an env with performance warnings disabled, should not warn
        env = dict(os.environ)
        env['NUMBA_DISABLE_PERFORMANCE_WARNINGS'] = "1"
        popen = subprocess.Popen([sys.executable, "-c", parallel_code], env=env)
        out, err = popen.communicate()
        self.assertEqual(popen.returncode, not_found_ret_code)

    def test_filter_deprecation_warnings(self):
        # Filter on base classes of deprecation warnings should apply to Numba's
        # deprecation warnings
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            warnings.simplefilter('ignore', category=DeprecationWarning)
            warnings.simplefilter('ignore', category=PendingDeprecationWarning)
            warnings.warn(DeprecationWarning("this is ignored"))
            warnings.warn(PendingDeprecationWarning("this is ignored"))
            warnings.warn(NumbaDeprecationWarning("this is ignored"))
            warnings.warn(NumbaPendingDeprecationWarning("this is ignored"))
            with self.assertRaises(NumbaWarning):
                warnings.warn(NumbaWarning("this is not ignored"))

    def test_filter_ignore_numba_deprecation_only(self):
        # Make a filter that ignores Numba's deprecation warnings but raises on
        # other deprecation warnings
        with warnings.catch_warnings():
            warnings.simplefilter('error', category=DeprecationWarning)
            warnings.simplefilter('error', category=PendingDeprecationWarning)
            warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
            warnings.simplefilter('ignore',
                                  category=NumbaPendingDeprecationWarning)

            with self.assertRaises(DeprecationWarning):
                warnings.warn(DeprecationWarning("this is not ignored"))
            with self.assertRaises(PendingDeprecationWarning):
                warnings.warn(PendingDeprecationWarning("this is not ignored"))

            warnings.warn(NumbaDeprecationWarning("this is ignored"))
            warnings.warn(NumbaPendingDeprecationWarning("this is ignored"))

            # now make it so that Numba deprecation warnings are raising
            warnings.simplefilter('error', category=NumbaDeprecationWarning)
            warnings.simplefilter('error',
                                  category=NumbaPendingDeprecationWarning)

            with self.assertRaises(DeprecationWarning):
                warnings.warn(NumbaDeprecationWarning("this is not ignored"))
            with self.assertRaises(PendingDeprecationWarning):
                warnings.warn(NumbaPendingDeprecationWarning(
                    "this is not ignored"))


if __name__ == '__main__':
    unittest.main()
