import os.path
import subprocess
import sys

import numpy as np

from numba.tests.support import skip_parfors_unsupported
from .test_caching import DispatcherCacheUsecasesTest


@skip_parfors_unsupported
class TestParforsCache(DispatcherCacheUsecasesTest):
    here = os.path.dirname(__file__)
    usecases_file = os.path.join(here, "parfors_cache_usecases.py")
    modname = "parfors_caching_test_fodder"

    def run_test(self, fname, num_funcs=1):
        mod = self.import_module()
        self.check_pycache(0)
        f = getattr(mod, fname)
        ary = np.ones(10)
        # The result of these functions is derived from e.g. out of order
        # accumulation so allclose should be fine.
        np.testing.assert_allclose(f(ary), f.py_func(ary))

        dynamic_globals = [cres.library.has_dynamic_globals
                           for cres in f.overloads.values()]
        [cres] = f.overloads.values()
        self.assertEqual(dynamic_globals, [False])
        # For each cached func, there're 2 entries (index + data)
        self.check_pycache(num_funcs * 2)

        self.run_in_separate_process()

    def test_arrayexprs(self):
        f = 'arrayexprs_case'
        self.run_test(f)

    def test_prange(self):
        f = 'prange_case'
        self.run_test(f)

    def test_caller(self):
        f = 'caller_case'
        # num_funcs=3 because, there's the `caller_case()` which calls
        # the `prange_case()` and `arrayexprs_case()`
        self.run_test(f, num_funcs=3)


@skip_parfors_unsupported
class TestParforsCacheChangingThreads(DispatcherCacheUsecasesTest):
    # NOTE: This test is checking issue #7518, that thread counts are not
    # baked into cached objects.

    here = os.path.dirname(__file__)
    usecases_file = os.path.join(here, "parfors_cache_usecases.py")
    modname = "parfors_caching_test_fodder"

    def run_in_separate_process(self, thread_count):
        # Cached functions can be run from a distinct process.
        code = """if 1:
            import sys

            sys.path.insert(0, %(tempdir)r)
            mod = __import__(%(modname)r)
            mod.self_run()
            """ % dict(tempdir=self.tempdir, modname=self.modname)

        new_env = {**os.environ, "NUMBA_NUM_THREADS" : str(thread_count)}
        popen = subprocess.Popen([sys.executable, "-c", code],
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 env=new_env)
        out, err = popen.communicate()
        if popen.returncode != 0:
            raise AssertionError(f"process failed with code {popen.returncode}:"
                                 f"stderr follows\n{err.decode()}\n")

    def test_caching(self):
        self.check_pycache(0)
        self.run_in_separate_process(1)
        self.check_pycache(3 * 2) # ran 3 functions, 2 entries each
        self.run_in_separate_process(2)
        self.check_pycache(3 * 2) # ran 3 functions, 2 entries each
