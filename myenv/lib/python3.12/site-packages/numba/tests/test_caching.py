import inspect
import llvmlite.binding as ll
import multiprocessing
import numpy as np
import os
import stat
import shutil
import subprocess
import sys
import traceback
import unittest
import warnings
from numba import njit
from numba.core import codegen
from numba.core.caching import _UserWideCacheLocator
from numba.core.errors import NumbaWarning
from numba.parfors import parfor
from numba.tests.support import (
    TestCase,
    SerialMixin,
    capture_cache_log,
    import_dynamic,
    override_config,
    run_in_new_process_caching,
    skip_if_typeguard,
    skip_parfors_unsupported,
    temp_directory,
)

try:
    import ipykernel
except ImportError:
    ipykernel = None


def check_access_is_preventable():
    # This exists to check whether it is possible to prevent access to
    # a file/directory through the use of `chmod 500`. If a user has
    # elevated rights (e.g. root) then writes are likely to be possible
    # anyway. Tests that require functioning access prevention are
    # therefore skipped based on the result of this check.
    tempdir = temp_directory('test_cache')
    test_dir = (os.path.join(tempdir, 'writable_test'))
    os.mkdir(test_dir)
    # check a write is possible
    with open(os.path.join(test_dir, 'write_ok'), 'wt') as f:
        f.write('check1')
    # now forbid access
    os.chmod(test_dir, 0o500)
    try:
        with open(os.path.join(test_dir, 'write_forbidden'), 'wt') as f:
            f.write('check2')
        # access prevention is not possible
        return False
    except PermissionError:
        # Check that the cause of the exception is due to access/permission
        # as per
        # https://github.com/conda/conda/blob/4.5.0/conda/gateways/disk/permissions.py#L35-L37  # noqa: E501
        # errno reports access/perm fail so access prevention via
        # `chmod 500` works for this user.
        return True
    finally:
        os.chmod(test_dir, 0o775)
        shutil.rmtree(test_dir)


_access_preventable = check_access_is_preventable()
_access_msg = "Cannot create a directory to which writes are preventable"
skip_bad_access = unittest.skipUnless(_access_preventable, _access_msg)


def constant_unicode_cache():
    c = "abcd"
    return hash(c), c


def check_constant_unicode_cache():
    pyfunc = constant_unicode_cache
    cfunc = njit(cache=True)(pyfunc)
    exp_hv, exp_str = pyfunc()
    got_hv, got_str = cfunc()
    assert exp_hv == got_hv
    assert exp_str == got_str


def dict_cache():
    return {'a': 1, 'b': 2}


def check_dict_cache():
    pyfunc = dict_cache
    cfunc = njit(cache=True)(pyfunc)
    exp = pyfunc()
    got = cfunc()
    assert exp == got


def generator_cache():
    for v in (1, 2, 3):
        yield v


def check_generator_cache():
    pyfunc = generator_cache
    cfunc = njit(cache=True)(pyfunc)
    exp = list(pyfunc())
    got = list(cfunc())
    assert exp == got


class TestCaching(SerialMixin, TestCase):
    def run_test(self, func):
        func()
        res = run_in_new_process_caching(func)
        self.assertEqual(res['exitcode'], 0)

    def test_constant_unicode_cache(self):
        self.run_test(check_constant_unicode_cache)

    def test_dict_cache(self):
        self.run_test(check_dict_cache)

    def test_generator_cache(self):
        self.run_test(check_generator_cache)

    def test_omitted(self):

        # Test in a new directory
        cache_dir = temp_directory(self.__class__.__name__)
        ctx = multiprocessing.get_context()
        result_queue = ctx.Queue()
        proc = ctx.Process(
            target=omitted_child_test_wrapper,
            args=(result_queue, cache_dir, False),
        )
        proc.start()
        proc.join()
        success, output = result_queue.get()

        # Ensure the child process is completed before checking its output
        if not success:
            self.fail(output)

        self.assertEqual(
            output,
            1000,
            "Omitted function returned an incorrect output"
        )

        proc = ctx.Process(
            target=omitted_child_test_wrapper,
            args=(result_queue, cache_dir, True)
        )
        proc.start()
        proc.join()
        success, output = result_queue.get()

        # Ensure the child process is completed before checking its output
        if not success:
            self.fail(output)

        self.assertEqual(
            output,
            1000,
            "Omitted function returned an incorrect output"
        )


def omitted_child_test_wrapper(result_queue, cache_dir, second_call):
    with override_config("CACHE_DIR", cache_dir):
        @njit(cache=True)
        def test(num=1000):
            return num

        try:
            output = test()
            # If we have a second call, we should have a cache hit.
            # Otherwise, we expect a cache miss.
            if second_call:
                assert test._cache_hits[test.signatures[0]] == 1, \
                    "Cache did not hit as expected"
                assert test._cache_misses[test.signatures[0]] == 0, \
                    "Cache has an unexpected miss"
            else:
                assert test._cache_misses[test.signatures[0]] == 1, \
                    "Cache did not miss as expected"
                assert test._cache_hits[test.signatures[0]] == 0, \
                    "Cache has an unexpected hit"
            success = True
        # Catch anything raised so it can be propagated
        except: # noqa: E722
            output = traceback.format_exc()
            success = False
        result_queue.put((success, output))


class BaseCacheTest(TestCase):
    # The source file that will be copied
    usecases_file = None
    # Make sure this doesn't conflict with another module
    modname = None

    def setUp(self):
        self.tempdir = temp_directory('test_cache')
        sys.path.insert(0, self.tempdir)
        self.modfile = os.path.join(self.tempdir, self.modname + ".py")
        self.cache_dir = os.path.join(self.tempdir, "__pycache__")
        shutil.copy(self.usecases_file, self.modfile)
        os.chmod(self.modfile, stat.S_IREAD | stat.S_IWRITE)
        self.maxDiff = None

    def tearDown(self):
        sys.modules.pop(self.modname, None)
        sys.path.remove(self.tempdir)

    def import_module(self):
        # Import a fresh version of the test module.  All jitted functions
        # in the test module will start anew and load overloads from
        # the on-disk cache if possible.
        old = sys.modules.pop(self.modname, None)
        if old is not None:
            # Make sure cached bytecode is removed
            cached = [old.__cached__]
            for fn in cached:
                try:
                    os.unlink(fn)
                except FileNotFoundError:
                    pass
        mod = import_dynamic(self.modname)
        self.assertEqual(mod.__file__.rstrip('co'), self.modfile)
        return mod

    def cache_contents(self):
        try:
            return [fn for fn in os.listdir(self.cache_dir)
                    if not fn.endswith(('.pyc', ".pyo"))]
        except FileNotFoundError:
            return []

    def get_cache_mtimes(self):
        return dict((fn, os.path.getmtime(os.path.join(self.cache_dir, fn)))
                    for fn in sorted(self.cache_contents()))

    def check_pycache(self, n):
        c = self.cache_contents()
        self.assertEqual(len(c), n, c)

    def dummy_test(self):
        pass


class DispatcherCacheUsecasesTest(BaseCacheTest):
    here = os.path.dirname(__file__)
    usecases_file = os.path.join(here, "cache_usecases.py")
    modname = "dispatcher_caching_test_fodder"

    def run_in_separate_process(self, *, envvars={}):
        # Cached functions can be run from a distinct process.
        # Also stresses issue #1603: uncached function calling cached function
        # shouldn't fail compiling.
        code = """if 1:
            import sys

            sys.path.insert(0, %(tempdir)r)
            mod = __import__(%(modname)r)
            mod.self_test()
            """ % dict(tempdir=self.tempdir, modname=self.modname)

        subp_env = os.environ.copy()
        subp_env.update(envvars)
        popen = subprocess.Popen([sys.executable, "-c", code],
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 env=subp_env)
        out, err = popen.communicate()
        if popen.returncode != 0:
            raise AssertionError(
                "process failed with code %s: \n"
                "stdout follows\n%s\n"
                "stderr follows\n%s\n"
                % (popen.returncode, out.decode(), err.decode()),
            )

    def check_hits(self, func, hits, misses=None):
        st = func.stats
        self.assertEqual(sum(st.cache_hits.values()), hits, st.cache_hits)
        if misses is not None:
            self.assertEqual(sum(st.cache_misses.values()), misses,
                             st.cache_misses)


class TestCache(DispatcherCacheUsecasesTest):

    def test_caching(self):
        self.check_pycache(0)
        mod = self.import_module()
        self.check_pycache(0)

        f = mod.add_usecase
        self.assertPreciseEqual(f(2, 3), 6)
        self.check_pycache(2)  # 1 index, 1 data
        self.assertPreciseEqual(f(2.5, 3), 6.5)
        self.check_pycache(3)  # 1 index, 2 data
        self.check_hits(f, 0, 2)

        f = mod.add_objmode_usecase
        self.assertPreciseEqual(f(2, 3), 6)
        self.check_pycache(5)  # 2 index, 3 data
        self.assertPreciseEqual(f(2.5, 3), 6.5)
        self.check_pycache(6)  # 2 index, 4 data
        self.check_hits(f, 0, 2)

        f = mod.record_return
        rec = f(mod.aligned_arr, 1)
        self.assertPreciseEqual(tuple(rec), (2, 43.5))
        rec = f(mod.packed_arr, 1)
        self.assertPreciseEqual(tuple(rec), (2, 43.5))
        self.check_pycache(9)  # 3 index, 6 data
        self.check_hits(f, 0, 2)

        # Check the code runs ok from another process
        self.run_in_separate_process()

    def test_caching_nrt_pruned(self):
        self.check_pycache(0)
        mod = self.import_module()
        self.check_pycache(0)

        f = mod.add_usecase
        self.assertPreciseEqual(f(2, 3), 6)
        self.check_pycache(2)  # 1 index, 1 data
        # NRT pruning may affect cache
        self.assertPreciseEqual(f(2, np.arange(3)), 2 + np.arange(3) + 1)
        self.check_pycache(3)  # 1 index, 2 data
        self.check_hits(f, 0, 2)

    def test_inner_then_outer(self):
        # Caching inner then outer function is ok
        mod = self.import_module()
        self.assertPreciseEqual(mod.inner(3, 2), 6)
        self.check_pycache(2)  # 1 index, 1 data
        # Uncached outer function shouldn't fail (issue #1603)
        f = mod.outer_uncached
        self.assertPreciseEqual(f(3, 2), 2)
        self.check_pycache(2)  # 1 index, 1 data
        mod = self.import_module()
        f = mod.outer_uncached
        self.assertPreciseEqual(f(3, 2), 2)
        self.check_pycache(2)  # 1 index, 1 data
        # Cached outer will create new cache entries
        f = mod.outer
        self.assertPreciseEqual(f(3, 2), 2)
        self.check_pycache(4)  # 2 index, 2 data
        self.assertPreciseEqual(f(3.5, 2), 2.5)
        self.check_pycache(6)  # 2 index, 4 data

    def test_outer_then_inner(self):
        # Caching outer then inner function is ok
        mod = self.import_module()
        self.assertPreciseEqual(mod.outer(3, 2), 2)
        self.check_pycache(4)  # 2 index, 2 data
        self.assertPreciseEqual(mod.outer_uncached(3, 2), 2)
        self.check_pycache(4)  # same
        mod = self.import_module()
        f = mod.inner
        self.assertPreciseEqual(f(3, 2), 6)
        self.check_pycache(4)  # same
        self.assertPreciseEqual(f(3.5, 2), 6.5)
        self.check_pycache(5)  # 2 index, 3 data

    def test_no_caching(self):
        mod = self.import_module()

        f = mod.add_nocache_usecase
        self.assertPreciseEqual(f(2, 3), 6)
        self.check_pycache(0)

    def test_looplifted(self):
        # Loop-lifted functions can't be cached and raise a warning
        mod = self.import_module()

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', NumbaWarning)

            f = mod.looplifted
            self.assertPreciseEqual(f(4), 6)
            self.check_pycache(0)

        self.assertEqual(len(w), 1)
        self.assertIn('Cannot cache compiled function "looplifted" '
                      'as it uses lifted code', str(w[0].message))

    def test_big_array(self):
        # Code references big array globals cannot be cached
        mod = self.import_module()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', NumbaWarning)

            f = mod.use_big_array
            np.testing.assert_equal(f(), mod.biggie)
            self.check_pycache(0)

        self.assertEqual(len(w), 1)
        self.assertIn('Cannot cache compiled function "use_big_array" '
                      'as it uses dynamic globals', str(w[0].message))

    def test_ctypes(self):
        # Functions using a ctypes pointer can't be cached and raise
        # a warning.
        mod = self.import_module()

        for f in [mod.use_c_sin, mod.use_c_sin_nest1, mod.use_c_sin_nest2]:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter('always', NumbaWarning)

                self.assertPreciseEqual(f(0.0), 0.0)
                self.check_pycache(0)

            self.assertEqual(len(w), 1)
            self.assertIn(
                'Cannot cache compiled function "{}"'.format(f.__name__),
                str(w[0].message),
            )

    def test_closure(self):
        mod = self.import_module()

        with warnings.catch_warnings():
            warnings.simplefilter('error', NumbaWarning)

            f = mod.closure1
            self.assertPreciseEqual(f(3), 6) # 3 + 3 = 6
            f = mod.closure2
            self.assertPreciseEqual(f(3), 8) # 3 + 5 = 8
            f = mod.closure3
            self.assertPreciseEqual(f(3), 10) # 3 + 7 = 10
            f = mod.closure4
            self.assertPreciseEqual(f(3), 12) # 3 + 9 = 12
            self.check_pycache(5) # 1 nbi, 4 nbc

    def test_first_class_function(self):
        mod = self.import_module()
        f = mod.first_class_function_usecase
        self.assertEqual(f(mod.first_class_function_mul, 1), 1)
        self.assertEqual(f(mod.first_class_function_mul, 10), 100)
        self.assertEqual(f(mod.first_class_function_add, 1), 2)
        self.assertEqual(f(mod.first_class_function_add, 10), 20)
        # 1 + 1 + 1 nbi, 1 + 1 + 2 nbc - a separate cache for each call to `f`
        # with a different callback.
        self.check_pycache(7)

    def test_cache_reuse(self):
        mod = self.import_module()
        mod.add_usecase(2, 3)
        mod.add_usecase(2.5, 3.5)
        mod.add_objmode_usecase(2, 3)
        mod.outer_uncached(2, 3)
        mod.outer(2, 3)
        mod.record_return(mod.packed_arr, 0)
        mod.record_return(mod.aligned_arr, 1)
        mtimes = self.get_cache_mtimes()
        # Two signatures compiled
        self.check_hits(mod.add_usecase, 0, 2)

        mod2 = self.import_module()
        self.assertIsNot(mod, mod2)
        f = mod2.add_usecase
        f(2, 3)
        self.check_hits(f, 1, 0)
        f(2.5, 3.5)
        self.check_hits(f, 2, 0)
        f = mod2.add_objmode_usecase
        f(2, 3)
        self.check_hits(f, 1, 0)

        # The files haven't changed
        self.assertEqual(self.get_cache_mtimes(), mtimes)

        self.run_in_separate_process()
        self.assertEqual(self.get_cache_mtimes(), mtimes)

    def test_cache_invalidate(self):
        mod = self.import_module()
        f = mod.add_usecase
        self.assertPreciseEqual(f(2, 3), 6)

        # This should change the functions' results
        with open(self.modfile, "a") as f:
            f.write("\nZ = 10\n")

        mod = self.import_module()
        f = mod.add_usecase
        self.assertPreciseEqual(f(2, 3), 15)
        f = mod.add_objmode_usecase
        self.assertPreciseEqual(f(2, 3), 15)

    def test_recompile(self):
        # Explicit call to recompile() should overwrite the cache
        mod = self.import_module()
        f = mod.add_usecase
        self.assertPreciseEqual(f(2, 3), 6)

        mod = self.import_module()
        f = mod.add_usecase
        mod.Z = 10
        self.assertPreciseEqual(f(2, 3), 6)
        f.recompile()
        self.assertPreciseEqual(f(2, 3), 15)

        # Freshly recompiled version is re-used from other imports
        mod = self.import_module()
        f = mod.add_usecase
        self.assertPreciseEqual(f(2, 3), 15)

    def test_same_names(self):
        # Function with the same names should still disambiguate
        mod = self.import_module()
        f = mod.renamed_function1
        self.assertPreciseEqual(f(2), 4)
        f = mod.renamed_function2
        self.assertPreciseEqual(f(2), 8)

    def test_frozen(self):
        from .dummy_module import function
        old_code = function.__code__
        code_obj = compile('pass', 'tests/dummy_module.py', 'exec')
        try:
            function.__code__ = code_obj

            source = inspect.getfile(function)
            # doesn't return anything, since it cannot find the module
            # fails unless the executable is frozen
            locator = _UserWideCacheLocator.from_function(function, source)
            self.assertIsNone(locator)

            sys.frozen = True
            # returns a cache locator object, only works when the executable
            # is frozen
            locator = _UserWideCacheLocator.from_function(function, source)
            self.assertIsInstance(locator, _UserWideCacheLocator)

        finally:
            function.__code__ = old_code
            del sys.frozen

    def _test_pycache_fallback(self):
        """
        With a disabled __pycache__, test there is a working fallback
        (e.g. on the user-wide cache dir)
        """
        mod = self.import_module()
        f = mod.add_usecase
        # Remove this function's cache files at the end, to avoid accumulation
        # across test calls.
        self.addCleanup(shutil.rmtree, f.stats.cache_path, ignore_errors=True)

        self.assertPreciseEqual(f(2, 3), 6)
        # It's a cache miss since the file was copied to a new temp location
        self.check_hits(f, 0, 1)

        # Test re-use
        mod2 = self.import_module()
        f = mod2.add_usecase
        self.assertPreciseEqual(f(2, 3), 6)
        self.check_hits(f, 1, 0)

        # The __pycache__ is empty (otherwise the test's preconditions
        # wouldn't be met)
        self.check_pycache(0)

    @skip_bad_access
    @unittest.skipIf(os.name == "nt",
                     "cannot easily make a directory read-only on Windows")
    def test_non_creatable_pycache(self):
        # Make it impossible to create the __pycache__ directory
        old_perms = os.stat(self.tempdir).st_mode
        os.chmod(self.tempdir, 0o500)
        self.addCleanup(os.chmod, self.tempdir, old_perms)

        self._test_pycache_fallback()

    @skip_bad_access
    @unittest.skipIf(os.name == "nt",
                     "cannot easily make a directory read-only on Windows")
    def test_non_writable_pycache(self):
        # Make it impossible to write to the __pycache__ directory
        pycache = os.path.join(self.tempdir, '__pycache__')
        os.mkdir(pycache)
        old_perms = os.stat(pycache).st_mode
        os.chmod(pycache, 0o500)
        self.addCleanup(os.chmod, pycache, old_perms)

        self._test_pycache_fallback()

    def test_ipython(self):
        # Test caching in an IPython session
        base_cmd = [sys.executable, '-m', 'IPython']
        base_cmd += ['--quiet', '--quick', '--no-banner', '--colors=NoColor']
        try:
            ver = subprocess.check_output(base_cmd + ['--version'])
        except subprocess.CalledProcessError as e:
            self.skipTest("ipython not available: return code %d"
                          % e.returncode)
        ver = ver.strip().decode()
        # Create test input
        inputfn = os.path.join(self.tempdir, "ipython_cache_usecase.txt")
        with open(inputfn, "w") as f:
            f.write(r"""
                import os
                import sys

                from numba import jit

                # IPython 5 does not support multiline input if stdin isn't
                # a tty (https://github.com/ipython/ipython/issues/9752)
                f = jit(cache=True)(lambda: 42)

                res = f()
                # IPython writes on stdout, so use stderr instead
                sys.stderr.write(u"cache hits = %d\n" % f.stats.cache_hits[()])

                # IPython hijacks sys.exit(), bypass it
                sys.stdout.flush()
                sys.stderr.flush()
                os._exit(res)
                """)

        def execute_with_input():
            # Feed the test input as stdin, to execute it in REPL context
            with open(inputfn, "rb") as stdin:
                p = subprocess.Popen(base_cmd, stdin=stdin,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     universal_newlines=True)
                out, err = p.communicate()
                if p.returncode != 42:
                    self.fail("unexpected return code %d\n"
                              "-- stdout:\n%s\n"
                              "-- stderr:\n%s\n"
                              % (p.returncode, out, err))
                return err

        execute_with_input()
        # Run a second time and check caching
        err = execute_with_input()
        self.assertIn("cache hits = 1", err.strip())

    @unittest.skipIf((ipykernel is None) or (ipykernel.version_info[0] < 6),
                     "requires ipykernel >= 6")
    def test_ipykernel(self):
        # Test caching in an IPython session using ipykernel

        base_cmd = [sys.executable, '-m', 'IPython']
        base_cmd += ['--quiet', '--quick', '--no-banner', '--colors=NoColor']
        try:
            ver = subprocess.check_output(base_cmd + ['--version'])
        except subprocess.CalledProcessError as e:
            self.skipTest("ipython not available: return code %d"
                          % e.returncode)
        ver = ver.strip().decode()
        # Create test input
        from ipykernel import compiler
        inputfn = compiler.get_tmp_directory()
        with open(inputfn, "w") as f:
            f.write(r"""
                import os
                import sys

                from numba import jit

                # IPython 5 does not support multiline input if stdin isn't
                # a tty (https://github.com/ipython/ipython/issues/9752)
                f = jit(cache=True)(lambda: 42)

                res = f()
                # IPython writes on stdout, so use stderr instead
                sys.stderr.write(u"cache hits = %d\n" % f.stats.cache_hits[()])

                # IPython hijacks sys.exit(), bypass it
                sys.stdout.flush()
                sys.stderr.flush()
                os._exit(res)
                """)

        def execute_with_input():
            # Feed the test input as stdin, to execute it in REPL context
            with open(inputfn, "rb") as stdin:
                p = subprocess.Popen(base_cmd, stdin=stdin,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     universal_newlines=True)
                out, err = p.communicate()
                if p.returncode != 42:
                    self.fail("unexpected return code %d\n"
                              "-- stdout:\n%s\n"
                              "-- stderr:\n%s\n"
                              % (p.returncode, out, err))
                return err

        execute_with_input()
        # Run a second time and check caching
        err = execute_with_input()
        self.assertIn("cache hits = 1", err.strip())


@skip_parfors_unsupported
class TestSequentialParForsCache(DispatcherCacheUsecasesTest):
    def setUp(self):
        super(TestSequentialParForsCache, self).setUp()
        # Turn on sequential parfor lowering
        parfor.sequential_parfor_lowering = True

    def tearDown(self):
        super(TestSequentialParForsCache, self).tearDown()
        # Turn off sequential parfor lowering
        parfor.sequential_parfor_lowering = False

    def test_caching(self):
        mod = self.import_module()
        self.check_pycache(0)
        f = mod.parfor_usecase
        ary = np.ones(10)
        self.assertPreciseEqual(f(ary), ary * ary + ary)
        dynamic_globals = [cres.library.has_dynamic_globals
                           for cres in f.overloads.values()]
        self.assertEqual(dynamic_globals, [False])
        self.check_pycache(2)  # 1 index, 1 data


class TestCacheWithCpuSetting(DispatcherCacheUsecasesTest):
    # Disable parallel testing due to envvars modification
    _numba_parallel_test_ = False

    def check_later_mtimes(self, mtimes_old):
        match_count = 0
        for k, v in self.get_cache_mtimes().items():
            if k in mtimes_old:
                self.assertGreaterEqual(v, mtimes_old[k])
                match_count += 1
        self.assertGreater(match_count, 0,
                           msg='nothing to compare')

    def test_user_set_cpu_name(self):
        self.check_pycache(0)
        mod = self.import_module()
        mod.self_test()
        cache_size = len(self.cache_contents())

        mtimes = self.get_cache_mtimes()
        # Change CPU name to generic
        self.run_in_separate_process(envvars={'NUMBA_CPU_NAME': 'generic'})

        self.check_later_mtimes(mtimes)
        self.assertGreater(len(self.cache_contents()), cache_size)
        # Check cache index
        cache = mod.add_usecase._cache
        cache_file = cache._cache_file
        cache_index = cache_file._load_index()
        self.assertEqual(len(cache_index), 2)
        [key_a, key_b] = cache_index.keys()
        if key_a[1][1] == ll.get_host_cpu_name():
            key_host, key_generic = key_a, key_b
        else:
            key_host, key_generic = key_b, key_a
        self.assertEqual(key_host[1][1], ll.get_host_cpu_name())
        self.assertEqual(key_host[1][2], codegen.get_host_cpu_features())
        self.assertEqual(key_generic[1][1], 'generic')
        self.assertEqual(key_generic[1][2], '')

    def test_user_set_cpu_features(self):
        self.check_pycache(0)
        mod = self.import_module()
        mod.self_test()
        cache_size = len(self.cache_contents())

        mtimes = self.get_cache_mtimes()
        # Change CPU feature
        my_cpu_features = '-sse;-avx'

        system_features = codegen.get_host_cpu_features()

        self.assertNotEqual(system_features, my_cpu_features)
        self.run_in_separate_process(
            envvars={'NUMBA_CPU_FEATURES': my_cpu_features},
        )
        self.check_later_mtimes(mtimes)
        self.assertGreater(len(self.cache_contents()), cache_size)
        # Check cache index
        cache = mod.add_usecase._cache
        cache_file = cache._cache_file
        cache_index = cache_file._load_index()
        self.assertEqual(len(cache_index), 2)
        [key_a, key_b] = cache_index.keys()

        if key_a[1][2] == system_features:
            key_host, key_generic = key_a, key_b
        else:
            key_host, key_generic = key_b, key_a

        self.assertEqual(key_host[1][1], ll.get_host_cpu_name())
        self.assertEqual(key_host[1][2], system_features)
        self.assertEqual(key_generic[1][1], ll.get_host_cpu_name())
        self.assertEqual(key_generic[1][2], my_cpu_features)


class TestMultiprocessCache(BaseCacheTest):

    # Nested multiprocessing.Pool raises AssertionError:
    # "daemonic processes are not allowed to have children"
    _numba_parallel_test_ = False

    here = os.path.dirname(__file__)
    usecases_file = os.path.join(here, "cache_usecases.py")
    modname = "dispatcher_caching_test_fodder"

    def test_multiprocessing(self):
        # Check caching works from multiple processes at once (#2028)
        mod = self.import_module()
        # Calling a pure Python caller of the JIT-compiled function is
        # necessary to reproduce the issue.
        f = mod.simple_usecase_caller
        n = 3
        try:
            ctx = multiprocessing.get_context('spawn')
        except AttributeError:
            ctx = multiprocessing
        pool = ctx.Pool(n)
        try:
            res = sum(pool.imap(f, range(n)))
        finally:
            pool.close()
        self.assertEqual(res, n * (n - 1) // 2)


@skip_if_typeguard
class TestCacheFileCollision(unittest.TestCase):
    _numba_parallel_test_ = False

    here = os.path.dirname(__file__)
    usecases_file = os.path.join(here, "cache_usecases.py")
    modname = "caching_file_loc_fodder"
    source_text_1 = """
from numba import njit
@njit(cache=True)
def bar():
    return 123
"""
    source_text_2 = """
from numba import njit
@njit(cache=True)
def bar():
    return 321
"""

    def setUp(self):
        self.tempdir = temp_directory('test_cache_file_loc')
        sys.path.insert(0, self.tempdir)
        self.modname = 'module_name_that_is_unlikely'
        self.assertNotIn(self.modname, sys.modules)
        self.modname_bar1 = self.modname
        self.modname_bar2 = '.'.join([self.modname, 'foo'])
        foomod = os.path.join(self.tempdir, self.modname)
        os.mkdir(foomod)
        with open(os.path.join(foomod, '__init__.py'), 'w') as fout:
            print(self.source_text_1, file=fout)
        with open(os.path.join(foomod, 'foo.py'), 'w') as fout:
            print(self.source_text_2, file=fout)

    def tearDown(self):
        sys.modules.pop(self.modname_bar1, None)
        sys.modules.pop(self.modname_bar2, None)
        sys.path.remove(self.tempdir)

    def import_bar1(self):
        return import_dynamic(self.modname_bar1).bar

    def import_bar2(self):
        return import_dynamic(self.modname_bar2).bar

    def test_file_location(self):
        bar1 = self.import_bar1()
        bar2 = self.import_bar2()
        # Check that the cache file is named correctly
        idxname1 = bar1._cache._cache_file._index_name
        idxname2 = bar2._cache._cache_file._index_name
        self.assertNotEqual(idxname1, idxname2)
        self.assertTrue(idxname1.startswith("__init__.bar-3.py"))
        self.assertTrue(idxname2.startswith("foo.bar-3.py"))

    @unittest.skipUnless(hasattr(multiprocessing, 'get_context'),
                         'Test requires multiprocessing.get_context')
    def test_no_collision(self):
        bar1 = self.import_bar1()
        bar2 = self.import_bar2()
        with capture_cache_log() as buf:
            res1 = bar1()
        cachelog = buf.getvalue()
        # bar1 should save new index and data
        self.assertEqual(cachelog.count('index saved'), 1)
        self.assertEqual(cachelog.count('data saved'), 1)
        self.assertEqual(cachelog.count('index loaded'), 0)
        self.assertEqual(cachelog.count('data loaded'), 0)
        with capture_cache_log() as buf:
            res2 = bar2()
        cachelog = buf.getvalue()
        # bar2 should save new index and data
        self.assertEqual(cachelog.count('index saved'), 1)
        self.assertEqual(cachelog.count('data saved'), 1)
        self.assertEqual(cachelog.count('index loaded'), 0)
        self.assertEqual(cachelog.count('data loaded'), 0)
        self.assertNotEqual(res1, res2)

        try:
            # Make sure we can spawn new process without inheriting
            # the parent context.
            mp = multiprocessing.get_context('spawn')
        except ValueError:
            print("missing spawn context")

        q = mp.Queue()
        # Start new process that calls `cache_file_collision_tester`
        proc = mp.Process(target=cache_file_collision_tester,
                          args=(q, self.tempdir,
                                self.modname_bar1,
                                self.modname_bar2))
        proc.start()
        # Get results from the process
        log1 = q.get()
        got1 = q.get()
        log2 = q.get()
        got2 = q.get()
        proc.join()

        # The remote execution result of bar1() and bar2() should match
        # the one executed locally.
        self.assertEqual(got1, res1)
        self.assertEqual(got2, res2)

        # The remote should have loaded bar1 from cache
        self.assertEqual(log1.count('index saved'), 0)
        self.assertEqual(log1.count('data saved'), 0)
        self.assertEqual(log1.count('index loaded'), 1)
        self.assertEqual(log1.count('data loaded'), 1)

        # The remote should have loaded bar2 from cache
        self.assertEqual(log2.count('index saved'), 0)
        self.assertEqual(log2.count('data saved'), 0)
        self.assertEqual(log2.count('index loaded'), 1)
        self.assertEqual(log2.count('data loaded'), 1)


def cache_file_collision_tester(q, tempdir, modname_bar1, modname_bar2):
    sys.path.insert(0, tempdir)
    bar1 = import_dynamic(modname_bar1).bar
    bar2 = import_dynamic(modname_bar2).bar
    with capture_cache_log() as buf:
        r1 = bar1()
    q.put(buf.getvalue())
    q.put(r1)
    with capture_cache_log() as buf:
        r2 = bar2()
    q.put(buf.getvalue())
    q.put(r2)


class TestCacheMultipleFilesWithSignature(unittest.TestCase):
    # Regression test for https://github.com/numba/numba/issues/3658

    _numba_parallel_test_ = False

    source_text_file1 = """
from file2 import function2
"""
    source_text_file2 = """
from numba import njit

@njit('float64(float64)', cache=True)
def function1(x):
    return x

@njit('float64(float64)', cache=True)
def function2(x):
    return x
"""

    def setUp(self):
        self.tempdir = temp_directory('test_cache_file_loc')

        self.file1 = os.path.join(self.tempdir, 'file1.py')
        with open(self.file1, 'w') as fout:
            print(self.source_text_file1, file=fout)

        self.file2 = os.path.join(self.tempdir, 'file2.py')
        with open(self.file2, 'w') as fout:
            print(self.source_text_file2, file=fout)

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_caching_mutliple_files_with_signature(self):
        # Execute file1.py
        popen = subprocess.Popen([sys.executable, self.file1],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
        out, err = popen.communicate()
        msg = f"stdout:\n{out.decode()}\n\nstderr:\n{err.decode()}"
        self.assertEqual(popen.returncode, 0, msg=msg)

        # Execute file2.py
        popen = subprocess.Popen([sys.executable, self.file2],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
        out, err = popen.communicate()
        msg = f"stdout:\n{out.decode()}\n\nstderr:\n{err.decode()}"
        self.assertEqual(popen.returncode, 0, msg)


class TestCFuncCache(BaseCacheTest):

    here = os.path.dirname(__file__)
    usecases_file = os.path.join(here, "cfunc_cache_usecases.py")
    modname = "cfunc_caching_test_fodder"

    def run_in_separate_process(self):
        # Cached functions can be run from a distinct process.
        code = """if 1:
            import sys

            sys.path.insert(0, %(tempdir)r)
            mod = __import__(%(modname)r)
            mod.self_test()

            f = mod.add_usecase
            assert f.cache_hits == 1
            f = mod.outer
            assert f.cache_hits == 1
            f = mod.div_usecase
            assert f.cache_hits == 1
            """ % dict(tempdir=self.tempdir, modname=self.modname)

        popen = subprocess.Popen([sys.executable, "-c", code],
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = popen.communicate()
        if popen.returncode != 0:
            raise AssertionError(f"process failed with code {popen.returncode}:"
                                 f"stderr follows\n{err.decode()}\n")

    def check_module(self, mod):
        mod.self_test()

    def test_caching(self):
        self.check_pycache(0)
        mod = self.import_module()
        self.check_pycache(6)  # 3 index, 3 data

        self.assertEqual(mod.add_usecase.cache_hits, 0)
        self.assertEqual(mod.outer.cache_hits, 0)
        self.assertEqual(mod.add_nocache_usecase.cache_hits, 0)
        self.assertEqual(mod.div_usecase.cache_hits, 0)
        self.check_module(mod)

        # Reload module to hit the cache
        mod = self.import_module()
        self.check_pycache(6)  # 3 index, 3 data

        self.assertEqual(mod.add_usecase.cache_hits, 1)
        self.assertEqual(mod.outer.cache_hits, 1)
        self.assertEqual(mod.add_nocache_usecase.cache_hits, 0)
        self.assertEqual(mod.div_usecase.cache_hits, 1)
        self.check_module(mod)

        self.run_in_separate_process()


if __name__ == '__main__':
    unittest.main()
