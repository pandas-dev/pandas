# -*- coding: utf-8 -*-

"""
Tests the parallel backend
"""
import faulthandler
import itertools
import multiprocessing
import os
import random
import re
import subprocess
import sys
import textwrap
import threading
import unittest

import numpy as np

from numba import jit, vectorize, guvectorize, set_num_threads
from numba.tests.support import (temp_directory, override_config, TestCase, tag,
                                 skip_parfors_unsupported, linux_only)

import queue as t_queue
from numba.testing.main import _TIMEOUT as _RUNNER_TIMEOUT
from numba.core import config


_TEST_TIMEOUT = _RUNNER_TIMEOUT - 60.


# Check which backends are available
# TODO: Put this in a subprocess so the address space is kept clean
try:
    # Check it's a compatible TBB before loading it
    from numba.np.ufunc.parallel import _check_tbb_version_compatible
    _check_tbb_version_compatible()
    from numba.np.ufunc import tbbpool    # noqa: F401
    _HAVE_TBB_POOL = True
except ImportError:
    _HAVE_TBB_POOL = False

try:
    from numba.np.ufunc import omppool
    _HAVE_OMP_POOL = True
except ImportError:
    _HAVE_OMP_POOL = False

try:
    import scipy.linalg.cython_lapack    # noqa: F401
    _HAVE_LAPACK = True
except ImportError:
    _HAVE_LAPACK = False

# test skipping decorators
skip_no_omp = unittest.skipUnless(_HAVE_OMP_POOL, "OpenMP threadpool required")
skip_no_tbb = unittest.skipUnless(_HAVE_TBB_POOL, "TBB threadpool required")

_gnuomp = _HAVE_OMP_POOL and omppool.openmp_vendor == "GNU"
skip_unless_gnu_omp = unittest.skipUnless(_gnuomp, "GNU OpenMP only tests")

_windows = sys.platform.startswith('win')
_osx = sys.platform.startswith('darwin')
_32bit = sys.maxsize <= 2 ** 32
_parfors_unsupported = _32bit

_HAVE_OS_FORK = not _windows


# some functions to jit

def foo(n, v):
    return np.ones(n) + v


if _HAVE_LAPACK:
    def linalg(n, v):
        x = np.dot(np.ones((n, n)), np.ones((n, n)))
        return x + np.arange(n) + v
else:
    def linalg(n, v):
        # no way to trigger MKL without the lapack bindings.
        return np.arange(n) + v


def ufunc_foo(a, b):
    return a + b


def gufunc_foo(a, b, out):
    out[0] = a + b


class runnable(object):
    def __init__(self, **options):
        self._options = options


class jit_runner(runnable):

    def __call__(self):
        cfunc = jit(**self._options)(foo)
        a = 4
        b = 10
        expected = foo(a, b)
        got = cfunc(a, b)
        np.testing.assert_allclose(expected, got)


class mask_runner(object):
    def __init__(self, runner, mask, **options):
        self.runner = runner
        self.mask = mask

    def __call__(self):
        if self.mask:
            # Tests are all run in isolated subprocesses, so we
            # don't have to worry about this affecting other tests
            set_num_threads(self.mask)
        self.runner()


class linalg_runner(runnable):

    def __call__(self):
        cfunc = jit(**self._options)(linalg)
        a = 4
        b = 10
        expected = linalg(a, b)
        got = cfunc(a, b)
        np.testing.assert_allclose(expected, got)


class vectorize_runner(runnable):

    def __call__(self):
        cfunc = vectorize(['(f4, f4)'], **self._options)(ufunc_foo)
        a = b = np.random.random(10).astype(np.float32)
        expected = ufunc_foo(a, b)
        got = cfunc(a, b)
        np.testing.assert_allclose(expected, got)


class guvectorize_runner(runnable):

    def __call__(self):
        sig = ['(f4, f4, f4[:])']
        cfunc = guvectorize(sig, '(),()->()', **self._options)(gufunc_foo)
        a = b = np.random.random(10).astype(np.float32)
        expected = ufunc_foo(a, b)
        got = cfunc(a, b)
        np.testing.assert_allclose(expected, got)


def chooser(fnlist, **kwargs):
    q = kwargs.get('queue')
    try:
        faulthandler.enable()
        for _ in range(int(len(fnlist) * 1.5)):
            fn = random.choice(fnlist)
            fn()
    except Exception as e:
        q.put(e)


def compile_factory(parallel_class, queue_impl):
    def run_compile(fnlist):
        q = queue_impl()
        kws = {'queue': q}
        ths = [parallel_class(target=chooser, args=(fnlist,), kwargs=kws)
               for i in range(4)]
        for th in ths:
            th.start()
        for th in ths:
            th.join()
        if not q.empty():
            errors = []
            while not q.empty():
                errors.append(q.get(False))
            _msg = "Error(s) occurred in delegated runner:\n%s"
            raise RuntimeError(_msg % '\n'.join([repr(x) for x in errors]))
    return run_compile


# workers
_thread_class = threading.Thread


class _proc_class_impl(object):

    def __init__(self, method):
        self._method = method

    def __call__(self, *args, **kwargs):
        ctx = multiprocessing.get_context(self._method)
        return ctx.Process(*args, **kwargs)


def _get_mp_classes(method):
    if method == 'default':
        method = None
    ctx = multiprocessing.get_context(method)
    proc = _proc_class_impl(method)
    queue = ctx.Queue
    return proc, queue


thread_impl = compile_factory(_thread_class, t_queue.Queue)
spawn_proc_impl = compile_factory(*_get_mp_classes('spawn'))
if not _windows:
    fork_proc_impl = compile_factory(*_get_mp_classes('fork'))
    forkserver_proc_impl = compile_factory(*_get_mp_classes('forkserver'))

# this is duplication as Py27, linux uses fork, windows uses spawn, it however
# is kept like this so that when tests fail it's less confusing!
default_proc_impl = compile_factory(*_get_mp_classes('default'))


class TestParallelBackendBase(TestCase):
    """
    Base class for testing the parallel backends
    """

    all_impls = [
        jit_runner(nopython=True),
        jit_runner(nopython=True, cache=True),
        jit_runner(nopython=True, nogil=True),
        linalg_runner(nopython=True),
        linalg_runner(nopython=True, nogil=True),
        vectorize_runner(nopython=True),
        vectorize_runner(nopython=True, target='parallel'),
        vectorize_runner(nopython=True, target='parallel', cache=True),
        guvectorize_runner(nopython=True),
        guvectorize_runner(nopython=True, target='parallel'),
        guvectorize_runner(nopython=True, target='parallel', cache=True),
    ]

    if not _parfors_unsupported:
        parfor_impls = [
            jit_runner(nopython=True, parallel=True),
            jit_runner(nopython=True, parallel=True, cache=True),
            linalg_runner(nopython=True, parallel=True),
            linalg_runner(nopython=True, parallel=True, cache=True),
        ]
        all_impls.extend(parfor_impls)

    if config.NUMBA_NUM_THREADS < 2:
        # Not enough cores
        masks = []
    else:
        masks = [1, 2]

    mask_impls = []
    for impl in all_impls:
        for mask in masks:
            mask_impls.append(mask_runner(impl, mask))

    parallelism = ['threading', 'random']
    parallelism.append('multiprocessing_spawn')
    if _HAVE_OS_FORK:
        parallelism.append('multiprocessing_fork')
        parallelism.append('multiprocessing_forkserver')

    runners = {
        'concurrent_jit': [
            jit_runner(nopython=True, parallel=(not _parfors_unsupported)),
        ],
        'concurrent_vectorize': [
            vectorize_runner(nopython=True, target='parallel'),
        ],
        'concurrent_guvectorize': [
            guvectorize_runner(nopython=True, target='parallel'),
        ],
        'concurrent_mix_use': all_impls,
        'concurrent_mix_use_masks': mask_impls,
    }

    safe_backends = {'omp', 'tbb'}

    def run_compile(self, fnlist, parallelism='threading'):
        self._cache_dir = temp_directory(self.__class__.__name__)
        with override_config('CACHE_DIR', self._cache_dir):
            if parallelism == 'threading':
                thread_impl(fnlist)
            elif parallelism == 'multiprocessing_fork':
                fork_proc_impl(fnlist)
            elif parallelism == 'multiprocessing_forkserver':
                forkserver_proc_impl(fnlist)
            elif parallelism == 'multiprocessing_spawn':
                spawn_proc_impl(fnlist)
            elif parallelism == 'multiprocessing_default':
                default_proc_impl(fnlist)
            elif parallelism == 'random':
                ps = [thread_impl, spawn_proc_impl]
                if _HAVE_OS_FORK:
                    ps.append(fork_proc_impl)
                    ps.append(forkserver_proc_impl)

                random.shuffle(ps)
                for impl in ps:
                    impl(fnlist)
            else:
                raise ValueError(
                    'Unknown parallelism supplied %s' % parallelism)


_specific_backends = config.THREADING_LAYER in ('omp', 'tbb', 'workqueue')


@unittest.skipUnless(_specific_backends, "Threading layer not explicit")
class TestParallelBackend(TestParallelBackendBase):
    """ These are like the numba.tests.test_threadsafety tests but designed
    instead to torture the parallel backend.
    If a suitable backend is supplied via NUMBA_THREADING_LAYER these tests
    can be run directly. This test class cannot be run using the multiprocessing
    option to the test runner (i.e. `./runtests -m`) as daemon processes cannot
    have children.
    """

    # NOTE: All tests are generated based on what a platform supports concurrent
    # execution wise from Python, irrespective of whether the native libraries
    # can actually handle the behaviour present.
    @classmethod
    def generate(cls):
        for p in cls.parallelism:
            for name, impl in cls.runners.items():
                methname = "test_" + p + '_' + name

                def methgen(impl, p):
                    def test_method(self):
                        selfproc = multiprocessing.current_process()
                        # daemonized processes cannot have children
                        if selfproc.daemon:
                            _msg = 'daemonized processes cannot have children'
                            self.skipTest(_msg)
                        else:
                            self.run_compile(impl, parallelism=p)
                    return test_method
                fn = methgen(impl, p)
                fn.__name__ = methname
                setattr(cls, methname, fn)


TestParallelBackend.generate()


class TestInSubprocess(object):
    backends = {'tbb': skip_no_tbb,
                'omp': skip_no_omp,
                'workqueue': unittest.skipIf(False, '')}

    def run_cmd(self, cmdline, env):
        popen = subprocess.Popen(cmdline,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 env=env)
        # finish in _TEST_TIMEOUT seconds or kill it
        timeout = threading.Timer(_TEST_TIMEOUT, popen.kill)
        try:
            timeout.start()
            out, err = popen.communicate()
            if popen.returncode != 0:
                raise AssertionError(
                    "process failed with code %s: stderr follows\n%s\n" %
                    (popen.returncode, err.decode()))
            return out.decode(), err.decode()
        finally:
            timeout.cancel()
        return None, None

    def run_test_in_separate_process(self, test, threading_layer):
        env_copy = os.environ.copy()
        env_copy['NUMBA_THREADING_LAYER'] = str(threading_layer)
        cmdline = [sys.executable, "-m", "numba.runtests", test]
        return self.run_cmd(cmdline, env_copy)


class TestSpecificBackend(TestInSubprocess, TestParallelBackendBase):
    """
    This is quite contrived, for each test in the TestParallelBackend tests it
    generates a test that will run the TestParallelBackend test in a new python
    process with an environment modified to ensure a specific threadsafe backend
    is used. This is with view of testing the backends independently and in an
    isolated manner such that if they hang/crash/have issues, it doesn't kill
    the test suite.
    """
    _DEBUG = False

    @classmethod
    def _inject(cls, p, name, backend, backend_guard):
        themod = cls.__module__
        thecls = TestParallelBackend.__name__
        methname = "test_" + p + '_' + name
        injected_method = '%s.%s.%s' % (themod, thecls, methname)

        def test_template(self):
            o, e = self.run_test_in_separate_process(injected_method, backend)
            if self._DEBUG:
                print('stdout:\n "%s"\n stderr:\n "%s"' % (o, e))
            # If the test was skipped in the subprocess, then mark this as a
            # skipped test.
            m = re.search(r"\.\.\. skipped '(.*?)'", e)
            if m is not None:
                self.skipTest(m.group(1))
            self.assertIn('OK', e)
            self.assertTrue('FAIL' not in e)
            self.assertTrue('ERROR' not in e)
        injected_test = "test_%s_%s_%s" % (p, name, backend)
        # Mark as long_running
        setattr(cls, injected_test,
                tag('long_running')(backend_guard(test_template)))

    @classmethod
    def generate(cls):
        for backend, backend_guard in cls.backends.items():
            for p in cls.parallelism:
                for name in cls.runners.keys():
                    # handle known problem cases...

                    # GNU OpenMP is not fork safe
                    if (p in ('multiprocessing_fork', 'random') and
                        backend == 'omp' and
                            sys.platform.startswith('linux')):
                        continue

                    # workqueue is not thread safe
                    if (p in ('threading', 'random') and
                            backend == 'workqueue'):
                        continue

                    cls._inject(p, name, backend, backend_guard)


TestSpecificBackend.generate()


class ThreadLayerTestHelper(TestCase):
    """
    Helper class for running an isolated piece of code based on a template
    """
    # sys path injection and separate usecase module to make sure everything
    # is importable by children of multiprocessing
    _here = "%r" % os.path.dirname(__file__)

    template = """if 1:
    import sys
    sys.path.insert(0, "%(here)r")
    import multiprocessing
    import numpy as np
    from numba import njit
    import numba
    try:
        import threading_backend_usecases
    except ImportError as e:
        print("DEBUG:", sys.path)
        raise e
    import os

    sigterm_handler = threading_backend_usecases.sigterm_handler
    busy_func = threading_backend_usecases.busy_func

    def the_test():
        %%s

    if __name__ == "__main__":
        the_test()
    """ % {'here': _here}

    def run_cmd(self, cmdline, env=None):
        if env is None:
            env = os.environ.copy()
            env['NUMBA_THREADING_LAYER'] = str("omp")
        popen = subprocess.Popen(cmdline,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 env=env)
        # finish in _TEST_TIMEOUT seconds or kill it
        timeout = threading.Timer(_TEST_TIMEOUT, popen.kill)
        try:
            timeout.start()
            out, err = popen.communicate()
            if popen.returncode != 0:
                raise AssertionError(
                    "process failed with code %s: stderr follows\n%s\n" %
                    (popen.returncode, err.decode()))
        finally:
            timeout.cancel()
        return out.decode(), err.decode()


@skip_parfors_unsupported
class TestThreadingLayerSelection(ThreadLayerTestHelper):
    """
    Checks that numba.threading_layer() reports correctly.
    """
    _DEBUG = False

    backends = {'tbb': skip_no_tbb,
                'omp': skip_no_omp,
                'workqueue': unittest.skipIf(False, '')}

    @classmethod
    def _inject(cls, backend, backend_guard):

        def test_template(self):
            body = """if 1:
                X = np.arange(1000000.)
                Y = np.arange(1000000.)
                Z = busy_func(X, Y)
                assert numba.threading_layer() == '%s'
            """
            runme = self.template % (body % backend)
            cmdline = [sys.executable, '-c', runme]
            env = os.environ.copy()
            env['NUMBA_THREADING_LAYER'] = str(backend)
            out, err = self.run_cmd(cmdline, env=env)
            if self._DEBUG:
                print(out, err)
        injected_test = "test_threading_layer_selector_%s" % backend
        setattr(cls, injected_test,
                tag("important")(backend_guard(test_template)))

    @classmethod
    def generate(cls):
        for backend, backend_guard in cls.backends.items():
            cls._inject(backend, backend_guard)


TestThreadingLayerSelection.generate()


@skip_parfors_unsupported
class TestThreadingLayerPriority(ThreadLayerTestHelper):

    def each_env_var(self, env_var: str):
        """Test setting priority via env var NUMBA_THREADING_LAYER_PRIORITY.
        """
        env = os.environ.copy()
        env['NUMBA_THREADING_LAYER'] = 'default'
        env['NUMBA_THREADING_LAYER_PRIORITY'] = env_var

        code = f"""
                import numba

                # trigger threading layer decision
                # hence catching invalid THREADING_LAYER_PRIORITY
                @numba.jit(
                    'float64[::1](float64[::1], float64[::1])',
                    nopython=True,
                    parallel=True,
                )
                def plus(x, y):
                    return x + y

                captured_envvar = list("{env_var}".split())
                assert numba.config.THREADING_LAYER_PRIORITY == \
                    captured_envvar, "priority mismatch"
                assert numba.threading_layer() == captured_envvar[0],\
                    "selected backend mismatch"
                """
        cmd = [
            sys.executable,
            '-c',
            textwrap.dedent(code),
        ]
        self.run_cmd(cmd, env=env)

    @skip_no_omp
    @skip_no_tbb
    def test_valid_env_var(self):
        default = ['tbb', 'omp', 'workqueue']
        for p in itertools.permutations(default):
            env_var = ' '.join(p)
            self.each_env_var(env_var)

    @skip_no_omp
    @skip_no_tbb
    def test_invalid_env_var(self):
        env_var = 'tbb omp workqueue notvalidhere'
        with self.assertRaises(AssertionError) as raises:
            self.each_env_var(env_var)
        for msg in (
            "THREADING_LAYER_PRIORITY invalid:",
            "It must be a permutation of"
        ):
            self.assertIn(f"{msg}", str(raises.exception))

    @skip_no_omp
    def test_omp(self):
        for env_var in ("omp tbb workqueue", "omp workqueue tbb"):
            self.each_env_var(env_var)

    @skip_no_tbb
    def test_tbb(self):
        for env_var in ("tbb omp workqueue", "tbb workqueue omp"):
            self.each_env_var(env_var)

    def test_workqueue(self):
        for env_var in ("workqueue tbb omp", "workqueue omp tbb"):
            self.each_env_var(env_var)


@skip_parfors_unsupported
class TestMiscBackendIssues(ThreadLayerTestHelper):
    """
    Checks fixes for the issues with threading backends implementation
    """
    _DEBUG = False

    @skip_no_omp
    def test_omp_stack_overflow(self):
        """
        Tests that OMP does not overflow stack
        """
        runme = """if 1:
            from numba import vectorize, threading_layer
            import numpy as np

            @vectorize(['f4(f4,f4,f4,f4,f4,f4,f4,f4)'], target='parallel')
            def foo(a, b, c, d, e, f, g, h):
                return a+b+c+d+e+f+g+h

            x = np.ones(2**20, np.float32)
            foo(*([x]*8))
            assert threading_layer() == "omp", "omp not found"
        """
        cmdline = [sys.executable, '-c', runme]
        env = os.environ.copy()
        env['NUMBA_THREADING_LAYER'] = "omp"
        env['OMP_STACKSIZE'] = "100K"
        self.run_cmd(cmdline, env=env)

    @skip_no_tbb
    def test_single_thread_tbb(self):
        """
        Tests that TBB works well with single thread
        https://github.com/numba/numba/issues/3440
        """
        runme = """if 1:
            from numba import njit, prange, threading_layer

            @njit(parallel=True)
            def foo(n):
                acc = 0
                for i in prange(n):
                    acc += i
                return acc

            foo(100)
            assert threading_layer() == "tbb", "tbb not found"
        """
        cmdline = [sys.executable, '-c', runme]
        env = os.environ.copy()
        env['NUMBA_THREADING_LAYER'] = "tbb"
        env['NUMBA_NUM_THREADS'] = "1"
        self.run_cmd(cmdline, env=env)

    def test_workqueue_aborts_on_nested_parallelism(self):
        """
        Tests workqueue raises sigabrt if a nested parallel call is performed
        """
        runme = """if 1:
            from numba import njit, prange
            import numpy as np

            @njit(parallel=True)
            def nested(x):
                for i in prange(len(x)):
                    x[i] += 1


            @njit(parallel=True)
            def main():
                Z = np.zeros((5, 10))
                for i in prange(Z.shape[0]):
                    nested(Z[i])
                return Z

            main()
        """
        cmdline = [sys.executable, '-c', runme]
        env = os.environ.copy()
        env['NUMBA_THREADING_LAYER'] = "workqueue"
        env['NUMBA_NUM_THREADS'] = "4"

        try:
            out, err = self.run_cmd(cmdline, env=env)
        except AssertionError as e:
            if self._DEBUG:
                print(out, err)
            e_msg = str(e)
            self.assertIn("failed with code", e_msg)
            # raised a SIGABRT, but the value is platform specific so just check
            # the error message
            expected = ("Numba workqueue threading layer is terminating: "
                        "Concurrent access has been detected.")
            self.assertIn(expected, e_msg)

    @unittest.skipUnless(_HAVE_OS_FORK, "Test needs fork(2)")
    def test_workqueue_handles_fork_from_non_main_thread(self):
        # For context see #7872, but essentially the multiprocessing pool
        # implementation has a number of Python threads for handling the worker
        # processes, one of which calls fork(2), this results in a fork from a
        # non-main thread.

        runme = """if 1:
            from numba import njit, prange, threading_layer
            import numpy as np
            import multiprocessing

            if __name__ == "__main__":
                # Need for force fork context (OSX default is "spawn")
                multiprocessing.set_start_method('fork')

                @njit(parallel=True)
                def func(x):
                    return 10. * x

                arr = np.arange(2.)

                # run in single process to start Numba's thread pool
                np.testing.assert_allclose(func(arr), func.py_func(arr))

                # now run in a multiprocessing pool to get a fork from a
                # non-main thread
                with multiprocessing.Pool(10) as p:
                    result = p.map(func, [arr])
                np.testing.assert_allclose(result,
                                           func.py_func(np.expand_dims(arr, 0)))

                assert threading_layer() == "workqueue"
        """
        cmdline = [sys.executable, '-c', runme]
        env = os.environ.copy()
        env['NUMBA_THREADING_LAYER'] = "workqueue"
        env['NUMBA_NUM_THREADS'] = "4"

        self.run_cmd(cmdline, env=env)


# 32bit or windows py27 (not that this runs on windows)
@skip_parfors_unsupported
@skip_unless_gnu_omp
class TestForkSafetyIssues(ThreadLayerTestHelper):
    """
    Checks Numba's behaviour in various situations involving GNU OpenMP and fork
    """
    _DEBUG = False

    def test_check_threading_layer_is_gnu(self):
        runme = """if 1:
            from numba.np.ufunc import omppool
            assert omppool.openmp_vendor == 'GNU'
            """
        cmdline = [sys.executable, '-c', runme]
        out, err = self.run_cmd(cmdline)

    def test_par_parent_os_fork_par_child(self):
        """
        Whilst normally valid, this actually isn't for Numba invariant of OpenMP
        Checks SIGABRT is received.
        """
        body = """if 1:
            X = np.arange(1000000.)
            Y = np.arange(1000000.)
            Z = busy_func(X, Y)
            pid = os.fork()
            if pid  == 0:
                Z = busy_func(X, Y)
            else:
                os.wait()
        """
        runme = self.template % body
        cmdline = [sys.executable, '-c', runme]
        try:
            out, err = self.run_cmd(cmdline)
        except AssertionError as e:
            self.assertIn("failed with code -6", str(e))

    def test_par_parent_implicit_mp_fork_par_child(self):
        """
        Implicit use of multiprocessing fork context.
        Does this:
        1. Start with OpenMP
        2. Fork to processes using OpenMP (this is invalid)
        3. Joins fork
        4. Check the exception pushed onto the queue that is a result of
           catching SIGTERM coming from the C++ aborting on illegal fork
           pattern for GNU OpenMP
        """
        body = """if 1:
            mp = multiprocessing.get_context('fork')
            X = np.arange(1000000.)
            Y = np.arange(1000000.)
            q = mp.Queue()

            # Start OpenMP runtime on parent via parallel function
            Z = busy_func(X, Y, q)

            # fork() underneath with no exec, will abort
            proc = mp.Process(target = busy_func, args=(X, Y, q))
            proc.start()

            err = q.get()
            assert "Caught SIGTERM" in str(err)
        """
        runme = self.template % body
        cmdline = [sys.executable, '-c', runme]
        out, err = self.run_cmd(cmdline)
        if self._DEBUG:
            print(out, err)

    @linux_only
    def test_par_parent_explicit_mp_fork_par_child(self):
        """
        Explicit use of multiprocessing fork context.
        Does this:
        1. Start with OpenMP
        2. Fork to processes using OpenMP (this is invalid)
        3. Joins fork
        4. Check the exception pushed onto the queue that is a result of
           catching SIGTERM coming from the C++ aborting on illegal fork
           pattern for GNU OpenMP
        """
        body = """if 1:
            X = np.arange(1000000.)
            Y = np.arange(1000000.)
            ctx = multiprocessing.get_context('fork')
            q = ctx.Queue()

            # Start OpenMP runtime on parent via parallel function
            Z = busy_func(X, Y, q)

            # fork() underneath with no exec, will abort
            proc = ctx.Process(target = busy_func, args=(X, Y, q))
            proc.start()
            proc.join()

            err = q.get()
            assert "Caught SIGTERM" in str(err)
        """
        runme = self.template % body
        cmdline = [sys.executable, '-c', runme]
        out, err = self.run_cmd(cmdline)
        if self._DEBUG:
            print(out, err)

    def test_par_parent_mp_spawn_par_child_par_parent(self):
        """
        Explicit use of multiprocessing spawn, this is safe.
        Does this:
        1. Start with OpenMP
        2. Spawn to processes using OpenMP
        3. Join spawns
        4. Run some more OpenMP
        """
        body = """if 1:
            X = np.arange(1000000.)
            Y = np.arange(1000000.)
            ctx = multiprocessing.get_context('spawn')
            q = ctx.Queue()

            # Start OpenMP runtime and run on parent via parallel function
            Z = busy_func(X, Y, q)
            procs = []
            for x in range(20): # start a lot to try and get overlap
                ## fork() + exec() to run some OpenMP on children
                proc = ctx.Process(target = busy_func, args=(X, Y, q))
                procs.append(proc)
                sys.stdout.flush()
                sys.stderr.flush()
                proc.start()

            [p.join() for p in procs]

            try:
                q.get(False)
            except multiprocessing.queues.Empty:
                pass
            else:
                raise RuntimeError("Queue was not empty")

            # Run some more OpenMP on parent
            Z = busy_func(X, Y, q)
        """
        runme = self.template % body
        cmdline = [sys.executable, '-c', runme]
        out, err = self.run_cmd(cmdline)
        if self._DEBUG:
            print(out, err)

    def test_serial_parent_implicit_mp_fork_par_child_then_par_parent(self):
        """
        Implicit use of multiprocessing (will be fork, but cannot declare that
        in Py2.7 as there's no process launch context).
        Does this:
        1. Start with no OpenMP
        2. Fork to processes using OpenMP
        3. Join forks
        4. Run some OpenMP
        """
        body = """if 1:
            X = np.arange(1000000.)
            Y = np.arange(1000000.)
            q = multiprocessing.Queue()

            # this is ok
            procs = []
            for x in range(10):
                # fork() underneath with but no OpenMP in parent, this is ok
                proc = multiprocessing.Process(target = busy_func,
                                               args=(X, Y, q))
                procs.append(proc)
                proc.start()

            [p.join() for p in procs]

            # and this is still ok as the OpenMP happened in forks
            Z = busy_func(X, Y, q)
            try:
                q.get(False)
            except multiprocessing.queues.Empty:
                pass
            else:
                raise RuntimeError("Queue was not empty")
        """
        runme = self.template % body
        cmdline = [sys.executable, '-c', runme]
        out, err = self.run_cmd(cmdline)
        if self._DEBUG:
            print(out, err)

    @linux_only
    def test_serial_parent_explicit_mp_fork_par_child_then_par_parent(self):
        """
        Explicit use of multiprocessing 'fork'.
        Does this:
        1. Start with no OpenMP
        2. Fork to processes using OpenMP
        3. Join forks
        4. Run some OpenMP
        """
        body = """if 1:
            X = np.arange(1000000.)
            Y = np.arange(1000000.)
            ctx = multiprocessing.get_context('fork')
            q = ctx.Queue()

            # this is ok
            procs = []
            for x in range(10):
                # fork() underneath with but no OpenMP in parent, this is ok
                proc = ctx.Process(target = busy_func, args=(X, Y, q))
                procs.append(proc)
                proc.start()

            [p.join() for p in procs]

            # and this is still ok as the OpenMP happened in forks
            Z = busy_func(X, Y, q)
            try:
                q.get(False)
            except multiprocessing.queues.Empty:
                pass
            else:
                raise RuntimeError("Queue was not empty")
        """
        runme = self.template % body
        cmdline = [sys.executable, '-c', runme]
        out, err = self.run_cmd(cmdline)
        if self._DEBUG:
            print(out, err)


@skip_parfors_unsupported
@skip_no_tbb
class TestTBBSpecificIssues(ThreadLayerTestHelper):

    _DEBUG = False

    @linux_only # os.fork required.
    def test_fork_from_non_main_thread(self):
        # See issue #5973 and PR #6208 for original context.
        # See issue #6963 for context on the following comments:
        #
        # Important things to note:
        # 1. Compilation of code containing an objmode block will result in the
        #    use of and `ObjModeLiftedWith` as the dispatcher. This inherits
        #    from `LiftedCode` which handles the serialization. In that
        #    serialization is a call to uuid.uuid1() which causes a fork_exec in
        #    CPython internals.
        # 2. The selected parallel backend thread pool is started during the
        #    compilation of a function that has `parallel=True`.
        # 3. The TBB backend can handle forks from the main thread, it will
        #    safely reinitialise after so doing. If a fork occurs from a
        #    non-main thread it will warn and the state is invalid in the child
        #    process.
        #
        # Due to 1. and 2. the `obj_mode_func` function separated out and is
        # `njit` decorated. This means during type inference of `work` it will
        # trigger a standard compilation of the function and the thread pools
        # won't have started yet as the parallelisation compiler passes for
        # `work` won't yet have run. This mitigates the fork() call from 1.
        # occurring after 2. The result of this is that 3. can be tested using
        # the threading etc herein with the state being known as the above
        # described, i.e. the TBB threading layer has not experienced a fork().

        runme = """if 1:
            import threading
            import numba
            numba.config.THREADING_LAYER='tbb'
            from numba import njit, prange, objmode
            from numba.core.serialize import PickleCallableByPath
            import os

            e_running = threading.Event()
            e_proceed = threading.Event()

            def indirect_core():
                e_running.set()
                # wait for forker() to have forked
                while not e_proceed.isSet():
                    pass

            indirect = PickleCallableByPath(indirect_core)

            @njit
            def obj_mode_func():
                with objmode():
                    indirect()

            @njit(parallel=True, nogil=True)
            def work():
                acc = 0
                for x in prange(10):
                    acc += x
                obj_mode_func()
                return acc

            def runner():
                work()

            def forker():
                # wait for the jit function to say it's running
                while not e_running.isSet():
                    pass
                # then fork
                os.fork()
                # now fork is done signal the runner to proceed to exit
                e_proceed.set()

            numba_runner = threading.Thread(target=runner,)
            fork_runner =  threading.Thread(target=forker,)

            threads = (numba_runner, fork_runner)
            for t in threads:
                t.start()
            for t in threads:
                t.join()
        """

        cmdline = [sys.executable, '-c', runme]
        out, err = self.run_cmd(cmdline)
        # assert error message printed on stderr
        msg_head = "Attempted to fork from a non-main thread, the TBB library"
        self.assertIn(msg_head, err)

        if self._DEBUG:
            print("OUT:", out)
            print("ERR:", err)

    @linux_only # fork required.
    def test_lifetime_of_task_scheduler_handle(self):

        self.skip_if_no_external_compiler() # external compiler needed

        # See PR #7280 for context.
        BROKEN_COMPILERS = 'SKIP: COMPILATION FAILED'
        runme = """if 1:
            import ctypes
            import sys
            import multiprocessing as mp
            from tempfile import TemporaryDirectory, NamedTemporaryFile
            from numba.pycc.platform import Toolchain, external_compiler_works
            from numba import njit, prange, threading_layer
            import faulthandler
            faulthandler.enable()
            if not external_compiler_works():
                raise AssertionError('External compilers are not found.')
            with TemporaryDirectory() as tmpdir:
                with NamedTemporaryFile(dir=tmpdir) as tmpfile:
                    try:
                        src = \"\"\"
                        #define TBB_PREVIEW_WAITING_FOR_WORKERS 1
                        #include <tbb/tbb.h>
                        static tbb::task_scheduler_handle tsh;
                        extern "C"
                        {
                        void launch(void)
                        {
                            tsh = tbb::task_scheduler_handle::get();
                        }
                        }
                        \"\"\"
                        cxxfile = f"{tmpfile.name}.cxx"
                        with open(cxxfile, 'wt') as f:
                            f.write(src)
                        tc = Toolchain()
                        object_files = tc.compile_objects([cxxfile,],
                                                           output_dir=tmpdir)
                        dso_name = f"{tmpfile.name}.so"
                        tc.link_shared(dso_name, object_files,
                                       libraries=['tbb',],
                                       export_symbols=['launch'])
                        # Load into the process, it doesn't matter whether the
                        # DSO exists on disk once it's loaded in.
                        DLL = ctypes.CDLL(dso_name)
                    except Exception as e:
                        # Something is broken in compilation, could be one of
                        # many things including, but not limited to: missing tbb
                        # headers, incorrect permissions, compilers that don't
                        # work for the above
                        print(e)
                        print('BROKEN_COMPILERS')
                        sys.exit(0)

                    # Do the test, launch this library and also execute a
                    # function with the TBB threading layer.

                    DLL.launch()

                    @njit(parallel=True)
                    def foo(n):
                        acc = 0
                        for i in prange(n):
                            acc += i
                        return acc

                    foo(1)

            # Check the threading layer used was TBB
            assert threading_layer() == 'tbb'

            # Use mp context for a controlled version of fork, this triggers the
            # reported bug.

            ctx = mp.get_context('fork')
            def nowork():
                pass
            p = ctx.Process(target=nowork)
            p.start()
            p.join(10)
            print("SUCCESS")
            """.replace('BROKEN_COMPILERS', BROKEN_COMPILERS)

        cmdline = [sys.executable, '-c', runme]
        env = os.environ.copy()
        env['NUMBA_THREADING_LAYER'] = 'tbb'
        out, err = self.run_cmd(cmdline, env=env)

        if BROKEN_COMPILERS in out:
            self.skipTest("Compilation of DSO failed. Check output for details")
        else:
            self.assertIn("SUCCESS", out)

        if self._DEBUG:
            print("OUT:", out)
            print("ERR:", err)


@skip_parfors_unsupported
class TestInitSafetyIssues(TestCase):

    _DEBUG = False

    def run_cmd(self, cmdline):
        popen = subprocess.Popen(cmdline,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,)
        # finish in _TEST_TIMEOUT seconds or kill it
        timeout = threading.Timer(_TEST_TIMEOUT, popen.kill)
        try:
            timeout.start()
            out, err = popen.communicate()
            if popen.returncode != 0:
                raise AssertionError(
                    "process failed with code %s: stderr follows\n%s\n" %
                    (popen.returncode, err.decode()))
        finally:
            timeout.cancel()
        return out.decode(), err.decode()

    @linux_only # only linux can leak semaphores
    def test_orphaned_semaphore(self):
        # sys path injection and separate usecase module to make sure everything
        # is importable by children of multiprocessing

        test_file = os.path.join(os.path.dirname(__file__),
                                 "orphaned_semaphore_usecase.py")
        cmdline = [sys.executable, test_file]
        out, err = self.run_cmd(cmdline)

        # assert no semaphore leaks reported on stderr
        self.assertNotIn("leaked semaphore", err)

        if self._DEBUG:
            print("OUT:", out)
            print("ERR:", err)

    def test_lazy_lock_init(self):
        # checks based on https://github.com/numba/numba/pull/5724
        # looking for "lazy" process lock initialisation so as to avoid setting
        # a multiprocessing context as part of import.
        for meth in ('fork', 'spawn', 'forkserver'):
            # if a context is available on the host check it can be set as the
            # start method in a separate process
            try:
                multiprocessing.get_context(meth)
            except ValueError:
                continue
            cmd = ("import numba; import multiprocessing;"
                   "multiprocessing.set_start_method('{}');"
                   "print(multiprocessing.get_context().get_start_method())")
            cmdline = [sys.executable, "-c", cmd.format(meth)]
            out, err = self.run_cmd(cmdline)
            if self._DEBUG:
                print("OUT:", out)
                print("ERR:", err)
            self.assertIn(meth, out)


@skip_parfors_unsupported
@skip_no_omp
class TestOpenMPVendors(TestCase):

    def test_vendors(self):
        """
        Checks the OpenMP vendor strings are correct
        """
        expected = dict()
        expected['win32'] = "MS"
        expected['darwin'] = "Intel"
        expected['linux'] = "GNU"

        # only check OS that are supported, custom toolchains may well work as
        # may other OS
        for k in expected.keys():
            if sys.platform.startswith(k):
                self.assertEqual(expected[k], omppool.openmp_vendor)


if __name__ == '__main__':
    unittest.main()
