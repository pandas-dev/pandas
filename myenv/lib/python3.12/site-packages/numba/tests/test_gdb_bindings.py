"""
Tests gdb bindings
"""
import os
import platform
import re
import subprocess
import sys
import threading
from itertools import permutations

from numba import njit, gdb, gdb_init, gdb_breakpoint, prange
from numba.core import errors
from numba import jit

from numba.tests.support import (TestCase, captured_stdout, tag,
                                 skip_parfors_unsupported)
from numba.tests.gdb_support import needs_gdb
import unittest


_platform = sys.platform

_unix_like = (_platform.startswith('linux')
              or _platform.startswith('darwin')
              or ('bsd' in _platform))

unix_only = unittest.skipUnless(_unix_like, "unix-like OS is required")
not_unix = unittest.skipIf(_unix_like, "non unix-like OS is required")

_arch_name = platform.machine()
_is_arm = _arch_name in {'aarch64', 'armv7l'}
not_arm = unittest.skipIf(_is_arm, "testing disabled on ARM")

_gdb_cond = os.environ.get('GDB_TEST', None) == '1'
needs_gdb_harness = unittest.skipUnless(_gdb_cond, "needs gdb harness")

long_running = tag('long_running')

_dbg_njit = njit(debug=True)
_dbg_jit = jit(forceobj=True, debug=True)


def impl_gdb_call(a):
    gdb('-ex', 'set confirm off', '-ex', 'c', '-ex', 'q')
    b = a + 1
    c = a * 2.34
    d = (a, b, c)
    print(a, b, c, d)


def impl_gdb_call_w_bp(a):
    gdb_init('-ex', 'set confirm off', '-ex', 'c', '-ex', 'q')
    b = a + 1
    c = a * 2.34
    d = (a, b, c)
    gdb_breakpoint()
    print(a, b, c, d)


def impl_gdb_split_init_and_break_w_parallel(a):
    gdb_init('-ex', 'set confirm off', '-ex', 'c', '-ex', 'q')
    a += 3
    for i in prange(4):
        b = a + 1
        c = a * 2.34
        d = (a, b, c)
        gdb_breakpoint()
        print(a, b, c, d)


@not_arm
@unix_only
class TestGdbBindImpls(TestCase):
    """
    Contains unit test implementations for gdb binding testing. Test must be
    decorated with `@needs_gdb_harness` to prevent their running under normal
    test conditions, the test methods must also end with `_impl` to be
    considered for execution. The tests themselves are invoked by the
    `TestGdbBinding` test class through the parsing of this class for test
    methods and then running the discovered tests in a separate process. Test
    names not including the word `quick` will be tagged as @tag('long_running')
    """

    @needs_gdb_harness
    def test_gdb_cmd_lang_cpython_quick_impl(self):
        with captured_stdout():
            impl_gdb_call(10)

    @needs_gdb_harness
    def test_gdb_cmd_lang_nopython_quick_impl(self):
        with captured_stdout():
            _dbg_njit(impl_gdb_call)(10)

    @needs_gdb_harness
    def test_gdb_cmd_lang_objmode_quick_impl(self):
        with captured_stdout():
            _dbg_jit(impl_gdb_call)(10)

    @needs_gdb_harness
    def test_gdb_split_init_and_break_cpython_impl(self):
        with captured_stdout():
            impl_gdb_call_w_bp(10)

    @needs_gdb_harness
    def test_gdb_split_init_and_break_nopython_impl(self):
        with captured_stdout():
            _dbg_njit(impl_gdb_call_w_bp)(10)

    @needs_gdb_harness
    def test_gdb_split_init_and_break_objmode_impl(self):
        with captured_stdout():
            _dbg_jit(impl_gdb_call_w_bp)(10)

    @skip_parfors_unsupported
    @needs_gdb_harness
    def test_gdb_split_init_and_break_w_parallel_cpython_impl(self):
        with captured_stdout():
            impl_gdb_split_init_and_break_w_parallel(10)

    @skip_parfors_unsupported
    @needs_gdb_harness
    def test_gdb_split_init_and_break_w_parallel_nopython_impl(self):
        with captured_stdout():
            _dbg_njit(impl_gdb_split_init_and_break_w_parallel)(10)

    @skip_parfors_unsupported
    @needs_gdb_harness
    def test_gdb_split_init_and_break_w_parallel_objmode_impl(self):
        with captured_stdout():
            _dbg_jit(impl_gdb_split_init_and_break_w_parallel)(10)


@not_arm
@unix_only
@needs_gdb
class TestGdbBinding(TestCase):
    """
    This test class is used to generate tests which will run the test cases
    defined in TestGdbBindImpls in isolated subprocesses, this is for safety
    in case something goes awry.
    """

    # test mutates env
    _numba_parallel_test_ = False

    _DEBUG = True

    def run_cmd(self, cmdline, env, kill_is_ok=False):
        popen = subprocess.Popen(cmdline,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 env=env,
                                 shell=True)
        # finish in 20s or kill it, there's no work being done

        def kill():
            popen.stdout.flush()
            popen.stderr.flush()
            popen.kill()
        timeout = threading.Timer(20., kill)
        try:
            timeout.start()
            out, err = popen.communicate()
            retcode = popen.returncode
            if retcode != 0:
                raise AssertionError(
                    "process failed with code %s: "
                    "stderr follows\n%s\n"
                    "stdout :%s" % (retcode, err.decode(), out.decode()))
            return out.decode(), err.decode()
        finally:
            timeout.cancel()
        return None, None

    def run_test_in_separate_process(self, test, **kwargs):
        env_copy = os.environ.copy()
        env_copy['NUMBA_OPT'] = '1'
        # Set GDB_TEST to permit the execution of tests decorated with
        # @needs_gdb_harness
        env_copy['GDB_TEST'] = '1'
        cmdline = [sys.executable, "-m", "numba.runtests", test]
        return self.run_cmd(' '.join(cmdline), env_copy, **kwargs)

    @classmethod
    def _inject(cls, name):
        themod = TestGdbBindImpls.__module__
        thecls = TestGdbBindImpls.__name__
        # strip impl
        assert name.endswith('_impl')
        methname = name.replace('_impl', '')
        injected_method = '%s.%s.%s' % (themod, thecls, name)

        def test_template(self):
            o, e = self.run_test_in_separate_process(injected_method)
            dbgmsg = f'\nSTDOUT={o}\nSTDERR={e}\n'
            # If the test was skipped in the subprocess, then mark this as a
            # skipped test.
            m = re.search(r"\.\.\. skipped '(.*?)'", e)
            if m is not None:
                self.skipTest(m.group(1))
            self.assertIn('GNU gdb', o, msg=dbgmsg)
            self.assertIn('OK', e, msg=dbgmsg)
            self.assertNotIn('FAIL', e, msg=dbgmsg)
            self.assertNotIn('ERROR', e, msg=dbgmsg)
        if 'quick' in name:
            setattr(cls, methname, test_template)
        else:
            setattr(cls, methname, long_running(test_template))

    @classmethod
    def generate(cls):
        for name in dir(TestGdbBindImpls):
            if name.startswith('test_gdb'):
                cls._inject(name)


TestGdbBinding.generate()


@not_arm
@unix_only
@needs_gdb
class TestGdbMisc(TestCase):

    @long_running
    def test_call_gdb_twice(self):
        def gen(f1, f2):
            @njit
            def impl():
                a = 1
                f1()
                b = 2
                f2()
                return a + b
            return impl

        msg_head = "Calling either numba.gdb() or numba.gdb_init() more than"

        def check(func):
            with self.assertRaises(errors.UnsupportedError) as raises:
                func()
            self.assertIn(msg_head, str(raises.exception))

        for g1, g2 in permutations([gdb, gdb_init]):
            func = gen(g1, g2)
            check(func)

        @njit
        def use_globals():
            a = 1
            gdb()
            b = 2
            gdb_init()
            return a + b

        check(use_globals)


@not_unix
class TestGdbExceptions(TestCase):

    def test_call_gdb(self):
        def nop_compiler(x):
            return x
        for compiler in [nop_compiler, jit(forceobj=True), njit]:
            for meth in [gdb, gdb_init]:
                def python_func():
                    meth()
                with self.assertRaises(errors.TypingError) as raises:
                    compiler(python_func)()
                msg = "gdb support is only available on unix-like systems"
                self.assertIn(msg, str(raises.exception))


if __name__ == '__main__':
    unittest.main()
