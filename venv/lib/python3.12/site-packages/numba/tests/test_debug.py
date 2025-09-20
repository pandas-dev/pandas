import os
import platform
import re
import textwrap
import warnings

import numpy as np

from numba.tests.support import (TestCase, override_config, override_env_config,
                      captured_stdout, forbid_codegen, skip_parfors_unsupported,
                      needs_blas)
from numba import jit, njit
from numba.core import types, compiler, utils
from numba.core.errors import NumbaPerformanceWarning
from numba import prange
from numba.experimental import jitclass
import unittest


def simple_nopython(somearg):
    retval = somearg + 1
    return retval

def simple_gen(x, y):
    yield x
    yield y


class SimpleClass(object):
    def __init__(self):
        self.h = 5

simple_class_spec = [('h', types.int32)]

def simple_class_user(obj):
    return obj.h

def unsupported_parfor(a, b):
    return np.dot(a, b) # dot as gemm unsupported

def supported_parfor(n):
    a = np.ones(n)
    for i in prange(n):
        a[i] = a[i] + np.sin(i)
    return a

def unsupported_prange(n):
    a = np.ones(n)
    for i in prange(n):
        a[i] = a[i] + np.sin(i)
        assert i + 13 < 100000
    return a


class DebugTestBase(TestCase):

    all_dumps = set(['bytecode', 'cfg', 'ir', 'typeinfer', 'llvm',
                     'func_opt_llvm', 'optimized_llvm', 'assembly'])

    def assert_fails(self, *args, **kwargs):
        self.assertRaises(AssertionError, *args, **kwargs)

    def check_debug_output(self, out, dump_names):
        enabled_dumps = dict.fromkeys(self.all_dumps, False)
        for name in dump_names:
            assert name in enabled_dumps
            enabled_dumps[name] = True
        for name, enabled in sorted(enabled_dumps.items()):
            check_meth = getattr(self, '_check_dump_%s' % name)
            if enabled:
                check_meth(out)
            else:
                self.assert_fails(check_meth, out)

    def _check_dump_bytecode(self, out):
        if utils.PYVERSION in ((3, 11), (3, 12), (3, 13)):
            self.assertIn('BINARY_OP', out)
        elif utils.PYVERSION in ((3, 10),):
            self.assertIn('BINARY_ADD', out)
        else:
            raise NotImplementedError(utils.PYVERSION)

    def _check_dump_cfg(self, out):
        self.assertIn('CFG dominators', out)

    def _check_dump_ir(self, out):
        self.assertIn('--IR DUMP: %s--' % self.func_name, out)

    def _check_dump_typeinfer(self, out):
        self.assertIn('--propagate--', out)

    def _check_dump_llvm(self, out):
        self.assertIn('--LLVM DUMP', out)
        if compiler.Flags.options["auto_parallel"].default.enabled == False:
            self.assertRegex(out, r'store i64 %\"\.\d", i64\* %"retptr"', out)

    def _check_dump_func_opt_llvm(self, out):
        self.assertIn('--FUNCTION OPTIMIZED DUMP %s' % self.func_name, out)
        # allocas have been optimized away
        self.assertIn('add nsw i64 %arg.somearg, 1', out)

    def _check_dump_optimized_llvm(self, out):
        self.assertIn('--OPTIMIZED DUMP %s' % self.func_name, out)
        self.assertIn('add nsw i64 %arg.somearg, 1', out)

    def _check_dump_assembly(self, out):
        self.assertIn('--ASSEMBLY %s' % self.func_name, out)
        if platform.machine() in ('x86_64', 'AMD64', 'i386', 'i686'):
            self.assertIn('xorl', out)


class FunctionDebugTestBase(DebugTestBase):

    func_name = 'simple_nopython'

    def compile_simple_nopython(self):
        with captured_stdout() as out:
            cfunc = njit((types.int64,))(simple_nopython)
            # Sanity check compiled function
            self.assertPreciseEqual(cfunc(2), 3)
        return out.getvalue()


class TestFunctionDebugOutput(FunctionDebugTestBase):

    def test_dump_bytecode(self):
        with override_config('DUMP_BYTECODE', True):
            out = self.compile_simple_nopython()
        self.check_debug_output(out, ['bytecode'])

    def test_dump_ir(self):
        with override_config('DUMP_IR', True):
            out = self.compile_simple_nopython()
        self.check_debug_output(out, ['ir'])

    def test_dump_cfg(self):
        with override_config('DUMP_CFG', True):
            out = self.compile_simple_nopython()
        self.check_debug_output(out, ['cfg'])

    def test_dump_llvm(self):
        with override_config('DUMP_LLVM', True):
            out = self.compile_simple_nopython()
        self.check_debug_output(out, ['llvm'])

    def test_dump_func_opt_llvm(self):
        with override_config('DUMP_FUNC_OPT', True):
            out = self.compile_simple_nopython()
        self.check_debug_output(out, ['func_opt_llvm'])

    def test_dump_optimized_llvm(self):
        with override_config('DUMP_OPTIMIZED', True):
            out = self.compile_simple_nopython()
        self.check_debug_output(out, ['optimized_llvm'])

    def test_dump_assembly(self):
        with override_config('DUMP_ASSEMBLY', True):
            out = self.compile_simple_nopython()
        self.check_debug_output(out, ['assembly'])


class TestGeneratorDebugOutput(DebugTestBase):

    func_name = 'simple_gen'

    def compile_simple_gen(self):
        with captured_stdout() as out:
            cfunc = njit((types.int64, types.int64))(simple_gen)
            # Sanity check compiled function
            self.assertPreciseEqual(list(cfunc(2, 5)), [2, 5])
        return out.getvalue()

    def test_dump_ir_generator(self):
        with override_config('DUMP_IR', True):
            out = self.compile_simple_gen()
        self.check_debug_output(out, ['ir'])
        self.assertIn('--GENERATOR INFO: %s' % self.func_name, out)
        expected_gen_info = textwrap.dedent("""
            generator state variables: ['x', 'y']
            yield point #1: live variables = ['y'], weak live variables = ['x']
            yield point #2: live variables = [], weak live variables = ['y']
            """)
        self.assertIn(expected_gen_info, out)


class TestDisableJIT(DebugTestBase):
    """
    Test the NUMBA_DISABLE_JIT environment variable.
    """

    def test_jit(self):
        with override_config('DISABLE_JIT', True):
            with forbid_codegen():
                cfunc = jit(nopython=True)(simple_nopython)
                self.assertPreciseEqual(cfunc(2), 3)

    def test_jitclass(self):
        with override_config('DISABLE_JIT', True):
            with forbid_codegen():
                SimpleJITClass = jitclass(simple_class_spec)(SimpleClass)

                obj = SimpleJITClass()
                self.assertPreciseEqual(obj.h, 5)

                cfunc = jit(nopython=True)(simple_class_user)
                self.assertPreciseEqual(cfunc(obj), 5)


class TestEnvironmentOverride(FunctionDebugTestBase):
    """
    Test that environment variables are reloaded by Numba when modified.
    """

    # mutates env with os.environ so must be run serially
    _numba_parallel_test_ = False

    def test_debug(self):
        out = self.compile_simple_nopython()
        self.assertFalse(out)
        with override_env_config('NUMBA_DEBUG', '1'):
            out = self.compile_simple_nopython()
            # Note that all variables dependent on NUMBA_DEBUG are
            # updated too.
            self.check_debug_output(out, ['ir', 'typeinfer',
                                          'llvm', 'func_opt_llvm',
                                          'optimized_llvm', 'assembly'])
        out = self.compile_simple_nopython()
        self.assertFalse(out)

class TestParforsDebug(TestCase):
    """
    Tests debug options associated with parfors
    """

    # mutates env with os.environ so must be run serially
    _numba_parallel_test_ = False

    def check_parfors_warning(self, warn_list):
        msg = ("'parallel=True' was specified but no transformation for "
               "parallel execution was possible.")
        warning_found = False
        for w in warn_list:
            if msg in str(w.message):
                warning_found = True
                break
        self.assertTrue(warning_found, "Warning message should be found.")

    def check_parfors_unsupported_prange_warning(self, warn_list):
        msg = ("prange or pndindex loop will not be executed in parallel "
               "due to there being more than one entry to or exit from the "
               "loop (e.g., an assertion).")
        warning_found = False
        for w in warn_list:
            if msg in str(w.message):
                warning_found = True
                break
        self.assertTrue(warning_found, "Warning message should be found.")

    @needs_blas
    @skip_parfors_unsupported
    def test_warns(self):
        """
        Test that using parallel=True on a function that does not have parallel
        semantics warns.
        """
        arr_ty = types.Array(types.float64, 2, "C")
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", NumbaPerformanceWarning)
            njit((arr_ty, arr_ty), parallel=True)(unsupported_parfor)
        self.check_parfors_warning(w)

    @needs_blas
    @skip_parfors_unsupported
    def test_unsupported_prange_warns(self):
        """
        Test that prange with multiple exits issues a warning
        """
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", NumbaPerformanceWarning)
            njit((types.int64,), parallel=True)(unsupported_prange)
        self.check_parfors_unsupported_prange_warning(w)

    @skip_parfors_unsupported
    def test_array_debug_opt_stats(self):
        """
        Test that NUMBA_DEBUG_ARRAY_OPT_STATS produces valid output
        """
        # deliberately trigger a compilation loop to increment the
        # Parfor class state, this is to ensure the test works based
        # on indices computed based on this state and not hard coded
        # indices.
        njit((types.int64,), parallel=True)(supported_parfor)

        with override_env_config('NUMBA_DEBUG_ARRAY_OPT_STATS', '1'):
            with captured_stdout() as out:
                njit((types.int64,), parallel=True)(supported_parfor)

            # grab the various parts out the output
            output = out.getvalue().split('\n')
            parallel_loop_output = \
                [x for x in output if 'is produced from pattern' in x]
            fuse_output = \
                [x for x in output if 'is fused into' in x]
            after_fusion_output = \
                [x for x in output if 'After fusion, function' in x]

            # Parfor's have a shared state index, grab the current value
            # as it will be used as an offset for all loop messages
            parfor_state = int(re.compile(r'#([0-9]+)').search(
                parallel_loop_output[0]).group(1))
            bounds = range(parfor_state,
                           parfor_state + len(parallel_loop_output))

            # Check the Parallel for-loop <index> is produced from <pattern>
            # works first
            pattern = ("('ones function', 'NumPy mapping')",
                       ('prange', 'user', ''))
            fmt = 'Parallel for-loop #{} is produced from pattern \'{}\' at'
            for i, trials, lpattern in zip(bounds, parallel_loop_output,
                                           pattern):
                to_match = fmt.format(i, lpattern)
                self.assertIn(to_match, trials)

            # Check the fusion statements are correct
            pattern = (parfor_state + 1, parfor_state + 0)
            fmt = 'Parallel for-loop #{} is fused into for-loop #{}.'
            for trials in fuse_output:
                to_match = fmt.format(*pattern)
                self.assertIn(to_match, trials)

            # Check the post fusion statements are correct
            pattern = (supported_parfor.__name__, 1, set([parfor_state]))
            fmt = 'After fusion, function {} has {} parallel for-loop(s) #{}.'
            for trials in after_fusion_output:
                to_match = fmt.format(*pattern)
                self.assertIn(to_match, trials)


if __name__ == '__main__':
    unittest.main()
