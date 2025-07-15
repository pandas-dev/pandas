"""
Assorted utilities for use in tests.
"""

import cmath
import contextlib
from collections import defaultdict
import enum
import gc
import math
import platform
import os
import signal
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import io
import ctypes
import multiprocessing as mp
import warnings
import traceback
from contextlib import contextmanager
import uuid
import importlib
import types as pytypes
from functools import cached_property

import numpy as np

from numba import testing, types
from numba.core import errors, typing, utils, config, cpu
from numba.core.typing import cffi_utils
from numba.core.compiler import (compile_extra, Flags,
                                 DEFAULT_FLAGS, CompilerBase,
                                 DefaultPassBuilder)
from numba.core.typed_passes import IRLegalization
from numba.core.untyped_passes import PreserveIR
import unittest
from numba.core.runtime import rtsys
from numba.np import numpy_support
from numba.core.runtime import _nrt_python as _nrt
from numba.core.extending import (
    overload_method,
    typeof_impl,
    register_model,
    unbox,
    NativeValue,
    models,
)
from numba.core.datamodel.models import OpaqueModel

try:
    import scipy
except ImportError:
    scipy = None

# Make sure that coverage is set up.
try:
    import coverage
except ImportError:
    pass
else:
    coverage.process_startup()

enable_pyobj_flags = Flags()
enable_pyobj_flags.enable_pyobject = True

force_pyobj_flags = Flags()
force_pyobj_flags.force_pyobject = True

no_pyobj_flags = Flags()

nrt_flags = Flags()
nrt_flags.nrt = True


tag = testing.make_tag_decorator(['important', 'long_running', 'always_test'])

# Use to mark a test as a test that must always run when sharded
always_test = tag('always_test')

_32bit = sys.maxsize <= 2 ** 32
is_parfors_unsupported = _32bit
skip_parfors_unsupported = unittest.skipIf(
    is_parfors_unsupported,
    'parfors not supported',
)

skip_unless_py10_or_later = unittest.skipUnless(
    utils.PYVERSION >= (3, 10),
    "needs Python 3.10 or later"
)

skip_unless_py10 = unittest.skipUnless(
    utils.PYVERSION == (3, 10),
    "needs Python 3.10"
)

skip_unless_py312 = unittest.skipUnless(
    utils.PYVERSION == (3, 12),
    "needs Python 3.12"
)

skip_if_py313_on_windows = unittest.skipIf(
     utils.PYVERSION == (3, 13) and sys.platform.startswith('win'),
     "Not supported on Python 3.13 on Windows"
 )

skip_if_32bit = unittest.skipIf(_32bit, "Not supported on 32 bit")

IS_NUMPY_2 = numpy_support.numpy_version >= (2, 0)
skip_if_numpy_2 = unittest.skipIf(IS_NUMPY_2,
                                  "Not supported on numpy 2.0+")

def expected_failure_py311(fn):
    if utils.PYVERSION == (3, 11):
        return unittest.expectedFailure(fn)
    else:
        return fn


def expected_failure_py312(fn):
    if utils.PYVERSION == (3, 12):
        return unittest.expectedFailure(fn)
    else:
        return fn


def expected_failure_py313(fn):
    if utils.PYVERSION == (3, 13):
        return unittest.expectedFailure(fn)
    else:
        return fn


def expected_failure_np2(fn):
    if numpy_support.numpy_version == (2, 0):
        return unittest.expectedFailure(fn)
    else:
        return fn

_msg = "SciPy needed for test"
skip_unless_scipy = unittest.skipIf(scipy is None, _msg)

skip_unless_cffi = unittest.skipUnless(cffi_utils.SUPPORTED, 'requires cffi')

_lnx_reason = 'linux only test'
linux_only = unittest.skipIf(not sys.platform.startswith('linux'), _lnx_reason)

_win_reason = 'Windows-only test'
windows_only = unittest.skipIf(not sys.platform.startswith('win'), _win_reason)

_is_armv7l = platform.machine() == 'armv7l'

disabled_test = unittest.skipIf(True, 'Test disabled')

# See issue #4563, PPC64LE LLVM bug
skip_ppc64le_issue4563 = unittest.skipIf(platform.machine() == 'ppc64le',
                                         ("Hits: 'Parameter area must exist "
                                          "to pass an argument in memory'"))

# Typeguard
has_typeguard = bool(os.environ.get('NUMBA_USE_TYPEGUARD', 0))

skip_unless_typeguard = unittest.skipUnless(
    has_typeguard, "Typeguard is not enabled",
)

skip_if_typeguard = unittest.skipIf(
    has_typeguard, "Broken if Typeguard is enabled",
)

# See issue #6465, PPC64LE LLVM bug
skip_ppc64le_issue6465 = unittest.skipIf(platform.machine() == 'ppc64le',
                                         ("Hits: 'mismatch in size of "
                                          "parameter area' in "
                                          "LowerCall_64SVR4"))

# LLVM PPC issue.
# Sample error message:
#   Invalid PPC CTR loop!
#   UNREACHABLE executed at /llvm/lib/Target/PowerPC/PPCCTRLoops.cpp:179!
skip_ppc64le_invalid_ctr_loop = unittest.skipIf(
    platform.machine() == 'ppc64le',
    "Invalid PPC CTR loop")

# fenv.h on M1 may have various issues:
# https://github.com/numba/numba/issues/7822#issuecomment-1065356758
_uname = platform.uname()
IS_MACOS = _uname.system == 'Darwin'
skip_macos_fenv_errors = unittest.skipIf(IS_MACOS,
    "fenv.h-like functionality unreliable on macOS")
IS_MACOS_ARM64 = IS_MACOS and _uname.machine == 'arm64'

try:
    import scipy.linalg.cython_lapack
    has_lapack = True
except ImportError:
    has_lapack = False

needs_lapack = unittest.skipUnless(has_lapack,
                                   "LAPACK needs SciPy 1.0+")

try:
    import scipy.linalg.cython_blas
    has_blas = True
except ImportError:
    has_blas = False

needs_blas = unittest.skipUnless(has_blas, "BLAS needs SciPy 1.0+")

# Decorate a test with @needs_subprocess to ensure it doesn't run unless the
# `SUBPROC_TEST` environment variable is set. Use this in conjunction with:
# TestCase::subprocess_test_runner which will execute a given test in subprocess
# with this environment variable set.
_exec_cond = os.environ.get('SUBPROC_TEST', None) == '1'
needs_subprocess = unittest.skipUnless(_exec_cond, "needs subprocess harness")


try:
    import setuptools
    has_setuptools = True
except ImportError:
    has_setuptools = False


# decorator for a test that need setuptools
needs_setuptools = unittest.skipUnless(has_setuptools, 'Test needs setuptools')


def ignore_internal_warnings():
    """Use in testing within a ` warnings.catch_warnings` block to filter out
    warnings that are unrelated/internally generated by Numba.
    """
    # Filter out warnings from typeguard
    warnings.filterwarnings('ignore', module="typeguard")
    # Filter out warnings about TBB interface mismatch
    warnings.filterwarnings(action='ignore',
                            message=r".*TBB_INTERFACE_VERSION.*",
                            category=errors.NumbaWarning,
                            module=r'numba\.np\.ufunc\.parallel.*')


class TestCase(unittest.TestCase):

    longMessage = True

    # A random state yielding the same random numbers for any test case.
    # Use as `self.random.<method name>`
    @cached_property
    def random(self):
        return np.random.RandomState(42)

    def reset_module_warnings(self, module):
        """
        Reset the warnings registry of a module.  This can be necessary
        as the warnings module is buggy in that regard.
        See http://bugs.python.org/issue4180
        """
        if isinstance(module, str):
            module = sys.modules[module]
        try:
            del module.__warningregistry__
        except AttributeError:
            pass

    @contextlib.contextmanager
    def assertTypingError(self):
        """
        A context manager that asserts the enclosed code block fails
        compiling in nopython mode.
        """
        _accepted_errors = (errors.LoweringError, errors.TypingError,
                            TypeError, NotImplementedError)
        with self.assertRaises(_accepted_errors) as cm:
            yield cm

    @contextlib.contextmanager
    def assertRefCount(self, *objects):
        """
        A context manager that asserts the given objects have the
        same reference counts before and after executing the
        enclosed block.
        """
        old_refcounts = [sys.getrefcount(x) for x in objects]
        yield
        gc.collect()
        new_refcounts = [sys.getrefcount(x) for x in objects]
        for old, new, obj in zip(old_refcounts, new_refcounts, objects):
            if old != new:
                self.fail("Refcount changed from %d to %d for object: %r"
                          % (old, new, obj))

    def assertRefCountEqual(self, *objects):
        gc.collect()
        rc = [sys.getrefcount(x) for x in objects]
        rc_0 = rc[0]
        for i in range(len(objects))[1:]:
            rc_i = rc[i]
            if rc_0 != rc_i:
                self.fail(f"Refcount for objects does not match. "
                          f"#0({rc_0}) != #{i}({rc_i}) does not match.")

    @contextlib.contextmanager
    def assertNoNRTLeak(self):
        """
        A context manager that asserts no NRT leak was created during
        the execution of the enclosed block.
        """
        old = rtsys.get_allocation_stats()
        yield
        new = rtsys.get_allocation_stats()
        total_alloc = new.alloc - old.alloc
        total_free = new.free - old.free
        total_mi_alloc = new.mi_alloc - old.mi_alloc
        total_mi_free = new.mi_free - old.mi_free
        self.assertEqual(total_alloc, total_free,
                         "number of data allocs != number of data frees")
        self.assertEqual(total_mi_alloc, total_mi_free,
                         "number of meminfo allocs != number of meminfo frees")


    _bool_types = (bool, np.bool_)
    _exact_typesets = [_bool_types, (int,), (str,), (np.integer,),
                       (bytes, np.bytes_)]
    _approx_typesets = [(float,), (complex,), (np.inexact)]
    _sequence_typesets = [(tuple, list)]
    _float_types = (float, np.floating)
    _complex_types = (complex, np.complexfloating)

    def _detect_family(self, numeric_object):
        """
        This function returns a string description of the type family
        that the object in question belongs to.  Possible return values
        are: "exact", "complex", "approximate", "sequence", and "unknown"
        """
        if isinstance(numeric_object, np.ndarray):
            return "ndarray"

        if isinstance(numeric_object, enum.Enum):
            return "enum"

        for tp in self._sequence_typesets:
            if isinstance(numeric_object, tp):
                return "sequence"

        for tp in self._exact_typesets:
            if isinstance(numeric_object, tp):
                return "exact"

        for tp in self._complex_types:
            if isinstance(numeric_object, tp):
                return "complex"

        for tp in self._approx_typesets:
            if isinstance(numeric_object, tp):
                return "approximate"

        return "unknown"

    def _fix_dtype(self, dtype):
        """
        Fix the given *dtype* for comparison.
        """
        # Under 64-bit Windows, Numpy may return either int32 or int64
        # arrays depending on the function.
        if (sys.platform == 'win32' and sys.maxsize > 2**32 and
            dtype == np.dtype('int32')):
            return np.dtype('int64')
        else:
            return dtype

    def _fix_strides(self, arr):
        """
        Return the strides of the given array, fixed for comparison.
        Strides for 0- or 1-sized dimensions are ignored.
        """
        if arr.size == 0:
            return [0] * arr.ndim
        else:
            return [stride / arr.itemsize
                    for (stride, shape) in zip(arr.strides, arr.shape)
                    if shape > 1]

    def assertStridesEqual(self, first, second):
        """
        Test that two arrays have the same shape and strides.
        """
        self.assertEqual(first.shape, second.shape, "shapes differ")
        self.assertEqual(first.itemsize, second.itemsize, "itemsizes differ")
        self.assertEqual(self._fix_strides(first), self._fix_strides(second),
                         "strides differ")

    def assertPreciseEqual(self, first, second, prec='exact', ulps=1,
                           msg=None, ignore_sign_on_zero=False,
                           abs_tol=None
                           ):
        """
        Versatile equality testing function with more built-in checks than
        standard assertEqual().

        For arrays, test that layout, dtype, shape are identical, and
        recursively call assertPreciseEqual() on the contents.

        For other sequences, recursively call assertPreciseEqual() on
        the contents.

        For scalars, test that two scalars or have similar types and are
        equal up to a computed precision.
        If the scalars are instances of exact types or if *prec* is
        'exact', they are compared exactly.
        If the scalars are instances of inexact types (float, complex)
        and *prec* is not 'exact', then the number of significant bits
        is computed according to the value of *prec*: 53 bits if *prec*
        is 'double', 24 bits if *prec* is single.  This number of bits
        can be lowered by raising the *ulps* value.
        ignore_sign_on_zero can be set to True if zeros are to be considered
        equal regardless of their sign bit.
        abs_tol if this is set to a float value its value is used in the
        following. If, however, this is set to the string "eps" then machine
        precision of the type(first) is used in the following instead. This
        kwarg is used to check if the absolute difference in value between first
        and second is less than the value set, if so the numbers being compared
        are considered equal. (This is to handle small numbers typically of
        magnitude less than machine precision).

        Any value of *prec* other than 'exact', 'single' or 'double'
        will raise an error.
        """
        try:
            self._assertPreciseEqual(first, second, prec, ulps, msg,
                ignore_sign_on_zero, abs_tol)
        except AssertionError as exc:
            failure_msg = str(exc)
            # Fall off of the 'except' scope to avoid Python 3 exception
            # chaining.
        else:
            return
        # Decorate the failure message with more information
        self.fail("when comparing %s and %s: %s" % (first, second, failure_msg))

    def _assertPreciseEqual(self, first, second, prec='exact', ulps=1,
                            msg=None, ignore_sign_on_zero=False,
                            abs_tol=None):
        """Recursive workhorse for assertPreciseEqual()."""

        def _assertNumberEqual(first, second, delta=None):
            if (delta is None or first == second == 0.0
                or math.isinf(first) or math.isinf(second)):
                self.assertEqual(first, second, msg=msg)
                # For signed zeros
                if not ignore_sign_on_zero:
                    try:
                        if math.copysign(1, first) != math.copysign(1, second):
                            self.fail(
                                self._formatMessage(msg,
                                                    "%s != %s" %
                                                    (first, second)))
                    except TypeError:
                        pass
            else:
                self.assertAlmostEqual(first, second, delta=delta, msg=msg)

        first_family = self._detect_family(first)
        second_family = self._detect_family(second)

        assertion_message = "Type Family mismatch. (%s != %s)" % (first_family,
            second_family)
        if msg:
            assertion_message += ': %s' % (msg,)
        self.assertEqual(first_family, second_family, msg=assertion_message)

        # We now know they are in the same comparison family
        compare_family = first_family

        # For recognized sequences, recurse
        if compare_family == "ndarray":
            dtype = self._fix_dtype(first.dtype)
            self.assertEqual(dtype, self._fix_dtype(second.dtype))
            self.assertEqual(first.ndim, second.ndim,
                             "different number of dimensions")
            self.assertEqual(first.shape, second.shape,
                             "different shapes")
            self.assertEqual(first.flags.writeable, second.flags.writeable,
                             "different mutability")
            # itemsize is already checked by the dtype test above
            self.assertEqual(self._fix_strides(first),
                self._fix_strides(second), "different strides")
            if first.dtype != dtype:
                first = first.astype(dtype)
            if second.dtype != dtype:
                second = second.astype(dtype)
            for a, b in zip(first.flat, second.flat):
                self._assertPreciseEqual(a, b, prec, ulps, msg,
                                         ignore_sign_on_zero, abs_tol)
            return

        elif compare_family == "sequence":
            self.assertEqual(len(first), len(second), msg=msg)
            for a, b in zip(first, second):
                self._assertPreciseEqual(a, b, prec, ulps, msg,
                                         ignore_sign_on_zero, abs_tol)
            return

        elif compare_family == "exact":
            exact_comparison = True

        elif compare_family in ["complex", "approximate"]:
            exact_comparison = False

        elif compare_family == "enum":
            self.assertIs(first.__class__, second.__class__)
            self._assertPreciseEqual(first.value, second.value,
                                     prec, ulps, msg,
                                     ignore_sign_on_zero, abs_tol)
            return

        elif compare_family == "unknown":
            # Assume these are non-numeric types: we will fall back
            # on regular unittest comparison.
            self.assertIs(first.__class__, second.__class__)
            exact_comparison = True

        else:
            assert 0, "unexpected family"

        # If a Numpy scalar, check the dtype is exactly the same too
        # (required for datetime64 and timedelta64).
        if hasattr(first, 'dtype') and hasattr(second, 'dtype'):
            self.assertEqual(first.dtype, second.dtype)

        # Mixing bools and non-bools should always fail
        if (isinstance(first, self._bool_types) !=
            isinstance(second, self._bool_types)):
            assertion_message = ("Mismatching return types (%s vs. %s)"
                                 % (first.__class__, second.__class__))
            if msg:
                assertion_message += ': %s' % (msg,)
            self.fail(assertion_message)

        try:
            if cmath.isnan(first) and cmath.isnan(second):
                # The NaNs will compare unequal, skip regular comparison
                return
        except TypeError:
            # Not floats.
            pass

        # if absolute comparison is set, use it
        if abs_tol is not None:
            if abs_tol == "eps":
                rtol = np.finfo(type(first)).eps
            elif isinstance(abs_tol, float):
                rtol = abs_tol
            else:
                raise ValueError("abs_tol is not \"eps\" or a float, found %s"
                    % abs_tol)
            if abs(first - second) < rtol:
                return

        exact_comparison = exact_comparison or prec == 'exact'

        if not exact_comparison and prec != 'exact':
            if prec == 'single':
                bits = 24
            elif prec == 'double':
                bits = 53
            else:
                raise ValueError("unsupported precision %r" % (prec,))
            k = 2 ** (ulps - bits - 1)
            delta = k * (abs(first) + abs(second))
        else:
            delta = None
        if isinstance(first, self._complex_types):
            _assertNumberEqual(first.real, second.real, delta)
            _assertNumberEqual(first.imag, second.imag, delta)
        elif isinstance(first, (np.timedelta64, np.datetime64)):
            # Since Np 1.16 NaT == NaT is False, so special comparison needed
            if np.isnat(first):
                self.assertEqual(np.isnat(first), np.isnat(second))
            else:
                _assertNumberEqual(first, second, delta)
        else:
            _assertNumberEqual(first, second, delta)

    def subprocess_test_runner(self, test_module, test_class=None,
                               test_name=None, envvars=None, timeout=60,
                               flags=None, _subproc_test_env="1"):
        """
        Runs named unit test(s) as specified in the arguments as:
        test_module.test_class.test_name. test_module must always be supplied
        and if no further refinement is made with test_class and test_name then
        all tests in the module will be run. The tests will be run in a
        subprocess with environment variables specified in `envvars`.
        If given, envvars must be a map of form:
            environment variable name (str) -> value (str)
        If given, flags must be a map of form:
            flag including the `-` (str) -> value (str)
        It is most convenient to use this method in conjunction with
        @needs_subprocess as the decorator will cause the decorated test to be
        skipped unless the `SUBPROC_TEST` environment variable is set to
        the same value of ``_subproc_test_env``
        (this special environment variable is set by this method such that the
        specified test(s) will not be skipped in the subprocess).


        Following execution in the subprocess this method will check the test(s)
        executed without error. The timeout kwarg can be used to allow more time
        for longer running tests, it defaults to 60 seconds.
        """
        themod = self.__module__
        thecls = type(self).__name__
        parts = (test_module, test_class, test_name)
        fully_qualified_test = '.'.join(x for x in parts if x is not None)
        flags_args = []
        if flags is not None:
            for flag, value in flags.items():
                flags_args.append(f'{flag}')
                flags_args.append(f'{value}')
        cmd = [sys.executable, *flags_args, '-m', 'numba.runtests',
               fully_qualified_test]
        env_copy = os.environ.copy()
        env_copy['SUBPROC_TEST'] = _subproc_test_env
        try:
            env_copy['COVERAGE_PROCESS_START'] = os.environ['COVERAGE_RCFILE']
        except KeyError:
            pass   # ignored
        envvars = pytypes.MappingProxyType({} if envvars is None else envvars)
        env_copy.update(envvars)
        status = subprocess.run(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, timeout=timeout,
                                env=env_copy, universal_newlines=True)
        streams = (f'\ncaptured stdout: {status.stdout}\n'
                   f'captured stderr: {status.stderr}')
        self.assertEqual(status.returncode, 0, streams)
        # Python 3.12.1 report
        no_tests_ran = "NO TESTS RAN"
        if no_tests_ran in status.stderr:
            self.skipTest(no_tests_ran)
        else:
            self.assertIn('OK', status.stderr)
        return status

    def run_test_in_subprocess(maybefunc=None, timeout=60, envvars=None):
        """Runs the decorated test in a subprocess via invoking numba's test
        runner. kwargs timeout and envvars are passed through to
        subprocess_test_runner."""
        def wrapper(func):
            def inner(self, *args, **kwargs):
                if os.environ.get("SUBPROC_TEST", None) != func.__name__:
                    # Not in a subprocess test env, so stage the call to run the
                    # test in a subprocess which will set the env var.
                    class_name = self.__class__.__name__
                    self.subprocess_test_runner(
                        test_module=self.__module__,
                        test_class=class_name,
                        test_name=func.__name__,
                        timeout=timeout,
                        envvars=envvars,
                        _subproc_test_env=func.__name__,
                    )
                else:
                    # env var is set, so we're in the subprocess, run the
                    # actual test.
                    func(self)
            return inner

        if isinstance(maybefunc, pytypes.FunctionType):
            return wrapper(maybefunc)
        else:
            return wrapper

    def make_dummy_type(self):
        """Use to generate a dummy type unique to this test. Returns a python
        Dummy class and a corresponding Numba type DummyType."""

        # Use test_id to make sure no collision is possible.
        test_id = self.id()
        DummyType = type('DummyTypeFor{}'.format(test_id), (types.Opaque,), {})

        dummy_type = DummyType("my_dummy")
        register_model(DummyType)(OpaqueModel)

        class Dummy(object):
            pass

        @typeof_impl.register(Dummy)
        def typeof_dummy(val, c):
            return dummy_type

        @unbox(DummyType)
        def unbox_dummy(typ, obj, c):
            return NativeValue(c.context.get_dummy_value())

        return Dummy, DummyType

    def skip_if_no_external_compiler(self):
        """
        Call this to ensure the test is skipped if no suitable external compiler
        is found. This is a method on the TestCase opposed to a stand-alone
        decorator so as to make it "lazy" via runtime evaluation opposed to
        running at test-discovery time.
        """
        # This is a local import to avoid deprecation warnings being generated
        # through the use of the numba.pycc module.
        from numba.pycc.platform import external_compiler_works
        if not external_compiler_works():
            self.skipTest("No suitable external compiler was found.")


class SerialMixin(object):
    """Mixin to mark test for serial execution.
    """
    _numba_parallel_test_ = False


# Various helpers

@contextlib.contextmanager
def override_config(name, value):
    """
    Return a context manager that temporarily sets Numba config variable
    *name* to *value*.  *name* must be the name of an existing variable
    in numba.config.
    """
    old_value = getattr(config, name)
    setattr(config, name, value)
    try:
        yield
    finally:
        setattr(config, name, old_value)


@contextlib.contextmanager
def override_env_config(name, value):
    """
    Return a context manager that temporarily sets an Numba config environment
    *name* to *value*.
    """
    old = os.environ.get(name)
    os.environ[name] = value
    config.reload_config()

    try:
        yield
    finally:
        if old is None:
            # If it wasn't set originally, delete the environ var
            del os.environ[name]
        else:
            # Otherwise, restore to the old value
            os.environ[name] = old
        # Always reload config
        config.reload_config()


def compile_function(name, code, globs):
    """
    Given a *code* string, compile it with globals *globs* and return
    the function named *name*.
    """
    co = compile(code.rstrip(), "<string>", "single")
    ns = {}
    eval(co, globs, ns)
    return ns[name]


_trashcan_dir = 'numba-tests'

if os.name == 'nt':
    # Under Windows, gettempdir() points to the user-local temp dir
    _trashcan_dir = os.path.join(tempfile.gettempdir(), _trashcan_dir)
else:
    # Mix the UID into the directory name to allow different users to
    # run the test suite without permission errors (issue #1586)
    _trashcan_dir = os.path.join(tempfile.gettempdir(),
                                 "%s.%s" % (_trashcan_dir, os.getuid()))

# Stale temporary directories are deleted after they are older than this value.
# The test suite probably won't ever take longer than this...
_trashcan_timeout = 24 * 3600  # 1 day

def _create_trashcan_dir():
    try:
        os.mkdir(_trashcan_dir)
    except FileExistsError:
        pass

def _purge_trashcan_dir():
    freshness_threshold = time.time() - _trashcan_timeout
    for fn in sorted(os.listdir(_trashcan_dir)):
        fn = os.path.join(_trashcan_dir, fn)
        try:
            st = os.stat(fn)
            if st.st_mtime < freshness_threshold:
                shutil.rmtree(fn, ignore_errors=True)
        except OSError as e:
            # In parallel testing, several processes can attempt to
            # remove the same entry at once, ignore.
            pass

def _create_trashcan_subdir(prefix):
    _purge_trashcan_dir()
    path = tempfile.mkdtemp(prefix=prefix + '-', dir=_trashcan_dir)
    return path

def temp_directory(prefix):
    """
    Create a temporary directory with the given *prefix* that will survive
    at least as long as this process invocation.  The temporary directory
    will be eventually deleted when it becomes stale enough.

    This is necessary because a DLL file can't be deleted while in use
    under Windows.

    An interesting side-effect is to be able to inspect the test files
    shortly after a test suite run.
    """
    _create_trashcan_dir()
    return _create_trashcan_subdir(prefix)


def import_dynamic(modname):
    """
    Import and return a module of the given name.  Care is taken to
    avoid issues due to Python's internal directory caching.
    """
    import importlib
    importlib.invalidate_caches()
    __import__(modname)
    return sys.modules[modname]


# From CPython

@contextlib.contextmanager
def captured_output(stream_name):
    """Return a context manager used by captured_stdout/stdin/stderr
    that temporarily replaces the sys stream *stream_name* with a StringIO."""
    orig_stdout = getattr(sys, stream_name)
    setattr(sys, stream_name, io.StringIO())
    try:
        yield getattr(sys, stream_name)
    finally:
        setattr(sys, stream_name, orig_stdout)

def captured_stdout():
    """Capture the output of sys.stdout:

       with captured_stdout() as stdout:
           print("hello")
       self.assertEqual(stdout.getvalue(), "hello\n")
    """
    return captured_output("stdout")

def captured_stderr():
    """Capture the output of sys.stderr:

       with captured_stderr() as stderr:
           print("hello", file=sys.stderr)
       self.assertEqual(stderr.getvalue(), "hello\n")
    """
    return captured_output("stderr")


@contextlib.contextmanager
def capture_cache_log():
    with captured_stdout() as out:
        with override_config('DEBUG_CACHE', True):
            yield out


class EnableNRTStatsMixin(object):
    """Mixin to enable the NRT statistics counters."""

    def setUp(self):
        _nrt.memsys_enable_stats()

    def tearDown(self):
        _nrt.memsys_disable_stats()


class MemoryLeak(object):

    __enable_leak_check = True

    def memory_leak_setup(self):
        # Clean up any NRT-backed objects hanging in a dead reference cycle
        gc.collect()
        self.__init_stats = rtsys.get_allocation_stats()

    def memory_leak_teardown(self):
        if self.__enable_leak_check:
            self.assert_no_memory_leak()

    def assert_no_memory_leak(self):
        old = self.__init_stats
        new = rtsys.get_allocation_stats()
        total_alloc = new.alloc - old.alloc
        total_free = new.free - old.free
        total_mi_alloc = new.mi_alloc - old.mi_alloc
        total_mi_free = new.mi_free - old.mi_free
        self.assertEqual(total_alloc, total_free)
        self.assertEqual(total_mi_alloc, total_mi_free)

    def disable_leak_check(self):
        # For per-test use when MemoryLeakMixin is injected into a TestCase
        self.__enable_leak_check = False


class MemoryLeakMixin(EnableNRTStatsMixin, MemoryLeak):

    def setUp(self):
        super(MemoryLeakMixin, self).setUp()
        self.memory_leak_setup()

    def tearDown(self):
        gc.collect()
        self.memory_leak_teardown()
        super(MemoryLeakMixin, self).tearDown()


@contextlib.contextmanager
def forbid_codegen():
    """
    Forbid LLVM code generation during the execution of the context
    manager's enclosed block.

    If code generation is invoked, a RuntimeError is raised.
    """
    from numba.core import codegen
    patchpoints = ['CPUCodeLibrary._finalize_final_module']

    old = {}
    def fail(*args, **kwargs):
        raise RuntimeError("codegen forbidden by test case")
    try:
        # XXX use the mock library instead?
        for name in patchpoints:
            parts = name.split('.')
            obj = codegen
            for attrname in parts[:-1]:
                obj = getattr(obj, attrname)
            attrname = parts[-1]
            value = getattr(obj, attrname)
            assert callable(value), ("%r should be callable" % name)
            old[obj, attrname] = value
            setattr(obj, attrname, fail)
        yield
    finally:
        for (obj, attrname), value in old.items():
            setattr(obj, attrname, value)


# For details about redirection of file-descriptor, read
# https://eli.thegreenplace.net/2015/redirecting-all-kinds-of-stdout-in-python/

@contextlib.contextmanager
def redirect_fd(fd):
    """
    Temporarily redirect *fd* to a pipe's write end and return a file object
    wrapping the pipe's read end.
    """

    from numba import _helperlib
    libnumba = ctypes.CDLL(_helperlib.__file__)

    libnumba._numba_flush_stdout()
    save = os.dup(fd)
    r, w = os.pipe()
    try:
        os.dup2(w, fd)
        yield io.open(r, "r")
    finally:
        libnumba._numba_flush_stdout()
        os.close(w)
        os.dup2(save, fd)
        os.close(save)


def redirect_c_stdout():
    """Redirect C stdout
    """
    fd = sys.__stdout__.fileno()
    return redirect_fd(fd)


def redirect_c_stderr():
    """Redirect C stderr
    """
    fd = sys.__stderr__.fileno()
    return redirect_fd(fd)


def run_in_new_process_caching(func, cache_dir_prefix=__name__, verbose=True):
    """Spawn a new process to run `func` with a temporary cache directory.

    The childprocess's stdout and stderr will be captured and redirected to
    the current process's stdout and stderr.

    Returns
    -------
    ret : dict
        exitcode: 0 for success. 1 for exception-raised.
        stdout: str
        stderr: str
    """
    cache_dir = temp_directory(cache_dir_prefix)
    return run_in_new_process_in_cache_dir(func, cache_dir, verbose=verbose)


def run_in_new_process_in_cache_dir(func, cache_dir, verbose=True):
    """Spawn a new process to run `func` with a temporary cache directory.

    The childprocess's stdout and stderr will be captured and redirected to
    the current process's stdout and stderr.

    Similar to ``run_in_new_process_caching()`` but the ``cache_dir`` is a
    directory path instead of a name prefix for the directory path.

    Returns
    -------
    ret : dict
        exitcode: 0 for success. 1 for exception-raised.
        stdout: str
        stderr: str
    """
    ctx = mp.get_context('spawn')
    qout = ctx.Queue()
    with override_env_config('NUMBA_CACHE_DIR', cache_dir):
        proc = ctx.Process(target=_remote_runner, args=[func, qout])
        proc.start()
        proc.join()
        stdout = qout.get_nowait()
        stderr = qout.get_nowait()
        if verbose and stdout.strip():
            print()
            print('STDOUT'.center(80, '-'))
            print(stdout)
        if verbose and stderr.strip():
            print(file=sys.stderr)
            print('STDERR'.center(80, '-'), file=sys.stderr)
            print(stderr, file=sys.stderr)
    return {
        'exitcode': proc.exitcode,
        'stdout': stdout,
        'stderr': stderr,
    }


def _remote_runner(fn, qout):
    """Used by `run_in_new_process_caching()`
    """
    with captured_stderr() as stderr:
        with captured_stdout() as stdout:
            try:
                fn()
            except Exception:
                traceback.print_exc()
                exitcode = 1
            else:
                exitcode = 0
        qout.put(stdout.getvalue())
    qout.put(stderr.getvalue())
    sys.exit(exitcode)

class CheckWarningsMixin(object):
    @contextlib.contextmanager
    def check_warnings(self, messages, category=RuntimeWarning):
        with warnings.catch_warnings(record=True) as catch:
            warnings.simplefilter("always")
            yield
        found = 0
        for w in catch:
            for m in messages:
                if m in str(w.message):
                    self.assertEqual(w.category, category)
                    found += 1
        self.assertEqual(found, len(messages))


def _format_jit_options(**jit_options):
    if not jit_options:
        return ''
    out = []
    for key, value in jit_options.items():
        if isinstance(value, str):
            value = '"{}"'.format(value)
        out.append('{}={}'.format(key, value))
    return ', '.join(out)


@contextlib.contextmanager
def create_temp_module(source_lines, **jit_options):
    """A context manager that creates and imports a temporary module
    from sources provided in ``source_lines``.

    Optionally it is possible to provide jit options for ``jit_module`` if it
    is explicitly used in ``source_lines`` like ``jit_module({jit_options})``.
    """
    # Use try/finally so cleanup happens even when an exception is raised
    try:
        tempdir = temp_directory('test_temp_module')
        # Generate random module name
        temp_module_name = 'test_temp_module_{}'.format(
            str(uuid.uuid4()).replace('-', '_'))
        temp_module_path = os.path.join(tempdir, temp_module_name + '.py')

        jit_options = _format_jit_options(**jit_options)
        with open(temp_module_path, 'w') as f:
            lines = source_lines.format(jit_options=jit_options)
            f.write(lines)
        # Add test_module to sys.path so it can be imported
        sys.path.insert(0, tempdir)
        test_module = importlib.import_module(temp_module_name)
        yield test_module
    finally:
        sys.modules.pop(temp_module_name, None)
        sys.path.remove(tempdir)
        shutil.rmtree(tempdir)


def run_in_subprocess(code, flags=None, env=None, timeout=30):
    """Run a snippet of Python code in a subprocess with flags, if any are
    given. 'env' is passed to subprocess.Popen(). 'timeout' is passed to
    popen.communicate().

    Returns the stdout and stderr of the subprocess after its termination.
    """
    if flags is None:
        flags = []
    cmd = [sys.executable,] + flags + ["-c", code]
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, env=env)
    out, err = popen.communicate(timeout=timeout)
    if popen.returncode != 0:
        msg = "process failed with code %s: stderr follows\n%s\n"
        raise AssertionError(msg % (popen.returncode, err.decode()))
    return out, err


def strace(work, syscalls, timeout=10):
    """Runs strace whilst executing the function work() in the current process,
    captures the listed syscalls (list of strings). Takes an optional timeout in
    seconds, default is 10, if this is exceeded the process will be sent a
    SIGKILL. Returns a list of lines that are output by strace.
    """

    # Open a tmpfile for strace to write into.
    with tempfile.NamedTemporaryFile('w+t') as ntf:

        parent_pid = os.getpid()
        strace_binary = shutil.which('strace')
        if strace_binary is None:
            raise ValueError("No valid 'strace' binary could be found")
        cmd = [strace_binary, # strace
               '-q', # quietly (no attach/detach print out)
               '-p', str(parent_pid), # this PID
               '-e', ','.join(syscalls), # these syscalls
               '-o', ntf.name] # put output into this file

        # redirect stdout, stderr is handled by the `-o` flag to strace.
        popen = subprocess.Popen(cmd, stdout=subprocess.PIPE,)
        strace_pid = popen.pid
        thread_timeout = threading.Timer(timeout, popen.kill)
        thread_timeout.start()

        def check_return(problem=''):
            ret = popen.returncode
            if ret != 0:
                msg = ("strace exited non-zero, process return code was:"
                       f"{ret}. {problem}")
                raise RuntimeError(msg)
        try:
            # push the communication onto a thread so it doesn't block.
            # start comms thread
            thread_comms = threading.Thread(target=popen.communicate)
            thread_comms.start()

            # do work
            work()
            # Flush the output buffer file
            ntf.flush()
            # interrupt the strace process to stop it if it's still running
            if popen.poll() is None:
                os.kill(strace_pid, signal.SIGINT)
            else:
                # it's not running, probably an issue, raise
                problem="If this is SIGKILL, increase the timeout?"
                check_return(problem)
            # Make sure the return code is 0, SIGINT to detach is considered
            # a successful exit.
            popen.wait()
            check_return()
            # collect the data
            strace_data = ntf.readlines()
        finally:
            # join communication, should be stopped now as process has
            # exited
            thread_comms.join()
            # should be stopped already
            thread_timeout.cancel()

    return strace_data


def strace_supported():
    """Checks if strace is supported and working"""

    # Only support this on linux where the `strace` binary is likely to be the
    # strace needed.
    if not sys.platform.startswith('linux'):
        return False

    def force_clone(): # subprocess triggers a clone
        subprocess.run([sys.executable, '-c', 'exit()'])

    syscall = 'clone'
    try:
        trace = strace(force_clone, [syscall,])
    except Exception:
        return False
    return syscall in ''.join(trace)


class IRPreservingTestPipeline(CompilerBase):
    """ Same as the standard pipeline, but preserves the func_ir into the
    metadata store after legalisation, useful for testing IR changes"""

    def define_pipelines(self):
        pipeline = DefaultPassBuilder.define_nopython_pipeline(
            self.state, "ir_preserving_custom_pipe")
        # mangle the default pipeline and inject DCE and IR preservation ahead
        # of legalisation

        # TODO: add a way to not do this! un-finalizing is not a good idea
        pipeline._finalized = False
        pipeline.add_pass_after(PreserveIR, IRLegalization)

        pipeline.finalize()
        return [pipeline]


def print_azure_matrix():
    """This is a utility function that prints out the map of NumPy to Python
    versions and how many of that combination are being tested across all the
    declared config for azure-pipelines. It is useful to run when updating the
    azure-pipelines config to be able to quickly see what the coverage is."""
    import yaml
    from yaml import Loader
    base_path = os.path.dirname(os.path.abspath(__file__))
    azure_pipe = os.path.join(base_path, '..', '..', 'azure-pipelines.yml')
    if not os.path.isfile(azure_pipe):
        raise RuntimeError("'azure-pipelines.yml' is not available")
    with open(os.path.abspath(azure_pipe), 'rt') as f:
        data = f.read()
    pipe_yml = yaml.load(data, Loader=Loader)

    templates = pipe_yml['jobs']
    # first look at the items in the first two templates, this is osx/linux
    py2np_map = defaultdict(lambda: defaultdict(int))
    for tmplt in templates[:2]:
        matrix = tmplt['parameters']['matrix']
        for setup in matrix.values():
            py2np_map[setup['NUMPY']][setup['PYTHON']]+=1

    # next look at the items in the windows only template
    winpath = ['..', '..', 'buildscripts', 'azure', 'azure-windows.yml']
    azure_windows = os.path.join(base_path, *winpath)
    if not os.path.isfile(azure_windows):
        raise RuntimeError("'azure-windows.yml' is not available")
    with open(os.path.abspath(azure_windows), 'rt') as f:
        data = f.read()
    windows_yml = yaml.load(data, Loader=Loader)

    # There's only one template in windows and its keyed differently to the
    # above, get its matrix.
    matrix = windows_yml['jobs'][0]['strategy']['matrix']
    for setup in matrix.values():
        py2np_map[setup['NUMPY']][setup['PYTHON']]+=1

    print("NumPy | Python | Count")
    print("-----------------------")
    for npver, pys in sorted(py2np_map.items()):
        for pyver, count in pys.items():
            print(f" {npver} |  {pyver:<4}  |   {count}")

    # print the "reverse" map
    rev_map = defaultdict(lambda: defaultdict(int))
    for npver, pys in sorted(py2np_map.items()):
        for pyver, count in pys.items():
            rev_map[pyver][npver] = count
    print("\nPython | NumPy | Count")
    print("-----------------------")
    sorter = lambda x: int(x[0].split('.')[1])
    for pyver, nps in sorted(rev_map.items(), key=sorter):
        for npver, count in nps.items():
            print(f" {pyver:<4} |  {npver}  |   {count}")
