"""
Unspecified error handling tests
"""

import sys
import subprocess
import numpy as np
import os
import warnings

from numba import jit, njit, types
from numba.core import errors
from numba.experimental import structref
from numba.extending import (overload, intrinsic, overload_method,
                             overload_attribute)
from numba.core.compiler import CompilerBase
from numba.core.untyped_passes import (TranslateByteCode, FixupArgs,
                                       IRProcessing,)
from numba.core.typed_passes import (NopythonTypeInference, DeadCodeElimination,
                                     NoPythonBackend, NativeLowering)
from numba.core.compiler_machinery import PassManager
from numba.core.types.functions import _err_reasons as error_reasons

from numba.tests.support import (skip_parfors_unsupported, override_config,
                                 SerialMixin, skip_unless_cffi,
                                 skip_unless_scipy, TestCase)
import unittest


class TestErrorHandlingBeforeLowering(unittest.TestCase):

    def test_unsupported_make_function_return_inner_func(self):
        def func(x):
            """ return the closure """
            z = x + 1

            def inner(x):
                return x + z
            return inner

        for pipeline in jit, njit:
            with self.assertRaises(errors.TypingError) as raises:
                pipeline(func)(1)

            expected = "Cannot capture the non-constant value"
            self.assertIn(expected, str(raises.exception))


class TestUnsupportedReporting(unittest.TestCase):

    def test_unsupported_numpy_function(self):
        # np.asanyarray(list) currently unsupported
        @njit
        def func():
            np.asanyarray([1,2,3])

        with self.assertRaises(errors.TypingError) as raises:
            func()

        expected = "Use of unsupported NumPy function 'numpy.asanyarray'"
        self.assertIn(expected, str(raises.exception))


class TestMiscErrorHandling(unittest.TestCase):

    def test_use_of_exception_for_flow_control(self):
        # constant inference uses exceptions with no Loc specified to determine
        # flow control, this asserts that the construction of the lowering
        # error context handler works in the case of an exception with no Loc
        # specified. See issue #3135.
        @njit
        def fn(x):
            return 10**x

        a = np.array([1.0],dtype=np.float64)
        fn(a) # should not raise

    def test_commented_func_definition_is_not_a_definition(self):
        # See issue #4056, the commented def should not be found as the
        # definition for reporting purposes when creating the synthetic
        # traceback because it is commented! Use of def in docstring would also
        # cause this issue hence is tested.

        def foo_commented():
            #def commented_definition()
            raise Exception('test_string')

        def foo_docstring():
            """ def docstring containing def might match function definition!"""
            raise Exception('test_string')

        for func in (foo_commented, foo_docstring):
            with self.assertRaises(Exception) as raises:
                func()

            self.assertIn("test_string", str(raises.exception))

    def test_use_of_ir_unknown_loc(self):
        # for context see # 3390
        class TestPipeline(CompilerBase):
            def define_pipelines(self):
                name = 'bad_DCE_pipeline'
                pm = PassManager(name)
                pm.add_pass(TranslateByteCode, "analyzing bytecode")
                pm.add_pass(FixupArgs, "fix up args")
                pm.add_pass(IRProcessing, "processing IR")
                # remove dead before type inference so that the Arg node is
                # removed and the location of the arg cannot be found
                pm.add_pass(DeadCodeElimination, "DCE")
                # typing
                pm.add_pass(NopythonTypeInference, "nopython frontend")
                pm.add_pass(NativeLowering, "native lowering")
                pm.add_pass(NoPythonBackend, "nopython mode backend")
                pm.finalize()
                return [pm]

        @njit(pipeline_class=TestPipeline)
        def f(a):
            return 0

        with self.assertRaises(errors.TypingError) as raises:
            f(iter([1,2]))  # use a type that Numba doesn't recognize

        expected = 'File "unknown location", line 0:'
        self.assertIn(expected, str(raises.exception))

    def check_write_to_globals(self, func):
        with self.assertRaises(errors.TypingError) as raises:
            func()

        expected = ["The use of a", "in globals, is not supported as globals"]
        for ex in expected:
            self.assertIn(ex, str(raises.exception))

    def test_handling_of_write_to_reflected_global(self):
        from numba.tests.errorhandling_usecases import global_reflected_write
        self.check_write_to_globals(njit(global_reflected_write))

    def test_handling_of_write_to_typed_dict_global(self):
        from numba.tests.errorhandling_usecases import global_dict_write
        self.check_write_to_globals(njit(global_dict_write))

    @skip_parfors_unsupported
    def test_handling_forgotten_numba_internal_import(self):
        @njit(parallel=True)
        def foo():
            for i in prange(10): # noqa: F821 prange is not imported
                pass

        with self.assertRaises(errors.TypingError) as raises:
            foo()

        expected = ("'prange' looks like a Numba internal function, "
                    "has it been imported")
        self.assertIn(expected, str(raises.exception))

    def test_handling_unsupported_generator_expression(self):
        def foo():
            (x for x in range(10))

        expected = "The use of yield in a closure is unsupported."

        for dec in jit(forceobj=True), njit:
            with self.assertRaises(errors.UnsupportedError) as raises:
                dec(foo)()
            self.assertIn(expected, str(raises.exception))

    def test_handling_undefined_variable(self):
        @njit
        def foo():
            return a # noqa: F821

        expected = "NameError: name 'a' is not defined"

        with self.assertRaises(errors.TypingError) as raises:
            foo()
        self.assertIn(expected, str(raises.exception))


class TestErrorMessages(unittest.TestCase):

    def test_specific_error(self):

        given_reason = "specific_reason"

        def foo():
            pass

        @overload(foo)
        def ol_foo():
            raise errors.NumbaValueError(given_reason)

        @njit
        def call_foo():
            foo()

        with self.assertRaises(errors.TypingError) as raises:
            call_foo()

        excstr = str(raises.exception)
        self.assertIn(error_reasons['specific_error'].splitlines()[0], excstr)
        self.assertIn(given_reason, excstr)

    def test_no_match_error(self):

        def foo():
            pass

        @overload(foo)
        def ol_foo():
            return None # emulate no impl available for type

        @njit
        def call_foo():
            foo()

        with self.assertRaises(errors.TypingError) as raises:
            call_foo()

        excstr = str(raises.exception)
        self.assertIn("No match", excstr)

    @skip_unless_scipy
    def test_error_function_source_is_correct(self):
        """ Checks that the reported source location for an overload is the
        overload implementation source, not the actual function source from the
        target library."""

        @njit
        def foo():
            np.linalg.svd("chars")

        with self.assertRaises(errors.TypingError) as raises:
            foo()

        excstr = str(raises.exception)
        self.assertIn(error_reasons['specific_error'].splitlines()[0], excstr)
        expected_file = os.path.join("numba", "np", "linalg.py")
        expected = f"Overload in function 'svd_impl': File: {expected_file}:"
        self.assertIn(expected.format(expected_file), excstr)

    def test_concrete_template_source(self):
        # hits ConcreteTemplate
        @njit
        def foo():
            return 'a' + 1

        with self.assertRaises(errors.TypingError) as raises:
            foo()

        excstr = str(raises.exception)

        self.assertIn("Overload of function 'add'", excstr)
        # there'll be numerous matched templates that don't work but as they
        # are mostly "overload"s they'll just appear as "No match".
        self.assertIn("No match.", excstr)

    def test_abstract_template_source(self):
        # hits AbstractTemplate
        @njit
        def foo():
            return len(1)

        with self.assertRaises(errors.TypingError) as raises:
            foo()

        excstr = str(raises.exception)
        self.assertIn("Overload of function 'len'", excstr)

    def test_callable_template_source(self):
        # hits CallableTemplate
        @njit
        def foo():
            return np.angle(None)

        with self.assertRaises(errors.TypingError) as raises:
            foo()

        excstr = str(raises.exception)
        self.assertIn("No implementation of function Function(<function angle",
                      excstr)

    def test_overloadfunction_template_source(self):
        # hits _OverloadFunctionTemplate
        def bar(x):
            pass

        @overload(bar)
        def ol_bar(x):
            pass

        @njit
        def foo():
            return bar(1)

        with self.assertRaises(errors.TypingError) as raises:
            foo()

        excstr = str(raises.exception)
        # there will not be "numerous" matched templates, there's just one,
        # the one above, so assert it is reported
        self.assertNotIn("<numerous>", excstr)
        expected_file = os.path.join("numba", "tests",
                                     "test_errorhandling.py")
        expected_ol = f"Overload of function 'bar': File: {expected_file}:"
        self.assertIn(expected_ol.format(expected_file), excstr)
        self.assertIn("No match.", excstr)

    def test_intrinsic_template_source(self):
        # hits _IntrinsicTemplate
        given_reason1 = "x must be literal"
        given_reason2 = "array.ndim must be 1"

        @intrinsic
        def myintrin(typingctx, x, arr):
            if not isinstance(x, types.IntegerLiteral):
                raise errors.RequireLiteralValue(given_reason1)

            if arr.ndim != 1:
                raise errors.NumbaValueError(given_reason2)

            sig = types.intp(x, arr)

            def codegen(context, builder, signature, args):
                pass
            return sig, codegen

        @njit
        def call_intrin():
            arr = np.zeros((2, 2))
            myintrin(1, arr)

        with self.assertRaises(errors.TypingError) as raises:
            call_intrin()

        excstr = str(raises.exception)
        self.assertIn(error_reasons['specific_error'].splitlines()[0], excstr)
        self.assertIn(given_reason1, excstr)
        self.assertIn(given_reason2, excstr)
        self.assertIn("Intrinsic in function", excstr)

    def test_overloadmethod_template_source(self):
        # doesn't hit _OverloadMethodTemplate for source as it's a nested
        # exception
        @overload_method(types.UnicodeType, 'isnonsense')
        def ol_unicode_isnonsense(self):
            pass

        @njit
        def foo():
            "abc".isnonsense()

        with self.assertRaises(errors.TypingError) as raises:
            foo()

        excstr = str(raises.exception)
        self.assertIn("Overload of function 'ol_unicode_isnonsense'", excstr)

    def test_overloadattribute_template_source(self):
        # doesn't hit _OverloadMethodTemplate for source as it's a nested
        # exception
        @overload_attribute(types.UnicodeType, 'isnonsense')
        def ol_unicode_isnonsense(self):
            pass

        @njit
        def foo():
            "abc".isnonsense

        with self.assertRaises(errors.TypingError) as raises:
            foo()

        excstr = str(raises.exception)
        self.assertIn("Overload of function 'ol_unicode_isnonsense'", excstr)

    def test_external_function_pointer_template_source(self):
        from numba.tests.ctypes_usecases import c_cos

        @njit
        def foo():
            c_cos('a')

        with self.assertRaises(errors.TypingError) as raises:
            foo()

        excstr = str(raises.exception)
        self.assertIn("Type Restricted Function in function 'unknown'", excstr)

    @skip_unless_cffi
    def test_cffi_function_pointer_template_source(self):
        from numba.tests import cffi_usecases as mod
        mod.init()
        func = mod.cffi_cos

        @njit
        def foo():
            func('a')

        with self.assertRaises(errors.TypingError) as raises:
            foo()

        excstr = str(raises.exception)
        self.assertIn("Type Restricted Function in function 'unknown'", excstr)

    def test_missing_source(self):

        @structref.register
        class ParticleType(types.StructRef):
            pass

        class Particle(structref.StructRefProxy):
            def __new__(cls, pos, mass):
                return structref.StructRefProxy.__new__(cls, pos)
                # didn't provide the required mass argument ----^

        structref.define_proxy(Particle, ParticleType, ["pos", "mass"])

        with self.assertRaises(errors.TypingError) as raises:
            Particle(pos=1, mass=2)

        excstr = str(raises.exception)
        self.assertIn("missing a required argument: 'mass'", excstr)


class TestDeveloperSpecificErrorMessages(SerialMixin, unittest.TestCase):

    def test_bound_function_error_string(self):
        # See PR #5952
        def foo(x):
            x.max(-1)

        with override_config('DEVELOPER_MODE', 1):
            with self.assertRaises(errors.TypingError) as raises:
                njit("void(int64[:,:])")(foo)

        excstr = str(raises.exception)
        self.assertIn("too many positional arguments", excstr)


class TestCapturedErrorHandling(SerialMixin, TestCase):
    """Checks that the way errors are captured changes depending on the env
    var "NUMBA_CAPTURED_ERRORS".
    """

    def test_error_in_overload(self):

        def bar(x):
            pass

        @overload(bar)
        def ol_bar(x):
            x.some_invalid_attr # doesn't exist!

            def impl(x):
                pass
            return impl

        with warnings.catch_warnings():
            # Suppress error going into stdout
            warnings.simplefilter("ignore",
                                  errors.NumbaPendingDeprecationWarning)
            # Check both new_style and old_style
            for style, err_class in (('new_style', AttributeError),
                                     ('old_style', errors.TypingError)):
                with override_config('CAPTURED_ERRORS', style):
                    with self.assertRaises(err_class) as raises:

                        @njit('void(int64)')
                        def foo(x):
                            bar(x)
                    expected = "object has no attribute 'some_invalid_attr'"
                    self.assertIn(expected, str(raises.exception))

    def _run_in_separate_process(self, runcode, env):
        # Run code in separate process with -Wall and specific env-vars
        code = f"""if 1:
            {runcode}\n
            """
        # On windows, missing the base environment variable can cause
        # Fatal Python error: _Py_HashRandomization_Init: failed to get random
        #                     numbers to initialize Python
        proc_env = os.environ.copy()
        proc_env.update(env)
        popen = subprocess.Popen([sys.executable, "-Wall", "-c", code],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 env=proc_env)

        out, err = popen.communicate()
        if popen.returncode != 0:
            raise AssertionError("process failed with code %s: stderr follows"
                                 "\n%s\n" % (popen.returncode, err.decode()))
        return out, err

    def test_old_style_deprecation_on_import(self):
        from numba.core.config import _old_style_deprecation_msg

        code = """
        import numba
        """
        # Check that the deprecated message is shown
        # if NUMBA_CAPTURED_ERRORS=old_style
        env = {"NUMBA_CAPTURED_ERRORS": "old_style"}
        _out, err = self._run_in_separate_process(code, env)
        self.assertIn(_old_style_deprecation_msg, err.decode())

        # Check that the deprecated message is NOT shown
        # if NUMBA_CAPTURED_ERRORS is unset
        env = {"NUMBA_CAPTURED_ERRORS": ""}
        _out, err = self._run_in_separate_process(code, env)
        # Check that the deprecated message is not shown
        self.assertNotIn("NumbaPendingDeprecationWarning", err.decode())

        # Check that the deprecated message is NOT shown
        # if NUMBA_CAPTURED_ERRORS=new_style
        env = {"NUMBA_CAPTURED_ERRORS": "new_style"}
        _out, err = self._run_in_separate_process(code, env)
        # Check that the deprecated message is not shown
        self.assertNotIn("NumbaPendingDeprecationWarning", err.decode())

    def _test_old_style_deprecation(self):
        # Verify that old_style error raise the correct deprecation warning
        warnings.simplefilter("always", errors.NumbaDeprecationWarning)

        def bar(x):
            pass

        @overload(bar)
        def ol_bar(x):
            raise AttributeError("Invalid attribute")

        with self.assertWarns(errors.NumbaDeprecationWarning) as warns:
            with self.assertRaises(errors.TypingError):
                @njit('void(int64)')
                def foo(x):
                    bar(x)

            self.assertIn(
                "Code using Numba extension API maybe depending on 'old_style' "
                "error-capturing",
                str(warns.warnings[0].message),
            )

    # Check deprecation warning when NUMBA_CAPTURED_ERRORS=old_style
    test_old_style_deprecation = TestCase.run_test_in_subprocess(
        envvars={"NUMBA_CAPTURED_ERRORS": "old_style"},
    )(_test_old_style_deprecation)

    @TestCase.run_test_in_subprocess(
        envvars={"NUMBA_CAPTURED_ERRORS": "old_style"},
    )
    def test_old_style_no_deprecation(self):
        # Verify that old_style error with NumbaError does not raise warnings
        warnings.simplefilter("always", errors.NumbaDeprecationWarning)

        def bar(x):
            pass

        @overload(bar)
        def ol_bar(x):
            raise errors.TypingError("Invalid attribute")

        with warnings.catch_warnings(record=True) as warns:
            with self.assertRaises(errors.TypingError):
                @njit('void(int64)')
                def foo(x):
                    bar(x)

            self.assertEqual(len(warns), 0,
                             msg="There should not be any warnings")

    def _test_new_style_no_warnings(self):
        # Verify that new_style error raise no warnings
        warnings.simplefilter("always", errors.NumbaDeprecationWarning)

        def bar(x):
            pass

        @overload(bar)
        def ol_bar(x):
            raise AttributeError("Invalid attribute")

        with warnings.catch_warnings(record=True) as warns:
            with self.assertRaises(AttributeError):
                @njit('void(int64)')
                def foo(x):
                    bar(x)
            # There should not be any warnings
            self.assertEqual(len(warns), 0,
                             msg="There should not be any warnings")

    test_new_style_no_warnings = TestCase.run_test_in_subprocess(
        envvars={"NUMBA_CAPTURED_ERRORS": "new_style"},
    )(_test_new_style_no_warnings)

    # Check deprecation warning when NUMBA_CAPTURED_ERRORS=default
    # ("default" means "old_style")
    test_default_new_style_no_deprecation = TestCase.run_test_in_subprocess(
        envvars={"NUMBA_CAPTURED_ERRORS": "default"},
    )(_test_new_style_no_warnings)


if __name__ == '__main__':
    unittest.main()
