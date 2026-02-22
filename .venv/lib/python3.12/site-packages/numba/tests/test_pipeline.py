from numba.core.compiler import Compiler, DefaultPassBuilder
from numba.core.compiler_machinery import (FunctionPass, AnalysisPass,
                                           register_pass)
from numba.core.untyped_passes import InlineInlinables
from numba.core.typed_passes import IRLegalization
from numba import jit, objmode, njit, cfunc
from numba.core import types, postproc, errors
from numba.core.ir import FunctionIR
from numba.tests.support import TestCase


class TestCustomPipeline(TestCase):
    def setUp(self):
        super(TestCustomPipeline, self).setUp()

        # Define custom pipeline class
        class CustomPipeline(Compiler):
            custom_pipeline_cache = []

            def compile_extra(self, func):
                # Store the compiled function
                self.custom_pipeline_cache.append(func)
                return super(CustomPipeline, self).compile_extra(func)

            def compile_ir(self, func_ir, *args, **kwargs):
                # Store the compiled function
                self.custom_pipeline_cache.append(func_ir)
                return super(CustomPipeline, self).compile_ir(
                    func_ir, *args, **kwargs)

        self.pipeline_class = CustomPipeline

    def test_jit_custom_pipeline(self):
        self.assertListEqual(self.pipeline_class.custom_pipeline_cache, [])

        @jit(pipeline_class=self.pipeline_class)
        def foo(x):
            return x

        self.assertEqual(foo(4), 4)
        self.assertListEqual(self.pipeline_class.custom_pipeline_cache,
                             [foo.py_func])

    def test_cfunc_custom_pipeline(self):
        self.assertListEqual(self.pipeline_class.custom_pipeline_cache, [])

        @cfunc(types.int64(types.int64), pipeline_class=self.pipeline_class)
        def foo(x):
            return x

        self.assertEqual(foo(4), 4)
        self.assertListEqual(self.pipeline_class.custom_pipeline_cache,
                             [foo.__wrapped__])

    def test_objmode_custom_pipeline(self):
        self.assertListEqual(self.pipeline_class.custom_pipeline_cache, [])

        @jit(pipeline_class=self.pipeline_class)
        def foo(x):
            with objmode(x="intp"):
                x += int(0x1)
            return x

        arg = 123
        self.assertEqual(foo(arg), arg + 1)
        # Two items in the list.
        self.assertEqual(len(self.pipeline_class.custom_pipeline_cache), 2)
        # First item is the `foo` function
        first = self.pipeline_class.custom_pipeline_cache[0]
        self.assertIs(first, foo.py_func)
        # Second item is a FunctionIR of the obj-lifted function
        second = self.pipeline_class.custom_pipeline_cache[1]
        self.assertIsInstance(second, FunctionIR)


class TestPassManagerFunctionality(TestCase):

    def _create_pipeline_w_del(self, base=None, inject_after=None):
        """
        Creates a new compiler pipeline with the _InjectDelsPass injected after
        the pass supplied in kwarg 'inject_after'.
        """
        self.assertTrue(inject_after is not None)
        self.assertTrue(base is not None)

        @register_pass(mutates_CFG=False, analysis_only=False)
        class _InjectDelsPass(base):
            """
            This pass injects ir.Del nodes into the IR
            """
            _name = "inject_dels_%s" % str(base)

            def __init__(self):
                base.__init__(self)

            def run_pass(self, state):
                pp = postproc.PostProcessor(state.func_ir)
                pp.run(emit_dels=True)
                return True

        class TestCompiler(Compiler):

            def define_pipelines(self):
                pm = DefaultPassBuilder.define_nopython_pipeline(self.state)
                pm.add_pass_after(_InjectDelsPass, inject_after)
                pm.finalize()
                return [pm]

        return TestCompiler

    def test_compiler_error_on_ir_del_from_functionpass(self):
        new_compiler = self._create_pipeline_w_del(FunctionPass,
                                                   InlineInlinables)

        @njit(pipeline_class=new_compiler)
        def foo(x):
            return x + 1

        with self.assertRaises(errors.CompilerError) as raises:
            foo(10)

        errstr = str(raises.exception)

        self.assertIn("Illegal IR, del found at:", errstr)
        self.assertIn("del x", errstr)

    def test_no_compiler_error_on_ir_del_after_legalization(self):
        # Legalization should be the last FunctionPass to execute so it's fine
        # for it to emit ir.Del nodes as no further FunctionPasses will run and
        # therefore the checking routine in the PassManager won't execute.
        # This test adds a new pass that is an AnalysisPass into the pipeline
        # after legalisation, this pass will return with already existing dels
        # in the IR but by virtue of it being an AnalysisPass the checking
        # routine won't execute.

        new_compiler = self._create_pipeline_w_del(AnalysisPass,
                                                   IRLegalization)

        @njit(pipeline_class=new_compiler)
        def foo(x):
            return x + 1

        self.assertTrue(foo(10), foo.py_func(10))
