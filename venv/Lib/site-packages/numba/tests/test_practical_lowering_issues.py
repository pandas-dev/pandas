"""
Tests for practical lowering specific errors.
"""

import numpy as np
from numba import njit
from numba.core import types, ir
from numba.core.compiler import CompilerBase, DefaultPassBuilder
from numba.core.typed_passes import NopythonTypeInference
from numba.core.compiler_machinery import register_pass, FunctionPass

from numba.tests.support import MemoryLeakMixin, TestCase


def issue7507_lround(a):
    """Dummy function used in test"""
    pass


class TestLowering(MemoryLeakMixin, TestCase):
    def test_issue4156_loop_vars_leak(self):
        """Test issues with zero-filling of refct'ed variables inside loops.

        Before the fix, the in-loop variables are always zero-filled at their
        definition location. As a result, their state from the previous
        iteration is erased. No decref is applied. To fix this, the
        zero-filling must only happen once after the alloca at the function
        entry block. The loop variables are technically defined once per
        function (one alloca per definition per function), but semantically
        defined once per assignment. Semantically, their lifetime stop only
        when the variable is re-assigned or when the function ends.
        """
        @njit
        def udt(N):
            sum_vec = np.zeros(3)
            for n in range(N):
                if n >= 0:
                    # `vec` would leak without the fix.
                    vec = np.ones(1)
                if n >= 0:
                    sum_vec += vec[0]

            return sum_vec

        got = udt(4)
        expect = udt.py_func(4)
        self.assertPreciseEqual(got, expect)

    def test_issue4156_loop_vars_leak_variant1(self):
        """Variant of test_issue4156_loop_vars_leak.

        Adding an outer loop.
        """
        @njit
        def udt(N):
            sum_vec = np.zeros(3)
            for x in range(N):
                for y in range(N):
                    n = x + y
                    if n >= 0:
                        # `vec` would leak without the fix.
                        vec = np.ones(1)
                    if n >= 0:
                        sum_vec += vec[0]

            return sum_vec

        got = udt(4)
        expect = udt.py_func(4)
        self.assertPreciseEqual(got, expect)

    def test_issue4156_loop_vars_leak_variant2(self):
        """Variant of test_issue4156_loop_vars_leak.

        Adding deeper outer loop.
        """
        @njit
        def udt(N):
            sum_vec = np.zeros(3)
            for z in range(N):
                for x in range(N):
                    for y in range(N):
                        n = x + y + z
                        if n >= 0:
                            # `vec` would leak without the fix.
                            vec = np.ones(1)
                        if n >= 0:
                            sum_vec += vec[0]

            return sum_vec

        got = udt(4)
        expect = udt.py_func(4)
        self.assertPreciseEqual(got, expect)

    def test_issue4156_loop_vars_leak_variant3(self):
        """Variant of test_issue4156_loop_vars_leak.

        Adding inner loop around allocation
        """
        @njit
        def udt(N):
            sum_vec = np.zeros(3)
            for z in range(N):
                for x in range(N):
                    n = x + z
                    if n >= 0:
                        for y in range(N):
                            # `vec` would leak without the fix.
                            vec = np.ones(y)
                    if n >= 0:
                        sum_vec += vec[0]

            return sum_vec

        got = udt(4)
        expect = udt.py_func(4)
        self.assertPreciseEqual(got, expect)

    def test_issue4156_loop_vars_leak_variant4(self):
        """Variant of test_issue4156_loop_vars_leak.

        Interleaves loops and allocations
        """
        @njit
        def udt(N):
            sum_vec = 0

            for n in range(N):
                vec = np.zeros(7)
                for n in range(N):
                    z = np.zeros(7)
                sum_vec += vec[0] + z[0]

            return sum_vec

        got = udt(4)
        expect = udt.py_func(4)
        self.assertPreciseEqual(got, expect)

    def test_issue_with_literal_in_static_getitem(self):
        """Test an issue with literal type used as index of static_getitem
        """

        @register_pass(mutates_CFG=False, analysis_only=False)
        class ForceStaticGetitemLiteral(FunctionPass):

            _name = "force_static_getitem_literal"

            def __init__(self):
                FunctionPass.__init__(self)

            def run_pass(self, state):
                repl = {}
                # Force the static_getitem to have a literal type as
                # index to replicate the problem.
                for inst, sig in state.calltypes.items():
                    if (isinstance(inst, ir.Expr) and
                            inst.op == 'static_getitem'):
                        [obj, idx] = sig.args
                        new_sig = sig.replace(args=(obj,
                                                    types.literal(inst.index)))
                        repl[inst] = new_sig
                state.calltypes.update(repl)
                return True

        class CustomPipeline(CompilerBase):
            def define_pipelines(self):
                pm = DefaultPassBuilder.define_nopython_pipeline(self.state)
                pm.add_pass_after(ForceStaticGetitemLiteral,
                                  NopythonTypeInference)
                pm.finalize()
                return [pm]

        @njit(pipeline_class=CustomPipeline)
        def foo(arr):
            return arr[4]  # force static_getitem

        arr = np.arange(10)
        got = foo(arr)
        expect = foo.py_func(arr)
        self.assertEqual(got, expect)

    def test_issue7507(self):
        """
        Test a problem with BaseContext.get_function() because of changes
        related to the new style error handling.
        """
        from numba.core.typing.templates import AbstractTemplate, infer_global
        from numba.core.imputils import lower_builtin

        @infer_global(issue7507_lround)
        class lroundTemplate(AbstractTemplate):
            key = issue7507_lround

            def generic(self, args, kws):
                signature = types.int64(types.float64)

                # insert a new builtin during the compilation process
                @lower_builtin(issue7507_lround, types.float64)
                def codegen(context, builder, sig, args):
                    # Simply truncate with the cast to integer.
                    return context.cast(builder, args[0], sig.args[0],
                                        sig.return_type)

                return signature

        @njit('int64(float64)')
        def foo(a):
            return issue7507_lround(a)

        self.assertEqual(foo(3.4), 3)
