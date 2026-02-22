"""
Tests for SSA reconstruction
"""
import sys
import copy
import logging

import numpy as np

from numba import njit, jit, types
from numba.core import errors, ir
from numba.core.compiler_machinery import FunctionPass, register_pass
from numba.core.compiler import DefaultPassBuilder, CompilerBase
from numba.core.untyped_passes import ReconstructSSA, PreserveIR
from numba.core.typed_passes import NativeLowering
from numba.extending import overload
from numba.tests.support import MemoryLeakMixin, TestCase, override_config


_DEBUG = False

if _DEBUG:
    # Enable debug logger on SSA reconstruction
    ssa_logger = logging.getLogger("numba.core.ssa")
    ssa_logger.setLevel(level=logging.DEBUG)
    ssa_logger.addHandler(logging.StreamHandler(sys.stderr))


class SSABaseTest(TestCase):

    def check_func(self, func, *args):
        got = func(*copy.deepcopy(args))
        exp = func.py_func(*copy.deepcopy(args))
        self.assertEqual(got, exp)


class TestSSA(SSABaseTest):
    """
    Contains tests to help isolate problems in SSA
    """

    def test_argument_name_reused(self):
        @njit
        def foo(x):
            x += 1
            return x

        self.check_func(foo, 123)

    def test_if_else_redefine(self):
        @njit
        def foo(x, y):
            z = x * y
            if x < y:
                z = x
            else:
                z = y
            return z

        self.check_func(foo, 3, 2)
        self.check_func(foo, 2, 3)

    def test_sum_loop(self):
        @njit
        def foo(n):
            c = 0
            for i in range(n):
                c += i
            return c

        self.check_func(foo, 0)
        self.check_func(foo, 10)

    def test_sum_loop_2vars(self):
        @njit
        def foo(n):
            c = 0
            d = n
            for i in range(n):
                c += i
                d += n
            return c, d

        self.check_func(foo, 0)
        self.check_func(foo, 10)

    def test_sum_2d_loop(self):
        @njit
        def foo(n):
            c = 0
            for i in range(n):
                for j in range(n):
                    c += j
                c += i
            return c

        self.check_func(foo, 0)
        self.check_func(foo, 10)

    def check_undefined_var(self, should_warn):
        @njit
        def foo(n):
            if n:
                if n > 0:
                    c = 0
                return c
            else:
                # variable c is not defined in this branch
                c += 1
                return c

        if should_warn:
            with self.assertWarns(errors.NumbaWarning) as warns:
                # n=1 so we won't actually run the branch with the uninitialized
                self.check_func(foo, 1)
            self.assertIn("Detected uninitialized variable c",
                          str(warns.warning))
        else:
            self.check_func(foo, 1)

        with self.assertRaises(UnboundLocalError):
            foo.py_func(0)

    def test_undefined_var(self):
        with override_config('ALWAYS_WARN_UNINIT_VAR', 0):
            self.check_undefined_var(should_warn=False)
        with override_config('ALWAYS_WARN_UNINIT_VAR', 1):
            self.check_undefined_var(should_warn=True)

    def test_phi_propagation(self):
        @njit
        def foo(actions):
            n = 1

            i = 0
            ct = 0
            while n > 0 and i < len(actions):
                n -= 1

                while actions[i]:
                    if actions[i]:
                        if actions[i]:
                            n += 10
                        actions[i] -= 1
                    else:
                        if actions[i]:
                            n += 20
                        actions[i] += 1

                    ct += n
                ct += n
            return ct, n

        self.check_func(foo, np.array([1, 2]))

    def test_unhandled_undefined(self):
        def function1(arg1, arg2, arg3, arg4, arg5):
            # This function is auto-generated.
            if arg1:
                var1 = arg2
                var2 = arg3
                var3 = var2
                var4 = arg1
                return
            else:
                if arg2:
                    if arg4:
                        var5 = arg4         # noqa: F841
                        return
                    else:
                        var6 = var4
                        return
                    return var6
                else:
                    if arg5:
                        if var1:
                            if arg5:
                                var1 = var6
                                return
                            else:
                                var7 = arg2     # noqa: F841
                                return arg2
                            return
                        else:
                            if var2:
                                arg5 = arg2
                                return arg1
                            else:
                                var6 = var3
                                return var4
                            return
                        return
                    else:
                        var8 = var1
                        return
                    return var8
                var9 = var3         # noqa: F841
                var10 = arg5        # noqa: F841
                return var1

        # The argument values is not critical for re-creating the bug
        # because the bug is in compile-time.
        expect = function1(2, 3, 6, 0, 7)
        got = njit(function1)(2, 3, 6, 0, 7)
        self.assertEqual(expect, got)


class TestReportedSSAIssues(SSABaseTest):
    # Tests from issues
    # https://github.com/numba/numba/issues?q=is%3Aopen+is%3Aissue+label%3ASSA

    def test_issue2194(self):

        @njit
        def foo():
            V = np.empty(1)
            s = np.uint32(1)

            for i in range(s):
                V[i] = 1
            for i in range(s, 1):
                pass

        self.check_func(foo, )

    def test_issue3094(self):

        @njit
        def doit(x):
            return x

        @njit
        def foo(pred):
            if pred:
                x = True
            else:
                x = False
            # do something with x
            return doit(x)

        self.check_func(foo, False)

    def test_issue3931(self):

        @njit
        def foo(arr):
            for i in range(1):
                arr = arr.reshape(3 * 2)
                arr = arr.reshape(3, 2)
            return (arr)

        np.testing.assert_allclose(foo(np.zeros((3, 2))),
                                   foo.py_func(np.zeros((3, 2))))

    def test_issue3976(self):

        def overload_this(a):
            return 'dummy'

        @njit
        def foo(a):
            if a:
                s = 5
                s = overload_this(s)
            else:
                s = 'b'

            return s

        @overload(overload_this)
        def ol(a):
            return overload_this

        self.check_func(foo, True)

    def test_issue3979(self):

        @njit
        def foo(A, B):
            x = A[0]
            y = B[0]
            for i in A:
                x = i
            for i in B:
                y = i
            return x, y

        self.check_func(foo, (1, 2), ('A', 'B'))

    def test_issue5219(self):

        def overload_this(a, b=None):
            if isinstance(b, tuple):
                b = b[0]
            return b

        @overload(overload_this)
        def ol(a, b=None):
            b_is_tuple = isinstance(b, (types.Tuple, types.UniTuple))

            def impl(a, b=None):
                if b_is_tuple is True:
                    b = b[0]
                return b
            return impl

        @njit
        def test_tuple(a, b):
            overload_this(a, b)

        self.check_func(test_tuple, 1, (2, ))

    def test_issue5223(self):

        @njit
        def bar(x):
            if len(x) == 5:
                return x
            x = x.copy()
            for i in range(len(x)):
                x[i] += 1
            return x

        a = np.ones(5)
        a.flags.writeable = False

        np.testing.assert_allclose(bar(a), bar.py_func(a))

    def test_issue5243(self):

        @njit
        def foo(q):
            lin = np.array((0.1, 0.6, 0.3))
            stencil = np.zeros((3, 3))
            stencil[0, 0] = q[0, 0]
            return lin[0]

        self.check_func(foo, np.zeros((2, 2)))

    def test_issue5482_missing_variable_init(self):
        # Test error that lowering fails because variable is missing
        # a definition before use.
        @njit("(intp, intp, intp)")
        def foo(x, v, n):
            for i in range(n):
                if i == 0:
                    if i == x:
                        pass
                    else:
                        problematic = v
                else:
                    if i == x:
                        pass
                    else:
                        problematic = problematic + v
            return problematic

    def test_issue5482_objmode_expr_null_lowering(self):
        # Existing pipelines will not have the Expr.null in objmode.
        # We have to create a custom pipeline to force a SSA reconstruction
        # and stripping.
        from numba.core.compiler import CompilerBase, DefaultPassBuilder
        from numba.core.untyped_passes import ReconstructSSA, IRProcessing
        from numba.core.typed_passes import PreLowerStripPhis

        class CustomPipeline(CompilerBase):
            def define_pipelines(self):
                pm = DefaultPassBuilder.define_objectmode_pipeline(self.state)
                # Force SSA reconstruction and stripping
                pm.add_pass_after(ReconstructSSA, IRProcessing)
                pm.add_pass_after(PreLowerStripPhis, ReconstructSSA)
                pm.finalize()
                return [pm]

        @jit("(intp, intp, intp)", looplift=False,
             pipeline_class=CustomPipeline)
        def foo(x, v, n):
            for i in range(n):
                if i == n:
                    if i == x:
                        pass
                    else:
                        problematic = v
                else:
                    if i == x:
                        pass
                    else:
                        problematic = problematic + v
            return problematic

    def test_issue5493_unneeded_phi(self):
        # Test error that unneeded phi is inserted because variable does not
        # have a dominance definition.
        data = (np.ones(2), np.ones(2))
        A = np.ones(1)
        B = np.ones((1,1))

        def foo(m, n, data):
            if len(data) == 1:
                v0 = data[0]
            else:
                v0 = data[0]
                # Unneeded PHI node for `problematic` would be placed here
                for _ in range(1, len(data)):
                    v0 += A

            for t in range(1, m):
                for idx in range(n):
                    t = B

                    if idx == 0:
                        if idx == n - 1:
                            pass
                        else:
                            problematic = t
                    else:
                        if idx == n - 1:
                            pass
                        else:
                            problematic = problematic + t
            return problematic

        expect = foo(10, 10, data)
        res1 = njit(foo)(10, 10, data)
        res2 = jit(forceobj=True, looplift=False)(foo)(10, 10, data)
        np.testing.assert_array_equal(expect, res1)
        np.testing.assert_array_equal(expect, res2)

    def test_issue5623_equal_statements_in_same_bb(self):

        def foo(pred, stack):
            i = 0
            c = 1

            if pred is True:
                stack[i] = c
                i += 1
                stack[i] = c
                i += 1

        python = np.array([0, 666])
        foo(True, python)

        nb = np.array([0, 666])
        njit(foo)(True, nb)

        expect = np.array([1, 1])

        np.testing.assert_array_equal(python, expect)
        np.testing.assert_array_equal(nb, expect)

    def test_issue5678_non_minimal_phi(self):
        # There should be only one phi for variable "i"

        from numba.core.compiler import CompilerBase, DefaultPassBuilder
        from numba.core.untyped_passes import (
            ReconstructSSA, FunctionPass, register_pass,
        )

        phi_counter = []

        @register_pass(mutates_CFG=False, analysis_only=True)
        class CheckSSAMinimal(FunctionPass):
            # A custom pass to count the number of phis

            _name = self.__class__.__qualname__ + ".CheckSSAMinimal"

            def __init__(self):
                super().__init__(self)

            def run_pass(self, state):
                ct = 0
                for blk in state.func_ir.blocks.values():
                    ct += len(list(blk.find_exprs('phi')))
                phi_counter.append(ct)
                return True

        class CustomPipeline(CompilerBase):
            def define_pipelines(self):
                pm = DefaultPassBuilder.define_nopython_pipeline(self.state)
                pm.add_pass_after(CheckSSAMinimal, ReconstructSSA)
                pm.finalize()
                return [pm]

        @njit(pipeline_class=CustomPipeline)
        def while_for(n, max_iter=1):
            a = np.empty((n,n))
            i = 0
            while i <= max_iter:
                for j in range(len(a)):
                    for k in range(len(a)):
                        a[j,k] = j + k
                i += 1
            return a

        # Runs fine?
        self.assertPreciseEqual(while_for(10), while_for.py_func(10))
        # One phi?
        self.assertEqual(phi_counter, [1])

    def test_issue9242_use_not_dom_def(self):
        from numba.core.ir import FunctionIR
        from numba.core.compiler_machinery import (
            AnalysisPass,
            register_pass,
        )

        def check(fir: FunctionIR):
            [blk, *_] = fir.blocks.values()
            var = blk.scope.get("d")
            defn = fir.get_definition(var)
            self.assertEqual(defn.op, "phi")
            self.assertIn(ir.UNDEFINED, defn.incoming_values)

        @register_pass(mutates_CFG=False, analysis_only=True)
        class SSACheck(AnalysisPass):
            """
            Check SSA on variable `d`
            """

            _name = "SSA_Check"

            def __init__(self):
                AnalysisPass.__init__(self)

            def run_pass(self, state):
                check(state.func_ir)
                return False

        class SSACheckPipeline(CompilerBase):
            """Inject SSACheck pass into the default pipeline following the SSA
            pass
            """

            def define_pipelines(self):
                pipeline = DefaultPassBuilder.define_nopython_pipeline(
                    self.state, "ssa_check_custom_pipeline")

                pipeline._finalized = False
                pipeline.add_pass_after(SSACheck, ReconstructSSA)

                pipeline.finalize()
                return [pipeline]

        @njit(pipeline_class=SSACheckPipeline)
        def py_func(a):
            c = a > 0
            if c:
                d = a + 5  # d is only defined here; undef in the else branch

            return c and d > 0

        py_func(10)


class TestSROAIssues(MemoryLeakMixin, TestCase):
    # This tests issues related to the SROA optimization done in lowering, which
    # reduces time spent in the LLVM SROA pass. The optimization is related to
    # SSA and tries to reduce the number of alloca statements for variables with
    # only a single assignment.
    def test_issue7258_multiple_assignment_post_SSA(self):
        # This test adds a pass that will duplicate assignment statements to
        # variables named "foobar".
        # In the reported issue, the bug will cause a memory leak.
        cloned = []

        @register_pass(analysis_only=False, mutates_CFG=True)
        class CloneFoobarAssignments(FunctionPass):
            # A pass that clones variable assignments into "foobar"
            _name = "clone_foobar_assignments_pass"

            def __init__(self):
                FunctionPass.__init__(self)

            def run_pass(self, state):
                mutated = False
                for blk in state.func_ir.blocks.values():
                    to_clone = []
                    # find assignments to "foobar"
                    for assign in blk.find_insts(ir.Assign):
                        if assign.target.name == "foobar":
                            to_clone.append(assign)
                    # clone
                    for assign in to_clone:
                        clone = copy.deepcopy(assign)
                        blk.insert_after(clone, assign)
                        mutated = True
                        # keep track of cloned statements
                        cloned.append(clone)
                return mutated

        class CustomCompiler(CompilerBase):
            def define_pipelines(self):
                pm = DefaultPassBuilder.define_nopython_pipeline(
                    self.state, "custom_pipeline",
                )
                pm._finalized = False
                # Insert the cloning pass after SSA
                pm.add_pass_after(CloneFoobarAssignments, ReconstructSSA)
                # Capture IR post lowering
                pm.add_pass_after(PreserveIR, NativeLowering)
                pm.finalize()
                return [pm]

        @njit(pipeline_class=CustomCompiler)
        def udt(arr):
            foobar = arr + 1  # this assignment will be cloned
            return foobar

        arr = np.arange(10)
        # Verify that the function works as expected
        self.assertPreciseEqual(udt(arr), arr + 1)
        # Verify that the expected statement is cloned
        self.assertEqual(len(cloned), 1)
        self.assertEqual(cloned[0].target.name, "foobar")
        # Verify in the Numba IR that the expected statement is cloned
        nir = udt.overloads[udt.signatures[0]].metadata['preserved_ir']
        self.assertEqual(len(nir.blocks), 1,
                         "only one block")
        [blk] = nir.blocks.values()
        assigns = blk.find_insts(ir.Assign)
        foobar_assigns = [stmt for stmt in assigns
                          if stmt.target.name == "foobar"]
        self.assertEqual(
            len(foobar_assigns), 2,
            "expected two assignment statements into 'foobar'",
        )
        self.assertEqual(
            foobar_assigns[0], foobar_assigns[1],
            "expected the two assignment statements to be the same",
        )
