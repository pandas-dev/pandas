from collections import namedtuple
import numpy as np

from numba.tests.support import (TestCase, MemoryLeakMixin,
                                 skip_parfors_unsupported, captured_stdout)
from numba import njit, typed, literal_unroll, prange
from numba.core import types, errors, ir
from numba.testing import unittest
from numba.core.extending import overload
from numba.core.compiler_machinery import (PassManager, register_pass,
                                           FunctionPass, AnalysisPass)
from numba.core.compiler import CompilerBase
from numba.core.untyped_passes import (FixupArgs, TranslateByteCode,
                                       IRProcessing, InlineClosureLikes,
                                       SimplifyCFG, IterLoopCanonicalization,
                                       LiteralUnroll, PreserveIR)
from numba.core.typed_passes import (NopythonTypeInference, IRLegalization,
                                     NoPythonBackend, PartialTypeInference,
                                     NativeLowering)
from numba.core.ir_utils import (compute_cfg_from_blocks, flatten_labels)
from numba.core.types.functions import _header_lead

_X_GLOBAL = (10, 11)


class TestLiteralTupleInterpretation(MemoryLeakMixin, TestCase):

    def check(self, func, var):
        cres = func.overloads[func.signatures[0]]
        ty = cres.fndesc.typemap[var]
        self.assertTrue(isinstance(ty, types.Tuple))
        for subty in ty:
            self.assertTrue(isinstance(subty, types.Literal), "non literal")

    def test_homogeneous_literal(self):
        @njit
        def foo():
            x = (1, 2, 3)
            return x[1]

        self.assertEqual(foo(), foo.py_func())
        self.check(foo, 'x')

    def test_heterogeneous_literal(self):
        @njit
        def foo():
            x = (1, 2, 3, 'a')
            return x[3]

        self.assertEqual(foo(), foo.py_func())
        self.check(foo, 'x')

    def test_non_literal(self):
        @njit
        def foo():
            x = (1, 2, 3, 'a', 1j)
            return x[4]

        self.assertEqual(foo(), foo.py_func())
        with self.assertRaises(AssertionError) as e:
            self.check(foo, 'x')

        self.assertIn("non literal", str(e.exception))


@register_pass(mutates_CFG=False, analysis_only=False)
class ResetTypeInfo(FunctionPass):
    _name = "reset_the_type_information"

    def __init__(self):
        FunctionPass.__init__(self)

    def run_pass(self, state):
        state.typemap = None
        state.return_type = None
        state.calltypes = None
        return True


class TestLoopCanonicalisation(MemoryLeakMixin, TestCase):

    def get_pipeline(use_canonicaliser, use_partial_typing=False):
        class NewCompiler(CompilerBase):

            def define_pipelines(self):
                pm = PassManager("custom_pipeline")

                # untyped
                pm.add_pass(TranslateByteCode, "analyzing bytecode")
                pm.add_pass(IRProcessing, "processing IR")
                pm.add_pass(InlineClosureLikes,
                            "inline calls to locally defined closures")
                if use_partial_typing:
                    pm.add_pass(PartialTypeInference, "do partial typing")
                if use_canonicaliser:
                    pm.add_pass(IterLoopCanonicalization, "Canonicalise loops")
                pm.add_pass(SimplifyCFG, "Simplify the CFG")

                # typed
                if use_partial_typing:
                    pm.add_pass(ResetTypeInfo, "resets the type info state")

                pm.add_pass(NopythonTypeInference, "nopython frontend")

                # legalise
                pm.add_pass(IRLegalization, "ensure IR is legal")

                # preserve
                pm.add_pass(PreserveIR, "save IR for later inspection")

                # lower
                pm.add_pass(NativeLowering, "native lowering")
                pm.add_pass(NoPythonBackend, "nopython mode backend")

                # finalise the contents
                pm.finalize()

                return [pm]
        return NewCompiler

    # generate variants
    LoopIgnoringCompiler = get_pipeline(False)
    LoopCanonicalisingCompiler = get_pipeline(True)
    TypedLoopCanonicalisingCompiler = get_pipeline(True, True)

    def test_simple_loop_in_depth(self):
        """ This heavily checks a simple loop transform """

        def get_info(pipeline):
            @njit(pipeline_class=pipeline)
            def foo(tup):
                acc = 0
                for i in tup:
                    acc += i
                return acc

            x = (1, 2, 3)
            self.assertEqual(foo(x), foo.py_func(x))
            cres = foo.overloads[foo.signatures[0]]
            func_ir = cres.metadata['preserved_ir']
            return func_ir, cres.fndesc

        ignore_loops_ir, ignore_loops_fndesc = \
            get_info(self.LoopIgnoringCompiler)
        canonicalise_loops_ir, canonicalise_loops_fndesc = \
            get_info(self.LoopCanonicalisingCompiler)

        # check CFG is the same
        def compare_cfg(a, b):
            a_cfg = compute_cfg_from_blocks(flatten_labels(a.blocks))
            b_cfg = compute_cfg_from_blocks(flatten_labels(b.blocks))
            self.assertEqual(a_cfg, b_cfg)

        compare_cfg(ignore_loops_ir, canonicalise_loops_ir)

        # check there's three more call types in the canonicalised one:
        # len(tuple arg)
        # range(of the len() above)
        # getitem(tuple arg, index)
        self.assertEqual(len(ignore_loops_fndesc.calltypes) + 3,
                         len(canonicalise_loops_fndesc.calltypes))

        def find_getX(fd, op):
            return [x for x in fd.calltypes.keys()
                    if isinstance(x, ir.Expr) and x.op == op]

        il_getiters = find_getX(ignore_loops_fndesc, "getiter")
        self.assertEqual(len(il_getiters), 1)  # tuple iterator

        cl_getiters = find_getX(canonicalise_loops_fndesc, "getiter")
        self.assertEqual(len(cl_getiters), 1)  # loop range iterator

        cl_getitems = find_getX(canonicalise_loops_fndesc, "getitem")
        self.assertEqual(len(cl_getitems), 1)  # tuple getitem induced by loop

        # check the value of the untransformed IR getiter is now the value of
        # the transformed getitem
        self.assertEqual(il_getiters[0].value.name, cl_getitems[0].value.name)

        # check the type of the transformed IR getiter is a range iter
        range_inst = canonicalise_loops_fndesc.calltypes[cl_getiters[0]].args[0]
        self.assertTrue(isinstance(range_inst, types.RangeType))

    def test_transform_scope(self):
        """ This checks the transform, when there's no typemap, will happily
        transform a loop on something that's not tuple-like
        """
        def get_info(pipeline):
            @njit(pipeline_class=pipeline)
            def foo():
                acc = 0
                for i in [1, 2, 3]:
                    acc += i
                return acc

            self.assertEqual(foo(), foo.py_func())
            cres = foo.overloads[foo.signatures[0]]
            func_ir = cres.metadata['preserved_ir']
            return func_ir, cres.fndesc

        ignore_loops_ir, ignore_loops_fndesc = \
            get_info(self.LoopIgnoringCompiler)
        canonicalise_loops_ir, canonicalise_loops_fndesc = \
            get_info(self.LoopCanonicalisingCompiler)

        # check CFG is the same
        def compare_cfg(a, b):
            a_cfg = compute_cfg_from_blocks(flatten_labels(a.blocks))
            b_cfg = compute_cfg_from_blocks(flatten_labels(b.blocks))
            self.assertEqual(a_cfg, b_cfg)

        compare_cfg(ignore_loops_ir, canonicalise_loops_ir)

        # check there's three more call types in the canonicalised one:
        # len(literal list)
        # range(of the len() above)
        # getitem(literal list arg, index)
        self.assertEqual(len(ignore_loops_fndesc.calltypes) + 3,
                         len(canonicalise_loops_fndesc.calltypes))

        def find_getX(fd, op):
            return [x for x in fd.calltypes.keys()
                    if isinstance(x, ir.Expr) and x.op == op]

        il_getiters = find_getX(ignore_loops_fndesc, "getiter")
        self.assertEqual(len(il_getiters), 1)  # list iterator

        cl_getiters = find_getX(canonicalise_loops_fndesc, "getiter")
        self.assertEqual(len(cl_getiters), 1)  # loop range iterator

        cl_getitems = find_getX(canonicalise_loops_fndesc, "getitem")
        self.assertEqual(len(cl_getitems), 1)  # list getitem induced by loop

        # check the value of the untransformed IR getiter is now the value of
        # the transformed getitem
        self.assertEqual(il_getiters[0].value.name, cl_getitems[0].value.name)

        # check the type of the transformed IR getiter is a range iter
        range_inst = canonicalise_loops_fndesc.calltypes[cl_getiters[0]].args[0]
        self.assertTrue(isinstance(range_inst, types.RangeType))

    @unittest.skip("Waiting for pass to be enabled for all tuples")
    def test_influence_of_typed_transform(self):
        """ This heavily checks a typed transformation only impacts tuple
        induced loops"""

        def get_info(pipeline):
            @njit(pipeline_class=pipeline)
            def foo(tup):
                acc = 0
                for i in range(4):
                    for y in tup:
                        for j in range(3):
                            acc += 1
                return acc

            x = (1, 2, 3)
            self.assertEqual(foo(x), foo.py_func(x))
            cres = foo.overloads[foo.signatures[0]]
            func_ir = cres.metadata['func_ir']
            return func_ir, cres.fndesc

        ignore_loops_ir, ignore_loops_fndesc = \
            get_info(self.LoopIgnoringCompiler)
        canonicalise_loops_ir, canonicalise_loops_fndesc = \
            get_info(self.TypedLoopCanonicalisingCompiler)

        # check CFG is the same
        def compare_cfg(a, b):
            a_cfg = compute_cfg_from_blocks(flatten_labels(a.blocks))
            b_cfg = compute_cfg_from_blocks(flatten_labels(b.blocks))
            self.assertEqual(a_cfg, b_cfg)

        compare_cfg(ignore_loops_ir, canonicalise_loops_ir)

        # check there's three more call types in the canonicalised one:
        # len(tuple arg)
        # range(of the len() above)
        # getitem(tuple arg, index)
        self.assertEqual(len(ignore_loops_fndesc.calltypes) + 3,
                         len(canonicalise_loops_fndesc.calltypes))

        def find_getX(fd, op):
            return [x for x in fd.calltypes.keys()
                    if isinstance(x, ir.Expr) and x.op == op]

        il_getiters = find_getX(ignore_loops_fndesc, "getiter")
        self.assertEqual(len(il_getiters), 3)  # 1 * tuple + 2 * loop range

        cl_getiters = find_getX(canonicalise_loops_fndesc, "getiter")
        self.assertEqual(len(cl_getiters), 3)  # 3 * loop range iterator

        cl_getitems = find_getX(canonicalise_loops_fndesc, "getitem")
        self.assertEqual(len(cl_getitems), 1)  # tuple getitem induced by loop

        # check the value of the untransformed IR getiter is now the value of
        # the transformed getitem
        self.assertEqual(il_getiters[1].value.name, cl_getitems[0].value.name)

        # check the type of the transformed IR getiter's are all range iter
        for x in cl_getiters:
            range_inst = canonicalise_loops_fndesc.calltypes[x].args[0]
            self.assertTrue(isinstance(range_inst, types.RangeType))

    def test_influence_of_typed_transform_literal_unroll(self):
        """ This heavily checks a typed transformation only impacts loops with
        literal_unroll marker"""

        def get_info(pipeline):
            @njit(pipeline_class=pipeline)
            def foo(tup):
                acc = 0
                for i in range(4):
                    for y in literal_unroll(tup):
                        for j in range(3):
                            acc += 1
                return acc

            x = (1, 2, 3)
            self.assertEqual(foo(x), foo.py_func(x))
            cres = foo.overloads[foo.signatures[0]]
            func_ir = cres.metadata['preserved_ir']
            return func_ir, cres.fndesc

        ignore_loops_ir, ignore_loops_fndesc = \
            get_info(self.LoopIgnoringCompiler)
        canonicalise_loops_ir, canonicalise_loops_fndesc = \
            get_info(self.TypedLoopCanonicalisingCompiler)

        # check CFG is the same
        def compare_cfg(a, b):
            a_cfg = compute_cfg_from_blocks(flatten_labels(a.blocks))
            b_cfg = compute_cfg_from_blocks(flatten_labels(b.blocks))
            self.assertEqual(a_cfg, b_cfg)

        compare_cfg(ignore_loops_ir, canonicalise_loops_ir)

        # check there's three more call types in the canonicalised one:
        # len(tuple arg)
        # range(of the len() above)
        # getitem(tuple arg, index)
        self.assertEqual(len(ignore_loops_fndesc.calltypes) + 3,
                         len(canonicalise_loops_fndesc.calltypes))

        def find_getX(fd, op):
            return [x for x in fd.calltypes.keys()
                    if isinstance(x, ir.Expr) and x.op == op]

        il_getiters = find_getX(ignore_loops_fndesc, "getiter")
        self.assertEqual(len(il_getiters), 3)  # 1 * tuple + 2 * loop range

        cl_getiters = find_getX(canonicalise_loops_fndesc, "getiter")
        self.assertEqual(len(cl_getiters), 3)  # 3 * loop range iterator

        cl_getitems = find_getX(canonicalise_loops_fndesc, "getitem")
        self.assertEqual(len(cl_getitems), 1)  # tuple getitem induced by loop

        # check the value of the untransformed IR getiter is now the value of
        # the transformed getitem
        self.assertEqual(il_getiters[1].value.name, cl_getitems[0].value.name)

        # check the type of the transformed IR getiter's are all range iter
        for x in cl_getiters:
            range_inst = canonicalise_loops_fndesc.calltypes[x].args[0]
            self.assertTrue(isinstance(range_inst, types.RangeType))

    @unittest.skip("Waiting for pass to be enabled for all tuples")
    def test_lots_of_loops(self):
        """ This heavily checks a simple loop transform """

        def get_info(pipeline):
            @njit(pipeline_class=pipeline)
            def foo(tup):
                acc = 0
                for i in tup:
                    acc += i
                    for j in tup + (4, 5, 6):
                        acc += 1 - j
                        if j > 5:
                            break
                    else:
                        acc -= 2
                for i in tup:
                    acc -= i % 2

                return acc

            x = (1, 2, 3)
            self.assertEqual(foo(x), foo.py_func(x))
            cres = foo.overloads[foo.signatures[0]]
            func_ir = cres.metadata['preserved_ir']
            return func_ir, cres.fndesc

        ignore_loops_ir, ignore_loops_fndesc = \
            get_info(self.LoopIgnoringCompiler)
        canonicalise_loops_ir, canonicalise_loops_fndesc = \
            get_info(self.LoopCanonicalisingCompiler)

        # check CFG is the same
        def compare_cfg(a, b):
            a_cfg = compute_cfg_from_blocks(flatten_labels(a.blocks))
            b_cfg = compute_cfg_from_blocks(flatten_labels(b.blocks))
            self.assertEqual(a_cfg, b_cfg)

        compare_cfg(ignore_loops_ir, canonicalise_loops_ir)

        # check there's three * N more call types in the canonicalised one:
        # len(tuple arg)
        # range(of the len() above)
        # getitem(tuple arg, index)
        self.assertEqual(len(ignore_loops_fndesc.calltypes) + 3 * 3,
                         len(canonicalise_loops_fndesc.calltypes))

    def test_inlined_loops(self):
        """ Checks a loop appearing from a closure """

        def get_info(pipeline):
            @njit(pipeline_class=pipeline)
            def foo(tup):
                def bar(n):
                    acc = 0
                    for i in range(n):
                        acc += 1
                    return acc

                acc = 0
                for i in tup:
                    acc += i
                    acc += bar(i)

                return acc

            x = (1, 2, 3)
            self.assertEqual(foo(x), foo.py_func(x))
            cres = foo.overloads[foo.signatures[0]]
            func_ir = cres.metadata['preserved_ir']
            return func_ir, cres.fndesc

        ignore_loops_ir, ignore_loops_fndesc = \
            get_info(self.LoopIgnoringCompiler)
        canonicalise_loops_ir, canonicalise_loops_fndesc = \
            get_info(self.LoopCanonicalisingCompiler)

        # check CFG is the same
        def compare_cfg(a, b):
            a_cfg = compute_cfg_from_blocks(flatten_labels(a.blocks))
            b_cfg = compute_cfg_from_blocks(flatten_labels(b.blocks))
            self.assertEqual(a_cfg, b_cfg)

        compare_cfg(ignore_loops_ir, canonicalise_loops_ir)

        # check there's 2 * N - 1 more call types in the canonicalised one:
        # The -1 comes from the closure being inlined and and the call removed.
        # len(tuple arg)
        # range(of the len() above)
        # getitem(tuple arg, index)
        self.assertEqual(len(ignore_loops_fndesc.calltypes) + 5,
                         len(canonicalise_loops_fndesc.calltypes))


class TestMixedTupleUnroll(MemoryLeakMixin, TestCase):

    def test_01(self):
        # test a case which is already in loop canonical form
        @njit
        def foo(idx, z):
            a = (12, 12.7, 3j, 4, z, 2 * z)
            acc = 0
            for i in range(len(literal_unroll(a))):
                acc += a[i]
                if acc.real < 26:
                    acc -= 1
                else:
                    break
            return acc

        f = 9
        k = f

        self.assertEqual(foo(2, k), foo.py_func(2, k))

    def test_02(self):
        # same as test_1 but without the explicit loop canonicalisation

        @njit
        def foo(idx, z):
            x = (12, 12.7, 3j, 4, z, 2 * z)
            acc = 0
            for a in literal_unroll(x):
                acc += a
                if acc.real < 26:
                    acc -= 1
                else:
                    break
            return acc

        f = 9
        k = f

        self.assertEqual(foo(2, k), foo.py_func(2, k))

    def test_03(self):
        # two unrolls
        @njit
        def foo(idx, z):
            x = (12, 12.7, 3j, 4, z, 2 * z)
            y = ('foo', z, 2 * z)
            acc = 0
            for a in literal_unroll(x):
                acc += a
                if acc.real < 26:
                    acc -= 1
                else:
                    for t in literal_unroll(y):
                        acc += t is False
                    break
            return acc

        f = 9
        k = f

        self.assertEqual(foo(2, k), foo.py_func(2, k))

    def test_04(self):
        # mixed ref counted types
        @njit
        def foo(tup):
            acc = 0
            for a in literal_unroll(tup):
                acc += a.sum()
            return acc

        n = 10
        tup = (np.ones((n,)), np.ones((n, n)), np.ones((n, n, n)))
        self.assertEqual(foo(tup), foo.py_func(tup))

    def test_05(self):
        # mix unroll and static_getitem
        @njit
        def foo(tup1, tup2):
            acc = 0
            for a in literal_unroll(tup1):
                if a == 'a':
                    acc += tup2[0].sum()
                elif a == 'b':
                    acc += tup2[1].sum()
                elif a == 'c':
                    acc += tup2[2].sum()
                elif a == 12:
                    acc += tup2[3].sum()
                elif a == 3j:
                    acc += tup2[4].sum()
                else:
                    raise RuntimeError("Unreachable")
            return acc

        n = 10
        tup1 = ('a', 'b', 'c', 12, 3j,)
        tup2 = (np.ones((n,)), np.ones((n, n)), np.ones((n, n, n)),
                np.ones((n, n, n, n)), np.ones((n, n, n, n, n)))
        self.assertEqual(foo(tup1, tup2), foo.py_func(tup1, tup2))

    @unittest.skip("needs more clever branch prune")
    def test_06(self):
        # This wont work because both sides of the branch need typing as neither
        # can be pruned by the current pruner
        @njit
        def foo(tup):
            acc = 0
            str_buf = typed.List.empty_list(types.unicode_type)
            for a in literal_unroll(tup):
                if a == 'a':
                    str_buf.append(a)
                else:
                    acc += a
            return acc

        tup = ('a', 12)
        self.assertEqual(foo(tup), foo.py_func(tup))

    def test_07(self):
        # A mix bag of stuff as an arg to a function that unifies as `intp`.
        @njit
        def foo(tup):
            acc = 0
            for a in literal_unroll(tup):
                acc += len(a)
            return acc

        n = 10
        tup = (np.ones((n,)), np.ones((n, n)), "ABCDEFGHJI", (1, 2, 3),
               (1, 'foo', 2, 'bar'), {3, 4, 5, 6, 7})
        self.assertEqual(foo(tup), foo.py_func(tup))

    def test_08(self):
        # dispatch to functions

        @njit
        def foo(tup1, tup2):
            acc = 0
            for a in literal_unroll(tup1):
                if a == 'a':
                    acc += tup2[0]()
                elif a == 'b':
                    acc += tup2[1]()
                elif a == 'c':
                    acc += tup2[2]()
            return acc

        def gen(x):
            def impl():
                return x
            return njit(impl)

        tup1 = ('a', 'b', 'c', 12, 3j, ('f',))
        tup2 = (gen(1), gen(2), gen(3))
        self.assertEqual(foo(tup1, tup2), foo.py_func(tup1, tup2))

    def test_09(self):
        # illegal RHS, has a mixed tuple being index dynamically

        @njit
        def foo(tup1, tup2):
            acc = 0
            idx = 0
            for a in literal_unroll(tup1):
                if a == 'a':
                    acc += tup2[idx]
                elif a == 'b':
                    acc += tup2[idx]
                elif a == 'c':
                    acc += tup2[idx]
                idx += 1
            return idx, acc

        @njit
        def func1():
            return 1

        @njit
        def func2():
            return 2

        @njit
        def func3():
            return 3

        tup1 = ('a', 'b', 'c')
        tup2 = (1j, 1, 2)

        with self.assertRaises(errors.TypingError) as raises:
            foo(tup1, tup2)

        self.assertIn(_header_lead, str(raises.exception))

    def test_10(self):
        # dispatch on literals triggering @overload resolution

        def dt(value):
            if value == "apple":
                return 1
            elif value == "orange":
                return 2
            elif value == "banana":
                return 3
            elif value == 0xca11ab1e:
                return 0x5ca1ab1e + value

        @overload(dt, inline='always')
        def ol_dt(li):
            if isinstance(li, types.StringLiteral):
                value = li.literal_value
                if value == "apple":
                    def impl(li):
                        return 1
                elif value == "orange":
                    def impl(li):
                        return 2
                elif value == "banana":
                    def impl(li):
                        return 3
                return impl
            elif isinstance(li, types.IntegerLiteral):
                value = li.literal_value
                if value == 0xca11ab1e:
                    def impl(li):
                        # close over the dispatcher :)
                        return 0x5ca1ab1e + value
                    return impl

        @njit
        def foo():
            acc = 0
            for t in literal_unroll(('apple', 'orange', 'banana', 3390155550)):
                acc += dt(t)
            return acc

        self.assertEqual(foo(), foo.py_func())

    def test_11(self):

        @njit
        def foo():
            x = []
            z = ('apple', 'orange', 'banana')
            for i in range(len(literal_unroll(z))):
                t = z[i]
                if t == "apple":
                    x.append("0")
                elif t == "orange":
                    x.append(t)
                elif t == "banana":
                    x.append("2.0")
            return x

        self.assertEqual(foo(), foo.py_func())

    def test_11a(self):

        @njit
        def foo():
            x = typed.List()
            z = ('apple', 'orange', 'banana')
            for i in range(len(literal_unroll(z))):
                t = z[i]
                if t == "apple":
                    x.append("0")
                elif t == "orange":
                    x.append(t)
                elif t == "banana":
                    x.append("2.0")
            return x

        self.assertEqual(foo(), foo.py_func())

    def test_12(self):
        # unroll the same target twice
        @njit
        def foo(idx, z):
            a = (12, 12.7, 3j, 4, z, 2 * z)
            acc = 0
            for i in literal_unroll(a):
                acc += i
                if acc.real < 26:
                    acc -= 1
                else:
                    for x in literal_unroll(a):
                        acc += x
                    break
            if a[0] < 23:
                acc += 2
            return acc

        f = 9
        k = f

        self.assertEqual(foo(2, k), foo.py_func(2, k))

    def test_13(self):
        # nesting unrolls is illegal
        @njit
        def foo(idx, z):
            a = (12, 12.7, 3j, 4, z, 2 * z)
            acc = 0
            for i in literal_unroll(a):
                acc += i
                if acc.real < 26:
                    acc -= 1
                else:
                    for x in literal_unroll(a):
                        for j in literal_unroll(a):
                            acc += j
                        acc += x
                for x in literal_unroll(a):
                    acc += x
            for x in literal_unroll(a):
                acc += x
            if a[0] < 23:
                acc += 2
            return acc

        f = 9
        k = f

        with self.assertRaises(errors.UnsupportedError) as raises:
            foo(2, k)

        self.assertIn("Nesting of literal_unroll is unsupported",
                      str(raises.exception))

    def test_14(self):
        # unituple unroll can return derivative of the induction var

        @njit
        def foo():
            x = (1, 2, 3, 4)
            acc = 0
            for a in literal_unroll(x):
                acc += a
            return a

        self.assertEqual(foo(), foo.py_func())

    def test_15(self):
        # mixed tuple unroll cannot return derivative of the induction var

        @njit
        def foo(x):
            acc = 0
            for a in literal_unroll(x):
                acc += len(a)
            return a

        n = 5
        tup = (np.ones((n,)), np.ones((n, n)), "ABCDEFGHJI", (1, 2, 3),
               (1, 'foo', 2, 'bar'), {3, 4, 5, 6, 7})

        with self.assertRaises(errors.TypingError) as raises:
            foo(tup)

        self.assertIn("Cannot unify", str(raises.exception))

    def test_16(self):
        # unituple slice and unroll is ok

        def dt(value):
            if value == 1000:
                return "a"
            elif value == 2000:
                return "b"
            elif value == 3000:
                return "c"
            elif value == 4000:
                return "d"

        @overload(dt, inline='always')
        def ol_dt(li):
            if isinstance(li, types.IntegerLiteral):
                value = li.literal_value
                if value == 1000:
                    def impl(li):
                        return "a"
                elif value == 2000:
                    def impl(li):
                        return "b"
                elif value == 3000:
                    def impl(li):
                        return "c"
                elif value == 4000:
                    def impl(li):
                        return "d"
                return impl

        @njit
        def foo():
            x = (1000, 2000, 3000, 4000)
            acc = ""
            for a in literal_unroll(x[:2]):
                acc += dt(a)
            return acc

        self.assertEqual(foo(), foo.py_func())

    def test_17(self):
        # mixed tuple slice and unroll is ok

        def dt(value):
            if value == 1000:
                return "a"
            elif value == 2000:
                return "b"
            elif value == 3000:
                return "c"
            elif value == 4000:
                return "d"
            elif value == 'f':
                return "EFF"

        @overload(dt, inline='always')
        def ol_dt(li):
            if isinstance(li, types.IntegerLiteral):
                value = li.literal_value
                if value == 1000:
                    def impl(li):
                        return "a"
                elif value == 2000:
                    def impl(li):
                        return "b"
                elif value == 3000:
                    def impl(li):
                        return "c"
                elif value == 4000:
                    def impl(li):
                        return "d"
                return impl
            elif isinstance(li, types.StringLiteral):
                value = li.literal_value
                if value == 'f':
                    def impl(li):
                        return "EFF"
                    return impl

        @njit
        def foo():
            x = (1000, 2000, 3000, 'f')
            acc = ""
            for a in literal_unroll(x[1:]):
                acc += dt(a)
            return acc

        self.assertEqual(foo(), foo.py_func())

    def test_18(self):
        # unituple backwards slice
        @njit
        def foo():
            x = (1000, 2000, 3000, 4000, 5000, 6000)
            count = 0
            for a in literal_unroll(x[::-1]):
                count += 1
                if a < 3000:
                    break
            return count

        self.assertEqual(foo(), foo.py_func())

    def test_19(self):
        # mixed bag of refcounted
        @njit
        def foo():
            acc = 0
            l1 = [1, 2, 3, 4]
            l2 = [10, 20]
            tup = (l1, l2)
            a1 = np.arange(20)
            a2 = np.ones(5, dtype=np.complex128)
            tup = (l1, a1, l2, a2)
            for t in literal_unroll(tup):
                acc += len(t)
            return acc

        self.assertEqual(foo(), foo.py_func())

    def test_20(self):
        # testing partial type inference survives as the list append in the
        # unrolled version is full inferable
        @njit
        def foo():
            l = []
            a1 = np.arange(20)
            a2 = np.ones(5, dtype=np.complex128)
            tup = (a1, a2)
            for t in literal_unroll(tup):
                l.append(t.sum())
            return l

        self.assertEqual(foo(), foo.py_func())

    def test_21(self):
        # unroll in closure that gets inlined
        @njit
        def foo(z):
            b = (23, 23.9, 6j, 8)

            def bar():
                acc = 0
                for j in literal_unroll(b):
                    acc += j
                return acc
            outer_acc = 0
            for x in (1, 2, 3, 4):
                outer_acc += bar() + x

            return outer_acc

        f = 9
        k = f
        self.assertEqual(foo(k), foo.py_func(k))

    def test_22(self):
        # NOTE: This test "worked" as a side effect of a bug that was discovered
        # during work that added support for Python 3.12. In the literal_unroll
        # transform, there was an accidental overwrite of a dictionary key that
        # could occur if nested unrolls happened to have the conditions:
        #
        # 1. outer unroll induction varible not used in the induced loop.
        # 2. the `max_label` from numba.core.ir_utils was in a specific state
        #    such that a collision occurred.
        # 3. the bug existed in the algorithm that checked for getitem access.
        #
        # This exceedingly rare usecase is now considered illegal as a result,
        # by banning this behaviour the analysis is considerably more simple.
        # Further there is no loss of functionality as the outer loop can be
        # trivially replaced (in the following) by using a `range(len())` based
        # iteration over `a` as there's no need for the loop in `a` to be
        # versioned!

        @njit
        def foo(z):
            a = (12, 12.7, 3j, 4, z, 2 * z, 'a')
            b = (23, 23.9, 6j, 8)

            def bar():
                acc = 0
                for j in literal_unroll(b):
                    acc += j
                return acc
            acc = 0
            # this loop is induced in `x` but `x` is not used, there is a nest
            # here by virtue of inlining
            for x in literal_unroll(a):
                acc += bar()

            return acc

        f = 9
        k = f

        with self.assertRaises(errors.UnsupportedError) as raises:
            foo(k)

        self.assertIn("Nesting of literal_unroll is unsupported",
                      str(raises.exception))

    def test_23(self):
        # unroll from closure that ends up banned as it leads to nesting
        @njit
        def foo(z):
            b = (23, 23.9, 6j, 8)

            def bar():
                acc = 0
                for j in literal_unroll(b):
                    acc += j
                return acc
            outer_acc = 0
            # this drives an inlined literal_unroll loop but also has access to
            # the induction variable, this is a nested literal_unroll so is
            # banned
            for x in literal_unroll(b):
                outer_acc += bar() + x

            return outer_acc

        f = 9
        k = f

        with self.assertRaises(errors.UnsupportedError) as raises:
            foo(k)

        self.assertIn("Nesting of literal_unroll is unsupported",
                      str(raises.exception))

    def test_24(self):
        # unroll something unsupported
        @njit
        def foo():
            for x in literal_unroll("ABCDE"):
                print(x)

        with self.assertRaises(errors.UnsupportedError) as raises:
            foo()

        msg = "argument should be a tuple or a list of constant values"
        self.assertIn(msg, str(raises.exception))

    def test_25(self):
        # use unroll by reference/alias
        @njit
        def foo():
            val = literal_unroll(((1, 2, 3), (2j, 3j), [1, 2], "xyz"))
            alias1 = val
            alias2 = alias1
            lens = []
            for x in alias2:
                lens.append(len(x))
            return lens

        self.assertEqual(foo(), foo.py_func())

    def test_26(self):
        # var defined in unrolled body escapes
        # untouched variable is untouched
        # read only variable is only read
        # mutated is muted correctly
        @njit
        def foo(z):
            a = (12, 12.7, 3j, 4, z, 2 * z)
            acc = 0
            count = 0
            untouched = 54
            read_only = 17
            mutated = np.empty((len(a),), dtype=np.complex128)
            for x in literal_unroll(a):
                acc += x
                mutated[count] = x
                count += 1
                escape = count + read_only
            return escape, acc, untouched, read_only, mutated

        f = 9
        k = f

        self.assertPreciseEqual(foo(k), foo.py_func(k))

    @skip_parfors_unsupported
    def test_27(self):
        # parfors loop in unrolled loop
        @njit(parallel=True)
        def foo(z):
            a = (12, 12.7, 3j, 4, z, 2 * z)
            acc = 0
            for x in literal_unroll(a):
                for k in prange(10):
                    acc += 1
            return acc

        f = 9
        k = f

        self.assertEqual(foo(k), foo.py_func(k))

    @skip_parfors_unsupported
    def test_28(self):
        # parfors reducing on the unrolled induction var
        @njit(parallel=True)
        def foo(z):
            a = (12, 12.7, 3j, 4, z, 2 * z)
            acc = 0
            for x in literal_unroll(a):
                for k in prange(10):
                    acc += x
            return acc

        f = 9
        k = f

        # summation is unstable
        np.testing.assert_allclose(foo(k), foo.py_func(k))

    @skip_parfors_unsupported
    def test_29(self):
        # This "works" but parfors is not producing a parallel loop
        # TODO: fix
        @njit(parallel=True)
        def foo(z):
            a = (12, 12.7, 3j, 4, z, 2 * z)
            acc = 0
            for k in prange(10):
                for x in literal_unroll(a):
                    acc += x
            return acc

        f = 9
        k = f

        self.assertEqual(foo(k), foo.py_func(k))

    def test_30(self):
        # function escaping containing an unroll
        @njit
        def foo():
            const = 1234

            def bar(t):
                acc = 0
                a = (12, 12.7, 3j, 4)
                for x in literal_unroll(a):
                    acc += x + const
                return acc, t
            return [x for x in map(bar, (1, 2))]

        self.assertEqual(foo(), foo.py_func())

    def test_31(self):
        # this is testing that generators can survive partial typing
        # invalid function escaping, map uses zip which can't handle the mixed
        # tuple
        @njit
        def foo():
            const = 1234

            def bar(t):
                acc = 0
                a = (12, 12.7, 3j, 4)
                for x in literal_unroll(a):
                    acc += x + const
                return acc, t
            return [x for x in map(bar, (1, 2j))]

        with self.assertRaises(errors.TypingError) as raises:
            foo()

        self.assertIn(_header_lead, str(raises.exception))
        self.assertIn("zip", str(raises.exception))

    def test_32(self):
        # test yielding from an unroll
        @njit
        def gen(a):
            for x in literal_unroll(a):
                yield x

        @njit
        def foo():
            return [x for x in gen((1, 2.3, 4j,))]

        self.assertEqual(foo(), foo.py_func())

    def test_33(self):
        # test yielding from unroll in escaping function that is consumed and
        # yields

        @njit
        def consumer(func, arg):
            yield func(arg)

        def get(cons):
            @njit
            def foo():
                def gen(a):
                    for x in literal_unroll(a):
                        yield x
                return [next(x) for x in cons(gen, (1, 2.3, 4j,))]
            return foo

        cfunc = get(consumer)
        pyfunc = get(consumer.py_func).py_func

        self.assertEqual(cfunc(), pyfunc())

    def test_34(self):
        # mixed bag, redefinition of tuple
        @njit
        def foo():
            acc = 0
            l1 = [1, 2, 3, 4]
            l2 = [10, 20]
            if acc - 2 > 3:
                tup = (l1, l2)
            else:
                a1 = np.arange(20)
                a2 = np.ones(5, dtype=np.complex128)
                tup = (l1, a1, l2, a2)
            for t in literal_unroll(tup):
                acc += len(t)
            return acc

        with self.assertRaises(errors.UnsupportedError) as raises:
            foo()

        self.assertIn("Invalid use of", str(raises.exception))
        self.assertIn("found multiple definitions of variable",
                      str(raises.exception))


class TestConstListUnroll(MemoryLeakMixin, TestCase):

    def test_01(self):

        @njit
        def foo():
            a = [12, 12.7, 3j, 4]
            acc = 0
            for i in range(len(literal_unroll(a))):
                acc += a[i]
                if acc.real < 26:
                    acc -= 1
                else:
                    break
            return acc

        self.assertEqual(foo(), foo.py_func())

    def test_02(self):
        # same as test_1 but without the explicit loop canonicalisation

        @njit
        def foo():
            x = [12, 12.7, 3j, 4]
            acc = 0
            for a in literal_unroll(x):
                acc += a
                if acc.real < 26:
                    acc -= 1
                else:
                    break
            return acc

        self.assertEqual(foo(), foo.py_func())

    def test_03(self):
        # two unrolls
        @njit
        def foo():
            x = [12, 12.7, 3j, 4]
            y = ['foo', 8]
            acc = 0
            for a in literal_unroll(x):
                acc += a
                if acc.real < 26:
                    acc -= 1
                else:
                    for t in literal_unroll(y):
                        acc += t is False
                    break
            return acc

        self.assertEqual(foo(), foo.py_func())

    def test_04(self):
        # two unrolls, one is a const list, one is a tuple
        @njit
        def foo():
            x = [12, 12.7, 3j, 4]
            y = ('foo', 8)
            acc = 0
            for a in literal_unroll(x):
                acc += a
                if acc.real < 26:
                    acc -= 1
                else:
                    for t in literal_unroll(y):
                        acc += t is False
                    break
            return acc

        self.assertEqual(foo(), foo.py_func())

    def test_05(self):
        # illegal, list has to be const
        @njit
        def foo(tup1, tup2):
            acc = 0
            for a in literal_unroll(tup1):
                if a[0] > 1:
                    acc += tup2[0].sum()
            return acc

        n = 10
        tup1 = [np.zeros(10), np.zeros(10)]
        tup2 = (np.ones((n,)), np.ones((n, n)), np.ones((n, n, n)),
                np.ones((n, n, n, n)), np.ones((n, n, n, n, n)))

        with self.assertRaises(errors.UnsupportedError) as raises:
            foo(tup1, tup2)

        msg = "Invalid use of literal_unroll with a function argument"
        self.assertIn(msg, str(raises.exception))

    def test_06(self):
        # illegal: list containing non const
        @njit
        def foo():
            n = 10
            tup = [np.ones((n,)), np.ones((n, n)), "ABCDEFGHJI", (1, 2, 3),
                   (1, 'foo', 2, 'bar'), {3, 4, 5, 6, 7}]
            acc = 0
            for a in literal_unroll(tup):
                acc += len(a)
            return acc

        with self.assertRaises(errors.UnsupportedError) as raises:
            foo()

        self.assertIn("Found non-constant value at position 0",
                      str(raises.exception))

    def test_7(self):
        # dispatch on literals triggering @overload resolution

        def dt(value):
            if value == "apple":
                return 1
            elif value == "orange":
                return 2
            elif value == "banana":
                return 3
            elif value == 0xca11ab1e:
                return 0x5ca1ab1e + value

        @overload(dt, inline='always')
        def ol_dt(li):
            if isinstance(li, types.StringLiteral):
                value = li.literal_value
                if value == "apple":
                    def impl(li):
                        return 1
                elif value == "orange":
                    def impl(li):
                        return 2
                elif value == "banana":
                    def impl(li):
                        return 3
                return impl
            elif isinstance(li, types.IntegerLiteral):
                value = li.literal_value
                if value == 0xca11ab1e:
                    def impl(li):
                        # close over the dispatcher :)
                        return 0x5ca1ab1e + value
                    return impl

        @njit
        def foo():
            acc = 0
            for t in literal_unroll(['apple', 'orange', 'banana', 3390155550]):
                acc += dt(t)
            return acc

        self.assertEqual(foo(), foo.py_func())

    def test_8(self):

        @njit
        def foo():
            x = []
            z = ['apple', 'orange', 'banana']
            for i in range(len(literal_unroll(z))):
                t = z[i]
                if t == "apple":
                    x.append("0")
                elif t == "orange":
                    x.append(t)
                elif t == "banana":
                    x.append("2.0")
            return x

        self.assertEqual(foo(), foo.py_func())

    def test_9(self):
        # unroll the same target twice
        @njit
        def foo(idx, z):
            a = [12, 12.7, 3j, 4]
            acc = 0
            for i in literal_unroll(a):
                acc += i
                if acc.real < 26:
                    acc -= 1
                else:
                    for x in literal_unroll(a):
                        acc += x
                    break
            if a[0] < 23:
                acc += 2
            return acc

        f = 9
        k = f

        self.assertEqual(foo(2, k), foo.py_func(2, k))

    def test_10(self):
        # nesting unrolls is illegal
        @njit
        def foo(idx, z):
            a = (12, 12.7, 3j, 4, z, 2 * z)
            b = [12, 12.7, 3j, 4]
            acc = 0
            for i in literal_unroll(a):
                acc += i
                if acc.real < 26:
                    acc -= 1
                else:
                    for x in literal_unroll(a):
                        for j in literal_unroll(b):
                            acc += j
                        acc += x
                for x in literal_unroll(a):
                    acc += x
            for x in literal_unroll(a):
                acc += x
            if a[0] < 23:
                acc += 2
            return acc

        f = 9
        k = f

        with self.assertRaises(errors.UnsupportedError) as raises:
            foo(2, k)

        self.assertIn("Nesting of literal_unroll is unsupported",
                      str(raises.exception))

    def test_11(self):
        # homogeneous const list unroll can return derivative of the induction
        # var

        @njit
        def foo():
            x = [1, 2, 3, 4]
            acc = 0
            for a in literal_unroll(x):
                acc += a
            return a

        self.assertEqual(foo(), foo.py_func())

    def test_12(self):
        # mixed unroll cannot return derivative of the induction var
        @njit
        def foo():
            acc = 0
            x = [1, 2, 'a']
            for a in literal_unroll(x):
                acc += bool(a)
            return a

        with self.assertRaises(errors.TypingError) as raises:
            foo()

        self.assertIn("Cannot unify", str(raises.exception))

    def test_13(self):
        # list slice is illegal

        @njit
        def foo():
            x = [1000, 2000, 3000, 4000]
            acc = 0
            for a in literal_unroll(x[:2]):
                acc += a
            return acc

        with self.assertRaises(errors.UnsupportedError) as raises:
            foo()

        self.assertIn("Invalid use of literal_unroll", str(raises.exception))

    def test_14(self):
        # list mutate is illegal

        @njit
        def foo():
            x = [1000, 2000, 3000, 4000]
            acc = 0
            for a in literal_unroll(x):
                acc += a
            x.append(10)
            return acc

        with self.assertRaises(errors.TypingError) as raises:
            foo()

        self.assertIn("Unknown attribute 'append' of type Tuple",
                      str(raises.exception))


class TestMore(TestCase):
    def test_invalid_use_of_unroller(self):
        @njit
        def foo():
            x = (10, 20)
            r = 0
            for a in literal_unroll(x, x):
                r += a
            return r

        with self.assertRaises(errors.UnsupportedError) as raises:
            foo()
        self.assertIn(
            "literal_unroll takes one argument, found 2",
            str(raises.exception),
        )

    def test_non_constant_list(self):

        @njit
        def foo(y):
            x = [10, y]
            r = 0
            for a in literal_unroll(x):
                r += a
            return r

        with self.assertRaises(errors.UnsupportedError) as raises:
            foo(10)
        self.assertIn(
            ("Found non-constant value at position 1 in a list argument to "
             "literal_unroll"),
            str(raises.exception)
        )

    @unittest.skip("numba.literally not supported yet")
    def test_literally_constant_list(self):
        # FAIL. May need to consider it in a future PR
        from numba import literally

        @njit
        def foo(y):
            x = [10, literally(y)]
            r = 0
            for a in literal_unroll(x):
                r += a
            return r

        # Found non-constant value at position 1 in a list argument to
        # literal_unroll
        foo(12)

        @njit
        def bar():
            return foo(12)

        # Found non-constant value at position 1 in a list argument to
        # literal_unroll
        bar()

    @unittest.skip("inlining of foo doesn't have const prop so y isn't const")
    def test_inlined_unroll_list(self):
        @njit(inline='always')
        def foo(y):
            x = [10, y]
            r = 0
            for a in literal_unroll(x):
                r += a
            return r

        @njit
        def bar():
            return foo(12)

        self.assertEqual(bar(), 10 + 12)

    def test_unroll_tuple_arg(self):
        @njit
        def foo(y):
            x = (10, y)
            r = 0
            for a in literal_unroll(x):
                r += a
            return r

        self.assertEqual(foo(12), foo.py_func(12))
        self.assertEqual(foo(1.2), foo.py_func(1.2))

    def test_unroll_tuple_arg2(self):
        @njit
        def foo(x):
            r = 0
            for a in literal_unroll(x):
                r += a
            return r

        self.assertEqual(foo((12, 1.2)), foo.py_func((12, 1.2)))
        self.assertEqual(foo((12, 1.2)), foo.py_func((12, 1.2)))

    def test_unroll_tuple_alias(self):
        @njit
        def foo():
            x = (10, 1.2)
            out = 0
            for i in literal_unroll(x):
                j = i
                k = j
                out += j + k + i
            return out

        self.assertEqual(foo(), foo.py_func())

    def test_unroll_tuple_nested(self):

        @njit
        def foo():
            x = ((10, 1.2), (1j, 3.))
            out = 0
            for i in literal_unroll(x):
                for j in (i):
                    out += j
            return out

        with self.assertRaises(errors.TypingError) as raises:
            foo()

        self.assertIn("getiter", str(raises.exception))
        re = r".*Tuple\(int[0-9][0-9], float64\).*"
        self.assertRegex(str(raises.exception), re)

    def test_unroll_tuple_of_dict(self):

        @njit
        def foo():
            x = {}
            x["a"] = 1
            x["b"] = 2
            y = {}
            y[3] = "c"
            y[4] = "d"

            for it in literal_unroll((x, y)):
                for k, v in it.items():
                    print(k, v)

        with captured_stdout() as stdout:
            foo()
        lines = stdout.getvalue().splitlines()
        self.assertEqual(
            lines,
            ['a 1', 'b 2', '3 c', '4 d'],
        )

    def test_unroll_named_tuple(self):
        ABC = namedtuple('ABC', ['a', 'b', 'c'])

        @njit
        def foo():
            abc = ABC(1, 2j, 3.4)
            out = 0
            for i in literal_unroll(abc):
                out += i
            return out

        self.assertEqual(foo(), foo.py_func())

    def test_unroll_named_tuple_arg(self):
        ABC = namedtuple('ABC', ['a', 'b', 'c'])

        @njit
        def foo(x):
            out = 0
            for i in literal_unroll(x):
                out += i
            return out

        abc = ABC(1, 2j, 3.4)

        self.assertEqual(foo(abc), foo.py_func(abc))

    def test_unroll_named_unituple(self):
        ABC = namedtuple('ABC', ['a', 'b', 'c'])

        @njit
        def foo():
            abc = ABC(1, 2, 3)
            out = 0
            for i in literal_unroll(abc):
                out += i
            return out

        self.assertEqual(foo(), foo.py_func())

    def test_unroll_named_unituple_arg(self):
        ABC = namedtuple('ABC', ['a', 'b', 'c'])

        @njit
        def foo(x):
            out = 0
            for i in literal_unroll(x):
                out += i
            return out

        abc = ABC(1, 2, 3)

        self.assertEqual(foo(abc), foo.py_func(abc))

    def test_unroll_global_tuple(self):

        @njit
        def foo():
            out = 0
            for i in literal_unroll(_X_GLOBAL):
                out += i
            return out

        self.assertEqual(foo(), foo.py_func())

    def test_unroll_freevar_tuple(self):
        x = (10, 11)

        @njit
        def foo():
            out = 0
            for i in literal_unroll(x):
                out += i
            return out

        self.assertEqual(foo(), foo.py_func())

    def test_unroll_function_tuple(self):
        @njit
        def a():
            return 1

        @njit
        def b():
            return 2

        x = (a, b)

        @njit
        def foo():
            out = 0
            for f in literal_unroll(x):
                out += f()
            return out

        self.assertEqual(foo(), foo.py_func())

    def test_unroll_indexing_list(self):
        # See issue #5477
        @njit
        def foo(cont):
            i = 0
            acc = 0
            normal_list = [a for a in cont]
            heter_tuple = ('a', 25, 0.23, None)
            for item in literal_unroll(heter_tuple):
                acc += normal_list[i]
                i += 1
                print(item)
            return i, acc

        data = [j for j in range(4)]

        # send stdout to nowhere, just check return values
        with captured_stdout():
            self.assertEqual(foo(data), foo.py_func(data))

        # now capture stdout for jit function and check
        with captured_stdout() as stdout:
            foo(data)
        lines = stdout.getvalue().splitlines()
        self.assertEqual(
            lines,
            ['a', '25', '0.23', 'None'],
        )

    def test_unroller_as_freevar(self):
        mixed = (np.ones((1,)), np.ones((1, 1)), np.ones((1, 1, 1)))
        from numba import literal_unroll as freevar_unroll

        @njit
        def foo():
            out = 0
            for i in freevar_unroll(mixed):
                out += i.ndim
            return out

        self.assertEqual(foo(), foo.py_func())

    def test_unroll_with_non_conformant_loops_present(self):
        # See issue #8311

        @njit('(Tuple((int64, float64)),)')
        def foo(tup):
            for t in literal_unroll(tup):
                pass

            x = 1
            while x == 1:
                x = 0

    def test_literal_unroll_legalize_var_names01(self):
        # See issue #8939
        test = np.array([(1, 2), (2, 3)], dtype=[("a1", "f8"), ("a2", "f8")])
        fields = tuple(test.dtype.fields.keys())

        @njit
        def foo(arr):
            res = 0
            for k in literal_unroll(fields):
                res = res + np.abs(arr[k]).sum()
            return res

        self.assertEqual(foo(test), 8.0)

    def test_literal_unroll_legalize_var_names02(self):
        # See issue #8939
        test = np.array([(1, 2), (2, 3)],
                        dtype=[("a1[0]", "f8"), ("a2[1]", "f8")])
        fields = tuple(test.dtype.fields.keys())

        @njit
        def foo(arr):
            res = 0
            for k in literal_unroll(fields):
                res = res + np.abs(arr[k]).sum()
            return res

        self.assertEqual(foo(test), 8.0)


def capture(real_pass):
    """ Returns a compiler pass that captures the mutation state reported
    by the pass used in the argument"""
    @register_pass(mutates_CFG=False, analysis_only=True)
    class ResultCapturer(AnalysisPass):
        _name = "capture_%s" % real_pass._name
        _real_pass = real_pass

        def __init__(self):
            FunctionPass.__init__(self)

        def run_pass(self, state):
            result = real_pass().run_pass(state)
            mutation_results = state.metadata.setdefault('mutation_results', {})
            mutation_results[real_pass] = result
            return result

    return ResultCapturer


class CapturingCompiler(CompilerBase):
    """ Simple pipeline that wraps passes with the ResultCapturer pass"""

    def define_pipelines(self):
        pm = PassManager("Capturing Compiler")

        def add_pass(x, y):
            return pm.add_pass(capture(x), y)

        add_pass(TranslateByteCode, "analyzing bytecode")
        add_pass(FixupArgs, "fix up args")
        add_pass(IRProcessing, "processing IR")
        add_pass(LiteralUnroll, "handles literal_unroll")

        # typing
        add_pass(NopythonTypeInference, "nopython frontend")

        # legalise
        add_pass(IRLegalization,
                 "ensure IR is legal prior to lowering")

        # lower
        add_pass(NativeLowering, "native lowering")
        add_pass(NoPythonBackend, "nopython mode backend")
        pm.finalize()
        return [pm]


class TestLiteralUnrollPassTriggering(TestCase):

    def test_literal_unroll_not_invoked(self):
        @njit(pipeline_class=CapturingCompiler)
        def foo():
            acc = 0
            for i in (1, 2, 3):
                acc += i
            return acc

        foo()
        cres = foo.overloads[foo.signatures[0]]
        self.assertFalse(cres.metadata['mutation_results'][LiteralUnroll])

    def test_literal_unroll_is_invoked(self):
        @njit(pipeline_class=CapturingCompiler)
        def foo():
            acc = 0
            for i in literal_unroll((1, 2, 3)):
                acc += i
            return acc

        foo()
        cres = foo.overloads[foo.signatures[0]]
        self.assertTrue(cres.metadata['mutation_results'][LiteralUnroll])

    def test_literal_unroll_is_invoked_via_alias(self):
        alias = literal_unroll

        @njit(pipeline_class=CapturingCompiler)
        def foo():
            acc = 0
            for i in alias((1, 2, 3)):
                acc += i
            return acc

        foo()
        cres = foo.overloads[foo.signatures[0]]
        self.assertTrue(cres.metadata['mutation_results'][LiteralUnroll])

    def test_literal_unroll_assess_empty_function(self):
        @njit(pipeline_class=CapturingCompiler)
        def foo():
            pass

        foo()
        cres = foo.overloads[foo.signatures[0]]
        self.assertFalse(cres.metadata['mutation_results'][LiteralUnroll])

    def test_literal_unroll_not_in_globals(self):
        f = """def foo():\n\tpass"""
        l = {}
        exec(f, {}, l)
        foo = njit(pipeline_class=CapturingCompiler)(l['foo'])
        foo()
        cres = foo.overloads[foo.signatures[0]]
        self.assertFalse(cres.metadata['mutation_results'][LiteralUnroll])

    def test_literal_unroll_globals_and_locals(self):
        f = """def foo():\n\tfor x in literal_unroll((1,)):\n\t\tpass"""
        l = {}
        exec(f, {}, l)
        foo = njit(pipeline_class=CapturingCompiler)(l['foo'])
        with self.assertRaises(errors.TypingError) as raises:
            foo()
        self.assertIn("Untyped global name 'literal_unroll'",
                      str(raises.exception))

        # same as above but now add literal_unroll to globals
        l = {}
        exec(f, {'literal_unroll': literal_unroll}, l)
        foo = njit(pipeline_class=CapturingCompiler)(l['foo'])
        foo()
        cres = foo.overloads[foo.signatures[0]]
        self.assertTrue(cres.metadata['mutation_results'][LiteralUnroll])

        # same as above, but now with import
        from textwrap import dedent
        f = """
            def gen():
                from numba import literal_unroll
                def foo():
                    for x in literal_unroll((1,)):
                        pass
                return foo
            bar = gen()
            """
        l = {}
        exec(dedent(f), {}, l)
        foo = njit(pipeline_class=CapturingCompiler)(l['bar'])
        foo()
        cres = foo.overloads[foo.signatures[0]]
        self.assertTrue(cres.metadata['mutation_results'][LiteralUnroll])

        # same as above, but now with import as something else
        from textwrap import dedent
        f = """
            def gen():
                from numba import literal_unroll as something_else
                def foo():
                    for x in something_else((1,)):
                        pass
                return foo
            bar = gen()
            """
        l = {}
        exec(dedent(f), {}, l)
        foo = njit(pipeline_class=CapturingCompiler)(l['bar'])
        foo()
        cres = foo.overloads[foo.signatures[0]]
        self.assertTrue(cres.metadata['mutation_results'][LiteralUnroll])


if __name__ == '__main__':
    unittest.main()
