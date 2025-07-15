import unittest
from llvmlite import ir
from llvmlite import binding as llvm
from llvmlite.tests import TestCase

import llvmlite.tests.refprune_proto as proto

# TODO:: Get rid of Legacy tests once completely transitioned to NewPassManager

# FIXME: Remove me once typed pointers are no longer supported.
from llvmlite import opaque_pointers_enabled


def _iterate_cases(generate_test):
    def wrap(fn):
        def wrapped(self):
            return generate_test(self, fn)
        wrapped.__doc__ = f"generated test for {fn.__module__}.{fn.__name__}"
        return wrapped

    for k, case_fn in proto.__dict__.items():
        if k.startswith('case'):
            yield f'test_{k}', wrap(case_fn)


class PassManagerMixin():

    def pb(self):
        llvm.initialize_native_target()
        tm = llvm.Target.from_default_triple().create_target_machine()
        pto = llvm.create_pipeline_tuning_options(speed_level=0, size_level=0)
        return llvm.create_pass_builder(tm, pto)


class TestRefPrunePrototype(TestCase):
    """
    Test that the prototype is working.
    """
    def generate_test(self, case_gen):
        nodes, edges, expected = case_gen()
        got = proto.FanoutAlgorithm(nodes, edges).run()
        self.assertEqual(expected, got)

    # Generate tests
    for name, case in _iterate_cases(generate_test):
        locals()[name] = case


ptr_ty = ir.IntType(8).as_pointer()


class TestRefPrunePass(TestCase, PassManagerMixin):
    """
    Test that the C++ implementation matches the expected behavior as for
    the prototype.

    This generates a LLVM module for each test case, runs the pruner and checks
    that the expected results are achieved.
    """

    def make_incref(self, m):
        fnty = ir.FunctionType(ir.VoidType(), [ptr_ty])
        return ir.Function(m, fnty, name='NRT_incref')

    def make_decref(self, m):
        fnty = ir.FunctionType(ir.VoidType(), [ptr_ty])
        return ir.Function(m, fnty, name='NRT_decref')

    def make_switcher(self, m):
        fnty = ir.FunctionType(ir.IntType(32), ())
        return ir.Function(m, fnty, name='switcher')

    def make_brancher(self, m):
        fnty = ir.FunctionType(ir.IntType(1), ())
        return ir.Function(m, fnty, name='brancher')

    def generate_ir(self, nodes, edges):
        # Build LLVM module for the CFG
        m = ir.Module()

        incref_fn = self.make_incref(m)
        decref_fn = self.make_decref(m)
        switcher_fn = self.make_switcher(m)
        brancher_fn = self.make_brancher(m)

        fnty = ir.FunctionType(ir.VoidType(), [ptr_ty])
        fn = ir.Function(m, fnty, name='main')
        [ptr] = fn.args
        ptr.name = 'mem'
        # populate the BB nodes
        bbmap = {}
        for bb in edges:
            bbmap[bb] = fn.append_basic_block(bb)
        # populate the BB
        builder = ir.IRBuilder()
        for bb, jump_targets in edges.items():
            builder.position_at_end(bbmap[bb])
            # Insert increfs and decrefs
            for action in nodes[bb]:
                if action == 'incref':
                    builder.call(incref_fn, [ptr])
                elif action == 'decref':
                    builder.call(decref_fn, [ptr])
                else:
                    raise AssertionError('unreachable')

            # Insert the terminator.
            # Switch base on the number of jump targets.
            n_targets = len(jump_targets)
            if n_targets == 0:
                builder.ret_void()
            elif n_targets == 1:
                [dst] = jump_targets
                builder.branch(bbmap[dst])
            elif n_targets == 2:
                [left, right] = jump_targets
                sel = builder.call(brancher_fn, ())
                builder.cbranch(sel, bbmap[left], bbmap[right])
            elif n_targets > 2:
                sel = builder.call(switcher_fn, ())
                [head, *tail] = jump_targets

                sw = builder.switch(sel, default=bbmap[head])
                for i, dst in enumerate(tail):
                    sw.add_case(sel.type(i), bbmap[dst])
            else:
                raise AssertionError('unreachable')

        return m

    def apply_refprune(self, irmod):
        mod = llvm.parse_assembly(str(irmod))
        pb = self.pb()
        pm = pb.getModulePassManager()
        pm.add_refprune_pass()
        pm.run(mod, pb)
        return mod

    def apply_refprune_legacy(self, irmod):
        mod = llvm.parse_assembly(str(irmod))
        pm = llvm.ModulePassManager()
        pm.add_refprune_pass()
        pm.run(mod)
        return mod

    def check(self, mod, expected, nodes):
        # preprocess incref/decref locations
        d = {}
        for k, vs in nodes.items():
            n_incref = vs.count('incref')
            n_decref = vs.count('decref')
            d[k] = {'incref': n_incref, 'decref': n_decref}
        for k, stats in d.items():
            if expected.get(k):
                stats['incref'] -= 1
                for dec_bb in expected[k]:
                    d[dec_bb]['decref'] -= 1

        # find the main function
        for f in mod.functions:
            if f.name == 'main':
                break
        # check each BB
        for bb in f.blocks:
            stats = d[bb.name]
            text = str(bb)
            n_incref = text.count('NRT_incref')
            n_decref = text.count('NRT_decref')
            self.assertEqual(stats['incref'], n_incref, msg=f'BB {bb}')
            self.assertEqual(stats['decref'], n_decref, msg=f'BB {bb}')

    def generate_test(self, case_gen):
        nodes, edges, expected = case_gen()
        irmod = self.generate_ir(nodes, edges)
        outmod = self.apply_refprune(irmod)
        self.check(outmod, expected, nodes)

    def generate_test_legacy(self, case_gen):
        nodes, edges, expected = case_gen()
        irmod = self.generate_ir(nodes, edges)
        outmod = self.apply_refprune_legacy(irmod)
        self.check(outmod, expected, nodes)

    # Generate tests
    for name, case in _iterate_cases(generate_test):
        locals()[name] = case

    for name, case in _iterate_cases(generate_test_legacy):
        locals()[name + "_legacy"] = case


class BaseTestByIR(TestCase, PassManagerMixin):
    refprune_bitmask = 0

    prologue = r"""
declare void @NRT_incref(i8* %ptr)
declare void @NRT_decref(i8* %ptr)
"""

    def check(self, irmod, subgraph_limit=None):
        mod = llvm.parse_assembly(f"{self.prologue}\n{irmod}")
        pb = self.pb()
        pm = pb.getModulePassManager()
        if subgraph_limit is None:
            pm.add_refprune_pass(self.refprune_bitmask)
        else:
            pm.add_refprune_pass(self.refprune_bitmask,
                                 subgraph_limit=subgraph_limit)
        before = llvm.dump_refprune_stats()
        pm.run(mod, pb)
        after = llvm.dump_refprune_stats()
        return mod, after - before

    def check_legacy(self, irmod, subgraph_limit=None):
        mod = llvm.parse_assembly(f"{self.prologue}\n{irmod}")
        pm = llvm.ModulePassManager()
        if subgraph_limit is None:
            pm.add_refprune_pass(self.refprune_bitmask)
        else:
            pm.add_refprune_pass(self.refprune_bitmask,
                                 subgraph_limit=subgraph_limit)
        before = llvm.dump_refprune_stats()
        pm.run(mod)
        after = llvm.dump_refprune_stats()
        return mod, after - before


class TestPerBB(BaseTestByIR):
    refprune_bitmask = llvm.RefPruneSubpasses.PER_BB

    per_bb_ir_1 = r"""
define void @main(i8* %ptr) {
    call void @NRT_incref(i8* %ptr)
    call void @NRT_decref(i8* %ptr)
    ret void
}
"""

    def test_per_bb_1(self):
        mod, stats = self.check(self.per_bb_ir_1)
        self.assertEqual(stats.basicblock, 2)

    def test_per_bb_1_legacy(self):
        mod, stats = self.check_legacy(self.per_bb_ir_1)
        self.assertEqual(stats.basicblock, 2)

    per_bb_ir_2 = r"""
define void @main(i8* %ptr) {
    call void @NRT_incref(i8* %ptr)
    call void @NRT_incref(i8* %ptr)
    call void @NRT_incref(i8* %ptr)
    call void @NRT_decref(i8* %ptr)
    call void @NRT_decref(i8* %ptr)
    ret void
}
"""

    def test_per_bb_2(self):
        mod, stats = self.check(self.per_bb_ir_2)
        self.assertEqual(stats.basicblock, 4)
        # not pruned
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.assertIn("call void @NRT_incref(ptr %ptr)", str(mod))
        else:
            self.assertIn("call void @NRT_incref(i8* %ptr)", str(mod))

    def test_per_bb_2_legacy(self):
        mod, stats = self.check_legacy(self.per_bb_ir_2)
        self.assertEqual(stats.basicblock, 4)
        # not pruned
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.assertIn("call void @NRT_incref(ptr %ptr)", str(mod))
        else:
            self.assertIn("call void @NRT_incref(i8* %ptr)", str(mod))

    per_bb_ir_3 = r"""
define void @main(ptr %ptr, ptr %other) {
    call void @NRT_incref(ptr %ptr)
    call void @NRT_incref(ptr %ptr)
    call void @NRT_decref(ptr %ptr)
    call void @NRT_decref(ptr %other)
    ret void
}
""" if opaque_pointers_enabled else r"""
define void @main(i8* %ptr, i8* %other) {
    call void @NRT_incref(i8* %ptr)
    call void @NRT_incref(i8* %ptr)
    call void @NRT_decref(i8* %ptr)
    call void @NRT_decref(i8* %other)
    ret void
}
"""

    def test_per_bb_3(self):
        mod, stats = self.check(self.per_bb_ir_3)
        self.assertEqual(stats.basicblock, 2)
        # not pruned
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.assertIn("call void @NRT_decref(ptr %other)", str(mod))
        else:
            self.assertIn("call void @NRT_decref(i8* %other)", str(mod))

    def test_per_bb_3_legacy(self):
        mod, stats = self.check_legacy(self.per_bb_ir_3)
        self.assertEqual(stats.basicblock, 2)
        # not pruned
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.assertIn("call void @NRT_decref(ptr %other)", str(mod))
        else:
            self.assertIn("call void @NRT_decref(i8* %other)", str(mod))

    per_bb_ir_4 = r"""
; reordered
define void @main(ptr %ptr, ptr %other) {
    call void @NRT_incref(ptr %ptr)
    call void @NRT_decref(ptr %ptr)
    call void @NRT_decref(ptr %ptr)
    call void @NRT_decref(ptr %other)
    call void @NRT_incref(ptr %ptr)
    ret void
}
""" if opaque_pointers_enabled else r"""
; reordered
define void @main(i8* %ptr, i8* %other) {
    call void @NRT_incref(i8* %ptr)
    call void @NRT_decref(i8* %ptr)
    call void @NRT_decref(i8* %ptr)
    call void @NRT_decref(i8* %other)
    call void @NRT_incref(i8* %ptr)
    ret void
}
"""

    def test_per_bb_4(self):
        mod, stats = self.check(self.per_bb_ir_4)
        self.assertEqual(stats.basicblock, 4)
        # not pruned
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.assertIn("call void @NRT_decref(ptr %other)", str(mod))
        else:
            self.assertIn("call void @NRT_decref(i8* %other)", str(mod))

    def test_per_bb_4_legacy(self):
        mod, stats = self.check_legacy(self.per_bb_ir_4)
        self.assertEqual(stats.basicblock, 4)
        # not pruned
        # FIXME: Remove `else' once TP are no longer supported.
        if opaque_pointers_enabled:
            self.assertIn("call void @NRT_decref(ptr %other)", str(mod))
        else:
            self.assertIn("call void @NRT_decref(i8* %other)", str(mod))


class TestDiamond(BaseTestByIR):
    refprune_bitmask = llvm.RefPruneSubpasses.DIAMOND

    per_diamond_1 = r"""
define void @main(i8* %ptr) {
bb_A:
    call void @NRT_incref(i8* %ptr)
    br label %bb_B
bb_B:
    call void @NRT_decref(i8* %ptr)
    ret void
}
"""

    def test_per_diamond_1(self):
        mod, stats = self.check(self.per_diamond_1)
        self.assertEqual(stats.diamond, 2)

    def test_per_diamond_1_legacy(self):
        mod, stats = self.check_legacy(self.per_diamond_1)
        self.assertEqual(stats.diamond, 2)

    per_diamond_2 = r"""
define void @main(i8* %ptr, i1 %cond) {
bb_A:
    call void @NRT_incref(i8* %ptr)
    br i1 %cond, label %bb_B, label %bb_C
bb_B:
    br label %bb_D
bb_C:
    br label %bb_D
bb_D:
    call void @NRT_decref(i8* %ptr)
    ret void
}
"""

    def test_per_diamond_2(self):
        mod, stats = self.check(self.per_diamond_2)
        self.assertEqual(stats.diamond, 2)

    def test_per_diamond_2_legacy(self):
        mod, stats = self.check_legacy(self.per_diamond_2)
        self.assertEqual(stats.diamond, 2)

    per_diamond_3 = r"""
define void @main(i8* %ptr, i1 %cond) {
bb_A:
    call void @NRT_incref(i8* %ptr)
    br i1 %cond, label %bb_B, label %bb_C
bb_B:
    br label %bb_D
bb_C:
    call void @NRT_decref(i8* %ptr)  ; reject because of decref in diamond
    br label %bb_D
bb_D:
    call void @NRT_decref(i8* %ptr)
    ret void
}
"""

    def test_per_diamond_3(self):
        mod, stats = self.check(self.per_diamond_3)
        self.assertEqual(stats.diamond, 0)

    def test_per_diamond_3_legacy(self):
        mod, stats = self.check_legacy(self.per_diamond_3)
        self.assertEqual(stats.diamond, 0)

    per_diamond_4 = r"""
define void @main(i8* %ptr, i1 %cond) {
bb_A:
    call void @NRT_incref(i8* %ptr)
    br i1 %cond, label %bb_B, label %bb_C
bb_B:
    call void @NRT_incref(i8* %ptr)     ; extra incref will not affect prune
    br label %bb_D
bb_C:
    br label %bb_D
bb_D:
    call void @NRT_decref(i8* %ptr)
    ret void
}
"""

    def test_per_diamond_4(self):
        mod, stats = self.check(self.per_diamond_4)
        self.assertEqual(stats.diamond, 2)

    def test_per_diamond_4_legacy(self):
        mod, stats = self.check_legacy(self.per_diamond_4)
        self.assertEqual(stats.diamond, 2)

    per_diamond_5 = r"""
define void @main(i8* %ptr, i1 %cond) {
bb_A:
    call void @NRT_incref(i8* %ptr)
    call void @NRT_incref(i8* %ptr)
    br i1 %cond, label %bb_B, label %bb_C
bb_B:
    br label %bb_D
bb_C:
    br label %bb_D
bb_D:
    call void @NRT_decref(i8* %ptr)
    call void @NRT_decref(i8* %ptr)
    ret void
}
"""

    def test_per_diamond_5(self):
        mod, stats = self.check(self.per_diamond_5)
        self.assertEqual(stats.diamond, 4)

    def test_per_diamond_5_legacy(self):
        mod, stats = self.check_legacy(self.per_diamond_5)
        self.assertEqual(stats.diamond, 4)


class TestFanout(BaseTestByIR):
    """More complex cases are tested in TestRefPrunePass
    """

    refprune_bitmask = llvm.RefPruneSubpasses.FANOUT

    fanout_1 = r"""
define void @main(i8* %ptr, i1 %cond) {
bb_A:
    call void @NRT_incref(i8* %ptr)
    br i1 %cond, label %bb_B, label %bb_C
bb_B:
    call void @NRT_decref(i8* %ptr)
    ret void
bb_C:
    call void @NRT_decref(i8* %ptr)
    ret void
}
"""

    def test_fanout_1(self):
        mod, stats = self.check(self.fanout_1)
        self.assertEqual(stats.fanout, 3)

    def test_fanout_1_legacy(self):
        mod, stats = self.check_legacy(self.fanout_1)
        self.assertEqual(stats.fanout, 3)

    fanout_2 = r"""
define void @main(i8* %ptr, i1 %cond, i8** %excinfo) {
bb_A:
    call void @NRT_incref(i8* %ptr)
    br i1 %cond, label %bb_B, label %bb_C
bb_B:
    call void @NRT_decref(i8* %ptr)
    ret void
bb_C:
    call void @NRT_decref(i8* %ptr)
    br label %bb_B                      ; illegal jump to other decref
}
"""

    def test_fanout_2(self):
        mod, stats = self.check(self.fanout_2)
        self.assertEqual(stats.fanout, 0)

    def test_fanout_2_legacy(self):
        mod, stats = self.check_legacy(self.fanout_2)
        self.assertEqual(stats.fanout, 0)

    fanout_3 = r"""
define void @main(i8* %ptr, i1 %cond) {
bb_A:
    call void @NRT_incref(i8* %ptr)
    call void @NRT_incref(i8* %ptr)
    br i1 %cond, label %bb_B, label %bb_C
bb_B:
    call void @NRT_decref(i8* %ptr)
    call void @NRT_decref(i8* %ptr)
    call void @NRT_decref(i8* %ptr)
    ret void
bb_C:
    call void @NRT_decref(i8* %ptr)
    call void @NRT_decref(i8* %ptr)
    ret void
}
"""

    def test_fanout_3(self):
        mod, stats = self.check(self.fanout_3)
        self.assertEqual(stats.fanout, 6)

    def test_fanout_3_limited(self):
        # With subgraph limit at 1, it is essentially turning off the fanout
        # pruner.
        mod, stats = self.check(self.fanout_3, subgraph_limit=1)
        self.assertEqual(stats.fanout, 0)

    def test_fanout_3_legacy(self):
        mod, stats = self.check_legacy(self.fanout_3)
        self.assertEqual(stats.fanout, 6)

    def test_fanout_3_limited_legacy(self):
        # With subgraph limit at 1, it is essentially turning off the fanout
        # pruner.
        mod, stats = self.check_legacy(self.fanout_3, subgraph_limit=1)
        self.assertEqual(stats.fanout, 0)


class TestFanoutRaise(BaseTestByIR):
    refprune_bitmask = llvm.RefPruneSubpasses.FANOUT_RAISE

    fanout_raise_1 = r"""
define i32 @main(i8* %ptr, i1 %cond, i8** %excinfo) {
bb_A:
    call void @NRT_incref(i8* %ptr)
    br i1 %cond, label %bb_B, label %bb_C
bb_B:
    call void @NRT_decref(i8* %ptr)
    ret i32 0
bb_C:
    store i8* null, i8** %excinfo, !numba_exception_output !0
    ret i32 1
}
!0 = !{i1 true}
"""

    def test_fanout_raise_1(self):
        mod, stats = self.check(self.fanout_raise_1)
        self.assertEqual(stats.fanout_raise, 2)

    def test_fanout_raise_1_legacy(self):
        mod, stats = self.check_legacy(self.fanout_raise_1)
        self.assertEqual(stats.fanout_raise, 2)

    fanout_raise_2 = r"""
define i32 @main(i8* %ptr, i1 %cond, i8** %excinfo) {
bb_A:
    call void @NRT_incref(i8* %ptr)
    br i1 %cond, label %bb_B, label %bb_C
bb_B:
    call void @NRT_decref(i8* %ptr)
    ret i32 0
bb_C:
    store i8* null, i8** %excinfo, !numba_exception_typo !0      ; bad metadata
    ret i32 1
}

!0 = !{i1 true}
"""

    def test_fanout_raise_2(self):
        # This is ensuring that fanout_raise is not pruning when the metadata
        # is incorrectly named.
        mod, stats = self.check(self.fanout_raise_2)
        self.assertEqual(stats.fanout_raise, 0)

    def test_fanout_raise_2_legacy(self):
        # This is ensuring that fanout_raise is not pruning when the metadata
        # is incorrectly named.
        mod, stats = self.check_legacy(self.fanout_raise_2)
        self.assertEqual(stats.fanout_raise, 0)

    fanout_raise_3 = r"""
define i32 @main(i8* %ptr, i1 %cond, i8** %excinfo) {
bb_A:
    call void @NRT_incref(i8* %ptr)
    br i1 %cond, label %bb_B, label %bb_C
bb_B:
    call void @NRT_decref(i8* %ptr)
    ret i32 0
bb_C:
    store i8* null, i8** %excinfo, !numba_exception_output !0
    ret i32 1
}

!0 = !{i32 1}       ; ok; use i32
"""

    def test_fanout_raise_3(self):
        mod, stats = self.check(self.fanout_raise_3)
        self.assertEqual(stats.fanout_raise, 2)

    def test_fanout_raise_3_legacy(self):
        mod, stats = self.check_legacy(self.fanout_raise_3)
        self.assertEqual(stats.fanout_raise, 2)

    fanout_raise_4 = r"""
define i32 @main(i8* %ptr, i1 %cond, i8** %excinfo) {
bb_A:
    call void @NRT_incref(i8* %ptr)
    br i1 %cond, label %bb_B, label %bb_C
bb_B:
    ret i32 1    ; BAD; all tails are raising without decref
bb_C:
    ret i32 1    ; BAD; all tails are raising without decref
}

!0 = !{i1 1}
"""

    def test_fanout_raise_4(self):
        mod, stats = self.check(self.fanout_raise_4)
        self.assertEqual(stats.fanout_raise, 0)

    def test_fanout_raise_4_legacy(self):
        mod, stats = self.check_legacy(self.fanout_raise_4)
        self.assertEqual(stats.fanout_raise, 0)

    fanout_raise_5 = r"""
define i32 @main(i8* %ptr, i1 %cond, i8** %excinfo) {
bb_A:
    call void @NRT_incref(i8* %ptr)
    br i1 %cond, label %bb_B, label %bb_C
bb_B:
    call void @NRT_decref(i8* %ptr)
    br label %common.ret
bb_C:
    store i8* null, i8** %excinfo, !numba_exception_output !0
    br label %common.ret
common.ret:
    %common.ret.op = phi i32 [ 0, %bb_B ], [ 1, %bb_C ]
    ret i32 %common.ret.op
}
!0 = !{i1 1}
"""

    def test_fanout_raise_5(self):
        mod, stats = self.check(self.fanout_raise_5)
        self.assertEqual(stats.fanout_raise, 2)

    def test_fanout_raise_5_legacy(self):
        mod, stats = self.check_legacy(self.fanout_raise_5)
        self.assertEqual(stats.fanout_raise, 2)

    # test case 6 is from https://github.com/numba/llvmlite/issues/1023
    fanout_raise_6 = r"""
define i32 @main(i8* %ptr, i1 %cond1, i1 %cond2, i1 %cond3, i8** %excinfo) {
bb_A:
    call void @NRT_incref(i8* %ptr)
    call void @NRT_incref(i8* %ptr)
    br i1 %cond1, label %bb_B, label %bb_C
bb_B:
    call void @NRT_decref(i8* %ptr)
    br i1 %cond2, label %bb_D, label %bb_E
bb_C:
    store i8* null, i8** %excinfo, !numba_exception_output !0
    ret i32 1
bb_D:
    call void @NRT_decref(i8* %ptr)
    ret i32 0
bb_E:
    call void @NRT_incref(i8* %ptr)
    br i1 %cond3, label %bb_F, label %bb_C
bb_F:
    call void @NRT_decref(i8* %ptr)
    call void @NRT_decref(i8* %ptr)
    ret i32 0
}
!0 = !{i1 1}
"""

    def test_fanout_raise_6(self):
        mod, stats = self.check(self.fanout_raise_6)
        self.assertEqual(stats.fanout_raise, 7)

    def test_fanout_raise_6_legacy(self):
        mod, stats = self.check_legacy(self.fanout_raise_6)
        self.assertEqual(stats.fanout_raise, 7)


if __name__ == '__main__':
    unittest.main()
