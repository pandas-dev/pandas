import unittest
from unittest.case import TestCase
import warnings
import numpy as np

from numba import objmode
from numba.core import ir, compiler
from numba.core import errors
from numba.core.compiler import (
    CompilerBase,
    ReconstructSSA,
)
from numba.core.compiler_machinery import (
    FunctionPass,
    PassManager,
    register_pass,
)
from numba.core.untyped_passes import (
    TranslateByteCode,
    IRProcessing,
)
from numba import njit


class TestIR(unittest.TestCase):

    def test_IRScope(self):
        filename = "<?>"
        top = ir.Scope(parent=None, loc=ir.Loc(filename=filename, line=1))
        local = ir.Scope(parent=top, loc=ir.Loc(filename=filename, line=2))

        apple = local.define('apple', loc=ir.Loc(filename=filename, line=3))
        self.assertIs(local.get('apple'), apple)
        self.assertEqual(len(local.localvars), 1)

        orange = top.define('orange', loc=ir.Loc(filename=filename, line=4))
        self.assertEqual(len(local.localvars), 1)
        self.assertEqual(len(top.localvars), 1)
        self.assertIs(top.get('orange'), orange)
        self.assertIs(local.get('orange'), orange)

        more_orange = local.define('orange', loc=ir.Loc(filename=filename,
                                                        line=5))
        self.assertIs(top.get('orange'), orange)
        self.assertIsNot(local.get('orange'), not orange)
        self.assertIs(local.get('orange'), more_orange)

        try:
            local.define('orange', loc=ir.Loc(filename=filename, line=5))
        except ir.RedefinedError:
            pass
        else:
            self.fail("Expecting an %s" % ir.RedefinedError)


class CheckEquality(unittest.TestCase):

    var_a = ir.Var(None, 'a', ir.unknown_loc)
    var_b = ir.Var(None, 'b', ir.unknown_loc)
    var_c = ir.Var(None, 'c', ir.unknown_loc)
    var_d = ir.Var(None, 'd', ir.unknown_loc)
    var_e = ir.Var(None, 'e', ir.unknown_loc)
    loc1 = ir.Loc('mock', 1, 0)
    loc2 = ir.Loc('mock', 2, 0)
    loc3 = ir.Loc('mock', 3, 0)

    def check(self, base, same=[], different=[]):
        for s in same:
            self.assertTrue(base == s)
        for d in different:
            self.assertTrue(base != d)


class TestIRMeta(CheckEquality):
    """
    Tests IR node meta, like Loc and Scope
    """
    def test_loc(self):
        a = ir.Loc('file', 1, 0)
        b = ir.Loc('file', 1, 0)
        c = ir.Loc('pile', 1, 0)
        d = ir.Loc('file', 2, 0)
        e = ir.Loc('file', 1, 1)
        self.check(a, same=[b,], different=[c, d, e])

        f = ir.Loc('file', 1, 0, maybe_decorator=False)
        g = ir.Loc('file', 1, 0, maybe_decorator=True)
        self.check(a, same=[f, g])

    def test_scope(self):
        parent1 = ir.Scope(None, self.loc1)
        parent2 = ir.Scope(None, self.loc1)
        parent3 = ir.Scope(None, self.loc2)
        self.check(parent1, same=[parent2, parent3,])

        a = ir.Scope(parent1, self.loc1)
        b = ir.Scope(parent1, self.loc1)
        c = ir.Scope(parent1, self.loc2)
        d = ir.Scope(parent3, self.loc1)
        self.check(a, same=[b, c, d])

        # parent1 and parent2 are equal, so children referring to either parent
        # should be equal
        e = ir.Scope(parent2, self.loc1)
        self.check(a, same=[e,])


class TestIRNodes(CheckEquality):
    """
    Tests IR nodes
    """
    def test_terminator(self):
        # terminator base class inst should always be equal
        t1 = ir.Terminator()
        t2 = ir.Terminator()
        self.check(t1, same=[t2])

    def test_jump(self):
        a = ir.Jump(1, self.loc1)
        b = ir.Jump(1, self.loc1)
        c = ir.Jump(1, self.loc2)
        d = ir.Jump(2, self.loc1)
        self.check(a, same=[b, c], different=[d])

    def test_return(self):
        a = ir.Return(self.var_a, self.loc1)
        b = ir.Return(self.var_a, self.loc1)
        c = ir.Return(self.var_a, self.loc2)
        d = ir.Return(self.var_b, self.loc1)
        self.check(a, same=[b, c], different=[d])

    def test_raise(self):
        a = ir.Raise(self.var_a, self.loc1)
        b = ir.Raise(self.var_a, self.loc1)
        c = ir.Raise(self.var_a, self.loc2)
        d = ir.Raise(self.var_b, self.loc1)
        self.check(a, same=[b, c], different=[d])

    def test_staticraise(self):
        a = ir.StaticRaise(AssertionError, None, self.loc1)
        b = ir.StaticRaise(AssertionError, None, self.loc1)
        c = ir.StaticRaise(AssertionError, None, self.loc2)
        e = ir.StaticRaise(AssertionError, ("str",), self.loc1)
        d = ir.StaticRaise(RuntimeError, None, self.loc1)
        self.check(a, same=[b, c], different=[d, e])

    def test_branch(self):
        a = ir.Branch(self.var_a, 1, 2, self.loc1)
        b = ir.Branch(self.var_a, 1, 2, self.loc1)
        c = ir.Branch(self.var_a, 1, 2, self.loc2)
        d = ir.Branch(self.var_b, 1, 2, self.loc1)
        e = ir.Branch(self.var_a, 2, 2, self.loc1)
        f = ir.Branch(self.var_a, 1, 3, self.loc1)
        self.check(a, same=[b, c], different=[d, e, f])

    def test_expr(self):
        a = ir.Expr('some_op', self.loc1)
        b = ir.Expr('some_op', self.loc1)
        c = ir.Expr('some_op', self.loc2)
        d = ir.Expr('some_other_op', self.loc1)
        self.check(a, same=[b, c], different=[d])

    def test_setitem(self):
        a = ir.SetItem(self.var_a, self.var_b, self.var_c, self.loc1)
        b = ir.SetItem(self.var_a, self.var_b, self.var_c, self.loc1)
        c = ir.SetItem(self.var_a, self.var_b, self.var_c, self.loc2)
        d = ir.SetItem(self.var_d, self.var_b, self.var_c, self.loc1)
        e = ir.SetItem(self.var_a, self.var_d, self.var_c, self.loc1)
        f = ir.SetItem(self.var_a, self.var_b, self.var_d, self.loc1)
        self.check(a, same=[b, c], different=[d, e, f])

    def test_staticsetitem(self):
        a = ir.StaticSetItem(self.var_a, 1, self.var_b, self.var_c, self.loc1)
        b = ir.StaticSetItem(self.var_a, 1, self.var_b, self.var_c, self.loc1)
        c = ir.StaticSetItem(self.var_a, 1, self.var_b, self.var_c, self.loc2)
        d = ir.StaticSetItem(self.var_d, 1, self.var_b, self.var_c, self.loc1)
        e = ir.StaticSetItem(self.var_a, 2, self.var_b, self.var_c, self.loc1)
        f = ir.StaticSetItem(self.var_a, 1, self.var_d, self.var_c, self.loc1)
        g = ir.StaticSetItem(self.var_a, 1, self.var_b, self.var_d, self.loc1)
        self.check(a, same=[b, c], different=[d, e, f, g])

    def test_delitem(self):
        a = ir.DelItem(self.var_a, self.var_b, self.loc1)
        b = ir.DelItem(self.var_a, self.var_b, self.loc1)
        c = ir.DelItem(self.var_a, self.var_b, self.loc2)
        d = ir.DelItem(self.var_c, self.var_b, self.loc1)
        e = ir.DelItem(self.var_a, self.var_c, self.loc1)
        self.check(a, same=[b, c], different=[d, e])

    def test_del(self):
        a = ir.Del(self.var_a.name, self.loc1)
        b = ir.Del(self.var_a.name, self.loc1)
        c = ir.Del(self.var_a.name, self.loc2)
        d = ir.Del(self.var_b.name, self.loc1)
        self.check(a, same=[b, c], different=[d])

    def test_setattr(self):
        a = ir.SetAttr(self.var_a, 'foo', self.var_b, self.loc1)
        b = ir.SetAttr(self.var_a, 'foo', self.var_b, self.loc1)
        c = ir.SetAttr(self.var_a, 'foo', self.var_b, self.loc2)
        d = ir.SetAttr(self.var_c, 'foo', self.var_b, self.loc1)
        e = ir.SetAttr(self.var_a, 'bar', self.var_b, self.loc1)
        f = ir.SetAttr(self.var_a, 'foo', self.var_c, self.loc1)
        self.check(a, same=[b, c], different=[d, e, f])

    def test_delattr(self):
        a = ir.DelAttr(self.var_a, 'foo', self.loc1)
        b = ir.DelAttr(self.var_a, 'foo', self.loc1)
        c = ir.DelAttr(self.var_a, 'foo', self.loc2)
        d = ir.DelAttr(self.var_c, 'foo', self.loc1)
        e = ir.DelAttr(self.var_a, 'bar', self.loc1)
        self.check(a, same=[b, c], different=[d, e])

    def test_assign(self):
        a = ir.Assign(self.var_a, self.var_b, self.loc1)
        b = ir.Assign(self.var_a, self.var_b, self.loc1)
        c = ir.Assign(self.var_a, self.var_b, self.loc2)
        d = ir.Assign(self.var_c, self.var_b, self.loc1)
        e = ir.Assign(self.var_a, self.var_c, self.loc1)
        self.check(a, same=[b, c], different=[d, e])

    def test_print(self):
        a = ir.Print((self.var_a,), self.var_b, self.loc1)
        b = ir.Print((self.var_a,), self.var_b, self.loc1)
        c = ir.Print((self.var_a,), self.var_b, self.loc2)
        d = ir.Print((self.var_c,), self.var_b, self.loc1)
        e = ir.Print((self.var_a,), self.var_c, self.loc1)
        self.check(a, same=[b, c], different=[d, e])

    def test_storemap(self):
        a = ir.StoreMap(self.var_a, self.var_b, self.var_c, self.loc1)
        b = ir.StoreMap(self.var_a, self.var_b, self.var_c, self.loc1)
        c = ir.StoreMap(self.var_a, self.var_b, self.var_c, self.loc2)
        d = ir.StoreMap(self.var_d, self.var_b, self.var_c, self.loc1)
        e = ir.StoreMap(self.var_a, self.var_d, self.var_c, self.loc1)
        f = ir.StoreMap(self.var_a, self.var_b, self.var_d, self.loc1)
        self.check(a, same=[b, c], different=[d, e, f])

    def test_yield(self):
        a = ir.Yield(self.var_a, self.loc1, 0)
        b = ir.Yield(self.var_a, self.loc1, 0)
        c = ir.Yield(self.var_a, self.loc2, 0)
        d = ir.Yield(self.var_b, self.loc1, 0)
        e = ir.Yield(self.var_a, self.loc1, 1)
        self.check(a, same=[b, c], different=[d, e])

    def test_enterwith(self):
        a = ir.EnterWith(self.var_a, 0, 1, self.loc1)
        b = ir.EnterWith(self.var_a, 0, 1, self.loc1)
        c = ir.EnterWith(self.var_a, 0, 1, self.loc2)
        d = ir.EnterWith(self.var_b, 0, 1, self.loc1)
        e = ir.EnterWith(self.var_a, 1, 1, self.loc1)
        f = ir.EnterWith(self.var_a, 0, 2, self.loc1)
        self.check(a, same=[b, c], different=[d, e, f])

    def test_arg(self):
        a = ir.Arg('foo', 0, self.loc1)
        b = ir.Arg('foo', 0, self.loc1)
        c = ir.Arg('foo', 0, self.loc2)
        d = ir.Arg('bar', 0, self.loc1)
        e = ir.Arg('foo', 1, self.loc1)
        self.check(a, same=[b, c], different=[d, e])

    def test_const(self):
        a = ir.Const(1, self.loc1)
        b = ir.Const(1, self.loc1)
        c = ir.Const(1, self.loc2)
        d = ir.Const(2, self.loc1)
        self.check(a, same=[b, c], different=[d])

    def test_global(self):
        a = ir.Global('foo', 0, self.loc1)
        b = ir.Global('foo', 0, self.loc1)
        c = ir.Global('foo', 0, self.loc2)
        d = ir.Global('bar', 0, self.loc1)
        e = ir.Global('foo', 1, self.loc1)
        self.check(a, same=[b, c], different=[d, e])

    def test_var(self):
        a = ir.Var(None, 'foo', self.loc1)
        b = ir.Var(None, 'foo', self.loc1)
        c = ir.Var(None, 'foo', self.loc2)
        d = ir.Var(ir.Scope(None, ir.unknown_loc), 'foo', self.loc1)
        e = ir.Var(None, 'bar', self.loc1)
        self.check(a, same=[b, c, d], different=[e])

    def test_undefinedtype(self):
        a = ir.UndefinedType()
        b = ir.UndefinedType()
        self.check(a, same=[b])

    def test_loop(self):
        a = ir.Loop(1, 3)
        b = ir.Loop(1, 3)
        c = ir.Loop(2, 3)
        d = ir.Loop(1, 4)
        self.check(a, same=[b], different=[c, d])

    def test_with(self):
        a = ir.With(1, 3)
        b = ir.With(1, 3)
        c = ir.With(2, 3)
        d = ir.With(1, 4)
        self.check(a, same=[b], different=[c, d])


# used later
_GLOBAL = 1234


class TestIRCompounds(CheckEquality):
    """
    Tests IR concepts that have state
    """
    def test_varmap(self):
        a = ir.VarMap()
        a.define(self.var_a, 'foo')
        a.define(self.var_b, 'bar')

        b = ir.VarMap()
        b.define(self.var_a, 'foo')
        b.define(self.var_b, 'bar')

        c = ir.VarMap()
        c.define(self.var_a, 'foo')
        c.define(self.var_c, 'bar')

        self.check(a, same=[b], different=[c])

    def test_block(self):
        def gen_block():
            parent = ir.Scope(None, self.loc1)
            tmp = ir.Block(parent, self.loc2)
            assign1 = ir.Assign(self.var_a, self.var_b, self.loc3)
            assign2 = ir.Assign(self.var_a, self.var_c, self.loc3)
            assign3 = ir.Assign(self.var_c, self.var_b, self.loc3)
            tmp.append(assign1)
            tmp.append(assign2)
            tmp.append(assign3)
            return tmp

        a = gen_block()
        b = gen_block()
        c = gen_block().append(ir.Assign(self.var_a, self.var_b, self.loc3))

        self.check(a, same=[b], different=[c])

    def test_functionir(self):

        def run_frontend(x):
            return compiler.run_frontend(x, emit_dels=True)

        # this creates a function full of all sorts of things to ensure the IR
        # is pretty involved, it then compares two instances of the compiled
        # function IR to check the IR is the same invariant of objects, and then
        # a tiny mutation is made to the IR in the second function and detection
        # of this change is checked.
        def gen():
            _FREEVAR = 0xCAFE

            def foo(a, b, c=12, d=1j, e=None):
                f = a + b
                a += _FREEVAR
                g = np.zeros(c, dtype=np.complex64)
                h = f + g
                i = 1j / d
                if np.abs(i) > 0:
                    k = h / i
                    l = np.arange(1, c + 1)
                    with objmode():
                        print(e, k)
                    m = np.sqrt(l - g)
                    if np.abs(m[0]) < 1:
                        n = 0
                        for o in range(a):
                            n += 0
                            if np.abs(n) < 3:
                                break
                        n += m[2]
                    p = g / l
                    q = []
                    for r in range(len(p)):
                        q.append(p[r])
                        if r > 4 + 1:
                            with objmode(s='intp', t='complex128'):
                                s = 123
                                t = 5
                            if s > 122:
                                t += s
                        t += q[0] + _GLOBAL

                return f + o + r + t + r + a + n

            return foo

        x = gen()
        y = gen()
        x_ir = run_frontend(x)
        y_ir = run_frontend(y)

        self.assertTrue(x_ir.equal_ir(y_ir))

        def check_diffstr(string, pointing_at=[]):
            lines = string.splitlines()
            for item in pointing_at:
                for l in lines:
                    if l.startswith('->'):
                        if item in l:
                            break
                else:
                    raise AssertionError("Could not find %s " % item)

        self.assertIn("IR is considered equivalent", x_ir.diff_str(y_ir))

        # minor mutation, simply switch branch targets on last branch
        for label in reversed(list(y_ir.blocks.keys())):
            blk = y_ir.blocks[label]
            if isinstance(blk.body[-1], ir.Branch):
                ref = blk.body[-1]
                ref.truebr, ref.falsebr = ref.falsebr, ref.truebr
                break

        check_diffstr(x_ir.diff_str(y_ir), ['branch'])

        z = gen()
        self.assertFalse(x_ir.equal_ir(y_ir))

        z_ir = run_frontend(z)

        change_set = set()
        for label in reversed(list(z_ir.blocks.keys())):
            blk = z_ir.blocks[label]
            ref = blk.body[:-1]
            idx = None
            for i in range(len(ref) - 1):
                # look for two adjacent Del
                if (isinstance(ref[i], ir.Del) and
                        isinstance(ref[i + 1], ir.Del)):
                    idx = i
                    break
            if idx is not None:
                b = blk.body
                change_set.add(str(b[idx + 1]))
                change_set.add(str(b[idx]))
                b[idx], b[idx + 1] = b[idx + 1], b[idx]
                break

        # ensure that a mutation occurred.
        self.assertTrue(change_set)

        self.assertFalse(x_ir.equal_ir(z_ir))
        self.assertEqual(len(change_set), 2)
        for item in change_set:
            self.assertTrue(item.startswith('del '))
        check_diffstr(x_ir.diff_str(z_ir), change_set)

        def foo(a, b):
            c = a * 2
            d = c + b
            e = np.sqrt(d)
            return e

        def bar(a, b): # same as foo
            c = a * 2
            d = c + b
            e = np.sqrt(d)
            return e

        def baz(a, b):
            c = a * 2
            d = b + c
            e = np.sqrt(d + 1)
            return e

        foo_ir = run_frontend(foo)
        bar_ir = run_frontend(bar)
        self.assertTrue(foo_ir.equal_ir(bar_ir))
        self.assertIn("IR is considered equivalent", foo_ir.diff_str(bar_ir))

        baz_ir = run_frontend(baz)
        self.assertFalse(foo_ir.equal_ir(baz_ir))
        tmp = foo_ir.diff_str(baz_ir)
        self.assertIn("Other block contains more statements", tmp)
        check_diffstr(tmp, ["c + b", "b + c"])


class TestIRPedanticChecks(TestCase):
    def test_var_in_scope_assumption(self):
        # Create a pass that clears ir.Scope in ir.Block
        @register_pass(mutates_CFG=False, analysis_only=False)
        class RemoveVarInScope(FunctionPass):
            _name = "_remove_var_in_scope"

            def __init__(self):
                FunctionPass.__init__(self)

            # implement method to do the work, "state" is the internal compiler
            # state from the CompilerBase instance.
            def run_pass(self, state):
                func_ir = state.func_ir
                # walk the blocks
                for blk in func_ir.blocks.values():
                    oldscope = blk.scope
                    # put in an empty Scope
                    blk.scope = ir.Scope(parent=oldscope.parent,
                                         loc=oldscope.loc)
                return True

        # Create a pass that always fails, to stop the compiler
        @register_pass(mutates_CFG=False, analysis_only=False)
        class FailPass(FunctionPass):
            _name = "_fail"

            def __init__(self, *args, **kwargs):
                FunctionPass.__init__(self)

            def run_pass(self, state):
                # This is unreachable. SSA pass should have raised before this
                # pass when run with `error.NumbaPedanticWarning`s raised as
                # errors.
                raise AssertionError("unreachable")

        class MyCompiler(CompilerBase):
            def define_pipelines(self):
                pm = PassManager("testing pm")
                pm.add_pass(TranslateByteCode, "analyzing bytecode")
                pm.add_pass(IRProcessing, "processing IR")
                pm.add_pass(RemoveVarInScope, "_remove_var_in_scope")
                pm.add_pass(ReconstructSSA, "ssa")
                pm.add_pass(FailPass, "_fail")
                pm.finalize()
                return [pm]

        @njit(pipeline_class=MyCompiler)
        def dummy(x):
            # To trigger SSA and the pedantic check, this function must have
            # multiple assignments to the same variable in different blocks.
            a = 1
            b = 2
            if a < b:
                a = 2
            else:
                b = 3
            return a, b

        with warnings.catch_warnings():
            # Make NumbaPedanticWarning an error
            warnings.simplefilter("error", errors.NumbaPedanticWarning)
            # Catch NumbaIRAssumptionWarning
            with self.assertRaises(errors.NumbaIRAssumptionWarning) as raises:
                dummy(1)
            # Verify the error message
            self.assertRegex(
                str(raises.exception),
                r"variable '[a-z]' is not in scope",
            )


if __name__ == '__main__':
    unittest.main()
