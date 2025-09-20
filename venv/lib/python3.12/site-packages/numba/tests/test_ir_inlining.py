"""
This tests the inline kwarg to @jit and @overload etc, it has nothing to do with
LLVM or low level inlining.
"""

import operator
import warnings
from itertools import product
import numpy as np

from numba import njit, typeof, literally, prange
from numba.core import types, ir, ir_utils, cgutils, errors, utils
from numba.core.extending import (
    overload,
    overload_method,
    overload_attribute,
    register_model,
    models,
    make_attribute_wrapper,
    intrinsic,
    register_jitable,
)
from numba.core.cpu import InlineOptions
from numba.core.compiler import DefaultPassBuilder, CompilerBase
from numba.core.typed_passes import InlineOverloads
from numba.core.typing import signature
from numba.tests.support import (TestCase, unittest,
                                 MemoryLeakMixin, IRPreservingTestPipeline,
                                 skip_parfors_unsupported,
                                 ignore_internal_warnings)


# this global has the same name as the global in inlining_usecases.py, it
# is here to check that inlined functions bind to their own globals
_GLOBAL1 = -50


@njit(inline='always')
def _global_func(x):
    return x + 1


# to be overloaded
def _global_defn(x):
    return x + 1


@overload(_global_defn, inline='always')
def _global_overload(x):
    return _global_defn


class InliningBase(TestCase):

    _DEBUG = False

    inline_opt_as_bool = {'always': True, 'never': False}

    # --------------------------------------------------------------------------
    # Example cost model

    def sentinel_17_cost_model(self, func_ir):
        # sentinel 17 cost model, this is a fake cost model that will return
        # True (i.e. inline) if the ir.FreeVar(17) is found in the func_ir,
        for blk in func_ir.blocks.values():
            for stmt in blk.body:
                if isinstance(stmt, ir.Assign):
                    if isinstance(stmt.value, ir.FreeVar):
                        if stmt.value.value == 17:
                            return True
        return False

    # --------------------------------------------------------------------------

    def check(self, test_impl, *args, **kwargs):
        inline_expect = kwargs.pop('inline_expect', None)
        assert inline_expect
        block_count = kwargs.pop('block_count', 1)
        assert not kwargs
        for k, v in inline_expect.items():
            assert isinstance(k, str)
            assert isinstance(v, bool)

        j_func = njit(pipeline_class=IRPreservingTestPipeline)(test_impl)

        # check they produce the same answer first!
        self.assertEqual(test_impl(*args), j_func(*args))

        # make sure IR doesn't have branches
        fir = j_func.overloads[j_func.signatures[0]].metadata['preserved_ir']
        fir.blocks = ir_utils.simplify_CFG(fir.blocks)
        if self._DEBUG:
            print("FIR".center(80, "-"))
            fir.dump()
        if block_count != 'SKIP':
            self.assertEqual(len(fir.blocks), block_count)
        block = next(iter(fir.blocks.values()))

        # if we don't expect the function to be inlined then make sure there is
        # 'call' present still
        exprs = [x for x in block.find_exprs()]
        assert exprs
        for k, v in inline_expect.items():
            found = False
            for expr in exprs:
                if getattr(expr, 'op', False) == 'call':
                    func_defn = fir.get_definition(expr.func)
                    found |= func_defn.name == k
                elif ir_utils.is_operator_or_getitem(expr):
                    found |= expr.fn.__name__ == k
            self.assertFalse(found == v)

        return fir  # for use in further analysis


# used in _gen_involved
_GLOBAL = 1234


def _gen_involved():
    _FREEVAR = 0xCAFE

    def foo(a, b, c=12, d=1j, e=None):
        f = a + b
        a += _FREEVAR
        g = np.zeros(c, dtype=np.complex64)
        h = f + g
        i = 1j / d
        # For SSA, zero init, n and t
        n = 0
        t = 0
        if np.abs(i) > 0:
            k = h / i
            l = np.arange(1, c + 1)
            m = np.sqrt(l - g) + e * k
            if np.abs(m[0]) < 1:
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
                    s = 123
                    t = 5
                    if s > 122 - c:
                        t += s
                t += q[0] + _GLOBAL

        return f + o + r + t + r + a + n

    return foo


class TestFunctionInlining(MemoryLeakMixin, InliningBase):

    def test_basic_inline_never(self):
        @njit(inline='never')
        def foo():
            return

        def impl():
            return foo()
        self.check(impl, inline_expect={'foo': False})

    def test_basic_inline_always(self):
        @njit(inline='always')
        def foo():
            return

        def impl():
            return foo()
        self.check(impl, inline_expect={'foo': True})

    def test_basic_inline_combos(self):

        def impl():
            x = foo()
            y = bar()
            z = baz()
            return x, y, z

        opts = (('always'), ('never'))

        for inline_foo, inline_bar, inline_baz in product(opts, opts, opts):

            @njit(inline=inline_foo)
            def foo():
                return

            @njit(inline=inline_bar)
            def bar():
                return

            @njit(inline=inline_baz)
            def baz():
                return

            inline_expect = {'foo': self.inline_opt_as_bool[inline_foo],
                             'bar': self.inline_opt_as_bool[inline_bar],
                             'baz': self.inline_opt_as_bool[inline_baz]}
            self.check(impl, inline_expect=inline_expect)

    @unittest.skip("Need to work out how to prevent this")
    def test_recursive_inline(self):

        @njit(inline='always')
        def foo(x):
            if x == 0:
                return 12
            else:
                foo(x - 1)

        a = 3

        def impl():
            b = 0
            if a > 1:
                b += 1
            foo(5)
            if b < a:
                b -= 1

        self.check(impl, inline_expect={'foo': True})

    def test_freevar_bindings(self):

        def factory(inline, x, y):
            z = x + 12

            @njit(inline=inline)
            def func():
                return (x, y + 3, z)
            return func

        def impl():
            x = foo()
            y = bar()
            z = baz()
            return x, y, z

        opts = (('always'), ('never'))

        for inline_foo, inline_bar, inline_baz in product(opts, opts, opts):

            foo = factory(inline_foo, 10, 20)
            bar = factory(inline_bar, 30, 40)
            baz = factory(inline_baz, 50, 60)

            inline_expect = {'foo': self.inline_opt_as_bool[inline_foo],
                             'bar': self.inline_opt_as_bool[inline_bar],
                             'baz': self.inline_opt_as_bool[inline_baz]}
            self.check(impl, inline_expect=inline_expect)

    def test_global_binding(self):

        def impl():
            x = 19
            return _global_func(x)

        self.check(impl, inline_expect={'_global_func': True})

    def test_inline_from_another_module(self):

        from .inlining_usecases import bar

        def impl():
            z = _GLOBAL1 + 2
            return bar(), z

        self.check(impl, inline_expect={'bar': True})

    def test_inline_from_another_module_w_getattr(self):

        import numba.tests.inlining_usecases as iuc

        def impl():
            z = _GLOBAL1 + 2
            return iuc.bar(), z

        self.check(impl, inline_expect={'bar': True})

    def test_inline_from_another_module_w_2_getattr(self):

        import numba.tests.inlining_usecases  # noqa forces registration
        import numba.tests as nt

        def impl():
            z = _GLOBAL1 + 2
            return nt.inlining_usecases.bar(), z

        self.check(impl, inline_expect={'bar': True})

    def test_inline_from_another_module_as_freevar(self):

        def factory():
            from .inlining_usecases import bar

            @njit(inline='always')
            def tmp():
                return bar()
            return tmp

        baz = factory()

        def impl():
            z = _GLOBAL1 + 2
            return baz(), z

        self.check(impl, inline_expect={'bar': True})

    def test_inline_w_freevar_from_another_module(self):

        from .inlining_usecases import baz_factory

        def gen(a, b):
            bar = baz_factory(a)

            def impl():
                z = _GLOBAL1 + a * b
                return bar(), z, a
            return impl

        impl = gen(10, 20)
        self.check(impl, inline_expect={'bar': True})

    def test_inlining_models(self):

        def s17_caller_model(expr, caller_info, callee_info):
            self.assertIsInstance(expr, ir.Expr)
            self.assertEqual(expr.op, "call")
            return self.sentinel_17_cost_model(caller_info)

        def s17_callee_model(expr, caller_info, callee_info):
            self.assertIsInstance(expr, ir.Expr)
            self.assertEqual(expr.op, "call")
            return self.sentinel_17_cost_model(callee_info)

        # caller has sentinel
        for caller, callee in ((11, 17), (17, 11)):

            @njit(inline=s17_caller_model)
            def foo():
                return callee

            def impl(z):
                x = z + caller
                y = foo()
                return y + 3, x

            self.check(impl, 10, inline_expect={'foo': caller == 17})

        # callee has sentinel
        for caller, callee in ((11, 17), (17, 11)):

            @njit(inline=s17_callee_model)
            def bar():
                return callee

            def impl(z):
                x = z + caller
                y = bar()
                return y + 3, x

            self.check(impl, 10, inline_expect={'bar': callee == 17})

    def test_inline_inside_loop(self):
        @njit(inline='always')
        def foo():
            return 12

        def impl():
            acc = 0.0
            for i in range(5):
                acc += foo()
            return acc

        self.check(impl, inline_expect={'foo': True}, block_count=4)

    def test_inline_inside_closure_inside_loop(self):
        @njit(inline='always')
        def foo():
            return 12

        def impl():
            acc = 0.0
            for i in range(5):
                def bar():
                    return foo() + 7
                acc += bar()
            return acc

        self.check(impl, inline_expect={'foo': True}, block_count=4)

    def test_inline_closure_inside_inlinable_inside_closure(self):
        @njit(inline='always')
        def foo(a):
            def baz():
                return 12 + a
            return baz() + 8

        def impl():
            z = 9

            def bar(x):
                return foo(z) + 7 + x
            return bar(z + 2)

        self.check(impl, inline_expect={'foo': True}, block_count=1)

    def test_inline_involved(self):

        fortran = njit(inline='always')(_gen_involved())

        @njit(inline='always')
        def boz(j):
            acc = 0

            def biz(t):
                return t + acc
            for x in range(j):
                acc += biz(8 + acc) + fortran(2., acc, 1, 12j, biz(acc))
            return acc

        @njit(inline='always')
        def foo(a):
            acc = 0
            for p in range(12):
                tmp = fortran(1, 1, 1, 1, 1)

                def baz(x):
                    return 12 + a + x + tmp
                acc += baz(p) + 8 + boz(p) + tmp
            return acc + baz(2)

        def impl():
            z = 9

            def bar(x):
                return foo(z) + 7 + x
            return bar(z + 2)

        # block count changes with Python version due to bytecode differences.
        if utils.PYVERSION in ((3, 12), (3, 13)):
            bc = 39
        elif utils.PYVERSION in ((3, 10), (3, 11)):
            bc = 35
        else:
            raise NotImplementedError(utils.PYVERSION)

        self.check(impl, inline_expect={'foo': True, 'boz': True,
                                        'fortran': True}, block_count=bc)

    def test_inline_renaming_scheme(self):
        # See #7380, this checks that inlined variables have a name derived from
        # the function they were defined in.

        @njit(inline="always")
        def bar(z):
            x = 5
            y = 10
            return x + y + z

        @njit(pipeline_class=IRPreservingTestPipeline)
        def foo(a, b):
            return bar(a), bar(b)

        self.assertEqual(foo(10, 20), (25, 35))

        # check IR. Look for the `x = 5`... there should be
        # Two lots of `const(int, 5)`, one for each inline
        # The LHS of the assignment will have a name like:
        # TestFunctionInlining_test_inline_renaming_scheme__locals__bar_v2.x
        # Ensure that this is the case!
        func_ir = foo.overloads[foo.signatures[0]].metadata['preserved_ir']
        store = []
        for blk in func_ir.blocks.values():
            for stmt in blk.body:
                if isinstance(stmt, ir.Assign):
                    if isinstance(stmt.value, ir.Const):
                        if stmt.value.value == 5:
                            store.append(stmt)

        self.assertEqual(len(store), 2)
        for i in store:
            name = i.target.name
            basename = self.id().lstrip(self.__module__)
            regex = rf'{basename}__locals__bar_v[0-9]+.x'
            self.assertRegex(name, regex)


class TestRegisterJitableInlining(MemoryLeakMixin, InliningBase):

    def test_register_jitable_inlines(self):

        @register_jitable(inline='always')
        def foo():
            return 1

        def impl():
            foo()

        self.check(impl, inline_expect={'foo': True})


class TestOverloadInlining(MemoryLeakMixin, InliningBase):

    def test_basic_inline_never(self):
        def foo():
            pass

        @overload(foo, inline='never')
        def foo_overload():
            def foo_impl():
                pass
            return foo_impl

        def impl():
            return foo()

        self.check(impl, inline_expect={'foo': False})

    def test_basic_inline_always(self):

        def foo():
            pass

        @overload(foo, inline='always')
        def foo_overload():
            def impl():
                pass
            return impl

        def impl():
            return foo()

        self.check(impl, inline_expect={'foo': True})

    def test_inline_always_kw_no_default(self):
        # pass call arg by name that doesn't have default value
        def foo(a, b):
            return a + b

        @overload(foo, inline='always')
        def overload_foo(a, b):
            return lambda a, b: a + b

        def impl():
            return foo(3, b=4)

        self.check(impl, inline_expect={'foo': True})

    def test_inline_operators_unary(self):

        def impl_inline(x):
            return -x

        def impl_noinline(x):
            return +x

        dummy_unary_impl = lambda x: True
        Dummy, DummyType = self.make_dummy_type()
        setattr(Dummy, '__neg__', dummy_unary_impl)
        setattr(Dummy, '__pos__', dummy_unary_impl)

        @overload(operator.neg, inline='always')
        def overload_dummy_neg(x):
            if isinstance(x, DummyType):
                return dummy_unary_impl

        @overload(operator.pos, inline='never')
        def overload_dummy_pos(x):
            if isinstance(x, DummyType):
                return dummy_unary_impl

        self.check(impl_inline, Dummy(), inline_expect={'neg': True})
        self.check(impl_noinline, Dummy(), inline_expect={'pos': False})

    def test_inline_operators_binop(self):

        def impl_inline(x):
            return x == 1

        def impl_noinline(x):
            return x != 1

        Dummy, DummyType = self.make_dummy_type()

        dummy_binop_impl = lambda a, b: True
        setattr(Dummy, '__eq__', dummy_binop_impl)
        setattr(Dummy, '__ne__', dummy_binop_impl)

        @overload(operator.eq, inline='always')
        def overload_dummy_eq(a, b):
            if isinstance(a, DummyType):
                return dummy_binop_impl

        @overload(operator.ne, inline='never')
        def overload_dummy_ne(a, b):
            if isinstance(a, DummyType):
                return dummy_binop_impl

        self.check(impl_inline, Dummy(), inline_expect={'eq': True})
        self.check(impl_noinline, Dummy(), inline_expect={'ne': False})

    def test_inline_operators_inplace_binop(self):

        def impl_inline(x):
            x += 1

        def impl_noinline(x):
            x -= 1

        Dummy, DummyType = self.make_dummy_type()

        dummy_inplace_binop_impl = lambda a, b: True
        setattr(Dummy, '__iadd__', dummy_inplace_binop_impl)
        setattr(Dummy, '__isub__', dummy_inplace_binop_impl)

        @overload(operator.iadd, inline='always')
        def overload_dummy_iadd(a, b):
            if isinstance(a, DummyType):
                return dummy_inplace_binop_impl

        @overload(operator.isub, inline='never')
        def overload_dummy_isub(a, b):
            if isinstance(a, DummyType):
                return dummy_inplace_binop_impl

        # DummyType is not mutable, so lowering 'inplace_binop' Expr
        # re-uses (requires) copying function definition
        @overload(operator.add, inline='always')
        def overload_dummy_add(a, b):
            if isinstance(a, DummyType):
                return dummy_inplace_binop_impl

        @overload(operator.sub, inline='never')
        def overload_dummy_sub(a, b):
            if isinstance(a, DummyType):
                return dummy_inplace_binop_impl

        self.check(impl_inline, Dummy(), inline_expect={'iadd': True})
        self.check(impl_noinline, Dummy(), inline_expect={'isub': False})

    def test_inline_always_operators_getitem(self):

        def impl(x, idx):
            return x[idx]

        def impl_static_getitem(x):
            return x[1]

        Dummy, DummyType = self.make_dummy_type()

        dummy_getitem_impl = lambda obj, idx: None
        setattr(Dummy, '__getitem__', dummy_getitem_impl)

        @overload(operator.getitem, inline='always')
        def overload_dummy_getitem(obj, idx):
            if isinstance(obj, DummyType):
                return dummy_getitem_impl

        # note getitem and static_getitem Exprs refer to operator.getitem
        # hence they are checked using the same expected key
        self.check(impl, Dummy(), 1, inline_expect={'getitem': True})
        self.check(impl_static_getitem, Dummy(),
                   inline_expect={'getitem': True})

    def test_inline_never_operators_getitem(self):

        def impl(x, idx):
            return x[idx]

        def impl_static_getitem(x):
            return x[1]

        Dummy, DummyType = self.make_dummy_type()

        dummy_getitem_impl = lambda obj, idx: None
        setattr(Dummy, '__getitem__', dummy_getitem_impl)

        @overload(operator.getitem, inline='never')
        def overload_dummy_getitem(obj, idx):
            if isinstance(obj, DummyType):
                return dummy_getitem_impl

        # both getitem and static_getitem Exprs refer to operator.getitem
        # hence they are checked using the same expect key
        self.check(impl, Dummy(), 1, inline_expect={'getitem': False})
        self.check(impl_static_getitem, Dummy(),
                   inline_expect={'getitem': False})

    def test_inline_stararg_error(self):
        def foo(a, *b):
            return a + b[0]

        @overload(foo, inline='always')
        def overload_foo(a, *b):
            return lambda a, *b: a + b[0]

        def impl():
            return foo(3, 3, 5)

        with self.assertRaises(NotImplementedError) as e:
            self.check(impl, inline_expect={'foo': True})

        self.assertIn("Stararg not supported in inliner for arg 1 *b",
                      str(e.exception))

    def test_basic_inline_combos(self):

        def impl():
            x = foo()
            y = bar()
            z = baz()
            return x, y, z

        opts = (('always'), ('never'))

        for inline_foo, inline_bar, inline_baz in product(opts, opts, opts):

            def foo():
                pass

            def bar():
                pass

            def baz():
                pass

            @overload(foo, inline=inline_foo)
            def foo_overload():
                def impl():
                    return
                return impl

            @overload(bar, inline=inline_bar)
            def bar_overload():
                def impl():
                    return
                return impl

            @overload(baz, inline=inline_baz)
            def baz_overload():
                def impl():
                    return
                return impl

            inline_expect = {'foo': self.inline_opt_as_bool[inline_foo],
                             'bar': self.inline_opt_as_bool[inline_bar],
                             'baz': self.inline_opt_as_bool[inline_baz]}
            self.check(impl, inline_expect=inline_expect)

    def test_freevar_bindings(self):

        def impl():
            x = foo()
            y = bar()
            z = baz()
            return x, y, z

        opts = (('always'), ('never'))

        for inline_foo, inline_bar, inline_baz in product(opts, opts, opts):
            # need to repeatedly clobber definitions of foo, bar, baz so
            # @overload binds to the right instance WRT inlining

            def foo():
                x = 10
                y = 20
                z = x + 12
                return (x, y + 3, z)

            def bar():
                x = 30
                y = 40
                z = x + 12
                return (x, y + 3, z)

            def baz():
                x = 60
                y = 80
                z = x + 12
                return (x, y + 3, z)

            def factory(target, x, y, inline=None):
                z = x + 12

                @overload(target, inline=inline)
                def func():
                    def impl():
                        return (x, y + 3, z)
                    return impl

            factory(foo, 10, 20, inline=inline_foo)
            factory(bar, 30, 40, inline=inline_bar)
            factory(baz, 60, 80, inline=inline_baz)

            inline_expect = {'foo': self.inline_opt_as_bool[inline_foo],
                             'bar': self.inline_opt_as_bool[inline_bar],
                             'baz': self.inline_opt_as_bool[inline_baz]}

            self.check(impl, inline_expect=inline_expect)

    def test_global_overload_binding(self):

        def impl():
            z = 19
            return _global_defn(z)

        self.check(impl, inline_expect={'_global_defn': True})

    def test_inline_from_another_module(self):

        from .inlining_usecases import baz

        def impl():
            z = _GLOBAL1 + 2
            return baz(), z

        self.check(impl, inline_expect={'baz': True})

    def test_inline_from_another_module_w_getattr(self):

        import numba.tests.inlining_usecases as iuc

        def impl():
            z = _GLOBAL1 + 2
            return iuc.baz(), z

        self.check(impl, inline_expect={'baz': True})

    def test_inline_from_another_module_w_2_getattr(self):

        import numba.tests.inlining_usecases  # noqa forces registration
        import numba.tests as nt

        def impl():
            z = _GLOBAL1 + 2
            return nt.inlining_usecases.baz(), z

        self.check(impl, inline_expect={'baz': True})

    def test_inline_from_another_module_as_freevar(self):

        def factory():
            from .inlining_usecases import baz

            @njit(inline='always')
            def tmp():
                return baz()
            return tmp

        bop = factory()

        def impl():
            z = _GLOBAL1 + 2
            return bop(), z

        self.check(impl, inline_expect={'baz': True})

    def test_inline_w_freevar_from_another_module(self):

        from .inlining_usecases import bop_factory

        def gen(a, b):
            bar = bop_factory(a)

            def impl():
                z = _GLOBAL1 + a * b
                return bar(), z, a
            return impl

        impl = gen(10, 20)
        self.check(impl, inline_expect={'bar': True})

    def test_inlining_models(self):

        def s17_caller_model(expr, caller_info, callee_info):
            self.assertIsInstance(expr, ir.Expr)
            self.assertEqual(expr.op, "call")
            return self.sentinel_17_cost_model(caller_info.func_ir)

        def s17_callee_model(expr, caller_info, callee_info):
            self.assertIsInstance(expr, ir.Expr)
            self.assertEqual(expr.op, "call")
            return self.sentinel_17_cost_model(callee_info.func_ir)

        # caller has sentinel
        for caller, callee in ((10, 11), (17, 11)):

            def foo():
                return callee

            @overload(foo, inline=s17_caller_model)
            def foo_ol():
                def impl():
                    return callee
                return impl

            def impl(z):
                x = z + caller
                y = foo()
                return y + 3, x

            self.check(impl, 10, inline_expect={'foo': caller == 17})

        # callee has sentinel
        for caller, callee in ((11, 17), (11, 10)):

            def bar():
                return callee

            @overload(bar, inline=s17_callee_model)
            def bar_ol():
                def impl():
                    return callee
                return impl

            def impl(z):
                x = z + caller
                y = bar()
                return y + 3, x

            self.check(impl, 10, inline_expect={'bar': callee == 17})

    def test_multiple_overloads_with_different_inline_characteristics(self):
        # check that having different inlining options for different overloads
        # of the same function works ok

        # this is the Python equiv of the overloads below
        def bar(x):
            if isinstance(typeof(x), types.Float):
                return x + 1234
            else:
                return x + 1

        @overload(bar, inline='always')
        def bar_int_ol(x):
            if isinstance(x, types.Integer):
                def impl(x):
                    return x + 1
                return impl

        @overload(bar, inline='never')
        def bar_float_ol(x):
            if isinstance(x, types.Float):
                def impl(x):
                    return x + 1234
                return impl

        def always_inline_cost_model(*args):
            return True

        @overload(bar, inline=always_inline_cost_model)
        def bar_complex_ol(x):
            if isinstance(x, types.Complex):
                def impl(x):
                    return x + 1
                return impl

        def impl():
            a = bar(1)  # integer literal, should inline
            b = bar(2.3)  # float literal, should not inline
            # complex literal, should inline by virtue of cost model
            c = bar(3j)
            return a + b + c

        # there should still be a `bar` not inlined
        fir = self.check(impl, inline_expect={'bar': False}, block_count=1)

        # check there is one call left in the IR
        block = next(iter(fir.blocks.items()))[1]
        calls = [x for x in block.find_exprs(op='call')]
        self.assertTrue(len(calls) == 1)

        # check that the constant "1234" is not in the IR
        consts = [x.value for x in block.find_insts(ir.Assign)
                  if isinstance(getattr(x, 'value', None), ir.Const)]
        for val in consts:
            self.assertNotEqual(val.value, 1234)

    def test_overload_inline_always_with_literally_in_inlinee(self):
        # See issue #5887

        def foo_ovld(dtype):

            if not isinstance(dtype, types.StringLiteral):
                def foo_noop(dtype):
                    return literally(dtype)
                return foo_noop

            if dtype.literal_value == 'str':
                def foo_as_str_impl(dtype):
                    return 10
                return foo_as_str_impl

            if dtype.literal_value in ('int64', 'float64'):
                def foo_as_num_impl(dtype):
                    return 20
                return foo_as_num_impl

        # define foo for literal str 'str'
        def foo(dtype):
            return 10

        overload(foo, inline='always')(foo_ovld)

        def test_impl(dtype):
            return foo(dtype)

        # check literal dispatch on 'str'
        dtype = 'str'
        self.check(test_impl, dtype, inline_expect={'foo': True})

        # redefine foo to be correct for literal str 'int64'
        def foo(dtype):
            return 20
        overload(foo, inline='always')(foo_ovld)

        # check literal dispatch on 'int64'
        dtype = 'int64'
        self.check(test_impl, dtype, inline_expect={'foo': True})

    def test_inline_always_ssa(self):
        # Make sure IR inlining uses SSA properly. Test for #6721.

        dummy_true = True

        def foo(A):
            return True

        @overload(foo, inline="always")
        def foo_overload(A):

            def impl(A):
                s = dummy_true
                for i in range(len(A)):
                    dummy = dummy_true
                    if A[i]:
                        dummy = A[i]
                    s *= dummy
                return s
            return impl

        def impl():
            return foo(np.array([True, False, True]))

        self.check(impl, block_count='SKIP', inline_expect={'foo': True})

    def test_inline_always_ssa_scope_validity(self):
        # Make sure IR inlining correctly updates the scope(s). See #7802

        def bar():
            b = 5
            while b > 1:
                b //= 2

            return 10

        @overload(bar, inline="always")
        def bar_impl():
            return bar

        @njit
        def foo():
            bar()

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', errors.NumbaIRAssumptionWarning)
            ignore_internal_warnings()
            self.assertEqual(foo(), foo.py_func())

        # There should be no warnings as the IR scopes should be consistent with
        # the IR involved.
        self.assertEqual(len(w), 0)


class TestOverloadMethsAttrsInlining(InliningBase):
    def setUp(self):
        self.make_dummy_type()
        super(TestOverloadMethsAttrsInlining, self).setUp()

    def check_method(self, test_impl, args, expected, block_count,
                     expects_inlined=True):
        j_func = njit(pipeline_class=IRPreservingTestPipeline)(test_impl)
        # check they produce the same answer first!
        self.assertEqual(j_func(*args), expected)

        # make sure IR doesn't have branches
        fir = j_func.overloads[j_func.signatures[0]].metadata['preserved_ir']
        fir.blocks = fir.blocks
        self.assertEqual(len(fir.blocks), block_count)
        if expects_inlined:
            # assert no calls
            for block in fir.blocks.values():
                calls = list(block.find_exprs('call'))
                self.assertFalse(calls)
        else:
            # assert has call
            allcalls = []
            for block in fir.blocks.values():
                allcalls += list(block.find_exprs('call'))
            self.assertTrue(allcalls)

    def check_getattr(self, test_impl, args, expected, block_count,
                      expects_inlined=True):
        j_func = njit(pipeline_class=IRPreservingTestPipeline)(test_impl)
        # check they produce the same answer first!
        self.assertEqual(j_func(*args), expected)

        # make sure IR doesn't have branches
        fir = j_func.overloads[j_func.signatures[0]].metadata['preserved_ir']
        fir.blocks = fir.blocks
        self.assertEqual(len(fir.blocks), block_count)
        if expects_inlined:
            # assert no getattr
            for block in fir.blocks.values():
                getattrs = list(block.find_exprs('getattr'))
                self.assertFalse(getattrs)
        else:
            # assert has getattr
            allgetattrs = []
            for block in fir.blocks.values():
                allgetattrs += list(block.find_exprs('getattr'))
            self.assertTrue(allgetattrs)

    def test_overload_method_default_args_always(self):
        Dummy, DummyType = self.make_dummy_type()

        @overload_method(DummyType, "inline_method", inline='always')
        def _get_inlined_method(obj, val=None, val2=None):
            def get(obj, val=None, val2=None):
                return ("THIS IS INLINED", val, val2)
            return get

        def foo(obj):
            return obj.inline_method(123), obj.inline_method(val2=321)

        self.check_method(
            test_impl=foo,
            args=[Dummy()],
            expected=(("THIS IS INLINED", 123, None),
                      ("THIS IS INLINED", None, 321)),
            block_count=1,
        )

    def make_overload_method_test(self, costmodel, should_inline):
        def costmodel(*args):
            return should_inline

        Dummy, DummyType = self.make_dummy_type()

        @overload_method(DummyType, "inline_method", inline=costmodel)
        def _get_inlined_method(obj, val):
            def get(obj, val):
                return ("THIS IS INLINED!!!", val)
            return get

        def foo(obj):
            return obj.inline_method(123)

        self.check_method(
            test_impl=foo,
            args=[Dummy()],
            expected=("THIS IS INLINED!!!", 123),
            block_count=1,
            expects_inlined=should_inline,
        )

    def test_overload_method_cost_driven_always(self):
        self.make_overload_method_test(
            costmodel='always',
            should_inline=True,
        )

    def test_overload_method_cost_driven_never(self):
        self.make_overload_method_test(
            costmodel='never',
            should_inline=False,
        )

    def test_overload_method_cost_driven_must_inline(self):
        self.make_overload_method_test(
            costmodel=lambda *args: True,
            should_inline=True,
        )

    def test_overload_method_cost_driven_no_inline(self):
        self.make_overload_method_test(
            costmodel=lambda *args: False,
            should_inline=False,
        )

    def make_overload_attribute_test(self, costmodel, should_inline):
        Dummy, DummyType = self.make_dummy_type()

        @overload_attribute(DummyType, "inlineme", inline=costmodel)
        def _get_inlineme(obj):
            def get(obj):
                return "MY INLINED ATTRS"
            return get

        def foo(obj):
            return obj.inlineme

        self.check_getattr(
            test_impl=foo,
            args=[Dummy()],
            expected="MY INLINED ATTRS",
            block_count=1,
            expects_inlined=should_inline,
        )

    def test_overload_attribute_always(self):
        self.make_overload_attribute_test(
            costmodel='always',
            should_inline=True,
        )

    def test_overload_attribute_never(self):
        self.make_overload_attribute_test(
            costmodel='never',
            should_inline=False,
        )

    def test_overload_attribute_costmodel_must_inline(self):
        self.make_overload_attribute_test(
            costmodel=lambda *args: True,
            should_inline=True,
        )

    def test_overload_attribute_costmodel_no_inline(self):
        self.make_overload_attribute_test(
            costmodel=lambda *args: False,
            should_inline=False,
        )


class TestGeneralInlining(MemoryLeakMixin, InliningBase):

    def test_with_inlined_and_noninlined_variants(self):
        # This test is contrived and was to demonstrate fixing a bug in the
        # template walking logic where inlinable and non-inlinable definitions
        # would not mix.

        @overload(len, inline='always')
        def overload_len(A):
            if False:
                return lambda A: 10

        def impl():
            return len([2, 3, 4])

        # len(list) won't be inlined because the overload above doesn't apply
        self.check(impl, inline_expect={'len': False})

    def test_with_kwargs(self):

        def foo(a, b=3, c=5):
            return a + b + c

        @overload(foo, inline='always')
        def overload_foo(a, b=3, c=5):
            def impl(a, b=3, c=5):
                return a + b + c
            return impl

        def impl():
            return foo(3, c=10)

        self.check(impl, inline_expect={'foo': True})

    def test_with_kwargs2(self):

        @njit(inline='always')
        def bar(a, b=12, c=9):
            return a + b

        def impl(a, b=7, c=5):
            return bar(a + b, c=19)

        self.check(impl, 3, 4, inline_expect={'bar': True})

    def test_inlining_optional_constant(self):
        # This testcase causes `b` to be a Optional(bool) constant once it is
        # inlined into foo().
        @njit(inline='always')
        def bar(a=None, b=None):
            if b is None:
                b = 123     # this changes the type of `b` due to lack of SSA
            return (a, b)

        def impl():
            return bar(), bar(123), bar(b=321)

        self.check(impl, block_count='SKIP', inline_expect={'bar': True})


class TestInlineOptions(TestCase):

    def test_basic(self):
        always = InlineOptions('always')
        self.assertTrue(always.is_always_inline)
        self.assertFalse(always.is_never_inline)
        self.assertFalse(always.has_cost_model)
        self.assertEqual(always.value, 'always')

        never = InlineOptions('never')
        self.assertFalse(never.is_always_inline)
        self.assertTrue(never.is_never_inline)
        self.assertFalse(never.has_cost_model)
        self.assertEqual(never.value, 'never')

        def cost_model(x):
            return x
        model = InlineOptions(cost_model)
        self.assertFalse(model.is_always_inline)
        self.assertFalse(model.is_never_inline)
        self.assertTrue(model.has_cost_model)
        self.assertIs(model.value, cost_model)


class TestInlineMiscIssues(TestCase):

    def test_issue4691(self):
        def output_factory(array, dtype):
            pass

        @overload(output_factory, inline='always')
        def ol_output_factory(array, dtype):
            if isinstance(array, types.npytypes.Array):
                def impl(array, dtype):
                    shape = array.shape[3:]
                    return np.zeros(shape, dtype=dtype)

                return impl

        @njit(nogil=True)
        def fn(array):
            out = output_factory(array, array.dtype)
            return out

        @njit(nogil=True)
        def fn2(array):
            return np.zeros(array.shape[3:], dtype=array.dtype)

        fn(np.ones((10, 20, 30, 40, 50)))
        fn2(np.ones((10, 20, 30, 40, 50)))

    def test_issue4693(self):

        @njit(inline='always')
        def inlining(array):
            if array.ndim != 1:
                raise ValueError("Invalid number of dimensions")

            return array

        @njit
        def fn(array):
            return inlining(array)

        fn(np.zeros(10))

    def test_issue5476(self):
        # Actual issue has the ValueError passed as an arg to `inlining` so is
        # a constant inference error
        @njit(inline='always')
        def inlining():
            msg = 'Something happened'
            raise ValueError(msg)

        @njit
        def fn():
            return inlining()

        with self.assertRaises(ValueError) as raises:
            fn()

        self.assertIn("Something happened", str(raises.exception))

    def test_issue5792(self):
        # Issue is that overloads cache their IR and closure inliner was
        # manipulating the cached IR in a way that broke repeated inlines.

        class Dummy:
            def __init__(self, data):
                self.data = data

            def div(self, other):
                return data / other.data

        class DummyType(types.Type):
            def __init__(self, data):
                self.data = data
                super().__init__(name=f'Dummy({self.data})')

        @register_model(DummyType)
        class DummyTypeModel(models.StructModel):
            def __init__(self, dmm, fe_type):
                members = [
                    ('data', fe_type.data),
                ]
                super().__init__(dmm, fe_type, members)

        make_attribute_wrapper(DummyType, 'data', '_data')

        @intrinsic
        def init_dummy(typingctx, data):
            def codegen(context, builder, sig, args):
                typ = sig.return_type
                data, = args
                dummy = cgutils.create_struct_proxy(typ)(context, builder)
                dummy.data = data

                if context.enable_nrt:
                    context.nrt.incref(builder, sig.args[0], data)

                return dummy._getvalue()

            ret_typ = DummyType(data)
            sig = signature(ret_typ, data)

            return sig, codegen

        @overload(Dummy, inline='always')
        def dummy_overload(data):
            def ctor(data):
                return init_dummy(data)

            return ctor

        @overload_method(DummyType, 'div', inline='always')
        def div_overload(self, other):
            def impl(self, other):
                return self._data / other._data

            return impl

        @njit
        def test_impl(data, other_data):
            dummy = Dummy(data) # ctor inlined once
            other = Dummy(other_data)  # ctor inlined again

            return dummy.div(other)

        data = 1.
        other_data = 2.
        res = test_impl(data, other_data)
        self.assertEqual(res, data / other_data)

    def test_issue5824(self):
        """ Similar to the above test_issue5792, checks mutation of the inlinee
        IR is local only"""

        class CustomCompiler(CompilerBase):

            def define_pipelines(self):
                pm = DefaultPassBuilder.define_nopython_pipeline(self.state)
                # Run the inliner twice!
                pm.add_pass_after(InlineOverloads, InlineOverloads)
                pm.finalize()
                return [pm]

        def bar(x):
            ...

        @overload(bar, inline='always')
        def ol_bar(x):
            if isinstance(x, types.Integer):
                def impl(x):
                    return x + 1.3
                return impl

        @njit(pipeline_class=CustomCompiler)
        def foo(z):
            return bar(z), bar(z)

        self.assertEqual(foo(10), (11.3, 11.3))

    @skip_parfors_unsupported
    def test_issue7380(self):
        # This checks that inlining a function containing a loop into another
        # loop where the induction variable in both loops is the same doesn't
        # end up with a name collision. Parfors can detect this so it is used.
        # See: https://github.com/numba/numba/issues/7380

        # Check Numba inlined function passes

        @njit(inline="always")
        def bar(x):
            for i in range(x.size):
                x[i] += 1

        @njit(parallel=True)
        def foo(a):
            for i in prange(a.shape[0]):
                bar(a[i])

        a = np.ones((10, 10))
        foo(a) # run
        # check mutation of data is correct
        self.assertPreciseEqual(a, 2 * np.ones_like(a))

        # Check manually inlined equivalent function fails
        @njit(parallel=True)
        def foo_bad(a):
            for i in prange(a.shape[0]):
                x = a[i]
                for i in range(x.size):
                    x[i] += 1

        with self.assertRaises(errors.UnsupportedRewriteError) as e:
            foo_bad(a)

        self.assertIn("Overwrite of parallel loop index", str(e.exception))


if __name__ == '__main__':
    unittest.main()
