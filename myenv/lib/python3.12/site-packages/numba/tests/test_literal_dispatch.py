import numpy as np

import numba
import unittest
from numba.tests.support import TestCase
from numba import njit
from numba.core import types, errors, cgutils
from numba.core.typing import signature
from numba.core.datamodel import models
from numba.core.extending import (
    overload, SentryLiteralArgs, overload_method, register_model, intrinsic,
)
from numba.misc.special import literally


class TestLiteralDispatch(TestCase):
    def check_literal_basic(self, literal_args):
        @njit
        def foo(x):
            return literally(x)

        # Test with int
        for lit in literal_args:
            self.assertEqual(foo(lit), lit)

        for lit, sig in zip(literal_args, foo.signatures):
            self.assertEqual(sig[0].literal_value, lit)

    def test_literal_basic(self):
        self.check_literal_basic([123, 321])
        self.check_literal_basic(["abc", "cb123"])

    def test_literal_nested(self):
        @njit
        def foo(x):
            return literally(x) * 2

        @njit
        def bar(y, x):
            return foo(y) + x

        y, x = 3, 7
        self.assertEqual(bar(y, x), y * 2 + x)
        [foo_sig] = foo.signatures
        self.assertEqual(foo_sig[0], types.literal(y))
        [bar_sig] = bar.signatures
        self.assertEqual(bar_sig[0], types.literal(y))
        self.assertNotIsInstance(bar_sig[1], types.Literal)

    def test_literally_freevar(self):
        # Try referring to numba.literally not in the globals
        import numba

        @njit
        def foo(x):
            return numba.literally(x)

        self.assertEqual(foo(123), 123)
        self.assertEqual(foo.signatures[0][0], types.literal(123))

    def test_mutual_recursion_literal(self):
        def get_functions(decor):
            @decor
            def outer_fac(n, value):
                if n < 1:
                    return value
                return n * inner_fac(n - 1, value)

            @decor
            def inner_fac(n, value):
                if n < 1:
                    return literally(value)
                return n * outer_fac(n - 1, value)

            return outer_fac, inner_fac

        ref_outer_fac, ref_inner_fac = get_functions(lambda x: x)
        outer_fac, inner_fac = get_functions(njit)

        self.assertEqual(outer_fac(10, 12), ref_outer_fac(10, 12))
        self.assertEqual(outer_fac.signatures[0][1].literal_value, 12)
        self.assertEqual(inner_fac.signatures[0][1].literal_value, 12)

        self.assertEqual(inner_fac(11, 13), ref_inner_fac(11, 13))
        self.assertEqual(outer_fac.signatures[1][1].literal_value, 13)
        self.assertEqual(inner_fac.signatures[1][1].literal_value, 13)

    def test_literal_nested_multi_arg(self):
        @njit
        def foo(a, b, c):
            return inner(a, c)

        @njit
        def inner(x, y):
            return x + literally(y)

        kwargs = dict(a=1, b=2, c=3)
        got = foo(**kwargs)
        expect = (lambda a, b, c: a + c)(**kwargs)
        self.assertEqual(got, expect)
        [foo_sig] = foo.signatures
        self.assertEqual(foo_sig[2], types.literal(3))

    def test_unsupported_literal_type(self):
        @njit
        def foo(a, b, c):
            return inner(a, c)

        @njit
        def inner(x, y):
            return x + literally(y)

        arr = np.arange(10)
        with self.assertRaises(errors.LiteralTypingError) as raises:
            foo(a=1, b=2, c=arr)
        self.assertIn("numpy.ndarray", str(raises.exception))

    def test_biliteral(self):
        # Test usecase with more than one literal call
        @njit
        def foo(a, b, c):
            return inner(a, b) + inner(b, c)

        @njit
        def inner(x, y):
            return x + literally(y)

        kwargs = dict(a=1, b=2, c=3)
        got = foo(**kwargs)
        expect = (lambda a, b, c: a + b + b + c)(**kwargs)
        self.assertEqual(got, expect)
        [(type_a, type_b, type_c)] = foo.signatures
        self.assertNotIsInstance(type_a, types.Literal)
        self.assertIsInstance(type_b, types.Literal)
        self.assertEqual(type_b.literal_value, 2)
        self.assertIsInstance(type_c, types.Literal)
        self.assertEqual(type_c.literal_value, 3)

    def test_literally_varargs(self):
        @njit
        def foo(a, *args):
            return literally(args)

        with self.assertRaises(errors.LiteralTypingError):
            foo(1, 2, 3)

        @njit
        def bar(a, b):
            foo(a, b)

        with self.assertRaises(errors.TypingError) as raises:
            bar(1, 2)
        self.assertIn(
            "Cannot request literal type",
            str(raises.exception),
        )

    @unittest.expectedFailure
    def test_literally_defaults(self):
        # Problem with OmittedArg
        @njit
        def foo(a, b=1):
            return (a, literally(b))
        foo(a=1)

    @unittest.expectedFailure
    def test_literally_defaults_inner(self):
        # Problem with Omitted
        @njit
        def foo(a, b=1):
            return (a, literally(b))

        @njit
        def bar(a):
            return foo(a) + 1

        bar(1)

    def test_literally_from_module(self):
        # Problem with Omitted
        @njit
        def foo(x):
            return numba.literally(x)

        got = foo(123)
        self.assertEqual(got, foo.py_func(123))
        self.assertIsInstance(foo.signatures[0][0], types.Literal)

    def test_non_literal(self):
        @njit
        def foo(a, b):
            return literally(1 + a)

        with self.assertRaises(errors.TypingError) as raises:
            foo(1, 2)
        self.assertIn(
            "Invalid use of non-Literal type",
            str(raises.exception),
        )

    def test_inlined_literal(self):
        # Check that literally accepts inlined literal
        @njit
        def foo(a, b):
            v = 1000
            return a + literally(v) + literally(b)

        got = foo(1, 2)
        self.assertEqual(got, foo.py_func(1, 2))

        @njit
        def bar():
            a = 100
            b = 9
            return foo(a=b, b=a)

        got = bar()
        self.assertEqual(got, bar.py_func())

    def test_aliased_variable(self):
        @njit
        def foo(a, b, c):
            def closure(d):
                return literally(d) + 10 * inner(a, b)
            # The inlining of the closure will create an alias to c
            return closure(c)

        @njit
        def inner(x, y):
            return x + literally(y)

        kwargs = dict(a=1, b=2, c=3)
        got = foo(**kwargs)
        expect = (lambda a, b, c: c + 10 * (a + b))(**kwargs)
        self.assertEqual(got, expect)
        [(type_a, type_b, type_c)] = foo.signatures
        self.assertNotIsInstance(type_a, types.Literal)
        self.assertIsInstance(type_b, types.Literal)
        self.assertEqual(type_b.literal_value, 2)
        self.assertIsInstance(type_c, types.Literal)
        self.assertEqual(type_c.literal_value, 3)

    def test_overload_explicit(self):
        # This test represents a more controlled usage with ensuring literal
        # typing for an argument.
        def do_this(x, y):
            return x + y

        @overload(do_this)
        def ov_do_this(x, y):
            SentryLiteralArgs(['x']).for_function(ov_do_this).bind(x, y)
            return lambda x, y: x + y

        @njit
        def foo(a, b):
            return do_this(a, b)

        a = 123
        b = 321
        r = foo(a, b)
        self.assertEqual(r, a + b)
        [type_a, type_b] = foo.signatures[0]
        self.assertIsInstance(type_a, types.Literal)
        self.assertEqual(type_a.literal_value, a)
        self.assertNotIsInstance(type_b, types.Literal)

    def test_overload_implicit(self):
        # This test represents the preferred usage style for using literally
        # in overload. Here, literally() is used inside the "implementation"
        # function of the overload.
        def do_this(x, y):
            return x + y

        @njit
        def hidden(x, y):
            return literally(x) + y

        @overload(do_this)
        def ov_do_this(x, y):
            if isinstance(x, types.Integer):
                # At this point, `x` can be a literal or not
                return lambda x, y: hidden(x, y)

        @njit
        def foo(a, b):
            return do_this(a, b)

        a = 123
        b = 321
        r = foo(a, b)
        self.assertEqual(r, a + b)
        [type_a, type_b] = foo.signatures[0]
        self.assertIsInstance(type_a, types.Literal)
        self.assertEqual(type_a.literal_value, a)
        self.assertNotIsInstance(type_b, types.Literal)

    def test_overload_error_loop(self):
        # Test a case where a infinite compiling loop is caused because a
        # literal type is requested but an error would raise for the
        # literal-ized code path. This causes the overload resolution to
        # retry by "de-literal-izing" the values.
        def do_this(x, y):
            return x + y

        @njit
        def hidden(x, y):
            return literally(x) + y

        @overload(do_this)
        def ov_do_this(x, y):
            if isinstance(y, types.IntegerLiteral):
                # This error is however suppressed because a non-literal
                # version is valid.
                raise errors.NumbaValueError("oops")
            else:
                def impl(x, y):
                    return hidden(x, y)
                return impl

        @njit
        def foo(a, b):
            return do_this(a, literally(b))

        # Expect raising CompilerError to stop re-compiling with duplicated
        # literal typing request.
        with self.assertRaises(errors.CompilerError) as raises:
            foo(a=123, b=321)
        self.assertIn("Repeated literal typing request",
                      str(raises.exception))


class TestLiteralDispatchWithCustomType(TestCase):
    def make_dummy_type(self):
        class Dummy(object):
            def lit(self, a):
                return a

        class DummyType(types.Type):
            def __init__(self):
                super(DummyType, self).__init__(name="dummy")

        @register_model(DummyType)
        class DummyTypeModel(models.StructModel):
            def __init__(self, dmm, fe_type):
                members = []
                super(DummyTypeModel, self).__init__(dmm, fe_type, members)

        @intrinsic
        def init_dummy(typingctx):
            def codegen(context, builder, signature, args):
                dummy = cgutils.create_struct_proxy(
                    signature.return_type)(context, builder)

                return dummy._getvalue()

            sig = signature(DummyType())
            return sig, codegen

        @overload(Dummy)
        def dummy_overload():
            def ctor():
                return init_dummy()

            return ctor

        return (DummyType, Dummy)

    def test_overload_method(self):
        # from issue #5011
        DummyType, Dummy = self.make_dummy_type()

        @overload_method(DummyType, 'lit')
        def lit_overload(self, a):
            def impl(self, a):
                return literally(a)  # <-- using literally here

            return impl

        @njit
        def test_impl(a):
            d = Dummy()

            return d.lit(a)

        # Successful case
        self.assertEqual(test_impl(5), 5)

        # Failing case
        @njit
        def inside(a):
            return test_impl(a + 1)

        with self.assertRaises(errors.TypingError) as raises:
            inside(4)

        self.assertIn("Cannot request literal type.", str(raises.exception))


if __name__ == '__main__':
    unittest.main()
