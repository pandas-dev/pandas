from numba import njit
from numba.core import errors
from numba.core.extending import overload
import numpy as np

import unittest


@njit
def consumer(func, *args):
    return func(*args)


@njit
def consumer2arg(func1, func2):
    return func2(func1)


_global = 123


class TestMakeFunctionToJITFunction(unittest.TestCase):
    """
    This tests the pass that converts ir.Expr.op == make_function (i.e. closure)
    into a JIT function.
    """
    # NOTE: testing this is a bit tricky. The function receiving a JIT'd closure
    # must also be under JIT control so as to handle the JIT'd closure
    # correctly, however, in the case of running the test implementations in the
    # interpreter, the receiving function cannot be JIT'd else it will receive
    # the Python closure and then complain about pyobjects as arguments.
    # The way around this is to use a factory function to close over either the
    # jitted or standard python function as the consumer depending on context.

    def test_escape(self):

        def impl_factory(consumer_func):
            def impl():
                def inner():
                    return 10
                return consumer_func(inner)
            return impl

        cfunc = njit(impl_factory(consumer))
        impl = impl_factory(consumer.py_func)

        self.assertEqual(impl(), cfunc())

    def test_nested_escape(self):

        def impl_factory(consumer_func):
            def impl():
                def inner():
                    return 10

                def innerinner(x):
                    return x()
                return consumer_func(inner, innerinner)
            return impl

        cfunc = njit(impl_factory(consumer2arg))
        impl = impl_factory(consumer2arg.py_func)

        self.assertEqual(impl(), cfunc())

    def test_closure_in_escaper(self):

        def impl_factory(consumer_func):
            def impl():
                def callinner():
                    def inner():
                        return 10
                    return inner()
                return consumer_func(callinner)
            return impl

        cfunc = njit(impl_factory(consumer))
        impl = impl_factory(consumer.py_func)

        self.assertEqual(impl(), cfunc())

    def test_close_over_consts(self):

        def impl_factory(consumer_func):
            def impl():
                y = 10

                def callinner(z):
                    return y + z + _global
                return consumer_func(callinner, 6)
            return impl

        cfunc = njit(impl_factory(consumer))
        impl = impl_factory(consumer.py_func)

        self.assertEqual(impl(), cfunc())

    def test_close_over_consts_w_args(self):

        def impl_factory(consumer_func):
            def impl(x):
                y = 10

                def callinner(z):
                    return y + z + _global
                return consumer_func(callinner, x)
            return impl

        cfunc = njit(impl_factory(consumer))
        impl = impl_factory(consumer.py_func)

        a = 5
        self.assertEqual(impl(a), cfunc(a))

    def test_with_overload(self):

        def foo(func, *args):
            nargs = len(args)
            if nargs == 1:
                return func(*args)
            elif nargs == 2:
                return func(func(*args))

        @overload(foo)
        def foo_ol(func, *args):
            # specialise on the number of args, as per `foo`
            nargs = len(args)
            if nargs == 1:
                def impl(func, *args):
                    return func(*args)
                return impl
            elif nargs == 2:
                def impl(func, *args):
                    return func(func(*args))
                return impl

        def impl_factory(consumer_func):
            def impl(x):
                y = 10

                def callinner(*z):
                    return y + np.sum(np.asarray(z)) + _global
                # run both specialisations, 1 arg, and 2 arg.
                return foo(callinner, x), foo(callinner, x, x)
            return impl

        cfunc = njit(impl_factory(consumer))
        impl = impl_factory(consumer.py_func)

        a = 5
        self.assertEqual(impl(a), cfunc(a))

    def test_basic_apply_like_case(self):
        def apply(array, func):
            return func(array)

        @overload(apply)
        def ov_apply(array, func):
            return lambda array, func: func(array)

        def impl(array):
            def mul10(x):
                return x * 10
            return apply(array, mul10)

        cfunc = njit(impl)

        a = np.arange(10)
        np.testing.assert_allclose(impl(a), cfunc(a))

    @unittest.skip("Needs option/flag inheritance to work")
    def test_jit_option_inheritance(self):

        def impl_factory(consumer_func):
            def impl(x):
                def inner(val):
                    return 1 / val
                return consumer_func(inner, x)
            return impl

        cfunc = njit(error_model='numpy')(impl_factory(consumer))
        impl = impl_factory(consumer.py_func)

        a = 0
        self.assertEqual(impl(a), cfunc(a))

    # this needs true SSA to be able to work correctly, check error for now
    def test_multiply_defined_freevar(self):

        @njit
        def impl(c):
            if c:
                x = 3

                def inner(y):
                    return y + x

                r = consumer(inner, 1)
            else:
                x = 6

                def inner(y):
                    return y + x

                r = consumer(inner, 2)
            return r

        with self.assertRaises(errors.TypingError) as e:
            impl(1)

        self.assertIn("Cannot capture a constant value for variable",
                      str(e.exception))

    def test_non_const_in_escapee(self):

        @njit
        def impl(x):
            z = np.arange(x)

            def inner(val):
                return 1 + z + val  # z is non-const freevar
            return consumer(inner, x)

        with self.assertRaises(errors.TypingError) as e:
            impl(1)

        self.assertIn("Cannot capture the non-constant value associated",
                      str(e.exception))

    def test_escape_with_kwargs(self):

        def impl_factory(consumer_func):
            def impl():
                t = 12

                def inner(a, b, c, mydefault1=123, mydefault2=456):
                    z = 4
                    return mydefault1 + mydefault2 + z + t + a + b + c
                # this is awkward, top and tail closure inlining with a escapees
                # in the middle that do/don't have defaults.
                return (inner(1, 2, 5, 91, 53),
                        consumer_func(inner, 1, 2, 3, 73),
                        consumer_func(inner, 1, 2, 3,),
                        inner(1, 2, 4))
            return impl

        cfunc = njit(impl_factory(consumer))
        impl = impl_factory(consumer.py_func)

        np.testing.assert_allclose(impl(), cfunc())

    def test_escape_with_kwargs_override_kwargs(self):

        @njit
        def specialised_consumer(func, *args):
            x, y, z = args  # unpack to avoid `CALL_FUNCTION_EX`
            a = func(x, y, z, mydefault1=1000)
            b = func(x, y, z, mydefault2=1000)
            c = func(x, y, z, mydefault1=1000, mydefault2=1000)
            return a + b + c

        def impl_factory(consumer_func):
            def impl():
                t = 12

                def inner(a, b, c, mydefault1=123, mydefault2=456):
                    z = 4
                    return mydefault1 + mydefault2 + z + t + a + b + c
                # this is awkward, top and tail closure inlining with a escapees
                # in the middle that get defaults specified in the consumer
                return (inner(1, 2, 5, 91, 53),
                        consumer_func(inner, 1, 2, 11),
                        consumer_func(inner, 1, 2, 3,),
                        inner(1, 2, 4))
            return impl

        cfunc = njit(impl_factory(specialised_consumer))
        impl = impl_factory(specialised_consumer.py_func)

        np.testing.assert_allclose(impl(), cfunc())


if __name__ == '__main__':
    unittest.main()
