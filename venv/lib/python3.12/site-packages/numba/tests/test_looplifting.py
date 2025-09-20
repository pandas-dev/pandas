from io import StringIO
import numpy as np

from numba.core import types
from numba.core.compiler import compile_extra, Flags
from numba.tests.support import TestCase, tag, MemoryLeakMixin
import unittest


looplift_flags = Flags()
looplift_flags.force_pyobject = True
looplift_flags.enable_looplift = True

pyobject_looplift_flags = looplift_flags.copy()
pyobject_looplift_flags.enable_pyobject_looplift = True


def compile_isolated(pyfunc, argtypes, **kwargs):
    from numba.core.registry import cpu_target

    kwargs.setdefault('return_type', None)
    kwargs.setdefault('locals', {})
    return compile_extra(
        cpu_target.typing_context,
        cpu_target.target_context,
        pyfunc,
        argtypes,
        **kwargs,
    )


def lift1(x):
    # Outer needs object mode because of np.empty()
    a = np.empty(3)
    for i in range(a.size):
        # Inner is nopython-compliant
        a[i] = x
    return a


def lift2(x):
    # Outer needs object mode because of np.empty()
    a = np.empty((3, 4))
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            # Inner is nopython-compliant
            a[i, j] = x
    return a


def lift3(x):
    # Output variable from the loop
    _ = object()
    a = np.arange(5, dtype=np.int64)
    c = 0
    for i in range(a.shape[0]):
        c += a[i] * x
    return c

def lift4(x):
    # Output two variables from the loop
    _ = object()
    a = np.arange(5, dtype=np.int64)
    c = 0
    d = 0
    for i in range(a.shape[0]):
        c += a[i] * x
        d += c
    return c + d

def lift5(x):
    _ = object()
    a = np.arange(4)
    for i in range(a.shape[0]):
        # Inner has a break statement
        if i > 2:
            break
    return a

def lift_gen1(x):
    # Outer needs object mode because of np.empty()
    a = np.empty(3)
    yield 0
    for i in range(a.size):
        # Inner is nopython-compliant
        a[i] = x
    yield np.sum(a)

def lift_issue2561():
    np.empty(1)   # This forces objectmode because no nrt
    for i in range(10):
        for j in range(10):
            return 1
    return 2

def reject1(x):
    a = np.arange(4)
    for i in range(a.shape[0]):
        # Inner returns a variable from outer "scope" => cannot loop-lift
        return a
    return a


def reject_gen1(x):
    _ = object()
    a = np.arange(4)
    for i in range(a.shape[0]):
        # Inner is a generator => cannot loop-lift
        yield a[i]

def reject_gen2(x):
    _ = object()
    a = np.arange(3)
    for i in range(a.size):
        # Middle has a yield => cannot loop-lift
        res = a[i] + x
        for j in range(i):
            # Inner is nopython-compliant, but the current algorithm isn't
            # able to separate it.
            res = res ** 2
        yield res

def reject_npm1(x):
    a = np.empty(3, dtype=np.int32)
    for i in range(a.size):
        # Inner uses object() => cannot loop-lift
        _ = object()
        a[i] = np.arange(i + 1)[i]

    return a


class TestLoopLifting(MemoryLeakMixin, TestCase):

    def try_lift(self, pyfunc, argtypes):
        from numba.core.registry import cpu_target

        cres = compile_extra(
            cpu_target.typing_context,
            cpu_target.target_context,
            pyfunc, argtypes,
            return_type=None, flags=looplift_flags, locals={},
        )
        # One lifted loop
        self.assertEqual(len(cres.lifted), 1)
        return cres

    def assert_lifted_native(self, cres):
        # Check if we have lifted in nopython mode
        jitloop = cres.lifted[0]
        [loopcres] = jitloop.overloads.values()
        self.assertTrue(loopcres.fndesc.native)  # Lifted function is native

    def check_lift_ok(self, pyfunc, argtypes, args):
        """
        Check that pyfunc can loop-lift even in nopython mode.
        """
        cres = self.try_lift(pyfunc, argtypes)
        expected = pyfunc(*args)
        got = cres.entry_point(*args)
        self.assert_lifted_native(cres)
        # Check return values
        self.assertPreciseEqual(expected, got)

    def check_lift_generator_ok(self, pyfunc, argtypes, args):
        """
        Check that pyfunc (a generator function) can loop-lift even in
        nopython mode.
        """
        cres = self.try_lift(pyfunc, argtypes)
        expected = list(pyfunc(*args))
        got = list(cres.entry_point(*args))
        self.assert_lifted_native(cres)
        # Check return values
        self.assertPreciseEqual(expected, got)

    def check_no_lift(self, pyfunc, argtypes, args):
        """
        Check that pyfunc can't loop-lift.
        """
        cres = compile_isolated(pyfunc, argtypes,
                                flags=looplift_flags)
        self.assertFalse(cres.lifted)
        expected = pyfunc(*args)
        got = cres.entry_point(*args)
        # Check return values
        self.assertPreciseEqual(expected, got)

    def check_no_lift_generator(self, pyfunc, argtypes, args):
        """
        Check that pyfunc (a generator function) can't loop-lift.
        """
        cres = compile_isolated(pyfunc, argtypes,
                                flags=looplift_flags)
        self.assertFalse(cres.lifted)
        expected = list(pyfunc(*args))
        got = list(cres.entry_point(*args))
        self.assertPreciseEqual(expected, got)

    def test_lift1(self):
        self.check_lift_ok(lift1, (types.intp,), (123,))

    def test_lift2(self):
        self.check_lift_ok(lift2, (types.intp,), (123,))

    def test_lift3(self):
        self.check_lift_ok(lift3, (types.intp,), (123,))

    def test_lift4(self):
        self.check_lift_ok(lift4, (types.intp,), (123,))

    def test_lift5(self):
        self.check_lift_ok(lift5, (types.intp,), (123,))

    def test_lift_issue2561(self):
        self.check_lift_ok(lift_issue2561, (), ())

    def test_lift_gen1(self):
        self.check_lift_generator_ok(lift_gen1, (types.intp,), (123,))

    def test_reject1(self):
        self.check_no_lift(reject1, (types.intp,), (123,))

    def test_reject_gen1(self):
        self.check_no_lift_generator(reject_gen1, (types.intp,), (123,))

    def test_reject_gen2(self):
        self.check_no_lift_generator(reject_gen2, (types.intp,), (123,))


class TestLoopLiftingAnnotate(TestCase):
    def test_annotate_1(self):
        """
        Verify that annotation works as expected with one lifted loop
        """
        from numba import jit

        # dummy function to force objmode
        def bar():
            pass

        def foo(x):
            bar()  # force obj
            for i in range(x.size):
                x[i] += 1

            return x

        cfoo = jit(forceobj=True)(foo)

        x = np.arange(10)
        xcopy = x.copy()
        r = cfoo(x)
        np.testing.assert_equal(r, xcopy + 1)

        buf = StringIO()
        cfoo.inspect_types(file=buf)
        annotation = buf.getvalue()
        buf.close()

        self.assertIn("The function contains lifted loops", annotation)
        line = foo.__code__.co_firstlineno + 2  # 2 lines down from func head
        self.assertIn("Loop at line {line}".format(line=line), annotation)
        self.assertIn("Has 1 overloads", annotation)

    def test_annotate_2(self):
        """
        Verify that annotation works as expected with two lifted loops
        """
        from numba import jit

        # dummy function to force objmode
        def bar():
            pass

        def foo(x):
            bar()  # force obj
            # first lifted loop
            for i in range(x.size):
                x[i] += 1
            # second lifted loop
            for j in range(x.size):
                x[j] *= 2
            return x

        cfoo = jit(forceobj=True)(foo)

        x = np.arange(10)
        xcopy = x.copy()
        r = cfoo(x)
        np.testing.assert_equal(r, (xcopy + 1) * 2)

        buf = StringIO()
        cfoo.inspect_types(file=buf)
        annotation = buf.getvalue()
        buf.close()

        self.assertIn("The function contains lifted loops", annotation)
        line1 = foo.__code__.co_firstlineno + 3  # 3 lines down from func head
        line2 = foo.__code__.co_firstlineno + 6  # 6 lines down from func head
        self.assertIn("Loop at line {line}".format(line=line1), annotation)
        self.assertIn("Loop at line {line}".format(line=line2), annotation)


class TestLoopLiftingInAction(MemoryLeakMixin, TestCase):
    def assert_has_lifted(self, jitted, loopcount):
        lifted = jitted.overloads[jitted.signatures[0]].lifted
        self.assertEqual(len(lifted), loopcount)

    def test_issue_734(self):
        from numba import jit, void, int32, double

        @jit(void(int32, double[:]), forceobj=True)
        def forloop_with_if(u, a):
            if u == 0:
                for i in range(a.shape[0]):
                    a[i] = a[i] * 2.0
            else:
                for i in range(a.shape[0]):
                    a[i] = a[i] + 1.0

        for u in (0, 1):
            nb_a = np.arange(10, dtype='int32')
            np_a = np.arange(10, dtype='int32')
            forloop_with_if(u, nb_a)
            forloop_with_if.py_func(u, np_a)
            self.assertPreciseEqual(nb_a, np_a)

    def test_issue_812(self):
        from numba import jit

        @jit('f8[:](f8[:])', forceobj=True)
        def test(x):
            res = np.zeros(len(x))
            ind = 0
            for ii in range(len(x)):
                ind += 1
                res[ind] = x[ind]
                if x[ind] >= 10:
                    break

            # Invalid loopjitting will miss the usage of `ind` in the
            # following loop.
            for ii in range(ind + 1, len(x)):
                res[ii] = 0
            return res

        x = np.array([1., 4, 2, -3, 5, 2, 10, 5, 2, 6])
        np.testing.assert_equal(test.py_func(x), test(x))

    def test_issue_2368(self):
        from numba import jit

        def lift_issue2368(a, b):
            s = 0
            for e in a:
                s += e
            h = b.__hash__()
            return s, h

        a = np.ones(10)
        b = object()
        jitted = jit(forceobj=True)(lift_issue2368)

        expected = lift_issue2368(a, b)
        got = jitted(a, b)

        self.assertEqual(expected[0], got[0])
        self.assertEqual(expected[1], got[1])

        jitloop = jitted.overloads[jitted.signatures[0]].lifted[0]
        [loopcres] = jitloop.overloads.values()
        # assert lifted function is native
        self.assertTrue(loopcres.fndesc.native)

    def test_no_iteration_w_redef(self):
        # redefinition of res in the loop with no use of res should not
        # prevent lifting
        from numba import jit

        @jit(forceobj=True)
        def test(n):
            res = 0
            for i in range(n):
                res = i
            return res

        # loop count = 1, loop lift but loop body not execute
        self.assertEqual(test.py_func(-1), test(-1))
        self.assert_has_lifted(test, loopcount=1)
        # loop count = 1, loop will lift and will execute
        self.assertEqual(test.py_func(1), test(1))
        self.assert_has_lifted(test, loopcount=1)

    def test_no_iteration(self):
        from numba import jit

        @jit(forceobj=True)
        def test(n):
            res = 0
            for i in range(n):
                res += i
            return res

        # loop count = 1
        self.assertEqual(test.py_func(-1), test(-1))
        self.assert_has_lifted(test, loopcount=1)
        # loop count = 1
        self.assertEqual(test.py_func(1), test(1))
        self.assert_has_lifted(test, loopcount=1)

    def test_define_in_loop_body(self):
        # tests a definition in a loop that leaves the loop is liftable
        from numba import jit

        @jit(forceobj=True)
        def test(n):
            for i in range(n):
                res = i
            return res

        # loop count = 1
        self.assertEqual(test.py_func(1), test(1))
        self.assert_has_lifted(test, loopcount=1)

    def test_invalid_argument(self):
        """Test a problem caused by invalid discovery of loop argument
        when a variable is used afterwards but not before.

        Before the fix, this will result in::

        numba.ir.NotDefinedError: 'i' is not defined
        """
        from numba import jit

        @jit(forceobj=True)
        def test(arg):
            if type(arg) == np.ndarray: # force object mode
                if arg.ndim == 1:
                    result = 0.0
                    j = 0
                    for i in range(arg.shape[0]):
                        pass
                else:
                    raise Exception
            else:
                result = 0.0
                i, j = 0, 0
                return result

        arg = np.arange(10)
        self.assertEqual(test.py_func(arg), test(arg))

    def test_conditionally_defined_in_loop(self):
        from numba import jit
        @jit(forceobj=True)
        def test():
            x = 5
            y = 0
            for i in range(2):
                if i > 0:
                   x = 6
                y += x
            return y, x

        self.assertEqual(test.py_func(), test())
        self.assert_has_lifted(test, loopcount=1)

    def test_stack_offset_error_when_has_no_return(self):
        from numba import jit
        import warnings

        def pyfunc(a):
            if a:
                for i in range(10):
                    pass

        with warnings.catch_warnings():
            warnings.simplefilter("error")

            cfunc = jit(forceobj=True)(pyfunc)
            self.assertEqual(pyfunc(True), cfunc(True))

    def test_variable_scope_bug(self):
        """
        https://github.com/numba/numba/issues/2179

        Looplifting transformation is using the wrong version of variable `h`.
        """
        from numba import jit

        def bar(x):
            return x

        def foo(x):
            h = 0.
            for k in range(x):
                h = h + k
            h = h - bar(x)
            return h

        cfoo = jit(forceobj=True)(foo)
        self.assertEqual(foo(10), cfoo(10))

    def test_recompilation_loop(self):
        """
        https://github.com/numba/numba/issues/2481
        """
        from numba import jit

        def foo(x, y):
            # slicing to make array `x` into different layout
            # to cause a new compilation of the lifted loop
            A = x[::y]
            c = 1
            for k in range(A.size):
                object()  # to force objectmode and looplifting
                c = c * A[::-1][k]   # the slice that is failing in static_getitem
            return c

        cfoo = jit(forceobj=True)(foo)
        # First run just works
        args = np.arange(10), 1
        self.assertEqual(foo(*args), cfoo(*args))
        # Exactly 1 lifted loop so far
        self.assertEqual(len(cfoo.overloads[cfoo.signatures[0]].lifted), 1)
        lifted = cfoo.overloads[cfoo.signatures[0]].lifted[0]
        # The lifted loop has 1 signature
        self.assertEqual(len(lifted.signatures), 1)
        # Use different argument to trigger a new compilation of the lifted loop
        args = np.arange(10), -1
        self.assertEqual(foo(*args), cfoo(*args))
        # Ensure that is really a new overload for the lifted loop
        self.assertEqual(len(lifted.signatures), 2)


    def test_lift_objectmode_issue_4223(self):
        from numba import jit

        @jit(forceobj=True)
        def foo(a, b, c, d, x0, y0, n):
            xs, ys = np.zeros(n), np.zeros(n)
            xs[0], ys[0] = x0, y0
            for i in np.arange(n-1):
                xs[i+1] = np.sin(a * ys[i]) + c * np.cos(a * xs[i])
                ys[i+1] = np.sin(b * xs[i]) + d * np.cos(b * ys[i])
            object() # ensure object mode
            return xs, ys

        kwargs = dict(a=1.7, b=1.7, c=0.6, d=1.2, x0=0, y0=0, n=200)
        got = foo(**kwargs)
        expected = foo.py_func(**kwargs)
        self.assertPreciseEqual(got[0], expected[0])
        self .assertPreciseEqual(got[1], expected[1])
        [lifted] = foo.overloads[foo.signatures[0]].lifted
        self.assertEqual(len(lifted.nopython_signatures), 1)

    def test_lift_zip_and_enumerate(self):
        # From issue https://github.com/numba/numba/issues/10076
        from numba import jit

        @jit(forceobj=True)
        def udt_zip(X, Y):
            i = 0
            for x, y in zip(X, Y):
                i = i + x
            return i


        @jit(forceobj=True)
        def udt_enumerate(X, Y):
            i = 0
            for n, x in enumerate(X):
                i = i + x
            return i


        @jit(forceobj=True)
        def udt_enumerate_zip(X, Y):
            i = 0
            for n, (x, y) in enumerate(zip(X, Y)):
                i = i + x
            return i

        X = np.ones(5)
        Y = np.ones(5)

        self.assertEqual(udt_zip(X, Y), udt_zip.py_func(X, Y))
        self.assertEqual(udt_enumerate(X, Y), udt_enumerate.py_func(X, Y))
        self.assertEqual(udt_enumerate_zip(X, Y),
                         udt_enumerate_zip.py_func(X, Y))


if __name__ == '__main__':
    unittest.main()
