import gc
import weakref

from numba import jit
from numba.core import types
from numba.tests.support import TestCase
import unittest


class Dummy(object):

    def __add__(self, other):
        return other + 5


def global_usecase1(x):
    return x + 1

def global_usecase2():
    return global_obj + 1


class TestFuncLifetime(TestCase):
    """
    Test the lifetime of compiled function objects and their dependencies.
    """

    def get_impl(self, dispatcher):
        """
        Get the single implementation (a C function object) of a dispatcher.
        """
        self.assertEqual(len(dispatcher.overloads), 1)
        cres = list(dispatcher.overloads.values())[0]
        return cres.entry_point

    def check_local_func_lifetime(self, **jitargs):
        def f(x):
            return x + 1

        c_f = jit('int32(int32)', **jitargs)(f)
        self.assertPreciseEqual(c_f(1), 2)

        cfunc = self.get_impl(c_f)

        # Since we can't take a weakref to a C function object
        # (see http://bugs.python.org/issue22116), ensure it's
        # collected by taking a weakref to its __self__ instead
        # (a _dynfunc._Closure object).
        refs = [weakref.ref(obj) for obj in (f, c_f, cfunc.__self__)]
        obj = f = c_f = cfunc = None
        gc.collect()
        self.assertEqual([wr() for wr in refs], [None] * len(refs))

    def test_local_func_lifetime(self):
        self.check_local_func_lifetime(forceobj=True)

    def test_local_func_lifetime_npm(self):
        self.check_local_func_lifetime(nopython=True)

    def check_global_func_lifetime(self, **jitargs):
        c_f = jit(**jitargs)(global_usecase1)
        self.assertPreciseEqual(c_f(1), 2)

        cfunc = self.get_impl(c_f)

        wr = weakref.ref(c_f)
        refs = [weakref.ref(obj) for obj in (c_f, cfunc.__self__)]
        obj = c_f = cfunc = None
        gc.collect()
        self.assertEqual([wr() for wr in refs], [None] * len(refs))

    def test_global_func_lifetime(self):
        self.check_global_func_lifetime(forceobj=True)

    def test_global_func_lifetime_npm(self):
        self.check_global_func_lifetime(nopython=True)

    def check_global_obj_lifetime(self, **jitargs):
        # Since global objects can be recorded for typing purposes,
        # check that they are not kept around after they are removed
        # from the globals.
        global global_obj
        global_obj = Dummy()

        c_f = jit(**jitargs)(global_usecase2)
        self.assertPreciseEqual(c_f(), 6)

        refs = [weakref.ref(obj) for obj in (c_f, global_obj)]
        obj = c_f = global_obj = None
        gc.collect()
        self.assertEqual([wr() for wr in refs], [None] * len(refs))

    def test_global_obj_lifetime(self):
        self.check_global_obj_lifetime(forceobj=True)

    def check_inner_function_lifetime(self, **jitargs):
        """
        When a jitted function calls into another jitted function, check
        that everything is collected as desired.
        """
        def mult_10(a):
            return a * 10

        c_mult_10 = jit('intp(intp)', **jitargs)(mult_10)
        c_mult_10.disable_compile()

        def do_math(x):
            return c_mult_10(x + 4)

        c_do_math = jit('intp(intp)', **jitargs)(do_math)
        c_do_math.disable_compile()

        self.assertEqual(c_do_math(1), 50)

        wrs = [weakref.ref(obj) for obj in
               (mult_10, c_mult_10, do_math, c_do_math,
                self.get_impl(c_mult_10).__self__,
                self.get_impl(c_do_math).__self__,
                )]
        obj = mult_10 = c_mult_10 = do_math = c_do_math = None
        gc.collect()
        self.assertEqual([w() for w in wrs], [None] * len(wrs))

    def test_inner_function_lifetime(self):
        self.check_inner_function_lifetime(forceobj=True)

    def test_inner_function_lifetime_npm(self):
        self.check_inner_function_lifetime(nopython=True)


class TestLifeTimeIssue(TestCase):
    def test_double_free(self):
        from numba import njit
        import numpy as np

        # This is the function that causes the crash

        @njit
        def is_point_in_polygons(point, polygons):
            num_polygons = polygons.shape[0]
            if num_polygons != 0:
                # An extra decref is inserted in this block
                intentionally_unused_variable = polygons[0]
            return 0

        # This function creates some NRT objects for the previous function
        # to corrupt.

        @njit
        def dummy():
            return np.empty(10, dtype=np.int64)

        polygons = np.array([[[0, 1]]])
        points = np.array([[-1.5, 0.5]])
        a = dummy()
        is_point_in_polygons(points[0], polygons)
        b = dummy()
        # Crash happens at second call
        is_point_in_polygons(points[0], polygons)
        c = dummy()


if __name__ == '__main__':
    unittest.main()
