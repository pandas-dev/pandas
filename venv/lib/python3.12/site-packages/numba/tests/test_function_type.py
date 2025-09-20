import re
import sys
import unittest
import types as pytypes
from numba import jit, njit, cfunc, types, int64, float64, float32, errors
from numba import literal_unroll, typeof
from numba.core.config import IS_WIN32
import ctypes
import warnings

from .support import TestCase, MemoryLeakMixin
from .support import redirect_c_stderr, captured_stderr

import numpy as np


def dump(foo):  # FOR DEBUGGING, TO BE REMOVED
    from numba.core import function
    foo_type = function.fromobject(foo)
    foo_sig = foo_type.signature()
    foo.compile(foo_sig)
    print('{" LLVM IR OF "+foo.__name__+" ":*^70}')
    print(foo.inspect_llvm(foo_sig.args))
    print('{"":*^70}')


# Decorators for transforming a Python function to different kinds of
# functions:

def mk_cfunc_func(sig):
    def cfunc_func(func):
        assert isinstance(func, pytypes.FunctionType), repr(func)
        f = cfunc(sig)(func)
        f.pyfunc = func
        return f
    return cfunc_func


def njit_func(func):
    assert isinstance(func, pytypes.FunctionType), repr(func)
    f = jit(nopython=True)(func)
    f.pyfunc = func
    return f


def mk_njit_with_sig_func(sig):
    def njit_with_sig_func(func):
        assert isinstance(func, pytypes.FunctionType), repr(func)
        f = jit(sig, nopython=True)(func)
        f.pyfunc = func
        return f
    return njit_with_sig_func


def mk_ctypes_func(sig):
    def ctypes_func(func, sig=int64(int64)):
        assert isinstance(func, pytypes.FunctionType), repr(func)
        cfunc = mk_cfunc_func(sig)(func)
        addr = cfunc._wrapper_address
        if sig == int64(int64):
            f = ctypes.CFUNCTYPE(ctypes.c_int64)(addr)
            f.pyfunc = func
            return f
        raise NotImplementedError(
            f'ctypes decorator for {func} with signature {sig}')
    return ctypes_func


class WAP(types.WrapperAddressProtocol):
    """An example implementation of wrapper address protocol.

    """
    def __init__(self, func, sig):
        self.pyfunc = func
        self.cfunc = cfunc(sig)(func)
        self.sig = sig

    def __wrapper_address__(self):
        return self.cfunc._wrapper_address

    def signature(self):
        return self.sig

    def __call__(self, *args, **kwargs):
        return self.pyfunc(*args, **kwargs)


def mk_wap_func(sig):
    def wap_func(func):
        return WAP(func, sig)
    return wap_func


class TestFunctionType(TestCase):
    """Test first-class functions in the context of a Numba jit compiled
    function.

    """

    def test_in__(self):
        """Function is passed in as an argument.
        """

        def a(i):
            return i + 1

        def foo(f):
            return 0

        sig = int64(int64)

        for decor in [mk_cfunc_func(sig),
                      njit_func,
                      mk_njit_with_sig_func(sig),
                      mk_ctypes_func(sig),
                      mk_wap_func(sig)]:
            for jit_opts in [dict(nopython=True), dict(forceobj=True)]:
                jit_ = jit(**jit_opts)
                with self.subTest(decor=decor.__name__, jit=jit_opts):
                    a_ = decor(a)
                    self.assertEqual(jit_(foo)(a_), foo(a))

    def test_in_call__(self):
        """Function is passed in as an argument and called.
        Also test different return values.
        """

        def a_i64(i):
            return i + 1234567

        def a_f64(i):
            return i + 1.5

        def a_str(i):
            return "abc"

        def foo(f):
            return f(123)

        for f, sig in [(a_i64, int64(int64)), (a_f64, float64(int64))]:
            for decor in [mk_cfunc_func(sig), njit_func,
                          mk_njit_with_sig_func(sig),
                          mk_wap_func(sig)]:
                for jit_opts in [dict(nopython=True), dict(forceobj=True)]:
                    jit_ = jit(**jit_opts)
                    with self.subTest(
                            sig=sig, decor=decor.__name__, jit=jit_opts):
                        f_ = decor(f)
                        self.assertEqual(jit_(foo)(f_), foo(f))

    def test_in_call_out(self):
        """Function is passed in as an argument, called, and returned.
        """

        def a(i):
            return i + 1

        def foo(f):
            f(123)
            return f

        sig = int64(int64)

        for decor in [mk_cfunc_func(sig), njit_func,
                      mk_njit_with_sig_func(sig), mk_wap_func(sig)]:
            for jit_opts in [dict(nopython=True), dict(forceobj=True)]:
                jit_ = jit(**jit_opts)
                with self.subTest(decor=decor.__name__):
                    a_ = decor(a)
                    r1 = jit_(foo)(a_).pyfunc
                    r2 = foo(a)
                    self.assertEqual(r1, r2)

    def test_in_seq_call(self):
        """Functions are passed in as arguments, used as tuple items, and
        called.

        """
        def a(i):
            return i + 1

        def b(i):
            return i + 2

        def foo(f, g):
            r = 0
            for f_ in (f, g):
                r = r + f_(r)
            return r

        sig = int64(int64)

        for decor in [mk_cfunc_func(sig), mk_wap_func(sig),
                      mk_njit_with_sig_func(sig)]:
            for jit_opts in [dict(nopython=True), dict(forceobj=True)]:
                jit_ = jit(**jit_opts)
                with self.subTest(decor=decor.__name__):
                    a_ = decor(a)
                    b_ = decor(b)
                    self.assertEqual(jit_(foo)(a_, b_), foo(a, b))

    def test_in_ns_seq_call(self):
        """Functions are passed in as an argument and via namespace scoping
        (mixed pathways), used as tuple items, and called.

        """

        def a(i):
            return i + 1

        def b(i):
            return i + 2

        def mkfoo(b_):
            def foo(f):
                r = 0
                for f_ in (f, b_):
                    r = r + f_(r)
                return r
            return foo

        sig = int64(int64)

        for decor in [mk_cfunc_func(sig),
                      mk_njit_with_sig_func(sig), mk_wap_func(sig),
                      mk_ctypes_func(sig)][:-1]:
            for jit_opts in [dict(nopython=True), dict(forceobj=True)]:
                jit_ = jit(**jit_opts)
                with self.subTest(decor=decor.__name__):
                    a_ = decor(a)
                    b_ = decor(b)
                    self.assertEqual(jit_(mkfoo(b_))(a_), mkfoo(b)(a))

    def test_ns_call(self):
        """Function is passed in via namespace scoping and called.

        """

        def a(i):
            return i + 1

        def mkfoo(a_):
            def foo():
                return a_(123)
            return foo

        sig = int64(int64)

        for decor in [mk_cfunc_func(sig), njit_func,
                      mk_njit_with_sig_func(sig), mk_wap_func(sig)]:
            for jit_opts in [dict(nopython=True), dict(forceobj=True)]:
                jit_ = jit(**jit_opts)
                with self.subTest(decor=decor.__name__):
                    a_ = decor(a)
                    self.assertEqual(jit_(mkfoo(a_))(), mkfoo(a)())

    def test_ns_out(self):
        """Function is passed in via namespace scoping and returned.

        """
        def a(i):
            return i + 1

        def mkfoo(a_):
            def foo():
                return a_
            return foo

        sig = int64(int64)

        for decor in [mk_cfunc_func(sig), njit_func,
                      mk_njit_with_sig_func(sig), mk_wap_func(sig),
                      mk_ctypes_func(sig)][:-1]:
            for jit_opts in [dict(nopython=True), dict(forceobj=True)]:
                jit_ = jit(**jit_opts)
                with self.subTest(decor=decor.__name__):
                    a_ = decor(a)
                    self.assertEqual(jit_(mkfoo(a_))().pyfunc, mkfoo(a)())

    def test_ns_call_out(self):
        """Function is passed in via namespace scoping, called, and then
        returned.

        """
        def a(i):
            return i + 1

        def mkfoo(a_):
            def foo():
                a_(123)
                return a_
            return foo

        sig = int64(int64)

        for decor in [mk_cfunc_func(sig), njit_func,
                      mk_njit_with_sig_func(sig), mk_wap_func(sig),
                      mk_ctypes_func(sig)]:
            for jit_opts in [dict(nopython=True), dict(forceobj=True)]:
                jit_ = jit(**jit_opts)
            with self.subTest(decor=decor.__name__):
                a_ = decor(a)
                self.assertEqual(jit_(mkfoo(a_))().pyfunc, mkfoo(a)())

    def test_in_overload(self):
        """Function is passed in as an argument and called with different
        argument types.

        """
        def a(i):
            return i + 1

        def foo(f):
            r1 = f(123)
            r2 = f(123.45)
            return (r1, r2)

        for decor in [njit_func]:
            for jit_opts in [dict(nopython=True), dict(forceobj=True)]:
                jit_ = jit(**jit_opts)
                with self.subTest(decor=decor.__name__):
                    a_ = decor(a)
                    self.assertEqual(jit_(foo)(a_), foo(a))

    def test_ns_overload(self):
        """Function is passed in via namespace scoping and called with
        different argument types.

        """
        def a(i):
            return i + 1

        def mkfoo(a_):
            def foo():
                r1 = a_(123)
                r2 = a_(123.45)
                return (r1, r2)
            return foo

        for decor in [njit_func]:
            for jit_opts in [dict(nopython=True), dict(forceobj=True)]:
                jit_ = jit(**jit_opts)
                with self.subTest(decor=decor.__name__):
                    a_ = decor(a)
                    self.assertEqual(jit_(mkfoo(a_))(), mkfoo(a)())

    def test_in_choose(self):
        """Functions are passed in as arguments and called conditionally.

        """
        def a(i):
            return i + 1

        def b(i):
            return i + 2

        def foo(a, b, choose_left):
            if choose_left:
                r = a(1)
            else:
                r = b(2)
            return r

        sig = int64(int64)

        for decor in [mk_cfunc_func(sig), njit_func,
                      mk_njit_with_sig_func(sig), mk_wap_func(sig)]:
            for jit_opts in [dict(nopython=True), dict(forceobj=True)]:
                jit_ = jit(**jit_opts)
                with self.subTest(decor=decor.__name__):
                    a_ = decor(a)
                    b_ = decor(b)
                    self.assertEqual(jit_(foo)(a_, b_, True), foo(a, b, True))
                    self.assertEqual(jit_(foo)(a_, b_, False),
                                     foo(a, b, False))
                    self.assertNotEqual(jit_(foo)(a_, b_, True),
                                        foo(a, b, False))

    def test_ns_choose(self):
        """Functions are passed in via namespace scoping and called
        conditionally.

        """
        def a(i):
            return i + 1

        def b(i):
            return i + 2

        def mkfoo(a_, b_):
            def foo(choose_left):
                if choose_left:
                    r = a_(1)
                else:
                    r = b_(2)
                return r
            return foo

        sig = int64(int64)

        for decor in [mk_cfunc_func(sig), njit_func,
                      mk_njit_with_sig_func(sig), mk_wap_func(sig)]:
            for jit_opts in [dict(nopython=True), dict(forceobj=True)]:
                jit_ = jit(**jit_opts)
                with self.subTest(decor=decor.__name__):
                    a_ = decor(a)
                    b_ = decor(b)
                    self.assertEqual(jit_(mkfoo(a_, b_))(True),
                                     mkfoo(a, b)(True))
                    self.assertEqual(jit_(mkfoo(a_, b_))(False),
                                     mkfoo(a, b)(False))
                    self.assertNotEqual(jit_(mkfoo(a_, b_))(True),
                                        mkfoo(a, b)(False))

    def test_in_choose_out(self):
        """Functions are passed in as arguments and returned conditionally.

        """
        def a(i):
            return i + 1

        def b(i):
            return i + 2

        def foo(a, b, choose_left):
            if choose_left:
                return a
            else:
                return b

        sig = int64(int64)

        for decor in [mk_cfunc_func(sig), njit_func,
                      mk_njit_with_sig_func(sig), mk_wap_func(sig)]:
            for jit_opts in [dict(nopython=True), dict(forceobj=True)]:
                jit_ = jit(**jit_opts)
                with self.subTest(decor=decor.__name__):
                    a_ = decor(a)
                    b_ = decor(b)
                    self.assertEqual(jit_(foo)(a_, b_, True).pyfunc,
                                     foo(a, b, True))
                    self.assertEqual(jit_(foo)(a_, b_, False).pyfunc,
                                     foo(a, b, False))
                    self.assertNotEqual(jit_(foo)(a_, b_, True).pyfunc,
                                        foo(a, b, False))

    def test_in_choose_func_value(self):
        """Functions are passed in as arguments, selected conditionally and
        called.

        """
        def a(i):
            return i + 1

        def b(i):
            return i + 2

        def foo(a, b, choose_left):
            if choose_left:
                f = a
            else:
                f = b
            return f(1)

        sig = int64(int64)

        for decor in [mk_cfunc_func(sig), mk_wap_func(sig), njit_func,
                      mk_njit_with_sig_func(sig)]:
            for jit_opts in [dict(nopython=True), dict(forceobj=True)]:
                jit_ = jit(**jit_opts)
                with self.subTest(decor=decor.__name__):
                    a_ = decor(a)
                    b_ = decor(b)
                    self.assertEqual(jit_(foo)(a_, b_, True), foo(a, b, True))
                    self.assertEqual(jit_(foo)(a_, b_, False),
                                     foo(a, b, False))
                    self.assertNotEqual(jit_(foo)(a_, b_, True),
                                        foo(a, b, False))

    def test_in_pick_func_call(self):
        """Functions are passed in as items of tuple argument, retrieved via
        indexing, and called.

        """
        def a(i):
            return i + 1

        def b(i):
            return i + 2

        def foo(funcs, i):
            f = funcs[i]
            r = f(123)
            return r

        sig = int64(int64)

        for decor in [mk_cfunc_func(sig), mk_wap_func(sig),
                      mk_njit_with_sig_func(sig)]:
            for jit_opts in [dict(nopython=True), dict(forceobj=True)]:
                jit_ = jit(**jit_opts)
                with self.subTest(decor=decor.__name__):
                    a_ = decor(a)
                    b_ = decor(b)
                    self.assertEqual(jit_(foo)((a_, b_), 0), foo((a, b), 0))
                    self.assertEqual(jit_(foo)((a_, b_), 1), foo((a, b), 1))
                    self.assertNotEqual(jit_(foo)((a_, b_), 0), foo((a, b), 1))

    def test_in_iter_func_call(self):
        """Functions are passed in as items of tuple argument, retrieved via
        indexing, and called within a variable for-loop.

        """
        def a(i):
            return i + 1

        def b(i):
            return i + 2

        def foo(funcs, n):
            r = 0
            for i in range(n):
                f = funcs[i]
                r = r + f(r)
            return r

        sig = int64(int64)

        for decor in [mk_cfunc_func(sig), mk_wap_func(sig),
                      mk_njit_with_sig_func(sig)]:
            for jit_opts in [dict(nopython=True), dict(forceobj=True)]:
                jit_ = jit(**jit_opts)
                with self.subTest(decor=decor.__name__):
                    a_ = decor(a)
                    b_ = decor(b)
                    self.assertEqual(jit_(foo)((a_, b_), 2), foo((a, b), 2))

    def test_experimental_feature_warning(self):
        @jit(nopython=True)
        def more(x):
            return x + 1

        @jit(nopython=True)
        def less(x):
            return x - 1

        @jit(nopython=True)
        def foo(sel, x):
            fn = more if sel else less
            return fn(x)

        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter("always")
            res = foo(True, 10)

        self.assertEqual(res, 11)
        self.assertEqual(foo(False, 10), 9)

        self.assertGreaterEqual(len(ws), 1)
        pat = "First-class function type feature is experimental"
        for w in ws:
            if pat in str(w.message):
                break
        else:
            self.fail("missing warning")


class TestFunctionTypeExtensions(TestCase):
    """Test calling external library functions within Numba jit compiled
    functions.

    """

    def test_wrapper_address_protocol_libm(self):
        """Call cos and sinf from standard math library.

        """
        import ctypes.util

        class LibM(types.WrapperAddressProtocol):

            def __init__(self, fname):
                if IS_WIN32:
                    lib = ctypes.cdll.msvcrt
                else:
                    libpath = ctypes.util.find_library('m')
                    lib = ctypes.cdll.LoadLibrary(libpath)
                self.lib = lib
                self._name = fname
                if fname == 'cos':
                    # test for double-precision math function
                    addr = ctypes.cast(self.lib.cos, ctypes.c_voidp).value
                    signature = float64(float64)
                elif fname == 'sinf':
                    # test for single-precision math function
                    # Other 32/64 bit platforms define sinf as the
                    # single-precision sin function
                    addr = ctypes.cast(self.lib.sinf, ctypes.c_voidp).value
                    signature = float32(float32)
                else:
                    raise NotImplementedError(
                        f'wrapper address of `{fname}`'
                        f' with signature `{signature}`')
                self._signature = signature
                self._address = addr

            def __repr__(self):
                return f'{type(self).__name__}({self._name!r})'

            def __wrapper_address__(self):
                return self._address

            def signature(self):
                return self._signature

        mycos = LibM('cos')
        mysin = LibM('sinf')

        def myeval(f, x):
            return f(x)

        # Not testing forceobj=True as it requires implementing
        # LibM.__call__ using ctypes which would be out-of-scope here.
        for jit_opts in [dict(nopython=True)]:
            jit_ = jit(**jit_opts)
            with self.subTest(jit=jit_opts):
                if mycos.signature() is not None:
                    self.assertEqual(jit_(myeval)(mycos, 0.0), 1.0)
                if mysin.signature() is not None:
                    self.assertEqual(jit_(myeval)(mysin, float32(0.0)), 0.0)

    def test_compilation_results(self):
        """Turn the existing compilation results of a dispatcher instance to
        first-class functions with precise types.
        """

        @jit(nopython=True)
        def add_template(x, y):
            return x + y

        # Trigger compilations
        self.assertEqual(add_template(1, 2), 3)
        self.assertEqual(add_template(1.2, 3.4), 4.6)

        cres1, cres2 = add_template.overloads.values()

        # Turn compilation results into first-class functions
        iadd = types.CompileResultWAP(cres1)
        fadd = types.CompileResultWAP(cres2)

        @jit(nopython=True)
        def foo(add, x, y):
            return add(x, y)

        @jit(forceobj=True)
        def foo_obj(add, x, y):
            return add(x, y)

        self.assertEqual(foo(iadd, 3, 4), 7)
        self.assertEqual(foo(fadd, 3.4, 4.5), 7.9)

        self.assertEqual(foo_obj(iadd, 3, 4), 7)
        self.assertEqual(foo_obj(fadd, 3.4, 4.5), 7.9)


class TestMiscIssues(TestCase):
    """Test issues of using first-class functions in the context of Numba
    jit compiled functions.

    """

    def test_issue_3405_using_cfunc(self):

        @cfunc('int64()')
        def a():
            return 2

        @cfunc('int64()')
        def b():
            return 3

        def g(arg):
            if arg:
                f = a
            else:
                f = b
            return f()

        self.assertEqual(jit(nopython=True)(g)(True), 2)
        self.assertEqual(jit(nopython=True)(g)(False), 3)

    def test_issue_3405_using_njit(self):

        @jit(nopython=True)
        def a():
            return 2

        @jit(nopython=True)
        def b():
            return 3

        def g(arg):
            if not arg:
                f = b
            else:
                f = a
            return f()

        self.assertEqual(jit(nopython=True)(g)(True), 2)
        self.assertEqual(jit(nopython=True)(g)(False), 3)

    def test_pr4967_example(self):

        @cfunc('int64(int64)')
        def a(i):
            return i + 1

        @cfunc('int64(int64)')
        def b(i):
            return i + 2

        @jit(nopython=True)
        def foo(f, g):
            i = f(2)
            seq = (f, g)
            for fun in seq:
                i += fun(i)
            return i

        a_ = a._pyfunc
        b_ = b._pyfunc
        self.assertEqual(foo(a, b),
                         a_(2) + a_(a_(2)) + b_(a_(2) + a_(a_(2))))

    def test_pr4967_array(self):
        import numpy as np

        @cfunc("intp(intp[:], float64[:])")
        def foo1(x, y):
            return x[0] + y[0]

        @cfunc("intp(intp[:], float64[:])")
        def foo2(x, y):
            return x[0] - y[0]

        def bar(fx, fy, i):
            a = np.array([10], dtype=np.intp)
            b = np.array([12], dtype=np.float64)
            if i == 0:
                f = fx
            elif i == 1:
                f = fy
            else:
                return
            return f(a, b)

        r = jit(nopython=True, no_cfunc_wrapper=True)(bar)(foo1, foo2, 0)
        self.assertEqual(r, bar(foo1, foo2, 0))
        self.assertNotEqual(r, bar(foo1, foo2, 1))

    def test_reference_example(self):
        import numba

        @numba.njit
        def composition(funcs, x):
            r = x
            for f in funcs[::-1]:
                r = f(r)
            return r

        @numba.cfunc("double(double)")
        def a(x):
            return x + 1.0

        @numba.njit()
        def b(x):
            return x * x

        r = composition((a, b, b, a), 0.5)
        self.assertEqual(r, (0.5 + 1.0) ** 4 + 1.0)

        r = composition((b, a, b, b, a), 0.5)
        self.assertEqual(r, ((0.5 + 1.0) ** 4 + 1.0) ** 2)

    def test_apply_function_in_function(self):

        def foo(f, f_inner):
            return f(f_inner)

        @cfunc('int64(float64)')
        def f_inner(i):
            return int64(i * 3)

        @cfunc(int64(types.FunctionType(f_inner._sig)))
        def f(f_inner):
            return f_inner(123.4)

        self.assertEqual(jit(nopython=True)(foo)(f, f_inner),
                         foo(f._pyfunc, f_inner._pyfunc))

    def test_function_with_none_argument(self):

        @cfunc(int64(types.none))
        def a(i):
            return 1

        @jit(nopython=True)
        def foo(f):
            return f(None)

        self.assertEqual(foo(a), 1)

    def test_constant_functions(self):

        @jit(nopython=True)
        def a():
            return 123

        @jit(nopython=True)
        def b():
            return 456

        @jit(nopython=True)
        def foo():
            return a() + b()

        r = foo()
        if r != 123 + 456:
            print(foo.overloads[()].library.get_llvm_str())
        self.assertEqual(r, 123 + 456)

    def test_generators(self):

        @jit(forceobj=True)
        def gen(xs):
            for x in xs:
                x += 1
                yield x

        @jit(forceobj=True)
        def con(gen_fn, xs):
            return [it for it in gen_fn(xs)]

        self.assertEqual(con(gen, (1, 2, 3)), [2, 3, 4])

        @jit(nopython=True)
        def gen_(xs):
            for x in xs:
                x += 1
                yield x
        self.assertEqual(con(gen_, (1, 2, 3)), [2, 3, 4])

    def test_jit_support(self):

        @jit(nopython=True)
        def foo(f, x):
            return f(x)

        @jit()
        def a(x):
            return x + 1

        @jit()
        def a2(x):
            return x - 1

        @jit()
        def b(x):
            return x + 1.5

        self.assertEqual(foo(a, 1), 2)
        a2(5)  # pre-compile
        self.assertEqual(foo(a2, 2), 1)
        self.assertEqual(foo(a2, 3), 2)
        self.assertEqual(foo(a, 2), 3)
        self.assertEqual(foo(a, 1.5), 2.5)
        self.assertEqual(foo(a2, 1), 0)
        self.assertEqual(foo(a, 2.5), 3.5)
        self.assertEqual(foo(b, 1.5), 3.0)
        self.assertEqual(foo(b, 1), 2.5)

    def test_signature_mismatch(self):
        @jit(nopython=True)
        def f1(x):
            return x

        @jit(nopython=True)
        def f2(x):
            return x

        @jit(nopython=True)
        def foo(disp1, disp2, sel):
            if sel == 1:
                fn = disp1
            else:
                fn = disp2
            return fn([1]), fn(2)

        with self.assertRaises(errors.UnsupportedError) as cm:
            foo(f1, f2, sel=1)
        self.assertRegex(
            str(cm.exception), 'mismatch of function types:')

        # this works because `sel == 1` condition is optimized away:
        self.assertEqual(foo(f1, f1, sel=1), ([1], 2))

    def test_unique_dispatcher(self):
        # In general, the type of a dispatcher instance is imprecise
        # and when used as an input to type-inference, the typing will
        # likely fail. However, if a dispatcher instance contains
        # exactly one overload and compilation is disabled for the dispatcher,
        # then the type of dispatcher instance is interpreted as precise
        # and is transformed to a FunctionType instance with the defined
        # signature of the single overload.

        def foo_template(funcs, x):
            r = x
            for f in funcs:
                r = f(r)
            return r

        a = jit(nopython=True)(lambda x: x + 1)
        b = jit(nopython=True)(lambda x: x + 2)
        foo = jit(nopython=True)(foo_template)

        # compiling and disabling compilation for `a` is sufficient,
        # `b` will inherit its type from the container Tuple type
        a(0)
        a.disable_compile()

        r = foo((a, b), 0)
        self.assertEqual(r, 3)
        # the Tuple type of foo's first argument is a precise FunctionType:
        self.assertEqual(foo.signatures[0][0].dtype.is_precise(), True)

    def test_zero_address(self):

        sig = int64()

        @cfunc(sig)
        def test():
            return 123

        class Good(types.WrapperAddressProtocol):
            """A first-class function type with valid address.
            """

            def __wrapper_address__(self):
                return test.address

            def signature(self):
                return sig

        class Bad(types.WrapperAddressProtocol):
            """A first-class function type with invalid 0 address.
            """

            def __wrapper_address__(self):
                return 0

            def signature(self):
                return sig

        class BadToGood(types.WrapperAddressProtocol):
            """A first-class function type with invalid address that is
            recovered to a valid address.
            """

            counter = -1

            def __wrapper_address__(self):
                self.counter += 1
                return test.address * min(1, self.counter)

            def signature(self):
                return sig

        good = Good()
        bad = Bad()
        bad2good = BadToGood()

        @jit(int64(sig.as_type()))
        def foo(func):
            return func()

        @jit(int64())
        def foo_good():
            return good()

        @jit(int64())
        def foo_bad():
            return bad()

        @jit(int64())
        def foo_bad2good():
            return bad2good()

        self.assertEqual(foo(good), 123)

        self.assertEqual(foo_good(), 123)

        with self.assertRaises(ValueError) as cm:
            foo(bad)
        self.assertRegex(
            str(cm.exception),
            'wrapper address of <.*> instance must be a positive')

        with self.assertRaises(RuntimeError) as cm:
            foo_bad()
        self.assertRegex(
            str(cm.exception), r'.* function address is null')

        self.assertEqual(foo_bad2good(), 123)

    def test_issue_5470(self):

        @njit()
        def foo1():
            return 10

        @njit()
        def foo2():
            return 20

        formulae_foo = (foo1, foo1)

        @njit()
        def bar_scalar(f1, f2):
            return f1() + f2()

        @njit()
        def bar():
            return bar_scalar(*formulae_foo)

        self.assertEqual(bar(), 20)

        formulae_foo = (foo1, foo2)

        @njit()
        def bar():
            return bar_scalar(*formulae_foo)

        self.assertEqual(bar(), 30)

    def test_issue_5540(self):

        @njit(types.int64(types.int64))
        def foo(x):
            return x + 1

        @njit
        def bar_bad(foos):
            f = foos[0]
            return f(x=1)

        @njit
        def bar_good(foos):
            f = foos[0]
            return f(1)

        self.assertEqual(bar_good((foo, )), 2)

        # new/old style error handling
        with self.assertRaises((errors.UnsupportedError,
                                errors.TypingError)) as cm:
            bar_bad((foo, ))

        self.assertRegex(
            str(cm.exception),
            r'.*first-class function call cannot use keyword arguments')

    def test_issue_5615(self):

        @njit
        def foo1(x):
            return x + 1

        @njit
        def foo2(x):
            return x + 2

        @njit
        def bar(fcs):
            x = 0
            a = 10
            i, j = fcs[0]
            x += i(j(a))
            for t in literal_unroll(fcs):
                i, j = t
                x += i(j(a))
            return x

        tup = ((foo1, foo2), (foo2, foo1))

        self.assertEqual(bar(tup), 39)

    def test_issue_5685(self):

        @njit
        def foo1():
            return 1

        @njit
        def foo2(x):
            return x + 1

        @njit
        def foo3(x):
            return x + 2

        @njit
        def bar(fcs):
            r = 0
            for pair in literal_unroll(fcs):
                f1, f2 = pair
                r += f1() + f2(2)
            return r

        self.assertEqual(bar(((foo1, foo2),)), 4)
        self.assertEqual(bar(((foo1, foo2), (foo1, foo3))), 9)  # reproducer


class TestBasicSubtyping(TestCase):
    def test_basic(self):
        """
        Test that a dispatcher object *with* a pre-compiled overload
        can be used as input to another function with locked-down signature
        """
        a = 1

        @njit
        def foo(x):
            return x + 1

        foo(a)
        int_int_fc = types.FunctionType(types.int64(types.int64,))

        @njit(types.int64(int_int_fc))
        def bar(fc):
            return fc(a)

        self.assertEqual(bar(foo), foo(a))

    def test_basic2(self):
        """
        Test that a dispatcher object *without* a pre-compiled overload
        can be used as input to another function with locked-down signature
        """
        a = 1

        @njit
        def foo(x):
            return x + 1

        int_int_fc = types.FunctionType(types.int64(types.int64,))

        @njit(types.int64(int_int_fc))
        def bar(fc):
            return fc(a)

        self.assertEqual(bar(foo), foo(a))

    def test_basic3(self):
        """
        Test that a dispatcher object *without* a pre-compiled overload
        can be used as input to another function with locked-down signature and
        that it behaves as a truly generic function (foo1 does not get locked)
        """
        a = 1

        @njit
        def foo1(x):
            return x + 1

        @njit
        def foo2(x):
            return x + 2

        int_int_fc = types.FunctionType(types.int64(types.int64,))

        @njit(types.int64(int_int_fc))
        def bar(fc):
            return fc(a)

        self.assertEqual(bar(foo1) + 1, bar(foo2))

    def test_basic4(self):
        """
        Test that a dispatcher object can be used as input to another
         function with signature as part of a tuple
        """
        a = 1

        @njit
        def foo1(x):
            return x + 1

        @njit
        def foo2(x):
            return x + 2

        tup = (foo1, foo2)
        int_int_fc = types.FunctionType(types.int64(types.int64,))

        @njit(types.int64(types.UniTuple(int_int_fc, 2)))
        def bar(fcs):
            x = 0
            for i in range(2):
                x += fcs[i](a)
            return x
        self.assertEqual(bar(tup), foo1(a) + foo2(a))

    def test_basic5(self):
        a = 1

        @njit
        def foo1(x):
            return x + 1

        @njit
        def foo2(x):
            return x + 2

        @njit
        def bar1(x):
            return x / 10

        @njit
        def bar2(x):
            return x / 1000

        tup = (foo1, foo2)
        tup_bar = (bar1, bar2)
        int_int_fc = types.FunctionType(types.int64(types.int64,))

        flt_flt_fc = types.FunctionType(types.float64(types.float64,))

        @njit((types.UniTuple(int_int_fc, 2), types.UniTuple(flt_flt_fc, 2)))
        def bar(fcs, ffs):
            x = 0
            for i in range(2):
                x += fcs[i](a)
            for fn in ffs:
                x += fn(a)
            return x

        got = bar(tup, tup_bar)
        expected = foo1(a) + foo2(a) + bar1(a) + bar2(a)
        self.assertEqual(got, expected)


class TestMultiFunctionType(MemoryLeakMixin, TestCase):

    def test_base(self):
        # The test is adapted from https://github.com/numba/numba/issues/9071
        nb_array = typeof(np.ones(2))
        callee_int_type = types.FunctionType(int64(int64))
        sig_int = int64(callee_int_type, int64)
        callee_array_type = types.FunctionType(float64(nb_array))
        sig_array = float64(callee_array_type, nb_array)

        @njit([sig_int, sig_array])
        def caller(callee, a):
            return callee(a)

        @njit
        def callee_int(b):
            return b

        @njit
        def callee_array(c):
            return c.sum()

        b = 1
        c = np.ones(2)

        self.assertEqual(caller(callee_int, b), b)
        self.assertEqual(caller(callee_array, c), c.sum())


class TestInliningFunctionType(MemoryLeakMixin, TestCase):
    def count_num_bb_in_cfg(self, dispatcher):
        dot = dispatcher.inspect_cfg(dispatcher.signatures[0]).dot
        num_of_nodes = re.findall(r"Node0x[0-9a-z]+", dot)
        return len(num_of_nodes)

    def test_inlining_global_dispatcher(self):
        @njit
        def add(x, y):
            return x + y

        fnty = types.FunctionType(int64(int64, int64))

        @njit(int64(fnty, int64, int64))
        def callme(fn, x, y):
            c = 0
            for i in range(100):
                c += fn(x, y)
            return c

        @njit
        def bar(x, y):
            return callme(add, x, y)

        res = bar(123, 321)
        self.assertEqual(100 * (123 + 321), res)
        # There's only one BB because LLVM will be able to fully optimize
        # the reduce-add loop if add() is properly inlined.
        self.assertEqual(self.count_num_bb_in_cfg(bar), 1)

    def test_not_inlining_dispatcher_args(self):
        @njit
        def add(x, y):
            return x + y

        fnty = types.FunctionType(int64(int64, int64))

        @njit(int64(fnty, int64, int64))
        def callme(fn, x, y):
            c = 0
            for i in range(100):
                c += fn(x, y)
            return c

        res = callme(add, 123, 321)

        self.assertEqual(100 * (123 + 321), res)
        # Since add() is not inline-able. The number of BB will be greater
        # than 1. See test_inlining_global_dispatcher().
        self.assertGreater(self.count_num_bb_in_cfg(callme), 1)


class TestExceptionInFunctionType(MemoryLeakMixin, TestCase):
    def test_exception_raising(self):
        class MyError(Exception):
            pass

        @njit
        def add(x, y):
            res = x + y
            if res > 100:
                raise MyError(res)
            return res

        fnty = types.FunctionType(int64(int64, int64))

        @njit(int64(fnty))
        def callme(fn):
            c = 0
            for i in range(100):
                c = fn(c, i)
            return c

        @njit
        def bar():
            return callme(add)

        # Pass Dispatcher as a global reference
        with self.assertRaises(MyError) as exc:
            bar()
        self.assertEqual(exc.exception.args, (105,))

        # Pass Dispatcher by argument
        with self.assertRaises(MyError) as exc:
            callme(add)
        self.assertEqual(exc.exception.args, (105,))

    def test_exception_ignored_in_cfunc(self):
        class MyError(Exception):
            pass

        @njit
        def add(x, y):
            res = x + y
            if res > 100:
                raise MyError(res)
            return res

        fnty = types.FunctionType(int64(int64, int64))

        @njit(int64(fnty, int64, int64))
        def callme(fn, x, y):
            return fn(x, y)

        # Cfunc as argument will ignore raised exception
        @cfunc(int64(int64, int64))
        def c_add(x, y):
            return add(x, y)

        self.assertEqual(callme(c_add, 12, 32), 44)

        # If unittest is buffering (-b), the message goes to Python level stderr
        # otherwise, it goes to C stderr.
        with redirect_c_stderr() as c_stderr, captured_stderr() as stderr:
            # raise ignored and result is garbage
            callme(c_add, 100, 1)
            sys.stderr.flush()

        err = c_stderr.read()
        if not err:
            err = stderr.getvalue()

        self.assertIn("Exception ignored in:", err)
        self.assertIn(str(MyError(101)), err)


if __name__ == '__main__':
    unittest.main()
