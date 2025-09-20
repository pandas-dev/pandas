import warnings
import dis
from itertools import product

import numpy as np

from numba import njit, typed, objmode, prange
from numba.core import ir_utils, ir
from numba.core.errors import (
    CompilerError, NumbaPerformanceWarning, TypingError,
    UnsupportedBytecodeError,
)
from numba.tests.support import (
    TestCase, unittest, captured_stdout, MemoryLeakMixin,
    skip_parfors_unsupported, skip_unless_scipy, expected_failure_py311,
    expected_failure_py312, expected_failure_py313,
)


class MyError(Exception):
    pass


class TestTryBareExcept(TestCase):
    """Test the following pattern:

        try:
            <body>
        except:
            <handling>
    """
    def test_try_inner_raise(self):
        @njit
        def inner(x):
            if x:
                raise MyError

        @njit
        def udt(x):
            try:
                inner(x)
                return "not raised"
            except:             # noqa: E722
                return "caught"

        self.assertEqual(udt(False), "not raised")
        self.assertEqual(udt(True), "caught")

    def test_try_state_reset(self):
        @njit
        def inner(x):
            if x == 1:
                raise MyError("one")
            elif x == 2:
                raise MyError("two")

        @njit
        def udt(x):
            try:
                inner(x)
                res = "not raised"
            except:             # noqa: E722
                res = "caught"
            if x == 0:
                inner(2)
            return res

        with self.assertRaises(MyError) as raises:
            udt(0)
        self.assertEqual(str(raises.exception), "two")
        self.assertEqual(udt(1), "caught")
        self.assertEqual(udt(-1), "not raised")

    def _multi_inner(self):
        @njit
        def inner(x):
            if x == 1:
                print("call_one")
                raise MyError("one")
            elif x == 2:
                print("call_two")
                raise MyError("two")
            elif x == 3:
                print("call_three")
                raise MyError("three")
            else:
                print("call_other")

        return inner

    def test_nested_try(self):
        inner = self._multi_inner()

        @njit
        def udt(x, y, z):
            try:
                try:
                    print("A")
                    inner(x)
                    print("B")
                except:         # noqa: E722
                    print("C")
                    inner(y)
                    print("D")
            except:             # noqa: E722
                print("E")
                inner(z)
                print("F")

        # case 1
        with self.assertRaises(MyError) as raises:
            with captured_stdout() as stdout:
                udt(1, 2, 3)
        self.assertEqual(
            stdout.getvalue().split(),
            ["A", "call_one", "C", "call_two", "E", "call_three"],
        )
        self.assertEqual(str(raises.exception), "three")

        # case 2
        with captured_stdout() as stdout:
            udt(1, 0, 3)
        self.assertEqual(
            stdout.getvalue().split(),
            ["A", "call_one", "C", "call_other", "D"],
        )

        # case 3
        with captured_stdout() as stdout:
            udt(1, 2, 0)
        self.assertEqual(
            stdout.getvalue().split(),
            ["A", "call_one", "C", "call_two", "E", "call_other", "F"],
        )

    def test_loop_in_try(self):
        inner = self._multi_inner()

        @njit
        def udt(x, n):
            try:
                print("A")
                for i in range(n):
                    print(i)
                    if i == x:
                        inner(i)
            except:             # noqa: E722
                print("B")
            return i

        # case 1
        with captured_stdout() as stdout:
            res = udt(3, 5)
        self.assertEqual(
            stdout.getvalue().split(),
            ["A", "0", "1", "2", "3", "call_three", "B"],
        )
        self.assertEqual(res, 3)

        # case 2
        with captured_stdout() as stdout:
            res = udt(1, 3)
        self.assertEqual(
            stdout.getvalue().split(),
            ["A", "0", "1", "call_one", "B"],
        )
        self.assertEqual(res, 1)

        # case 3
        with captured_stdout() as stdout:
            res = udt(0, 3)
        self.assertEqual(
            stdout.getvalue().split(),
            ["A", "0", "call_other", "1", "2"],
        )
        self.assertEqual(res, 2)

    def test_raise_in_try(self):
        @njit
        def udt(x):
            try:
                print("A")
                if x:
                    raise MyError("my_error")
                print("B")
            except:             # noqa: E722
                print("C")
                return 321
            return 123

        # case 1
        with captured_stdout() as stdout:
            res = udt(True)

        self.assertEqual(
            stdout.getvalue().split(),
            ["A", "C"],
        )
        self.assertEqual(res, 321)

        # case 2
        with captured_stdout() as stdout:
            res = udt(False)

        self.assertEqual(
            stdout.getvalue().split(),
            ["A", "B"],
        )
        self.assertEqual(res, 123)

    def test_recursion(self):
        @njit
        def foo(x):
            if x > 0:
                try:
                    foo(x - 1)
                except:   # noqa: E722
                    print("CAUGHT")
                    return 12
            if x == 1:
                raise ValueError("exception")

        with captured_stdout() as stdout:
            res = foo(10)

        self.assertIsNone(res)
        self.assertEqual(stdout.getvalue().split(), ["CAUGHT",],)

    def test_yield(self):
        @njit
        def foo(x):
            if x > 0:
                try:
                    yield 7
                    raise ValueError("exception")   # never hit
                except Exception:
                    print("CAUGHT")

        @njit
        def bar(z):
            return next(foo(z))

        with captured_stdout() as stdout:
            res = bar(10)

        self.assertEqual(res, 7)
        self.assertEqual(stdout.getvalue().split(), [])

    def test_closure2(self):
        @njit
        def foo(x):
            def bar():
                try:
                    raise ValueError("exception")
                except:  # noqa: E722
                    print("CAUGHT")
                    return 12
            bar()

        with captured_stdout() as stdout:
            foo(10)

        self.assertEqual(stdout.getvalue().split(), ["CAUGHT",],)

    def test_closure3(self):
        @njit
        def foo(x):
            def bar(z):
                try:
                    raise ValueError("exception")
                except:  # noqa: E722
                    print("CAUGHT")
                    return z
            return [x for x in map(bar, [1, 2, 3])]

        with captured_stdout() as stdout:
            res = foo(10)

        self.assertEqual(res, [1, 2, 3])
        self.assertEqual(stdout.getvalue().split(), ["CAUGHT",] * 3,)

    def test_closure4(self):
        @njit
        def foo(x):
            def bar(z):
                if z < 0:
                    raise ValueError("exception")
                return z

            try:
                return [x for x in map(bar, [1, 2, 3, x])]
            except:  # noqa: E722
                print("CAUGHT")

        with captured_stdout() as stdout:
            res = foo(-1)

        self.assertEqual(stdout.getvalue().strip(), "CAUGHT")
        self.assertIsNone(res)

        with captured_stdout() as stdout:
            res = foo(4)
        self.assertEqual(stdout.getvalue(), "")
        self.assertEqual(res, [1, 2, 3, 4])

    @skip_unless_scipy
    def test_real_problem(self):
        @njit
        def foo():
            a = np.zeros((4, 4))
            try:
                chol = np.linalg.cholesky(a)
            except:  # noqa: E722
                print("CAUGHT")
                return chol

        with captured_stdout() as stdout:
            foo()

        self.assertEqual(stdout.getvalue().split(), ["CAUGHT",])

    def test_for_loop(self):
        @njit
        def foo(n):

            for i in range(n):
                try:
                    if i > 5:
                        raise ValueError
                except:  # noqa: E722
                    print("CAUGHT")
            else:
                try:
                    try:
                        try:
                            if i > 5:
                                raise ValueError
                        except:  # noqa: E722
                            print("CAUGHT1")
                            raise ValueError
                    except:  # noqa: E722
                        print("CAUGHT2")
                        raise ValueError
                except:  # noqa: E722
                    print("CAUGHT3")

        with captured_stdout() as stdout:
            foo(10)

        self.assertEqual(
            stdout.getvalue().split(),
            ["CAUGHT",] * 4 + ["CAUGHT%s" % i for i in range(1, 4)],
        )

    def test_try_pass(self):
        @njit
        def foo(x):
            try:
                pass
            except:     # noqa: E722
                pass
            return x

        res = foo(123)
        self.assertEqual(res, 123)

    def test_try_except_reraise(self):
        @njit
        def udt():
            try:
                raise ValueError("ERROR")
            except:    # noqa: E722
                raise

        with self.assertRaises(UnsupportedBytecodeError) as raises:
            udt()
        self.assertIn(
            "The re-raising of an exception is not yet supported.",
            str(raises.exception),
        )


class TestTryExceptCaught(TestCase):
    def test_catch_exception(self):
        @njit
        def udt(x):
            try:
                print("A")
                if x:
                    raise ZeroDivisionError("321")
                print("B")
            except Exception:
                print("C")
            print("D")

        # case 1
        with captured_stdout() as stdout:
            udt(True)

        self.assertEqual(
            stdout.getvalue().split(),
            ["A", "C", "D"],
        )

        # case 2
        with captured_stdout() as stdout:
            udt(False)

        self.assertEqual(
            stdout.getvalue().split(),
            ["A", "B", "D"],
        )

    def test_return_in_catch(self):
        @njit
        def udt(x):
            try:
                print("A")
                if x:
                    raise ZeroDivisionError
                print("B")
                r = 123
            except Exception:
                print("C")
                r = 321
                return r
            print("D")
            return r

        # case 1
        with captured_stdout() as stdout:
            res = udt(True)

        self.assertEqual(
            stdout.getvalue().split(),
            ["A", "C"],
        )
        self.assertEqual(res, 321)

        # case 2
        with captured_stdout() as stdout:
            res = udt(False)

        self.assertEqual(
            stdout.getvalue().split(),
            ["A", "B", "D"],
        )
        self.assertEqual(res, 123)

    def test_save_caught(self):
        @njit
        def udt(x):
            try:
                if x:
                    raise ZeroDivisionError
                r = 123
            except Exception as e:  # noqa: F841
                r = 321
                return r
            return r

        with self.assertRaises(UnsupportedBytecodeError) as raises:
            udt(True)
        self.assertIn(
            "Exception object cannot be stored into variable (e)",
            str(raises.exception)
        )

    def test_try_except_reraise(self):
        @njit
        def udt():
            try:
                raise ValueError("ERROR")
            except Exception:
                raise

        with self.assertRaises(UnsupportedBytecodeError) as raises:
            udt()
        self.assertIn(
            "The re-raising of an exception is not yet supported.",
            str(raises.exception),
        )

    def test_try_except_reraise_chain(self):
        @njit
        def udt():
            try:
                raise ValueError("ERROR")
            except Exception:
                try:
                    raise
                except Exception:
                    raise

        with self.assertRaises(UnsupportedBytecodeError) as raises:
            udt()
        self.assertIn(
            "The re-raising of an exception is not yet supported.",
            str(raises.exception),
        )

    def test_division_operator(self):
        # This test that old-style implementation propagate exception
        # to the exception handler properly.
        @njit
        def udt(y):
            try:
                1 / y
            except Exception:
                return 0xdead
            else:
                return 1 / y

        self.assertEqual(udt(0), 0xdead)
        self.assertEqual(udt(2), 0.5)


class TestTryExceptNested(TestCase):
    "Tests for complicated nesting"

    def check_compare(self, cfunc, pyfunc, *args, **kwargs):
        with captured_stdout() as stdout:
            pyfunc(*args, **kwargs)
        expect = stdout.getvalue()

        with captured_stdout() as stdout:
            cfunc(*args, **kwargs)
        got = stdout.getvalue()
        self.assertEqual(
            expect, got,
            msg="args={} kwargs={}".format(args, kwargs)
        )

    def test_try_except_else(self):
        @njit
        def udt(x, y, z, p):
            print('A')
            if x:
                print('B')
                try:
                    print('C')
                    if y:
                        print('D')
                        raise MyError("y")
                    print('E')
                except Exception: # noqa: F722
                    print('F')
                    try:
                        print('H')
                        try:
                            print('I')
                            if z:
                                print('J')
                                raise MyError('z')
                            print('K')
                        except Exception:
                            print('L')
                        else:
                            print('M')
                    except Exception:
                        print('N')
                    else:
                        print('O')
                    print('P')
                else:
                    print('G')
                print('Q')
            print('R')

        cases = list(product([True, False], repeat=4))
        self.assertTrue(cases)
        for x, y, z, p in cases:
            self.check_compare(
                udt, udt.py_func,
                x=x, y=y, z=z, p=p,
            )

    def test_try_except_finally(self):
        @njit
        def udt(p, q):
            try:
                print('A')
                if p:
                    print('B')
                    raise MyError
                print('C')
            except:             # noqa: E722
                print('D')
            finally:
                try:
                    print('E')
                    if q:
                        print('F')
                        raise MyError
                except Exception:
                    print('G')
                else:
                    print('H')
                finally:
                    print('I')

        cases = list(product([True, False], repeat=2))
        self.assertTrue(cases)
        for p, q in cases:
            self.check_compare(
                udt, udt.py_func,
                p=p, q=q,
            )


class TestTryExceptRefct(MemoryLeakMixin, TestCase):
    def test_list_direct_raise(self):
        @njit
        def udt(n, raise_at):
            lst = typed.List()
            try:
                for i in range(n):
                    if i == raise_at:
                        raise IndexError
                    lst.append(i)
            except Exception:
                return lst
            else:
                return lst

        out = udt(10, raise_at=5)
        self.assertEqual(list(out), list(range(5)))
        out = udt(10, raise_at=10)
        self.assertEqual(list(out), list(range(10)))

    def test_list_indirect_raise(self):
        @njit
        def appender(lst, n, raise_at):
            for i in range(n):
                if i == raise_at:
                    raise IndexError
                lst.append(i)
            return lst

        @njit
        def udt(n, raise_at):
            lst = typed.List()
            lst.append(0xbe11)
            try:
                appender(lst, n, raise_at)
            except Exception:
                return lst
            else:
                return lst

        out = udt(10, raise_at=5)
        self.assertEqual(list(out), [0xbe11] + list(range(5)))
        out = udt(10, raise_at=10)
        self.assertEqual(list(out), [0xbe11] + list(range(10)))

    def test_incompatible_refinement(self):
        @njit
        def udt():
            try:
                lst = typed.List()
                print("A")
                lst.append(0)
                print("B")
                lst.append("fda") # invalid type will cause typing error
                print("C")
                return lst
            except Exception:
                print("D")

        with self.assertRaises(TypingError) as raises:
            udt()
        self.assertRegex(
            str(raises.exception),
            r"Cannot refine type|cannot safely cast unicode_type to int(32|64)"
        )


class TestTryExceptOtherControlFlow(TestCase):
    def test_yield(self):
        @njit
        def udt(n, x):
            for i in range(n):
                try:
                    if i == x:
                        raise ValueError
                    yield i
                except Exception:
                    return

        self.assertEqual(list(udt(10, 5)), list(range(5)))
        self.assertEqual(list(udt(10, 10)), list(range(10)))

    @expected_failure_py311
    @expected_failure_py312
    @expected_failure_py313
    def test_objmode(self):
        @njit
        def udt():
            try:
                with objmode():
                    print(object())
            except Exception:
                return

        with self.assertRaises(CompilerError) as raises:
            udt()
        msg = ("unsupported control flow: with-context contains branches "
               "(i.e. break/return/raise) that can leave the block ")
        self.assertIn(
            msg,
            str(raises.exception),
        )

    @expected_failure_py311
    @expected_failure_py312
    @expected_failure_py313
    def test_objmode_output_type(self):
        def bar(x):
            return np.asarray(list(reversed(x.tolist())))

        @njit
        def test_objmode():
            x = np.arange(5)
            y = np.zeros_like(x)
            try:
                with objmode(y='intp[:]'):  # annotate return type
                    # this region is executed by object-mode.
                    y += bar(x)
            except Exception:
                pass
            return y

        with self.assertRaises(CompilerError) as raises:
            test_objmode()
        msg = ("unsupported control flow: with-context contains branches "
               "(i.e. break/return/raise) that can leave the block ")
        self.assertIn(
            msg,
            str(raises.exception),
        )

    def test_reraise_opcode_unreachable(self):
        # The opcode RERAISE was added in python 3.9, there should be no
        # supported way to actually reach it. This test just checks that an
        # exception is present to deal with if it is reached in a case known
        # to produce this opcode.
        def pyfunc():
            try:
                raise Exception
            except Exception:
                raise ValueError("ERROR")
        for inst in dis.get_instructions(pyfunc):
            if inst.opname == 'RERAISE':
                break
        else:
            self.fail("expected RERAISE opcode not found")
        func_ir = ir_utils.get_ir_of_code({}, pyfunc.__code__)
        found = False
        for lbl, blk in func_ir.blocks.items():
            for stmt in blk.find_insts(ir.StaticRaise):
                # don't worry about guarding this strongly, if the exec_args[0]
                # is a string it'll either be "ERROR" or the guard message
                # saying unreachable has been reached
                msg = "Unreachable condition reached (op code RERAISE executed)"
                if stmt.exc_args and msg in stmt.exc_args[0]:
                    found = True
        if not found:
            self.fail("expected RERAISE unreachable message not found")


@skip_parfors_unsupported
class TestTryExceptParfors(TestCase):

    def test_try_in_prange_reduction(self):
        # The try-except is transformed basically into chains of if-else
        def udt(n):
            c = 0
            for i in prange(n):
                try:
                    c += 1
                except Exception:
                    c += 1
            return c

        args = [10]
        expect = udt(*args)
        self.assertEqual(njit(parallel=False)(udt)(*args), expect)
        self.assertEqual(njit(parallel=True)(udt)(*args), expect)

    def test_try_outside_prange_reduction(self):
        # The try-except is transformed basically into chains of if-else
        def udt(n):
            c = 0
            try:
                for i in prange(n):
                    c += 1
            except Exception:
                return 0xdead
            else:
                return c

        args = [10]
        expect = udt(*args)
        self.assertEqual(njit(parallel=False)(udt)(*args), expect)
        # Parfors transformation didn't happen
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', NumbaPerformanceWarning)
            self.assertEqual(njit(parallel=True)(udt)(*args), expect)
        self.assertEqual(len(w), 1)
        self.assertIn("no transformation for parallel execution was possible",
                      str(w[0]))

    def test_try_in_prange_map(self):
        def udt(arr, x):
            out = arr.copy()
            for i in prange(arr.size):
                try:
                    if i == x:
                        raise ValueError
                    out[i] = arr[i] + i
                except Exception:
                    out[i] = -1
            return out

        args = [np.arange(10), 6]
        expect = udt(*args)
        self.assertPreciseEqual(njit(parallel=False)(udt)(*args), expect)
        self.assertPreciseEqual(njit(parallel=True)(udt)(*args), expect)

    def test_try_outside_prange_map(self):
        def udt(arr, x):
            out = arr.copy()
            try:
                for i in prange(arr.size):
                    if i == x:
                        raise ValueError
                    out[i] = arr[i] + i
            except Exception:
                out[i] = -1
            return out

        args = [np.arange(10), 6]
        expect = udt(*args)
        self.assertPreciseEqual(njit(parallel=False)(udt)(*args), expect)
        self.assertPreciseEqual(njit(parallel=True)(udt)(*args), expect)


if __name__ == '__main__':
    unittest.main()
