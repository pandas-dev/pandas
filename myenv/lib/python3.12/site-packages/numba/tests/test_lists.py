from collections import namedtuple
import contextlib
import itertools
import math
import sys
import ctypes as ct
import numpy as np

from numba import jit, typeof, njit, literal_unroll, literally
import unittest
from numba.core import types, errors
from numba.tests.support import TestCase, MemoryLeakMixin
from numba.experimental import jitclass
from numba.core.extending import overload


Point = namedtuple('Point', ('a', 'b'))


def noop(x):
    pass

def unbox_usecase(x):
    """
    Expect a list of numbers
    """
    res = 0
    for v in x:
        res += v
    return res

def unbox_usecase2(x):
    """
    Expect a list of tuples
    """
    res = 0
    for v in x:
        res += len(v)
    return res

def unbox_usecase3(x):
    """
    Expect a (number, list of numbers) tuple.
    """
    a, b = x
    res = a
    for v in b:
        res += v
    return res

def unbox_usecase4(x):
    """
    Expect a (number, list of tuples) tuple.
    """
    a, b = x
    res = a
    for v in b:
        res += len(v)
    return res


def create_list(x, y, z):
    return [x, y, z]

def create_nested_list(x, y, z, a, b, c):
    return [[x, y, z], [a, b, c]]

def list_comprehension1():
    return sum([x**2 for x in range(10)])

def list_comprehension2():
    return sum([x for x in range(10) if x % 2 == 0])

def list_comprehension3():
    return sum([math.pow(x, 2) for x in range(10)])

def list_comprehension4():
    return sum([x * y for x in range(10) for y in range(10)])

def list_comprehension5():
    return [x * 2 for x in range(10)]

def list_comprehension6():
    return [[x for x in range(y)] for y in range(3)]


def list_constructor(n):
    return list(range(n))

def list_constructor_empty():
    # cannot be typed, list is empty and no typing information is present to
    # infer a type
    return list()

def list_constructor_empty_but_typeable(n):
    # can be typed, list is empty but later append has typing info that allows
    # for inference
    y = list()
    return y.append(n)

def list_append(n):
    l = []
    l.append(42)
    for i in range(n):
        l.append(i)
    return l

def list_append_heterogeneous(n):
    l = []
    l.append(42.0)
    for i in range(n):
        l.append(i)
    return l

def list_extend(n):
    l = []
    # A non-list iterable and a list
    l.extend(range(n))
    l.extend(l[:-1])
    l.extend(range(n, 0, -1))
    return l

def list_extend_heterogeneous(n):
    l = []
    # Extend with various iterables, including lists, with different types
    l.extend(range(n))
    l.extend(l[:-1])
    l.extend((5, 42))
    l.extend([123.0])
    return l

def list_pop0(n):
    l = list(range(n))
    res = 0
    while len(l) > 0:
        res += len(l) * l.pop()
    return res

def list_pop1(n, i):
    l = list(range(n))
    x = l.pop(i)
    return x, l

def list_len(n):
    l = list(range(n))
    return len(l)

def list_getitem(n):
    l = list(range(n))
    res = 0
    # Positive indices
    for i in range(len(l)):
        res += i * l[i]
    # Negative indices
    for i in range(-len(l), 0):
        res -= i * l[i]
    return res

def list_setitem(n):
    l = list(range(n))
    res = 0
    # Positive indices
    for i in range(len(l)):
        l[i] = i * l[i]
    # Negative indices
    for i in range(-len(l), 0):
        l[i] = i * l[i]
    for i in range(len(l)):
        res += l[i]
    return res

def list_getslice2(n, start, stop):
    l = list(range(n))
    return l[start:stop]

def list_getslice3(n, start, stop, step):
    l = list(range(n))
    return l[start:stop:step]

def list_setslice2(n, n_source, start, stop):
    # Generic setslice with size change
    l = list(range(n))
    v = list(range(100, 100 + n_source))
    l[start:stop] = v
    return l

def list_setslice3(n, start, stop, step):
    l = list(range(n))
    v = l[start:stop:step]
    for i in range(len(v)):
        v[i] += 100
    l[start:stop:step] = v
    return l

def list_setslice3_arbitrary(n, n_src, start, stop, step):
    l = list(range(n))
    l[start:stop:step] = list(range(100, 100 + n_src))
    return l

def list_delslice0(n):
    l = list(range(n))
    del l[:]
    return l

def list_delslice1(n, start, stop):
    l = list(range(n))
    del l[start:]
    del l[:stop]
    return l

def list_delslice2(n, start, stop):
    l = list(range(n))
    del l[start:stop]
    return l

def list_clear(n):
    l = list(range(n))
    l.clear()
    return l

def list_copy(n):
    l = list(range(n))
    ll = l.copy()
    l.append(42)
    return l, ll

def list_iteration(n):
    l = list(range(n))
    res = 0
    for i, v in enumerate(l):
        res += i * v
    return res

def list_contains(n):
    l = list(range(n))
    return (0 in l, 1 in l, n - 1 in l, n in l,
            0 not in l, 1 not in l, n - 1 not in l, n not in l,
            )

def list_index1(n, v):
    l = list(range(n, 0, -1))
    return l.index(v)

def list_index2(n, v, start):
    l = list(range(n, 0, -1))
    return l.index(v, start)

def list_index3(n, v, start, stop):
    l = list(range(n, 0, -1))
    return l.index(v, start, stop)

def list_remove(n, v):
    l = list(range(n - 1, -1, -1))
    l.remove(v)
    return l

def list_insert(n, pos, v):
    l = list(range(0, n))
    l.insert(pos, v)
    return l

def list_count(n, v):
    l = []
    for x in range(n):
        l.append(x & 3)
    return l.count(v)

def list_reverse(n):
    l = list(range(n))
    l.reverse()
    return l

def list_add(m, n):
    a = list(range(0, m))
    b = list(range(100, 100 + n))
    res = a + b
    res.append(42)   # check result is a copy
    return a, b, res

def list_add_heterogeneous():
    a = [1]
    b = [2.0]
    c = a + b
    d = b + a
    # check result is a copy
    a.append(3)
    b.append(4.0)
    return a, b, c, d

def list_add_inplace(m, n):
    a = list(range(0, m))
    b = list(range(100, 100 + n))
    a += b
    return a, b

def list_add_inplace_heterogeneous():
    a = [1]
    b = [2.0]
    a += b
    b += a
    return a, b

def list_mul(n, v):
    a = list(range(n))
    return a * v

def list_mul2(n, v):
    a = list(range(n))
    return v * a

def list_mul_inplace(n, v):
    a = list(range(n))
    a *= v
    return a

def list_bool(n):
    a = list(range(n))
    return bool(a), (True if a else False)

def eq_usecase(a, b):
    return list(a) == list(b)

def ne_usecase(a, b):
    return list(a) != list(b)

def gt_usecase(a, b):
    return list(a) > list(b)

def ge_usecase(a, b):
    return list(a) >= list(b)

def lt_usecase(a, b):
    return list(a) < list(b)

def le_usecase(a, b):
    return list(a) <= list(b)

def identity_usecase(n):
    a = list(range(n))
    b = a
    c = a[:]
    return (a is b), (a is not b), (a is c), (a is not c)

def bool_list_usecase():
    # Exercise getitem, setitem, iteration with bool values (issue #1373)
    l = [False]
    l[0] = True
    x = False
    for v in l:
        x = x ^ v
    return l, x

def reflect_simple(l, ll):
    x = l.pop()
    y = l.pop()
    l[0] = 42.
    l.extend(ll)
    return l, x, y

def reflect_conditional(l, ll):
    # `l` may or may not actually reflect a Python list
    if ll[0]:
        l = [11., 22., 33., 44.]
    x = l.pop()
    y = l.pop()
    l[0] = 42.
    l.extend(ll)
    return l, x, y

def reflect_exception(l):
    l.append(42)
    raise ZeroDivisionError

def reflect_dual(l, ll):
    l.append(ll.pop())
    return l is ll


class TestLists(MemoryLeakMixin, TestCase):

    def test_create_list(self):
        pyfunc = create_list
        cfunc = njit((types.int32, types.int32, types.int32))(pyfunc)
        self.assertEqual(cfunc(1, 2, 3), pyfunc(1, 2, 3))

    def test_create_nested_list(self):
        pyfunc = create_nested_list
        cfunc = njit((types.int32, types.int32, types.int32,
                      types.int32, types.int32, types.int32))(pyfunc)
        self.assertEqual(cfunc(1, 2, 3, 4, 5, 6), pyfunc(1, 2, 3, 4, 5, 6))

    def check_unary_with_size(self, pyfunc, precise=True):
        cfunc = jit(nopython=True)(pyfunc)
        # Use various sizes, to stress the allocation algorithm
        for n in [0, 3, 16, 70, 400]:
            eq = self.assertPreciseEqual if precise else self.assertEqual
            eq(cfunc(n), pyfunc(n))

    def test_constructor(self):
        self.check_unary_with_size(list_constructor)

    def test_constructor_empty(self):
        self.disable_leak_check()
        cfunc = jit(nopython=True)(list_constructor_empty)
        with self.assertRaises(errors.TypingError) as raises:
            cfunc()
        errmsg = str(raises.exception)
        self.assertIn("Cannot infer the type of variable", errmsg)
        self.assertIn("list(undefined)", errmsg)
        # check error message went in
        self.assertIn("For Numba to be able to compile a list", errmsg)

    def test_constructor_empty_but_typeable(self):
        args = [np.int32(1), 10., 1 + 3j, [7], [17., 14.], np.array([10])]
        pyfunc = list_constructor_empty_but_typeable
        for arg in args:
            cfunc = jit(nopython=True)(pyfunc)
            expected = pyfunc(arg)
            got = cfunc(arg)
            self.assertPreciseEqual(got, expected)

    def test_append(self):
        self.check_unary_with_size(list_append)

    def test_append_heterogeneous(self):
        self.check_unary_with_size(list_append_heterogeneous, precise=False)

    def test_extend(self):
        self.check_unary_with_size(list_extend)

    def test_extend_heterogeneous(self):
        self.check_unary_with_size(list_extend_heterogeneous, precise=False)

    def test_pop0(self):
        self.check_unary_with_size(list_pop0)

    def test_pop1(self):
        pyfunc = list_pop1
        cfunc = jit(nopython=True)(pyfunc)
        for n in [5, 40]:
            for i in [0, 1, n - 2, n - 1, -1, -2, -n + 3, -n + 1]:
                expected = pyfunc(n, i)
                self.assertPreciseEqual(cfunc(n, i), expected)

    def test_pop_errors(self):
        # XXX References are leaked when an exception is raised
        self.disable_leak_check()
        cfunc = jit(nopython=True)(list_pop1)
        with self.assertRaises(IndexError) as cm:
            cfunc(0, 5)
        self.assertEqual(str(cm.exception), "pop from empty list")
        with self.assertRaises(IndexError) as cm:
            cfunc(1, 5)
        self.assertEqual(str(cm.exception), "pop index out of range")

    def test_insert(self):
        pyfunc = list_insert
        cfunc = jit(nopython=True)(pyfunc)
        for n in [5, 40]:
            indices = [0, 1, n - 2, n - 1, n + 1, -1, -2, -n + 3, -n - 1]
            for i in indices:
                expected = pyfunc(n, i, 42)
                self.assertPreciseEqual(cfunc(n, i, 42), expected)

    def test_len(self):
        self.check_unary_with_size(list_len)

    def test_getitem(self):
        self.check_unary_with_size(list_getitem)

    def test_setitem(self):
        self.check_unary_with_size(list_setitem)

    def check_slicing2(self, pyfunc):
        cfunc = jit(nopython=True)(pyfunc)
        sizes = [5, 40]
        for n in sizes:
            indices = [0, 1, n - 2, -1, -2, -n + 3, -n - 1, -n]
            for start, stop in itertools.product(indices, indices):
                expected = pyfunc(n, start, stop)
                self.assertPreciseEqual(cfunc(n, start, stop), expected)

    def test_getslice2(self):
        self.check_slicing2(list_getslice2)

    def test_setslice2(self):
        pyfunc = list_setslice2
        cfunc = jit(nopython=True)(pyfunc)
        sizes = [5, 40]
        for n, n_src in itertools.product(sizes, sizes):
            indices = [0, 1, n - 2, -1, -2, -n + 3, -n - 1, -n]
            for start, stop in itertools.product(indices, indices):
                expected = pyfunc(n, n_src, start, stop)
                self.assertPreciseEqual(cfunc(n, n_src, start, stop), expected)

    def test_getslice3(self):
        pyfunc = list_getslice3
        cfunc = jit(nopython=True)(pyfunc)
        for n in [10]:
            indices = [0, 1, n - 2, -1, -2, -n + 3, -n - 1, -n]
            steps = [4, 1, -1, 2, -3]
            for start, stop, step in itertools.product(indices, indices, steps):
                expected = pyfunc(n, start, stop, step)
                self.assertPreciseEqual(cfunc(n, start, stop, step), expected)

    def test_setslice3(self):
        pyfunc = list_setslice3
        cfunc = jit(nopython=True)(pyfunc)
        for n in [10]:
            indices = [0, 1, n - 2, -1, -2, -n + 3, -n - 1, -n]
            steps = [4, 1, -1, 2, -3]
            for start, stop, step in itertools.product(indices, indices, steps):
                expected = pyfunc(n, start, stop, step)
                self.assertPreciseEqual(cfunc(n, start, stop, step), expected)

    def test_setslice3_resize(self):
        # XXX References are leaked when an exception is raised
        self.disable_leak_check()
        pyfunc = list_setslice3_arbitrary
        cfunc = jit(nopython=True)(pyfunc)
        # step == 1 => can resize
        cfunc(5, 10, 0, 2, 1)
        # step != 1 => cannot resize
        with self.assertRaises(ValueError) as cm:
            cfunc(5, 100, 0, 3, 2)
        self.assertIn("cannot resize", str(cm.exception))

    def test_delslice0(self):
        self.check_unary_with_size(list_delslice0)

    def test_delslice1(self):
        self.check_slicing2(list_delslice1)

    def test_delslice2(self):
        self.check_slicing2(list_delslice2)

    def test_invalid_slice(self):
        self.disable_leak_check()
        pyfunc = list_getslice3
        cfunc = jit(nopython=True)(pyfunc)
        with self.assertRaises(ValueError) as cm:
            cfunc(10, 1, 2, 0)
        self.assertEqual(str(cm.exception), "slice step cannot be zero")

    def test_iteration(self):
        self.check_unary_with_size(list_iteration)

    def test_reverse(self):
        self.check_unary_with_size(list_reverse)

    def test_contains(self):
        self.check_unary_with_size(list_contains)

    def check_index_result(self, pyfunc, cfunc, args):
        try:
            expected = pyfunc(*args)
        except ValueError:
            with self.assertRaises(ValueError):
                cfunc(*args)
        else:
            self.assertPreciseEqual(cfunc(*args), expected)

    def test_index1(self):
        self.disable_leak_check()
        pyfunc = list_index1
        cfunc = jit(nopython=True)(pyfunc)
        for v in (0, 1, 5, 10, 99999999):
            self.check_index_result(pyfunc, cfunc, (16, v))

    def test_index2(self):
        self.disable_leak_check()
        pyfunc = list_index2
        cfunc = jit(nopython=True)(pyfunc)
        n = 16
        for v in (0, 1, 5, 10, 99999999):
            indices = [0, 1, n - 2, n - 1, n + 1, -1, -2, -n + 3, -n - 1]
            for start in indices:
                self.check_index_result(pyfunc, cfunc, (16, v, start))

    def test_index3(self):
        self.disable_leak_check()
        pyfunc = list_index3
        cfunc = jit(nopython=True)(pyfunc)
        n = 16
        for v in (0, 1, 5, 10, 99999999):
            indices = [0, 1, n - 2, n - 1, n + 1, -1, -2, -n + 3, -n - 1]
            for start, stop in itertools.product(indices, indices):
                self.check_index_result(pyfunc, cfunc, (16, v, start, stop))

    def test_index_exception1(self):
        pyfunc = list_index3
        cfunc = jit(nopython=True)(pyfunc)
        msg = 'arg "start" must be an Integer.'
        with self.assertRaisesRegex(errors.TypingError, msg):
            cfunc(10, 0, 'invalid', 5)

    def test_index_exception2(self):
        pyfunc = list_index3
        cfunc = jit(nopython=True)(pyfunc)
        msg = 'arg "stop" must be an Integer.'
        with self.assertRaisesRegex(errors.TypingError, msg):
            cfunc(10, 0, 0, 'invalid')

    def test_remove(self):
        pyfunc = list_remove
        cfunc = jit(nopython=True)(pyfunc)
        n = 16
        for v in (0, 1, 5, 15):
            expected = pyfunc(n, v)
            self.assertPreciseEqual(cfunc(n, v), expected)

    def test_remove_error(self):
        self.disable_leak_check()
        pyfunc = list_remove
        cfunc = jit(nopython=True)(pyfunc)
        with self.assertRaises(ValueError) as cm:
            cfunc(10, 42)
        self.assertEqual(str(cm.exception), "list.remove(x): x not in list")

    def test_count(self):
        pyfunc = list_count
        cfunc = jit(nopython=True)(pyfunc)
        for v in range(5):
            self.assertPreciseEqual(cfunc(18, v), pyfunc(18, v))

    def test_clear(self):
        self.check_unary_with_size(list_clear)

    def test_copy(self):
        self.check_unary_with_size(list_copy)

    def check_add(self, pyfunc):
        cfunc = jit(nopython=True)(pyfunc)
        sizes = [0, 3, 50, 300]
        for m, n in itertools.product(sizes, sizes):
            expected = pyfunc(m, n)
            self.assertPreciseEqual(cfunc(m, n), expected)

    def test_add(self):
        self.check_add(list_add)

    def test_add_heterogeneous(self):
        pyfunc = list_add_heterogeneous
        cfunc = jit(nopython=True)(pyfunc)
        expected = pyfunc()
        self.assertEqual(cfunc(), expected)

    def test_add_inplace(self):
        self.check_add(list_add_inplace)

    def test_add_inplace_heterogeneous(self):
        pyfunc = list_add_inplace_heterogeneous
        cfunc = jit(nopython=True)(pyfunc)
        expected = pyfunc()
        self.assertEqual(cfunc(), expected)

    def check_mul(self, pyfunc):
        cfunc = jit(nopython=True)(pyfunc)
        for n in [0, 3, 50, 300]:
            for v in [1, 2, 3, 0, -1, -42]:
                expected = pyfunc(n, v)
                self.assertPreciseEqual(cfunc(n, v), expected)

    def test_mul(self):
        self.check_mul(list_mul)

    def test_mul2(self):
        self.check_mul(list_mul2)

    def test_mul_inplace(self):
        self.check_mul(list_mul_inplace)

    @unittest.skipUnless(sys.maxsize >= 2**32,
                         "need a 64-bit system to test for MemoryError")
    def test_mul_error(self):
        self.disable_leak_check()
        pyfunc = list_mul
        cfunc = jit(nopython=True)(pyfunc)
        # Fail in malloc()
        with self.assertRaises(MemoryError):
            cfunc(1, 2**58)
        if sys.platform.startswith('darwin'):
            libc = ct.CDLL('libc.dylib')
            libc.printf("###Please ignore the above error message i.e. \
can't allocate region. It is in fact the purpose of this test to \
request more memory than can be provided###\n".encode("UTF-8"))
        # Overflow size computation when multiplying by item size
        with self.assertRaises(MemoryError):
            cfunc(1, 2**62)

    def test_bool(self):
        pyfunc = list_bool
        cfunc = jit(nopython=True)(pyfunc)
        for n in [0, 1, 3]:
            expected = pyfunc(n)
            self.assertPreciseEqual(cfunc(n), expected)

    def test_list_passing(self):
        # Check one can pass a list from a Numba function to another
        @jit(nopython=True)
        def inner(lst):
            return len(lst), lst[-1]

        @jit(nopython=True)
        def outer(n):
            l = list(range(n))
            return inner(l)

        self.assertPreciseEqual(outer(5), (5, 4))

    def _test_compare(self, pyfunc):
        def eq(args):
            self.assertIs(cfunc(*args), pyfunc(*args),
                          "mismatch for arguments %s" % (args,))

        cfunc = jit(nopython=True)(pyfunc)
        eq(((1, 2), (1, 2)))
        eq(((1, 2, 3), (1, 2)))
        eq(((1, 2), (1, 2, 3)))
        eq(((1, 2, 4), (1, 2, 3)))
        eq(((1.0, 2.0, 3.0), (1, 2, 3)))
        eq(((1.0, 2.0, 3.5), (1, 2, 3)))

    def test_eq(self):
        self._test_compare(eq_usecase)

    def test_ne(self):
        self._test_compare(ne_usecase)

    def test_le(self):
        self._test_compare(le_usecase)

    def test_lt(self):
        self._test_compare(lt_usecase)

    def test_ge(self):
        self._test_compare(ge_usecase)

    def test_gt(self):
        self._test_compare(gt_usecase)

    def test_identity(self):
        pyfunc = identity_usecase
        cfunc = jit(nopython=True)(pyfunc)
        self.assertPreciseEqual(cfunc(3), pyfunc(3))

    def test_bool_list(self):
        # Check lists of bools compile and run successfully
        pyfunc = bool_list_usecase
        cfunc = jit(nopython=True)(pyfunc)
        self.assertPreciseEqual(cfunc(), pyfunc())


class TestUnboxing(MemoryLeakMixin, TestCase):
    """
    Test unboxing of Python lists into native Numba lists.
    """

    @contextlib.contextmanager
    def assert_type_error(self, msg):
        with self.assertRaises(TypeError) as raises:
            yield
        if msg is not None:
            self.assertRegex(str(raises.exception), msg)

    def check_unary(self, pyfunc):
        cfunc = jit(nopython=True)(pyfunc)
        def check(arg):
            expected = pyfunc(arg)
            got = cfunc(arg)
            self.assertPreciseEqual(got, expected)
        return check

    def test_numbers(self):
        check = self.check_unary(unbox_usecase)
        check([1, 2])
        check([1j, 2.5j])

    def test_tuples(self):
        check = self.check_unary(unbox_usecase2)
        check([(1, 2), (3, 4)])
        check([(1, 2j), (3, 4j)])
        check([(), (), ()])

    def test_list_inside_tuple(self):
        check = self.check_unary(unbox_usecase3)
        check((1, [2, 3, 4]))

    def test_list_of_tuples_inside_tuple(self):
        check = self.check_unary(unbox_usecase4)
        check((1, [(2,), (3,)]))

    def test_errors(self):
        # See #1545 and #1594: error checking should ensure the list is
        # homogeneous
        msg = "can't unbox heterogeneous list"
        pyfunc = noop
        cfunc = jit(nopython=True)(pyfunc)
        lst = [1, 2.5]
        with self.assert_type_error(msg):
            cfunc(lst)
        # The list hasn't been changed (bogus reflecting)
        self.assertEqual(lst, [1, 2.5])
        with self.assert_type_error(msg):
            cfunc([1, 2j])
        # Same when the list is nested in a tuple or namedtuple
        with self.assert_type_error(msg):
            cfunc((1, [1, 2j]))
        with self.assert_type_error(msg):
            cfunc(Point(1, [1, 2j]))
        # Issue #1638: tuples of different size.
        # Note the check is really on the tuple side.
        lst = [(1,), (2, 3)]
        with self.assertRaises(TypeError) as raises:
            cfunc(lst)
        msg = ("can't unbox heterogeneous list: "
                "UniTuple({0} x 1) != UniTuple({0} x 2)")
        self.assertEqual(str(raises.exception), msg.format(types.intp))



class TestListReflection(MemoryLeakMixin, TestCase):
    """
    Test reflection of native Numba lists on Python list objects.
    """

    def check_reflection(self, pyfunc):
        cfunc = jit(nopython=True)(pyfunc)
        samples = [([1., 2., 3., 4.], [0.]),
                   ([1., 2., 3., 4.], [5., 6., 7., 8., 9.]),
                   ]
        for dest, src in samples:
            expected = list(dest)
            got = list(dest)
            pyres = pyfunc(expected, src)
            with self.assertRefCount(got, src):
                cres = cfunc(got, src)
                self.assertPreciseEqual(cres, pyres)
                self.assertPreciseEqual(expected, got)
                self.assertEqual(pyres[0] is expected, cres[0] is got)
                del pyres, cres

    def test_reflect_simple(self):
        self.check_reflection(reflect_simple)

    def test_reflect_conditional(self):
        self.check_reflection(reflect_conditional)

    def test_reflect_exception(self):
        """
        When the function exits with an exception, lists should still be
        reflected.
        """
        pyfunc = reflect_exception
        cfunc = jit(nopython=True)(pyfunc)
        l = [1, 2, 3]
        with self.assertRefCount(l):
            with self.assertRaises(ZeroDivisionError):
                cfunc(l)
            self.assertPreciseEqual(l, [1, 2, 3, 42])

    def test_reflect_same_list(self):
        """
        When the same list object is reflected twice, behaviour should
        be consistent.
        """
        pyfunc = reflect_dual
        cfunc = jit(nopython=True)(pyfunc)
        pylist = [1, 2, 3]
        clist = pylist[:]
        expected = pyfunc(pylist, pylist)
        got = cfunc(clist, clist)
        self.assertPreciseEqual(expected, got)
        self.assertPreciseEqual(pylist, clist)
        self.assertRefCountEqual(pylist, clist)

    def test_reflect_clean(self):
        """
        When the list wasn't mutated, no reflection should take place.
        """
        cfunc = jit(nopython=True)(noop)
        # Use a complex, as Python integers can be cached
        l = [12.5j]
        ids = [id(x) for x in l]
        cfunc(l)
        self.assertEqual([id(x) for x in l], ids)


class ManagedListTestCase(MemoryLeakMixin, TestCase):

    def assert_list_element_precise_equal(self, expect, got):
        self.assertEqual(len(expect), len(got))
        for a, b in zip(expect, got):
            self.assertPreciseEqual(a, b)


class TestListManagedElements(ManagedListTestCase):
    "Test list containing objects that need refct"

    def _check_element_equal(self, pyfunc):
        cfunc = jit(nopython=True)(pyfunc)
        con = [np.arange(3).astype(np.intp), np.arange(5).astype(np.intp)]
        expect = list(con)
        pyfunc(expect)
        got = list(con)
        cfunc(got)
        self.assert_list_element_precise_equal(
            expect=expect, got=got
            )

    def test_reflect_passthru(self):
        def pyfunc(con):
            pass
        self._check_element_equal(pyfunc)

    def test_reflect_appended(self):
        def pyfunc(con):
            con.append(np.arange(10).astype(np.intp))

        self._check_element_equal(pyfunc)

    def test_reflect_setitem(self):
        def pyfunc(con):
            con[1] = np.arange(10)

        self._check_element_equal(pyfunc)

    def test_reflect_popped(self):
        def pyfunc(con):
            con.pop()

        self._check_element_equal(pyfunc)

    def test_reflect_insert(self):
        """make sure list.insert() doesn't crash for refcounted objects (see #7553)
        """

        def pyfunc(con):
            con.insert(1, np.arange(4).astype(np.intp))

        self._check_element_equal(pyfunc)

    def test_append(self):
        def pyfunc():
            con = []
            for i in range(300):
                con.append(np.arange(i, ).astype(np.intp))
            return con

        cfunc = jit(nopython=True)(pyfunc)
        expect = pyfunc()
        got = cfunc()

        self.assert_list_element_precise_equal(
            expect=expect, got=got
            )

    def test_append_noret(self):
        # This test make sure local dtor works
        def pyfunc():
            con = []
            for i in range(300):
                con.append(np.arange(i))
            c = 0.0
            for arr in con:
                c += arr.sum() / (1 + arr.size)
            return c

        cfunc = jit(nopython=True)(pyfunc)
        expect = pyfunc()
        got = cfunc()

        self.assertEqual(expect, got)

    def test_reassign_refct(self):
        def pyfunc():
            con = []
            for i in range(5):
                con.append(np.arange(2))
            con[0] = np.arange(4)
            return con

        cfunc = jit(nopython=True)(pyfunc)
        expect = pyfunc()
        got = cfunc()

        self.assert_list_element_precise_equal(
            expect=expect, got=got
            )

    def test_get_slice(self):
        def pyfunc():
            con = []
            for i in range(5):
                con.append(np.arange(2))
            return con[2:4]

        cfunc = jit(nopython=True)(pyfunc)
        expect = pyfunc()
        got = cfunc()

        self.assert_list_element_precise_equal(
            expect=expect, got=got
            )

    def test_set_slice(self):
        def pyfunc():
            con = []
            for i in range(5):
                con.append(np.arange(2))
            con[1:3] = con[2:4]
            return con

        cfunc = jit(nopython=True)(pyfunc)
        expect = pyfunc()
        got = cfunc()

        self.assert_list_element_precise_equal(
            expect=expect, got=got
            )

    def test_pop(self):
        def pyfunc():
            con = []
            for i in range(20):
                con.append(np.arange(i + 1))
            while len(con) > 2:
                con.pop()
            return con

        cfunc = jit(nopython=True)(pyfunc)
        expect = pyfunc()
        got = cfunc()

        self.assert_list_element_precise_equal(
            expect=expect, got=got
            )

    def test_pop_loc(self):
        def pyfunc():
            con = []
            for i in range(1000):
                con.append(np.arange(i + 1))
            while len(con) > 2:
                con.pop(1)
            return con

        cfunc = jit(nopython=True)(pyfunc)
        expect = pyfunc()
        got = cfunc()

        self.assert_list_element_precise_equal(
            expect=expect, got=got
            )

    def test_del_range(self):
        def pyfunc():
            con = []
            for i in range(20):
                con.append(np.arange(i + 1))
            del con[3:10]
            return con

        cfunc = jit(nopython=True)(pyfunc)
        expect = pyfunc()
        got = cfunc()

        self.assert_list_element_precise_equal(
            expect=expect, got=got
            )

    def test_list_of_list(self):
        def pyfunc():
            con = []
            for i in range(10):
                con.append([0] * i)
            return con

        cfunc = jit(nopython=True)(pyfunc)
        expect = pyfunc()
        got = cfunc()

        self.assertEqual(expect, got)



def expect_reflection_failure(fn):
    def wrapped(self, *args, **kwargs):
        self.disable_leak_check()
        with self.assertRaises(TypeError) as raises:
            fn(self, *args, **kwargs)
        expect_msg = 'cannot reflect element of reflected container'
        self.assertIn(expect_msg, str(raises.exception))

    return wrapped


class TestListOfList(ManagedListTestCase):

    def compile_and_test(self, pyfunc, *args):
        from copy import deepcopy
        expect_args = deepcopy(args)
        expect = pyfunc(*expect_args)

        njit_args = deepcopy(args)
        cfunc = jit(nopython=True)(pyfunc)
        got = cfunc(*njit_args)

        self.assert_list_element_precise_equal(
            expect=expect, got=got
            )
        # Check reflection
        self.assert_list_element_precise_equal(
            expect=expect_args, got=njit_args
            )

    def test_returning_list_of_list(self):
        def pyfunc():
            a = [[np.arange(i)] for i in range(4)]
            return a

        self.compile_and_test(pyfunc)

    @expect_reflection_failure
    def test_heterogeneous_list_error(self):
        def pyfunc(x):
            return x[1]

        cfunc = jit(nopython=True)(pyfunc)
        l2 = [[np.zeros(i) for i in range(5)],
              [np.ones(i)+1j for i in range(5)]]
        l3 = [[np.zeros(i) for i in range(5)], [(1,)]]
        l4 = [[1], [{1}]]
        l5 = [[1], [{'a': 1}]]

        # TODO: this triggers a reflection error.
        # Remove this line when nested reflection is supported
        cfunc(l2)

        # error_cases
        with self.assertRaises(TypeError) as raises:
            cfunc(l2)

        self.assertIn(
            ("reflected list(array(float64, 1d, C)) != "
             "reflected list(array(complex128, 1d, C))"),
            str(raises.exception)
            )

        with self.assertRaises(TypeError) as raises:
            cfunc(l3)

        self.assertIn(
            ("reflected list(array(float64, 1d, C)) != "
             "reflected list((int64 x 1))"),
            str(raises.exception)
            )

        with self.assertRaises(TypeError) as raises:
            cfunc(l4)
        self.assertIn(
            "reflected list(int64) != reflected list(reflected set(int64))",
            str(raises.exception)
            )

        with self.assertRaises(ValueError) as raises:
            cfunc(l5)
        self.assertIn(
            "Cannot type list element of <class 'dict'>",
            str(raises.exception)
            )

    @expect_reflection_failure
    def test_list_of_list_reflected(self):
        def pyfunc(l1, l2):
            l1.append(l2)
            l1[-1].append(123)

        cfunc = jit(nopython=True)(pyfunc)
        l1 = [[0, 1], [2, 3]]
        l2 = [4, 5]
        expect = list(l1), list(l2)
        got = list(l1), list(l2)
        pyfunc(*expect)
        cfunc(*got)
        self.assertEqual(expect, got)

    @expect_reflection_failure
    def test_heterogeneous_list(self):
        def pyfunc(x):
            return x[1]

        l1 = [[np.zeros(i) for i in range(5)], [np.ones(i) for i in range(5)]]

        cfunc = jit(nopython=True)(pyfunc)
        l1_got = cfunc(l1)
        self.assertPreciseEqual(pyfunc(l1), l1_got)

    @expect_reflection_failure
    def test_c01(self):
        def bar(x):
            return x.pop()

        r = [[np.zeros(0)], [np.zeros(10)*1j]]
        # TODO: this triggers a reflection error.
        # Remove this line when nested reflection is supported
        self.compile_and_test(bar, r)

        with self.assertRaises(TypeError) as raises:
            self.compile_and_test(bar, r)
        self.assertIn(
            ("reflected list(array(float64, 1d, C)) != "
             "reflected list(array(complex128, 1d, C))"),
            str(raises.exception),
            )

    def test_c02(self):
        def bar(x):
            x.append(x)
            return x

        r = [[np.zeros(0)]]

        with self.assertRaises(errors.TypingError) as raises:
            self.compile_and_test(bar, r)
        self.assertIn(
            "Invalid use of BoundFunction(list.append",
            str(raises.exception),
            )

    def test_c03(self):
        def bar(x):
            f = x
            f[0] = 1
            return f

        r = [[np.arange(3)]]

        with self.assertRaises(errors.TypingError) as raises:
            self.compile_and_test(bar, r)
        self.assertIn(
            "invalid setitem with value of {} to element of {}".format(
                typeof(1),
                typeof(r[0]),
                ),
            str(raises.exception),
        )

    def test_c04(self):
        def bar(x):
            f = x
            f[0][0] = 10
            return f

        r = [[np.arange(3)]]
        with self.assertRaises(errors.TypingError) as raises:
            self.compile_and_test(bar, r)
        self.assertIn(
            "invalid setitem with value of {} to element of {}".format(
                typeof(10),
                typeof(r[0][0]),
                ),
            str(raises.exception),
            )

    @expect_reflection_failure
    def test_c05(self):
        def bar(x):
            f = x
            f[0][0] = np.array([x for x in np.arange(10).astype(np.intp)])
            return f

        r = [[np.arange(3).astype(np.intp)]]
        self.compile_and_test(bar, r)

    def test_c06(self):
        def bar(x):
            f = x
            f[0][0] = np.array([x + 1j for x in np.arange(10)])
            return f

        r = [[np.arange(3)]]
        with self.assertRaises(errors.TypingError) as raises:
            self.compile_and_test(bar, r)
        self.assertIn("invalid setitem with value", str(raises.exception))

    @expect_reflection_failure
    def test_c07(self):
        self.disable_leak_check()

        def bar(x):
            return x[-7]

        r = [[np.arange(3)]]
        cfunc = jit(nopython=True)(bar)
        with self.assertRaises(IndexError) as raises:
            cfunc(r)
        self.assertIn("getitem out of range", str(raises.exception))

    def test_c08(self):
        self.disable_leak_check()

        def bar(x):
            x[5] = 7
            return x

        r = [1, 2, 3]
        cfunc = jit(nopython=True)(bar)
        with self.assertRaises(IndexError) as raises:
            cfunc(r)
        self.assertIn("setitem out of range", str(raises.exception))

    def test_c09(self):
        def bar(x):
            x[-2] = 7j
            return x

        r = [1, 2, 3]
        with self.assertRaises(errors.TypingError) as raises:
            self.compile_and_test(bar, r)
        self.assertIn("invalid setitem with value", str(raises.exception))

    @expect_reflection_failure
    def test_c10(self):
        def bar(x):
            x[0], x[1] = x[1], x[0]
            return x

        r = [[1, 2, 3], [4, 5, 6]]
        self.compile_and_test(bar, r)

    @expect_reflection_failure
    def test_c11(self):
        def bar(x):
            x[:] = x[::-1]
            return x

        r = [[1, 2, 3], [4, 5, 6]]
        self.compile_and_test(bar, r)

    def test_c12(self):
        def bar(x):
            del x[-1]
            return x

        r = [x for x in range(10)]
        self.compile_and_test(bar, r)


class Item(object):
    def __init__(self, many, scalar):
        self.many = many
        self.scalar = scalar


class Container(object):
    def __init__(self, n):
        self.data = [[np.arange(i).astype(np.float64)] for i in range(n)]

    def more(self, n):
        for i in range(n):
            self.data.append([np.arange(i).astype(np.float64)])


class TestListAndJitClasses(ManagedListTestCase):
    def make_jitclass_element(self):
        spec = [
            ('many', types.float64[:]),
            ('scalar', types.float64),
        ]
        JCItem = jitclass(spec)(Item)
        return JCItem

    def make_jitclass_container(self):
        spec = {
            'data': types.List(dtype=types.List(types.float64[::1])),
        }
        JCContainer = jitclass(spec)(Container)
        return JCContainer

    def assert_list_element_with_tester(self, tester, expect, got):
        for x, y in zip(expect, got):
            tester(x, y)

    def test_jitclass_instance_elements(self):
        JCItem = self.make_jitclass_element()

        def pyfunc(xs):
            xs[1], xs[0] = xs[0], xs[1]
            return xs

        def eq(x, y):
            self.assertPreciseEqual(x.many, y.many)
            self.assertPreciseEqual(x.scalar, y.scalar)

        cfunc = jit(nopython=True)(pyfunc)

        arg = [JCItem(many=np.random.random(n + 1), scalar=n * 1.2)
               for n in range(5)]

        expect_arg = list(arg)
        got_arg = list(arg)

        expect_res = pyfunc(expect_arg)
        got_res = cfunc(got_arg)

        self.assert_list_element_with_tester(eq, expect_arg, got_arg)
        self.assert_list_element_with_tester(eq, expect_res, got_res)

    def test_jitclass_containing_list(self):
        JCContainer = self.make_jitclass_container()

        expect = Container(n=4)
        got = JCContainer(n=4)
        self.assert_list_element_precise_equal(got.data, expect.data)
        expect.more(3)
        got.more(3)
        self.assert_list_element_precise_equal(got.data, expect.data)


class TestListInitialValues(MemoryLeakMixin, TestCase):
    """Tests that lists carry their initial value if present"""

    def test_homogeneous_and_literal(self):
        def bar(l):
            ...

        @overload(bar)
        def ol_bar(l):
            if l.initial_value is None:
                return lambda l: literally(l)
            self.assertTrue(isinstance(l, types.List))
            self.assertEqual(l.initial_value, [1, 2, 3])
            self.assertEqual(hasattr(l, 'literal_value'), False)
            return lambda l: l

        @njit
        def foo():
            # keys and values all have literal representation
            x = [1, 2, 3]
            bar(x)

        foo()

    def test_heterogeneous_but_castable_to_homogeneous(self):
        def bar(l):
            ...

        @overload(bar)
        def ol_bar(l):
            self.assertTrue(isinstance(l, types.List))
            self.assertEqual(l.initial_value, None)
            self.assertEqual(hasattr(l, 'literal_value'), False)
            return lambda l: l

        @njit
        def foo():
            # This list will be typed based on 1j, i.e. complex128
            # as the values are not all literals, there's no "initial_value"
            # available irrespective of whether it's possible to rip this
            # information out of the bytecode.
            x = [1j, 2, 3]
            bar(x)

        foo()

    def test_mutation_not_carried(self):
        def bar(d):
            ...

        @overload(bar)
        def ol_bar(d):
            if d.initial_value is None:
                return lambda d: literally(d)
            self.assertTrue(isinstance(d, types.List))
            self.assertEqual(d.initial_value, [1, 2, 3])
            return lambda d: d

        @njit
        def foo():
            # This list is mutated, check the initial_value carries
            # correctly and is not mutated
            x = [1, 2, 3]
            x.append(4)
            bar(x)

        foo()

    def test_mutation_not_carried_single_function(self):
        # this is another pattern for using literally

        @njit
        def nop(*args):
            pass

        for fn, iv in (nop, None), (literally, [1, 2, 3]):
            @njit
            def baz(x):
                pass

            def bar(z):
                pass

            @overload(bar)
            def ol_bar(z):
                def impl(z):
                    fn(z)
                    baz(z)
                return impl

            @njit
            def foo():
                x = [1, 2, 3]
                bar(x)
                x.append(2)
                return x

            foo()
            # baz should be specialised based on literally being invoked and
            # the literal/unliteral arriving at the call site
            larg = baz.signatures[0][0]
            self.assertEqual(larg.initial_value, iv)

    def test_list_of_list_ctor(self):
        # see issue 6082
        @njit
        def bar(x):
            pass

        @njit
        def foo():
            x = [[1, 2, 3, 4, 5], [1, 2, 3, 4, 6]]
            bar(x)

        foo()
        larg = bar.signatures[0][0]
        self.assertEqual(larg.initial_value, None)
        self.assertEqual(larg.dtype.initial_value, None)


class TestLiteralLists(MemoryLeakMixin, TestCase):

    def test_basic_compile(self):
        @njit
        def foo():
            l = [1, 'a']

        foo()

    def test_literal_value_passthrough(self):

        def bar(x):
            pass

        @overload(bar)
        def ol_bar(x):
            self.assertTrue(isinstance(x, types.LiteralList))
            lv = x.literal_value
            self.assertTrue(isinstance(lv, list))
            self.assertEqual(lv[0], types.literal(1))
            self.assertEqual(lv[1], types.literal('a'))
            self.assertEqual(lv[2], types.Array(types.float64, 1, 'C'))
            self.assertEqual(lv[3], types.List(types.intp, reflected=False,
                                               initial_value=[1, 2, 3]))
            self.assertTrue(isinstance(lv[4], types.LiteralList))
            self.assertEqual(lv[4].literal_value[0], types.literal('cat'))
            self.assertEqual(lv[4].literal_value[1], types.literal(10))
            return lambda x: x

        @njit
        def foo():
            otherhomogeneouslist = [1, 2, 3]
            otherheterogeneouslist = ['cat', 10]
            zeros = np.zeros(5)
            l = [1, 'a', zeros, otherhomogeneouslist, otherheterogeneouslist]
            bar(l)

        foo()

    def test_literal_value_involved_passthrough(self):

        def bar(x):
            pass

        @overload(bar)
        def ol_bar(x):
            self.assertTrue(isinstance(x, types.LiteralStrKeyDict))
            dlv = x.literal_value
            inner_literal = {types.literal('g'): types.literal('h'),
                             types.literal('i'):
                                 types.Array(types.float64, 1, 'C')}
            inner_dict = types.LiteralStrKeyDict(inner_literal)
            outer_literal = {types.literal('a'):
                            types.LiteralList([types.literal(1),
                                               types.literal('a'),
                                               types.DictType(
                                                   types.unicode_type,
                                                   types.intp,
                                                   initial_value={'f': 1}),
                                               inner_dict]),
                        types.literal('b'): types.literal(2),
                        types.literal('c'): types.List(types.complex128,
                                                       reflected=False)}
            def check_same(a, b):
                if (isinstance(a, types.LiteralList) and
                    isinstance(b, types.LiteralList)):
                    for i, j in zip(a.literal_value, b.literal_value):
                        check_same(a.literal_value, b.literal_value)
                elif (isinstance(a, list) and
                     isinstance(b, list)):
                    for i, j in zip(a, b):
                        check_same(i, j)
                elif (isinstance(a, types.LiteralStrKeyDict) and
                      isinstance(b, types.LiteralStrKeyDict)):
                    for (ki, vi), (kj, vj) in zip(a.literal_value.items(),
                                                  b.literal_value.items()):
                        check_same(ki, kj)
                        check_same(vi, vj)
                elif (isinstance(a, dict) and
                      isinstance(b, dict)):
                    for (ki, vi), (kj, vj) in zip(a.items(), b.items()):
                        check_same(ki, kj)
                        check_same(vi, vj)
                else:
                    self.assertEqual(a, b)
            check_same(dlv, outer_literal)
            return lambda x: x

        @njit
        def foo():
            # this stretches what's possible with LiteralStrKeyDict and
            # LiteralList, it's nested and contains typed.Dict and np arrays
            # as well as more constant type things.
            l = {'a': [1, 'a', {'f': 1}, {'g': 'h', 'i': np.zeros(5)}], 'b': 2,
                 'c': [1j, 2j, 3j]}
            bar(l)

        foo()

    def test_mutation_failure(self):

        def staticsetitem():
            l = ['a', 1]
            l[0] = 'b'

        def delitem():
            l = ['a', 1]
            del l[0]

        def append():
            l = ['a', 1]
            l.append(2j)

        def extend():
            l = ['a', 1]
            l.extend([2j, 3j])

        def insert():
            l = ['a', 1]
            l.insert(0, 2j)

        def remove():
            l = ['a', 1]
            l.remove('a')

        def pop():
            l = ['a', 1]
            l.pop()

        def clear():
            l = ['a', 1]
            l.clear()

        def sort():
            l = ['a', 1]
            l.sort()

        def reverse():
            l = ['a', 1]
            l.reverse()

        illegals = (staticsetitem, delitem, append, extend, insert, remove, pop,
                    clear, sort, reverse)

        for test in illegals:
            with self.subTest(test.__name__):
                with self.assertRaises(errors.TypingError) as raises:
                    njit(test)()
                expect = "Cannot mutate a literal list"
                self.assertIn(expect, str(raises.exception))

    def test_count(self):
        @njit
        def foo():
            l = ['a', 1, 'a', 2, 'a', 3, 'b', 4, 'b', 5, 'c']
            r = []
            for x in 'abc':
                r.append(l.count(x))
            return r

        self.assertEqual(foo.py_func(), foo())

    def test_len(self):
        @njit
        def foo():
            l = ['a', 1, 'a', 2, 'a', 3, 'b', 4, 'b', 5, 'c']
            return len(l)

        self.assertEqual(foo.py_func(), foo())

    def test_contains(self):
        @njit
        def foo():
            l = ['a', 1, 'a', 2, 'a', 3, 'b', 4, 'b', 5, 'c']
            r = []
            for x in literal_unroll(('a', 'd', 2, 6)):
                r.append(x in l)
            return r

        self.assertEqual(foo.py_func(), foo())

    def test_getitem(self):

        @njit
        def foo(x):
            l = ['a', 1]
            return l[x]

        with self.assertRaises(errors.TypingError) as raises:
            foo(0)
        expect = "Cannot __getitem__ on a literal list"
        self.assertIn(expect, str(raises.exception))

    def test_staticgetitem(self):

        @njit
        def foo():
            l = ['a', 1]
            return l[0], l[1]

        self.assertEqual(foo.py_func(), foo())

    def test_staticgetitem_slice(self):
        # this is forbidden by typing as there's no way to serialize a list of
        # any kind as required by returning a (static) slice of a LiteralList
        @njit
        def foo():
            l = ['a', 'b', 1]
            return l[:2]

        with self.assertRaises(errors.TypingError) as raises:
            foo()
        expect = "Cannot __getitem__ on a literal list"
        self.assertIn(expect, str(raises.exception))

    def test_setitem(self):

        @njit
        def foo(x):
            l = ['a', 1]
            l[x] = 'b'

        with self.assertRaises(errors.TypingError) as raises:
            foo(0)
        expect = "Cannot mutate a literal list"
        self.assertIn(expect, str(raises.exception))

    def test_unify(self):

        @njit
        def foo(x):
            if x + 1 > 3:
                l = ['a', 1]
            else:
                l = ['b', 2]
            return l[0]

        for x in (-100, 100):
            self.assertEqual(foo.py_func(x), foo(x))

    def test_not_unify(self):

        @njit
        def foo(x):
            if x + 1 > 3:
                l = ['a', 1, 2j]
            else:
                l = ['b', 2]
            return l[0], l[1], l[0], l[1] # defeat py310 inliner

        with self.assertRaises(errors.TypingError) as raises:
            foo(100)
        expect = "Cannot unify LiteralList"
        self.assertIn(expect, str(raises.exception))

    def test_index(self):
        @njit
        def foo():
            l = ['a', 1]
            l.index('a')

        with self.assertRaises(errors.TypingError) as raises:
            foo()
        expect = "list.index is unsupported for literal lists"
        self.assertIn(expect, str(raises.exception))

    def test_copy(self):
        @njit
        def foo():
            l = ['a', 1].copy()
            return l[0], l[1]

        self.assertEqual(foo(), foo.py_func())

    def test_tuple_not_in_mro(self):
        # Related to #6094, make sure that LiteralList does not inherit from
        # types.BaseTuple as this breaks isinstance checks.
        def bar(x):
            pass

        @overload(bar)
        def ol_bar(x):
            self.assertFalse(isinstance(x, types.BaseTuple))
            self.assertTrue(isinstance(x, types.LiteralList))
            return lambda x: ...

        @njit
        def foo():
            l = ['a', 1]
            bar(l)

        foo()


if __name__ == '__main__':
    unittest.main()
