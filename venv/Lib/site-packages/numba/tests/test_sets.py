import unittest

from collections import namedtuple
import contextlib
import itertools
import random
from numba.core.errors import TypingError

import numpy as np

from numba import jit, njit
from numba.tests.support import (TestCase, enable_pyobj_flags, MemoryLeakMixin,
                                 compile_function)


Point = namedtuple('Point', ('a', 'b'))


def _build_set_literal_usecase(code, args):
    code = code % {'initializer': ', '.join(repr(arg) for arg in args)}
    return compile_function('build_set', code, globals())

def set_literal_return_usecase(args):
    code = """if 1:
    def build_set():
        return {%(initializer)s}
    """
    return _build_set_literal_usecase(code, args)

def set_literal_convert_usecase(args):
    code = """if 1:
    def build_set():
        my_set = {%(initializer)s}
        return list(my_set)
    """
    return _build_set_literal_usecase(code, args)


def empty_constructor_usecase():
    s = set()
    s.add(1)
    return len(s)

def constructor_usecase(arg):
    s = set(arg)
    return len(s)

def iterator_usecase(arg):
    s = set(arg)
    l = []
    for v in s:
        l.append(v)
    return l

def update_usecase(a, b, c):
    s = set()
    s.update(a)
    s.update(b)
    s.update(c)
    return list(s)

def bool_usecase(arg):
    # Remove one element to allow for empty sets.
    s = set(arg[1:])
    return bool(s)

def remove_usecase(a, b):
    s = set(a)
    for v in b:
        s.remove(v)
    return list(s)

def discard_usecase(a, b):
    s = set(a)
    for v in b:
        s.discard(v)
    return list(s)

def add_discard_usecase(a, u, v):
    s = set(a)
    for i in range(1000):
        s.add(u)
        s.discard(v)
    return list(s)

def pop_usecase(a):
    s = set(a)
    l = []
    while len(s) > 0:
        l.append(s.pop())
    return l

def contains_usecase(a, b):
    s = set(a)
    l = []
    for v in b:
        l.append(v in s)
    return l

def difference_update_usecase(a, b):
    s = set(a)
    s.difference_update(set(b))
    return list(s)

def intersection_update_usecase(a, b):
    s = set(a)
    s.intersection_update(set(b))
    return list(s)

def symmetric_difference_update_usecase(a, b):
    s = set(a)
    s.symmetric_difference_update(set(b))
    return list(s)

def isdisjoint_usecase(a, b):
    return set(a).isdisjoint(set(b))

def issubset_usecase(a, b):
    return set(a).issubset(set(b))

def issuperset_usecase(a, b):
    return set(a).issuperset(set(b))

def clear_usecase(a):
    s = set(a)
    s.clear()
    return len(s), list(s)

def copy_usecase(a):
    s = set(a)
    ss = s.copy()
    s.pop()
    return len(ss), list(ss)

def copy_usecase_empty(a):
    s = set(a)
    s.clear()
    ss = s.copy()
    s.add(a[0])
    return len(ss), list(ss)

def copy_usecase_deleted(a, b):
    s = set(a)
    s.remove(b)
    ss = s.copy()
    s.pop()
    return len(ss), list(ss)

def difference_usecase(a, b):
    sa = set(a)
    s = sa.difference(set(b))
    return list(s)

def intersection_usecase(a, b):
    sa = set(a)
    s = sa.intersection(set(b))
    return list(s)

def symmetric_difference_usecase(a, b):
    sa = set(a)
    s = sa.symmetric_difference(set(b))
    return list(s)

def union_usecase(a, b):
    sa = set(a)
    s = sa.union(set(b))
    return list(s)

def set_return_usecase(a):
    s = set(a)
    return s


def noop(x):
    pass

def unbox_usecase(x):
    """
    Expect a set of numbers
    """
    res = 0
    for v in x:
        res += v
    return res

def unbox_usecase2(x):
    """
    Expect a set of tuples
    """
    res = 0
    for v in x:
        res += len(v)
    return res

def unbox_usecase3(x):
    """
    Expect a (number, set of numbers) tuple.
    """
    a, b = x
    res = a
    for v in b:
        res += v
    return res

def unbox_usecase4(x):
    """
    Expect a (number, set of tuples) tuple.
    """
    a, b = x
    res = a
    for v in b:
        res += len(v)
    return res


def reflect_simple(sa, sb):
    sa.add(42)
    sa.update(sb)
    return sa, len(sa), len(sb)

def reflect_conditional(sa, sb):
    # `sa` may or may not actually reflect a Python set
    if len(sb) > 1:
        sa = set((11., 22., 33., 44.))
    sa.add(42.)
    sa.update(sb)
    # Combine with a non-reflected set (to check method typing)
    sc = set((55., 66.))
    sa.symmetric_difference_update(sc)
    return sa, len(sa), len(sb)

def reflect_exception(s):
    s.add(42)
    raise ZeroDivisionError

def reflect_dual(sa, sb):
    sa.add(sb.pop())
    return sa is sb


def unique_usecase(src):
    seen = set()
    res = []
    for v in src:
        if v not in seen:
            seen.add(v)
            res.append(v)
    return res


class BaseTest(MemoryLeakMixin, TestCase):

    def setUp(self):
        super(BaseTest, self).setUp()
        self.rnd = random.Random(42)

    def _range(self, stop):
        return np.arange(int(stop)).tolist()

    def _random_choice(self, seq, n):
        """
        Choose *n* possibly duplicate items from sequence.
        """
        l = [self.rnd.choice(list(seq)) for i in range(n)]
        if isinstance(seq, np.ndarray):
            return np.array(l, dtype=seq.dtype)
        else:
            return l

    def duplicates_array(self, n):
        """
        Get a 1d array with many duplicate values.
        """
        a = self._range(np.sqrt(n))
        return self._random_choice(a, n)

    def sparse_array(self, n):
        """
        Get a 1d array with values spread around.
        """
        # Note two calls to sparse_array() should generate reasonable overlap
        a = self._range(n ** 1.3)
        return self._random_choice(a, n)

    def _assert_equal_unordered(self, a, b):
        if isinstance(a, tuple):
            self.assertIsInstance(b, tuple)
            for u, v in zip(a, b):
                self._assert_equal_unordered(u, v)
        elif isinstance(a, list):
            self.assertIsInstance(b, list)
            self.assertPreciseEqual(sorted(a), sorted(b))
        else:
            self.assertPreciseEqual(a, b)

    def unordered_checker(self, pyfunc):
        cfunc = jit(nopython=True)(pyfunc)
        def check(*args):
            expected = pyfunc(*args)
            got = cfunc(*args)
            self._assert_equal_unordered(expected, got)
        return check


class TestSetLiterals(BaseTest):

    def check(self, pyfunc):
        cfunc = njit(pyfunc)
        expected = pyfunc()
        got = cfunc()
        self.assertPreciseEqual(expected, got)
        return got, expected

    def test_build_set(self):
        pyfunc = set_literal_return_usecase((1, 2, 3, 2))
        self.check(pyfunc)

    def test_build_heterogeneous_set(self, flags=enable_pyobj_flags):
        pyfunc = set_literal_return_usecase((1, 2.0, 3j, 2))
        self.check(pyfunc)
        pyfunc = set_literal_return_usecase((2.0, 2))
        got, expected = self.check(pyfunc)
        self.assertIs(type(got.pop()), type(expected.pop()))

    def test_build_set_nopython(self):
        arg = list(self.sparse_array(50))
        pyfunc = set_literal_convert_usecase(arg)
        cfunc = jit(nopython=True)(pyfunc)

        expected = pyfunc()
        got = cfunc()
        self.assertPreciseEqual(sorted(expected), sorted(got))


class TestSets(BaseTest):

    def test_constructor(self):
        pyfunc = empty_constructor_usecase
        cfunc = jit(nopython=True)(pyfunc)
        self.assertPreciseEqual(cfunc(), pyfunc())

        pyfunc = constructor_usecase
        cfunc = jit(nopython=True)(pyfunc)
        def check(arg):
            self.assertPreciseEqual(pyfunc(arg), cfunc(arg))

        check(self.duplicates_array(200))
        check(self.sparse_array(200))

    def test_set_return(self):
        pyfunc = set_return_usecase
        cfunc = jit(nopython=True)(pyfunc)

        arg = self.duplicates_array(200)
        self.assertEqual(cfunc(arg), set(arg))

    def test_iterator(self):
        pyfunc = iterator_usecase
        check = self.unordered_checker(pyfunc)

        check(self.duplicates_array(200))
        check(self.sparse_array(200))

    def test_update(self):
        pyfunc = update_usecase
        check = self.unordered_checker(pyfunc)

        a = self.sparse_array(50)
        b = self.duplicates_array(50)
        c = self.sparse_array(50)
        check(a, b, c)

    def test_remove(self):
        pyfunc = remove_usecase
        check = self.unordered_checker(pyfunc)

        a = self.sparse_array(50)
        b = a[::10]
        check(a, b)

    def test_remove_error(self):
        # References are leaked on exception
        self.disable_leak_check()

        pyfunc = remove_usecase
        cfunc = jit(nopython=True)(pyfunc)

        # ensure that there will be a key error
        items = tuple(set(self.sparse_array(3)))
        a = items[1:]
        b = (items[0],)
        with self.assertRaises(KeyError):
            cfunc(a, b)

    def test_discard(self):
        pyfunc = discard_usecase
        check = self.unordered_checker(pyfunc)

        a = self.sparse_array(50)
        b = self.sparse_array(50)
        check(a, b)

    def test_add_discard(self):
        """
        Check that the insertion logic does not create an infinite lookup
        chain with deleted entries (insertion should happen at the first
        deleted entry, not at the free entry at the end of the chain).
        See issue #1913.
        """
        pyfunc = add_discard_usecase
        check = self.unordered_checker(pyfunc)

        # ensure a and b are different
        a = b = None
        while a == b:
            a, b = self.sparse_array(2)
        check((a,), b, b)

    def test_pop(self):
        pyfunc = pop_usecase
        check = self.unordered_checker(pyfunc)

        check(self.sparse_array(50))

    def test_contains(self):
        pyfunc = contains_usecase
        cfunc = jit(nopython=True)(pyfunc)
        def check(a, b):
            self.assertPreciseEqual(pyfunc(a, b), cfunc(a, b))

        a = self.sparse_array(50)
        b = self.sparse_array(50)
        check(a, b)

    def _test_xxx_update(self, pyfunc):
        check = self.unordered_checker(pyfunc)

        sizes = (1, 50, 500)
        for na, nb in itertools.product(sizes, sizes):
            a = self.sparse_array(na)
            b = self.sparse_array(nb)
            check(a, b)

    def test_difference_update(self):
        self._test_xxx_update(difference_update_usecase)

    def test_intersection_update(self):
        self._test_xxx_update(intersection_update_usecase)

    def test_symmetric_difference_update(self):
        self._test_xxx_update(symmetric_difference_update_usecase)

    def _test_comparator(self, pyfunc):
        cfunc = jit(nopython=True)(pyfunc)
        def check(a, b):
            self.assertPreciseEqual(pyfunc(a, b), cfunc(a, b))

        a, b = map(set, [self.sparse_array(10), self.sparse_array(15)])
        args = [a & b, a - b, a | b, a ^ b]
        args = [tuple(x) for x in args]
        for a, b in itertools.product(args, args):
            check(a, b)

    def test_isdisjoint(self):
        self._test_comparator(isdisjoint_usecase)

    def test_issubset(self):
        self._test_comparator(issubset_usecase)

    def test_issuperset(self):
        self._test_comparator(issuperset_usecase)

    def test_clear(self):
        pyfunc = clear_usecase
        check = self.unordered_checker(pyfunc)

        check(self.sparse_array(50))

    def test_copy(self):
        # Source set doesn't have any deleted entries
        pyfunc = copy_usecase
        check = self.unordered_checker(pyfunc)
        check(self.sparse_array(50))

        pyfunc = copy_usecase_empty
        check = self.unordered_checker(pyfunc)
        a = self.sparse_array(1)
        check(a)

        # Source set has deleted entries
        pyfunc = copy_usecase_deleted
        check = self.unordered_checker(pyfunc)
        check((1, 2, 4, 11), 2)
        a = self.sparse_array(50)
        check(a, a[len(a) // 2])

    def test_bool(self):
        pyfunc = bool_usecase
        check = self.unordered_checker(pyfunc)

        check(self.sparse_array(1))
        check(self.sparse_array(2))

    def _test_set_operator(self, pyfunc):
        check = self.unordered_checker(pyfunc)

        a, b = (1, 2, 4, 11), (2, 3, 5, 11, 42)
        check(a, b)

        sizes = (1, 50, 500)
        for na, nb in itertools.product(sizes, sizes):
            a = self.sparse_array(na)
            b = self.sparse_array(nb)
            check(a, b)

    def make_operator_usecase(self, op):
        code = """if 1:
        def operator_usecase(a, b):
            s = set(a) %(op)s set(b)
            return list(s)
        """ % dict(op=op)
        return compile_function('operator_usecase', code, globals())

    def make_inplace_operator_usecase(self, op):
        code = """if 1:
        def inplace_operator_usecase(a, b):
            sa = set(a)
            sb = set(b)
            sc = sa
            sc %(op)s sb
            return list(sc), list(sa)
        """ % dict(op=op)
        return compile_function('inplace_operator_usecase', code, globals())

    def make_comparison_usecase(self, op):
        code = """if 1:
        def comparison_usecase(a, b):
            return set(a) %(op)s set(b)
        """ % dict(op=op)
        return compile_function('comparison_usecase', code, globals())

    def test_difference(self):
        self._test_set_operator(difference_usecase)

    def test_intersection(self):
        self._test_set_operator(intersection_usecase)

    def test_symmetric_difference(self):
        self._test_set_operator(symmetric_difference_usecase)

    def test_union(self):
        self._test_set_operator(union_usecase)

    def test_and(self):
        self._test_set_operator(self.make_operator_usecase('&'))

    def test_or(self):
        self._test_set_operator(self.make_operator_usecase('|'))

    def test_sub(self):
        self._test_set_operator(self.make_operator_usecase('-'))

    def test_xor(self):
        self._test_set_operator(self.make_operator_usecase('^'))

    def test_eq(self):
        self._test_set_operator(self.make_comparison_usecase('=='))

    def test_ne(self):
        self._test_set_operator(self.make_comparison_usecase('!='))

    def test_le(self):
        self._test_set_operator(self.make_comparison_usecase('<='))

    def test_lt(self):
        self._test_set_operator(self.make_comparison_usecase('<'))

    def test_ge(self):
        self._test_set_operator(self.make_comparison_usecase('>='))

    def test_gt(self):
        self._test_set_operator(self.make_comparison_usecase('>'))

    def test_iand(self):
        self._test_set_operator(self.make_inplace_operator_usecase('&='))

    def test_ior(self):
        self._test_set_operator(self.make_inplace_operator_usecase('|='))

    def test_isub(self):
        self._test_set_operator(self.make_inplace_operator_usecase('-='))

    def test_ixor(self):
        self._test_set_operator(self.make_inplace_operator_usecase('^='))


class TestFloatSets(TestSets):
    """
    Test sets with floating-point keys.
    """
    # Only a few basic tests here, as the sanity of most operations doesn't
    # depend on the key type.

    def _range(self, stop):
        return np.arange(stop, dtype=np.float32) * np.float32(0.1)


class TestTupleSets(TestSets):
    """
    Test sets with tuple keys.
    """
    def _range(self, stop):
        a = np.arange(stop, dtype=np.int64)
        b = a & 0x5555555555555555
        c = (a & 0xaaaaaaaa).astype(np.int32)
        d = ((a >> 32) & 1).astype(np.bool_)
        return list(zip(b, c, d))


class TestUnicodeSets(TestSets):
    """
    Test sets with unicode keys. For the purpose of testing refcounted sets.
    """
    def _range(self, stop):
        return ['A{}'.format(i) for i in range(int(stop))]


class TestSetsInvalidDtype(TestSets):

    def _test_set_operator(self, pyfunc):
        # it is invalid to apply some set operations on
        # sets with different dtype
        cfunc = jit(nopython=True)(pyfunc)

        a = set([1, 2, 4, 11])
        b = set(['a', 'b', 'c'])
        msg = 'All Sets must be of the same type'
        with self.assertRaisesRegex(TypingError, msg):
            cfunc(a, b)


class TestSetsInvalid(TestSets):

    def symmetric_difference_usecase(a, b):
        s = a.symmetric_difference(b)
        return list(s)

    def difference_usecase(a, b):
        s = a.difference(b)
        return list(s)

    def intersection_usecase(a, b):
        s = a.intersection(b)
        return list(s)

    def union_usecase(a, b):
        s = a.union(b)
        return list(s)

    def _test_set_operator(self, pyfunc):
        # it is invalid to apply some set operations on
        # sets with different dtype
        cfunc = jit(nopython=True)(pyfunc)

        a = set([1, 2, 4, 11])
        b = (1, 2, 3)
        msg = 'All arguments must be Sets'
        with self.assertRaisesRegex(TypingError, msg):
            cfunc(a, b)

    def test_difference(self):
        self._test_set_operator(TestSetsInvalid.difference_usecase)

    def test_intersection(self):
        self._test_set_operator(TestSetsInvalid.intersection_usecase)

    def test_symmetric_difference(self):
        self._test_set_operator(TestSetsInvalid.symmetric_difference_usecase)

    def test_union(self):
        self._test_set_operator(TestSetsInvalid.union_usecase)

    def make_operator_usecase(self, op):
        code = """if 1:
        def operator_usecase(a, b):
            s = a %(op)s b
            return list(s)
        """ % dict(op=op)
        return compile_function('operator_usecase', code, globals())

    def make_inplace_operator_usecase(self, op):
        code = """if 1:
        def inplace_operator_usecase(a, b):
            sa = a
            sb = b
            sc = sa
            sc %(op)s sb
            return list(sc), list(sa)
        """ % dict(op=op)
        return compile_function('inplace_operator_usecase', code, globals())

    def make_comparison_usecase(self, op):
        code = """if 1:
        def comparison_usecase(a, b):
            return set(a) %(op)s b
        """ % dict(op=op)
        return compile_function('comparison_usecase', code, globals())


class TestUnboxing(BaseTest):
    """
    Test unboxing of Python sets into native Numba sets.
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
        check(set([1, 2]))
        check(set([1j, 2.5j]))
        # Check allocation and sizing
        check(set(range(100)))

    def test_tuples(self):
        check = self.check_unary(unbox_usecase2)
        check(set([(1, 2), (3, 4)]))
        check(set([(1, 2j), (3, 4j)]))

    def test_set_inside_tuple(self):
        check = self.check_unary(unbox_usecase3)
        check((1, set([2, 3, 4])))

    def test_set_of_tuples_inside_tuple(self):
        check = self.check_unary(unbox_usecase4)
        check((1, set([(2,), (3,)])))

    def test_errors(self):
        # Error checking should ensure the set is homogeneous
        msg = "can't unbox heterogeneous set"
        pyfunc = noop
        cfunc = jit(nopython=True)(pyfunc)
        val = set([1, 2.5])
        with self.assert_type_error(msg):
            cfunc(val)
        # The set hasn't been changed (bogus reflecting)
        self.assertEqual(val, set([1, 2.5]))
        with self.assert_type_error(msg):
            cfunc(set([1, 2j]))
        # Same when the set is nested in a tuple or namedtuple
        with self.assert_type_error(msg):
            cfunc((1, set([1, 2j])))
        with self.assert_type_error(msg):
            cfunc(Point(1, set([1, 2j])))
        # Tuples of different size.
        # Note the check is really on the tuple side.
        lst = set([(1,), (2, 3)])
        # Depending on which tuple is examined first, we could get
        # a IndexError or a ValueError.
        with self.assertRaises((IndexError, ValueError)) as raises:
            cfunc(lst)


class TestSetReflection(BaseTest):
    """
    Test reflection of native Numba sets on Python set objects.
    """

    def check_reflection(self, pyfunc):
        cfunc = jit(nopython=True)(pyfunc)
        samples = [(set([1., 2., 3., 4.]), set([0.])),
                   (set([1., 2., 3., 4.]), set([5., 6., 7., 8., 9.])),
                   ]
        for dest, src in samples:
            expected = set(dest)
            got = set(dest)
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
        When the function exits with an exception, sets should still be
        reflected.
        """
        pyfunc = reflect_exception
        cfunc = jit(nopython=True)(pyfunc)
        s = set([1, 2, 3])
        with self.assertRefCount(s):
            with self.assertRaises(ZeroDivisionError):
                cfunc(s)
            self.assertPreciseEqual(s, set([1, 2, 3, 42]))

    def test_reflect_same_set(self):
        """
        When the same set object is reflected twice, behaviour should
        be consistent.
        """
        pyfunc = reflect_dual
        cfunc = jit(nopython=True)(pyfunc)
        pyset = set([1, 2, 3])
        cset = pyset.copy()
        expected = pyfunc(pyset, pyset)
        got = cfunc(cset, cset)
        self.assertPreciseEqual(expected, got)
        self.assertPreciseEqual(pyset, cset)
        self.assertRefCountEqual(pyset, cset)

    def test_reflect_clean(self):
        """
        When the set wasn't mutated, no reflection should take place.
        """
        cfunc = jit(nopython=True)(noop)
        # Use a complex, as Python integers can be cached
        s = set([12.5j])
        ids = [id(x) for x in s]
        cfunc(s)
        self.assertEqual([id(x) for x in s], ids)


class TestExamples(BaseTest):
    """
    Examples of using sets.
    """

    def test_unique(self):
        pyfunc = unique_usecase
        check = self.unordered_checker(pyfunc)

        check(self.duplicates_array(200))
        check(self.sparse_array(200))

    def test_type_coercion_from_update(self):
        # see issue #6621
        def impl():
            i = np.uint64(1)
            R = set()
            R.update({1, 2, 3})
            R.add(i)
            return R
        check = self.unordered_checker(impl)
        check()


if __name__ == '__main__':
    unittest.main()
