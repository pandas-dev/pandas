""" Tests for the compiler components of the Numba typed-list.

The tests here should exercise everything within an `@njit` context.
Importantly, the tests should not return a typed list from within such a
context as this would require code from numba/typed/typedlist.py (this is
tested separately).  Tests in this file build on each other in the order of
writing. For example, the first test, tests the creation, append and len of the
list. These are the barebones to do anything useful with a list. The subsequent
test for getitem assumes makes use of these three operations and therefore
assumes that they work.

"""

from textwrap import dedent

from numba import njit
from numba import int32
from numba.extending import register_jitable
from numba.core import types
from numba.core.errors import TypingError
from numba.tests.support import (TestCase, MemoryLeakMixin, override_config,
                                 forbid_codegen)
from numba.typed import listobject, List


class TestCreateAppendLength(MemoryLeakMixin, TestCase):
    """Test list creation, append and len. """

    def test_list_create(self):
        @njit
        def foo(n):
            l = listobject.new_list(int32)
            for i in range(n):
                l.append(i)
            return len(l)

        for i in (0, 1, 2, 100):
            self.assertEqual(foo(i), i)

    def test_list_create_no_jit(self):
        with override_config('DISABLE_JIT', True):
            with forbid_codegen():
                l = listobject.new_list(int32)
                self.assertEqual(type(l), list)

    def test_nonempty_list_create_no_jit(self):
        # See Issue #6001: https://github.com/numba/numba/issues/6001
        with override_config('DISABLE_JIT', True):
            with forbid_codegen():
                l = List([1, 2, 3])
                self.assertEqual(type(l), list)
                self.assertEqual(l, [1, 2, 3])


class TestBool(MemoryLeakMixin, TestCase):
    """Test list bool."""

    def test_list_bool(self):
        @njit
        def foo(n):
            l = listobject.new_list(int32)
            for i in range(n):
                l.append(i)
            return bool(l)

        for i in (0, 1, 2, 100):
            self.assertEqual(foo(i), i > 0)


class TestAllocation(MemoryLeakMixin, TestCase):

    def test_list_allocation(self):
        @njit
        def foo_kwarg(n):
            l = listobject.new_list(int32, allocated=n)
            return l._allocated()

        for i in range(16):
            self.assertEqual(foo_kwarg(i), i)

        @njit
        def foo_posarg(n):
            l = listobject.new_list(int32, n)
            return l._allocated()
        for i in range(16):
            self.assertEqual(foo_posarg(i), i)

    def test_list_allocation_negative(self):
        @njit
        def foo():
            l = listobject.new_list(int32, -1)
            return l._allocated()

        with self.assertRaises(RuntimeError) as raises:
            self.assertEqual(foo(), -1)
        self.assertIn(
            "expecting *allocated* to be >= 0",
            str(raises.exception),
        )


class TestToFromMeminfo(MemoryLeakMixin, TestCase):

    def test_list_to_from_meminfo(self):
        """
        Exercise listobject.{_as_meminfo, _from_meminfo}
        """

        @njit
        def boxer():
            l = listobject.new_list(int32)
            for i in range(10, 20):
                l.append(i)
            return listobject._as_meminfo(l)

        lsttype = types.ListType(int32)

        @njit
        def unboxer(mi):
            l = listobject._from_meminfo(mi, lsttype)
            return l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9]

        mi = boxer()
        self.assertEqual(mi.refcount, 1)

        received = list(unboxer(mi))
        expected = list(range(10, 20))
        self.assertEqual(received, expected)


class TestGetitem(MemoryLeakMixin, TestCase):
    """Test list getitem. """

    def test_list_getitem_singleton(self):
        @njit
        def foo(n):
            l = listobject.new_list(int32)
            l.append(n)
            return l[0]

        self.assertEqual(foo(0), 0)

    def test_list_getitem_singleton_negtive_index(self):
        @njit
        def foo(n):
            l = listobject.new_list(int32)
            l.append(n)
            return l[-1]

        self.assertEqual(foo(0), 0)

    def test_list_getitem_multiple(self):
        @njit
        def foo(i):
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            return l[i]

        for i,j in ((0, 10), (9, 19), (4, 14), (-5, 15), (-1, 19), (-10, 10)):
            self.assertEqual(foo(i), j)

    def test_list_getitem_empty_index_error(self):
        self.disable_leak_check()

        @njit
        def foo(i):
            l = listobject.new_list(int32)
            return l[i]

        for i in (1, 0, -1):
            with self.assertRaises(IndexError) as raises:
                foo(i)
            self.assertIn(
                "list index out of range",
                str(raises.exception),
            )

    def test_list_getitem_multiple_index_error(self):
        self.disable_leak_check()

        @njit
        def foo(i):
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            return l[i]

        for i in (10, -11):
            with self.assertRaises(IndexError) as raises:
                foo(i)
            self.assertIn(
                "list index out of range",
                str(raises.exception),
            )

    def test_list_getitem_empty_typing_error(self):
        self.disable_leak_check()

        @njit
        def foo(i):
            l = listobject.new_list(int32)
            return l[i]

        for i in "xyz", 1.0, 1j:
            with self.assertRaises(TypingError) as raises:
                foo(i)
            self.assertIn(
                "list indices must be integers or slices",
                str(raises.exception),
            )

    def test_list_getitem_integer_types_as_index(self):

        @njit
        def foo(i):
            l = listobject.new_list(int32)
            l.append(0)
            return l[i]

        # try all signed integers and make sure they are cast
        for t in (types.signed_domain
                  ):
            self.assertEqual(foo((t(0))), 0)

    def test_list_getitem_different_sized_uint_index(self):
        # Checks that the index type cast and ext/trunc to the
        # type of the length is correct, both wraparound and
        # direct index is tested via -1/0.

        for ty in types.unsigned_domain:
            @njit
            def foo():
                l = listobject.new_list(int32)
                l.append(7)
                return l[ty(0)]

            self.assertEqual(foo(), 7)

    def test_list_getitem_different_sized_int_index(self):
        # Checks that the index type cast and ext/trunc to the
        # type of the length is correct, both wraparound and
        # direct index is tested via -1/0.

        for ty in types.signed_domain:
            @njit
            def foo():
                l = listobject.new_list(int32)
                l.append(7)
                return l[ty(0)], l[ty(-1)]

            self.assertEqual(foo(), (7, 7))


class TestGetitemSlice(MemoryLeakMixin, TestCase):
    """Test list getitem when indexing with slices. """

    def test_list_getitem_empty_slice_defaults(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            n = l[:]
            return len(n)

        self.assertEqual(foo(), 0)

    def test_list_getitem_singleton_slice_defaults(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            l.append(0)
            n = l[:]
            return len(n)

        self.assertEqual(foo(), 1)

    def test_list_getitem_multiple_slice_defaults(self):
        @njit
        def foo(i):
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            n = l[:]
            return n[i]

        for i,j in ((0, 10), (9, 19), (4, 14), (-5, 15), (-1, 19), (-10, 10)):
            self.assertEqual(foo(i), j)

    def test_list_getitem_multiple_slice_pos_start(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            n = l[5:]
            return len(n), (n[0], n[1], n[2], n[3], n[4])

        length, items = foo()
        self.assertEqual(length, 5)
        self.assertEqual(items, (15, 16, 17, 18, 19))

    def test_list_getitem_multiple_slice_pos_stop(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            n = l[:5]
            return len(n), (n[0], n[1], n[2], n[3], n[4])

        length, items = foo()
        self.assertEqual(length, 5)
        self.assertEqual(items, (10, 11, 12, 13, 14))

    def test_list_getitem_multiple_slice_pos_start_pos_stop(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            n = l[2:7]
            return len(n), (n[0], n[1], n[2], n[3], n[4])

        length, items = foo()
        self.assertEqual(length, 5)
        self.assertEqual(items, (12, 13, 14, 15, 16))

    def test_list_getitem_multiple_slice_pos_start_pos_stop_pos_step(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            n = l[1:9:2]
            return len(n), (n[0], n[1], n[2], n[3])

        length, items = foo()
        self.assertEqual(length, 4)
        self.assertEqual(items, (11, 13, 15, 17))

    def test_list_getitem_multiple_slice_neg_start(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            n = l[-5:]
            return len(n), (n[0], n[1], n[2], n[3], n[4])

        length, items = foo()
        self.assertEqual(length, 5)
        self.assertEqual(items, (15, 16, 17, 18, 19))

    def test_list_getitem_multiple_slice_neg_stop(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            n = l[:-5]
            return len(n), (n[0], n[1], n[2], n[3], n[4])

        length, items = foo()
        self.assertEqual(length, 5)
        self.assertEqual(items, (10, 11, 12, 13, 14))

    def test_list_getitem_multiple_slice_neg_step(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            n = l[::-2]
            return len(n), (n[0], n[1], n[2], n[3], n[4])

        length, items = foo()
        self.assertEqual(length, 5)
        self.assertEqual(items, (19, 17, 15, 13, 11))

    def test_list_getitem_multiple_slice_pos_start_neg_step(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            n = l[4::-1]
            return len(n), (n[0], n[1], n[2], n[3], n[4])

        length, items = foo()
        self.assertEqual(length, 5)
        self.assertEqual(items, (14, 13, 12, 11, 10))

    def test_list_getitem_multiple_slice_neg_start_neg_step(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            n = l[-6::-1]
            return len(n), (n[0], n[1], n[2], n[3], n[4])

        length, items = foo()
        self.assertEqual(length, 5)
        self.assertEqual(items, (14, 13, 12, 11, 10))

    def test_list_getitem_multiple_slice_pos_stop_neg_step(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            n = l[:4:-1]
            return len(n), (n[0], n[1], n[2], n[3], n[4])

        length, items = foo()
        self.assertEqual(length, 5)
        self.assertEqual(items, (19, 18, 17, 16, 15))

    def test_list_getitem_multiple_slice_neg_stop_neg_step(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            n = l[:-6:-1]
            return len(n), (n[0], n[1], n[2], n[3], n[4])

        length, items = foo()
        self.assertEqual(length, 5)
        self.assertEqual(items, (19, 18, 17, 16, 15))

    def test_list_getitem_multiple_slice_pos_start_pos_stop_neg_step(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            n = l[8:3:-1]
            return len(n), (n[0], n[1], n[2], n[3], n[4])

        length, items = foo()
        self.assertEqual(length, 5)
        self.assertEqual(items, (18, 17, 16, 15, 14))

    def test_list_getitem_multiple_slice_neg_start_neg_stop_neg_step(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            n = l[-2:-7:-1]
            return len(n), (n[0], n[1], n[2], n[3], n[4])

        length, items = foo()
        self.assertEqual(length, 5)
        self.assertEqual(items, (18, 17, 16, 15, 14))

    def test_list_getitem_multiple_slice_start_out_of_range(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            n = l[10:]
            return len(n)

        self.assertEqual(foo(), 0)

    def test_list_getitem_multiple_slice_stop_zero(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            n = l[:0]
            return len(n)

        self.assertEqual(foo(), 0)

    def test_list_getitem_multiple_slice_zero_step_index_error(self):
        self.disable_leak_check()

        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            l[::0]

        with self.assertRaises(ValueError) as raises:
            foo()
        self.assertIn(
            "slice step cannot be zero",
            str(raises.exception),
        )


class TestSetitem(MemoryLeakMixin, TestCase):
    """Test list setitem. """

    def test_list_setitem_singleton(self):
        @njit
        def foo(n):
            l = listobject.new_list(int32)
            l.append(0)
            l[0] = n
            return l[0]

        for i in (0, 1, 2, 100):
            self.assertEqual(foo(i), i)

    def test_list_setitem_singleton_negative_index(self):
        @njit
        def foo(n):
            l = listobject.new_list(int32)
            l.append(0)
            l[0] = n
            return l[-1]

        for i in (0, 1, 2, 100):
            self.assertEqual(foo(i), i)

    def test_list_setitem_singleton_index_error(self):
        self.disable_leak_check()

        @njit
        def foo(i):
            l = listobject.new_list(int32)
            l.append(0)
            l[i] = 1

        with self.assertRaises(IndexError):
            foo(1)

        with self.assertRaises(IndexError):
            foo(-2)

    def test_list_setitem_multiple(self):

        @njit
        def foo(i, n):
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            l[i] = n
            return l[i]

        for i,n in zip(range(0,10), range(20,30)):
            self.assertEqual(foo(i, n), n)

    def test_list_setitem_multiple_index_error(self):
        self.disable_leak_check()

        @njit
        def foo(i):
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            l[i] = 0

        with self.assertRaises(IndexError):
            foo(10)

        with self.assertRaises(IndexError):
            foo(-11)

    def test_list_setitem_singleton_typing_error_on_index(self):
        self.disable_leak_check()

        @njit
        def foo(i):
            l = listobject.new_list(int32)
            l.append(0)
            # slice with a non-{integer,slice}
            l[i] = 1

        for i in "xyz", 1.0, 1j:
            with self.assertRaises(TypingError) as raises:
                foo(i)
            self.assertIn(
                "list indices must be integers or slices",
                str(raises.exception),
            )

    def test_list_setitem_singleton_typing_error_on_item(self):
        self.disable_leak_check()

        @njit
        def foo():
            l = listobject.new_list(int32)
            l.append(0)
            # assign a non-iterable to a slice
            l[:] = 1

        with self.assertRaises(TypingError) as raises:
            foo()
        self.assertIn(
            "can only assign an iterable when using a slice "
            "with assignment/setitem",
            str(raises.exception),
        )

    def test_list_setitem_integer_types_as_index(self):

        @njit
        def foo(i):
            l = listobject.new_list(int32)
            l.append(0)
            l[i] = 1
            return l[i]

        # try all signed integers and make sure they are cast
        for t in (types.signed_domain
                  ):
            self.assertEqual(foo((t(0))), 1)


class TestPop(MemoryLeakMixin, TestCase):
    """Test list pop. """

    def test_list_pop_singleton(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            l.append(0)
            return l.pop(), len(l)

        self.assertEqual(foo(), (0, 0))

    def test_list_pop_singleton_index(self):
        @njit
        def foo(i):
            l = listobject.new_list(int32)
            l.append(0)
            return l.pop(i), len(l)

        self.assertEqual(foo(0), (0, 0))
        self.assertEqual(foo(-1), (0, 0))

    def test_list_pop_multiple(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in (10, 11, 12):
                l.append(j)
            return l.pop(), len(l)

        self.assertEqual(foo(), (12, 2))

    def test_list_pop_multiple_index(self):
        @njit
        def foo(i):
            l = listobject.new_list(int32)
            for j in (10, 11, 12):
                l.append(j)
            return l.pop(i), len(l)

        for i, n in ((0, 10), (1, 11), (2, 12)):
            self.assertEqual(foo(i), (n, 2))

        for i, n in ((-3, 10), (-2, 11), (-1, 12)):
            self.assertEqual(foo(i), (n, 2))

    def test_list_pop_integer_types_as_index(self):

        @njit
        def foo(i):
            l = listobject.new_list(int32)
            l.append(0)
            return l.pop(i)

        # try all signed integers and make sure they are cast
        for t in (types.signed_domain
                  ):
            self.assertEqual(foo((t(0))), 0)

    def test_list_pop_empty_index_error_no_index(self):
        self.disable_leak_check()

        @njit
        def foo():
            l = listobject.new_list(int32)
            l.pop()

        with self.assertRaises(IndexError) as raises:
            foo()
        self.assertIn(
            "pop from empty list",
            str(raises.exception),
        )

    def test_list_pop_empty_index_error_with_index(self):
        self.disable_leak_check()

        @njit
        def foo(i):
            l = listobject.new_list(int32)
            l.pop(i)

        with self.assertRaises(IndexError) as raises:
            foo(-1)
        self.assertIn(
            "pop from empty list",
            str(raises.exception),
        )

        with self.assertRaises(IndexError) as raises:
            foo(0)
        self.assertIn(
            "pop from empty list",
            str(raises.exception),
        )

        with self.assertRaises(IndexError) as raises:
            foo(1)
        self.assertIn(
            "pop from empty list",
            str(raises.exception),
        )

    def test_list_pop_mutiple_index_error_with_index(self):
        self.disable_leak_check()

        @njit
        def foo(i):
            l = listobject.new_list(int32)
            for j in (10, 11, 12):
                l.append(j)
            l.pop(i)

        with self.assertRaises(IndexError) as raises:
            foo(-4)
        self.assertIn(
            "list index out of range",
            str(raises.exception),
        )

        with self.assertRaises(IndexError) as raises:
            foo(3)
        self.assertIn(
            "list index out of range",
            str(raises.exception),
        )

    def test_list_pop_singleton_typing_error_on_index(self):
        self.disable_leak_check()

        @njit
        def foo(i):
            l = listobject.new_list(int32)
            l.append(0)
            # slice with a non-{integer,slice}
            return l.pop(i)

        for i in "xyz", 1.0, 1j:
            with self.assertRaises(TypingError) as raises:
                foo(i)
            self.assertIn(
                "argument for pop must be an integer",
                str(raises.exception),
            )


class TestListObjectDelitem(MemoryLeakMixin, TestCase):
    """Test list delitem.
    """

    def test_list_singleton_delitem_index(self):

        @njit
        def foo():
            l = listobject.new_list(int32)
            l.append(0)
            del l[0]
            return len(l)
        self.assertEqual(foo(), 0)

    def test_list_singleton_delitem_slice_defaults(self):

        @njit
        def foo():
            l = listobject.new_list(int32)
            l.append(0)
            del l[:]
            return len(l)
        self.assertEqual(foo(), 0)

    def test_list_singleton_delitem_slice_start(self):

        @njit
        def foo():
            l = listobject.new_list(int32)
            l.append(0)
            del l[0:]
            return len(l)
        self.assertEqual(foo(), 0)

    def test_list_singleton_delitem_slice_stop(self):

        @njit
        def foo():
            l = listobject.new_list(int32)
            l.append(0)
            del l[:1]
            return len(l)
        self.assertEqual(foo(), 0)

    def test_list_singleton_delitem_slice_start_stop(self):

        @njit
        def foo():
            l = listobject.new_list(int32)
            l.append(0)
            del l[0:1]
            return len(l)
        self.assertEqual(foo(), 0)

    def test_list_singleton_delitem_slice_start_step(self):

        @njit
        def foo():
            l = listobject.new_list(int32)
            l.append(0)
            del l[0::1]
            return len(l)
        self.assertEqual(foo(), 0)

    def test_list_singleton_delitem_slice_start_stop_step(self):

        @njit
        def foo():
            l = listobject.new_list(int32)
            l.append(0)
            del l[0:1:1]
            return len(l)
        self.assertEqual(foo(), 0)

    def test_list_multiple_delitem(self):

        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in (10, 11, 12):
                l.append(j)
            del l[0]
            return len(l), l[0], l[1]
        self.assertEqual(foo(), (2, 11, 12))

    def test_list_multiple_delitem_slice(self):

        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in (10, 11, 12):
                l.append(j)
            del l[:]
            return len(l)
        self.assertEqual(foo(), 0)

    def test_list_multiple_delitem_off_by_one(self):
        # this was exposing a nasty off-by-one error, leaving it in to detect
        # and regressions
        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            k = listobject.new_list(int32)
            for j in range(10, 20):
                k.append(j)
            # should be a no-op
            del l[-9:-20]
            return k == l
        self.assertTrue(foo())


class TestContains(MemoryLeakMixin, TestCase):
    """Test list contains. """

    def test_list_contains_empty(self):
        @njit
        def foo(i):
            l = listobject.new_list(int32)
            return i in l

        self.assertFalse(foo(0))
        self.assertFalse(foo(1))

    def test_list_contains_singleton(self):
        @njit
        def foo(i):
            l = listobject.new_list(int32)
            l.append(0)
            return i in l

        self.assertTrue(foo(0))
        self.assertFalse(foo(1))

    def test_list_contains_multiple(self):
        @njit
        def foo(i):
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            return i in l

        for i in range(10, 20):
            self.assertTrue(foo(i))

        for i in range(20, 30):
            self.assertFalse(foo(i))


class TestCount(MemoryLeakMixin, TestCase):
    """Test list count. """

    def test_list_count_empty(self):
        @njit
        def foo(i):
            l = listobject.new_list(int32)
            return l.count(i)

        self.assertEqual(foo(10), 0)

    def test_list_count_singleton(self):
        @njit
        def foo(i):
            l = listobject.new_list(int32)
            l.append(10)
            return l.count(i)

        self.assertEqual(foo(1), 0)
        self.assertEqual(foo(10), 1)

    def test_list_count_mutiple(self):
        @njit
        def foo(i):
            l = listobject.new_list(int32)
            for j in [11, 12, 12, 13, 13, 13]:
                l.append(j)
            return l.count(i)

        self.assertEqual(foo(10), 0)
        self.assertEqual(foo(11), 1)
        self.assertEqual(foo(12), 2)
        self.assertEqual(foo(13), 3)


class TestExtend(MemoryLeakMixin, TestCase):
    """Test list extend. """

    def test_list_extend_empty(self):
        @njit
        def foo(items):
            l = listobject.new_list(int32)
            l.extend(items)
            return len(l)

        self.assertEqual(foo((1,)), 1)
        self.assertEqual(foo((1,2)), 2)
        self.assertEqual(foo((1,2,3)), 3)

    def test_list_extend_typing_error_non_iterable(self):
        self.disable_leak_check()

        @njit
        def foo():
            l = listobject.new_list(int32)
            l.extend(1)

        with self.assertRaises(TypingError) as raises:
            foo()
        self.assertIn(
            "extend argument must be iterable",
            str(raises.exception),
        )


class TestInsert(MemoryLeakMixin, TestCase):
    """Test list insert. """

    def test_list_insert_empty(self):
        @njit
        def foo(i):
            l = listobject.new_list(int32)
            l.insert(i, 1)
            return len(l), l[0]

        for i in (-10, -5, -1, 0, 1, 4, 9):
            self.assertEqual(foo(i), (1, 1))

    def test_list_insert_singleton(self):
        @njit
        def foo(i):
            l = listobject.new_list(int32)
            l.append(0)
            l.insert(i, 1)
            return len(l), l[0], l[1]

        # insert before
        for i in (-10, -3, -2, -1, 0):
            self.assertEqual(foo(i), (2, 1, 0))

        # insert after
        for i in (1, 2, 3, 10):
            self.assertEqual(foo(i), (2, 0, 1))

    def test_list_insert_multiple(self):
        @njit
        def foo(i):
            l = listobject.new_list(int32)
            for j in range(10):
                l.append(0)
            l.insert(i, 1)
            return len(l), l[i]

        for i in (0, 4, 9):
            self.assertEqual(foo(i), (11, 1))

    def test_list_insert_multiple_before(self):
        @njit
        def foo(i):
            l = listobject.new_list(int32)
            for j in range(10):
                l.append(0)
            l.insert(i, 1)
            return len(l), l[0]

        for i in (-12, -11, -10, 0):
            self.assertEqual(foo(i), (11, 1))

    def test_list_insert_multiple_after(self):
        @njit
        def foo(i):
            l = listobject.new_list(int32)
            for j in range(10):
                l.append(0)
            l.insert(i, 1)
            return len(l), l[10]

        for i in (10, 11, 12):
            self.assertEqual(foo(i), (11, 1))

    def test_list_insert_typing_error(self):
        self.disable_leak_check()

        @njit
        def foo():
            l = listobject.new_list(int32)
            l.insert("a", 0)

        with self.assertRaises(TypingError) as raises:
            foo()
        self.assertIn(
            "list insert indices must be integers",
            str(raises.exception),
        )


class TestRemove(MemoryLeakMixin, TestCase):
    """Test list remove. """

    def test_list_remove_empty(self):
        self.disable_leak_check()

        @njit
        def foo():
            l = listobject.new_list(int32)
            l.remove(0)

        with self.assertRaises(ValueError):
            foo()

    def test_list_remove_singleton(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            l.append(0)
            l.remove(0)
            return len(l)

        self.assertEqual(foo(), 0)

    def test_list_remove_singleton_value_error(self):
        self.disable_leak_check()

        @njit
        def foo():
            l = listobject.new_list(int32)
            l.append(1)
            l.remove(0)

        with self.assertRaises(ValueError):
            foo()

    def test_list_remove_multiple(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            l.remove(13)
            l.remove(19)
            return len(l)

        self.assertEqual(foo(), 8)

    def test_list_remove_multiple_value_error(self):
        self.disable_leak_check()

        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            l.remove(23)

        with self.assertRaises(ValueError):
            foo()


class TestClear(MemoryLeakMixin, TestCase):
    """Test list clear. """

    def test_list_clear_empty(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            l.clear()
            return len(l)

        self.assertEqual(foo(), 0)

    def test_list_clear_singleton(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            l.append(0)
            l.clear()
            return len(l)

        self.assertEqual(foo(), 0)

    def test_list_clear_multiple(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10):
                l.append(0)
            l.clear()
            return len(l)
        self.assertEqual(foo(), 0)


class TestReverse(MemoryLeakMixin, TestCase):
    """Test list reverse. """

    def test_list_reverse_empty(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            l.reverse()
            return len(l)

        self.assertEqual(foo(), 0)

    def test_list_reverse_singleton(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            l.append(0)
            l.reverse()
            return len(l), l[0]

        self.assertEqual(foo(), (1, 0))

    def test_list_reverse_multiple(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10, 13):
                l.append(j)
            l.reverse()
            return len(l), l[0], l[1], l[2]
        self.assertEqual(foo(), (3, 12, 11, 10))


class TestCopy(MemoryLeakMixin, TestCase):
    """Test list copy. """

    def test_list_copy_empty(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            n = l.copy()
            return len(l), len(n)

        self.assertEqual(foo(), (0, 0))

    def test_list_copy_singleton(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            l.append(0)
            n = l.copy()
            return len(l), len(n), l[0], n[0]

        self.assertEqual(foo(), (1, 1, 0, 0))

    def test_list_copy_multiple(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10, 13):
                l.append(j)
            n = l.copy()
            return len(l), len(n), l[0], l[1], l[2], l[0], l[1], l[2]

        self.assertEqual(foo(), (3, 3, 10, 11, 12, 10, 11, 12))


class TestIndex(MemoryLeakMixin, TestCase):

    def test_index_singleton(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            l.append(1)
            return l.index(1)

        self.assertEqual(foo(), 0)

    def test_index_multiple(self):
        @njit
        def foo(i):
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            return l.index(i)

        for i,v in zip(range(10), range(10,20)):
            self.assertEqual(foo(v), i)

    def test_index_duplicate(self):
        @njit
        def foo():
            l = listobject.new_list(int32)
            for _ in range(10, 20):
                l.append(1)
            return l.index(1)

        self.assertEqual(foo(), 0)

    def test_index_duplicate_with_start(self):
        @njit
        def foo(start):
            l = listobject.new_list(int32)
            for _ in range(10, 20):
                l.append(1)
            return l.index(1, start)

        for i in range(10):
            self.assertEqual(foo(i), i)

    def test_index_singleton_value_error(self):
        self.disable_leak_check()

        @njit
        def foo():
            l = listobject.new_list(int32)
            l.append(0)
            return l.index(1)

        with self.assertRaises(ValueError) as raises:
            foo()
        self.assertIn(
            "item not in list",
            str(raises.exception),
        )

    def test_index_multiple_value_error(self):
        self.disable_leak_check()

        @njit
        def foo():
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            return l.index(23)

        with self.assertRaises(ValueError) as raises:
            foo()
        self.assertIn(
            "item not in list",
            str(raises.exception),
        )

    def test_index_multiple_value_error_start(self):
        self.disable_leak_check()

        @njit
        def foo(start):
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            return l.index(10, start)

        self.assertEqual(foo(0), 0)
        for i in range(1,10):
            with self.assertRaises(ValueError) as raises:
                foo(i)
            self.assertIn(
                "item not in list",
                str(raises.exception),
            )

    def test_index_multiple_value_error_end(self):
        self.disable_leak_check()

        @njit
        def foo(end):
            l = listobject.new_list(int32)
            for j in range(10, 20):
                l.append(j)
            return l.index(19, 0, end)

        self.assertEqual(foo(10), 9)
        for i in range(0,9):
            with self.assertRaises(ValueError) as raises:
                foo(i)
            self.assertIn(
                "item not in list",
                str(raises.exception),
            )

    def test_index_typing_error_start(self):
        self.disable_leak_check()

        @njit
        def foo():
            l = listobject.new_list(int32)
            l.append(0)
            return l.index(0, start="a")

        with self.assertRaises(TypingError) as raises:
            foo()
        self.assertIn(
            "start argument for index must be an integer",
            str(raises.exception),
        )

    def test_index_typing_error_end(self):
        self.disable_leak_check()

        @njit
        def foo():
            l = listobject.new_list(int32)
            l.append(0)
            return l.index(0, end="a")

        with self.assertRaises(TypingError) as raises:
            foo()
        self.assertIn(
            "end argument for index must be an integer",
            str(raises.exception),
        )


class TestEqualNotEqual(MemoryLeakMixin, TestCase):
    """Test list equal and not equal. """

    def test_list_empty_equal(self):
        @njit
        def foo():
            t = listobject.new_list(int32)
            o = listobject.new_list(int32)
            return t == o, t != o

        self.assertEqual(foo(), (True, False))

    def test_list_singleton_equal(self):
        @njit
        def foo():
            t = listobject.new_list(int32)
            t.append(0)
            o = listobject.new_list(int32)
            o.append(0)
            return t == o, t != o

        self.assertEqual(foo(), (True, False))

    def test_list_singleton_not_equal(self):
        @njit
        def foo():
            t = listobject.new_list(int32)
            t.append(0)
            o = listobject.new_list(int32)
            o.append(1)
            return t == o, t != o

        self.assertEqual(foo(), (False, True))

    def test_list_length_mismatch(self):
        @njit
        def foo():
            t = listobject.new_list(int32)
            t.append(0)
            o = listobject.new_list(int32)
            return t == o, t != o

        self.assertEqual(foo(), (False, True))

    def test_list_multiple_equal(self):
        @njit
        def foo():
            t = listobject.new_list(int32)
            o = listobject.new_list(int32)
            for i in range(10):
                t.append(i)
                o.append(i)
            return t == o, t != o

        self.assertEqual(foo(), (True, False))

    def test_list_multiple_not_equal(self):
        @njit
        def foo():
            t = listobject.new_list(int32)
            o = listobject.new_list(int32)
            for i in range(10):
                t.append(i)
                o.append(i)
            o[-1] = 42
            return t == o, t != o

        self.assertEqual(foo(), (False, True))


class TestIter(MemoryLeakMixin, TestCase):
    """Test list iter. """

    def test_list_iter(self):
        @njit
        def foo(items):
            l = listobject.new_list(int32)
            l.extend(items)
            # use a simple sum to check this w/o having to return a list
            r = 0
            for j in l:
                r += j
            return r

        items = (1, 2, 3, 4)

        self.assertEqual(
            foo(items),
            sum(items)
        )

    def test_list_iter_self_mutation(self):
        self.disable_leak_check()

        @njit
        def foo():
            l = listobject.new_list(int32)
            l.extend((1, 2, 3, 4))
            for i in l:
                l.append(i)

        with self.assertRaises(RuntimeError) as raises:
            foo()
        self.assertIn(
            'list was mutated during iteration'.format(**locals()),
            str(raises.exception),
        )


class TestStringItem(MemoryLeakMixin, TestCase):
    """Test list can take strings as items. """

    def test_string_item(self):
        @njit
        def foo():
            l = listobject.new_list(types.unicode_type)
            l.append('a')
            l.append('b')
            l.append('c')
            l.append('d')
            return l[0], l[1], l[2], l[3]

        items = foo()
        self.assertEqual(['a', 'b', 'c', 'd'], list(items))


class TestItemCasting(TestCase):

    @njit
    def foo(fromty, toty):
        l = listobject.new_list(toty)
        l.append(fromty(0))

    def check_good(self, fromty, toty):
        TestItemCasting.foo(fromty, toty)

    def check_bad(self, fromty, toty):
        with self.assertRaises(TypingError) as raises:
            TestItemCasting.foo(fromty, toty)
        self.assertIn(
            'cannot safely cast {fromty} to {toty}'.format(**locals()),
            str(raises.exception),
        )

    def test_cast_int_to(self):
        self.check_good(types.int32, types.float32)
        self.check_good(types.int32, types.float64)
        self.check_good(types.int32, types.complex128)
        self.check_good(types.int64, types.complex128)
        self.check_bad(types.int32, types.complex64)
        self.check_good(types.int8, types.complex64)

    def test_cast_float_to(self):
        self.check_good(types.float32, types.float64)
        self.check_good(types.float32, types.complex64)
        self.check_good(types.float64, types.complex128)

    def test_cast_bool_to(self):
        self.check_good(types.boolean, types.int32)
        self.check_good(types.boolean, types.float64)
        self.check_good(types.boolean, types.complex128)

    def test_cast_fail_unicode_int(self):

        @njit
        def foo():
            l = listobject.new_list(int32)
            l.append("xyz")

        with self.assertRaises(TypingError) as raises:
            foo()
        self.assertIn(
            'cannot safely cast unicode_type to int32',
            str(raises.exception),
        )

    def test_cast_fail_int_unicode(self):

        @njit
        def foo():
            l = listobject.new_list(types.unicode_type)
            l.append(int32(0))

        with self.assertRaises(TypingError) as raises:
            foo()
        self.assertIn(
            'Cannot cast int32 to unicode_type',
            str(raises.exception),
        )


@register_jitable
def make_test_list():
    l = listobject.new_list(int32)
    l.append(int32(1))
    return l


class TestImmutable(MemoryLeakMixin, TestCase):

    def test_is_immutable(self):
        @njit
        def foo():
            l = make_test_list()
            return l._is_mutable()
        self.assertTrue(foo())

    def test_make_immutable_is_immutable(self):
        @njit
        def foo():
            l = make_test_list()
            l._make_immutable()
            return l._is_mutable()
        self.assertFalse(foo())

    def test_length_still_works_when_immutable(self):
        @njit
        def foo():
            l = make_test_list()
            l._make_immutable()
            return len(l),l._is_mutable()
        length, mutable = foo()
        self.assertEqual(length, 1)
        self.assertFalse(mutable)

    def test_getitem_still_works_when_immutable(self):
        @njit
        def foo():
            l = make_test_list()
            l._make_immutable()
            return l[0], l._is_mutable()
        test_item, mutable = foo()
        self.assertEqual(test_item, 1)
        self.assertFalse(mutable)

    def test_append_fails(self):
        self.disable_leak_check()

        @njit
        def foo():
            l = make_test_list()
            l._make_immutable()
            l.append(int32(1))
        with self.assertRaises(ValueError) as raises:
            foo()
        self.assertIn(
            'list is immutable',
            str(raises.exception),
        )

    def test_mutation_fails(self):
        """ Test that any attempt to mutate an immutable typed list fails. """
        self.disable_leak_check()

        def generate_function(line):
            context = {}
            exec(dedent("""
                from numba.typed import listobject
                from numba import int32
                def bar():
                    lst = listobject.new_list(int32)
                    lst.append(int32(1))
                    lst._make_immutable()
                    zero = int32(0)
                    {}
                """.format(line)), context)
            return njit(context["bar"])
        for line in ("lst.append(zero)",
                     "lst[0] = zero",
                     "lst.pop()",
                     "del lst[0]",
                     "lst.extend((zero,))",
                     "lst.insert(0, zero)",
                     "lst.clear()",
                     "lst.reverse()",
                     "lst.sort()",
                     ):
            foo = generate_function(line)
            with self.assertRaises(ValueError) as raises:
                foo()
            self.assertIn(
                "list is immutable",
                str(raises.exception),
            )
