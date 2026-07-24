"""
Testing numba implementation of the numba set.

The tests here only check that the numba typing and codegen are working
correctly.  Detailed testing of the underlying set operations is done
in test_setimpl.py.
"""

from numba import njit
from numba import int32
from numba.typed import Set, setobject
from numba.tests.support import TestCase, MemoryLeakMixin


class TestSetObject(MemoryLeakMixin, TestCase):
    def test_set_create(self):
        """
        Test set creation, insertion and len
        """
        @njit
        def foo(n):
            s = setobject.new_set(int32)
            for i in range(n):
                s.add(i + 1)
            return len(s)

        # Insert nothing
        self.assertEqual(foo(n=0), 0)
        # Insert 1 entry
        self.assertEqual(foo(n=1), 1)
        # Insert 2 entries
        self.assertEqual(foo(n=2), 2)
        # Insert 100 entries
        self.assertEqual(foo(n=100), 100)

        @njit
        def foo(n):
            s = Set.empty(int32)
            for i in range(n):
                s.add(i + 1)
            return len(s)

        # Insert nothing
        self.assertEqual(foo(n=0), 0)
        # Insert 1 entry
        self.assertEqual(foo(n=1), 1)
        # Insert 2 entries
        self.assertEqual(foo(n=2), 2)
        # Insert 100 entries
        self.assertEqual(foo(n=100), 100)

    def test_box_unbox(self):
        """
        Test set boxing and unboxing
        """
        s = Set.empty(int32)

        self.assertEqual(len(s), 0)

        @njit
        def foo(s):
            return s

        s = foo(s)
        self.assertEqual(len(s), 0)

    def test_set_add(self):
        """
        Test set add method
        """
        s = Set.empty(int32)
        # Add elements to a set
        for i in range(4):
            s.add(i + 1)

        # length of set should equal number of elements added
        self.assertEqual(len(s), 4)

        @njit
        def foo(s):
            # Add extra elements to the set
            for i in range(4):
                s.add(-i)
            return s

        s = foo(s)
        # Length of set should equal total number of elements added
        self.assertEqual(len(s), 8)

    def test_set_add_duplicates(self):
        """
        Test set add method
        """
        s = Set.empty(int32)
        # Add elements to a set
        for i in range(4):
            s.add(i + 1)

        # length of set should equal number of elements added
        self.assertEqual(len(s), 4)

        @njit
        def foo(s):
            # Add some overlapping elements to the set
            for i in range(8):
                s.add(i + 1)
            return s

        s = foo(s)
        # Length of set should equal total number of
        # non overlapping elements added
        self.assertEqual(len(s), 8)

    def test_contains(self):
        """
        Test set contains method
        """
        s = Set.empty(int32)
        num_keys = 10

        for key in range(num_keys):
            s.add(key)

        self.assertEqual(len(s), num_keys)

        @njit
        def foo(set, key):
            return bool(key in set)

        for key in range(num_keys):
            self.assertEqual(bool(key in s), True)
            self.assertEqual(foo(s, key), True)

    def test_discard(self):
        """
        Test set discard method
        """
        s = Set.empty(int32)
        # Add elements to a set
        for i in range(4):
            s.add(i + 1)
        self.assertEqual(len(s), 4)

        # Discard an element from the set
        s.discard(1) # In the set
        s.discard(5) # Not in the set
        self.assertEqual(len(s), 3)

        @njit
        def foo(s):
            # Discard another element from the set
            # from within JIT compiled function
            s.discard(1)
            s.discard(2)
            return s

        s = foo(s)
        # Length of set should be number of remaining elements
        self.assertEqual(len(s), 2)
        # Both the discarded elements shouldn't be present in the set
        self.assertEqual(bool(1 in s), False)
        self.assertEqual(bool(2 in s), False)

    def test_iter(self):
        # Test iteration
        s = Set.empty(int32)
        nmax = 1000

        # Add elements to the set
        for i in range(nmax):
            s.add(i)

        @njit
        def foo(set):
            keys = []
            # Iterate over the set
            for key in set:
                keys.append(key)
            return keys

        keys = foo.py_func(s)
        # Check every key was present in the iteration
        self.assertEqual(len(s), len(keys))
        for i in range(nmax):
            self.assertIn(i, keys)

        keys = foo(s)
        # Check every key was present in the iteration
        self.assertEqual(len(s), len(keys))
        for i in range(nmax):
            self.assertIn(i, keys)

    def test_equality(self):
        s1 = Set.empty(int32)
        s2 = Set.empty(int32)

        for i in range(10):
            s1.add(i)
            s2.add(i)

        self.assertEqual(s1, s2)

        @njit
        def foo(s1, s2):
            return s1 == s2

        self.assertEqual(foo(s1, s2), True)
