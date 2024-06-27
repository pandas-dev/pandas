# Contents in this file are referenced from the sphinx-generated docs.
# "magictoken" is used for markers as beginning and ending of example text.

import unittest
from numba.tests.support import TestCase


class DocsJitclassUsageTest(TestCase):

    def test_ex_jitclass(self):
        # magictoken.ex_jitclass.begin
        import numpy as np
        from numba import int32, float32    # import the types
        from numba.experimental import jitclass

        spec = [
            ('value', int32),               # a simple scalar field
            ('array', float32[:]),          # an array field
        ]

        @jitclass(spec)
        class Bag(object):
            def __init__(self, value):
                self.value = value
                self.array = np.zeros(value, dtype=np.float32)

            @property
            def size(self):
                return self.array.size

            def increment(self, val):
                for i in range(self.size):
                    self.array[i] += val
                return self.array

            @staticmethod
            def add(x, y):
                return x + y

        n = 21
        mybag = Bag(n)
        # magictoken.ex_jitclass.end

        self.assertTrue(isinstance(mybag, Bag))
        self.assertPreciseEqual(mybag.value, n)
        np.testing.assert_allclose(mybag.array, np.zeros(n, dtype=np.float32))
        self.assertPreciseEqual(mybag.size, n)
        np.testing.assert_allclose(mybag.increment(3),
                                   3 * np.ones(n, dtype=np.float32))
        np.testing.assert_allclose(mybag.increment(6),
                                   9 * np.ones(n, dtype=np.float32))
        self.assertPreciseEqual(mybag.add(1, 1), 2)
        self.assertPreciseEqual(Bag.add(1, 2), 3)

    def test_ex_jitclass_type_hints(self):
        # magictoken.ex_jitclass_type_hints.begin
        from typing import List
        from numba.experimental import jitclass
        from numba.typed import List as NumbaList

        @jitclass
        class Counter:
            value: int

            def __init__(self):
                self.value = 0

            def get(self) -> int:
                ret = self.value
                self.value += 1
                return ret

        @jitclass
        class ListLoopIterator:
            counter: Counter
            items: List[float]

            def __init__(self, items: List[float]):
                self.items = items
                self.counter = Counter()

            def get(self) -> float:
                idx = self.counter.get() % len(self.items)
                return self.items[idx]

        items = NumbaList([3.14, 2.718, 0.123, -4.])
        loop_itr = ListLoopIterator(items)
        # magictoken.ex_jitclass_type_hints.end

        for idx in range(10):
            self.assertEqual(loop_itr.counter.value, idx)
            self.assertAlmostEqual(loop_itr.get(), items[idx % len(items)])
            self.assertEqual(loop_itr.counter.value, idx + 1)


if __name__ == '__main__':
    unittest.main()
