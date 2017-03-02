import numpy as np
from pandas.util import testing as tm
from pandas.tests.test_base import CheckImmutable, CheckStringMixin
from pandas.indexes.frozen import FrozenList, FrozenNDArray
from pandas.compat import u


class TestFrozenList(CheckImmutable, CheckStringMixin, tm.TestCase):
    mutable_methods = ('extend', 'pop', 'remove', 'insert')
    unicode_container = FrozenList([u("\u05d0"), u("\u05d1"), "c"])

    def setUp(self):
        self.lst = [1, 2, 3, 4, 5]
        self.container = FrozenList(self.lst)
        self.klass = FrozenList

    def test_add(self):
        q = FrozenList([1])
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            q = q + [2, 3]
        expected = FrozenList([1, 2, 3])
        self.check_result(q, expected)

    def test_iadd(self):
        q = FrozenList([1])
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            q += [2, 3]
        expected = FrozenList([1, 2, 3])
        self.check_result(q, expected)

    def test_union(self):
        result = self.container.union((1, 2, 3))
        expected = FrozenList(self.lst + [1, 2, 3])
        self.check_result(result, expected)

    def test_difference(self):
        result = self.container.difference([2])
        expected = FrozenList([1, 3, 4, 5])
        self.check_result(result, expected)

    def test_difference_dupe(self):
        result = FrozenList([1, 2, 3, 2]).difference([2])
        expected = FrozenList([1, 3])
        self.check_result(result, expected)


class TestFrozenNDArray(CheckImmutable, CheckStringMixin, tm.TestCase):
    mutable_methods = ('put', 'itemset', 'fill')
    unicode_container = FrozenNDArray([u("\u05d0"), u("\u05d1"), "c"])

    def setUp(self):
        self.lst = [3, 5, 7, -2]
        self.container = FrozenNDArray(self.lst)
        self.klass = FrozenNDArray

    def test_shallow_copying(self):
        original = self.container.copy()
        self.assertIsInstance(self.container.view(), FrozenNDArray)
        self.assertFalse(isinstance(
            self.container.view(np.ndarray), FrozenNDArray))
        self.assertIsNot(self.container.view(), self.container)
        self.assert_numpy_array_equal(self.container, original)
        # shallow copy should be the same too
        self.assertIsInstance(self.container._shallow_copy(), FrozenNDArray)

        # setting should not be allowed
        def testit(container):
            container[0] = 16

        self.check_mutable_error(testit, self.container)

    def test_values(self):
        original = self.container.view(np.ndarray).copy()
        n = original[0] + 15
        vals = self.container.values()
        self.assert_numpy_array_equal(original, vals)
        self.assertIsNot(original, vals)
        vals[0] = n
        self.assertIsInstance(self.container, FrozenNDArray)
        self.assert_numpy_array_equal(self.container.values(), original)
        self.assertEqual(vals[0], n)
