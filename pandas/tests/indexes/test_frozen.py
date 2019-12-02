import pytest

from pandas.core.indexes.frozen import FrozenList
from pandas.tests.test_base import CheckImmutable, CheckStringMixin


class TestFrozenList(CheckImmutable, CheckStringMixin):
    mutable_methods = ("extend", "pop", "remove", "insert")
    unicode_container = FrozenList(["\u05d0", "\u05d1", "c"])

    def setup_method(self, _):
        self.lst = [1, 2, 3, 4, 5]
        self.container = FrozenList(self.lst)
        self.klass = FrozenList

    def test_add(self):
        result = self.container + (1, 2, 3)
        expected = FrozenList(self.lst + [1, 2, 3])
        self.check_result(result, expected)

        result = (1, 2, 3) + self.container
        expected = FrozenList([1, 2, 3] + self.lst)
        self.check_result(result, expected)

    def test_iadd(self):
        q = r = self.container

        q += [5]
        self.check_result(q, self.lst + [5])

        # Other shouldn't be mutated.
        self.check_result(r, self.lst)

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

    def test_tricky_container_to_bytes_raises(self):
        # GH 26447
        msg = "^'str' object cannot be interpreted as an integer$"
        with pytest.raises(TypeError, match=msg):
            bytes(self.unicode_container)
