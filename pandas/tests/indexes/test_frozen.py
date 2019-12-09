import re

import pytest

from pandas.core.indexes.frozen import FrozenList


class CheckImmutableMixin:
    mutable_regex = re.compile("does not support mutable operations")

    def check_mutable_error(self, *args, **kwargs):
        # Pass whatever function you normally would to pytest.raises
        # (after the Exception kind).
        with pytest.raises(TypeError):
            self.mutable_regex(*args, **kwargs)

    def test_no_mutable_funcs(self):
        def setitem():
            self.container[0] = 5

        self.check_mutable_error(setitem)

        def setslice():
            self.container[1:2] = 3

        self.check_mutable_error(setslice)

        def delitem():
            del self.container[0]

        self.check_mutable_error(delitem)

        def delslice():
            del self.container[0:3]

        self.check_mutable_error(delslice)
        mutable_methods = getattr(self, "mutable_methods", [])

        for meth in mutable_methods:
            self.check_mutable_error(getattr(self.container, meth))

    def test_slicing_maintains_type(self):
        result = self.container[1:2]
        expected = self.lst[1:2]
        self.check_result(result, expected)

    def check_result(self, result, expected, klass=None):
        klass = klass or self.klass
        assert isinstance(result, klass)
        assert result == expected


class CheckStringMixin:
    def test_string_methods_dont_fail(self):
        repr(self.container)
        str(self.container)
        bytes(self.container)

    def test_tricky_container(self):
        if not hasattr(self, "unicode_container"):
            pytest.skip("Need unicode_container to test with this")
        repr(self.unicode_container)
        str(self.unicode_container)


class TestFrozenList(CheckImmutableMixin, CheckStringMixin):
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
