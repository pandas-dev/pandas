"""
Testing that functions from compat work as expected
"""
import builtins
import re

from pandas.compat import lmap, lrange, lzip, re_type


class TestBuiltinIterators:

    @classmethod
    def check_results(cls, results, expecteds, lengths):
        for result, expected, length in zip(results, expecteds, lengths):
            assert isinstance(result, list)
            assert len(result) == length
            assert result == expected

    def test_lrange(self):
        results = lrange(10),
        expecteds = list(builtins.range(10)),
        lengths = 10,

        results += lrange(1, 10, 2),
        lengths += 5,
        expecteds += list(builtins.range(1, 10, 2)),
        self.check_results(results, expecteds, lengths)

    def test_lmap(self):
        func = lambda x, y, z: x + y + z
        lst = [builtins.range(10), builtins.range(10), builtins.range(10)]
        results = lmap(func, *lst),
        expecteds = list(builtins.map(func, *lst)),
        lengths = 10,
        self.check_results(results, expecteds, lengths)

    def test_lzip(self):
        lst = [builtins.range(10), builtins.range(10), builtins.range(10)]
        results = lzip(*lst),
        expecteds = list(builtins.zip(*lst)),
        lengths = 10,
        self.check_results(results, expecteds, lengths)


def test_re_type():
    assert isinstance(re.compile(''), re_type)
