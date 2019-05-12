"""
Testing that functions from compat work as expected
"""
import builtins

from pandas.compat import lrange


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
