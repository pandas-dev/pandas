import numpy as np

from pandas import DataFrame, Series
import pandas.util.testing as tm

from pandas import cut


class TestCut(object):

    def test_cut_duplicates_drop(self):
        values = Series(np.array([1, 3, 5, 7, 9]),index=["a","b","c","d","e"])
        results = cut(values, [0,2,4,6,10,10], labels=False, right=False, duplicates="drop")
        expected = DataFrame({"a": 0,
                              "b": 1,
                              "c": 2,
                              "d": 3,
                              "e": 3})
        assert_frame_equal(result, expected)

    def test_cut_duplicates_raise(self):
        values = Series(np.array([1, 3, 5, 7, 9]),index=["a","b","c","d","e"])
        assertRaises(ValueError, cut, values, [0,2,4,6,10,10], duplicates='raise')
