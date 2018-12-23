# -*- coding: utf-8 -*-

from pandas import Categorical
from pandas.api.types import CategoricalDtype
import pandas.util.testing as tm


class TestCategoricalSubclassing(object):

    def test_constructor(self):
        subclassed = tm.SubclassedCategorical(['a', 'b', 'c'])
        assert isinstance(subclassed, tm.SubclassedCategorical)
        tm.assert_categorical_equal(subclassed, Categorical(['a', 'b', 'c']))

    def test_from_codes(self):
        dtype = CategoricalDtype(['a', 'b', 'c'])
        subclassed = tm.SubclassedCategorical.from_codes([1, 0, 2],
                                                         dtype=dtype)
        assert isinstance(subclassed, tm.SubclassedCategorical)

        expected = Categorical.from_codes([1, 0, 2], dtype=dtype)
        tm.assert_categorical_equal(subclassed, expected)

    def test_map(self):
        subclassed = tm.SubclassedCategorical(['a', 'b', 'c'])
        result = subclassed.map(lambda x: x.upper())
        assert isinstance(result, tm.SubclassedCategorical)
        expected = Categorical(['A', 'B', 'C'])
        tm.assert_categorical_equal(result, expected)
