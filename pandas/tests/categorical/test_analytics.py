# -*- coding: utf-8 -*-

import pytest

import numpy as np

import pandas.util.testing as tm
from pandas import Categorical


class TestCategoricalAnalytics(object):

    def test_min_max(self):

        # unordered cats have no min/max
        cat = Categorical(["a", "b", "c", "d"], ordered=False)
        pytest.raises(TypeError, lambda: cat.min())
        pytest.raises(TypeError, lambda: cat.max())
        cat = Categorical(["a", "b", "c", "d"], ordered=True)
        _min = cat.min()
        _max = cat.max()
        assert _min == "a"
        assert _max == "d"
        cat = Categorical(["a", "b", "c", "d"],
                          categories=['d', 'c', 'b', 'a'], ordered=True)
        _min = cat.min()
        _max = cat.max()
        assert _min == "d"
        assert _max == "a"
        cat = Categorical([np.nan, "b", "c", np.nan],
                          categories=['d', 'c', 'b', 'a'], ordered=True)
        _min = cat.min()
        _max = cat.max()
        assert np.isnan(_min)
        assert _max == "b"

        _min = cat.min(numeric_only=True)
        assert _min == "c"
        _max = cat.max(numeric_only=True)
        assert _max == "b"

        cat = Categorical([np.nan, 1, 2, np.nan], categories=[5, 4, 3, 2, 1],
                          ordered=True)
        _min = cat.min()
        _max = cat.max()
        assert np.isnan(_min)
        assert _max == 1

        _min = cat.min(numeric_only=True)
        assert _min == 2
        _max = cat.max(numeric_only=True)
        assert _max == 1

    def test_mode(self):
        s = Categorical([1, 1, 2, 4, 5, 5, 5], categories=[5, 4, 3, 2, 1],
                        ordered=True)
        res = s.mode()
        exp = Categorical([5], categories=[5, 4, 3, 2, 1], ordered=True)
        tm.assert_categorical_equal(res, exp)
        s = Categorical([1, 1, 1, 4, 5, 5, 5], categories=[5, 4, 3, 2, 1],
                        ordered=True)
        res = s.mode()
        exp = Categorical([5, 1], categories=[5, 4, 3, 2, 1], ordered=True)
        tm.assert_categorical_equal(res, exp)
        s = Categorical([1, 2, 3, 4, 5], categories=[5, 4, 3, 2, 1],
                        ordered=True)
        res = s.mode()
        exp = Categorical([5, 4, 3, 2, 1],
                          categories=[5, 4, 3, 2, 1], ordered=True)
        tm.assert_categorical_equal(res, exp)
        # NaN should not become the mode!
        s = Categorical([np.nan, np.nan, np.nan, 4, 5],
                        categories=[5, 4, 3, 2, 1], ordered=True)
        res = s.mode()
        exp = Categorical([5, 4], categories=[5, 4, 3, 2, 1], ordered=True)
        tm.assert_categorical_equal(res, exp)
        s = Categorical([np.nan, np.nan, np.nan, 4, 5, 4],
                        categories=[5, 4, 3, 2, 1], ordered=True)
        res = s.mode()
        exp = Categorical([4], categories=[5, 4, 3, 2, 1], ordered=True)
        tm.assert_categorical_equal(res, exp)
        s = Categorical([np.nan, np.nan, 4, 5, 4], categories=[5, 4, 3, 2, 1],
                        ordered=True)
        res = s.mode()
        exp = Categorical([4], categories=[5, 4, 3, 2, 1], ordered=True)
        tm.assert_categorical_equal(res, exp)
