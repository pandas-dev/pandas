import string

import pytest
import numpy as np

import pandas as pd
import pandas.util.testing as tm
from pandas.api.types import CategoricalDtype
from pandas import Categorical
from .base import BaseArrayTests, BaseDtypeTests


class TestCategoricalDtype(BaseDtypeTests):
    @pytest.fixture
    def dtype(self):
        return CategoricalDtype()


def make_data():
    return np.random.choice(list(string.ascii_letters), size=100)


class TestCategoricalArray(BaseArrayTests):

    @pytest.fixture
    def data(self):
        """Length-100 PeriodArray for semantics test."""
        return Categorical(make_data())

    @pytest.fixture
    def data_missing(self):
        """Length 2 array with [NA, Valid]"""
        return Categorical([np.nan, 'A'])

    @pytest.mark.skip(reason="Memory usage doesn't match")
    def test_memory_usage(self):
        # Is this deliberate?
        pass

    @pytest.mark.skip(reason="Backwards compatability")
    def test_getitem_scalar(self):
        # CategoricalDtype.type isn't "correct" since it should
        # be a parent of the elements (object). But don't want
        # to break things by changing.
        pass

    def test_align(self, data):
        # Override to pass through dtype
        a = data[:3]
        b = data[2:5]
        r1, r2 = pd.Series(a).align(pd.Series(b, index=[1, 2, 3]))

        e1 = pd.Series(type(data)(list(a) + [data._fill_value],
                                  dtype=data.dtype))
        e2 = pd.Series(type(data)([data._fill_value] + list(b),
                                  dtype=data.dtype))
        tm.assert_series_equal(r1, e1)
        tm.assert_series_equal(r2, e2)

    @pytest.mark.skip(reason="Different value_counts semantics.")
    def test_value_counts(self, all_data, dropna):
        pass
