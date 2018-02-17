import string

import pytest
import numpy as np

import pandas as pd
import pandas.util.testing as tm
from pandas.api.types import CategoricalDtype
from pandas import Categorical
from . import base


def make_data():
    return np.random.choice(list(string.ascii_letters), size=100)


@pytest.fixture
def dtype():
    return CategoricalDtype()


@pytest.fixture
def data():
    """Length-100 PeriodArray for semantics test."""
    return Categorical(make_data())


@pytest.fixture
def data_missing():
    """Length 2 array with [NA, Valid]"""
    return Categorical([np.nan, 'A'])


class TestDtype(base.BaseDtypeTests):
    pass


class TestInterface(base.BaseInterfaceTests):
    @pytest.mark.skip(reason="Memory usage doesn't match")
    def test_memory_usage(self):
        # Is this deliberate?
        pass


class TestConstructors(base.BaseConstructorsTests):
    pass


class TestReshaping(base.BaseReshapingTests):
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


class TestGetitem(base.BaseGetitemTests):
    @pytest.mark.skip(reason="Backwards compatability")
    def test_getitem_scalar(self):
        # CategoricalDtype.type isn't "correct" since it should
        # be a parent of the elements (object). But don't want
        # to break things by changing.
        pass


class TestMissing(base.BaseMissingTests):
    pass


class TestMethods(base.BaseMethodsTests):
    pass

    @pytest.mark.skip(reason="Different value_counts semantics.")
    def test_value_counts(self, all_data, dropna):
        pass
