import string

import pytest
import numpy as np

from pandas.api.types import CategoricalDtype
from pandas import Categorical
from pandas.tests.extension import base


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


@pytest.fixture
def data_for_sorting():
    return Categorical(['A', 'B', 'C'], categories=['C', 'A', 'B'],
                       ordered=True)


@pytest.fixture
def data_missing_for_sorting():
    return Categorical(['A', None, 'B'], categories=['B', 'A'],
                       ordered=True)


@pytest.fixture
def na_value():
    return np.nan


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
    @pytest.mark.skip(reason="Unobserved categories preseved in concat.")
    def test_align(self, data, na_value):
        pass

    @pytest.mark.skip(reason="Unobserved categories preseved in concat.")
    def test_align_frame(self, data, na_value):
        pass


class TestGetitem(base.BaseGetitemTests):
    @pytest.mark.skip(reason="Backwards compatibility")
    def test_getitem_scalar(self):
        # CategoricalDtype.type isn't "correct" since it should
        # be a parent of the elements (object). But don't want
        # to break things by changing.
        pass


class TestMissing(base.BaseMissingTests):

    @pytest.mark.skip(reason="Not implemented")
    def test_fillna_limit_pad(self):
        pass

    @pytest.mark.skip(reason="Not implemented")
    def test_fillna_limit_backfill(self):
        pass


class TestMethods(base.BaseMethodsTests):
    pass

    @pytest.mark.skip(reason="Unobserved categories included")
    def test_value_counts(self, all_data, dropna):
        pass


class TestCasting(base.BaseCastingTests):
    pass
