import random
import string

import numpy as np
import pytest

import pandas as pd
from pandas.core.arrays.string_ import StringArray, StringDtype
from pandas.tests.extension import base


@pytest.fixture
def dtype():
    return StringDtype()


@pytest.fixture
def data():
    # strings = random.choices(string.ascii_letters, k=100)
    strings = np.random.choice(list(string.ascii_letters), size=100)
    while strings[0] == strings[1]:
        strings = random.choices(string.ascii_letters, k=100)

    return StringArray._from_sequence(strings)


@pytest.fixture
def data_missing():
    """Length 2 array with [NA, Valid]"""
    return StringArray._from_sequence([np.nan, "A"])


@pytest.fixture
def data_for_sorting():
    return StringArray._from_sequence(["B", "C", "A"])


@pytest.fixture
def data_missing_for_sorting():
    return StringArray._from_sequence(["B", np.nan, "A"])


@pytest.fixture
def na_value():
    return np.nan


@pytest.fixture
def data_for_grouping():
    return StringArray._from_sequence(["B", "B", np.nan, np.nan, "A", "A", "B", "C"])


class TestDtype(base.BaseDtypeTests):
    pass


class TestInterface(base.BaseInterfaceTests):
    pass


class TestConstructors(base.BaseConstructorsTests):
    pass


class TestReshaping(base.BaseReshapingTests):
    pass


class TestGetitem(base.BaseGetitemTests):
    pass


class TestSetitem(base.BaseSetitemTests):
    pass


class TestMissing(base.BaseMissingTests):
    pass


class TestReduce(base.BaseNoReduceTests):
    pass


class TestMethods(base.BaseMethodsTests):
    pass


class TestCasting(base.BaseCastingTests):
    pass


class TestComparisonOps(base.BaseComparisonOpsTests):
    def _compare_other(self, s, data, op_name, other):
        result = getattr(s, op_name)(other)
        expected = getattr(s.astype(object), op_name)(other)
        self.assert_series_equal(result, expected)

    def test_compare_scalar(self, data, all_compare_operators):
        op_name = all_compare_operators
        s = pd.Series(data)
        self._compare_other(s, data, op_name, "abc")


class TestParsing(base.BaseParsingTests):
    pass
