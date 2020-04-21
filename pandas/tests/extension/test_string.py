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
    strings = np.random.choice(list(string.ascii_letters), size=100)
    while strings[0] == strings[1]:
        strings = np.random.choice(list(string.ascii_letters), size=100)

    return StringArray._from_sequence(strings)


@pytest.fixture
def data_missing():
    """Length 2 array with [NA, Valid]"""
    return StringArray._from_sequence([pd.NA, "A"])


@pytest.fixture
def data_for_sorting():
    return StringArray._from_sequence(["B", "C", "A"])


@pytest.fixture
def data_missing_for_sorting():
    return StringArray._from_sequence(["B", pd.NA, "A"])


@pytest.fixture
def na_value():
    return pd.NA


@pytest.fixture
def data_for_grouping():
    return StringArray._from_sequence(["B", "B", pd.NA, pd.NA, "A", "A", "B", "C"])


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


class TestNoReduce(base.BaseNoReduceTests):
    pass


class TestMethods(base.BaseMethodsTests):
    pass


class TestCasting(base.BaseCastingTests):
    pass


class TestComparisonOps(base.BaseComparisonOpsTests):
    def _compare_other(self, s, data, op_name, other):
        result = getattr(s, op_name)(other)
        expected = getattr(s.astype(object), op_name)(other).astype("boolean")
        self.assert_series_equal(result, expected)

    def test_compare_scalar(self, data, all_compare_operators):
        op_name = all_compare_operators
        s = pd.Series(data)
        self._compare_other(s, data, op_name, "abc")


class TestParsing(base.BaseParsingTests):
    pass


class TestPrinting(base.BasePrintingTests):
    pass


class TestGroupBy(base.BaseGroupbyTests):
    pass
