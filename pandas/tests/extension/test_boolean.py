"""
This file contains a minimal set of tests for compliance with the extension
array interface test suite, and should contain no other tests.
The test suite for the full functionality of the array is located in
`pandas/tests/arrays/`.

The tests in this file are inherited from the BaseExtensionTests, and only
minimal tweaks should be applied to get the tests passing (by overwriting a
parent method).

Additional tests should either be added to one of the BaseExtensionTests
classes (if they are relevant for the extension interface for all dtypes), or
be added to the array-specific tests in `pandas/tests/arrays/`.

"""
import numpy as np
import pytest

import pandas as pd
from pandas.core.arrays.boolean import BooleanDtype
from pandas.tests.extension import base
import pandas.util.testing as tm


def make_data():
    return [True, False] * 4 + [np.nan] + [True, False] * 44 + [np.nan] + [True, False]


@pytest.fixture
def dtype():
    return BooleanDtype()


@pytest.fixture
def data(dtype):
    return pd.array(make_data(), dtype=dtype)


@pytest.fixture
def data_for_twos(dtype):
    return pd.array(np.ones(100), dtype=dtype)


@pytest.fixture
def data_missing(dtype):
    return pd.array([np.nan, True], dtype=dtype)


@pytest.fixture
def data_for_sorting(dtype):
    return pd.array([True, True, False], dtype=dtype)


@pytest.fixture
def data_missing_for_sorting(dtype):
    return pd.array([True, np.nan, False], dtype=dtype)


@pytest.fixture
def na_cmp():
    # we are np.nan
    return lambda x, y: np.isnan(x) and np.isnan(y)


@pytest.fixture
def na_value():
    return np.nan


@pytest.fixture
def data_for_grouping(dtype):
    b = True
    a = False
    na = np.nan
    return pd.array([b, b, na, na, a, a, b], dtype=dtype)


class TestDtype(base.BaseDtypeTests):
    pass


class TestInterface(base.BaseInterfaceTests):
    pass


class TestConstructors(base.BaseConstructorsTests):
    pass


class TestGetitem(base.BaseGetitemTests):
    pass


class TestSetitem(base.BaseSetitemTests):
    pass


class TestMissing(base.BaseMissingTests):
    pass


# class TestArithmeticOps(base.BaseArithmeticOpsTests):
#    pass


# class TestComparisonOps(base.BaseComparisonOpsTests):
#     pass


class TestReshaping(base.BaseReshapingTests):
    pass


class TestMethods(base.BaseMethodsTests):
    @pytest.mark.parametrize("na_sentinel", [-1, -2])
    def test_factorize(self, data_for_grouping, na_sentinel):
        # override because we only have 2 unique values
        labels, uniques = pd.factorize(data_for_grouping, na_sentinel=na_sentinel)
        expected_labels = np.array(
            [0, 0, na_sentinel, na_sentinel, 1, 1, 0], dtype=np.intp
        )
        expected_uniques = data_for_grouping.take([0, 4])

        tm.assert_numpy_array_equal(labels, expected_labels)
        self.assert_extension_array_equal(uniques, expected_uniques)

    def test_combine_le(self, data_repeated):
        # override because expected needs to be boolean instead of bool dtype
        orig_data1, orig_data2 = data_repeated(2)
        s1 = pd.Series(orig_data1)
        s2 = pd.Series(orig_data2)
        result = s1.combine(s2, lambda x1, x2: x1 <= x2)
        expected = pd.Series(
            [a <= b for (a, b) in zip(list(orig_data1), list(orig_data2))],
            dtype="boolean",
        )
        self.assert_series_equal(result, expected)

        val = s1.iloc[0]
        result = s1.combine(val, lambda x1, x2: x1 <= x2)
        expected = pd.Series([a <= val for a in list(orig_data1)], dtype="boolean")
        self.assert_series_equal(result, expected)

    def test_searchsorted(self, data_for_sorting, as_series):
        # override because we only have 2 unique values
        data_for_sorting = pd.array([True, False], dtype="boolean")
        b, a = data_for_sorting
        arr = type(data_for_sorting)._from_sequence([a, b])

        if as_series:
            arr = pd.Series(arr)
        assert arr.searchsorted(a) == 0
        assert arr.searchsorted(a, side="right") == 1

        assert arr.searchsorted(b) == 1
        assert arr.searchsorted(b, side="right") == 2

        result = arr.searchsorted(arr.take([0, 1]))
        expected = np.array([0, 1], dtype=np.intp)

        tm.assert_numpy_array_equal(result, expected)

        # sorter
        sorter = np.array([1, 0])
        assert data_for_sorting.searchsorted(a, sorter=sorter) == 0


class TestCasting(base.BaseCastingTests):
    pass


class TestGroupby(base.BaseGroupbyTests):
    """
    Groupby-specific tests are overridden because boolean only has 2
    unique values, base tests uses 3 groups.
    """

    def test_grouping_grouper(self, data_for_grouping):
        df = pd.DataFrame(
            {"A": ["B", "B", None, None, "A", "A", "B"], "B": data_for_grouping}
        )
        gr1 = df.groupby("A").grouper.groupings[0]
        gr2 = df.groupby("B").grouper.groupings[0]

        tm.assert_numpy_array_equal(gr1.grouper, df.A.values)
        tm.assert_extension_array_equal(gr2.grouper, data_for_grouping)

    @pytest.mark.parametrize("as_index", [True, False])
    def test_groupby_extension_agg(self, as_index, data_for_grouping):
        df = pd.DataFrame({"A": [1, 1, 2, 2, 3, 3, 1], "B": data_for_grouping})
        result = df.groupby("B", as_index=as_index).A.mean()
        _, index = pd.factorize(data_for_grouping, sort=True)

        index = pd.Index(index, name="B")
        expected = pd.Series([3, 1], index=index, name="A")
        if as_index:
            self.assert_series_equal(result, expected)
        else:
            expected = expected.reset_index()
            self.assert_frame_equal(result, expected)

    def test_groupby_extension_no_sort(self, data_for_grouping):
        df = pd.DataFrame({"A": [1, 1, 2, 2, 3, 3, 1], "B": data_for_grouping})
        result = df.groupby("B", sort=False).A.mean()
        _, index = pd.factorize(data_for_grouping, sort=False)

        index = pd.Index(index, name="B")
        expected = pd.Series([1, 3], index=index, name="A")
        self.assert_series_equal(result, expected)

    def test_groupby_extension_transform(self, data_for_grouping):
        valid = data_for_grouping[~data_for_grouping.isna()]
        df = pd.DataFrame({"A": [1, 1, 3, 3, 1], "B": valid})

        result = df.groupby("B").A.transform(len)
        expected = pd.Series([3, 3, 2, 2, 3], name="A")

        self.assert_series_equal(result, expected)

    def test_groupby_extension_apply(self, data_for_grouping, groupby_apply_op):
        df = pd.DataFrame({"A": [1, 1, 2, 2, 3, 3, 1], "B": data_for_grouping})
        df.groupby("B").apply(groupby_apply_op)
        df.groupby("B").A.apply(groupby_apply_op)
        df.groupby("A").apply(groupby_apply_op)
        df.groupby("A").B.apply(groupby_apply_op)

    def test_groupby_apply_identity(self, data_for_grouping):
        df = pd.DataFrame({"A": [1, 1, 2, 2, 3, 3, 1], "B": data_for_grouping})
        result = df.groupby("A").B.apply(lambda x: x.array)
        expected = pd.Series(
            [
                df.B.iloc[[0, 1, 6]].array,
                df.B.iloc[[2, 3]].array,
                df.B.iloc[[4, 5]].array,
            ],
            index=pd.Index([1, 2, 3], name="A"),
            name="B",
        )
        self.assert_series_equal(result, expected)

    def test_in_numeric_groupby(self, data_for_grouping):
        df = pd.DataFrame(
            {
                "A": [1, 1, 2, 2, 3, 3, 1],
                "B": data_for_grouping,
                "C": [1, 1, 1, 1, 1, 1, 1],
            }
        )
        result = df.groupby("A").sum().columns

        if data_for_grouping.dtype._is_numeric:
            expected = pd.Index(["B", "C"])
        else:
            expected = pd.Index(["C"])

        tm.assert_index_equal(result, expected)


class TestNumericReduce(base.BaseNumericReduceTests):
    pass


class TestBooleanReduce(base.BaseBooleanReduceTests):
    pass


class TestPrinting(base.BasePrintingTests):
    pass


# TODO parsing not yet supported
# class TestParsing(base.BaseParsingTests):
#     pass
