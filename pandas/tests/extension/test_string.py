import string

import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
from pandas.core.arrays.string_ import StringDtype
from pandas.core.arrays.string_arrow import ArrowStringDtype
from pandas.tests.extension import base


@pytest.fixture(
    params=[
        StringDtype,
        pytest.param(
            ArrowStringDtype, marks=td.skip_if_no("pyarrow", min_version="1.0.0")
        ),
    ]
)
def dtype(request):
    return request.param()


@pytest.fixture
def data(dtype):
    strings = np.random.choice(list(string.ascii_letters), size=100)
    while strings[0] == strings[1]:
        strings = np.random.choice(list(string.ascii_letters), size=100)

    return dtype.construct_array_type()._from_sequence(strings)


@pytest.fixture
def data_missing(dtype):
    """Length 2 array with [NA, Valid]"""
    return dtype.construct_array_type()._from_sequence([pd.NA, "A"])


@pytest.fixture
def data_for_sorting(dtype):
    return dtype.construct_array_type()._from_sequence(["B", "C", "A"])


@pytest.fixture
def data_missing_for_sorting(dtype):
    return dtype.construct_array_type()._from_sequence(["B", pd.NA, "A"])


@pytest.fixture
def na_value():
    return pd.NA


@pytest.fixture
def data_for_grouping(dtype):
    return dtype.construct_array_type()._from_sequence(
        ["B", "B", pd.NA, pd.NA, "A", "A", "B", "C"]
    )


class TestDtype(base.BaseDtypeTests):
    pass


class TestInterface(base.BaseInterfaceTests):
    def test_view(self, data, request):
        if isinstance(data.dtype, ArrowStringDtype):
            mark = pytest.mark.xfail(reason="not implemented")
            request.node.add_marker(mark)
        super().test_view(data)


class TestConstructors(base.BaseConstructorsTests):
    pass


class TestReshaping(base.BaseReshapingTests):
    def test_transpose(self, data, dtype, request):
        if isinstance(dtype, ArrowStringDtype):
            mark = pytest.mark.xfail(reason="not implemented")
            request.node.add_marker(mark)
        super().test_transpose(data)


class TestGetitem(base.BaseGetitemTests):
    pass


class TestSetitem(base.BaseSetitemTests):
    def test_setitem_preserves_views(self, data, dtype, request):
        if isinstance(dtype, ArrowStringDtype):
            mark = pytest.mark.xfail(reason="not implemented")
            request.node.add_marker(mark)
        super().test_setitem_preserves_views(data)


class TestMissing(base.BaseMissingTests):
    pass


class TestNoReduce(base.BaseNoReduceTests):
    @pytest.mark.parametrize("skipna", [True, False])
    def test_reduce_series_numeric(self, data, all_numeric_reductions, skipna):
        op_name = all_numeric_reductions

        if op_name in ["min", "max"]:
            return None

        s = pd.Series(data)
        with pytest.raises(TypeError):
            getattr(s, op_name)(skipna=skipna)


class TestMethods(base.BaseMethodsTests):
    @pytest.mark.skip(reason="returns nullable")
    def test_value_counts(self, all_data, dropna):
        return super().test_value_counts(all_data, dropna)

    @pytest.mark.skip(reason="returns nullable")
    def test_value_counts_with_normalize(self, data):
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
