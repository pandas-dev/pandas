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
import string

import numpy as np
import pytest

from pandas.errors import PerformanceWarning

import pandas as pd
import pandas._testing as tm
from pandas.api.types import is_string_dtype
from pandas.core.arrays import ArrowStringArray
from pandas.core.arrays.string_ import StringDtype
from pandas.tests.extension import base


def split_array(arr):
    if arr.dtype.storage != "pyarrow":
        pytest.skip("only applicable for pyarrow chunked array n/a")

    def _split_array(arr):
        import pyarrow as pa

        arrow_array = arr._data
        split = len(arrow_array) // 2
        arrow_array = pa.chunked_array(
            [*arrow_array[:split].chunks, *arrow_array[split:].chunks]
        )
        assert arrow_array.num_chunks == 2
        return type(arr)(arrow_array)

    return _split_array(arr)


@pytest.fixture(params=[True, False])
def chunked(request):
    return request.param


@pytest.fixture
def dtype(string_storage):
    return StringDtype(storage=string_storage)


@pytest.fixture
def data(dtype, chunked):
    strings = np.random.choice(list(string.ascii_letters), size=100)
    while strings[0] == strings[1]:
        strings = np.random.choice(list(string.ascii_letters), size=100)

    arr = dtype.construct_array_type()._from_sequence(strings)
    return split_array(arr) if chunked else arr


@pytest.fixture
def data_missing(dtype, chunked):
    """Length 2 array with [NA, Valid]"""
    arr = dtype.construct_array_type()._from_sequence([pd.NA, "A"])
    return split_array(arr) if chunked else arr


@pytest.fixture
def data_for_sorting(dtype, chunked):
    arr = dtype.construct_array_type()._from_sequence(["B", "C", "A"])
    return split_array(arr) if chunked else arr


@pytest.fixture
def data_missing_for_sorting(dtype, chunked):
    arr = dtype.construct_array_type()._from_sequence(["B", pd.NA, "A"])
    return split_array(arr) if chunked else arr


@pytest.fixture
def na_value():
    return pd.NA


@pytest.fixture
def data_for_grouping(dtype, chunked):
    arr = dtype.construct_array_type()._from_sequence(
        ["B", "B", pd.NA, pd.NA, "A", "A", "B", "C"]
    )
    return split_array(arr) if chunked else arr


class TestDtype(base.BaseDtypeTests):
    def test_eq_with_str(self, dtype):
        assert dtype == f"string[{dtype.storage}]"
        super().test_eq_with_str(dtype)

    def test_is_not_string_type(self, dtype):
        # Different from BaseDtypeTests.test_is_not_string_type
        # because StringDtype is a string type
        assert is_string_dtype(dtype)


class TestInterface(base.BaseInterfaceTests):
    def test_view(self, data, request):
        if data.dtype.storage == "pyarrow":
            pytest.skip(reason="2D support not implemented for ArrowStringArray")
        super().test_view(data)


class TestConstructors(base.BaseConstructorsTests):
    def test_from_dtype(self, data):
        # base test uses string representation of dtype
        pass

    def test_constructor_from_list(self):
        # GH 27673
        pytest.importorskip("pyarrow", minversion="1.0.0")
        result = pd.Series(["E"], dtype=StringDtype(storage="pyarrow"))
        assert isinstance(result.dtype, StringDtype)
        assert result.dtype.storage == "pyarrow"


class TestReshaping(base.BaseReshapingTests):
    def test_transpose(self, data, request):
        if data.dtype.storage == "pyarrow":
            pytest.skip(reason="2D support not implemented for ArrowStringArray")
        super().test_transpose(data)


class TestGetitem(base.BaseGetitemTests):
    pass


class TestSetitem(base.BaseSetitemTests):
    def test_setitem_preserves_views(self, data, request):
        if data.dtype.storage == "pyarrow":
            pytest.skip(reason="2D support not implemented for ArrowStringArray")
        super().test_setitem_preserves_views(data)


class TestIndex(base.BaseIndexTests):
    pass


class TestMissing(base.BaseMissingTests):
    def test_dropna_array(self, data_missing):
        result = data_missing.dropna()
        expected = data_missing[[1]]
        self.assert_extension_array_equal(result, expected)

    def test_fillna_no_op_returns_copy(self, data):
        data = data[~data.isna()]

        valid = data[0]
        result = data.fillna(valid)
        assert result is not data
        self.assert_extension_array_equal(result, data)

        with tm.maybe_produces_warning(
            PerformanceWarning, data.dtype.storage == "pyarrow"
        ):
            result = data.fillna(method="backfill")
        assert result is not data
        self.assert_extension_array_equal(result, data)

    def test_fillna_series_method(self, data_missing, fillna_method):
        with tm.maybe_produces_warning(
            PerformanceWarning,
            fillna_method is not None and data_missing.dtype.storage == "pyarrow",
            check_stacklevel=False,
        ):
            super().test_fillna_series_method(data_missing, fillna_method)


class TestNoReduce(base.BaseNoReduceTests):
    @pytest.mark.parametrize("skipna", [True, False])
    def test_reduce_series_numeric(self, data, all_numeric_reductions, skipna):
        op_name = all_numeric_reductions

        if op_name in ["min", "max"]:
            return None

        ser = pd.Series(data)
        with pytest.raises(TypeError):
            getattr(ser, op_name)(skipna=skipna)


class TestMethods(base.BaseMethodsTests):
    def test_value_counts_with_normalize(self, data):
        data = data[:10].unique()
        values = np.array(data[~data.isna()])
        ser = pd.Series(data, dtype=data.dtype)

        result = ser.value_counts(normalize=True).sort_index()

        expected = pd.Series(
            [1 / len(values)] * len(values), index=result.index, name="proportion"
        )
        if getattr(data.dtype, "storage", "") == "pyarrow":
            expected = expected.astype("double[pyarrow]")
        else:
            expected = expected.astype("Float64")

        self.assert_series_equal(result, expected)


class TestCasting(base.BaseCastingTests):
    pass


class TestComparisonOps(base.BaseComparisonOpsTests):
    def _compare_other(self, ser, data, op, other):
        op_name = f"__{op.__name__}__"
        result = getattr(ser, op_name)(other)
        dtype = "boolean[pyarrow]" if ser.dtype.storage == "pyarrow" else "boolean"
        expected = getattr(ser.astype(object), op_name)(other).astype(dtype)
        self.assert_series_equal(result, expected)

    def test_compare_scalar(self, data, comparison_op):
        ser = pd.Series(data)
        self._compare_other(ser, data, comparison_op, "abc")


class TestParsing(base.BaseParsingTests):
    pass


class TestPrinting(base.BasePrintingTests):
    pass


class TestGroupBy(base.BaseGroupbyTests):
    @pytest.mark.parametrize("as_index", [True, False])
    def test_groupby_extension_agg(self, as_index, data_for_grouping):
        df = pd.DataFrame({"A": [1, 1, 2, 2, 3, 3, 1, 4], "B": data_for_grouping})
        result = df.groupby("B", as_index=as_index).A.mean()
        _, uniques = pd.factorize(data_for_grouping, sort=True)

        if as_index:
            index = pd.Index(uniques, name="B")
            expected = pd.Series([3.0, 1.0, 4.0], index=index, name="A")
            self.assert_series_equal(result, expected)
        else:
            expected = pd.DataFrame({"B": uniques, "A": [3.0, 1.0, 4.0]})
            self.assert_frame_equal(result, expected)

    @pytest.mark.filterwarnings("ignore:Falling back:pandas.errors.PerformanceWarning")
    def test_groupby_extension_apply(self, data_for_grouping, groupby_apply_op):
        super().test_groupby_extension_apply(data_for_grouping, groupby_apply_op)


class Test2DCompat(base.Dim2CompatTests):
    @pytest.fixture(autouse=True)
    def arrow_not_supported(self, data, request):
        if isinstance(data, ArrowStringArray):
            pytest.skip(reason="2D support not implemented for ArrowStringArray")


def test_searchsorted_with_na_raises(data_for_sorting, as_series):
    # GH50447
    b, c, a = data_for_sorting
    arr = data_for_sorting.take([2, 0, 1])  # to get [a, b, c]
    arr[-1] = pd.NA

    if as_series:
        arr = pd.Series(arr)

    msg = (
        "searchsorted requires array to be sorted, "
        "which is impossible with NAs present."
    )
    with pytest.raises(ValueError, match=msg):
        arr.searchsorted(b)
