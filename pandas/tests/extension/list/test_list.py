import pyarrow as pa
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays.list_ import (
    ListArray,
    ListDtype,
)
from pandas.tests.extension.base.accumulate import BaseAccumulateTests
from pandas.tests.extension.base.casting import BaseCastingTests
from pandas.tests.extension.base.constructors import BaseConstructorsTests
from pandas.tests.extension.base.dim2 import (  # noqa: F401
    Dim2CompatTests,
    NDArrayBacked2DTests,
)
from pandas.tests.extension.base.dtype import BaseDtypeTests
from pandas.tests.extension.base.groupby import BaseGroupbyTests
from pandas.tests.extension.base.index import BaseIndexTests
from pandas.tests.extension.base.missing import BaseMissingTests
from pandas.tests.extension.base.ops import (  # noqa: F401
    BaseArithmeticOpsTests,
    BaseComparisonOpsTests,
    BaseOpsUtil,
    BaseUnaryOpsTests,
)
from pandas.tests.extension.base.printing import BasePrintingTests
from pandas.tests.extension.base.reduce import BaseReduceTests

# TODO(wayd): This is copied from string tests - is it required here?
# @pytest.fixture(params=[True, False])
# def chunked(request):
#     return request.param


@pytest.fixture
def dtype():
    return ListDtype(pa.large_string())


@pytest.fixture
def data():
    """Length-100 ListArray for semantics test."""
    # TODO: make better random data
    data = [list("a"), list("ab"), list("abc")] * 33 + [None]
    return ListArray(data)


@pytest.fixture
def data_missing(dtype):
    """Length 2 array with [NA, Valid]"""
    arr = dtype.construct_array_type()._from_sequence([pd.NA, [1, 2, 3]], dtype=dtype)
    return arr


@pytest.fixture
def data_for_grouping(dtype):
    A = ["a"]
    B = ["a", "b"]
    NA = None
    C = ["a", "b", "c"]
    return ListArray([B, B, NA, NA, A, A, B, C])


class TestListArray(
    BaseAccumulateTests,
    BaseCastingTests,
    BaseConstructorsTests,
    BaseDtypeTests,
    # BaseGetitemTests,
    BaseGroupbyTests,
    BaseIndexTests,
    # BaseInterfaceTests,
    # BaseParsingTests,
    # BaseMethodsTests,
    BaseMissingTests,
    # BaseArithmeticOpsTests,
    # BaseComparisonOpsTests,
    # BaseUnaryOpsTests,
    BasePrintingTests,
    BaseReduceTests,
    # BaseReshapingTests,
    # BaseSetitemTests,
    Dim2CompatTests,
):
    # TODO(wayd): The tests here are copied from test_arrow.py
    # It appears the TestArrowArray class has different expectations around
    # when copies should be made then the base.ExtensionTests
    # Assuming intentional, maybe in the long term this should just
    # inherit from TestArrowArray
    def test_fillna_no_op_returns_copy(self, data):
        data = data[~data.isna()]

        valid = data[0]
        result = data.fillna(valid)
        assert result is not data
        tm.assert_extension_array_equal(result, data)

    def test_kind(self, dtype):
        assert dtype.kind == "+L"

    @pytest.mark.parametrize("as_index", [True, False])
    def test_groupby_extension_agg(self, as_index, data_for_grouping):
        pytest.skip(reason="ListArray does not implement mean")

    def test_groupby_extension_no_sort(self, data_for_grouping):
        pytest.skip(reason="ListArray does not implement mean")

    def test_groupby_extension_transform(self, data_for_grouping):
        pytest.skip(reason="ListArray does not implement dictionary_encode")

    def test_groupby_extension_apply(self, data_for_grouping, groupby_apply_op):
        pytest.skip(reason="ListArray does not implement dictionary_encode")


def test_to_csv(data):
    # https://github.com/pandas-dev/pandas/issues/28840
    # array with list-likes fail when doing astype(str) on the numpy array
    # which was done in get_values_for_csv
    df = pd.DataFrame({"a": data})
    res = df.to_csv()
    assert str(data[0]) in res
