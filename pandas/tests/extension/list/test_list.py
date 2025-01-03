import itertools
import operator

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
from pandas.tests.extension.base.interface import BaseInterfaceTests
from pandas.tests.extension.base.io import BaseParsingTests
from pandas.tests.extension.base.missing import BaseMissingTests
from pandas.tests.extension.base.ops import (  # noqa: F401
    BaseArithmeticOpsTests,
    BaseComparisonOpsTests,
    BaseOpsUtil,
    BaseUnaryOpsTests,
)
from pandas.tests.extension.base.printing import BasePrintingTests
from pandas.tests.extension.base.reduce import BaseReduceTests
from pandas.tests.extension.base.reshaping import BaseReshapingTests
from pandas.tests.extension.base.setitem import BaseSetitemTests

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
    BaseInterfaceTests,
    BaseParsingTests,
    # BaseMethodsTests,
    BaseMissingTests,
    BaseArithmeticOpsTests,
    BaseComparisonOpsTests,
    BaseUnaryOpsTests,
    BasePrintingTests,
    BaseReduceTests,
    BaseReshapingTests,
    BaseSetitemTests,
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

    def test_array_interface(self, data):
        pytest.skip(reason="ListArrayScalar does not compare to numpy object-dtype")

    @pytest.mark.parametrize("engine", ["c", "python"])
    def test_EA_types(self, engine, data, request):
        pytest.skip(reason="ListArray has not implemented parsing from string")

    def test_arith_series_with_scalar(self, data, all_arithmetic_operators):
        if all_arithmetic_operators in ("__mod__", "__rmod__"):
            pytest.skip("ListArray does not implement __mod__ or __rmod__")

        super().test_arith_series_with_scalar(data, all_arithmetic_operators)

    def test_arith_series_with_array(self, data, all_arithmetic_operators, request):
        if all_arithmetic_operators in ("__mod__", "__rmod__"):
            pytest.skip("ListArray does not implement __mod__ or __rmod__")

        super().test_arith_series_with_array(data, all_arithmetic_operators)

    def test_arith_frame_with_scalar(self, data, all_arithmetic_operators):
        if all_arithmetic_operators in ("__mod__", "__rmod__"):
            pytest.skip("ListArray does not implement __mod__ or __rmod__")

        super().test_arith_frame_with_scalar(data, all_arithmetic_operators)

    def test_divmod(self, data):
        pytest.skip("ListArray does not implement divmod")

    def test_compare_scalar(self, data, comparison_op):
        if comparison_op in (operator.eq, operator.ne):
            pytest.skip("Series.combine does not properly handle missing values")

        super().test_compare_scalar(data, comparison_op)

    def test_compare_array(self, data, comparison_op):
        if comparison_op in (operator.eq, operator.ne):
            pytest.skip("Series.combine does not properly handle missing values")

        super().test_compare_array(data, comparison_op)

    def test_invert(self, data):
        pytest.skip("ListArray does not implement invert")

    def test_merge_on_extension_array(self, data):
        pytest.skip("ListArray cannot be factorized")

    def test_merge_on_extension_array_duplicates(self, data):
        pytest.skip("ListArray cannot be factorized")

    @pytest.mark.parametrize(
        "index",
        [
            # Two levels, uniform.
            pd.MultiIndex.from_product(([["A", "B"], ["a", "b"]]), names=["a", "b"]),
            # non-uniform
            pd.MultiIndex.from_tuples([("A", "a"), ("A", "b"), ("B", "b")]),
            # three levels, non-uniform
            pd.MultiIndex.from_product([("A", "B"), ("a", "b", "c"), (0, 1, 2)]),
            pd.MultiIndex.from_tuples(
                [
                    ("A", "a", 1),
                    ("A", "b", 0),
                    ("A", "a", 0),
                    ("B", "a", 0),
                    ("B", "c", 1),
                ]
            ),
        ],
    )
    @pytest.mark.parametrize("obj", ["series", "frame"])
    def test_unstack(self, data, index, obj):
        # TODO: the base class test casts everything to object
        # If you remove the object casts, these tests pass...
        # Check if still needed in base class
        data = data[: len(index)]
        if obj == "series":
            ser = pd.Series(data, index=index)
        else:
            ser = pd.DataFrame({"A": data, "B": data}, index=index)

        n = index.nlevels
        levels = list(range(n))
        # [0, 1, 2]
        # [(0,), (1,), (2,), (0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)]
        combinations = itertools.chain.from_iterable(
            itertools.permutations(levels, i) for i in range(1, n)
        )

        for level in combinations:
            result = ser.unstack(level=level)
            assert all(
                isinstance(result[col].array, type(data)) for col in result.columns
            )

            if obj == "series":
                # We should get the same result with to_frame+unstack+droplevel
                df = ser.to_frame()

                alt = df.unstack(level=level).droplevel(0, axis=1)
                tm.assert_frame_equal(result, alt)

            # obj_ser = ser.astype(object)

            expected = ser.unstack(level=level, fill_value=data.dtype.na_value)
            # if obj == "series":
            #    assert (expected.dtypes == object).all()

            # result = result.astype(object)
            tm.assert_frame_equal(result, expected)


def test_to_csv(data):
    # https://github.com/pandas-dev/pandas/issues/28840
    # array with list-likes fail when doing astype(str) on the numpy array
    # which was done in get_values_for_csv
    df = pd.DataFrame({"a": data})
    res = df.to_csv()
    assert str(data[0]) in res
