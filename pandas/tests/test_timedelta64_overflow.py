"""
Check overflow behavior of operations on timedelta-valued
ExtensionArrays/Indexes/Series/DataFrames.
"""

from contextlib import (
    AbstractContextManager,
    nullcontext,
)
from functools import partial
from typing import (
    List,
    Type,
    Union,
)

import pytest

from pandas._libs.lib import is_list_like
from pandas.errors import OutOfBoundsDatetime

from pandas import (
    DataFrame,
    DatetimeIndex,
    Index,
    NaT,
    Series,
    Timedelta,
    TimedeltaIndex,
    Timestamp,
    array,
)
import pandas._testing as tm
from pandas.core.arrays import (
    DatetimeArray,
    ExtensionArray,
    TimedeltaArray,
)

td64_types = (Timedelta, TimedeltaArray, TimedeltaIndex, Series, DataFrame)
td64_box_types = td64_types[slice(1, None)]
td64_arraylike_types = td64_types[slice(1, 4)]
dt64_types = (Timestamp, DatetimeArray, DatetimeIndex, Series, DataFrame)

TD64_TYPE = Union[Timedelta, TimedeltaArray, TimedeltaIndex, Series, DataFrame]
TD64_BOX_TYPE = Union[TimedeltaArray, TimedeltaIndex, Series, DataFrame]
TD64_ARRAYLIKE_TYPE = Union[TimedeltaArray, TimedeltaIndex, Series]
DT64_TYPE = Union[Timestamp, DatetimeArray, DatetimeIndex, Series, DataFrame]


TD64_VALUE_ERROR_MSG = "overflow in timedelta operation"
TD64_OVERFLOW_MSG = "|".join(
    [
        "int too big to convert",
        "Python int too large to convert to C long",
        "Overflow in int64 addition",
    ]
)

does_not_raise = nullcontext
raises_overflow_error = partial(pytest.raises, OverflowError, match=TD64_OVERFLOW_MSG)
raises_value_error = partial(pytest.raises, ValueError, match=TD64_VALUE_ERROR_MSG)


def wrap_value(value, type_):
    """
    Return value wrapped in a container of given type_, or as-is if type_ is a scalar.
    """
    if issubclass(type_, (Timedelta, Timestamp)):
        return type_(value)

    if issubclass(type_, ExtensionArray):
        box_cls = array
    elif issubclass(type_, Index):
        box_cls = Index
    else:
        box_cls = type_

    if not is_list_like(value):
        value = [value]
    return tm.box_expected(value, box_cls, transpose=False)


@pytest.fixture(name="td64_type", params=td64_types, scope="module")
def fixture_td64_type(request) -> Type[TD64_TYPE]:
    return request.param


@pytest.fixture(name="td64_arraylike_type", params=td64_arraylike_types, scope="module")
def fixture_td64_arraylike_type(request) -> Type[TD64_ARRAYLIKE_TYPE]:
    return request.param


@pytest.fixture(name="td64_box_type", params=td64_box_types, scope="module")
def fixture_td64_box_type(request) -> Type[TD64_BOX_TYPE]:
    return request.param


@pytest.fixture(name="dt64_type", params=dt64_types, scope="module")
def fixture_dt64_type(request) -> Type[DT64_TYPE]:
    return request.param


@pytest.fixture(name="max_td64")
def fixture_max_td64(td64_box_type: Type[TD64_BOX_TYPE]) -> TD64_BOX_TYPE:
    """
    A 1-elem ExtensionArray/Index/Series, or 2x1 DataFrame, w/ all elements set to
    Timestamp.max.
    """
    return wrap_value(Timedelta.max, td64_box_type)


@pytest.fixture(
    name="positive_td64",
    params=[Timedelta(1), Timedelta(1024), Timedelta.max],
    ids=["1ns", "1024ns", "td_max"],
)
def fixture_positive_td64(request, td64_type: Type[TD64_TYPE]) -> TD64_TYPE:
    """
    A scalar, 1-elem ExtensionArray/Index/Series, or 2x1 DataFrame.
    """
    value = request.param
    return wrap_value(value, td64_type)


class TestBoxReductionMethods:
    """
    For timedelta64-valued ExtensionArrays/Indexes/Series/DataFrames.
    """

    @pytest.mark.parametrize(
        ["value", "expected_exs"],
        [
            (Timedelta.min, does_not_raise()),
            (Timedelta.min + Timedelta(511), does_not_raise()),
            (Timedelta.max - Timedelta(511), raises_overflow_error()),
            (Timedelta.max, raises_overflow_error()),
        ],
    )
    def test_arraylike_sum_fails_with_large_single_elem(
        self,
        value: Timedelta,
        expected_exs: AbstractContextManager,
        td64_arraylike_type: Type[TD64_ARRAYLIKE_TYPE],
    ):
        td64_arraylike = wrap_value(value, td64_arraylike_type)
        with expected_exs:
            result = td64_arraylike.sum()
            # for large negative values, sum() doesn't raise but does return NaT
            assert result is NaT

    @pytest.mark.parametrize(
        ("values", "expected_exs"),
        (
            ([Timedelta.min] * 2, raises_value_error()),
            ([Timedelta.min, Timedelta(-1025)], raises_value_error()),
            ([Timedelta.min, Timedelta(-1024)], does_not_raise()),
            ([Timedelta.min, Timedelta(-1)], does_not_raise()),
            ([Timedelta.max, Timedelta(1)], raises_overflow_error()),
            ([Timedelta.max, Timedelta(1024)], raises_overflow_error()),
            ([Timedelta.max, Timedelta(1025)], raises_value_error()),
            ([Timedelta.max] * 2, raises_value_error()),
        ),
    )
    def test_arraylike_sum_usually_raises_for_overflow(
        self,
        values: List[Timedelta],
        expected_exs: AbstractContextManager,
        td64_arraylike_type: Type[TD64_ARRAYLIKE_TYPE],
    ):
        td64_arraylike = wrap_value(values, td64_arraylike_type)
        with expected_exs:
            result = td64_arraylike.sum()
            # for small negative overflows, sum() doesn't raise but does return NaT
            assert result is NaT

    @pytest.mark.parametrize(
        "values",
        (
            [Timedelta.min] * 2,
            [Timedelta.min, Timedelta(-1)],
            [Timedelta.max, Timedelta(1)],
            [Timedelta.max] * 2,
        ),
        ids=["double_td_min", "over_by_-1ns", "over_by_1ns", "double_td_max"],
    )
    def test_df_sum_returns_nat_for_all_overflows(self, values: List[Timedelta]):
        td64_df = wrap_value(values, DataFrame)
        result = td64_df.sum()
        expected = Series(NaT, index=[0], dtype="timedelta64[ns]")

        tm.assert_series_equal(result, expected)


class TestBinaryOps:
    """
    Operations between timedelta64-valued ExtensionArrays/Indexes/Series/DataFrames, and
    a numeric/timelike scalar or timelike-valued ExtensionArray/Index/Series/DataFrame.
    """

    def test_add_raises_if_result_would_overflow(
        self,
        max_td64: TD64_TYPE,
        positive_td64: TD64_BOX_TYPE,
    ):
        with raises_overflow_error():
            max_td64 + positive_td64

        with raises_overflow_error():
            positive_td64 + max_td64

    @pytest.mark.parametrize(
        ["rval", "expected_exs"],
        [
            (Timedelta(1), does_not_raise()),
            (Timedelta(2), raises_overflow_error()),
            (Timedelta.max, raises_overflow_error()),
        ],
    )
    def test_sub_raises_if_result_would_overflow(
        self,
        max_td64: TD64_BOX_TYPE,
        rval: Timedelta,
        expected_exs: AbstractContextManager,
        td64_type: Type[TD64_TYPE],
    ):
        rvalue = wrap_value(rval, td64_type)
        min_td64 = -1 * max_td64

        with expected_exs:
            min_td64 - rvalue

        with expected_exs:
            -1 * rvalue - max_td64

    def test_add_dt64_raises_if_result_would_overflow(
        self,
        max_td64: TD64_BOX_TYPE,
        dt64_type: Type[DT64_TYPE],
    ):
        max_dt64 = wrap_value(Timestamp.max, dt64_type)
        ex = (OutOfBoundsDatetime, OverflowError)
        msg = TD64_OVERFLOW_MSG + "|Out of bounds nanosecond timestamp"

        with pytest.raises(ex, match=msg):
            max_td64 + max_dt64

        with pytest.raises(ex, match=msg):
            max_dt64 + max_td64

    def test_sub_td64_raises_if_result_would_overflow(
        self,
        max_td64: TD64_BOX_TYPE,
        dt64_type: Type[DT64_TYPE],
    ):
        min_dt64 = wrap_value(Timestamp.min, dt64_type)
        ex = (OutOfBoundsDatetime, OverflowError)
        msg = TD64_OVERFLOW_MSG + "|Out of bounds nanosecond timestamp"

        with pytest.raises(ex, match=msg):
            min_dt64 - max_td64

    @pytest.mark.xfail(reason="Not implemented")
    def test_scalar_mul_raises_if_result_would_overflow(self, max_td64: TD64_BOX_TYPE):
        with raises_overflow_error():
            max_td64 * 1.01

        with raises_overflow_error():
            1.01 * max_td64
