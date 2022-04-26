"""
Check overflow behavior of operations on timedelta-valued
DataFrames/Series/Indexes/ExtensionArrays.
"""

from itertools import (
    chain,
    combinations_with_replacement,
    product,
)
from operator import attrgetter
from typing import (
    NamedTuple,
    Type,
    Union,
)

from hypothesis import given
import hypothesis.strategies as st
import pytest

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
    to_timedelta,
)
import pandas._testing as tm
from pandas.core.arrays import (
    DatetimeArray,
    ExtensionArray,
    TimedeltaArray,
)

timedelta_types = (Timedelta, TimedeltaArray, TimedeltaIndex, Series, DataFrame)
timestamp_types = (Timestamp, DatetimeArray, DatetimeIndex, Series, DataFrame)
containers = slice(1, None)
get_item_names = lambda t: "-".join(map(attrgetter("__name__"), t))


# TODO: consolidate remaining td64 overflow tests here?
#   - tests/tslibs/test_conversion.py::test_ensure_timedelta64ns_overflows()
#   - tests/tslibs/test_timedeltas.py::test_huge_nanoseconds_overflow()
#   - tests/scalar/timedelta/test_timedelta.py::test_mul_preserves_reso()
#   - tests/scalar/timedelta/test_constructors.py::test_construct_from_td64_with_unit(),
#     test_overflow_on_construction()


class BinaryOpTypes(NamedTuple):
    """
    The expected operand and result types for a binary operation.
    """

    left: Type
    right: Type
    result: Type

    def __str__(self) -> str:
        return get_item_names(self)

    def __repr__(self) -> str:
        return f"BinaryOpTypes({self})"


positive_tds = st.integers(min_value=1, max_value=Timedelta.max.value).map(Timedelta)

xfail_no_overflow_check = pytest.mark.xfail(reason="No overflow check")


@given(
    st.integers(
        min_value=Timedelta.max.value - 2**9 + 1,
        max_value=Timedelta.max.value,
    ).map(Timedelta)
)
def test_td64_summation_raises_spurious_overflow_error_for_single_elem_series(
    value: Timedelta,
):
    s = Series(value)

    msg = "int too big to convert|Python int too large to convert to C long"
    with pytest.raises(OverflowError, match=msg):
        s.sum()


@given(st.integers(min_value=1, max_value=2**10).map(Timedelta))
def test_td64_summation_raises_overflow_error_for_small_overflows(value: Timedelta):
    s = Series([Timedelta.max, value])

    msg = "int too big to convert|Python int too large to convert to C long"
    with pytest.raises(OverflowError, match=msg):
        s.sum()


@given(
    st.integers(
        min_value=2**10 + 1,
        max_value=Timedelta.max.value,
    ).map(Timedelta)
)
def test_td64_summation_raises_value_error_for_most_overflows(value: Timedelta):
    s = Series([Timedelta.max, value])

    msg = "overflow in timedelta operation"
    with pytest.raises(ValueError, match=msg):
        s.sum()


@pytest.fixture(
    name="add_sub_types",
    scope="module",
    params=tuple(combinations_with_replacement(timedelta_types, 2)),
    ids=get_item_names,
)
def fixture_add_sub_types(request) -> BinaryOpTypes:
    """
    Expected types when adding, subtracting Timedeltas.
    """
    return_type = max(request.param, key=lambda t: timedelta_types.index(t))
    return BinaryOpTypes(request.param[0], request.param[1], return_type)


@pytest.fixture(
    name="ts_add_sub_types",
    scope="module",
    params=tuple(product(timedelta_types, timestamp_types)),
    ids=get_item_names,
)
def fixture_ts_add_sub_types(request) -> BinaryOpTypes:
    """
    Expected types when adding, subtracting Timedeltas and Timestamps.
    """
    type_hierarchy = {
        name: i
        for i, name in chain(enumerate(timedelta_types), enumerate(timestamp_types))
    }
    return_type = timestamp_types[max(type_hierarchy[t] for t in request.param)]

    return BinaryOpTypes(request.param[0], request.param[1], return_type)


def wrap_value(value: Union[Timestamp, Timedelta], type_):
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

    return type_(tm.box_expected([value], box_cls))


@given(positive_td=positive_tds)
def test_add_raises_expected_error_if_result_would_overflow(
    add_sub_types: BinaryOpTypes,
    positive_td: Timedelta,
):
    left = wrap_value(Timedelta.max, add_sub_types.left)
    right = wrap_value(positive_td, add_sub_types.right)

    if add_sub_types.result is Timedelta:
        msg = "|".join(
            [
                "int too big to convert",
                "Python int too large to convert to C long",
            ]
        )
    else:
        msg = "Overflow in int64 addition"

    with pytest.raises(OverflowError, match=msg):
        left + right

    with pytest.raises(OverflowError, match=msg):
        right + left


@xfail_no_overflow_check
@given(positive_td=positive_tds)
def test_sub_raises_expected_error_if_result_would_overflow(
    add_sub_types: BinaryOpTypes,
    positive_td: Timedelta,
):
    left = wrap_value(Timedelta.min, add_sub_types.left)
    right = wrap_value(positive_td, add_sub_types.right)

    msg = "Overflow in int64 addition"
    with pytest.raises(OverflowError, match=msg):
        left - right

    with pytest.raises(OverflowError, match=msg):
        (-1 * right) - abs(left)


@given(td_value=positive_tds)
def test_add_timestamp_raises_expected_error_if_result_would_overflow(
    ts_add_sub_types: BinaryOpTypes,
    td_value: Timedelta,
):
    left = wrap_value(td_value, ts_add_sub_types.left)
    right = wrap_value(Timestamp.max, ts_add_sub_types.right)

    ex = (OutOfBoundsDatetime, OverflowError)
    msg = "|".join(["Out of bounds nanosecond timestamp", "Overflow in int64 addition"])

    with pytest.raises(ex, match=msg):
        left + right

    with pytest.raises(ex, match=msg):
        right + left


@xfail_no_overflow_check
@given(td_value=positive_tds)
def test_sub_timestamp_raises_expected_error_if_result_would_overflow(
    ts_add_sub_types: BinaryOpTypes,
    td_value: Timedelta,
):
    right = wrap_value(td_value, ts_add_sub_types[0])
    left = wrap_value(Timestamp.min, ts_add_sub_types[1])

    ex = (OutOfBoundsDatetime, OverflowError)
    msg = "|".join(["Out of bounds nanosecond timestamp", "Overflow in int64 addition"])

    with pytest.raises(ex, match=msg):
        left - right


@given(value=st.floats().filter(lambda f: abs(f) > 1))
def test_scalar_multiplication_raises_expected_error_if_result_would_overflow(
    value: float,
):
    td = Timedelta.max

    msg = "|".join(
        [
            "cannot convert float infinity to integer",
            "Python int too large to convert to C long",
            "int too big to convert",
        ]
    )
    with pytest.raises(OverflowError, match=msg):
        td * value

    with pytest.raises(OverflowError, match=msg):
        value * td


@xfail_no_overflow_check
@given(value=st.floats().filter(lambda f: abs(f) > 1))
@pytest.mark.parametrize(
    argnames="td_type",
    argvalues=timedelta_types[containers],  # type: ignore[arg-type]
    ids=attrgetter("__name__"),
)
def test_container_scalar_multiplication_raises_expected_error_if_result_would_overflow(
    value: float,
    td_type: Type,
):
    td = wrap_value(Timedelta.max, td_type)

    msg = "Overflow in int64 addition"
    with pytest.raises(OverflowError, match=msg):
        td * value

    with pytest.raises(OverflowError, match=msg):
        value * td


class TestAddSubNaTMasking:
    # TODO: parametrize over boxes

    @pytest.mark.parametrize("str_ts", ["1950-01-01", "1980-01-01"])
    def test_tdarr_add_timestamp_nat_masking(self, box_with_array, str_ts):
        # GH#17991 checking for overflow-masking with NaT
        tdinat = to_timedelta(["24658 days 11:15:00", "NaT"])
        tdobj = tm.box_expected(tdinat, box_with_array)

        ts = Timestamp(str_ts)
        ts_variants = [
            ts,
            ts.to_pydatetime(),
            ts.to_datetime64().astype("datetime64[ns]"),
            ts.to_datetime64().astype("datetime64[D]"),
        ]

        for variant in ts_variants:
            res = tdobj + variant
            if box_with_array is DataFrame:
                assert res.iloc[1, 1] is NaT
            else:
                assert res[1] is NaT

    def test_tdi_add_overflow(self):
        # These should not overflow!
        exp = TimedeltaIndex([NaT])
        result = to_timedelta([NaT]) - Timedelta("1 days")
        tm.assert_index_equal(result, exp)

        exp = TimedeltaIndex(["4 days", NaT])
        result = to_timedelta(["5 days", NaT]) - Timedelta("1 days")
        tm.assert_index_equal(result, exp)

        exp = TimedeltaIndex([NaT, NaT, "5 hours"])
        result = to_timedelta([NaT, "5 days", "1 hours"]) + to_timedelta(
            ["7 seconds", NaT, "4 hours"]
        )
        tm.assert_index_equal(result, exp)
