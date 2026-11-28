"""
Run the base extension tests against an EA whose ``dtype.na_value`` is None.

See GH#44602.
"""

from __future__ import annotations

import decimal
import functools

import numpy as np
import pytest

from pandas.api.extensions import (
    register_extension_dtype,
    take,
)
from pandas.api.types import (
    is_float,
    is_integer,
    is_list_like,
    is_scalar,
)
from pandas.core.indexers import check_array_indexer
from pandas.tests.extension import base
from pandas.tests.extension.decimal import test_decimal
from pandas.tests.extension.decimal.array import (
    DecimalArray,
    DecimalDtype,
    make_data,
)


def _to_decimal_or_none(val):
    # reject float NaN so the tests catch pandas filling NaN instead of None
    if val is None or isinstance(val, decimal.Decimal):
        return val
    if is_integer(val) or (is_float(val) and not np.isnan(val)):
        return decimal.Decimal(val)
    raise TypeError(f"Cannot convert {val!r} to Decimal or None")


@register_extension_dtype
class DecimalNoneDtype(DecimalDtype):
    name = "decimal_none"
    na_value = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        return f"DecimalNoneDtype(context={self.context})"

    def construct_array_type(self):
        return DecimalNoneArray


class DecimalNoneArray(DecimalArray):
    def __init__(self, values, dtype=None, copy=False, context=None) -> None:
        values = np.asarray(values, dtype=object)
        for i, val in enumerate(values):
            new_val = _to_decimal_or_none(val)
            if new_val is not val:
                values[i] = new_val

        self._data = self._items = self.data = values
        self._dtype = DecimalNoneDtype(context)

    def __setitem__(self, key, value) -> None:
        if self._readonly:
            raise ValueError("Cannot modify read-only array")

        if is_list_like(value):
            if is_scalar(key):
                raise ValueError("setting an array element with a sequence.")
            value = [_to_decimal_or_none(v) for v in value]
        else:
            value = _to_decimal_or_none(value)

        key = check_array_indexer(self, key)
        self._data[key] = value

    def take(self, indexer, allow_fill=False, fill_value=None):
        # take treats fill_value=None as "use NaN", so fill by hand
        indexer = np.asarray(indexer, dtype=np.intp)
        result = take(self._data, indexer, allow_fill=allow_fill, fill_value=fill_value)
        if allow_fill and fill_value is None:
            result[indexer == -1] = None
        return type(self)(result, context=self.dtype.context)

    def __contains__(self, item) -> bool | np.bool_:
        if item is None:
            return self.isna().any()
        if isinstance(item, decimal.Decimal) and item.is_nan():
            return False
        return super().__contains__(item)

    def isna(self):
        return np.array([x is None for x in self._data], dtype=bool)

    @property
    def _na_value(self):
        return None

    @classmethod
    def _create_arithmetic_method(cls, op):
        # propagate None instead of raising on e.g. Decimal - None
        @functools.wraps(op)
        def na_op(left, right):
            if left is None or right is None:
                return (None, None) if op.__name__ in ["divmod", "rdivmod"] else None
            return op(left, right)

        return cls._create_method(na_op)


DecimalNoneArray._add_arithmetic_ops()


@pytest.fixture
def dtype():
    return DecimalNoneDtype()


@pytest.fixture
def data():
    return DecimalNoneArray(make_data(10))


@pytest.fixture
def data_for_twos():
    return DecimalNoneArray([decimal.Decimal(2) for _ in range(10)])


@pytest.fixture
def data_missing():
    return DecimalNoneArray([None, decimal.Decimal(1)])


@pytest.fixture
def data_for_sorting():
    return DecimalNoneArray(
        [decimal.Decimal("1"), decimal.Decimal("2"), decimal.Decimal("0")]
    )


@pytest.fixture
def data_missing_for_sorting():
    return DecimalNoneArray([decimal.Decimal("1"), None, decimal.Decimal("0")])


@pytest.fixture
def data_for_grouping():
    b = decimal.Decimal("1.0")
    a = decimal.Decimal("0.0")
    c = decimal.Decimal("2.0")
    return DecimalNoneArray([b, b, None, None, a, a, b, c])


@pytest.mark.filterwarnings(
    "ignore:DecimalNoneArray uses the default:pandas.errors.PerformanceWarning"
)
class TestDecimalNoneArray(test_decimal.TestDecimalArray):
    def test_fillna_with_none(self, data_missing):
        # None is the NA value, so unlike in the parent class this is a no-op
        base.BaseMissingTests.test_fillna_with_none(self, data_missing)
