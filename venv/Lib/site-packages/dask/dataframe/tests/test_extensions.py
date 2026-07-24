from __future__ import annotations

from decimal import Decimal

import pytest

import dask.dataframe as dd
from dask.dataframe.utils import assert_eq

pd = pytest.importorskip("pandas")

from pandas.tests.extension.decimal.array import DecimalArray, DecimalDtype

from dask.dataframe.extensions import make_array_nonempty, make_scalar


@make_array_nonempty.register(DecimalDtype)
def _(dtype):
    return DecimalArray._from_sequence([Decimal("0"), Decimal("NaN")], dtype=dtype)


@make_scalar.register(Decimal)
def _(x):
    return Decimal("1")


def test_register_extension_type():
    arr = DecimalArray._from_sequence([Decimal("1.0")] * 10)
    ser = pd.Series(arr)
    dser = dd.from_pandas(ser, 2)
    assert_eq(ser, dser)

    df = pd.DataFrame({"A": ser})
    ddf = dd.from_pandas(df, 2)
    assert_eq(df, ddf)


def test_reduction():
    ser = pd.Series(DecimalArray._from_sequence([Decimal("0"), Decimal("1")]))
    dser = dd.from_pandas(ser, 2)
    assert_eq(ser.mean(skipna=False), dser.mean(skipna=False))

    # It's unclear whether this can be reliably provided, at least with the current
    # implementation, which uses pandas.DataFrame.sum(), returning a (homogeneous)
    # series which has potentially cast values.

    # assert_eq(ser.to_frame().mean(skipna=False), dser.to_frame().mean(skipna=False))


def test_scalar():
    result = dd.utils.make_meta(Decimal("1.0"), parent_meta=pd.DataFrame())
    assert result == Decimal("1.0")
