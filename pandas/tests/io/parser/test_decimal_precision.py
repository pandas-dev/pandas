import io
import decimal
import pandas as pd
import numpy as np
import pytest

CSV = """a,b
0.123456789012345678901234567890,1.23456789012345678901234567890
"""

def test_high_precision_default_float():
    df = pd.read_csv(io.StringIO(CSV))
    # Default: float64, precision truncated
    assert df["a"].dtype == np.float64
    assert str(df["a"][0]).startswith("0.1234567890123456")
    assert not str(df["a"][0]).startswith("0.12345678901234567890")

def test_high_precision_object():
    df = pd.read_csv(io.StringIO(CSV), dtype="object")
    # Should preserve full string
    assert df["a"].dtype == object
    assert df["a"][0] == "0.123456789012345678901234567890"

def test_high_precision_decimal_dtype():
    df = pd.read_csv(io.StringIO(CSV), dtype={"a": decimal.Decimal, "b": decimal.Decimal}, converters=None)
    assert isinstance(df["a"][0], decimal.Decimal)
    assert df["a"][0] == decimal.Decimal("0.123456789012345678901234567890")

def test_high_precision_decimal_converter():
    df = pd.read_csv(io.StringIO(CSV), converters={"a": decimal.Decimal, "b": decimal.Decimal})
    assert isinstance(df["a"][0], decimal.Decimal)
    assert df["a"][0] == decimal.Decimal("0.123456789012345678901234567890")

def test_roundtrip_object():
    df = pd.read_csv(io.StringIO(CSV), dtype="object")
    out = io.StringIO()
    df.to_csv(out, index=False)
    out.seek(0)
    df2 = pd.read_csv(out, dtype="object")
    assert df2.equals(df)
