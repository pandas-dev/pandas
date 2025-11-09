import pytest
import pandas as pd
import numpy as np
import datetime


def test_min_with_mixed_nan_and_dates():
    """Test min() with mix of NaN and datetime.date objects"""
    data = {
        "dates": [
            np.nan,
            np.nan,
            datetime.date(2025, 1, 3),
            datetime.date(2025, 1, 4),
        ],
    }
    df = pd.DataFrame(data)
    result = df.min()
    expected = pd.Series({"dates": datetime.date(2025, 1, 3)}, dtype="datetime64[ns]")
    pd.testing.assert_series_equal(result, expected)


def test_min_with_all_nan():
    """Test min() with all NaN values (should return NaN)"""
    data = {
        "dates": [
            np.nan,
            np.nan,
            np.nan,
        ],
    }
    df = pd.DataFrame(data)
    result = df.min()
    expected = pd.Series({"dates": np.nan}, dtype="datetime64[ns]")
    pd.testing.assert_series_equal(result, expected)


def test_min_with_only_dates():
    """Test min() with only datetime.date objects (no NaN)"""
    data = {
        "dates": [
            datetime.date(2025, 1, 3),
            datetime.date(2025, 1, 4),
            datetime.date(2025, 1, 1),
        ],
    }
    df = pd.DataFrame(data)
    result = df.min()
    expected = pd.Series({"dates": datetime.date(2025, 1, 1)}, dtype="datetime64[ns]")
    pd.testing.assert_series_equal(result, expected)


def test_min_with_mixed_types_error():
    """Test min() with incompatible types (should raise TypeError)"""
    data = {
        "mixed": [
            np.nan,
            "string",
            datetime.date(2025, 1, 3),
        ],
    }
    df = pd.DataFrame(data)
    with pytest.raises(TypeError):
        df.max()
