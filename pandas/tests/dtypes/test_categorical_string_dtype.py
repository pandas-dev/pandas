import numpy as np

from pandas.core.dtypes.common import is_string_dtype
from pandas.core.dtypes.dtypes import CategoricalDtype

import pandas as pd


def test_is_string_dtype_categorical_consistency():
    """Test that is_string_dtype returns consistent results for
    Categorical series and dtype."""
    # Test with CategoricalDtype directly
    categorical_dtype = CategoricalDtype()
    assert not is_string_dtype(categorical_dtype)

    # Test with Series containing Categorical
    categorical_series = pd.Series(pd.Categorical(["a", "b", "c"]))
    assert not is_string_dtype(categorical_series)

    # Test with ordered CategoricalDtype
    ordered_categorical_dtype = CategoricalDtype(ordered=True)
    assert not is_string_dtype(ordered_categorical_dtype)

    # Test with Series containing ordered Categorical
    ordered_categorical_series = pd.Series(
        pd.Categorical(["a", "b", "c"], ordered=True)
    )
    assert not is_string_dtype(ordered_categorical_series)

    # Test with CategoricalDtype with specific categories
    specific_categorical_dtype = CategoricalDtype(categories=["x", "y", "z"])
    assert not is_string_dtype(specific_categorical_dtype)

    # Test with Series containing Categorical with specific categories
    specific_categorical_series = pd.Series(
        pd.Categorical(["x", "y", "z"], categories=["x", "y", "z"])
    )
    assert not is_string_dtype(specific_categorical_series)

    # Test with empty Categorical
    empty_categorical = pd.Series(pd.Categorical([]))
    assert not is_string_dtype(empty_categorical)

    # Test with Categorical containing NaN values
    nan_categorical = pd.Series(pd.Categorical([np.nan, "a", "b"]))
    assert not is_string_dtype(nan_categorical)

    # Test with numeric Categorical
    numeric_categorical = pd.Series(pd.Categorical([1, 2, 3]))
    assert not is_string_dtype(numeric_categorical)
