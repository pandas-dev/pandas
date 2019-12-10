from typing import Any

import numpy as np

from pandas import Index, Series
import pandas.util.testing as tm


def allow_na_ops(obj: Any) -> bool:
    """Whether to skip test cases including NaN"""
    is_bool_index = isinstance(obj, Index) and obj.is_boolean()
    return not is_bool_index and obj._can_hold_na


def check_ops_properties_valid(obj: Any, op: str) -> None:
    """ Validates that certain properties are available """
    if isinstance(obj, Series):
        expected = Series(getattr(obj.index, op), index=obj.index, name="a")
    else:
        expected = getattr(obj, op)

    result = getattr(obj, op)

    # these could be series, arrays or scalars
    if isinstance(result, Series) and isinstance(expected, Series):
        tm.assert_series_equal(result, expected)
    elif isinstance(result, Index) and isinstance(expected, Index):
        tm.assert_index_equal(result, expected)
    elif isinstance(result, np.ndarray) and isinstance(expected, np.ndarray):
        tm.assert_numpy_array_equal(result, expected)
    else:
        assert result == expected
