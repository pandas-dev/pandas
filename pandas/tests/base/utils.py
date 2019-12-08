from typing import Any, Callable, List, Type

import numpy as np
import pytest

from pandas import Index, Series
from pandas.core.indexes.datetimelike import DatetimeIndexOpsMixin
import pandas.util.testing as tm


def allow_na_ops(obj: Any) -> bool:
    """Whether to skip test cases including NaN"""
    is_bool_index = isinstance(obj, Index) and obj.is_boolean()
    return not is_bool_index and obj._can_hold_na


def check_ops_properties_valid(obj: Any, props: List[str], filter: Callable) -> None:
    """ Validates that certain properties are available """
    for op in props:
        # if a filter, skip if it doesn't match
        index = obj.index if isinstance(obj, Series) else obj
        if not filter(index):
            continue

        try:
            if isinstance(obj, Series):
                expected = Series(getattr(obj.index, op), index=obj.index, name="a")
            else:
                expected = getattr(obj, op)
        except AttributeError:
            continue

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


def check_ops_properties_invalid(obj: Any, props: List[str]) -> None:
    """ Validates that certain properties are not available """
    for op in props:
        # freq raises AttributeError on an Int64Index because its not
        # defined we mostly care about Series here anyhow

        # an object that is datetimelike will raise a TypeError,
        # otherwise an AttributeError
        err: Type[Exception] = AttributeError
        if issubclass(type(obj), DatetimeIndexOpsMixin):
            err = TypeError

        with pytest.raises(err):
            getattr(obj, op)
