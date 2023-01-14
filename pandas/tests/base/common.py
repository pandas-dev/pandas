from typing import Any

from pandas import Index
from pandas.api.types import is_bool_dtype


def allow_na_ops(obj: Any) -> bool:
    """Whether to skip test cases including NaN"""
    is_bool_index = isinstance(obj, Index) and is_bool_dtype(obj)
    return not is_bool_index and obj._can_hold_na
