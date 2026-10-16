from typing import Any

import pandas as pd


def allow_na_ops(obj: Any) -> bool:
    """Whether to skip test cases including NaN"""
    is_bool_index = isinstance(obj, pd.Index) and obj.inferred_type == "boolean"
    return not is_bool_index and obj._can_hold_na
