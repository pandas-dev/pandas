from datetime import tzinfo
from typing import Optional

import numpy as np

from pandas import Timestamp

def _test_parse_iso8601(ts: str) -> Timestamp: ...


def format_array_from_datetime(
    values: np.ndarray,   # np.ndarray[np.int64]
    tz: Optional[tzinfo] = ...,
    format: Optional[str] = ...,
    na_rep: object = ...
) -> np.ndarray: ...  # np.ndarray[object]


def array_with_unit_to_datetime(
    values: np.ndarray,
    unit: str,
    errors: str = ...,
) -> tuple[np.ndarray, Optional[tzinfo]]: ...


def array_to_datetime(
    values: np.ndarray,  # np.ndarray[object]
    errors: str = ...,
    dayfirst: bool = ...,
    yearfirst: bool = ...,
    utc: bool = ...,
    require_iso8601: bool = ...,
    allow_mixed: bool = ...,
) -> tuple[np.ndarray, Optional[tzinfo]]: ...
# returned ndarray may be object dtype or datetime64[ns]
