from collections.abc import Mapping
from typing import Any

import numpy as np

from pandas import Index

def assert_dict_equal(
    a: Mapping[Any, Any], b: Mapping[Any, Any], compare_keys: bool = ...
) -> bool: ...
def assert_almost_equal(
    a: object,
    b: object,
    rtol: float = ...,
    atol: float = ...,
    check_dtype: bool = ...,
    obj: object = ...,
    lobj: object = ...,
    robj: object = ...,
    index_values: Index | np.ndarray | None = ...,
) -> bool: ...
