from typing import (
    Any,
    Callable,
)

import numpy as np

binop = Callable[[Any, Any], Any]
bool_op = Callable[[Any, Any], bool]


def scalar_compare(
    values: np.ndarray,  # object[:]
    val: object,
    op: bool_op,          # {operator.eq, operator.ne, ...}
) -> np.ndarray: ...     # np.ndarray[bool]

def vec_compare(
    left: np.ndarray,   # np.ndarray[object]
    right: np.ndarray,  # np.ndarray[object]
    op: bool_op,         # {operator.eq, operator.ne, ...}
) -> np.ndarray: ...    # np.ndarray[bool]


def scalar_binop(
    values: np.ndarray,   # object[:]
    val: object,
    op: binop,           # binary operator
) -> np.ndarray: ...


def vec_binop(
    left: np.ndarray,   # object[:]
    right: np.ndarray,  # object[:]
    op: binop,         # binary operator
) -> np.ndarray: ...


def maybe_convert_bool(
    arr: np.ndarray,  # np.ndarray[object]
    true_values=None,
    false_values=None
) -> np.ndarray: ...
