"""
Templates for invalid operations.
"""

from __future__ import annotations

import operator
from typing import (
    TYPE_CHECKING,
    Any,
    NoReturn,
)

import numpy as np

from pandas._libs import lib

from pandas.core.dtypes.dtypes import (
    ArrowDtype,
    CategoricalDtype,
    PeriodDtype,
)

if TYPE_CHECKING:
    from collections.abc import Callable

    from pandas._typing import (
        ArrayLike,
        Scalar,
        npt,
    )


def invalid_comparison(
    left: ArrayLike,
    right: ArrayLike | list | range | Scalar,
    op: Callable[[Any, Any], bool],
) -> npt.NDArray[np.bool_]:
    """
    If a comparison has mismatched types and is not necessarily meaningful,
    follow python3 conventions by:

        - returning all-False for equality
        - returning all-True for inequality
        - raising TypeError otherwise

    Parameters
    ----------
    left : array-like
    right : scalar, array-like
    op : operator.{eq, ne, lt, le, gt}

    Raises
    ------
    TypeError : on inequality comparisons
    """
    if op is operator.eq:
        res_values = np.zeros(left.shape, dtype=bool)
    elif op is operator.ne:
        res_values = np.ones(left.shape, dtype=bool)
    else:
        typ = type(right).__name__
        raise TypeError(f"Invalid comparison between dtype={left.dtype} and {typ}")
    return res_values


def _is_datetimelike_array(obj: object) -> bool:
    if lib.is_scalar(obj):
        return False
    dtype = getattr(obj, "dtype", None)
    if dtype is None:
        return False
    if dtype.kind in "biufc":
        # no numeric or bool dtype is datetimelike, and logical_op is hot
        return False
    if isinstance(dtype, CategoricalDtype) and dtype.categories is not None:
        # a Categorical hides its categories behind kind "O"
        dtype = dtype.categories.dtype
    if isinstance(dtype, ArrowDtype):
        import pyarrow as pa

        pa_type = dtype.pyarrow_dtype
        if pa.types.is_dictionary(pa_type) or pa.types.is_run_end_encoded(pa_type):
            # the arrow spellings of the CategoricalDtype unwrap above
            pa_type = pa_type.value_type
        # kind is "O" for those two, for time32/time64 and for month_day_nano_interval
        return pa.types.is_temporal(pa_type)
    # kind covers numpy M8/m8 and DatetimeTZDtype; PeriodDtype has kind "O" and needs
    #  naming.  Plain object dtype is excluded on purpose, see
    #  test_logical_op_object_dtype_still_truthy
    return dtype.kind in "mM" or isinstance(dtype, PeriodDtype)


def _operand_repr(obj: object) -> str:
    dtype = getattr(obj, "dtype", None)
    if dtype is None:
        return f"object of type {type(obj).__name__}"
    return f"dtype '{dtype}'"


def disallow_datetimelike_logical_op(
    left: ArrayLike, right: Any, op: Callable[[Any, Any], Any]
) -> None:
    """
    Raise TypeError if either operand of a logical op is a datetimelike array.

    Parameters
    ----------
    left : array-like
    right : array-like or scalar
    op : operator.{and_, or_, xor}
        Or one of the reversed variants from roperator.

    Raises
    ------
    TypeError : if either operand is a datetimelike array, including a categorical
        of one
    """
    # GH#68452 these have no truth value; casting them to bool would make every
    #  entry, NaT included, True
    if _is_datetimelike_array(left) or _is_datetimelike_array(right):
        raise TypeError(
            f"operation '{op.__name__}' not supported for "
            f"{_operand_repr(left)} with {_operand_repr(right)}"
        )


def make_invalid_op(name: str) -> Callable[..., NoReturn]:
    """
    Return a binary method that always raises a TypeError.

    Parameters
    ----------
    name : str

    Returns
    -------
    invalid_op : function
    """

    def invalid_op(self: object, other: object = None) -> NoReturn:
        typ = type(self).__name__
        raise TypeError(f"cannot perform {name} with this index type: {typ}")

    invalid_op.__name__ = name
    return invalid_op
