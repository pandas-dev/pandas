"""
Boilerplate functions used in defining binary operations.
"""

from __future__ import annotations

from functools import wraps
from typing import (
    TYPE_CHECKING,
    Any,
)

import numpy as np

from pandas._libs.lib import item_from_zerodim
from pandas._libs.missing import is_matching_na

from pandas.core.dtypes.generic import (
    ABCExtensionArray,
    ABCIndex,
    ABCSeries,
)

if TYPE_CHECKING:
    from collections.abc import Callable

    from pandas._typing import (
        ArrayLike,
        F,
    )


def get_shape_exception_message(left: ArrayLike, right: ArrayLike) -> str:
    """
    Find the standardized exception message to give for operations between
    arrays of mismatched length or shape.
    """
    if left.ndim == right.ndim == 1:
        return "Lengths must match"
    else:
        return "Shapes must match"


def get_op_exception_message(op_name: str, left: ArrayLike, right: Any) -> str:
    """
    Find the standardized exception message to give for op(left, right).
    """
    if isinstance(right, (np.ndarray, ABCExtensionArray)):
        msg = (
            f"Cannot perform operation '{op_name}' between object "
            f"with dtype '{left.dtype}' and "
            f"dtype '{right.dtype}'"
        )
    else:
        msg = (
            f"Cannot perform operation '{op_name}' between object "
            f"with dtype '{left.dtype}' and "
            f"type '{type(right).__name__}'"
        )
    return msg


def unpack_zerodim_and_defer(name: str) -> Callable[[F], F]:
    """
    Boilerplate for pandas conventions in arithmetic and comparison methods.

    Parameters
    ----------
    name : str

    Returns
    -------
    decorator
    """

    def wrapper(method: F) -> F:
        return _unpack_zerodim_and_defer(method, name)

    return wrapper


def _unpack_zerodim_and_defer(method: F, name: str) -> F:
    """
    Boilerplate for pandas conventions in arithmetic and comparison methods.

    Ensure method returns NotImplemented when operating against "senior"
    classes.  Ensure zero-dimensional ndarrays are always unpacked.

    Parameters
    ----------
    method : binary method
    name : str

    Returns
    -------
    method
    """

    @wraps(method)
    def new_method(self, other):
        prio = getattr(other, "__pandas_priority__", None)
        if prio is not None:
            if prio > self.__pandas_priority__:
                # e.g. other is DataFrame while self is Index/Series/EA
                return NotImplemented

        other = item_from_zerodim(other)

        if isinstance(self, ABCExtensionArray):
            if isinstance(other, (np.ndarray, ABCExtensionArray)):
                if not self._supports_array_op(other, name):
                    msg = get_op_exception_message(name, self, other)
                    raise TypeError(msg)

                if other.shape != self.shape:
                    msg = get_shape_exception_message(self, other)
                    raise ValueError(msg)
            else:
                if not self._supports_scalar_op(other, name):
                    msg = get_op_exception_message(name, self, other)
                    raise TypeError(msg)

        return method(self, other)

    # error: Incompatible return value type (got "Callable[[Any, Any], Any]",
    # expected "F")
    return new_method  # type: ignore[return-value]


def get_op_result_name(left, right):
    """
    Find the appropriate name to pin to an operation result.  This result
    should always be either an Index or a Series.

    Parameters
    ----------
    left : {Series, Index}
    right : object

    Returns
    -------
    name : object
        Usually a string
    """
    if isinstance(right, (ABCSeries, ABCIndex)):
        name = _maybe_match_name(left, right)
    else:
        name = left.name
    return name


def _maybe_match_name(a, b):
    """
    Try to find a name to attach to the result of an operation between
    a and b.  If only one of these has a `name` attribute, return that
    name.  Otherwise return a consensus name if they match or None if
    they have different names.

    Parameters
    ----------
    a : object
    b : object

    Returns
    -------
    name : str or None

    See Also
    --------
    pandas.core.common.consensus_name_attr
    """
    a_has = hasattr(a, "name")
    b_has = hasattr(b, "name")
    if a_has and b_has:
        try:
            if a.name == b.name:
                return a.name
            elif is_matching_na(a.name, b.name):
                # e.g. both are np.nan
                return a.name
            else:
                return None
        except TypeError:
            # pd.NA
            if is_matching_na(a.name, b.name):
                return a.name
            return None
        except ValueError:
            # e.g. np.int64(1) vs (np.int64(1), np.int64(2))
            return None
    elif a_has:
        return a.name
    elif b_has:
        return b.name
    return None
