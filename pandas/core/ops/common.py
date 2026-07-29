"""
Boilerplate functions used in defining binary operations.
"""

from __future__ import annotations

from functools import wraps
from typing import TYPE_CHECKING
import warnings

import numpy as np

from pandas._libs.lib import (
    is_list_like,
    item_from_zerodim,
)
from pandas._libs.missing import is_matching_na
from pandas.errors import Pandas4Warning
from pandas.util._exceptions import find_stack_level

from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCExtensionArray,
    ABCIndex,
    ABCSeries,
)

from pandas.core.construction import (
    ensure_wrapped_if_datetimelike,
    sanitize_array,
)

if TYPE_CHECKING:
    from collections.abc import Callable

    from pandas._typing import F


def has_castable_attr(obj) -> bool:
    attrs = ["__array__", "__dlpack__", "__arrow_c_array__", "__arrow_c_stream__"]
    return any(hasattr(obj, name) for name in attrs)


def maybe_warn_listlike(other) -> None:
    """
    Warn when operating against a list-like that is neither a standard container
    (``list``, ``np.ndarray``) nor a pandas object nor array-castable.

    Such operations (e.g. with ``tuple``, ``range``, ``deque``) are deprecated
    (GH#62423) and will treat ``other`` as scalar-like in a future version.
    """
    if (
        is_list_like(other)
        and not isinstance(
            other,
            (list, np.ndarray, ABCExtensionArray, ABCIndex, ABCSeries, ABCDataFrame),
        )
        and not has_castable_attr(other)
    ):
        warnings.warn(
            f"Operation with {type(other).__name__} are deprecated. "
            "In a future version these will be treated as scalar-like. "
            "To retain the old behavior, explicitly wrap in a Series "
            "instead.",
            Pandas4Warning,
            stacklevel=find_stack_level(),
        )


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
    is_logical = name.strip("_") in ["or", "xor", "and", "ror", "rxor", "rand"]

    @wraps(method)
    def new_method(self, other):
        prio = getattr(other, "__pandas_priority__", None)
        if prio is not None:
            if prio > self.__pandas_priority__:
                # e.g. other is DataFrame while self is Index/Series/EA
                return NotImplemented

        other = item_from_zerodim(other)
        if (
            isinstance(self, ABCExtensionArray)
            and isinstance(other, list)
            and not is_logical
        ):
            # See GH#62423
            other = sanitize_array(other, None)
            other = ensure_wrapped_if_datetimelike(other)

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
