"""
Boilerplate functions used in defining binary operations.
"""
from functools import wraps
from typing import Any, Callable

import numpy as np

from pandas._libs.lib import item_from_zerodim
from pandas._libs.ops_dispatch import maybe_dispatch_ufunc_to_dunder_op
from pandas._typing import F

from pandas.core.dtypes.generic import ABCDataFrame, ABCIndexClass, ABCSeries


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


def _unpack_zerodim_and_defer(method, name: str):
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
    is_cmp = name.strip("__") in {"eq", "ne", "lt", "le", "gt", "ge"}

    @wraps(method)
    def new_method(self, other):

        if is_cmp and isinstance(self, ABCIndexClass) and isinstance(other, ABCSeries):
            # For comparison ops, Index does *not* defer to Series
            pass
        else:
            for cls in [ABCDataFrame, ABCSeries, ABCIndexClass]:
                if isinstance(self, cls):
                    break
                if isinstance(other, cls):
                    return NotImplemented

        other = item_from_zerodim(other)

        return method(self, other)

    return new_method


def defer_or_dispatch_ufunc(meth):
    """
    Boilerplate for pandas conventions in arithmetic and comparison methods.

    Ensure method returns NotImplemented when operating against "senior"
    classes.  Ensure zero-dimensional ndarrays are always unpacked.

    Parameters
    ----------
    method : binary method

    Returns
    -------
    method
    """

    @wraps(meth)
    def new_method(self, ufunc: Callable, method: str, *inputs: Any, **kwargs: Any):
        cls = type(self)

        # for binary ops, use our custom dunder methods
        result = maybe_dispatch_ufunc_to_dunder_op(
            self, ufunc, method, *inputs, **kwargs
        )
        if result is not NotImplemented:
            return result

        # Determine if we should defer.
        no_defer = (np.ndarray.__array_ufunc__, cls.__array_ufunc__)

        for item in inputs:
            higher_priority = (
                hasattr(item, "__array_priority__")
                and item.__array_priority__ > self.__array_priority__
            )
            has_array_ufunc = (
                hasattr(item, "__array_ufunc__")
                and type(item).__array_ufunc__ not in no_defer
                and not isinstance(item, self._HANDLED_TYPES)
            )
            if higher_priority or has_array_ufunc:
                return NotImplemented

        return meth(self, ufunc, method, *inputs, **kwargs)

    return new_method
