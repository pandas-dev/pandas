"""
Boilerplate functions used in defining binary operations.
"""
from functools import wraps

from pandas._libs.lib import item_from_zerodim

from pandas.core.dtypes.generic import ABCDataFrame, ABCIndexClass, ABCSeries


def unpack_zerodim_and_defer(name: str):
    """
    Boilerplate for pandas conventions in arithmetic and comparison methods.

    Parameters
    ----------
    name : str

    Returns
    -------
    decorator
    """

    def wrapper(method):
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
