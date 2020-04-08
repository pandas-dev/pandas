"""
Functions for defining unary operations.
"""
from typing import Any

from pandas._typing import ArrayLike

from pandas.core.dtypes.common import (
    is_datetime64_dtype,
    is_integer_dtype,
    is_object_dtype,
    is_timedelta64_dtype,
)
from pandas.core.dtypes.generic import ABCExtensionArray


def should_extension_dispatch(left: ArrayLike, right: Any) -> bool:
    """
    Identify cases where Series operation should dispatch to ExtensionArray method.

    Parameters
    ----------
    left : np.ndarray or ExtensionArray
    right : object

    Returns
    -------
    bool
    """
    return isinstance(left, ABCExtensionArray) or isinstance(right, ABCExtensionArray)


def should_series_dispatch(left, right, op):
    """
    Identify cases where a DataFrame operation should dispatch to its
    Series counterpart.

    Parameters
    ----------
    left : DataFrame
    right : DataFrame or Series
    op : binary operator

    Returns
    -------
    override : bool
    """
    if left._is_mixed_type or right._is_mixed_type:
        return True

    if op.__name__.strip("_") in ["and", "or", "xor", "rand", "ror", "rxor"]:
        # TODO: GH references for what this fixes
        # Note: this check must come before the check for nonempty columns.
        return True

    if right.ndim == 1:
        # operating with Series, short-circuit checks that would fail
        #  with AttributeError.
        return False

    if not len(left.columns) or not len(right.columns):
        # ensure obj.dtypes[0] exists for each obj
        return False

    ldtype = left.dtypes.iloc[0]
    rdtype = right.dtypes.iloc[0]

    if (is_timedelta64_dtype(ldtype) and is_integer_dtype(rdtype)) or (
        is_timedelta64_dtype(rdtype) and is_integer_dtype(ldtype)
    ):
        # numpy integer dtypes as timedelta64 dtypes in this scenario
        return True

    if is_datetime64_dtype(ldtype) and is_object_dtype(rdtype):
        # in particular case where right is an array of DateOffsets
        return True

    return False
