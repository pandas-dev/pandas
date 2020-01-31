"""
Functions for defining unary operations.
"""
from typing import Any

import numpy as np

from pandas._typing import ArrayLike

from pandas.core.dtypes.common import (
    is_datetime64_dtype,
    is_extension_array_dtype,
    is_integer_dtype,
    is_object_dtype,
    is_scalar,
    is_timedelta64_dtype,
)
from pandas.core.dtypes.generic import ABCSeries

from pandas.core.construction import array


def should_extension_dispatch(left: ABCSeries, right: Any) -> bool:
    """
    Identify cases where Series operation should use dispatch_to_extension_op.

    Parameters
    ----------
    left : Series
    right : object

    Returns
    -------
    bool
    """
    if (
        is_extension_array_dtype(left.dtype)
        or is_datetime64_dtype(left.dtype)
        or is_timedelta64_dtype(left.dtype)
    ):
        return True

    if not is_scalar(right) and is_extension_array_dtype(right):
        # GH#22378 disallow scalar to exclude e.g. "category", "Int64"
        return True

    return False


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


def dispatch_to_extension_op(op, left: ArrayLike, right: Any):
    """
    Assume that left or right is a Series backed by an ExtensionArray,
    apply the operator defined by op.

    Parameters
    ----------
    op : binary operator
    left : ExtensionArray or np.ndarray
    right : object

    Returns
    -------
    ExtensionArray or np.ndarray
        2-tuple of these if op is divmod or rdivmod
    """
    # NB: left and right should already be unboxed, so neither should be
    #  a Series or Index.

    if left.dtype.kind in "mM" and isinstance(left, np.ndarray):
        # We need to cast datetime64 and timedelta64 ndarrays to
        #  DatetimeArray/TimedeltaArray.  But we avoid wrapping others in
        #  PandasArray as that behaves poorly with e.g. IntegerArray.
        left = array(left)

    # The op calls will raise TypeError if the op is not defined
    # on the ExtensionArray
    res_values = op(left, right)
    return res_values
