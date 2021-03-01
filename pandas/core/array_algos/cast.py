import inspect

import numpy as np

from pandas._typing import (
    ArrayLike,
    Dtype,
    DtypeObj,
)

from pandas.core.dtypes.cast import (
    astype_dt64_to_dt64tz,
    astype_nansafe,
)
from pandas.core.dtypes.common import (
    is_datetime64_dtype,
    is_datetime64tz_dtype,
    is_dtype_equal,
    pandas_dtype,
)
from pandas.core.dtypes.dtypes import ExtensionDtype

from pandas.core.arrays import ExtensionArray


def astype_array(values: ArrayLike, dtype: DtypeObj, copy: bool = False):
    """
    Cast array to the new dtype.

    Parameters
    ----------
    values : ndarray or ExtensionArray
    dtype : str, dtype convertible
    copy : bool, default False
        copy if indicated

    Returns
    -------
    ndarray or ExtensionArray
    """
    if (
        values.dtype.kind in ["m", "M"]
        and dtype.kind in ["i", "u"]
        and isinstance(dtype, np.dtype)
        and dtype.itemsize != 8
    ):
        # TODO(2.0) remove special case once deprecation on DTA/TDA is enforced
        msg = rf"cannot astype a datetimelike from [{values.dtype}] to [{dtype}]"
        raise TypeError(msg)

    if is_datetime64tz_dtype(dtype) and is_datetime64_dtype(values.dtype):
        return astype_dt64_to_dt64tz(values, dtype, copy, via_utc=True)

    if is_dtype_equal(values.dtype, dtype):
        if copy:
            return values.copy()
        return values

    if isinstance(values, ExtensionArray):
        values = values.astype(dtype, copy=copy)

    else:
        values = astype_nansafe(values, dtype, copy=copy)

    # now in ObjectBlock._maybe_coerce_values(cls, values):
    if isinstance(dtype, np.dtype) and issubclass(values.dtype.type, str):
        values = np.array(values, dtype=object)

    return values


def astype_array_safe(values, dtype: Dtype, copy: bool = False, errors: str = "raise"):
    """
    Cast array to the new dtype.

    Parameters
    ----------
    values : ndarray or ExtensionArray
    dtype : str, dtype convertible
    copy : bool, default False
        copy if indicated
    errors : str, {'raise', 'ignore'}, default 'raise'
        - ``raise`` : allow exceptions to be raised
        - ``ignore`` : suppress exceptions. On error return original object

    Returns
    -------
    ndarray or ExtensionArray
    """
    errors_legal_values = ("raise", "ignore")

    if errors not in errors_legal_values:
        invalid_arg = (
            "Expected value of kwarg 'errors' to be one of "
            f"{list(errors_legal_values)}. Supplied value is '{errors}'"
        )
        raise ValueError(invalid_arg)

    if inspect.isclass(dtype) and issubclass(dtype, ExtensionDtype):
        msg = (
            f"Expected an instance of {dtype.__name__}, "
            "but got the class instead. Try instantiating 'dtype'."
        )
        raise TypeError(msg)

    dtype = pandas_dtype(dtype)

    try:
        new_values = astype_array(values, dtype, copy=copy)
    except (ValueError, TypeError):
        # e.g. astype_nansafe can fail on object-dtype of strings
        #  trying to convert to float
        if errors == "ignore":
            new_values = values
        else:
            raise

    return new_values
