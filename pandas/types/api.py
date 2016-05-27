# flake8: noqa

import numpy as np
from pandas.compat import string_types

from .dtypes import (CategoricalDtype, CategoricalDtypeType,
                     DatetimeTZDtype, DatetimeTZDtypeType)
from .generic import (ABCIndex, ABCInt64Index, ABCRangeIndex,
                      ABCFloat64Index, ABCMultiIndex,
                      ABCDatetimeIndex,
                      ABCTimedeltaIndex, ABCPeriodIndex,
                      ABCCategoricalIndex,
                      ABCIndexClass,
                      ABCSeries, ABCDataFrame, ABCPanel,
                      ABCSparseSeries, ABCSparseArray,
                      ABCCategorical, ABCPeriod,
                      ABCGeneric)

def pandas_dtype(dtype):
    """
    Converts input into a pandas only dtype object or a numpy dtype object.

    Parameters
    ----------
    dtype : object to be converted

    Returns
    -------
    np.dtype or a pandas dtype
    """
    if isinstance(dtype, DatetimeTZDtype):
        return dtype
    elif isinstance(dtype, CategoricalDtype):
        return dtype
    elif isinstance(dtype, string_types):
        try:
            return DatetimeTZDtype.construct_from_string(dtype)
        except TypeError:
            pass

        try:
            return CategoricalDtype.construct_from_string(dtype)
        except TypeError:
            pass

    return np.dtype(dtype)

def na_value_for_dtype(dtype):
    """
    Return a dtype compat na value

    Parameters
    ----------
    dtype : string / dtype

    Returns
    -------
    dtype compat na value
    """

    from pandas.core import common as com
    from pandas import NaT
    dtype = pandas_dtype(dtype)

    if (com.is_datetime64_dtype(dtype) or
        com.is_datetime64tz_dtype(dtype) or
        com.is_timedelta64_dtype(dtype)):
        return NaT
    elif com.is_float_dtype(dtype):
        return np.nan
    elif com.is_integer_dtype(dtype):
        return 0
    elif com.is_bool_dtype(dtype):
        return False
    return np.nan
