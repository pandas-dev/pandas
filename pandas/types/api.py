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
    if isinstance(dtype, string_types):
        try:
            return DatetimeTZDtype.construct_from_string(dtype)
        except TypeError:
            pass

        try:
            return CategoricalDtype.construct_from_string(dtype)
        except TypeError:
            pass

    return np.dtype(dtype)
