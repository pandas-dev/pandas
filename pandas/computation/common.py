import numpy as np
import pandas as pd
from pandas.compat import reduce


def _ensure_decoded(s):
    """ if we have bytes, decode them to unicode """
    if isinstance(s, (np.bytes_, bytes)):
        s = s.decode(pd.get_option('display.encoding'))
    return s


def _result_type_many(*arrays_and_dtypes):
    """ wrapper around numpy.result_type which overcomes the NPY_MAXARGS (32)
    argument limit """
    try:
        return np.result_type(*arrays_and_dtypes)
    except ValueError:
        # length 0 or length > NPY_MAXARGS both throw a ValueError, so check
        # which one we're dealing with
        if len(arrays_and_dtypes) == 0:
            raise ValueError('at least one array or dtype is required')
        return reduce(np.result_type, arrays_and_dtypes)


class NameResolutionError(NameError):
    pass
