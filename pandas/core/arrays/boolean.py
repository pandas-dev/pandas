import numpy as np

from pandas.core.arrays.integer import IntegerArray, _IntegerDtype
from pandas.core.dtypes.dtypes import registry
from pandas.util._decorators import cache_readonly


class BooleanDtype(_IntegerDtype):
    name = "EABool"
    type = np.bool
    na_value = np.nan

    @cache_readonly
    def is_signed_integer(self):
        return False

    @cache_readonly
    def is_unsigned_integer(self):
        return False

    @classmethod
    def construct_array_type(cls):
        return BooleanArray


class BooleanArray(IntegerArray):

    @cache_readonly
    def dtype(self):
        return BooleanDtype()


def to_boolean_array(values):
    """
    Infer and return an integer array of the values.

    Parameters
    ----------
    values : 1D list-like

    Returns
    -------
    IntegerArray

    Raises
    ------
    TypeError if incompatible types
    """
    return BooleanArray(values, dtype='uint8', copy=False)


registry.register(BooleanDtype)
