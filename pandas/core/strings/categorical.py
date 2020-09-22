from typing import cast

import numpy as np

from pandas.core.algorithms import take_1d
from pandas.core.strings.object_array import ObjectStringArray
from pandas.core.arrays import Categorical


class CategoricalStringMethods(Categorical, ObjectStringArray):
    """
    Extension array implementing _str methods for Categorical.

    We implement this just to avoid mucking up the inheritance chain
    for Categorical. This inherits from ObjectStringArray just for
    convenience.
    """
    # Probably lots of room for improvement here.
    def _map(self, f, na_value=np.nan, dtype=np.dtype(object)):
        from pandas import Categorical

        arr = cast(Categorical, self._array)

        categories = arr.categories
        codes = arr.codes
        result = ObjectStringArray(categories)._map(f, na_value, dtype)
        return take_1d(result, codes, fill_value=na_value)

    def _str_get_dummies(self, sep="|"):
        # sep may not be in categories. Just bail on this.
        return ObjectStringArray(self.astype(str))._str_get_dummies(sep)
