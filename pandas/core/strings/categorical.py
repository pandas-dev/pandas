from typing import cast

import numpy as np

from pandas.core.algorithms import take_1d
from pandas.core.strings.object_array import ObjectArrayMethods


class CategoricalStringMethods(ObjectArrayMethods):
    def _map(self, f, na_value=np.nan, dtype=np.dtype(object)):
        from pandas import Categorical

        arr = cast(Categorical, self._array)

        categories = arr.categories
        codes = arr.codes
        result = ObjectArrayMethods(categories)._map(f, na_value, dtype)
        return take_1d(result, codes, fill_value=na_value)

    def get_dummies(self, sep="|"):
        # sep may not be in categories. Just bail on this.
        return ObjectArrayMethods(self._array.astype(str)).get_dummies(sep)
