import numpy as np

from pandas import (
    Categorical,
    CategoricalDtype,
    array,
)
import pandas._testing as tm


class TestAstype:
    def test_astype_str_int_categories_to_nullable_int(self):
        # GH#39616
        dtype = CategoricalDtype([str(i) for i in range(5)])
        arr = Categorical.from_codes(np.random.randint(5, size=20), dtype=dtype)

        res = arr.astype("Int64")
        expected = array(arr.astype("int64"), dtype="Int64")
        tm.assert_extension_array_equal(res, expected)
