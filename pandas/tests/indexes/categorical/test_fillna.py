import numpy as np
import pytest

from pandas import CategoricalIndex
import pandas._testing as tm


class TestFillNA:
    def test_fillna_categorical(self):
        # GH#11343
        idx = CategoricalIndex([1.0, np.nan, 3.0, 1.0], name="x")
        # fill by value in categories
        exp = CategoricalIndex([1.0, 1.0, 3.0, 1.0], name="x")
        tm.assert_index_equal(idx.fillna(1.0), exp)

        # fill by value not in categories raises ValueError
        msg = "fill value must be in categories"
        with pytest.raises(ValueError, match=msg):
            idx.fillna(2.0)
