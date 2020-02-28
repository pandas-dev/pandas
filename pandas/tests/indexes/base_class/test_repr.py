import numpy as np

from pandas import Index


class TestIndexRepr:
    # Tests for the Index representation,
    # specifically for the case that includes bools and NANs

    def test_index_repr_bool_nan(self):
        # GH32146
        arr = Index([True, False, np.nan], dtype=object)
        exp = arr.format()
        out = ["True", "False", "NaN"]
        assert out == exp
