import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "input_dict,expected",
    [
        ({0: 0}, np.array([[0]], dtype=int)),
        ({"a": "a"}, np.array([["a"]], dtype=object)),
        ({1: 1}, np.array([[1]], dtype=int)),
    ],
)
def test_numpy_array(input_dict, expected):
    result = np.array([pd.Series(input_dict)])
    tm.assert_numpy_array_equal(result, expected)
