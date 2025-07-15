import numpy as np
import pandas as pd
import pytest

@pytest.mark.parametrize("op", [lambda x, y: x & y, lambda x, y: x | y])
def test_logical_op_2d_extension_array(op):
    df = pd.DataFrame(np.arange(50).reshape(10, 5)).notna().values

    EA_array = pd.array([i for i in range(10)], dtype="Int64").reshape(10, 1)
    NP_array = np.arange(10).reshape(10, 1)

    expected = op(df, NP_array)
    result = op(df, EA_array)

    assert np.array_equal(result, expected)
