"""
Tests involving custom Index subclasses
"""

import numpy as np

import pandas as pd
import pandas._testing as tm


class CustomIndex(pd.Index):
    def __new__(cls, data, name=None):
        # assert that this index class cannot hold strings
        if any(isinstance(val, str) for val in data):
            raise TypeError("CustomIndex cannot hold strings")

        if name is None and hasattr(data, "name"):
            name = data.name
        data = np.array(data, dtype="O")

        return cls._simple_new(data, name)


def test_insert_fallback_to_base_index():
    # https://github.com/pandas-dev/pandas/issues/47071

    idx = CustomIndex([1, 2, 3])
    result = idx.insert(0, "string")
    expected = pd.Index(["string", 1, 2, 3], dtype=object)
    tm.assert_index_equal(result, expected)

    df = pd.DataFrame(
        np.random.default_rng(2).standard_normal((2, 3)),
        columns=idx,
        index=pd.Index([1, 2], name="string"),
    )
    result = df.reset_index()
    tm.assert_index_equal(result.columns, expected)
