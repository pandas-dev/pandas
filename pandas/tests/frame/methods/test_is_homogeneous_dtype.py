import numpy as np
import pytest

import pandas as pd


@pytest.mark.parametrize(
    "data, expected",
    [
        # empty
        (pd.DataFrame(), True),
        # multi-same
        (pd.DataFrame({"A": [1, 2], "B": [1, 2]}), True),
        # multi-object
        (
            pd.DataFrame(
                {
                    "A": np.array([1, 2], dtype=object),
                    "B": np.array(["a", "b"], dtype=object),
                },
                dtype="object",
            ),
            True,
        ),
        # multi-extension
        (
            pd.DataFrame(
                {"A": pd.Categorical(["a", "b"]), "B": pd.Categorical(["a", "b"])}
            ),
            True,
        ),
        # differ types
        (pd.DataFrame({"A": [1, 2], "B": [1.0, 2.0]}), False),
        # differ sizes
        (
            pd.DataFrame(
                {
                    "A": np.array([1, 2], dtype=np.int32),
                    "B": np.array([1, 2], dtype=np.int64),
                }
            ),
            False,
        ),
        # multi-extension differ
        (
            pd.DataFrame(
                {"A": pd.Categorical(["a", "b"]), "B": pd.Categorical(["b", "c"])}
            ),
            False,
        ),
    ],
)
def test_is_homogeneous_type(data, expected):
    assert data._is_homogeneous_type is expected
