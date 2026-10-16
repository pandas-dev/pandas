import numpy as np
import pytest

import pandas as pd


@pytest.mark.parametrize(
    "data,exp_size",
    [
        # see gh-16362.
        ([[pd.NaT, "a", "b", 0], [pd.NaT, "b", "c", 1]], 8),
        ([[pd.NaT, "a", 0], [pd.NaT, "b", 1]], 6),
    ],
)
def test_maybe_infer_to_datetimelike_df_construct(data, exp_size):
    result = pd.DataFrame(np.array(data))
    assert result.size == exp_size


def test_maybe_infer_to_datetimelike_ser_construct():
    # see gh-19671.
    result = pd.Series(["M1701", pd.Timestamp("20130101")])
    assert result.dtype.kind == "O"
