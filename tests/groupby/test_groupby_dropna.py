import numpy as np
import pandas as pd
import pandas._testing as tm

def test_groupby_dataframe_dropna_false_preserves_nan_group():
    # Ensure DataFrame.groupby(..., dropna=False) preserves NA entries as a single group
    # Tests-only addition to lock current behavior (GHxxxx)
    data = {"group": ["g1", np.nan, "g1", "g2", np.nan], "val": [0, 1, 2, 3, 4]}
    df = pd.DataFrame(data)

    gb_keepna = df.groupby("group", dropna=False)
    result = gb_keepna.indices

    # expected: g1 -> [0,2], g2 -> [3], NaN -> [1,4]
    expected = {
        "g1": np.array([0, 2], dtype=np.intp),
        "g2": np.array([3], dtype=np.intp),
        np.nan: np.array([1, 4], dtype=np.intp),
    }

    # Compare group indices allowing for np.nan key
    for res_vals, exp_vals in zip(result.values(), expected.values()):
        tm.assert_numpy_array_equal(res_vals, exp_vals)
    # check there is an NaN key present
    assert any(pd.isna(k) for k in result.keys())


def test_groupby_series_dropna_false_preserves_nan_group():
    # Verify Series.groupby(..., dropna=False) also preserves NA groups
    s = pd.Series([1, 2, 3, 4], index=["a", np.nan, "a", np.nan], name="s")
    gb = s.groupby(level=0, dropna=False)
    res = gb.indices

    expected = {
        "a": np.array([0, 2], dtype=np.intp),
        np.nan: np.array([1, 3], dtype=np.intp),
    }

    for res_vals, exp_vals in zip(res.values(), expected.values()):
        tm.assert_numpy_array_equal(res_vals, exp_vals)
    assert any(pd.isna(k) for k in res.keys())
