import numpy as np

import pandas as pd


def test_size():
    # test for Series object
    s = pd.Series({"a": 1, "b": 2, "c": 3})
    assert s.size == 3
    assert isinstance(s.size, int)

    # test for DataFrame object
    df = pd.DataFrame({"col1": [1, 2], "col2": [3, 4]})
    assert df.size == 4
    assert isinstance(df.size, int)

    # test for empty DataFrame object
    empty_df = pd.DataFrame()
    assert empty_df.size == 0
    assert isinstance(empty_df.size, int)

    # test for DataFrame with missing values
    df_with_missing = pd.DataFrame({"col1": [1, np.nan], "col2": [3, 4]})
    assert df_with_missing.size == 4
    assert isinstance(df_with_missing.size, int)

    # test for MultiIndex DataFrame
    multi_df = pd.DataFrame(
        {"col1": [1, 2], "col2": [3, 4]}, index=[["a", "b"], [1, 2]]
    )
    assert multi_df.size == 4
    assert isinstance(multi_df.size, int)
