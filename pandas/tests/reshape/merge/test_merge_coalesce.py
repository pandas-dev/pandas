import numpy as np

import pandas as pd
import pandas._testing as tm


def test_merge_with_and_without_coalesce_keys():
    # Create two example DataFrames to merge
    df1 = pd.DataFrame({"id": [1, 2, 3], "value1": ["A", "B", "C"]})

    df2 = pd.DataFrame({"id": [2, 3, 4], "value2": ["X", "Y", "Z"]})

    # Merge with coalesce_keys=True (default behavior)
    result_coalesced = pd.merge(
        df1,
        df2,
        on="id",
        how="outer",
        coalesce_keys=True,
        suffixes=("", "_right"),
    )

    # Expected result: keys from both sides are coalesced into a single 'id' column
    expected_coalesced = pd.DataFrame(
        {
            "id": [1.0, 2.0, 3.0, 4.0],
            "value1": ["A", "B", "C", np.nan],
            "value2": [np.nan, "X", "Y", "Z"],
        }
    )

    # Use assert_frame_equal to validate the result
    tm.assert_frame_equal(
        result_coalesced.sort_index(axis=1),
        expected_coalesced.sort_index(axis=1),
        check_dtype=False,
    )

    # Merge with coalesce_keys=False (new functionality)
    result_separated = pd.merge(
        df1,
        df2,
        on="id",
        how="outer",
        coalesce_keys=False,
        suffixes=("", "_right"),
    )

    # Expected result: original left 'id' and right 'id' (renamed to 'id_right')
    # are preserved
    expected_separated = pd.DataFrame(
        {
            "id": [1.0, 2.0, 3.0, np.nan],
            "value1": ["A", "B", "C", np.nan],
            "id_right": [np.nan, 2.0, 3.0, 4.0],
            "value2": [np.nan, "X", "Y", "Z"],
        }
    )

    # Check that no columns have value None (should be np.nan)
    assert None not in result_separated.columns

    # Check shape of the result
    assert result_separated.shape == (4, 4)

    # Validate content
    tm.assert_frame_equal(
        result_separated.sort_index(axis=1),
        expected_separated.sort_index(axis=1),
        check_dtype=False,
    )
