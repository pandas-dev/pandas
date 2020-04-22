import numpy as np

import pandas as pd
import pandas._testing as tm
from pandas.core.construction import create_series_with_explicit_index


def test_data_frame_value_counts_unsorted():
    df = pd.DataFrame(
        {"num_legs": [2, 4, 4, 6], "num_wings": [2, 0, 0, 0]},
        index=["falcon", "dog", "cat", "ant"],
    )

    result = df.value_counts(sort=False)
    expected = pd.Series(
        data=[1, 2, 1],
        index=pd.MultiIndex.from_arrays(
            [(2, 4, 6), (2, 0, 0)], names=["num_legs", "num_wings"]
        ),
    )

    tm.assert_series_equal(result, expected)


def test_data_frame_value_counts_ascending():
    df = pd.DataFrame(
        {"num_legs": [2, 4, 4, 6], "num_wings": [2, 0, 0, 0]},
        index=["falcon", "dog", "cat", "ant"],
    )

    result = df.value_counts(ascending=True)
    expected = pd.Series(
        data=[1, 1, 2],
        index=pd.MultiIndex.from_arrays(
            [(2, 6, 4), (2, 0, 0)], names=["num_legs", "num_wings"]
        ),
    )

    tm.assert_series_equal(result, expected)


def test_data_frame_value_counts_default():
    df = pd.DataFrame(
        {"num_legs": [2, 4, 4, 6], "num_wings": [2, 0, 0, 0]},
        index=["falcon", "dog", "cat", "ant"],
    )

    result = df.value_counts()
    expected = pd.Series(
        data=[2, 1, 1],
        index=pd.MultiIndex.from_arrays(
            [(4, 6, 2), (0, 0, 2)], names=["num_legs", "num_wings"]
        ),
    )

    tm.assert_series_equal(result, expected)


def test_data_frame_value_counts_normalize():
    df = pd.DataFrame(
        {"num_legs": [2, 4, 4, 6], "num_wings": [2, 0, 0, 0]},
        index=["falcon", "dog", "cat", "ant"],
    )

    result = df.value_counts(normalize=True)
    expected = pd.Series(
        data=[0.5, 0.25, 0.25],
        index=pd.MultiIndex.from_arrays(
            [(4, 6, 2), (0, 0, 2)], names=["num_legs", "num_wings"]
        ),
    )

    tm.assert_series_equal(result, expected)


def test_data_frame_value_counts_single_col_default():
    df = pd.DataFrame({"num_legs": [2, 4, 4, 6]})

    result = df.value_counts()
    expected = pd.Series(
        data=[2, 1, 1],
        index=pd.MultiIndex.from_arrays([[4, 6, 2]], names=["num_legs"]),
    )

    tm.assert_series_equal(result, expected)


def test_data_frame_value_counts_empty():
    df_no_cols = pd.DataFrame()

    result = df_no_cols.value_counts()
    expected = create_series_with_explicit_index([], dtype=np.int64)

    tm.assert_series_equal(result, expected)


def test_data_frame_value_counts_empty_normalize():
    df_no_cols = pd.DataFrame()

    result = df_no_cols.value_counts(normalize=True)
    expected = create_series_with_explicit_index([], dtype=np.float64)

    tm.assert_series_equal(result, expected)
