import numpy as np

import pandas as pd
import pandas._testing as tm


def test_data_frame_value_counts_unsorted():
    df = pd.DataFrame(
        {"num_legs": [2, 4, 4, 6], "num_wings": [2, 0, 0, 0]},
        index=["falcon", "dog", "cat", "ant"],
    )

    with tm.assert_produces_warning(FutureWarning, match="In pandas 2.0.0, the name"):
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

    with tm.assert_produces_warning(FutureWarning, match="In pandas 2.0.0, the name"):
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

    with tm.assert_produces_warning(FutureWarning, match="In pandas 2.0.0, the name"):
        result = df.value_counts()
    expected = pd.Series(
        data=[2, 1, 1],
        index=pd.MultiIndex.from_arrays(
            [(4, 2, 6), (0, 2, 0)], names=["num_legs", "num_wings"]
        ),
    )

    tm.assert_series_equal(result, expected)


def test_data_frame_value_counts_normalize():
    df = pd.DataFrame(
        {"num_legs": [2, 4, 4, 6], "num_wings": [2, 0, 0, 0]},
        index=["falcon", "dog", "cat", "ant"],
    )

    with tm.assert_produces_warning(FutureWarning, match="In pandas 2.0.0, the name"):
        result = df.value_counts(normalize=True)
    expected = pd.Series(
        data=[0.5, 0.25, 0.25],
        index=pd.MultiIndex.from_arrays(
            [(4, 2, 6), (0, 2, 0)], names=["num_legs", "num_wings"]
        ),
    )

    tm.assert_series_equal(result, expected)


def test_data_frame_value_counts_single_col_default():
    df = pd.DataFrame({"num_legs": [2, 4, 4, 6]})

    with tm.assert_produces_warning(FutureWarning, match="In pandas 2.0.0, the name"):
        result = df.value_counts()
    expected = pd.Series(
        data=[2, 1, 1],
        index=pd.MultiIndex.from_arrays([[4, 2, 6]], names=["num_legs"]),
    )

    tm.assert_series_equal(result, expected)


def test_data_frame_value_counts_empty():
    df_no_cols = pd.DataFrame()

    with tm.assert_produces_warning(FutureWarning, match="In pandas 2.0.0, the name"):
        result = df_no_cols.value_counts()
    expected = pd.Series([], dtype=np.int64)

    tm.assert_series_equal(result, expected)


def test_data_frame_value_counts_empty_normalize():
    df_no_cols = pd.DataFrame()

    with tm.assert_produces_warning(FutureWarning, match="In pandas 2.0.0, the name"):
        result = df_no_cols.value_counts(normalize=True)
    expected = pd.Series([], dtype=np.float64)

    tm.assert_series_equal(result, expected)


def test_data_frame_value_counts_dropna_true(nulls_fixture):
    # GH 41334
    df = pd.DataFrame(
        {
            "first_name": ["John", "Anne", "John", "Beth"],
            "middle_name": ["Smith", nulls_fixture, nulls_fixture, "Louise"],
        },
    )
    with tm.assert_produces_warning(FutureWarning, match="In pandas 2.0.0, the name"):
        result = df.value_counts()
    expected = pd.Series(
        data=[1, 1],
        index=pd.MultiIndex.from_arrays(
            [("Beth", "John"), ("Louise", "Smith")], names=["first_name", "middle_name"]
        ),
    )

    tm.assert_series_equal(result, expected)


def test_data_frame_value_counts_dropna_false(nulls_fixture):
    # GH 41334
    df = pd.DataFrame(
        {
            "first_name": ["John", "Anne", "John", "Beth"],
            "middle_name": ["Smith", nulls_fixture, nulls_fixture, "Louise"],
        },
    )

    with tm.assert_produces_warning(FutureWarning, match="In pandas 2.0.0, the name"):
        result = df.value_counts(dropna=False)
    expected = pd.Series(
        data=[1, 1, 1, 1],
        index=pd.MultiIndex(
            levels=[
                pd.Index(["Anne", "Beth", "John"]),
                pd.Index(["Louise", "Smith", nulls_fixture]),
            ],
            codes=[[0, 1, 2, 2], [2, 0, 1, 2]],
            names=["first_name", "middle_name"],
        ),
    )

    tm.assert_series_equal(result, expected)


def test_value_counts_name():
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    result = df.value_counts(name="counts")
    expected_idx = pd.MultiIndex.from_arrays([[1, 2, 3], [4, 5, 6]], names=["a", "b"])
    expected = pd.Series([1, 1, 1], name="counts", index=expected_idx)
    tm.assert_series_equal(result, expected)
