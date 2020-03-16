import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "msg,labels,level",
    [
        (r"labels \[4\] not found in level", 4, "a"),
        (r"labels \[7\] not found in level", 7, "b"),
    ],
)
def test_drop_raise_exception_if_labels_not_in_level(msg, labels, level):
    # GH 8594
    mi = pd.MultiIndex.from_arrays([[1, 2, 3], [4, 5, 6]], names=["a", "b"])
    s = pd.Series([10, 20, 30], index=mi)
    df = pd.DataFrame([10, 20, 30], index=mi)

    with pytest.raises(KeyError, match=msg):
        s.drop(labels, level=level)
    with pytest.raises(KeyError, match=msg):
        df.drop(labels, level=level)


@pytest.mark.parametrize("labels,level", [(4, "a"), (7, "b")])
def test_drop_errors_ignore(labels, level):
    # GH 8594
    mi = pd.MultiIndex.from_arrays([[1, 2, 3], [4, 5, 6]], names=["a", "b"])
    s = pd.Series([10, 20, 30], index=mi)
    df = pd.DataFrame([10, 20, 30], index=mi)

    expected_s = s.drop(labels, level=level, errors="ignore")
    tm.assert_series_equal(s, expected_s)

    expected_df = df.drop(labels, level=level, errors="ignore")
    tm.assert_frame_equal(df, expected_df)


def test_drop_with_non_unique_datetime_index_and_invalid_keys():
    # GH 30399

    # define dataframe with unique datetime index
    df = pd.DataFrame(
        np.random.randn(5, 3),
        columns=["a", "b", "c"],
        index=pd.date_range("2012", freq="H", periods=5),
    )
    # create dataframe with non-unique datetime index
    df = df.iloc[[0, 2, 2, 3]].copy()

    with pytest.raises(KeyError, match="not found in axis"):
        df.drop(["a", "b"])  # Dropping with labels not exist in the index
