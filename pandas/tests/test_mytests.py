import re

import numpy as np
import pytest

from pandas import DataFrame
import pandas._testing as tm


def test_dup_with_mismatched_column_lengths():
    # dup across dtypes
    df = DataFrame(
        np.random.randn(3, 4),
        columns=["foo", "bar", "foo", "hello"],
    )
    # check(df)

    msg = "Errored123"
    with pytest.raises(ValueError, match=msg):
        df["foo"] = np.random.randn(3, 3)

    df["foo"] = np.random.randn(3, 2)
    df["foo"] = np.random.randn(3, 1)


@pytest.mark.parametrize(
    argnames="arr",
    argvalues=[(4, 3), (4, 4), (4, 10), (4, 20), (4, 30)],
)
def test_setitem_size_incompatible_ndarray(arr):
    # GH#40827
    # Assigning a dataframe column to an ndarray with more than one columns
    # should raise an exception.
    data = DataFrame(np.zeros((4, 2)), columns=["A", "B"])
    msg = ""
    msg = "Errored123"
    with pytest.raises(ValueError, match=msg):
        data["A"] = np.random.randn(arr[0], arr[1])


# @pytest.mark.parametrize(
#     argnames="arr",
#     argvalues=[
#         # (4, 3),
#          (4, 4), (4, 10), (4, 20), (4, 30)
#     ],
# )
# def test_setitem_size_incompatible_ndarray2(arr):
#     # GH#40827
#     # Assigning a dataframe column to an ndarray with more than one columns
#     # should raise an exception.
#     data = DataFrame(
#         [[1, "A", 1.0], [2, "B", 2.0], [3, "C", 3.0], [4, "D", 4.0]],
#         columns=["A", "A", "B"],
#     )
#     msg = "Errored123"
#     with pytest.raises(ValueError, match=msg):
#         data["A"] = np.random.randn(arr[0], arr[1])

#     data = DataFrame(
#         [[1, 1.0], [2, 2.0], [3, 3.0], [4, 4.0]],
#         columns=["A", "B"],
#     )
#     msg = "Errored123"
#     # TODO: This should fail with pytest.raises()
#     with pytest.raises(ValueError, match=msg):
#         data["A"] = np.random.randn(arr[0], arr[1])


def test_setitem_size_compatible_ndarray():
    data = DataFrame(np.zeros(4), columns=["A"])

    data_to_set = np.random.randn(4, 1)

    expected = DataFrame(data=data_to_set, columns=["A"])

    data["A"] = data_to_set
    tm.assert_frame_equal(data, expected)


@pytest.mark.parametrize(
    argnames="arr", argvalues=[(4, 2), (4, 3), (4, 4), (4, 10), (4, 20), (4, 30)],
)
def test_loc_setitem_with_size_incompatible_ndarray(arr):
    # GH#40827
    # Assigning a dataframe column to an ndarray with more than one columns
    # should raise an exception.
    data = DataFrame(np.zeros(4), columns=["A"])
    msg = re.escape(
        f"could not broadcast input array from shape (4,{arr[1]}) into shape (4,1)"
    )
    with pytest.raises(Exception, match=msg):
        data.iloc[:] = np.random.randn(arr[0], arr[1])


def test_loc_setitem_with_size_compatible_ndarray():
    data = DataFrame(np.zeros(4), columns=["A"])

    data_to_set = np.random.randn(4, 1)

    expected = DataFrame(data=data_to_set, columns=["A"])

    data.iloc[:] = data_to_set
    tm.assert_frame_equal(data, expected)


def test_iloc_setitem_with_size_compatible_ndarray():
    data = DataFrame(np.zeros(4), columns=["A"])

    data_to_set = np.random.randn(4, 1)

    expected = DataFrame(data=data_to_set, columns=["A"])

    data.iloc[:] = data_to_set
    tm.assert_frame_equal(data, expected)


@pytest.mark.parametrize(
    argnames="arr", argvalues=[(4, 2), (4, 3), (4, 4), (4, 10), (4, 20), (4, 30)]
)
def test_iloc_setitem_with_size_incompatible_ndarray(arr):
    # GH#40827
    # Assigning a dataframe column to an ndarray with more than one columns
    # should raise an exception.
    data = DataFrame(np.zeros(4), columns=["A"])
    msg = re.escape(
        f"could not broadcast input array from shape (4,{arr[1]}) " "into shape (4,1)"
    )
    with pytest.raises(Exception, match=msg):
        data.iloc[:] = np.random.randn(arr[0], arr[1])
