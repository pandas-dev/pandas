import numpy as np
import pytest

from pandas import MultiIndex, Series


@pytest.mark.parametrize(
    "arr, idx",
    [
        ([1, 2, 3, 4], [0, 2, 1, 3]),
        ([1, np.nan, 3, np.nan], [0, 2, 1, 3]),
        (
            [1, np.nan, 3, np.nan],
            MultiIndex.from_tuples([(0, "a"), (1, "b"), (2, "c"), (3, "c")]),
        ),
    ],
)
def test_equals(arr, idx):
    s1 = Series(arr, index=idx)
    s2 = s1.copy()
    assert s1.equals(s2)

    s1[1] = 9
    assert not s1.equals(s2)


def test_equals_list_array():
    # GH20676 Verify equals operator for list of Numpy arrays
    arr = np.array([1, 2])
    s1 = Series([arr, arr])
    s2 = s1.copy()
    assert s1.equals(s2)

    # TODO: Series equals should also work between single value and list
    # s1[1] = 9
    # assert not s1.equals(s2)


def test_equals_false_negative():
    # GH8437 Verify false negative behavior of equals function for dtype object
    arr = [False, np.nan]
    s1 = Series(arr)
    s2 = s1.copy()
    s3 = Series(index=range(2), dtype=object)
    s4 = s3.copy()
    s5 = s3.copy()
    s6 = s3.copy()

    s3[:-1] = s4[:-1] = s5[0] = s6[0] = False
    assert s1.equals(s1)
    assert s1.equals(s2)
    assert s1.equals(s3)
    assert s1.equals(s4)
    assert s1.equals(s5)
    assert s5.equals(s6)
