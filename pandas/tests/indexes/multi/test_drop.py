import numpy as np
import pytest

from pandas.errors import PerformanceWarning

import pandas as pd
from pandas import Index, MultiIndex
import pandas._testing as tm


def test_drop(idx):
    dropped = idx.drop([("foo", "two"), ("qux", "one")])

    index = MultiIndex.from_tuples([("foo", "two"), ("qux", "one")])
    dropped2 = idx.drop(index)

    expected = idx[[0, 2, 3, 5]]
    tm.assert_index_equal(dropped, expected)
    tm.assert_index_equal(dropped2, expected)

    dropped = idx.drop(["bar"])
    expected = idx[[0, 1, 3, 4, 5]]
    tm.assert_index_equal(dropped, expected)

    dropped = idx.drop("foo")
    expected = idx[[2, 3, 4, 5]]
    tm.assert_index_equal(dropped, expected)

    index = MultiIndex.from_tuples([("bar", "two")])
    with pytest.raises(KeyError, match=r"^10$"):
        idx.drop([("bar", "two")])
    with pytest.raises(KeyError, match=r"^10$"):
        idx.drop(index)
    with pytest.raises(KeyError, match=r"^'two'$"):
        idx.drop(["foo", "two"])

    # partially correct argument
    mixed_index = MultiIndex.from_tuples([("qux", "one"), ("bar", "two")])
    with pytest.raises(KeyError, match=r"^10$"):
        idx.drop(mixed_index)

    # error='ignore'
    dropped = idx.drop(index, errors="ignore")
    expected = idx[[0, 1, 2, 3, 4, 5]]
    tm.assert_index_equal(dropped, expected)

    dropped = idx.drop(mixed_index, errors="ignore")
    expected = idx[[0, 1, 2, 3, 5]]
    tm.assert_index_equal(dropped, expected)

    dropped = idx.drop(["foo", "two"], errors="ignore")
    expected = idx[[2, 3, 4, 5]]
    tm.assert_index_equal(dropped, expected)

    # mixed partial / full drop
    dropped = idx.drop(["foo", ("qux", "one")])
    expected = idx[[2, 3, 5]]
    tm.assert_index_equal(dropped, expected)

    # mixed partial / full drop / error='ignore'
    mixed_index = ["foo", ("qux", "one"), "two"]
    with pytest.raises(KeyError, match=r"^'two'$"):
        idx.drop(mixed_index)
    dropped = idx.drop(mixed_index, errors="ignore")
    expected = idx[[2, 3, 5]]
    tm.assert_index_equal(dropped, expected)


def test_droplevel_with_names(idx):
    index = idx[idx.get_loc("foo")]
    dropped = index.droplevel(0)
    assert dropped.name == "second"

    index = MultiIndex(
        levels=[Index(range(4)), Index(range(4)), Index(range(4))],
        codes=[
            np.array([0, 0, 1, 2, 2, 2, 3, 3]),
            np.array([0, 1, 0, 0, 0, 1, 0, 1]),
            np.array([1, 0, 1, 1, 0, 0, 1, 0]),
        ],
        names=["one", "two", "three"],
    )
    dropped = index.droplevel(0)
    assert dropped.names == ("two", "three")

    dropped = index.droplevel("two")
    expected = index.droplevel(1)
    assert dropped.equals(expected)


def test_droplevel_list():
    index = MultiIndex(
        levels=[Index(range(4)), Index(range(4)), Index(range(4))],
        codes=[
            np.array([0, 0, 1, 2, 2, 2, 3, 3]),
            np.array([0, 1, 0, 0, 0, 1, 0, 1]),
            np.array([1, 0, 1, 1, 0, 0, 1, 0]),
        ],
        names=["one", "two", "three"],
    )

    dropped = index[:2].droplevel(["three", "one"])
    expected = index[:2].droplevel(2).droplevel(0)
    assert dropped.equals(expected)

    dropped = index[:2].droplevel([])
    expected = index[:2]
    assert dropped.equals(expected)

    msg = (
        "Cannot remove 3 levels from an index with 3 levels: at least one"
        " level must be left"
    )
    with pytest.raises(ValueError, match=msg):
        index[:2].droplevel(["one", "two", "three"])

    with pytest.raises(KeyError, match="'Level four not found'"):
        index[:2].droplevel(["one", "four"])


def test_drop_not_lexsorted():
    # GH 12078

    # define the lexsorted version of the multi-index
    tuples = [("a", ""), ("b1", "c1"), ("b2", "c2")]
    lexsorted_mi = MultiIndex.from_tuples(tuples, names=["b", "c"])
    assert lexsorted_mi.is_lexsorted()

    # and the not-lexsorted version
    df = pd.DataFrame(
        columns=["a", "b", "c", "d"], data=[[1, "b1", "c1", 3], [1, "b2", "c2", 4]]
    )
    df = df.pivot_table(index="a", columns=["b", "c"], values="d")
    df = df.reset_index()
    not_lexsorted_mi = df.columns
    assert not not_lexsorted_mi.is_lexsorted()

    # compare the results
    tm.assert_index_equal(lexsorted_mi, not_lexsorted_mi)
    with tm.assert_produces_warning(PerformanceWarning):
        tm.assert_index_equal(lexsorted_mi.drop("a"), not_lexsorted_mi.drop("a"))


@pytest.mark.parametrize(
    "msg,labels,level",
    [
        (r"labels \[4\] not found in level", 4, "a"),
        (r"labels \[7\] not found in level", 7, "b"),
    ],
)
def test_drop_raise_exception_if_labels_not_in_level(msg, labels, level):
    # GH 8594
    mi = MultiIndex.from_arrays([[1, 2, 3], [4, 5, 6]], names=["a", "b"])
    s = pd.Series([10, 20, 30], index=mi)
    df = pd.DataFrame([10, 20, 30], index=mi)

    with pytest.raises(KeyError, match=msg):
        s.drop(labels, level=level)
    with pytest.raises(KeyError, match=msg):
        df.drop(labels, level=level)


@pytest.mark.parametrize("labels,level", [(4, "a"), (7, "b")])
def test_drop_errors_ignore(labels, level):
    # GH 8594
    mi = MultiIndex.from_arrays([[1, 2, 3], [4, 5, 6]], names=["a", "b"])
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
