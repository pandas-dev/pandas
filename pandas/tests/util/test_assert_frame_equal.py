# -*- coding: utf-8 -*-

import pytest

from pandas import DataFrame
from pandas.util.testing import assert_frame_equal


@pytest.fixture(params=[True, False])
def by_blocks(request):
    return request.param


def _assert_frame_equal_both(a, b, **kwargs):
    """
    Check that two DataFrame equal.

    This check is performed commutatively.

    Parameters
    ----------
    a : DataFrame
        The first DataFrame to compare.
    b : DataFrame
        The second DataFrame to compare.
    kwargs : dict
        The arguments passed to `assert_frame_equal`.
    """
    assert_frame_equal(a, b, **kwargs)
    assert_frame_equal(b, a, **kwargs)


def _assert_not_frame_equal(a, b, **kwargs):
    """
    Check that two DataFrame are not equal.

    Parameters
    ----------
    a : DataFrame
        The first DataFrame to compare.
    b : DataFrame
        The second DataFrame to compare.
    kwargs : dict
        The arguments passed to `assert_frame_equal`.
    """
    try:
        assert_frame_equal(a, b, **kwargs)
        msg = "The two DataFrames were equal when they shouldn't have been"

        pytest.fail(msg=msg)
    except AssertionError:
        pass


def _assert_not_frame_equal_both(a, b, **kwargs):
    """
    Check that two DataFrame are not equal.

    This check is performed commutatively.

    Parameters
    ----------
    a : DataFrame
        The first DataFrame to compare.
    b : DataFrame
        The second DataFrame to compare.
    kwargs : dict
        The arguments passed to `assert_frame_equal`.
    """
    _assert_not_frame_equal(a, b, **kwargs)
    _assert_not_frame_equal(b, a, **kwargs)


@pytest.mark.parametrize("check_like", [True, False])
def test_frame_equal_row_order_mismatch(check_like):
    df1 = DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]},
                    index=["a", "b", "c"])
    df2 = DataFrame({"A": [3, 2, 1], "B": [6, 5, 4]},
                    index=["c", "b", "a"])

    if not check_like:  # Do not ignore row-column orderings.
        msg = "DataFrame.index are different"
        with pytest.raises(AssertionError, match=msg):
            assert_frame_equal(df1, df2, check_like=check_like)
    else:
        _assert_frame_equal_both(df1, df2, check_like=check_like)


@pytest.mark.parametrize("df1,df2", [
    (DataFrame({"A": [1, 2, 3]}), DataFrame({"A": [1, 2, 3, 4]})),
    (DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}), DataFrame({"A": [1, 2, 3]})),
])
def test_frame_equal_shape_mismatch(df1, df2):
    msg = "DataFrame are different"

    with pytest.raises(AssertionError, match=msg):
        assert_frame_equal(df1, df2)


@pytest.mark.parametrize("df1,df2,msg", [
    # Index
    (DataFrame.from_records({"a": [1, 2],
                             "c": ["l1", "l2"]}, index=["a"]),
     DataFrame.from_records({"a": [1.0, 2.0],
                             "c": ["l1", "l2"]}, index=["a"]),
     "DataFrame\\.index are different"),

    # MultiIndex
    (DataFrame.from_records({"a": [1, 2], "b": [2.1, 1.5],
                             "c": ["l1", "l2"]}, index=["a", "b"]),
     DataFrame.from_records({"a": [1.0, 2.0], "b": [2.1, 1.5],
                             "c": ["l1", "l2"]}, index=["a", "b"]),
     "MultiIndex level \\[0\\] are different")
])
def test_frame_equal_index_dtype_mismatch(df1, df2, msg, check_index_type):
    kwargs = dict(check_index_type=check_index_type)

    if check_index_type:
        with pytest.raises(AssertionError, match=msg):
            assert_frame_equal(df1, df2, **kwargs)
    else:
        assert_frame_equal(df1, df2, **kwargs)


def test_empty_dtypes(check_dtype):
    columns = ["col1", "col2"]
    df1 = DataFrame(columns=columns)
    df2 = DataFrame(columns=columns)

    kwargs = dict(check_dtype=check_dtype)
    df1["col1"] = df1["col1"].astype("int64")

    if check_dtype:
        msg = "Attributes are different"
        with pytest.raises(AssertionError, match=msg):
            assert_frame_equal(df1, df2, **kwargs)
    else:
        assert_frame_equal(df1, df2, **kwargs)


def test_frame_equal_index_mismatch():
    msg = """DataFrame\\.index are different

DataFrame\\.index values are different \\(33\\.33333 %\\)
\\[left\\]:  Index\\(\\[u?'a', u?'b', u?'c'\\], dtype='object'\\)
\\[right\\]: Index\\(\\[u?'a', u?'b', u?'d'\\], dtype='object'\\)"""

    df1 = DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]},
                    index=["a", "b", "c"])
    df2 = DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]},
                    index=["a", "b", "d"])

    with pytest.raises(AssertionError, match=msg):
        assert_frame_equal(df1, df2)


def test_frame_equal_columns_mismatch():
    msg = """DataFrame\\.columns are different

DataFrame\\.columns values are different \\(50\\.0 %\\)
\\[left\\]:  Index\\(\\[u?'A', u?'B'\\], dtype='object'\\)
\\[right\\]: Index\\(\\[u?'A', u?'b'\\], dtype='object'\\)"""

    df1 = DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]},
                    index=["a", "b", "c"])
    df2 = DataFrame({"A": [1, 2, 3], "b": [4, 5, 6]},
                    index=["a", "b", "c"])

    with pytest.raises(AssertionError, match=msg):
        assert_frame_equal(df1, df2)


def test_frame_equal_block_mismatch(by_blocks):
    msg = """DataFrame\\.iloc\\[:, 1\\] are different

DataFrame\\.iloc\\[:, 1\\] values are different \\(33\\.33333 %\\)
\\[left\\]:  \\[4, 5, 6\\]
\\[right\\]: \\[4, 5, 7\\]"""

    df1 = DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
    df2 = DataFrame({"A": [1, 2, 3], "B": [4, 5, 7]})

    with pytest.raises(AssertionError, match=msg):
        assert_frame_equal(df1, df2, by_blocks=by_blocks)


@pytest.mark.parametrize("df1,df2,msg", [
    (DataFrame({"A": ["á", "à", "ä"], "E": ["é", "è", "ë"]}),
     DataFrame({"A": ["á", "à", "ä"], "E": ["é", "è", "e̊"]}),
     """DataFrame\\.iloc\\[:, 1\\] are different

DataFrame\\.iloc\\[:, 1\\] values are different \\(33\\.33333 %\\)
\\[left\\]:  \\[é, è, ë\\]
\\[right\\]: \\[é, è, e̊\\]"""),
    (DataFrame({"A": ["á", "à", "ä"], "E": ["é", "è", "ë"]}),
     DataFrame({"A": ["a", "a", "a"], "E": ["e", "e", "e"]}),
     """DataFrame\\.iloc\\[:, 0\\] are different

DataFrame\\.iloc\\[:, 0\\] values are different \\(100\\.0 %\\)
\\[left\\]:  \\[á, à, ä\\]
\\[right\\]: \\[a, a, a\\]"""),
])
def test_frame_equal_unicode(df1, df2, msg, by_blocks):
    # see gh-20503
    #
    # Test ensures that `assert_frame_equals` raises the right exception
    # when comparing DataFrames containing differing unicode objects.
    with pytest.raises(AssertionError, match=msg):
        assert_frame_equal(df1, df2, by_blocks=by_blocks)
