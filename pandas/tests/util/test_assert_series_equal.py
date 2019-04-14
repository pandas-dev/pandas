# -*- coding: utf-8 -*-

import pytest

from pandas import Categorical, DataFrame, Series
from pandas.util.testing import assert_series_equal


def _assert_series_equal_both(a, b, **kwargs):
    """
    Check that two Series equal.

    This check is performed commutatively.

    Parameters
    ----------
    a : Series
        The first Series to compare.
    b : Series
        The second Series to compare.
    kwargs : dict
        The arguments passed to `assert_series_equal`.
    """
    assert_series_equal(a, b, **kwargs)
    assert_series_equal(b, a, **kwargs)


def _assert_not_series_equal(a, b, **kwargs):
    """
    Check that two Series are not equal.

    Parameters
    ----------
    a : Series
        The first Series to compare.
    b : Series
        The second Series to compare.
    kwargs : dict
        The arguments passed to `assert_series_equal`.
    """
    try:
        assert_series_equal(a, b, **kwargs)
        msg = "The two Series were equal when they shouldn't have been"

        pytest.fail(msg=msg)
    except AssertionError:
        pass


def _assert_not_series_equal_both(a, b, **kwargs):
    """
    Check that two Series are not equal.

    This check is performed commutatively.

    Parameters
    ----------
    a : Series
        The first Series to compare.
    b : Series
        The second Series to compare.
    kwargs : dict
        The arguments passed to `assert_series_equal`.
    """
    _assert_not_series_equal(a, b, **kwargs)
    _assert_not_series_equal(b, a, **kwargs)


@pytest.mark.parametrize("data", [
    range(3), list("abc"), list("áàä"),
])
def test_series_equal(data):
    _assert_series_equal_both(Series(data), Series(data))


@pytest.mark.parametrize("data1,data2", [
    (range(3), range(1, 4)),
    (list("abc"), list("xyz")),
    (list("áàä"), list("éèë")),
    (list("áàä"), list(b"aaa")),
    (range(3), range(4)),
])
def test_series_not_equal_value_mismatch(data1, data2):
    _assert_not_series_equal_both(Series(data1), Series(data2))


@pytest.mark.parametrize("kwargs", [
    dict(dtype="float64"),  # dtype mismatch
    dict(index=[1, 2, 4]),  # index mismatch
    dict(name="foo"),       # name mismatch
])
def test_series_not_equal_metadata_mismatch(kwargs):
    data = range(3)
    s1 = Series(data)

    s2 = Series(data, **kwargs)
    _assert_not_series_equal_both(s1, s2)


@pytest.mark.parametrize("data1,data2", [(0.12345, 0.12346), (0.1235, 0.1236)])
@pytest.mark.parametrize("dtype", ["float32", "float64"])
@pytest.mark.parametrize("check_less_precise", [False, True, 0, 1, 2, 3, 10])
def test_less_precise(data1, data2, dtype, check_less_precise):
    s1 = Series([data1], dtype=dtype)
    s2 = Series([data2], dtype=dtype)

    kwargs = dict(check_less_precise=check_less_precise)

    if ((check_less_precise is False or check_less_precise == 10) or
            ((check_less_precise is True or check_less_precise >= 3) and
             abs(data1 - data2) >= 0.0001)):
        msg = "Series values are different"
        with pytest.raises(AssertionError, match=msg):
            assert_series_equal(s1, s2, **kwargs)
    else:
        _assert_series_equal_both(s1, s2, **kwargs)


@pytest.mark.parametrize("s1,s2,msg", [
    # Index
    (Series(["l1", "l2"], index=[1, 2]),
     Series(["l1", "l2"], index=[1., 2.]),
     "Series\\.index are different"),

    # MultiIndex
    (DataFrame.from_records({"a": [1, 2], "b": [2.1, 1.5],
                             "c": ["l1", "l2"]}, index=["a", "b"]).c,
     DataFrame.from_records({"a": [1., 2.], "b": [2.1, 1.5],
                             "c": ["l1", "l2"]}, index=["a", "b"]).c,
     "MultiIndex level \\[0\\] are different")
])
def test_series_equal_index_dtype(s1, s2, msg, check_index_type):
    kwargs = dict(check_index_type=check_index_type)

    if check_index_type:
        with pytest.raises(AssertionError, match=msg):
            assert_series_equal(s1, s2, **kwargs)
    else:
        assert_series_equal(s1, s2, **kwargs)


def test_series_equal_length_mismatch(check_less_precise):
    msg = """Series are different

Series length are different
\\[left\\]:  3, RangeIndex\\(start=0, stop=3, step=1\\)
\\[right\\]: 4, RangeIndex\\(start=0, stop=4, step=1\\)"""

    s1 = Series([1, 2, 3])
    s2 = Series([1, 2, 3, 4])

    with pytest.raises(AssertionError, match=msg):
        assert_series_equal(s1, s2, check_less_precise=check_less_precise)


def test_series_equal_values_mismatch(check_less_precise):
    msg = """Series are different

Series values are different \\(33\\.33333 %\\)
\\[left\\]:  \\[1, 2, 3\\]
\\[right\\]: \\[1, 2, 4\\]"""

    s1 = Series([1, 2, 3])
    s2 = Series([1, 2, 4])

    with pytest.raises(AssertionError, match=msg):
        assert_series_equal(s1, s2, check_less_precise=check_less_precise)


def test_series_equal_categorical_mismatch(check_categorical):
    msg = """Attributes are different

Attribute "dtype" are different
\\[left\\]:  CategoricalDtype\\(categories=\\[u?'a', u?'b'\\], ordered=False\\)
\\[right\\]: CategoricalDtype\\(categories=\\[u?'a', u?'b', u?'c'\\], \
ordered=False\\)"""

    s1 = Series(Categorical(["a", "b"]))
    s2 = Series(Categorical(["a", "b"], categories=list("abc")))

    if check_categorical:
        with pytest.raises(AssertionError, match=msg):
            assert_series_equal(s1, s2, check_categorical=check_categorical)
    else:
        _assert_series_equal_both(s1, s2, check_categorical=check_categorical)
