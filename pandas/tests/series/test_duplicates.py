# coding=utf-8

import numpy as np
import pytest

import pandas.util.testing as tm
from pandas import Categorical, Series


def test_value_counts_nunique():
    # basics.rst doc example
    series = Series(np.random.randn(500))
    series[20:500] = np.nan
    series[10:20] = 5000
    result = series.nunique()
    assert result == 11

    # GH 18051
    s = Series(Categorical([]))
    assert s.nunique() == 0
    s = Series(Categorical([np.nan]))
    assert s.nunique() == 0


def test_unique():
    # GH714 also, dtype=float
    s = Series([1.2345] * 100)
    s[::2] = np.nan
    result = s.unique()
    assert len(result) == 2

    s = Series([1.2345] * 100, dtype='f4')
    s[::2] = np.nan
    result = s.unique()
    assert len(result) == 2

    # NAs in object arrays #714
    s = Series(['foo'] * 100, dtype='O')
    s[::2] = np.nan
    result = s.unique()
    assert len(result) == 2

    # decision about None
    s = Series([1, 2, 3, None, None, None], dtype=object)
    result = s.unique()
    expected = np.array([1, 2, 3, None], dtype=object)
    tm.assert_numpy_array_equal(result, expected)

    # GH 18051
    s = Series(Categorical([]))
    tm.assert_categorical_equal(s.unique(), Categorical([]), check_dtype=False)
    s = Series(Categorical([np.nan]))
    tm.assert_categorical_equal(s.unique(), Categorical([np.nan]),
                                check_dtype=False)


def test_unique_data_ownership():
    # it works! #1807
    Series(Series(["a", "c", "b"]).unique()).sort_values()


def test_is_unique():
    # GH11946
    s = Series(np.random.randint(0, 10, size=1000))
    assert s.is_unique is False
    s = Series(np.arange(1000))
    assert s.is_unique is True


def test_is_unique_class_ne(capsys):
    # GH 20661
    class Foo(object):
        def __init__(self, val):
            self._value = val

        def __ne__(self, other):
            raise Exception("NEQ not supported")

    li = [Foo(i) for i in range(5)]
    s = Series(li, index=[i for i in range(5)])
    _, err = capsys.readouterr()
    s.is_unique
    _, err = capsys.readouterr()
    assert len(err) == 0


@pytest.mark.parametrize(
    'keep, expected',
    [
        ('first', Series([False, False, False, False, True, True, False])),
        ('last', Series([False, True, True, False, False, False, False])),
        (False, Series([False, True, True, False, True, True, False]))
    ])
def test_drop_duplicates_non_bool(any_numpy_dtype, keep, expected):
    tc = Series([1, 2, 3, 5, 3, 2, 4], dtype=np.dtype(any_numpy_dtype))

    tm.assert_series_equal(tc.duplicated(keep=keep), expected)
    tm.assert_series_equal(tc.drop_duplicates(keep=keep), tc[~expected])
    sc = tc.copy()
    sc.drop_duplicates(keep=keep, inplace=True)
    tm.assert_series_equal(sc, tc[~expected])


@pytest.mark.parametrize('keep, expected',
                         [('first', Series([False, False, True, True])),
                          ('last', Series([True, True, False, False])),
                          (False, Series([True, True, True, True]))])
def test_drop_duplicates_bool(keep, expected):
    tc = Series([True, False, True, False])

    tm.assert_series_equal(tc.duplicated(keep=keep), expected)
    tm.assert_series_equal(tc.drop_duplicates(keep=keep), tc[~expected])
    sc = tc.copy()
    sc.drop_duplicates(keep=keep, inplace=True)
    tm.assert_series_equal(sc, tc[~expected])


@pytest.mark.parametrize('keep, expected', [
    ('first', Series([False, False, True, False, True], name='name')),
    ('last', Series([True, True, False, False, False], name='name')),
    (False, Series([True, True, True, False, True], name='name'))
])
def test_duplicated_keep(keep, expected):
    s = Series(['a', 'b', 'b', 'c', 'a'], name='name')

    result = s.duplicated(keep=keep)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize('keep, expected', [
    ('first', Series([False, False, True, False, True])),
    ('last', Series([True, True, False, False, False])),
    (False, Series([True, True, True, False, True]))
])
def test_duplicated_nan_none(keep, expected):
    s = Series([np.nan, 3, 3, None, np.nan], dtype=object)

    result = s.duplicated(keep=keep)
    tm.assert_series_equal(result, expected)
