# coding=utf-8

import pytest

import numpy as np

from pandas import Series, Categorical
import pandas.util.testing as tm


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
    assert not s.is_unique
    s = Series(np.arange(1000))
    assert s.is_unique


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


@pytest.mark.parametrize('keep, expected_inv_values', [
    ('first', [1, 4, 4, 16, 1]),
    ('last', [25, 9, 9, 16, 25])
])
def test_duplicated_inverse(keep, expected_inv_values):
    # GH 21357
    # check that return_inverse kwarg does not affect outcome;
    # index of inverse must be correctly transformed as well
    idx = [1, 4, 9, 16, 25]
    s = Series(['a', 'b', 'b', 'c', 'a'], index=idx)

    expected_isdup = s.duplicated(keep=keep)
    expected_inv = Series(expected_inv_values, index=idx)
    result_isdup, result_inv = s.duplicated(keep=keep,
                                            return_inverse=True)
    tm.assert_series_equal(result_isdup, expected_isdup)
    tm.assert_series_equal(result_inv, expected_inv)

    # test that result_inv works (and fits together with expected_isdup)
    unique = s.loc[~expected_isdup]
    reconstr = unique.reindex(result_inv)
    # Series has no set_index (GH21684)
    reconstr.index = result_inv.index
    tm.assert_series_equal(reconstr, s)


def test_duplicated_inverse_raises():
    s = Series(['a', 'b', 'b', 'c', 'a'])

    rgx = 'The parameters return_inverse=True and keep=False cannot be.*'
    with tm.assert_raises_regex(ValueError, rgx):
        s.duplicated(keep=False, return_inverse=True)


@pytest.mark.parametrize('keep', ['first', 'last'])
def test_duplicated_inverse_large(keep):
    # unsorted index important to check 'first'/'last' functionality
    s = Series(np.random.randint(0, 1000, 10000)).sample(5000)

    expected_isdup = s.duplicated(keep=keep)
    result_isdup, result_inv = s.duplicated(keep=keep, return_inverse=True)
    tm.assert_series_equal(result_isdup, expected_isdup)

    # test that result_inv works (and fits together with expected_isdup)
    unique = s.loc[~expected_isdup]
    reconstr = unique.reindex(result_inv)
    # Series has no set_index (GH21684)
    reconstr.index = result_inv.index
    tm.assert_series_equal(reconstr, s)


@pytest.mark.parametrize('keep', ['first', 'last'])
def test_duplicated_inverse_fastpath(keep):
    s = Series(range(10))  # no duplicates

    expected_isdup = s.duplicated(keep=keep)
    result_isdup, result_inv = s.duplicated(keep=keep,
                                            return_inverse=True)
    tm.assert_series_equal(result_isdup, expected_isdup)

    expected_inv = Series(range(10))
    tm.assert_series_equal(result_inv, expected_inv)
