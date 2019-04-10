# -*- coding: utf-8 -*-

from __future__ import print_function

import numpy as np
import pytest

from pandas.compat import lrange, string_types

from pandas import DataFrame, Series
import pandas.util.testing as tm


@pytest.mark.parametrize('subset', ['a', ['a'], ['a', 'B']])
def test_duplicated_with_misspelled_column_name(subset):
    # GH 19730
    df = DataFrame({'A': [0, 0, 1],
                    'B': [0, 0, 1],
                    'C': [0, 0, 1]})

    with pytest.raises(KeyError):
        df.duplicated(subset)

    with pytest.raises(KeyError):
        df.drop_duplicates(subset)


@pytest.mark.slow
def test_duplicated_do_not_fail_on_wide_dataframes():
    # gh-21524
    # Given the wide dataframe with a lot of columns
    # with different (important!) values
    data = {'col_{0:02d}'.format(i): np.random.randint(0, 1000, 30000)
            for i in range(100)}
    df = DataFrame(data).T
    result = df.duplicated()

    # Then duplicates produce the bool Series as a result and don't fail during
    # calculation. Actual values doesn't matter here, though usually it's all
    # False in this case
    assert isinstance(result, Series)
    assert result.dtype == np.bool


@pytest.mark.parametrize('keep, expected', [
    ('first', Series([False, False, True, False, True])),
    ('last', Series([True, True, False, False, False])),
    (False, Series([True, True, True, False, True]))
])
def test_duplicated_keep(keep, expected):
    df = DataFrame({'A': [0, 1, 1, 2, 0], 'B': ['a', 'b', 'b', 'c', 'a']})

    result = df.duplicated(keep=keep)
    tm.assert_series_equal(result, expected)


@pytest.mark.xfail(reason="GH#21720; nan/None falsely considered equal")
@pytest.mark.parametrize('keep, expected', [
    ('first', Series([False, False, True, False, True])),
    ('last', Series([True, True, False, False, False])),
    (False, Series([True, True, True, False, True]))
])
def test_duplicated_nan_none(keep, expected):
    df = DataFrame({'C': [np.nan, 3, 3, None, np.nan]}, dtype=object)

    result = df.duplicated(keep=keep)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize('keep', ['first', 'last', False])
@pytest.mark.parametrize('subset', [None, ['A', 'B'], 'A'])
def test_duplicated_subset(subset, keep):
    df = DataFrame({'A': [0, 1, 1, 2, 0],
                    'B': ['a', 'b', 'b', 'c', 'a'],
                    'C': [np.nan, 3, 3, None, np.nan]})

    if subset is None:
        subset = list(df.columns)
    elif isinstance(subset, string_types):
        # need to have a DataFrame, not a Series
        # -> select columns with singleton list, not string
        subset = [subset]

    expected = df[subset].duplicated(keep=keep)
    result = df.duplicated(keep=keep, subset=subset)
    tm.assert_series_equal(result, expected)


def test_drop_duplicates():
    df = DataFrame({'AAA': ['foo', 'bar', 'foo', 'bar',
                            'foo', 'bar', 'bar', 'foo'],
                    'B': ['one', 'one', 'two', 'two',
                          'two', 'two', 'one', 'two'],
                    'C': [1, 1, 2, 2, 2, 2, 1, 2],
                    'D': lrange(8)})

    # single column
    result = df.drop_duplicates('AAA')
    expected = df[:2]
    tm.assert_frame_equal(result, expected)

    result = df.drop_duplicates('AAA', keep='last')
    expected = df.loc[[6, 7]]
    tm.assert_frame_equal(result, expected)

    result = df.drop_duplicates('AAA', keep=False)
    expected = df.loc[[]]
    tm.assert_frame_equal(result, expected)
    assert len(result) == 0

    # multi column
    expected = df.loc[[0, 1, 2, 3]]
    result = df.drop_duplicates(np.array(['AAA', 'B']))
    tm.assert_frame_equal(result, expected)
    result = df.drop_duplicates(['AAA', 'B'])
    tm.assert_frame_equal(result, expected)

    result = df.drop_duplicates(('AAA', 'B'), keep='last')
    expected = df.loc[[0, 5, 6, 7]]
    tm.assert_frame_equal(result, expected)

    result = df.drop_duplicates(('AAA', 'B'), keep=False)
    expected = df.loc[[0]]
    tm.assert_frame_equal(result, expected)

    # consider everything
    df2 = df.loc[:, ['AAA', 'B', 'C']]

    result = df2.drop_duplicates()
    # in this case only
    expected = df2.drop_duplicates(['AAA', 'B'])
    tm.assert_frame_equal(result, expected)

    result = df2.drop_duplicates(keep='last')
    expected = df2.drop_duplicates(['AAA', 'B'], keep='last')
    tm.assert_frame_equal(result, expected)

    result = df2.drop_duplicates(keep=False)
    expected = df2.drop_duplicates(['AAA', 'B'], keep=False)
    tm.assert_frame_equal(result, expected)

    # integers
    result = df.drop_duplicates('C')
    expected = df.iloc[[0, 2]]
    tm.assert_frame_equal(result, expected)
    result = df.drop_duplicates('C', keep='last')
    expected = df.iloc[[-2, -1]]
    tm.assert_frame_equal(result, expected)

    df['E'] = df['C'].astype('int8')
    result = df.drop_duplicates('E')
    expected = df.iloc[[0, 2]]
    tm.assert_frame_equal(result, expected)
    result = df.drop_duplicates('E', keep='last')
    expected = df.iloc[[-2, -1]]
    tm.assert_frame_equal(result, expected)

    # GH 11376
    df = DataFrame({'x': [7, 6, 3, 3, 4, 8, 0],
                    'y': [0, 6, 5, 5, 9, 1, 2]})
    expected = df.loc[df.index != 3]
    tm.assert_frame_equal(df.drop_duplicates(), expected)

    df = DataFrame([[1, 0], [0, 2]])
    tm.assert_frame_equal(df.drop_duplicates(), df)

    df = DataFrame([[-2, 0], [0, -4]])
    tm.assert_frame_equal(df.drop_duplicates(), df)

    x = np.iinfo(np.int64).max / 3 * 2
    df = DataFrame([[-x, x], [0, x + 4]])
    tm.assert_frame_equal(df.drop_duplicates(), df)

    df = DataFrame([[-x, x], [x, x + 4]])
    tm.assert_frame_equal(df.drop_duplicates(), df)

    # GH 11864
    df = DataFrame([i] * 9 for i in range(16))
    df = df.append([[1] + [0] * 8], ignore_index=True)

    for keep in ['first', 'last', False]:
        assert df.duplicated(keep=keep).sum() == 0


def test_duplicated_on_empty_frame():
    # GH 25184

    df = DataFrame(columns=['a', 'b'])
    dupes = df.duplicated('a')

    result = df[dupes]
    expected = df.copy()
    tm.assert_frame_equal(result, expected)


def test_drop_duplicates_with_duplicate_column_names():
    # GH17836
    df = DataFrame([
        [1, 2, 5],
        [3, 4, 6],
        [3, 4, 7]
    ], columns=['a', 'a', 'b'])

    result0 = df.drop_duplicates()
    tm.assert_frame_equal(result0, df)

    result1 = df.drop_duplicates('a')
    expected1 = df[:2]
    tm.assert_frame_equal(result1, expected1)


def test_drop_duplicates_for_take_all():
    df = DataFrame({'AAA': ['foo', 'bar', 'baz', 'bar',
                            'foo', 'bar', 'qux', 'foo'],
                    'B': ['one', 'one', 'two', 'two',
                          'two', 'two', 'one', 'two'],
                    'C': [1, 1, 2, 2, 2, 2, 1, 2],
                    'D': lrange(8)})

    # single column
    result = df.drop_duplicates('AAA')
    expected = df.iloc[[0, 1, 2, 6]]
    tm.assert_frame_equal(result, expected)

    result = df.drop_duplicates('AAA', keep='last')
    expected = df.iloc[[2, 5, 6, 7]]
    tm.assert_frame_equal(result, expected)

    result = df.drop_duplicates('AAA', keep=False)
    expected = df.iloc[[2, 6]]
    tm.assert_frame_equal(result, expected)

    # multiple columns
    result = df.drop_duplicates(['AAA', 'B'])
    expected = df.iloc[[0, 1, 2, 3, 4, 6]]
    tm.assert_frame_equal(result, expected)

    result = df.drop_duplicates(['AAA', 'B'], keep='last')
    expected = df.iloc[[0, 1, 2, 5, 6, 7]]
    tm.assert_frame_equal(result, expected)

    result = df.drop_duplicates(['AAA', 'B'], keep=False)
    expected = df.iloc[[0, 1, 2, 6]]
    tm.assert_frame_equal(result, expected)


def test_drop_duplicates_tuple():
    df = DataFrame({('AA', 'AB'): ['foo', 'bar', 'foo', 'bar',
                                   'foo', 'bar', 'bar', 'foo'],
                    'B': ['one', 'one', 'two', 'two',
                          'two', 'two', 'one', 'two'],
                    'C': [1, 1, 2, 2, 2, 2, 1, 2],
                    'D': lrange(8)})

    # single column
    result = df.drop_duplicates(('AA', 'AB'))
    expected = df[:2]
    tm.assert_frame_equal(result, expected)

    result = df.drop_duplicates(('AA', 'AB'), keep='last')
    expected = df.loc[[6, 7]]
    tm.assert_frame_equal(result, expected)

    result = df.drop_duplicates(('AA', 'AB'), keep=False)
    expected = df.loc[[]]  # empty df
    assert len(result) == 0
    tm.assert_frame_equal(result, expected)

    # multi column
    expected = df.loc[[0, 1, 2, 3]]
    result = df.drop_duplicates((('AA', 'AB'), 'B'))
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize('df', [
    DataFrame(),
    DataFrame(columns=[]),
    DataFrame(columns=['A', 'B', 'C']),
    DataFrame(index=[]),
    DataFrame(index=['A', 'B', 'C'])
])
def test_drop_duplicates_empty(df):
    # GH 20516
    result = df.drop_duplicates()
    tm.assert_frame_equal(result, df)

    result = df.copy()
    result.drop_duplicates(inplace=True)
    tm.assert_frame_equal(result, df)


def test_drop_duplicates_NA():
    # none
    df = DataFrame({'A': [None, None, 'foo', 'bar',
                          'foo', 'bar', 'bar', 'foo'],
                    'B': ['one', 'one', 'two', 'two',
                          'two', 'two', 'one', 'two'],
                    'C': [1.0, np.nan, np.nan, np.nan, 1., 1., 1, 1.],
                    'D': lrange(8)})

    # single column
    result = df.drop_duplicates('A')
    expected = df.loc[[0, 2, 3]]
    tm.assert_frame_equal(result, expected)

    result = df.drop_duplicates('A', keep='last')
    expected = df.loc[[1, 6, 7]]
    tm.assert_frame_equal(result, expected)

    result = df.drop_duplicates('A', keep=False)
    expected = df.loc[[]]  # empty df
    tm.assert_frame_equal(result, expected)
    assert len(result) == 0

    # multi column
    result = df.drop_duplicates(['A', 'B'])
    expected = df.loc[[0, 2, 3, 6]]
    tm.assert_frame_equal(result, expected)

    result = df.drop_duplicates(['A', 'B'], keep='last')
    expected = df.loc[[1, 5, 6, 7]]
    tm.assert_frame_equal(result, expected)

    result = df.drop_duplicates(['A', 'B'], keep=False)
    expected = df.loc[[6]]
    tm.assert_frame_equal(result, expected)

    # nan
    df = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                          'foo', 'bar', 'bar', 'foo'],
                    'B': ['one', 'one', 'two', 'two',
                          'two', 'two', 'one', 'two'],
                    'C': [1.0, np.nan, np.nan, np.nan, 1., 1., 1, 1.],
                    'D': lrange(8)})

    # single column
    result = df.drop_duplicates('C')
    expected = df[:2]
    tm.assert_frame_equal(result, expected)

    result = df.drop_duplicates('C', keep='last')
    expected = df.loc[[3, 7]]
    tm.assert_frame_equal(result, expected)

    result = df.drop_duplicates('C', keep=False)
    expected = df.loc[[]]  # empty df
    tm.assert_frame_equal(result, expected)
    assert len(result) == 0

    # multi column
    result = df.drop_duplicates(['C', 'B'])
    expected = df.loc[[0, 1, 2, 4]]
    tm.assert_frame_equal(result, expected)

    result = df.drop_duplicates(['C', 'B'], keep='last')
    expected = df.loc[[1, 3, 6, 7]]
    tm.assert_frame_equal(result, expected)

    result = df.drop_duplicates(['C', 'B'], keep=False)
    expected = df.loc[[1]]
    tm.assert_frame_equal(result, expected)


def test_drop_duplicates_NA_for_take_all():
    # none
    df = DataFrame({'A': [None, None, 'foo', 'bar',
                          'foo', 'baz', 'bar', 'qux'],
                    'C': [1.0, np.nan, np.nan, np.nan, 1., 2., 3, 1.]})

    # single column
    result = df.drop_duplicates('A')
    expected = df.iloc[[0, 2, 3, 5, 7]]
    tm.assert_frame_equal(result, expected)

    result = df.drop_duplicates('A', keep='last')
    expected = df.iloc[[1, 4, 5, 6, 7]]
    tm.assert_frame_equal(result, expected)

    result = df.drop_duplicates('A', keep=False)
    expected = df.iloc[[5, 7]]
    tm.assert_frame_equal(result, expected)

    # nan

    # single column
    result = df.drop_duplicates('C')
    expected = df.iloc[[0, 1, 5, 6]]
    tm.assert_frame_equal(result, expected)

    result = df.drop_duplicates('C', keep='last')
    expected = df.iloc[[3, 5, 6, 7]]
    tm.assert_frame_equal(result, expected)

    result = df.drop_duplicates('C', keep=False)
    expected = df.iloc[[5, 6]]
    tm.assert_frame_equal(result, expected)


def test_drop_duplicates_inplace():
    orig = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                            'foo', 'bar', 'bar', 'foo'],
                      'B': ['one', 'one', 'two', 'two',
                            'two', 'two', 'one', 'two'],
                      'C': [1, 1, 2, 2, 2, 2, 1, 2],
                      'D': lrange(8)})

    # single column
    df = orig.copy()
    df.drop_duplicates('A', inplace=True)
    expected = orig[:2]
    result = df
    tm.assert_frame_equal(result, expected)

    df = orig.copy()
    df.drop_duplicates('A', keep='last', inplace=True)
    expected = orig.loc[[6, 7]]
    result = df
    tm.assert_frame_equal(result, expected)

    df = orig.copy()
    df.drop_duplicates('A', keep=False, inplace=True)
    expected = orig.loc[[]]
    result = df
    tm.assert_frame_equal(result, expected)
    assert len(df) == 0

    # multi column
    df = orig.copy()
    df.drop_duplicates(['A', 'B'], inplace=True)
    expected = orig.loc[[0, 1, 2, 3]]
    result = df
    tm.assert_frame_equal(result, expected)

    df = orig.copy()
    df.drop_duplicates(['A', 'B'], keep='last', inplace=True)
    expected = orig.loc[[0, 5, 6, 7]]
    result = df
    tm.assert_frame_equal(result, expected)

    df = orig.copy()
    df.drop_duplicates(['A', 'B'], keep=False, inplace=True)
    expected = orig.loc[[0]]
    result = df
    tm.assert_frame_equal(result, expected)

    # consider everything
    orig2 = orig.loc[:, ['A', 'B', 'C']].copy()

    df2 = orig2.copy()
    df2.drop_duplicates(inplace=True)
    # in this case only
    expected = orig2.drop_duplicates(['A', 'B'])
    result = df2
    tm.assert_frame_equal(result, expected)

    df2 = orig2.copy()
    df2.drop_duplicates(keep='last', inplace=True)
    expected = orig2.drop_duplicates(['A', 'B'], keep='last')
    result = df2
    tm.assert_frame_equal(result, expected)

    df2 = orig2.copy()
    df2.drop_duplicates(keep=False, inplace=True)
    expected = orig2.drop_duplicates(['A', 'B'], keep=False)
    result = df2
    tm.assert_frame_equal(result, expected)
