import numpy as np
import pytest

from pandas.compat import lrange, product as cart_product

from pandas import DataFrame, Index, MultiIndex, Series, concat, date_range
import pandas.core.common as com
from pandas.util import testing as tm


@pytest.fixture
def four_level_index_dataframe():
    arr = np.array([[-0.5109, -2.3358, -0.4645, 0.05076, 0.364],
                    [0.4473, 1.4152, 0.2834, 1.00661, 0.1744],
                    [-0.6662, -0.5243, -0.358, 0.89145, 2.5838]])
    index = MultiIndex(
        levels=[['a', 'x'], ['b', 'q'], [10.0032, 20.0, 30.0], [3, 4, 5]],
        codes=[[0, 0, 1], [0, 1, 1], [0, 1, 2], [2, 1, 0]],
        names=['one', 'two', 'three', 'four'])
    return DataFrame(arr, index=index, columns=list('ABCDE'))


@pytest.mark.parametrize('key, level, exp_arr, exp_index', [
    ('a', 'lvl0', lambda x: x[:, 0:2], Index(['bar', 'foo'], name='lvl1')),
    ('foo', 'lvl1', lambda x: x[:, 1:2], Index(['a'], name='lvl0'))
])
def test_xs_named_levels_axis_eq_1(key, level, exp_arr, exp_index):
    # see gh-2903
    arr = np.random.randn(4, 4)
    index = MultiIndex(levels=[['a', 'b'], ['bar', 'foo', 'hello', 'world']],
                       codes=[[0, 0, 1, 1], [0, 1, 2, 3]],
                       names=['lvl0', 'lvl1'])
    df = DataFrame(arr, columns=index)
    result = df.xs(key, level=level, axis=1)
    expected = DataFrame(exp_arr(arr), columns=exp_index)
    tm.assert_frame_equal(result, expected)


def test_xs_values(multiindex_dataframe_random_data):
    df = multiindex_dataframe_random_data
    result = df.xs(('bar', 'two')).values
    expected = df.values[4]
    tm.assert_almost_equal(result, expected)


def test_xs_loc_equality(multiindex_dataframe_random_data):
    df = multiindex_dataframe_random_data
    result = df.xs(('bar', 'two'))
    expected = df.loc[('bar', 'two')]
    tm.assert_series_equal(result, expected)


def test_xs_missing_values_in_index():
    # see gh-6574
    # missing values in returned index should be preserrved
    acc = [
        ('a', 'abcde', 1),
        ('b', 'bbcde', 2),
        ('y', 'yzcde', 25),
        ('z', 'xbcde', 24),
        ('z', None, 26),
        ('z', 'zbcde', 25),
        ('z', 'ybcde', 26),
    ]
    df = DataFrame(acc,
                   columns=['a1', 'a2', 'cnt']).set_index(['a1', 'a2'])
    expected = DataFrame({'cnt': [24, 26, 25, 26]}, index=Index(
        ['xbcde', np.nan, 'zbcde', 'ybcde'], name='a2'))

    result = df.xs('z', level='a1')
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize('key, level', [
    ('one', 'second'),
    (['one'], ['second'])
])
def test_xs_with_duplicates(key, level, multiindex_dataframe_random_data):
    # see gh-13719
    frame = multiindex_dataframe_random_data
    df = concat([frame] * 2)
    assert df.index.is_unique is False
    expected = concat([frame.xs('one', level='second')] * 2)

    result = df.xs(key, level=level)
    tm.assert_frame_equal(result, expected)


def test_xs_level(multiindex_dataframe_random_data):
    df = multiindex_dataframe_random_data
    result = df.xs('two', level='second')
    expected = df[df.index.get_level_values(1) == 'two']
    expected.index = Index(['foo', 'bar', 'baz', 'qux'], name='first')
    tm.assert_frame_equal(result, expected)


def test_xs_level_eq_2():
    arr = np.random.randn(3, 5)
    index = MultiIndex(
        levels=[['a', 'p', 'x'], ['b', 'q', 'y'], ['c', 'r', 'z']],
        codes=[[2, 0, 1], [2, 0, 1], [2, 0, 1]])
    df = DataFrame(arr, index=index)
    expected = DataFrame(arr[1:2], index=[['a'], ['b']])
    result = df.xs('c', level=2)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize('indexer', [
    lambda df: df.xs(('a', 4), level=['one', 'four']),
    lambda df: df.xs('a').xs(4, level='four')
])
def test_xs_level_multiple(indexer, four_level_index_dataframe):
    df = four_level_index_dataframe
    expected_values = [[0.4473, 1.4152, 0.2834, 1.00661, 0.1744]]
    expected_index = MultiIndex(
        levels=[['q'], [20.0]],
        codes=[[0], [0]],
        names=['two', 'three'])
    expected = DataFrame(
        expected_values, index=expected_index, columns=list('ABCDE'))
    result = indexer(df)
    tm.assert_frame_equal(result, expected)


def test_xs_setting_with_copy_error(multiindex_dataframe_random_data):
    # this is a copy in 0.14
    df = multiindex_dataframe_random_data
    result = df.xs('two', level='second')

    # setting this will give a SettingWithCopyError
    # as we are trying to write a view
    msg = 'A value is trying to be set on a copy of a slice from a DataFrame'
    with pytest.raises(com.SettingWithCopyError, match=msg):
        result[:] = 10


def test_xs_setting_with_copy_error_multiple(four_level_index_dataframe):
    # this is a copy in 0.14
    df = four_level_index_dataframe
    result = df.xs(('a', 4), level=['one', 'four'])

    # setting this will give a SettingWithCopyError
    # as we are trying to write a view
    msg = 'A value is trying to be set on a copy of a slice from a DataFrame'
    with pytest.raises(com.SettingWithCopyError, match=msg):
        result[:] = 10


def test_xs_integer_key():
    # see gh-2107
    dates = lrange(20111201, 20111205)
    ids = 'abcde'
    index = MultiIndex.from_tuples(
        [x for x in cart_product(dates, ids)],
        names=['date', 'secid'])
    df = DataFrame(
        np.random.randn(len(index), 3), index, ['X', 'Y', 'Z'])

    result = df.xs(20111201, level='date')
    expected = df.loc[20111201, :]
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize('indexer', [
    lambda df: df.xs('a', level=0),
    lambda df: df.xs('a')
])
def test_xs_level0(indexer, four_level_index_dataframe):
    df = four_level_index_dataframe
    expected_values = [[-0.5109, -2.3358, -0.4645, 0.05076, 0.364],
                       [0.4473, 1.4152, 0.2834, 1.00661, 0.1744]]
    expected_index = MultiIndex(
        levels=[['b', 'q'], [10.0032, 20.0], [4, 5]],
        codes=[[0, 1], [0, 1], [1, 0]],
        names=['two', 'three', 'four'])
    expected = DataFrame(
        expected_values, index=expected_index, columns=list('ABCDE'))

    result = indexer(df)
    tm.assert_frame_equal(result, expected)


def test_xs_level_series(multiindex_dataframe_random_data):
    # this test is not explicitly testing .xs functionality
    # TODO: move to another module or refactor
    df = multiindex_dataframe_random_data
    s = df['A']
    result = s[:, 'two']
    expected = df.xs('two', level=1)['A']
    tm.assert_series_equal(result, expected)


def test_xs_level_series_ymd(multiindex_year_month_day_dataframe_random_data):
    # this test is not explicitly testing .xs functionality
    # TODO: move to another module or refactor
    df = multiindex_year_month_day_dataframe_random_data
    s = df['A']
    result = s[2000, 5]
    expected = df.loc[2000, 5]['A']
    tm.assert_series_equal(result, expected)


def test_xs_level_series_slice_not_implemented(
        multiindex_year_month_day_dataframe_random_data):
    # this test is not explicitly testing .xs functionality
    # TODO: move to another module or refactor
    # not implementing this for now
    df = multiindex_year_month_day_dataframe_random_data
    s = df['A']

    msg = r'\(2000, slice\(3, 4, None\)\)'
    with pytest.raises(TypeError, match=msg):
        s[2000, 3:4]


def test_series_getitem_multiindex_xs():
    # GH6258
    dt = list(date_range('20130903', periods=3))
    idx = MultiIndex.from_product([list('AB'), dt])
    s = Series([1, 3, 4, 1, 3, 4], index=idx)
    expected = Series([1, 1], index=list('AB'))

    result = s.xs('20130903', level=1)
    tm.assert_series_equal(result, expected)


def test_series_getitem_multiindex_xs_by_label():
    # GH5684
    idx = MultiIndex.from_tuples([('a', 'one'), ('a', 'two'), ('b', 'one'),
                                  ('b', 'two')])
    s = Series([1, 2, 3, 4], index=idx)
    s.index.set_names(['L1', 'L2'], inplace=True)
    expected = Series([1, 3], index=['a', 'b'])
    expected.index.set_names(['L1'], inplace=True)

    result = s.xs('one', level='L2')
    tm.assert_series_equal(result, expected)
