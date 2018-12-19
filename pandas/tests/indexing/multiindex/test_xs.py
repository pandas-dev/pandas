import numpy as np
import pytest

from pandas.compat import StringIO, lrange, product as cart_product

from pandas import DataFrame, Index, MultiIndex, concat, read_csv
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


def test_xs_level0():
    text = """                      A       B       C       D        E
one two three   four
a   b   10.0032 5    -0.5109 -2.3358 -0.4645  0.05076  0.3640
a   q   20      4     0.4473  1.4152  0.2834  1.00661  0.1744
x   q   30      3    -0.6662 -0.5243 -0.3580  0.89145  2.5838"""

    df = read_csv(StringIO(text), sep=r'\s+', engine='python')

    result = df.xs('a', level=0)
    expected = df.xs('a')
    assert len(result) == 2
    tm.assert_frame_equal(result, expected)


def test_xs_level_series(multiindex_dataframe_random_data,
                         multiindex_year_month_day_dataframe_random_data):
    frame = multiindex_dataframe_random_data
    ymd = multiindex_year_month_day_dataframe_random_data
    s = frame['A']
    result = s[:, 'two']
    expected = frame.xs('two', level=1)['A']
    tm.assert_series_equal(result, expected)

    s = ymd['A']
    result = s[2000, 5]
    expected = ymd.loc[2000, 5]['A']
    tm.assert_series_equal(result, expected)

    # not implementing this for now

    pytest.raises(TypeError, s.__getitem__, (2000, slice(3, 4)))

    # result = s[2000, 3:4]
    # lv =s.index.get_level_values(1)
    # expected = s[(lv == 3) | (lv == 4)]
    # expected.index = expected.index.droplevel(0)
    # tm.assert_series_equal(result, expected)

    # can do this though
