# pylint: disable=E1103

from collections import OrderedDict

import numpy as np
import pytest
from numpy import nan
from numpy.random import randn

import pandas as pd
import pandas.util.testing as tm
from pandas import (DataFrame, Index, MultiIndex, Series)
from pandas.core.reshape.concat import concat
from pandas.core.reshape.merge import merge


@pytest.fixture
def left():
    """left dataframe (not multi-indexed) for multi-index join tests"""
    # a little relevant example with NAs
    key1 = ['bar', 'bar', 'bar', 'foo', 'foo', 'baz', 'baz', 'qux',
            'qux', 'snap']
    key2 = ['two', 'one', 'three', 'one', 'two', 'one', 'two', 'two',
            'three', 'one']

    data = np.random.randn(len(key1))
    return DataFrame({'key1': key1, 'key2': key2, 'data': data})


@pytest.fixture
def right():
    """right dataframe (multi-indexed) for multi-index join tests"""
    index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                               ['one', 'two', 'three']],
                       labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                               [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                       names=['key1', 'key2'])

    return DataFrame(np.random.randn(10, 3), index=index,
                     columns=['j_one', 'j_two', 'j_three'])


class TestMergeMulti(object):

    def test_merge_on_multikey(self, left, right, join_type):
        on_cols = ['key1', 'key2']
        result = (left.join(right, on=on_cols, how=join_type)
                  .reset_index(drop=True))

        expected = pd.merge(left, right.reset_index(),
                            on=on_cols, how=join_type)

        tm.assert_frame_equal(result, expected)

        result = (left.join(right, on=on_cols, how=join_type, sort=True)
                  .reset_index(drop=True))

        expected = pd.merge(left, right.reset_index(),
                            on=on_cols, how=join_type, sort=True)

        tm.assert_frame_equal(result, expected)

    """
    @pytest.mark.parametrize("sort", [False, True])
    def test_left_join_multi_index(self, sort):

        icols = ['1st', '2nd', '3rd']

        def bind_cols(df):
            def iord(a):
                return 0 if a != a else ord(a)

            def f(ts):
                return ts.map(iord) - ord('a')

            return (f(df['1st']) + f(df['3rd']) * 1e2 +
                    df['2nd'].fillna(0) * 1e4)

        def run_asserts(left, right, sort):
            res = left.join(right, on=icols, how='left', sort=sort)

            assert len(left) < len(res) + 1
            assert not res['4th'].isna().any()
            assert not res['5th'].isna().any()

            tm.assert_series_equal(res['4th'], - res['5th'],
                                   check_names=False)
            result = bind_cols(res.iloc[:, :-2])

            tm.assert_series_equal(res['4th'], result, check_names=False)
            assert result.name is None

            if sort:
                tm.assert_frame_equal(
                    res, res.sort_values(icols, kind='mergesort'))

            out = merge(left, right.reset_index(), on=icols,
                        sort=sort, how='left')

            res.index = np.arange(len(res))
            tm.assert_frame_equal(out, res)

        lc = list(map(chr, np.arange(ord('a'), ord('z') + 1)))
        left = DataFrame(np.random.choice(lc, (5000, 2)),
                         columns=['1st', '3rd'])
        left.insert(1, '2nd', np.random.randint(0, 1000, len(left)))

        i = np.random.permutation(len(left))
        right = left.iloc[i].copy()

        left['4th'] = bind_cols(left)
        right['5th'] = - bind_cols(right)
        right.set_index(icols, inplace=True)

        run_asserts(left, right)

        # inject some nulls
        left.loc[1::23, '1st'] = np.nan
        left.loc[2::37, '2nd'] = np.nan
        left.loc[3::43, '3rd'] = np.nan
        left['4th'] = bind_cols(left)

        i = np.random.permutation(len(left))
        right = left.iloc[i, :-1]
        right['5th'] = - bind_cols(right)
        right.set_index(icols, inplace=True)

        run_asserts(left, right)
    """

    @pytest.mark.parametrize("sort", [False, True])
    def test_merge_right_vs_left(self, left, right, sort):
        # compare left vs right merge with multikey
        on_cols = ['key1', 'key2']
        merged_left_right = left.merge(right,
                                       left_on=on_cols, right_index=True,
                                       how='left', sort=sort)

        merge_right_left = right.merge(left,
                                       right_on=on_cols, left_index=True,
                                       how='right', sort=sort)

        # Reorder columns
        merge_right_left = merge_right_left[merged_left_right.columns]

        tm.assert_frame_equal(merged_left_right, merge_right_left)

    def test_compress_group_combinations(self):

        # ~ 40000000 possible unique groups
        key1 = tm.rands_array(10, 10000)
        key1 = np.tile(key1, 2)
        key2 = key1[::-1]

        df = DataFrame({'key1': key1, 'key2': key2,
                        'value1': np.random.randn(20000)})

        df2 = DataFrame({'key1': key1[::2], 'key2': key2[::2],
                         'value2': np.random.randn(10000)})

        # just to hit the label compression code path
        merge(df, df2, how='outer')

    def test_left_join_index_preserve_order(self):

        on_cols = ['k1', 'k2']
        left = DataFrame({'k1': [0, 1, 2] * 8,
                          'k2': ['foo', 'bar'] * 12,
                          'v': np.array(np.arange(24), dtype=np.int64)})

        index = MultiIndex.from_tuples([(2, 'bar'), (1, 'foo')])
        right = DataFrame({'v2': [5, 7]}, index=index)

        result = left.join(right, on=on_cols)

        expected = left.copy()
        expected['v2'] = np.nan
        expected.loc[(expected.k1 == 2) & (expected.k2 == 'bar'), 'v2'] = 5
        expected.loc[(expected.k1 == 1) & (expected.k2 == 'foo'), 'v2'] = 7

        tm.assert_frame_equal(result, expected)

        result.sort_values(on_cols, kind='mergesort', inplace=True)
        expected = left.join(right, on=on_cols, sort=True)

        tm.assert_frame_equal(result, expected)

        # test join with multi dtypes blocks
        left = DataFrame({'k1': [0, 1, 2] * 8,
                          'k2': ['foo', 'bar'] * 12,
                          'k3': np.array([0, 1, 2] * 8, dtype=np.float32),
                          'v': np.array(np.arange(24), dtype=np.int32)})

        index = MultiIndex.from_tuples([(2, 'bar'), (1, 'foo')])
        right = DataFrame({'v2': [5, 7]}, index=index)

        result = left.join(right, on=on_cols)

        expected = left.copy()
        expected['v2'] = np.nan
        expected.loc[(expected.k1 == 2) & (expected.k2 == 'bar'), 'v2'] = 5
        expected.loc[(expected.k1 == 1) & (expected.k2 == 'foo'), 'v2'] = 7

        tm.assert_frame_equal(result, expected)

        result = result.sort_values(on_cols, kind='mergesort')
        expected = left.join(right, on=on_cols, sort=True)

        tm.assert_frame_equal(result, expected)

    def test_left_join_index_multi_match_multiindex(self):
        left = DataFrame([
            ['X', 'Y', 'C', 'a'],
            ['W', 'Y', 'C', 'e'],
            ['V', 'Q', 'A', 'h'],
            ['V', 'R', 'D', 'i'],
            ['X', 'Y', 'D', 'b'],
            ['X', 'Y', 'A', 'c'],
            ['W', 'Q', 'B', 'f'],
            ['W', 'R', 'C', 'g'],
            ['V', 'Y', 'C', 'j'],
            ['X', 'Y', 'B', 'd']],
            columns=['cola', 'colb', 'colc', 'tag'],
            index=[3, 2, 0, 1, 7, 6, 4, 5, 9, 8])

        right = DataFrame([
            ['W', 'R', 'C', 0],
            ['W', 'Q', 'B', 3],
            ['W', 'Q', 'B', 8],
            ['X', 'Y', 'A', 1],
            ['X', 'Y', 'A', 4],
            ['X', 'Y', 'B', 5],
            ['X', 'Y', 'C', 6],
            ['X', 'Y', 'C', 9],
            ['X', 'Q', 'C', -6],
            ['X', 'R', 'C', -9],
            ['V', 'Y', 'C', 7],
            ['V', 'R', 'D', 2],
            ['V', 'R', 'D', -1],
            ['V', 'Q', 'A', -3]],
            columns=['col1', 'col2', 'col3', 'val'])
        right.set_index(['col1', 'col2', 'col3'], inplace=True)

        result = left.join(right, on=['cola', 'colb', 'colc'], how='left')

        expected = DataFrame([
            ['X', 'Y', 'C', 'a', 6],
            ['X', 'Y', 'C', 'a', 9],
            ['W', 'Y', 'C', 'e', nan],
            ['V', 'Q', 'A', 'h', -3],
            ['V', 'R', 'D', 'i', 2],
            ['V', 'R', 'D', 'i', -1],
            ['X', 'Y', 'D', 'b', nan],
            ['X', 'Y', 'A', 'c', 1],
            ['X', 'Y', 'A', 'c', 4],
            ['W', 'Q', 'B', 'f', 3],
            ['W', 'Q', 'B', 'f', 8],
            ['W', 'R', 'C', 'g', 0],
            ['V', 'Y', 'C', 'j', 7],
            ['X', 'Y', 'B', 'd', 5]],
            columns=['cola', 'colb', 'colc', 'tag', 'val'],
            index=[3, 3, 2, 0, 1, 1, 7, 6, 6, 4, 4, 5, 9, 8])

        tm.assert_frame_equal(result, expected)

        result = left.join(right, on=['cola', 'colb', 'colc'],
                           how='left', sort=True)

        expected = expected.sort_values(['cola', 'colb', 'colc'],
                                        kind='mergesort')

        tm.assert_frame_equal(result, expected)

    def test_left_join_index_multi_match(self):
        left = DataFrame([
            ['c', 0],
            ['b', 1],
            ['a', 2],
            ['b', 3]],
            columns=['tag', 'val'],
            index=[2, 0, 1, 3])

        right = DataFrame([
            ['a', 'v'],
            ['c', 'w'],
            ['c', 'x'],
            ['d', 'y'],
            ['a', 'z'],
            ['c', 'r'],
            ['e', 'q'],
            ['c', 's']],
            columns=['tag', 'char'])

        right.set_index('tag', inplace=True)
        result = left.join(right, on='tag', how='left')

        expected = DataFrame([
            ['c', 0, 'w'],
            ['c', 0, 'x'],
            ['c', 0, 'r'],
            ['c', 0, 's'],
            ['b', 1, nan],
            ['a', 2, 'v'],
            ['a', 2, 'z'],
            ['b', 3, nan]],
            columns=['tag', 'val', 'char'],
            index=[2, 2, 2, 2, 0, 1, 1, 3])

        tm.assert_frame_equal(result, expected)

        result = left.join(right, on='tag', how='left', sort=True)
        tm.assert_frame_equal(
            result, expected.sort_values('tag', kind='mergesort'))

        # GH7331 - maintain left frame order in left merge
        result = merge(left, right.reset_index(), how='left', on='tag')
        expected.index = np.arange(len(expected))
        tm.assert_frame_equal(result, expected)

    def test_left_merge_na_buglet(self):
        left = DataFrame({'id': list('abcde'), 'v1': randn(5),
                          'v2': randn(5), 'dummy': list('abcde'),
                          'v3': randn(5)},
                         columns=['id', 'v1', 'v2', 'dummy', 'v3'])
        right = DataFrame({'id': ['a', 'b', np.nan, np.nan, np.nan],
                           'sv3': [1.234, 5.678, np.nan, np.nan, np.nan]})

        result = merge(left, right, on='id', how='left')

        rdf = right.drop(['id'], axis=1)
        expected = left.join(rdf)
        tm.assert_frame_equal(result, expected)

    def test_merge_na_keys(self):
        data = [[1950, "A", 1.5],
                [1950, "B", 1.5],
                [1955, "B", 1.5],
                [1960, "B", np.nan],
                [1970, "B", 4.],
                [1950, "C", 4.],
                [1960, "C", np.nan],
                [1965, "C", 3.],
                [1970, "C", 4.]]

        frame = DataFrame(data, columns=["year", "panel", "data"])

        other_data = [[1960, 'A', np.nan],
                      [1970, 'A', np.nan],
                      [1955, 'A', np.nan],
                      [1965, 'A', np.nan],
                      [1965, 'B', np.nan],
                      [1955, 'C', np.nan]]
        other = DataFrame(other_data, columns=['year', 'panel', 'data'])

        result = frame.merge(other, how='outer')

        expected = frame.fillna(-999).merge(other.fillna(-999), how='outer')
        expected = expected.replace(-999, np.nan)

        tm.assert_frame_equal(result, expected)

    def test_join_multi_levels(self):

        # GH 3662
        # merge multi-levels
        household = (
            DataFrame(
                dict(household_id=[1, 2, 3],
                     male=[0, 1, 0],
                     wealth=[196087.3, 316478.7, 294750]),
                columns=['household_id', 'male', 'wealth'])
            .set_index('household_id'))
        portfolio = (
            DataFrame(
                dict(household_id=[1, 2, 2, 3, 3, 3, 4],
                     asset_id=["nl0000301109", "nl0000289783", "gb00b03mlx29",
                               "gb00b03mlx29", "lu0197800237", "nl0000289965",
                               np.nan],
                     name=["ABN Amro", "Robeco", "Royal Dutch Shell",
                           "Royal Dutch Shell",
                           "AAB Eastern Europe Equity Fund",
                           "Postbank BioTech Fonds", np.nan],
                     share=[1.0, 0.4, 0.6, 0.15, 0.6, 0.25, 1.0]),
                columns=['household_id', 'asset_id', 'name', 'share'])
            .set_index(['household_id', 'asset_id']))
        result = household.join(portfolio, how='inner')
        expected = (
            DataFrame(
                dict(male=[0, 1, 1, 0, 0, 0],
                     wealth=[196087.3, 316478.7, 316478.7,
                             294750.0, 294750.0, 294750.0],
                     name=['ABN Amro', 'Robeco', 'Royal Dutch Shell',
                           'Royal Dutch Shell',
                           'AAB Eastern Europe Equity Fund',
                           'Postbank BioTech Fonds'],
                     share=[1.00, 0.40, 0.60, 0.15, 0.60, 0.25],
                     household_id=[1, 2, 2, 3, 3, 3],
                     asset_id=['nl0000301109', 'nl0000289783', 'gb00b03mlx29',
                               'gb00b03mlx29', 'lu0197800237',
                               'nl0000289965']))
            .set_index(['household_id', 'asset_id'])
            .reindex(columns=['male', 'wealth', 'name', 'share']))
        tm.assert_frame_equal(result, expected)

        # equivalency
        result = (merge(household.reset_index(), portfolio.reset_index(),
                        on=['household_id'], how='inner')
                  .set_index(['household_id', 'asset_id']))
        tm.assert_frame_equal(result, expected)

        result = household.join(portfolio, how='outer')
        expected = (concat([
            expected,
            (DataFrame(
                dict(share=[1.00]),
                index=MultiIndex.from_tuples(
                    [(4, np.nan)],
                    names=['household_id', 'asset_id'])))
        ], axis=0, sort=True).reindex(columns=expected.columns))
        tm.assert_frame_equal(result, expected)

        # invalid cases
        household.index.name = 'foo'

        def f():
            household.join(portfolio, how='inner')

        pytest.raises(ValueError, f)

        portfolio2 = portfolio.copy()
        portfolio2.index.set_names(['household_id', 'foo'])

        def f():
            portfolio2.join(portfolio, how='inner')

        pytest.raises(ValueError, f)

    def test_join_multi_levels2(self):

        # some more advanced merges
        # GH6360
        household = (
            DataFrame(
                dict(household_id=[1, 2, 2, 3, 3, 3, 4],
                     asset_id=["nl0000301109", "nl0000301109", "gb00b03mlx29",
                               "gb00b03mlx29", "lu0197800237", "nl0000289965",
                               np.nan],
                     share=[1.0, 0.4, 0.6, 0.15, 0.6, 0.25, 1.0]),
                columns=['household_id', 'asset_id', 'share'])
            .set_index(['household_id', 'asset_id']))

        log_return = DataFrame(dict(
            asset_id=["gb00b03mlx29", "gb00b03mlx29",
                      "gb00b03mlx29", "lu0197800237", "lu0197800237"],
            t=[233, 234, 235, 180, 181],
            log_return=[.09604978, -.06524096, .03532373, .03025441, .036997]
        )).set_index(["asset_id", "t"])

        expected = (
            DataFrame(dict(
                household_id=[2, 2, 2, 3, 3, 3, 3, 3],
                asset_id=["gb00b03mlx29", "gb00b03mlx29",
                          "gb00b03mlx29", "gb00b03mlx29",
                          "gb00b03mlx29", "gb00b03mlx29",
                          "lu0197800237", "lu0197800237"],
                t=[233, 234, 235, 233, 234, 235, 180, 181],
                share=[0.6, 0.6, 0.6, 0.15, 0.15, 0.15, 0.6, 0.6],
                log_return=[.09604978, -.06524096, .03532373,
                            .09604978, -.06524096, .03532373,
                            .03025441, .036997]
            ))
            .set_index(["household_id", "asset_id", "t"])
            .reindex(columns=['share', 'log_return']))

        # this is the equivalency
        result = (merge(household.reset_index(), log_return.reset_index(),
                        on=['asset_id'], how='inner')
                  .set_index(['household_id', 'asset_id', 't']))
        tm.assert_frame_equal(result, expected)

        expected = (
            DataFrame(dict(
                household_id=[1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4],
                asset_id=["nl0000301109", "nl0000301109", "gb00b03mlx29",
                          "gb00b03mlx29", "gb00b03mlx29",
                          "gb00b03mlx29", "gb00b03mlx29", "gb00b03mlx29",
                          "lu0197800237", "lu0197800237",
                          "nl0000289965", None],
                t=[None, None, 233, 234, 235, 233, 234,
                   235, 180, 181, None, None],
                share=[1.0, 0.4, 0.6, 0.6, 0.6, 0.15,
                       0.15, 0.15, 0.6, 0.6, 0.25, 1.0],
                log_return=[None, None, .09604978, -.06524096, .03532373,
                            .09604978, -.06524096, .03532373,
                            .03025441, .036997, None, None]
            ))
            .set_index(["household_id", "asset_id", "t"])
            .reindex(columns=['share', 'log_return']))

        result = (merge(household.reset_index(), log_return.reset_index(),
                  on=['asset_id'], how='outer')
                  .set_index(['household_id', 'asset_id', 't']))

        tm.assert_frame_equal(result, expected)


@pytest.fixture
def left_multi():
    return (
        DataFrame(
            dict(Origin=['A', 'A', 'B', 'B', 'C'],
                 Destination=['A', 'B', 'A', 'C', 'A'],
                 Period=['AM', 'AM', 'IP', 'AM', 'OP'],
                 TripPurp=['hbw', 'nhb', 'hbo', 'nhb', 'hbw'],
                 Trips=[1987, 3647, 2470, 4296, 4444]),
            columns=['Origin', 'Destination', 'Period',
                     'TripPurp', 'Trips'])
        .set_index(['Origin', 'Destination', 'Period', 'TripPurp']))


@pytest.fixture
def right_multi():
    return (
        DataFrame(
            dict(Origin=['A', 'A', 'B', 'B', 'C', 'C', 'E'],
                 Destination=['A', 'B', 'A', 'B', 'A', 'B', 'F'],
                 Period=['AM', 'AM', 'IP', 'AM', 'OP', 'IP', 'AM'],
                 LinkType=['a', 'b', 'c', 'b', 'a', 'b', 'a'],
                 Distance=[100, 80, 90, 80, 75, 35, 55]),
            columns=['Origin', 'Destination', 'Period',
                     'LinkType', 'Distance'])
        .set_index(['Origin', 'Destination', 'Period', 'LinkType']))


@pytest.fixture
def on_cols():
    return ['Origin', 'Destination', 'Period']


@pytest.fixture
def idx_cols():
    return ['Origin', 'Destination', 'Period', 'TripPurp', 'LinkType']


class TestJoinMultiMulti(object):

    def test_join_multi_multi(self, left_multi, right_multi, join_type,
                              on_cols, idx_cols):
        # Multi-index join tests
        expected = (pd.merge(left_multi.reset_index(),
                             right_multi.reset_index(),
                             how=join_type, on=on_cols).set_index(idx_cols)
                    .sort_index())

        result = left_multi.join(right_multi, how=join_type).sort_index()
        tm.assert_frame_equal(result, expected)

    def test_join_multi_empty_frames(self, left_multi, right_multi, join_type,
                                     on_cols, idx_cols):

        left_multi = left_multi.drop(columns=left_multi.columns)
        right_multi = right_multi.drop(columns=right_multi.columns)

        expected = (pd.merge(left_multi.reset_index(),
                             right_multi.reset_index(),
                             how=join_type, on=on_cols).set_index(idx_cols)
                    .sort_index())

        result = left_multi.join(right_multi, how=join_type).sort_index()
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("box", [None, np.asarray, Series, Index])
    def test_merge_datetime_index(self, box):
        # see gh-19038
        df = DataFrame([1, 2, 3],
                       ["2016-01-01", "2017-01-01", "2018-01-01"],
                       columns=["a"])
        df.index = pd.to_datetime(df.index)
        on_vector = df.index.year

        if box is not None:
            on_vector = box(on_vector)

        expected = DataFrame(
            OrderedDict([
                ("a", [1, 2, 3]),
                ("key_1", [2016, 2017, 2018]),
            ])
        )

        result = df.merge(df, on=["a", on_vector], how="inner")
        tm.assert_frame_equal(result, expected)

        expected = DataFrame(
            OrderedDict([
                ("key_0", [2016, 2017, 2018]),
                ("a_x", [1, 2, 3]),
                ("a_y", [1, 2, 3]),
            ])
        )

        result = df.merge(df, on=[df.index.year], how="inner")
        tm.assert_frame_equal(result, expected)

    def test_single_common_level(self):
        index_left = pd.MultiIndex.from_tuples([('K0', 'X0'), ('K0', 'X1'),
                                                ('K1', 'X2')],
                                               names=['key', 'X'])

        left = pd.DataFrame({'A': ['A0', 'A1', 'A2'],
                             'B': ['B0', 'B1', 'B2']},
                            index=index_left)

        index_right = pd.MultiIndex.from_tuples([('K0', 'Y0'), ('K1', 'Y1'),
                                                 ('K2', 'Y2'), ('K2', 'Y3')],
                                                names=['key', 'Y'])

        right = pd.DataFrame({'C': ['C0', 'C1', 'C2', 'C3'],
                              'D': ['D0', 'D1', 'D2', 'D3']},
                             index=index_right)

        result = left.join(right)
        expected = (pd.merge(left.reset_index(), right.reset_index(),
                             on=['key'], how='inner')
                    .set_index(['key', 'X', 'Y']))

        tm.assert_frame_equal(result, expected)
