# -*- coding: utf-8 -*-

import numpy as np
import pytest

from pandas._libs.tslib import iNaT

import pandas as pd
from pandas import Int64Index, MultiIndex, PeriodIndex, UInt64Index
from pandas.core.indexes.datetimelike import DatetimeIndexOpsMixin
import pandas.util.testing as tm


def test_fillna(idx):
    # GH 11343

    # TODO: Remove or Refactor.  Not Implemented for MultiIndex
    for name, index in [('idx', idx), ]:
        if len(index) == 0:
            pass
        elif isinstance(index, MultiIndex):
            idx = index.copy()
            msg = "isna is not defined for MultiIndex"
            with pytest.raises(NotImplementedError, match=msg):
                idx.fillna(idx[0])
        else:
            idx = index.copy()
            result = idx.fillna(idx[0])
            tm.assert_index_equal(result, idx)
            assert result is not idx

            msg = "'value' must be a scalar, passed: "
            with pytest.raises(TypeError, match=msg):
                idx.fillna([idx[0]])

            idx = index.copy()
            values = idx.values

            if isinstance(index, DatetimeIndexOpsMixin):
                values[1] = iNaT
            elif isinstance(index, (Int64Index, UInt64Index)):
                continue
            else:
                values[1] = np.nan

            if isinstance(index, PeriodIndex):
                idx = index.__class__(values, freq=index.freq)
            else:
                idx = index.__class__(values)

            expected = np.array([False] * len(idx), dtype=bool)
            expected[1] = True
            tm.assert_numpy_array_equal(idx._isnan, expected)
            assert idx.hasnans is True


def test_dropna():
    # GH 6194
    idx = pd.MultiIndex.from_arrays([[1, np.nan, 3, np.nan, 5],
                                     [1, 2, np.nan, np.nan, 5],
                                     ['a', 'b', 'c', np.nan, 'e']])

    exp = pd.MultiIndex.from_arrays([[1, 5],
                                     [1, 5],
                                     ['a', 'e']])
    tm.assert_index_equal(idx.dropna(), exp)
    tm.assert_index_equal(idx.dropna(how='any'), exp)

    exp = pd.MultiIndex.from_arrays([[1, np.nan, 3, 5],
                                     [1, 2, np.nan, 5],
                                     ['a', 'b', 'c', 'e']])
    tm.assert_index_equal(idx.dropna(how='all'), exp)

    msg = "invalid how option: xxx"
    with pytest.raises(ValueError, match=msg):
        idx.dropna(how='xxx')


def test_nulls(idx):
    # this is really a smoke test for the methods
    # as these are adequately tested for function elsewhere

    msg = "isna is not defined for MultiIndex"
    with pytest.raises(NotImplementedError, match=msg):
        idx.isna()


@pytest.mark.xfail
def test_hasnans_isnans(idx):
    # GH 11343, added tests for hasnans / isnans
    index = idx.copy()

    # cases in indices doesn't include NaN
    expected = np.array([False] * len(index), dtype=bool)
    tm.assert_numpy_array_equal(index._isnan, expected)
    assert index.hasnans is False

    index = idx.copy()
    values = index.values
    values[1] = np.nan

    index = idx.__class__(values)

    expected = np.array([False] * len(index), dtype=bool)
    expected[1] = True
    tm.assert_numpy_array_equal(index._isnan, expected)
    assert index.hasnans is True


def test_nan_stays_float():

    # GH 7031
    idx0 = pd.MultiIndex(levels=[["A", "B"], []],
                         codes=[[1, 0], [-1, -1]],
                         names=[0, 1])
    idx1 = pd.MultiIndex(levels=[["C"], ["D"]],
                         codes=[[0], [0]],
                         names=[0, 1])
    idxm = idx0.join(idx1, how='outer')
    assert pd.isna(idx0.get_level_values(1)).all()
    # the following failed in 0.14.1
    assert pd.isna(idxm.get_level_values(1)[:-1]).all()

    df0 = pd.DataFrame([[1, 2]], index=idx0)
    df1 = pd.DataFrame([[3, 4]], index=idx1)
    dfm = df0 - df1
    assert pd.isna(df0.index.get_level_values(1)).all()
    # the following failed in 0.14.1
    assert pd.isna(dfm.index.get_level_values(1)[:-1]).all()


def test_nan_multi_index():
    # GH 22247
    # When using the MultiIndex features of pandas, when an `np.nan`
    # is in the index when new values are added to the DF then the
    # values are not `np.nan`, but copied from the `np.nan` row.
    df = pd.DataFrame(
        [
            ['A', np.nan, 1.23, 4.56],
            ['A', 'G', 1.23, 4.56],
            ['A', 'D', 9.87, 10.54],
        ],
        columns=['pivot_0', 'pivot_1', 'col_1', 'col_2'],
    )
    df.set_index(['pivot_0', 'pivot_1'], inplace=True)
    pivot_0 = 'A'
    pivot_1_values = ['D', 'E', 'F']
    for value in pivot_1_values:
        if value not in df.index.get_level_values('pivot_1').tolist():
            df.at[(pivot_0, value), 'col_2'] = 0.0

    assert df.loc[('A', 'F')]['col_2'] == 0.0  # Pass
    # Fails: value of 1.23 from the first row in the df is copied
    # This behavior shows for all versions v0.23.x, however is fine for 0.22.0.
    assert pd.isna(df.loc[('A', 'F')]['col_1'])


def test_nan_set_value_multi_index():
    # GH 22247
    # When using the MultiIndex features of pandas, when an `np.nan`
    # is in the index when new values are added to the DF then the
    # values are not `np.nan`, but copied from the `np.nan` row.
    df = pd.DataFrame(
        [
            ['A', 'G', 1.23, 4.56],
            ['A', 'D', 9.87, 10.54],
        ],
        columns=['pivot_0', 'pivot_1', 'col_1', 'col_2'],
    )
    df.set_index(['pivot_0', 'pivot_1'], inplace=True)
    df.at[('A', 'E'), 'col_2'] = 0.0
    df.at[('A', 'F'), 'col_2'] = 0.0
    # Fails: raise exception
    # This behavior shows for all versions v0.23.x, however is fine for 0.22.0.
    df.at[('A', np.nan), 'col_2'] = 0.0

    assert df.loc[('A', np.nan)]['col_2'] == 0.0
    assert pd.isna(df.loc[('A', np.nan)]['col_1'])


def test_nan_sigle_index():
    # GH 22247
    df = pd.DataFrame(
        [
            [np.nan, 1.23, 4.56],
            ['G', 1.23, 4.56],
            ['D', 9.87, 10.54],
        ],
        columns=['pivot_0', 'col_1', 'col_2'],
    )
    df.set_index(['pivot_0'], inplace=True)

    pivot_0_values = ['D', 'E', 'F']
    for value in pivot_0_values:
        if value not in df.index.get_level_values('pivot_0').tolist():
            df.at[(value), 'col_2'] = 0.0

    assert df.loc[('F')]['col_2'] == 0.0
    assert pd.isna(df.loc[('F')]['col_1'])
