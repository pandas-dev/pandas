# -*- coding: utf-8 -*-

import pandas as pd
from pandas import MultiIndex
import pytest
from pandas.compat import PY3, PYPY, lrange, lzip, range, u
import numpy as np
import pandas.util.testing as tm


def test_contains_top_level():
    midx = MultiIndex.from_product([['A', 'B'], [1, 2]])
    assert 'A' in midx
    assert 'A' not in midx._engine


def test_contains_with_nat():
    # MI with a NaT
    mi = MultiIndex(levels=[['C'],
                            pd.date_range('2012-01-01', periods=5)],
                    labels=[[0, 0, 0, 0, 0, 0], [-1, 0, 1, 2, 3, 4]],
                    names=[None, 'B'])
    assert ('C', pd.Timestamp('2012-01-01')) in mi
    for val in mi.values:
        assert val in mi


def test_contains(_index):
    assert ('foo', 'two') in _index
    assert ('bar', 'two') not in _index
    assert None not in _index


@pytest.mark.skipif(not PYPY, reason="tuples cmp recursively on PyPy")
def test_isin_nan_pypy():
    idx = MultiIndex.from_arrays([['foo', 'bar'], [1.0, np.nan]])
    tm.assert_numpy_array_equal(idx.isin([('bar', np.nan)]),
                                np.array([False, True]))
    tm.assert_numpy_array_equal(idx.isin([('bar', float('nan'))]),
                                np.array([False, True]))


def test_isin():
    values = [('foo', 2), ('bar', 3), ('quux', 4)]

    idx = MultiIndex.from_arrays([['qux', 'baz', 'foo', 'bar'], np.arange(
        4)])
    result = idx.isin(values)
    expected = np.array([False, False, True, True])
    tm.assert_numpy_array_equal(result, expected)

    # empty, return dtype bool
    idx = MultiIndex.from_arrays([[], []])
    result = idx.isin(values)
    assert len(result) == 0
    assert result.dtype == np.bool_


@pytest.mark.skipif(PYPY, reason="tuples cmp recursively on PyPy")
def test_isin_nan_not_pypy():
    idx = MultiIndex.from_arrays([['foo', 'bar'], [1.0, np.nan]])
    tm.assert_numpy_array_equal(idx.isin([('bar', np.nan)]),
                                np.array([False, False]))
    tm.assert_numpy_array_equal(idx.isin([('bar', float('nan'))]),
                                np.array([False, False]))


def test_isin_level_kwarg():
    idx = MultiIndex.from_arrays([['qux', 'baz', 'foo', 'bar'], np.arange(
        4)])

    vals_0 = ['foo', 'bar', 'quux']
    vals_1 = [2, 3, 10]

    expected = np.array([False, False, True, True])
    tm.assert_numpy_array_equal(expected, idx.isin(vals_0, level=0))
    tm.assert_numpy_array_equal(expected, idx.isin(vals_0, level=-2))

    tm.assert_numpy_array_equal(expected, idx.isin(vals_1, level=1))
    tm.assert_numpy_array_equal(expected, idx.isin(vals_1, level=-1))

    pytest.raises(IndexError, idx.isin, vals_0, level=5)
    pytest.raises(IndexError, idx.isin, vals_0, level=-5)

    pytest.raises(KeyError, idx.isin, vals_0, level=1.0)
    pytest.raises(KeyError, idx.isin, vals_1, level=-1.0)
    pytest.raises(KeyError, idx.isin, vals_1, level='A')

    idx.names = ['A', 'B']
    tm.assert_numpy_array_equal(expected, idx.isin(vals_0, level='A'))
    tm.assert_numpy_array_equal(expected, idx.isin(vals_1, level='B'))

    pytest.raises(KeyError, idx.isin, vals_1, level='C')


def test_hasnans_isnans(named_index):
    # GH 11343, added tests for hasnans / isnans
    for name, index in named_index.items():
        if isinstance(index, MultiIndex):
            pass
        else:
            idx = index.copy()

            # cases in indices doesn't include NaN
            expected = np.array([False] * len(idx), dtype=bool)
            tm.assert_numpy_array_equal(idx._isnan, expected)
            assert not idx.hasnans

            idx = index.copy()
            values = idx.values

            if len(index) == 0:
                continue
            elif isinstance(index, DatetimeIndexOpsMixin):
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
            assert idx.hasnans
