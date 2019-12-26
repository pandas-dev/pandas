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
    for name, index in [("idx", idx)]:
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
                idx = type(index)(values, freq=index.freq)
            else:
                idx = type(index)(values)

            expected = np.array([False] * len(idx), dtype=bool)
            expected[1] = True
            tm.assert_numpy_array_equal(idx._isnan, expected)
            assert idx.hasnans is True


def test_dropna():
    # GH 6194
    idx = pd.MultiIndex.from_arrays(
        [
            [1, np.nan, 3, np.nan, 5],
            [1, 2, np.nan, np.nan, 5],
            ["a", "b", "c", np.nan, "e"],
        ]
    )

    exp = pd.MultiIndex.from_arrays([[1, 5], [1, 5], ["a", "e"]])
    tm.assert_index_equal(idx.dropna(), exp)
    tm.assert_index_equal(idx.dropna(how="any"), exp)

    exp = pd.MultiIndex.from_arrays(
        [[1, np.nan, 3, 5], [1, 2, np.nan, 5], ["a", "b", "c", "e"]]
    )
    tm.assert_index_equal(idx.dropna(how="all"), exp)

    msg = "invalid how option: xxx"
    with pytest.raises(ValueError, match=msg):
        idx.dropna(how="xxx")

    # GH26408
    # test if missing values are dropped for multiindex constructed
    # from codes and values
    idx = MultiIndex(
        levels=[[np.nan, None, pd.NaT, "128", 2], [np.nan, None, pd.NaT, "128", 2]],
        codes=[[0, -1, 1, 2, 3, 4], [0, -1, 3, 3, 3, 4]],
    )
    expected = MultiIndex.from_arrays([["128", 2], ["128", 2]])
    tm.assert_index_equal(idx.dropna(), expected)
    tm.assert_index_equal(idx.dropna(how="any"), expected)

    expected = MultiIndex.from_arrays(
        [[np.nan, np.nan, "128", 2], ["128", "128", "128", 2]]
    )
    tm.assert_index_equal(idx.dropna(how="all"), expected)


def test_nulls(idx):
    # this is really a smoke test for the methods
    # as these are adequately tested for function elsewhere

    msg = "isna is not defined for MultiIndex"
    with pytest.raises(NotImplementedError, match=msg):
        idx.isna()


@pytest.mark.xfail(reason="isna is not defined for MultiIndex")
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

    index = type(idx)(values)

    expected = np.array([False] * len(index), dtype=bool)
    expected[1] = True
    tm.assert_numpy_array_equal(index._isnan, expected)
    assert index.hasnans is True


def test_nan_stays_float():

    # GH 7031
    idx0 = pd.MultiIndex(
        levels=[["A", "B"], []], codes=[[1, 0], [-1, -1]], names=[0, 1]
    )
    idx1 = pd.MultiIndex(levels=[["C"], ["D"]], codes=[[0], [0]], names=[0, 1])
    idxm = idx0.join(idx1, how="outer")
    assert pd.isna(idx0.get_level_values(1)).all()
    # the following failed in 0.14.1
    assert pd.isna(idxm.get_level_values(1)[:-1]).all()

    df0 = pd.DataFrame([[1, 2]], index=idx0)
    df1 = pd.DataFrame([[3, 4]], index=idx1)
    dfm = df0 - df1
    assert pd.isna(df0.index.get_level_values(1)).all()
    # the following failed in 0.14.1
    assert pd.isna(dfm.index.get_level_values(1)[:-1]).all()
