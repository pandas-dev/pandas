import numpy as np

from pandas import DataFrame
import pandas.util.testing as tm
from .common import _check_roundtrip, _check_double_roundtrip


def test_sparse_series():
    s = tm.makeStringSeries()
    s.iloc[3:5] = np.nan
    ss = s.to_sparse()
    _check_roundtrip(ss, tm.assert_series_equal, check_series_type=True)

    ss2 = s.to_sparse(kind='integer')
    _check_roundtrip(ss2, tm.assert_series_equal, check_series_type=True)

    ss3 = s.to_sparse(fill_value=0)
    _check_roundtrip(ss3, tm.assert_series_equal, check_series_type=True)


def test_sparse_frame():
    s = tm.makeDataFrame()
    s.iloc[3:5, 1:3] = np.nan
    s.iloc[8:10, -2] = np.nan
    ss = s.to_sparse()

    _check_double_roundtrip(ss, tm.assert_frame_equal, check_frame_type=True)

    ss2 = s.to_sparse(kind='integer')
    _check_double_roundtrip(ss2, tm.assert_frame_equal, check_frame_type=True)

    ss3 = s.to_sparse(fill_value=0)
    _check_double_roundtrip(ss3, tm.assert_frame_equal, check_frame_type=True)


def test_sparse_with_compression():
    # GH 2931
    # make sparse dataframe
    arr = np.random.binomial(n=1, p=.01, size=(1000, 10))
    df = DataFrame(arr).to_sparse(fill_value=0)

    # case 1: store uncompressed
    _check_double_roundtrip(df, tm.assert_frame_equal,
                            compression=False,
                            check_frame_type=True)

    # case 2: store compressed (works)
    _check_double_roundtrip(df, tm.assert_frame_equal,
                            compression='zlib',
                            check_frame_type=True)

    # set one series to be completely sparse
    df[0] = np.zeros(1000)

    # case 3: store df with completely sparse series uncompressed
    _check_double_roundtrip(df, tm.assert_frame_equal,
                            compression=False,
                            check_frame_type=True)

    # case 4: try storing df with completely sparse series compressed
    # (fails)
    _check_double_roundtrip(df, tm.assert_frame_equal,
                            compression='zlib',
                            check_frame_type=True)
