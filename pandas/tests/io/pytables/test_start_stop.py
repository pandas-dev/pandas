import pytest

import numpy as np
import pandas as pd
from pandas import DataFrame
from pandas.tests.io.pytables.common import ensure_clean_store
import pandas.util.testing as tm


def test_start_stop_table():
    with ensure_clean_store() as store:
        # table
        df = DataFrame(dict(A=np.random.rand(20), B=np.random.rand(20)))
        store.append('df', df)

        result = store.select(
            'df', "columns=['A']", start=0, stop=5)
        expected = df.loc[0:4, ['A']]
        tm.assert_frame_equal(result, expected)

        # out of range
        result = store.select(
            'df', "columns=['A']", start=30, stop=40)
        assert len(result) == 0
        expected = df.loc[30:40, ['A']]
        tm.assert_frame_equal(result, expected)


def test_start_stop_multiple():
    # GH 16209
    with ensure_clean_store() as store:
        df = DataFrame({"foo": [1, 2], "bar": [1, 2]})

        store.append_to_multiple({'selector': ['foo'], 'data': None}, df,
                                 selector='selector')
        result = store.select_as_multiple(['selector', 'data'],
                                          selector='selector', start=0,
                                          stop=1)
        expected = df.loc[[0], ['foo', 'bar']]
        tm.assert_frame_equal(result, expected)


def test_start_stop_fixed():
    with ensure_clean_store() as store:
        # fixed, GH 8287
        df = DataFrame(dict(A=np.random.rand(20),
                            B=np.random.rand(20)),
                       index=pd.date_range('20130101', periods=20))
        store.put('df', df)

        result = store.select(
            'df', start=0, stop=5)
        expected = df.iloc[0:5, :]
        tm.assert_frame_equal(result, expected)

        result = store.select(
            'df', start=5, stop=10)
        expected = df.iloc[5:10, :]
        tm.assert_frame_equal(result, expected)

        # out of range
        result = store.select(
            'df', start=30, stop=40)
        expected = df.iloc[30:40, :]
        tm.assert_frame_equal(result, expected)

        # series
        s = df.A
        store.put('s', s)
        result = store.select('s', start=0, stop=5)
        expected = s.iloc[0:5]
        tm.assert_series_equal(result, expected)

        result = store.select('s', start=5, stop=10)
        expected = s.iloc[5:10]
        tm.assert_series_equal(result, expected)

        # sparse; not implemented
        df = tm.makeDataFrame()
        df.iloc[3:5, 1:3] = np.nan
        df.iloc[8:10, -2] = np.nan
        dfs = df.to_sparse()
        store.put('dfs', dfs)
        with pytest.raises(NotImplementedError):
            store.select('dfs', start=0, stop=5)
