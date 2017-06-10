""" test ipc compat """

import pytest
pyarrow = pytest.importorskip('pyarrow')

from distutils.version import LooseVersion
import numpy as np
import pandas as pd
from pandas import Series, Index, DataFrame
from pandas.io.ipc import (to_ipc, read_ipc,
                           _to_pickle, _to_pyarrow,
                           _read_pickle, _read_pyarrow)

import pandas.util.testing as tm

_HAVE_LATEST_PYARROW = LooseVersion(pyarrow.__version__) > '0.2.0'


@pytest.fixture(
    params=[('pickle', _to_pickle, _read_pickle),
            pytest.mark.skipif(not _HAVE_LATEST_PYARROW,
                               reason='need newer pyarrow version')(
                'pyarrow', _to_pyarrow, _read_pyarrow)],
    ids=lambda x: x[0])
def engine(request):
    return request.param


@pytest.fixture
def pa():
    if not _HAVE_LATEST_PYARROW:
        pytest.skip("need newer pyarrow")


def make_mixed_frame(N):
    return DataFrame(
        {'A': np.arange(N),
         'B': np.random.randn(N),
         'C': 'foo',
         'D': tm.makeStringIndex(N),
         'E': pd.Categorical.from_codes(np.repeat([0, 1], N // 2),
                                        categories=['foo', 'bar']),
         'F': pd.date_range('20130101', freq='s', periods=N)})


class TestIPC(object):

    def check_error_on_write(self, df, exc):
        # check that we are raising the exception
        # on writing

        with pytest.raises(exc):
            to_ipc(df)

    def check_round_trip(self, df, engine=None):

        if engine is None:
            writer = to_ipc
            reader = read_ipc
            b = writer(df)
        else:
            _, writer, reader = engine
            b = writer(df)

            # we are calling a lower-level routine
            b = b['data']

        result = reader(b)
        tm.assert_frame_equal(result, df)

    def test_error(self):
        for obj in [1, 'foo', pd.Timestamp('20130101'),
                    np.array([1, 2, 3])]:
            self.check_error_on_write(obj, ValueError)

    def test_with_small_size(self, engine):

        N = 100
        df = make_mixed_frame(N)
        self.check_round_trip(df, engine)

    def test_with_med_size(self, engine):

        # large size
        N = 10000
        df = make_mixed_frame(N)
        self.check_round_trip(df, engine)

    def test_with_large_size(self, engine):

        # large size
        N = 1000000
        df = make_mixed_frame(N)
        self.check_round_trip(df, engine)

    def test_non_dataframe(self):

        i = Index(['foo', 'bar'])
        b = to_ipc(i)
        result = read_ipc(b)
        tm.assert_index_equal(result, i)

        s = Series(['foo', 'bar'])
        b = to_ipc(s)
        result = read_ipc(b)
        tm.assert_series_equal(result, s)

    def test_basic(self, pa):

        df = pd.DataFrame({
            'string': list('abc'),
            'int': list(range(1, 4)),
            'uint': np.arange(3, 6).astype('u1'),
            'float': np.arange(4.0, 7.0, dtype='float64'),
            'bool': [True, False, True],
            'bool_with_nan': [True, None, True],
            'cat': pd.Categorical(list('abc')),
            'date_range': pd.date_range('20130101', periods=3),
            'date_range_tz': pd.date_range('20130101', periods=3,
                                           tz='US/Eastern'),
            'timedelta': pd.timedelta_range('1 day', periods=3)})

        # should work both on pickle & pyarrow
        # TODO: how to assure this?
        self.check_round_trip(df)

    def test_pickle_only(self):

        # period
        df = pd.DataFrame({'a': pd.period_range('2013', freq='M', periods=3)})
        self.check_round_trip(df)

        # non-strings
        df = pd.DataFrame({'a': ['a', 1, 2.0]})
        self.check_round_trip(df)

    def test_duplicate_columns(self, pa):

        df = pd.DataFrame(np.arange(12).reshape(4, 3),
                          columns=list('aaa')).copy()
        self.check_round_trip(df)

    def test_stringify_columns(self, pa):

        df = pd.DataFrame(np.arange(12).reshape(4, 3)).copy()
        self.check_round_trip(df)
