""" test parquet compat """

import pytest
import datetime
from distutils.version import LooseVersion
from warnings import catch_warnings

import numpy as np
import pandas as pd
from pandas.compat import PY3
from pandas.io.parquet import (to_parquet, read_parquet, get_engine,
                               PyArrowImpl, FastParquetImpl)
from pandas.util import testing as tm

try:
    import pyarrow  # noqa
    _HAVE_PYARROW = True
except ImportError:
    _HAVE_PYARROW = False

try:
    import fastparquet  # noqa
    _HAVE_FASTPARQUET = True
except ImportError:
    _HAVE_FASTPARQUET = False


# setup engines & skips
@pytest.fixture(params=[
    pytest.param('fastparquet',
                 marks=pytest.mark.skipif(not _HAVE_FASTPARQUET,
                                          reason='fastparquet is '
                                                 'not installed')),
    pytest.param('pyarrow',
                 marks=pytest.mark.skipif(not _HAVE_PYARROW,
                                          reason='pyarrow is '
                                                 'not installed'))])
def engine(request):
    return request.param


@pytest.fixture
def pa():
    if not _HAVE_PYARROW:
        pytest.skip("pyarrow is not installed")
    return 'pyarrow'


@pytest.fixture
def pa_lt_070():
    if not _HAVE_PYARROW:
        pytest.skip("pyarrow is not installed")
    if LooseVersion(pyarrow.__version__) >= '0.7.0':
        pytest.skip("pyarrow is >= 0.7.0")
    return 'pyarrow'


@pytest.fixture
def pa_ge_070():
    if not _HAVE_PYARROW:
        pytest.skip("pyarrow is not installed")
    if LooseVersion(pyarrow.__version__) < '0.7.0':
        pytest.skip("pyarrow is < 0.7.0")
    return 'pyarrow'


@pytest.fixture
def fp():
    if not _HAVE_FASTPARQUET:
        pytest.skip("fastparquet is not installed")
    return 'fastparquet'


@pytest.fixture
def df_compat():
    return pd.DataFrame({'A': [1, 2, 3], 'B': 'foo'})


@pytest.fixture
def df_cross_compat():
    df = pd.DataFrame({'a': list('abc'),
                       'b': list(range(1, 4)),
                       'c': np.arange(3, 6).astype('u1'),
                       'd': np.arange(4.0, 7.0, dtype='float64'),
                       'e': [True, False, True],
                       'f': pd.date_range('20130101', periods=3),
                       'g': pd.date_range('20130101', periods=3,
                                          tz='US/Eastern'),
                       'h': pd.date_range('20130101', periods=3, freq='ns')})
    return df


def test_invalid_engine(df_compat):

    with pytest.raises(ValueError):
        df_compat.to_parquet('foo', 'bar')


def test_options_py(df_compat, pa):
    # use the set option

    df = df_compat
    with tm.ensure_clean() as path:

        with pd.option_context('io.parquet.engine', 'pyarrow'):
            df.to_parquet(path)

            result = read_parquet(path, compression=None)
            tm.assert_frame_equal(result, df)


def test_options_fp(df_compat, fp):
    # use the set option

    df = df_compat
    with tm.ensure_clean() as path:

        with pd.option_context('io.parquet.engine', 'fastparquet'):
            df.to_parquet(path, compression=None)

            result = read_parquet(path, compression=None)
            tm.assert_frame_equal(result, df)


def test_options_auto(df_compat, fp, pa):

    df = df_compat
    with tm.ensure_clean() as path:

        with pd.option_context('io.parquet.engine', 'auto'):
            df.to_parquet(path)

            result = read_parquet(path, compression=None)
            tm.assert_frame_equal(result, df)


def test_options_get_engine(fp, pa):
    assert isinstance(get_engine('pyarrow'), PyArrowImpl)
    assert isinstance(get_engine('fastparquet'), FastParquetImpl)

    with pd.option_context('io.parquet.engine', 'pyarrow'):
        assert isinstance(get_engine('auto'), PyArrowImpl)
        assert isinstance(get_engine('pyarrow'), PyArrowImpl)
        assert isinstance(get_engine('fastparquet'), FastParquetImpl)

    with pd.option_context('io.parquet.engine', 'fastparquet'):
        assert isinstance(get_engine('auto'), FastParquetImpl)
        assert isinstance(get_engine('pyarrow'), PyArrowImpl)
        assert isinstance(get_engine('fastparquet'), FastParquetImpl)

    with pd.option_context('io.parquet.engine', 'auto'):
        assert isinstance(get_engine('auto'), PyArrowImpl)
        assert isinstance(get_engine('pyarrow'), PyArrowImpl)
        assert isinstance(get_engine('fastparquet'), FastParquetImpl)


@pytest.mark.xfail(reason="fp does not ignore pa index __index_level_0__")
def test_cross_engine_pa_fp(df_cross_compat, pa, fp):
    # cross-compat with differing reading/writing engines

    df = df_cross_compat
    with tm.ensure_clean() as path:
        df.to_parquet(path, engine=pa, compression=None)

        result = read_parquet(path, engine=fp, compression=None)
        tm.assert_frame_equal(result, df)


@pytest.mark.xfail(reason="pyarrow reading fp in some cases")
def test_cross_engine_fp_pa(df_cross_compat, pa, fp):
    # cross-compat with differing reading/writing engines

    df = df_cross_compat
    with tm.ensure_clean() as path:
        df.to_parquet(path, engine=fp, compression=None)

        result = read_parquet(path, engine=pa, compression=None)
        tm.assert_frame_equal(result, df)


class Base(object):

    def check_error_on_write(self, df, engine, exc):
        # check that we are raising the exception
        # on writing

        with pytest.raises(exc):
            with tm.ensure_clean() as path:
                to_parquet(df, path, engine, compression=None)

    def check_round_trip(self, df, engine, expected=None, **kwargs):

        with tm.ensure_clean() as path:
            df.to_parquet(path, engine, **kwargs)
            result = read_parquet(path, engine)

            if expected is None:
                expected = df
            tm.assert_frame_equal(result, expected)

            # repeat
            to_parquet(df, path, engine, **kwargs)
            result = pd.read_parquet(path, engine)

            if expected is None:
                expected = df
            tm.assert_frame_equal(result, expected)


class TestBasic(Base):

    def test_error(self, engine):

        for obj in [pd.Series([1, 2, 3]), 1, 'foo', pd.Timestamp('20130101'),
                    np.array([1, 2, 3])]:
            self.check_error_on_write(obj, engine, ValueError)

    def test_columns_dtypes(self, engine):

        df = pd.DataFrame({'string': list('abc'),
                           'int': list(range(1, 4))})

        # unicode
        df.columns = [u'foo', u'bar']
        self.check_round_trip(df, engine, compression=None)

    def test_columns_dtypes_invalid(self, engine):

        df = pd.DataFrame({'string': list('abc'),
                           'int': list(range(1, 4))})

        # numeric
        df.columns = [0, 1]
        self.check_error_on_write(df, engine, ValueError)

        if PY3:
            # bytes on PY3, on PY2 these are str
            df.columns = [b'foo', b'bar']
            self.check_error_on_write(df, engine, ValueError)

        # python object
        df.columns = [datetime.datetime(2011, 1, 1, 0, 0),
                      datetime.datetime(2011, 1, 1, 1, 1)]
        self.check_error_on_write(df, engine, ValueError)

    def test_write_with_index(self, engine):

        df = pd.DataFrame({'A': [1, 2, 3]})
        self.check_round_trip(df, engine, compression=None)

        # non-default index
        for index in [[2, 3, 4],
                      pd.date_range('20130101', periods=3),
                      list('abc'),
                      [1, 3, 4],
                      pd.MultiIndex.from_tuples([('a', 1), ('a', 2),
                                                 ('b', 1)]),
                      ]:

            df.index = index
            self.check_error_on_write(df, engine, ValueError)

        # index with meta-data
        df.index = [0, 1, 2]
        df.index.name = 'foo'
        self.check_error_on_write(df, engine, ValueError)

        # column multi-index
        df.index = [0, 1, 2]
        df.columns = pd.MultiIndex.from_tuples([('a', 1), ('a', 2), ('b', 1)]),
        self.check_error_on_write(df, engine, ValueError)

    @pytest.mark.parametrize('compression', [None, 'gzip', 'snappy', 'brotli'])
    def test_compression(self, engine, compression):

        if compression == 'snappy':
            pytest.importorskip('snappy')

        elif compression == 'brotli':
            pytest.importorskip('brotli')

        df = pd.DataFrame({'A': [1, 2, 3]})
        self.check_round_trip(df, engine, compression=compression)


class TestParquetPyArrow(Base):

    def test_basic(self, pa):

        df = pd.DataFrame({'string': list('abc'),
                           'string_with_nan': ['a', np.nan, 'c'],
                           'string_with_none': ['a', None, 'c'],
                           'bytes': [b'foo', b'bar', b'baz'],
                           'unicode': [u'foo', u'bar', u'baz'],
                           'int': list(range(1, 4)),
                           'uint': np.arange(3, 6).astype('u1'),
                           'float': np.arange(4.0, 7.0, dtype='float64'),
                           'float_with_nan': [2., np.nan, 3.],
                           'bool': [True, False, True],
                           'bool_with_none': [True, None, True],
                           'datetime_ns': pd.date_range('20130101', periods=3),
                           'datetime_with_nat': [pd.Timestamp('20130101'),
                                                 pd.NaT,
                                                 pd.Timestamp('20130103')]
                           })

        self.check_round_trip(df, pa)

    def test_duplicate_columns(self, pa):

        # not currently able to handle duplicate columns
        df = pd.DataFrame(np.arange(12).reshape(4, 3),
                          columns=list('aaa')).copy()
        self.check_error_on_write(df, pa, ValueError)

    def test_unsupported(self, pa):

        # period
        df = pd.DataFrame({'a': pd.period_range('2013', freq='M', periods=3)})
        self.check_error_on_write(df, pa, ValueError)

        # timedelta
        df = pd.DataFrame({'a': pd.timedelta_range('1 day',
                                                   periods=3)})
        self.check_error_on_write(df, pa, NotImplementedError)

        # mixed python objects
        df = pd.DataFrame({'a': ['a', 1, 2.0]})
        self.check_error_on_write(df, pa, ValueError)

    def test_categorical(self, pa_ge_070):
        pa = pa_ge_070

        # supported in >= 0.7.0
        df = pd.DataFrame({'a': pd.Categorical(list('abc'))})

        # de-serialized as object
        expected = df.assign(a=df.a.astype(object))
        self.check_round_trip(df, pa, expected)

    def test_categorical_unsupported(self, pa_lt_070):
        pa = pa_lt_070

        # supported in >= 0.7.0
        df = pd.DataFrame({'a': pd.Categorical(list('abc'))})
        self.check_error_on_write(df, pa, NotImplementedError)


class TestParquetFastParquet(Base):

    def test_basic(self, fp):

        df = pd.DataFrame(
            {'string': list('abc'),
             'string_with_nan': ['a', np.nan, 'c'],
             'string_with_none': ['a', None, 'c'],
             'bytes': [b'foo', b'bar', b'baz'],
             'unicode': [u'foo', u'bar', u'baz'],
             'int': list(range(1, 4)),
             'uint': np.arange(3, 6).astype('u1'),
             'float': np.arange(4.0, 7.0, dtype='float64'),
             'float_with_nan': [2., np.nan, 3.],
             'bool': [True, False, True],
             'datetime': pd.date_range('20130101', periods=3),
             'datetime_with_nat': [pd.Timestamp('20130101'),
                                   pd.NaT,
                                   pd.Timestamp('20130103')],
             'timedelta': pd.timedelta_range('1 day', periods=3),
             })

        self.check_round_trip(df, fp, compression=None)

    @pytest.mark.skip(reason="not supported")
    def test_duplicate_columns(self, fp):

        # not currently able to handle duplicate columns
        df = pd.DataFrame(np.arange(12).reshape(4, 3),
                          columns=list('aaa')).copy()
        self.check_error_on_write(df, fp, ValueError)

    def test_bool_with_none(self, fp):
        df = pd.DataFrame({'a': [True, None, False]})
        expected = pd.DataFrame({'a': [1.0, np.nan, 0.0]}, dtype='float16')
        self.check_round_trip(df, fp, expected=expected, compression=None)

    def test_unsupported(self, fp):

        # period
        df = pd.DataFrame({'a': pd.period_range('2013', freq='M', periods=3)})
        self.check_error_on_write(df, fp, ValueError)

        # mixed
        df = pd.DataFrame({'a': ['a', 1, 2.0]})
        self.check_error_on_write(df, fp, ValueError)

    def test_categorical(self, fp):
        if LooseVersion(fastparquet.__version__) < LooseVersion("0.1.3"):
            pytest.skip("CategoricalDtype not supported for older fp")
        df = pd.DataFrame({'a': pd.Categorical(list('abc'))})
        self.check_round_trip(df, fp, compression=None)

    def test_datetime_tz(self, fp):
        # doesn't preserve tz
        df = pd.DataFrame({'a': pd.date_range('20130101', periods=3,
                                              tz='US/Eastern')})

        # warns on the coercion
        with catch_warnings(record=True):
            self.check_round_trip(df, fp, df.astype('datetime64[ns]'),
                                  compression=None)
