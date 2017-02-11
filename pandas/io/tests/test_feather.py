""" test feather-format compat """

import pytest
feather = pytest.importorskip('feather')

import numpy as np
import pandas as pd
from pandas.io.feather_format import to_feather, read_feather

from feather import FeatherError
import pandas.util.testing as tm
from pandas.util.testing import assert_frame_equal, ensure_clean


class TestFeather(tm.TestCase):

    def setUp(self):
        pass

    def check_error_on_write(self, df, exc):
        # check that we are raising the exception
        # on writing

        def f():
            with ensure_clean() as path:
                to_feather(df, path)
        self.assertRaises(exc, f)

    def check_round_trip(self, df):

        with ensure_clean() as path:
            to_feather(df, path)
            result = read_feather(path)
            assert_frame_equal(result, df)

    def test_error(self):

        for obj in [pd.Series([1, 2, 3]), 1, 'foo', pd.Timestamp('20130101'),
                    np.array([1, 2, 3])]:
            self.check_error_on_write(obj, ValueError)

    def test_basic(self):

        df = pd.DataFrame({'a': list('abc'),
                           'b': list(range(1, 4)),
                           'c': np.arange(3, 6).astype('u1'),
                           'd': np.arange(4.0, 7.0, dtype='float64'),
                           'e': [True, False, True],
                           'f': pd.Categorical(list('abc')),
                           'g': pd.date_range('20130101', periods=3),
                           'h': pd.date_range('20130101', periods=3,
                                              tz='US/Eastern'),
                           'i': pd.date_range('20130101', periods=3,
                                              freq='ns')})

        self.check_round_trip(df)

    def test_strided_data_issues(self):

        # strided data issuehttps://github.com/wesm/feather/issues/97
        df = pd.DataFrame(np.arange(12).reshape(4, 3), columns=list('abc'))
        self.check_error_on_write(df, FeatherError)

    def test_duplicate_columns(self):

        # https://github.com/wesm/feather/issues/53
        # not currently able to handle duplicate columns
        df = pd.DataFrame(np.arange(12).reshape(4, 3),
                          columns=list('aaa')).copy()
        self.check_error_on_write(df, ValueError)

    def test_stringify_columns(self):

        df = pd.DataFrame(np.arange(12).reshape(4, 3)).copy()
        self.check_error_on_write(df, ValueError)

    def test_unsupported(self):

        # period
        df = pd.DataFrame({'a': pd.period_range('2013', freq='M', periods=3)})
        self.check_error_on_write(df, ValueError)

        # non-strings
        df = pd.DataFrame({'a': ['a', 1, 2.0]})
        self.check_error_on_write(df, ValueError)

    def test_write_with_index(self):

        df = pd.DataFrame({'A': [1, 2, 3]})
        self.check_round_trip(df)

        # non-default index
        for index in [[2, 3, 4],
                      pd.date_range('20130101', periods=3),
                      list('abc'),
                      [1, 3, 4],
                      pd.MultiIndex.from_tuples([('a', 1), ('a', 2),
                                                 ('b', 1)]),
                      ]:

            df.index = index
            self.check_error_on_write(df, ValueError)

        # index with meta-data
        df.index = [0, 1, 2]
        df.index.name = 'foo'
        self.check_error_on_write(df, ValueError)

        # column multi-index
        df.index = [0, 1, 2]
        df.columns = pd.MultiIndex.from_tuples([('a', 1), ('a', 2), ('b', 1)]),
        self.check_error_on_write(df, ValueError)
