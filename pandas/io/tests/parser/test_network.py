# -*- coding: utf-8 -*-

"""
Tests parsers ability to read and parse non-local files
and hence require a network connection to be read.
"""

import os
import nose

import pandas.util.testing as tm
from pandas import DataFrame
from pandas import compat
from pandas.io.parsers import read_csv, read_table


class TestUrlGz(tm.TestCase):

    def setUp(self):
        dirpath = tm.get_data_path()
        localtable = os.path.join(dirpath, 'salary.table.csv')
        self.local_table = read_table(localtable)

    @tm.network
    def test_url_gz(self):
        url = ('https://raw.github.com/pydata/pandas/'
               'master/pandas/io/tests/parser/data/salary.table.gz')
        url_table = read_table(url, compression="gzip", engine="python")
        tm.assert_frame_equal(url_table, self.local_table)

    @tm.network
    def test_url_gz_infer(self):
        url = 'https://s3.amazonaws.com/pandas-test/salary.table.gz'
        url_table = read_table(url, compression="infer", engine="python")
        tm.assert_frame_equal(url_table, self.local_table)


class TestS3(tm.TestCase):

    def setUp(self):
        try:
            import boto  # noqa
        except ImportError:
            raise nose.SkipTest("boto not installed")

    @tm.network
    def test_parse_public_s3_bucket(self):
        for ext, comp in [('', None), ('.gz', 'gzip'), ('.bz2', 'bz2')]:
            if comp == 'bz2' and compat.PY2:
                # The Python 2 C parser can't read bz2 from S3.
                self.assertRaises(ValueError, read_csv,
                                  's3://pandas-test/tips.csv' + ext,
                                  compression=comp)
            else:
                df = read_csv('s3://pandas-test/tips.csv' +
                              ext, compression=comp)
                self.assertTrue(isinstance(df, DataFrame))
                self.assertFalse(df.empty)
                tm.assert_frame_equal(read_csv(
                    tm.get_data_path('tips.csv')), df)

        # Read public file from bucket with not-public contents
        df = read_csv('s3://cant_get_it/tips.csv')
        self.assertTrue(isinstance(df, DataFrame))
        self.assertFalse(df.empty)
        tm.assert_frame_equal(read_csv(tm.get_data_path('tips.csv')), df)

    @tm.network
    def test_parse_public_s3n_bucket(self):
        # Read from AWS s3 as "s3n" URL
        df = read_csv('s3n://pandas-test/tips.csv', nrows=10)
        self.assertTrue(isinstance(df, DataFrame))
        self.assertFalse(df.empty)
        tm.assert_frame_equal(read_csv(
            tm.get_data_path('tips.csv')).iloc[:10], df)

    @tm.network
    def test_parse_public_s3a_bucket(self):
        # Read from AWS s3 as "s3a" URL
        df = read_csv('s3a://pandas-test/tips.csv', nrows=10)
        self.assertTrue(isinstance(df, DataFrame))
        self.assertFalse(df.empty)
        tm.assert_frame_equal(read_csv(
            tm.get_data_path('tips.csv')).iloc[:10], df)

    @tm.network
    def test_parse_public_s3_bucket_nrows(self):
        for ext, comp in [('', None), ('.gz', 'gzip'), ('.bz2', 'bz2')]:
            if comp == 'bz2' and compat.PY2:
                # The Python 2 C parser can't read bz2 from S3.
                self.assertRaises(ValueError, read_csv,
                                  's3://pandas-test/tips.csv' + ext,
                                  compression=comp)
            else:
                df = read_csv('s3://pandas-test/tips.csv' +
                              ext, nrows=10, compression=comp)
                self.assertTrue(isinstance(df, DataFrame))
                self.assertFalse(df.empty)
                tm.assert_frame_equal(read_csv(
                    tm.get_data_path('tips.csv')).iloc[:10], df)

    @tm.network
    def test_parse_public_s3_bucket_chunked(self):
        # Read with a chunksize
        chunksize = 5
        local_tips = read_csv(tm.get_data_path('tips.csv'))
        for ext, comp in [('', None), ('.gz', 'gzip'), ('.bz2', 'bz2')]:
            if comp == 'bz2' and compat.PY2:
                # The Python 2 C parser can't read bz2 from S3.
                self.assertRaises(ValueError, read_csv,
                                  's3://pandas-test/tips.csv' + ext,
                                  compression=comp)
            else:
                df_reader = read_csv('s3://pandas-test/tips.csv' + ext,
                                     chunksize=chunksize, compression=comp)
                self.assertEqual(df_reader.chunksize, chunksize)
                for i_chunk in [0, 1, 2]:
                    # Read a couple of chunks and make sure we see them
                    # properly.
                    df = df_reader.get_chunk()
                    self.assertTrue(isinstance(df, DataFrame))
                    self.assertFalse(df.empty)
                    true_df = local_tips.iloc[
                        chunksize * i_chunk: chunksize * (i_chunk + 1)]
                    # Chunking doesn't preserve row numbering
                    true_df = true_df.reset_index().drop('index', axis=1)
                    tm.assert_frame_equal(true_df, df)

    @tm.network
    def test_parse_public_s3_bucket_chunked_python(self):
        # Read with a chunksize using the Python parser
        chunksize = 5
        local_tips = read_csv(tm.get_data_path('tips.csv'))
        for ext, comp in [('', None), ('.gz', 'gzip'), ('.bz2', 'bz2')]:
            df_reader = read_csv('s3://pandas-test/tips.csv' + ext,
                                 chunksize=chunksize, compression=comp,
                                 engine='python')
            self.assertEqual(df_reader.chunksize, chunksize)
            for i_chunk in [0, 1, 2]:
                # Read a couple of chunks and make sure we see them properly.
                df = df_reader.get_chunk()
                self.assertTrue(isinstance(df, DataFrame))
                self.assertFalse(df.empty)
                true_df = local_tips.iloc[
                    chunksize * i_chunk: chunksize * (i_chunk + 1)]
                # Chunking doesn't preserve row numbering
                true_df = true_df.reset_index().drop('index', axis=1)
                tm.assert_frame_equal(true_df, df)

    @tm.network
    def test_parse_public_s3_bucket_python(self):
        for ext, comp in [('', None), ('.gz', 'gzip'), ('.bz2', 'bz2')]:
            df = read_csv('s3://pandas-test/tips.csv' + ext, engine='python',
                          compression=comp)
            self.assertTrue(isinstance(df, DataFrame))
            self.assertFalse(df.empty)
            tm.assert_frame_equal(read_csv(
                tm.get_data_path('tips.csv')), df)

    @tm.network
    def test_infer_s3_compression(self):
        for ext in ['', '.gz', '.bz2']:
            df = read_csv('s3://pandas-test/tips.csv' + ext,
                          engine='python', compression='infer')
            self.assertTrue(isinstance(df, DataFrame))
            self.assertFalse(df.empty)
            tm.assert_frame_equal(read_csv(
                tm.get_data_path('tips.csv')), df)

    @tm.network
    def test_parse_public_s3_bucket_nrows_python(self):
        for ext, comp in [('', None), ('.gz', 'gzip'), ('.bz2', 'bz2')]:
            df = read_csv('s3://pandas-test/tips.csv' + ext, engine='python',
                          nrows=10, compression=comp)
            self.assertTrue(isinstance(df, DataFrame))
            self.assertFalse(df.empty)
            tm.assert_frame_equal(read_csv(
                tm.get_data_path('tips.csv')).iloc[:10], df)

    @tm.network
    def test_s3_fails(self):
        import boto
        with tm.assertRaisesRegexp(boto.exception.S3ResponseError,
                                   'S3ResponseError: 404 Not Found'):
            read_csv('s3://nyqpug/asdf.csv')

        # Receive a permission error when trying to read a private bucket.
        # It's irrelevant here that this isn't actually a table.
        with tm.assertRaisesRegexp(boto.exception.S3ResponseError,
                                   'S3ResponseError: 403 Forbidden'):
            read_csv('s3://cant_get_it/')

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
