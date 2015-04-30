import ast
import datetime
import json
import nose
import os
import pytz
import shutil
import subprocess
import sys
import platform
from time import sleep

import numpy as np

from distutils.version import LooseVersion
from pandas import compat

from pandas import NaT
from pandas.compat import u
from pandas.core.frame import DataFrame
import pandas.io.gbq as gbq
import pandas.util.testing as tm

PROJECT_ID = None

VERSION = platform.python_version()

_IMPORTS = False
_GOOGLE_API_CLIENT_INSTALLED = False
_GOOGLE_API_CLIENT_VALID_VERSION = False
_HTTPLIB2_INSTALLED = False
_SETUPTOOLS_INSTALLED = False

def missing_bq():
    try:
        subprocess.call('bq')
        return False
    except OSError:
        return True

def _test_imports():
    if not compat.PY3:

        global _GOOGLE_API_CLIENT_INSTALLED, _GOOGLE_API_CLIENT_VALID_VERSION, \
               _HTTPLIB2_INSTALLED, _SETUPTOOLS_INSTALLED

        try:
            import pkg_resources
            _SETUPTOOLS_INSTALLED = True
        except ImportError:
            _SETUPTOOLS_INSTALLED = False

        if _SETUPTOOLS_INSTALLED:
            try:
                from apiclient.discovery import build
                from apiclient.errors import HttpError

                from oauth2client.client import OAuth2WebServerFlow
                from oauth2client.client import AccessTokenRefreshError

                from oauth2client.file import Storage
                from oauth2client.tools import run_flow
                _GOOGLE_API_CLIENT_INSTALLED=True
                _GOOGLE_API_CLIENT_VERSION = pkg_resources.get_distribution('google-api-python-client').version

                if LooseVersion(_GOOGLE_API_CLIENT_VERSION) >= '1.2.0':
                    _GOOGLE_API_CLIENT_VALID_VERSION = True

            except ImportError:
                _GOOGLE_API_CLIENT_INSTALLED = False


            try:
                import httplib2
                _HTTPLIB2_INSTALLED = True
            except ImportError:
                _HTTPLIB2_INSTALLED = False
    

    if compat.PY3:
        raise NotImplementedError("Google's libraries do not support Python 3 yet")

    if not _SETUPTOOLS_INSTALLED:
        raise ImportError('Could not import pkg_resources (setuptools).')

    if not _GOOGLE_API_CLIENT_INSTALLED:
        raise ImportError('Could not import Google API Client.')

    if not _GOOGLE_API_CLIENT_VALID_VERSION:
        raise ImportError("pandas requires google-api-python-client >= 1.2.0 for Google "
                          "BigQuery support, current version " + _GOOGLE_API_CLIENT_VERSION)

    if not _HTTPLIB2_INSTALLED:
        raise ImportError("pandas requires httplib2 for Google BigQuery support")

def test_requirements():
    try:
        _test_imports()
    except (ImportError, NotImplementedError) as import_exception:
        raise nose.SkipTest(import_exception)

class TestGBQConnectorIntegration(tm.TestCase):
    def setUp(self):
        test_requirements()

        if not PROJECT_ID:
            raise nose.SkipTest("Cannot run integration tests without a project id")
        
        self.sut = gbq.GbqConnector(PROJECT_ID)

    def test_should_be_able_to_make_a_connector(self):
        self.assertTrue(self.sut is not None, 'Could not create a GbqConnector')

    def test_should_be_able_to_get_valid_credentials(self):
        credentials = self.sut.get_credentials()
        self.assertFalse(credentials.invalid, 'Returned credentials invalid')

    def test_should_be_able_to_get_a_bigquery_service(self):
        credentials = self.sut.get_credentials()
        bigquery_service = self.sut.get_service(credentials)
        self.assertTrue(bigquery_service is not None, 'No service returned')

    def test_should_be_able_to_get_schema_from_query(self):
        schema, pages = self.sut.run_query('SELECT 1')
        self.assertTrue(schema is not None)

    def test_should_be_able_to_get_results_from_query(self):
        schema, pages = self.sut.run_query('SELECT 1')
        self.assertTrue(pages is not None)

class TestReadGBQUnitTests(tm.TestCase):
    def setUp(self):
        test_requirements()

    def test_should_return_bigquery_integers_as_python_floats(self):
        result = gbq._parse_entry(1, 'INTEGER')
        tm.assert_equal(result, float(1))

    def test_should_return_bigquery_floats_as_python_floats(self):
        result = gbq._parse_entry(1, 'FLOAT')
        tm.assert_equal(result, float(1))

    def test_should_return_bigquery_timestamps_as_numpy_datetime(self):
        result = gbq._parse_entry('0e9', 'TIMESTAMP')
        tm.assert_equal(result, np.datetime64('1970-01-01T00:00:00Z'))

    def test_should_return_bigquery_booleans_as_python_booleans(self):
        result = gbq._parse_entry('false', 'BOOLEAN')
        tm.assert_equal(result, False)

    def test_should_return_bigquery_strings_as_python_strings(self):
        result = gbq._parse_entry('STRING', 'STRING')
        tm.assert_equal(result, 'STRING')

    def test_to_gbq_should_fail_if_invalid_table_name_passed(self):
        with tm.assertRaises(gbq.NotFoundException):
            gbq.to_gbq(DataFrame(), 'invalid_table_name', project_id="1234")

    def test_to_gbq_with_no_project_id_given_should_fail(self):
        with tm.assertRaises(TypeError):
            gbq.to_gbq(DataFrame(), 'dataset.tablename')

    def test_read_gbq_with_no_project_id_given_should_fail(self):
        with tm.assertRaises(TypeError):
            gbq.read_gbq('SELECT "1" as NUMBER_1')

    def test_that_parse_data_works_properly(self):
        test_schema = {'fields': [{'mode': 'NULLABLE',
                  'name': 'VALID_STRING',
                  'type': 'STRING'}]}
        test_page = [{'f': [{'v': 'PI'}]}]

        test_output = gbq._parse_data(test_schema, test_page)
        correct_output = DataFrame({'VALID_STRING' : ['PI']})
        tm.assert_frame_equal(test_output, correct_output)

class TestReadGBQIntegration(tm.TestCase):
    def setUp(self):
        test_requirements()

        if not PROJECT_ID:
            raise nose.SkipTest("Cannot run integration tests without a project id")

    def test_should_properly_handle_valid_strings(self):
        query = 'SELECT "PI" as VALID_STRING'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'VALID_STRING' : ['PI']}))

    def test_should_properly_handle_empty_strings(self):
        query = 'SELECT "" as EMPTY_STRING'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'EMPTY_STRING' : [""]}))

    def test_should_properly_handle_null_strings(self):
        query = 'SELECT STRING(NULL) as NULL_STRING'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'NULL_STRING' : [None]}))

    def test_should_properly_handle_valid_integers(self):
        query = 'SELECT INTEGER(3) as VALID_INTEGER'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'VALID_INTEGER' : [3]}))

    def test_should_properly_handle_null_integers(self):
        query = 'SELECT INTEGER(NULL) as NULL_INTEGER'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'NULL_INTEGER' : [np.nan]}))

    def test_should_properly_handle_valid_floats(self):
        query = 'SELECT PI() as VALID_FLOAT'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'VALID_FLOAT' : [3.141592653589793]}))

    def test_should_properly_handle_null_floats(self):
        query = 'SELECT FLOAT(NULL) as NULL_FLOAT'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'NULL_FLOAT' : [np.nan]}))

    def test_should_properly_handle_timestamp_unix_epoch(self):
        query = 'SELECT TIMESTAMP("1970-01-01 00:00:00") as UNIX_EPOCH'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'UNIX_EPOCH' : [np.datetime64('1970-01-01T00:00:00.000000Z')]}))

    def test_should_properly_handle_arbitrary_timestamp(self):
        query = 'SELECT TIMESTAMP("2004-09-15 05:00:00") as VALID_TIMESTAMP'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'VALID_TIMESTAMP' : [np.datetime64('2004-09-15T05:00:00.000000Z')]}))

    def test_should_properly_handle_null_timestamp(self):
        query = 'SELECT TIMESTAMP(NULL) as NULL_TIMESTAMP'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'NULL_TIMESTAMP' :[NaT]}))

    def test_should_properly_handle_true_boolean(self):
        query = 'SELECT BOOLEAN(TRUE) as TRUE_BOOLEAN'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'TRUE_BOOLEAN' : [True]}))

    def test_should_properly_handle_false_boolean(self):
        query = 'SELECT BOOLEAN(FALSE) as FALSE_BOOLEAN'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'FALSE_BOOLEAN' : [False]}))

    def test_should_properly_handle_null_boolean(self):
        query = 'SELECT BOOLEAN(NULL) as NULL_BOOLEAN'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'NULL_BOOLEAN' : [None]}))

    def test_unicode_string_conversion_and_normalization(self):
        correct_test_datatype = DataFrame(
            {'UNICODE_STRING' : [u("\xe9\xfc")]}
        )

        query = 'SELECT "\xc3\xa9\xc3\xbc" as UNICODE_STRING'

        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, correct_test_datatype)

    def test_index_column(self):
        query = "SELECT 'a' as STRING_1, 'b' as STRING_2"
        result_frame = gbq.read_gbq(query, project_id=PROJECT_ID, index_col="STRING_1")
        correct_frame = DataFrame({'STRING_1' : ['a'], 'STRING_2' : ['b']}).set_index("STRING_1")
        tm.assert_equal(result_frame.index.name, correct_frame.index.name)

    def test_column_order(self):
        query = "SELECT 'a' as STRING_1, 'b' as STRING_2, 'c' as STRING_3"
        col_order = ['STRING_3', 'STRING_1', 'STRING_2']
        result_frame = gbq.read_gbq(query, project_id=PROJECT_ID, col_order=col_order)
        correct_frame = DataFrame({'STRING_1' : ['a'], 'STRING_2' : ['b'], 'STRING_3' : ['c']})[col_order]
        tm.assert_frame_equal(result_frame, correct_frame)

    def test_column_order_plus_index(self):
        query = "SELECT 'a' as STRING_1, 'b' as STRING_2, 'c' as STRING_3"
        col_order = ['STRING_3', 'STRING_2']
        result_frame = gbq.read_gbq(query, project_id=PROJECT_ID, index_col='STRING_1', col_order=col_order)
        correct_frame = DataFrame({'STRING_1' : ['a'], 'STRING_2' : ['b'], 'STRING_3' : ['c']})
        correct_frame.set_index('STRING_1', inplace=True)
        correct_frame = correct_frame[col_order]
        tm.assert_frame_equal(result_frame, correct_frame)

    def test_malformed_query(self):
        with tm.assertRaises(gbq.InvalidQueryException):
            gbq.read_gbq("SELCET * FORM [publicdata:samples.shakespeare]", project_id=PROJECT_ID)

    def test_bad_project_id(self):
        with tm.assertRaises(gbq.NotFoundException):
            gbq.read_gbq("SELECT 1", project_id='001')

    def test_bad_table_name(self):
        with tm.assertRaises(gbq.NotFoundException):
            gbq.read_gbq("SELECT * FROM [publicdata:samples.nope]", project_id=PROJECT_ID)

    def test_download_dataset_larger_than_200k_rows(self):
        # Test for known BigQuery bug in datasets larger than 100k rows
        # http://stackoverflow.com/questions/19145587/bq-py-not-paging-results
        df = gbq.read_gbq("SELECT id FROM [publicdata:samples.wikipedia] GROUP EACH BY id ORDER BY id ASC LIMIT 200005", project_id=PROJECT_ID)
        self.assertEqual(len(df.drop_duplicates()), 200005)

class TestToGBQIntegration(tm.TestCase):
    # This class requires bq.py to be installed for setup/teardown. 
    # It will also need to be preconfigured with a default dataset,
    # so, be sure to `bq init` in terminal before running.

    def setUp(self):
        test_requirements()

        if not PROJECT_ID:
            raise nose.SkipTest("Cannot run integration tests without a project id")
        if missing_bq():
            raise nose.SkipTest("Cannot run to_gbq tests without bq command line client")

    @classmethod
    def setUpClass(cls):
        if PROJECT_ID and not missing_bq():
            subprocess.call(['bq','mk','pydata_pandas_bq_testing'])
            subprocess.call(['bq','mk','pydata_pandas_bq_testing.new_test','bools:BOOLEAN,flts:FLOAT,ints:INTEGER,strs:STRING,times:TIMESTAMP'])

    def test_upload_data(self):
        test_size = 1000001
        #create df to test for all BQ datatypes except RECORD
        bools = np.random.randint(2, size=(1,test_size)).astype(bool)
        flts = np.random.randn(1,test_size)
        ints = np.random.randint(1,10, size=(1,test_size))
        strs = np.random.randint(1,10, size=(1,test_size)).astype(str)
        times = [datetime.datetime.now(pytz.timezone('US/Arizona')) for t in xrange(test_size)]
        df = DataFrame({'bools':bools[0], 'flts':flts[0], 'ints':ints[0], 'strs':strs[0], 'times':times[0]}, index=range(test_size))
        gbq.to_gbq(df,"pydata_pandas_bq_testing.new_test", project_id=PROJECT_ID, chunksize=10000)
        sleep(60) # <- Curses Google!!!

        result = gbq.read_gbq("SELECT COUNT(*) as NUM_ROWS FROM pydata_pandas_bq_testing.new_test", project_id=PROJECT_ID)
        self.assertEqual(result['NUM_ROWS'][0], test_size)

    def test_google_upload_errors_should_raise_exception(self):
        test_timestamp = datetime.datetime.now(pytz.timezone('US/Arizona'))
        bad_df = DataFrame( {'bools': [False, False],
                             'flts': [0.0,1.0],
                             'ints': [0,'1'],
                             'strs': ['a', 1],
                             'times': [test_timestamp, test_timestamp]
                             }, index=range(2))
        with tm.assertRaises(gbq.UnknownGBQException):
            gbq.to_gbq(bad_df, 'pydata_pandas_bq_testing.new_test', project_id = PROJECT_ID)

    def test_generate_bq_schema(self):

        df = tm.makeMixedDataFrame()
        schema = gbq.generate_bq_schema(df)

        test_schema = {'fields': [{'name': 'A', 'type': 'FLOAT'},
                                  {'name': 'B', 'type': 'FLOAT'},
                                  {'name': 'C', 'type': 'STRING'},
                                  {'name': 'D', 'type': 'TIMESTAMP'}]}

        self.assertEqual(schema, test_schema)

    @classmethod
    def tearDownClass(cls):
        if PROJECT_ID and not missing_bq():
            subprocess.call(['bq','rm','-f','pydata_pandas_bq_testing.new_test'])
            subprocess.call(['bq','rm','-f','pydata_pandas_bq_testing'])

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)

