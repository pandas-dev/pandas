from datetime import datetime
import nose
import pytz
import platform
from time import sleep

import numpy as np

from distutils.version import StrictVersion
from pandas import compat

from pandas import NaT
from pandas.compat import u, range
from pandas.core.frame import DataFrame
import pandas.io.gbq as gbq
import pandas.util.testing as tm

PROJECT_ID = None
DATASET_ID = 'pydata_pandas_bq_testing'
TABLE_ID = 'new_test'
DESTINATION_TABLE = "{0}.{1}".format(DATASET_ID + "1", TABLE_ID)

VERSION = platform.python_version()

_IMPORTS = False
_GOOGLE_API_CLIENT_INSTALLED = False
_GOOGLE_API_CLIENT_VALID_VERSION = False
_HTTPLIB2_INSTALLED = False
_SETUPTOOLS_INSTALLED = False


def _test_imports():
    global _GOOGLE_API_CLIENT_INSTALLED, _GOOGLE_API_CLIENT_VALID_VERSION, \
           _HTTPLIB2_INSTALLED, _SETUPTOOLS_INSTALLED

    try:
        import pkg_resources
        _SETUPTOOLS_INSTALLED = True
    except ImportError:
        _SETUPTOOLS_INSTALLED = False

    if compat.PY3:
        google_api_minimum_version = '1.4.1'
    else:
        google_api_minimum_version = '1.2.0'

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

            if StrictVersion(_GOOGLE_API_CLIENT_VERSION) >= StrictVersion(google_api_minimum_version):
                _GOOGLE_API_CLIENT_VALID_VERSION = True

        except ImportError:
            _GOOGLE_API_CLIENT_INSTALLED = False

        try:
            import httplib2
            _HTTPLIB2_INSTALLED = True
        except ImportError:
            _HTTPLIB2_INSTALLED = False

    if not _SETUPTOOLS_INSTALLED:
        raise ImportError('Could not import pkg_resources (setuptools).')

    if not _GOOGLE_API_CLIENT_INSTALLED:
        raise ImportError('Could not import Google API Client.')

    if not _GOOGLE_API_CLIENT_VALID_VERSION:
        raise ImportError("pandas requires google-api-python-client >= {0} for Google BigQuery support, "
                          "current version {1}".format(google_api_minimum_version, _GOOGLE_API_CLIENT_VERSION))

    if not _HTTPLIB2_INSTALLED:
        raise ImportError("pandas requires httplib2 for Google BigQuery support")


def test_requirements():
    try:
        _test_imports()
    except (ImportError, NotImplementedError) as import_exception:
        raise nose.SkipTest(import_exception)


def clean_gbq_environment():
    dataset = gbq._Dataset(PROJECT_ID)

    for i in range(1, 10):
        if DATASET_ID + str(i) in dataset.datasets():
            dataset_id = DATASET_ID + str(i)
            table = gbq._Table(PROJECT_ID, dataset_id)
            for j in range(1, 20):
                if TABLE_ID + str(j) in dataset.tables(dataset_id):
                    table.delete(TABLE_ID + str(j))

            dataset.delete(dataset_id)


def make_mixed_dataframe_v2(test_size):
    # create df to test for all BQ datatypes except RECORD
    bools = np.random.randint(2, size=(1, test_size)).astype(bool)
    flts = np.random.randn(1, test_size)
    ints = np.random.randint(1, 10, size=(1, test_size))
    strs = np.random.randint(1, 10, size=(1, test_size)).astype(str)
    times = [datetime.now(pytz.timezone('US/Arizona')) for t in range(test_size)]
    return DataFrame({'bools': bools[0],
                      'flts': flts[0],
                      'ints': ints[0],
                      'strs': strs[0],
                      'times': times[0]},
                      index=range(test_size))


def test_generate_bq_schema_deprecated():
    # 11121 Deprecation of generate_bq_schema
    with tm.assert_produces_warning(FutureWarning):
        df = make_mixed_dataframe_v2(10)
        gbq.generate_bq_schema(df)

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
        test_schema = {'fields': [{'mode': 'NULLABLE', 'name': 'VALID_STRING', 'type': 'STRING'}]}
        test_page = [{'f': [{'v': 'PI'}]}]

        test_output = gbq._parse_data(test_schema, test_page)
        correct_output = DataFrame({'VALID_STRING': ['PI']})
        tm.assert_frame_equal(test_output, correct_output)


class TestReadGBQIntegration(tm.TestCase):
    @classmethod
    def setUpClass(cls):
        # - GLOBAL CLASS FIXTURES -
        #   put here any instruction you want to execute only *ONCE* *BEFORE* executing *ALL* tests
        #   described below.

        if not PROJECT_ID:
            raise nose.SkipTest("Cannot run integration tests without a project id")

        test_requirements()

    def setUp(self):
        # - PER-TEST FIXTURES -
        #   put here any instruction you want to be run *BEFORE* *EVERY* test is executed.
        pass

    @classmethod
    def tearDownClass(cls):
        # - GLOBAL CLASS FIXTURES -
        #   put here any instruction you want to execute only *ONCE* *AFTER* executing all tests.
        pass

    def tearDown(self):
        # - PER-TEST FIXTURES -
        #   put here any instructions you want to be run *AFTER* *EVERY* test is executed.
        pass

    def test_should_properly_handle_valid_strings(self):
        query = 'SELECT "PI" as VALID_STRING'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'VALID_STRING': ['PI']}))

    def test_should_properly_handle_empty_strings(self):
        query = 'SELECT "" as EMPTY_STRING'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'EMPTY_STRING': [""]}))

    def test_should_properly_handle_null_strings(self):
        query = 'SELECT STRING(NULL) as NULL_STRING'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'NULL_STRING': [None]}))

    def test_should_properly_handle_valid_integers(self):
        query = 'SELECT INTEGER(3) as VALID_INTEGER'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'VALID_INTEGER': [3]}))

    def test_should_properly_handle_null_integers(self):
        query = 'SELECT INTEGER(NULL) as NULL_INTEGER'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'NULL_INTEGER': [np.nan]}))

    def test_should_properly_handle_valid_floats(self):
        query = 'SELECT PI() as VALID_FLOAT'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'VALID_FLOAT': [3.141592653589793]}))

    def test_should_properly_handle_null_floats(self):
        query = 'SELECT FLOAT(NULL) as NULL_FLOAT'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'NULL_FLOAT': [np.nan]}))

    def test_should_properly_handle_timestamp_unix_epoch(self):
        query = 'SELECT TIMESTAMP("1970-01-01 00:00:00") as UNIX_EPOCH'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'UNIX_EPOCH': [np.datetime64('1970-01-01T00:00:00.000000Z')]}))

    def test_should_properly_handle_arbitrary_timestamp(self):
        query = 'SELECT TIMESTAMP("2004-09-15 05:00:00") as VALID_TIMESTAMP'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'VALID_TIMESTAMP': [np.datetime64('2004-09-15T05:00:00.000000Z')]}))

    def test_should_properly_handle_null_timestamp(self):
        query = 'SELECT TIMESTAMP(NULL) as NULL_TIMESTAMP'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'NULL_TIMESTAMP': [NaT]}))

    def test_should_properly_handle_true_boolean(self):
        query = 'SELECT BOOLEAN(TRUE) as TRUE_BOOLEAN'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'TRUE_BOOLEAN': [True]}))

    def test_should_properly_handle_false_boolean(self):
        query = 'SELECT BOOLEAN(FALSE) as FALSE_BOOLEAN'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'FALSE_BOOLEAN': [False]}))

    def test_should_properly_handle_null_boolean(self):
        query = 'SELECT BOOLEAN(NULL) as NULL_BOOLEAN'
        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, DataFrame({'NULL_BOOLEAN': [None]}))

    def test_unicode_string_conversion_and_normalization(self):
        correct_test_datatype = DataFrame(
            {'UNICODE_STRING': [u("\xe9\xfc")]}
        )

        unicode_string = "\xc3\xa9\xc3\xbc"

        if compat.PY3:
            unicode_string = unicode_string.encode('latin-1').decode('utf8')

        query = 'SELECT "{0}" as UNICODE_STRING'.format(unicode_string)

        df = gbq.read_gbq(query, project_id=PROJECT_ID)
        tm.assert_frame_equal(df, correct_test_datatype)

    def test_index_column(self):
        query = "SELECT 'a' as STRING_1, 'b' as STRING_2"
        result_frame = gbq.read_gbq(query, project_id=PROJECT_ID, index_col="STRING_1")
        correct_frame = DataFrame({'STRING_1': ['a'], 'STRING_2': ['b']}).set_index("STRING_1")
        tm.assert_equal(result_frame.index.name, correct_frame.index.name)

    def test_column_order(self):
        query = "SELECT 'a' as STRING_1, 'b' as STRING_2, 'c' as STRING_3"
        col_order = ['STRING_3', 'STRING_1', 'STRING_2']
        result_frame = gbq.read_gbq(query, project_id=PROJECT_ID, col_order=col_order)
        correct_frame = DataFrame({'STRING_1': ['a'], 'STRING_2': ['b'], 'STRING_3': ['c']})[col_order]
        tm.assert_frame_equal(result_frame, correct_frame)

    def test_column_order_plus_index(self):
        query = "SELECT 'a' as STRING_1, 'b' as STRING_2, 'c' as STRING_3"
        col_order = ['STRING_3', 'STRING_2']
        result_frame = gbq.read_gbq(query, project_id=PROJECT_ID, index_col='STRING_1', col_order=col_order)
        correct_frame = DataFrame({'STRING_1': ['a'], 'STRING_2': ['b'], 'STRING_3': ['c']})
        correct_frame.set_index('STRING_1', inplace=True)
        correct_frame = correct_frame[col_order]
        tm.assert_frame_equal(result_frame, correct_frame)

    def test_malformed_query(self):
        with tm.assertRaises(gbq.GenericGBQException):
            gbq.read_gbq("SELCET * FORM [publicdata:samples.shakespeare]", project_id=PROJECT_ID)

    def test_bad_project_id(self):
        with tm.assertRaises(gbq.GenericGBQException):
            gbq.read_gbq("SELECT 1", project_id='001')

    def test_bad_table_name(self):
        with tm.assertRaises(gbq.GenericGBQException):
            gbq.read_gbq("SELECT * FROM [publicdata:samples.nope]", project_id=PROJECT_ID)

    def test_download_dataset_larger_than_200k_rows(self):
        test_size = 200005
        # Test for known BigQuery bug in datasets larger than 100k rows
        # http://stackoverflow.com/questions/19145587/bq-py-not-paging-results
        df = gbq.read_gbq("SELECT id FROM [publicdata:samples.wikipedia] GROUP EACH BY id ORDER BY id ASC LIMIT {0}".format(test_size),
                          project_id=PROJECT_ID)
        self.assertEqual(len(df.drop_duplicates()), test_size)

    def test_zero_rows(self):
        # Bug fix for https://github.com/pydata/pandas/issues/10273
        df = gbq.read_gbq("SELECT title, language  FROM [publicdata:samples.wikipedia] where timestamp=-9999999",
                          project_id=PROJECT_ID)
        expected_result = DataFrame(columns=['title', 'language'])
        self.assert_frame_equal(df, expected_result)


class TestToGBQIntegration(tm.TestCase):
    # Changes to BigQuery table schema may take up to 2 minutes as of May 2015
    # As a workaround to this issue, each test should use a unique table name.
    # Make sure to modify the for loop range in the tearDownClass when a new test is added
    # See `Issue 191 <https://code.google.com/p/google-bigquery/issues/detail?id=191>`__

    @classmethod
    def setUpClass(cls):
        # - GLOBAL CLASS FIXTURES -
        # put here any instruction you want to execute only *ONCE* *BEFORE* executing *ALL* tests
        # described below.

        if not PROJECT_ID:
            raise nose.SkipTest("Cannot run integration tests without a project id")

        test_requirements()
        clean_gbq_environment()

        gbq._Dataset(PROJECT_ID).create(DATASET_ID + "1")

    def setUp(self):
        # - PER-TEST FIXTURES -
        # put here any instruction you want to be run *BEFORE* *EVERY* test is executed.

        self.dataset = gbq._Dataset(PROJECT_ID)
        self.table = gbq._Table(PROJECT_ID, DATASET_ID + "1")

    @classmethod
    def tearDownClass(cls):
        # - GLOBAL CLASS FIXTURES -
        #   put here any instruction you want to execute only *ONCE* *AFTER* executing all tests.

        clean_gbq_environment()

    def tearDown(self):
        # - PER-TEST FIXTURES -
        # put here any instructions you want to be run *AFTER* *EVERY* test is executed.
        pass

    def test_upload_data(self):
        destination_table = DESTINATION_TABLE + "1"

        test_size = 1000001
        df = make_mixed_dataframe_v2(test_size)

        gbq.to_gbq(df, destination_table, PROJECT_ID, chunksize=10000)

        sleep(60)  # <- Curses Google!!!

        result = gbq.read_gbq("SELECT COUNT(*) as NUM_ROWS FROM {0}".format(destination_table),
                              project_id=PROJECT_ID)
        self.assertEqual(result['NUM_ROWS'][0], test_size)

    def test_upload_data_if_table_exists_fail(self):
        destination_table = DESTINATION_TABLE + "2"

        test_size = 10
        df = make_mixed_dataframe_v2(test_size)
        self.table.create(TABLE_ID + "2", gbq._generate_bq_schema(df))

        # Test the default value of if_exists is 'fail'
        with tm.assertRaises(gbq.TableCreationError):
            gbq.to_gbq(df, destination_table, PROJECT_ID)

        # Test the if_exists parameter with value 'fail'
        with tm.assertRaises(gbq.TableCreationError):
            gbq.to_gbq(df, destination_table, PROJECT_ID, if_exists='fail')

    def test_upload_data_if_table_exists_append(self):
        destination_table = DESTINATION_TABLE + "3"

        test_size = 10
        df = make_mixed_dataframe_v2(test_size)
        df_different_schema = tm.makeMixedDataFrame()

        # Initialize table with sample data
        gbq.to_gbq(df, destination_table, PROJECT_ID, chunksize=10000)

        # Test the if_exists parameter with value 'append'
        gbq.to_gbq(df, destination_table, PROJECT_ID, if_exists='append')

        sleep(60)  # <- Curses Google!!!

        result = gbq.read_gbq("SELECT COUNT(*) as NUM_ROWS FROM {0}".format(destination_table), project_id=PROJECT_ID)
        self.assertEqual(result['NUM_ROWS'][0], test_size * 2)

        # Try inserting with a different schema, confirm failure
        with tm.assertRaises(gbq.InvalidSchema):
            gbq.to_gbq(df_different_schema, destination_table, PROJECT_ID, if_exists='append')

    def test_upload_data_if_table_exists_replace(self):
        destination_table = DESTINATION_TABLE + "4"

        test_size = 10
        df = make_mixed_dataframe_v2(test_size)
        df_different_schema = tm.makeMixedDataFrame()

        # Initialize table with sample data
        gbq.to_gbq(df, destination_table, PROJECT_ID, chunksize=10000)

        # Test the if_exists parameter with the value 'replace'.
        gbq.to_gbq(df_different_schema, destination_table, PROJECT_ID, if_exists='replace')

        sleep(60)  # <- Curses Google!!!

        result = gbq.read_gbq("SELECT COUNT(*) as NUM_ROWS FROM {0}".format(destination_table), project_id=PROJECT_ID)
        self.assertEqual(result['NUM_ROWS'][0], 5)

    def test_google_upload_errors_should_raise_exception(self):
        destination_table = DESTINATION_TABLE + "5"

        test_timestamp = datetime.now(pytz.timezone('US/Arizona'))
        bad_df = DataFrame({'bools': [False, False], 'flts': [0.0, 1.0], 'ints': [0, '1'], 'strs': ['a', 1],
                            'times': [test_timestamp, test_timestamp]}, index=range(2))

        with tm.assertRaises(gbq.StreamingInsertError):
            gbq.to_gbq(bad_df, destination_table, PROJECT_ID, verbose=True)

    def test_generate_schema(self):
        df = tm.makeMixedDataFrame()
        schema = gbq._generate_bq_schema(df)

        test_schema = {'fields': [{'name': 'A', 'type': 'FLOAT'},
                                  {'name': 'B', 'type': 'FLOAT'},
                                  {'name': 'C', 'type': 'STRING'},
                                  {'name': 'D', 'type': 'TIMESTAMP'}]}

        self.assertEqual(schema, test_schema)

    def test_create_table(self):
        destination_table = TABLE_ID + "6"
        test_schema = {'fields': [{'name': 'A', 'type': 'FLOAT'}, {'name': 'B', 'type': 'FLOAT'},
                                  {'name': 'C', 'type': 'STRING'}, {'name': 'D', 'type': 'TIMESTAMP'}]}
        self.table.create(destination_table, test_schema)
        self.assertTrue(self.table.exists(destination_table), 'Expected table to exist')

    def test_table_does_not_exist(self):
        self.assertTrue(not self.table.exists(TABLE_ID + "7"), 'Expected table not to exist')

    def test_delete_table(self):
        destination_table = TABLE_ID + "8"
        test_schema = {'fields': [{'name': 'A', 'type': 'FLOAT'}, {'name': 'B', 'type': 'FLOAT'},
                                  {'name': 'C', 'type': 'STRING'}, {'name': 'D', 'type': 'TIMESTAMP'}]}
        self.table.create(destination_table, test_schema)
        self.table.delete(destination_table)
        self.assertTrue(not self.table.exists(destination_table), 'Expected table not to exist')

    def test_list_table(self):
        destination_table = TABLE_ID + "9"
        test_schema = {'fields': [{'name': 'A', 'type': 'FLOAT'}, {'name': 'B', 'type': 'FLOAT'},
                                  {'name': 'C', 'type': 'STRING'}, {'name': 'D', 'type': 'TIMESTAMP'}]}
        self.table.create(destination_table, test_schema)
        self.assertTrue(destination_table in self.dataset.tables(DATASET_ID + "1"),
                        'Expected table list to contain table {0}'.format(destination_table))

    def test_list_dataset(self):
        dataset_id = DATASET_ID + "1"
        self.assertTrue(dataset_id in self.dataset.datasets(),
                        'Expected dataset list to contain dataset {0}'.format(dataset_id))

    def test_list_table_zero_results(self):
        dataset_id = DATASET_ID + "2"
        self.dataset.create(dataset_id)
        table_list = gbq._Dataset(PROJECT_ID).tables(dataset_id)
        self.assertEqual(len(table_list), 0, 'Expected gbq.list_table() to return 0')

    def test_create_dataset(self):
        dataset_id = DATASET_ID + "3"
        self.dataset.create(dataset_id)
        self.assertTrue(dataset_id in self.dataset.datasets(), 'Expected dataset to exist')

    def test_delete_dataset(self):
        dataset_id = DATASET_ID + "4"
        self.dataset.create(dataset_id)
        self.dataset.delete(dataset_id)
        self.assertTrue(dataset_id not in self.dataset.datasets(), 'Expected dataset not to exist')

    def test_dataset_exists(self):
        dataset_id = DATASET_ID + "5"
        self.dataset.create(dataset_id)
        self.assertTrue(self.dataset.exists(dataset_id), 'Expected dataset to exist')

    def create_table_data_dataset_does_not_exist(self):
        dataset_id = DATASET_ID + "6"
        table_id = TABLE_ID + "1"
        table_with_new_dataset = gbq._Table(PROJECT_ID, dataset_id)
        df = make_mixed_dataframe_v2(10)
        table_with_new_dataset.create(table_id, gbq._generate_bq_schema(df))
        self.assertTrue(self.dataset.exists(dataset_id), 'Expected dataset to exist')
        self.assertTrue(table_with_new_dataset.exists(table_id), 'Expected dataset to exist')

    def test_dataset_does_not_exist(self):
        self.assertTrue(not self.dataset.exists(DATASET_ID + "_not_found"), 'Expected dataset not to exist')

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
