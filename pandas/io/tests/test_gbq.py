import ast
import nose
import os
import shutil
import subprocess

import numpy as np

import pandas.io.gbq as gbq
import pandas.util.testing as tm

from pandas.core.frame import DataFrame
from pandas.util.testing import with_connectivity_check
from pandas import NaT


try:
    import bq
    import bigquery_client
    import gflags as flags
except ImportError:
    raise nose.SkipTest

####################################################################################
# Fake Google BigQuery Client

class FakeClient:
    def __init__(self):
        self.apiclient = FakeApiClient()
    def GetTableSchema(self,table_dict):
        retval =  {'fields': [
                    {'type': 'STRING', 'name': 'corpus', 'mode': 'NULLABLE'},
                    {'type': 'INTEGER', 'name': 'corpus_date', 'mode': 'NULLABLE'},
                    {'type': 'STRING', 'name': 'word', 'mode': 'NULLABLE'},
                    {'type': 'INTEGER', 'name': 'word_count', 'mode': 'NULLABLE'}
                  ]}
        return retval

# Fake Google BigQuery API Client
class FakeApiClient:
    def __init__(self):
        self._fakejobs = FakeJobs()


    def jobs(self):
        return self._fakejobs

class FakeJobs:
    def __init__(self):
        self._fakequeryresults = FakeResults()

    def getQueryResults(self, job_id=None, project_id=None,
                        max_results=None, timeout_ms=None, **kwargs):
        return self._fakequeryresults

class FakeResults:
    def execute(self):
        return {'rows': [  {'f': [{'v': 'othello'}, {'v': '1603'}, {'v': 'brave'}, {'v': '3'}]},
                            {'f': [{'v': 'othello'}, {'v': '1603'}, {'v': 'attended'}, {'v': '1'}]},
                            {'f': [{'v': 'othello'}, {'v': '1603'}, {'v': 'treason'}, {'v': '1'}]},
                            {'f': [{'v': 'othello'}, {'v': '1603'}, {'v': 'islanders'}, {'v': '1'}]},
                            {'f': [{'v': 'othello'}, {'v': '1603'}, {'v': 'heed'}, {'v': '3'}]},
                            {'f': [{'v': 'othello'}, {'v': '1603'}, {'v': 'alehouse'}, {'v': '1'}]},
                            {'f': [{'v': 'othello'}, {'v': '1603'}, {'v': 'corrigible'}, {'v': '1'}]},
                            {'f': [{'v': 'othello'}, {'v': '1603'}, {'v': 'brawl'}, {'v': '2'}]},
                            {'f': [{'v': 'othello'}, {'v': '1603'}, {'v': "'"}, {'v': '17'}]},
                            {'f': [{'v': 'othello'}, {'v': '1603'}, {'v': 'troubled'}, {'v': '1'}]}
                        ],
                'kind': 'bigquery#tableDataList',
                'etag': '"4PTsVxg68bQkQs1RJ1Ndewqkgg4/hoRHzb4qfhJAIa2mEewC-jhs9Bg"',
                'totalRows': '10',
                'jobComplete' : True}

####################################################################################

class TestGbq(tm.TestCase):
    def setUp(self):
        with open(self.fake_job_path, 'r') as fin:
            self.fake_job = ast.literal_eval(fin.read())

        self.test_data_small = [{'f': [{'v': 'othello'}, {'v': '1603'}, {'v': 'brave'}, {'v': '3'}]},
                {'f': [{'v': 'othello'}, {'v': '1603'}, {'v': 'attended'}, {'v': '1'}]},
                {'f': [{'v': 'othello'}, {'v': '1603'}, {'v': 'treason'}, {'v': '1'}]},
                {'f': [{'v': 'othello'}, {'v': '1603'}, {'v': 'islanders'}, {'v': '1'}]},
                {'f': [{'v': 'othello'}, {'v': '1603'}, {'v': 'heed'}, {'v': '3'}]},
                {'f': [{'v': 'othello'}, {'v': '1603'}, {'v': 'alehouse'}, {'v': '1'}]},
                {'f': [{'v': 'othello'}, {'v': '1603'}, {'v': 'corrigible'}, {'v': '1'}]},
                {'f': [{'v': 'othello'}, {'v': '1603'}, {'v': 'brawl'}, {'v': '2'}]},
                {'f': [{'v': 'othello'}, {'v': '1603'}, {'v': "'"}, {'v': '17'}]},
                {'f': [{'v': 'othello'}, {'v': '1603'}, {'v': 'troubled'},
                {'v': '1'}]}]

        self.correct_data_small = np.array(
            [('othello', 1603, 'brave', 3),
             ('othello', 1603, 'attended', 1),
             ('othello', 1603, 'treason', 1),
             ('othello', 1603, 'islanders', 1),
             ('othello', 1603, 'heed', 3),
             ('othello', 1603, 'alehouse', 1),
             ('othello', 1603, 'corrigible', 1),
             ('othello', 1603, 'brawl', 2),
             ('othello', 1603, "'", 17),
             ('othello', 1603, 'troubled', 1)
            ],
            dtype=[('corpus', 'S16'),
                   ('corpus_date', '<i8'),
                   ('word', 'S16'),
                   ('word_count', '<i8')]
        )

        self.correct_test_datatype = DataFrame(
            {'VALID_STRING' : ['PI'],
             'EMPTY_STRING' :  [""],
             'NULL_STRING' : [None],
             'VALID_INTEGER' : [3],
             'NULL_INTEGER' : [np.nan],
             'VALID_FLOAT' : [3.141592653589793],
             'NULL_FLOAT' : [np.nan],
             'UNIX_EPOCH' : [np.datetime64('1970-01-01T00:00:00.000000Z')],
             'VALID_TIMESTAMP' : [np.datetime64('2004-09-15T05:00:00.000000Z')],
             'NULL_TIMESTAMP' :[NaT],
             'TRUE_BOOLEAN' : [True],
             'FALSE_BOOLEAN' : [False],
             'NULL_BOOLEAN' : [None]
            }
        )[['VALID_STRING',
          'EMPTY_STRING',
          'NULL_STRING',
          'VALID_INTEGER',
          'NULL_INTEGER',
          'VALID_FLOAT',
          'NULL_FLOAT',
          'UNIX_EPOCH',
          'VALID_TIMESTAMP',
          'NULL_TIMESTAMP',
          'TRUE_BOOLEAN',
          'FALSE_BOOLEAN',
          'NULL_BOOLEAN']]

    @classmethod
    def setUpClass(cls):
        # Integration tests require a valid bigquery token
        # be present in the user's home directory. This
        # can be generated with 'bq init' in the command line
        super(TestGbq, cls).setUpClass()
        cls.dirpath = tm.get_data_path()
        home = os.path.expanduser("~")
        cls.bq_token = os.path.join(home, '.bigquery.v2.token')
        cls.fake_job_path = os.path.join(cls.dirpath, 'gbq_fake_job.txt')

        # If we're using a valid token, make a test dataset
        # Note, dataset functionality is beyond the scope
        # of the module under test, so we rely on the command
        # line utility for this.
        if os.path.exists(cls.bq_token):
            subprocess.call(['bq','mk', '-d', 'pandas_testing_dataset'])

    @classmethod
    def tearDownClass(cls):
        super(TestGbq, cls).tearDownClass()

        # If we're using a valid token, remove the test dataset
        # created.
        if os.path.exists(cls.bq_token):
            subprocess.call(['bq', 'rm', '-r', '-f', '-d', 'pandas_testing_dataset'])

    @with_connectivity_check
    def test_valid_authentication(self):
        # If the user has a token file, they should recieve a client from gbq._authenticate
        if not os.path.exists(self.bq_token):
            raise nose.SkipTest('Skipped because authentication information is not available.')

        self.assertTrue(gbq._authenticate is not None, 'Authentication To GBQ Failed')

    @with_connectivity_check
    def test_malformed_query(self):
        # If the user has a connection file, performing an invalid query should raise an error
        if not os.path.exists(self.bq_token):
            raise nose.SkipTest('Skipped because authentication information is not available.')
        else:
            self.assertRaises(bigquery_client.BigqueryInvalidQueryError,
                              gbq.read_gbq, "SELCET * FORM [publicdata:samples.shakespeare]")

    def test_type_conversion(self):
        # All BigQuery Types should be cast into appropriate numpy types
        sample_input = [('1.095292800E9', 'TIMESTAMP'),
                 ('false', 'BOOLEAN'),
                 ('2', 'INTEGER'),
                 ('3.14159', 'FLOAT'),
                 ('Hello World', 'STRING')]
        actual_output = [gbq._parse_entry(result[0],result[1]) for result in sample_input]
        sample_output = [np.datetime64('2004-09-16T00:00:00.000000Z'),
                  np.bool(False),
                  np.int('2'),
                  np.float('3.14159'),
                  'Hello World']
        self.assertEqual(actual_output, sample_output, 'A format conversion failed')

    def test_data_small(self):
        # Parsing a fixed page of data should return the proper fixed np.array()
        result_frame = gbq._parse_page(self.test_data_small,
                                        ['corpus','corpus_date','word','word_count'],
                                        ['STRING','INTEGER','STRING','INTEGER'],
                                        [object,np.dtype(int),object,np.dtype(int)]
                                      )
        tm.assert_frame_equal(DataFrame(result_frame), DataFrame(self.correct_data_small),
                        'An element in the result DataFrame didn\'t match the sample set')

    def test_index_column(self):
        # A user should be able to specify an index column for return
        result_frame = gbq._parse_data(FakeClient(), self.fake_job, index_col='word')
        correct_frame = DataFrame(self.correct_data_small)
        correct_frame.set_index('word', inplace=True)
        self.assertTrue(result_frame.index.name == correct_frame.index.name)

    def test_column_order(self):
        # A User should be able to specify the order in which columns are returned in the dataframe
        col_order = ['corpus_date', 'word_count', 'corpus', 'word']
        result_frame = gbq._parse_data(FakeClient(), self.fake_job, col_order=col_order)
        tm.assert_index_equal(result_frame.columns, DataFrame(self.correct_data_small)[col_order].columns)

    def test_column_order_plus_index(self):
        # A User should be able to specify an index and the order of THE REMAINING columns.. they should be notified
        # if they screw up
        col_order = ['corpus_date', 'word', 'corpus']
        result_frame = gbq._parse_data(FakeClient(), self.fake_job, index_col='word_count', col_order=col_order)
        correct_frame_small = DataFrame(self.correct_data_small)
        correct_frame_small.set_index('word_count',inplace=True)
        correct_frame_small = DataFrame(correct_frame_small)[col_order]
        tm.assert_index_equal(result_frame.columns, correct_frame_small.columns)

    @with_connectivity_check
    def test_download_dataset_larger_than_100k_rows(self):
        # Test for known BigQuery bug in datasets larger than 100k rows
        # http://stackoverflow.com/questions/19145587/bq-py-not-paging-results
        if not os.path.exists(self.bq_token):
            raise nose.SkipTest('Skipped because authentication information is not available.')

        client = gbq._authenticate()
        a = gbq.read_gbq("SELECT id, FROM [publicdata:samples.wikipedia] LIMIT 100005")
        self.assertTrue(len(a) == 100005)

    @with_connectivity_check
    def test_download_all_data_types(self):
        # Test that all available data types from BigQuery (as of now)
        # are handled properly
        if not os.path.exists(self.bq_token):
            raise nose.SkipTest('Skipped because authentication information is not available.')

        query = """SELECT "PI" as VALID_STRING,
                   "" as EMPTY_STRING,
                   STRING(NULL) as NULL_STRING,
                   INTEGER(3) as VALID_INTEGER,
                   INTEGER(NULL) as NULL_INTEGER,
                   PI() as VALID_FLOAT,
                   FLOAT(NULL) as NULL_FLOAT,
                   TIMESTAMP("1970-01-01 00:00:00") as UNIX_EPOCH,
                   TIMESTAMP("2004-09-15 05:00:00") as VALID_TIMESTAMP,
                   TIMESTAMP(NULL) as NULL_TIMESTAMP,
                   BOOLEAN(TRUE) as TRUE_BOOLEAN,
                   BOOLEAN(FALSE) as FALSE_BOOLEAN,
                   BOOLEAN(NULL) as NULL_BOOLEAN"""

        client = gbq._authenticate()
        a = gbq.read_gbq(query, col_order = ['VALID_STRING',
                                             'EMPTY_STRING',
                                             'NULL_STRING',
                                             'VALID_INTEGER',
                                             'NULL_INTEGER',
                                             'VALID_FLOAT',
                                             'NULL_FLOAT',
                                             'UNIX_EPOCH',
                                             'VALID_TIMESTAMP',
                                             'NULL_TIMESTAMP',
                                             'TRUE_BOOLEAN',
                                             'FALSE_BOOLEAN',
                                             'NULL_BOOLEAN'])

        tm.assert_frame_equal(a, self.correct_test_datatype)

    @with_connectivity_check
    def test_table_exists(self):
        # Given a table name in the format {dataset}.{tablename}, if a table exists,
        # the GetTableReference should accurately indicate this.
        # This could possibly change in future implementations of bq,
        # but it is the simplest way to provide users with appropriate
        # error messages regarding schemas.
        if not os.path.exists(self.bq_token):
            raise nose.SkipTest('Skipped because authentication information is not available.')

        client = gbq._authenticate()
        table_reference = client.GetTableReference("publicdata:samples.shakespeare")
        self.assertTrue(client.TableExists(table_reference))

    @with_connectivity_check
    def test_table__not_exists(self):
        # Test the inverse of `test_table_exists`
        if not os.path.exists(self.bq_token):
            raise nose.SkipTest('Skipped because authentication information is not available.')

        client = gbq._authenticate()
        table_reference = client.GetTableReference("publicdata:samples.does_not_exist")
        self.assertFalse(client.TableExists(table_reference))

    @with_connectivity_check
    def test_upload_new_table_schema_error(self):
        # Attempting to upload to a non-existent table without a schema should fail
        if not os.path.exists(self.bq_token):
            raise nose.SkipTest('Skipped because authentication information is not available.')

        df = DataFrame(self.correct_data_small)
        with self.assertRaises(gbq.SchemaMissing):
            gbq.to_gbq(df, 'pandas_testing_dataset.test_database', schema=None, col_order=None, if_exists='fail')

    @with_connectivity_check
    def test_upload_replace_schema_error(self):
        # Attempting to replace an existing table without specifying a schema should fail
        if not os.path.exists(self.bq_token):
            raise nose.SkipTest('Skipped because authentication information is not available.')

        df = DataFrame(self.correct_data_small)
        with self.assertRaises(gbq.SchemaMissing):
            gbq.to_gbq(df, 'pandas_testing_dataset.test_database', schema=None, col_order=None, if_exists='replace')

    @with_connectivity_check
    def test_upload_public_data_error(self):
        # Attempting to upload to a public, read-only, dataset should fail
        if not os.path.exists(self.bq_token):
            raise nose.SkipTest('Skipped because authentication information is not available.')

        array = [['TESTING_GBQ', 999999999, 'hi', 0, True, 9999999999, '00.000.00.000', 1, 'hola',
                 99999999, False, False, 1, 'Jedi', 11210]]
        df = DataFrame(array)
        with self.assertRaises(bigquery_client.BigqueryServiceError):
            gbq.to_gbq(df, 'publicdata:samples.wikipedia', schema=None, col_order=None, if_exists='append')

    @with_connectivity_check
    def test_upload_new_table(self):
        # Attempting to upload to a new table with valid data and a valid schema should succeed
        if not os.path.exists(self.bq_token):
            raise nose.SkipTest('Skipped because authentication information is not available.')

        schema = ['STRING', 'INTEGER', 'STRING', 'INTEGER', 'BOOLEAN',
                  'INTEGER', 'STRING', 'INTEGER',
                  'STRING', 'INTEGER', 'BOOLEAN', 'BOOLEAN',
                  'INTEGER', 'STRING', 'INTEGER']

        array = [['TESTING_GBQ', 999999999, 'hi', 0, True, 9999999999, '00.000.00.000', 1, 'hola',
                 99999999, False, False, 1, 'Jedi', 11210]]
        df = DataFrame(array, columns=['title','id','language','wp_namespace','is_redirect','revision_id',
                                       'contributor_ip','contributor_id','contributor_username','timestamp',
                                       'is_minor','is_bot','reversion_id','comment','num_characters'])
        gbq.to_gbq(df, 'pandas_testing_dataset.test_data2', schema=schema, col_order=None, if_exists='append')
        a = gbq.read_gbq("SELECT * FROM pandas_testing_dataset.test_data2")
        self.assertTrue((a == df).all().all())

    @with_connectivity_check
    def test_upload_bad_data_table(self):
        # Attempting to upload data that does not match schema should fail
        if not os.path.exists(self.bq_token):
            raise nose.SkipTest('Skipped because authentication information is not available.')

        schema = ['STRING', 'INTEGER', 'STRING', 'INTEGER', 'BOOLEAN',
                  'INTEGER', 'STRING', 'INTEGER',
                  'STRING', 'INTEGER', 'BOOLEAN', 'BOOLEAN',
                  'INTEGER', 'STRING', 'INTEGER']

        array = [['TESTING_GBQ\',', False, 'hi', 0, True, 'STRING IN INTEGER', '00.000.00.000', 1, 'hola',
                 99999999, -100, 1000, 1, 'Jedi', 11210]]
        df = DataFrame(array, columns=['title','id','language','wp_namespace','is_redirect','revision_id',
                                       'contributor_ip','contributor_id','contributor_username','timestamp',
                                       'is_minor','is_bot','reversion_id','comment','num_characters'])
        with self.assertRaises(bigquery_client.BigqueryServiceError):
            gbq.to_gbq(df, 'pandas_testing_dataset.test_data1', schema=schema, col_order=None, if_exists='append')

    @with_connectivity_check
    def test_invalid_column_name_schema(self):
        # Specifying a schema that contains an invalid column name should fail
        if not os.path.exists(self.bq_token):
            raise nose.SkipTest('Skipped because authentication information is not available.')

        schema = ['INCORRECT']
        df = DataFrame([[1]],columns=['fake'])
        with self.assertRaises(gbq.InvalidSchema):
            gbq.to_gbq(df, 'pandas_testing_dataset.test_data', schema=schema, col_order=None, if_exists='append')

    @with_connectivity_check
    def test_invalid_number_of_columns_schema(self):
        # Specifying a schema that does not have same shape as dataframe should fail
        if not os.path.exists(self.bq_token):
            raise nose.SkipTest('Skipped because authentication information is not available.')

        schema = ['INTEGER']
        df = DataFrame([[1, 'STRING']],columns=['fake','fake'])
        with self.assertRaises(gbq.InvalidSchema):
            gbq.to_gbq(df, 'pandas_testing_dataset.test_data4', schema=schema, col_order=None, if_exists='append')

    @with_connectivity_check
    def test_upload_fail_if_exists(self):
        # Attempting to upload to a new table with valid data and a valid schema should succeed
        if not os.path.exists(self.bq_token):
            raise nose.SkipTest('Skipped because authentication information is not available.')

        schema = ['STRING', 'INTEGER', 'STRING', 'INTEGER', 'BOOLEAN',
                  'INTEGER', 'STRING', 'INTEGER',
                  'STRING', 'INTEGER', 'BOOLEAN', 'BOOLEAN',
                  'INTEGER', 'STRING', 'INTEGER']

        array = [['TESTING_GBQ', 999999999, 'hi', 0, True, 9999999999, '00.000.00.000', 1, 'hola',
                 99999999, False, False, 1, 'Jedi', 11210]]
        df = DataFrame(array, columns=['title','id','language','wp_namespace','is_redirect','revision_id',
                                       'contributor_ip','contributor_id','contributor_username','timestamp',
                                       'is_minor','is_bot','reversion_id','comment','num_characters'])
        gbq.to_gbq(df, 'pandas_testing_dataset.test_data3', schema=schema, col_order=None, if_exists='fail')

        with self.assertRaises(gbq.TableExistsFail):
            gbq.to_gbq(df, 'pandas_testing_dataset.test_data3', schema=schema, col_order=None, if_exists='fail')

    @with_connectivity_check
    def test_upload_replace(self):
        # Attempting to overwrite an existing table with valid data and a valid schema should succeed
        if not os.path.exists(self.bq_token):
            raise nose.SkipTest('Skipped because authentication information is not available.')

        schema = ['STRING', 'INTEGER', 'STRING', 'INTEGER', 'BOOLEAN',
                  'INTEGER', 'STRING', 'INTEGER',
                  'STRING', 'INTEGER', 'BOOLEAN', 'BOOLEAN',
                  'INTEGER', 'STRING', 'INTEGER']

        # Setup an existing table
        array1 = [['', 1, '', 1, False, 1, '00.111.00.111', 1, 'hola',
                 1, True, True, 1, 'Sith', 1]]
        df1 = DataFrame(array1, columns=['title','id','language','wp_namespace','is_redirect','revision_id',
                                       'contributor_ip','contributor_id','contributor_username','timestamp',
                                       'is_minor','is_bot','reversion_id','comment','num_characters'])
        gbq.to_gbq(df1, 'pandas_testing_dataset.test_data5', schema=schema, col_order=None, if_exists='fail')

        array2 = [['TESTING_GBQ', 999999999, 'hi', 0, True, 9999999999, '00.000.00.000', 1, 'hola',
                 99999999, False, False, 1, 'Jedi', 11210]]

        # Overwrite the existing table with different data
        df2 = DataFrame(array2, columns=['title','id','language','wp_namespace','is_redirect','revision_id',
                                       'contributor_ip','contributor_id','contributor_username','timestamp',
                                       'is_minor','is_bot','reversion_id','comment','num_characters'])
        gbq.to_gbq(df2, 'pandas_testing_dataset.test_data5', schema=schema, col_order=None, if_exists='replace')

        # Read the table and confirm the new data is all that is there
        a = gbq.read_gbq("SELECT * FROM pandas_testing_dataset.test_data5")
        self.assertTrue((a == df2).all().all())


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
