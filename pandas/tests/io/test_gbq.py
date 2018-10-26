import pytest
from datetime import datetime
import pytz
import platform
import os

try:
    from unittest import mock
except ImportError:
    mock = pytest.importorskip("mock")

import numpy as np
import pandas as pd
from pandas import compat, DataFrame
from pandas.compat import range
import pandas.util.testing as tm


api_exceptions = pytest.importorskip("google.api_core.exceptions")
bigquery = pytest.importorskip("google.cloud.bigquery")
service_account = pytest.importorskip("google.oauth2.service_account")
pandas_gbq = pytest.importorskip('pandas_gbq')

PROJECT_ID = None
PRIVATE_KEY_JSON_PATH = None
PRIVATE_KEY_JSON_CONTENTS = None

if compat.PY3:
    DATASET_ID = 'pydata_pandas_bq_testing_py3'
else:
    DATASET_ID = 'pydata_pandas_bq_testing_py2'

TABLE_ID = 'new_test'
DESTINATION_TABLE = "{0}.{1}".format(DATASET_ID + "1", TABLE_ID)

VERSION = platform.python_version()


def _skip_if_no_project_id():
    if not _get_project_id():
        pytest.skip(
            "Cannot run integration tests without a project id")


def _skip_if_no_private_key_path():
    if not _get_private_key_path():
        pytest.skip("Cannot run integration tests without a "
                    "private key json file path")


def _in_travis_environment():
    return 'TRAVIS_BUILD_DIR' in os.environ and \
           'GBQ_PROJECT_ID' in os.environ


def _get_project_id():
    if _in_travis_environment():
        return os.environ.get('GBQ_PROJECT_ID')
    return PROJECT_ID or os.environ.get('GBQ_PROJECT_ID')


def _get_private_key_path():
    if _in_travis_environment():
        return os.path.join(*[os.environ.get('TRAVIS_BUILD_DIR'), 'ci',
                              'travis_gbq.json'])

    private_key_path = PRIVATE_KEY_JSON_PATH
    if not private_key_path:
        private_key_path = os.environ.get('GBQ_GOOGLE_APPLICATION_CREDENTIALS')
    return private_key_path


def _get_client():
    project_id = _get_project_id()
    credentials = None

    private_key_path = _get_private_key_path()
    if private_key_path:
        credentials = service_account.Credentials.from_service_account_file(
            private_key_path)

    return bigquery.Client(project=project_id, credentials=credentials)


def make_mixed_dataframe_v2(test_size):
    # create df to test for all BQ datatypes except RECORD
    bools = np.random.randint(2, size=(1, test_size)).astype(bool)
    flts = np.random.randn(1, test_size)
    ints = np.random.randint(1, 10, size=(1, test_size))
    strs = np.random.randint(1, 10, size=(1, test_size)).astype(str)
    times = [datetime.now(pytz.timezone('US/Arizona'))
             for t in range(test_size)]
    return DataFrame({'bools': bools[0],
                      'flts': flts[0],
                      'ints': ints[0],
                      'strs': strs[0],
                      'times': times[0]},
                     index=range(test_size))


def test_read_gbq_without_dialect_warns_future_change(monkeypatch):
    # Default dialect is changing to standard SQL. See:
    # https://github.com/pydata/pandas-gbq/issues/195
    mock_read_gbq = mock.Mock()
    mock_read_gbq.return_value = DataFrame([[1.0]])
    monkeypatch.setattr(pandas_gbq, 'read_gbq', mock_read_gbq)
    with tm.assert_produces_warning(FutureWarning):
        pd.read_gbq("SELECT 1")


@pytest.mark.single
class TestToGBQIntegrationWithServiceAccountKeyPath(object):

    @classmethod
    def setup_class(cls):
        # - GLOBAL CLASS FIXTURES -
        # put here any instruction you want to execute only *ONCE* *BEFORE*
        # executing *ALL* tests described below.

        _skip_if_no_project_id()
        _skip_if_no_private_key_path()

        cls.client = _get_client()
        cls.dataset = cls.client.dataset(DATASET_ID + "1")
        try:
            # Clean-up previous test runs.
            cls.client.delete_dataset(cls.dataset, delete_contents=True)
        except api_exceptions.NotFound:
            pass  # It's OK if the dataset doesn't already exist.

        cls.client.create_dataset(bigquery.Dataset(cls.dataset))

    @classmethod
    def teardown_class(cls):
        # - GLOBAL CLASS FIXTURES -
        # put here any instruction you want to execute only *ONCE* *AFTER*
        # executing all tests.
        cls.client.delete_dataset(cls.dataset, delete_contents=True)

    def test_roundtrip(self):
        destination_table = DESTINATION_TABLE + "1"

        test_size = 20001
        df = make_mixed_dataframe_v2(test_size)

        df.to_gbq(destination_table, _get_project_id(), chunksize=None,
                  private_key=_get_private_key_path())

        result = pd.read_gbq("SELECT COUNT(*) AS num_rows FROM {0}"
                             .format(destination_table),
                             project_id=_get_project_id(),
                             private_key=_get_private_key_path(),
                             dialect="standard")
        assert result['num_rows'][0] == test_size
