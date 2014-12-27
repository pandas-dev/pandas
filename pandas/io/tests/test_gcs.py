import pandas.io.gcs as gcs
import pandas.util.testing as tm
import pandas as pd
from time import sleep

import subprocess

import nose

PROJECT_ID = "eighth-physics-623"

def missing_gsutil():
    try:
        subprocess.call(['which', 'gsutil'])
        return False
    except OSError:
        return True

def test_requirements():
    try:
        gcs._test_imports()
    except (ImportError, NotImplementedError) as import_exception:
        raise nose.SkipTest(import_exception)

class TestGCSConnectorIntegration(tm.TestCase):
    def setUp(self):
        test_requirements()

        if not PROJECT_ID:
            raise nose.SkipTest("Cannot run integration tests without a project id")

        self.sut = gcs.GCSConnector(PROJECT_ID)

    def test_should_be_able_to_make_a_connector(self):
        self.assertTrue(self.sut is not None, 'Could not create a GCSConnector')

    def test_should_be_able_to_get_valid_credentials(self):
        credentials = self.sut.get_credentials()
        self.assertFalse(credentials.invalid, 'Returned credentials invalid')

    def test_should_be_able_to_get_a_cloudstorage_service(self):
        credentials = self.sut.get_credentials()
        cloudstorage_service = self.sut.get_service(credentials)
        self.assertTrue(cloudstorage_service is not None, 'No service returned')

class TestGCSReadWrite(tm.TestCase):

    def setUp(self):
        test_requirements()

        if not PROJECT_ID:
            raise nose.SkipTest("Cannot run write tests without a project id")
        if missing_gsutil():
            raise nose.SkipTest("Cannot run write tests without a gsutil.")

    @classmethod
    def setUpClass(cls):
        if PROJECT_ID and not missing_gsutil():
            subprocess.call(["gsutil", "mb", "gs://pandas-write-test"])

    @classmethod
    def tearDownClass(cls):
        if PROJECT_ID and not missing_gsutil():
            subprocess.call(['gsutil', 'rb', '-f', 'gs://pandas-write-test'])

    def test_upload_data(self):
        #test that we can upload data and that from_csv_kwargs works correctly
        fake_df = pd.DataFrame({'a': [1, 2, 3]}, index=pd.Index([1, 2, 3], name='test_idx'))
        gcs.to_gcs(fake_df, "pandas-write-test", "new_test.csv", project_id=PROJECT_ID)

        sleep(60)

        result = gcs.from_gcs("pandas-write-test", "new_test.csv", project_id=PROJECT_ID,
                              from_csv_kwargs={'index_col': 'test_idx'})

        tm.assert_frame_equal(fake_df, result)

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x'], exit=False)
