from pandas.util import testing as tm

from pandas.io.common import _is_s3_url


class TestS3URL(tm.TestCase):

    def test_is_s3_url(self):
        self.assertTrue(_is_s3_url("s3://pandas/somethingelse.com"))
        self.assertFalse(_is_s3_url("s4://pandas/somethingelse.com"))
