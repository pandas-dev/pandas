from pandas.io.common import _is_hdfs_url


class TestHDFSURL(object):

    def test_is_hdfs_url(self):
        assert _is_hdfs_url("hdfs://pandas/somethingelse.com")
        assert not _is_hdfs_url("hdf://pandas/somethingelse.com")
