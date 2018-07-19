import numpy as np
import pytest

from pandas import DataFrame, date_range, read_csv
from pandas.compat import StringIO
from pandas.io.common import is_gcs_url
from pandas.util import _test_decorators as td
from pandas.util.testing import assert_frame_equal


def test_is_gcs_url():
    assert is_gcs_url("gcs://pandas/somethingelse.com")
    assert is_gcs_url("gs://pandas/somethingelse.com")
    assert not is_gcs_url("s3://pandas/somethingelse.com")


@td.skip_if_no('gcsfs')
def test_read_csv_gcs(mock):
    df1 = DataFrame({'int': [1, 3], 'float': [2.0, np.nan], 'str': ['t', 's'],
                     'dt': date_range('2018-06-18', periods=2)})
    with mock.patch('gcsfs.GCSFileSystem') as MockFileSystem:
        instance = MockFileSystem.return_value
        instance.open.return_value = StringIO(df1.to_csv(index=False))
        df2 = read_csv('gs://test/test.csv', parse_dates=['dt'])

    assert_frame_equal(df1, df2)


@td.skip_if_no('gcsfs')
def test_gcs_get_filepath_or_buffer(mock):
    df1 = DataFrame({'int': [1, 3], 'float': [2.0, np.nan], 'str': ['t', 's'],
                     'dt': date_range('2018-06-18', periods=2)})
    with mock.patch('pandas.io.gcs.get_filepath_or_buffer') as MockGetFilepath:
        MockGetFilepath.return_value = (StringIO(df1.to_csv(index=False)),
                                        None, None, False)
        df2 = read_csv('gs://test/test.csv', parse_dates=['dt'])

    assert_frame_equal(df1, df2)
    assert MockGetFilepath.called


@pytest.mark.skipif(td.safe_import('gcsfs'),
                    reason='Only check when gcsfs not installed')
def test_gcs_not_present_exception():
    with pytest.raises(ImportError) as e:
        read_csv('gs://test/test.csv')
        assert 'gcsfs library is required' in str(e.value)
