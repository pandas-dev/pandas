import pytest
import pandas
import pandas.util._test_decorators as td


@pytest.fixture(params=[None, 'gzip', 'bz2',
                        pytest.param('xz', marks=td.skip_if_no_lzma)])
def compression(request):
    """
    Fixture for trying common compression types in compression tests
    """
    return request.param
