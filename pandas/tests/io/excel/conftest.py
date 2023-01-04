import pytest

from pandas.compat import is_platform_windows
import pandas.util._test_decorators as td

import pandas._testing as tm

from pandas.io.parsers import read_csv


@pytest.fixture(name="frame")
def fixture_frame(float_frame):
    """
    Returns the first ten items in fixture "float_frame".
    """
    return float_frame[:10]


@pytest.fixture(name="tsframe")
def fixture_tsframe():
    return tm.makeTimeDataFrame()[:5]


@pytest.fixture(name="merge_cells", params=[True, False])
def fixture_merge_cells(request):
    return request.param


@pytest.fixture(name="df_ref")
def fixture_df_ref(datapath):
    """
    Obtain the reference data from read_csv with the Python engine.
    """
    filepath = datapath("io", "data", "csv", "test1.csv")
    df_ref = read_csv(filepath, index_col=0, parse_dates=True, engine="python")
    return df_ref


@pytest.fixture(name="read_ext", params=[".xls", ".xlsx", ".xlsm", ".ods", ".xlsb"])
def fixture_read_ext(request):
    """
    Valid extensions for reading Excel files.
    """
    return request.param


# Checking for file leaks can hang on Windows CI
@pytest.fixture(name="check_for_file_leaks", autouse=not is_platform_windows())
def fixture_check_for_file_leaks():
    """
    Fixture to run around every test to ensure that we are not leaking files.

    See also
    --------
    _test_decorators.check_file_leaks
    """
    # GH#30162
    psutil = td.safe_import("psutil")
    if not psutil:
        yield

    else:
        proc = psutil.Process()
        flist = proc.open_files()
        yield
        flist2 = proc.open_files()
        assert flist == flist2
