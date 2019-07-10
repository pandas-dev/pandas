import pytest

import pandas.util.testing as tm

from pandas.io.parsers import read_csv


@pytest.fixture
def frame(float_frame):
    return float_frame[:10]


@pytest.fixture
def tsframe():
    return tm.makeTimeDataFrame()[:5]


@pytest.fixture(params=[True, False])
def merge_cells(request):
    return request.param


@pytest.fixture
def df_ref():
    """
    Obtain the reference data from read_csv with the Python engine.
    """
    df_ref = read_csv("test1.csv", index_col=0, parse_dates=True, engine="python")
    return df_ref


@pytest.fixture(params=[".xls", ".xlsx", ".xlsm", ".ods"])
def read_ext(request):
    """
    Valid extensions for reading Excel files.
    """
    return request.param
