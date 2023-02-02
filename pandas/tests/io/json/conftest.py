import pytest


@pytest.fixture(params=["split", "records", "index", "columns", "values"])
def orient(request):
    """
    Fixture for orients excluding the table format.
    """
    return request.param


@pytest.fixture
def json_dir_path(datapath):
    """
    The directory path to the data files needed for parser tests.
    """
    return datapath("io", "json", "data")


@pytest.fixture(params=["ujson", "pyarrow"])
def engine(request):
    if request.param == "pyarrow":
        pytest.importorskip("pyarrow.json")
        return request.param
    else:
        return request.param
