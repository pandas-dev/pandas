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
    return request.param


@pytest.fixture
def json_engine_pyarrow_xfail(request):
    """
    Fixture that xfails a test if the engine is pyarrow.
    """
    engine = request.getfixturevalue("engine")
    if engine == "pyarrow":
        mark = pytest.mark.xfail(reason="pyarrow doesn't support this.")
        request.node.add_marker(mark)
