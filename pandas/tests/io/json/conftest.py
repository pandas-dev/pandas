import pytest


@pytest.fixture(params=["split", "records", "index", "columns", "values"])
def orient(request):
    """
    Fixture for orients excluding the table format.
    """
    return request.param


@pytest.fixture(params=["ujson", "orjson"])
def json_engines_no_ndjson(request):
    """
    Fixture for json decoders that won't read newline delimited JSON.
    """
    return request.param
