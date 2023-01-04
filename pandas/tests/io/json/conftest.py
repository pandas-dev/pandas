import pytest


@pytest.fixture(
    name="orient", params=["split", "records", "index", "columns", "values"]
)
def fixture_orient(request):
    """
    Fixture for orients excluding the table format.
    """
    return request.param
