import pytest

@pytest.fixture(params=["split", "records", "index", "columns", "values"])
def df_orient(request):
    """
     Fixture for orients applicable to a DataFrame, excluding the table format.
     """
    return request.param
