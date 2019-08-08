import pytest

@pytest.fixture(params=["split", "records", "index", "columns", "values", "table"])
def df_orient(request):
    """
     Fixture for orients applicable to a DataFrame.
     """
    return request.param
