import pytest


@pytest.fixture(name="sort", params=[True, False])
def fixture_sort(request):
    """Boolean sort keyword for concat and DataFrame.append."""
    return request.param
