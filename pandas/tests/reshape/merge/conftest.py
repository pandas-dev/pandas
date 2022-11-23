import pytest


@pytest.fixture(params=["left", "right", "inner", "outer"])
def how(request):
    """Boolean sort keyword for concat and DataFrame.append."""
    return request.param
