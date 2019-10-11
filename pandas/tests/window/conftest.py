import pytest


@pytest.fixture(params=[True, False])
def raw(request):
    return request.param


@pytest.fixture(
    params=[
        "triang",
        "blackman",
        "hamming",
        "bartlett",
        "bohman",
        "blackmanharris",
        "nuttall",
        "barthann",
    ]
)
def win_types(request):
    return request.param


@pytest.fixture(params=["kaiser", "gaussian", "general_gaussian", "exponential"])
def win_types_special(request):
    return request.param


@pytest.fixture(
    params=["sum", "mean", "median", "max", "min", "var", "std", "kurt", "skew"]
)
def arithmetic_win_operators(request):
    return request.param


@pytest.fixture(params=["right", "left", "both", "neither"])
def closed(request):
    return request.param


@pytest.fixture(params=[True, False])
def center(request):
    return request.param


@pytest.fixture(params=[None, 1])
def min_periods(request):
    return request.param
