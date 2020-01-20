import pytest

import pandas.util._test_decorators as td


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


@pytest.fixture(params=[True, False])
def parallel(request):
    """parallel keyword argument for numba.jit"""
    return request.param


@pytest.fixture(params=[True, False])
def nogil(request):
    """nogil keyword argument for numba.jit"""
    return request.param


@pytest.fixture(params=[True, False])
def nopython(request):
    """nopython keyword argument for numba.jit"""
    return request.param


@pytest.fixture(
    params=[pytest.param("numba", marks=td.skip_if_no("numba", "0.46.0")), "cython"]
)
def engine(request):
    """engine keyword argument for rolling.apply"""
    return request.param


@pytest.fixture(
    params=[
        pytest.param(("numba", True), marks=td.skip_if_no("numba", "0.46.0")),
        ("cython", True),
        ("cython", False),
    ]
)
def engine_and_raw(request):
    """engine and raw keyword arguments for rolling.apply"""
    return request.param
