import pytest


@pytest.fixture(name="check_dtype", params=[True, False])
def fixture_check_dtype(request):
    return request.param


@pytest.fixture(name="check_exact", params=[True, False])
def fixture_check_exact(request):
    return request.param


@pytest.fixture(name="check_index_type", params=[True, False])
def fixture_check_index_type(request):
    return request.param


@pytest.fixture(name="rtol", params=[0.5e-3, 0.5e-5])
def fixture_rtol(request):
    return request.param


@pytest.fixture(name="check_categorical", params=[True, False])
def fixture_check_categorical(request):
    return request.param
