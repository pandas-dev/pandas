import pytest


@pytest.fixture(params=[True, False])
def check_dtype(request):
    return request.param


@pytest.fixture(params=[True, False])
def check_exact(request):
    return request.param


@pytest.fixture(params=[True, False])
def check_index_type(request):
    return request.param


@pytest.fixture(params=[True, False])
def check_less_precise(request):
    return request.param


@pytest.fixture(params=[True, False])
def check_categorical(request):
    return request.param
