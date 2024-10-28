import pytest


@pytest.fixture(params=[True, False])
def check_dtype(request):
    """
    Fixture returning `True` or `False`, determining whether to check
    if the `dtype` is identical or not, when comparing two data structures,
    e.g. `Series`, `SparseArray` or `DataFrame`.
    """
    return request.param


@pytest.fixture(params=[True, False])
def check_exact(request):
    """
    Fixture returning `True` or `False`, determining whether to
    compare floating point numbers exactly or not.
    """
    return request.param


@pytest.fixture(params=[True, False])
def check_index_type(request):
    """
    Fixture returning `True` or `False`, determining whether to check
    if the `Index` types are identical or not.
    """
    return request.param


@pytest.fixture(params=[0.5e-3, 0.5e-5])
def rtol(request):
    """
    Fixture returning 0.5e-3 or 0.5e-5. Those values are used as relative tolerance.
    """
    return request.param


@pytest.fixture(params=[True, False])
def check_categorical(request):
    """
    Fixture returning `True` or `False`, determining whether to
    compare internal `Categorical` exactly or not.
    """
    return request.param
