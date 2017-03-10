import pytest

import pandas.util.testing as tm


@pytest.fixture(params=['bsr', 'coo', 'csc', 'csr', 'dia', 'dok', 'lil'])
def spmatrix(request):
    tm._skip_if_no_scipy()
    from scipy import sparse
    return getattr(sparse, request.param + '_matrix')
