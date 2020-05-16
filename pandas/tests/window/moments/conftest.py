import numpy as np
from numpy.random import randn
import pytest

from pandas import Series


@pytest.fixture
def binary_ew_data():
    A = Series(randn(50), index=np.arange(50))
    B = A[2:] + randn(48)

    A[:10] = np.NaN
    B[-10:] = np.NaN
    return A, B


@pytest.fixture(params=[0, 1, 2])
def min_periods(request):
    return request.param
