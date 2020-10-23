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


base_functions_list = [
    (lambda v: Series(v).count(), None, "count"),
    (lambda v: Series(v).max(), None, "max"),
    (lambda v: Series(v).min(), None, "min"),
    (lambda v: Series(v).sum(), None, "sum"),
    (lambda v: Series(v).mean(), None, "mean"),
    (lambda v: Series(v).std(), 1, "std"),
    (lambda v: Series(v).cov(Series(v)), None, "cov"),
    (lambda v: Series(v).corr(Series(v)), None, "corr"),
    (lambda v: Series(v).var(), 1, "var"),
    # restore once GH 8086 is fixed
    # lambda v: Series(v).skew(), 3, 'skew'),
    # (lambda v: Series(v).kurt(), 4, 'kurt'),
    # restore once GH 8084 is fixed
    # lambda v: Series(v).quantile(0.3), None, 'quantile'),
    (lambda v: Series(v).median(), None, "median"),
    (np.nanmax, 1, "max"),
    (np.nanmin, 1, "min"),
    (np.nansum, 1, "sum"),
    (np.nanmean, 1, "mean"),
    (lambda v: np.nanstd(v, ddof=1), 1, "std"),
    (lambda v: np.nanvar(v, ddof=1), 1, "var"),
    (np.nanmedian, 1, "median"),
]

no_nan_functions_list = [
    (np.max, None, "max"),
    (np.min, None, "min"),
    (np.sum, None, "sum"),
    (np.mean, None, "mean"),
    (lambda v: np.std(v, ddof=1), 1, "std"),
    (lambda v: np.var(v, ddof=1), 1, "var"),
    (np.median, None, "median"),
]


@pytest.fixture(scope="session")
def base_functions():
    """Fixture for base functions.

    Returns
    -------
    List of tuples: (applied function, require_min_periods, name of applied function)
    """
    return base_functions_list


@pytest.fixture(scope="session")
def no_nan_functions():
    """Fixture for no nan functions.

    Returns
    -------
    List of tuples: (applied function, require_min_periods, name of applied function)
    """
    return no_nan_functions_list
