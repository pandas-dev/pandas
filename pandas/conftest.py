import os
import importlib

import pytest

import pandas
import numpy as np
import pandas as pd
from pandas.compat import PY3
import pandas.util._test_decorators as td
import hypothesis


hypothesis.settings.register_profile(
    "ci",
    suppress_health_check=(hypothesis.HealthCheck.too_slow,)
)
hypothesis.settings.load_profile("ci")


def pytest_addoption(parser):
    parser.addoption("--skip-slow", action="store_true",
                     help="skip slow tests")
    parser.addoption("--skip-network", action="store_true",
                     help="skip network tests")
    parser.addoption("--run-high-memory", action="store_true",
                     help="run high memory tests")
    parser.addoption("--only-slow", action="store_true",
                     help="run only slow tests")
    parser.addoption("--strict-data-files", action="store_true",
                     help="Fail if a test is skipped for missing data file.")


def pytest_runtest_setup(item):
    if 'slow' in item.keywords and item.config.getoption("--skip-slow"):
        pytest.skip("skipping due to --skip-slow")

    if 'slow' not in item.keywords and item.config.getoption("--only-slow"):
        pytest.skip("skipping due to --only-slow")

    if 'network' in item.keywords and item.config.getoption("--skip-network"):
        pytest.skip("skipping due to --skip-network")

    if 'high_memory' in item.keywords and not item.config.getoption(
            "--run-high-memory"):
        pytest.skip(
            "skipping high memory test since --run-high-memory was not set")


# Configurations for all tests and all test modules

@pytest.fixture(autouse=True)
def configure_tests():
    pd.set_option('chained_assignment', 'raise')


# For running doctests: make np and pd names available

@pytest.fixture(autouse=True)
def add_imports(doctest_namespace):
    doctest_namespace['np'] = np
    doctest_namespace['pd'] = pd


@pytest.fixture(params=['bsr', 'coo', 'csc', 'csr', 'dia', 'dok', 'lil'])
def spmatrix(request):
    from scipy import sparse
    return getattr(sparse, request.param + '_matrix')


@pytest.fixture(params=[0, 1, 'index', 'columns'],
                ids=lambda x: "axis {!r}".format(x))
def axis(request):
    """
     Fixture for returning the axis numbers of a DataFrame.
     """
    return request.param


axis_frame = axis


@pytest.fixture(params=[0, 'index'], ids=lambda x: "axis {!r}".format(x))
def axis_series(request):
    """
     Fixture for returning the axis numbers of a Series.
     """
    return request.param


@pytest.fixture
def ip():
    """
    Get an instance of IPython.InteractiveShell.

    Will raise a skip if IPython is not installed.
    """

    pytest.importorskip('IPython', minversion="6.0.0")
    from IPython.core.interactiveshell import InteractiveShell
    return InteractiveShell()


@pytest.fixture(params=[True, False, None])
def observed(request):
    """ pass in the observed keyword to groupby for [True, False]
    This indicates whether categoricals should return values for
    values which are not in the grouper [False / None], or only values which
    appear in the grouper [True]. [None] is supported for future compatiblity
    if we decide to change the default (and would need to warn if this
    parameter is not passed)"""
    return request.param


_all_arithmetic_operators = ['__add__', '__radd__',
                             '__sub__', '__rsub__',
                             '__mul__', '__rmul__',
                             '__floordiv__', '__rfloordiv__',
                             '__truediv__', '__rtruediv__',
                             '__pow__', '__rpow__',
                             '__mod__', '__rmod__']
if not PY3:
    _all_arithmetic_operators.extend(['__div__', '__rdiv__'])


@pytest.fixture(params=_all_arithmetic_operators)
def all_arithmetic_operators(request):
    """
    Fixture for dunder names for common arithmetic operations
    """
    return request.param


_all_numeric_reductions = ['sum', 'max', 'min',
                           'mean', 'prod', 'std', 'var', 'median']


@pytest.fixture(params=_all_numeric_reductions)
def all_numeric_reductions(request):
    """
    Fixture for numeric reduction names
    """
    return request.param


_all_boolean_reductions = ['all', 'any']


@pytest.fixture(params=_all_boolean_reductions)
def all_boolean_reductions(request):
    """
    Fixture for boolean reduction names
    """
    return request.param


_cython_table = pd.core.base.SelectionMixin._cython_table.items()


@pytest.fixture(params=list(_cython_table))
def cython_table_items(request):
    return request.param


def _get_cython_table_params(ndframe, func_names_and_expected):
    """combine frame, functions from SelectionMixin._cython_table
    keys and expected result.

    Parameters
    ----------
    ndframe : DataFrame or Series
    func_names_and_expected : Sequence of two items
        The first item is a name of a NDFrame method ('sum', 'prod') etc.
        The second item is the expected return value

    Returns
    -------
    results : list
        List of three items (DataFrame, function, expected result)
    """
    results = []
    for func_name, expected in func_names_and_expected:
        results.append((ndframe, func_name, expected))
        results += [(ndframe, func, expected) for func, name in _cython_table
                    if name == func_name]
    return results


@pytest.fixture(params=['__eq__', '__ne__', '__le__',
                        '__lt__', '__ge__', '__gt__'])
def all_compare_operators(request):
    """
    Fixture for dunder names for common compare operations

    * >=
    * >
    * ==
    * !=
    * <
    * <=
    """
    return request.param


@pytest.fixture(params=[None, 'gzip', 'bz2', 'zip',
                        pytest.param('xz', marks=td.skip_if_no_lzma)])
def compression(request):
    """
    Fixture for trying common compression types in compression tests
    """
    return request.param


@pytest.fixture(params=['gzip', 'bz2', 'zip',
                        pytest.param('xz', marks=td.skip_if_no_lzma)])
def compression_only(request):
    """
    Fixture for trying common compression types in compression tests excluding
    uncompressed case
    """
    return request.param


@pytest.fixture(params=[True, False])
def writable(request):
    """
    Fixture that an array is writable
    """
    return request.param


@pytest.fixture(scope='module')
def datetime_tz_utc():
    from datetime import timezone
    return timezone.utc


@pytest.fixture(params=['inner', 'outer', 'left', 'right'])
def join_type(request):
    """
    Fixture for trying all types of join operations
    """
    return request.param


@pytest.fixture
def datapath(request):
    """Get the path to a data file.

    Parameters
    ----------
    path : str
        Path to the file, relative to ``pandas/tests/``

    Returns
    -------
    path : path including ``pandas/tests``.

    Raises
    ------
    ValueError
        If the path doesn't exist and the --strict-data-files option is set.
    """
    BASE_PATH = os.path.join(os.path.dirname(__file__), 'tests')

    def deco(*args):
        path = os.path.join(BASE_PATH, *args)
        if not os.path.exists(path):
            if request.config.getoption("--strict-data-files"):
                msg = "Could not find file {} and --strict-data-files is set."
                raise ValueError(msg.format(path))
            else:
                msg = "Could not find {}."
                pytest.skip(msg.format(path))
        return path
    return deco


@pytest.fixture
def iris(datapath):
    """The iris dataset as a DataFrame."""
    return pandas.read_csv(datapath('data', 'iris.csv'))


@pytest.fixture(params=['nlargest', 'nsmallest'])
def nselect_method(request):
    """
    Fixture for trying all nselect methods
    """
    return request.param


@pytest.fixture(params=['left', 'right', 'both', 'neither'])
def closed(request):
    """
    Fixture for trying all interval closed parameters
    """
    return request.param


@pytest.fixture(params=[None, np.nan, pd.NaT, float('nan'), np.float('NaN')])
def nulls_fixture(request):
    """
    Fixture for each null type in pandas
    """
    return request.param


nulls_fixture2 = nulls_fixture  # Generate cartesian product of nulls_fixture


@pytest.fixture(params=[None, np.nan, pd.NaT])
def unique_nulls_fixture(request):
    """
    Fixture for each null type in pandas, each null type exactly once
    """
    return request.param


# Generate cartesian product of unique_nulls_fixture:
unique_nulls_fixture2 = unique_nulls_fixture


TIMEZONES = [None, 'UTC', 'US/Eastern', 'Asia/Tokyo', 'dateutil/US/Pacific',
             'dateutil/Asia/Singapore']


@td.parametrize_fixture_doc(str(TIMEZONES))
@pytest.fixture(params=TIMEZONES)
def tz_naive_fixture(request):
    """
    Fixture for trying timezones including default (None): {0}
    """
    return request.param


@td.parametrize_fixture_doc(str(TIMEZONES[1:]))
@pytest.fixture(params=TIMEZONES[1:])
def tz_aware_fixture(request):
    """
    Fixture for trying explicit timezones: {0}
    """
    return request.param


UNSIGNED_INT_DTYPES = ["uint8", "uint16", "uint32", "uint64"]
SIGNED_INT_DTYPES = [int, "int8", "int16", "int32", "int64"]
ALL_INT_DTYPES = UNSIGNED_INT_DTYPES + SIGNED_INT_DTYPES

FLOAT_DTYPES = [float, "float32", "float64"]
COMPLEX_DTYPES = [complex, "complex64", "complex128"]
STRING_DTYPES = [str, 'str', 'U']

ALL_REAL_DTYPES = FLOAT_DTYPES + ALL_INT_DTYPES
ALL_NUMPY_DTYPES = ALL_REAL_DTYPES + COMPLEX_DTYPES + STRING_DTYPES


@pytest.fixture(params=STRING_DTYPES)
def string_dtype(request):
    """Parametrized fixture for string dtypes.

    * str
    * 'str'
    * 'U'
    """
    return request.param


@pytest.fixture(params=FLOAT_DTYPES)
def float_dtype(request):
    """
    Parameterized fixture for float dtypes.

    * float32
    * float64
    """

    return request.param


@pytest.fixture(params=COMPLEX_DTYPES)
def complex_dtype(request):
    """
    Parameterized fixture for complex dtypes.

    * complex64
    * complex128
    """

    return request.param


@pytest.fixture(params=SIGNED_INT_DTYPES)
def sint_dtype(request):
    """
    Parameterized fixture for signed integer dtypes.

    * int8
    * int16
    * int32
    * int64
    """

    return request.param


@pytest.fixture(params=UNSIGNED_INT_DTYPES)
def uint_dtype(request):
    """
    Parameterized fixture for unsigned integer dtypes.

    * uint8
    * uint16
    * uint32
    * uint64
    """

    return request.param


@pytest.fixture(params=ALL_INT_DTYPES)
def any_int_dtype(request):
    """
    Parameterized fixture for any integer dtypes.

    * int8
    * uint8
    * int16
    * uint16
    * int32
    * uint32
    * int64
    * uint64
    """

    return request.param


@pytest.fixture(params=ALL_REAL_DTYPES)
def any_real_dtype(request):
    """
    Parameterized fixture for any (purely) real numeric dtypes.

    * int8
    * uint8
    * int16
    * uint16
    * int32
    * uint32
    * int64
    * uint64
    * float32
    * float64
    """

    return request.param


@pytest.fixture(params=ALL_NUMPY_DTYPES)
def any_numpy_dtype(request):
    """
    Parameterized fixture for all numpy dtypes.

    * int8
    * uint8
    * int16
    * uint16
    * int32
    * uint32
    * int64
    * uint64
    * float32
    * float64
    * complex64
    * complex128
    * str
    * 'str'
    * 'U'
    """

    return request.param


@pytest.fixture
def mock():
    """
    Fixture providing the 'mock' module.

    Uses 'unittest.mock' for Python 3. Attempts to import the 3rd party 'mock'
    package for Python 2, skipping if not present.
    """
    if PY3:
        return importlib.import_module("unittest.mock")
    else:
        return pytest.importorskip("mock")


# ----------------------------------------------------------------
# Global setup for tests using Hypothesis

from hypothesis import strategies as st

# Registering these strategies makes them globally available via st.from_type,
# which is use for offsets in tests/tseries/offsets/test_offsets_properties.py
for name in 'MonthBegin MonthEnd BMonthBegin BMonthEnd'.split():
    cls = getattr(pd.tseries.offsets, name)
    st.register_type_strategy(cls, st.builds(
        cls,
        n=st.integers(-99, 99),
        normalize=st.booleans(),
    ))

for name in 'YearBegin YearEnd BYearBegin BYearEnd'.split():
    cls = getattr(pd.tseries.offsets, name)
    st.register_type_strategy(cls, st.builds(
        cls,
        n=st.integers(-5, 5),
        normalize=st.booleans(),
        month=st.integers(min_value=1, max_value=12),
    ))

for name in 'QuarterBegin QuarterEnd BQuarterBegin BQuarterEnd'.split():
    cls = getattr(pd.tseries.offsets, name)
    st.register_type_strategy(cls, st.builds(
        cls,
        n=st.integers(-24, 24),
        normalize=st.booleans(),
        startingMonth=st.integers(min_value=1, max_value=12)
    ))
