from collections import abc
from datetime import date, time, timedelta, timezone
from decimal import Decimal
import operator
import os

from dateutil.tz import tzlocal, tzutc
import hypothesis
from hypothesis import strategies as st
import numpy as np
import pytest
from pytz import FixedOffset, utc

import pandas.util._test_decorators as td

import pandas as pd
from pandas import DataFrame
import pandas._testing as tm
from pandas.core import ops
from pandas.core.indexes.api import Index, MultiIndex

hypothesis.settings.register_profile(
    "ci",
    # Hypothesis timing checks are tuned for scalars by default, so we bump
    # them from 200ms to 500ms per test case as the global default.  If this
    # is too short for a specific test, (a) try to make it faster, and (b)
    # if it really is slow add `@settings(deadline=...)` with a working value,
    # or `deadline=None` to entirely disable timeouts for that test.
    deadline=500,
    suppress_health_check=(hypothesis.HealthCheck.too_slow,),
)
hypothesis.settings.load_profile("ci")


def pytest_addoption(parser):
    parser.addoption("--skip-slow", action="store_true", help="skip slow tests")
    parser.addoption("--skip-network", action="store_true", help="skip network tests")
    parser.addoption("--skip-db", action="store_true", help="skip db tests")
    parser.addoption(
        "--run-high-memory", action="store_true", help="run high memory tests"
    )
    parser.addoption("--only-slow", action="store_true", help="run only slow tests")
    parser.addoption(
        "--strict-data-files",
        action="store_true",
        help="Fail if a test is skipped for missing data file.",
    )


def pytest_runtest_setup(item):
    if "slow" in item.keywords and item.config.getoption("--skip-slow"):
        pytest.skip("skipping due to --skip-slow")

    if "slow" not in item.keywords and item.config.getoption("--only-slow"):
        pytest.skip("skipping due to --only-slow")

    if "network" in item.keywords and item.config.getoption("--skip-network"):
        pytest.skip("skipping due to --skip-network")

    if "db" in item.keywords and item.config.getoption("--skip-db"):
        pytest.skip("skipping due to --skip-db")

    if "high_memory" in item.keywords and not item.config.getoption(
        "--run-high-memory"
    ):
        pytest.skip("skipping high memory test since --run-high-memory was not set")


@pytest.fixture(autouse=True)
def configure_tests():
    """
    Configure settings for all tests and test modules.
    """
    pd.set_option("chained_assignment", "raise")


@pytest.fixture(autouse=True)
def add_imports(doctest_namespace):
    """
    Make `np` and `pd` names available for doctests.
    """
    doctest_namespace["np"] = np
    doctest_namespace["pd"] = pd


@pytest.fixture(params=["bsr", "coo", "csc", "csr", "dia", "dok", "lil"])
def spmatrix(request):
    """
    Yields scipy sparse matrix classes.
    """
    from scipy import sparse

    return getattr(sparse, request.param + "_matrix")


@pytest.fixture(params=[0, 1, "index", "columns"], ids=lambda x: f"axis {repr(x)}")
def axis(request):
    """
    Fixture for returning the axis numbers of a DataFrame.
    """
    return request.param


axis_frame = axis


@pytest.fixture(params=[0, "index"], ids=lambda x: f"axis {repr(x)}")
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
    pytest.importorskip("IPython", minversion="6.0.0")
    from IPython.core.interactiveshell import InteractiveShell

    return InteractiveShell()


@pytest.fixture(params=[True, False, None])
def observed(request):
    """
    Pass in the observed keyword to groupby for [True, False]
    This indicates whether categoricals should return values for
    values which are not in the grouper [False / None], or only values which
    appear in the grouper [True]. [None] is supported for future compatibility
    if we decide to change the default (and would need to warn if this
    parameter is not passed).
    """
    return request.param


@pytest.fixture(params=[True, False, None])
def ordered_fixture(request):
    """
    Boolean 'ordered' parameter for Categorical.
    """
    return request.param


_all_arithmetic_operators = [
    "__add__",
    "__radd__",
    "__sub__",
    "__rsub__",
    "__mul__",
    "__rmul__",
    "__floordiv__",
    "__rfloordiv__",
    "__truediv__",
    "__rtruediv__",
    "__pow__",
    "__rpow__",
    "__mod__",
    "__rmod__",
]


@pytest.fixture(params=_all_arithmetic_operators)
def all_arithmetic_operators(request):
    """
    Fixture for dunder names for common arithmetic operations.
    """
    return request.param


@pytest.fixture(
    params=[
        operator.add,
        ops.radd,
        operator.sub,
        ops.rsub,
        operator.mul,
        ops.rmul,
        operator.truediv,
        ops.rtruediv,
        operator.floordiv,
        ops.rfloordiv,
        operator.mod,
        ops.rmod,
        operator.pow,
        ops.rpow,
    ]
)
def all_arithmetic_functions(request):
    """
    Fixture for operator and roperator arithmetic functions.

    Notes
    -----
    This includes divmod and rdivmod, whereas all_arithmetic_operators
    does not.
    """
    return request.param


_all_numeric_reductions = [
    "sum",
    "max",
    "min",
    "mean",
    "prod",
    "std",
    "var",
    "median",
    "kurt",
    "skew",
]


@pytest.fixture(params=_all_numeric_reductions)
def all_numeric_reductions(request):
    """
    Fixture for numeric reduction names.
    """
    return request.param


_all_boolean_reductions = ["all", "any"]


@pytest.fixture(params=_all_boolean_reductions)
def all_boolean_reductions(request):
    """
    Fixture for boolean reduction names.
    """
    return request.param


_cython_table = pd.core.base.SelectionMixin._cython_table.items()


@pytest.fixture(params=list(_cython_table))
def cython_table_items(request):
    """
    Yields a tuple of a function and its corresponding name. Correspond to
    the list of aggregator "Cython functions" used on selected table items.
    """
    return request.param


def _get_cython_table_params(ndframe, func_names_and_expected):
    """
    Combine frame, functions from SelectionMixin._cython_table
    keys and expected result.

    Parameters
    ----------
    ndframe : DataFrame or Series
    func_names_and_expected : Sequence of two items
        The first item is a name of a NDFrame method ('sum', 'prod') etc.
        The second item is the expected return value.

    Returns
    -------
    list
        List of three items (DataFrame, function, expected result)
    """
    results = []
    for func_name, expected in func_names_and_expected:
        results.append((ndframe, func_name, expected))
        results += [
            (ndframe, func, expected)
            for func, name in _cython_table
            if name == func_name
        ]
    return results


@pytest.fixture(params=["__eq__", "__ne__", "__le__", "__lt__", "__ge__", "__gt__"])
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


@pytest.fixture(params=["__le__", "__lt__", "__ge__", "__gt__"])
def compare_operators_no_eq_ne(request):
    """
    Fixture for dunder names for compare operations except == and !=

    * >=
    * >
    * <
    * <=
    """
    return request.param


@pytest.fixture(
    params=["__and__", "__rand__", "__or__", "__ror__", "__xor__", "__rxor__"]
)
def all_logical_operators(request):
    """
    Fixture for dunder names for common logical operations

    * |
    * &
    * ^
    """
    return request.param


@pytest.fixture(params=[None, "gzip", "bz2", "zip", "xz"])
def compression(request):
    """
    Fixture for trying common compression types in compression tests.
    """
    return request.param


@pytest.fixture(params=["gzip", "bz2", "zip", "xz"])
def compression_only(request):
    """
    Fixture for trying common compression types in compression tests excluding
    uncompressed case.
    """
    return request.param


@pytest.fixture(params=[True, False])
def writable(request):
    """
    Fixture that an array is writable.
    """
    return request.param


@pytest.fixture(scope="module")
def datetime_tz_utc():
    """
    Yields the UTC timezone object from the datetime module.
    """
    return timezone.utc


@pytest.fixture(params=["utc", "dateutil/UTC", utc, tzutc(), timezone.utc])
def utc_fixture(request):
    """
    Fixture to provide variants of UTC timezone strings and tzinfo objects.
    """
    return request.param


@pytest.fixture(params=["inner", "outer", "left", "right"])
def join_type(request):
    """
    Fixture for trying all types of join operations.
    """
    return request.param


@pytest.fixture
def strict_data_files(pytestconfig):
    """
    Returns the configuration for the test setting `--strict-data-files`.
    """
    return pytestconfig.getoption("--strict-data-files")


@pytest.fixture
def datapath(strict_data_files):
    """
    Get the path to a data file.

    Parameters
    ----------
    path : str
        Path to the file, relative to ``pandas/tests/``

    Returns
    -------
    path including ``pandas/tests``.

    Raises
    ------
    ValueError
        If the path doesn't exist and the --strict-data-files option is set.
    """
    BASE_PATH = os.path.join(os.path.dirname(__file__), "tests")

    def deco(*args):
        path = os.path.join(BASE_PATH, *args)
        if not os.path.exists(path):
            if strict_data_files:
                raise ValueError(
                    f"Could not find file {path} and --strict-data-files is set."
                )
            else:
                pytest.skip(f"Could not find {path}.")
        return path

    return deco


@pytest.fixture
def iris(datapath):
    """
    The iris dataset as a DataFrame.
    """
    return pd.read_csv(datapath("data", "iris.csv"))


@pytest.fixture(params=["nlargest", "nsmallest"])
def nselect_method(request):
    """
    Fixture for trying all nselect methods.
    """
    return request.param


@pytest.fixture(params=["first", "last", False])
def keep(request):
    """
    Valid values for the 'keep' parameter used in
    .duplicated or .drop_duplicates
    """
    return request.param


@pytest.fixture(params=["left", "right", "both", "neither"])
def closed(request):
    """
    Fixture for trying all interval closed parameters.
    """
    return request.param


@pytest.fixture(params=["left", "right", "both", "neither"])
def other_closed(request):
    """
    Secondary closed fixture to allow parametrizing over all pairs of closed.
    """
    return request.param


@pytest.fixture(params=[None, np.nan, pd.NaT, float("nan"), np.float("NaN"), pd.NA])
def nulls_fixture(request):
    """
    Fixture for each null type in pandas.
    """
    return request.param


nulls_fixture2 = nulls_fixture  # Generate cartesian product of nulls_fixture


@pytest.fixture(params=[None, np.nan, pd.NaT])
def unique_nulls_fixture(request):
    """
    Fixture for each null type in pandas, each null type exactly once.
    """
    return request.param


# Generate cartesian product of unique_nulls_fixture:
unique_nulls_fixture2 = unique_nulls_fixture


TIMEZONES = [
    None,
    "UTC",
    "US/Eastern",
    "Asia/Tokyo",
    "dateutil/US/Pacific",
    "dateutil/Asia/Singapore",
    tzutc(),
    tzlocal(),
    FixedOffset(300),
    FixedOffset(0),
    FixedOffset(-300),
    timezone.utc,
    timezone(timedelta(hours=1)),
    timezone(timedelta(hours=-1), name="foo"),
]
TIMEZONE_IDS = [repr(i) for i in TIMEZONES]


@td.parametrize_fixture_doc(str(TIMEZONE_IDS))
@pytest.fixture(params=TIMEZONES, ids=TIMEZONE_IDS)
def tz_naive_fixture(request):
    """
    Fixture for trying timezones including default (None): {0}
    """
    return request.param


@td.parametrize_fixture_doc(str(TIMEZONE_IDS[1:]))
@pytest.fixture(params=TIMEZONES[1:], ids=TIMEZONE_IDS[1:])
def tz_aware_fixture(request):
    """
    Fixture for trying explicit timezones: {0}
    """
    return request.param


# Generate cartesian product of tz_aware_fixture:
tz_aware_fixture2 = tz_aware_fixture


# ----------------------------------------------------------------
# Dtypes
# ----------------------------------------------------------------

UNSIGNED_INT_DTYPES = ["uint8", "uint16", "uint32", "uint64"]
UNSIGNED_EA_INT_DTYPES = ["UInt8", "UInt16", "UInt32", "UInt64"]
SIGNED_INT_DTYPES = [int, "int8", "int16", "int32", "int64"]
SIGNED_EA_INT_DTYPES = ["Int8", "Int16", "Int32", "Int64"]
ALL_INT_DTYPES = UNSIGNED_INT_DTYPES + SIGNED_INT_DTYPES
ALL_EA_INT_DTYPES = UNSIGNED_EA_INT_DTYPES + SIGNED_EA_INT_DTYPES

FLOAT_DTYPES = [float, "float32", "float64"]
COMPLEX_DTYPES = [complex, "complex64", "complex128"]
STRING_DTYPES = [str, "str", "U"]

DATETIME64_DTYPES = ["datetime64[ns]", "M8[ns]"]
TIMEDELTA64_DTYPES = ["timedelta64[ns]", "m8[ns]"]

BOOL_DTYPES = [bool, "bool"]
BYTES_DTYPES = [bytes, "bytes"]
OBJECT_DTYPES = [object, "object"]

ALL_REAL_DTYPES = FLOAT_DTYPES + ALL_INT_DTYPES
ALL_NUMPY_DTYPES = (
    ALL_REAL_DTYPES
    + COMPLEX_DTYPES
    + STRING_DTYPES
    + DATETIME64_DTYPES
    + TIMEDELTA64_DTYPES
    + BOOL_DTYPES
    + OBJECT_DTYPES
    + BYTES_DTYPES
)


@pytest.fixture(params=STRING_DTYPES)
def string_dtype(request):
    """
    Parametrized fixture for string dtypes.

    * str
    * 'str'
    * 'U'
    """
    return request.param


@pytest.fixture(params=BYTES_DTYPES)
def bytes_dtype(request):
    """
    Parametrized fixture for bytes dtypes.

    * bytes
    * 'bytes'
    """
    return request.param


@pytest.fixture(params=OBJECT_DTYPES)
def object_dtype(request):
    """
    Parametrized fixture for object dtypes.

    * object
    * 'object'
    """
    return request.param


@pytest.fixture(params=DATETIME64_DTYPES)
def datetime64_dtype(request):
    """
    Parametrized fixture for datetime64 dtypes.

    * 'datetime64[ns]'
    * 'M8[ns]'
    """
    return request.param


@pytest.fixture(params=TIMEDELTA64_DTYPES)
def timedelta64_dtype(request):
    """
    Parametrized fixture for timedelta64 dtypes.

    * 'timedelta64[ns]'
    * 'm8[ns]'
    """
    return request.param


@pytest.fixture(params=FLOAT_DTYPES)
def float_dtype(request):
    """
    Parameterized fixture for float dtypes.

    * float
    * 'float32'
    * 'float64'
    """
    return request.param


@pytest.fixture(params=COMPLEX_DTYPES)
def complex_dtype(request):
    """
    Parameterized fixture for complex dtypes.

    * complex
    * 'complex64'
    * 'complex128'
    """
    return request.param


@pytest.fixture(params=SIGNED_INT_DTYPES)
def sint_dtype(request):
    """
    Parameterized fixture for signed integer dtypes.

    * int
    * 'int8'
    * 'int16'
    * 'int32'
    * 'int64'
    """
    return request.param


@pytest.fixture(params=UNSIGNED_INT_DTYPES)
def uint_dtype(request):
    """
    Parameterized fixture for unsigned integer dtypes.

    * 'uint8'
    * 'uint16'
    * 'uint32'
    * 'uint64'
    """
    return request.param


@pytest.fixture(params=ALL_INT_DTYPES)
def any_int_dtype(request):
    """
    Parameterized fixture for any integer dtype.

    * int
    * 'int8'
    * 'uint8'
    * 'int16'
    * 'uint16'
    * 'int32'
    * 'uint32'
    * 'int64'
    * 'uint64'
    """
    return request.param


@pytest.fixture(params=ALL_EA_INT_DTYPES)
def any_nullable_int_dtype(request):
    """
    Parameterized fixture for any nullable integer dtype.

    * 'UInt8'
    * 'Int8'
    * 'UInt16'
    * 'Int16'
    * 'UInt32'
    * 'Int32'
    * 'UInt64'
    * 'Int64'
    """
    return request.param


@pytest.fixture(params=ALL_REAL_DTYPES)
def any_real_dtype(request):
    """
    Parameterized fixture for any (purely) real numeric dtype.

    * int
    * 'int8'
    * 'uint8'
    * 'int16'
    * 'uint16'
    * 'int32'
    * 'uint32'
    * 'int64'
    * 'uint64'
    * float
    * 'float32'
    * 'float64'
    """
    return request.param


@pytest.fixture(params=ALL_NUMPY_DTYPES)
def any_numpy_dtype(request):
    """
    Parameterized fixture for all numpy dtypes.

    * bool
    * 'bool'
    * int
    * 'int8'
    * 'uint8'
    * 'int16'
    * 'uint16'
    * 'int32'
    * 'uint32'
    * 'int64'
    * 'uint64'
    * float
    * 'float32'
    * 'float64'
    * complex
    * 'complex64'
    * 'complex128'
    * str
    * 'str'
    * 'U'
    * bytes
    * 'bytes'
    * 'datetime64[ns]'
    * 'M8[ns]'
    * 'timedelta64[ns]'
    * 'm8[ns]'
    * object
    * 'object'
    """
    return request.param


# categoricals are handled separately
_any_skipna_inferred_dtype = [
    ("string", ["a", np.nan, "c"]),
    ("string", ["a", pd.NA, "c"]),
    ("bytes", [b"a", np.nan, b"c"]),
    ("empty", [np.nan, np.nan, np.nan]),
    ("empty", []),
    ("mixed-integer", ["a", np.nan, 2]),
    ("mixed", ["a", np.nan, 2.0]),
    ("floating", [1.0, np.nan, 2.0]),
    ("integer", [1, np.nan, 2]),
    ("mixed-integer-float", [1, np.nan, 2.0]),
    ("decimal", [Decimal(1), np.nan, Decimal(2)]),
    ("boolean", [True, np.nan, False]),
    ("boolean", [True, pd.NA, False]),
    ("datetime64", [np.datetime64("2013-01-01"), np.nan, np.datetime64("2018-01-01")]),
    ("datetime", [pd.Timestamp("20130101"), np.nan, pd.Timestamp("20180101")]),
    ("date", [date(2013, 1, 1), np.nan, date(2018, 1, 1)]),
    # The following two dtypes are commented out due to GH 23554
    # ('complex', [1 + 1j, np.nan, 2 + 2j]),
    # ('timedelta64', [np.timedelta64(1, 'D'),
    #                  np.nan, np.timedelta64(2, 'D')]),
    ("timedelta", [timedelta(1), np.nan, timedelta(2)]),
    ("time", [time(1), np.nan, time(2)]),
    ("period", [pd.Period(2013), pd.NaT, pd.Period(2018)]),
    ("interval", [pd.Interval(0, 1), np.nan, pd.Interval(0, 2)]),
]
ids, _ = zip(*_any_skipna_inferred_dtype)  # use inferred type as fixture-id


@pytest.fixture(params=_any_skipna_inferred_dtype, ids=ids)
def any_skipna_inferred_dtype(request):
    """
    Fixture for all inferred dtypes from _libs.lib.infer_dtype

    The covered (inferred) types are:
    * 'string'
    * 'empty'
    * 'bytes'
    * 'mixed'
    * 'mixed-integer'
    * 'mixed-integer-float'
    * 'floating'
    * 'integer'
    * 'decimal'
    * 'boolean'
    * 'datetime64'
    * 'datetime'
    * 'date'
    * 'timedelta'
    * 'time'
    * 'period'
    * 'interval'

    Returns
    -------
    inferred_dtype : str
        The string for the inferred dtype from _libs.lib.infer_dtype
    values : np.ndarray
        An array of object dtype that will be inferred to have
        `inferred_dtype`

    Examples
    --------
    >>> import pandas._libs.lib as lib
    >>>
    >>> def test_something(any_skipna_inferred_dtype):
    ...     inferred_dtype, values = any_skipna_inferred_dtype
    ...     # will pass
    ...     assert lib.infer_dtype(values, skipna=True) == inferred_dtype
    """
    inferred_dtype, values = request.param
    values = np.array(values, dtype=object)  # object dtype to avoid casting

    # correctness of inference tested in tests/dtypes/test_inference.py
    return inferred_dtype, values


@pytest.fixture(
    params=[
        getattr(pd.offsets, o)
        for o in pd.offsets.__all__
        if issubclass(getattr(pd.offsets, o), pd.offsets.Tick)
    ]
)
def tick_classes(request):
    """
    Fixture for Tick based datetime offsets available for a time series.
    """
    return request.param


# ----------------------------------------------------------------
# Global setup for tests using Hypothesis


# Registering these strategies makes them globally available via st.from_type,
# which is use for offsets in tests/tseries/offsets/test_offsets_properties.py
for name in "MonthBegin MonthEnd BMonthBegin BMonthEnd".split():
    cls = getattr(pd.tseries.offsets, name)
    st.register_type_strategy(
        cls, st.builds(cls, n=st.integers(-99, 99), normalize=st.booleans())
    )

for name in "YearBegin YearEnd BYearBegin BYearEnd".split():
    cls = getattr(pd.tseries.offsets, name)
    st.register_type_strategy(
        cls,
        st.builds(
            cls,
            n=st.integers(-5, 5),
            normalize=st.booleans(),
            month=st.integers(min_value=1, max_value=12),
        ),
    )

for name in "QuarterBegin QuarterEnd BQuarterBegin BQuarterEnd".split():
    cls = getattr(pd.tseries.offsets, name)
    st.register_type_strategy(
        cls,
        st.builds(
            cls,
            n=st.integers(-24, 24),
            normalize=st.booleans(),
            startingMonth=st.integers(min_value=1, max_value=12),
        ),
    )


@pytest.fixture
def empty_series():
    return pd.Series([], index=[], dtype=np.float64)


@pytest.fixture
def datetime_series():
    """
    Fixture for Series of floats with DatetimeIndex
    """
    s = tm.makeTimeSeries()
    s.name = "ts"
    return s


@pytest.fixture
def string_series():
    """
    Fixture for Series of floats with Index of unique strings
    """
    s = tm.makeStringSeries()
    s.name = "series"
    return s


@pytest.fixture
def object_series():
    """
    Fixture for Series of dtype object with Index of unique strings
    """
    s = tm.makeObjectSeries()
    s.name = "objects"
    return s


@pytest.fixture
def float_frame():
    """
    Fixture for DataFrame of floats with index of unique strings

    Columns are ['A', 'B', 'C', 'D'].

                       A         B         C         D
    P7GACiRnxd -0.465578 -0.361863  0.886172 -0.053465
    qZKh6afn8n -0.466693 -0.373773  0.266873  1.673901
    tkp0r6Qble  0.148691 -0.059051  0.174817  1.598433
    wP70WOCtv8  0.133045 -0.581994 -0.992240  0.261651
    M2AeYQMnCz -1.207959 -0.185775  0.588206  0.563938
    QEPzyGDYDo -0.381843 -0.758281  0.502575 -0.565053
    r78Jwns6dn -0.653707  0.883127  0.682199  0.206159
    ...              ...       ...       ...       ...
    IHEGx9NO0T -0.277360  0.113021 -1.018314  0.196316
    lPMj8K27FA -1.313667 -0.604776 -1.305618 -0.863999
    qa66YMWQa5  1.110525  0.475310 -0.747865  0.032121
    yOa0ATsmcE -0.431457  0.067094  0.096567 -0.264962
    65znX3uRNG  1.528446  0.160416 -0.109635 -0.032987
    eCOBvKqf3e  0.235281  1.622222  0.781255  0.392871
    xSucinXxuV -1.263557  0.252799 -0.552247  0.400426

    [30 rows x 4 columns]
    """
    return DataFrame(tm.getSeriesData())


@pytest.fixture(params=[pd.Index, pd.Series], ids=["index", "series"])
def index_or_series(request):
    """
    Fixture to parametrize over Index and Series, made necessary by a mypy
    bug, giving an error:

    List item 0 has incompatible type "Type[Series]"; expected "Type[PandasObject]"

    See GH#29725
    """
    return request.param


@pytest.fixture
def dict_subclass():
    """
    Fixture for a dictionary subclass.
    """

    class TestSubDict(dict):
        def __init__(self, *args, **kwargs):
            dict.__init__(self, *args, **kwargs)

    return TestSubDict


@pytest.fixture
def non_mapping_dict_subclass():
    """
    Fixture for a non-mapping dictionary subclass.
    """

    class TestNonDictMapping(abc.Mapping):
        def __init__(self, underlying_dict):
            self._data = underlying_dict

        def __getitem__(self, key):
            return self._data.__getitem__(key)

        def __iter__(self):
            return self._data.__iter__()

        def __len__(self):
            return self._data.__len__()

    return TestNonDictMapping


def _gen_mi():
    # a MultiIndex used to test the general functionality of this object

    # See Also: tests.multi.conftest.idx
    major_axis = Index(["foo", "bar", "baz", "qux"])
    minor_axis = Index(["one", "two"])

    major_codes = np.array([0, 0, 1, 2, 3, 3])
    minor_codes = np.array([0, 1, 0, 1, 0, 1])
    index_names = ["first", "second"]
    mi = MultiIndex(
        levels=[major_axis, minor_axis],
        codes=[major_codes, minor_codes],
        names=index_names,
        verify_integrity=False,
    )
    return mi


indices_dict = {
    "unicode": tm.makeUnicodeIndex(100),
    "string": tm.makeStringIndex(100),
    "datetime": tm.makeDateIndex(100),
    "datetime-tz": tm.makeDateIndex(100, tz="US/Pacific"),
    "period": tm.makePeriodIndex(100),
    "timedelta": tm.makeTimedeltaIndex(100),
    "int": tm.makeIntIndex(100),
    "uint": tm.makeUIntIndex(100),
    "range": tm.makeRangeIndex(100),
    "float": tm.makeFloatIndex(100),
    "bool": tm.makeBoolIndex(10),
    "categorical": tm.makeCategoricalIndex(100),
    "interval": tm.makeIntervalIndex(100),
    "empty": Index([]),
    "tuples": MultiIndex.from_tuples(zip(["foo", "bar", "baz"], [1, 2, 3])),
    "multi": _gen_mi(),
    "repeats": Index([0, 0, 1, 1, 2, 2]),
}


@pytest.fixture(params=indices_dict.keys())
def indices(request):
    """
    Fixture for many "simple" kinds of indices.

    These indices are unlikely to cover corner cases, e.g.
        - no names
        - no NaTs/NaNs
        - no values near implementation bounds
        - ...
    """
    # copy to avoid mutation, e.g. setting .name
    return indices_dict[request.param].copy()


def _create_series(index):
    """ Helper for the _series dict """
    size = len(index)
    data = np.random.randn(size)
    return pd.Series(data, index=index, name="a")


_series = {
    f"series-with-{index_id}-index": _create_series(index)
    for index_id, index in indices_dict.items()
}


@pytest.fixture
def series_with_simple_index(indices):
    """
    Fixture for tests on series with changing types of indices.
    """
    return _create_series(indices)


_narrow_dtypes = [
    np.float16,
    np.float32,
    np.int8,
    np.int16,
    np.int32,
    np.uint8,
    np.uint16,
    np.uint32,
]
_narrow_series = {
    f"{dtype.__name__}-series": tm.makeFloatSeries(name="a").astype(dtype)
    for dtype in _narrow_dtypes
}


@pytest.fixture(params=_narrow_series.keys())
def narrow_series(request):
    """
    Fixture for Series with low precision data types
    """
    # copy to avoid mutation, e.g. setting .name
    return _narrow_series[request.param].copy()


_index_or_series_objs = {**indices_dict, **_series, **_narrow_series}


@pytest.fixture(params=_index_or_series_objs.keys())
def index_or_series_obj(request):
    """
    Fixture for tests on indexes, series and series with a narrow dtype
    copy to avoid mutation, e.g. setting .name
    """
    return _index_or_series_objs[request.param].copy(deep=True)


@pytest.fixture
def multiindex_year_month_day_dataframe_random_data():
    """
    DataFrame with 3 level MultiIndex (year, month, day) covering
    first 100 business days from 2000-01-01 with random data
    """
    tdf = tm.makeTimeDataFrame(100)
    ymd = tdf.groupby([lambda x: x.year, lambda x: x.month, lambda x: x.day]).sum()
    # use Int64Index, to make sure things work
    ymd.index.set_levels([lev.astype("i8") for lev in ymd.index.levels], inplace=True)
    ymd.index.set_names(["year", "month", "day"], inplace=True)
    return ymd
