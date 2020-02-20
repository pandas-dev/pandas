from functools import wraps
import os
from typing import Callable, List, Type
import warnings

import numpy as np
from numpy.random import rand

from pandas._config.localization import (  # noqa:F401
    can_set_locale,
    get_locales,
    set_locale,
)

from pandas.compat import _get_lzma_file, _import_lzma

from pandas.core.dtypes.common import (
    is_datetime64_dtype,
    is_datetime64tz_dtype,
    is_period_dtype,
    is_timedelta64_dtype,
)

import pandas as pd
from pandas import Categorical, DataFrame, Series
from pandas.core.arrays import DatetimeArray, PeriodArray, TimedeltaArray, period_array

from .asserters import (  # noqa:F401
    assert_almost_equal,
    assert_attr_equal,
    assert_categorical_equal,
    assert_class_equal,
    assert_contains_all,
    assert_copy,
    assert_datetime_array_equal,
    assert_dict_equal,
    assert_equal,
    assert_extension_array_equal,
    assert_frame_equal,
    assert_index_equal,
    assert_interval_array_equal,
    assert_is_sorted,
    assert_numpy_array_equal,
    assert_period_array_equal,
    assert_series_equal,
    assert_sp_array_equal,
    assert_timedelta_array_equal,
    raise_assert_detail,
)
from .contexts import (  # noqa:F401
    RNGContext,
    assert_produces_warning,
    decompress_file,
    ensure_clean,
    ensure_clean_dir,
    ensure_safe_environment_variables,
    set_timezone,
    use_numexpr,
    with_csv_dialect,
)
from .graphics import assert_is_valid_plot_return_object, close  # noqa:F401
from .makers import (  # noqa:F401
    _create_missing_idx,
    _make_timeseries,
    getMixedTypeDict,
    getPeriodData,
    getSeriesData,
    getTimeSeriesData,
    makeBoolIndex,
    makeCategoricalIndex,
    makeCustomDataframe,
    makeCustomIndex,
    makeDataFrame,
    makeDateIndex,
    makeFloatIndex,
    makeFloatSeries,
    makeIntervalIndex,
    makeIntIndex,
    makeMissingCustomDataframe,
    makeMissingDataframe,
    makeMixedDataFrame,
    makeMultiIndex,
    makeObjectSeries,
    makePeriodFrame,
    makePeriodIndex,
    makePeriodSeries,
    makeRangeIndex,
    makeStringIndex,
    makeStringSeries,
    makeTimeDataFrame,
    makeTimedeltaIndex,
    makeTimeSeries,
    makeUIntIndex,
    makeUnicodeIndex,
    rands,
    rands_array,
    randu,
    randu_array,
)
from .network import (  # noqa:F401
    can_connect,
    network,
    optional_args,
    with_connectivity_check,
)
from .roundtrips import (  # noqa:F401
    round_trip_localpath,
    round_trip_pathlib,
    round_trip_pickle,
)

lzma = _import_lzma()


# set testing_mode
_testing_mode_warnings = (DeprecationWarning, ResourceWarning)


def set_testing_mode():
    # set the testing mode filters
    testing_mode = os.environ.get("PANDAS_TESTING_MODE", "None")
    if "deprecate" in testing_mode:
        warnings.simplefilter("always", _testing_mode_warnings)


def reset_testing_mode():
    # reset the testing mode filters
    testing_mode = os.environ.get("PANDAS_TESTING_MODE", "None")
    if "deprecate" in testing_mode:
        warnings.simplefilter("ignore", _testing_mode_warnings)


set_testing_mode()


def reset_display_options():
    """
    Reset the display options for printing and representing objects.
    """
    pd.reset_option("^display.", silent=True)


def write_to_compressed(compression, path, data, dest="test"):
    """
    Write data to a compressed file.

    Parameters
    ----------
    compression : {'gzip', 'bz2', 'zip', 'xz'}
        The compression type to use.
    path : str
        The file path to write the data.
    data : str
        The data to write.
    dest : str, default "test"
        The destination file (for ZIP only)

    Raises
    ------
    ValueError : An invalid compression value was passed in.
    """
    if compression == "zip":
        import zipfile

        compress_method = zipfile.ZipFile
    elif compression == "gzip":
        import gzip

        compress_method = gzip.GzipFile
    elif compression == "bz2":
        import bz2

        compress_method = bz2.BZ2File
    elif compression == "xz":
        compress_method = _get_lzma_file(lzma)
    else:
        raise ValueError(f"Unrecognized compression type: {compression}")

    if compression == "zip":
        mode = "w"
        args = (dest, data)
        method = "writestr"
    else:
        mode = "wb"
        args = (data,)
        method = "write"

    with compress_method(path, mode=mode) as f:
        getattr(f, method)(*args)


def randbool(size=(), p: float = 0.5):
    return rand(*size) <= p


# -----------------------------------------------------------------------------
# Comparators


def equalContents(arr1, arr2) -> bool:
    """
    Checks if the set of unique elements of arr1 and arr2 are equivalent.
    """
    return frozenset(arr1) == frozenset(arr2)


def isiterable(obj):
    return hasattr(obj, "__iter__")


def box_expected(expected, box_cls, transpose=True):
    """
    Helper function to wrap the expected output of a test in a given box_class.

    Parameters
    ----------
    expected : np.ndarray, Index, Series
    box_cls : {Index, Series, DataFrame}

    Returns
    -------
    subclass of box_cls
    """
    if box_cls is pd.Index:
        expected = pd.Index(expected)
    elif box_cls is pd.Series:
        expected = pd.Series(expected)
    elif box_cls is pd.DataFrame:
        expected = pd.Series(expected).to_frame()
        if transpose:
            # for vector operations, we we need a DataFrame to be a single-row,
            #  not a single-column, in order to operate against non-DataFrame
            #  vectors of the same length.
            expected = expected.T
    elif box_cls is PeriodArray:
        # the PeriodArray constructor is not as flexible as period_array
        expected = period_array(expected)
    elif box_cls is DatetimeArray:
        expected = DatetimeArray(expected)
    elif box_cls is TimedeltaArray:
        expected = TimedeltaArray(expected)
    elif box_cls is np.ndarray:
        expected = np.array(expected)
    elif box_cls is to_array:
        expected = to_array(expected)
    else:
        raise NotImplementedError(box_cls)
    return expected


def to_array(obj):
    # temporary implementation until we get pd.array in place
    if is_period_dtype(obj):
        return period_array(obj)
    elif is_datetime64_dtype(obj) or is_datetime64tz_dtype(obj):
        return DatetimeArray._from_sequence(obj)
    elif is_timedelta64_dtype(obj):
        return TimedeltaArray._from_sequence(obj)
    else:
        return np.array(obj)


# -----------------------------------------------------------------------------
# Others


def all_index_generator(k=10):
    """
    Generator which can be iterated over to get instances of all the various
    index classes.

    Parameters
    ----------
    k: length of each of the index instances
    """
    all_make_index_funcs = [
        makeIntIndex,
        makeFloatIndex,
        makeStringIndex,
        makeUnicodeIndex,
        makeDateIndex,
        makePeriodIndex,
        makeTimedeltaIndex,
        makeBoolIndex,
        makeRangeIndex,
        makeIntervalIndex,
        makeCategoricalIndex,
    ]
    for make_index_func in all_make_index_funcs:
        yield make_index_func(k=k)


def index_subclass_makers_generator():
    make_index_funcs = [
        makeDateIndex,
        makePeriodIndex,
        makeTimedeltaIndex,
        makeRangeIndex,
        makeIntervalIndex,
        makeCategoricalIndex,
        makeMultiIndex,
    ]
    for make_index_func in make_index_funcs:
        yield make_index_func


def all_timeseries_index_generator(k=10):
    """
    Generator which can be iterated over to get instances of all the classes
    which represent time-series.

    Parameters
    ----------
    k: length of each of the index instances
    """
    make_index_funcs = [makeDateIndex, makePeriodIndex, makeTimedeltaIndex]
    for make_index_func in make_index_funcs:
        yield make_index_func(k=k)


def test_parallel(num_threads=2, kwargs_list=None):
    """
    Decorator to run the same function multiple times in parallel.

    Parameters
    ----------
    num_threads : int, optional
        The number of times the function is run in parallel.
    kwargs_list : list of dicts, optional
        The list of kwargs to update original
        function kwargs on different threads.

    Notes
    -----
    This decorator does not pass the return value of the decorated function.

    Original from scikit-image:

    https://github.com/scikit-image/scikit-image/pull/1519

    """
    assert num_threads > 0
    has_kwargs_list = kwargs_list is not None
    if has_kwargs_list:
        assert len(kwargs_list) == num_threads
    import threading

    def wrapper(func):
        @wraps(func)
        def inner(*args, **kwargs):
            if has_kwargs_list:
                update_kwargs = lambda i: dict(kwargs, **kwargs_list[i])
            else:
                update_kwargs = lambda i: kwargs
            threads = []
            for i in range(num_threads):
                updated_kwargs = update_kwargs(i)
                thread = threading.Thread(target=func, args=args, kwargs=updated_kwargs)
                threads.append(thread)
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join()

        return inner

    return wrapper


class SubclassedSeries(Series):
    _metadata = ["testattr", "name"]

    @property
    def _constructor(self):
        return SubclassedSeries

    @property
    def _constructor_expanddim(self):
        return SubclassedDataFrame


class SubclassedDataFrame(DataFrame):
    _metadata = ["testattr"]

    @property
    def _constructor(self):
        return SubclassedDataFrame

    @property
    def _constructor_sliced(self):
        return SubclassedSeries


class SubclassedCategorical(Categorical):
    @property
    def _constructor(self):
        return SubclassedCategorical


def _make_skipna_wrapper(alternative, skipna_alternative=None):
    """
    Create a function for calling on an array.

    Parameters
    ----------
    alternative : function
        The function to be called on the array with no NaNs.
        Only used when 'skipna_alternative' is None.
    skipna_alternative : function
        The function to be called on the original array

    Returns
    -------
    function
    """
    if skipna_alternative:

        def skipna_wrapper(x):
            return skipna_alternative(x.values)

    else:

        def skipna_wrapper(x):
            nona = x.dropna()
            if len(nona) == 0:
                return np.nan
            return alternative(nona)

    return skipna_wrapper


def convert_rows_list_to_csv_str(rows_list: List[str]):
    """
    Convert list of CSV rows to single CSV-formatted string for current OS.

    This method is used for creating expected value of to_csv() method.

    Parameters
    ----------
    rows_list : List[str]
        Each element represents the row of csv.

    Returns
    -------
    str
        Expected output of to_csv() in current OS.
    """
    sep = os.linesep
    expected = sep.join(rows_list) + sep
    return expected


def external_error_raised(
    expected_exception: Type[Exception],
) -> Callable[[Type[Exception], None], None]:
    """
    Helper function to mark pytest.raises that have an external error message.

    Parameters
    ----------
    expected_exception : Exception
        Expected error to raise.

    Returns
    -------
    Callable
        Regular `pytest.raises` function with `match` equal to `None`.
    """
    import pytest

    return pytest.raises(expected_exception, match=None)
