import bz2
from collections import Counter
from contextlib import contextmanager
from datetime import datetime
from functools import wraps
import gzip
import operator
import os
import re
import string
from typing import (
    Any,
    Callable,
    ContextManager,
    List,
    Optional,
    Sequence,
    Tuple,
    Type,
    Union,
    cast,
)
import warnings
import zipfile

import numpy as np

from pandas._config.localization import (  # noqa:F401
    can_set_locale,
    get_locales,
    set_locale,
)

from pandas._typing import Dtype, FilePathOrBuffer, FrameOrSeries
from pandas.compat import get_lzma_file, import_lzma

from pandas.core.dtypes.common import (
    is_datetime64_dtype,
    is_datetime64tz_dtype,
    is_period_dtype,
    is_sequence,
    is_timedelta64_dtype,
)

import pandas as pd
from pandas import (
    Categorical,
    CategoricalIndex,
    DataFrame,
    DatetimeIndex,
    Index,
    IntervalIndex,
    MultiIndex,
    RangeIndex,
    Series,
    bdate_range,
)
from pandas._testing.asserters import (  # noqa:F401
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
    assert_is_valid_plot_return_object,
    assert_numpy_array_equal,
    assert_period_array_equal,
    assert_series_equal,
    assert_sp_array_equal,
    assert_timedelta_array_equal,
    raise_assert_detail,
)
from pandas._testing.contexts import (  # noqa:F401
    decompress_file,
    ensure_clean,
    ensure_clean_dir,
    ensure_safe_environment_variables,
    set_timezone,
    use_numexpr,
    with_csv_dialect,
)
from pandas.core.arrays import DatetimeArray, PeriodArray, TimedeltaArray, period_array

from pandas.io.common import urlopen

lzma = import_lzma()

_N = 30
_K = 4
_RAISE_NETWORK_ERROR_DEFAULT = False

UNSIGNED_INT_DTYPES: List[Dtype] = ["uint8", "uint16", "uint32", "uint64"]
UNSIGNED_EA_INT_DTYPES: List[Dtype] = ["UInt8", "UInt16", "UInt32", "UInt64"]
SIGNED_INT_DTYPES: List[Dtype] = [int, "int8", "int16", "int32", "int64"]
SIGNED_EA_INT_DTYPES: List[Dtype] = ["Int8", "Int16", "Int32", "Int64"]
ALL_INT_DTYPES = UNSIGNED_INT_DTYPES + SIGNED_INT_DTYPES
ALL_EA_INT_DTYPES = UNSIGNED_EA_INT_DTYPES + SIGNED_EA_INT_DTYPES

FLOAT_DTYPES: List[Dtype] = [float, "float32", "float64"]
FLOAT_EA_DTYPES: List[Dtype] = ["Float32", "Float64"]
COMPLEX_DTYPES: List[Dtype] = [complex, "complex64", "complex128"]
STRING_DTYPES: List[Dtype] = [str, "str", "U"]

DATETIME64_DTYPES: List[Dtype] = ["datetime64[ns]", "M8[ns]"]
TIMEDELTA64_DTYPES: List[Dtype] = ["timedelta64[ns]", "m8[ns]"]

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

NULL_OBJECTS = [None, np.nan, pd.NaT, float("nan"), pd.NA]


# set testing_mode
_testing_mode_warnings = (DeprecationWarning, ResourceWarning)


def set_testing_mode():
    # set the testing mode filters
    testing_mode = os.environ.get("PANDAS_TESTING_MODE", "None")
    if "deprecate" in testing_mode:
        # pandas\_testing.py:119: error: Argument 2 to "simplefilter" has
        # incompatible type "Tuple[Type[DeprecationWarning],
        # Type[ResourceWarning]]"; expected "Type[Warning]"
        warnings.simplefilter(
            "always", _testing_mode_warnings  # type: ignore[arg-type]
        )


def reset_testing_mode():
    # reset the testing mode filters
    testing_mode = os.environ.get("PANDAS_TESTING_MODE", "None")
    if "deprecate" in testing_mode:
        # pandas\_testing.py:126: error: Argument 2 to "simplefilter" has
        # incompatible type "Tuple[Type[DeprecationWarning],
        # Type[ResourceWarning]]"; expected "Type[Warning]"
        warnings.simplefilter(
            "ignore", _testing_mode_warnings  # type: ignore[arg-type]
        )


set_testing_mode()


def reset_display_options():
    """
    Reset the display options for printing and representing objects.
    """
    pd.reset_option("^display.", silent=True)


def round_trip_pickle(
    obj: Any, path: Optional[FilePathOrBuffer] = None
) -> FrameOrSeries:
    """
    Pickle an object and then read it again.

    Parameters
    ----------
    obj : any object
        The object to pickle and then re-read.
    path : str, path object or file-like object, default None
        The path where the pickled object is written and then read.

    Returns
    -------
    pandas object
        The original object that was pickled and then re-read.
    """
    _path = path
    if _path is None:
        _path = f"__{rands(10)}__.pickle"
    with ensure_clean(_path) as temp_path:
        pd.to_pickle(obj, temp_path)
        return pd.read_pickle(temp_path)


def round_trip_pathlib(writer, reader, path: Optional[str] = None):
    """
    Write an object to file specified by a pathlib.Path and read it back

    Parameters
    ----------
    writer : callable bound to pandas object
        IO writing function (e.g. DataFrame.to_csv )
    reader : callable
        IO reading function (e.g. pd.read_csv )
    path : str, default None
        The path where the object is written and then read.

    Returns
    -------
    pandas object
        The original object that was serialized and then re-read.
    """
    import pytest

    Path = pytest.importorskip("pathlib").Path
    if path is None:
        path = "___pathlib___"
    with ensure_clean(path) as path:
        writer(Path(path))
        obj = reader(Path(path))
    return obj


def round_trip_localpath(writer, reader, path: Optional[str] = None):
    """
    Write an object to file specified by a py.path LocalPath and read it back.

    Parameters
    ----------
    writer : callable bound to pandas object
        IO writing function (e.g. DataFrame.to_csv )
    reader : callable
        IO reading function (e.g. pd.read_csv )
    path : str, default None
        The path where the object is written and then read.

    Returns
    -------
    pandas object
        The original object that was serialized and then re-read.
    """
    import pytest

    LocalPath = pytest.importorskip("py.path").local
    if path is None:
        path = "___localpath___"
    with ensure_clean(path) as path:
        writer(LocalPath(path))
        obj = reader(LocalPath(path))
    return obj


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
    args: Tuple[Any, ...] = (data,)
    mode = "wb"
    method = "write"
    compress_method: Callable

    if compression == "zip":
        compress_method = zipfile.ZipFile
        mode = "w"
        args = (dest, data)
        method = "writestr"
    elif compression == "gzip":
        compress_method = gzip.GzipFile
    elif compression == "bz2":
        compress_method = bz2.BZ2File
    elif compression == "xz":
        compress_method = get_lzma_file(lzma)
    else:
        raise ValueError(f"Unrecognized compression type: {compression}")

    with compress_method(path, mode=mode) as f:
        getattr(f, method)(*args)


def randbool(size=(), p: float = 0.5):
    return np.random.rand(*size) <= p


RANDS_CHARS = np.array(list(string.ascii_letters + string.digits), dtype=(np.str_, 1))
RANDU_CHARS = np.array(
    list("".join(map(chr, range(1488, 1488 + 26))) + string.digits),
    dtype=(np.unicode_, 1),
)


def rands_array(nchars, size, dtype="O"):
    """
    Generate an array of byte strings.
    """
    retval = (
        np.random.choice(RANDS_CHARS, size=nchars * np.prod(size))
        .view((np.str_, nchars))
        .reshape(size)
    )
    return retval.astype(dtype)


def randu_array(nchars, size, dtype="O"):
    """
    Generate an array of unicode strings.
    """
    retval = (
        np.random.choice(RANDU_CHARS, size=nchars * np.prod(size))
        .view((np.unicode_, nchars))
        .reshape(size)
    )
    return retval.astype(dtype)


def rands(nchars):
    """
    Generate one random byte string.

    See `rands_array` if you want to create an array of random strings.

    """
    return "".join(np.random.choice(RANDS_CHARS, nchars))


def close(fignum=None):
    from matplotlib.pyplot import close as _close, get_fignums

    if fignum is None:
        for fignum in get_fignums():
            _close(fignum)
    else:
        _close(fignum)


# -----------------------------------------------------------------------------
# Comparators


def equalContents(arr1, arr2) -> bool:
    """
    Checks if the set of unique elements of arr1 and arr2 are equivalent.
    """
    return frozenset(arr1) == frozenset(arr2)


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
    if box_cls is pd.array:
        expected = pd.array(expected)
    elif box_cls is pd.Index:
        expected = pd.Index(expected)
    elif box_cls is pd.Series:
        expected = pd.Series(expected)
    elif box_cls is pd.DataFrame:
        expected = pd.Series(expected).to_frame()
        if transpose:
            # for vector operations, we need a DataFrame to be a single-row,
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
    dtype = getattr(obj, "dtype", None)

    if is_period_dtype(dtype):
        return period_array(obj)
    elif is_datetime64_dtype(dtype) or is_datetime64tz_dtype(dtype):
        return DatetimeArray._from_sequence(obj)
    elif is_timedelta64_dtype(dtype):
        return TimedeltaArray._from_sequence(obj)
    else:
        return np.array(obj)


# -----------------------------------------------------------------------------
# Others


def getCols(k):
    return string.ascii_uppercase[:k]


# make index
def makeStringIndex(k=10, name=None):
    return Index(rands_array(nchars=10, size=k), name=name)


def makeUnicodeIndex(k=10, name=None):
    return Index(randu_array(nchars=10, size=k), name=name)


def makeCategoricalIndex(k=10, n=3, name=None, **kwargs):
    """ make a length k index or n categories """
    x = rands_array(nchars=4, size=n)
    return CategoricalIndex(
        Categorical.from_codes(np.arange(k) % n, categories=x), name=name, **kwargs
    )


def makeIntervalIndex(k=10, name=None, **kwargs):
    """ make a length k IntervalIndex """
    x = np.linspace(0, 100, num=(k + 1))
    return IntervalIndex.from_breaks(x, name=name, **kwargs)


def makeBoolIndex(k=10, name=None):
    if k == 1:
        return Index([True], name=name)
    elif k == 2:
        return Index([False, True], name=name)
    return Index([False, True] + [False] * (k - 2), name=name)


def makeIntIndex(k=10, name=None):
    return Index(list(range(k)), name=name)


def makeUIntIndex(k=10, name=None):
    return Index([2 ** 63 + i for i in range(k)], name=name)


def makeRangeIndex(k=10, name=None, **kwargs):
    return RangeIndex(0, k, 1, name=name, **kwargs)


def makeFloatIndex(k=10, name=None):
    values = sorted(np.random.random_sample(k)) - np.random.random_sample(1)
    return Index(values * (10 ** np.random.randint(0, 9)), name=name)


def makeDateIndex(k=10, freq="B", name=None, **kwargs):
    dt = datetime(2000, 1, 1)
    dr = bdate_range(dt, periods=k, freq=freq, name=name)
    return DatetimeIndex(dr, name=name, **kwargs)


def makeTimedeltaIndex(k=10, freq="D", name=None, **kwargs):
    return pd.timedelta_range(start="1 day", periods=k, freq=freq, name=name, **kwargs)


def makePeriodIndex(k=10, name=None, **kwargs):
    dt = datetime(2000, 1, 1)
    return pd.period_range(start=dt, periods=k, freq="B", name=name, **kwargs)


def makeMultiIndex(k=10, names=None, **kwargs):
    return MultiIndex.from_product((("foo", "bar"), (1, 2)), names=names, **kwargs)


_names = [
    "Alice",
    "Bob",
    "Charlie",
    "Dan",
    "Edith",
    "Frank",
    "George",
    "Hannah",
    "Ingrid",
    "Jerry",
    "Kevin",
    "Laura",
    "Michael",
    "Norbert",
    "Oliver",
    "Patricia",
    "Quinn",
    "Ray",
    "Sarah",
    "Tim",
    "Ursula",
    "Victor",
    "Wendy",
    "Xavier",
    "Yvonne",
    "Zelda",
]


def _make_timeseries(start="2000-01-01", end="2000-12-31", freq="1D", seed=None):
    """
    Make a DataFrame with a DatetimeIndex

    Parameters
    ----------
    start : str or Timestamp, default "2000-01-01"
        The start of the index. Passed to date_range with `freq`.
    end : str or Timestamp, default "2000-12-31"
        The end of the index. Passed to date_range with `freq`.
    freq : str or Freq
        The frequency to use for the DatetimeIndex
    seed : int, optional
        The random state seed.

        * name : object dtype with string names
        * id : int dtype with
        * x, y : float dtype

    Examples
    --------
    >>> _make_timeseries()
                  id    name         x         y
    timestamp
    2000-01-01   982   Frank  0.031261  0.986727
    2000-01-02  1025   Edith -0.086358 -0.032920
    2000-01-03   982   Edith  0.473177  0.298654
    2000-01-04  1009   Sarah  0.534344 -0.750377
    2000-01-05   963   Zelda -0.271573  0.054424
    ...          ...     ...       ...       ...
    2000-12-27   980  Ingrid -0.132333 -0.422195
    2000-12-28   972   Frank -0.376007 -0.298687
    2000-12-29  1009  Ursula -0.865047 -0.503133
    2000-12-30  1000  Hannah -0.063757 -0.507336
    2000-12-31   972     Tim -0.869120  0.531685
    """
    index = pd.date_range(start=start, end=end, freq=freq, name="timestamp")
    n = len(index)
    state = np.random.RandomState(seed)
    columns = {
        "name": state.choice(_names, size=n),
        "id": state.poisson(1000, size=n),
        "x": state.rand(n) * 2 - 1,
        "y": state.rand(n) * 2 - 1,
    }
    df = pd.DataFrame(columns, index=index, columns=sorted(columns))
    if df.index[-1] == end:
        df = df.iloc[:-1]
    return df


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
    yield from make_index_funcs


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
        # pandas\_testing.py:1986: error: Cannot call function of unknown type
        yield make_index_func(k=k)  # type: ignore[operator]


# make series
def makeFloatSeries(name=None):
    index = makeStringIndex(_N)
    return Series(np.random.randn(_N), index=index, name=name)


def makeStringSeries(name=None):
    index = makeStringIndex(_N)
    return Series(np.random.randn(_N), index=index, name=name)


def makeObjectSeries(name=None):
    data = makeStringIndex(_N)
    data = Index(data, dtype=object)
    index = makeStringIndex(_N)
    return Series(data, index=index, name=name)


def getSeriesData():
    index = makeStringIndex(_N)
    return {c: Series(np.random.randn(_N), index=index) for c in getCols(_K)}


def makeTimeSeries(nper=None, freq="B", name=None):
    if nper is None:
        nper = _N
    return Series(
        np.random.randn(nper), index=makeDateIndex(nper, freq=freq), name=name
    )


def makePeriodSeries(nper=None, name=None):
    if nper is None:
        nper = _N
    return Series(np.random.randn(nper), index=makePeriodIndex(nper), name=name)


def getTimeSeriesData(nper=None, freq="B"):
    return {c: makeTimeSeries(nper, freq) for c in getCols(_K)}


def getPeriodData(nper=None):
    return {c: makePeriodSeries(nper) for c in getCols(_K)}


# make frame
def makeTimeDataFrame(nper=None, freq="B"):
    data = getTimeSeriesData(nper, freq)
    return DataFrame(data)


def makeDataFrame():
    data = getSeriesData()
    return DataFrame(data)


def getMixedTypeDict():
    index = Index(["a", "b", "c", "d", "e"])

    data = {
        "A": [0.0, 1.0, 2.0, 3.0, 4.0],
        "B": [0.0, 1.0, 0.0, 1.0, 0.0],
        "C": ["foo1", "foo2", "foo3", "foo4", "foo5"],
        "D": bdate_range("1/1/2009", periods=5),
    }

    return index, data


def makeMixedDataFrame():
    return DataFrame(getMixedTypeDict()[1])


def makePeriodFrame(nper=None):
    data = getPeriodData(nper)
    return DataFrame(data)


def makeCustomIndex(
    nentries, nlevels, prefix="#", names=False, ndupe_l=None, idx_type=None
):
    """
    Create an index/multindex with given dimensions, levels, names, etc'

    nentries - number of entries in index
    nlevels - number of levels (> 1 produces multindex)
    prefix - a string prefix for labels
    names - (Optional), bool or list of strings. if True will use default
       names, if false will use no names, if a list is given, the name of
       each level in the index will be taken from the list.
    ndupe_l - (Optional), list of ints, the number of rows for which the
       label will repeated at the corresponding level, you can specify just
       the first few, the rest will use the default ndupe_l of 1.
       len(ndupe_l) <= nlevels.
    idx_type - "i"/"f"/"s"/"u"/"dt"/"p"/"td".
       If idx_type is not None, `idx_nlevels` must be 1.
       "i"/"f" creates an integer/float index,
       "s"/"u" creates a string/unicode index
       "dt" create a datetime index.
       "td" create a datetime index.

        if unspecified, string labels will be generated.
    """
    if ndupe_l is None:
        ndupe_l = [1] * nlevels
    assert is_sequence(ndupe_l) and len(ndupe_l) <= nlevels
    assert names is None or names is False or names is True or len(names) is nlevels
    assert idx_type is None or (
        idx_type in ("i", "f", "s", "u", "dt", "p", "td") and nlevels == 1
    )

    if names is True:
        # build default names
        names = [prefix + str(i) for i in range(nlevels)]
    if names is False:
        # pass None to index constructor for no name
        names = None

    # make singleton case uniform
    if isinstance(names, str) and nlevels == 1:
        names = [names]

    # specific 1D index type requested?
    idx_func = {
        "i": makeIntIndex,
        "f": makeFloatIndex,
        "s": makeStringIndex,
        "u": makeUnicodeIndex,
        "dt": makeDateIndex,
        "td": makeTimedeltaIndex,
        "p": makePeriodIndex,
    }.get(idx_type)
    if idx_func:
        # pandas\_testing.py:2120: error: Cannot call function of unknown type
        idx = idx_func(nentries)  # type: ignore[operator]
        # but we need to fill in the name
        if names:
            idx.name = names[0]
        return idx
    elif idx_type is not None:
        raise ValueError(
            f"{repr(idx_type)} is not a legal value for `idx_type`, "
            "use  'i'/'f'/'s'/'u'/'dt'/'p'/'td'."
        )

    if len(ndupe_l) < nlevels:
        ndupe_l.extend([1] * (nlevels - len(ndupe_l)))
    assert len(ndupe_l) == nlevels

    assert all(x > 0 for x in ndupe_l)

    tuples = []
    for i in range(nlevels):

        def keyfunc(x):
            import re

            numeric_tuple = re.sub(r"[^\d_]_?", "", x).split("_")
            return [int(num) for num in numeric_tuple]

        # build a list of lists to create the index from
        div_factor = nentries // ndupe_l[i] + 1
        # pandas\_testing.py:2148: error: Need type annotation for 'cnt'
        cnt = Counter()  # type: ignore[var-annotated]
        for j in range(div_factor):
            label = f"{prefix}_l{i}_g{j}"
            cnt[label] = ndupe_l[i]
        # cute Counter trick
        result = sorted(cnt.elements(), key=keyfunc)[:nentries]
        tuples.append(result)

    tuples = list(zip(*tuples))

    # convert tuples to index
    if nentries == 1:
        # we have a single level of tuples, i.e. a regular Index
        index = Index(tuples[0], name=names[0])
    elif nlevels == 1:
        name = None if names is None else names[0]
        index = Index((x[0] for x in tuples), name=name)
    else:
        index = MultiIndex.from_tuples(tuples, names=names)
    return index


def makeCustomDataframe(
    nrows,
    ncols,
    c_idx_names=True,
    r_idx_names=True,
    c_idx_nlevels=1,
    r_idx_nlevels=1,
    data_gen_f=None,
    c_ndupe_l=None,
    r_ndupe_l=None,
    dtype=None,
    c_idx_type=None,
    r_idx_type=None,
):
    """
    Create a DataFrame using supplied parameters.

    Parameters
    ----------
    nrows,  ncols - number of data rows/cols
    c_idx_names, idx_names  - False/True/list of strings,  yields No names ,
            default names or uses the provided names for the levels of the
            corresponding index. You can provide a single string when
            c_idx_nlevels ==1.
    c_idx_nlevels - number of levels in columns index. > 1 will yield MultiIndex
    r_idx_nlevels - number of levels in rows index. > 1 will yield MultiIndex
    data_gen_f - a function f(row,col) which return the data value
            at that position, the default generator used yields values of the form
            "RxCy" based on position.
    c_ndupe_l, r_ndupe_l - list of integers, determines the number
            of duplicates for each label at a given level of the corresponding
            index. The default `None` value produces a multiplicity of 1 across
            all levels, i.e. a unique index. Will accept a partial list of length
            N < idx_nlevels, for just the first N levels. If ndupe doesn't divide
            nrows/ncol, the last label might have lower multiplicity.
    dtype - passed to the DataFrame constructor as is, in case you wish to
            have more control in conjunction with a custom `data_gen_f`
    r_idx_type, c_idx_type -  "i"/"f"/"s"/"u"/"dt"/"td".
        If idx_type is not None, `idx_nlevels` must be 1.
        "i"/"f" creates an integer/float index,
        "s"/"u" creates a string/unicode index
        "dt" create a datetime index.
        "td" create a timedelta index.

            if unspecified, string labels will be generated.

    Examples
    --------
    # 5 row, 3 columns, default names on both, single index on both axis
    >> makeCustomDataframe(5,3)

    # make the data a random int between 1 and 100
    >> mkdf(5,3,data_gen_f=lambda r,c:randint(1,100))

    # 2-level multiindex on rows with each label duplicated
    # twice on first level, default names on both axis, single
    # index on both axis
    >> a=makeCustomDataframe(5,3,r_idx_nlevels=2,r_ndupe_l=[2])

    # DatetimeIndex on row, index with unicode labels on columns
    # no names on either axis
    >> a=makeCustomDataframe(5,3,c_idx_names=False,r_idx_names=False,
                             r_idx_type="dt",c_idx_type="u")

    # 4-level multindex on rows with names provided, 2-level multindex
    # on columns with default labels and default names.
    >> a=makeCustomDataframe(5,3,r_idx_nlevels=4,
                             r_idx_names=["FEE","FI","FO","FAM"],
                             c_idx_nlevels=2)

    >> a=mkdf(5,3,r_idx_nlevels=2,c_idx_nlevels=4)
    """
    assert c_idx_nlevels > 0
    assert r_idx_nlevels > 0
    assert r_idx_type is None or (
        r_idx_type in ("i", "f", "s", "u", "dt", "p", "td") and r_idx_nlevels == 1
    )
    assert c_idx_type is None or (
        c_idx_type in ("i", "f", "s", "u", "dt", "p", "td") and c_idx_nlevels == 1
    )

    columns = makeCustomIndex(
        ncols,
        nlevels=c_idx_nlevels,
        prefix="C",
        names=c_idx_names,
        ndupe_l=c_ndupe_l,
        idx_type=c_idx_type,
    )
    index = makeCustomIndex(
        nrows,
        nlevels=r_idx_nlevels,
        prefix="R",
        names=r_idx_names,
        ndupe_l=r_ndupe_l,
        idx_type=r_idx_type,
    )

    # by default, generate data based on location
    if data_gen_f is None:
        data_gen_f = lambda r, c: f"R{r}C{c}"

    data = [[data_gen_f(r, c) for c in range(ncols)] for r in range(nrows)]

    return DataFrame(data, index, columns, dtype=dtype)


def _create_missing_idx(nrows, ncols, density, random_state=None):
    if random_state is None:
        random_state = np.random
    else:
        random_state = np.random.RandomState(random_state)

    # below is cribbed from scipy.sparse
    size = int(np.round((1 - density) * nrows * ncols))
    # generate a few more to ensure unique values
    min_rows = 5
    fac = 1.02
    extra_size = min(size + min_rows, fac * size)

    def _gen_unique_rand(rng, _extra_size):
        ind = rng.rand(int(_extra_size))
        return np.unique(np.floor(ind * nrows * ncols))[:size]

    ind = _gen_unique_rand(random_state, extra_size)
    while ind.size < size:
        extra_size *= 1.05
        ind = _gen_unique_rand(random_state, extra_size)

    j = np.floor(ind * 1.0 / nrows).astype(int)
    i = (ind - j * nrows).astype(int)
    return i.tolist(), j.tolist()


def makeMissingDataframe(density=0.9, random_state=None):
    df = makeDataFrame()
    # pandas\_testing.py:2306: error: "_create_missing_idx" gets multiple
    # values for keyword argument "density"  [misc]

    # pandas\_testing.py:2306: error: "_create_missing_idx" gets multiple
    # values for keyword argument "random_state"  [misc]
    i, j = _create_missing_idx(  # type: ignore[misc]
        *df.shape, density=density, random_state=random_state
    )
    df.values[i, j] = np.nan
    return df


def optional_args(decorator):
    """
    allows a decorator to take optional positional and keyword arguments.
    Assumes that taking a single, callable, positional argument means that
    it is decorating a function, i.e. something like this::

        @my_decorator
        def function(): pass

    Calls decorator with decorator(f, *args, **kwargs)
    """

    @wraps(decorator)
    def wrapper(*args, **kwargs):
        def dec(f):
            return decorator(f, *args, **kwargs)

        is_decorating = not kwargs and len(args) == 1 and callable(args[0])
        if is_decorating:
            f = args[0]
            # pandas\_testing.py:2331: error: Incompatible types in assignment
            # (expression has type "List[<nothing>]", variable has type
            # "Tuple[Any, ...]")
            args = []  # type: ignore[assignment]
            return dec(f)
        else:
            return dec

    return wrapper


# skip tests on exceptions with this message
_network_error_messages = (
    # 'urlopen error timed out',
    # 'timeout: timed out',
    # 'socket.timeout: timed out',
    "timed out",
    "Server Hangup",
    "HTTP Error 503: Service Unavailable",
    "502: Proxy Error",
    "HTTP Error 502: internal error",
    "HTTP Error 502",
    "HTTP Error 503",
    "HTTP Error 403",
    "HTTP Error 400",
    "Temporary failure in name resolution",
    "Name or service not known",
    "Connection refused",
    "certificate verify",
)

# or this e.errno/e.reason.errno
_network_errno_vals = (
    101,  # Network is unreachable
    111,  # Connection refused
    110,  # Connection timed out
    104,  # Connection reset Error
    54,  # Connection reset by peer
    60,  # urllib.error.URLError: [Errno 60] Connection timed out
)

# Both of the above shouldn't mask real issues such as 404's
# or refused connections (changed DNS).
# But some tests (test_data yahoo) contact incredibly flakey
# servers.

# and conditionally raise on exception types in _get_default_network_errors


def _get_default_network_errors():
    # Lazy import for http.client because it imports many things from the stdlib
    import http.client

    return (IOError, http.client.HTTPException, TimeoutError)


def can_connect(url, error_classes=None):
    """
    Try to connect to the given url. True if succeeds, False if IOError
    raised

    Parameters
    ----------
    url : basestring
        The URL to try to connect to

    Returns
    -------
    connectable : bool
        Return True if no IOError (unable to connect) or URLError (bad url) was
        raised
    """
    if error_classes is None:
        error_classes = _get_default_network_errors()

    try:
        with urlopen(url):
            pass
    except error_classes:
        return False
    else:
        return True


@optional_args
def network(
    t,
    url="https://www.google.com",
    raise_on_error=_RAISE_NETWORK_ERROR_DEFAULT,
    check_before_test=False,
    error_classes=None,
    skip_errnos=_network_errno_vals,
    _skip_on_messages=_network_error_messages,
):
    """
    Label a test as requiring network connection and, if an error is
    encountered, only raise if it does not find a network connection.

    In comparison to ``network``, this assumes an added contract to your test:
    you must assert that, under normal conditions, your test will ONLY fail if
    it does not have network connectivity.

    You can call this in 3 ways: as a standard decorator, with keyword
    arguments, or with a positional argument that is the url to check.

    Parameters
    ----------
    t : callable
        The test requiring network connectivity.
    url : path
        The url to test via ``pandas.io.common.urlopen`` to check
        for connectivity. Defaults to 'https://www.google.com'.
    raise_on_error : bool
        If True, never catches errors.
    check_before_test : bool
        If True, checks connectivity before running the test case.
    error_classes : tuple or Exception
        error classes to ignore. If not in ``error_classes``, raises the error.
        defaults to IOError. Be careful about changing the error classes here.
    skip_errnos : iterable of int
        Any exception that has .errno or .reason.erno set to one
        of these values will be skipped with an appropriate
        message.
    _skip_on_messages: iterable of string
        any exception e for which one of the strings is
        a substring of str(e) will be skipped with an appropriate
        message. Intended to suppress errors where an errno isn't available.

    Notes
    -----
    * ``raise_on_error`` supersedes ``check_before_test``

    Returns
    -------
    t : callable
        The decorated test ``t``, with checks for connectivity errors.

    Example
    -------

    Tests decorated with @network will fail if it's possible to make a network
    connection to another URL (defaults to google.com)::

      >>> from pandas._testing import network
      >>> from pandas.io.common import urlopen
      >>> @network
      ... def test_network():
      ...     with urlopen("rabbit://bonanza.com"):
      ...         pass
      Traceback
         ...
      URLError: <urlopen error unknown url type: rabit>

      You can specify alternative URLs::

        >>> @network("https://www.yahoo.com")
        ... def test_something_with_yahoo():
        ...    raise IOError("Failure Message")
        >>> test_something_with_yahoo()
        Traceback (most recent call last):
            ...
        IOError: Failure Message

    If you set check_before_test, it will check the url first and not run the
    test on failure::

        >>> @network("failing://url.blaher", check_before_test=True)
        ... def test_something():
        ...     print("I ran!")
        ...     raise ValueError("Failure")
        >>> test_something()
        Traceback (most recent call last):
            ...

    Errors not related to networking will always be raised.
    """
    from pytest import skip

    if error_classes is None:
        error_classes = _get_default_network_errors()

    t.network = True

    @wraps(t)
    def wrapper(*args, **kwargs):
        if (
            check_before_test
            and not raise_on_error
            and not can_connect(url, error_classes)
        ):
            skip()
        try:
            return t(*args, **kwargs)
        except Exception as err:
            errno = getattr(err, "errno", None)
            if not errno and hasattr(errno, "reason"):
                # pandas\_testing.py:2521: error: "Exception" has no attribute
                # "reason"
                errno = getattr(err.reason, "errno", None)  # type: ignore[attr-defined]

            if errno in skip_errnos:
                skip(f"Skipping test due to known errno and error {err}")

            e_str = str(err)

            if any(m.lower() in e_str.lower() for m in _skip_on_messages):
                skip(
                    f"Skipping test because exception message is known and error {err}"
                )

            if not isinstance(err, error_classes):
                raise

            if raise_on_error or can_connect(url, error_classes):
                raise
            else:
                skip(f"Skipping test due to lack of connectivity and error {err}")

    return wrapper


with_connectivity_check = network


@contextmanager
def assert_produces_warning(
    expected_warning: Optional[Union[Type[Warning], bool]] = Warning,
    filter_level="always",
    check_stacklevel: bool = True,
    raise_on_extra_warnings: bool = True,
    match: Optional[str] = None,
):
    """
    Context manager for running code expected to either raise a specific
    warning, or not raise any warnings. Verifies that the code raises the
    expected warning, and that it does not raise any other unexpected
    warnings. It is basically a wrapper around ``warnings.catch_warnings``.

    Parameters
    ----------
    expected_warning : {Warning, False, None}, default Warning
        The type of Exception raised. ``exception.Warning`` is the base
        class for all warnings. To check that no warning is returned,
        specify ``False`` or ``None``.
    filter_level : str or None, default "always"
        Specifies whether warnings are ignored, displayed, or turned
        into errors.
        Valid values are:

        * "error" - turns matching warnings into exceptions
        * "ignore" - discard the warning
        * "always" - always emit a warning
        * "default" - print the warning the first time it is generated
          from each location
        * "module" - print the warning the first time it is generated
          from each module
        * "once" - print the warning the first time it is generated

    check_stacklevel : bool, default True
        If True, displays the line that called the function containing
        the warning to show were the function is called. Otherwise, the
        line that implements the function is displayed.
    raise_on_extra_warnings : bool, default True
        Whether extra warnings not of the type `expected_warning` should
        cause the test to fail.
    match : str, optional
        Match warning message.

    Examples
    --------
    >>> import warnings
    >>> with assert_produces_warning():
    ...     warnings.warn(UserWarning())
    ...
    >>> with assert_produces_warning(False):
    ...     warnings.warn(RuntimeWarning())
    ...
    Traceback (most recent call last):
        ...
    AssertionError: Caused unexpected warning(s): ['RuntimeWarning'].
    >>> with assert_produces_warning(UserWarning):
    ...     warnings.warn(RuntimeWarning())
    Traceback (most recent call last):
        ...
    AssertionError: Did not see expected warning of class 'UserWarning'.

    ..warn:: This is *not* thread-safe.
    """
    __tracebackhide__ = True

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter(filter_level)
        yield w

        if expected_warning:
            expected_warning = cast(Type[Warning], expected_warning)
            _assert_caught_expected_warning(
                caught_warnings=w,
                expected_warning=expected_warning,
                match=match,
                check_stacklevel=check_stacklevel,
            )

        if raise_on_extra_warnings:
            _assert_caught_no_extra_warnings(
                caught_warnings=w,
                expected_warning=expected_warning,
            )


def _assert_caught_expected_warning(
    *,
    caught_warnings: Sequence[warnings.WarningMessage],
    expected_warning: Type[Warning],
    match: Optional[str],
    check_stacklevel: bool,
) -> None:
    """Assert that there was the expected warning among the caught warnings."""
    saw_warning = False
    matched_message = False

    for actual_warning in caught_warnings:
        if issubclass(actual_warning.category, expected_warning):
            saw_warning = True

            if check_stacklevel and issubclass(
                actual_warning.category, (FutureWarning, DeprecationWarning)
            ):
                _assert_raised_with_correct_stacklevel(actual_warning)

            if match is not None and re.search(match, str(actual_warning.message)):
                matched_message = True

    if not saw_warning:
        raise AssertionError(
            f"Did not see expected warning of class "
            f"{repr(expected_warning.__name__)}"
        )

    if match and not matched_message:
        raise AssertionError(
            f"Did not see warning {repr(expected_warning.__name__)} "
            f"matching {match}"
        )


def _assert_caught_no_extra_warnings(
    *,
    caught_warnings: Sequence[warnings.WarningMessage],
    expected_warning: Optional[Union[Type[Warning], bool]],
) -> None:
    """Assert that no extra warnings apart from the expected ones are caught."""
    extra_warnings = []

    for actual_warning in caught_warnings:
        if _is_unexpected_warning(actual_warning, expected_warning):
            extra_warnings.append(
                (
                    actual_warning.category.__name__,
                    actual_warning.message,
                    actual_warning.filename,
                    actual_warning.lineno,
                )
            )

    if extra_warnings:
        raise AssertionError(f"Caused unexpected warning(s): {repr(extra_warnings)}")


def _is_unexpected_warning(
    actual_warning: warnings.WarningMessage,
    expected_warning: Optional[Union[Type[Warning], bool]],
) -> bool:
    """Check if the actual warning issued is unexpected."""
    if actual_warning and not expected_warning:
        return True
    expected_warning = cast(Type[Warning], expected_warning)
    return bool(not issubclass(actual_warning.category, expected_warning))


def _assert_raised_with_correct_stacklevel(
    actual_warning: warnings.WarningMessage,
) -> None:
    from inspect import getframeinfo, stack

    caller = getframeinfo(stack()[4][0])
    msg = (
        "Warning not set with correct stacklevel. "
        f"File where warning is raised: {actual_warning.filename} != "
        f"{caller.filename}. Warning message: {actual_warning.message}"
    )
    assert actual_warning.filename == caller.filename, msg


class RNGContext:
    """
    Context manager to set the numpy random number generator speed. Returns
    to the original value upon exiting the context manager.

    Parameters
    ----------
    seed : int
        Seed for numpy.random.seed

    Examples
    --------
    with RNGContext(42):
        np.random.randn()
    """

    def __init__(self, seed):
        self.seed = seed

    def __enter__(self):

        self.start_state = np.random.get_state()
        np.random.seed(self.seed)

    def __exit__(self, exc_type, exc_value, traceback):

        np.random.set_state(self.start_state)


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
    return sep.join(rows_list) + sep


def external_error_raised(expected_exception: Type[Exception]) -> ContextManager:
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


cython_table = pd.core.base.SelectionMixin._cython_table.items()


def get_cython_table_params(ndframe, func_names_and_expected):
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
            for func, name in cython_table
            if name == func_name
        ]
    return results


def get_op_from_name(op_name: str) -> Callable:
    """
    The operator function for a given op name.

    Parameters
    ----------
    op_name : string
        The op name, in form of "add" or "__add__".

    Returns
    -------
    function
        A function performing the operation.
    """
    short_opname = op_name.strip("_")
    try:
        op = getattr(operator, short_opname)
    except AttributeError:
        # Assume it is the reverse operator
        rop = getattr(operator, short_opname[1:])
        op = lambda x, y: rop(y, x)

    return op


# -----------------------------------------------------------------------------
# Indexing test helpers


def getitem(x):
    return x


def setitem(x):
    return x


def loc(x):
    return x.loc


def iloc(x):
    return x.iloc
