import bz2
from collections import Counter
from contextlib import contextmanager
from datetime import datetime
from functools import wraps
import gzip
import os
from shutil import rmtree
import string
import tempfile
from typing import Any, Callable, List, Optional, Type, Union, cast
import warnings
import zipfile

import numpy as np
from numpy.random import rand, randn

from pandas._config.localization import (  # noqa:F401
    can_set_locale,
    get_locales,
    set_locale,
)

import pandas._libs.testing as _testing
from pandas._typing import FilePathOrBuffer, FrameOrSeries
from pandas.compat import _get_lzma_file, _import_lzma

from pandas.core.dtypes.common import (
    is_bool,
    is_categorical_dtype,
    is_datetime64_dtype,
    is_datetime64tz_dtype,
    is_extension_array_dtype,
    is_interval_dtype,
    is_list_like,
    is_number,
    is_period_dtype,
    is_sequence,
    is_timedelta64_dtype,
    needs_i8_conversion,
)
from pandas.core.dtypes.missing import array_equivalent

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
from pandas.core.algorithms import take_1d
from pandas.core.arrays import (
    DatetimeArray,
    ExtensionArray,
    IntervalArray,
    PeriodArray,
    TimedeltaArray,
    period_array,
)

from pandas.io.common import urlopen
from pandas.io.formats.printing import pprint_thing

lzma = _import_lzma()

N = 30
K = 4
_RAISE_NETWORK_ERROR_DEFAULT = False

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


@contextmanager
def decompress_file(path, compression):
    """
    Open a compressed file and return a file object.

    Parameters
    ----------
    path : str
        The path where the file is read from.

    compression : {'gzip', 'bz2', 'zip', 'xz', None}
        Name of the decompression to use

    Returns
    -------
    file object
    """
    if compression is None:
        f = open(path, "rb")
    elif compression == "gzip":
        f = gzip.open(path, "rb")
    elif compression == "bz2":
        f = bz2.BZ2File(path, "rb")
    elif compression == "xz":
        f = _get_lzma_file(lzma)(path, "rb")
    elif compression == "zip":
        zip_file = zipfile.ZipFile(path)
        zip_names = zip_file.namelist()
        if len(zip_names) == 1:
            f = zip_file.open(zip_names.pop())
        else:
            raise ValueError(f"ZIP file {path} error. Only one file per ZIP.")
    else:
        raise ValueError(f"Unrecognized compression type: {compression}")

    try:
        yield f
    finally:
        f.close()
        if compression == "zip":
            zip_file.close()


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


def assert_almost_equal(
    left,
    right,
    check_dtype: Union[bool, str] = "equiv",
    check_less_precise: Union[bool, int] = False,
    **kwargs,
):
    """
    Check that the left and right objects are approximately equal.

    By approximately equal, we refer to objects that are numbers or that
    contain numbers which may be equivalent to specific levels of precision.

    Parameters
    ----------
    left : object
    right : object
    check_dtype : bool or {'equiv'}, default 'equiv'
        Check dtype if both a and b are the same type. If 'equiv' is passed in,
        then `RangeIndex` and `Int64Index` are also considered equivalent
        when doing type checking.
    check_less_precise : bool or int, default False
        Specify comparison precision. 5 digits (False) or 3 digits (True)
        after decimal points are compared. If int, then specify the number
        of digits to compare.

        When comparing two numbers, if the first number has magnitude less
        than 1e-5, we compare the two numbers directly and check whether
        they are equivalent within the specified precision. Otherwise, we
        compare the **ratio** of the second number to the first number and
        check whether it is equivalent to 1 within the specified precision.
    """
    if isinstance(left, pd.Index):
        assert_index_equal(
            left,
            right,
            check_exact=False,
            exact=check_dtype,
            check_less_precise=check_less_precise,
            **kwargs,
        )

    elif isinstance(left, pd.Series):
        assert_series_equal(
            left,
            right,
            check_exact=False,
            check_dtype=check_dtype,
            check_less_precise=check_less_precise,
            **kwargs,
        )

    elif isinstance(left, pd.DataFrame):
        assert_frame_equal(
            left,
            right,
            check_exact=False,
            check_dtype=check_dtype,
            check_less_precise=check_less_precise,
            **kwargs,
        )

    else:
        # Other sequences.
        if check_dtype:
            if is_number(left) and is_number(right):
                # Do not compare numeric classes, like np.float64 and float.
                pass
            elif is_bool(left) and is_bool(right):
                # Do not compare bool classes, like np.bool_ and bool.
                pass
            else:
                if isinstance(left, np.ndarray) or isinstance(right, np.ndarray):
                    obj = "numpy array"
                else:
                    obj = "Input"
                assert_class_equal(left, right, obj=obj)
        _testing.assert_almost_equal(
            left,
            right,
            check_dtype=check_dtype,
            check_less_precise=check_less_precise,
            **kwargs,
        )


def _check_isinstance(left, right, cls):
    """
    Helper method for our assert_* methods that ensures that
    the two objects being compared have the right type before
    proceeding with the comparison.

    Parameters
    ----------
    left : The first object being compared.
    right : The second object being compared.
    cls : The class type to check against.

    Raises
    ------
    AssertionError : Either `left` or `right` is not an instance of `cls`.
    """
    cls_name = cls.__name__

    if not isinstance(left, cls):
        raise AssertionError(
            f"{cls_name} Expected type {cls}, found {type(left)} instead"
        )
    if not isinstance(right, cls):
        raise AssertionError(
            f"{cls_name} Expected type {cls}, found {type(right)} instead"
        )


def assert_dict_equal(left, right, compare_keys: bool = True):

    _check_isinstance(left, right, dict)
    _testing.assert_dict_equal(left, right, compare_keys=compare_keys)


def randbool(size=(), p: float = 0.5):
    return rand(*size) <= p


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
    if dtype is None:
        return retval
    else:
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
    if dtype is None:
        return retval
    else:
        return retval.astype(dtype)


def rands(nchars):
    """
    Generate one random byte string.

    See `rands_array` if you want to create an array of random strings.

    """
    return "".join(np.random.choice(RANDS_CHARS, nchars))


def randu(nchars):
    """
    Generate one random unicode string.

    See `randu_array` if you want to create an array of random unicode strings.

    """
    return "".join(np.random.choice(RANDU_CHARS, nchars))


def close(fignum=None):
    from matplotlib.pyplot import get_fignums, close as _close

    if fignum is None:
        for fignum in get_fignums():
            _close(fignum)
    else:
        _close(fignum)


# -----------------------------------------------------------------------------
# contextmanager to ensure the file cleanup


@contextmanager
def ensure_clean(filename=None, return_filelike=False, **kwargs):
    """
    Gets a temporary path and agrees to remove on close.

    Parameters
    ----------
    filename : str (optional)
        if None, creates a temporary file which is then removed when out of
        scope. if passed, creates temporary file with filename as ending.
    return_filelike : bool (default False)
        if True, returns a file-like which is *always* cleaned. Necessary for
        savefig and other functions which want to append extensions.
    **kwargs
        Additional keywords passed in for creating a temporary file.
        :meth:`tempFile.TemporaryFile` is used when `return_filelike` is ``True``.
        :meth:`tempfile.mkstemp` is used when `return_filelike` is ``False``.
        Note that the `filename` parameter will be passed in as the `suffix`
        argument to either function.

    See Also
    --------
    tempfile.TemporaryFile
    tempfile.mkstemp
    """
    filename = filename or ""
    fd = None

    kwargs["suffix"] = filename

    if return_filelike:
        f = tempfile.TemporaryFile(**kwargs)

        try:
            yield f
        finally:
            f.close()
    else:
        # Don't generate tempfile if using a path with directory specified.
        if len(os.path.dirname(filename)):
            raise ValueError("Can't pass a qualified name to ensure_clean()")

        try:
            fd, filename = tempfile.mkstemp(**kwargs)
        except UnicodeEncodeError:
            import pytest

            pytest.skip("no unicode file names on this system")

        try:
            yield filename
        finally:
            try:
                os.close(fd)
            except OSError:
                print(f"Couldn't close file descriptor: {fd} (file: {filename})")
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except OSError as e:
                print(f"Exception on removing file: {e}")


@contextmanager
def ensure_clean_dir():
    """
    Get a temporary directory path and agrees to remove on close.

    Yields
    ------
    Temporary directory path
    """
    directory_name = tempfile.mkdtemp(suffix="")
    try:
        yield directory_name
    finally:
        try:
            rmtree(directory_name)
        except OSError:
            pass


@contextmanager
def ensure_safe_environment_variables():
    """
    Get a context manager to safely set environment variables

    All changes will be undone on close, hence environment variables set
    within this contextmanager will neither persist nor change global state.
    """
    saved_environ = dict(os.environ)
    try:
        yield
    finally:
        os.environ.clear()
        os.environ.update(saved_environ)


# -----------------------------------------------------------------------------
# Comparators


def equalContents(arr1, arr2) -> bool:
    """
    Checks if the set of unique elements of arr1 and arr2 are equivalent.
    """
    return frozenset(arr1) == frozenset(arr2)


def assert_index_equal(
    left: Index,
    right: Index,
    exact: Union[bool, str] = "equiv",
    check_names: bool = True,
    check_less_precise: Union[bool, int] = False,
    check_exact: bool = True,
    check_categorical: bool = True,
    obj: str = "Index",
) -> None:
    """
    Check that left and right Index are equal.

    Parameters
    ----------
    left : Index
    right : Index
    exact : bool or {'equiv'}, default 'equiv'
        Whether to check the Index class, dtype and inferred_type
        are identical. If 'equiv', then RangeIndex can be substituted for
        Int64Index as well.
    check_names : bool, default True
        Whether to check the names attribute.
    check_less_precise : bool or int, default False
        Specify comparison precision. Only used when check_exact is False.
        5 digits (False) or 3 digits (True) after decimal points are compared.
        If int, then specify the digits to compare.
    check_exact : bool, default True
        Whether to compare number exactly.
    check_categorical : bool, default True
        Whether to compare internal Categorical exactly.
    obj : str, default 'Index'
        Specify object name being compared, internally used to show appropriate
        assertion message.
    """
    __tracebackhide__ = True

    def _check_types(l, r, obj="Index"):
        if exact:
            assert_class_equal(l, r, exact=exact, obj=obj)

            # Skip exact dtype checking when `check_categorical` is False
            if check_categorical:
                assert_attr_equal("dtype", l, r, obj=obj)

            # allow string-like to have different inferred_types
            if l.inferred_type in ("string"):
                assert r.inferred_type in ("string")
            else:
                assert_attr_equal("inferred_type", l, r, obj=obj)

    def _get_ilevel_values(index, level):
        # accept level number only
        unique = index.levels[level]
        level_codes = index.codes[level]
        filled = take_1d(unique._values, level_codes, fill_value=unique._na_value)
        values = unique._shallow_copy(filled, name=index.names[level])
        return values

    # instance validation
    _check_isinstance(left, right, Index)

    # class / dtype comparison
    _check_types(left, right, obj=obj)

    # level comparison
    if left.nlevels != right.nlevels:
        msg1 = f"{obj} levels are different"
        msg2 = f"{left.nlevels}, {left}"
        msg3 = f"{right.nlevels}, {right}"
        raise_assert_detail(obj, msg1, msg2, msg3)

    # length comparison
    if len(left) != len(right):
        msg1 = f"{obj} length are different"
        msg2 = f"{len(left)}, {left}"
        msg3 = f"{len(right)}, {right}"
        raise_assert_detail(obj, msg1, msg2, msg3)

    # MultiIndex special comparison for little-friendly error messages
    if left.nlevels > 1:
        left = cast(MultiIndex, left)
        right = cast(MultiIndex, right)

        for level in range(left.nlevels):
            # cannot use get_level_values here because it can change dtype
            llevel = _get_ilevel_values(left, level)
            rlevel = _get_ilevel_values(right, level)

            lobj = f"MultiIndex level [{level}]"
            assert_index_equal(
                llevel,
                rlevel,
                exact=exact,
                check_names=check_names,
                check_less_precise=check_less_precise,
                check_exact=check_exact,
                obj=lobj,
            )
            # get_level_values may change dtype
            _check_types(left.levels[level], right.levels[level], obj=obj)

    # skip exact index checking when `check_categorical` is False
    if check_exact and check_categorical:
        if not left.equals(right):
            diff = np.sum((left.values != right.values).astype(int)) * 100.0 / len(left)
            msg = f"{obj} values are different ({np.round(diff, 5)} %)"
            raise_assert_detail(obj, msg, left, right)
    else:
        _testing.assert_almost_equal(
            left.values,
            right.values,
            check_less_precise=check_less_precise,
            check_dtype=exact,
            obj=obj,
            lobj=left,
            robj=right,
        )

    # metadata comparison
    if check_names:
        assert_attr_equal("names", left, right, obj=obj)
    if isinstance(left, pd.PeriodIndex) or isinstance(right, pd.PeriodIndex):
        assert_attr_equal("freq", left, right, obj=obj)
    if isinstance(left, pd.IntervalIndex) or isinstance(right, pd.IntervalIndex):
        assert_interval_array_equal(left.values, right.values)

    if check_categorical:
        if is_categorical_dtype(left) or is_categorical_dtype(right):
            assert_categorical_equal(left.values, right.values, obj=f"{obj} category")


def assert_class_equal(left, right, exact: Union[bool, str] = True, obj="Input"):
    """
    Checks classes are equal.
    """
    __tracebackhide__ = True

    def repr_class(x):
        if isinstance(x, Index):
            # return Index as it is to include values in the error message
            return x

        try:
            return type(x).__name__
        except AttributeError:
            return repr(type(x))

    if exact == "equiv":
        if type(left) != type(right):
            # allow equivalence of Int64Index/RangeIndex
            types = {type(left).__name__, type(right).__name__}
            if len(types - {"Int64Index", "RangeIndex"}):
                msg = f"{obj} classes are not equivalent"
                raise_assert_detail(obj, msg, repr_class(left), repr_class(right))
    elif exact:
        if type(left) != type(right):
            msg = f"{obj} classes are different"
            raise_assert_detail(obj, msg, repr_class(left), repr_class(right))


def assert_attr_equal(attr, left, right, obj="Attributes"):
    """
    checks attributes are equal. Both objects must have attribute.

    Parameters
    ----------
    attr : str
        Attribute name being compared.
    left : object
    right : object
    obj : str, default 'Attributes'
        Specify object name being compared, internally used to show appropriate
        assertion message
    """
    __tracebackhide__ = True

    left_attr = getattr(left, attr)
    right_attr = getattr(right, attr)

    if left_attr is right_attr:
        return True
    elif (
        is_number(left_attr)
        and np.isnan(left_attr)
        and is_number(right_attr)
        and np.isnan(right_attr)
    ):
        # np.nan
        return True

    try:
        result = left_attr == right_attr
    except TypeError:
        # datetimetz on rhs may raise TypeError
        result = False
    if not isinstance(result, bool):
        result = result.all()

    if result:
        return True
    else:
        msg = f'Attribute "{attr}" are different'
        raise_assert_detail(obj, msg, left_attr, right_attr)


def assert_is_valid_plot_return_object(objs):
    import matplotlib.pyplot as plt

    if isinstance(objs, (pd.Series, np.ndarray)):
        for el in objs.ravel():
            msg = (
                "one of 'objs' is not a matplotlib Axes instance, "
                f"type encountered {repr(type(el).__name__)}"
            )
            assert isinstance(el, (plt.Axes, dict)), msg
    else:
        msg = (
            "objs is neither an ndarray of Artist instances nor a single "
            "ArtistArtist instance, tuple, or dict, 'objs' is a "
            f"{repr(type(objs).__name__)}"
        )
        assert isinstance(objs, (plt.Artist, tuple, dict)), msg


def assert_is_sorted(seq):
    """Assert that the sequence is sorted."""
    if isinstance(seq, (Index, Series)):
        seq = seq.values
    # sorting does not change precisions
    assert_numpy_array_equal(seq, np.sort(np.array(seq)))


def assert_categorical_equal(
    left, right, check_dtype=True, check_category_order=True, obj="Categorical"
):
    """
    Test that Categoricals are equivalent.

    Parameters
    ----------
    left : Categorical
    right : Categorical
    check_dtype : bool, default True
        Check that integer dtype of the codes are the same
    check_category_order : bool, default True
        Whether the order of the categories should be compared, which
        implies identical integer codes.  If False, only the resulting
        values are compared.  The ordered attribute is
        checked regardless.
    obj : str, default 'Categorical'
        Specify object name being compared, internally used to show appropriate
        assertion message
    """
    _check_isinstance(left, right, Categorical)

    if check_category_order:
        assert_index_equal(left.categories, right.categories, obj=f"{obj}.categories")
        assert_numpy_array_equal(
            left.codes, right.codes, check_dtype=check_dtype, obj=f"{obj}.codes",
        )
    else:
        assert_index_equal(
            left.categories.sort_values(),
            right.categories.sort_values(),
            obj=f"{obj}.categories",
        )
        assert_index_equal(
            left.categories.take(left.codes),
            right.categories.take(right.codes),
            obj=f"{obj}.values",
        )

    assert_attr_equal("ordered", left, right, obj=obj)


def assert_interval_array_equal(left, right, exact="equiv", obj="IntervalArray"):
    """
    Test that two IntervalArrays are equivalent.

    Parameters
    ----------
    left, right : IntervalArray
        The IntervalArrays to compare.
    exact : bool or {'equiv'}, default 'equiv'
        Whether to check the Index class, dtype and inferred_type
        are identical. If 'equiv', then RangeIndex can be substituted for
        Int64Index as well.
    obj : str, default 'IntervalArray'
        Specify object name being compared, internally used to show appropriate
        assertion message
    """
    _check_isinstance(left, right, IntervalArray)

    assert_index_equal(left.left, right.left, exact=exact, obj=f"{obj}.left")
    assert_index_equal(left.right, right.right, exact=exact, obj=f"{obj}.left")
    assert_attr_equal("closed", left, right, obj=obj)


def assert_period_array_equal(left, right, obj="PeriodArray"):
    _check_isinstance(left, right, PeriodArray)

    assert_numpy_array_equal(left._data, right._data, obj=f"{obj}.values")
    assert_attr_equal("freq", left, right, obj=obj)


def assert_datetime_array_equal(left, right, obj="DatetimeArray"):
    __tracebackhide__ = True
    _check_isinstance(left, right, DatetimeArray)

    assert_numpy_array_equal(left._data, right._data, obj=f"{obj}._data")
    assert_attr_equal("freq", left, right, obj=obj)
    assert_attr_equal("tz", left, right, obj=obj)


def assert_timedelta_array_equal(left, right, obj="TimedeltaArray"):
    __tracebackhide__ = True
    _check_isinstance(left, right, TimedeltaArray)
    assert_numpy_array_equal(left._data, right._data, obj=f"{obj}._data")
    assert_attr_equal("freq", left, right, obj=obj)


def raise_assert_detail(obj, message, left, right, diff=None):
    __tracebackhide__ = True

    if isinstance(left, np.ndarray):
        left = pprint_thing(left)
    elif is_categorical_dtype(left):
        left = repr(left)

    if isinstance(right, np.ndarray):
        right = pprint_thing(right)
    elif is_categorical_dtype(right):
        right = repr(right)

    msg = f"""{obj} are different

{message}
[left]:  {left}
[right]: {right}"""

    if diff is not None:
        msg += f"\n[diff]: {diff}"

    raise AssertionError(msg)


def assert_numpy_array_equal(
    left,
    right,
    strict_nan=False,
    check_dtype=True,
    err_msg=None,
    check_same=None,
    obj="numpy array",
):
    """
    Check that 'np.ndarray' is equivalent.

    Parameters
    ----------
    left, right : numpy.ndarray or iterable
        The two arrays to be compared.
    strict_nan : bool, default False
        If True, consider NaN and None to be different.
    check_dtype : bool, default True
        Check dtype if both a and b are np.ndarray.
    err_msg : str, default None
        If provided, used as assertion message.
    check_same : None|'copy'|'same', default None
        Ensure left and right refer/do not refer to the same memory area.
    obj : str, default 'numpy array'
        Specify object name being compared, internally used to show appropriate
        assertion message.
    """
    __tracebackhide__ = True

    # instance validation
    # Show a detailed error message when classes are different
    assert_class_equal(left, right, obj=obj)
    # both classes must be an np.ndarray
    _check_isinstance(left, right, np.ndarray)

    def _get_base(obj):
        return obj.base if getattr(obj, "base", None) is not None else obj

    left_base = _get_base(left)
    right_base = _get_base(right)

    if check_same == "same":
        if left_base is not right_base:
            raise AssertionError(f"{repr(left_base)} is not {repr(right_base)}")
    elif check_same == "copy":
        if left_base is right_base:
            raise AssertionError(f"{repr(left_base)} is {repr(right_base)}")

    def _raise(left, right, err_msg):
        if err_msg is None:
            if left.shape != right.shape:
                raise_assert_detail(
                    obj, f"{obj} shapes are different", left.shape, right.shape,
                )

            diff = 0
            for l, r in zip(left, right):
                # count up differences
                if not array_equivalent(l, r, strict_nan=strict_nan):
                    diff += 1

            diff = diff * 100.0 / left.size
            msg = f"{obj} values are different ({np.round(diff, 5)} %)"
            raise_assert_detail(obj, msg, left, right)

        raise AssertionError(err_msg)

    # compare shape and values
    if not array_equivalent(left, right, strict_nan=strict_nan):
        _raise(left, right, err_msg)

    if check_dtype:
        if isinstance(left, np.ndarray) and isinstance(right, np.ndarray):
            assert_attr_equal("dtype", left, right, obj=obj)


def assert_extension_array_equal(
    left, right, check_dtype=True, check_less_precise=False, check_exact=False
):
    """
    Check that left and right ExtensionArrays are equal.

    Parameters
    ----------
    left, right : ExtensionArray
        The two arrays to compare.
    check_dtype : bool, default True
        Whether to check if the ExtensionArray dtypes are identical.
    check_less_precise : bool or int, default False
        Specify comparison precision. Only used when check_exact is False.
        5 digits (False) or 3 digits (True) after decimal points are compared.
        If int, then specify the digits to compare.
    check_exact : bool, default False
        Whether to compare number exactly.

    Notes
    -----
    Missing values are checked separately from valid values.
    A mask of missing values is computed for each and checked to match.
    The remaining all-valid values are cast to object dtype and checked.
    """
    assert isinstance(left, ExtensionArray), "left is not an ExtensionArray"
    assert isinstance(right, ExtensionArray), "right is not an ExtensionArray"
    if check_dtype:
        assert_attr_equal("dtype", left, right, obj="ExtensionArray")

    if hasattr(left, "asi8") and type(right) == type(left):
        # Avoid slow object-dtype comparisons
        assert_numpy_array_equal(left.asi8, right.asi8)
        return

    left_na = np.asarray(left.isna())
    right_na = np.asarray(right.isna())
    assert_numpy_array_equal(left_na, right_na, obj="ExtensionArray NA mask")

    left_valid = np.asarray(left[~left_na].astype(object))
    right_valid = np.asarray(right[~right_na].astype(object))
    if check_exact:
        assert_numpy_array_equal(left_valid, right_valid, obj="ExtensionArray")
    else:
        _testing.assert_almost_equal(
            left_valid,
            right_valid,
            check_dtype=check_dtype,
            check_less_precise=check_less_precise,
            obj="ExtensionArray",
        )


# This could be refactored to use the NDFrame.equals method
def assert_series_equal(
    left,
    right,
    check_dtype=True,
    check_index_type="equiv",
    check_series_type=True,
    check_less_precise=False,
    check_names=True,
    check_exact=False,
    check_datetimelike_compat=False,
    check_categorical=True,
    check_category_order=True,
    obj="Series",
):
    """
    Check that left and right Series are equal.

    Parameters
    ----------
    left : Series
    right : Series
    check_dtype : bool, default True
        Whether to check the Series dtype is identical.
    check_index_type : bool or {'equiv'}, default 'equiv'
        Whether to check the Index class, dtype and inferred_type
        are identical.
    check_series_type : bool, default True
        Whether to check the Series class is identical.
    check_less_precise : bool or int, default False
        Specify comparison precision. Only used when check_exact is False.
        5 digits (False) or 3 digits (True) after decimal points are compared.
        If int, then specify the digits to compare.

        When comparing two numbers, if the first number has magnitude less
        than 1e-5, we compare the two numbers directly and check whether
        they are equivalent within the specified precision. Otherwise, we
        compare the **ratio** of the second number to the first number and
        check whether it is equivalent to 1 within the specified precision.
    check_names : bool, default True
        Whether to check the Series and Index names attribute.
    check_exact : bool, default False
        Whether to compare number exactly.
    check_datetimelike_compat : bool, default False
        Compare datetime-like which is comparable ignoring dtype.
    check_categorical : bool, default True
        Whether to compare internal Categorical exactly.
    check_category_order : bool, default True
        Whether to compare category order of internal Categoricals

        .. versionadded:: 1.0.2
    obj : str, default 'Series'
        Specify object name being compared, internally used to show appropriate
        assertion message.
    """
    __tracebackhide__ = True

    # instance validation
    _check_isinstance(left, right, Series)

    if check_series_type:
        # ToDo: There are some tests using rhs is sparse
        # lhs is dense. Should use assert_class_equal in future
        assert isinstance(left, type(right))
        # assert_class_equal(left, right, obj=obj)

    # length comparison
    if len(left) != len(right):
        msg1 = f"{len(left)}, {left.index}"
        msg2 = f"{len(right)}, {right.index}"
        raise_assert_detail(obj, "Series length are different", msg1, msg2)

    # index comparison
    assert_index_equal(
        left.index,
        right.index,
        exact=check_index_type,
        check_names=check_names,
        check_less_precise=check_less_precise,
        check_exact=check_exact,
        check_categorical=check_categorical,
        obj=f"{obj}.index",
    )

    if check_dtype:
        # We want to skip exact dtype checking when `check_categorical`
        # is False. We'll still raise if only one is a `Categorical`,
        # regardless of `check_categorical`
        if (
            is_categorical_dtype(left)
            and is_categorical_dtype(right)
            and not check_categorical
        ):
            pass
        else:
            assert_attr_equal("dtype", left, right, obj=f"Attributes of {obj}")

    if check_exact:
        assert_numpy_array_equal(
            left._internal_get_values(),
            right._internal_get_values(),
            check_dtype=check_dtype,
            obj=str(obj),
        )
    elif check_datetimelike_compat:
        # we want to check only if we have compat dtypes
        # e.g. integer and M|m are NOT compat, but we can simply check
        # the values in that case
        if needs_i8_conversion(left) or needs_i8_conversion(right):

            # datetimelike may have different objects (e.g. datetime.datetime
            # vs Timestamp) but will compare equal
            if not Index(left.values).equals(Index(right.values)):
                msg = (
                    f"[datetimelike_compat=True] {left.values} "
                    f"is not equal to {right.values}."
                )
                raise AssertionError(msg)
        else:
            assert_numpy_array_equal(
                left._internal_get_values(),
                right._internal_get_values(),
                check_dtype=check_dtype,
            )
    elif is_interval_dtype(left) or is_interval_dtype(right):
        assert_interval_array_equal(left.array, right.array)
    elif is_extension_array_dtype(left.dtype) and is_datetime64tz_dtype(left.dtype):
        # .values is an ndarray, but ._values is the ExtensionArray.
        # TODO: Use .array
        assert is_extension_array_dtype(right.dtype)
        assert_extension_array_equal(left._values, right._values)
    elif (
        is_extension_array_dtype(left)
        and not is_categorical_dtype(left)
        and is_extension_array_dtype(right)
        and not is_categorical_dtype(right)
    ):
        assert_extension_array_equal(left.array, right.array)
    else:
        _testing.assert_almost_equal(
            left._internal_get_values(),
            right._internal_get_values(),
            check_less_precise=check_less_precise,
            check_dtype=check_dtype,
            obj=str(obj),
        )

    # metadata comparison
    if check_names:
        assert_attr_equal("name", left, right, obj=obj)

    if check_categorical:
        if is_categorical_dtype(left) or is_categorical_dtype(right):
            assert_categorical_equal(
                left.values,
                right.values,
                obj=f"{obj} category",
                check_category_order=check_category_order,
            )


# This could be refactored to use the NDFrame.equals method
def assert_frame_equal(
    left,
    right,
    check_dtype=True,
    check_index_type="equiv",
    check_column_type="equiv",
    check_frame_type=True,
    check_less_precise=False,
    check_names=True,
    by_blocks=False,
    check_exact=False,
    check_datetimelike_compat=False,
    check_categorical=True,
    check_like=False,
    obj="DataFrame",
):
    """
    Check that left and right DataFrame are equal.

    This function is intended to compare two DataFrames and output any
    differences. Is is mostly intended for use in unit tests.
    Additional parameters allow varying the strictness of the
    equality checks performed.

    Parameters
    ----------
    left : DataFrame
        First DataFrame to compare.
    right : DataFrame
        Second DataFrame to compare.
    check_dtype : bool, default True
        Whether to check the DataFrame dtype is identical.
    check_index_type : bool or {'equiv'}, default 'equiv'
        Whether to check the Index class, dtype and inferred_type
        are identical.
    check_column_type : bool or {'equiv'}, default 'equiv'
        Whether to check the columns class, dtype and inferred_type
        are identical. Is passed as the ``exact`` argument of
        :func:`assert_index_equal`.
    check_frame_type : bool, default True
        Whether to check the DataFrame class is identical.
    check_less_precise : bool or int, default False
        Specify comparison precision. Only used when check_exact is False.
        5 digits (False) or 3 digits (True) after decimal points are compared.
        If int, then specify the digits to compare.

        When comparing two numbers, if the first number has magnitude less
        than 1e-5, we compare the two numbers directly and check whether
        they are equivalent within the specified precision. Otherwise, we
        compare the **ratio** of the second number to the first number and
        check whether it is equivalent to 1 within the specified precision.
    check_names : bool, default True
        Whether to check that the `names` attribute for both the `index`
        and `column` attributes of the DataFrame is identical.
    by_blocks : bool, default False
        Specify how to compare internal data. If False, compare by columns.
        If True, compare by blocks.
    check_exact : bool, default False
        Whether to compare number exactly.
    check_datetimelike_compat : bool, default False
        Compare datetime-like which is comparable ignoring dtype.
    check_categorical : bool, default True
        Whether to compare internal Categorical exactly.
    check_like : bool, default False
        If True, ignore the order of index & columns.
        Note: index labels must match their respective rows
        (same as in columns) - same labels must be with the same data.
    obj : str, default 'DataFrame'
        Specify object name being compared, internally used to show appropriate
        assertion message.

    See Also
    --------
    assert_series_equal : Equivalent method for asserting Series equality.
    DataFrame.equals : Check DataFrame equality.

    Examples
    --------
    This example shows comparing two DataFrames that are equal
    but with columns of differing dtypes.

    >>> from pandas._testing import assert_frame_equal
    >>> df1 = pd.DataFrame({'a': [1, 2], 'b': [3, 4]})
    >>> df2 = pd.DataFrame({'a': [1, 2], 'b': [3.0, 4.0]})

    df1 equals itself.

    >>> assert_frame_equal(df1, df1)

    df1 differs from df2 as column 'b' is of a different type.

    >>> assert_frame_equal(df1, df2)
    Traceback (most recent call last):
    ...
    AssertionError: Attributes of DataFrame.iloc[:, 1] (column name="b") are different

    Attribute "dtype" are different
    [left]:  int64
    [right]: float64

    Ignore differing dtypes in columns with check_dtype.

    >>> assert_frame_equal(df1, df2, check_dtype=False)
    """
    __tracebackhide__ = True

    # instance validation
    _check_isinstance(left, right, DataFrame)

    if check_frame_type:
        assert isinstance(left, type(right))
        # assert_class_equal(left, right, obj=obj)

    # shape comparison
    if left.shape != right.shape:
        raise_assert_detail(
            obj, f"{obj} shape mismatch", f"{repr(left.shape)}", f"{repr(right.shape)}",
        )

    if check_like:
        left, right = left.reindex_like(right), right

    # index comparison
    assert_index_equal(
        left.index,
        right.index,
        exact=check_index_type,
        check_names=check_names,
        check_less_precise=check_less_precise,
        check_exact=check_exact,
        check_categorical=check_categorical,
        obj=f"{obj}.index",
    )

    # column comparison
    assert_index_equal(
        left.columns,
        right.columns,
        exact=check_column_type,
        check_names=check_names,
        check_less_precise=check_less_precise,
        check_exact=check_exact,
        check_categorical=check_categorical,
        obj=f"{obj}.columns",
    )

    # compare by blocks
    if by_blocks:
        rblocks = right._to_dict_of_blocks()
        lblocks = left._to_dict_of_blocks()
        for dtype in list(set(list(lblocks.keys()) + list(rblocks.keys()))):
            assert dtype in lblocks
            assert dtype in rblocks
            assert_frame_equal(
                lblocks[dtype], rblocks[dtype], check_dtype=check_dtype, obj=obj
            )

    # compare by columns
    else:
        for i, col in enumerate(left.columns):
            assert col in right
            lcol = left.iloc[:, i]
            rcol = right.iloc[:, i]
            assert_series_equal(
                lcol,
                rcol,
                check_dtype=check_dtype,
                check_index_type=check_index_type,
                check_less_precise=check_less_precise,
                check_exact=check_exact,
                check_names=check_names,
                check_datetimelike_compat=check_datetimelike_compat,
                check_categorical=check_categorical,
                obj=f'{obj}.iloc[:, {i}] (column name="{col}")',
            )


def assert_equal(left, right, **kwargs):
    """
    Wrapper for tm.assert_*_equal to dispatch to the appropriate test function.

    Parameters
    ----------
    left, right : Index, Series, DataFrame, ExtensionArray, or np.ndarray
        The two items to be compared.
    **kwargs
        All keyword arguments are passed through to the underlying assert method.
    """
    __tracebackhide__ = True

    if isinstance(left, pd.Index):
        assert_index_equal(left, right, **kwargs)
    elif isinstance(left, pd.Series):
        assert_series_equal(left, right, **kwargs)
    elif isinstance(left, pd.DataFrame):
        assert_frame_equal(left, right, **kwargs)
    elif isinstance(left, IntervalArray):
        assert_interval_array_equal(left, right, **kwargs)
    elif isinstance(left, PeriodArray):
        assert_period_array_equal(left, right, **kwargs)
    elif isinstance(left, DatetimeArray):
        assert_datetime_array_equal(left, right, **kwargs)
    elif isinstance(left, TimedeltaArray):
        assert_timedelta_array_equal(left, right, **kwargs)
    elif isinstance(left, ExtensionArray):
        assert_extension_array_equal(left, right, **kwargs)
    elif isinstance(left, np.ndarray):
        assert_numpy_array_equal(left, right, **kwargs)
    elif isinstance(left, str):
        assert kwargs == {}
        assert left == right
    else:
        raise NotImplementedError(type(left))


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
# Sparse


def assert_sp_array_equal(
    left,
    right,
    check_dtype=True,
    check_kind=True,
    check_fill_value=True,
    consolidate_block_indices=False,
):
    """
    Check that the left and right SparseArray are equal.

    Parameters
    ----------
    left : SparseArray
    right : SparseArray
    check_dtype : bool, default True
        Whether to check the data dtype is identical.
    check_kind : bool, default True
        Whether to just the kind of the sparse index for each column.
    check_fill_value : bool, default True
        Whether to check that left.fill_value matches right.fill_value
    consolidate_block_indices : bool, default False
        Whether to consolidate contiguous blocks for sparse arrays with
        a BlockIndex. Some operations, e.g. concat, will end up with
        block indices that could be consolidated. Setting this to true will
        create a new BlockIndex for that array, with consolidated
        block indices.
    """
    _check_isinstance(left, right, pd.arrays.SparseArray)

    assert_numpy_array_equal(left.sp_values, right.sp_values, check_dtype=check_dtype)

    # SparseIndex comparison
    assert isinstance(left.sp_index, pd._libs.sparse.SparseIndex)
    assert isinstance(right.sp_index, pd._libs.sparse.SparseIndex)

    if not check_kind:
        left_index = left.sp_index.to_block_index()
        right_index = right.sp_index.to_block_index()
    else:
        left_index = left.sp_index
        right_index = right.sp_index

    if consolidate_block_indices and left.kind == "block":
        # we'll probably remove this hack...
        left_index = left_index.to_int_index().to_block_index()
        right_index = right_index.to_int_index().to_block_index()

    if not left_index.equals(right_index):
        raise_assert_detail(
            "SparseArray.index", "index are not equal", left_index, right_index
        )
    else:
        # Just ensure a
        pass

    if check_fill_value:
        assert_attr_equal("fill_value", left, right)
    if check_dtype:
        assert_attr_equal("dtype", left, right)
    assert_numpy_array_equal(left.to_dense(), right.to_dense(), check_dtype=check_dtype)


# -----------------------------------------------------------------------------
# Others


def assert_contains_all(iterable, dic):
    for k in iterable:
        assert k in dic, f"Did not contain item: {repr(k)}"


def assert_copy(iter1, iter2, **eql_kwargs):
    """
    iter1, iter2: iterables that produce elements
    comparable with assert_almost_equal

    Checks that the elements are equal, but not
    the same object. (Does not check that items
    in sequences are also not the same object)
    """
    for elem1, elem2 in zip(iter1, iter2):
        assert_almost_equal(elem1, elem2, **eql_kwargs)
        msg = (
            f"Expected object {repr(type(elem1))} and object {repr(type(elem2))} to be "
            "different objects, but they were the same object."
        )
        assert elem1 is not elem2, msg


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
    dr = pd.period_range(start=dt, periods=k, freq="B", name=name, **kwargs)
    return dr


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


# make series
def makeFloatSeries(name=None):
    index = makeStringIndex(N)
    return Series(randn(N), index=index, name=name)


def makeStringSeries(name=None):
    index = makeStringIndex(N)
    return Series(randn(N), index=index, name=name)


def makeObjectSeries(name=None):
    data = makeStringIndex(N)
    data = Index(data, dtype=object)
    index = makeStringIndex(N)
    return Series(data, index=index, name=name)


def getSeriesData():
    index = makeStringIndex(N)
    return {c: Series(randn(N), index=index) for c in getCols(K)}


def makeTimeSeries(nper=None, freq="B", name=None):
    if nper is None:
        nper = N
    return Series(randn(nper), index=makeDateIndex(nper, freq=freq), name=name)


def makePeriodSeries(nper=None, name=None):
    if nper is None:
        nper = N
    return Series(randn(nper), index=makePeriodIndex(nper), name=name)


def getTimeSeriesData(nper=None, freq="B"):
    return {c: makeTimeSeries(nper, freq) for c in getCols(K)}


def getPeriodData(nper=None):
    return {c: makePeriodSeries(nper) for c in getCols(K)}


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
    idx_func = dict(
        i=makeIntIndex,
        f=makeFloatIndex,
        s=makeStringIndex,
        u=makeUnicodeIndex,
        dt=makeDateIndex,
        td=makeTimedeltaIndex,
        p=makePeriodIndex,
    ).get(idx_type)
    if idx_func:
        idx = idx_func(nentries)
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
        cnt = Counter()
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


def makeMissingCustomDataframe(
    nrows,
    ncols,
    density=0.9,
    random_state=None,
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
    Parameters
    ----------
    Density : float, optional
        Float in (0, 1) that gives the percentage of non-missing numbers in
        the DataFrame.
    random_state : {np.random.RandomState, int}, optional
        Random number generator or random seed.

    See makeCustomDataframe for descriptions of the rest of the parameters.
    """
    df = makeCustomDataframe(
        nrows,
        ncols,
        c_idx_names=c_idx_names,
        r_idx_names=r_idx_names,
        c_idx_nlevels=c_idx_nlevels,
        r_idx_nlevels=r_idx_nlevels,
        data_gen_f=data_gen_f,
        c_ndupe_l=c_ndupe_l,
        r_ndupe_l=r_ndupe_l,
        dtype=dtype,
        c_idx_type=c_idx_type,
        r_idx_type=r_idx_type,
    )

    i, j = _create_missing_idx(nrows, ncols, density, random_state)
    df.values[i, j] = np.nan
    return df


def makeMissingDataframe(density=0.9, random_state=None):
    df = makeDataFrame()
    i, j = _create_missing_idx(*df.shape, density=density, random_state=random_state)
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
            args = []
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
    url="http://www.google.com",
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
        for connectivity. Defaults to 'http://www.google.com'.
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
    * ``raise_on_error`` supercedes ``check_before_test``

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

        >>> @network("http://www.yahoo.com")
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
        if check_before_test and not raise_on_error:
            if not can_connect(url, error_classes):
                skip()
        try:
            return t(*args, **kwargs)
        except Exception as err:
            errno = getattr(err, "errno", None)
            if not errno and hasattr(errno, "reason"):
                errno = getattr(err.reason, "errno", None)

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
    expected_warning=Warning,
    filter_level="always",
    clear=None,
    check_stacklevel=True,
    raise_on_extra_warnings=True,
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

    clear : str, default None
        If not ``None`` then remove any previously raised warnings from
        the ``__warningsregistry__`` to ensure that no warning messages are
        suppressed by this context manager. If ``None`` is specified,
        the ``__warningsregistry__`` keeps track of which warnings have been
        shown, and does not show them again.
    check_stacklevel : bool, default True
        If True, displays the line that called the function containing
        the warning to show were the function is called. Otherwise, the
        line that implements the function is displayed.
    raise_on_extra_warnings : bool, default True
        Whether extra warnings not of the type `expected_warning` should
        cause the test to fail.

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

        if clear is not None:
            # make sure that we are clearing these warnings
            # if they have happened before
            # to guarantee that we will catch them
            if not is_list_like(clear):
                clear = [clear]
            for m in clear:
                try:
                    m.__warningregistry__.clear()
                except AttributeError:
                    # module may not have __warningregistry__
                    pass

        saw_warning = False
        warnings.simplefilter(filter_level)
        yield w
        extra_warnings = []

        for actual_warning in w:
            if expected_warning and issubclass(
                actual_warning.category, expected_warning
            ):
                saw_warning = True

                if check_stacklevel and issubclass(
                    actual_warning.category, (FutureWarning, DeprecationWarning)
                ):
                    from inspect import getframeinfo, stack

                    caller = getframeinfo(stack()[2][0])
                    msg = (
                        "Warning not set with correct stacklevel. "
                        f"File where warning is raised: {actual_warning.filename} != "
                        f"{caller.filename}. Warning message: {actual_warning.message}"
                    )
                    assert actual_warning.filename == caller.filename, msg
            else:
                extra_warnings.append(
                    (
                        actual_warning.category.__name__,
                        actual_warning.message,
                        actual_warning.filename,
                        actual_warning.lineno,
                    )
                )
        if expected_warning:
            msg = (
                f"Did not see expected warning of class "
                f"{repr(expected_warning.__name__)}"
            )
            assert saw_warning, msg
        if raise_on_extra_warnings and extra_warnings:
            raise AssertionError(
                f"Caused unexpected warning(s): {repr(extra_warnings)}"
            )


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


@contextmanager
def with_csv_dialect(name, **kwargs):
    """
    Context manager to temporarily register a CSV dialect for parsing CSV.

    Parameters
    ----------
    name : str
        The name of the dialect.
    kwargs : mapping
        The parameters for the dialect.

    Raises
    ------
    ValueError : the name of the dialect conflicts with a builtin one.

    See Also
    --------
    csv : Python's CSV library.
    """
    import csv

    _BUILTIN_DIALECTS = {"excel", "excel-tab", "unix"}

    if name in _BUILTIN_DIALECTS:
        raise ValueError("Cannot override builtin dialect.")

    csv.register_dialect(name, **kwargs)
    yield
    csv.unregister_dialect(name)


@contextmanager
def use_numexpr(use, min_elements=None):
    from pandas.core.computation import expressions as expr

    if min_elements is None:
        min_elements = expr._MIN_ELEMENTS

    olduse = expr._USE_NUMEXPR
    oldmin = expr._MIN_ELEMENTS
    expr.set_use_numexpr(use)
    expr._MIN_ELEMENTS = min_elements
    yield
    expr._MIN_ELEMENTS = oldmin
    expr.set_use_numexpr(olduse)


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


@contextmanager
def set_timezone(tz: str):
    """
    Context manager for temporarily setting a timezone.

    Parameters
    ----------
    tz : str
        A string representing a valid timezone.

    Examples
    --------
    >>> from datetime import datetime
    >>> from dateutil.tz import tzlocal
    >>> tzlocal().tzname(datetime.now())
    'IST'

    >>> with set_timezone('US/Eastern'):
    ...     tzlocal().tzname(datetime.now())
    ...
    'EDT'
    """
    import os
    import time

    def setTZ(tz):
        if tz is None:
            try:
                del os.environ["TZ"]
            except KeyError:
                pass
        else:
            os.environ["TZ"] = tz
            time.tzset()

    orig_tz = os.environ.get("TZ")
    setTZ(tz)
    try:
        yield
    finally:
        setTZ(orig_tz)


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
