from __future__ import division

from contextlib import contextmanager
from functools import wraps
import locale
import os
import re
from shutil import rmtree
import string
import subprocess
import sys
import tempfile
import traceback
import warnings

import numpy as np
from numpy.random import rand

import pandas.compat as compat
from pandas.compat import (
    PY2, PY3, callable, filter, httplib, lrange, map, range, u, unichr, zip)

from pandas.core.dtypes.common import (
    is_datetime64_dtype, is_datetime64tz_dtype, is_period_dtype,
    is_timedelta64_dtype)

import pandas as pd
from pandas import Categorical, DataFrame, Index, Series
from pandas.core.arrays import (
    DatetimeArrayMixin as DatetimeArray, PeriodArray,
    TimedeltaArrayMixin as TimedeltaArray, period_array)
import pandas.core.common as com

from pandas.io.common import urlopen

from .asserters import (  # noqa:F401
    assert_almost_equal, assert_attr_equal, assert_categorical_equal,
    assert_class_equal, assert_datetime_array_equal, assert_dict_equal,
    assert_equal, assert_extension_array_equal, assert_frame_equal,
    assert_index_equal, assert_interval_array_equal,
    assert_is_valid_plot_return_object, assert_numpy_array_equal,
    assert_panel_equal, assert_period_array_equal, assert_produces_warning,
    assert_raises_regex, assert_series_equal, assert_sp_array_equal,
    assert_sp_frame_equal, assert_sp_series_equal,
    assert_timedelta_array_equal, raise_assert_detail)
from .strategies import (  # noqa:F401
    getMixedTypeDict, getPeriodData, getSeriesData, getTimeSeriesData,
    makeBoolIndex, makeCategoricalIndex, makeCustomDataframe, makeCustomIndex,
    makeDataFrame, makeDateIndex, makeFloatIndex, makeFloatSeries,
    makeIntervalIndex, makeIntIndex, makeMissingCustomDataframe,
    makeMissingDataframe, makeMixedDataFrame, makeMultiIndex, makeObjectSeries,
    makePanel, makePeriodFrame, makePeriodIndex, makePeriodPanel,
    makePeriodSeries, makeRangeIndex, makeStringIndex, makeStringSeries,
    makeTimeDataFrame, makeTimedeltaIndex, makeTimeSeries, makeUIntIndex,
    makeUnicodeIndex)

N = 30
K = 4
_RAISE_NETWORK_ERROR_DEFAULT = False

# set testing_mode
_testing_mode_warnings = (DeprecationWarning, compat.ResourceWarning)


def set_testing_mode():
    # set the testing mode filters
    testing_mode = os.environ.get('PANDAS_TESTING_MODE', 'None')
    if 'deprecate' in testing_mode:
        warnings.simplefilter('always', _testing_mode_warnings)


def reset_testing_mode():
    # reset the testing mode filters
    testing_mode = os.environ.get('PANDAS_TESTING_MODE', 'None')
    if 'deprecate' in testing_mode:
        warnings.simplefilter('ignore', _testing_mode_warnings)


set_testing_mode()


def reset_display_options():
    """
    Reset the display options for printing and representing objects.
    """

    pd.reset_option('^display.', silent=True)


def round_trip_pickle(obj, path=None):
    """
    Pickle an object and then read it again.

    Parameters
    ----------
    obj : pandas object
        The object to pickle and then re-read.
    path : str, default None
        The path where the pickled object is written and then read.

    Returns
    -------
    round_trip_pickled_object : pandas object
        The original object that was pickled and then re-read.
    """

    if path is None:
        path = u('__{random_bytes}__.pickle'.format(random_bytes=rands(10)))
    with ensure_clean(path) as path:
        pd.to_pickle(obj, path)
        return pd.read_pickle(path)


def round_trip_pathlib(writer, reader, path=None):
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
    round_trip_object : pandas object
        The original object that was serialized and then re-read.
    """

    import pytest
    Path = pytest.importorskip('pathlib').Path
    if path is None:
        path = '___pathlib___'
    with ensure_clean(path) as path:
        writer(Path(path))
        obj = reader(Path(path))
    return obj


def round_trip_localpath(writer, reader, path=None):
    """
    Write an object to file specified by a py.path LocalPath and read it back

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
    round_trip_object : pandas object
        The original object that was serialized and then re-read.
    """
    import pytest
    LocalPath = pytest.importorskip('py.path').local
    if path is None:
        path = '___localpath___'
    with ensure_clean(path) as path:
        writer(LocalPath(path))
        obj = reader(LocalPath(path))
    return obj


@contextmanager
def decompress_file(path, compression):
    """
    Open a compressed file and return a file object

    Parameters
    ----------
    path : str
        The path where the file is read from

    compression : {'gzip', 'bz2', 'zip', 'xz', None}
        Name of the decompression to use

    Returns
    -------
    f : file object
    """

    if compression is None:
        f = open(path, 'rb')
    elif compression == 'gzip':
        import gzip
        f = gzip.open(path, 'rb')
    elif compression == 'bz2':
        import bz2
        f = bz2.BZ2File(path, 'rb')
    elif compression == 'xz':
        lzma = compat.import_lzma()
        f = lzma.LZMAFile(path, 'rb')
    elif compression == 'zip':
        import zipfile
        zip_file = zipfile.ZipFile(path)
        zip_names = zip_file.namelist()
        if len(zip_names) == 1:
            f = zip_file.open(zip_names.pop())
        else:
            raise ValueError('ZIP file {} error. Only one file per ZIP.'
                             .format(path))
    else:
        msg = 'Unrecognized compression type: {}'.format(compression)
        raise ValueError(msg)

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
        lzma = compat.import_lzma()
        compress_method = lzma.LZMAFile
    else:
        msg = "Unrecognized compression type: {}".format(compression)
        raise ValueError(msg)

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


def randbool(size=(), p=0.5):
    return rand(*size) <= p


RANDS_CHARS = np.array(list(string.ascii_letters + string.digits),
                       dtype=(np.str_, 1))
RANDU_CHARS = np.array(list(u("").join(map(unichr, lrange(1488, 1488 + 26))) +
                            string.digits), dtype=(np.unicode_, 1))


def rands_array(nchars, size, dtype='O'):
    """Generate an array of byte strings."""
    retval = (np.random.choice(RANDS_CHARS, size=nchars * np.prod(size))
              .view((np.str_, nchars)).reshape(size))
    if dtype is None:
        return retval
    else:
        return retval.astype(dtype)


def randu_array(nchars, size, dtype='O'):
    """Generate an array of unicode strings."""
    retval = (np.random.choice(RANDU_CHARS, size=nchars * np.prod(size))
              .view((np.unicode_, nchars)).reshape(size))
    if dtype is None:
        return retval
    else:
        return retval.astype(dtype)


def rands(nchars):
    """
    Generate one random byte string.

    See `rands_array` if you want to create an array of random strings.

    """
    return ''.join(np.random.choice(RANDS_CHARS, nchars))


def randu(nchars):
    """
    Generate one random unicode string.

    See `randu_array` if you want to create an array of random unicode strings.

    """
    return ''.join(np.random.choice(RANDU_CHARS, nchars))


def close(fignum=None):
    from matplotlib.pyplot import get_fignums, close as _close

    if fignum is None:
        for fignum in get_fignums():
            _close(fignum)
    else:
        _close(fignum)


# -----------------------------------------------------------------------------
# locale utilities


def check_output(*popenargs, **kwargs):
    # shamelessly taken from Python 2.7 source
    r"""Run command with arguments and return its output as a byte string.

    If the exit code was non-zero it raises a CalledProcessError.  The
    CalledProcessError object will have the return code in the returncode
    attribute and output in the output attribute.

    The arguments are the same as for the Popen constructor.  Example:

    >>> check_output(["ls", "-l", "/dev/null"])
    'crw-rw-rw- 1 root root 1, 3 Oct 18  2007 /dev/null\n'

    The stdout argument is not allowed as it is used internally.
    To capture standard error in the result, use stderr=STDOUT.

    >>> check_output(["/bin/sh", "-c",
    ...               "ls -l non_existent_file ; exit 0"],
    ...              stderr=STDOUT)
    'ls: non_existent_file: No such file or directory\n'
    """
    if 'stdout' in kwargs:
        raise ValueError('stdout argument not allowed, it will be overridden.')
    process = subprocess.Popen(stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               *popenargs, **kwargs)
    output, unused_err = process.communicate()
    retcode = process.poll()
    if retcode:
        cmd = kwargs.get("args")
        if cmd is None:
            cmd = popenargs[0]
        raise subprocess.CalledProcessError(retcode, cmd, output=output)
    return output


def _default_locale_getter():
    try:
        raw_locales = check_output(['locale -a'], shell=True)
    except subprocess.CalledProcessError as e:
        raise type(e)("{exception}, the 'locale -a' command cannot be found "
                      "on your system".format(exception=e))
    return raw_locales


def get_locales(prefix=None, normalize=True,
                locale_getter=_default_locale_getter):
    """Get all the locales that are available on the system.

    Parameters
    ----------
    prefix : str
        If not ``None`` then return only those locales with the prefix
        provided. For example to get all English language locales (those that
        start with ``"en"``), pass ``prefix="en"``.
    normalize : bool
        Call ``locale.normalize`` on the resulting list of available locales.
        If ``True``, only locales that can be set without throwing an
        ``Exception`` are returned.
    locale_getter : callable
        The function to use to retrieve the current locales. This should return
        a string with each locale separated by a newline character.

    Returns
    -------
    locales : list of strings
        A list of locale strings that can be set with ``locale.setlocale()``.
        For example::

            locale.setlocale(locale.LC_ALL, locale_string)

    On error will return None (no locale available, e.g. Windows)

    """
    try:
        raw_locales = locale_getter()
    except Exception:
        return None

    try:
        # raw_locales is "\n" separated list of locales
        # it may contain non-decodable parts, so split
        # extract what we can and then rejoin.
        raw_locales = raw_locales.split(b'\n')
        out_locales = []
        for x in raw_locales:
            if PY3:
                out_locales.append(str(
                    x, encoding=pd.options.display.encoding))
            else:
                out_locales.append(str(x))

    except TypeError:
        pass

    if prefix is None:
        return _valid_locales(out_locales, normalize)

    pattern = re.compile('{prefix}.*'.format(prefix=prefix))
    found = pattern.findall('\n'.join(out_locales))
    return _valid_locales(found, normalize)


@contextmanager
def set_locale(new_locale, lc_var=locale.LC_ALL):
    """Context manager for temporarily setting a locale.

    Parameters
    ----------
    new_locale : str or tuple
        A string of the form <language_country>.<encoding>. For example to set
        the current locale to US English with a UTF8 encoding, you would pass
        "en_US.UTF-8".
    lc_var : int, default `locale.LC_ALL`
        The category of the locale being set.

    Notes
    -----
    This is useful when you want to run a particular block of code under a
    particular locale, without globally setting the locale. This probably isn't
    thread-safe.
    """
    current_locale = locale.getlocale()

    try:
        locale.setlocale(lc_var, new_locale)
        normalized_locale = locale.getlocale()
        if com._all_not_none(*normalized_locale):
            yield '.'.join(normalized_locale)
        else:
            yield new_locale
    finally:
        locale.setlocale(lc_var, current_locale)


def can_set_locale(lc, lc_var=locale.LC_ALL):
    """
    Check to see if we can set a locale, and subsequently get the locale,
    without raising an Exception.

    Parameters
    ----------
    lc : str
        The locale to attempt to set.
    lc_var : int, default `locale.LC_ALL`
        The category of the locale being set.

    Returns
    -------
    is_valid : bool
        Whether the passed locale can be set
    """

    try:
        with set_locale(lc, lc_var=lc_var):
            pass
    except (ValueError,
            locale.Error):  # horrible name for a Exception subclass
        return False
    else:
        return True


def _valid_locales(locales, normalize):
    """Return a list of normalized locales that do not throw an ``Exception``
    when set.

    Parameters
    ----------
    locales : str
        A string where each locale is separated by a newline.
    normalize : bool
        Whether to call ``locale.normalize`` on each locale.

    Returns
    -------
    valid_locales : list
        A list of valid locales.
    """
    if normalize:
        normalizer = lambda x: locale.normalize(x.strip())
    else:
        normalizer = lambda x: x.strip()

    return list(filter(can_set_locale, map(normalizer, locales)))

# -----------------------------------------------------------------------------
# Stdout / stderr decorators


@contextmanager
def set_defaultencoding(encoding):
    """
    Set default encoding (as given by sys.getdefaultencoding()) to the given
    encoding; restore on exit.

    Parameters
    ----------
    encoding : str
    """
    if not PY2:
        raise ValueError("set_defaultencoding context is only available "
                         "in Python 2.")
    orig = sys.getdefaultencoding()
    reload(sys)  # noqa:F821
    sys.setdefaultencoding(encoding)
    try:
        yield
    finally:
        sys.setdefaultencoding(orig)


# -----------------------------------------------------------------------------
# Console debugging tools


def debug(f, *args, **kwargs):
    from pdb import Pdb as OldPdb
    try:
        from IPython.core.debugger import Pdb
        kw = dict(color_scheme='Linux')
    except ImportError:
        Pdb = OldPdb
        kw = {}
    pdb = Pdb(**kw)
    return pdb.runcall(f, *args, **kwargs)


def pudebug(f, *args, **kwargs):
    import pudb
    return pudb.runcall(f, *args, **kwargs)


def set_trace():
    from IPython.core.debugger import Pdb
    try:
        Pdb(color_scheme='Linux').set_trace(sys._getframe().f_back)
    except Exception:
        from pdb import Pdb as OldPdb
        OldPdb().set_trace(sys._getframe().f_back)

# -----------------------------------------------------------------------------
# contextmanager to ensure the file cleanup


@contextmanager
def ensure_clean(filename=None, return_filelike=False):
    """Gets a temporary path and agrees to remove on close.

    Parameters
    ----------
    filename : str (optional)
        if None, creates a temporary file which is then removed when out of
        scope. if passed, creates temporary file with filename as ending.
    return_filelike : bool (default False)
        if True, returns a file-like which is *always* cleaned. Necessary for
        savefig and other functions which want to append extensions.
    """
    filename = filename or ''
    fd = None

    if return_filelike:
        f = tempfile.TemporaryFile(suffix=filename)
        try:
            yield f
        finally:
            f.close()
    else:
        # don't generate tempfile if using a path with directory specified
        if len(os.path.dirname(filename)):
            raise ValueError("Can't pass a qualified name to ensure_clean()")

        try:
            fd, filename = tempfile.mkstemp(suffix=filename)
        except UnicodeEncodeError:
            import pytest
            pytest.skip('no unicode file names on this system')

        try:
            yield filename
        finally:
            try:
                os.close(fd)
            except Exception:
                print("Couldn't close file descriptor: {fdesc} (file: {fname})"
                      .format(fdesc=fd, fname=filename))
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception as e:
                print("Exception on removing file: {error}".format(error=e))


@contextmanager
def ensure_clean_dir():
    """
    Get a temporary directory path and agrees to remove on close.

    Yields
    ------
    Temporary directory path
    """
    directory_name = tempfile.mkdtemp(suffix='')
    try:
        yield directory_name
    finally:
        try:
            rmtree(directory_name)
        except Exception:
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


def equalContents(arr1, arr2):
    """Checks if the set of unique elements of arr1 and arr2 are equivalent.
    """
    return frozenset(arr1) == frozenset(arr2)


def isiterable(obj):
    return hasattr(obj, '__iter__')


def is_sorted(seq):
    if isinstance(seq, (Index, Series)):
        seq = seq.values
    # sorting does not change precisions
    return assert_numpy_array_equal(seq, np.sort(np.array(seq)))


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


def assert_contains_all(iterable, dic):
    for k in iterable:
        assert k in dic, "Did not contain item: '{key!r}'".format(key=k)


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
        msg = ("Expected object {obj1!r} and object {obj2!r} to be "
               "different objects, but they were the same object."
               ).format(obj1=type(elem1), obj2=type(elem2))
        assert elem1 is not elem2, msg


def getArangeMat():
    return np.arange(N * K).reshape((N, K))


def all_index_generator(k=10):
    """Generator which can be iterated over to get instances of all the various
    index classes.

    Parameters
    ----------
    k: length of each of the index instances
    """
    all_make_index_funcs = [makeIntIndex, makeFloatIndex, makeStringIndex,
                            makeUnicodeIndex, makeDateIndex, makePeriodIndex,
                            makeTimedeltaIndex, makeBoolIndex, makeRangeIndex,
                            makeIntervalIndex,
                            makeCategoricalIndex]
    for make_index_func in all_make_index_funcs:
        yield make_index_func(k=k)


def index_subclass_makers_generator():
    make_index_funcs = [
        makeDateIndex, makePeriodIndex,
        makeTimedeltaIndex, makeRangeIndex,
        makeIntervalIndex, makeCategoricalIndex,
        makeMultiIndex
    ]
    for make_index_func in make_index_funcs:
        yield make_index_func


def all_timeseries_index_generator(k=10):
    """Generator which can be iterated over to get instances of all the classes
    which represent time-seires.

    Parameters
    ----------
    k: length of each of the index instances
    """
    make_index_funcs = [makeDateIndex, makePeriodIndex, makeTimedeltaIndex]
    for make_index_func in make_index_funcs:
        yield make_index_func(k=k)


def add_nans(panel):
    I, J, N = panel.shape
    for i, item in enumerate(panel.items):
        dm = panel[item]
        for j, col in enumerate(dm.columns):
            dm[col][:i + j] = np.NaN
    return panel


def add_nans_panel4d(panel4d):
    for l, label in enumerate(panel4d.labels):
        panel = panel4d[label]
        add_nans(panel)
    return panel4d


class TestSubDict(dict):

    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)


def optional_args(decorator):
    """allows a decorator to take optional positional and keyword arguments.
    Assumes that taking a single, callable, positional argument means that
    it is decorating a function, i.e. something like this::

        @my_decorator
        def function(): pass

    Calls decorator with decorator(f, *args, **kwargs)"""

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
    'timed out',
    'Server Hangup',
    'HTTP Error 503: Service Unavailable',
    '502: Proxy Error',
    'HTTP Error 502: internal error',
    'HTTP Error 502',
    'HTTP Error 503',
    'HTTP Error 403',
    'HTTP Error 400',
    'Temporary failure in name resolution',
    'Name or service not known',
    'Connection refused',
    'certificate verify',
)

# or this e.errno/e.reason.errno
_network_errno_vals = (
    101,  # Network is unreachable
    111,  # Connection refused
    110,  # Connection timed out
    104,  # Connection reset Error
    54,   # Connection reset by peer
    60,   # urllib.error.URLError: [Errno 60] Connection timed out
)

# Both of the above shouldn't mask real issues such as 404's
# or refused connections (changed DNS).
# But some tests (test_data yahoo) contact incredibly flakey
# servers.

# and conditionally raise on these exception types
_network_error_classes = (IOError, httplib.HTTPException)

if PY3:
    _network_error_classes += (TimeoutError,)  # noqa


def can_connect(url, error_classes=_network_error_classes):
    """Try to connect to the given url. True if succeeds, False if IOError
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
    try:
        with urlopen(url):
            pass
    except error_classes:
        return False
    else:
        return True


@optional_args
def network(t, url="http://www.google.com",
            raise_on_error=_RAISE_NETWORK_ERROR_DEFAULT,
            check_before_test=False,
            error_classes=_network_error_classes,
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

      >>> from pandas.util.testing import network
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
    t.network = True

    @compat.wraps(t)
    def wrapper(*args, **kwargs):
        if check_before_test and not raise_on_error:
            if not can_connect(url, error_classes):
                skip()
        try:
            return t(*args, **kwargs)
        except Exception as e:
            errno = getattr(e, 'errno', None)
            if not errno and hasattr(errno, "reason"):
                errno = getattr(e.reason, 'errno', None)

            if errno in skip_errnos:
                skip("Skipping test due to known errno"
                     " and error {error}".format(error=e))

            try:
                e_str = traceback.format_exc(e)
            except Exception:
                e_str = str(e)

            if any(m.lower() in e_str.lower() for m in _skip_on_messages):
                skip("Skipping test because exception "
                     "message is known and error {error}".format(error=e))

            if not isinstance(e, error_classes):
                raise

            if raise_on_error or can_connect(url, error_classes):
                raise
            else:
                skip("Skipping test due to lack of connectivity"
                     " and error {error}".format(error=e))

    return wrapper


with_connectivity_check = network


class RNGContext(object):
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
    """Decorator to run the same function multiple times in parallel.

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
                thread = threading.Thread(target=func, args=args,
                                          kwargs=updated_kwargs)
                threads.append(thread)
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join()
        return inner
    return wrapper


class SubclassedSeries(Series):
    _metadata = ['testattr', 'name']

    @property
    def _constructor(self):
        return SubclassedSeries

    @property
    def _constructor_expanddim(self):
        return SubclassedDataFrame


class SubclassedDataFrame(DataFrame):
    _metadata = ['testattr']

    @property
    def _constructor(self):
        return SubclassedDataFrame

    @property
    def _constructor_sliced(self):
        return SubclassedSeries


class SubclassedSparseSeries(pd.SparseSeries):
    _metadata = ['testattr']

    @property
    def _constructor(self):
        return SubclassedSparseSeries

    @property
    def _constructor_expanddim(self):
        return SubclassedSparseDataFrame


class SubclassedSparseDataFrame(pd.SparseDataFrame):
    _metadata = ['testattr']

    @property
    def _constructor(self):
        return SubclassedSparseDataFrame

    @property
    def _constructor_sliced(self):
        return SubclassedSparseSeries


class SubclassedCategorical(Categorical):

    @property
    def _constructor(self):
        return SubclassedCategorical


@contextmanager
def set_timezone(tz):
    """Context manager for temporarily setting a timezone.

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
                del os.environ['TZ']
            except KeyError:
                pass
        else:
            os.environ['TZ'] = tz
            time.tzset()

    orig_tz = os.environ.get('TZ')
    setTZ(tz)
    try:
        yield
    finally:
        setTZ(orig_tz)


def _make_skipna_wrapper(alternative, skipna_alternative=None):
    """Create a function for calling on an array.

    Parameters
    ----------
    alternative : function
        The function to be called on the array with no NaNs.
        Only used when 'skipna_alternative' is None.
    skipna_alternative : function
        The function to be called on the original array

    Returns
    -------
    skipna_wrapper : function
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


def convert_rows_list_to_csv_str(rows_list):
    """
    Convert list of CSV rows to single CSV-formatted string for current OS.

    This method is used for creating expected value of to_csv() method.

    Parameters
    ----------
    rows_list : list
        The list of string. Each element represents the row of csv.

    Returns
    -------
    expected : string
        Expected output of to_csv() in current OS
    """
    sep = os.linesep
    expected = sep.join(rows_list) + sep
    return expected
