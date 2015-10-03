from __future__ import division
# pylint: disable-msg=W0402

import random
import re
import string
import sys
import tempfile
import warnings
import inspect
import os
import subprocess
import locale
import unittest
import traceback

from datetime import datetime
from functools import wraps, partial
from contextlib import contextmanager
from distutils.version import LooseVersion

from numpy.random import randn, rand
import numpy as np

import pandas as pd
from pandas.core.common import (is_sequence, array_equivalent, is_list_like, is_number,
                                is_datetimelike_v_numeric, is_datetimelike_v_object,
                                is_number, pprint_thing, take_1d,
                                needs_i8_conversion)

import pandas.compat as compat
from pandas.compat import(
    filter, map, zip, range, unichr, lrange, lmap, lzip, u, callable, Counter,
    raise_with_traceback, httplib, is_platform_windows, is_platform_32bit
)

from pandas.computation import expressions as expr

from pandas import (bdate_range, CategoricalIndex, DatetimeIndex, TimedeltaIndex, PeriodIndex,
                    Index, MultiIndex, Series, DataFrame, Panel, Panel4D)
from pandas.util.decorators import deprecate
from pandas import _testing
from pandas.io.common import urlopen

N = 30
K = 4
_RAISE_NETWORK_ERROR_DEFAULT = False

# set testing_mode
def set_testing_mode():
    # set the testing mode filters
    testing_mode = os.environ.get('PANDAS_TESTING_MODE','None')
    if 'deprecate' in testing_mode:
        warnings.simplefilter('always', DeprecationWarning)

def reset_testing_mode():
    # reset the testing mode filters
    testing_mode = os.environ.get('PANDAS_TESTING_MODE','None')
    if 'deprecate' in testing_mode:
        warnings.simplefilter('ignore', DeprecationWarning)


set_testing_mode()

class TestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pd.set_option('chained_assignment', 'raise')

    @classmethod
    def tearDownClass(cls):
        pass

    def reset_display_options(self):
        # reset the display options
        pd.reset_option('^display.', silent=True)

    def round_trip_pickle(self, obj, path=None):
        if path is None:
            path = u('__%s__.pickle' % rands(10))
        with ensure_clean(path) as path:
            pd.to_pickle(obj, path)
            return pd.read_pickle(path)

    # https://docs.python.org/3/library/unittest.html#deprecated-aliases
    def assertEquals(self, *args, **kwargs):
        return deprecate('assertEquals', self.assertEqual)(*args, **kwargs)

    def assertNotEquals(self, *args, **kwargs):
        return deprecate('assertNotEquals', self.assertNotEqual)(*args, **kwargs)

    def assert_(self, *args, **kwargs):
        return deprecate('assert_', self.assertTrue)(*args, **kwargs)

    def assertAlmostEquals(self, *args, **kwargs):
        return deprecate('assertAlmostEquals', self.assertAlmostEqual)(*args, **kwargs)

    def assertNotAlmostEquals(self, *args, **kwargs):
        return deprecate('assertNotAlmostEquals', self.assertNotAlmostEqual)(*args, **kwargs)


# NOTE: don't pass an NDFrame or index to this function - may not handle it
# well.
assert_almost_equal = _testing.assert_almost_equal

assert_dict_equal = _testing.assert_dict_equal

def randbool(size=(), p=0.5):
    return rand(*size) <= p


RANDS_CHARS = np.array(list(string.ascii_letters + string.digits),
                       dtype=(np.str_, 1))
RANDU_CHARS = np.array(list(u("").join(map(unichr, lrange(1488, 1488 + 26))) +
                            string.digits), dtype=(np.unicode_, 1))


def rands_array(nchars, size, dtype='O'):
    """Generate an array of byte strings."""
    retval = (choice(RANDS_CHARS, size=nchars * np.prod(size))
              .view((np.str_, nchars)).reshape(size))
    if dtype is None:
        return retval
    else:
        return retval.astype(dtype)


def randu_array(nchars, size, dtype='O'):
    """Generate an array of unicode strings."""
    retval = (choice(RANDU_CHARS, size=nchars * np.prod(size))
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
    return ''.join(choice(RANDS_CHARS, nchars))


def randu(nchars):
    """
    Generate one random unicode string.

    See `randu_array` if you want to create an array of random unicode strings.

    """
    return ''.join(choice(RANDU_CHARS, nchars))


def choice(x, size=10):
    """sample with replacement; uniform over the input"""
    try:
        return np.random.choice(x, size=size)
    except AttributeError:
        return np.random.randint(len(x), size=size).choose(x)


def close(fignum=None):
    from matplotlib.pyplot import get_fignums, close as _close

    if fignum is None:
        for fignum in get_fignums():
            _close(fignum)
    else:
        _close(fignum)


def _skip_if_32bit():
    import nose
    if is_platform_32bit():
        raise nose.SkipTest("skipping for 32 bit")

def mplskip(cls):
    """Skip a TestCase instance if matplotlib isn't installed"""

    @classmethod
    def setUpClass(cls):
        try:
            import matplotlib as mpl
            mpl.use("Agg", warn=False)
        except ImportError:
            import nose
            raise nose.SkipTest("matplotlib not installed")

    cls.setUpClass = setUpClass
    return cls


def _skip_if_mpl_1_5():
    import matplotlib
    v = matplotlib.__version__
    if v > LooseVersion('1.4.3') or v[0] == '0':
        import nose
        raise nose.SkipTest("matplotlib 1.5")


def _skip_if_no_scipy():
    try:
        import scipy.stats
    except ImportError:
        import nose
        raise nose.SkipTest("no scipy.stats module")
    try:
        import scipy.interpolate
    except ImportError:
        import nose
        raise nose.SkipTest('scipy.interpolate missing')


def _skip_if_no_pytz():
    try:
        import pytz
    except ImportError:
        import nose
        raise nose.SkipTest("pytz not installed")


def _skip_if_no_dateutil():
    try:
        import dateutil
    except ImportError:
        import nose
        raise nose.SkipTest("dateutil not installed")


def _skip_if_windows_python_3():
    if compat.PY3 and is_platform_windows():
        import nose
        raise nose.SkipTest("not used on python 3/win32")

def _skip_if_windows():
    if is_platform_windows():
        import nose
        raise nose.SkipTest("Running on Windows")


def _skip_if_no_cday():
    from pandas.core.datetools import cday
    if cday is None:
        import nose
        raise nose.SkipTest("CustomBusinessDay not available.")


def _skip_if_python26():
    if sys.version_info[:2] == (2, 6):
        import nose
        raise nose.SkipTest("skipping on python2.6")

def _incompat_bottleneck_version(method):
    """ skip if we have bottleneck installed
    and its >= 1.0
    as we don't match the nansum/nanprod behavior for all-nan
    ops, see GH9422
    """
    if method not in ['sum','prod']:
        return False
    try:
        import bottleneck as bn
        return bn.__version__ >= LooseVersion('1.0')
    except ImportError:
        return False

def skip_if_no_ne(engine='numexpr'):
    import nose
    _USE_NUMEXPR = pd.computation.expressions._USE_NUMEXPR

    if engine == 'numexpr':
        try:
            import numexpr as ne
        except ImportError:
            raise nose.SkipTest("numexpr not installed")

        if not _USE_NUMEXPR:
            raise nose.SkipTest("numexpr disabled")

        if ne.__version__ < LooseVersion('2.0'):
            raise nose.SkipTest("numexpr version too low: "
                                "%s" % ne.__version__)



#------------------------------------------------------------------------------
# locale utilities

def check_output(*popenargs, **kwargs):  # shamelessly taken from Python 2.7 source
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
    process = subprocess.Popen(stdout=subprocess.PIPE,stderr=subprocess.PIPE,
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
        raise type(e)("%s, the 'locale -a' command cannot be found on your "
                      "system" % e)
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
    except:
        return None

    try:
        # raw_locales is "\n" seperated list of locales
        # it may contain non-decodable parts, so split
        # extract what we can and then rejoin.
        raw_locales = raw_locales.split(b'\n')
        out_locales = []
        for x in raw_locales:
            if compat.PY3:
                out_locales.append(str(x, encoding=pd.options.display.encoding))
            else:
                out_locales.append(str(x))

    except TypeError:
        pass

    if prefix is None:
        return _valid_locales(out_locales, normalize)

    found = re.compile('%s.*' % prefix).findall('\n'.join(out_locales))
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

    Notes
    -----
    This is useful when you want to run a particular block of code under a
    particular locale, without globally setting the locale. This probably isn't
    thread-safe.
    """
    current_locale = locale.getlocale()

    try:
        locale.setlocale(lc_var, new_locale)

        try:
            normalized_locale = locale.getlocale()
        except ValueError:
            yield new_locale
        else:
            if all(lc is not None for lc in normalized_locale):
                yield '.'.join(normalized_locale)
            else:
                yield new_locale
    finally:
        locale.setlocale(lc_var, current_locale)


def _can_set_locale(lc):
    """Check to see if we can set a locale without throwing an exception.

    Parameters
    ----------
    lc : str
        The locale to attempt to set.

    Returns
    -------
    isvalid : bool
        Whether the passed locale can be set
    """
    try:
        with set_locale(lc):
            pass
    except locale.Error:  # horrible name for a Exception subclass
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

    return list(filter(_can_set_locale, map(normalizer, locales)))


#------------------------------------------------------------------------------
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
    except:
        from pdb import Pdb as OldPdb
        OldPdb().set_trace(sys._getframe().f_back)

#------------------------------------------------------------------------------
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
            import nose
            raise nose.SkipTest('no unicode file names on this system')

        try:
            yield filename
        finally:
            try:
                os.close(fd)
            except Exception as e:
                print("Couldn't close file descriptor: %d (file: %s)" %
                    (fd, filename))
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception as e:
                print("Exception on removing file: %s" % e)

def get_data_path(f=''):
    """Return the path of a data file, these are relative to the current test
    directory.
    """
    # get our callers file
    _, filename, _, _, _, _ = inspect.getouterframes(inspect.currentframe())[1]
    base_dir = os.path.abspath(os.path.dirname(filename))
    return os.path.join(base_dir, 'data', f)

#------------------------------------------------------------------------------
# Comparators


def equalContents(arr1, arr2):
    """Checks if the set of unique elements of arr1 and arr2 are equivalent.
    """
    return frozenset(arr1) == frozenset(arr2)


def assert_equal(a, b, msg=""):
    """asserts that a equals b, like nose's assert_equal, but allows custom message to start.
    Passes a and b to format string as well. So you can use '{0}' and '{1}' to display a and b.

    Examples
    --------
    >>> assert_equal(2, 2, "apples")
    >>> assert_equal(5.2, 1.2, "{0} was really a dead parrot")
    Traceback (most recent call last):
        ...
    AssertionError: 5.2 was really a dead parrot: 5.2 != 1.2
    """
    assert a == b, "%s: %r != %r" % (msg.format(a,b), a, b)


def assert_index_equal(left, right, exact=False, check_names=True,
                       check_less_precise=False, check_exact=True, obj='Index'):
    """Check that left and right Index are equal.

    Parameters
    ----------
    left : Index
    right : Index
    exact : bool, default False
        Whether to check the Index class, dtype and inferred_type are identical.
    check_names : bool, default True
        Whether to check the names attribute.
    check_less_precise : bool, default False
        Specify comparison precision. Only used when check_exact is False.
        5 digits (False) or 3 digits (True) after decimal points are compared.
    check_exact : bool, default True
        Whether to compare number exactly.
    obj : str, default 'Index'
        Specify object name being compared, internally used to show appropriate
        assertion message
    """

    def _check_types(l, r, obj='Index'):
        if exact:
            if type(l) != type(r):
                msg = '{0} classes are different'.format(obj)
                raise_assert_detail(obj, msg, l, r)
            assert_attr_equal('dtype', l, r, obj=obj)
            assert_attr_equal('inferred_type', l, r, obj=obj)

    def _get_ilevel_values(index, level):
        # accept level number only
        unique = index.levels[level]
        labels = index.labels[level]
        filled = take_1d(unique.values, labels, fill_value=unique._na_value)
        values = unique._simple_new(filled, index.names[level],
                                    freq=getattr(unique, 'freq', None),
                                    tz=getattr(unique, 'tz', None))
        return values

    # instance validation
    assertIsInstance(left, Index, '[index] ')
    assertIsInstance(right, Index, '[index] ')

    # class / dtype comparison
    _check_types(left, right)

    # level comparison
    if left.nlevels != right.nlevels:
        raise_assert_detail(obj, '{0} levels are different'.format(obj),
                            '{0}, {1}'.format(left.nlevels, left),
                            '{0}, {1}'.format(right.nlevels, right))

    # length comparison
    if len(left) != len(right):
        raise_assert_detail(obj, '{0} length are different'.format(obj),
                           '{0}, {1}'.format(len(left), left),
                           '{0}, {1}'.format(len(right), right))

    # MultiIndex special comparison for little-friendly error messages
    if left.nlevels > 1:
        for level in range(left.nlevels):
            # cannot use get_level_values here because it can change dtype
            llevel = _get_ilevel_values(left, level)
            rlevel = _get_ilevel_values(right, level)

            lobj = 'MultiIndex level [{0}]'.format(level)
            assert_index_equal(llevel, rlevel,
                               exact=exact, check_names=check_names,
                               check_less_precise=check_less_precise,
                               check_exact=check_exact, obj=lobj)
            # get_level_values may change dtype
            _check_types(left.levels[level], right.levels[level], obj=obj)

    if check_exact:
        if not left.equals(right):
            diff = np.sum((left.values != right.values).astype(int)) * 100.0 / len(left)
            msg = '{0} values are different ({1} %)'.format(obj, np.round(diff, 5))
            raise_assert_detail(obj, msg, left, right)
    else:
        assert_almost_equal(left.values, right.values,
                            check_less_precise=check_less_precise,
                            obj=obj, lobj=left, robj=right)

    # metadata comparison
    if check_names:
        assert_attr_equal('names', left, right, obj=obj)


def assert_attr_equal(attr, left, right, obj='Attributes'):
    """checks attributes are equal. Both objects must have attribute.

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

    left_attr = getattr(left, attr)
    right_attr = getattr(right, attr)

    if left_attr is right_attr:
        return True
    elif (is_number(left_attr) and np.isnan(left_attr) and
          is_number(right_attr) and np.isnan(right_attr)):
        # np.nan
        return True

    result = left_attr == right_attr
    if not isinstance(result, bool):
        result = result.all()

    if result:
        return True
    else:
        raise_assert_detail(obj, 'Attribute "{0}" are different'.format(attr),
                        left_attr, right_attr)


def isiterable(obj):
    return hasattr(obj, '__iter__')

def is_sorted(seq):
    return assert_almost_equal(seq, np.sort(np.array(seq)))


def assertIs(first, second, msg=''):
    """Checks that 'first' is 'second'"""
    a, b = first, second
    assert a is b, "%s: %r is not %r" % (msg.format(a, b), a, b)


def assertIsNot(first, second, msg=''):
    """Checks that 'first' is not 'second'"""
    a, b = first, second
    assert a is not b, "%s: %r is %r" % (msg.format(a, b), a, b)


def assertIn(first, second, msg=''):
    """Checks that 'first' is in 'second'"""
    a, b = first, second
    assert a in b, "%s: %r is not in %r" % (msg.format(a, b), a, b)


def assertNotIn(first, second, msg=''):
    """Checks that 'first' is not in 'second'"""
    a, b = first, second
    assert a not in b, "%s: %r is in %r" % (msg.format(a, b), a, b)


def assertIsNone(expr, msg=''):
    """Checks that 'expr' is None"""
    return assertIs(expr, None, msg)


def assertIsNotNone(expr, msg=''):
    """Checks that 'expr' is not None"""
    return assertIsNot(expr, None, msg)


def assertIsInstance(obj, cls, msg=''):
    """Test that obj is an instance of cls
    (which can be a class or a tuple of classes,
    as supported by isinstance())."""
    assert isinstance(obj, cls), (
        "%sExpected object to be of type %r, found %r instead" % (
            msg, cls, type(obj)))

def assert_isinstance(obj, class_type_or_tuple, msg=''):
    return deprecate('assert_isinstance', assertIsInstance)(obj, class_type_or_tuple, msg=msg)


def assertNotIsInstance(obj, cls, msg=''):
    """Test that obj is not an instance of cls
    (which can be a class or a tuple of classes,
    as supported by isinstance())."""
    assert not isinstance(obj, cls), (
        "%sExpected object to be of type %r, found %r instead" % (
            msg, cls, type(obj)))


def assert_categorical_equal(res, exp):

    if not array_equivalent(res.categories, exp.categories):
        raise AssertionError(
            'categories not equivalent: {0} vs {1}.'.format(res.categories,
                                                            exp.categories))
    if not array_equivalent(res.codes, exp.codes):
        raise AssertionError(
            'codes not equivalent: {0} vs {1}.'.format(res.codes, exp.codes))

    if res.ordered != exp.ordered:
        raise AssertionError("ordered not the same")


def raise_assert_detail(obj, message, left, right):
    if isinstance(left, np.ndarray):
        left = pprint_thing(left)
    if isinstance(right, np.ndarray):
        right = pprint_thing(right)

    msg = """{0} are different

{1}
[left]:  {2}
[right]: {3}""".format(obj, message, left, right)
    raise AssertionError(msg)


def assert_numpy_array_equal(left, right,
                             strict_nan=False, err_msg=None,
                             obj='numpy array'):
    """Checks that 'np_array' is equivalent to 'assert_equal'.

    This is similar to ``numpy.testing.assert_array_equal``, but can
    check equality including ``np.nan``. Two numpy arrays are regarded as
    equivalent if the arrays have equal non-NaN elements,
    and `np.nan` in corresponding locations.
    """

    # compare shape and values
    if array_equivalent(left, right, strict_nan=strict_nan):
        return

    if err_msg is None:
        # show detailed error

        if np.isscalar(left) and np.isscalar(right):
            # show scalar comparison error
            assert_equal(left, right)
        elif is_list_like(left) and is_list_like(right):
            # some test cases pass list
            left = np.asarray(left)
            right = np.array(right)

            if left.shape != right.shape:
                raise_assert_detail(obj, '{0} shapes are different'.format(obj),
                                    left.shape, right.shape)

            diff = 0
            for l, r in zip(left, right):
                # count up differences
                if not array_equivalent(l, r, strict_nan=strict_nan):
                    diff += 1

            diff = diff * 100.0 / left.size
            msg = '{0} values are different ({1} %)'.format(obj, np.round(diff, 5))
            raise_assert_detail(obj, msg, left, right)
        elif is_list_like(left):
            msg = "First object is iterable, second isn't"
            raise_assert_detail(obj, msg, left, right)
        else:
            msg = "Second object is iterable, first isn't"
            raise_assert_detail(obj, msg, left, right)

    raise AssertionError(err_msg)


# This could be refactored to use the NDFrame.equals method
def assert_series_equal(left, right, check_dtype=True,
                        check_index_type=False,
                        check_series_type=False,
                        check_less_precise=False,
                        check_names=True,
                        check_exact=False,
                        check_datetimelike_compat=False,
                        obj='Series'):

    """Check that left and right Series are equal.

    Parameters
    ----------
    left : Series
    right : Series
    check_dtype : bool, default True
        Whether to check the Series dtype is identical.
    check_index_type : bool, default False
        Whether to check the Index class, dtype and inferred_type are identical.
    check_series_type : bool, default False
        Whether to check the Series class is identical.
    check_less_precise : bool, default False
        Specify comparison precision. Only used when check_exact is False.
        5 digits (False) or 3 digits (True) after decimal points are compared.
    check_exact : bool, default False
        Whether to compare number exactly.
    check_names : bool, default True
        Whether to check the Series and Index names attribute.
    check_dateteimelike_compat : bool, default False
        Compare datetime-like which is comparable ignoring dtype.
    obj : str, default 'Series'
        Specify object name being compared, internally used to show appropriate
        assertion message
    """

    # instance validation
    assertIsInstance(left, Series, '[Series] ')
    assertIsInstance(right, Series, '[Series] ')

    if check_series_type:
        assertIsInstance(left, type(right))

    # length comparison
    if len(left) != len(right):
        raise_assert_detail(obj, 'Series length are different',
                            '{0}, {1}'.format(len(left), left.index),
                            '{0}, {1}'.format(len(right), right.index))

    # index comparison
    assert_index_equal(left.index, right.index, exact=check_index_type,
                       check_names=check_names,
                       check_less_precise=check_less_precise, check_exact=check_exact,
                       obj='{0}.index'.format(obj))

    if check_dtype:
        assert_attr_equal('dtype', left, right)

    if check_exact:
        assert_numpy_array_equal(left.get_values(), right.get_values(),
                                 obj='{0}'.format(obj))
    elif check_datetimelike_compat:
        # we want to check only if we have compat dtypes
        # e.g. integer and M|m are NOT compat, but we can simply check the values in that case
        if is_datetimelike_v_numeric(left, right) or is_datetimelike_v_object(left, right) or needs_i8_conversion(left) or needs_i8_conversion(right):

            # datetimelike may have different objects (e.g. datetime.datetime vs Timestamp) but will compare equal
            if not Index(left.values).equals(Index(right.values)):
                raise AssertionError(
                    '[datetimelike_compat=True] {0} is not equal to {1}.'.format(left.values,
                                                                                 right.values))
        else:
            assert_numpy_array_equal(left.values, right.values)
    else:
        assert_almost_equal(left.get_values(), right.get_values(),
                            check_less_precise, obj='{0}'.format(obj))

    # metadata comparison
    if check_names:
        assert_attr_equal('name', left, right, obj=obj)


# This could be refactored to use the NDFrame.equals method
def assert_frame_equal(left, right, check_dtype=True,
                       check_index_type=False,
                       check_column_type=False,
                       check_frame_type=False,
                       check_less_precise=False,
                       check_names=True,
                       by_blocks=False,
                       check_exact=False,
                       check_datetimelike_compat=False,
                       obj='DataFrame'):

    """Check that left and right DataFrame are equal.

    Parameters
    ----------
    left : DataFrame
    right : DataFrame
    check_dtype : bool, default True
        Whether to check the DataFrame dtype is identical.
    check_index_type : bool, default False
        Whether to check the Index class, dtype and inferred_type are identical.
    check_column_type : bool, default False
        Whether to check the columns class, dtype and inferred_type are identical.
    check_frame_type : bool, default False
        Whether to check the DataFrame class is identical.
    check_less_precise : bool, default False
        Specify comparison precision. Only used when check_exact is False.
        5 digits (False) or 3 digits (True) after decimal points are compared.
    check_names : bool, default True
        Whether to check the Index names attribute.
    by_blocks : bool, default False
        Specify how to compare internal data. If False, compare by columns.
        If True, compare by blocks.
    check_exact : bool, default False
        Whether to compare number exactly.
    check_dateteimelike_compat : bool, default False
        Compare datetime-like which is comparable ignoring dtype.
    obj : str, default 'DataFrame'
        Specify object name being compared, internally used to show appropriate
        assertion message
    """

    # instance validation
    assertIsInstance(left, DataFrame, '[DataFrame] ')
    assertIsInstance(right, DataFrame, '[DataFrame] ')

    if check_frame_type:
        assertIsInstance(left, type(right))

    # shape comparison (row)
    if left.shape[0] != right.shape[0]:
        raise_assert_detail(obj, 'DataFrame shape (number of rows) are different',
                            '{0}, {1}'.format(left.shape[0], left.index),
                            '{0}, {1}'.format(right.shape[0], right.index))
    # shape comparison (columns)
    if left.shape[1] != right.shape[1]:
        raise_assert_detail(obj, 'DataFrame shape (number of columns) are different',
                            '{0}, {1}'.format(left.shape[1], left.columns),
                            '{0}, {1}'.format(right.shape[1], right.columns))

    # index comparison
    assert_index_equal(left.index, right.index, exact=check_index_type,
                       check_names=check_names,
                       check_less_precise=check_less_precise, check_exact=check_exact,
                       obj='{0}.index'.format(obj))

    # column comparison
    assert_index_equal(left.columns, right.columns, exact=check_column_type,
                       check_names=check_names,
                       check_less_precise=check_less_precise, check_exact=check_exact,
                       obj='{0}.columns'.format(obj))

    # compare by blocks
    if by_blocks:
        rblocks = right.blocks
        lblocks = left.blocks
        for dtype in list(set(list(lblocks.keys()) + list(rblocks.keys()))):
            assert dtype in lblocks
            assert dtype in rblocks
            assert_frame_equal(lblocks[dtype], rblocks[dtype],
                               check_dtype=check_dtype, obj='DataFrame.blocks')

    # compare by columns
    else:
        for i, col in enumerate(left.columns):
            assert col in right
            lcol = left.iloc[:, i]
            rcol = right.iloc[:, i]
            assert_series_equal(lcol, rcol,
                                check_dtype=check_dtype,
                                check_index_type=check_index_type,
                                check_less_precise=check_less_precise,
                                check_exact=check_exact,
                                check_names=check_names,
                                check_datetimelike_compat=check_datetimelike_compat,
                                obj='DataFrame.iloc[:, {0}]'.format(i))


def assert_panelnd_equal(left, right,
                         check_panel_type=False,
                         check_less_precise=False,
                         assert_func=assert_frame_equal,
                         check_names=False):
    if check_panel_type:
        assertIsInstance(left, type(right))

    for axis in ['items', 'major_axis', 'minor_axis']:
        left_ind = getattr(left, axis)
        right_ind = getattr(right, axis)
        assert_index_equal(left_ind, right_ind, check_names=check_names)

    for i, item in enumerate(left._get_axis(0)):
        assert item in right, "non-matching item (right) '%s'" % item
        litem = left.iloc[i]
        ritem = right.iloc[i]
        assert_func(litem, ritem, check_less_precise=check_less_precise)

    for i, item in enumerate(right._get_axis(0)):
        assert item in left, "non-matching item (left) '%s'" % item

# TODO: strangely check_names fails in py3 ?
_panel_frame_equal = partial(assert_frame_equal, check_names=False)
assert_panel_equal = partial(assert_panelnd_equal,
                             assert_func=_panel_frame_equal)
assert_panel4d_equal = partial(assert_panelnd_equal,
                               assert_func=assert_panel_equal)


def assert_contains_all(iterable, dic):
    for k in iterable:
        assert k in dic, "Did not contain item: '%r'" % k


def assert_copy(iter1, iter2, **eql_kwargs):
    """
    iter1, iter2: iterables that produce elements comparable with assert_almost_equal

    Checks that the elements are equal, but not the same object. (Does not
    check that items in sequences are also not the same object)
    """
    for elem1, elem2 in zip(iter1, iter2):
        assert_almost_equal(elem1, elem2, **eql_kwargs)
        assert elem1 is not elem2, "Expected object %r and object %r to be different objects, were same." % (
                                    type(elem1), type(elem2))


def getCols(k):
    return string.ascii_uppercase[:k]

def getArangeMat():
    return np.arange(N * K).reshape((N, K))


# make index
def makeStringIndex(k=10, name=None):
    return Index(rands_array(nchars=10, size=k), name=name)


def makeUnicodeIndex(k=10, name=None):
    return Index(randu_array(nchars=10, size=k))

def makeCategoricalIndex(k=10, n=3, name=None):
    """ make a length k index or n categories """
    x = rands_array(nchars=4, size=n)
    return CategoricalIndex(np.random.choice(x,k), name=name)

def makeBoolIndex(k=10, name=None):
    if k == 1:
        return Index([True], name=name)
    elif k == 2:
        return Index([False,True], name=name)
    return Index([False,True] + [False]*(k-2), name=name)

def makeIntIndex(k=10, name=None):
    return Index(lrange(k), name=name)

def makeFloatIndex(k=10, name=None):
    values = sorted(np.random.random_sample(k)) - np.random.random_sample(1)
    return Index(values * (10 ** np.random.randint(0, 9)), name=name)

def makeDateIndex(k=10, freq='B', name=None):
    dt = datetime(2000, 1, 1)
    dr = bdate_range(dt, periods=k, freq=freq, name=name)
    return DatetimeIndex(dr, name=name)

def makeTimedeltaIndex(k=10, freq='D', name=None):
    return TimedeltaIndex(start='1 day', periods=k, freq=freq, name=name)

def makePeriodIndex(k=10, name=None):
    dt = datetime(2000, 1, 1)
    dr = PeriodIndex(start=dt, periods=k, freq='B', name=name)
    return dr

def all_index_generator(k=10):
    """Generator which can be iterated over to get instances of all the various
    index classes.

    Parameters
    ----------
    k: length of each of the index instances
    """
    all_make_index_funcs = [makeIntIndex, makeFloatIndex, makeStringIndex,
                            makeUnicodeIndex, makeDateIndex, makePeriodIndex,
                            makeTimedeltaIndex, makeBoolIndex,
                            makeCategoricalIndex]
    for make_index_func in all_make_index_funcs:
        yield make_index_func(k=k)

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


# make series
def makeFloatSeries(name=None):
    index = makeStringIndex(N)
    return Series(randn(N), index=index, name=name)


def makeStringSeries(name=None):
    index = makeStringIndex(N)
    return Series(randn(N), index=index, name=name)


def makeObjectSeries(name=None):
    dateIndex = makeDateIndex(N)
    dateIndex = Index(dateIndex, dtype=object)
    index = makeStringIndex(N)
    return Series(dateIndex, index=index, name=name)


def getSeriesData():
    index = makeStringIndex(N)
    return dict((c, Series(randn(N), index=index)) for c in getCols(K))


def makeTimeSeries(nper=None, freq='B', name=None):
    if nper is None:
        nper = N
    return Series(randn(nper), index=makeDateIndex(nper, freq=freq), name=name)


def makePeriodSeries(nper=None, name=None):
    if nper is None:
        nper = N
    return Series(randn(nper), index=makePeriodIndex(nper), name=name)


def getTimeSeriesData(nper=None, freq='B'):
    return dict((c, makeTimeSeries(nper, freq)) for c in getCols(K))


def getPeriodData(nper=None):
    return dict((c, makePeriodSeries(nper)) for c in getCols(K))

# make frame
def makeTimeDataFrame(nper=None, freq='B'):
    data = getTimeSeriesData(nper, freq)
    return DataFrame(data)


def makeDataFrame():
    data = getSeriesData()
    return DataFrame(data)


def getMixedTypeDict():
    index = Index(['a', 'b', 'c', 'd', 'e'])

    data = {
        'A': [0., 1., 2., 3., 4.],
        'B': [0., 1., 0., 1., 0.],
        'C': ['foo1', 'foo2', 'foo3', 'foo4', 'foo5'],
        'D': bdate_range('1/1/2009', periods=5)
    }

    return index, data

def makeMixedDataFrame():
    return DataFrame(getMixedTypeDict()[1])

def makePeriodFrame(nper=None):
    data = getPeriodData(nper)
    return DataFrame(data)


def makePanel(nper=None):
    cols = ['Item' + c for c in string.ascii_uppercase[:K - 1]]
    data = dict((c, makeTimeDataFrame(nper)) for c in cols)
    return Panel.fromDict(data)


def makePeriodPanel(nper=None):
    cols = ['Item' + c for c in string.ascii_uppercase[:K - 1]]
    data = dict((c, makePeriodFrame(nper)) for c in cols)
    return Panel.fromDict(data)


def makePanel4D(nper=None):
    return Panel4D(dict(l1=makePanel(nper), l2=makePanel(nper),
                        l3=makePanel(nper)))


def makeCustomIndex(nentries, nlevels, prefix='#', names=False, ndupe_l=None,
                    idx_type=None):
    """Create an index/multindex with given dimensions, levels, names, etc'

    nentries - number of entries in index
    nlevels - number of levels (> 1 produces multindex)
    prefix - a string prefix for labels
    names - (Optional), bool or list of strings. if True will use default names,
       if false will use no names, if a list is given,  the name of each level
       in the index will be taken from the list.
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
    assert (is_sequence(ndupe_l) and len(ndupe_l) <= nlevels)
    assert (names is None or names is False
            or names is True or len(names) is nlevels)
    assert idx_type is None or \
        (idx_type in ('i', 'f', 's', 'u', 'dt', 'p', 'td') and nlevels == 1)

    if names is True:
        # build default names
        names = [prefix + str(i) for i in range(nlevels)]
    if names is False:
        # pass None to index constructor for no name
        names = None

    # make singelton case uniform
    if isinstance(names, compat.string_types) and nlevels == 1:
        names = [names]

    # specific 1D index type requested?
    idx_func = dict(i=makeIntIndex, f=makeFloatIndex, s=makeStringIndex,
                    u=makeUnicodeIndex, dt=makeDateIndex, td=makeTimedeltaIndex,
                    p=makePeriodIndex).get(idx_type)
    if idx_func:
        idx = idx_func(nentries)
        # but we need to fill in the name
        if names:
            idx.name = names[0]
        return idx
    elif idx_type is not None:
        raise ValueError('"%s" is not a legal value for `idx_type`, use  '
                         '"i"/"f"/"s"/"u"/"dt/"p"/"td".' % idx_type)

    if len(ndupe_l) < nlevels:
        ndupe_l.extend([1] * (nlevels - len(ndupe_l)))
    assert len(ndupe_l) == nlevels

    assert all([x > 0 for x in ndupe_l])

    tuples = []
    for i in range(nlevels):
        def keyfunc(x):
            import re
            numeric_tuple = re.sub("[^\d_]_?", "", x).split("_")
            return lmap(int, numeric_tuple)

        # build a list of lists to create the index from
        div_factor = nentries // ndupe_l[i] + 1
        cnt = Counter()
        for j in range(div_factor):
            label = prefix + '_l%d_g' % i + str(j)
            cnt[label] = ndupe_l[i]
        # cute Counter trick
        result = list(sorted(cnt.elements(), key=keyfunc))[:nentries]
        tuples.append(result)

    tuples = lzip(*tuples)

    # convert tuples to index
    if nentries == 1:
        index = Index(tuples[0], name=names[0])
    else:
        index = MultiIndex.from_tuples(tuples, names=names)
    return index


def makeCustomDataframe(nrows, ncols, c_idx_names=True, r_idx_names=True,
                        c_idx_nlevels=1, r_idx_nlevels=1, data_gen_f=None,
                        c_ndupe_l=None, r_ndupe_l=None, dtype=None,
                        c_idx_type=None, r_idx_type=None):
    """
   nrows,  ncols - number of data rows/cols
   c_idx_names, idx_names  - False/True/list of strings,  yields No names ,
        default names or  uses the provided names for the levels of the
        corresponding  index. You can provide a single string when
        c_idx_nlevels ==1.
   c_idx_nlevels - number of levels in columns index. > 1 will yield MultiIndex
   r_idx_nlevels - number of levels in rows index. > 1 will yield MultiIndex
   data_gen_f - a function f(row,col) which return the data value at that position,
        the default generator used yields values of the form "RxCy" based on position.
   c_ndupe_l, r_ndupe_l - list of integers, determines the number
        of duplicates for each label at a given level of the corresponding index.
        The default `None` value produces a multiplicity of 1 across
        all levels, i.e. a unique index. Will accept a partial list of
        length N < idx_nlevels, for just the first N levels. If ndupe
        doesn't divide nrows/ncol, the last label might have lower multiplicity.
   dtype - passed to the DataFrame constructor as is, in case you wish to
        have more control in conjuncion with a custom `data_gen_f`
   r_idx_type, c_idx_type -  "i"/"f"/"s"/"u"/"dt"/"td".
       If idx_type is not None, `idx_nlevels` must be 1.
       "i"/"f" creates an integer/float index,
       "s"/"u" creates a string/unicode index
       "dt" create a datetime index.
       "td" create a timedelta index.

        if unspecified, string labels will be generated.

    Examples:

    # 5 row, 3 columns, default names on both, single index on both axis
    >> makeCustomDataframe(5,3)

    # make the data a random int between 1 and 100
    >> mkdf(5,3,data_gen_f=lambda r,c:randint(1,100))

    # 2-level multiindex on rows with each label duplicated twice on first level,
    # default names on both axis, single index on both axis
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
    assert r_idx_type is None or \
        (r_idx_type in ('i', 'f', 's', 'u', 'dt', 'p', 'td') and r_idx_nlevels == 1)
    assert c_idx_type is None or \
        (c_idx_type in ('i', 'f', 's', 'u', 'dt', 'p', 'td') and c_idx_nlevels == 1)

    columns = makeCustomIndex(ncols, nlevels=c_idx_nlevels, prefix='C',
                              names=c_idx_names, ndupe_l=c_ndupe_l,
                              idx_type=c_idx_type)
    index = makeCustomIndex(nrows, nlevels=r_idx_nlevels, prefix='R',
                            names=r_idx_names, ndupe_l=r_ndupe_l,
                            idx_type=r_idx_type)

    # by default, generate data based on location
    if data_gen_f is None:
        data_gen_f = lambda r, c: "R%dC%d" % (r, c)

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

    j = np.floor(ind * 1. / nrows).astype(int)
    i = (ind - j * nrows).astype(int)
    return i.tolist(), j.tolist()


def makeMissingCustomDataframe(nrows, ncols, density=.9, random_state=None,
                               c_idx_names=True, r_idx_names=True,
                               c_idx_nlevels=1, r_idx_nlevels=1,
                               data_gen_f=None,
                               c_ndupe_l=None, r_ndupe_l=None, dtype=None,
                               c_idx_type=None, r_idx_type=None):
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
    df = makeCustomDataframe(nrows, ncols, c_idx_names=c_idx_names,
                             r_idx_names=r_idx_names,
                             c_idx_nlevels=c_idx_nlevels,
                             r_idx_nlevels=r_idx_nlevels,
                             data_gen_f=data_gen_f,
                             c_ndupe_l=c_ndupe_l, r_ndupe_l=r_ndupe_l,
                             dtype=dtype, c_idx_type=c_idx_type,
                             r_idx_type=r_idx_type)

    i, j = _create_missing_idx(nrows, ncols, density, random_state)
    df.values[i, j] = np.nan
    return df


def makeMissingDataframe(density=.9, random_state=None):
    df = makeDataFrame()
    i, j = _create_missing_idx(*df.shape, density=density,
                               random_state=random_state)
    df.values[i, j] = np.nan
    return df


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


# Dependency checks.  Copied this from Nipy/Nipype (Copyright of
# respective developers, license: BSD-3)
def package_check(pkg_name, version=None, app='pandas', checker=LooseVersion,
                  exc_failed_import=ImportError,
                  exc_failed_check=RuntimeError):
    """Check that the minimal version of the required package is installed.

    Parameters
    ----------
    pkg_name : string
        Name of the required package.
    version : string, optional
        Minimal version number for required package.
    app : string, optional
        Application that is performing the check.  For instance, the
        name of the tutorial being executed that depends on specific
        packages.
    checker : object, optional
        The class that will perform the version checking.  Default is
        distutils.version.LooseVersion.
    exc_failed_import : Exception, optional
        Class of the exception to be thrown if import failed.
    exc_failed_check : Exception, optional
        Class of the exception to be thrown if version check failed.

    Examples
    --------
    package_check('numpy', '1.3')
    package_check('networkx', '1.0', 'tutorial1')

    """

    if app:
        msg = '%s requires %s' % (app, pkg_name)
    else:
        msg = 'module requires %s' % pkg_name
    if version:
        msg += ' with version >= %s' % (version,)
    try:
        mod = __import__(pkg_name)
    except ImportError:
        raise exc_failed_import(msg)
    if not version:
        return
    try:
        have_version = mod.__version__
    except AttributeError:
        raise exc_failed_check('Cannot find version for %s' % pkg_name)
    if checker(have_version) < checker(version):
        raise exc_failed_check(msg)


def skip_if_no_package(*args, **kwargs):
    """Raise SkipTest if package_check fails

    Parameters
    ----------
    *args Positional parameters passed to `package_check`
    *kwargs Keyword parameters passed to `package_check`
    """
    from nose import SkipTest
    package_check(exc_failed_import=SkipTest,
                  exc_failed_check=SkipTest,
                  *args, **kwargs)

#
# Additional tags decorators for nose
#


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
    'Temporary failure in name resolution',
    'Name or service not known',
    'Connection refused',
)

# or this e.errno/e.reason.errno
_network_errno_vals = (
    101, # Network is unreachable
    111, # Connection refused
    110, # Connection timed out
    104, # Connection reset Error
    54,  # Connection reset by peer
    60,  # urllib.error.URLError: [Errno 60] Connection timed out
    )

# Both of the above shouldn't mask real issues such as 404's
# or refused connections (changed DNS).
# But some tests (test_data yahoo) contact incredibly flakey
# servers.

# and conditionally raise on these exception types
_network_error_classes = (IOError, httplib.HTTPException)

if sys.version_info >= (3, 3):
    _network_error_classes += (TimeoutError,)

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
        The url to test via ``pandas.io.common.urlopen`` to check for connectivity.
        Defaults to 'http://www.google.com'.
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
        message. Intended to supress errors where an errno isn't available.

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
        SkipTest

    Errors not related to networking will always be raised.
    """
    from nose import SkipTest
    t.network = True

    @wraps(t)
    def wrapper(*args, **kwargs):
        if check_before_test and not raise_on_error:
            if not can_connect(url, error_classes):
                raise SkipTest
        try:
            return t(*args, **kwargs)
        except Exception as e:
            errno = getattr(e, 'errno', None)
            if not errno and hasattr(errno, "reason"):
                errno = getattr(e.reason, 'errno', None)

            if errno in skip_errnos:
                raise SkipTest("Skipping test due to known errno"
                               " and error %s" % e)

            try:
                e_str = traceback.format_exc(e)
            except:
                e_str = str(e)

            if any([m.lower() in e_str.lower() for m in _skip_on_messages]):
                raise SkipTest("Skipping test because exception message is known"
                               " and error %s" % e)

            if not isinstance(e, error_classes):
                raise

            if raise_on_error or can_connect(url, error_classes):
                raise
            else:
                raise SkipTest("Skipping test due to lack of connectivity"
                               " and error %s" % e)

    return wrapper


with_connectivity_check = network


class SimpleMock(object):

    """
    Poor man's mocking object

    Note: only works for new-style classes, assumes  __getattribute__ exists.

    >>> a = type("Duck",(),{})
    >>> a.attr1,a.attr2 ="fizz","buzz"
    >>> b = SimpleMock(a,"attr1","bar")
    >>> b.attr1 == "bar" and b.attr2 == "buzz"
    True
    >>> a.attr1 == "fizz" and a.attr2 == "buzz"
    True
    """

    def __init__(self, obj, *args, **kwds):
        assert(len(args) % 2 == 0)
        attrs = kwds.get("attrs", {})
        for k, v in zip(args[::2], args[1::2]):
            # dict comprehensions break 2.6
            attrs[k] = v
        self.attrs = attrs
        self.obj = obj

    def __getattribute__(self, name):
        attrs = object.__getattribute__(self, "attrs")
        obj = object.__getattribute__(self, "obj")
        return attrs.get(name, type(obj).__getattribute__(obj, name))


@contextmanager
def stdin_encoding(encoding=None):
    """
    Context manager for running bits of code while emulating an arbitrary
    stdin encoding.

    >>> import sys
    >>> _encoding = sys.stdin.encoding
    >>> with stdin_encoding('AES'): sys.stdin.encoding
    'AES'
    >>> sys.stdin.encoding==_encoding
    True

    """
    import sys

    _stdin = sys.stdin
    sys.stdin = SimpleMock(sys.stdin, "encoding", encoding)
    yield
    sys.stdin = _stdin


def assertRaises(_exception, _callable=None, *args, **kwargs):
    """assertRaises that is usable as context manager or in a with statement

    Exceptions that don't match the given Exception type fall through::

    >>> with assertRaises(ValueError):
    ...     raise TypeError("banana")
    ...
    Traceback (most recent call last):
        ...
    TypeError: banana

    If it raises the given Exception type, the test passes
    >>> with assertRaises(KeyError):
    ...     dct = dict()
    ...     dct["apple"]

    If the expected error doesn't occur, it raises an error.
    >>> with assertRaises(KeyError):
    ...     dct = {'apple':True}
    ...     dct["apple"]
    Traceback (most recent call last):
        ...
    AssertionError: KeyError not raised.

    In addition to using it as a contextmanager, you can also use it as a
    function, just like the normal assertRaises

    >>> assertRaises(TypeError, ",".join, [1, 3, 5]);
    """
    manager = _AssertRaisesContextmanager(exception=_exception)
    # don't return anything if used in function form
    if _callable is not None:
        with manager:
            _callable(*args, **kwargs)
    else:
        return manager

def assertRaisesRegexp(_exception, _regexp, _callable=None, *args, **kwargs):
    """ Port of assertRaisesRegexp from unittest in Python 2.7 - used in with statement.

    Explanation from standard library:
        Like assertRaises() but also tests that regexp matches on the string
        representation of the raised exception. regexp may be a regular expression
        object or a string containing a regular expression suitable for use by
        re.search().

    You can pass either a regular expression or a compiled regular expression object.
    >>> assertRaisesRegexp(ValueError, 'invalid literal for.*XYZ',
    ...                                int, 'XYZ');
    >>> import re
    >>> assertRaisesRegexp(ValueError, re.compile('literal'), int, 'XYZ');

    If an exception of a different type is raised, it bubbles up.

    >>> assertRaisesRegexp(TypeError, 'literal', int, 'XYZ');
    Traceback (most recent call last):
        ...
    ValueError: invalid literal for int() with base 10: 'XYZ'
    >>> dct = dict()
    >>> assertRaisesRegexp(KeyError, 'pear', dct.__getitem__, 'apple');
    Traceback (most recent call last):
        ...
    AssertionError: "pear" does not match "'apple'"

    You can also use this in a with statement.
    >>> with assertRaisesRegexp(TypeError, 'unsupported operand type\(s\)'):
    ...     1 + {}
    >>> with assertRaisesRegexp(TypeError, 'banana'):
    ...     'apple'[0] = 'b'
    Traceback (most recent call last):
        ...
    AssertionError: "banana" does not match "'str' object does not support \
item assignment"
    """
    manager = _AssertRaisesContextmanager(exception=_exception, regexp=_regexp)
    if _callable is not None:
        with manager:
            _callable(*args, **kwargs)
    else:
        return manager


class _AssertRaisesContextmanager(object):
    """handles the behind the scenes work for assertRaises and assertRaisesRegexp"""
    def __init__(self, exception, regexp=None, *args, **kwargs):
        self.exception = exception
        if regexp is not None and not hasattr(regexp, "search"):
            regexp = re.compile(regexp, re.DOTALL)
        self.regexp = regexp

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        expected = self.exception
        if not exc_type:
            name = getattr(expected, "__name__", str(expected))
            raise AssertionError("{0} not raised.".format(name))
        if issubclass(exc_type, expected):
            return self.handle_success(exc_type, exc_value, traceback)
        return self.handle_failure(exc_type, exc_value, traceback)

    def handle_failure(*args, **kwargs):
        # Failed, so allow Exception to bubble up
        return False

    def handle_success(self, exc_type, exc_value, traceback):
        if self.regexp is not None:
            val = str(exc_value)
            if not self.regexp.search(val):
                e = AssertionError('"%s" does not match "%s"' %
                                   (self.regexp.pattern, str(val)))
                raise_with_traceback(e, traceback)
        return True


@contextmanager
def assert_produces_warning(expected_warning=Warning, filter_level="always",
                            clear=None, check_stacklevel=True):
    """
    Context manager for running code that expects to raise (or not raise)
    warnings.  Checks that code raises the expected warning and only the
    expected warning. Pass ``False`` or ``None`` to check that it does *not*
    raise a warning. Defaults to ``exception.Warning``, baseclass of all
    Warnings. (basically a wrapper around ``warnings.catch_warnings``).

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
    with warnings.catch_warnings(record=True) as w:

        if clear is not None:
            # make sure that we are clearning these warnings
            # if they have happened before
            # to guarantee that we will catch them
            if not is_list_like(clear):
                clear = [ clear ]
            for m in clear:
                try:
                    m.__warningregistry__.clear()
                except:
                    pass

        saw_warning = False
        warnings.simplefilter(filter_level)
        yield w
        extra_warnings = []
        for actual_warning in w:
            if (expected_warning and issubclass(actual_warning.category,
                                                expected_warning)):
                saw_warning = True

                if check_stacklevel and issubclass(actual_warning.category,
                                                   (FutureWarning, DeprecationWarning)):
                    from inspect import getframeinfo, stack
                    caller = getframeinfo(stack()[2][0])
                    msg = ("Warning not set with correct stacklevel. File were warning"
                           " is raised: {0} != {1}. Warning message: {2}".format(
                               actual_warning.filename, caller.filename,
                               actual_warning.message))
                    assert actual_warning.filename == caller.filename, msg
            else:
                extra_warnings.append(actual_warning.category.__name__)
        if expected_warning:
            assert saw_warning, ("Did not see expected warning of class %r."
                                 % expected_warning.__name__)
        assert not extra_warnings, ("Caused unexpected warning(s): %r."
                                    % extra_warnings)


def disabled(t):
    t.disabled = True
    return t


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
def use_numexpr(use, min_elements=expr._MIN_ELEMENTS):
    olduse = expr._USE_NUMEXPR
    oldmin = expr._MIN_ELEMENTS
    expr.set_use_numexpr(use)
    expr._MIN_ELEMENTS = min_elements
    yield
    expr._MIN_ELEMENTS = oldmin
    expr.set_use_numexpr(olduse)


# Also provide all assert_* functions in the TestCase class
for name, obj in inspect.getmembers(sys.modules[__name__]):
    if inspect.isfunction(obj) and name.startswith('assert'):
        setattr(TestCase, name, staticmethod(obj))


def test_parallel(num_threads=2, kwargs_list=None):
    """Decorator to run the same function multiple times in parallel.

    Parameters
    ----------
    num_threads : int, optional
        The number of times the function is run in parallel.
    kwargs_list : list of dicts, optional
        The list of kwargs to update original function kwargs on different threads.
    Notes
    -----
    This decorator does not pass the return value of the decorated function.

    Original from scikit-image: https://github.com/scikit-image/scikit-image/pull/1519

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


class SubclassedDataFrame(DataFrame):
    _metadata = ['testattr']

    @property
    def _constructor(self):
        return SubclassedDataFrame
