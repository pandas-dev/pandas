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
from pandas.core.common import _is_sequence
import pandas.core.index as index
import pandas.core.series as series
import pandas.core.frame as frame
import pandas.core.panel as panel
import pandas.core.panel4d as panel4d
import pandas.compat as compat
from pandas.compat import(
    filter, map, zip, range, unichr, lrange, lmap, lzip, u, callable, Counter,
    raise_with_traceback, httplib
)

from pandas import bdate_range
from pandas.tseries.index import DatetimeIndex
from pandas.tseries.period import PeriodIndex

from pandas import _testing

from pandas.io.common import urlopen

Index = index.Index
MultiIndex = index.MultiIndex
Series = series.Series
DataFrame = frame.DataFrame
Panel = panel.Panel
Panel4D = panel4d.Panel4D

N = 30
K = 4
_RAISE_NETWORK_ERROR_DEFAULT = False

class TestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pd.set_option('chained_assignment','raise')
        #print("setting up: {0}".format(cls))

    @classmethod
    def tearDownClass(cls):
        #print("tearing down up: {0}".format(cls))
        pass

# NOTE: don't pass an NDFrame or index to this function - may not handle it
# well.
assert_almost_equal = _testing.assert_almost_equal

assert_dict_equal = _testing.assert_dict_equal

def randbool(size=(), p=0.5):
    return rand(*size) <= p


def rands(n):
    choices = string.ascii_letters + string.digits
    return ''.join(random.choice(choices) for _ in range(n))


def randu(n):
    choices = u("").join(map(unichr, lrange(1488, 1488 + 26)))
    choices += string.digits
    return ''.join([random.choice(choices) for _ in range(n)])


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
        locales = raw_locales.split(b'\n')
        raw_locales = []
        for x in raw_locales:
            try:
                raw_locales.append(str(x, encoding=pd.options.display.encoding))
            except:
                pass
    except TypeError:
        pass

    if prefix is None:
        return _valid_locales(raw_locales, normalize)

    found = re.compile('%s.*' % prefix).findall('\n'.join(raw_locales))
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


def assert_isinstance(obj, class_type_or_tuple, msg=''):
    """asserts that obj is an instance of class_type_or_tuple"""
    assert isinstance(obj, class_type_or_tuple), (
        "%sExpected object to be of type %r, found %r instead" % (
            msg, class_type_or_tuple, type(obj)))


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


def assert_index_equal(left, right):
    assert_isinstance(left, Index, '[index] ')
    assert_isinstance(right, Index, '[index] ')
    if not left.equals(right):
        raise AssertionError("[index] left [{0} {1}], right [{2} {3}]".format(left.dtype,
                                                                              left,
                                                                              right,
                                                                              right.dtype))


def assert_attr_equal(attr, left, right):
    """checks attributes are equal. Both objects must have attribute."""
    left_attr = getattr(left, attr)
    right_attr = getattr(right, attr)
    assert_equal(left_attr,right_attr,"attr is not equal [{0}]" .format(attr))


def isiterable(obj):
    return hasattr(obj, '__iter__')

def is_sorted(seq):
    return assert_almost_equal(seq, np.sort(np.array(seq)))

# This could be refactored to use the NDFrame.equals method
def assert_series_equal(left, right, check_dtype=True,
                        check_index_type=False,
                        check_series_type=False,
                        check_less_precise=False):
    if check_series_type:
        assert_isinstance(left, type(right))
    if check_dtype:
        assert_attr_equal('dtype', left, right)
    assert_almost_equal(left.values, right.values, check_less_precise)
    if check_less_precise:
        assert_almost_equal(
            left.index.values, right.index.values, check_less_precise)
    else:
        assert_index_equal(left.index, right.index)
    if check_index_type:
        assert_isinstance(left.index, type(right.index))
        assert_attr_equal('dtype', left.index, right.index)
        assert_attr_equal('inferred_type', left.index, right.index)

# This could be refactored to use the NDFrame.equals method
def assert_frame_equal(left, right, check_dtype=True,
                       check_index_type=False,
                       check_column_type=False,
                       check_frame_type=False,
                       check_less_precise=False,
                       check_names=True,
                       by_blocks=False):
    if check_frame_type:
        assert_isinstance(left, type(right))
    assert_isinstance(left, DataFrame)
    assert_isinstance(right, DataFrame)

    if check_less_precise:
        if not by_blocks:
            assert_almost_equal(left.columns, right.columns)
        assert_almost_equal(left.index, right.index)
    else:
        if not by_blocks:
            assert_index_equal(left.columns, right.columns)
        assert_index_equal(left.index, right.index)

    # compare by blocks
    if by_blocks:
        rblocks = right.blocks
        lblocks = left.blocks
        for dtype in list(set(list(lblocks.keys()) + list(rblocks.keys()))):
            assert dtype in lblocks
            assert dtype in rblocks
            assert_frame_equal(lblocks[dtype],rblocks[dtype],check_dtype=check_dtype)

    # compare by columns
    else:
        for i, col in enumerate(left.columns):
            assert col in right
            lcol = left.icol(i)
            rcol = right.icol(i)
            assert_series_equal(lcol, rcol,
                                check_dtype=check_dtype,
                                check_index_type=check_index_type,
                                check_less_precise=check_less_precise)

    if check_index_type:
        assert_isinstance(left.index, type(right.index))
        assert_attr_equal('dtype', left.index, right.index)
        assert_attr_equal('inferred_type', left.index, right.index)
    if check_column_type:
        assert_isinstance(left.columns, type(right.columns))
        assert_attr_equal('dtype', left.columns, right.columns)
        assert_attr_equal('inferred_type', left.columns, right.columns)
    if check_names:
        assert_attr_equal('names', left.index, right.index)
        assert_attr_equal('names', left.columns, right.columns)


def assert_panelnd_equal(left, right,
                         check_panel_type=False,
                         check_less_precise=False,
                         assert_func=assert_frame_equal):
    if check_panel_type:
        assert_isinstance(left, type(right))

    for axis in ['items', 'major_axis', 'minor_axis']:
        left_ind = getattr(left, axis)
        right_ind = getattr(right, axis)
        assert_index_equal(left_ind, right_ind)

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
def makeStringIndex(k=10):
    return Index([rands(10) for _ in range(k)])


def makeUnicodeIndex(k=10):
    return Index([randu(10) for _ in range(k)])


def makeIntIndex(k=10):
    return Index(lrange(k))


def makeFloatIndex(k=10):
    values = sorted(np.random.random_sample(k)) - np.random.random_sample(1)
    return Index(values * (10 ** np.random.randint(0, 9)))


def makeDateIndex(k=10):
    dt = datetime(2000, 1, 1)
    dr = bdate_range(dt, periods=k)
    return DatetimeIndex(dr)


def makePeriodIndex(k=10):
    dt = datetime(2000, 1, 1)
    dr = PeriodIndex(start=dt, periods=k, freq='B')
    return dr


# make series
def makeFloatSeries():
    index = makeStringIndex(N)
    return Series(randn(N), index=index)


def makeStringSeries():
    index = makeStringIndex(N)
    return Series(randn(N), index=index)


def makeObjectSeries():
    dateIndex = makeDateIndex(N)
    dateIndex = Index(dateIndex, dtype=object)
    index = makeStringIndex(N)
    return Series(dateIndex, index=index)


def getSeriesData():
    index = makeStringIndex(N)
    return dict((c, Series(randn(N), index=index)) for c in getCols(K))


def makeTimeSeries(nper=None):
    if nper is None:
        nper = N
    return Series(randn(nper), index=makeDateIndex(nper))


def makePeriodSeries(nper=None):
    if nper is None:
        nper = N
    return Series(randn(nper), index=makePeriodIndex(nper))


def getTimeSeriesData(nper=None):
    return dict((c, makeTimeSeries(nper)) for c in getCols(K))


def getPeriodData(nper=None):
    return dict((c, makePeriodSeries(nper)) for c in getCols(K))

# make frame
def makeTimeDataFrame(nper=None):
    data = getTimeSeriesData(nper)
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
    idx_type - "i"/"f"/"s"/"u"/"dt/"p".
       If idx_type is not None, `idx_nlevels` must be 1.
       "i"/"f" creates an integer/float index,
       "s"/"u" creates a string/unicode index
       "dt" create a datetime index.

        if unspecified, string labels will be generated.
    """

    if ndupe_l is None:
        ndupe_l = [1] * nlevels
    assert (_is_sequence(ndupe_l) and len(ndupe_l) <= nlevels)
    assert (names is None or names is False
            or names is True or len(names) is nlevels)
    assert idx_type is None or \
        (idx_type in ('i', 'f', 's', 'u', 'dt', 'p') and nlevels == 1)

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
                    u=makeUnicodeIndex, dt=makeDateIndex, p=makePeriodIndex).get(idx_type)
    if idx_func:
        idx = idx_func(nentries)
        # but we need to fill in the name
        if names:
            idx.name = names[0]
        return idx
    elif idx_type is not None:
        raise ValueError('"%s" is not a legal value for `idx_type`, use  '
                         '"i"/"f"/"s"/"u"/"dt/"p".' % idx_type)

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
   r_idx_type, c_idx_type -  "i"/"f"/"s"/"u"/"dt".
       If idx_type is not None, `idx_nlevels` must be 1.
       "i"/"f" creates an integer/float index,
       "s"/"u" creates a string/unicode index
       "dt" create a datetime index.

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
        (r_idx_type in ('i', 'f', 's', 'u', 'dt', 'p') and r_idx_nlevels == 1)
    assert c_idx_type is None or \
        (c_idx_type in ('i', 'f', 's', 'u', 'dt', 'p') and c_idx_nlevels == 1)

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
    )

# or this e.errno/e.reason.errno
_network_errno_vals = (
    101, # Network is unreachable
    110, # Connection timed out
    104, # Connection reset Error
    54,  # Connection reset by peer
    60,  # urllib.error.URLError: [Errno 60] Connection timed out
    )

# Both of the above shouldn't mask reasl issues such as 404's
# or refused connections (changed DNS).
# But some tests (test_data yahoo) contact incredibly flakey
# servers.

# and conditionally raise on these exception types
_network_error_classes = (IOError, httplib.HTTPException)

if sys.version_info[:2] >= (3,3):
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
            regexp = re.compile(regexp)
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
def assert_produces_warning(expected_warning=Warning, filter_level="always"):
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
        saw_warning = False
        warnings.simplefilter(filter_level)
        yield w
        extra_warnings = []
        for actual_warning in w:
            if (expected_warning and issubclass(actual_warning.category,
                                                expected_warning)):
                saw_warning = True
            else:
                extra_warnings.append(actual_warning.category.__name__)
        if expected_warning:
            assert saw_warning, ("Did not see expected warning of class %r."
                                 % expected_warning.__name__)
        assert not extra_warnings, ("Caused unexpected warning(s): %r."
                                    % extra_warnings)


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
