from __future__ import division

# pylint: disable-msg=W0402

import random
import string
import sys
import tempfile
import warnings
import inspect
import os

from datetime import datetime
from functools import wraps
from contextlib import contextmanager, closing
from httplib import HTTPException
from urllib2 import urlopen
from distutils.version import LooseVersion

from numpy.random import randn
import numpy as np

from pandas.core.common import isnull, _is_sequence
import pandas.core.index as index
import pandas.core.series as series
import pandas.core.frame as frame
import pandas.core.panel as panel
import pandas.core.panel4d as panel4d

from pandas import bdate_range
from pandas.tseries.index import DatetimeIndex
from pandas.tseries.period import PeriodIndex

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


def rands(n):
    choices = string.ascii_letters + string.digits
    return ''.join(random.choice(choices) for _ in xrange(n))


def randu(n):
    choices = u"".join(map(unichr, range(1488, 1488 + 26))) + string.digits
    return ''.join([random.choice(choices) for _ in xrange(n)])

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
def ensure_clean(filename = None):
    # if we are not passed a filename, generate a temporary
    if filename is None:
        filename = tempfile.mkstemp()[1]

    try:
        yield filename
    finally:
        try:
            os.remove(filename)
        except:
            pass


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


def isiterable(obj):
    return hasattr(obj, '__iter__')


def assert_almost_equal(a, b, check_less_precise = False):
    if isinstance(a, dict) or isinstance(b, dict):
        return assert_dict_equal(a, b)

    if isinstance(a, basestring):
        assert a == b, "%s != %s" % (a, b)
        return True

    if isiterable(a):
        np.testing.assert_(isiterable(b))
        na, nb = len(a), len(b)
        assert na == nb, "%s != %s" % (na, nb)

        if np.array_equal(a, b):
            return True
        else:
            for i in xrange(na):
                assert_almost_equal(a[i], b[i], check_less_precise)
        return True

    err_msg = lambda a, b: 'expected %.5f but got %.5f' % (b, a)

    if isnull(a):
        np.testing.assert_(isnull(b))
        return

    if isinstance(a, (bool, float, int, np.float32)):
        decimal = 5

        # deal with differing dtypes
        if check_less_precise:
            dtype_a = np.dtype(type(a))
            dtype_b = np.dtype(type(b))
            if dtype_a.kind == 'f' and dtype_b == 'f':
                if dtype_a.itemsize <= 4 and dtype_b.itemsize <= 4:
                    decimal = 3

        if np.isinf(a):
            assert np.isinf(b), err_msg(a, b)

        # case for zero
        elif abs(a) < 1e-5:
            np.testing.assert_almost_equal(
                a, b, decimal=decimal, err_msg=err_msg(a, b), verbose=False)
        else:
            np.testing.assert_almost_equal(
                1, a / b, decimal=decimal, err_msg=err_msg(a, b), verbose=False)
    else:
        assert a == b, "%s != %s" % (a, b)


def is_sorted(seq):
    return assert_almost_equal(seq, np.sort(np.array(seq)))


def assert_dict_equal(a, b, compare_keys=True):
    a_keys = frozenset(a.keys())
    b_keys = frozenset(b.keys())

    if compare_keys:
        assert(a_keys == b_keys)

    for k in a_keys:
        assert_almost_equal(a[k], b[k])


def assert_series_equal(left, right, check_dtype=True,
                        check_index_type=False,
                        check_index_freq=False,
                        check_series_type=False,
                        check_less_precise=False):
    if check_series_type:
        assert(type(left) == type(right))
    assert_almost_equal(left.values, right.values, check_less_precise)
    if check_dtype:
        assert(left.dtype == right.dtype)
    if check_less_precise:
        assert_almost_equal(left.index.values, right.index.values, check_less_precise)
    else:
        assert(left.index.equals(right.index))
    if check_index_type:
        assert(type(left.index) == type(right.index))
        assert(left.index.dtype == right.index.dtype)
        assert(left.index.inferred_type == right.index.inferred_type)
    if check_index_freq:
        assert(getattr(left, 'freqstr', None) ==
               getattr(right, 'freqstr', None))


def assert_frame_equal(left, right, check_dtype=True,
                       check_index_type=False,
                       check_column_type=False,
                       check_frame_type=False,
                       check_less_precise=False,
                       check_names=True):
    if check_frame_type:
        assert(type(left) == type(right))
    assert(isinstance(left, DataFrame))
    assert(isinstance(right, DataFrame))

    if check_less_precise:
        assert_almost_equal(left.columns,right.columns)
        assert_almost_equal(left.index,right.index)
    else:
        assert(left.columns.equals(right.columns))
        assert(left.index.equals(right.index))

    for i, col in enumerate(left.columns):
        assert(col in right)
        lcol = left.icol(i)
        rcol = right.icol(i)
        assert_series_equal(lcol, rcol,
                            check_dtype=check_dtype,
                            check_index_type=check_index_type,
                            check_less_precise=check_less_precise)

    if check_index_type:
        assert(type(left.index) == type(right.index))
        assert(left.index.dtype == right.index.dtype)
        assert(left.index.inferred_type == right.index.inferred_type)
    if check_column_type:
        assert(type(left.columns) == type(right.columns))
        assert(left.columns.dtype == right.columns.dtype)
        assert(left.columns.inferred_type == right.columns.inferred_type)
    if check_names:
        assert(left.index.names == right.index.names)
        assert(left.columns.names == right.columns.names)


def assert_panel_equal(left, right,
                       check_panel_type=False,
                       check_less_precise=False):
    if check_panel_type:
        assert(type(left) == type(right))

    assert(left.items.equals(right.items))
    assert(left.major_axis.equals(right.major_axis))
    assert(left.minor_axis.equals(right.minor_axis))

    for col, series in left.iterkv():
        assert(col in right)
        assert_frame_equal(series, right[col], check_less_precise=check_less_precise, check_names=False)  # TODO strangely check_names fails in py3 ?

    for col in right:
        assert(col in left)


def assert_panel4d_equal(left, right,
                         check_less_precise=False):
    assert(left.labels.equals(right.labels))
    assert(left.items.equals(right.items))
    assert(left.major_axis.equals(right.major_axis))
    assert(left.minor_axis.equals(right.minor_axis))

    for col, series in left.iterkv():
        assert(col in right)
        assert_panel_equal(series, right[col], check_less_precise=check_less_precise)

    for col in right:
        assert(col in left)


def assert_contains_all(iterable, dic):
    for k in iterable:
        assert(k in dic)


def getCols(k):
    return string.ascii_uppercase[:k]


def makeStringIndex(k):
    return Index([rands(10) for _ in xrange(k)])


def makeUnicodeIndex(k):
    return Index([randu(10) for _ in xrange(k)])


def makeIntIndex(k):
    return Index(range(k))


def makeFloatIndex(k):
    values = sorted(np.random.random_sample(k)) - np.random.random_sample(1)
    return Index(values * (10 ** np.random.randint(0, 9)))


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


def makeDataFrame():
    data = getSeriesData()
    return DataFrame(data)


def getArangeMat():
    return np.arange(N * K).reshape((N, K))


def getMixedTypeDict():
    index = Index(['a', 'b', 'c', 'd', 'e'])

    data = {
        'A': [0., 1., 2., 3., 4.],
        'B': [0., 1., 0., 1., 0.],
        'C': ['foo1', 'foo2', 'foo3', 'foo4', 'foo5'],
        'D': bdate_range('1/1/2009', periods=5)
    }

    return index, data


def makeDateIndex(k):
    dt = datetime(2000, 1, 1)
    dr = bdate_range(dt, periods=k)
    return DatetimeIndex(dr)


def makePeriodIndex(k):
    dt = datetime(2000, 1, 1)
    dr = PeriodIndex(start=dt, periods=k, freq='B')
    return dr


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


def makeTimeDataFrame(nper=None):
    data = getTimeSeriesData(nper)
    return DataFrame(data)


def getPeriodData():
    return dict((c, makePeriodSeries()) for c in getCols(K))


def makePeriodFrame():
    data = getPeriodData()
    return DataFrame(data)


def makePanel(nper=None):
    cols = ['Item' + c for c in string.ascii_uppercase[:K - 1]]
    data = dict((c, makeTimeDataFrame(nper)) for c in cols)
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

    from pandas.util.compat import Counter
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
    if isinstance(names, basestring) and nlevels == 1:
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
            numeric_tuple = re.sub("[^\d_]_?","",x).split("_")
            return map(int,numeric_tuple)

        # build a list of lists to create the index from
        div_factor = nentries // ndupe_l[i] + 1
        cnt = Counter()
        for j in range(div_factor):
            label = prefix + '_l%d_g' % i + str(j)
            cnt[label] = ndupe_l[i]
        # cute Counter trick
        result = list(sorted(cnt.elements(), key=keyfunc))[:nentries]
        tuples.append(result)

    tuples = zip(*tuples)

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


def add_nans_panel4d(panel4d):
    for l, label in enumerate(panel4d.labels):
        panel = panel4d[label]
        add_nans(panel)


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


_network_error_classes = IOError, HTTPException

@optional_args
def network(t, raise_on_error=_RAISE_NETWORK_ERROR_DEFAULT,
            error_classes=_network_error_classes):
    """
    Label a test as requiring network connection and skip test if it encounters a ``URLError``.

    In some cases it is not possible to assume network presence (e.g. Debian
    build hosts).

    You can pass an optional ``raise_on_error`` argument to the decorator, in
    which case it will always raise an error even if it's not a subclass of
    ``error_classes``.

    Parameters
    ----------
    t : callable
        The test requiring network connectivity.
    raise_on_error : bool
        If True, never catches errors.
    error_classes : iterable
        error classes to ignore. If not in ``error_classes``, raises the error.
        defaults to URLError. Be careful about changing the error classes here,
        it may result in undefined behavior.

    Returns
    -------
    t : callable
        The decorated test `t`.

    Examples
    --------
    A test can be decorated as requiring network like this::

      >>> from pandas.util.testing import network
      >>> import urllib2
      >>> import nose
      >>> @network
      ... def test_network():
      ...   urllib2.urlopen("rabbit://bonanza.com")
      ...
      >>> try:
      ...   test_network()
      ... except nose.SkipTest:
      ...   print "SKIPPING!"
      ...
      SKIPPING!

    Alternatively, you can use set ``raise_on_error`` in order to get
    the error to bubble up, e.g.::

      >>> @network(raise_on_error=True)
      ... def test_network():
      ...   urllib2.urlopen("complaint://deadparrot.com")
      ...
      >>> test_network()
      Traceback (most recent call last):
        ...
      URLError: <urlopen error unknown url type: complaint>

    And use ``nosetests -a '!network'`` to exclude running tests requiring
    network connectivity. ``_RAISE_NETWORK_ERROR_DEFAULT`` in
    ``pandas/util/testing.py`` sets the default behavior (currently False).
    """
    from nose import SkipTest
    t.network = True

    @wraps(t)
    def network_wrapper(*args, **kwargs):
        if raise_on_error:
            return t(*args, **kwargs)
        else:
            try:
                return t(*args, **kwargs)
            except error_classes as e:
                raise SkipTest("Skipping test %s" % e)

    return network_wrapper


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
def with_connectivity_check(t, url="http://www.google.com",
                            raise_on_error=_RAISE_NETWORK_ERROR_DEFAULT,
                            check_before_test=False,
                            error_classes=_network_error_classes):
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
        The url to test via ``urllib2.urlopen`` to check for connectivity.
        Defaults to 'http://www.google.com'.
    raise_on_error : bool
        If True, never catches errors.
    check_before_test : bool
        If True, checks connectivity before running the test case.
    error_classes : tuple or Exception
        error classes to ignore. If not in ``error_classes``, raises the error.
        defaults to IOError. Be careful about changing the error classes here.

    Notes
    -----
    * ``raise_on_error`` supercedes ``check_before_test``

    Returns
    -------
    t : callable
        The decorated test ``t``, with checks for connectivity errors.

    Example
    -------

    In this example, you see how it will raise the error if it can connect to
    the url::
        >>> @with_connectivity_check("http://www.yahoo.com")
        ... def test_something_with_yahoo():
        ...    raise IOError("Failure Message")
        >>> test_something_with_yahoo()
        Traceback (most recent call last):
            ...
        IOError: Failure Message

    I you set check_before_test, it will check the url first and not run the test on failure::
        >>> @with_connectivity_check("failing://url.blaher", check_before_test=True)
        ... def test_something():
        ...     print("I ran!")
        ...     raise ValueError("Failure")
        >>> test_something()
        Traceback (most recent call last):
            ...
        SkipTest
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
        except error_classes as e:
            if raise_on_error or can_connect(url, error_classes):
                raise
            else:
                raise SkipTest("Skipping test due to lack of connectivity"
                               " and error %s" % e)

    return wrapper


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


def assertRaisesRegexp(exception, regexp, callable, *args, **kwargs):
    """ Port of assertRaisesRegexp from unittest in Python 2.7 - used in with statement.

    Explanation from standard library:
        Like assertRaises() but also tests that regexp matches on the string
        representation of the raised exception. regexp may be a regular expression
        object or a string containing a regular expression suitable for use by
        re.search().

    You can pass either a regular expression or a compiled regular expression object.
    >>> assertRaisesRegexp(ValueError, 'invalid literal for.*XYZ',
    ...                                int, 'XYZ')
    >>> import re
    >>> assertRaisesRegexp(ValueError, re.compile('literal'), int, 'XYZ')

    If an exception of a different type is raised, it bubbles up.

    >>> assertRaisesRegexp(TypeError, 'literal', int, 'XYZ')
    Traceback (most recent call last):
        ...
    ValueError: invalid literal for int() with base 10: 'XYZ'
    >>> dct = {}
    >>> assertRaisesRegexp(KeyError, 'pear', dct.__getitem__, 'apple')
    Traceback (most recent call last):
        ...
    AssertionError: "pear" does not match "'apple'"
    >>> assertRaisesRegexp(KeyError, 'apple', dct.__getitem__, 'apple')
    >>> assertRaisesRegexp(Exception, 'operand type.*int.*dict', lambda : 2 + {})
    """

    import re

    try:
        callable(*args, **kwargs)
    except Exception as e:
        if not issubclass(e.__class__, exception):
            # mimics behavior of unittest
            raise
        # don't recompile
        if hasattr(regexp, "search"):
            expected_regexp = regexp
        else:
            expected_regexp = re.compile(regexp)
        if not expected_regexp.search(str(e)):
            raise AssertionError('"%s" does not match "%s"' %
                                 (expected_regexp.pattern, str(e)))
    else:
        # Apparently some exceptions don't have a __name__ attribute? Just aping unittest library here
        name = getattr(exception, "__name__", str(exception))
        raise AssertionError("{0} not raised".format(name))


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
