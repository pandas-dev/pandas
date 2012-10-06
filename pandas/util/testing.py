from __future__ import division

# pylint: disable-msg=W0402

from datetime import datetime
import random
import string
import sys

from contextlib import contextmanager  # contextlib is available since 2.5

from distutils.version import LooseVersion

from numpy.random import randn
import numpy as np

from pandas.core.common import isnull
import pandas.core.index as index
import pandas.core.series as series
import pandas.core.frame as frame
import pandas.core.panel as panel

from pandas import bdate_range
from pandas.tseries.index import DatetimeIndex
from pandas.tseries.period import PeriodIndex
from pandas.tseries.interval import IntervalIndex

Index = index.Index
Series = series.Series
DataFrame = frame.DataFrame
Panel = panel.Panel

N = 30
K = 4

def rands(n):
    choices = string.ascii_letters + string.digits
    return ''.join([random.choice(choices) for _ in xrange(n)])

def randu(n):
    choices = u"".join(map(unichr,range(1488,1488+26))) + string.digits
    return ''.join([random.choice(choices) for _ in xrange(n)])

#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
# Comparators

def equalContents(arr1, arr2):
    """Checks if the set of unique elements of arr1 and arr2 are equivalent.
    """
    return frozenset(arr1) == frozenset(arr2)

def isiterable(obj):
    return hasattr(obj, '__iter__')

def assert_almost_equal(a, b):
    if isinstance(a, dict) or isinstance(b, dict):
        return assert_dict_equal(a, b)

    if isinstance(a, basestring):
        assert a == b, (a, b)
        return True

    if isiterable(a):
        np.testing.assert_(isiterable(b))
        np.testing.assert_equal(len(a), len(b))
        if np.array_equal(a, b):
            return True
        else:
            for i in xrange(len(a)):
                assert_almost_equal(a[i], b[i])
        return True

    err_msg = lambda a, b: 'expected %.5f but got %.5f' % (a, b)

    if isnull(a):
        np.testing.assert_(isnull(b))
        return

    if isinstance(a, (bool, float, int)):
        # case for zero
        if abs(a) < 1e-5:
            np.testing.assert_almost_equal(
                a, b, decimal=5, err_msg=err_msg(a, b), verbose=False)
        else:
            np.testing.assert_almost_equal(
                1, a/b, decimal=5, err_msg=err_msg(a, b), verbose=False)
    else:
        assert(a == b)

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
                        check_series_type=False):
    if check_series_type:
        assert(type(left) == type(right))
    assert_almost_equal(left.values, right.values)
    if check_dtype:
        assert(left.dtype == right.dtype)
    assert(left.index.equals(right.index))
    if check_index_type:
        assert(type(left.index) == type(right.index))
        assert(left.index.dtype == right.index.dtype)
        assert(left.index.inferred_type == right.index.inferred_type)
    if check_index_freq:
        assert(getattr(left, 'freqstr', None) ==
               getattr(right, 'freqstr', None))

def assert_frame_equal(left, right, check_index_type=False,
                       check_column_type=False,
                       check_frame_type=False):
    if check_frame_type:
        assert(type(left) == type(right))
    assert(isinstance(left, DataFrame))
    assert(isinstance(right, DataFrame))
    for col, series in left.iterkv():
        assert(col in right)
        assert_series_equal(series, right[col])
    for col in right:
        assert(col in left)
    assert(left.index.equals(right.index))
    assert(left.columns.equals(right.columns))
    if check_index_type:
        assert(type(left.index) == type(right.index))
        assert(left.index.dtype == right.index.dtype)
        assert(left.index.inferred_type == right.index.inferred_type)
    if check_column_type:
        assert(type(left.columns) == type(right.columns))
        assert(left.columns.dtype == right.columns.dtype)
        assert(left.columns.inferred_type == right.columns.inferred_type)

def assert_panel_equal(left, right, check_panel_type=False):
    if check_panel_type:
        assert(type(left) == type(right))

    assert(left.items.equals(right.items))
    assert(left.major_axis.equals(right.major_axis))
    assert(left.minor_axis.equals(right.minor_axis))

    for col, series in left.iterkv():
        assert(col in right)
        assert_frame_equal(series, right[col])

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
        'A' : [0., 1., 2., 3., 4.],
        'B' : [0., 1., 0., 1., 0.],
        'C' : ['foo1', 'foo2', 'foo3', 'foo4', 'foo5'],
        'D' : bdate_range('1/1/2009', periods=5)
    }

    return index, data

def makeDateIndex(k):
    dt = datetime(2000,1,1)
    dr = bdate_range(dt, periods=k)
    return DatetimeIndex(dr)

def makePeriodIndex(k):
    dt = datetime(2000,1,1)
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

def getTimeSeriesData():
    return dict((c, makeTimeSeries()) for c in getCols(K))

def makeTimeDataFrame():
    data = getTimeSeriesData()
    return DataFrame(data)

def getPeriodData():
    return dict((c, makePeriodSeries()) for c in getCols(K))

def makePeriodFrame():
    data = getPeriodData()
    return DataFrame(data)

def makePanel():
    cols = ['Item' + c for c in string.ascii_uppercase[:K - 1]]
    data = dict((c, makeTimeDataFrame()) for c in cols)
    return Panel.fromDict(data)

def add_nans(panel):
    I, J, N = panel.shape
    for i, item in enumerate(panel.items):
        dm = panel[item]
        for j, col in enumerate(dm.columns):
            dm[col][:i + j] = np.NaN

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
def network(t):
    """
    Label a test as requiring network connection.

    In some cases it is not possible to assume network presence (e.g. Debian
    build hosts).

    Parameters
    ----------
    t : callable
        The test requiring network connectivity.

    Returns
    -------
    t : callable
        The decorated test `t`.

    Examples
    --------
    A test can be decorated as requiring network like this::

      from pandas.util.testing import *

      @network
      def test_network(self):
          print 'Fetch the stars from http://'

    And use ``nosetests -a '!network'`` to exclude running tests requiring
    network connectivity.
    """

    t.network = True
    return t


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
            attrs[k]=v
        self.attrs = attrs
        self.obj = obj

    def __getattribute__(self,name):
        attrs = object.__getattribute__(self, "attrs")
        obj = object.__getattribute__(self, "obj")
        return attrs.get(name, type(obj).__getattribute__(obj,name))

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
