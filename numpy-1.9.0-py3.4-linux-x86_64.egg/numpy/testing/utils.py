"""
Utility function to facilitate testing.

"""
from __future__ import division, absolute_import, print_function

import os
import sys
import re
import operator
import warnings
from functools import partial
import shutil
import contextlib
from tempfile import mkdtemp
from .nosetester import import_nose
from numpy.core import float32, empty, arange, array_repr, ndarray

if sys.version_info[0] >= 3:
    from io import StringIO
else:
    from StringIO import StringIO

__all__ = ['assert_equal', 'assert_almost_equal', 'assert_approx_equal',
           'assert_array_equal', 'assert_array_less', 'assert_string_equal',
           'assert_array_almost_equal', 'assert_raises', 'build_err_msg',
           'decorate_methods', 'jiffies', 'memusage', 'print_assert_equal',
           'raises', 'rand', 'rundocs', 'runstring', 'verbose', 'measure',
           'assert_', 'assert_array_almost_equal_nulp', 'assert_raises_regex',
           'assert_array_max_ulp', 'assert_warns', 'assert_no_warnings',
           'assert_allclose', 'IgnoreException']


verbose = 0


def assert_(val, msg='') :
    """
    Assert that works in release mode.
    Accepts callable msg to allow deferring evaluation until failure.

    The Python built-in ``assert`` does not work when executing code in
    optimized mode (the ``-O`` flag) - no byte-code is generated for it.

    For documentation on usage, refer to the Python documentation.

    """
    if not val :
        try:
            smsg = msg()
        except TypeError:
            smsg = msg
        raise AssertionError(smsg)

def gisnan(x):
    """like isnan, but always raise an error if type not supported instead of
    returning a TypeError object.

    Notes
    -----
    isnan and other ufunc sometimes return a NotImplementedType object instead
    of raising any exception. This function is a wrapper to make sure an
    exception is always raised.

    This should be removed once this problem is solved at the Ufunc level."""
    from numpy.core import isnan
    st = isnan(x)
    if isinstance(st, type(NotImplemented)):
        raise TypeError("isnan not supported for this type")
    return st

def gisfinite(x):
    """like isfinite, but always raise an error if type not supported instead of
    returning a TypeError object.

    Notes
    -----
    isfinite and other ufunc sometimes return a NotImplementedType object instead
    of raising any exception. This function is a wrapper to make sure an
    exception is always raised.

    This should be removed once this problem is solved at the Ufunc level."""
    from numpy.core import isfinite, errstate
    with errstate(invalid='ignore'):
        st = isfinite(x)
        if isinstance(st, type(NotImplemented)):
            raise TypeError("isfinite not supported for this type")
    return st

def gisinf(x):
    """like isinf, but always raise an error if type not supported instead of
    returning a TypeError object.

    Notes
    -----
    isinf and other ufunc sometimes return a NotImplementedType object instead
    of raising any exception. This function is a wrapper to make sure an
    exception is always raised.

    This should be removed once this problem is solved at the Ufunc level."""
    from numpy.core import isinf, errstate
    with errstate(invalid='ignore'):
        st = isinf(x)
        if isinstance(st, type(NotImplemented)):
            raise TypeError("isinf not supported for this type")
    return st

def rand(*args):
    """Returns an array of random numbers with the given shape.

    This only uses the standard library, so it is useful for testing purposes.
    """
    import random
    from numpy.core import zeros, float64
    results = zeros(args, float64)
    f = results.flat
    for i in range(len(f)):
        f[i] = random.random()
    return results

if sys.platform[:5]=='linux':
    def jiffies(_proc_pid_stat = '/proc/%s/stat'%(os.getpid()),
                _load_time=[]):
        """ Return number of jiffies (1/100ths of a second) that this
    process has been scheduled in user mode. See man 5 proc. """
        import time
        if not _load_time:
            _load_time.append(time.time())
        try:
            f=open(_proc_pid_stat, 'r')
            l = f.readline().split(' ')
            f.close()
            return int(l[13])
        except:
            return int(100*(time.time()-_load_time[0]))

    def memusage(_proc_pid_stat = '/proc/%s/stat'%(os.getpid())):
        """ Return virtual memory size in bytes of the running python.
        """
        try:
            f=open(_proc_pid_stat, 'r')
            l = f.readline().split(' ')
            f.close()
            return int(l[22])
        except:
            return
else:
    # os.getpid is not in all platforms available.
    # Using time is safe but inaccurate, especially when process
    # was suspended or sleeping.
    def jiffies(_load_time=[]):
        """ Return number of jiffies (1/100ths of a second) that this
    process has been scheduled in user mode. [Emulation with time.time]. """
        import time
        if not _load_time:
            _load_time.append(time.time())
        return int(100*(time.time()-_load_time[0]))
    def memusage():
        """ Return memory usage of running python. [Not implemented]"""
        raise NotImplementedError

if os.name=='nt' and sys.version[:3] > '2.3':
    # Code "stolen" from enthought/debug/memusage.py
    def GetPerformanceAttributes(object, counter, instance = None,
                                 inum=-1, format = None, machine=None):
        # NOTE: Many counters require 2 samples to give accurate results,
        # including "% Processor Time" (as by definition, at any instant, a
        # thread's CPU usage is either 0 or 100).  To read counters like this,
        # you should copy this function, but keep the counter open, and call
        # CollectQueryData() each time you need to know.
        # See http://msdn.microsoft.com/library/en-us/dnperfmo/html/perfmonpt2.asp
        # My older explanation for this was that the "AddCounter" process forced
        # the CPU to 100%, but the above makes more sense :)
        import win32pdh
        if format is None: format = win32pdh.PDH_FMT_LONG
        path = win32pdh.MakeCounterPath( (machine, object, instance, None, inum, counter) )
        hq = win32pdh.OpenQuery()
        try:
            hc = win32pdh.AddCounter(hq, path)
            try:
                win32pdh.CollectQueryData(hq)
                type, val = win32pdh.GetFormattedCounterValue(hc, format)
                return val
            finally:
                win32pdh.RemoveCounter(hc)
        finally:
            win32pdh.CloseQuery(hq)

    def memusage(processName="python", instance=0):
        # from win32pdhutil, part of the win32all package
        import win32pdh
        return GetPerformanceAttributes("Process", "Virtual Bytes",
                                        processName, instance,
                                        win32pdh.PDH_FMT_LONG, None)

def build_err_msg(arrays, err_msg, header='Items are not equal:',
                  verbose=True, names=('ACTUAL', 'DESIRED'), precision=8):
    msg = ['\n' + header]
    if err_msg:
        if err_msg.find('\n') == -1 and len(err_msg) < 79-len(header):
            msg = [msg[0] + ' ' + err_msg]
        else:
            msg.append(err_msg)
    if verbose:
        for i, a in enumerate(arrays):

            if isinstance(a, ndarray):
                # precision argument is only needed if the objects are ndarrays
                r_func = partial(array_repr, precision=precision)
            else:
                r_func = repr

            try:
                r = r_func(a)
            except:
                r = '[repr failed]'
            if r.count('\n') > 3:
                r = '\n'.join(r.splitlines()[:3])
                r += '...'
            msg.append(' %s: %s' % (names[i], r))
    return '\n'.join(msg)

def assert_equal(actual,desired,err_msg='',verbose=True):
    """
    Raises an AssertionError if two objects are not equal.

    Given two objects (scalars, lists, tuples, dictionaries or numpy arrays),
    check that all elements of these objects are equal. An exception is raised
    at the first conflicting values.

    Parameters
    ----------
    actual : array_like
        The object to check.
    desired : array_like
        The expected object.
    err_msg : str, optional
        The error message to be printed in case of failure.
    verbose : bool, optional
        If True, the conflicting values are appended to the error message.

    Raises
    ------
    AssertionError
        If actual and desired are not equal.

    Examples
    --------
    >>> np.testing.assert_equal([4,5], [4,6])
    ...
    <type 'exceptions.AssertionError'>:
    Items are not equal:
    item=1
     ACTUAL: 5
     DESIRED: 6

    """
    if isinstance(desired, dict):
        if not isinstance(actual, dict) :
            raise AssertionError(repr(type(actual)))
        assert_equal(len(actual), len(desired), err_msg, verbose)
        for k, i in desired.items():
            if k not in actual :
                raise AssertionError(repr(k))
            assert_equal(actual[k], desired[k], 'key=%r\n%s' % (k, err_msg), verbose)
        return
    if isinstance(desired, (list, tuple)) and isinstance(actual, (list, tuple)):
        assert_equal(len(actual), len(desired), err_msg, verbose)
        for k in range(len(desired)):
            assert_equal(actual[k], desired[k], 'item=%r\n%s' % (k, err_msg), verbose)
        return
    from numpy.core import ndarray, isscalar, signbit
    from numpy.lib import iscomplexobj, real, imag
    if isinstance(actual, ndarray) or isinstance(desired, ndarray):
        return assert_array_equal(actual, desired, err_msg, verbose)
    msg = build_err_msg([actual, desired], err_msg, verbose=verbose)

    # Handle complex numbers: separate into real/imag to handle
    # nan/inf/negative zero correctly
    # XXX: catch ValueError for subclasses of ndarray where iscomplex fail
    try:
        usecomplex = iscomplexobj(actual) or iscomplexobj(desired)
    except ValueError:
        usecomplex = False

    if usecomplex:
        if iscomplexobj(actual):
            actualr = real(actual)
            actuali = imag(actual)
        else:
            actualr = actual
            actuali = 0
        if iscomplexobj(desired):
            desiredr = real(desired)
            desiredi = imag(desired)
        else:
            desiredr = desired
            desiredi = 0
        try:
            assert_equal(actualr, desiredr)
            assert_equal(actuali, desiredi)
        except AssertionError:
            raise AssertionError(msg)

    # Inf/nan/negative zero handling
    try:
        # isscalar test to check cases such as [np.nan] != np.nan
        if isscalar(desired) != isscalar(actual):
            raise AssertionError(msg)

        # If one of desired/actual is not finite, handle it specially here:
        # check that both are nan if any is a nan, and test for equality
        # otherwise
        if not (gisfinite(desired) and gisfinite(actual)):
            isdesnan = gisnan(desired)
            isactnan = gisnan(actual)
            if isdesnan or isactnan:
                if not (isdesnan and isactnan):
                    raise AssertionError(msg)
            else:
                if not desired == actual:
                    raise AssertionError(msg)
            return
        elif desired == 0 and actual == 0:
            if not signbit(desired) == signbit(actual):
                raise AssertionError(msg)
    # If TypeError or ValueError raised while using isnan and co, just handle
    # as before
    except (TypeError, ValueError, NotImplementedError):
        pass

    # Explicitly use __eq__ for comparison, ticket #2552
    if not (desired == actual):
        raise AssertionError(msg)

def print_assert_equal(test_string, actual, desired):
    """
    Test if two objects are equal, and print an error message if test fails.

    The test is performed with ``actual == desired``.

    Parameters
    ----------
    test_string : str
        The message supplied to AssertionError.
    actual : object
        The object to test for equality against `desired`.
    desired : object
        The expected result.

    Examples
    --------
    >>> np.testing.print_assert_equal('Test XYZ of func xyz', [0, 1], [0, 1])
    >>> np.testing.print_assert_equal('Test XYZ of func xyz', [0, 1], [0, 2])
    Traceback (most recent call last):
    ...
    AssertionError: Test XYZ of func xyz failed
    ACTUAL:
    [0, 1]
    DESIRED:
    [0, 2]

    """
    import pprint

    if not (actual == desired):
        msg = StringIO()
        msg.write(test_string)
        msg.write(' failed\nACTUAL: \n')
        pprint.pprint(actual, msg)
        msg.write('DESIRED: \n')
        pprint.pprint(desired, msg)
        raise AssertionError(msg.getvalue())

def assert_almost_equal(actual,desired,decimal=7,err_msg='',verbose=True):
    """
    Raises an AssertionError if two items are not equal up to desired
    precision.

    .. note:: It is recommended to use one of `assert_allclose`,
              `assert_array_almost_equal_nulp` or `assert_array_max_ulp`
              instead of this function for more consistent floating point
              comparisons.

    The test is equivalent to ``abs(desired-actual) < 0.5 * 10**(-decimal)``.

    Given two objects (numbers or ndarrays), check that all elements of these
    objects are almost equal. An exception is raised at conflicting values.
    For ndarrays this delegates to assert_array_almost_equal

    Parameters
    ----------
    actual : array_like
        The object to check.
    desired : array_like
        The expected object.
    decimal : int, optional
        Desired precision, default is 7.
    err_msg : str, optional
        The error message to be printed in case of failure.
    verbose : bool, optional
        If True, the conflicting values are appended to the error message.

    Raises
    ------
    AssertionError
      If actual and desired are not equal up to specified precision.

    See Also
    --------
    assert_allclose: Compare two array_like objects for equality with desired
                     relative and/or absolute precision.
    assert_array_almost_equal_nulp, assert_array_max_ulp, assert_equal

    Examples
    --------
    >>> import numpy.testing as npt
    >>> npt.assert_almost_equal(2.3333333333333, 2.33333334)
    >>> npt.assert_almost_equal(2.3333333333333, 2.33333334, decimal=10)
    ...
    <type 'exceptions.AssertionError'>:
    Items are not equal:
     ACTUAL: 2.3333333333333002
     DESIRED: 2.3333333399999998

    >>> npt.assert_almost_equal(np.array([1.0,2.3333333333333]),
    ...                         np.array([1.0,2.33333334]), decimal=9)
    ...
    <type 'exceptions.AssertionError'>:
    Arrays are not almost equal
    <BLANKLINE>
    (mismatch 50.0%)
     x: array([ 1.        ,  2.33333333])
     y: array([ 1.        ,  2.33333334])

    """
    from numpy.core import ndarray
    from numpy.lib import iscomplexobj, real, imag

    # Handle complex numbers: separate into real/imag to handle
    # nan/inf/negative zero correctly
    # XXX: catch ValueError for subclasses of ndarray where iscomplex fail
    try:
        usecomplex = iscomplexobj(actual) or iscomplexobj(desired)
    except ValueError:
        usecomplex = False

    def _build_err_msg():
        header = ('Arrays are not almost equal to %d decimals' % decimal)
        return build_err_msg([actual, desired], err_msg, verbose=verbose,
                             header=header)

    if usecomplex:
        if iscomplexobj(actual):
            actualr = real(actual)
            actuali = imag(actual)
        else:
            actualr = actual
            actuali = 0
        if iscomplexobj(desired):
            desiredr = real(desired)
            desiredi = imag(desired)
        else:
            desiredr = desired
            desiredi = 0
        try:
            assert_almost_equal(actualr, desiredr, decimal=decimal)
            assert_almost_equal(actuali, desiredi, decimal=decimal)
        except AssertionError:
            raise AssertionError(_build_err_msg())

    if isinstance(actual, (ndarray, tuple, list)) \
            or isinstance(desired, (ndarray, tuple, list)):
        return assert_array_almost_equal(actual, desired, decimal, err_msg)
    try:
        # If one of desired/actual is not finite, handle it specially here:
        # check that both are nan if any is a nan, and test for equality
        # otherwise
        if not (gisfinite(desired) and gisfinite(actual)):
            if gisnan(desired) or gisnan(actual):
                if not (gisnan(desired) and gisnan(actual)):
                    raise AssertionError(_build_err_msg())
            else:
                if not desired == actual:
                    raise AssertionError(_build_err_msg())
            return
    except (NotImplementedError, TypeError):
        pass
    if round(abs(desired - actual), decimal) != 0 :
        raise AssertionError(_build_err_msg())


def assert_approx_equal(actual,desired,significant=7,err_msg='',verbose=True):
    """
    Raises an AssertionError if two items are not equal up to significant
    digits.

    .. note:: It is recommended to use one of `assert_allclose`,
              `assert_array_almost_equal_nulp` or `assert_array_max_ulp`
              instead of this function for more consistent floating point
              comparisons.

    Given two numbers, check that they are approximately equal.
    Approximately equal is defined as the number of significant digits
    that agree.

    Parameters
    ----------
    actual : scalar
        The object to check.
    desired : scalar
        The expected object.
    significant : int, optional
        Desired precision, default is 7.
    err_msg : str, optional
        The error message to be printed in case of failure.
    verbose : bool, optional
        If True, the conflicting values are appended to the error message.

    Raises
    ------
    AssertionError
      If actual and desired are not equal up to specified precision.

    See Also
    --------
    assert_allclose: Compare two array_like objects for equality with desired
                     relative and/or absolute precision.
    assert_array_almost_equal_nulp, assert_array_max_ulp, assert_equal

    Examples
    --------
    >>> np.testing.assert_approx_equal(0.12345677777777e-20, 0.1234567e-20)
    >>> np.testing.assert_approx_equal(0.12345670e-20, 0.12345671e-20,
                                       significant=8)
    >>> np.testing.assert_approx_equal(0.12345670e-20, 0.12345672e-20,
                                       significant=8)
    ...
    <type 'exceptions.AssertionError'>:
    Items are not equal to 8 significant digits:
     ACTUAL: 1.234567e-021
     DESIRED: 1.2345672000000001e-021

    the evaluated condition that raises the exception is

    >>> abs(0.12345670e-20/1e-21 - 0.12345672e-20/1e-21) >= 10**-(8-1)
    True

    """
    import numpy as np

    (actual, desired) = map(float, (actual, desired))
    if desired==actual:
        return
    # Normalized the numbers to be in range (-10.0,10.0)
    # scale = float(pow(10,math.floor(math.log10(0.5*(abs(desired)+abs(actual))))))
    with np.errstate(invalid='ignore'):
        scale = 0.5*(np.abs(desired) + np.abs(actual))
        scale = np.power(10, np.floor(np.log10(scale)))
    try:
        sc_desired = desired/scale
    except ZeroDivisionError:
        sc_desired = 0.0
    try:
        sc_actual = actual/scale
    except ZeroDivisionError:
        sc_actual = 0.0
    msg = build_err_msg([actual, desired], err_msg,
                header='Items are not equal to %d significant digits:' %
                                 significant,
                verbose=verbose)
    try:
        # If one of desired/actual is not finite, handle it specially here:
        # check that both are nan if any is a nan, and test for equality
        # otherwise
        if not (gisfinite(desired) and gisfinite(actual)):
            if gisnan(desired) or gisnan(actual):
                if not (gisnan(desired) and gisnan(actual)):
                    raise AssertionError(msg)
            else:
                if not desired == actual:
                    raise AssertionError(msg)
            return
    except (TypeError, NotImplementedError):
        pass
    if np.abs(sc_desired - sc_actual) >= np.power(10., -(significant-1)) :
        raise AssertionError(msg)

def assert_array_compare(comparison, x, y, err_msg='', verbose=True,
                         header='', precision=6):
    from numpy.core import array, isnan, isinf, any, all, inf
    x = array(x, copy=False, subok=True)
    y = array(y, copy=False, subok=True)

    def isnumber(x):
        return x.dtype.char in '?bhilqpBHILQPefdgFDG'

    def chk_same_position(x_id, y_id, hasval='nan'):
        """Handling nan/inf: check that x and y have the nan/inf at the same
        locations."""
        try:
            assert_array_equal(x_id, y_id)
        except AssertionError:
            msg = build_err_msg([x, y],
                                err_msg + '\nx and y %s location mismatch:' \
                                % (hasval), verbose=verbose, header=header,
                                names=('x', 'y'), precision=precision)
            raise AssertionError(msg)

    try:
        cond = (x.shape==() or y.shape==()) or x.shape == y.shape
        if not cond:
            msg = build_err_msg([x, y],
                                err_msg
                                + '\n(shapes %s, %s mismatch)' % (x.shape,
                                                                  y.shape),
                                verbose=verbose, header=header,
                                names=('x', 'y'), precision=precision)
            if not cond :
                raise AssertionError(msg)

        if isnumber(x) and isnumber(y):
            x_isnan, y_isnan = isnan(x), isnan(y)
            x_isinf, y_isinf = isinf(x), isinf(y)

            # Validate that the special values are in the same place
            if any(x_isnan) or any(y_isnan):
                chk_same_position(x_isnan, y_isnan, hasval='nan')
            if any(x_isinf) or any(y_isinf):
                # Check +inf and -inf separately, since they are different
                chk_same_position(x == +inf, y == +inf, hasval='+inf')
                chk_same_position(x == -inf, y == -inf, hasval='-inf')

            # Combine all the special values
            x_id, y_id = x_isnan, y_isnan
            x_id |= x_isinf
            y_id |= y_isinf

            # Only do the comparison if actual values are left
            if all(x_id):
                return

            if any(x_id):
                val = comparison(x[~x_id], y[~y_id])
            else:
                val = comparison(x, y)
        else:
            val = comparison(x, y)

        if isinstance(val, bool):
            cond = val
            reduced = [0]
        else:
            reduced = val.ravel()
            cond = reduced.all()
            reduced = reduced.tolist()
        if not cond:
            match = 100-100.0*reduced.count(1)/len(reduced)
            msg = build_err_msg([x, y],
                                err_msg
                                + '\n(mismatch %s%%)' % (match,),
                                verbose=verbose, header=header,
                                names=('x', 'y'), precision=precision)
            if not cond :
                raise AssertionError(msg)
    except ValueError as e:
        import traceback
        efmt = traceback.format_exc()
        header = 'error during assertion:\n\n%s\n\n%s' % (efmt, header)

        msg = build_err_msg([x, y], err_msg, verbose=verbose, header=header,
                            names=('x', 'y'), precision=precision)
        raise ValueError(msg)

def assert_array_equal(x, y, err_msg='', verbose=True):
    """
    Raises an AssertionError if two array_like objects are not equal.

    Given two array_like objects, check that the shape is equal and all
    elements of these objects are equal. An exception is raised at
    shape mismatch or conflicting values. In contrast to the standard usage
    in numpy, NaNs are compared like numbers, no assertion is raised if
    both objects have NaNs in the same positions.

    The usual caution for verifying equality with floating point numbers is
    advised.

    Parameters
    ----------
    x : array_like
        The actual object to check.
    y : array_like
        The desired, expected object.
    err_msg : str, optional
        The error message to be printed in case of failure.
    verbose : bool, optional
        If True, the conflicting values are appended to the error message.

    Raises
    ------
    AssertionError
        If actual and desired objects are not equal.

    See Also
    --------
    assert_allclose: Compare two array_like objects for equality with desired
                     relative and/or absolute precision.
    assert_array_almost_equal_nulp, assert_array_max_ulp, assert_equal

    Examples
    --------
    The first assert does not raise an exception:

    >>> np.testing.assert_array_equal([1.0,2.33333,np.nan],
    ...                               [np.exp(0),2.33333, np.nan])

    Assert fails with numerical inprecision with floats:

    >>> np.testing.assert_array_equal([1.0,np.pi,np.nan],
    ...                               [1, np.sqrt(np.pi)**2, np.nan])
    ...
    <type 'exceptions.ValueError'>:
    AssertionError:
    Arrays are not equal
    <BLANKLINE>
    (mismatch 50.0%)
     x: array([ 1.        ,  3.14159265,         NaN])
     y: array([ 1.        ,  3.14159265,         NaN])

    Use `assert_allclose` or one of the nulp (number of floating point values)
    functions for these cases instead:

    >>> np.testing.assert_allclose([1.0,np.pi,np.nan],
    ...                            [1, np.sqrt(np.pi)**2, np.nan],
    ...                            rtol=1e-10, atol=0)

    """
    assert_array_compare(operator.__eq__, x, y, err_msg=err_msg,
                         verbose=verbose, header='Arrays are not equal')

def assert_array_almost_equal(x, y, decimal=6, err_msg='', verbose=True):
    """
    Raises an AssertionError if two objects are not equal up to desired
    precision.

    .. note:: It is recommended to use one of `assert_allclose`,
              `assert_array_almost_equal_nulp` or `assert_array_max_ulp`
              instead of this function for more consistent floating point
              comparisons.

    The test verifies identical shapes and verifies values with
    ``abs(desired-actual) < 0.5 * 10**(-decimal)``.

    Given two array_like objects, check that the shape is equal and all
    elements of these objects are almost equal. An exception is raised at
    shape mismatch or conflicting values. In contrast to the standard usage
    in numpy, NaNs are compared like numbers, no assertion is raised if
    both objects have NaNs in the same positions.

    Parameters
    ----------
    x : array_like
        The actual object to check.
    y : array_like
        The desired, expected object.
    decimal : int, optional
        Desired precision, default is 6.
    err_msg : str, optional
      The error message to be printed in case of failure.
    verbose : bool, optional
        If True, the conflicting values are appended to the error message.

    Raises
    ------
    AssertionError
        If actual and desired are not equal up to specified precision.

    See Also
    --------
    assert_allclose: Compare two array_like objects for equality with desired
                     relative and/or absolute precision.
    assert_array_almost_equal_nulp, assert_array_max_ulp, assert_equal

    Examples
    --------
    the first assert does not raise an exception

    >>> np.testing.assert_array_almost_equal([1.0,2.333,np.nan],
                                             [1.0,2.333,np.nan])

    >>> np.testing.assert_array_almost_equal([1.0,2.33333,np.nan],
    ...                                      [1.0,2.33339,np.nan], decimal=5)
    ...
    <type 'exceptions.AssertionError'>:
    AssertionError:
    Arrays are not almost equal
    <BLANKLINE>
    (mismatch 50.0%)
     x: array([ 1.     ,  2.33333,      NaN])
     y: array([ 1.     ,  2.33339,      NaN])

    >>> np.testing.assert_array_almost_equal([1.0,2.33333,np.nan],
    ...                                      [1.0,2.33333, 5], decimal=5)
    <type 'exceptions.ValueError'>:
    ValueError:
    Arrays are not almost equal
     x: array([ 1.     ,  2.33333,      NaN])
     y: array([ 1.     ,  2.33333,  5.     ])

    """
    from numpy.core import around, number, float_, result_type, array
    from numpy.core.numerictypes import issubdtype
    from numpy.core.fromnumeric import any as npany
    def compare(x, y):
        try:
            if npany(gisinf(x)) or npany( gisinf(y)):
                xinfid = gisinf(x)
                yinfid = gisinf(y)
                if not xinfid == yinfid:
                    return False
                # if one item, x and y is +- inf
                if x.size == y.size == 1:
                    return x == y
                x = x[~xinfid]
                y = y[~yinfid]
        except (TypeError, NotImplementedError):
            pass

        # make sure y is an inexact type to avoid abs(MIN_INT); will cause
        # casting of x later.
        dtype = result_type(y, 1.)
        y = array(y, dtype=dtype, copy=False, subok=True)
        z = abs(x-y)

        if not issubdtype(z.dtype, number):
            z = z.astype(float_) # handle object arrays

        return around(z, decimal) <= 10.0**(-decimal)

    assert_array_compare(compare, x, y, err_msg=err_msg, verbose=verbose,
             header=('Arrays are not almost equal to %d decimals' % decimal),
             precision=decimal)


def assert_array_less(x, y, err_msg='', verbose=True):
    """
    Raises an AssertionError if two array_like objects are not ordered by less
    than.

    Given two array_like objects, check that the shape is equal and all
    elements of the first object are strictly smaller than those of the
    second object. An exception is raised at shape mismatch or incorrectly
    ordered values. Shape mismatch does not raise if an object has zero
    dimension. In contrast to the standard usage in numpy, NaNs are
    compared, no assertion is raised if both objects have NaNs in the same
    positions.



    Parameters
    ----------
    x : array_like
      The smaller object to check.
    y : array_like
      The larger object to compare.
    err_msg : string
      The error message to be printed in case of failure.
    verbose : bool
        If True, the conflicting values are appended to the error message.

    Raises
    ------
    AssertionError
      If actual and desired objects are not equal.

    See Also
    --------
    assert_array_equal: tests objects for equality
    assert_array_almost_equal: test objects for equality up to precision



    Examples
    --------
    >>> np.testing.assert_array_less([1.0, 1.0, np.nan], [1.1, 2.0, np.nan])
    >>> np.testing.assert_array_less([1.0, 1.0, np.nan], [1, 2.0, np.nan])
    ...
    <type 'exceptions.ValueError'>:
    Arrays are not less-ordered
    (mismatch 50.0%)
     x: array([  1.,   1.,  NaN])
     y: array([  1.,   2.,  NaN])

    >>> np.testing.assert_array_less([1.0, 4.0], 3)
    ...
    <type 'exceptions.ValueError'>:
    Arrays are not less-ordered
    (mismatch 50.0%)
     x: array([ 1.,  4.])
     y: array(3)

    >>> np.testing.assert_array_less([1.0, 2.0, 3.0], [4])
    ...
    <type 'exceptions.ValueError'>:
    Arrays are not less-ordered
    (shapes (3,), (1,) mismatch)
     x: array([ 1.,  2.,  3.])
     y: array([4])

    """
    assert_array_compare(operator.__lt__, x, y, err_msg=err_msg,
                         verbose=verbose,
                         header='Arrays are not less-ordered')

def runstring(astr, dict):
    exec(astr, dict)

def assert_string_equal(actual, desired):
    """
    Test if two strings are equal.

    If the given strings are equal, `assert_string_equal` does nothing.
    If they are not equal, an AssertionError is raised, and the diff
    between the strings is shown.

    Parameters
    ----------
    actual : str
        The string to test for equality against the expected string.
    desired : str
        The expected string.

    Examples
    --------
    >>> np.testing.assert_string_equal('abc', 'abc')
    >>> np.testing.assert_string_equal('abc', 'abcd')
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    ...
    AssertionError: Differences in strings:
    - abc+ abcd?    +

    """
    # delay import of difflib to reduce startup time
    import difflib

    if not isinstance(actual, str) :
        raise AssertionError(repr(type(actual)))
    if not isinstance(desired, str):
        raise AssertionError(repr(type(desired)))
    if re.match(r'\A'+desired+r'\Z', actual, re.M):
        return

    diff = list(difflib.Differ().compare(actual.splitlines(1), desired.splitlines(1)))
    diff_list = []
    while diff:
        d1 = diff.pop(0)
        if d1.startswith('  '):
            continue
        if d1.startswith('- '):
            l = [d1]
            d2 = diff.pop(0)
            if d2.startswith('? '):
                l.append(d2)
                d2 = diff.pop(0)
            if not d2.startswith('+ ') :
                raise AssertionError(repr(d2))
            l.append(d2)
            d3 = diff.pop(0)
            if d3.startswith('? '):
                l.append(d3)
            else:
                diff.insert(0, d3)
            if re.match(r'\A'+d2[2:]+r'\Z', d1[2:]):
                continue
            diff_list.extend(l)
            continue
        raise AssertionError(repr(d1))
    if not diff_list:
        return
    msg = 'Differences in strings:\n%s' % (''.join(diff_list)).rstrip()
    if actual != desired :
        raise AssertionError(msg)


def rundocs(filename=None, raise_on_error=True):
    """
    Run doctests found in the given file.

    By default `rundocs` raises an AssertionError on failure.

    Parameters
    ----------
    filename : str
        The path to the file for which the doctests are run.
    raise_on_error : bool
        Whether to raise an AssertionError when a doctest fails. Default is
        True.

    Notes
    -----
    The doctests can be run by the user/developer by adding the ``doctests``
    argument to the ``test()`` call. For example, to run all tests (including
    doctests) for `numpy.lib`:

    >>> np.lib.test(doctests=True) #doctest: +SKIP
    """
    import doctest, imp
    if filename is None:
        f = sys._getframe(1)
        filename = f.f_globals['__file__']
    name = os.path.splitext(os.path.basename(filename))[0]
    path = [os.path.dirname(filename)]
    file, pathname, description = imp.find_module(name, path)
    try:
        m = imp.load_module(name, file, pathname, description)
    finally:
        file.close()

    tests = doctest.DocTestFinder().find(m)
    runner = doctest.DocTestRunner(verbose=False)

    msg = []
    if raise_on_error:
        out = lambda s: msg.append(s)
    else:
        out = None

    for test in tests:
        runner.run(test, out=out)

    if runner.failures > 0 and raise_on_error:
        raise AssertionError("Some doctests failed:\n%s" % "\n".join(msg))


def raises(*args,**kwargs):
    nose = import_nose()
    return nose.tools.raises(*args,**kwargs)


def assert_raises(*args,**kwargs):
    """
    assert_raises(exception_class, callable, *args, **kwargs)

    Fail unless an exception of class exception_class is thrown
    by callable when invoked with arguments args and keyword
    arguments kwargs. If a different type of exception is
    thrown, it will not be caught, and the test case will be
    deemed to have suffered an error, exactly as for an
    unexpected exception.

    """
    nose = import_nose()
    return nose.tools.assert_raises(*args,**kwargs)


assert_raises_regex_impl = None


def assert_raises_regex(exception_class, expected_regexp,
                        callable_obj=None, *args, **kwargs):
    """
    Fail unless an exception of class exception_class and with message that
    matches expected_regexp is thrown by callable when invoked with arguments
    args and keyword arguments kwargs.

    Name of this function adheres to Python 3.2+ reference, but should work in
    all versions down to 2.6.

    """
    nose = import_nose()

    global assert_raises_regex_impl
    if assert_raises_regex_impl is None:
        try:
            # Python 3.2+
            assert_raises_regex_impl = nose.tools.assert_raises_regex
        except AttributeError:
            try:
                # 2.7+
                assert_raises_regex_impl = nose.tools.assert_raises_regexp
            except AttributeError:
                # 2.6

                # This class is copied from Python2.7 stdlib almost verbatim
                class _AssertRaisesContext(object):
                    """A context manager used to implement TestCase.assertRaises* methods."""

                    def __init__(self, expected, expected_regexp=None):
                        self.expected = expected
                        self.expected_regexp = expected_regexp

                    def failureException(self, msg):
                        return AssertionError(msg)

                    def __enter__(self):
                        return self

                    def __exit__(self, exc_type, exc_value, tb):
                        if exc_type is None:
                            try:
                                exc_name = self.expected.__name__
                            except AttributeError:
                                exc_name = str(self.expected)
                            raise self.failureException(
                                "{0} not raised".format(exc_name))
                        if not issubclass(exc_type, self.expected):
                            # let unexpected exceptions pass through
                            return False
                        self.exception = exc_value  # store for later retrieval
                        if self.expected_regexp is None:
                            return True

                        expected_regexp = self.expected_regexp
                        if isinstance(expected_regexp, basestring):
                            expected_regexp = re.compile(expected_regexp)
                        if not expected_regexp.search(str(exc_value)):
                            raise self.failureException(
                                '"%s" does not match "%s"' %
                                (expected_regexp.pattern, str(exc_value)))
                        return True

                def impl(cls, regex, callable_obj, *a, **kw):
                    mgr = _AssertRaisesContext(cls, regex)
                    if callable_obj is None:
                        return mgr
                    with mgr:
                        callable_obj(*a, **kw)
                assert_raises_regex_impl = impl

    return assert_raises_regex_impl(exception_class, expected_regexp,
                                    callable_obj, *args, **kwargs)


def decorate_methods(cls, decorator, testmatch=None):
    """
    Apply a decorator to all methods in a class matching a regular expression.

    The given decorator is applied to all public methods of `cls` that are
    matched by the regular expression `testmatch`
    (``testmatch.search(methodname)``). Methods that are private, i.e. start
    with an underscore, are ignored.

    Parameters
    ----------
    cls : class
        Class whose methods to decorate.
    decorator : function
        Decorator to apply to methods
    testmatch : compiled regexp or str, optional
        The regular expression. Default value is None, in which case the
        nose default (``re.compile(r'(?:^|[\\b_\\.%s-])[Tt]est' % os.sep)``)
        is used.
        If `testmatch` is a string, it is compiled to a regular expression
        first.

    """
    if testmatch is None:
        testmatch = re.compile(r'(?:^|[\\b_\\.%s-])[Tt]est' % os.sep)
    else:
        testmatch = re.compile(testmatch)
    cls_attr = cls.__dict__

    # delayed import to reduce startup time
    from inspect import isfunction

    methods = [_m for _m in cls_attr.values() if isfunction(_m)]
    for function in methods:
        try:
            if hasattr(function, 'compat_func_name'):
                funcname = function.compat_func_name
            else:
                funcname = function.__name__
        except AttributeError:
            # not a function
            continue
        if testmatch.search(funcname) and not funcname.startswith('_'):
            setattr(cls, funcname, decorator(function))
    return


def measure(code_str,times=1,label=None):
    """
    Return elapsed time for executing code in the namespace of the caller.

    The supplied code string is compiled with the Python builtin ``compile``.
    The precision of the timing is 10 milli-seconds. If the code will execute
    fast on this timescale, it can be executed many times to get reasonable
    timing accuracy.

    Parameters
    ----------
    code_str : str
        The code to be timed.
    times : int, optional
        The number of times the code is executed. Default is 1. The code is
        only compiled once.
    label : str, optional
        A label to identify `code_str` with. This is passed into ``compile``
        as the second argument (for run-time error messages).

    Returns
    -------
    elapsed : float
        Total elapsed time in seconds for executing `code_str` `times` times.

    Examples
    --------
    >>> etime = np.testing.measure('for i in range(1000): np.sqrt(i**2)',
    ...                            times=times)
    >>> print "Time for a single execution : ", etime / times, "s"
    Time for a single execution :  0.005 s

    """
    frame = sys._getframe(1)
    locs, globs = frame.f_locals, frame.f_globals

    code = compile(code_str,
                   'Test name: %s ' % label,
                   'exec')
    i = 0
    elapsed = jiffies()
    while i < times:
        i += 1
        exec(code, globs, locs)
    elapsed = jiffies() - elapsed
    return 0.01*elapsed

def _assert_valid_refcount(op):
    """
    Check that ufuncs don't mishandle refcount of object `1`.
    Used in a few regression tests.
    """
    import numpy as np
    a = np.arange(100 * 100)
    b = np.arange(100*100).reshape(100, 100)
    c = b

    i = 1

    rc = sys.getrefcount(i)
    for j in range(15):
        d = op(b, c)

    assert_(sys.getrefcount(i) >= rc)

def assert_allclose(actual, desired, rtol=1e-7, atol=0,
                    err_msg='', verbose=True):
    """
    Raises an AssertionError if two objects are not equal up to desired
    tolerance.

    The test is equivalent to ``allclose(actual, desired, rtol, atol)``.
    It compares the difference between `actual` and `desired` to
    ``atol + rtol * abs(desired)``.

    .. versionadded:: 1.5.0

    Parameters
    ----------
    actual : array_like
        Array obtained.
    desired : array_like
        Array desired.
    rtol : float, optional
        Relative tolerance.
    atol : float, optional
        Absolute tolerance.
    err_msg : str, optional
        The error message to be printed in case of failure.
    verbose : bool, optional
        If True, the conflicting values are appended to the error message.

    Raises
    ------
    AssertionError
        If actual and desired are not equal up to specified precision.

    See Also
    --------
    assert_array_almost_equal_nulp, assert_array_max_ulp

    Examples
    --------
    >>> x = [1e-5, 1e-3, 1e-1]
    >>> y = np.arccos(np.cos(x))
    >>> assert_allclose(x, y, rtol=1e-5, atol=0)

    """
    import numpy as np
    def compare(x, y):
        return np.allclose(x, y, rtol=rtol, atol=atol)

    actual, desired = np.asanyarray(actual), np.asanyarray(desired)
    header = 'Not equal to tolerance rtol=%g, atol=%g' % (rtol, atol)
    assert_array_compare(compare, actual, desired, err_msg=str(err_msg),
                         verbose=verbose, header=header)

def assert_array_almost_equal_nulp(x, y, nulp=1):
    """
    Compare two arrays relatively to their spacing.

    This is a relatively robust method to compare two arrays whose amplitude
    is variable.

    Parameters
    ----------
    x, y : array_like
        Input arrays.
    nulp : int, optional
        The maximum number of unit in the last place for tolerance (see Notes).
        Default is 1.

    Returns
    -------
    None

    Raises
    ------
    AssertionError
        If the spacing between `x` and `y` for one or more elements is larger
        than `nulp`.

    See Also
    --------
    assert_array_max_ulp : Check that all items of arrays differ in at most
        N Units in the Last Place.
    spacing : Return the distance between x and the nearest adjacent number.

    Notes
    -----
    An assertion is raised if the following condition is not met::

        abs(x - y) <= nulps * spacing(max(abs(x), abs(y)))

    Examples
    --------
    >>> x = np.array([1., 1e-10, 1e-20])
    >>> eps = np.finfo(x.dtype).eps
    >>> np.testing.assert_array_almost_equal_nulp(x, x*eps/2 + x)

    >>> np.testing.assert_array_almost_equal_nulp(x, x*eps + x)
    Traceback (most recent call last):
      ...
    AssertionError: X and Y are not equal to 1 ULP (max is 2)

    """
    import numpy as np
    ax = np.abs(x)
    ay = np.abs(y)
    ref = nulp * np.spacing(np.where(ax > ay, ax, ay))
    if not np.all(np.abs(x-y) <= ref):
        if np.iscomplexobj(x) or np.iscomplexobj(y):
            msg = "X and Y are not equal to %d ULP" % nulp
        else:
            max_nulp = np.max(nulp_diff(x, y))
            msg = "X and Y are not equal to %d ULP (max is %g)" % (nulp, max_nulp)
        raise AssertionError(msg)

def assert_array_max_ulp(a, b, maxulp=1, dtype=None):
    """
    Check that all items of arrays differ in at most N Units in the Last Place.

    Parameters
    ----------
    a, b : array_like
        Input arrays to be compared.
    maxulp : int, optional
        The maximum number of units in the last place that elements of `a` and
        `b` can differ. Default is 1.
    dtype : dtype, optional
        Data-type to convert `a` and `b` to if given. Default is None.

    Returns
    -------
    ret : ndarray
        Array containing number of representable floating point numbers between
        items in `a` and `b`.

    Raises
    ------
    AssertionError
        If one or more elements differ by more than `maxulp`.

    See Also
    --------
    assert_array_almost_equal_nulp : Compare two arrays relatively to their
        spacing.

    Examples
    --------
    >>> a = np.linspace(0., 1., 100)
    >>> res = np.testing.assert_array_max_ulp(a, np.arcsin(np.sin(a)))

    """
    import numpy as np
    ret = nulp_diff(a, b, dtype)
    if not np.all(ret <= maxulp):
        raise AssertionError("Arrays are not almost equal up to %g ULP" % \
                             maxulp)
    return ret

def nulp_diff(x, y, dtype=None):
    """For each item in x and y, return the number of representable floating
    points between them.

    Parameters
    ----------
    x : array_like
        first input array
    y : array_like
        second input array

    Returns
    -------
    nulp : array_like
        number of representable floating point numbers between each item in x
        and y.

    Examples
    --------
    # By definition, epsilon is the smallest number such as 1 + eps != 1, so
    # there should be exactly one ULP between 1 and 1 + eps
    >>> nulp_diff(1, 1 + np.finfo(x.dtype).eps)
    1.0
    """
    import numpy as np
    if dtype:
        x = np.array(x, dtype=dtype)
        y = np.array(y, dtype=dtype)
    else:
        x = np.array(x)
        y = np.array(y)

    t = np.common_type(x, y)
    if np.iscomplexobj(x) or np.iscomplexobj(y):
        raise NotImplementedError("_nulp not implemented for complex array")

    x = np.array(x, dtype=t)
    y = np.array(y, dtype=t)

    if not x.shape == y.shape:
        raise ValueError("x and y do not have the same shape: %s - %s" % \
                         (x.shape, y.shape))

    def _diff(rx, ry, vdt):
        diff = np.array(rx-ry, dtype=vdt)
        return np.abs(diff)

    rx = integer_repr(x)
    ry = integer_repr(y)
    return _diff(rx, ry, t)

def _integer_repr(x, vdt, comp):
    # Reinterpret binary representation of the float as sign-magnitude:
    # take into account two-complement representation
    # See also
    # http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
    rx = x.view(vdt)
    if not (rx.size == 1):
        rx[rx < 0] = comp - rx[rx<0]
    else:
        if rx < 0:
            rx = comp - rx

    return rx

def integer_repr(x):
    """Return the signed-magnitude interpretation of the binary representation of
    x."""
    import numpy as np
    if x.dtype == np.float32:
        return _integer_repr(x, np.int32, np.int32(-2**31))
    elif x.dtype == np.float64:
        return _integer_repr(x, np.int64, np.int64(-2**63))
    else:
        raise ValueError("Unsupported dtype %s" % x.dtype)

# The following two classes are copied from python 2.6 warnings module (context
# manager)
class WarningMessage(object):

    """
    Holds the result of a single showwarning() call.

    Deprecated in 1.8.0

    Notes
    -----
    `WarningMessage` is copied from the Python 2.6 warnings module,
    so it can be used in NumPy with older Python versions.

    """

    _WARNING_DETAILS = ("message", "category", "filename", "lineno", "file",
                        "line")

    def __init__(self, message, category, filename, lineno, file=None,
                    line=None):
        local_values = locals()
        for attr in self._WARNING_DETAILS:
            setattr(self, attr, local_values[attr])
        if category:
            self._category_name = category.__name__
        else:
            self._category_name = None

    def __str__(self):
        return ("{message : %r, category : %r, filename : %r, lineno : %s, "
                    "line : %r}" % (self.message, self._category_name,
                                    self.filename, self.lineno, self.line))

class WarningManager(object):
    """
    A context manager that copies and restores the warnings filter upon
    exiting the context.

    The 'record' argument specifies whether warnings should be captured by a
    custom implementation of ``warnings.showwarning()`` and be appended to a
    list returned by the context manager. Otherwise None is returned by the
    context manager. The objects appended to the list are arguments whose
    attributes mirror the arguments to ``showwarning()``.

    The 'module' argument is to specify an alternative module to the module
    named 'warnings' and imported under that name. This argument is only useful
    when testing the warnings module itself.

    Deprecated in 1.8.0

    Notes
    -----
    `WarningManager` is a copy of the ``catch_warnings`` context manager
    from the Python 2.6 warnings module, with slight modifications.
    It is copied so it can be used in NumPy with older Python versions.

    """
    def __init__(self, record=False, module=None):
        self._record = record
        if module is None:
            self._module = sys.modules['warnings']
        else:
            self._module = module
        self._entered = False

    def __enter__(self):
        if self._entered:
            raise RuntimeError("Cannot enter %r twice" % self)
        self._entered = True
        self._filters = self._module.filters
        self._module.filters = self._filters[:]
        self._showwarning = self._module.showwarning
        if self._record:
            log = []
            def showwarning(*args, **kwargs):
                log.append(WarningMessage(*args, **kwargs))
            self._module.showwarning = showwarning
            return log
        else:
            return None

    def __exit__(self):
        if not self._entered:
            raise RuntimeError("Cannot exit %r without entering first" % self)
        self._module.filters = self._filters
        self._module.showwarning = self._showwarning


def assert_warns(warning_class, func, *args, **kw):
    """
    Fail unless the given callable throws the specified warning.

    A warning of class warning_class should be thrown by the callable when
    invoked with arguments args and keyword arguments kwargs.
    If a different type of warning is thrown, it will not be caught, and the
    test case will be deemed to have suffered an error.

    .. versionadded:: 1.4.0

    Parameters
    ----------
    warning_class : class
        The class defining the warning that `func` is expected to throw.
    func : callable
        The callable to test.
    \\*args : Arguments
        Arguments passed to `func`.
    \\*\\*kwargs : Kwargs
        Keyword arguments passed to `func`.

    Returns
    -------
    The value returned by `func`.

    """
    with warnings.catch_warnings(record=True) as l:
        warnings.simplefilter('always')
        result = func(*args, **kw)
        if not len(l) > 0:
            raise AssertionError("No warning raised when calling %s"
                    % func.__name__)
        if not l[0].category is warning_class:
            raise AssertionError("First warning for %s is not a " \
                    "%s( is %s)" % (func.__name__, warning_class, l[0]))
    return result

def assert_no_warnings(func, *args, **kw):
    """
    Fail if the given callable produces any warnings.

    .. versionadded:: 1.7.0

    Parameters
    ----------
    func : callable
        The callable to test.
    \\*args : Arguments
        Arguments passed to `func`.
    \\*\\*kwargs : Kwargs
        Keyword arguments passed to `func`.

    Returns
    -------
    The value returned by `func`.

    """
    with warnings.catch_warnings(record=True) as l:
        warnings.simplefilter('always')
        result = func(*args, **kw)
        if len(l) > 0:
            raise AssertionError("Got warnings when calling %s: %s"
                    % (func.__name__, l))
    return result


def _gen_alignment_data(dtype=float32, type='binary', max_size=24):
    """
    generator producing data with different alignment and offsets
    to test simd vectorization

    Parameters
    ----------
    dtype : dtype
        data type to produce
    type : string
        'unary': create data for unary operations, creates one input
                 and output array
        'binary': create data for unary operations, creates two input
                 and output array
    max_size : integer
        maximum size of data to produce

    Returns
    -------
    if type is 'unary' yields one output, one input array and a message
    containing information on the data
    if type is 'binary' yields one output array, two input array and a message
    containing information on the data

    """
    ufmt = 'unary offset=(%d, %d), size=%d, dtype=%r, %s'
    bfmt = 'binary offset=(%d, %d, %d), size=%d, dtype=%r, %s'
    for o in range(3):
        for s in range(o + 2, max(o + 3, max_size)):
            if type == 'unary':
                inp = lambda : arange(s, dtype=dtype)[o:]
                out = empty((s,), dtype=dtype)[o:]
                yield out, inp(), ufmt % (o, o, s, dtype, 'out of place')
                yield inp(), inp(), ufmt % (o, o, s, dtype, 'in place')
                yield out[1:], inp()[:-1], ufmt % \
                    (o + 1, o, s - 1, dtype, 'out of place')
                yield out[:-1], inp()[1:], ufmt % \
                    (o, o + 1, s - 1, dtype, 'out of place')
                yield inp()[:-1], inp()[1:], ufmt % \
                    (o, o + 1, s - 1, dtype, 'aliased')
                yield inp()[1:], inp()[:-1], ufmt % \
                    (o + 1, o, s - 1, dtype, 'aliased')
            if type == 'binary':
                inp1 = lambda :arange(s, dtype=dtype)[o:]
                inp2 = lambda :arange(s, dtype=dtype)[o:]
                out = empty((s,), dtype=dtype)[o:]
                yield out, inp1(), inp2(),  bfmt % \
                    (o, o, o, s, dtype, 'out of place')
                yield inp1(), inp1(), inp2(), bfmt % \
                    (o, o, o, s, dtype, 'in place1')
                yield inp2(), inp1(), inp2(), bfmt % \
                    (o, o, o, s, dtype, 'in place2')
                yield out[1:], inp1()[:-1], inp2()[:-1], bfmt % \
                    (o + 1, o, o, s - 1, dtype, 'out of place')
                yield out[:-1], inp1()[1:], inp2()[:-1], bfmt % \
                    (o, o + 1, o, s - 1, dtype, 'out of place')
                yield out[:-1], inp1()[:-1], inp2()[1:], bfmt % \
                    (o, o, o + 1, s - 1, dtype, 'out of place')
                yield inp1()[1:], inp1()[:-1], inp2()[:-1], bfmt % \
                    (o + 1, o, o, s - 1, dtype, 'aliased')
                yield inp1()[:-1], inp1()[1:], inp2()[:-1], bfmt % \
                    (o, o + 1, o, s - 1, dtype, 'aliased')
                yield inp1()[:-1], inp1()[:-1], inp2()[1:], bfmt % \
                    (o, o, o + 1, s - 1, dtype, 'aliased')


class IgnoreException(Exception):
    "Ignoring this exception due to disabled feature"


@contextlib.contextmanager
def tempdir(*args, **kwargs):
    """Context manager to provide a temporary test folder.

    All arguments are passed as this to the underlying tempfile.mkdtemp
    function.

    """
    tmpdir = mkdtemp(*args, **kwargs)
    yield tmpdir
    shutil.rmtree(tmpdir)
