"""
Module of functions that are like ufuncs in acting on arrays and optionally
storing results in an output array.

"""
from __future__ import division, absolute_import, print_function

__all__ = ['fix', 'isneginf', 'isposinf']

import numpy.core.numeric as nx

def fix(x, y=None):
    """
    Round to nearest integer towards zero.

    Round an array of floats element-wise to nearest integer towards zero.
    The rounded values are returned as floats.

    Parameters
    ----------
    x : array_like
        An array of floats to be rounded
    y : ndarray, optional
        Output array

    Returns
    -------
    out : ndarray of floats
        The array of rounded numbers

    See Also
    --------
    trunc, floor, ceil
    around : Round to given number of decimals

    Examples
    --------
    >>> np.fix(3.14)
    3.0
    >>> np.fix(3)
    3.0
    >>> np.fix([2.1, 2.9, -2.1, -2.9])
    array([ 2.,  2., -2., -2.])

    """
    x = nx.asanyarray(x)
    y1 = nx.floor(x)
    y2 = nx.ceil(x)
    if y is None:
        y = nx.asanyarray(y1)
    y[...] = nx.where(x >= 0, y1, y2)
    return y

def isposinf(x, y=None):
    """
    Test element-wise for positive infinity, return result as bool array.

    Parameters
    ----------
    x : array_like
        The input array.
    y : array_like, optional
        A boolean array with the same shape as `x` to store the result.

    Returns
    -------
    y : ndarray
        A boolean array with the same dimensions as the input.
        If second argument is not supplied then a boolean array is returned
        with values True where the corresponding element of the input is
        positive infinity and values False where the element of the input is
        not positive infinity.

        If a second argument is supplied the result is stored there. If the
        type of that array is a numeric type the result is represented as zeros
        and ones, if the type is boolean then as False and True.
        The return value `y` is then a reference to that array.

    See Also
    --------
    isinf, isneginf, isfinite, isnan

    Notes
    -----
    Numpy uses the IEEE Standard for Binary Floating-Point for Arithmetic
    (IEEE 754).

    Errors result if the second argument is also supplied when `x` is a
    scalar input, or if first and second arguments have different shapes.

    Examples
    --------
    >>> np.isposinf(np.PINF)
    array(True, dtype=bool)
    >>> np.isposinf(np.inf)
    array(True, dtype=bool)
    >>> np.isposinf(np.NINF)
    array(False, dtype=bool)
    >>> np.isposinf([-np.inf, 0., np.inf])
    array([False, False,  True], dtype=bool)

    >>> x = np.array([-np.inf, 0., np.inf])
    >>> y = np.array([2, 2, 2])
    >>> np.isposinf(x, y)
    array([0, 0, 1])
    >>> y
    array([0, 0, 1])

    """
    if y is None:
        x = nx.asarray(x)
        y = nx.empty(x.shape, dtype=nx.bool_)
    nx.logical_and(nx.isinf(x), ~nx.signbit(x), y)
    return y

def isneginf(x, y=None):
    """
    Test element-wise for negative infinity, return result as bool array.

    Parameters
    ----------
    x : array_like
        The input array.
    y : array_like, optional
        A boolean array with the same shape and type as `x` to store the
        result.

    Returns
    -------
    y : ndarray
        A boolean array with the same dimensions as the input.
        If second argument is not supplied then a numpy boolean array is
        returned with values True where the corresponding element of the
        input is negative infinity and values False where the element of
        the input is not negative infinity.

        If a second argument is supplied the result is stored there. If the
        type of that array is a numeric type the result is represented as
        zeros and ones, if the type is boolean then as False and True. The
        return value `y` is then a reference to that array.

    See Also
    --------
    isinf, isposinf, isnan, isfinite

    Notes
    -----
    Numpy uses the IEEE Standard for Binary Floating-Point for Arithmetic
    (IEEE 754).

    Errors result if the second argument is also supplied when x is a scalar
    input, or if first and second arguments have different shapes.

    Examples
    --------
    >>> np.isneginf(np.NINF)
    array(True, dtype=bool)
    >>> np.isneginf(np.inf)
    array(False, dtype=bool)
    >>> np.isneginf(np.PINF)
    array(False, dtype=bool)
    >>> np.isneginf([-np.inf, 0., np.inf])
    array([ True, False, False], dtype=bool)

    >>> x = np.array([-np.inf, 0., np.inf])
    >>> y = np.array([2, 2, 2])
    >>> np.isneginf(x, y)
    array([1, 0, 0])
    >>> y
    array([1, 0, 0])

    """
    if y is None:
        x = nx.asarray(x)
        y = nx.empty(x.shape, dtype=nx.bool_)
    nx.logical_and(nx.isinf(x), nx.signbit(x), y)
    return y
