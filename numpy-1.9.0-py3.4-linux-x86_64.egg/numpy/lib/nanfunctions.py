"""
Functions that ignore NaN.

Functions
---------

- `nanmin` -- minimum non-NaN value
- `nanmax` -- maximum non-NaN value
- `nanargmin` -- index of minimum non-NaN value
- `nanargmax` -- index of maximum non-NaN value
- `nansum` -- sum of non-NaN values
- `nanmean` -- mean of non-NaN values
- `nanvar` -- variance of non-NaN values
- `nanstd` -- standard deviation of non-NaN values

"""
from __future__ import division, absolute_import, print_function

import warnings
import numpy as np
from numpy.lib.function_base import _ureduce as _ureduce

__all__ = [
    'nansum', 'nanmax', 'nanmin', 'nanargmax', 'nanargmin', 'nanmean',
    'nanmedian', 'nanpercentile', 'nanvar', 'nanstd'
    ]


def _replace_nan(a, val):
    """
    If `a` is of inexact type, make a copy of `a`, replace NaNs with
    the `val` value, and return the copy together with a boolean mask
    marking the locations where NaNs were present. If `a` is not of
    inexact type, do nothing and return `a` together with a mask of None.

    Parameters
    ----------
    a : array-like
        Input array.
    val : float
        NaN values are set to val before doing the operation.

    Returns
    -------
    y : ndarray
        If `a` is of inexact type, return a copy of `a` with the NaNs
        replaced by the fill value, otherwise return `a`.
    mask: {bool, None}
        If `a` is of inexact type, return a boolean mask marking locations of
        NaNs, otherwise return None.

    """
    is_new = not isinstance(a, np.ndarray)
    if is_new:
        a = np.array(a)
    if not issubclass(a.dtype.type, np.inexact):
        return a, None
    if not is_new:
        # need copy
        a = np.array(a, subok=True)

    mask = np.isnan(a)
    np.copyto(a, val, where=mask)
    return a, mask


def _copyto(a, val, mask):
    """
    Replace values in `a` with NaN where `mask` is True.  This differs from
    copyto in that it will deal with the case where `a` is a numpy scalar.

    Parameters
    ----------
    a : ndarray or numpy scalar
        Array or numpy scalar some of whose values are to be replaced
        by val.
    val : numpy scalar
        Value used a replacement.
    mask : ndarray, scalar
        Boolean array. Where True the corresponding element of `a` is
        replaced by `val`. Broadcasts.

    Returns
    -------
    res : ndarray, scalar
        Array with elements replaced or scalar `val`.

    """
    if isinstance(a, np.ndarray):
        np.copyto(a, val, where=mask, casting='unsafe')
    else:
        a = a.dtype.type(val)
    return a


def _divide_by_count(a, b, out=None):
    """
    Compute a/b ignoring invalid results. If `a` is an array the division
    is done in place. If `a` is a scalar, then its type is preserved in the
    output. If out is None, then then a is used instead so that the
    division is in place. Note that this is only called with `a` an inexact
    type.

    Parameters
    ----------
    a : {ndarray, numpy scalar}
        Numerator. Expected to be of inexact type but not checked.
    b : {ndarray, numpy scalar}
        Denominator.
    out : ndarray, optional
        Alternate output array in which to place the result.  The default
        is ``None``; if provided, it must have the same shape as the
        expected output, but the type will be cast if necessary.

    Returns
    -------
    ret : {ndarray, numpy scalar}
        The return value is a/b. If `a` was an ndarray the division is done
        in place. If `a` is a numpy scalar, the division preserves its type.

    """
    with np.errstate(invalid='ignore'):
        if isinstance(a, np.ndarray):
            if out is None:
                return np.divide(a, b, out=a, casting='unsafe')
            else:
                return np.divide(a, b, out=out, casting='unsafe')
        else:
            if out is None:
                return a.dtype.type(a / b)
            else:
                # This is questionable, but currently a numpy scalar can
                # be output to a zero dimensional array.
                return np.divide(a, b, out=out, casting='unsafe')


def nanmin(a, axis=None, out=None, keepdims=False):
    """
    Return minimum of an array or minimum along an axis, ignoring any NaNs.
    When all-NaN slices are encountered a ``RuntimeWarning`` is raised and
    Nan is returned for that slice.

    Parameters
    ----------
    a : array_like
        Array containing numbers whose minimum is desired. If `a` is not an
        array, a conversion is attempted.
    axis : int, optional
        Axis along which the minimum is computed. The default is to compute
        the minimum of the flattened array.
    out : ndarray, optional
        Alternate output array in which to place the result.  The default
        is ``None``; if provided, it must have the same shape as the
        expected output, but the type will be cast if necessary.  See
        `doc.ufuncs` for details.

        .. versionadded:: 1.8.0
    keepdims : bool, optional
        If this is set to True, the axes which are reduced are left in the
        result as dimensions with size one. With this option, the result
        will broadcast correctly against the original `a`.

        .. versionadded:: 1.8.0

    Returns
    -------
    nanmin : ndarray
        An array with the same shape as `a`, with the specified axis
        removed.  If `a` is a 0-d array, or if axis is None, an ndarray
        scalar is returned.  The same dtype as `a` is returned.

    See Also
    --------
    nanmax :
        The maximum value of an array along a given axis, ignoring any NaNs.
    amin :
        The minimum value of an array along a given axis, propagating any NaNs.
    fmin :
        Element-wise minimum of two arrays, ignoring any NaNs.
    minimum :
        Element-wise minimum of two arrays, propagating any NaNs.
    isnan :
        Shows which elements are Not a Number (NaN).
    isfinite:
        Shows which elements are neither NaN nor infinity.

    amax, fmax, maximum

    Notes
    -----
    Numpy uses the IEEE Standard for Binary Floating-Point for Arithmetic
    (IEEE 754). This means that Not a Number is not equivalent to infinity.
    Positive infinity is treated as a very large number and negative
    infinity is treated as a very small (i.e. negative) number.

    If the input has a integer type the function is equivalent to np.min.

    Examples
    --------
    >>> a = np.array([[1, 2], [3, np.nan]])
    >>> np.nanmin(a)
    1.0
    >>> np.nanmin(a, axis=0)
    array([ 1.,  2.])
    >>> np.nanmin(a, axis=1)
    array([ 1.,  3.])

    When positive infinity and negative infinity are present:

    >>> np.nanmin([1, 2, np.nan, np.inf])
    1.0
    >>> np.nanmin([1, 2, np.nan, np.NINF])
    -inf

    """
    if not isinstance(a, np.ndarray) or type(a) is np.ndarray:
        # Fast, but not safe for subclasses of ndarray
        res = np.fmin.reduce(a, axis=axis, out=out, keepdims=keepdims)
        if np.isnan(res).any():
            warnings.warn("All-NaN axis encountered", RuntimeWarning)
    else:
        # Slow, but safe for subclasses of ndarray
        a, mask = _replace_nan(a, +np.inf)
        res = np.amin(a, axis=axis, out=out, keepdims=keepdims)
        if mask is None:
            return res

        # Check for all-NaN axis
        mask = np.all(mask, axis=axis, keepdims=keepdims)
        if np.any(mask):
            res = _copyto(res, np.nan, mask)
            warnings.warn("All-NaN axis encountered", RuntimeWarning)
    return res


def nanmax(a, axis=None, out=None, keepdims=False):
    """
    Return the maximum of an array or maximum along an axis, ignoring any
    NaNs.  When all-NaN slices are encountered a ``RuntimeWarning`` is
    raised and NaN is returned for that slice.

    Parameters
    ----------
    a : array_like
        Array containing numbers whose maximum is desired. If `a` is not an
        array, a conversion is attempted.
    axis : int, optional
        Axis along which the maximum is computed. The default is to compute
        the maximum of the flattened array.
    out : ndarray, optional
        Alternate output array in which to place the result.  The default
        is ``None``; if provided, it must have the same shape as the
        expected output, but the type will be cast if necessary.  See
        `doc.ufuncs` for details.

        .. versionadded:: 1.8.0
    keepdims : bool, optional
        If this is set to True, the axes which are reduced are left in the
        result as dimensions with size one. With this option, the result
        will broadcast correctly against the original `a`.

        .. versionadded:: 1.8.0

    Returns
    -------
    nanmax : ndarray
        An array with the same shape as `a`, with the specified axis removed.
        If `a` is a 0-d array, or if axis is None, an ndarray scalar is
        returned.  The same dtype as `a` is returned.

    See Also
    --------
    nanmin :
        The minimum value of an array along a given axis, ignoring any NaNs.
    amax :
        The maximum value of an array along a given axis, propagating any NaNs.
    fmax :
        Element-wise maximum of two arrays, ignoring any NaNs.
    maximum :
        Element-wise maximum of two arrays, propagating any NaNs.
    isnan :
        Shows which elements are Not a Number (NaN).
    isfinite:
        Shows which elements are neither NaN nor infinity.

    amin, fmin, minimum

    Notes
    -----
    Numpy uses the IEEE Standard for Binary Floating-Point for Arithmetic
    (IEEE 754). This means that Not a Number is not equivalent to infinity.
    Positive infinity is treated as a very large number and negative
    infinity is treated as a very small (i.e. negative) number.

    If the input has a integer type the function is equivalent to np.max.

    Examples
    --------
    >>> a = np.array([[1, 2], [3, np.nan]])
    >>> np.nanmax(a)
    3.0
    >>> np.nanmax(a, axis=0)
    array([ 3.,  2.])
    >>> np.nanmax(a, axis=1)
    array([ 2.,  3.])

    When positive infinity and negative infinity are present:

    >>> np.nanmax([1, 2, np.nan, np.NINF])
    2.0
    >>> np.nanmax([1, 2, np.nan, np.inf])
    inf

    """
    if not isinstance(a, np.ndarray) or type(a) is np.ndarray:
        # Fast, but not safe for subclasses of ndarray
        res = np.fmax.reduce(a, axis=axis, out=out, keepdims=keepdims)
        if np.isnan(res).any():
            warnings.warn("All-NaN slice encountered", RuntimeWarning)
    else:
        # Slow, but safe for subclasses of ndarray
        a, mask = _replace_nan(a, -np.inf)
        res = np.amax(a, axis=axis, out=out, keepdims=keepdims)
        if mask is None:
            return res

        # Check for all-NaN axis
        mask = np.all(mask, axis=axis, keepdims=keepdims)
        if np.any(mask):
            res = _copyto(res, np.nan, mask)
            warnings.warn("All-NaN axis encountered", RuntimeWarning)
    return res


def nanargmin(a, axis=None):
    """
    Return the indices of the minimum values in the specified axis ignoring
    NaNs. For all-NaN slices ``ValueError`` is raised. Warning: the results
    cannot be trusted if a slice contains only NaNs and Infs.

    Parameters
    ----------
    a : array_like
        Input data.
    axis : int, optional
        Axis along which to operate.  By default flattened input is used.

    Returns
    -------
    index_array : ndarray
        An array of indices or a single index value.

    See Also
    --------
    argmin, nanargmax

    Examples
    --------
    >>> a = np.array([[np.nan, 4], [2, 3]])
    >>> np.argmin(a)
    0
    >>> np.nanargmin(a)
    2
    >>> np.nanargmin(a, axis=0)
    array([1, 1])
    >>> np.nanargmin(a, axis=1)
    array([1, 0])

    """
    a, mask = _replace_nan(a, np.inf)
    res = np.argmin(a, axis=axis)
    if mask is not None:
        mask = np.all(mask, axis=axis)
        if np.any(mask):
            raise ValueError("All-NaN slice encountered")
    return res


def nanargmax(a, axis=None):
    """
    Return the indices of the maximum values in the specified axis ignoring
    NaNs. For all-NaN slices ``ValueError`` is raised. Warning: the
    results cannot be trusted if a slice contains only NaNs and -Infs.


    Parameters
    ----------
    a : array_like
        Input data.
    axis : int, optional
        Axis along which to operate.  By default flattened input is used.

    Returns
    -------
    index_array : ndarray
        An array of indices or a single index value.

    See Also
    --------
    argmax, nanargmin

    Examples
    --------
    >>> a = np.array([[np.nan, 4], [2, 3]])
    >>> np.argmax(a)
    0
    >>> np.nanargmax(a)
    1
    >>> np.nanargmax(a, axis=0)
    array([1, 0])
    >>> np.nanargmax(a, axis=1)
    array([1, 1])

    """
    a, mask = _replace_nan(a, -np.inf)
    res = np.argmax(a, axis=axis)
    if mask is not None:
        mask = np.all(mask, axis=axis)
        if np.any(mask):
            raise ValueError("All-NaN slice encountered")
    return res


def nansum(a, axis=None, dtype=None, out=None, keepdims=0):
    """
    Return the sum of array elements over a given axis treating Not a
    Numbers (NaNs) as zero.

    In Numpy versions <= 1.8 Nan is returned for slices that are all-NaN or
    empty. In later versions zero is returned.

    Parameters
    ----------
    a : array_like
        Array containing numbers whose sum is desired. If `a` is not an
        array, a conversion is attempted.
    axis : int, optional
        Axis along which the sum is computed. The default is to compute the
        sum of the flattened array.
    dtype : data-type, optional
        The type of the returned array and of the accumulator in which the
        elements are summed.  By default, the dtype of `a` is used.  An
        exception is when `a` has an integer type with less precision than
        the platform (u)intp. In that case, the default will be either
        (u)int32 or (u)int64 depending on whether the platform is 32 or 64
        bits. For inexact inputs, dtype must be inexact.

        .. versionadded:: 1.8.0
    out : ndarray, optional
        Alternate output array in which to place the result.  The default
        is ``None``. If provided, it must have the same shape as the
        expected output, but the type will be cast if necessary.  See
        `doc.ufuncs` for details. The casting of NaN to integer can yield
        unexpected results.

        .. versionadded:: 1.8.0
    keepdims : bool, optional
        If True, the axes which are reduced are left in the result as
        dimensions with size one. With this option, the result will
        broadcast correctly against the original `arr`.

        .. versionadded:: 1.8.0

    Returns
    -------
    y : ndarray or numpy scalar

    See Also
    --------
    numpy.sum : Sum across array propagating NaNs.
    isnan : Show which elements are NaN.
    isfinite: Show which elements are not NaN or +/-inf.

    Notes
    -----
    If both positive and negative infinity are present, the sum will be Not
    A Number (NaN).

    Numpy integer arithmetic is modular. If the size of a sum exceeds the
    size of an integer accumulator, its value will wrap around and the
    result will be incorrect. Specifying ``dtype=double`` can alleviate
    that problem.

    Examples
    --------
    >>> np.nansum(1)
    1
    >>> np.nansum([1])
    1
    >>> np.nansum([1, np.nan])
    1.0
    >>> a = np.array([[1, 1], [1, np.nan]])
    >>> np.nansum(a)
    3.0
    >>> np.nansum(a, axis=0)
    array([ 2.,  1.])
    >>> np.nansum([1, np.nan, np.inf])
    inf
    >>> np.nansum([1, np.nan, np.NINF])
    -inf
    >>> np.nansum([1, np.nan, np.inf, -np.inf]) # both +/- infinity present
    nan

    """
    a, mask = _replace_nan(a, 0)
    return np.sum(a, axis=axis, dtype=dtype, out=out, keepdims=keepdims)


def nanmean(a, axis=None, dtype=None, out=None, keepdims=False):
    """
    Compute the arithmetic mean along the specified axis, ignoring NaNs.

    Returns the average of the array elements.  The average is taken over
    the flattened array by default, otherwise over the specified axis.
    `float64` intermediate and return values are used for integer inputs.

    For all-NaN slices, NaN is returned and a `RuntimeWarning` is raised.

    .. versionadded:: 1.8.0

    Parameters
    ----------
    a : array_like
        Array containing numbers whose mean is desired. If `a` is not an
        array, a conversion is attempted.
    axis : int, optional
        Axis along which the means are computed. The default is to compute
        the mean of the flattened array.
    dtype : data-type, optional
        Type to use in computing the mean.  For integer inputs, the default
        is `float64`; for inexact inputs, it is the same as the input
        dtype.
    out : ndarray, optional
        Alternate output array in which to place the result.  The default
        is ``None``; if provided, it must have the same shape as the
        expected output, but the type will be cast if necessary.  See
        `doc.ufuncs` for details.
    keepdims : bool, optional
        If this is set to True, the axes which are reduced are left in the
        result as dimensions with size one. With this option, the result
        will broadcast correctly against the original `arr`.

    Returns
    -------
    m : ndarray, see dtype parameter above
        If `out=None`, returns a new array containing the mean values,
        otherwise a reference to the output array is returned. Nan is
        returned for slices that contain only NaNs.

    See Also
    --------
    average : Weighted average
    mean : Arithmetic mean taken while not ignoring NaNs
    var, nanvar

    Notes
    -----
    The arithmetic mean is the sum of the non-NaN elements along the axis
    divided by the number of non-NaN elements.

    Note that for floating-point input, the mean is computed using the same
    precision the input has.  Depending on the input data, this can cause
    the results to be inaccurate, especially for `float32`.  Specifying a
    higher-precision accumulator using the `dtype` keyword can alleviate
    this issue.

    Examples
    --------
    >>> a = np.array([[1, np.nan], [3, 4]])
    >>> np.nanmean(a)
    2.6666666666666665
    >>> np.nanmean(a, axis=0)
    array([ 2.,  4.])
    >>> np.nanmean(a, axis=1)
    array([ 1.,  3.5])

    """
    arr, mask = _replace_nan(a, 0)
    if mask is None:
        return np.mean(arr, axis=axis, dtype=dtype, out=out, keepdims=keepdims)

    if dtype is not None:
        dtype = np.dtype(dtype)
    if dtype is not None and not issubclass(dtype.type, np.inexact):
        raise TypeError("If a is inexact, then dtype must be inexact")
    if out is not None and not issubclass(out.dtype.type, np.inexact):
        raise TypeError("If a is inexact, then out must be inexact")

    # The warning context speeds things up.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        cnt = np.sum(~mask, axis=axis, dtype=np.intp, keepdims=keepdims)
        tot = np.sum(arr, axis=axis, dtype=dtype, out=out, keepdims=keepdims)
        avg = _divide_by_count(tot, cnt, out=out)

    isbad = (cnt == 0)
    if isbad.any():
        warnings.warn("Mean of empty slice", RuntimeWarning)
        # NaN is the only possible bad value, so no further
        # action is needed to handle bad results.
    return avg


def _nanmedian1d(arr1d, overwrite_input=False):
    """
    Private function for rank 1 arrays. Compute the median ignoring NaNs.
    See nanmedian for parameter usage
    """
    c = np.isnan(arr1d)
    s = np.where(c)[0]
    if s.size == arr1d.size:
        warnings.warn("All-NaN slice encountered", RuntimeWarning)
        return np.nan
    elif s.size == 0:
        return np.median(arr1d, overwrite_input=overwrite_input)
    else:
        if overwrite_input:
            x = arr1d
        else:
            x = arr1d.copy()
        # select non-nans at end of array
        enonan = arr1d[-s.size:][~c[-s.size:]]
        # fill nans in beginning of array with non-nans of end
        x[s[:enonan.size]] = enonan
        # slice nans away
        return np.median(x[:-s.size], overwrite_input=True)


def _nanmedian(a, axis=None, out=None, overwrite_input=False):
    """
    Private function that doesn't support extended axis or keepdims.
    These methods are extended to this function using _ureduce
    See nanmedian for parameter usage

    """
    if axis is None or a.ndim == 1:
        part = a.ravel()
        if out is None:
            return _nanmedian1d(part, overwrite_input)
        else:
            out[...] = _nanmedian1d(part, overwrite_input)
            return out
    else:
        # for small medians use sort + indexing which is still faster than
        # apply_along_axis
        if a.shape[axis] < 400:
            return _nanmedian_small(a, axis, out, overwrite_input)
        result = np.apply_along_axis(_nanmedian1d, axis, a, overwrite_input)
        if out is not None:
            out[...] = result
        return result

def _nanmedian_small(a, axis=None, out=None, overwrite_input=False):
    """
    sort + indexing median, faster for small medians along multiple dimensions
    due to the high overhead of apply_along_axis
    see nanmedian for parameter usage
    """
    a = np.ma.masked_array(a, np.isnan(a))
    m = np.ma.median(a, axis=axis, overwrite_input=overwrite_input)
    for i in range(np.count_nonzero(m.mask.ravel())):
        warnings.warn("All-NaN slice encountered", RuntimeWarning)
    if out is not None:
        out[...] = m.filled(np.nan)
        return out
    return m.filled(np.nan)

def nanmedian(a, axis=None, out=None, overwrite_input=False, keepdims=False):
    """
    Compute the median along the specified axis, while ignoring NaNs.

    Returns the median of the array elements.

    .. versionadded:: 1.9.0

    Parameters
    ----------
    a : array_like
        Input array or object that can be converted to an array.
    axis : int, optional
        Axis along which the medians are computed. The default (axis=None)
        is to compute the median along a flattened version of the array.
        A sequence of axes is supported since version 1.9.0.
    out : ndarray, optional
        Alternative output array in which to place the result. It must have
        the same shape and buffer length as the expected output, but the
        type (of the output) will be cast if necessary.
    overwrite_input : bool, optional
       If True, then allow use of memory of input array (a) for
       calculations. The input array will be modified by the call to
       median. This will save memory when you do not need to preserve
       the contents of the input array. Treat the input as undefined,
       but it will probably be fully or partially sorted. Default is
       False. Note that, if `overwrite_input` is True and the input
       is not already an ndarray, an error will be raised.
    keepdims : bool, optional
        If this is set to True, the axes which are reduced are left
        in the result as dimensions with size one. With this option,
        the result will broadcast correctly against the original `arr`.



    Returns
    -------
    median : ndarray
        A new array holding the result. If the input contains integers, or
        floats of smaller precision than 64, then the output data-type is
        float64.  Otherwise, the output data-type is the same as that of the
        input.

    See Also
    --------
    mean, median, percentile

    Notes
    -----
    Given a vector V of length N, the median of V is the middle value of
    a sorted copy of V, ``V_sorted`` - i.e., ``V_sorted[(N-1)/2]``, when N is
    odd.  When N is even, it is the average of the two middle values of
    ``V_sorted``.

    Examples
    --------
    >>> a = np.array([[10.0, 7, 4], [3, 2, 1]])
    >>> a[0, 1] = np.nan
    >>> a
    array([[ 10.,  nan,   4.],
       [  3.,   2.,   1.]])
    >>> np.median(a)
    nan
    >>> np.nanmedian(a)
    3.0
    >>> np.nanmedian(a, axis=0)
    array([ 6.5,  2.,  2.5])
    >>> np.median(a, axis=1)
    array([ 7.,  2.])
    >>> b = a.copy()
    >>> np.nanmedian(b, axis=1, overwrite_input=True)
    array([ 7.,  2.])
    >>> assert not np.all(a==b)
    >>> b = a.copy()
    >>> np.nanmedian(b, axis=None, overwrite_input=True)
    3.0
    >>> assert not np.all(a==b)

    """
    a = np.asanyarray(a)
    # apply_along_axis in _nanmedian doesn't handle empty arrays well,
    # so deal them upfront
    if a.size == 0:
        return np.nanmean(a, axis, out=out, keepdims=keepdims)

    r, k = _ureduce(a, func=_nanmedian, axis=axis, out=out,
                    overwrite_input=overwrite_input)
    if keepdims:
        return r.reshape(k)
    else:
        return r


def nanpercentile(a, q, axis=None, out=None, overwrite_input=False,
                  interpolation='linear', keepdims=False):
    """
    Compute the qth percentile of the data along the specified axis, while
    ignoring nan values.

    Returns the qth percentile of the array elements.

    Parameters
    ----------
    a : array_like
        Input array or object that can be converted to an array.
    q : float in range of [0,100] (or sequence of floats)
        Percentile to compute which must be between 0 and 100 inclusive.
    axis : int or sequence of int, optional
        Axis along which the percentiles are computed. The default (None)
        is to compute the percentiles along a flattened version of the array.
        A sequence of axes is supported since version 1.9.0.
    out : ndarray, optional
        Alternative output array in which to place the result. It must
        have the same shape and buffer length as the expected output,
        but the type (of the output) will be cast if necessary.
    overwrite_input : bool, optional
        If True, then allow use of memory of input array `a` for
        calculations. The input array will be modified by the call to
        percentile. This will save memory when you do not need to preserve
        the contents of the input array. In this case you should not make
        any assumptions about the content of the passed in array `a` after
        this function completes -- treat it as undefined. Default is False.
        Note that, if the `a` input is not already an array this parameter
        will have no effect, `a` will be converted to an array internally
        regardless of the value of this parameter.
    interpolation : {'linear', 'lower', 'higher', 'midpoint', 'nearest'}
        This optional parameter specifies the interpolation method to use,
        when the desired quantile lies between two data points `i` and `j`:
            * linear: `i + (j - i) * fraction`, where `fraction` is the
              fractional part of the index surrounded by `i` and `j`.
            * lower: `i`.
            * higher: `j`.
            * nearest: `i` or `j` whichever is nearest.
            * midpoint: (`i` + `j`) / 2.

    keepdims : bool, optional
        If this is set to True, the axes which are reduced are left
        in the result as dimensions with size one. With this option,
        the result will broadcast correctly against the original `arr`.


    Returns
    -------
    nanpercentile : scalar or ndarray
        If a single percentile `q` is given and axis=None a scalar is
        returned.  If multiple percentiles `q` are given an array holding
        the result is returned. The results are listed in the first axis.
        (If `out` is specified, in which case that array is returned
        instead).  If the input contains integers, or floats of smaller
        precision than 64, then the output data-type is float64. Otherwise,
        the output data-type is the same as that of the input.

    See Also
    --------
    nanmean, nanmedian, percentile, median, mean

    Notes
    -----
    Given a vector V of length N, the q-th percentile of V is the q-th ranked
    value in a sorted copy of V.  The values and distances of the two
    nearest neighbors as well as the `interpolation` parameter will
    determine the percentile if the normalized ranking does not match q
    exactly. This function is the same as the median if ``q=50``, the same
    as the minimum if ``q=0``and the same as the maximum if ``q=100``.

    Examples
    --------
    >>> a = np.array([[10., 7., 4.], [3., 2., 1.]])
    >>> a[0][1] = np.nan
    >>> a
    array([[ 10.,  nan,   4.],
       [  3.,   2.,   1.]])
    >>> np.percentile(a, 50)
    nan
    >>> np.nanpercentile(a, 50)
    3.5
    >>> np.nanpercentile(a, 50, axis=0)
    array([[ 6.5,  4.5,  2.5]])
    >>> np.nanpercentile(a, 50, axis=1)
    array([[ 7.],
           [ 2.]])
    >>> m = np.nanpercentile(a, 50, axis=0)
    >>> out = np.zeros_like(m)
    >>> np.nanpercentile(a, 50, axis=0, out=m)
    array([[ 6.5,  4.5,  2.5]])
    >>> m
    array([[ 6.5,  4.5,  2.5]])
    >>> b = a.copy()
    >>> np.nanpercentile(b, 50, axis=1, overwrite_input=True)
    array([[ 7.],
           [ 2.]])
    >>> assert not np.all(a==b)
    >>> b = a.copy()
    >>> np.nanpercentile(b, 50, axis=None, overwrite_input=True)
    array([ 3.5])

    """

    a = np.asanyarray(a)
    q = np.asanyarray(q)
    # apply_along_axis in _nanpercentile doesn't handle empty arrays well,
    # so deal them upfront
    if a.size == 0:
        return np.nanmean(a, axis, out=out, keepdims=keepdims)

    r, k = _ureduce(a, func=_nanpercentile, q=q, axis=axis, out=out,
                    overwrite_input=overwrite_input,
                    interpolation=interpolation)
    if keepdims:
        if q.ndim == 0:
            return r.reshape(k)
        else:
            return r.reshape([len(q)] + k)
    else:
        return r


def _nanpercentile(a, q, axis=None, out=None, overwrite_input=False,
                   interpolation='linear', keepdims=False):
    """
    Private function that doesn't support extended axis or keepdims.
    These methods are extended to this function using _ureduce
    See nanpercentile for parameter usage

    """
    if axis is None:
        part = a.ravel()
        result = _nanpercentile1d(part, q, overwrite_input, interpolation)
    else:
        result = np.apply_along_axis(_nanpercentile1d, axis, a, q,
                                     overwrite_input, interpolation)

    if out is not None:
        out[...] = result
    return result


def _nanpercentile1d(arr1d, q, overwrite_input=False, interpolation='linear'):
    """
    Private function for rank 1 arrays. Compute percentile ignoring NaNs.
    See nanpercentile for parameter usage

    """
    c = np.isnan(arr1d)
    s = np.where(c)[0]
    if s.size == arr1d.size:
        warnings.warn("All-NaN slice encountered", RuntimeWarning)
        return np.nan
    elif s.size == 0:
        return np.percentile(arr1d, q, overwrite_input=overwrite_input,
                             interpolation=interpolation)
    else:
        if overwrite_input:
            x = arr1d
        else:
            x = arr1d.copy()
        # select non-nans at end of array
        enonan = arr1d[-s.size:][~c[-s.size:]]
        # fill nans in beginning of array with non-nans of end
        x[s[:enonan.size]] = enonan
        # slice nans away
        return np.percentile(x[:-s.size], q, overwrite_input=True,
                             interpolation=interpolation)


def nanvar(a, axis=None, dtype=None, out=None, ddof=0, keepdims=False):
    """
    Compute the variance along the specified axis, while ignoring NaNs.

    Returns the variance of the array elements, a measure of the spread of
    a distribution.  The variance is computed for the flattened array by
    default, otherwise over the specified axis.

    For all-NaN slices or slices with zero degrees of freedom, NaN is
    returned and a `RuntimeWarning` is raised.

    .. versionadded:: 1.8.0

    Parameters
    ----------
    a : array_like
        Array containing numbers whose variance is desired.  If `a` is not an
        array, a conversion is attempted.
    axis : int, optional
        Axis along which the variance is computed.  The default is to compute
        the variance of the flattened array.
    dtype : data-type, optional
        Type to use in computing the variance.  For arrays of integer type
        the default is `float32`; for arrays of float types it is the same as
        the array type.
    out : ndarray, optional
        Alternate output array in which to place the result.  It must have
        the same shape as the expected output, but the type is cast if
        necessary.
    ddof : int, optional
        "Delta Degrees of Freedom": the divisor used in the calculation is
        ``N - ddof``, where ``N`` represents the number of non-NaN
        elements. By default `ddof` is zero.
    keepdims : bool, optional
        If this is set to True, the axes which are reduced are left
        in the result as dimensions with size one. With this option,
        the result will broadcast correctly against the original `arr`.

    Returns
    -------
    variance : ndarray, see dtype parameter above
        If `out` is None, return a new array containing the variance,
        otherwise return a reference to the output array. If ddof is >= the
        number of non-NaN elements in a slice or the slice contains only
        NaNs, then the result for that slice is NaN.

    See Also
    --------
    std : Standard deviation
    mean : Average
    var : Variance while not ignoring NaNs
    nanstd, nanmean
    numpy.doc.ufuncs : Section "Output arguments"

    Notes
    -----
    The variance is the average of the squared deviations from the mean,
    i.e.,  ``var = mean(abs(x - x.mean())**2)``.

    The mean is normally calculated as ``x.sum() / N``, where ``N = len(x)``.
    If, however, `ddof` is specified, the divisor ``N - ddof`` is used
    instead.  In standard statistical practice, ``ddof=1`` provides an
    unbiased estimator of the variance of a hypothetical infinite
    population.  ``ddof=0`` provides a maximum likelihood estimate of the
    variance for normally distributed variables.

    Note that for complex numbers, the absolute value is taken before
    squaring, so that the result is always real and nonnegative.

    For floating-point input, the variance is computed using the same
    precision the input has.  Depending on the input data, this can cause
    the results to be inaccurate, especially for `float32` (see example
    below).  Specifying a higher-accuracy accumulator using the ``dtype``
    keyword can alleviate this issue.

    Examples
    --------
    >>> a = np.array([[1, np.nan], [3, 4]])
    >>> np.var(a)
    1.5555555555555554
    >>> np.nanvar(a, axis=0)
    array([ 1.,  0.])
    >>> np.nanvar(a, axis=1)
    array([ 0.,  0.25])

    """
    arr, mask = _replace_nan(a, 0)
    if mask is None:
        return np.var(arr, axis=axis, dtype=dtype, out=out, ddof=ddof,
                      keepdims=keepdims)

    if dtype is not None:
        dtype = np.dtype(dtype)
    if dtype is not None and not issubclass(dtype.type, np.inexact):
        raise TypeError("If a is inexact, then dtype must be inexact")
    if out is not None and not issubclass(out.dtype.type, np.inexact):
        raise TypeError("If a is inexact, then out must be inexact")

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')

        # Compute mean
        cnt = np.sum(~mask, axis=axis, dtype=np.intp, keepdims=True)
        avg = np.sum(arr, axis=axis, dtype=dtype, keepdims=True)
        avg = _divide_by_count(avg, cnt)

        # Compute squared deviation from mean.
        arr -= avg
        arr = _copyto(arr, 0, mask)
        if issubclass(arr.dtype.type, np.complexfloating):
            sqr = np.multiply(arr, arr.conj(), out=arr).real
        else:
            sqr = np.multiply(arr, arr, out=arr)

        # Compute variance.
        var = np.sum(sqr, axis=axis, dtype=dtype, out=out, keepdims=keepdims)
        if var.ndim < cnt.ndim:
            # Subclasses of ndarray may ignore keepdims, so check here.
            cnt = cnt.squeeze(axis)
        dof = cnt - ddof
        var = _divide_by_count(var, dof)

    isbad = (dof <= 0)
    if np.any(isbad):
        warnings.warn("Degrees of freedom <= 0 for slice.", RuntimeWarning)
        # NaN, inf, or negative numbers are all possible bad
        # values, so explicitly replace them with NaN.
        var = _copyto(var, np.nan, isbad)
    return var


def nanstd(a, axis=None, dtype=None, out=None, ddof=0, keepdims=False):
    """
    Compute the standard deviation along the specified axis, while
    ignoring NaNs.

    Returns the standard deviation, a measure of the spread of a
    distribution, of the non-NaN array elements. The standard deviation is
    computed for the flattened array by default, otherwise over the
    specified axis.

    For all-NaN slices or slices with zero degrees of freedom, NaN is
    returned and a `RuntimeWarning` is raised.

    .. versionadded:: 1.8.0

    Parameters
    ----------
    a : array_like
        Calculate the standard deviation of the non-NaN values.
    axis : int, optional
        Axis along which the standard deviation is computed. The default is
        to compute the standard deviation of the flattened array.
    dtype : dtype, optional
        Type to use in computing the standard deviation. For arrays of
        integer type the default is float64, for arrays of float types it
        is the same as the array type.
    out : ndarray, optional
        Alternative output array in which to place the result. It must have
        the same shape as the expected output but the type (of the
        calculated values) will be cast if necessary.
    ddof : int, optional
        Means Delta Degrees of Freedom.  The divisor used in calculations
        is ``N - ddof``, where ``N`` represents the number of non-NaN
        elements.  By default `ddof` is zero.
    keepdims : bool, optional
        If this is set to True, the axes which are reduced are left
        in the result as dimensions with size one. With this option,
        the result will broadcast correctly against the original `arr`.

    Returns
    -------
    standard_deviation : ndarray, see dtype parameter above.
        If `out` is None, return a new array containing the standard
        deviation, otherwise return a reference to the output array. If
        ddof is >= the number of non-NaN elements in a slice or the slice
        contains only NaNs, then the result for that slice is NaN.

    See Also
    --------
    var, mean, std
    nanvar, nanmean
    numpy.doc.ufuncs : Section "Output arguments"

    Notes
    -----
    The standard deviation is the square root of the average of the squared
    deviations from the mean: ``std = sqrt(mean(abs(x - x.mean())**2))``.

    The average squared deviation is normally calculated as
    ``x.sum() / N``, where ``N = len(x)``.  If, however, `ddof` is
    specified, the divisor ``N - ddof`` is used instead. In standard
    statistical practice, ``ddof=1`` provides an unbiased estimator of the
    variance of the infinite population. ``ddof=0`` provides a maximum
    likelihood estimate of the variance for normally distributed variables.
    The standard deviation computed in this function is the square root of
    the estimated variance, so even with ``ddof=1``, it will not be an
    unbiased estimate of the standard deviation per se.

    Note that, for complex numbers, `std` takes the absolute value before
    squaring, so that the result is always real and nonnegative.

    For floating-point input, the *std* is computed using the same
    precision the input has. Depending on the input data, this can cause
    the results to be inaccurate, especially for float32 (see example
    below).  Specifying a higher-accuracy accumulator using the `dtype`
    keyword can alleviate this issue.

    Examples
    --------
    >>> a = np.array([[1, np.nan], [3, 4]])
    >>> np.nanstd(a)
    1.247219128924647
    >>> np.nanstd(a, axis=0)
    array([ 1.,  0.])
    >>> np.nanstd(a, axis=1)
    array([ 0.,  0.5])

    """
    var = nanvar(a, axis=axis, dtype=dtype, out=out, ddof=ddof,
                 keepdims=keepdims)
    if isinstance(var, np.ndarray):
        std = np.sqrt(var, out=var)
    else:
        std = var.dtype.type(np.sqrt(var))
    return std
