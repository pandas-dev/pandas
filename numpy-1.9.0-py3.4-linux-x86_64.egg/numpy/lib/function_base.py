from __future__ import division, absolute_import, print_function

import warnings
import sys
import collections
import operator

import numpy as np
import numpy.core.numeric as _nx
from numpy.core import linspace, atleast_1d, atleast_2d
from numpy.core.numeric import (
    ones, zeros, arange, concatenate, array, asarray, asanyarray, empty,
    empty_like, ndarray, around, floor, ceil, take, dot, where, intp,
    integer, isscalar
    )
from numpy.core.umath import (
    pi, multiply, add, arctan2, frompyfunc, cos, less_equal, sqrt, sin,
    mod, exp, log10
    )
from numpy.core.fromnumeric import (
    ravel, nonzero, sort, partition, mean
    )
from numpy.core.numerictypes import typecodes, number
from numpy.lib.twodim_base import diag
from .utils import deprecate
from ._compiled_base import _insert, add_docstring
from ._compiled_base import digitize, bincount, interp as compiled_interp
from ._compiled_base import add_newdoc_ufunc
from numpy.compat import long

# Force range to be a generator, for np.delete's usage.
if sys.version_info[0] < 3:
    range = xrange


__all__ = [
    'select', 'piecewise', 'trim_zeros', 'copy', 'iterable', 'percentile',
    'diff', 'gradient', 'angle', 'unwrap', 'sort_complex', 'disp',
    'extract', 'place', 'vectorize', 'asarray_chkfinite', 'average',
    'histogram', 'histogramdd', 'bincount', 'digitize', 'cov', 'corrcoef',
    'msort', 'median', 'sinc', 'hamming', 'hanning', 'bartlett',
    'blackman', 'kaiser', 'trapz', 'i0', 'add_newdoc', 'add_docstring',
    'meshgrid', 'delete', 'insert', 'append', 'interp', 'add_newdoc_ufunc'
    ]


def iterable(y):
    """
    Check whether or not an object can be iterated over.

    Parameters
    ----------
    y : object
      Input object.

    Returns
    -------
    b : {0, 1}
      Return 1 if the object has an iterator method or is a sequence,
      and 0 otherwise.


    Examples
    --------
    >>> np.iterable([1, 2, 3])
    1
    >>> np.iterable(2)
    0

    """
    try:
        iter(y)
    except:
        return 0
    return 1


def histogram(a, bins=10, range=None, normed=False, weights=None,
              density=None):
    """
    Compute the histogram of a set of data.

    Parameters
    ----------
    a : array_like
        Input data. The histogram is computed over the flattened array.
    bins : int or sequence of scalars, optional
        If `bins` is an int, it defines the number of equal-width
        bins in the given range (10, by default). If `bins` is a sequence,
        it defines the bin edges, including the rightmost edge, allowing
        for non-uniform bin widths.
    range : (float, float), optional
        The lower and upper range of the bins.  If not provided, range
        is simply ``(a.min(), a.max())``.  Values outside the range are
        ignored.
    normed : bool, optional
        This keyword is deprecated in Numpy 1.6 due to confusing/buggy
        behavior. It will be removed in Numpy 2.0. Use the density keyword
        instead.
        If False, the result will contain the number of samples
        in each bin.  If True, the result is the value of the
        probability *density* function at the bin, normalized such that
        the *integral* over the range is 1. Note that this latter behavior is
        known to be buggy with unequal bin widths; use `density` instead.
    weights : array_like, optional
        An array of weights, of the same shape as `a`.  Each value in `a`
        only contributes its associated weight towards the bin count
        (instead of 1).  If `normed` is True, the weights are normalized,
        so that the integral of the density over the range remains 1
    density : bool, optional
        If False, the result will contain the number of samples
        in each bin.  If True, the result is the value of the
        probability *density* function at the bin, normalized such that
        the *integral* over the range is 1. Note that the sum of the
        histogram values will not be equal to 1 unless bins of unity
        width are chosen; it is not a probability *mass* function.
        Overrides the `normed` keyword if given.

    Returns
    -------
    hist : array
        The values of the histogram. See `normed` and `weights` for a
        description of the possible semantics.
    bin_edges : array of dtype float
        Return the bin edges ``(length(hist)+1)``.


    See Also
    --------
    histogramdd, bincount, searchsorted, digitize

    Notes
    -----
    All but the last (righthand-most) bin is half-open.  In other words, if
    `bins` is::

      [1, 2, 3, 4]

    then the first bin is ``[1, 2)`` (including 1, but excluding 2) and the
    second ``[2, 3)``.  The last bin, however, is ``[3, 4]``, which *includes*
    4.

    Examples
    --------
    >>> np.histogram([1, 2, 1], bins=[0, 1, 2, 3])
    (array([0, 2, 1]), array([0, 1, 2, 3]))
    >>> np.histogram(np.arange(4), bins=np.arange(5), density=True)
    (array([ 0.25,  0.25,  0.25,  0.25]), array([0, 1, 2, 3, 4]))
    >>> np.histogram([[1, 2, 1], [1, 0, 1]], bins=[0,1,2,3])
    (array([1, 4, 1]), array([0, 1, 2, 3]))

    >>> a = np.arange(5)
    >>> hist, bin_edges = np.histogram(a, density=True)
    >>> hist
    array([ 0.5,  0. ,  0.5,  0. ,  0. ,  0.5,  0. ,  0.5,  0. ,  0.5])
    >>> hist.sum()
    2.4999999999999996
    >>> np.sum(hist*np.diff(bin_edges))
    1.0

    """

    a = asarray(a)
    if weights is not None:
        weights = asarray(weights)
        if np.any(weights.shape != a.shape):
            raise ValueError(
                'weights should have the same shape as a.')
        weights = weights.ravel()
    a = a.ravel()

    if (range is not None):
        mn, mx = range
        if (mn > mx):
            raise AttributeError(
                'max must be larger than min in range parameter.')

    if not iterable(bins):
        if np.isscalar(bins) and bins < 1:
            raise ValueError(
                '`bins` should be a positive integer.')
        if range is None:
            if a.size == 0:
                # handle empty arrays. Can't determine range, so use 0-1.
                range = (0, 1)
            else:
                range = (a.min(), a.max())
        mn, mx = [mi + 0.0 for mi in range]
        if mn == mx:
            mn -= 0.5
            mx += 0.5
        bins = linspace(mn, mx, bins + 1, endpoint=True)
    else:
        bins = asarray(bins)
        if (np.diff(bins) < 0).any():
            raise AttributeError(
                'bins must increase monotonically.')

    # Histogram is an integer or a float array depending on the weights.
    if weights is None:
        ntype = int
    else:
        ntype = weights.dtype
    n = np.zeros(bins.shape, ntype)

    block = 65536
    if weights is None:
        for i in arange(0, len(a), block):
            sa = sort(a[i:i+block])
            n += np.r_[sa.searchsorted(bins[:-1], 'left'),
                       sa.searchsorted(bins[-1], 'right')]
    else:
        zero = array(0, dtype=ntype)
        for i in arange(0, len(a), block):
            tmp_a = a[i:i+block]
            tmp_w = weights[i:i+block]
            sorting_index = np.argsort(tmp_a)
            sa = tmp_a[sorting_index]
            sw = tmp_w[sorting_index]
            cw = np.concatenate(([zero, ], sw.cumsum()))
            bin_index = np.r_[sa.searchsorted(bins[:-1], 'left'),
                              sa.searchsorted(bins[-1], 'right')]
            n += cw[bin_index]

    n = np.diff(n)

    if density is not None:
        if density:
            db = array(np.diff(bins), float)
            return n/db/n.sum(), bins
        else:
            return n, bins
    else:
        # deprecated, buggy behavior. Remove for Numpy 2.0
        if normed:
            db = array(np.diff(bins), float)
            return n/(n*db).sum(), bins
        else:
            return n, bins


def histogramdd(sample, bins=10, range=None, normed=False, weights=None):
    """
    Compute the multidimensional histogram of some data.

    Parameters
    ----------
    sample : array_like
        The data to be histogrammed. It must be an (N,D) array or data
        that can be converted to such. The rows of the resulting array
        are the coordinates of points in a D dimensional polytope.
    bins : sequence or int, optional
        The bin specification:

        * A sequence of arrays describing the bin edges along each dimension.
        * The number of bins for each dimension (nx, ny, ... =bins)
        * The number of bins for all dimensions (nx=ny=...=bins).

    range : sequence, optional
        A sequence of lower and upper bin edges to be used if the edges are
        not given explicitly in `bins`. Defaults to the minimum and maximum
        values along each dimension.
    normed : bool, optional
        If False, returns the number of samples in each bin. If True,
        returns the bin density ``bin_count / sample_count / bin_volume``.
    weights : array_like (N,), optional
        An array of values `w_i` weighing each sample `(x_i, y_i, z_i, ...)`.
        Weights are normalized to 1 if normed is True. If normed is False,
        the values of the returned histogram are equal to the sum of the
        weights belonging to the samples falling into each bin.

    Returns
    -------
    H : ndarray
        The multidimensional histogram of sample x. See normed and weights
        for the different possible semantics.
    edges : list
        A list of D arrays describing the bin edges for each dimension.

    See Also
    --------
    histogram: 1-D histogram
    histogram2d: 2-D histogram

    Examples
    --------
    >>> r = np.random.randn(100,3)
    >>> H, edges = np.histogramdd(r, bins = (5, 8, 4))
    >>> H.shape, edges[0].size, edges[1].size, edges[2].size
    ((5, 8, 4), 6, 9, 5)

    """

    try:
        # Sample is an ND-array.
        N, D = sample.shape
    except (AttributeError, ValueError):
        # Sample is a sequence of 1D arrays.
        sample = atleast_2d(sample).T
        N, D = sample.shape

    nbin = empty(D, int)
    edges = D*[None]
    dedges = D*[None]
    if weights is not None:
        weights = asarray(weights)

    try:
        M = len(bins)
        if M != D:
            raise AttributeError(
                'The dimension of bins must be equal to the dimension of the '
                ' sample x.')
    except TypeError:
        # bins is an integer
        bins = D*[bins]

    # Select range for each dimension
    # Used only if number of bins is given.
    if range is None:
        # Handle empty input. Range can't be determined in that case, use 0-1.
        if N == 0:
            smin = zeros(D)
            smax = ones(D)
        else:
            smin = atleast_1d(array(sample.min(0), float))
            smax = atleast_1d(array(sample.max(0), float))
    else:
        smin = zeros(D)
        smax = zeros(D)
        for i in arange(D):
            smin[i], smax[i] = range[i]

    # Make sure the bins have a finite width.
    for i in arange(len(smin)):
        if smin[i] == smax[i]:
            smin[i] = smin[i] - .5
            smax[i] = smax[i] + .5

    # avoid rounding issues for comparisons when dealing with inexact types
    if np.issubdtype(sample.dtype, np.inexact):
        edge_dt = sample.dtype
    else:
        edge_dt = float
    # Create edge arrays
    for i in arange(D):
        if isscalar(bins[i]):
            if bins[i] < 1:
                raise ValueError(
                    "Element at index %s in `bins` should be a positive "
                    "integer." % i)
            nbin[i] = bins[i] + 2  # +2 for outlier bins
            edges[i] = linspace(smin[i], smax[i], nbin[i]-1, dtype=edge_dt)
        else:
            edges[i] = asarray(bins[i], edge_dt)
            nbin[i] = len(edges[i]) + 1  # +1 for outlier bins
        dedges[i] = diff(edges[i])
        if np.any(np.asarray(dedges[i]) <= 0):
            raise ValueError(
                "Found bin edge of size <= 0. Did you specify `bins` with"
                "non-monotonic sequence?")

    nbin = asarray(nbin)

    # Handle empty input.
    if N == 0:
        return np.zeros(nbin-2), edges

    # Compute the bin number each sample falls into.
    Ncount = {}
    for i in arange(D):
        Ncount[i] = digitize(sample[:, i], edges[i])

    # Using digitize, values that fall on an edge are put in the right bin.
    # For the rightmost bin, we want values equal to the right edge to be
    # counted in the last bin, and not as an outlier.
    for i in arange(D):
        # Rounding precision
        mindiff = dedges[i].min()
        if not np.isinf(mindiff):
            decimal = int(-log10(mindiff)) + 6
            # Find which points are on the rightmost edge.
            not_smaller_than_edge = (sample[:, i] >= edges[i][-1])
            on_edge = (around(sample[:, i], decimal) ==
                       around(edges[i][-1], decimal))
            # Shift these points one bin to the left.
            Ncount[i][where(on_edge & not_smaller_than_edge)[0]] -= 1

    # Flattened histogram matrix (1D)
    # Reshape is used so that overlarge arrays
    # will raise an error.
    hist = zeros(nbin, float).reshape(-1)

    # Compute the sample indices in the flattened histogram matrix.
    ni = nbin.argsort()
    xy = zeros(N, int)
    for i in arange(0, D-1):
        xy += Ncount[ni[i]] * nbin[ni[i+1:]].prod()
    xy += Ncount[ni[-1]]

    # Compute the number of repetitions in xy and assign it to the
    # flattened histmat.
    if len(xy) == 0:
        return zeros(nbin-2, int), edges

    flatcount = bincount(xy, weights)
    a = arange(len(flatcount))
    hist[a] = flatcount

    # Shape into a proper matrix
    hist = hist.reshape(sort(nbin))
    for i in arange(nbin.size):
        j = ni.argsort()[i]
        hist = hist.swapaxes(i, j)
        ni[i], ni[j] = ni[j], ni[i]

    # Remove outliers (indices 0 and -1 for each dimension).
    core = D*[slice(1, -1)]
    hist = hist[core]

    # Normalize if normed is True
    if normed:
        s = hist.sum()
        for i in arange(D):
            shape = ones(D, int)
            shape[i] = nbin[i] - 2
            hist = hist / dedges[i].reshape(shape)
        hist /= s

    if (hist.shape != nbin - 2).any():
        raise RuntimeError(
            "Internal Shape Error")
    return hist, edges


def average(a, axis=None, weights=None, returned=False):
    """
    Compute the weighted average along the specified axis.

    Parameters
    ----------
    a : array_like
        Array containing data to be averaged. If `a` is not an array, a
        conversion is attempted.
    axis : int, optional
        Axis along which to average `a`. If `None`, averaging is done over
        the flattened array.
    weights : array_like, optional
        An array of weights associated with the values in `a`. Each value in
        `a` contributes to the average according to its associated weight.
        The weights array can either be 1-D (in which case its length must be
        the size of `a` along the given axis) or of the same shape as `a`.
        If `weights=None`, then all data in `a` are assumed to have a
        weight equal to one.
    returned : bool, optional
        Default is `False`. If `True`, the tuple (`average`, `sum_of_weights`)
        is returned, otherwise only the average is returned.
        If `weights=None`, `sum_of_weights` is equivalent to the number of
        elements over which the average is taken.


    Returns
    -------
    average, [sum_of_weights] : {array_type, double}
        Return the average along the specified axis. When returned is `True`,
        return a tuple with the average as the first element and the sum
        of the weights as the second element. The return type is `Float`
        if `a` is of integer type, otherwise it is of the same type as `a`.
        `sum_of_weights` is of the same type as `average`.

    Raises
    ------
    ZeroDivisionError
        When all weights along axis are zero. See `numpy.ma.average` for a
        version robust to this type of error.
    TypeError
        When the length of 1D `weights` is not the same as the shape of `a`
        along axis.

    See Also
    --------
    mean

    ma.average : average for masked arrays -- useful if your data contains
                 "missing" values

    Examples
    --------
    >>> data = range(1,5)
    >>> data
    [1, 2, 3, 4]
    >>> np.average(data)
    2.5
    >>> np.average(range(1,11), weights=range(10,0,-1))
    4.0

    >>> data = np.arange(6).reshape((3,2))
    >>> data
    array([[0, 1],
           [2, 3],
           [4, 5]])
    >>> np.average(data, axis=1, weights=[1./4, 3./4])
    array([ 0.75,  2.75,  4.75])
    >>> np.average(data, weights=[1./4, 3./4])
    Traceback (most recent call last):
    ...
    TypeError: Axis must be specified when shapes of a and weights differ.

    """
    if not isinstance(a, np.matrix):
        a = np.asarray(a)

    if weights is None:
        avg = a.mean(axis)
        scl = avg.dtype.type(a.size/avg.size)
    else:
        a = a + 0.0
        wgt = np.array(weights, dtype=a.dtype, copy=0)

        # Sanity checks
        if a.shape != wgt.shape:
            if axis is None:
                raise TypeError(
                    "Axis must be specified when shapes of a and weights "
                    "differ.")
            if wgt.ndim != 1:
                raise TypeError(
                    "1D weights expected when shapes of a and weights differ.")
            if wgt.shape[0] != a.shape[axis]:
                raise ValueError(
                    "Length of weights not compatible with specified axis.")

            # setup wgt to broadcast along axis
            wgt = np.array(wgt, copy=0, ndmin=a.ndim).swapaxes(-1, axis)

        scl = wgt.sum(axis=axis)
        if (scl == 0.0).any():
            raise ZeroDivisionError(
                "Weights sum to zero, can't be normalized")

        avg = np.multiply(a, wgt).sum(axis)/scl

    if returned:
        scl = np.multiply(avg, 0) + scl
        return avg, scl
    else:
        return avg


def asarray_chkfinite(a, dtype=None, order=None):
    """
    Convert the input to an array, checking for NaNs or Infs.

    Parameters
    ----------
    a : array_like
        Input data, in any form that can be converted to an array.  This
        includes lists, lists of tuples, tuples, tuples of tuples, tuples
        of lists and ndarrays.  Success requires no NaNs or Infs.
    dtype : data-type, optional
        By default, the data-type is inferred from the input data.
    order : {'C', 'F'}, optional
        Whether to use row-major ('C') or column-major ('FORTRAN') memory
        representation.  Defaults to 'C'.

    Returns
    -------
    out : ndarray
        Array interpretation of `a`.  No copy is performed if the input
        is already an ndarray.  If `a` is a subclass of ndarray, a base
        class ndarray is returned.

    Raises
    ------
    ValueError
        Raises ValueError if `a` contains NaN (Not a Number) or Inf (Infinity).

    See Also
    --------
    asarray : Create and array.
    asanyarray : Similar function which passes through subclasses.
    ascontiguousarray : Convert input to a contiguous array.
    asfarray : Convert input to a floating point ndarray.
    asfortranarray : Convert input to an ndarray with column-major
                     memory order.
    fromiter : Create an array from an iterator.
    fromfunction : Construct an array by executing a function on grid
                   positions.

    Examples
    --------
    Convert a list into an array.  If all elements are finite
    ``asarray_chkfinite`` is identical to ``asarray``.

    >>> a = [1, 2]
    >>> np.asarray_chkfinite(a, dtype=float)
    array([1., 2.])

    Raises ValueError if array_like contains Nans or Infs.

    >>> a = [1, 2, np.inf]
    >>> try:
    ...     np.asarray_chkfinite(a)
    ... except ValueError:
    ...     print 'ValueError'
    ...
    ValueError

    """
    a = asarray(a, dtype=dtype, order=order)
    if a.dtype.char in typecodes['AllFloat'] and not np.isfinite(a).all():
        raise ValueError(
            "array must not contain infs or NaNs")
    return a


def piecewise(x, condlist, funclist, *args, **kw):
    """
    Evaluate a piecewise-defined function.

    Given a set of conditions and corresponding functions, evaluate each
    function on the input data wherever its condition is true.

    Parameters
    ----------
    x : ndarray
        The input domain.
    condlist : list of bool arrays
        Each boolean array corresponds to a function in `funclist`.  Wherever
        `condlist[i]` is True, `funclist[i](x)` is used as the output value.

        Each boolean array in `condlist` selects a piece of `x`,
        and should therefore be of the same shape as `x`.

        The length of `condlist` must correspond to that of `funclist`.
        If one extra function is given, i.e. if
        ``len(funclist) - len(condlist) == 1``, then that extra function
        is the default value, used wherever all conditions are false.
    funclist : list of callables, f(x,*args,**kw), or scalars
        Each function is evaluated over `x` wherever its corresponding
        condition is True.  It should take an array as input and give an array
        or a scalar value as output.  If, instead of a callable,
        a scalar is provided then a constant function (``lambda x: scalar``) is
        assumed.
    args : tuple, optional
        Any further arguments given to `piecewise` are passed to the functions
        upon execution, i.e., if called ``piecewise(..., ..., 1, 'a')``, then
        each function is called as ``f(x, 1, 'a')``.
    kw : dict, optional
        Keyword arguments used in calling `piecewise` are passed to the
        functions upon execution, i.e., if called
        ``piecewise(..., ..., lambda=1)``, then each function is called as
        ``f(x, lambda=1)``.

    Returns
    -------
    out : ndarray
        The output is the same shape and type as x and is found by
        calling the functions in `funclist` on the appropriate portions of `x`,
        as defined by the boolean arrays in `condlist`.  Portions not covered
        by any condition have a default value of 0.


    See Also
    --------
    choose, select, where

    Notes
    -----
    This is similar to choose or select, except that functions are
    evaluated on elements of `x` that satisfy the corresponding condition from
    `condlist`.

    The result is::

            |--
            |funclist[0](x[condlist[0]])
      out = |funclist[1](x[condlist[1]])
            |...
            |funclist[n2](x[condlist[n2]])
            |--

    Examples
    --------
    Define the sigma function, which is -1 for ``x < 0`` and +1 for ``x >= 0``.

    >>> x = np.linspace(-2.5, 2.5, 6)
    >>> np.piecewise(x, [x < 0, x >= 0], [-1, 1])
    array([-1., -1., -1.,  1.,  1.,  1.])

    Define the absolute value, which is ``-x`` for ``x <0`` and ``x`` for
    ``x >= 0``.

    >>> np.piecewise(x, [x < 0, x >= 0], [lambda x: -x, lambda x: x])
    array([ 2.5,  1.5,  0.5,  0.5,  1.5,  2.5])

    """
    x = asanyarray(x)
    n2 = len(funclist)
    if (isscalar(condlist) or not (isinstance(condlist[0], list) or
                                   isinstance(condlist[0], ndarray))):
        condlist = [condlist]
    condlist = array(condlist, dtype=bool)
    n = len(condlist)
    # This is a hack to work around problems with NumPy's
    #  handling of 0-d arrays and boolean indexing with
    #  numpy.bool_ scalars
    zerod = False
    if x.ndim == 0:
        x = x[None]
        zerod = True
        if condlist.shape[-1] != 1:
            condlist = condlist.T
    if n == n2 - 1:  # compute the "otherwise" condition.
        totlist = np.logical_or.reduce(condlist, axis=0)
        condlist = np.vstack([condlist, ~totlist])
        n += 1
    if (n != n2):
        raise ValueError(
                "function list and condition list must be the same")

    y = zeros(x.shape, x.dtype)
    for k in range(n):
        item = funclist[k]
        if not isinstance(item, collections.Callable):
            y[condlist[k]] = item
        else:
            vals = x[condlist[k]]
            if vals.size > 0:
                y[condlist[k]] = item(vals, *args, **kw)
    if zerod:
        y = y.squeeze()
    return y


def select(condlist, choicelist, default=0):
    """
    Return an array drawn from elements in choicelist, depending on conditions.

    Parameters
    ----------
    condlist : list of bool ndarrays
        The list of conditions which determine from which array in `choicelist`
        the output elements are taken. When multiple conditions are satisfied,
        the first one encountered in `condlist` is used.
    choicelist : list of ndarrays
        The list of arrays from which the output elements are taken. It has
        to be of the same length as `condlist`.
    default : scalar, optional
        The element inserted in `output` when all conditions evaluate to False.

    Returns
    -------
    output : ndarray
        The output at position m is the m-th element of the array in
        `choicelist` where the m-th element of the corresponding array in
        `condlist` is True.

    See Also
    --------
    where : Return elements from one of two arrays depending on condition.
    take, choose, compress, diag, diagonal

    Examples
    --------
    >>> x = np.arange(10)
    >>> condlist = [x<3, x>5]
    >>> choicelist = [x, x**2]
    >>> np.select(condlist, choicelist)
    array([ 0,  1,  2,  0,  0,  0, 36, 49, 64, 81])

    """
    # Check the size of condlist and choicelist are the same, or abort.
    if len(condlist) != len(choicelist):
        raise ValueError(
            'list of cases must be same length as list of conditions')

    # Now that the dtype is known, handle the deprecated select([], []) case
    if len(condlist) == 0:
        warnings.warn("select with an empty condition list is not possible"
                      "and will be deprecated",
                      DeprecationWarning)
        return np.asarray(default)[()]

    choicelist = [np.asarray(choice) for choice in choicelist]
    choicelist.append(np.asarray(default))

    # need to get the result type before broadcasting for correct scalar
    # behaviour
    dtype = np.result_type(*choicelist)

    # Convert conditions to arrays and broadcast conditions and choices
    # as the shape is needed for the result. Doing it seperatly optimizes
    # for example when all choices are scalars.
    condlist = np.broadcast_arrays(*condlist)
    choicelist = np.broadcast_arrays(*choicelist)

    # If cond array is not an ndarray in boolean format or scalar bool, abort.
    deprecated_ints = False
    for i in range(len(condlist)):
        cond = condlist[i]
        if cond.dtype.type is not np.bool_:
            if np.issubdtype(cond.dtype, np.integer):
                # A previous implementation accepted int ndarrays accidentally.
                # Supported here deliberately, but deprecated.
                condlist[i] = condlist[i].astype(bool)
                deprecated_ints = True
            else:
                raise ValueError(
                    'invalid entry in choicelist: should be boolean ndarray')

    if deprecated_ints:
        msg = "select condlists containing integer ndarrays is deprecated " \
            "and will be removed in the future. Use `.astype(bool)` to " \
            "convert to bools."
        warnings.warn(msg, DeprecationWarning)

    if choicelist[0].ndim == 0:
        # This may be common, so avoid the call.
        result_shape = condlist[0].shape
    else:
        result_shape = np.broadcast_arrays(condlist[0], choicelist[0])[0].shape

    result = np.full(result_shape, choicelist[-1], dtype)

    # Use np.copyto to burn each choicelist array onto result, using the
    # corresponding condlist as a boolean mask. This is done in reverse
    # order since the first choice should take precedence.
    choicelist = choicelist[-2::-1]
    condlist = condlist[::-1]
    for choice, cond in zip(choicelist, condlist):
        np.copyto(result, choice, where=cond)

    return result


def copy(a, order='K'):
    """
    Return an array copy of the given object.

    Parameters
    ----------
    a : array_like
        Input data.
    order : {'C', 'F', 'A', 'K'}, optional
        Controls the memory layout of the copy. 'C' means C-order,
        'F' means F-order, 'A' means 'F' if `a` is Fortran contiguous,
        'C' otherwise. 'K' means match the layout of `a` as closely
        as possible. (Note that this function and :meth:ndarray.copy are very
        similar, but have different default values for their order=
        arguments.)

    Returns
    -------
    arr : ndarray
        Array interpretation of `a`.

    Notes
    -----
    This is equivalent to

    >>> np.array(a, copy=True)                              #doctest: +SKIP

    Examples
    --------
    Create an array x, with a reference y and a copy z:

    >>> x = np.array([1, 2, 3])
    >>> y = x
    >>> z = np.copy(x)

    Note that, when we modify x, y changes, but not z:

    >>> x[0] = 10
    >>> x[0] == y[0]
    True
    >>> x[0] == z[0]
    False

    """
    return array(a, order=order, copy=True)

# Basic operations


def gradient(f, *varargs):
    """
    Return the gradient of an N-dimensional array.

    The gradient is computed using second order accurate central differences
    in the interior and second order accurate one-sides (forward or backwards)
    differences at the boundaries. The returned gradient hence has the same
    shape as the input array.

    Parameters
    ----------
    f : array_like
      An N-dimensional array containing samples of a scalar function.
    `*varargs` : scalars
      0, 1, or N scalars specifying the sample distances in each direction,
      that is: `dx`, `dy`, `dz`, ... The default distance is 1.

    Returns
    -------
    gradient : ndarray
      N arrays of the same shape as `f` giving the derivative of `f` with
      respect to each dimension.

    Examples
    --------
    >>> x = np.array([1, 2, 4, 7, 11, 16], dtype=np.float)
    >>> np.gradient(x)
    array([ 1. ,  1.5,  2.5,  3.5,  4.5,  5. ])
    >>> np.gradient(x, 2)
    array([ 0.5 ,  0.75,  1.25,  1.75,  2.25,  2.5 ])

    >>> np.gradient(np.array([[1, 2, 6], [3, 4, 5]], dtype=np.float))
    [array([[ 2.,  2., -1.],
           [ 2.,  2., -1.]]),
    array([[ 1. ,  2.5,  4. ],
           [ 1. ,  1. ,  1. ]])]

    >>> x = np.array([0,1,2,3,4])
    >>> dx = gradient(x)
    >>> y = x**2
    >>> gradient(y,dx)
    array([0.,  2.,  4.,  6.,  8.])
    """
    f = np.asanyarray(f)
    N = len(f.shape)  # number of dimensions
    n = len(varargs)
    if n == 0:
        dx = [1.0]*N
    elif n == 1:
        dx = [varargs[0]]*N
    elif n == N:
        dx = list(varargs)
    else:
        raise SyntaxError(
            "invalid number of arguments")

    # use central differences on interior and one-sided differences on the
    # endpoints. This preserves second order-accuracy over the full domain.

    outvals = []

    # create slice objects --- initially all are [:, :, ..., :]
    slice1 = [slice(None)]*N
    slice2 = [slice(None)]*N
    slice3 = [slice(None)]*N
    slice4 = [slice(None)]*N

    otype = f.dtype.char
    if otype not in ['f', 'd', 'F', 'D', 'm', 'M']:
        otype = 'd'

    # Difference of datetime64 elements results in timedelta64
    if otype == 'M':
        # Need to use the full dtype name because it contains unit information
        otype = f.dtype.name.replace('datetime', 'timedelta')
    elif otype == 'm':
        # Needs to keep the specific units, can't be a general unit
        otype = f.dtype

    # Convert datetime64 data into ints. Make dummy variable `y`
    # that is a view of ints if the data is datetime64, otherwise
    # just set y equal to the the array `f`.
    if f.dtype.char in ["M", "m"]:
        y = f.view('int64')
    else:
        y = f

    for axis in range(N):

        if y.shape[axis] < 2:
            raise ValueError(
                "Shape of array too small to calculate a numerical gradient, "
                "at least two elements are required.")

        # Numerical differentiation: 1st order edges, 2nd order interior
        if y.shape[axis] == 2:
            # Use first order differences for time data
            out = np.empty_like(y, dtype=otype)

            slice1[axis] = slice(1, -1)
            slice2[axis] = slice(2, None)
            slice3[axis] = slice(None, -2)
            # 1D equivalent -- out[1:-1] = (y[2:] - y[:-2])/2.0
            out[slice1] = (y[slice2] - y[slice3])/2.0

            slice1[axis] = 0
            slice2[axis] = 1
            slice3[axis] = 0
            # 1D equivalent -- out[0] = (y[1] - y[0])
            out[slice1] = (y[slice2] - y[slice3])

            slice1[axis] = -1
            slice2[axis] = -1
            slice3[axis] = -2
            # 1D equivalent -- out[-1] = (y[-1] - y[-2])
            out[slice1] = (y[slice2] - y[slice3])

        # Numerical differentiation: 2st order edges, 2nd order interior
        else:
            # Use second order differences where possible
            out = np.empty_like(y, dtype=otype)

            slice1[axis] = slice(1, -1)
            slice2[axis] = slice(2, None)
            slice3[axis] = slice(None, -2)
            # 1D equivalent -- out[1:-1] = (y[2:] - y[:-2])/2.0
            out[slice1] = (y[slice2] - y[slice3])/2.0

            slice1[axis] = 0
            slice2[axis] = 0
            slice3[axis] = 1
            slice4[axis] = 2
            # 1D equivalent -- out[0] = -(3*y[0] - 4*y[1] + y[2]) / 2.0
            out[slice1] = -(3.0*y[slice2] - 4.0*y[slice3] + y[slice4])/2.0

            slice1[axis] = -1
            slice2[axis] = -1
            slice3[axis] = -2
            slice4[axis] = -3
            # 1D equivalent -- out[-1] = (3*y[-1] - 4*y[-2] + y[-3])
            out[slice1] = (3.0*y[slice2] - 4.0*y[slice3] + y[slice4])/2.0

        # divide by step size
        outvals.append(out / dx[axis])

        # reset the slice object in this dimension to ":"
        slice1[axis] = slice(None)
        slice2[axis] = slice(None)
        slice3[axis] = slice(None)
        slice4[axis] = slice(None)

    if N == 1:
        return outvals[0]
    else:
        return outvals


def diff(a, n=1, axis=-1):
    """
    Calculate the n-th order discrete difference along given axis.

    The first order difference is given by ``out[n] = a[n+1] - a[n]`` along
    the given axis, higher order differences are calculated by using `diff`
    recursively.

    Parameters
    ----------
    a : array_like
        Input array
    n : int, optional
        The number of times values are differenced.
    axis : int, optional
        The axis along which the difference is taken, default is the last axis.

    Returns
    -------
    diff : ndarray
        The `n` order differences. The shape of the output is the same as `a`
        except along `axis` where the dimension is smaller by `n`.

    See Also
    --------
    gradient, ediff1d, cumsum

    Examples
    --------
    >>> x = np.array([1, 2, 4, 7, 0])
    >>> np.diff(x)
    array([ 1,  2,  3, -7])
    >>> np.diff(x, n=2)
    array([  1,   1, -10])

    >>> x = np.array([[1, 3, 6, 10], [0, 5, 6, 8]])
    >>> np.diff(x)
    array([[2, 3, 4],
           [5, 1, 2]])
    >>> np.diff(x, axis=0)
    array([[-1,  2,  0, -2]])

    """
    if n == 0:
        return a
    if n < 0:
        raise ValueError(
            "order must be non-negative but got " + repr(n))
    a = asanyarray(a)
    nd = len(a.shape)
    slice1 = [slice(None)]*nd
    slice2 = [slice(None)]*nd
    slice1[axis] = slice(1, None)
    slice2[axis] = slice(None, -1)
    slice1 = tuple(slice1)
    slice2 = tuple(slice2)
    if n > 1:
        return diff(a[slice1]-a[slice2], n-1, axis=axis)
    else:
        return a[slice1]-a[slice2]


def interp(x, xp, fp, left=None, right=None):
    """
    One-dimensional linear interpolation.

    Returns the one-dimensional piecewise linear interpolant to a function
    with given values at discrete data-points.

    Parameters
    ----------
    x : array_like
        The x-coordinates of the interpolated values.

    xp : 1-D sequence of floats
        The x-coordinates of the data points, must be increasing.

    fp : 1-D sequence of floats
        The y-coordinates of the data points, same length as `xp`.

    left : float, optional
        Value to return for `x < xp[0]`, default is `fp[0]`.

    right : float, optional
        Value to return for `x > xp[-1]`, default is `fp[-1]`.

    Returns
    -------
    y : {float, ndarray}
        The interpolated values, same shape as `x`.

    Raises
    ------
    ValueError
        If `xp` and `fp` have different length

    Notes
    -----
    Does not check that the x-coordinate sequence `xp` is increasing.
    If `xp` is not increasing, the results are nonsense.
    A simple check for increasing is::

        np.all(np.diff(xp) > 0)


    Examples
    --------
    >>> xp = [1, 2, 3]
    >>> fp = [3, 2, 0]
    >>> np.interp(2.5, xp, fp)
    1.0
    >>> np.interp([0, 1, 1.5, 2.72, 3.14], xp, fp)
    array([ 3. ,  3. ,  2.5 ,  0.56,  0. ])
    >>> UNDEF = -99.0
    >>> np.interp(3.14, xp, fp, right=UNDEF)
    -99.0

    Plot an interpolant to the sine function:

    >>> x = np.linspace(0, 2*np.pi, 10)
    >>> y = np.sin(x)
    >>> xvals = np.linspace(0, 2*np.pi, 50)
    >>> yinterp = np.interp(xvals, x, y)
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(x, y, 'o')
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.plot(xvals, yinterp, '-x')
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.show()

    """
    if isinstance(x, (float, int, number)):
        return compiled_interp([x], xp, fp, left, right).item()
    elif isinstance(x, np.ndarray) and x.ndim == 0:
        return compiled_interp([x], xp, fp, left, right).item()
    else:
        return compiled_interp(x, xp, fp, left, right)


def angle(z, deg=0):
    """
    Return the angle of the complex argument.

    Parameters
    ----------
    z : array_like
        A complex number or sequence of complex numbers.
    deg : bool, optional
        Return angle in degrees if True, radians if False (default).

    Returns
    -------
    angle : {ndarray, scalar}
        The counterclockwise angle from the positive real axis on
        the complex plane, with dtype as numpy.float64.

    See Also
    --------
    arctan2
    absolute



    Examples
    --------
    >>> np.angle([1.0, 1.0j, 1+1j])               # in radians
    array([ 0.        ,  1.57079633,  0.78539816])
    >>> np.angle(1+1j, deg=True)                  # in degrees
    45.0

    """
    if deg:
        fact = 180/pi
    else:
        fact = 1.0
    z = asarray(z)
    if (issubclass(z.dtype.type, _nx.complexfloating)):
        zimag = z.imag
        zreal = z.real
    else:
        zimag = 0
        zreal = z
    return arctan2(zimag, zreal) * fact


def unwrap(p, discont=pi, axis=-1):
    """
    Unwrap by changing deltas between values to 2*pi complement.

    Unwrap radian phase `p` by changing absolute jumps greater than
    `discont` to their 2*pi complement along the given axis.

    Parameters
    ----------
    p : array_like
        Input array.
    discont : float, optional
        Maximum discontinuity between values, default is ``pi``.
    axis : int, optional
        Axis along which unwrap will operate, default is the last axis.

    Returns
    -------
    out : ndarray
        Output array.

    See Also
    --------
    rad2deg, deg2rad

    Notes
    -----
    If the discontinuity in `p` is smaller than ``pi``, but larger than
    `discont`, no unwrapping is done because taking the 2*pi complement
    would only make the discontinuity larger.

    Examples
    --------
    >>> phase = np.linspace(0, np.pi, num=5)
    >>> phase[3:] += np.pi
    >>> phase
    array([ 0.        ,  0.78539816,  1.57079633,  5.49778714,  6.28318531])
    >>> np.unwrap(phase)
    array([ 0.        ,  0.78539816,  1.57079633, -0.78539816,  0.        ])

    """
    p = asarray(p)
    nd = len(p.shape)
    dd = diff(p, axis=axis)
    slice1 = [slice(None, None)]*nd     # full slices
    slice1[axis] = slice(1, None)
    ddmod = mod(dd + pi, 2*pi) - pi
    _nx.copyto(ddmod, pi, where=(ddmod == -pi) & (dd > 0))
    ph_correct = ddmod - dd
    _nx.copyto(ph_correct, 0, where=abs(dd) < discont)
    up = array(p, copy=True, dtype='d')
    up[slice1] = p[slice1] + ph_correct.cumsum(axis)
    return up


def sort_complex(a):
    """
    Sort a complex array using the real part first, then the imaginary part.

    Parameters
    ----------
    a : array_like
        Input array

    Returns
    -------
    out : complex ndarray
        Always returns a sorted complex array.

    Examples
    --------
    >>> np.sort_complex([5, 3, 6, 2, 1])
    array([ 1.+0.j,  2.+0.j,  3.+0.j,  5.+0.j,  6.+0.j])

    >>> np.sort_complex([1 + 2j, 2 - 1j, 3 - 2j, 3 - 3j, 3 + 5j])
    array([ 1.+2.j,  2.-1.j,  3.-3.j,  3.-2.j,  3.+5.j])

    """
    b = array(a, copy=True)
    b.sort()
    if not issubclass(b.dtype.type, _nx.complexfloating):
        if b.dtype.char in 'bhBH':
            return b.astype('F')
        elif b.dtype.char == 'g':
            return b.astype('G')
        else:
            return b.astype('D')
    else:
        return b


def trim_zeros(filt, trim='fb'):
    """
    Trim the leading and/or trailing zeros from a 1-D array or sequence.

    Parameters
    ----------
    filt : 1-D array or sequence
        Input array.
    trim : str, optional
        A string with 'f' representing trim from front and 'b' to trim from
        back. Default is 'fb', trim zeros from both front and back of the
        array.

    Returns
    -------
    trimmed : 1-D array or sequence
        The result of trimming the input. The input data type is preserved.

    Examples
    --------
    >>> a = np.array((0, 0, 0, 1, 2, 3, 0, 2, 1, 0))
    >>> np.trim_zeros(a)
    array([1, 2, 3, 0, 2, 1])

    >>> np.trim_zeros(a, 'b')
    array([0, 0, 0, 1, 2, 3, 0, 2, 1])

    The input data type is preserved, list/tuple in means list/tuple out.

    >>> np.trim_zeros([0, 1, 2, 0])
    [1, 2]

    """
    first = 0
    trim = trim.upper()
    if 'F' in trim:
        for i in filt:
            if i != 0.:
                break
            else:
                first = first + 1
    last = len(filt)
    if 'B' in trim:
        for i in filt[::-1]:
            if i != 0.:
                break
            else:
                last = last - 1
    return filt[first:last]


@deprecate
def unique(x):
    """
    This function is deprecated.  Use numpy.lib.arraysetops.unique()
    instead.
    """
    try:
        tmp = x.flatten()
        if tmp.size == 0:
            return tmp
        tmp.sort()
        idx = concatenate(([True], tmp[1:] != tmp[:-1]))
        return tmp[idx]
    except AttributeError:
        items = sorted(set(x))
        return asarray(items)


def extract(condition, arr):
    """
    Return the elements of an array that satisfy some condition.

    This is equivalent to ``np.compress(ravel(condition), ravel(arr))``.  If
    `condition` is boolean ``np.extract`` is equivalent to ``arr[condition]``.

    Parameters
    ----------
    condition : array_like
        An array whose nonzero or True entries indicate the elements of `arr`
        to extract.
    arr : array_like
        Input array of the same size as `condition`.

    Returns
    -------
    extract : ndarray
        Rank 1 array of values from `arr` where `condition` is True.

    See Also
    --------
    take, put, copyto, compress

    Examples
    --------
    >>> arr = np.arange(12).reshape((3, 4))
    >>> arr
    array([[ 0,  1,  2,  3],
           [ 4,  5,  6,  7],
           [ 8,  9, 10, 11]])
    >>> condition = np.mod(arr, 3)==0
    >>> condition
    array([[ True, False, False,  True],
           [False, False,  True, False],
           [False,  True, False, False]], dtype=bool)
    >>> np.extract(condition, arr)
    array([0, 3, 6, 9])


    If `condition` is boolean:

    >>> arr[condition]
    array([0, 3, 6, 9])

    """
    return _nx.take(ravel(arr), nonzero(ravel(condition))[0])


def place(arr, mask, vals):
    """
    Change elements of an array based on conditional and input values.

    Similar to ``np.copyto(arr, vals, where=mask)``, the difference is that
    `place` uses the first N elements of `vals`, where N is the number of
    True values in `mask`, while `copyto` uses the elements where `mask`
    is True.

    Note that `extract` does the exact opposite of `place`.

    Parameters
    ----------
    arr : array_like
        Array to put data into.
    mask : array_like
        Boolean mask array. Must have the same size as `a`.
    vals : 1-D sequence
        Values to put into `a`. Only the first N elements are used, where
        N is the number of True values in `mask`. If `vals` is smaller
        than N it will be repeated.

    See Also
    --------
    copyto, put, take, extract

    Examples
    --------
    >>> arr = np.arange(6).reshape(2, 3)
    >>> np.place(arr, arr>2, [44, 55])
    >>> arr
    array([[ 0,  1,  2],
           [44, 55, 44]])

    """
    return _insert(arr, mask, vals)


def disp(mesg, device=None, linefeed=True):
    """
    Display a message on a device.

    Parameters
    ----------
    mesg : str
        Message to display.
    device : object
        Device to write message. If None, defaults to ``sys.stdout`` which is
        very similar to ``print``. `device` needs to have ``write()`` and
        ``flush()`` methods.
    linefeed : bool, optional
        Option whether to print a line feed or not. Defaults to True.

    Raises
    ------
    AttributeError
        If `device` does not have a ``write()`` or ``flush()`` method.

    Examples
    --------
    Besides ``sys.stdout``, a file-like object can also be used as it has
    both required methods:

    >>> from StringIO import StringIO
    >>> buf = StringIO()
    >>> np.disp('"Display" in a file', device=buf)
    >>> buf.getvalue()
    '"Display" in a file\\n'

    """
    if device is None:
        device = sys.stdout
    if linefeed:
        device.write('%s\n' % mesg)
    else:
        device.write('%s' % mesg)
    device.flush()
    return


class vectorize(object):
    """
    vectorize(pyfunc, otypes='', doc=None, excluded=None, cache=False)

    Generalized function class.

    Define a vectorized function which takes a nested sequence
    of objects or numpy arrays as inputs and returns a
    numpy array as output. The vectorized function evaluates `pyfunc` over
    successive tuples of the input arrays like the python map function,
    except it uses the broadcasting rules of numpy.

    The data type of the output of `vectorized` is determined by calling
    the function with the first element of the input.  This can be avoided
    by specifying the `otypes` argument.

    Parameters
    ----------
    pyfunc : callable
        A python function or method.
    otypes : str or list of dtypes, optional
        The output data type. It must be specified as either a string of
        typecode characters or a list of data type specifiers. There should
        be one data type specifier for each output.
    doc : str, optional
        The docstring for the function. If `None`, the docstring will be the
        ``pyfunc.__doc__``.
    excluded : set, optional
        Set of strings or integers representing the positional or keyword
        arguments for which the function will not be vectorized.  These will be
        passed directly to `pyfunc` unmodified.

        .. versionadded:: 1.7.0

    cache : bool, optional
       If `True`, then cache the first function call that determines the number
       of outputs if `otypes` is not provided.

        .. versionadded:: 1.7.0

    Returns
    -------
    vectorized : callable
        Vectorized function.

    Examples
    --------
    >>> def myfunc(a, b):
    ...     "Return a-b if a>b, otherwise return a+b"
    ...     if a > b:
    ...         return a - b
    ...     else:
    ...         return a + b

    >>> vfunc = np.vectorize(myfunc)
    >>> vfunc([1, 2, 3, 4], 2)
    array([3, 4, 1, 2])

    The docstring is taken from the input function to `vectorize` unless it
    is specified

    >>> vfunc.__doc__
    'Return a-b if a>b, otherwise return a+b'
    >>> vfunc = np.vectorize(myfunc, doc='Vectorized `myfunc`')
    >>> vfunc.__doc__
    'Vectorized `myfunc`'

    The output type is determined by evaluating the first element of the input,
    unless it is specified

    >>> out = vfunc([1, 2, 3, 4], 2)
    >>> type(out[0])
    <type 'numpy.int32'>
    >>> vfunc = np.vectorize(myfunc, otypes=[np.float])
    >>> out = vfunc([1, 2, 3, 4], 2)
    >>> type(out[0])
    <type 'numpy.float64'>

    The `excluded` argument can be used to prevent vectorizing over certain
    arguments.  This can be useful for array-like arguments of a fixed length
    such as the coefficients for a polynomial as in `polyval`:

    >>> def mypolyval(p, x):
    ...     _p = list(p)
    ...     res = _p.pop(0)
    ...     while _p:
    ...         res = res*x + _p.pop(0)
    ...     return res
    >>> vpolyval = np.vectorize(mypolyval, excluded=['p'])
    >>> vpolyval(p=[1, 2, 3], x=[0, 1])
    array([3, 6])

    Positional arguments may also be excluded by specifying their position:

    >>> vpolyval.excluded.add(0)
    >>> vpolyval([1, 2, 3], x=[0, 1])
    array([3, 6])

    Notes
    -----
    The `vectorize` function is provided primarily for convenience, not for
    performance. The implementation is essentially a for loop.

    If `otypes` is not specified, then a call to the function with the
    first argument will be used to determine the number of outputs.  The
    results of this call will be cached if `cache` is `True` to prevent
    calling the function twice.  However, to implement the cache, the
    original function must be wrapped which will slow down subsequent
    calls, so only do this if your function is expensive.

    The new keyword argument interface and `excluded` argument support
    further degrades performance.

    """

    def __init__(self, pyfunc, otypes='', doc=None, excluded=None,
                 cache=False):
        self.pyfunc = pyfunc
        self.cache = cache
        self._ufunc = None    # Caching to improve default performance

        if doc is None:
            self.__doc__ = pyfunc.__doc__
        else:
            self.__doc__ = doc

        if isinstance(otypes, str):
            self.otypes = otypes
            for char in self.otypes:
                if char not in typecodes['All']:
                    raise ValueError(
                        "Invalid otype specified: %s" % (char,))
        elif iterable(otypes):
            self.otypes = ''.join([_nx.dtype(x).char for x in otypes])
        else:
            raise ValueError(
                "Invalid otype specification")

        # Excluded variable support
        if excluded is None:
            excluded = set()
        self.excluded = set(excluded)

    def __call__(self, *args, **kwargs):
        """
        Return arrays with the results of `pyfunc` broadcast (vectorized) over
        `args` and `kwargs` not in `excluded`.
        """
        excluded = self.excluded
        if not kwargs and not excluded:
            func = self.pyfunc
            vargs = args
        else:
            # The wrapper accepts only positional arguments: we use `names` and
            # `inds` to mutate `the_args` and `kwargs` to pass to the original
            # function.
            nargs = len(args)

            names = [_n for _n in kwargs if _n not in excluded]
            inds = [_i for _i in range(nargs) if _i not in excluded]
            the_args = list(args)

            def func(*vargs):
                for _n, _i in enumerate(inds):
                    the_args[_i] = vargs[_n]
                kwargs.update(zip(names, vargs[len(inds):]))
                return self.pyfunc(*the_args, **kwargs)

            vargs = [args[_i] for _i in inds]
            vargs.extend([kwargs[_n] for _n in names])

        return self._vectorize_call(func=func, args=vargs)

    def _get_ufunc_and_otypes(self, func, args):
        """Return (ufunc, otypes)."""
        # frompyfunc will fail if args is empty
        if not args:
            raise ValueError('args can not be empty')

        if self.otypes:
            otypes = self.otypes
            nout = len(otypes)

            # Note logic here: We only *use* self._ufunc if func is self.pyfunc
            # even though we set self._ufunc regardless.
            if func is self.pyfunc and self._ufunc is not None:
                ufunc = self._ufunc
            else:
                ufunc = self._ufunc = frompyfunc(func, len(args), nout)
        else:
            # Get number of outputs and output types by calling the function on
            # the first entries of args.  We also cache the result to prevent
            # the subsequent call when the ufunc is evaluated.
            # Assumes that ufunc first evaluates the 0th elements in the input
            # arrays (the input values are not checked to ensure this)
            inputs = [asarray(_a).flat[0] for _a in args]
            outputs = func(*inputs)

            # Performance note: profiling indicates that -- for simple
            # functions at least -- this wrapping can almost double the
            # execution time.
            # Hence we make it optional.
            if self.cache:
                _cache = [outputs]

                def _func(*vargs):
                    if _cache:
                        return _cache.pop()
                    else:
                        return func(*vargs)
            else:
                _func = func

            if isinstance(outputs, tuple):
                nout = len(outputs)
            else:
                nout = 1
                outputs = (outputs,)

            otypes = ''.join([asarray(outputs[_k]).dtype.char
                              for _k in range(nout)])

            # Performance note: profiling indicates that creating the ufunc is
            # not a significant cost compared with wrapping so it seems not
            # worth trying to cache this.
            ufunc = frompyfunc(_func, len(args), nout)

        return ufunc, otypes

    def _vectorize_call(self, func, args):
        """Vectorized call to `func` over positional `args`."""
        if not args:
            _res = func()
        else:
            ufunc, otypes = self._get_ufunc_and_otypes(func=func, args=args)

            # Convert args to object arrays first
            inputs = [array(_a, copy=False, subok=True, dtype=object)
                      for _a in args]

            outputs = ufunc(*inputs)

            if ufunc.nout == 1:
                _res = array(outputs,
                             copy=False, subok=True, dtype=otypes[0])
            else:
                _res = tuple([array(_x, copy=False, subok=True, dtype=_t)
                              for _x, _t in zip(outputs, otypes)])
        return _res


def cov(m, y=None, rowvar=1, bias=0, ddof=None):
    """
    Estimate a covariance matrix, given data.

    Covariance indicates the level to which two variables vary together.
    If we examine N-dimensional samples, :math:`X = [x_1, x_2, ... x_N]^T`,
    then the covariance matrix element :math:`C_{ij}` is the covariance of
    :math:`x_i` and :math:`x_j`. The element :math:`C_{ii}` is the variance
    of :math:`x_i`.

    Parameters
    ----------
    m : array_like
        A 1-D or 2-D array containing multiple variables and observations.
        Each row of `m` represents a variable, and each column a single
        observation of all those variables. Also see `rowvar` below.
    y : array_like, optional
        An additional set of variables and observations. `y` has the same
        form as that of `m`.
    rowvar : int, optional
        If `rowvar` is non-zero (default), then each row represents a
        variable, with observations in the columns. Otherwise, the relationship
        is transposed: each column represents a variable, while the rows
        contain observations.
    bias : int, optional
        Default normalization is by ``(N - 1)``, where ``N`` is the number of
        observations given (unbiased estimate). If `bias` is 1, then
        normalization is by ``N``. These values can be overridden by using
        the keyword ``ddof`` in numpy versions >= 1.5.
    ddof : int, optional
        .. versionadded:: 1.5
        If not ``None`` normalization is by ``(N - ddof)``, where ``N`` is
        the number of observations; this overrides the value implied by
        ``bias``. The default value is ``None``.

    Returns
    -------
    out : ndarray
        The covariance matrix of the variables.

    See Also
    --------
    corrcoef : Normalized covariance matrix

    Examples
    --------
    Consider two variables, :math:`x_0` and :math:`x_1`, which
    correlate perfectly, but in opposite directions:

    >>> x = np.array([[0, 2], [1, 1], [2, 0]]).T
    >>> x
    array([[0, 1, 2],
           [2, 1, 0]])

    Note how :math:`x_0` increases while :math:`x_1` decreases. The covariance
    matrix shows this clearly:

    >>> np.cov(x)
    array([[ 1., -1.],
           [-1.,  1.]])

    Note that element :math:`C_{0,1}`, which shows the correlation between
    :math:`x_0` and :math:`x_1`, is negative.

    Further, note how `x` and `y` are combined:

    >>> x = [-2.1, -1,  4.3]
    >>> y = [3,  1.1,  0.12]
    >>> X = np.vstack((x,y))
    >>> print np.cov(X)
    [[ 11.71        -4.286     ]
     [ -4.286        2.14413333]]
    >>> print np.cov(x, y)
    [[ 11.71        -4.286     ]
     [ -4.286        2.14413333]]
    >>> print np.cov(x)
    11.71

    """
    # Check inputs
    if ddof is not None and ddof != int(ddof):
        raise ValueError(
            "ddof must be integer")

    # Handles complex arrays too
    m = np.asarray(m)
    if y is None:
        dtype = np.result_type(m, np.float64)
    else:
        y = np.asarray(y)
        dtype = np.result_type(m, y, np.float64)
    X = array(m, ndmin=2, dtype=dtype)

    if X.shape[0] == 1:
        rowvar = 1
    if rowvar:
        N = X.shape[1]
        axis = 0
    else:
        N = X.shape[0]
        axis = 1

    # check ddof
    if ddof is None:
        if bias == 0:
            ddof = 1
        else:
            ddof = 0
    fact = float(N - ddof)
    if fact <= 0:
        warnings.warn("Degrees of freedom <= 0 for slice", RuntimeWarning)
        fact = 0.0

    if y is not None:
        y = array(y, copy=False, ndmin=2, dtype=dtype)
        X = concatenate((X, y), axis)

    X -= X.mean(axis=1-axis, keepdims=True)
    if not rowvar:
        return (dot(X.T, X.conj()) / fact).squeeze()
    else:
        return (dot(X, X.T.conj()) / fact).squeeze()


def corrcoef(x, y=None, rowvar=1, bias=0, ddof=None):
    """
    Return correlation coefficients.

    Please refer to the documentation for `cov` for more detail.  The
    relationship between the correlation coefficient matrix, `P`, and the
    covariance matrix, `C`, is

    .. math:: P_{ij} = \\frac{ C_{ij} } { \\sqrt{ C_{ii} * C_{jj} } }

    The values of `P` are between -1 and 1, inclusive.

    Parameters
    ----------
    x : array_like
        A 1-D or 2-D array containing multiple variables and observations.
        Each row of `m` represents a variable, and each column a single
        observation of all those variables. Also see `rowvar` below.
    y : array_like, optional
        An additional set of variables and observations. `y` has the same
        shape as `m`.
    rowvar : int, optional
        If `rowvar` is non-zero (default), then each row represents a
        variable, with observations in the columns. Otherwise, the relationship
        is transposed: each column represents a variable, while the rows
        contain observations.
    bias : int, optional
        Default normalization is by ``(N - 1)``, where ``N`` is the number of
        observations (unbiased estimate). If `bias` is 1, then
        normalization is by ``N``. These values can be overridden by using
        the keyword ``ddof`` in numpy versions >= 1.5.
    ddof : {None, int}, optional
        .. versionadded:: 1.5
        If not ``None`` normalization is by ``(N - ddof)``, where ``N`` is
        the number of observations; this overrides the value implied by
        ``bias``. The default value is ``None``.

    Returns
    -------
    out : ndarray
        The correlation coefficient matrix of the variables.

    See Also
    --------
    cov : Covariance matrix

    """
    c = cov(x, y, rowvar, bias, ddof)
    try:
        d = diag(c)
    except ValueError:  # scalar covariance
        # nan if incorrect value (nan, inf, 0), 1 otherwise
        return c / c
    return c / sqrt(multiply.outer(d, d))


def blackman(M):
    """
    Return the Blackman window.

    The Blackman window is a taper formed by using the first three
    terms of a summation of cosines. It was designed to have close to the
    minimal leakage possible.  It is close to optimal, only slightly worse
    than a Kaiser window.

    Parameters
    ----------
    M : int
        Number of points in the output window. If zero or less, an empty
        array is returned.

    Returns
    -------
    out : ndarray
        The window, with the maximum value normalized to one (the value one
        appears only if the number of samples is odd).

    See Also
    --------
    bartlett, hamming, hanning, kaiser

    Notes
    -----
    The Blackman window is defined as

    .. math::  w(n) = 0.42 - 0.5 \\cos(2\\pi n/M) + 0.08 \\cos(4\\pi n/M)

    Most references to the Blackman window come from the signal processing
    literature, where it is used as one of many windowing functions for
    smoothing values.  It is also known as an apodization (which means
    "removing the foot", i.e. smoothing discontinuities at the beginning
    and end of the sampled signal) or tapering function. It is known as a
    "near optimal" tapering function, almost as good (by some measures)
    as the kaiser window.

    References
    ----------
    Blackman, R.B. and Tukey, J.W., (1958) The measurement of power spectra,
    Dover Publications, New York.

    Oppenheim, A.V., and R.W. Schafer. Discrete-Time Signal Processing.
    Upper Saddle River, NJ: Prentice-Hall, 1999, pp. 468-471.

    Examples
    --------
    >>> np.blackman(12)
    array([ -1.38777878e-17,   3.26064346e-02,   1.59903635e-01,
             4.14397981e-01,   7.36045180e-01,   9.67046769e-01,
             9.67046769e-01,   7.36045180e-01,   4.14397981e-01,
             1.59903635e-01,   3.26064346e-02,  -1.38777878e-17])


    Plot the window and the frequency response:

    >>> from numpy.fft import fft, fftshift
    >>> window = np.blackman(51)
    >>> plt.plot(window)
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.title("Blackman window")
    <matplotlib.text.Text object at 0x...>
    >>> plt.ylabel("Amplitude")
    <matplotlib.text.Text object at 0x...>
    >>> plt.xlabel("Sample")
    <matplotlib.text.Text object at 0x...>
    >>> plt.show()

    >>> plt.figure()
    <matplotlib.figure.Figure object at 0x...>
    >>> A = fft(window, 2048) / 25.5
    >>> mag = np.abs(fftshift(A))
    >>> freq = np.linspace(-0.5, 0.5, len(A))
    >>> response = 20 * np.log10(mag)
    >>> response = np.clip(response, -100, 100)
    >>> plt.plot(freq, response)
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.title("Frequency response of Blackman window")
    <matplotlib.text.Text object at 0x...>
    >>> plt.ylabel("Magnitude [dB]")
    <matplotlib.text.Text object at 0x...>
    >>> plt.xlabel("Normalized frequency [cycles per sample]")
    <matplotlib.text.Text object at 0x...>
    >>> plt.axis('tight')
    (-0.5, 0.5, -100.0, ...)
    >>> plt.show()

    """
    if M < 1:
        return array([])
    if M == 1:
        return ones(1, float)
    n = arange(0, M)
    return 0.42 - 0.5*cos(2.0*pi*n/(M-1)) + 0.08*cos(4.0*pi*n/(M-1))


def bartlett(M):
    """
    Return the Bartlett window.

    The Bartlett window is very similar to a triangular window, except
    that the end points are at zero.  It is often used in signal
    processing for tapering a signal, without generating too much
    ripple in the frequency domain.

    Parameters
    ----------
    M : int
        Number of points in the output window. If zero or less, an
        empty array is returned.

    Returns
    -------
    out : array
        The triangular window, with the maximum value normalized to one
        (the value one appears only if the number of samples is odd), with
        the first and last samples equal to zero.

    See Also
    --------
    blackman, hamming, hanning, kaiser

    Notes
    -----
    The Bartlett window is defined as

    .. math:: w(n) = \\frac{2}{M-1} \\left(
              \\frac{M-1}{2} - \\left|n - \\frac{M-1}{2}\\right|
              \\right)

    Most references to the Bartlett window come from the signal
    processing literature, where it is used as one of many windowing
    functions for smoothing values.  Note that convolution with this
    window produces linear interpolation.  It is also known as an
    apodization (which means"removing the foot", i.e. smoothing
    discontinuities at the beginning and end of the sampled signal) or
    tapering function. The fourier transform of the Bartlett is the product
    of two sinc functions.
    Note the excellent discussion in Kanasewich.

    References
    ----------
    .. [1] M.S. Bartlett, "Periodogram Analysis and Continuous Spectra",
           Biometrika 37, 1-16, 1950.
    .. [2] E.R. Kanasewich, "Time Sequence Analysis in Geophysics",
           The University of Alberta Press, 1975, pp. 109-110.
    .. [3] A.V. Oppenheim and R.W. Schafer, "Discrete-Time Signal
           Processing", Prentice-Hall, 1999, pp. 468-471.
    .. [4] Wikipedia, "Window function",
           http://en.wikipedia.org/wiki/Window_function
    .. [5] W.H. Press,  B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling,
           "Numerical Recipes", Cambridge University Press, 1986, page 429.


    Examples
    --------
    >>> np.bartlett(12)
    array([ 0.        ,  0.18181818,  0.36363636,  0.54545455,  0.72727273,
            0.90909091,  0.90909091,  0.72727273,  0.54545455,  0.36363636,
            0.18181818,  0.        ])

    Plot the window and its frequency response (requires SciPy and matplotlib):

    >>> from numpy.fft import fft, fftshift
    >>> window = np.bartlett(51)
    >>> plt.plot(window)
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.title("Bartlett window")
    <matplotlib.text.Text object at 0x...>
    >>> plt.ylabel("Amplitude")
    <matplotlib.text.Text object at 0x...>
    >>> plt.xlabel("Sample")
    <matplotlib.text.Text object at 0x...>
    >>> plt.show()

    >>> plt.figure()
    <matplotlib.figure.Figure object at 0x...>
    >>> A = fft(window, 2048) / 25.5
    >>> mag = np.abs(fftshift(A))
    >>> freq = np.linspace(-0.5, 0.5, len(A))
    >>> response = 20 * np.log10(mag)
    >>> response = np.clip(response, -100, 100)
    >>> plt.plot(freq, response)
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.title("Frequency response of Bartlett window")
    <matplotlib.text.Text object at 0x...>
    >>> plt.ylabel("Magnitude [dB]")
    <matplotlib.text.Text object at 0x...>
    >>> plt.xlabel("Normalized frequency [cycles per sample]")
    <matplotlib.text.Text object at 0x...>
    >>> plt.axis('tight')
    (-0.5, 0.5, -100.0, ...)
    >>> plt.show()

    """
    if M < 1:
        return array([])
    if M == 1:
        return ones(1, float)
    n = arange(0, M)
    return where(less_equal(n, (M-1)/2.0), 2.0*n/(M-1), 2.0 - 2.0*n/(M-1))


def hanning(M):
    """
    Return the Hanning window.

    The Hanning window is a taper formed by using a weighted cosine.

    Parameters
    ----------
    M : int
        Number of points in the output window. If zero or less, an
        empty array is returned.

    Returns
    -------
    out : ndarray, shape(M,)
        The window, with the maximum value normalized to one (the value
        one appears only if `M` is odd).

    See Also
    --------
    bartlett, blackman, hamming, kaiser

    Notes
    -----
    The Hanning window is defined as

    .. math::  w(n) = 0.5 - 0.5cos\\left(\\frac{2\\pi{n}}{M-1}\\right)
               \\qquad 0 \\leq n \\leq M-1

    The Hanning was named for Julius van Hann, an Austrian meteorologist.
    It is also known as the Cosine Bell. Some authors prefer that it be
    called a Hann window, to help avoid confusion with the very similar
    Hamming window.

    Most references to the Hanning window come from the signal processing
    literature, where it is used as one of many windowing functions for
    smoothing values.  It is also known as an apodization (which means
    "removing the foot", i.e. smoothing discontinuities at the beginning
    and end of the sampled signal) or tapering function.

    References
    ----------
    .. [1] Blackman, R.B. and Tukey, J.W., (1958) The measurement of power
           spectra, Dover Publications, New York.
    .. [2] E.R. Kanasewich, "Time Sequence Analysis in Geophysics",
           The University of Alberta Press, 1975, pp. 106-108.
    .. [3] Wikipedia, "Window function",
           http://en.wikipedia.org/wiki/Window_function
    .. [4] W.H. Press,  B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling,
           "Numerical Recipes", Cambridge University Press, 1986, page 425.

    Examples
    --------
    >>> np.hanning(12)
    array([ 0.        ,  0.07937323,  0.29229249,  0.57115742,  0.82743037,
            0.97974649,  0.97974649,  0.82743037,  0.57115742,  0.29229249,
            0.07937323,  0.        ])

    Plot the window and its frequency response:

    >>> from numpy.fft import fft, fftshift
    >>> window = np.hanning(51)
    >>> plt.plot(window)
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.title("Hann window")
    <matplotlib.text.Text object at 0x...>
    >>> plt.ylabel("Amplitude")
    <matplotlib.text.Text object at 0x...>
    >>> plt.xlabel("Sample")
    <matplotlib.text.Text object at 0x...>
    >>> plt.show()

    >>> plt.figure()
    <matplotlib.figure.Figure object at 0x...>
    >>> A = fft(window, 2048) / 25.5
    >>> mag = np.abs(fftshift(A))
    >>> freq = np.linspace(-0.5, 0.5, len(A))
    >>> response = 20 * np.log10(mag)
    >>> response = np.clip(response, -100, 100)
    >>> plt.plot(freq, response)
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.title("Frequency response of the Hann window")
    <matplotlib.text.Text object at 0x...>
    >>> plt.ylabel("Magnitude [dB]")
    <matplotlib.text.Text object at 0x...>
    >>> plt.xlabel("Normalized frequency [cycles per sample]")
    <matplotlib.text.Text object at 0x...>
    >>> plt.axis('tight')
    (-0.5, 0.5, -100.0, ...)
    >>> plt.show()

    """
    if M < 1:
        return array([])
    if M == 1:
        return ones(1, float)
    n = arange(0, M)
    return 0.5 - 0.5*cos(2.0*pi*n/(M-1))


def hamming(M):
    """
    Return the Hamming window.

    The Hamming window is a taper formed by using a weighted cosine.

    Parameters
    ----------
    M : int
        Number of points in the output window. If zero or less, an
        empty array is returned.

    Returns
    -------
    out : ndarray
        The window, with the maximum value normalized to one (the value
        one appears only if the number of samples is odd).

    See Also
    --------
    bartlett, blackman, hanning, kaiser

    Notes
    -----
    The Hamming window is defined as

    .. math::  w(n) = 0.54 - 0.46cos\\left(\\frac{2\\pi{n}}{M-1}\\right)
               \\qquad 0 \\leq n \\leq M-1

    The Hamming was named for R. W. Hamming, an associate of J. W. Tukey
    and is described in Blackman and Tukey. It was recommended for
    smoothing the truncated autocovariance function in the time domain.
    Most references to the Hamming window come from the signal processing
    literature, where it is used as one of many windowing functions for
    smoothing values.  It is also known as an apodization (which means
    "removing the foot", i.e. smoothing discontinuities at the beginning
    and end of the sampled signal) or tapering function.

    References
    ----------
    .. [1] Blackman, R.B. and Tukey, J.W., (1958) The measurement of power
           spectra, Dover Publications, New York.
    .. [2] E.R. Kanasewich, "Time Sequence Analysis in Geophysics", The
           University of Alberta Press, 1975, pp. 109-110.
    .. [3] Wikipedia, "Window function",
           http://en.wikipedia.org/wiki/Window_function
    .. [4] W.H. Press,  B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling,
           "Numerical Recipes", Cambridge University Press, 1986, page 425.

    Examples
    --------
    >>> np.hamming(12)
    array([ 0.08      ,  0.15302337,  0.34890909,  0.60546483,  0.84123594,
            0.98136677,  0.98136677,  0.84123594,  0.60546483,  0.34890909,
            0.15302337,  0.08      ])

    Plot the window and the frequency response:

    >>> from numpy.fft import fft, fftshift
    >>> window = np.hamming(51)
    >>> plt.plot(window)
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.title("Hamming window")
    <matplotlib.text.Text object at 0x...>
    >>> plt.ylabel("Amplitude")
    <matplotlib.text.Text object at 0x...>
    >>> plt.xlabel("Sample")
    <matplotlib.text.Text object at 0x...>
    >>> plt.show()

    >>> plt.figure()
    <matplotlib.figure.Figure object at 0x...>
    >>> A = fft(window, 2048) / 25.5
    >>> mag = np.abs(fftshift(A))
    >>> freq = np.linspace(-0.5, 0.5, len(A))
    >>> response = 20 * np.log10(mag)
    >>> response = np.clip(response, -100, 100)
    >>> plt.plot(freq, response)
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.title("Frequency response of Hamming window")
    <matplotlib.text.Text object at 0x...>
    >>> plt.ylabel("Magnitude [dB]")
    <matplotlib.text.Text object at 0x...>
    >>> plt.xlabel("Normalized frequency [cycles per sample]")
    <matplotlib.text.Text object at 0x...>
    >>> plt.axis('tight')
    (-0.5, 0.5, -100.0, ...)
    >>> plt.show()

    """
    if M < 1:
        return array([])
    if M == 1:
        return ones(1, float)
    n = arange(0, M)
    return 0.54 - 0.46*cos(2.0*pi*n/(M-1))

## Code from cephes for i0

_i0A = [
    -4.41534164647933937950E-18,
    3.33079451882223809783E-17,
    -2.43127984654795469359E-16,
    1.71539128555513303061E-15,
    -1.16853328779934516808E-14,
    7.67618549860493561688E-14,
    -4.85644678311192946090E-13,
    2.95505266312963983461E-12,
    -1.72682629144155570723E-11,
    9.67580903537323691224E-11,
    -5.18979560163526290666E-10,
    2.65982372468238665035E-9,
    -1.30002500998624804212E-8,
    6.04699502254191894932E-8,
    -2.67079385394061173391E-7,
    1.11738753912010371815E-6,
    -4.41673835845875056359E-6,
    1.64484480707288970893E-5,
    -5.75419501008210370398E-5,
    1.88502885095841655729E-4,
    -5.76375574538582365885E-4,
    1.63947561694133579842E-3,
    -4.32430999505057594430E-3,
    1.05464603945949983183E-2,
    -2.37374148058994688156E-2,
    4.93052842396707084878E-2,
    -9.49010970480476444210E-2,
    1.71620901522208775349E-1,
    -3.04682672343198398683E-1,
    6.76795274409476084995E-1
    ]

_i0B = [
    -7.23318048787475395456E-18,
    -4.83050448594418207126E-18,
    4.46562142029675999901E-17,
    3.46122286769746109310E-17,
    -2.82762398051658348494E-16,
    -3.42548561967721913462E-16,
    1.77256013305652638360E-15,
    3.81168066935262242075E-15,
    -9.55484669882830764870E-15,
    -4.15056934728722208663E-14,
    1.54008621752140982691E-14,
    3.85277838274214270114E-13,
    7.18012445138366623367E-13,
    -1.79417853150680611778E-12,
    -1.32158118404477131188E-11,
    -3.14991652796324136454E-11,
    1.18891471078464383424E-11,
    4.94060238822496958910E-10,
    3.39623202570838634515E-9,
    2.26666899049817806459E-8,
    2.04891858946906374183E-7,
    2.89137052083475648297E-6,
    6.88975834691682398426E-5,
    3.36911647825569408990E-3,
    8.04490411014108831608E-1
    ]


def _chbevl(x, vals):
    b0 = vals[0]
    b1 = 0.0

    for i in range(1, len(vals)):
        b2 = b1
        b1 = b0
        b0 = x*b1 - b2 + vals[i]

    return 0.5*(b0 - b2)


def _i0_1(x):
    return exp(x) * _chbevl(x/2.0-2, _i0A)


def _i0_2(x):
    return exp(x) * _chbevl(32.0/x - 2.0, _i0B) / sqrt(x)


def i0(x):
    """
    Modified Bessel function of the first kind, order 0.

    Usually denoted :math:`I_0`.  This function does broadcast, but will *not*
    "up-cast" int dtype arguments unless accompanied by at least one float or
    complex dtype argument (see Raises below).

    Parameters
    ----------
    x : array_like, dtype float or complex
        Argument of the Bessel function.

    Returns
    -------
    out : ndarray, shape = x.shape, dtype = x.dtype
        The modified Bessel function evaluated at each of the elements of `x`.

    Raises
    ------
    TypeError: array cannot be safely cast to required type
        If argument consists exclusively of int dtypes.

    See Also
    --------
    scipy.special.iv, scipy.special.ive

    Notes
    -----
    We use the algorithm published by Clenshaw [1]_ and referenced by
    Abramowitz and Stegun [2]_, for which the function domain is
    partitioned into the two intervals [0,8] and (8,inf), and Chebyshev
    polynomial expansions are employed in each interval. Relative error on
    the domain [0,30] using IEEE arithmetic is documented [3]_ as having a
    peak of 5.8e-16 with an rms of 1.4e-16 (n = 30000).

    References
    ----------
    .. [1] C. W. Clenshaw, "Chebyshev series for mathematical functions", in
           *National Physical Laboratory Mathematical Tables*, vol. 5, London:
           Her Majesty's Stationery Office, 1962.
    .. [2] M. Abramowitz and I. A. Stegun, *Handbook of Mathematical
           Functions*, 10th printing, New York: Dover, 1964, pp. 379.
           http://www.math.sfu.ca/~cbm/aands/page_379.htm
    .. [3] http://kobesearch.cpan.org/htdocs/Math-Cephes/Math/Cephes.html

    Examples
    --------
    >>> np.i0([0.])
    array(1.0)
    >>> np.i0([0., 1. + 2j])
    array([ 1.00000000+0.j        ,  0.18785373+0.64616944j])

    """
    x = atleast_1d(x).copy()
    y = empty_like(x)
    ind = (x < 0)
    x[ind] = -x[ind]
    ind = (x <= 8.0)
    y[ind] = _i0_1(x[ind])
    ind2 = ~ind
    y[ind2] = _i0_2(x[ind2])
    return y.squeeze()

## End of cephes code for i0


def kaiser(M, beta):
    """
    Return the Kaiser window.

    The Kaiser window is a taper formed by using a Bessel function.

    Parameters
    ----------
    M : int
        Number of points in the output window. If zero or less, an
        empty array is returned.
    beta : float
        Shape parameter for window.

    Returns
    -------
    out : array
        The window, with the maximum value normalized to one (the value
        one appears only if the number of samples is odd).

    See Also
    --------
    bartlett, blackman, hamming, hanning

    Notes
    -----
    The Kaiser window is defined as

    .. math::  w(n) = I_0\\left( \\beta \\sqrt{1-\\frac{4n^2}{(M-1)^2}}
               \\right)/I_0(\\beta)

    with

    .. math:: \\quad -\\frac{M-1}{2} \\leq n \\leq \\frac{M-1}{2},

    where :math:`I_0` is the modified zeroth-order Bessel function.

    The Kaiser was named for Jim Kaiser, who discovered a simple
    approximation to the DPSS window based on Bessel functions.  The Kaiser
    window is a very good approximation to the Digital Prolate Spheroidal
    Sequence, or Slepian window, which is the transform which maximizes the
    energy in the main lobe of the window relative to total energy.

    The Kaiser can approximate many other windows by varying the beta
    parameter.

    ====  =======================
    beta  Window shape
    ====  =======================
    0     Rectangular
    5     Similar to a Hamming
    6     Similar to a Hanning
    8.6   Similar to a Blackman
    ====  =======================

    A beta value of 14 is probably a good starting point. Note that as beta
    gets large, the window narrows, and so the number of samples needs to be
    large enough to sample the increasingly narrow spike, otherwise NaNs will
    get returned.

    Most references to the Kaiser window come from the signal processing
    literature, where it is used as one of many windowing functions for
    smoothing values.  It is also known as an apodization (which means
    "removing the foot", i.e. smoothing discontinuities at the beginning
    and end of the sampled signal) or tapering function.

    References
    ----------
    .. [1] J. F. Kaiser, "Digital Filters" - Ch 7 in "Systems analysis by
           digital computer", Editors: F.F. Kuo and J.F. Kaiser, p 218-285.
           John Wiley and Sons, New York, (1966).
    .. [2] E.R. Kanasewich, "Time Sequence Analysis in Geophysics", The
           University of Alberta Press, 1975, pp. 177-178.
    .. [3] Wikipedia, "Window function",
           http://en.wikipedia.org/wiki/Window_function

    Examples
    --------
    >>> np.kaiser(12, 14)
    array([  7.72686684e-06,   3.46009194e-03,   4.65200189e-02,
             2.29737120e-01,   5.99885316e-01,   9.45674898e-01,
             9.45674898e-01,   5.99885316e-01,   2.29737120e-01,
             4.65200189e-02,   3.46009194e-03,   7.72686684e-06])


    Plot the window and the frequency response:

    >>> from numpy.fft import fft, fftshift
    >>> window = np.kaiser(51, 14)
    >>> plt.plot(window)
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.title("Kaiser window")
    <matplotlib.text.Text object at 0x...>
    >>> plt.ylabel("Amplitude")
    <matplotlib.text.Text object at 0x...>
    >>> plt.xlabel("Sample")
    <matplotlib.text.Text object at 0x...>
    >>> plt.show()

    >>> plt.figure()
    <matplotlib.figure.Figure object at 0x...>
    >>> A = fft(window, 2048) / 25.5
    >>> mag = np.abs(fftshift(A))
    >>> freq = np.linspace(-0.5, 0.5, len(A))
    >>> response = 20 * np.log10(mag)
    >>> response = np.clip(response, -100, 100)
    >>> plt.plot(freq, response)
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.title("Frequency response of Kaiser window")
    <matplotlib.text.Text object at 0x...>
    >>> plt.ylabel("Magnitude [dB]")
    <matplotlib.text.Text object at 0x...>
    >>> plt.xlabel("Normalized frequency [cycles per sample]")
    <matplotlib.text.Text object at 0x...>
    >>> plt.axis('tight')
    (-0.5, 0.5, -100.0, ...)
    >>> plt.show()

    """
    from numpy.dual import i0
    if M == 1:
        return np.array([1.])
    n = arange(0, M)
    alpha = (M-1)/2.0
    return i0(beta * sqrt(1-((n-alpha)/alpha)**2.0))/i0(float(beta))


def sinc(x):
    """
    Return the sinc function.

    The sinc function is :math:`\\sin(\\pi x)/(\\pi x)`.

    Parameters
    ----------
    x : ndarray
        Array (possibly multi-dimensional) of values for which to to
        calculate ``sinc(x)``.

    Returns
    -------
    out : ndarray
        ``sinc(x)``, which has the same shape as the input.

    Notes
    -----
    ``sinc(0)`` is the limit value 1.

    The name sinc is short for "sine cardinal" or "sinus cardinalis".

    The sinc function is used in various signal processing applications,
    including in anti-aliasing, in the construction of a Lanczos resampling
    filter, and in interpolation.

    For bandlimited interpolation of discrete-time signals, the ideal
    interpolation kernel is proportional to the sinc function.

    References
    ----------
    .. [1] Weisstein, Eric W. "Sinc Function." From MathWorld--A Wolfram Web
           Resource. http://mathworld.wolfram.com/SincFunction.html
    .. [2] Wikipedia, "Sinc function",
           http://en.wikipedia.org/wiki/Sinc_function

    Examples
    --------
    >>> x = np.linspace(-4, 4, 41)
    >>> np.sinc(x)
    array([ -3.89804309e-17,  -4.92362781e-02,  -8.40918587e-02,
            -8.90384387e-02,  -5.84680802e-02,   3.89804309e-17,
             6.68206631e-02,   1.16434881e-01,   1.26137788e-01,
             8.50444803e-02,  -3.89804309e-17,  -1.03943254e-01,
            -1.89206682e-01,  -2.16236208e-01,  -1.55914881e-01,
             3.89804309e-17,   2.33872321e-01,   5.04551152e-01,
             7.56826729e-01,   9.35489284e-01,   1.00000000e+00,
             9.35489284e-01,   7.56826729e-01,   5.04551152e-01,
             2.33872321e-01,   3.89804309e-17,  -1.55914881e-01,
            -2.16236208e-01,  -1.89206682e-01,  -1.03943254e-01,
            -3.89804309e-17,   8.50444803e-02,   1.26137788e-01,
             1.16434881e-01,   6.68206631e-02,   3.89804309e-17,
            -5.84680802e-02,  -8.90384387e-02,  -8.40918587e-02,
            -4.92362781e-02,  -3.89804309e-17])

    >>> plt.plot(x, np.sinc(x))
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.title("Sinc Function")
    <matplotlib.text.Text object at 0x...>
    >>> plt.ylabel("Amplitude")
    <matplotlib.text.Text object at 0x...>
    >>> plt.xlabel("X")
    <matplotlib.text.Text object at 0x...>
    >>> plt.show()

    It works in 2-D as well:

    >>> x = np.linspace(-4, 4, 401)
    >>> xx = np.outer(x, x)
    >>> plt.imshow(np.sinc(xx))
    <matplotlib.image.AxesImage object at 0x...>

    """
    x = np.asanyarray(x)
    y = pi * where(x == 0, 1.0e-20, x)
    return sin(y)/y


def msort(a):
    """
    Return a copy of an array sorted along the first axis.

    Parameters
    ----------
    a : array_like
        Array to be sorted.

    Returns
    -------
    sorted_array : ndarray
        Array of the same type and shape as `a`.

    See Also
    --------
    sort

    Notes
    -----
    ``np.msort(a)`` is equivalent to  ``np.sort(a, axis=0)``.

    """
    b = array(a, subok=True, copy=True)
    b.sort(0)
    return b


def _ureduce(a, func, **kwargs):
    """
    Internal Function.
    Call `func` with `a` as first argument swapping the axes to use extended
    axis on functions that don't support it natively.

    Returns result and a.shape with axis dims set to 1.

    Parameters
    ----------
    a : array_like
        Input array or object that can be converted to an array.
    func : callable
        Reduction function Kapable of receiving an axis argument.
        It is is called with `a` as first argument followed by `kwargs`.
     kwargs : keyword arguments
        additional keyword arguments to pass to `func`.

    Returns
    -------
    result : tuple
        Result of func(a, **kwargs) and a.shape with axis dims set to 1
        which can be used to reshape the result to the same shape a ufunc with
        keepdims=True would produce.

    """
    a = np.asanyarray(a)
    axis = kwargs.get('axis', None)
    if axis is not None:
        keepdim = list(a.shape)
        nd = a.ndim
        try:
            axis = operator.index(axis)
            if axis >= nd or axis < -nd:
                raise IndexError("axis %d out of bounds (%d)" % (axis, a.ndim))
            keepdim[axis] = 1
        except TypeError:
            sax = set()
            for x in axis:
                if x >= nd or x < -nd:
                    raise IndexError("axis %d out of bounds (%d)" % (x, nd))
                if x in sax:
                    raise ValueError("duplicate value in axis")
                sax.add(x % nd)
                keepdim[x] = 1
            keep = sax.symmetric_difference(frozenset(range(nd)))
            nkeep = len(keep)
            # swap axis that should not be reduced to front
            for i, s in enumerate(sorted(keep)):
                a = a.swapaxes(i, s)
            # merge reduced axis
            a = a.reshape(a.shape[:nkeep] + (-1,))
            kwargs['axis'] = -1
    else:
        keepdim = [1] * a.ndim

    r = func(a, **kwargs)
    return r, keepdim


def median(a, axis=None, out=None, overwrite_input=False, keepdims=False):
    """
    Compute the median along the specified axis.

    Returns the median of the array elements.

    Parameters
    ----------
    a : array_like
        Input array or object that can be converted to an array.
    axis : int or sequence of int, optional
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
       median. This will save memory when you do not need to preserve the
       contents of the input array. Treat the input as undefined, but it
       will probably be fully or partially sorted. Default is False. Note
       that, if `overwrite_input` is True and the input is not already an
       ndarray, an error will be raised.
    keepdims : bool, optional
        If this is set to True, the axes which are reduced are left
        in the result as dimensions with size one. With this option,
        the result will broadcast correctly against the original `arr`.

        .. versionadded:: 1.9.0


    Returns
    -------
    median : ndarray
        A new array holding the result (unless `out` is specified, in which
        case that array is returned instead).  If the input contains
        integers, or floats of smaller precision than 64, then the output
        data-type is float64.  Otherwise, the output data-type is the same
        as that of the input.

    See Also
    --------
    mean, percentile

    Notes
    -----
    Given a vector V of length N, the median of V is the middle value of
    a sorted copy of V, ``V_sorted`` - i.e., ``V_sorted[(N-1)/2]``, when N is
    odd.  When N is even, it is the average of the two middle values of
    ``V_sorted``.

    Examples
    --------
    >>> a = np.array([[10, 7, 4], [3, 2, 1]])
    >>> a
    array([[10,  7,  4],
           [ 3,  2,  1]])
    >>> np.median(a)
    3.5
    >>> np.median(a, axis=0)
    array([ 6.5,  4.5,  2.5])
    >>> np.median(a, axis=1)
    array([ 7.,  2.])
    >>> m = np.median(a, axis=0)
    >>> out = np.zeros_like(m)
    >>> np.median(a, axis=0, out=m)
    array([ 6.5,  4.5,  2.5])
    >>> m
    array([ 6.5,  4.5,  2.5])
    >>> b = a.copy()
    >>> np.median(b, axis=1, overwrite_input=True)
    array([ 7.,  2.])
    >>> assert not np.all(a==b)
    >>> b = a.copy()
    >>> np.median(b, axis=None, overwrite_input=True)
    3.5
    >>> assert not np.all(a==b)

    """
    r, k = _ureduce(a, func=_median, axis=axis, out=out,
                    overwrite_input=overwrite_input)
    if keepdims:
        return r.reshape(k)
    else:
        return r

def _median(a, axis=None, out=None, overwrite_input=False):
    # can't be reasonably be implemented in terms of percentile as we have to
    # call mean to not break astropy
    a = np.asanyarray(a)
    if axis is not None and axis >= a.ndim:
        raise IndexError(
            "axis %d out of bounds (%d)" % (axis, a.ndim))

    if overwrite_input:
        if axis is None:
            part = a.ravel()
            sz = part.size
            if sz % 2 == 0:
                szh = sz // 2
                part.partition((szh - 1, szh))
            else:
                part.partition((sz - 1) // 2)
        else:
            sz = a.shape[axis]
            if sz % 2 == 0:
                szh = sz // 2
                a.partition((szh - 1, szh), axis=axis)
            else:
                a.partition((sz - 1) // 2, axis=axis)
            part = a
    else:
        if axis is None:
            sz = a.size
        else:
            sz = a.shape[axis]
        if sz % 2 == 0:
            part = partition(a, ((sz // 2) - 1, sz // 2), axis=axis)
        else:
            part = partition(a, (sz - 1) // 2, axis=axis)
    if part.shape == ():
        # make 0-D arrays work
        return part.item()
    if axis is None:
        axis = 0
    indexer = [slice(None)] * part.ndim
    index = part.shape[axis] // 2
    if part.shape[axis] % 2 == 1:
        # index with slice to allow mean (below) to work
        indexer[axis] = slice(index, index+1)
    else:
        indexer[axis] = slice(index-1, index+1)
    # Use mean in odd and even case to coerce data type
    # and check, use out array.
    return mean(part[indexer], axis=axis, out=out)


def percentile(a, q, axis=None, out=None,
               overwrite_input=False, interpolation='linear', keepdims=False):
    """
    Compute the qth percentile of the data along the specified axis.

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

        .. versionadded:: 1.9.0
    keepdims : bool, optional
        If this is set to True, the axes which are reduced are left
        in the result as dimensions with size one. With this option,
        the result will broadcast correctly against the original `arr`.

        .. versionadded:: 1.9.0

    Returns
    -------
    percentile : scalar or ndarray
        If a single percentile `q` is given and axis=None a scalar is
        returned.  If multiple percentiles `q` are given an array holding
        the result is returned. The results are listed in the first axis.
        (If `out` is specified, in which case that array is returned
        instead).  If the input contains integers, or floats of smaller
        precision than 64, then the output data-type is float64. Otherwise,
        the output data-type is the same as that of the input.

    See Also
    --------
    mean, median

    Notes
    -----
    Given a vector V of length N, the q-th percentile of V is the q-th ranked
    value in a sorted copy of V.  The values and distances of the two
    nearest neighbors as well as the `interpolation` parameter will
    determine the percentile if the normalized ranking does not match q
    exactly. This function is the same as the median if ``q=50``, the same
    as the minimum if ``q=0`` and the same as the maximum if ``q=100``.

    Examples
    --------
    >>> a = np.array([[10, 7, 4], [3, 2, 1]])
    >>> a
    array([[10,  7,  4],
           [ 3,  2,  1]])
    >>> np.percentile(a, 50)
    array([ 3.5])
    >>> np.percentile(a, 50, axis=0)
    array([[ 6.5,  4.5,  2.5]])
    >>> np.percentile(a, 50, axis=1)
    array([[ 7.],
           [ 2.]])

    >>> m = np.percentile(a, 50, axis=0)
    >>> out = np.zeros_like(m)
    >>> np.percentile(a, 50, axis=0, out=m)
    array([[ 6.5,  4.5,  2.5]])
    >>> m
    array([[ 6.5,  4.5,  2.5]])

    >>> b = a.copy()
    >>> np.percentile(b, 50, axis=1, overwrite_input=True)
    array([[ 7.],
           [ 2.]])
    >>> assert not np.all(a==b)
    >>> b = a.copy()
    >>> np.percentile(b, 50, axis=None, overwrite_input=True)
    array([ 3.5])

    """
    q = array(q, dtype=np.float64, copy=True)
    r, k = _ureduce(a, func=_percentile, q=q, axis=axis, out=out,
                    overwrite_input=overwrite_input,
                    interpolation=interpolation)
    if keepdims:
        if q.ndim == 0:
            return r.reshape(k)
        else:
            return r.reshape([len(q)] + k)
    else:
        return r


def _percentile(a, q, axis=None, out=None,
                overwrite_input=False, interpolation='linear', keepdims=False):
    a = asarray(a)
    if q.ndim == 0:
        # Do not allow 0-d arrays because following code fails for scalar
        zerod = True
        q = q[None]
    else:
        zerod = False

    # avoid expensive reductions, relevant for arrays with < O(1000) elements
    if q.size < 10:
        for i in range(q.size):
            if q[i] < 0. or q[i] > 100.:
                raise ValueError("Percentiles must be in the range [0,100]")
            q[i] /= 100.
    else:
        # faster than any()
        if np.count_nonzero(q < 0.) or np.count_nonzero(q > 100.):
            raise ValueError("Percentiles must be in the range [0,100]")
        q /= 100.

    # prepare a for partioning
    if overwrite_input:
        if axis is None:
            ap = a.ravel()
        else:
            ap = a
    else:
        if axis is None:
            ap = a.flatten()
        else:
            ap = a.copy()

    if axis is None:
        axis = 0

    Nx = ap.shape[axis]
    indices = q * (Nx - 1)

    # round fractional indices according to interpolation method
    if interpolation == 'lower':
        indices = floor(indices).astype(intp)
    elif interpolation == 'higher':
        indices = ceil(indices).astype(intp)
    elif interpolation == 'midpoint':
        indices = floor(indices) + 0.5
    elif interpolation == 'nearest':
        indices = around(indices).astype(intp)
    elif interpolation == 'linear':
        pass  # keep index as fraction and interpolate
    else:
        raise ValueError(
            "interpolation can only be 'linear', 'lower' 'higher', "
            "'midpoint', or 'nearest'")

    if indices.dtype == intp:  # take the points along axis
        ap.partition(indices, axis=axis)
        # ensure axis with qth is first
        ap = np.rollaxis(ap, axis, 0)
        axis = 0

        if zerod:
            indices = indices[0]
        r = take(ap, indices, axis=axis, out=out)
    else:  # weight the points above and below the indices
        indices_below = floor(indices).astype(intp)
        indices_above = indices_below + 1
        indices_above[indices_above > Nx - 1] = Nx - 1

        weights_above = indices - indices_below
        weights_below = 1.0 - weights_above

        weights_shape = [1, ] * ap.ndim
        weights_shape[axis] = len(indices)
        weights_below.shape = weights_shape
        weights_above.shape = weights_shape

        ap.partition(concatenate((indices_below, indices_above)), axis=axis)
        x1 = take(ap, indices_below, axis=axis) * weights_below
        x2 = take(ap, indices_above, axis=axis) * weights_above

        # ensure axis with qth is first
        x1 = np.rollaxis(x1, axis, 0)
        x2 = np.rollaxis(x2, axis, 0)

        if zerod:
            x1 = x1.squeeze(0)
            x2 = x2.squeeze(0)

        if out is not None:
            r = add(x1, x2, out=out)
        else:
            r = add(x1, x2)

    return r


def trapz(y, x=None, dx=1.0, axis=-1):
    """
    Integrate along the given axis using the composite trapezoidal rule.

    Integrate `y` (`x`) along given axis.

    Parameters
    ----------
    y : array_like
        Input array to integrate.
    x : array_like, optional
        If `x` is None, then spacing between all `y` elements is `dx`.
    dx : scalar, optional
        If `x` is None, spacing given by `dx` is assumed. Default is 1.
    axis : int, optional
        Specify the axis.

    Returns
    -------
    trapz : float
        Definite integral as approximated by trapezoidal rule.

    See Also
    --------
    sum, cumsum

    Notes
    -----
    Image [2]_ illustrates trapezoidal rule -- y-axis locations of points
    will be taken from `y` array, by default x-axis distances between
    points will be 1.0, alternatively they can be provided with `x` array
    or with `dx` scalar.  Return value will be equal to combined area under
    the red lines.


    References
    ----------
    .. [1] Wikipedia page: http://en.wikipedia.org/wiki/Trapezoidal_rule

    .. [2] Illustration image:
           http://en.wikipedia.org/wiki/File:Composite_trapezoidal_rule_illustration.png

    Examples
    --------
    >>> np.trapz([1,2,3])
    4.0
    >>> np.trapz([1,2,3], x=[4,6,8])
    8.0
    >>> np.trapz([1,2,3], dx=2)
    8.0
    >>> a = np.arange(6).reshape(2, 3)
    >>> a
    array([[0, 1, 2],
           [3, 4, 5]])
    >>> np.trapz(a, axis=0)
    array([ 1.5,  2.5,  3.5])
    >>> np.trapz(a, axis=1)
    array([ 2.,  8.])

    """
    y = asanyarray(y)
    if x is None:
        d = dx
    else:
        x = asanyarray(x)
        if x.ndim == 1:
            d = diff(x)
            # reshape to correct shape
            shape = [1]*y.ndim
            shape[axis] = d.shape[0]
            d = d.reshape(shape)
        else:
            d = diff(x, axis=axis)
    nd = len(y.shape)
    slice1 = [slice(None)]*nd
    slice2 = [slice(None)]*nd
    slice1[axis] = slice(1, None)
    slice2[axis] = slice(None, -1)
    try:
        ret = (d * (y[slice1] + y[slice2]) / 2.0).sum(axis)
    except ValueError:
        # Operations didn't work, cast to ndarray
        d = np.asarray(d)
        y = np.asarray(y)
        ret = add.reduce(d * (y[slice1]+y[slice2])/2.0, axis)
    return ret


#always succeed
def add_newdoc(place, obj, doc):
    """Adds documentation to obj which is in module place.

    If doc is a string add it to obj as a docstring

    If doc is a tuple, then the first element is interpreted as
       an attribute of obj and the second as the docstring
          (method, docstring)

    If doc is a list, then each element of the list should be a
       sequence of length two --> [(method1, docstring1),
       (method2, docstring2), ...]

    This routine never raises an error.

    This routine cannot modify read-only docstrings, as appear
    in new-style classes or built-in functions. Because this
    routine never raises an error the caller must check manually
    that the docstrings were changed.
       """
    try:
        new = getattr(__import__(place, globals(), {}, [obj]), obj)
        if isinstance(doc, str):
            add_docstring(new, doc.strip())
        elif isinstance(doc, tuple):
            add_docstring(getattr(new, doc[0]), doc[1].strip())
        elif isinstance(doc, list):
            for val in doc:
                add_docstring(getattr(new, val[0]), val[1].strip())
    except:
        pass


# Based on scitools meshgrid
def meshgrid(*xi, **kwargs):
    """
    Return coordinate matrices from coordinate vectors.

    Make N-D coordinate arrays for vectorized evaluations of
    N-D scalar/vector fields over N-D grids, given
    one-dimensional coordinate arrays x1, x2,..., xn.

    .. versionchanged:: 1.9
       1-D and 0-D cases are allowed.

    Parameters
    ----------
    x1, x2,..., xn : array_like
        1-D arrays representing the coordinates of a grid.
    indexing : {'xy', 'ij'}, optional
        Cartesian ('xy', default) or matrix ('ij') indexing of output.
        See Notes for more details.

        .. versionadded:: 1.7.0
    sparse : bool, optional
        If True a sparse grid is returned in order to conserve memory.
        Default is False.

        .. versionadded:: 1.7.0
    copy : bool, optional
        If False, a view into the original arrays are returned in order to
        conserve memory.  Default is True.  Please note that
        ``sparse=False, copy=False`` will likely return non-contiguous
        arrays.  Furthermore, more than one element of a broadcast array
        may refer to a single memory location.  If you need to write to the
        arrays, make copies first.

        .. versionadded:: 1.7.0

    Returns
    -------
    X1, X2,..., XN : ndarray
        For vectors `x1`, `x2`,..., 'xn' with lengths ``Ni=len(xi)`` ,
        return ``(N1, N2, N3,...Nn)`` shaped arrays if indexing='ij'
        or ``(N2, N1, N3,...Nn)`` shaped arrays if indexing='xy'
        with the elements of `xi` repeated to fill the matrix along
        the first dimension for `x1`, the second for `x2` and so on.

    Notes
    -----
    This function supports both indexing conventions through the indexing
    keyword argument.  Giving the string 'ij' returns a meshgrid with
    matrix indexing, while 'xy' returns a meshgrid with Cartesian indexing.
    In the 2-D case with inputs of length M and N, the outputs are of shape
    (N, M) for 'xy' indexing and (M, N) for 'ij' indexing.  In the 3-D case
    with inputs of length M, N and P, outputs are of shape (N, M, P) for
    'xy' indexing and (M, N, P) for 'ij' indexing.  The difference is
    illustrated by the following code snippet::

        xv, yv = meshgrid(x, y, sparse=False, indexing='ij')
        for i in range(nx):
            for j in range(ny):
                # treat xv[i,j], yv[i,j]

        xv, yv = meshgrid(x, y, sparse=False, indexing='xy')
        for i in range(nx):
            for j in range(ny):
                # treat xv[j,i], yv[j,i]

    In the 1-D and 0-D case, the indexing and sparse keywords have no effect.

    See Also
    --------
    index_tricks.mgrid : Construct a multi-dimensional "meshgrid"
                     using indexing notation.
    index_tricks.ogrid : Construct an open multi-dimensional "meshgrid"
                     using indexing notation.

    Examples
    --------
    >>> nx, ny = (3, 2)
    >>> x = np.linspace(0, 1, nx)
    >>> y = np.linspace(0, 1, ny)
    >>> xv, yv = meshgrid(x, y)
    >>> xv
    array([[ 0. ,  0.5,  1. ],
           [ 0. ,  0.5,  1. ]])
    >>> yv
    array([[ 0.,  0.,  0.],
           [ 1.,  1.,  1.]])
    >>> xv, yv = meshgrid(x, y, sparse=True)  # make sparse output arrays
    >>> xv
    array([[ 0. ,  0.5,  1. ]])
    >>> yv
    array([[ 0.],
           [ 1.]])

    `meshgrid` is very useful to evaluate functions on a grid.

    >>> x = np.arange(-5, 5, 0.1)
    >>> y = np.arange(-5, 5, 0.1)
    >>> xx, yy = meshgrid(x, y, sparse=True)
    >>> z = np.sin(xx**2 + yy**2) / (xx**2 + yy**2)
    >>> h = plt.contourf(x,y,z)

    """
    ndim = len(xi)

    copy_ = kwargs.pop('copy', True)
    sparse = kwargs.pop('sparse', False)
    indexing = kwargs.pop('indexing', 'xy')

    if kwargs:
        raise TypeError("meshgrid() got an unexpected keyword argument '%s'"
                        % (list(kwargs)[0],))

    if indexing not in ['xy', 'ij']:
        raise ValueError(
            "Valid values for `indexing` are 'xy' and 'ij'.")

    s0 = (1,) * ndim
    output = [np.asanyarray(x).reshape(s0[:i] + (-1,) + s0[i + 1::])
              for i, x in enumerate(xi)]

    shape = [x.size for x in output]

    if indexing == 'xy' and ndim > 1:
        # switch first and second axis
        output[0].shape = (1, -1) + (1,)*(ndim - 2)
        output[1].shape = (-1, 1) + (1,)*(ndim - 2)
        shape[0], shape[1] = shape[1], shape[0]

    if sparse:
        if copy_:
            return [x.copy() for x in output]
        else:
            return output
    else:
        # Return the full N-D matrix (not only the 1-D vector)
        if copy_:
            mult_fact = np.ones(shape, dtype=int)
            return [x * mult_fact for x in output]
        else:
            return np.broadcast_arrays(*output)


def delete(arr, obj, axis=None):
    """
    Return a new array with sub-arrays along an axis deleted. For a one
    dimensional array, this returns those entries not returned by
    `arr[obj]`.

    Parameters
    ----------
    arr : array_like
      Input array.
    obj : slice, int or array of ints
      Indicate which sub-arrays to remove.
    axis : int, optional
      The axis along which to delete the subarray defined by `obj`.
      If `axis` is None, `obj` is applied to the flattened array.

    Returns
    -------
    out : ndarray
        A copy of `arr` with the elements specified by `obj` removed. Note
        that `delete` does not occur in-place. If `axis` is None, `out` is
        a flattened array.

    See Also
    --------
    insert : Insert elements into an array.
    append : Append elements at the end of an array.

    Notes
    -----
    Often it is preferable to use a boolean mask. For example:

    >>> mask = np.ones(len(arr), dtype=bool)
    >>> mask[[0,2,4]] = False
    >>> result = arr[mask,...]

    Is equivalent to `np.delete(arr, [0,2,4], axis=0)`, but allows further
    use of `mask`.

    Examples
    --------
    >>> arr = np.array([[1,2,3,4], [5,6,7,8], [9,10,11,12]])
    >>> arr
    array([[ 1,  2,  3,  4],
           [ 5,  6,  7,  8],
           [ 9, 10, 11, 12]])
    >>> np.delete(arr, 1, 0)
    array([[ 1,  2,  3,  4],
           [ 9, 10, 11, 12]])

    >>> np.delete(arr, np.s_[::2], 1)
    array([[ 2,  4],
           [ 6,  8],
           [10, 12]])
    >>> np.delete(arr, [1,3,5], None)
    array([ 1,  3,  5,  7,  8,  9, 10, 11, 12])

    """
    wrap = None
    if type(arr) is not ndarray:
        try:
            wrap = arr.__array_wrap__
        except AttributeError:
            pass

    arr = asarray(arr)
    ndim = arr.ndim
    if axis is None:
        if ndim != 1:
            arr = arr.ravel()
        ndim = arr.ndim
        axis = ndim - 1
    if ndim == 0:
        warnings.warn(
            "in the future the special handling of scalars will be removed "
            "from delete and raise an error", DeprecationWarning)
        if wrap:
            return wrap(arr)
        else:
            return arr.copy()

    slobj = [slice(None)]*ndim
    N = arr.shape[axis]
    newshape = list(arr.shape)

    if isinstance(obj, slice):
        start, stop, step = obj.indices(N)
        xr = range(start, stop, step)
        numtodel = len(xr)

        if numtodel <= 0:
            if wrap:
                return wrap(arr.copy())
            else:
                return arr.copy()

        # Invert if step is negative:
        if step < 0:
            step = -step
            start = xr[-1]
            stop = xr[0] + 1

        newshape[axis] -= numtodel
        new = empty(newshape, arr.dtype, arr.flags.fnc)
        # copy initial chunk
        if start == 0:
            pass
        else:
            slobj[axis] = slice(None, start)
            new[slobj] = arr[slobj]
        # copy end chunck
        if stop == N:
            pass
        else:
            slobj[axis] = slice(stop-numtodel, None)
            slobj2 = [slice(None)]*ndim
            slobj2[axis] = slice(stop, None)
            new[slobj] = arr[slobj2]
        # copy middle pieces
        if step == 1:
            pass
        else:  # use array indexing.
            keep = ones(stop-start, dtype=bool)
            keep[:stop-start:step] = False
            slobj[axis] = slice(start, stop-numtodel)
            slobj2 = [slice(None)]*ndim
            slobj2[axis] = slice(start, stop)
            arr = arr[slobj2]
            slobj2[axis] = keep
            new[slobj] = arr[slobj2]
        if wrap:
            return wrap(new)
        else:
            return new

    _obj = obj
    obj = np.asarray(obj)
    # After removing the special handling of booleans and out of
    # bounds values, the conversion to the array can be removed.
    if obj.dtype == bool:
        warnings.warn(
            "in the future insert will treat boolean arrays and array-likes "
            "as boolean index instead of casting it to integer", FutureWarning)
        obj = obj.astype(intp)
    if isinstance(_obj, (int, long, integer)):
        # optimization for a single value
        obj = obj.item()
        if (obj < -N or obj >= N):
            raise IndexError(
                "index %i is out of bounds for axis %i with "
                "size %i" % (obj, axis, N))
        if (obj < 0):
            obj += N
        newshape[axis] -= 1
        new = empty(newshape, arr.dtype, arr.flags.fnc)
        slobj[axis] = slice(None, obj)
        new[slobj] = arr[slobj]
        slobj[axis] = slice(obj, None)
        slobj2 = [slice(None)]*ndim
        slobj2[axis] = slice(obj+1, None)
        new[slobj] = arr[slobj2]
    else:
        if obj.size == 0 and not isinstance(_obj, np.ndarray):
            obj = obj.astype(intp)
        if not np.can_cast(obj, intp, 'same_kind'):
            # obj.size = 1 special case always failed and would just
            # give superfluous warnings.
            warnings.warn(
                "using a non-integer array as obj in delete will result in an "
                "error in the future", DeprecationWarning)
            obj = obj.astype(intp)
        keep = ones(N, dtype=bool)

        # Test if there are out of bound indices, this is deprecated
        inside_bounds = (obj < N) & (obj >= -N)
        if not inside_bounds.all():
            warnings.warn(
                "in the future out of bounds indices will raise an error "
                "instead of being ignored by `numpy.delete`.",
                DeprecationWarning)
            obj = obj[inside_bounds]
        positive_indices = obj >= 0
        if not positive_indices.all():
            warnings.warn(
                "in the future negative indices will not be ignored by "
                "`numpy.delete`.", FutureWarning)
            obj = obj[positive_indices]

        keep[obj, ] = False
        slobj[axis] = keep
        new = arr[slobj]

    if wrap:
        return wrap(new)
    else:
        return new


def insert(arr, obj, values, axis=None):
    """
    Insert values along the given axis before the given indices.

    Parameters
    ----------
    arr : array_like
        Input array.
    obj : int, slice or sequence of ints
        Object that defines the index or indices before which `values` is
        inserted.

        .. versionadded:: 1.8.0

        Support for multiple insertions when `obj` is a single scalar or a
        sequence with one element (similar to calling insert multiple
        times).
    values : array_like
        Values to insert into `arr`. If the type of `values` is different
        from that of `arr`, `values` is converted to the type of `arr`.
        `values` should be shaped so that ``arr[...,obj,...] = values``
        is legal.
    axis : int, optional
        Axis along which to insert `values`.  If `axis` is None then `arr`
        is flattened first.

    Returns
    -------
    out : ndarray
        A copy of `arr` with `values` inserted.  Note that `insert`
        does not occur in-place: a new array is returned. If
        `axis` is None, `out` is a flattened array.

    See Also
    --------
    append : Append elements at the end of an array.
    concatenate : Join a sequence of arrays together.
    delete : Delete elements from an array.

    Notes
    -----
    Note that for higher dimensional inserts `obj=0` behaves very different
    from `obj=[0]` just like `arr[:,0,:] = values` is different from
    `arr[:,[0],:] = values`.

    Examples
    --------
    >>> a = np.array([[1, 1], [2, 2], [3, 3]])
    >>> a
    array([[1, 1],
           [2, 2],
           [3, 3]])
    >>> np.insert(a, 1, 5)
    array([1, 5, 1, 2, 2, 3, 3])
    >>> np.insert(a, 1, 5, axis=1)
    array([[1, 5, 1],
           [2, 5, 2],
           [3, 5, 3]])

    Difference between sequence and scalars:
    >>> np.insert(a, [1], [[1],[2],[3]], axis=1)
    array([[1, 1, 1],
           [2, 2, 2],
           [3, 3, 3]])
    >>> np.array_equal(np.insert(a, 1, [1, 2, 3], axis=1),
    ...                np.insert(a, [1], [[1],[2],[3]], axis=1))
    True

    >>> b = a.flatten()
    >>> b
    array([1, 1, 2, 2, 3, 3])
    >>> np.insert(b, [2, 2], [5, 6])
    array([1, 1, 5, 6, 2, 2, 3, 3])

    >>> np.insert(b, slice(2, 4), [5, 6])
    array([1, 1, 5, 2, 6, 2, 3, 3])

    >>> np.insert(b, [2, 2], [7.13, False]) # type casting
    array([1, 1, 7, 0, 2, 2, 3, 3])

    >>> x = np.arange(8).reshape(2, 4)
    >>> idx = (1, 3)
    >>> np.insert(x, idx, 999, axis=1)
    array([[  0, 999,   1,   2, 999,   3],
           [  4, 999,   5,   6, 999,   7]])

    """
    wrap = None
    if type(arr) is not ndarray:
        try:
            wrap = arr.__array_wrap__
        except AttributeError:
            pass

    arr = asarray(arr)
    ndim = arr.ndim
    if axis is None:
        if ndim != 1:
            arr = arr.ravel()
        ndim = arr.ndim
        axis = ndim - 1
    else:
        if ndim > 0 and (axis < -ndim or axis >= ndim):
            raise IndexError(
                "axis %i is out of bounds for an array of "
                "dimension %i" % (axis, ndim))
        if (axis < 0):
            axis += ndim
    if (ndim == 0):
        warnings.warn(
            "in the future the special handling of scalars will be removed "
            "from insert and raise an error", DeprecationWarning)
        arr = arr.copy()
        arr[...] = values
        if wrap:
            return wrap(arr)
        else:
            return arr
    slobj = [slice(None)]*ndim
    N = arr.shape[axis]
    newshape = list(arr.shape)

    if isinstance(obj, slice):
        # turn it into a range object
        indices = arange(*obj.indices(N), **{'dtype': intp})
    else:
        # need to copy obj, because indices will be changed in-place
        indices = np.array(obj)
        if indices.dtype == bool:
            # See also delete
            warnings.warn(
                "in the future insert will treat boolean arrays and "
                "array-likes as a boolean index instead of casting it to "
                "integer", FutureWarning)
            indices = indices.astype(intp)
            # Code after warning period:
            #if obj.ndim != 1:
            #    raise ValueError('boolean array argument obj to insert '
            #                     'must be one dimensional')
            #indices = np.flatnonzero(obj)
        elif indices.ndim > 1:
            raise ValueError(
                "index array argument obj to insert must be one dimensional "
                "or scalar")
    if indices.size == 1:
        index = indices.item()
        if index < -N or index > N:
            raise IndexError(
                "index %i is out of bounds for axis %i with "
                "size %i" % (obj, axis, N))
        if (index < 0):
            index += N

        # There are some object array corner cases here, but we cannot avoid
        # that:
        values = array(values, copy=False, ndmin=arr.ndim, dtype=arr.dtype)
        if indices.ndim == 0:
            # broadcasting is very different here, since a[:,0,:] = ... behaves
            # very different from a[:,[0],:] = ...! This changes values so that
            # it works likes the second case. (here a[:,0:1,:])
            values = np.rollaxis(values, 0, (axis % values.ndim) + 1)
        numnew = values.shape[axis]
        newshape[axis] += numnew
        new = empty(newshape, arr.dtype, arr.flags.fnc)
        slobj[axis] = slice(None, index)
        new[slobj] = arr[slobj]
        slobj[axis] = slice(index, index+numnew)
        new[slobj] = values
        slobj[axis] = slice(index+numnew, None)
        slobj2 = [slice(None)] * ndim
        slobj2[axis] = slice(index, None)
        new[slobj] = arr[slobj2]
        if wrap:
            return wrap(new)
        return new
    elif indices.size == 0 and not isinstance(obj, np.ndarray):
        # Can safely cast the empty list to intp
        indices = indices.astype(intp)

    if not np.can_cast(indices, intp, 'same_kind'):
        warnings.warn(
            "using a non-integer array as obj in insert will result in an "
            "error in the future", DeprecationWarning)
        indices = indices.astype(intp)

    indices[indices < 0] += N

    numnew = len(indices)
    order = indices.argsort(kind='mergesort')   # stable sort
    indices[order] += np.arange(numnew)

    newshape[axis] += numnew
    old_mask = ones(newshape[axis], dtype=bool)
    old_mask[indices] = False

    new = empty(newshape, arr.dtype, arr.flags.fnc)
    slobj2 = [slice(None)]*ndim
    slobj[axis] = indices
    slobj2[axis] = old_mask
    new[slobj] = values
    new[slobj2] = arr

    if wrap:
        return wrap(new)
    return new


def append(arr, values, axis=None):
    """
    Append values to the end of an array.

    Parameters
    ----------
    arr : array_like
        Values are appended to a copy of this array.
    values : array_like
        These values are appended to a copy of `arr`.  It must be of the
        correct shape (the same shape as `arr`, excluding `axis`).  If
        `axis` is not specified, `values` can be any shape and will be
        flattened before use.
    axis : int, optional
        The axis along which `values` are appended.  If `axis` is not
        given, both `arr` and `values` are flattened before use.

    Returns
    -------
    append : ndarray
        A copy of `arr` with `values` appended to `axis`.  Note that
        `append` does not occur in-place: a new array is allocated and
        filled.  If `axis` is None, `out` is a flattened array.

    See Also
    --------
    insert : Insert elements into an array.
    delete : Delete elements from an array.

    Examples
    --------
    >>> np.append([1, 2, 3], [[4, 5, 6], [7, 8, 9]])
    array([1, 2, 3, 4, 5, 6, 7, 8, 9])

    When `axis` is specified, `values` must have the correct shape.

    >>> np.append([[1, 2, 3], [4, 5, 6]], [[7, 8, 9]], axis=0)
    array([[1, 2, 3],
           [4, 5, 6],
           [7, 8, 9]])
    >>> np.append([[1, 2, 3], [4, 5, 6]], [7, 8, 9], axis=0)
    Traceback (most recent call last):
    ...
    ValueError: arrays must have same number of dimensions

    """
    arr = asanyarray(arr)
    if axis is None:
        if arr.ndim != 1:
            arr = arr.ravel()
        values = ravel(values)
        axis = arr.ndim-1
    return concatenate((arr, values), axis=axis)
