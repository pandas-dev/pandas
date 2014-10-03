"""
Masked arrays add-ons.

A collection of utilities for `numpy.ma`.

:author: Pierre Gerard-Marchant
:contact: pierregm_at_uga_dot_edu
:version: $Id: extras.py 3473 2007-10-29 15:18:13Z jarrod.millman $

"""
from __future__ import division, absolute_import, print_function

__author__ = "Pierre GF Gerard-Marchant ($Author: jarrod.millman $)"
__version__ = '1.0'
__revision__ = "$Revision: 3473 $"
__date__ = '$Date: 2007-10-29 17:18:13 +0200 (Mon, 29 Oct 2007) $'

__all__ = ['apply_along_axis', 'apply_over_axes', 'atleast_1d', 'atleast_2d',
           'atleast_3d', 'average',
           'clump_masked', 'clump_unmasked', 'column_stack', 'compress_cols',
           'compress_rowcols', 'compress_rows', 'count_masked', 'corrcoef',
           'cov',
           'diagflat', 'dot', 'dstack',
           'ediff1d',
           'flatnotmasked_contiguous', 'flatnotmasked_edges',
           'hsplit', 'hstack',
           'in1d', 'intersect1d',
           'mask_cols', 'mask_rowcols', 'mask_rows', 'masked_all',
           'masked_all_like', 'median', 'mr_',
           'notmasked_contiguous', 'notmasked_edges',
           'polyfit',
           'row_stack',
           'setdiff1d', 'setxor1d',
           'unique', 'union1d',
           'vander', 'vstack',
           ]

import itertools
import warnings

from . import core as ma
from .core import MaskedArray, MAError, add, array, asarray, concatenate, count, \
    filled, getmask, getmaskarray, make_mask_descr, masked, masked_array, \
    mask_or, nomask, ones, sort, zeros
#from core import *

import numpy as np
from numpy import ndarray, array as nxarray
import numpy.core.umath as umath
from numpy.lib.index_tricks import AxisConcatenator
from numpy.linalg import lstsq


#...............................................................................
def issequence(seq):
    """Is seq a sequence (ndarray, list or tuple)?"""
    if isinstance(seq, (ndarray, tuple, list)):
        return True
    return False

def count_masked(arr, axis=None):
    """
    Count the number of masked elements along the given axis.

    Parameters
    ----------
    arr : array_like
        An array with (possibly) masked elements.
    axis : int, optional
        Axis along which to count. If None (default), a flattened
        version of the array is used.

    Returns
    -------
    count : int, ndarray
        The total number of masked elements (axis=None) or the number
        of masked elements along each slice of the given axis.

    See Also
    --------
    MaskedArray.count : Count non-masked elements.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = np.arange(9).reshape((3,3))
    >>> a = ma.array(a)
    >>> a[1, 0] = ma.masked
    >>> a[1, 2] = ma.masked
    >>> a[2, 1] = ma.masked
    >>> a
    masked_array(data =
     [[0 1 2]
     [-- 4 --]
     [6 -- 8]],
          mask =
     [[False False False]
     [ True False  True]
     [False  True False]],
          fill_value=999999)
    >>> ma.count_masked(a)
    3

    When the `axis` keyword is used an array is returned.

    >>> ma.count_masked(a, axis=0)
    array([1, 1, 1])
    >>> ma.count_masked(a, axis=1)
    array([0, 2, 1])

    """
    m = getmaskarray(arr)
    return m.sum(axis)

def masked_all(shape, dtype=float):
    """
    Empty masked array with all elements masked.

    Return an empty masked array of the given shape and dtype, where all the
    data are masked.

    Parameters
    ----------
    shape : tuple
        Shape of the required MaskedArray.
    dtype : dtype, optional
        Data type of the output.

    Returns
    -------
    a : MaskedArray
        A masked array with all data masked.

    See Also
    --------
    masked_all_like : Empty masked array modelled on an existing array.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> ma.masked_all((3, 3))
    masked_array(data =
     [[-- -- --]
     [-- -- --]
     [-- -- --]],
          mask =
     [[ True  True  True]
     [ True  True  True]
     [ True  True  True]],
          fill_value=1e+20)

    The `dtype` parameter defines the underlying data type.

    >>> a = ma.masked_all((3, 3))
    >>> a.dtype
    dtype('float64')
    >>> a = ma.masked_all((3, 3), dtype=np.int32)
    >>> a.dtype
    dtype('int32')

    """
    a = masked_array(np.empty(shape, dtype),
                     mask=np.ones(shape, make_mask_descr(dtype)))
    return a

def masked_all_like(arr):
    """
    Empty masked array with the properties of an existing array.

    Return an empty masked array of the same shape and dtype as
    the array `arr`, where all the data are masked.

    Parameters
    ----------
    arr : ndarray
        An array describing the shape and dtype of the required MaskedArray.

    Returns
    -------
    a : MaskedArray
        A masked array with all data masked.

    Raises
    ------
    AttributeError
        If `arr` doesn't have a shape attribute (i.e. not an ndarray)

    See Also
    --------
    masked_all : Empty masked array with all elements masked.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> arr = np.zeros((2, 3), dtype=np.float32)
    >>> arr
    array([[ 0.,  0.,  0.],
           [ 0.,  0.,  0.]], dtype=float32)
    >>> ma.masked_all_like(arr)
    masked_array(data =
     [[-- -- --]
     [-- -- --]],
          mask =
     [[ True  True  True]
     [ True  True  True]],
          fill_value=1e+20)

    The dtype of the masked array matches the dtype of `arr`.

    >>> arr.dtype
    dtype('float32')
    >>> ma.masked_all_like(arr).dtype
    dtype('float32')

    """
    a = np.empty_like(arr).view(MaskedArray)
    a._mask = np.ones(a.shape, dtype=make_mask_descr(a.dtype))
    return a


#####--------------------------------------------------------------------------
#---- --- Standard functions ---
#####--------------------------------------------------------------------------
class _fromnxfunction:
    """
    Defines a wrapper to adapt NumPy functions to masked arrays.


    An instance of `_fromnxfunction` can be called with the same parameters
    as the wrapped NumPy function. The docstring of `newfunc` is adapted from
    the wrapped function as well, see `getdoc`.

    Parameters
    ----------
    funcname : str
        The name of the function to be adapted. The function should be
        in the NumPy namespace (i.e. ``np.funcname``).

    """

    def __init__(self, funcname):
        self.__name__ = funcname
        self.__doc__ = self.getdoc()

    def getdoc(self):
        """
        Retrieve the docstring and signature from the function.

        The ``__doc__`` attribute of the function is used as the docstring for
        the new masked array version of the function. A note on application
        of the function to the mask is appended.

        .. warning::
          If the function docstring already contained a Notes section, the
          new docstring will have two Notes sections instead of appending a note
          to the existing section.

        Parameters
        ----------
        None

        """
        npfunc = getattr(np, self.__name__, None)
        doc = getattr(npfunc, '__doc__', None)
        if doc:
            sig = self.__name__ + ma.get_object_signature(npfunc)
            locdoc = "Notes\n-----\nThe function is applied to both the _data"\
                     " and the _mask, if any."
            return '\n'.join((sig, doc, locdoc))
        return


    def __call__(self, *args, **params):
        func = getattr(np, self.__name__)
        if len(args) == 1:
            x = args[0]
            if isinstance(x, ndarray):
                _d = func(x.__array__(), **params)
                _m = func(getmaskarray(x), **params)
                return masked_array(_d, mask=_m)
            elif isinstance(x, tuple) or isinstance(x, list):
                _d = func(tuple([np.asarray(a) for a in x]), **params)
                _m = func(tuple([getmaskarray(a) for a in x]), **params)
                return masked_array(_d, mask=_m)
        else:
            arrays = []
            args = list(args)
            while len(args) > 0 and issequence(args[0]):
                arrays.append(args.pop(0))
            res = []
            for x in arrays:
                _d = func(np.asarray(x), *args, **params)
                _m = func(getmaskarray(x), *args, **params)
                res.append(masked_array(_d, mask=_m))
            return res

atleast_1d = _fromnxfunction('atleast_1d')
atleast_2d = _fromnxfunction('atleast_2d')
atleast_3d = _fromnxfunction('atleast_3d')
#atleast_1d = np.atleast_1d
#atleast_2d = np.atleast_2d
#atleast_3d = np.atleast_3d

vstack = row_stack = _fromnxfunction('vstack')
hstack = _fromnxfunction('hstack')
column_stack = _fromnxfunction('column_stack')
dstack = _fromnxfunction('dstack')

hsplit = _fromnxfunction('hsplit')

diagflat = _fromnxfunction('diagflat')


#####--------------------------------------------------------------------------
#----
#####--------------------------------------------------------------------------
def flatten_inplace(seq):
    """Flatten a sequence in place."""
    k = 0
    while (k != len(seq)):
        while hasattr(seq[k], '__iter__'):
            seq[k:(k + 1)] = seq[k]
        k += 1
    return seq


def apply_along_axis(func1d, axis, arr, *args, **kwargs):
    """
    (This docstring should be overwritten)
    """
    arr = array(arr, copy=False, subok=True)
    nd = arr.ndim
    if axis < 0:
        axis += nd
    if (axis >= nd):
        raise ValueError("axis must be less than arr.ndim; axis=%d, rank=%d."
            % (axis, nd))
    ind = [0] * (nd - 1)
    i = np.zeros(nd, 'O')
    indlist = list(range(nd))
    indlist.remove(axis)
    i[axis] = slice(None, None)
    outshape = np.asarray(arr.shape).take(indlist)
    i.put(indlist, ind)
    j = i.copy()
    res = func1d(arr[tuple(i.tolist())], *args, **kwargs)
    #  if res is a number, then we have a smaller output array
    asscalar = np.isscalar(res)
    if not asscalar:
        try:
            len(res)
        except TypeError:
            asscalar = True
    # Note: we shouldn't set the dtype of the output from the first result...
    #...so we force the type to object, and build a list of dtypes
    #...we'll just take the largest, to avoid some downcasting
    dtypes = []
    if asscalar:
        dtypes.append(np.asarray(res).dtype)
        outarr = zeros(outshape, object)
        outarr[tuple(ind)] = res
        Ntot = np.product(outshape)
        k = 1
        while k < Ntot:
            # increment the index
            ind[-1] += 1
            n = -1
            while (ind[n] >= outshape[n]) and (n > (1 - nd)):
                ind[n - 1] += 1
                ind[n] = 0
                n -= 1
            i.put(indlist, ind)
            res = func1d(arr[tuple(i.tolist())], *args, **kwargs)
            outarr[tuple(ind)] = res
            dtypes.append(asarray(res).dtype)
            k += 1
    else:
        res = array(res, copy=False, subok=True)
        j = i.copy()
        j[axis] = ([slice(None, None)] * res.ndim)
        j.put(indlist, ind)
        Ntot = np.product(outshape)
        holdshape = outshape
        outshape = list(arr.shape)
        outshape[axis] = res.shape
        dtypes.append(asarray(res).dtype)
        outshape = flatten_inplace(outshape)
        outarr = zeros(outshape, object)
        outarr[tuple(flatten_inplace(j.tolist()))] = res
        k = 1
        while k < Ntot:
            # increment the index
            ind[-1] += 1
            n = -1
            while (ind[n] >= holdshape[n]) and (n > (1 - nd)):
                ind[n - 1] += 1
                ind[n] = 0
                n -= 1
            i.put(indlist, ind)
            j.put(indlist, ind)
            res = func1d(arr[tuple(i.tolist())], *args, **kwargs)
            outarr[tuple(flatten_inplace(j.tolist()))] = res
            dtypes.append(asarray(res).dtype)
            k += 1
    max_dtypes = np.dtype(np.asarray(dtypes).max())
    if not hasattr(arr, '_mask'):
        result = np.asarray(outarr, dtype=max_dtypes)
    else:
        result = asarray(outarr, dtype=max_dtypes)
        result.fill_value = ma.default_fill_value(result)
    return result
apply_along_axis.__doc__ = np.apply_along_axis.__doc__


def apply_over_axes(func, a, axes):
    """
    (This docstring will be overwritten)
    """
    val = asarray(a)
    N = a.ndim
    if array(axes).ndim == 0:
        axes = (axes,)
    for axis in axes:
        if axis < 0: axis = N + axis
        args = (val, axis)
        res = func(*args)
        if res.ndim == val.ndim:
            val = res
        else:
            res = ma.expand_dims(res, axis)
            if res.ndim == val.ndim:
                val = res
            else:
                raise ValueError("function is not returning "
                        "an array of the correct shape")
    return val
apply_over_axes.__doc__ = np.apply_over_axes.__doc__[
    :np.apply_over_axes.__doc__.find('Notes')].rstrip() + \
    """

    Examples
    --------
    >>> a = ma.arange(24).reshape(2,3,4)
    >>> a[:,0,1] = ma.masked
    >>> a[:,1,:] = ma.masked
    >>> print a
    [[[0 -- 2 3]
      [-- -- -- --]
      [8 9 10 11]]

     [[12 -- 14 15]
      [-- -- -- --]
      [20 21 22 23]]]
    >>> print ma.apply_over_axes(ma.sum, a, [0,2])
    [[[46]
      [--]
      [124]]]

    Tuple axis arguments to ufuncs are equivalent:

    >>> print ma.sum(a, axis=(0,2)).reshape((1,-1,1))
    [[[46]
      [--]
      [124]]]
"""


def average(a, axis=None, weights=None, returned=False):
    """
    Return the weighted average of array over the given axis.

    Parameters
    ----------
    a : array_like
        Data to be averaged.
        Masked entries are not taken into account in the computation.
    axis : int, optional
        Axis along which the average is computed. The default is to compute
        the average of the flattened array.
    weights : array_like, optional
        The importance that each element has in the computation of the average.
        The weights array can either be 1-D (in which case its length must be
        the size of `a` along the given axis) or of the same shape as `a`.
        If ``weights=None``, then all data in `a` are assumed to have a
        weight equal to one.   If `weights` is complex, the imaginary parts
        are ignored.
    returned : bool, optional
        Flag indicating whether a tuple ``(result, sum of weights)``
        should be returned as output (True), or just the result (False).
        Default is False.

    Returns
    -------
    average, [sum_of_weights] : (tuple of) scalar or MaskedArray
        The average along the specified axis. When returned is `True`,
        return a tuple with the average as the first element and the sum
        of the weights as the second element. The return type is `np.float64`
        if `a` is of integer type, otherwise it is of the same type as `a`.
        If returned, `sum_of_weights` is of the same type as `average`.

    Examples
    --------
    >>> a = np.ma.array([1., 2., 3., 4.], mask=[False, False, True, True])
    >>> np.ma.average(a, weights=[3, 1, 0, 0])
    1.25

    >>> x = np.ma.arange(6.).reshape(3, 2)
    >>> print x
    [[ 0.  1.]
     [ 2.  3.]
     [ 4.  5.]]
    >>> avg, sumweights = np.ma.average(x, axis=0, weights=[1, 2, 3],
    ...                                 returned=True)
    >>> print avg
    [2.66666666667 3.66666666667]

    """
    a = asarray(a)
    mask = a.mask
    ash = a.shape
    if ash == ():
        ash = (1,)
    if axis is None:
        if mask is nomask:
            if weights is None:
                n = a.sum(axis=None)
                d = float(a.size)
            else:
                w = filled(weights, 0.0).ravel()
                n = umath.add.reduce(a._data.ravel() * w)
                d = umath.add.reduce(w)
                del w
        else:
            if weights is None:
                n = a.filled(0).sum(axis=None)
                d = float(umath.add.reduce((~mask).ravel()))
            else:
                w = array(filled(weights, 0.0), float, mask=mask).ravel()
                n = add.reduce(a.ravel() * w)
                d = add.reduce(w)
                del w
    else:
        if mask is nomask:
            if weights is None:
                d = ash[axis] * 1.0
                n = add.reduce(a._data, axis)
            else:
                w = filled(weights, 0.0)
                wsh = w.shape
                if wsh == ():
                    wsh = (1,)
                if wsh == ash:
                    w = np.array(w, float, copy=0)
                    n = add.reduce(a * w, axis)
                    d = add.reduce(w, axis)
                    del w
                elif wsh == (ash[axis],):
                    ni = ash[axis]
                    r = [None] * len(ash)
                    r[axis] = slice(None, None, 1)
                    w = eval ("w[" + repr(tuple(r)) + "] * ones(ash, float)")
                    n = add.reduce(a * w, axis)
                    d = add.reduce(w, axis, dtype=float)
                    del w, r
                else:
                    raise ValueError('average: weights wrong shape.')
        else:
            if weights is None:
                n = add.reduce(a, axis)
                d = umath.add.reduce((~mask), axis=axis, dtype=float)
            else:
                w = filled(weights, 0.0)
                wsh = w.shape
                if wsh == ():
                    wsh = (1,)
                if wsh == ash:
                    w = array(w, dtype=float, mask=mask, copy=0)
                    n = add.reduce(a * w, axis)
                    d = add.reduce(w, axis, dtype=float)
                elif wsh == (ash[axis],):
                    ni = ash[axis]
                    r = [None] * len(ash)
                    r[axis] = slice(None, None, 1)
                    w = eval ("w[" + repr(tuple(r)) + \
                              "] * masked_array(ones(ash, float), mask)")
                    n = add.reduce(a * w, axis)
                    d = add.reduce(w, axis, dtype=float)
                else:
                    raise ValueError('average: weights wrong shape.')
                del w
    if n is masked or d is masked:
        return masked
    result = n / d
    del n

    if isinstance(result, MaskedArray):
        if ((axis is None) or (axis == 0 and a.ndim == 1)) and \
           (result.mask is nomask):
            result = result._data
        if returned:
            if not isinstance(d, MaskedArray):
                d = masked_array(d)
            if isinstance(d, ndarray) and (not d.shape == result.shape):
                d = ones(result.shape, dtype=float) * d
    if returned:
        return result, d
    else:
        return result


def median(a, axis=None, out=None, overwrite_input=False):
    """
    Compute the median along the specified axis.

    Returns the median of the array elements.

    Parameters
    ----------
    a : array_like
        Input array or object that can be converted to an array.
    axis : int, optional
        Axis along which the medians are computed. The default (None) is
        to compute the median along a flattened version of the array.
    out : ndarray, optional
        Alternative output array in which to place the result. It must
        have the same shape and buffer length as the expected output
        but the type will be cast if necessary.
    overwrite_input : bool, optional
        If True, then allow use of memory of input array (a) for
        calculations. The input array will be modified by the call to
        median. This will save memory when you do not need to preserve
        the contents of the input array. Treat the input as undefined,
        but it will probably be fully or partially sorted. Default is
        False. Note that, if `overwrite_input` is True, and the input
        is not already an `ndarray`, an error will be raised.

    Returns
    -------
    median : ndarray
        A new array holding the result is returned unless out is
        specified, in which case a reference to out is returned.
        Return data-type is `float64` for integers and floats smaller than
        `float64`, or the input data-type, otherwise.

    See Also
    --------
    mean

    Notes
    -----
    Given a vector ``V`` with ``N`` non masked values, the median of ``V``
    is the middle value of a sorted copy of ``V`` (``Vs``) - i.e.
    ``Vs[(N-1)/2]``, when ``N`` is odd, or ``{Vs[N/2 - 1] + Vs[N/2]}/2``
    when ``N`` is even.

    Examples
    --------
    >>> x = np.ma.array(np.arange(8), mask=[0]*4 + [1]*4)
    >>> np.ma.extras.median(x)
    1.5

    >>> x = np.ma.array(np.arange(10).reshape(2, 5), mask=[0]*6 + [1]*4)
    >>> np.ma.extras.median(x)
    2.5
    >>> np.ma.extras.median(x, axis=-1, overwrite_input=True)
    masked_array(data = [ 2.  5.],
                 mask = False,
           fill_value = 1e+20)

    """
    if not hasattr(a, 'mask') or np.count_nonzero(a.mask) == 0:
        return masked_array(np.median(a, axis=axis, out=out,
                                  overwrite_input=overwrite_input), copy=False)
    if overwrite_input:
        if axis is None:
            asorted = a.ravel()
            asorted.sort()
        else:
            a.sort(axis=axis)
            asorted = a
    else:
        asorted = sort(a, axis=axis)
    if axis is None:
        axis = 0
    elif axis < 0:
        axis += a.ndim

    counts = asorted.shape[axis] - (asorted.mask).sum(axis=axis)
    h = counts // 2
    # create indexing mesh grid for all but reduced axis
    axes_grid = [np.arange(x) for i, x in enumerate(asorted.shape)
                 if i != axis]
    ind = np.meshgrid(*axes_grid, sparse=True, indexing='ij')
    # insert indices of low and high median
    ind.insert(axis, h - 1)
    low = asorted[ind]
    ind[axis] = h
    high = asorted[ind]
    # duplicate high if odd number of elements so mean does nothing
    odd = counts % 2 == 1
    if asorted.ndim == 1:
        if odd:
            low = high
    else:
        low[odd] = high[odd]
    return np.ma.mean([low, high], axis=0, out=out)


#..............................................................................
def compress_rowcols(x, axis=None):
    """
    Suppress the rows and/or columns of a 2-D array that contain
    masked values.

    The suppression behavior is selected with the `axis` parameter.

    - If axis is None, both rows and columns are suppressed.
    - If axis is 0, only rows are suppressed.
    - If axis is 1 or -1, only columns are suppressed.

    Parameters
    ----------
    axis : int, optional
        Axis along which to perform the operation. Default is None.

    Returns
    -------
    compressed_array : ndarray
        The compressed array.

    Examples
    --------
    >>> x = np.ma.array(np.arange(9).reshape(3, 3), mask=[[1, 0, 0],
    ...                                                   [1, 0, 0],
    ...                                                   [0, 0, 0]])
    >>> x
    masked_array(data =
     [[-- 1 2]
     [-- 4 5]
     [6 7 8]],
                 mask =
     [[ True False False]
     [ True False False]
     [False False False]],
           fill_value = 999999)

    >>> np.ma.extras.compress_rowcols(x)
    array([[7, 8]])
    >>> np.ma.extras.compress_rowcols(x, 0)
    array([[6, 7, 8]])
    >>> np.ma.extras.compress_rowcols(x, 1)
    array([[1, 2],
           [4, 5],
           [7, 8]])

    """
    x = asarray(x)
    if x.ndim != 2:
        raise NotImplementedError("compress2d works for 2D arrays only.")
    m = getmask(x)
    # Nothing is masked: return x
    if m is nomask or not m.any():
        return x._data
    # All is masked: return empty
    if m.all():
        return nxarray([])
    # Builds a list of rows/columns indices
    (idxr, idxc) = (list(range(len(x))), list(range(x.shape[1])))
    masked = m.nonzero()
    if not axis:
        for i in np.unique(masked[0]):
            idxr.remove(i)
    if axis in [None, 1, -1]:
        for j in np.unique(masked[1]):
            idxc.remove(j)
    return x._data[idxr][:, idxc]

def compress_rows(a):
    """
    Suppress whole rows of a 2-D array that contain masked values.

    This is equivalent to ``np.ma.extras.compress_rowcols(a, 0)``, see
    `extras.compress_rowcols` for details.

    See Also
    --------
    extras.compress_rowcols

    """
    return compress_rowcols(a, 0)

def compress_cols(a):
    """
    Suppress whole columns of a 2-D array that contain masked values.

    This is equivalent to ``np.ma.extras.compress_rowcols(a, 1)``, see
    `extras.compress_rowcols` for details.

    See Also
    --------
    extras.compress_rowcols

    """
    return compress_rowcols(a, 1)

def mask_rowcols(a, axis=None):
    """
    Mask rows and/or columns of a 2D array that contain masked values.

    Mask whole rows and/or columns of a 2D array that contain
    masked values.  The masking behavior is selected using the
    `axis` parameter.

      - If `axis` is None, rows *and* columns are masked.
      - If `axis` is 0, only rows are masked.
      - If `axis` is 1 or -1, only columns are masked.

    Parameters
    ----------
    a : array_like, MaskedArray
        The array to mask.  If not a MaskedArray instance (or if no array
        elements are masked).  The result is a MaskedArray with `mask` set
        to `nomask` (False). Must be a 2D array.
    axis : int, optional
        Axis along which to perform the operation. If None, applies to a
        flattened version of the array.

    Returns
    -------
    a : MaskedArray
        A modified version of the input array, masked depending on the value
        of the `axis` parameter.

    Raises
    ------
    NotImplementedError
        If input array `a` is not 2D.

    See Also
    --------
    mask_rows : Mask rows of a 2D array that contain masked values.
    mask_cols : Mask cols of a 2D array that contain masked values.
    masked_where : Mask where a condition is met.

    Notes
    -----
    The input array's mask is modified by this function.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = np.zeros((3, 3), dtype=np.int)
    >>> a[1, 1] = 1
    >>> a
    array([[0, 0, 0],
           [0, 1, 0],
           [0, 0, 0]])
    >>> a = ma.masked_equal(a, 1)
    >>> a
    masked_array(data =
     [[0 0 0]
     [0 -- 0]
     [0 0 0]],
          mask =
     [[False False False]
     [False  True False]
     [False False False]],
          fill_value=999999)
    >>> ma.mask_rowcols(a)
    masked_array(data =
     [[0 -- 0]
     [-- -- --]
     [0 -- 0]],
          mask =
     [[False  True False]
     [ True  True  True]
     [False  True False]],
          fill_value=999999)

    """
    a = array(a, subok=False)
    if a.ndim != 2:
        raise NotImplementedError("mask_rowcols works for 2D arrays only.")
    m = getmask(a)
    # Nothing is masked: return a
    if m is nomask or not m.any():
        return a
    maskedval = m.nonzero()
    a._mask = a._mask.copy()
    if not axis:
        a[np.unique(maskedval[0])] = masked
    if axis in [None, 1, -1]:
        a[:, np.unique(maskedval[1])] = masked
    return a

def mask_rows(a, axis=None):
    """
    Mask rows of a 2D array that contain masked values.

    This function is a shortcut to ``mask_rowcols`` with `axis` equal to 0.

    See Also
    --------
    mask_rowcols : Mask rows and/or columns of a 2D array.
    masked_where : Mask where a condition is met.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = np.zeros((3, 3), dtype=np.int)
    >>> a[1, 1] = 1
    >>> a
    array([[0, 0, 0],
           [0, 1, 0],
           [0, 0, 0]])
    >>> a = ma.masked_equal(a, 1)
    >>> a
    masked_array(data =
     [[0 0 0]
     [0 -- 0]
     [0 0 0]],
          mask =
     [[False False False]
     [False  True False]
     [False False False]],
          fill_value=999999)
    >>> ma.mask_rows(a)
    masked_array(data =
     [[0 0 0]
     [-- -- --]
     [0 0 0]],
          mask =
     [[False False False]
     [ True  True  True]
     [False False False]],
          fill_value=999999)

    """
    return mask_rowcols(a, 0)

def mask_cols(a, axis=None):
    """
    Mask columns of a 2D array that contain masked values.

    This function is a shortcut to ``mask_rowcols`` with `axis` equal to 1.

    See Also
    --------
    mask_rowcols : Mask rows and/or columns of a 2D array.
    masked_where : Mask where a condition is met.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = np.zeros((3, 3), dtype=np.int)
    >>> a[1, 1] = 1
    >>> a
    array([[0, 0, 0],
           [0, 1, 0],
           [0, 0, 0]])
    >>> a = ma.masked_equal(a, 1)
    >>> a
    masked_array(data =
     [[0 0 0]
     [0 -- 0]
     [0 0 0]],
          mask =
     [[False False False]
     [False  True False]
     [False False False]],
          fill_value=999999)
    >>> ma.mask_cols(a)
    masked_array(data =
     [[0 -- 0]
     [0 -- 0]
     [0 -- 0]],
          mask =
     [[False  True False]
     [False  True False]
     [False  True False]],
          fill_value=999999)

    """
    return mask_rowcols(a, 1)


def dot(a, b, strict=False):
    """
    Return the dot product of two arrays.

    .. note::
      Works only with 2-D arrays at the moment.

    This function is the equivalent of `numpy.dot` that takes masked values
    into account, see `numpy.dot` for details.

    Parameters
    ----------
    a, b : ndarray
        Inputs arrays.
    strict : bool, optional
        Whether masked data are propagated (True) or set to 0 (False) for the
        computation. Default is False.
        Propagating the mask means that if a masked value appears in a row or
        column, the whole row or column is considered masked.

    See Also
    --------
    numpy.dot : Equivalent function for ndarrays.

    Examples
    --------
    >>> a = ma.array([[1, 2, 3], [4, 5, 6]], mask=[[1, 0, 0], [0, 0, 0]])
    >>> b = ma.array([[1, 2], [3, 4], [5, 6]], mask=[[1, 0], [0, 0], [0, 0]])
    >>> np.ma.dot(a, b)
    masked_array(data =
     [[21 26]
     [45 64]],
                 mask =
     [[False False]
     [False False]],
           fill_value = 999999)
    >>> np.ma.dot(a, b, strict=True)
    masked_array(data =
     [[-- --]
     [-- 64]],
                 mask =
     [[ True  True]
     [ True False]],
           fill_value = 999999)

    """
    #!!!: Works only with 2D arrays. There should be a way to get it to run with higher dimension
    if strict and (a.ndim == 2) and (b.ndim == 2):
        a = mask_rows(a)
        b = mask_cols(b)
    #
    d = np.dot(filled(a, 0), filled(b, 0))
    #
    am = (~getmaskarray(a))
    bm = (~getmaskarray(b))
    m = ~np.dot(am, bm)
    return masked_array(d, mask=m)

#####--------------------------------------------------------------------------
#---- --- arraysetops ---
#####--------------------------------------------------------------------------

def ediff1d(arr, to_end=None, to_begin=None):
    """
    Compute the differences between consecutive elements of an array.

    This function is the equivalent of `numpy.ediff1d` that takes masked
    values into account, see `numpy.ediff1d` for details.

    See Also
    --------
    numpy.ediff1d : Equivalent function for ndarrays.

    """
    arr = ma.asanyarray(arr).flat
    ed = arr[1:] - arr[:-1]
    arrays = [ed]
    #
    if to_begin is not None:
        arrays.insert(0, to_begin)
    if to_end is not None:
        arrays.append(to_end)
    #
    if len(arrays) != 1:
        # We'll save ourselves a copy of a potentially large array in the common
        # case where neither to_begin or to_end was given.
        ed = hstack(arrays)
    #
    return ed


def unique(ar1, return_index=False, return_inverse=False):
    """
    Finds the unique elements of an array.

    Masked values are considered the same element (masked). The output array
    is always a masked array. See `numpy.unique` for more details.

    See Also
    --------
    numpy.unique : Equivalent function for ndarrays.

    """
    output = np.unique(ar1,
                       return_index=return_index,
                       return_inverse=return_inverse)
    if isinstance(output, tuple):
        output = list(output)
        output[0] = output[0].view(MaskedArray)
        output = tuple(output)
    else:
        output = output.view(MaskedArray)
    return output


def intersect1d(ar1, ar2, assume_unique=False):
    """
    Returns the unique elements common to both arrays.

    Masked values are considered equal one to the other.
    The output is always a masked array.

    See `numpy.intersect1d` for more details.

    See Also
    --------
    numpy.intersect1d : Equivalent function for ndarrays.

    Examples
    --------
    >>> x = array([1, 3, 3, 3], mask=[0, 0, 0, 1])
    >>> y = array([3, 1, 1, 1], mask=[0, 0, 0, 1])
    >>> intersect1d(x, y)
    masked_array(data = [1 3 --],
                 mask = [False False  True],
           fill_value = 999999)

    """
    if assume_unique:
        aux = ma.concatenate((ar1, ar2))
    else:
        # Might be faster than unique( intersect1d( ar1, ar2 ) )?
        aux = ma.concatenate((unique(ar1), unique(ar2)))
    aux.sort()
    return aux[:-1][aux[1:] == aux[:-1]]


def setxor1d(ar1, ar2, assume_unique=False):
    """
    Set exclusive-or of 1-D arrays with unique elements.

    The output is always a masked array. See `numpy.setxor1d` for more details.

    See Also
    --------
    numpy.setxor1d : Equivalent function for ndarrays.

    """
    if not assume_unique:
        ar1 = unique(ar1)
        ar2 = unique(ar2)

    aux = ma.concatenate((ar1, ar2))
    if aux.size == 0:
        return aux
    aux.sort()
    auxf = aux.filled()
#    flag = ediff1d( aux, to_end = 1, to_begin = 1 ) == 0
    flag = ma.concatenate(([True], (auxf[1:] != auxf[:-1]), [True]))
#    flag2 = ediff1d( flag ) == 0
    flag2 = (flag[1:] == flag[:-1])
    return aux[flag2]

def in1d(ar1, ar2, assume_unique=False, invert=False):
    """
    Test whether each element of an array is also present in a second
    array.

    The output is always a masked array. See `numpy.in1d` for more details.

    See Also
    --------
    numpy.in1d : Equivalent function for ndarrays.

    Notes
    -----
    .. versionadded:: 1.4.0

    """
    if not assume_unique:
        ar1, rev_idx = unique(ar1, return_inverse=True)
        ar2 = unique(ar2)

    ar = ma.concatenate((ar1, ar2))
    # We need this to be a stable sort, so always use 'mergesort'
    # here. The values from the first array should always come before
    # the values from the second array.
    order = ar.argsort(kind='mergesort')
    sar = ar[order]
    if invert:
        bool_ar = (sar[1:] != sar[:-1])
    else:
        bool_ar = (sar[1:] == sar[:-1])
    flag = ma.concatenate((bool_ar, [invert]))
    indx = order.argsort(kind='mergesort')[:len(ar1)]

    if assume_unique:
        return flag[indx]
    else:
        return flag[indx][rev_idx]


def union1d(ar1, ar2):
    """
    Union of two arrays.

    The output is always a masked array. See `numpy.union1d` for more details.

    See also
    --------
    numpy.union1d : Equivalent function for ndarrays.

    """
    return unique(ma.concatenate((ar1, ar2)))


def setdiff1d(ar1, ar2, assume_unique=False):
    """
    Set difference of 1D arrays with unique elements.

    The output is always a masked array. See `numpy.setdiff1d` for more
    details.

    See Also
    --------
    numpy.setdiff1d : Equivalent function for ndarrays.

    Examples
    --------
    >>> x = np.ma.array([1, 2, 3, 4], mask=[0, 1, 0, 1])
    >>> np.ma.extras.setdiff1d(x, [1, 2])
    masked_array(data = [3 --],
                 mask = [False  True],
           fill_value = 999999)

    """
    if not assume_unique:
        ar1 = unique(ar1)
        ar2 = unique(ar2)
    aux = in1d(ar1, ar2, assume_unique=True)
    if aux.size == 0:
        return aux
    else:
        return ma.asarray(ar1)[aux == 0]


#####--------------------------------------------------------------------------
#---- --- Covariance ---
#####--------------------------------------------------------------------------




def _covhelper(x, y=None, rowvar=True, allow_masked=True):
    """
    Private function for the computation of covariance and correlation
    coefficients.

    """
    x = ma.array(x, ndmin=2, copy=True, dtype=float)
    xmask = ma.getmaskarray(x)
    # Quick exit if we can't process masked data
    if not allow_masked and xmask.any():
        raise ValueError("Cannot process masked data...")
    #
    if x.shape[0] == 1:
        rowvar = True
    # Make sure that rowvar is either 0 or 1
    rowvar = int(bool(rowvar))
    axis = 1 - rowvar
    if rowvar:
        tup = (slice(None), None)
    else:
        tup = (None, slice(None))
    #
    if y is None:
        xnotmask = np.logical_not(xmask).astype(int)
    else:
        y = array(y, copy=False, ndmin=2, dtype=float)
        ymask = ma.getmaskarray(y)
        if not allow_masked and ymask.any():
            raise ValueError("Cannot process masked data...")
        if xmask.any() or ymask.any():
            if y.shape == x.shape:
                # Define some common mask
                common_mask = np.logical_or(xmask, ymask)
                if common_mask is not nomask:
                    x.unshare_mask()
                    y.unshare_mask()
                    xmask = x._mask = y._mask = ymask = common_mask
        x = ma.concatenate((x, y), axis)
        xnotmask = np.logical_not(np.concatenate((xmask, ymask), axis)).astype(int)
    x -= x.mean(axis=rowvar)[tup]
    return (x, xnotmask, rowvar)


def cov(x, y=None, rowvar=True, bias=False, allow_masked=True, ddof=None):
    """
    Estimate the covariance matrix.

    Except for the handling of missing data this function does the same as
    `numpy.cov`. For more details and examples, see `numpy.cov`.

    By default, masked values are recognized as such. If `x` and `y` have the
    same shape, a common mask is allocated: if ``x[i,j]`` is masked, then
    ``y[i,j]`` will also be masked.
    Setting `allow_masked` to False will raise an exception if values are
    missing in either of the input arrays.

    Parameters
    ----------
    x : array_like
        A 1-D or 2-D array containing multiple variables and observations.
        Each row of `x` represents a variable, and each column a single
        observation of all those variables. Also see `rowvar` below.
    y : array_like, optional
        An additional set of variables and observations. `y` has the same
        form as `x`.
    rowvar : bool, optional
        If `rowvar` is True (default), then each row represents a
        variable, with observations in the columns. Otherwise, the relationship
        is transposed: each column represents a variable, while the rows
        contain observations.
    bias : bool, optional
        Default normalization (False) is by ``(N-1)``, where ``N`` is the
        number of observations given (unbiased estimate). If `bias` is True,
        then normalization is by ``N``. This keyword can be overridden by
        the keyword ``ddof`` in numpy versions >= 1.5.
    allow_masked : bool, optional
        If True, masked values are propagated pair-wise: if a value is masked
        in `x`, the corresponding value is masked in `y`.
        If False, raises a `ValueError` exception when some values are missing.
    ddof : {None, int}, optional
        .. versionadded:: 1.5
        If not ``None`` normalization is by ``(N - ddof)``, where ``N`` is
        the number of observations; this overrides the value implied by
        ``bias``. The default value is ``None``.


    Raises
    ------
    ValueError
        Raised if some values are missing and `allow_masked` is False.

    See Also
    --------
    numpy.cov

    """
    # Check inputs
    if ddof is not None and ddof != int(ddof):
        raise ValueError("ddof must be an integer")
    # Set up ddof
    if ddof is None:
        if bias:
            ddof = 0
        else:
            ddof = 1

    (x, xnotmask, rowvar) = _covhelper(x, y, rowvar, allow_masked)
    if not rowvar:
        fact = np.dot(xnotmask.T, xnotmask) * 1. - ddof
        result = (dot(x.T, x.conj(), strict=False) / fact).squeeze()
    else:
        fact = np.dot(xnotmask, xnotmask.T) * 1. - ddof
        result = (dot(x, x.T.conj(), strict=False) / fact).squeeze()
    return result


def corrcoef(x, y=None, rowvar=True, bias=False, allow_masked=True, ddof=None):
    """
    Return correlation coefficients of the input array.

    Except for the handling of missing data this function does the same as
    `numpy.corrcoef`. For more details and examples, see `numpy.corrcoef`.

    Parameters
    ----------
    x : array_like
        A 1-D or 2-D array containing multiple variables and observations.
        Each row of `x` represents a variable, and each column a single
        observation of all those variables. Also see `rowvar` below.
    y : array_like, optional
        An additional set of variables and observations. `y` has the same
        shape as `x`.
    rowvar : bool, optional
        If `rowvar` is True (default), then each row represents a
        variable, with observations in the columns. Otherwise, the relationship
        is transposed: each column represents a variable, while the rows
        contain observations.
    bias : bool, optional
        Default normalization (False) is by ``(N-1)``, where ``N`` is the
        number of observations given (unbiased estimate). If `bias` is 1,
        then normalization is by ``N``. This keyword can be overridden by
        the keyword ``ddof`` in numpy versions >= 1.5.
    allow_masked : bool, optional
        If True, masked values are propagated pair-wise: if a value is masked
        in `x`, the corresponding value is masked in `y`.
        If False, raises an exception.
    ddof : {None, int}, optional
        .. versionadded:: 1.5
        If not ``None`` normalization is by ``(N - ddof)``, where ``N`` is
        the number of observations; this overrides the value implied by
        ``bias``. The default value is ``None``.

    See Also
    --------
    numpy.corrcoef : Equivalent function in top-level NumPy module.
    cov : Estimate the covariance matrix.

    """
    # Check inputs
    if ddof is not None and ddof != int(ddof):
        raise ValueError("ddof must be an integer")
    # Set up ddof
    if ddof is None:
        if bias:
            ddof = 0
        else:
            ddof = 1

    # Get the data
    (x, xnotmask, rowvar) = _covhelper(x, y, rowvar, allow_masked)
    # Compute the covariance matrix
    if not rowvar:
        fact = np.dot(xnotmask.T, xnotmask) * 1. - ddof
        c = (dot(x.T, x.conj(), strict=False) / fact).squeeze()
    else:
        fact = np.dot(xnotmask, xnotmask.T) * 1. - ddof
        c = (dot(x, x.T.conj(), strict=False) / fact).squeeze()
    # Check whether we have a scalar
    try:
        diag = ma.diagonal(c)
    except ValueError:
        return 1
    #
    if xnotmask.all():
        _denom = ma.sqrt(ma.multiply.outer(diag, diag))
    else:
        _denom = diagflat(diag)
        n = x.shape[1 - rowvar]
        if rowvar:
            for i in range(n - 1):
                for j in range(i + 1, n):
                    _x = mask_cols(vstack((x[i], x[j]))).var(axis=1, ddof=ddof)
                    _denom[i, j] = _denom[j, i] = ma.sqrt(ma.multiply.reduce(_x))
        else:
            for i in range(n - 1):
                for j in range(i + 1, n):
                    _x = mask_cols(
                            vstack((x[:, i], x[:, j]))).var(axis=1, ddof=ddof)
                    _denom[i, j] = _denom[j, i] = ma.sqrt(ma.multiply.reduce(_x))
    return c / _denom

#####--------------------------------------------------------------------------
#---- --- Concatenation helpers ---
#####--------------------------------------------------------------------------

class MAxisConcatenator(AxisConcatenator):
    """
    Translate slice objects to concatenation along an axis.

    For documentation on usage, see `mr_class`.

    See Also
    --------
    mr_class

    """

    def __init__(self, axis=0):
        AxisConcatenator.__init__(self, axis, matrix=False)

    def __getitem__(self, key):
        if isinstance(key, str):
            raise MAError("Unavailable for masked array.")
        if not isinstance(key, tuple):
            key = (key,)
        objs = []
        scalars = []
        final_dtypedescr = None
        for k in range(len(key)):
            scalar = False
            if isinstance(key[k], slice):
                step = key[k].step
                start = key[k].start
                stop = key[k].stop
                if start is None:
                    start = 0
                if step is None:
                    step = 1
                if isinstance(step, complex):
                    size = int(abs(step))
                    newobj = np.linspace(start, stop, num=size)
                else:
                    newobj = np.arange(start, stop, step)
            elif isinstance(key[k], str):
                if (key[k] in 'rc'):
                    self.matrix = True
                    self.col = (key[k] == 'c')
                    continue
                try:
                    self.axis = int(key[k])
                    continue
                except (ValueError, TypeError):
                    raise ValueError("Unknown special directive")
            elif type(key[k]) in np.ScalarType:
                newobj = asarray([key[k]])
                scalars.append(k)
                scalar = True
            else:
                newobj = key[k]
            objs.append(newobj)
            if isinstance(newobj, ndarray) and not scalar:
                if final_dtypedescr is None:
                    final_dtypedescr = newobj.dtype
                elif newobj.dtype > final_dtypedescr:
                    final_dtypedescr = newobj.dtype
        if final_dtypedescr is not None:
            for k in scalars:
                objs[k] = objs[k].astype(final_dtypedescr)
        res = concatenate(tuple(objs), axis=self.axis)
        return self._retval(res)

class mr_class(MAxisConcatenator):
    """
    Translate slice objects to concatenation along the first axis.

    This is the masked array version of `lib.index_tricks.RClass`.

    See Also
    --------
    lib.index_tricks.RClass

    Examples
    --------
    >>> np.ma.mr_[np.ma.array([1,2,3]), 0, 0, np.ma.array([4,5,6])]
    array([1, 2, 3, 0, 0, 4, 5, 6])

    """
    def __init__(self):
        MAxisConcatenator.__init__(self, 0)

mr_ = mr_class()

#####--------------------------------------------------------------------------
#---- Find unmasked data ---
#####--------------------------------------------------------------------------

def flatnotmasked_edges(a):
    """
    Find the indices of the first and last unmasked values.

    Expects a 1-D `MaskedArray`, returns None if all values are masked.

    Parameters
    ----------
    arr : array_like
        Input 1-D `MaskedArray`

    Returns
    -------
    edges : ndarray or None
        The indices of first and last non-masked value in the array.
        Returns None if all values are masked.

    See Also
    --------
    flatnotmasked_contiguous, notmasked_contiguous, notmasked_edges,
    clump_masked, clump_unmasked

    Notes
    -----
    Only accepts 1-D arrays.

    Examples
    --------
    >>> a = np.ma.arange(10)
    >>> flatnotmasked_edges(a)
    [0,-1]

    >>> mask = (a < 3) | (a > 8) | (a == 5)
    >>> a[mask] = np.ma.masked
    >>> np.array(a[~a.mask])
    array([3, 4, 6, 7, 8])

    >>> flatnotmasked_edges(a)
    array([3, 8])

    >>> a[:] = np.ma.masked
    >>> print flatnotmasked_edges(ma)
    None

    """
    m = getmask(a)
    if m is nomask or not np.any(m):
        return np.array([0, a.size - 1])
    unmasked = np.flatnonzero(~m)
    if len(unmasked) > 0:
        return unmasked[[0, -1]]
    else:
        return None


def notmasked_edges(a, axis=None):
    """
    Find the indices of the first and last unmasked values along an axis.

    If all values are masked, return None.  Otherwise, return a list
    of two tuples, corresponding to the indices of the first and last
    unmasked values respectively.

    Parameters
    ----------
    a : array_like
        The input array.
    axis : int, optional
        Axis along which to perform the operation.
        If None (default), applies to a flattened version of the array.

    Returns
    -------
    edges : ndarray or list
        An array of start and end indexes if there are any masked data in
        the array. If there are no masked data in the array, `edges` is a
        list of the first and last index.

    See Also
    --------
    flatnotmasked_contiguous, flatnotmasked_edges, notmasked_contiguous,
    clump_masked, clump_unmasked

    Examples
    --------
    >>> a = np.arange(9).reshape((3, 3))
    >>> m = np.zeros_like(a)
    >>> m[1:, 1:] = 1

    >>> am = np.ma.array(a, mask=m)
    >>> np.array(am[~am.mask])
    array([0, 1, 2, 3, 6])

    >>> np.ma.extras.notmasked_edges(ma)
    array([0, 6])

    """
    a = asarray(a)
    if axis is None or a.ndim == 1:
        return flatnotmasked_edges(a)
    m = getmaskarray(a)
    idx = array(np.indices(a.shape), mask=np.asarray([m] * a.ndim))
    return [tuple([idx[i].min(axis).compressed() for i in range(a.ndim)]),
            tuple([idx[i].max(axis).compressed() for i in range(a.ndim)]), ]


def flatnotmasked_contiguous(a):
    """
    Find contiguous unmasked data in a masked array along the given axis.

    Parameters
    ----------
    a : narray
        The input array.

    Returns
    -------
    slice_list : list
        A sorted sequence of slices (start index, end index).

    See Also
    --------
    flatnotmasked_edges, notmasked_contiguous, notmasked_edges,
    clump_masked, clump_unmasked

    Notes
    -----
    Only accepts 2-D arrays at most.

    Examples
    --------
    >>> a = np.ma.arange(10)
    >>> np.ma.extras.flatnotmasked_contiguous(a)
    slice(0, 10, None)

    >>> mask = (a < 3) | (a > 8) | (a == 5)
    >>> a[mask] = np.ma.masked
    >>> np.array(a[~a.mask])
    array([3, 4, 6, 7, 8])

    >>> np.ma.extras.flatnotmasked_contiguous(a)
    [slice(3, 5, None), slice(6, 9, None)]
    >>> a[:] = np.ma.masked
    >>> print np.ma.extras.flatnotmasked_edges(a)
    None

    """
    m = getmask(a)
    if m is nomask:
        return slice(0, a.size, None)
    i = 0
    result = []
    for (k, g) in itertools.groupby(m.ravel()):
        n = len(list(g))
        if not k:
            result.append(slice(i, i + n))
        i += n
    return result or None

def notmasked_contiguous(a, axis=None):
    """
    Find contiguous unmasked data in a masked array along the given axis.

    Parameters
    ----------
    a : array_like
        The input array.
    axis : int, optional
        Axis along which to perform the operation.
        If None (default), applies to a flattened version of the array.

    Returns
    -------
    endpoints : list
        A list of slices (start and end indexes) of unmasked indexes
        in the array.

    See Also
    --------
    flatnotmasked_edges, flatnotmasked_contiguous, notmasked_edges,
    clump_masked, clump_unmasked

    Notes
    -----
    Only accepts 2-D arrays at most.

    Examples
    --------
    >>> a = np.arange(9).reshape((3, 3))
    >>> mask = np.zeros_like(a)
    >>> mask[1:, 1:] = 1

    >>> ma = np.ma.array(a, mask=mask)
    >>> np.array(ma[~ma.mask])
    array([0, 1, 2, 3, 6])

    >>> np.ma.extras.notmasked_contiguous(ma)
    [slice(0, 4, None), slice(6, 7, None)]

    """
    a = asarray(a)
    nd = a.ndim
    if nd > 2:
        raise NotImplementedError("Currently limited to atmost 2D array.")
    if axis is None or nd == 1:
        return flatnotmasked_contiguous(a)
    #
    result = []
    #
    other = (axis + 1) % 2
    idx = [0, 0]
    idx[axis] = slice(None, None)
    #
    for i in range(a.shape[other]):
        idx[other] = i
        result.append(flatnotmasked_contiguous(a[idx]) or None)
    return result


def _ezclump(mask):
    """
    Finds the clumps (groups of data with the same values) for a 1D bool array.

    Returns a series of slices.
    """
    #def clump_masked(a):
    if mask.ndim > 1:
        mask = mask.ravel()
    idx = (mask[1:] ^ mask[:-1]).nonzero()
    idx = idx[0] + 1
    slices = [slice(left, right)
              for (left, right) in zip(itertools.chain([0], idx),
                                       itertools.chain(idx, [len(mask)]),)]
    return slices


def clump_unmasked(a):
    """
    Return list of slices corresponding to the unmasked clumps of a 1-D array.
    (A "clump" is defined as a contiguous region of the array).

    Parameters
    ----------
    a : ndarray
        A one-dimensional masked array.

    Returns
    -------
    slices : list of slice
        The list of slices, one for each continuous region of unmasked
        elements in `a`.

    Notes
    -----
    .. versionadded:: 1.4.0

    See Also
    --------
    flatnotmasked_edges, flatnotmasked_contiguous, notmasked_edges,
    notmasked_contiguous, clump_masked

    Examples
    --------
    >>> a = np.ma.masked_array(np.arange(10))
    >>> a[[0, 1, 2, 6, 8, 9]] = np.ma.masked
    >>> np.ma.extras.clump_unmasked(a)
    [slice(3, 6, None), slice(7, 8, None)]

    """
    mask = getattr(a, '_mask', nomask)
    if mask is nomask:
        return [slice(0, a.size)]
    slices = _ezclump(mask)
    if a[0] is masked:
        result = slices[1::2]
    else:
        result = slices[::2]
    return result


def clump_masked(a):
    """
    Returns a list of slices corresponding to the masked clumps of a 1-D array.
    (A "clump" is defined as a contiguous region of the array).

    Parameters
    ----------
    a : ndarray
        A one-dimensional masked array.

    Returns
    -------
    slices : list of slice
        The list of slices, one for each continuous region of masked elements
        in `a`.

    Notes
    -----
    .. versionadded:: 1.4.0

    See Also
    --------
    flatnotmasked_edges, flatnotmasked_contiguous, notmasked_edges,
    notmasked_contiguous, clump_unmasked

    Examples
    --------
    >>> a = np.ma.masked_array(np.arange(10))
    >>> a[[0, 1, 2, 6, 8, 9]] = np.ma.masked
    >>> np.ma.extras.clump_masked(a)
    [slice(0, 3, None), slice(6, 7, None), slice(8, 10, None)]

    """
    mask = ma.getmask(a)
    if mask is nomask:
        return []
    slices = _ezclump(mask)
    if len(slices):
        if a[0] is masked:
            slices = slices[::2]
        else:
            slices = slices[1::2]
    return slices



#####--------------------------------------------------------------------------
#---- Polynomial fit ---
#####--------------------------------------------------------------------------

def vander(x, n=None):
    """
    Masked values in the input array result in rows of zeros.
    """
    _vander = np.vander(x, n)
    m = getmask(x)
    if m is not nomask:
        _vander[m] = 0
    return _vander
vander.__doc__ = ma.doc_note(np.vander.__doc__, vander.__doc__)


def polyfit(x, y, deg, rcond=None, full=False, w=None, cov=False):
    """
    Any masked values in x is propagated in y, and vice-versa.
    """
    x = asarray(x)
    y = asarray(y)

    m = getmask(x)
    if y.ndim == 1:
        m = mask_or(m, getmask(y))
    elif y.ndim == 2:
        my = getmask(mask_rows(y))
        if my is not nomask:
            m = mask_or(m, my[:, 0])
    else:
        raise TypeError("Expected a 1D or 2D array for y!")

    if w is not None:
        w = asarray(w)
        if w.ndim != 1:
            raise TypeError("expected a 1-d array for weights")
        if w.shape[0] != y.shape[0] :
            raise TypeError("expected w and y to have the same length")
        m = mask_or(m, getmask(w))

    if m is not nomask:
        if w is not None:
            w = ~m*w
        else:
            w = ~m

    return np.polyfit(x, y, deg, rcond, full, w, cov)

polyfit.__doc__ = ma.doc_note(np.polyfit.__doc__, polyfit.__doc__)

################################################################################
