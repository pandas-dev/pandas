from __future__ import division, absolute_import, print_function

import os
import sys
import warnings
import collections
from . import multiarray
from . import umath
from .umath import (invert, sin, UFUNC_BUFSIZE_DEFAULT, ERR_IGNORE,
                    ERR_WARN, ERR_RAISE, ERR_CALL, ERR_PRINT, ERR_LOG,
                    ERR_DEFAULT, PINF, NAN)
from . import numerictypes
from .numerictypes import longlong, intc, int_, float_, complex_, bool_

if sys.version_info[0] >= 3:
    import pickle
    basestring = str
else:
    import cPickle as pickle

loads = pickle.loads


__all__ = ['newaxis', 'ndarray', 'flatiter', 'nditer', 'nested_iters', 'ufunc',
           'arange', 'array', 'zeros', 'count_nonzero',
           'empty', 'broadcast', 'dtype', 'fromstring', 'fromfile',
           'frombuffer', 'int_asbuffer', 'where', 'argwhere', 'copyto',
           'concatenate', 'fastCopyAndTranspose', 'lexsort', 'set_numeric_ops',
           'can_cast', 'promote_types', 'min_scalar_type', 'result_type',
           'asarray', 'asanyarray', 'ascontiguousarray', 'asfortranarray',
           'isfortran', 'empty_like', 'zeros_like', 'ones_like',
           'correlate', 'convolve', 'inner', 'dot', 'einsum', 'outer', 'vdot',
           'alterdot', 'restoredot', 'roll', 'rollaxis', 'cross', 'tensordot',
           'array2string', 'get_printoptions', 'set_printoptions',
           'array_repr', 'array_str', 'set_string_function',
           'little_endian', 'require',
           'fromiter', 'array_equal', 'array_equiv',
           'indices', 'fromfunction', 'isclose',
           'load', 'loads', 'isscalar', 'binary_repr', 'base_repr',
           'ones', 'identity', 'allclose', 'compare_chararrays', 'putmask',
           'seterr', 'geterr', 'setbufsize', 'getbufsize',
           'seterrcall', 'geterrcall', 'errstate', 'flatnonzero',
           'Inf', 'inf', 'infty', 'Infinity',
           'nan', 'NaN', 'False_', 'True_', 'bitwise_not',
           'CLIP', 'RAISE', 'WRAP', 'MAXDIMS', 'BUFSIZE', 'ALLOW_THREADS',
           'ComplexWarning', 'may_share_memory', 'full', 'full_like']

if sys.version_info[0] < 3:
    __all__.extend(['getbuffer', 'newbuffer'])


class ComplexWarning(RuntimeWarning):
    """
    The warning raised when casting a complex dtype to a real dtype.

    As implemented, casting a complex number to a real discards its imaginary
    part, but this behavior may not be what the user actually wants.

    """
    pass

bitwise_not = invert

CLIP = multiarray.CLIP
WRAP = multiarray.WRAP
RAISE = multiarray.RAISE
MAXDIMS = multiarray.MAXDIMS
ALLOW_THREADS = multiarray.ALLOW_THREADS
BUFSIZE = multiarray.BUFSIZE

ndarray = multiarray.ndarray
flatiter = multiarray.flatiter
nditer = multiarray.nditer
nested_iters = multiarray.nested_iters
broadcast = multiarray.broadcast
dtype = multiarray.dtype
copyto = multiarray.copyto
ufunc = type(sin)


def zeros_like(a, dtype=None, order='K', subok=True):
    """
    Return an array of zeros with the same shape and type as a given array.

    Parameters
    ----------
    a : array_like
        The shape and data-type of `a` define these same attributes of
        the returned array.
    dtype : data-type, optional
        .. versionadded:: 1.6.0
        Overrides the data type of the result.
    order : {'C', 'F', 'A', or 'K'}, optional
        .. versionadded:: 1.6.0
        Overrides the memory layout of the result. 'C' means C-order,
        'F' means F-order, 'A' means 'F' if `a` is Fortran contiguous,
        'C' otherwise. 'K' means match the layout of `a` as closely
        as possible.
    subok : bool, optional.
        If True, then the newly created array will use the sub-class
        type of 'a', otherwise it will be a base-class array. Defaults
        to True.

    Returns
    -------
    out : ndarray
        Array of zeros with the same shape and type as `a`.

    See Also
    --------
    ones_like : Return an array of ones with shape and type of input.
    empty_like : Return an empty array with shape and type of input.
    zeros : Return a new array setting values to zero.
    ones : Return a new array setting values to one.
    empty : Return a new uninitialized array.

    Examples
    --------
    >>> x = np.arange(6)
    >>> x = x.reshape((2, 3))
    >>> x
    array([[0, 1, 2],
           [3, 4, 5]])
    >>> np.zeros_like(x)
    array([[0, 0, 0],
           [0, 0, 0]])

    >>> y = np.arange(3, dtype=np.float)
    >>> y
    array([ 0.,  1.,  2.])
    >>> np.zeros_like(y)
    array([ 0.,  0.,  0.])

    """
    res = empty_like(a, dtype=dtype, order=order, subok=subok)
    # needed instead of a 0 to get same result as zeros for for string dtypes
    z = zeros(1, dtype=res.dtype)
    multiarray.copyto(res, z, casting='unsafe')
    return res

def ones(shape, dtype=None, order='C'):
    """
    Return a new array of given shape and type, filled with ones.

    Parameters
    ----------
    shape : int or sequence of ints
        Shape of the new array, e.g., ``(2, 3)`` or ``2``.
    dtype : data-type, optional
        The desired data-type for the array, e.g., `numpy.int8`.  Default is
        `numpy.float64`.
    order : {'C', 'F'}, optional
        Whether to store multidimensional data in C- or Fortran-contiguous
        (row- or column-wise) order in memory.

    Returns
    -------
    out : ndarray
        Array of ones with the given shape, dtype, and order.

    See Also
    --------
    zeros, ones_like

    Examples
    --------
    >>> np.ones(5)
    array([ 1.,  1.,  1.,  1.,  1.])

    >>> np.ones((5,), dtype=np.int)
    array([1, 1, 1, 1, 1])

    >>> np.ones((2, 1))
    array([[ 1.],
           [ 1.]])

    >>> s = (2,2)
    >>> np.ones(s)
    array([[ 1.,  1.],
           [ 1.,  1.]])

    """
    a = empty(shape, dtype, order)
    multiarray.copyto(a, 1, casting='unsafe')
    return a

def ones_like(a, dtype=None, order='K', subok=True):
    """
    Return an array of ones with the same shape and type as a given array.

    Parameters
    ----------
    a : array_like
        The shape and data-type of `a` define these same attributes of
        the returned array.
    dtype : data-type, optional
        .. versionadded:: 1.6.0
        Overrides the data type of the result.
    order : {'C', 'F', 'A', or 'K'}, optional
        .. versionadded:: 1.6.0
        Overrides the memory layout of the result. 'C' means C-order,
        'F' means F-order, 'A' means 'F' if `a` is Fortran contiguous,
        'C' otherwise. 'K' means match the layout of `a` as closely
        as possible.
    subok : bool, optional.
        If True, then the newly created array will use the sub-class
        type of 'a', otherwise it will be a base-class array. Defaults
        to True.

    Returns
    -------
    out : ndarray
        Array of ones with the same shape and type as `a`.

    See Also
    --------
    zeros_like : Return an array of zeros with shape and type of input.
    empty_like : Return an empty array with shape and type of input.
    zeros : Return a new array setting values to zero.
    ones : Return a new array setting values to one.
    empty : Return a new uninitialized array.

    Examples
    --------
    >>> x = np.arange(6)
    >>> x = x.reshape((2, 3))
    >>> x
    array([[0, 1, 2],
           [3, 4, 5]])
    >>> np.ones_like(x)
    array([[1, 1, 1],
           [1, 1, 1]])

    >>> y = np.arange(3, dtype=np.float)
    >>> y
    array([ 0.,  1.,  2.])
    >>> np.ones_like(y)
    array([ 1.,  1.,  1.])

    """
    res = empty_like(a, dtype=dtype, order=order, subok=subok)
    multiarray.copyto(res, 1, casting='unsafe')
    return res

def full(shape, fill_value, dtype=None, order='C'):
    """
    Return a new array of given shape and type, filled with `fill_value`.

    Parameters
    ----------
    shape : int or sequence of ints
        Shape of the new array, e.g., ``(2, 3)`` or ``2``.
    fill_value : scalar
        Fill value.
    dtype : data-type, optional
        The desired data-type for the array, e.g., `numpy.int8`.  Default is
        is chosen as `np.array(fill_value).dtype`.
    order : {'C', 'F'}, optional
        Whether to store multidimensional data in C- or Fortran-contiguous
        (row- or column-wise) order in memory.

    Returns
    -------
    out : ndarray
        Array of `fill_value` with the given shape, dtype, and order.

    See Also
    --------
    zeros_like : Return an array of zeros with shape and type of input.
    ones_like : Return an array of ones with shape and type of input.
    empty_like : Return an empty array with shape and type of input.
    full_like : Fill an array with shape and type of input.
    zeros : Return a new array setting values to zero.
    ones : Return a new array setting values to one.
    empty : Return a new uninitialized array.

    Examples
    --------
    >>> np.full((2, 2), np.inf)
    array([[ inf,  inf],
           [ inf,  inf]])
    >>> np.full((2, 2), 10, dtype=np.int)
    array([[10, 10],
           [10, 10]])

    """
    a = empty(shape, dtype, order)
    multiarray.copyto(a, fill_value, casting='unsafe')
    return a

def full_like(a, fill_value, dtype=None, order='K', subok=True):
    """
    Return a full array with the same shape and type as a given array.

    Parameters
    ----------
    a : array_like
        The shape and data-type of `a` define these same attributes of
        the returned array.
    fill_value : scalar
        Fill value.
    dtype : data-type, optional
        Overrides the data type of the result.
    order : {'C', 'F', 'A', or 'K'}, optional
        Overrides the memory layout of the result. 'C' means C-order,
        'F' means F-order, 'A' means 'F' if `a` is Fortran contiguous,
        'C' otherwise. 'K' means match the layout of `a` as closely
        as possible.
    subok : bool, optional.
        If True, then the newly created array will use the sub-class
        type of 'a', otherwise it will be a base-class array. Defaults
        to True.

    Returns
    -------
    out : ndarray
        Array of `fill_value` with the same shape and type as `a`.

    See Also
    --------
    zeros_like : Return an array of zeros with shape and type of input.
    ones_like : Return an array of ones with shape and type of input.
    empty_like : Return an empty array with shape and type of input.
    zeros : Return a new array setting values to zero.
    ones : Return a new array setting values to one.
    empty : Return a new uninitialized array.
    full : Fill a new array.

    Examples
    --------
    >>> x = np.arange(6, dtype=np.int)
    >>> np.full_like(x, 1)
    array([1, 1, 1, 1, 1, 1])
    >>> np.full_like(x, 0.1)
    array([0, 0, 0, 0, 0, 0])
    >>> np.full_like(x, 0.1, dtype=np.double)
    array([ 0.1,  0.1,  0.1,  0.1,  0.1,  0.1])
    >>> np.full_like(x, np.nan, dtype=np.double)
    array([ nan,  nan,  nan,  nan,  nan,  nan])

    >>> y = np.arange(6, dtype=np.double)
    >>> np.full_like(y, 0.1)
    array([ 0.1,  0.1,  0.1,  0.1,  0.1,  0.1])

    """
    res = empty_like(a, dtype=dtype, order=order, subok=subok)
    multiarray.copyto(res, fill_value, casting='unsafe')
    return res


def extend_all(module):
    adict = {}
    for a in __all__:
        adict[a] = 1
    try:
        mall = getattr(module, '__all__')
    except AttributeError:
        mall = [k for k in module.__dict__.keys() if not k.startswith('_')]
    for a in mall:
        if a not in adict:
            __all__.append(a)

newaxis = None


arange = multiarray.arange
array = multiarray.array
zeros = multiarray.zeros
count_nonzero = multiarray.count_nonzero
empty = multiarray.empty
empty_like = multiarray.empty_like
fromstring = multiarray.fromstring
fromiter = multiarray.fromiter
fromfile = multiarray.fromfile
frombuffer = multiarray.frombuffer
may_share_memory = multiarray.may_share_memory
if sys.version_info[0] < 3:
    newbuffer = multiarray.newbuffer
    getbuffer = multiarray.getbuffer
int_asbuffer = multiarray.int_asbuffer
where = multiarray.where
concatenate = multiarray.concatenate
fastCopyAndTranspose = multiarray._fastCopyAndTranspose
set_numeric_ops = multiarray.set_numeric_ops
can_cast = multiarray.can_cast
promote_types = multiarray.promote_types
min_scalar_type = multiarray.min_scalar_type
result_type = multiarray.result_type
lexsort = multiarray.lexsort
compare_chararrays = multiarray.compare_chararrays
putmask = multiarray.putmask
einsum = multiarray.einsum

def asarray(a, dtype=None, order=None):
    """
    Convert the input to an array.

    Parameters
    ----------
    a : array_like
        Input data, in any form that can be converted to an array.  This
        includes lists, lists of tuples, tuples, tuples of tuples, tuples
        of lists and ndarrays.
    dtype : data-type, optional
        By default, the data-type is inferred from the input data.
    order : {'C', 'F'}, optional
        Whether to use row-major ('C') or column-major ('F' for FORTRAN)
        memory representation.  Defaults to 'C'.

    Returns
    -------
    out : ndarray
        Array interpretation of `a`.  No copy is performed if the input
        is already an ndarray.  If `a` is a subclass of ndarray, a base
        class ndarray is returned.

    See Also
    --------
    asanyarray : Similar function which passes through subclasses.
    ascontiguousarray : Convert input to a contiguous array.
    asfarray : Convert input to a floating point ndarray.
    asfortranarray : Convert input to an ndarray with column-major
                     memory order.
    asarray_chkfinite : Similar function which checks input for NaNs and Infs.
    fromiter : Create an array from an iterator.
    fromfunction : Construct an array by executing a function on grid
                   positions.

    Examples
    --------
    Convert a list into an array:

    >>> a = [1, 2]
    >>> np.asarray(a)
    array([1, 2])

    Existing arrays are not copied:

    >>> a = np.array([1, 2])
    >>> np.asarray(a) is a
    True

    If `dtype` is set, array is copied only if dtype does not match:

    >>> a = np.array([1, 2], dtype=np.float32)
    >>> np.asarray(a, dtype=np.float32) is a
    True
    >>> np.asarray(a, dtype=np.float64) is a
    False

    Contrary to `asanyarray`, ndarray subclasses are not passed through:

    >>> issubclass(np.matrix, np.ndarray)
    True
    >>> a = np.matrix([[1, 2]])
    >>> np.asarray(a) is a
    False
    >>> np.asanyarray(a) is a
    True

    """
    return array(a, dtype, copy=False, order=order)

def asanyarray(a, dtype=None, order=None):
    """
    Convert the input to an ndarray, but pass ndarray subclasses through.

    Parameters
    ----------
    a : array_like
        Input data, in any form that can be converted to an array.  This
        includes scalars, lists, lists of tuples, tuples, tuples of tuples,
        tuples of lists, and ndarrays.
    dtype : data-type, optional
        By default, the data-type is inferred from the input data.
    order : {'C', 'F'}, optional
        Whether to use row-major ('C') or column-major ('F') memory
        representation.  Defaults to 'C'.

    Returns
    -------
    out : ndarray or an ndarray subclass
        Array interpretation of `a`.  If `a` is an ndarray or a subclass
        of ndarray, it is returned as-is and no copy is performed.

    See Also
    --------
    asarray : Similar function which always returns ndarrays.
    ascontiguousarray : Convert input to a contiguous array.
    asfarray : Convert input to a floating point ndarray.
    asfortranarray : Convert input to an ndarray with column-major
                     memory order.
    asarray_chkfinite : Similar function which checks input for NaNs and
                        Infs.
    fromiter : Create an array from an iterator.
    fromfunction : Construct an array by executing a function on grid
                   positions.

    Examples
    --------
    Convert a list into an array:

    >>> a = [1, 2]
    >>> np.asanyarray(a)
    array([1, 2])

    Instances of `ndarray` subclasses are passed through as-is:

    >>> a = np.matrix([1, 2])
    >>> np.asanyarray(a) is a
    True

    """
    return array(a, dtype, copy=False, order=order, subok=True)

def ascontiguousarray(a, dtype=None):
    """
    Return a contiguous array in memory (C order).

    Parameters
    ----------
    a : array_like
        Input array.
    dtype : str or dtype object, optional
        Data-type of returned array.

    Returns
    -------
    out : ndarray
        Contiguous array of same shape and content as `a`, with type `dtype`
        if specified.

    See Also
    --------
    asfortranarray : Convert input to an ndarray with column-major
                     memory order.
    require : Return an ndarray that satisfies requirements.
    ndarray.flags : Information about the memory layout of the array.

    Examples
    --------
    >>> x = np.arange(6).reshape(2,3)
    >>> np.ascontiguousarray(x, dtype=np.float32)
    array([[ 0.,  1.,  2.],
           [ 3.,  4.,  5.]], dtype=float32)
    >>> x.flags['C_CONTIGUOUS']
    True

    """
    return array(a, dtype, copy=False, order='C', ndmin=1)

def asfortranarray(a, dtype=None):
    """
    Return an array laid out in Fortran order in memory.

    Parameters
    ----------
    a : array_like
        Input array.
    dtype : str or dtype object, optional
        By default, the data-type is inferred from the input data.

    Returns
    -------
    out : ndarray
        The input `a` in Fortran, or column-major, order.

    See Also
    --------
    ascontiguousarray : Convert input to a contiguous (C order) array.
    asanyarray : Convert input to an ndarray with either row or
        column-major memory order.
    require : Return an ndarray that satisfies requirements.
    ndarray.flags : Information about the memory layout of the array.

    Examples
    --------
    >>> x = np.arange(6).reshape(2,3)
    >>> y = np.asfortranarray(x)
    >>> x.flags['F_CONTIGUOUS']
    False
    >>> y.flags['F_CONTIGUOUS']
    True

    """
    return array(a, dtype, copy=False, order='F', ndmin=1)

def require(a, dtype=None, requirements=None):
    """
    Return an ndarray of the provided type that satisfies requirements.

    This function is useful to be sure that an array with the correct flags
    is returned for passing to compiled code (perhaps through ctypes).

    Parameters
    ----------
    a : array_like
       The object to be converted to a type-and-requirement-satisfying array.
    dtype : data-type
       The required data-type, the default data-type is float64).
    requirements : str or list of str
       The requirements list can be any of the following

       * 'F_CONTIGUOUS' ('F') - ensure a Fortran-contiguous array
       * 'C_CONTIGUOUS' ('C') - ensure a C-contiguous array
       * 'ALIGNED' ('A')      - ensure a data-type aligned array
       * 'WRITEABLE' ('W')    - ensure a writable array
       * 'OWNDATA' ('O')      - ensure an array that owns its own data

    See Also
    --------
    asarray : Convert input to an ndarray.
    asanyarray : Convert to an ndarray, but pass through ndarray subclasses.
    ascontiguousarray : Convert input to a contiguous array.
    asfortranarray : Convert input to an ndarray with column-major
                     memory order.
    ndarray.flags : Information about the memory layout of the array.

    Notes
    -----
    The returned array will be guaranteed to have the listed requirements
    by making a copy if needed.

    Examples
    --------
    >>> x = np.arange(6).reshape(2,3)
    >>> x.flags
      C_CONTIGUOUS : True
      F_CONTIGUOUS : False
      OWNDATA : False
      WRITEABLE : True
      ALIGNED : True
      UPDATEIFCOPY : False

    >>> y = np.require(x, dtype=np.float32, requirements=['A', 'O', 'W', 'F'])
    >>> y.flags
      C_CONTIGUOUS : False
      F_CONTIGUOUS : True
      OWNDATA : True
      WRITEABLE : True
      ALIGNED : True
      UPDATEIFCOPY : False

    """
    if requirements is None:
        requirements = []
    else:
        requirements = [x.upper() for x in requirements]

    if not requirements:
        return asanyarray(a, dtype=dtype)

    if 'ENSUREARRAY' in requirements or 'E' in requirements:
        subok = False
    else:
        subok = True

    arr = array(a, dtype=dtype, copy=False, subok=subok)

    copychar = 'A'
    if 'FORTRAN' in requirements or \
       'F_CONTIGUOUS' in requirements or \
       'F' in requirements:
        copychar = 'F'
    elif 'CONTIGUOUS' in requirements or \
         'C_CONTIGUOUS' in requirements or \
         'C' in requirements:
        copychar = 'C'

    for prop in requirements:
        if not arr.flags[prop]:
            arr = arr.copy(copychar)
            break
    return arr

def isfortran(a):
    """
    Returns True if array is arranged in Fortran-order in memory
    and not C-order.

    Parameters
    ----------
    a : ndarray
        Input array.


    Examples
    --------

    np.array allows to specify whether the array is written in C-contiguous
    order (last index varies the fastest), or FORTRAN-contiguous order in
    memory (first index varies the fastest).

    >>> a = np.array([[1, 2, 3], [4, 5, 6]], order='C')
    >>> a
    array([[1, 2, 3],
           [4, 5, 6]])
    >>> np.isfortran(a)
    False

    >>> b = np.array([[1, 2, 3], [4, 5, 6]], order='FORTRAN')
    >>> b
    array([[1, 2, 3],
           [4, 5, 6]])
    >>> np.isfortran(b)
    True


    The transpose of a C-ordered array is a FORTRAN-ordered array.

    >>> a = np.array([[1, 2, 3], [4, 5, 6]], order='C')
    >>> a
    array([[1, 2, 3],
           [4, 5, 6]])
    >>> np.isfortran(a)
    False
    >>> b = a.T
    >>> b
    array([[1, 4],
           [2, 5],
           [3, 6]])
    >>> np.isfortran(b)
    True

    C-ordered arrays evaluate as False even if they are also FORTRAN-ordered.

    >>> np.isfortran(np.array([1, 2], order='FORTRAN'))
    False

    """
    return a.flags.fnc

def argwhere(a):
    """
    Find the indices of array elements that are non-zero, grouped by element.

    Parameters
    ----------
    a : array_like
        Input data.

    Returns
    -------
    index_array : ndarray
        Indices of elements that are non-zero. Indices are grouped by element.

    See Also
    --------
    where, nonzero

    Notes
    -----
    ``np.argwhere(a)`` is the same as ``np.transpose(np.nonzero(a))``.

    The output of ``argwhere`` is not suitable for indexing arrays.
    For this purpose use ``where(a)`` instead.

    Examples
    --------
    >>> x = np.arange(6).reshape(2,3)
    >>> x
    array([[0, 1, 2],
           [3, 4, 5]])
    >>> np.argwhere(x>1)
    array([[0, 2],
           [1, 0],
           [1, 1],
           [1, 2]])

    """
    return transpose(nonzero(a))

def flatnonzero(a):
    """
    Return indices that are non-zero in the flattened version of a.

    This is equivalent to a.ravel().nonzero()[0].

    Parameters
    ----------
    a : ndarray
        Input array.

    Returns
    -------
    res : ndarray
        Output array, containing the indices of the elements of `a.ravel()`
        that are non-zero.

    See Also
    --------
    nonzero : Return the indices of the non-zero elements of the input array.
    ravel : Return a 1-D array containing the elements of the input array.

    Examples
    --------
    >>> x = np.arange(-2, 3)
    >>> x
    array([-2, -1,  0,  1,  2])
    >>> np.flatnonzero(x)
    array([0, 1, 3, 4])

    Use the indices of the non-zero elements as an index array to extract
    these elements:

    >>> x.ravel()[np.flatnonzero(x)]
    array([-2, -1,  1,  2])

    """
    return a.ravel().nonzero()[0]

_mode_from_name_dict = {'v': 0,
                        's' : 1,
                        'f' : 2}

def _mode_from_name(mode):
    if isinstance(mode, basestring):
        return _mode_from_name_dict[mode.lower()[0]]
    return mode

def correlate(a, v, mode='valid', old_behavior=False):
    """
    Cross-correlation of two 1-dimensional sequences.

    This function computes the correlation as generally defined in signal
    processing texts::

        c_{av}[k] = sum_n a[n+k] * conj(v[n])

    with a and v sequences being zero-padded where necessary and conj being
    the conjugate.

    Parameters
    ----------
    a, v : array_like
        Input sequences.
    mode : {'valid', 'same', 'full'}, optional
        Refer to the `convolve` docstring.  Note that the default
        is `valid`, unlike `convolve`, which uses `full`.
    old_behavior : bool
        If True, uses the old behavior from Numeric,
        (correlate(a,v) == correlate(v,a), and the conjugate is not taken
        for complex arrays). If False, uses the conventional signal
        processing definition.

    Returns
    -------
    out : ndarray
        Discrete cross-correlation of `a` and `v`.

    See Also
    --------
    convolve : Discrete, linear convolution of two one-dimensional sequences.

    Notes
    -----
    The definition of correlation above is not unique and sometimes correlation
    may be defined differently. Another common definition is::

        c'_{av}[k] = sum_n a[n] conj(v[n+k])

    which is related to ``c_{av}[k]`` by ``c'_{av}[k] = c_{av}[-k]``.

    Examples
    --------
    >>> np.correlate([1, 2, 3], [0, 1, 0.5])
    array([ 3.5])
    >>> np.correlate([1, 2, 3], [0, 1, 0.5], "same")
    array([ 2. ,  3.5,  3. ])
    >>> np.correlate([1, 2, 3], [0, 1, 0.5], "full")
    array([ 0.5,  2. ,  3.5,  3. ,  0. ])

    Using complex sequences:

    >>> np.correlate([1+1j, 2, 3-1j], [0, 1, 0.5j], 'full')
    array([ 0.5-0.5j,  1.0+0.j ,  1.5-1.5j,  3.0-1.j ,  0.0+0.j ])

    Note that you get the time reversed, complex conjugated result
    when the two input sequences change places, i.e.,
    ``c_{va}[k] = c^{*}_{av}[-k]``:

    >>> np.correlate([0, 1, 0.5j], [1+1j, 2, 3-1j], 'full')
    array([ 0.0+0.j ,  3.0+1.j ,  1.5+1.5j,  1.0+0.j ,  0.5+0.5j])

    """
    mode = _mode_from_name(mode)
# the old behavior should be made available under a different name, see thread
# http://thread.gmane.org/gmane.comp.python.numeric.general/12609/focus=12630
    if old_behavior:
        warnings.warn("""
The old behavior of correlate was deprecated for 1.4.0, and will be completely removed
for NumPy 2.0.

The new behavior fits the conventional definition of correlation: inputs are
never swapped, and the second argument is conjugated for complex arrays.""",
            DeprecationWarning)
        return multiarray.correlate(a, v, mode)
    else:
        return multiarray.correlate2(a, v, mode)

def convolve(a,v,mode='full'):
    """
    Returns the discrete, linear convolution of two one-dimensional sequences.

    The convolution operator is often seen in signal processing, where it
    models the effect of a linear time-invariant system on a signal [1]_.  In
    probability theory, the sum of two independent random variables is
    distributed according to the convolution of their individual
    distributions.

    If `v` is longer than `a`, the arrays are swapped before computation.

    Parameters
    ----------
    a : (N,) array_like
        First one-dimensional input array.
    v : (M,) array_like
        Second one-dimensional input array.
    mode : {'full', 'valid', 'same'}, optional
        'full':
          By default, mode is 'full'.  This returns the convolution
          at each point of overlap, with an output shape of (N+M-1,). At
          the end-points of the convolution, the signals do not overlap
          completely, and boundary effects may be seen.

        'same':
          Mode `same` returns output of length ``max(M, N)``.  Boundary
          effects are still visible.

        'valid':
          Mode `valid` returns output of length
          ``max(M, N) - min(M, N) + 1``.  The convolution product is only given
          for points where the signals overlap completely.  Values outside
          the signal boundary have no effect.

    Returns
    -------
    out : ndarray
        Discrete, linear convolution of `a` and `v`.

    See Also
    --------
    scipy.signal.fftconvolve : Convolve two arrays using the Fast Fourier
                               Transform.
    scipy.linalg.toeplitz : Used to construct the convolution operator.
    polymul : Polynomial multiplication. Same output as convolve, but also
              accepts poly1d objects as input.

    Notes
    -----
    The discrete convolution operation is defined as

    .. math:: (a * v)[n] = \\sum_{m = -\\infty}^{\\infty} a[m] v[n - m]

    It can be shown that a convolution :math:`x(t) * y(t)` in time/space
    is equivalent to the multiplication :math:`X(f) Y(f)` in the Fourier
    domain, after appropriate padding (padding is necessary to prevent
    circular convolution).  Since multiplication is more efficient (faster)
    than convolution, the function `scipy.signal.fftconvolve` exploits the
    FFT to calculate the convolution of large data-sets.

    References
    ----------
    .. [1] Wikipedia, "Convolution", http://en.wikipedia.org/wiki/Convolution.

    Examples
    --------
    Note how the convolution operator flips the second array
    before "sliding" the two across one another:

    >>> np.convolve([1, 2, 3], [0, 1, 0.5])
    array([ 0. ,  1. ,  2.5,  4. ,  1.5])

    Only return the middle values of the convolution.
    Contains boundary effects, where zeros are taken
    into account:

    >>> np.convolve([1,2,3],[0,1,0.5], 'same')
    array([ 1. ,  2.5,  4. ])

    The two arrays are of the same length, so there
    is only one position where they completely overlap:

    >>> np.convolve([1,2,3],[0,1,0.5], 'valid')
    array([ 2.5])

    """
    a, v = array(a, ndmin=1), array(v, ndmin=1)
    if (len(v) > len(a)):
        a, v = v, a
    if len(a) == 0 :
        raise ValueError('a cannot be empty')
    if len(v) == 0 :
        raise ValueError('v cannot be empty')
    mode = _mode_from_name(mode)
    return multiarray.correlate(a, v[::-1], mode)

def outer(a, b, out=None):
    """
    Compute the outer product of two vectors.

    Given two vectors, ``a = [a0, a1, ..., aM]`` and
    ``b = [b0, b1, ..., bN]``,
    the outer product [1]_ is::

      [[a0*b0  a0*b1 ... a0*bN ]
       [a1*b0    .
       [ ...          .
       [aM*b0            aM*bN ]]

    Parameters
    ----------
    a : (M,) array_like
        First input vector.  Input is flattened if
        not already 1-dimensional.
    b : (N,) array_like
        Second input vector.  Input is flattened if
        not already 1-dimensional.
    out : (M, N) ndarray, optional
          A location where the result is stored

        .. versionadded:: 1.9.0

    Returns
    -------
    out : (M, N) ndarray
        ``out[i, j] = a[i] * b[j]``

    See also
    --------
    inner, einsum

    References
    ----------
    .. [1] : G. H. Golub and C. F. van Loan, *Matrix Computations*, 3rd
             ed., Baltimore, MD, Johns Hopkins University Press, 1996,
             pg. 8.

    Examples
    --------
    Make a (*very* coarse) grid for computing a Mandelbrot set:

    >>> rl = np.outer(np.ones((5,)), np.linspace(-2, 2, 5))
    >>> rl
    array([[-2., -1.,  0.,  1.,  2.],
           [-2., -1.,  0.,  1.,  2.],
           [-2., -1.,  0.,  1.,  2.],
           [-2., -1.,  0.,  1.,  2.],
           [-2., -1.,  0.,  1.,  2.]])
    >>> im = np.outer(1j*np.linspace(2, -2, 5), np.ones((5,)))
    >>> im
    array([[ 0.+2.j,  0.+2.j,  0.+2.j,  0.+2.j,  0.+2.j],
           [ 0.+1.j,  0.+1.j,  0.+1.j,  0.+1.j,  0.+1.j],
           [ 0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j],
           [ 0.-1.j,  0.-1.j,  0.-1.j,  0.-1.j,  0.-1.j],
           [ 0.-2.j,  0.-2.j,  0.-2.j,  0.-2.j,  0.-2.j]])
    >>> grid = rl + im
    >>> grid
    array([[-2.+2.j, -1.+2.j,  0.+2.j,  1.+2.j,  2.+2.j],
           [-2.+1.j, -1.+1.j,  0.+1.j,  1.+1.j,  2.+1.j],
           [-2.+0.j, -1.+0.j,  0.+0.j,  1.+0.j,  2.+0.j],
           [-2.-1.j, -1.-1.j,  0.-1.j,  1.-1.j,  2.-1.j],
           [-2.-2.j, -1.-2.j,  0.-2.j,  1.-2.j,  2.-2.j]])

    An example using a "vector" of letters:

    >>> x = np.array(['a', 'b', 'c'], dtype=object)
    >>> np.outer(x, [1, 2, 3])
    array([[a, aa, aaa],
           [b, bb, bbb],
           [c, cc, ccc]], dtype=object)

    """
    a = asarray(a)
    b = asarray(b)
    return multiply(a.ravel()[:, newaxis], b.ravel()[newaxis,:], out)

# try to import blas optimized dot if available
envbak = os.environ.copy()
try:
    # importing this changes the dot function for basic 4 types
    # to blas-optimized versions.

    # disables openblas affinity setting of the main thread that limits
    # python threads or processes to one core
    if 'OPENBLAS_MAIN_FREE' not in os.environ:
        os.environ['OPENBLAS_MAIN_FREE'] = '1'
    if 'GOTOBLAS_MAIN_FREE' not in os.environ:
        os.environ['GOTOBLAS_MAIN_FREE'] = '1'
    from ._dotblas import dot, vdot, inner, alterdot, restoredot
except ImportError:
    # docstrings are in add_newdocs.py
    inner = multiarray.inner
    dot = multiarray.dot
    def vdot(a, b):
        return dot(asarray(a).ravel().conj(), asarray(b).ravel())
    def alterdot():
        pass
    def restoredot():
        pass
finally:
    os.environ.clear()
    os.environ.update(envbak)
    del envbak

def tensordot(a, b, axes=2):
    """
    Compute tensor dot product along specified axes for arrays >= 1-D.

    Given two tensors (arrays of dimension greater than or equal to one),
    `a` and `b`, and an array_like object containing two array_like
    objects, ``(a_axes, b_axes)``, sum the products of `a`'s and `b`'s
    elements (components) over the axes specified by ``a_axes`` and
    ``b_axes``. The third argument can be a single non-negative
    integer_like scalar, ``N``; if it is such, then the last ``N``
    dimensions of `a` and the first ``N`` dimensions of `b` are summed
    over.

    Parameters
    ----------
    a, b : array_like, len(shape) >= 1
        Tensors to "dot".
    axes : variable type
        * integer_like scalar
          Number of axes to sum over (applies to both arrays); or
        * (2,) array_like, both elements array_like of the same length
          List of axes to be summed over, first sequence applying to `a`,
          second to `b`.

    See Also
    --------
    dot, einsum

    Notes
    -----
    When there is more than one axis to sum over - and they are not the last
    (first) axes of `a` (`b`) - the argument `axes` should consist of
    two sequences of the same length, with the first axis to sum over given
    first in both sequences, the second axis second, and so forth.

    Examples
    --------
    A "traditional" example:

    >>> a = np.arange(60.).reshape(3,4,5)
    >>> b = np.arange(24.).reshape(4,3,2)
    >>> c = np.tensordot(a,b, axes=([1,0],[0,1]))
    >>> c.shape
    (5, 2)
    >>> c
    array([[ 4400.,  4730.],
           [ 4532.,  4874.],
           [ 4664.,  5018.],
           [ 4796.,  5162.],
           [ 4928.,  5306.]])
    >>> # A slower but equivalent way of computing the same...
    >>> d = np.zeros((5,2))
    >>> for i in range(5):
    ...   for j in range(2):
    ...     for k in range(3):
    ...       for n in range(4):
    ...         d[i,j] += a[k,n,i] * b[n,k,j]
    >>> c == d
    array([[ True,  True],
           [ True,  True],
           [ True,  True],
           [ True,  True],
           [ True,  True]], dtype=bool)

    An extended example taking advantage of the overloading of + and \\*:

    >>> a = np.array(range(1, 9))
    >>> a.shape = (2, 2, 2)
    >>> A = np.array(('a', 'b', 'c', 'd'), dtype=object)
    >>> A.shape = (2, 2)
    >>> a; A
    array([[[1, 2],
            [3, 4]],
           [[5, 6],
            [7, 8]]])
    array([[a, b],
           [c, d]], dtype=object)

    >>> np.tensordot(a, A) # third argument default is 2
    array([abbcccdddd, aaaaabbbbbbcccccccdddddddd], dtype=object)

    >>> np.tensordot(a, A, 1)
    array([[[acc, bdd],
            [aaacccc, bbbdddd]],
           [[aaaaacccccc, bbbbbdddddd],
            [aaaaaaacccccccc, bbbbbbbdddddddd]]], dtype=object)

    >>> np.tensordot(a, A, 0) # "Left for reader" (result too long to incl.)
    array([[[[[a, b],
              [c, d]],
              ...

    >>> np.tensordot(a, A, (0, 1))
    array([[[abbbbb, cddddd],
            [aabbbbbb, ccdddddd]],
           [[aaabbbbbbb, cccddddddd],
            [aaaabbbbbbbb, ccccdddddddd]]], dtype=object)

    >>> np.tensordot(a, A, (2, 1))
    array([[[abb, cdd],
            [aaabbbb, cccdddd]],
           [[aaaaabbbbbb, cccccdddddd],
            [aaaaaaabbbbbbbb, cccccccdddddddd]]], dtype=object)

    >>> np.tensordot(a, A, ((0, 1), (0, 1)))
    array([abbbcccccddddddd, aabbbbccccccdddddddd], dtype=object)

    >>> np.tensordot(a, A, ((2, 1), (1, 0)))
    array([acccbbdddd, aaaaacccccccbbbbbbdddddddd], dtype=object)

    """
    try:
        iter(axes)
    except:
        axes_a = list(range(-axes, 0))
        axes_b = list(range(0, axes))
    else:
        axes_a, axes_b = axes
    try:
        na = len(axes_a)
        axes_a = list(axes_a)
    except TypeError:
        axes_a = [axes_a]
        na = 1
    try:
        nb = len(axes_b)
        axes_b = list(axes_b)
    except TypeError:
        axes_b = [axes_b]
        nb = 1

    a, b = asarray(a), asarray(b)
    as_ = a.shape
    nda = len(a.shape)
    bs = b.shape
    ndb = len(b.shape)
    equal = True
    if (na != nb): equal = False
    else:
        for k in range(na):
            if as_[axes_a[k]] != bs[axes_b[k]]:
                equal = False
                break
            if axes_a[k] < 0:
                axes_a[k] += nda
            if axes_b[k] < 0:
                axes_b[k] += ndb
    if not equal:
        raise ValueError("shape-mismatch for sum")

    # Move the axes to sum over to the end of "a"
    # and to the front of "b"
    notin = [k for k in range(nda) if k not in axes_a]
    newaxes_a = notin + axes_a
    N2 = 1
    for axis in axes_a:
        N2 *= as_[axis]
    newshape_a = (-1, N2)
    olda = [as_[axis] for axis in notin]

    notin = [k for k in range(ndb) if k not in axes_b]
    newaxes_b = axes_b + notin
    N2 = 1
    for axis in axes_b:
        N2 *= bs[axis]
    newshape_b = (N2, -1)
    oldb = [bs[axis] for axis in notin]

    at = a.transpose(newaxes_a).reshape(newshape_a)
    bt = b.transpose(newaxes_b).reshape(newshape_b)
    res = dot(at, bt)
    return res.reshape(olda + oldb)

def roll(a, shift, axis=None):
    """
    Roll array elements along a given axis.

    Elements that roll beyond the last position are re-introduced at
    the first.

    Parameters
    ----------
    a : array_like
        Input array.
    shift : int
        The number of places by which elements are shifted.
    axis : int, optional
        The axis along which elements are shifted.  By default, the array
        is flattened before shifting, after which the original
        shape is restored.

    Returns
    -------
    res : ndarray
        Output array, with the same shape as `a`.

    See Also
    --------
    rollaxis : Roll the specified axis backwards, until it lies in a
               given position.

    Examples
    --------
    >>> x = np.arange(10)
    >>> np.roll(x, 2)
    array([8, 9, 0, 1, 2, 3, 4, 5, 6, 7])

    >>> x2 = np.reshape(x, (2,5))
    >>> x2
    array([[0, 1, 2, 3, 4],
           [5, 6, 7, 8, 9]])
    >>> np.roll(x2, 1)
    array([[9, 0, 1, 2, 3],
           [4, 5, 6, 7, 8]])
    >>> np.roll(x2, 1, axis=0)
    array([[5, 6, 7, 8, 9],
           [0, 1, 2, 3, 4]])
    >>> np.roll(x2, 1, axis=1)
    array([[4, 0, 1, 2, 3],
           [9, 5, 6, 7, 8]])

    """
    a = asanyarray(a)
    if axis is None:
        n = a.size
        reshape = True
    else:
        try:
            n = a.shape[axis]
        except IndexError:
            raise ValueError('axis must be >= 0 and < %d' % a.ndim)
        reshape = False
    if n == 0:
        return a
    shift %= n
    indexes = concatenate((arange(n - shift, n), arange(n - shift)))
    res = a.take(indexes, axis)
    if reshape:
        res = res.reshape(a.shape)
    return res

def rollaxis(a, axis, start=0):
    """
    Roll the specified axis backwards, until it lies in a given position.

    Parameters
    ----------
    a : ndarray
        Input array.
    axis : int
        The axis to roll backwards.  The positions of the other axes do not
        change relative to one another.
    start : int, optional
        The axis is rolled until it lies before this position.  The default,
        0, results in a "complete" roll.

    Returns
    -------
    res : ndarray
        Output array.

    See Also
    --------
    roll : Roll the elements of an array by a number of positions along a
        given axis.

    Examples
    --------
    >>> a = np.ones((3,4,5,6))
    >>> np.rollaxis(a, 3, 1).shape
    (3, 6, 4, 5)
    >>> np.rollaxis(a, 2).shape
    (5, 3, 4, 6)
    >>> np.rollaxis(a, 1, 4).shape
    (3, 5, 6, 4)

    """
    n = a.ndim
    if axis < 0:
        axis += n
    if start < 0:
        start += n
    msg = 'rollaxis: %s (%d) must be >=0 and < %d'
    if not (0 <= axis < n):
        raise ValueError(msg % ('axis', axis, n))
    if not (0 <= start < n+1):
        raise ValueError(msg % ('start', start, n+1))
    if (axis < start): # it's been removed
        start -= 1
    if axis==start:
        return a
    axes = list(range(0, n))
    axes.remove(axis)
    axes.insert(start, axis)
    return a.transpose(axes)

# fix hack in scipy which imports this function
def _move_axis_to_0(a, axis):
    return rollaxis(a, axis, 0)

def cross(a, b, axisa=-1, axisb=-1, axisc=-1, axis=None):
    """
    Return the cross product of two (arrays of) vectors.

    The cross product of `a` and `b` in :math:`R^3` is a vector perpendicular
    to both `a` and `b`.  If `a` and `b` are arrays of vectors, the vectors
    are defined by the last axis of `a` and `b` by default, and these axes
    can have dimensions 2 or 3.  Where the dimension of either `a` or `b` is
    2, the third component of the input vector is assumed to be zero and the
    cross product calculated accordingly.  In cases where both input vectors
    have dimension 2, the z-component of the cross product is returned.

    Parameters
    ----------
    a : array_like
        Components of the first vector(s).
    b : array_like
        Components of the second vector(s).
    axisa : int, optional
        Axis of `a` that defines the vector(s).  By default, the last axis.
    axisb : int, optional
        Axis of `b` that defines the vector(s).  By default, the last axis.
    axisc : int, optional
        Axis of `c` containing the cross product vector(s).  By default, the
        last axis.
    axis : int, optional
        If defined, the axis of `a`, `b` and `c` that defines the vector(s)
        and cross product(s).  Overrides `axisa`, `axisb` and `axisc`.

    Returns
    -------
    c : ndarray
        Vector cross product(s).

    Raises
    ------
    ValueError
        When the dimension of the vector(s) in `a` and/or `b` does not
        equal 2 or 3.

    See Also
    --------
    inner : Inner product
    outer : Outer product.
    ix_ : Construct index arrays.

    Notes
    -----
    .. versionadded:: 1.9.0
    Supports full broadcasting of the inputs.

    Examples
    --------
    Vector cross-product.

    >>> x = [1, 2, 3]
    >>> y = [4, 5, 6]
    >>> np.cross(x, y)
    array([-3,  6, -3])

    One vector with dimension 2.

    >>> x = [1, 2]
    >>> y = [4, 5, 6]
    >>> np.cross(x, y)
    array([12, -6, -3])

    Equivalently:

    >>> x = [1, 2, 0]
    >>> y = [4, 5, 6]
    >>> np.cross(x, y)
    array([12, -6, -3])

    Both vectors with dimension 2.

    >>> x = [1,2]
    >>> y = [4,5]
    >>> np.cross(x, y)
    -3

    Multiple vector cross-products. Note that the direction of the cross
    product vector is defined by the `right-hand rule`.

    >>> x = np.array([[1,2,3], [4,5,6]])
    >>> y = np.array([[4,5,6], [1,2,3]])
    >>> np.cross(x, y)
    array([[-3,  6, -3],
           [ 3, -6,  3]])

    The orientation of `c` can be changed using the `axisc` keyword.

    >>> np.cross(x, y, axisc=0)
    array([[-3,  3],
           [ 6, -6],
           [-3,  3]])

    Change the vector definition of `x` and `y` using `axisa` and `axisb`.

    >>> x = np.array([[1,2,3], [4,5,6], [7, 8, 9]])
    >>> y = np.array([[7, 8, 9], [4,5,6], [1,2,3]])
    >>> np.cross(x, y)
    array([[ -6,  12,  -6],
           [  0,   0,   0],
           [  6, -12,   6]])
    >>> np.cross(x, y, axisa=0, axisb=0)
    array([[-24,  48, -24],
           [-30,  60, -30],
           [-36,  72, -36]])

    """
    if axis is not None:
        axisa, axisb, axisc = (axis,) * 3
    a = asarray(a)
    b = asarray(b)
    # Move working axis to the end of the shape
    a = rollaxis(a, axisa, a.ndim)
    b = rollaxis(b, axisb, b.ndim)
    msg = ("incompatible dimensions for cross product\n"
           "(dimension must be 2 or 3)")
    if a.shape[-1] not in (2, 3) or b.shape[-1] not in (2, 3):
        raise ValueError(msg)

        # Create the output array
    shape = broadcast(a[..., 0], b[..., 0]).shape
    if a.shape[-1] == 3 or b.shape[-1] == 3:
        shape += (3,)
    dtype = promote_types(a.dtype, b.dtype)
    cp = empty(shape, dtype)

    # create local aliases for readability
    a0 = a[..., 0]
    a1 = a[..., 1]
    if a.shape[-1] == 3:
        a2 = a[..., 2]
    b0 = b[..., 0]
    b1 = b[..., 1]
    if b.shape[-1] == 3:
        b2 = b[..., 2]
    if cp.ndim != 0 and cp.shape[-1] == 3:
        cp0 = cp[..., 0]
        cp1 = cp[..., 1]
        cp2 = cp[..., 2]

    if a.shape[-1] == 2:
        if b.shape[-1] == 2:
            # a0 * b1 - a1 * b0
            multiply(a0, b1, out=cp)
            cp -= a1 * b0
            if cp.ndim == 0:
                return cp
            else:
                # This works because we are moving the last axis
                return rollaxis(cp, -1, axisc)
        else:
            # cp0 = a1 * b2 - 0  (a2 = 0)
            # cp1 = 0 - a0 * b2  (a2 = 0)
            # cp2 = a0 * b1 - a1 * b0
            multiply(a1, b2, out=cp0)
            multiply(a0, b2, out=cp1)
            negative(cp1, out=cp1)
            multiply(a0, b1, out=cp2)
            cp2 -= a1 * b0
    elif a.shape[-1] == 3:
        if b.shape[-1] == 3:
            # cp0 = a1 * b2 - a2 * b1
            # cp1 = a2 * b0 - a0 * b2
            # cp2 = a0 * b1 - a1 * b0
            multiply(a1, b2, out=cp0)
            tmp = array(a2 * b1)
            cp0 -= tmp
            multiply(a2, b0, out=cp1)
            multiply(a0, b2, out=tmp)
            cp1 -= tmp
            multiply(a0, b1, out=cp2)
            multiply(a1, b0, out=tmp)
            cp2 -= tmp
        else:
            # cp0 = 0 - a2 * b1  (b2 = 0)
            # cp1 = a2 * b0 - 0  (b2 = 0)
            # cp2 = a0 * b1 - a1 * b0
            multiply(a2, b1, out=cp0)
            negative(cp0, out=cp0)
            multiply(a2, b0, out=cp1)
            multiply(a0, b1, out=cp2)
            cp2 -= a1 * b0

    if cp.ndim == 1:
        return cp
    else:
        # This works because we are moving the last axis
        return rollaxis(cp, -1, axisc)

#Use numarray's printing function
from .arrayprint import array2string, get_printoptions, set_printoptions

_typelessdata = [int_, float_, complex_]
if issubclass(intc, int):
    _typelessdata.append(intc)

if issubclass(longlong, int):
    _typelessdata.append(longlong)

def array_repr(arr, max_line_width=None, precision=None, suppress_small=None):
    """
    Return the string representation of an array.

    Parameters
    ----------
    arr : ndarray
        Input array.
    max_line_width : int, optional
        The maximum number of columns the string should span. Newline
        characters split the string appropriately after array elements.
    precision : int, optional
        Floating point precision. Default is the current printing precision
        (usually 8), which can be altered using `set_printoptions`.
    suppress_small : bool, optional
        Represent very small numbers as zero, default is False. Very small
        is defined by `precision`, if the precision is 8 then
        numbers smaller than 5e-9 are represented as zero.

    Returns
    -------
    string : str
      The string representation of an array.

    See Also
    --------
    array_str, array2string, set_printoptions

    Examples
    --------
    >>> np.array_repr(np.array([1,2]))
    'array([1, 2])'
    >>> np.array_repr(np.ma.array([0.]))
    'MaskedArray([ 0.])'
    >>> np.array_repr(np.array([], np.int32))
    'array([], dtype=int32)'

    >>> x = np.array([1e-6, 4e-7, 2, 3])
    >>> np.array_repr(x, precision=6, suppress_small=True)
    'array([ 0.000001,  0.      ,  2.      ,  3.      ])'

    """
    if arr.size > 0 or arr.shape==(0,):
        lst = array2string(arr, max_line_width, precision, suppress_small,
                           ', ', "array(")
    else: # show zero-length shape unless it is (0,)
        lst = "[], shape=%s" % (repr(arr.shape),)

    if arr.__class__ is not ndarray:
        cName= arr.__class__.__name__
    else:
        cName = "array"

    skipdtype = (arr.dtype.type in _typelessdata) and arr.size > 0

    if skipdtype:
        return "%s(%s)" % (cName, lst)
    else:
        typename = arr.dtype.name
        # Quote typename in the output if it is "complex".
        if typename and not (typename[0].isalpha() and typename.isalnum()):
            typename = "'%s'" % typename

        lf = ''
        if issubclass(arr.dtype.type, flexible):
            if arr.dtype.names:
                typename = "%s" % str(arr.dtype)
            else:
                typename = "'%s'" % str(arr.dtype)
            lf = '\n'+' '*len("array(")
        return cName + "(%s, %sdtype=%s)" % (lst, lf, typename)

def array_str(a, max_line_width=None, precision=None, suppress_small=None):
    """
    Return a string representation of the data in an array.

    The data in the array is returned as a single string.  This function is
    similar to `array_repr`, the difference being that `array_repr` also
    returns information on the kind of array and its data type.

    Parameters
    ----------
    a : ndarray
        Input array.
    max_line_width : int, optional
        Inserts newlines if text is longer than `max_line_width`.  The
        default is, indirectly, 75.
    precision : int, optional
        Floating point precision.  Default is the current printing precision
        (usually 8), which can be altered using `set_printoptions`.
    suppress_small : bool, optional
        Represent numbers "very close" to zero as zero; default is False.
        Very close is defined by precision: if the precision is 8, e.g.,
        numbers smaller (in absolute value) than 5e-9 are represented as
        zero.

    See Also
    --------
    array2string, array_repr, set_printoptions

    Examples
    --------
    >>> np.array_str(np.arange(3))
    '[0 1 2]'

    """
    return array2string(a, max_line_width, precision, suppress_small, ' ', "", str)

def set_string_function(f, repr=True):
    """
    Set a Python function to be used when pretty printing arrays.

    Parameters
    ----------
    f : function or None
        Function to be used to pretty print arrays. The function should expect
        a single array argument and return a string of the representation of
        the array. If None, the function is reset to the default NumPy function
        to print arrays.
    repr : bool, optional
        If True (default), the function for pretty printing (``__repr__``)
        is set, if False the function that returns the default string
        representation (``__str__``) is set.

    See Also
    --------
    set_printoptions, get_printoptions

    Examples
    --------
    >>> def pprint(arr):
    ...     return 'HA! - What are you going to do now?'
    ...
    >>> np.set_string_function(pprint)
    >>> a = np.arange(10)
    >>> a
    HA! - What are you going to do now?
    >>> print a
    [0 1 2 3 4 5 6 7 8 9]

    We can reset the function to the default:

    >>> np.set_string_function(None)
    >>> a
    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

    `repr` affects either pretty printing or normal string representation.
    Note that ``__repr__`` is still affected by setting ``__str__``
    because the width of each array element in the returned string becomes
    equal to the length of the result of ``__str__()``.

    >>> x = np.arange(4)
    >>> np.set_string_function(lambda x:'random', repr=False)
    >>> x.__str__()
    'random'
    >>> x.__repr__()
    'array([     0,      1,      2,      3])'

    """
    if f is None:
        if repr:
            return multiarray.set_string_function(array_repr, 1)
        else:
            return multiarray.set_string_function(array_str, 0)
    else:
        return multiarray.set_string_function(f, repr)

set_string_function(array_str, 0)
set_string_function(array_repr, 1)

little_endian = (sys.byteorder == 'little')


def indices(dimensions, dtype=int):
    """
    Return an array representing the indices of a grid.

    Compute an array where the subarrays contain index values 0,1,...
    varying only along the corresponding axis.

    Parameters
    ----------
    dimensions : sequence of ints
        The shape of the grid.
    dtype : dtype, optional
        Data type of the result.

    Returns
    -------
    grid : ndarray
        The array of grid indices,
        ``grid.shape = (len(dimensions),) + tuple(dimensions)``.

    See Also
    --------
    mgrid, meshgrid

    Notes
    -----
    The output shape is obtained by prepending the number of dimensions
    in front of the tuple of dimensions, i.e. if `dimensions` is a tuple
    ``(r0, ..., rN-1)`` of length ``N``, the output shape is
    ``(N,r0,...,rN-1)``.

    The subarrays ``grid[k]`` contains the N-D array of indices along the
    ``k-th`` axis. Explicitly::

        grid[k,i0,i1,...,iN-1] = ik

    Examples
    --------
    >>> grid = np.indices((2, 3))
    >>> grid.shape
    (2, 2, 3)
    >>> grid[0]        # row indices
    array([[0, 0, 0],
           [1, 1, 1]])
    >>> grid[1]        # column indices
    array([[0, 1, 2],
           [0, 1, 2]])

    The indices can be used as an index into an array.

    >>> x = np.arange(20).reshape(5, 4)
    >>> row, col = np.indices((2, 3))
    >>> x[row, col]
    array([[0, 1, 2],
           [4, 5, 6]])

    Note that it would be more straightforward in the above example to
    extract the required elements directly with ``x[:2, :3]``.

    """
    dimensions = tuple(dimensions)
    N = len(dimensions)
    if N == 0:
        return array([], dtype=dtype)
    res = empty((N,)+dimensions, dtype=dtype)
    for i, dim in enumerate(dimensions):
        tmp = arange(dim, dtype=dtype)
        tmp.shape = (1,)*i + (dim,)+(1,)*(N-i-1)
        newdim = dimensions[:i] + (1,)+ dimensions[i+1:]
        val = zeros(newdim, dtype)
        add(tmp, val, res[i])
    return res

def fromfunction(function, shape, **kwargs):
    """
    Construct an array by executing a function over each coordinate.

    The resulting array therefore has a value ``fn(x, y, z)`` at
    coordinate ``(x, y, z)``.

    Parameters
    ----------
    function : callable
        The function is called with N parameters, where N is the rank of
        `shape`.  Each parameter represents the coordinates of the array
        varying along a specific axis.  For example, if `shape`
        were ``(2, 2)``, then the parameters in turn be (0, 0), (0, 1),
        (1, 0), (1, 1).
    shape : (N,) tuple of ints
        Shape of the output array, which also determines the shape of
        the coordinate arrays passed to `function`.
    dtype : data-type, optional
        Data-type of the coordinate arrays passed to `function`.
        By default, `dtype` is float.

    Returns
    -------
    fromfunction : any
        The result of the call to `function` is passed back directly.
        Therefore the shape of `fromfunction` is completely determined by
        `function`.  If `function` returns a scalar value, the shape of
        `fromfunction` would match the `shape` parameter.

    See Also
    --------
    indices, meshgrid

    Notes
    -----
    Keywords other than `dtype` are passed to `function`.

    Examples
    --------
    >>> np.fromfunction(lambda i, j: i == j, (3, 3), dtype=int)
    array([[ True, False, False],
           [False,  True, False],
           [False, False,  True]], dtype=bool)

    >>> np.fromfunction(lambda i, j: i + j, (3, 3), dtype=int)
    array([[0, 1, 2],
           [1, 2, 3],
           [2, 3, 4]])

    """
    dtype = kwargs.pop('dtype', float)
    args = indices(shape, dtype=dtype)
    return function(*args,**kwargs)

def isscalar(num):
    """
    Returns True if the type of `num` is a scalar type.

    Parameters
    ----------
    num : any
        Input argument, can be of any type and shape.

    Returns
    -------
    val : bool
        True if `num` is a scalar type, False if it is not.

    Examples
    --------
    >>> np.isscalar(3.1)
    True
    >>> np.isscalar([3.1])
    False
    >>> np.isscalar(False)
    True

    """
    if isinstance(num, generic):
        return True
    else:
        return type(num) in ScalarType

_lkup = {
    '0':'0000',
    '1':'0001',
    '2':'0010',
    '3':'0011',
    '4':'0100',
    '5':'0101',
    '6':'0110',
    '7':'0111',
    '8':'1000',
    '9':'1001',
    'a':'1010',
    'b':'1011',
    'c':'1100',
    'd':'1101',
    'e':'1110',
    'f':'1111',
    'A':'1010',
    'B':'1011',
    'C':'1100',
    'D':'1101',
    'E':'1110',
    'F':'1111',
    'L':''}

def binary_repr(num, width=None):
    """
    Return the binary representation of the input number as a string.

    For negative numbers, if width is not given, a minus sign is added to the
    front. If width is given, the two's complement of the number is
    returned, with respect to that width.

    In a two's-complement system negative numbers are represented by the two's
    complement of the absolute value. This is the most common method of
    representing signed integers on computers [1]_. A N-bit two's-complement
    system can represent every integer in the range
    :math:`-2^{N-1}` to :math:`+2^{N-1}-1`.

    Parameters
    ----------
    num : int
        Only an integer decimal number can be used.
    width : int, optional
        The length of the returned string if `num` is positive, the length of
        the two's complement if `num` is negative.

    Returns
    -------
    bin : str
        Binary representation of `num` or two's complement of `num`.

    See Also
    --------
    base_repr: Return a string representation of a number in the given base
               system.

    Notes
    -----
    `binary_repr` is equivalent to using `base_repr` with base 2, but about 25x
    faster.

    References
    ----------
    .. [1] Wikipedia, "Two's complement",
        http://en.wikipedia.org/wiki/Two's_complement

    Examples
    --------
    >>> np.binary_repr(3)
    '11'
    >>> np.binary_repr(-3)
    '-11'
    >>> np.binary_repr(3, width=4)
    '0011'

    The two's complement is returned when the input number is negative and
    width is specified:

    >>> np.binary_repr(-3, width=4)
    '1101'

    """
    # ' <-- unbreak Emacs fontification
    sign = ''
    if num < 0:
        if width is None:
            sign = '-'
            num = -num
        else:
            # replace num with its 2-complement
            num = 2**width + num
    elif num == 0:
        return '0'*(width or 1)
    ostr = hex(num)
    bin = ''.join([_lkup[ch] for ch in ostr[2:]])
    bin = bin.lstrip('0')
    if width is not None:
        bin = bin.zfill(width)
    return sign + bin

def base_repr(number, base=2, padding=0):
    """
    Return a string representation of a number in the given base system.

    Parameters
    ----------
    number : int
        The value to convert. Only positive values are handled.
    base : int, optional
        Convert `number` to the `base` number system. The valid range is 2-36,
        the default value is 2.
    padding : int, optional
        Number of zeros padded on the left. Default is 0 (no padding).

    Returns
    -------
    out : str
        String representation of `number` in `base` system.

    See Also
    --------
    binary_repr : Faster version of `base_repr` for base 2.

    Examples
    --------
    >>> np.base_repr(5)
    '101'
    >>> np.base_repr(6, 5)
    '11'
    >>> np.base_repr(7, base=5, padding=3)
    '00012'

    >>> np.base_repr(10, base=16)
    'A'
    >>> np.base_repr(32, base=16)
    '20'

    """
    digits = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    if base > len(digits):
        raise ValueError("Bases greater than 36 not handled in base_repr.")

    num = abs(number)
    res = []
    while num:
        res.append(digits[num % base])
        num //= base
    if padding:
        res.append('0' * padding)
    if number < 0:
        res.append('-')
    return ''.join(reversed(res or '0'))


def load(file):
    """
    Wrapper around cPickle.load which accepts either a file-like object or
    a filename.

    Note that the NumPy binary format is not based on pickle/cPickle anymore.
    For details on the preferred way of loading and saving files, see `load`
    and `save`.

    See Also
    --------
    load, save

    """
    if isinstance(file, type("")):
        file = open(file, "rb")
    return pickle.load(file)

# These are all essentially abbreviations
# These might wind up in a special abbreviations module

def _maketup(descr, val):
    dt = dtype(descr)
    # Place val in all scalar tuples:
    fields = dt.fields
    if fields is None:
        return val
    else:
        res = [_maketup(fields[name][0], val) for name in dt.names]
        return tuple(res)

def identity(n, dtype=None):
    """
    Return the identity array.

    The identity array is a square array with ones on
    the main diagonal.

    Parameters
    ----------
    n : int
        Number of rows (and columns) in `n` x `n` output.
    dtype : data-type, optional
        Data-type of the output.  Defaults to ``float``.

    Returns
    -------
    out : ndarray
        `n` x `n` array with its main diagonal set to one,
        and all other elements 0.

    Examples
    --------
    >>> np.identity(3)
    array([[ 1.,  0.,  0.],
           [ 0.,  1.,  0.],
           [ 0.,  0.,  1.]])

    """
    from numpy import eye
    return eye(n, dtype=dtype)

def allclose(a, b, rtol=1.e-5, atol=1.e-8):
    """
    Returns True if two arrays are element-wise equal within a tolerance.

    The tolerance values are positive, typically very small numbers.  The
    relative difference (`rtol` * abs(`b`)) and the absolute difference
    `atol` are added together to compare against the absolute difference
    between `a` and `b`.

    If either array contains one or more NaNs, False is returned.
    Infs are treated as equal if they are in the same place and of the same
    sign in both arrays.

    Parameters
    ----------
    a, b : array_like
        Input arrays to compare.
    rtol : float
        The relative tolerance parameter (see Notes).
    atol : float
        The absolute tolerance parameter (see Notes).

    Returns
    -------
    allclose : bool
        Returns True if the two arrays are equal within the given
        tolerance; False otherwise.

    See Also
    --------
    isclose, all, any

    Notes
    -----
    If the following equation is element-wise True, then allclose returns
    True.

     absolute(`a` - `b`) <= (`atol` + `rtol` * absolute(`b`))

    The above equation is not symmetric in `a` and `b`, so that
    `allclose(a, b)` might be different from `allclose(b, a)` in
    some rare cases.

    Examples
    --------
    >>> np.allclose([1e10,1e-7], [1.00001e10,1e-8])
    False
    >>> np.allclose([1e10,1e-8], [1.00001e10,1e-9])
    True
    >>> np.allclose([1e10,1e-8], [1.0001e10,1e-9])
    False
    >>> np.allclose([1.0, np.nan], [1.0, np.nan])
    False

    """
    x = array(a, copy=False, ndmin=1)
    y = array(b, copy=False, ndmin=1)

    # make sure y is an inexact type to avoid abs(MIN_INT); will cause
    # casting of x later.
    dtype = multiarray.result_type(y, 1.)
    y = array(y, dtype=dtype, copy=False)

    xinf = isinf(x)
    yinf = isinf(y)
    if any(xinf) or any(yinf):
        # Check that x and y have inf's only in the same positions
        if not all(xinf == yinf):
            return False
        # Check that sign of inf's in x and y is the same
        if not all(x[xinf] == y[xinf]):
            return False

        x = x[~xinf]
        y = y[~xinf]

    # ignore invalid fpe's
    with errstate(invalid='ignore'):
        r = all(less_equal(abs(x - y), atol + rtol * abs(y)))

    return r

def isclose(a, b, rtol=1.e-5, atol=1.e-8, equal_nan=False):
    """
    Returns a boolean array where two arrays are element-wise equal within a
    tolerance.

    The tolerance values are positive, typically very small numbers.  The
    relative difference (`rtol` * abs(`b`)) and the absolute difference
    `atol` are added together to compare against the absolute difference
    between `a` and `b`.

    Parameters
    ----------
    a, b : array_like
        Input arrays to compare.
    rtol : float
        The relative tolerance parameter (see Notes).
    atol : float
        The absolute tolerance parameter (see Notes).
    equal_nan : bool
        Whether to compare NaN's as equal.  If True, NaN's in `a` will be
        considered equal to NaN's in `b` in the output array.

    Returns
    -------
    y : array_like
        Returns a boolean array of where `a` and `b` are equal within the
        given tolerance. If both `a` and `b` are scalars, returns a single
        boolean value.

    See Also
    --------
    allclose

    Notes
    -----
    .. versionadded:: 1.7.0

    For finite values, isclose uses the following equation to test whether
    two floating point values are equivalent.

     absolute(`a` - `b`) <= (`atol` + `rtol` * absolute(`b`))

    The above equation is not symmetric in `a` and `b`, so that
    `isclose(a, b)` might be different from `isclose(b, a)` in
    some rare cases.

    Examples
    --------
    >>> np.isclose([1e10,1e-7], [1.00001e10,1e-8])
    array([True, False])
    >>> np.isclose([1e10,1e-8], [1.00001e10,1e-9])
    array([True, True])
    >>> np.isclose([1e10,1e-8], [1.0001e10,1e-9])
    array([False, True])
    >>> np.isclose([1.0, np.nan], [1.0, np.nan])
    array([True, False])
    >>> np.isclose([1.0, np.nan], [1.0, np.nan], equal_nan=True)
    array([True, True])
    """
    def within_tol(x, y, atol, rtol):
        with errstate(invalid='ignore'):
            result = less_equal(abs(x-y), atol + rtol * abs(y))
        if isscalar(a) and isscalar(b):
            result = bool(result)
        return result

    x = array(a, copy=False, subok=True, ndmin=1)
    y = array(b, copy=False, subok=True, ndmin=1)
    xfin = isfinite(x)
    yfin = isfinite(y)
    if all(xfin) and all(yfin):
        return within_tol(x, y, atol, rtol)
    else:
        finite = xfin & yfin
        cond = zeros_like(finite, subok=True)
        # Because we're using boolean indexing, x & y must be the same shape.
        # Ideally, we'd just do x, y = broadcast_arrays(x, y). It's in
        # lib.stride_tricks, though, so we can't import it here.
        x = x * ones_like(cond)
        y = y * ones_like(cond)
        # Avoid subtraction with infinite/nan values...
        cond[finite] = within_tol(x[finite], y[finite], atol, rtol)
        # Check for equality of infinite values...
        cond[~finite] = (x[~finite] == y[~finite])
        if equal_nan:
            # Make NaN == NaN
            both_nan = isnan(x) & isnan(y)
            cond[both_nan] = both_nan[both_nan]
        return cond

def array_equal(a1, a2):
    """
    True if two arrays have the same shape and elements, False otherwise.

    Parameters
    ----------
    a1, a2 : array_like
        Input arrays.

    Returns
    -------
    b : bool
        Returns True if the arrays are equal.

    See Also
    --------
    allclose: Returns True if two arrays are element-wise equal within a
              tolerance.
    array_equiv: Returns True if input arrays are shape consistent and all
                 elements equal.

    Examples
    --------
    >>> np.array_equal([1, 2], [1, 2])
    True
    >>> np.array_equal(np.array([1, 2]), np.array([1, 2]))
    True
    >>> np.array_equal([1, 2], [1, 2, 3])
    False
    >>> np.array_equal([1, 2], [1, 4])
    False

    """
    try:
        a1, a2 = asarray(a1), asarray(a2)
    except:
        return False
    if a1.shape != a2.shape:
        return False
    return bool(asarray(a1 == a2).all())

def array_equiv(a1, a2):
    """
    Returns True if input arrays are shape consistent and all elements equal.

    Shape consistent means they are either the same shape, or one input array
    can be broadcasted to create the same shape as the other one.

    Parameters
    ----------
    a1, a2 : array_like
        Input arrays.

    Returns
    -------
    out : bool
        True if equivalent, False otherwise.

    Examples
    --------
    >>> np.array_equiv([1, 2], [1, 2])
    True
    >>> np.array_equiv([1, 2], [1, 3])
    False

    Showing the shape equivalence:

    >>> np.array_equiv([1, 2], [[1, 2], [1, 2]])
    True
    >>> np.array_equiv([1, 2], [[1, 2, 1, 2], [1, 2, 1, 2]])
    False

    >>> np.array_equiv([1, 2], [[1, 2], [1, 3]])
    False

    """
    try:
        a1, a2 = asarray(a1), asarray(a2)
    except:
        return False
    try:
        multiarray.broadcast(a1, a2)
    except:
        return False

    return bool(asarray(a1 == a2).all())


_errdict = {"ignore":ERR_IGNORE,
            "warn":ERR_WARN,
            "raise":ERR_RAISE,
            "call":ERR_CALL,
            "print":ERR_PRINT,
            "log":ERR_LOG}

_errdict_rev = {}
for key in _errdict.keys():
    _errdict_rev[_errdict[key]] = key
del key

def seterr(all=None, divide=None, over=None, under=None, invalid=None):
    """
    Set how floating-point errors are handled.

    Note that operations on integer scalar types (such as `int16`) are
    handled like floating point, and are affected by these settings.

    Parameters
    ----------
    all : {'ignore', 'warn', 'raise', 'call', 'print', 'log'}, optional
        Set treatment for all types of floating-point errors at once:

        - ignore: Take no action when the exception occurs.
        - warn: Print a `RuntimeWarning` (via the Python `warnings` module).
        - raise: Raise a `FloatingPointError`.
        - call: Call a function specified using the `seterrcall` function.
        - print: Print a warning directly to ``stdout``.
        - log: Record error in a Log object specified by `seterrcall`.

        The default is not to change the current behavior.
    divide : {'ignore', 'warn', 'raise', 'call', 'print', 'log'}, optional
        Treatment for division by zero.
    over : {'ignore', 'warn', 'raise', 'call', 'print', 'log'}, optional
        Treatment for floating-point overflow.
    under : {'ignore', 'warn', 'raise', 'call', 'print', 'log'}, optional
        Treatment for floating-point underflow.
    invalid : {'ignore', 'warn', 'raise', 'call', 'print', 'log'}, optional
        Treatment for invalid floating-point operation.

    Returns
    -------
    old_settings : dict
        Dictionary containing the old settings.

    See also
    --------
    seterrcall : Set a callback function for the 'call' mode.
    geterr, geterrcall, errstate

    Notes
    -----
    The floating-point exceptions are defined in the IEEE 754 standard [1]:

    - Division by zero: infinite result obtained from finite numbers.
    - Overflow: result too large to be expressed.
    - Underflow: result so close to zero that some precision
      was lost.
    - Invalid operation: result is not an expressible number, typically
      indicates that a NaN was produced.

    .. [1] http://en.wikipedia.org/wiki/IEEE_754

    Examples
    --------
    >>> old_settings = np.seterr(all='ignore')  #seterr to known value
    >>> np.seterr(over='raise')
    {'over': 'ignore', 'divide': 'ignore', 'invalid': 'ignore',
     'under': 'ignore'}
    >>> np.seterr(**old_settings)  # reset to default
    {'over': 'raise', 'divide': 'ignore', 'invalid': 'ignore', 'under': 'ignore'}

    >>> np.int16(32000) * np.int16(3)
    30464
    >>> old_settings = np.seterr(all='warn', over='raise')
    >>> np.int16(32000) * np.int16(3)
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    FloatingPointError: overflow encountered in short_scalars

    >>> old_settings = np.seterr(all='print')
    >>> np.geterr()
    {'over': 'print', 'divide': 'print', 'invalid': 'print', 'under': 'print'}
    >>> np.int16(32000) * np.int16(3)
    Warning: overflow encountered in short_scalars
    30464

    """

    pyvals = umath.geterrobj()
    old = geterr()

    if divide is None: divide = all or old['divide']
    if over is None: over = all or old['over']
    if under is None: under = all or old['under']
    if invalid is None: invalid = all or old['invalid']

    maskvalue = ((_errdict[divide] << SHIFT_DIVIDEBYZERO) +
                 (_errdict[over] << SHIFT_OVERFLOW ) +
                 (_errdict[under] << SHIFT_UNDERFLOW) +
                 (_errdict[invalid] << SHIFT_INVALID))

    pyvals[1] = maskvalue
    umath.seterrobj(pyvals)
    return old


def geterr():
    """
    Get the current way of handling floating-point errors.

    Returns
    -------
    res : dict
        A dictionary with keys "divide", "over", "under", and "invalid",
        whose values are from the strings "ignore", "print", "log", "warn",
        "raise", and "call". The keys represent possible floating-point
        exceptions, and the values define how these exceptions are handled.

    See Also
    --------
    geterrcall, seterr, seterrcall

    Notes
    -----
    For complete documentation of the types of floating-point exceptions and
    treatment options, see `seterr`.

    Examples
    --------
    >>> np.geterr()
    {'over': 'warn', 'divide': 'warn', 'invalid': 'warn',
    'under': 'ignore'}
    >>> np.arange(3.) / np.arange(3.)
    array([ NaN,   1.,   1.])

    >>> oldsettings = np.seterr(all='warn', over='raise')
    >>> np.geterr()
    {'over': 'raise', 'divide': 'warn', 'invalid': 'warn', 'under': 'warn'}
    >>> np.arange(3.) / np.arange(3.)
    __main__:1: RuntimeWarning: invalid value encountered in divide
    array([ NaN,   1.,   1.])

    """
    maskvalue = umath.geterrobj()[1]
    mask = 7
    res = {}
    val = (maskvalue >> SHIFT_DIVIDEBYZERO) & mask
    res['divide'] = _errdict_rev[val]
    val = (maskvalue >> SHIFT_OVERFLOW) & mask
    res['over'] = _errdict_rev[val]
    val = (maskvalue >> SHIFT_UNDERFLOW) & mask
    res['under'] = _errdict_rev[val]
    val = (maskvalue >> SHIFT_INVALID) & mask
    res['invalid'] = _errdict_rev[val]
    return res

def setbufsize(size):
    """
    Set the size of the buffer used in ufuncs.

    Parameters
    ----------
    size : int
        Size of buffer.

    """
    if size > 10e6:
        raise ValueError("Buffer size, %s, is too big." % size)
    if size < 5:
        raise ValueError("Buffer size, %s, is too small." %size)
    if size % 16 != 0:
        raise ValueError("Buffer size, %s, is not a multiple of 16." %size)

    pyvals = umath.geterrobj()
    old = getbufsize()
    pyvals[0] = size
    umath.seterrobj(pyvals)
    return old

def getbufsize():
    """
    Return the size of the buffer used in ufuncs.

    Returns
    -------
    getbufsize : int
        Size of ufunc buffer in bytes.

    """
    return umath.geterrobj()[0]

def seterrcall(func):
    """
    Set the floating-point error callback function or log object.

    There are two ways to capture floating-point error messages.  The first
    is to set the error-handler to 'call', using `seterr`.  Then, set
    the function to call using this function.

    The second is to set the error-handler to 'log', using `seterr`.
    Floating-point errors then trigger a call to the 'write' method of
    the provided object.

    Parameters
    ----------
    func : callable f(err, flag) or object with write method
        Function to call upon floating-point errors ('call'-mode) or
        object whose 'write' method is used to log such message ('log'-mode).

        The call function takes two arguments. The first is the
        type of error (one of "divide", "over", "under", or "invalid"),
        and the second is the status flag.  The flag is a byte, whose
        least-significant bits indicate the status::

          [0 0 0 0 invalid over under invalid]

        In other words, ``flags = divide + 2*over + 4*under + 8*invalid``.

        If an object is provided, its write method should take one argument,
        a string.

    Returns
    -------
    h : callable, log instance or None
        The old error handler.

    See Also
    --------
    seterr, geterr, geterrcall

    Examples
    --------
    Callback upon error:

    >>> def err_handler(type, flag):
    ...     print "Floating point error (%s), with flag %s" % (type, flag)
    ...

    >>> saved_handler = np.seterrcall(err_handler)
    >>> save_err = np.seterr(all='call')

    >>> np.array([1, 2, 3]) / 0.0
    Floating point error (divide by zero), with flag 1
    array([ Inf,  Inf,  Inf])

    >>> np.seterrcall(saved_handler)
    <function err_handler at 0x...>
    >>> np.seterr(**save_err)
    {'over': 'call', 'divide': 'call', 'invalid': 'call', 'under': 'call'}

    Log error message:

    >>> class Log(object):
    ...     def write(self, msg):
    ...         print "LOG: %s" % msg
    ...

    >>> log = Log()
    >>> saved_handler = np.seterrcall(log)
    >>> save_err = np.seterr(all='log')

    >>> np.array([1, 2, 3]) / 0.0
    LOG: Warning: divide by zero encountered in divide
    <BLANKLINE>
    array([ Inf,  Inf,  Inf])

    >>> np.seterrcall(saved_handler)
    <__main__.Log object at 0x...>
    >>> np.seterr(**save_err)
    {'over': 'log', 'divide': 'log', 'invalid': 'log', 'under': 'log'}

    """
    if func is not None and not isinstance(func, collections.Callable):
        if not hasattr(func, 'write') or not isinstance(func.write, collections.Callable):
            raise ValueError("Only callable can be used as callback")
    pyvals = umath.geterrobj()
    old = geterrcall()
    pyvals[2] = func
    umath.seterrobj(pyvals)
    return old

def geterrcall():
    """
    Return the current callback function used on floating-point errors.

    When the error handling for a floating-point error (one of "divide",
    "over", "under", or "invalid") is set to 'call' or 'log', the function
    that is called or the log instance that is written to is returned by
    `geterrcall`. This function or log instance has been set with
    `seterrcall`.

    Returns
    -------
    errobj : callable, log instance or None
        The current error handler. If no handler was set through `seterrcall`,
        ``None`` is returned.

    See Also
    --------
    seterrcall, seterr, geterr

    Notes
    -----
    For complete documentation of the types of floating-point exceptions and
    treatment options, see `seterr`.

    Examples
    --------
    >>> np.geterrcall()  # we did not yet set a handler, returns None

    >>> oldsettings = np.seterr(all='call')
    >>> def err_handler(type, flag):
    ...     print "Floating point error (%s), with flag %s" % (type, flag)
    >>> oldhandler = np.seterrcall(err_handler)
    >>> np.array([1, 2, 3]) / 0.0
    Floating point error (divide by zero), with flag 1
    array([ Inf,  Inf,  Inf])

    >>> cur_handler = np.geterrcall()
    >>> cur_handler is err_handler
    True

    """
    return umath.geterrobj()[2]

class _unspecified(object):
    pass
_Unspecified = _unspecified()

class errstate(object):
    """
    errstate(**kwargs)

    Context manager for floating-point error handling.

    Using an instance of `errstate` as a context manager allows statements in
    that context to execute with a known error handling behavior. Upon entering
    the context the error handling is set with `seterr` and `seterrcall`, and
    upon exiting it is reset to what it was before.

    Parameters
    ----------
    kwargs : {divide, over, under, invalid}
        Keyword arguments. The valid keywords are the possible floating-point
        exceptions. Each keyword should have a string value that defines the
        treatment for the particular error. Possible values are
        {'ignore', 'warn', 'raise', 'call', 'print', 'log'}.

    See Also
    --------
    seterr, geterr, seterrcall, geterrcall

    Notes
    -----
    The ``with`` statement was introduced in Python 2.5, and can only be used
    there by importing it: ``from __future__ import with_statement``. In
    earlier Python versions the ``with`` statement is not available.

    For complete documentation of the types of floating-point exceptions and
    treatment options, see `seterr`.

    Examples
    --------
    >>> from __future__ import with_statement  # use 'with' in Python 2.5
    >>> olderr = np.seterr(all='ignore')  # Set error handling to known state.

    >>> np.arange(3) / 0.
    array([ NaN,  Inf,  Inf])
    >>> with np.errstate(divide='warn'):
    ...     np.arange(3) / 0.
    ...
    __main__:2: RuntimeWarning: divide by zero encountered in divide
    array([ NaN,  Inf,  Inf])

    >>> np.sqrt(-1)
    nan
    >>> with np.errstate(invalid='raise'):
    ...     np.sqrt(-1)
    Traceback (most recent call last):
      File "<stdin>", line 2, in <module>
    FloatingPointError: invalid value encountered in sqrt

    Outside the context the error handling behavior has not changed:

    >>> np.geterr()
    {'over': 'warn', 'divide': 'warn', 'invalid': 'warn',
    'under': 'ignore'}

    """
    # Note that we don't want to run the above doctests because they will fail
    # without a from __future__ import with_statement
    def __init__(self, **kwargs):
        self.call = kwargs.pop('call', _Unspecified)
        self.kwargs = kwargs

    def __enter__(self):
        self.oldstate = seterr(**self.kwargs)
        if self.call is not _Unspecified:
            self.oldcall = seterrcall(self.call)

    def __exit__(self, *exc_info):
        seterr(**self.oldstate)
        if self.call is not _Unspecified:
            seterrcall(self.oldcall)


def _setdef():
    defval = [UFUNC_BUFSIZE_DEFAULT, ERR_DEFAULT, None]
    umath.seterrobj(defval)

# set the default values
_setdef()

Inf = inf = infty = Infinity = PINF
nan = NaN = NAN
False_ = bool_(False)
True_ = bool_(True)

from .umath import *
from .numerictypes import *
from . import fromnumeric
from .fromnumeric import *
extend_all(fromnumeric)
extend_all(umath)
extend_all(numerictypes)
