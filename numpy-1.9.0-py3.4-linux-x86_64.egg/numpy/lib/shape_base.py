from __future__ import division, absolute_import, print_function

import warnings

import numpy.core.numeric as _nx
from numpy.core.numeric import (
    asarray, zeros, outer, concatenate, isscalar, array, asanyarray
    )
from numpy.core.fromnumeric import product, reshape
from numpy.core import vstack, atleast_3d


__all__ = [
    'column_stack', 'row_stack', 'dstack', 'array_split', 'split',
    'hsplit', 'vsplit', 'dsplit', 'apply_over_axes', 'expand_dims',
    'apply_along_axis', 'kron', 'tile', 'get_array_wrap'
    ]


def apply_along_axis(func1d, axis, arr, *args, **kwargs):
    """
    Apply a function to 1-D slices along the given axis.

    Execute `func1d(a, *args)` where `func1d` operates on 1-D arrays and `a`
    is a 1-D slice of `arr` along `axis`.

    Parameters
    ----------
    func1d : function
        This function should accept 1-D arrays. It is applied to 1-D
        slices of `arr` along the specified axis.
    axis : integer
        Axis along which `arr` is sliced.
    arr : ndarray
        Input array.
    args : any
        Additional arguments to `func1d`.
    kwargs: any
        Additional named arguments to `func1d`.

        .. versionadded:: 1.9.0


    Returns
    -------
    apply_along_axis : ndarray
        The output array. The shape of `outarr` is identical to the shape of
        `arr`, except along the `axis` dimension, where the length of `outarr`
        is equal to the size of the return value of `func1d`.  If `func1d`
        returns a scalar `outarr` will have one fewer dimensions than `arr`.

    See Also
    --------
    apply_over_axes : Apply a function repeatedly over multiple axes.

    Examples
    --------
    >>> def my_func(a):
    ...     \"\"\"Average first and last element of a 1-D array\"\"\"
    ...     return (a[0] + a[-1]) * 0.5
    >>> b = np.array([[1,2,3], [4,5,6], [7,8,9]])
    >>> np.apply_along_axis(my_func, 0, b)
    array([ 4.,  5.,  6.])
    >>> np.apply_along_axis(my_func, 1, b)
    array([ 2.,  5.,  8.])

    For a function that doesn't return a scalar, the number of dimensions in
    `outarr` is the same as `arr`.

    >>> b = np.array([[8,1,7], [4,3,9], [5,2,6]])
    >>> np.apply_along_axis(sorted, 1, b)
    array([[1, 7, 8],
           [3, 4, 9],
           [2, 5, 6]])

    """
    arr = asarray(arr)
    nd = arr.ndim
    if axis < 0:
        axis += nd
    if (axis >= nd):
        raise ValueError("axis must be less than arr.ndim; axis=%d, rank=%d."
            % (axis, nd))
    ind = [0]*(nd-1)
    i = zeros(nd, 'O')
    indlist = list(range(nd))
    indlist.remove(axis)
    i[axis] = slice(None, None)
    outshape = asarray(arr.shape).take(indlist)
    i.put(indlist, ind)
    res = func1d(arr[tuple(i.tolist())], *args, **kwargs)
    #  if res is a number, then we have a smaller output array
    if isscalar(res):
        outarr = zeros(outshape, asarray(res).dtype)
        outarr[tuple(ind)] = res
        Ntot = product(outshape)
        k = 1
        while k < Ntot:
            # increment the index
            ind[-1] += 1
            n = -1
            while (ind[n] >= outshape[n]) and (n > (1-nd)):
                ind[n-1] += 1
                ind[n] = 0
                n -= 1
            i.put(indlist, ind)
            res = func1d(arr[tuple(i.tolist())], *args, **kwargs)
            outarr[tuple(ind)] = res
            k += 1
        return outarr
    else:
        Ntot = product(outshape)
        holdshape = outshape
        outshape = list(arr.shape)
        outshape[axis] = len(res)
        outarr = zeros(outshape, asarray(res).dtype)
        outarr[tuple(i.tolist())] = res
        k = 1
        while k < Ntot:
            # increment the index
            ind[-1] += 1
            n = -1
            while (ind[n] >= holdshape[n]) and (n > (1-nd)):
                ind[n-1] += 1
                ind[n] = 0
                n -= 1
            i.put(indlist, ind)
            res = func1d(arr[tuple(i.tolist())], *args, **kwargs)
            outarr[tuple(i.tolist())] = res
            k += 1
        return outarr


def apply_over_axes(func, a, axes):
    """
    Apply a function repeatedly over multiple axes.

    `func` is called as `res = func(a, axis)`, where `axis` is the first
    element of `axes`.  The result `res` of the function call must have
    either the same dimensions as `a` or one less dimension.  If `res`
    has one less dimension than `a`, a dimension is inserted before
    `axis`.  The call to `func` is then repeated for each axis in `axes`,
    with `res` as the first argument.

    Parameters
    ----------
    func : function
        This function must take two arguments, `func(a, axis)`.
    a : array_like
        Input array.
    axes : array_like
        Axes over which `func` is applied; the elements must be integers.

    Returns
    -------
    apply_over_axis : ndarray
        The output array.  The number of dimensions is the same as `a`,
        but the shape can be different.  This depends on whether `func`
        changes the shape of its output with respect to its input.

    See Also
    --------
    apply_along_axis :
        Apply a function to 1-D slices of an array along the given axis.

    Notes
    ------
    This function is equivalent to tuple axis arguments to reorderable ufuncs
    with keepdims=True. Tuple axis arguments to ufuncs have been availabe since
    version 1.7.0.

    Examples
    --------
    >>> a = np.arange(24).reshape(2,3,4)
    >>> a
    array([[[ 0,  1,  2,  3],
            [ 4,  5,  6,  7],
            [ 8,  9, 10, 11]],
           [[12, 13, 14, 15],
            [16, 17, 18, 19],
            [20, 21, 22, 23]]])

    Sum over axes 0 and 2. The result has same number of dimensions
    as the original array:

    >>> np.apply_over_axes(np.sum, a, [0,2])
    array([[[ 60],
            [ 92],
            [124]]])

    Tuple axis arguments to ufuncs are equivalent:

    >>> np.sum(a, axis=(0,2), keepdims=True)
    array([[[ 60],
            [ 92],
            [124]]])

    """
    val = asarray(a)
    N = a.ndim
    if array(axes).ndim == 0:
        axes = (axes,)
    for axis in axes:
        if axis < 0:
            axis = N + axis
        args = (val, axis)
        res = func(*args)
        if res.ndim == val.ndim:
            val = res
        else:
            res = expand_dims(res, axis)
            if res.ndim == val.ndim:
                val = res
            else:
                raise ValueError("function is not returning "
                        "an array of the correct shape")
    return val

def expand_dims(a, axis):
    """
    Expand the shape of an array.

    Insert a new axis, corresponding to a given position in the array shape.

    Parameters
    ----------
    a : array_like
        Input array.
    axis : int
        Position (amongst axes) where new axis is to be inserted.

    Returns
    -------
    res : ndarray
        Output array. The number of dimensions is one greater than that of
        the input array.

    See Also
    --------
    doc.indexing, atleast_1d, atleast_2d, atleast_3d

    Examples
    --------
    >>> x = np.array([1,2])
    >>> x.shape
    (2,)

    The following is equivalent to ``x[np.newaxis,:]`` or ``x[np.newaxis]``:

    >>> y = np.expand_dims(x, axis=0)
    >>> y
    array([[1, 2]])
    >>> y.shape
    (1, 2)

    >>> y = np.expand_dims(x, axis=1)  # Equivalent to x[:,newaxis]
    >>> y
    array([[1],
           [2]])
    >>> y.shape
    (2, 1)

    Note that some examples may use ``None`` instead of ``np.newaxis``.  These
    are the same objects:

    >>> np.newaxis is None
    True

    """
    a = asarray(a)
    shape = a.shape
    if axis < 0:
        axis = axis + len(shape) + 1
    return a.reshape(shape[:axis] + (1,) + shape[axis:])

row_stack = vstack

def column_stack(tup):
    """
    Stack 1-D arrays as columns into a 2-D array.

    Take a sequence of 1-D arrays and stack them as columns
    to make a single 2-D array. 2-D arrays are stacked as-is,
    just like with `hstack`.  1-D arrays are turned into 2-D columns
    first.

    Parameters
    ----------
    tup : sequence of 1-D or 2-D arrays.
        Arrays to stack. All of them must have the same first dimension.

    Returns
    -------
    stacked : 2-D array
        The array formed by stacking the given arrays.

    See Also
    --------
    hstack, vstack, concatenate

    Examples
    --------
    >>> a = np.array((1,2,3))
    >>> b = np.array((2,3,4))
    >>> np.column_stack((a,b))
    array([[1, 2],
           [2, 3],
           [3, 4]])

    """
    arrays = []
    for v in tup:
        arr = array(v, copy=False, subok=True)
        if arr.ndim < 2:
            arr = array(arr, copy=False, subok=True, ndmin=2).T
        arrays.append(arr)
    return _nx.concatenate(arrays, 1)

def dstack(tup):
    """
    Stack arrays in sequence depth wise (along third axis).

    Takes a sequence of arrays and stack them along the third axis
    to make a single array. Rebuilds arrays divided by `dsplit`.
    This is a simple way to stack 2D arrays (images) into a single
    3D array for processing.

    Parameters
    ----------
    tup : sequence of arrays
        Arrays to stack. All of them must have the same shape along all
        but the third axis.

    Returns
    -------
    stacked : ndarray
        The array formed by stacking the given arrays.

    See Also
    --------
    vstack : Stack along first axis.
    hstack : Stack along second axis.
    concatenate : Join arrays.
    dsplit : Split array along third axis.

    Notes
    -----
    Equivalent to ``np.concatenate(tup, axis=2)``.

    Examples
    --------
    >>> a = np.array((1,2,3))
    >>> b = np.array((2,3,4))
    >>> np.dstack((a,b))
    array([[[1, 2],
            [2, 3],
            [3, 4]]])

    >>> a = np.array([[1],[2],[3]])
    >>> b = np.array([[2],[3],[4]])
    >>> np.dstack((a,b))
    array([[[1, 2]],
           [[2, 3]],
           [[3, 4]]])

    """
    return _nx.concatenate([atleast_3d(_m) for _m in tup], 2)

def _replace_zero_by_x_arrays(sub_arys):
    for i in range(len(sub_arys)):
        if len(_nx.shape(sub_arys[i])) == 0:
            sub_arys[i] = _nx.empty(0, dtype=sub_arys[i].dtype)
        elif _nx.sometrue(_nx.equal(_nx.shape(sub_arys[i]), 0)):
            sub_arys[i] = _nx.empty(0, dtype=sub_arys[i].dtype)
    return sub_arys

def array_split(ary, indices_or_sections, axis=0):
    """
    Split an array into multiple sub-arrays.

    Please refer to the ``split`` documentation.  The only difference
    between these functions is that ``array_split`` allows
    `indices_or_sections` to be an integer that does *not* equally
    divide the axis.

    See Also
    --------
    split : Split array into multiple sub-arrays of equal size.

    Examples
    --------
    >>> x = np.arange(8.0)
    >>> np.array_split(x, 3)
        [array([ 0.,  1.,  2.]), array([ 3.,  4.,  5.]), array([ 6.,  7.])]

    """
    try:
        Ntotal = ary.shape[axis]
    except AttributeError:
        Ntotal = len(ary)
    try:
        # handle scalar case.
        Nsections = len(indices_or_sections) + 1
        div_points = [0] + list(indices_or_sections) + [Ntotal]
    except TypeError:
        # indices_or_sections is a scalar, not an array.
        Nsections = int(indices_or_sections)
        if Nsections <= 0:
            raise ValueError('number sections must be larger than 0.')
        Neach_section, extras = divmod(Ntotal, Nsections)
        section_sizes = ([0] +
                         extras * [Neach_section+1] +
                         (Nsections-extras) * [Neach_section])
        div_points = _nx.array(section_sizes).cumsum()

    sub_arys = []
    sary = _nx.swapaxes(ary, axis, 0)
    for i in range(Nsections):
        st = div_points[i]
        end = div_points[i + 1]
        sub_arys.append(_nx.swapaxes(sary[st:end], axis, 0))

    # This "kludge" was introduced here to replace arrays shaped (0, 10)
    # or similar with an array shaped (0,).
    # There seems no need for this, so give a FutureWarning to remove later.
    if sub_arys[-1].size == 0 and sub_arys[-1].ndim != 1:
        warnings.warn("in the future np.array_split will retain the shape of "
                      "arrays with a zero size, instead of replacing them by "
                      "`array([])`, which always has a shape of (0,).",
                      FutureWarning)
        sub_arys = _replace_zero_by_x_arrays(sub_arys)

    return sub_arys

def split(ary,indices_or_sections,axis=0):
    """
    Split an array into multiple sub-arrays.

    Parameters
    ----------
    ary : ndarray
        Array to be divided into sub-arrays.
    indices_or_sections : int or 1-D array
        If `indices_or_sections` is an integer, N, the array will be divided
        into N equal arrays along `axis`.  If such a split is not possible,
        an error is raised.

        If `indices_or_sections` is a 1-D array of sorted integers, the entries
        indicate where along `axis` the array is split.  For example,
        ``[2, 3]`` would, for ``axis=0``, result in

          - ary[:2]
          - ary[2:3]
          - ary[3:]

        If an index exceeds the dimension of the array along `axis`,
        an empty sub-array is returned correspondingly.
    axis : int, optional
        The axis along which to split, default is 0.

    Returns
    -------
    sub-arrays : list of ndarrays
        A list of sub-arrays.

    Raises
    ------
    ValueError
        If `indices_or_sections` is given as an integer, but
        a split does not result in equal division.

    See Also
    --------
    array_split : Split an array into multiple sub-arrays of equal or
                  near-equal size.  Does not raise an exception if
                  an equal division cannot be made.
    hsplit : Split array into multiple sub-arrays horizontally (column-wise).
    vsplit : Split array into multiple sub-arrays vertically (row wise).
    dsplit : Split array into multiple sub-arrays along the 3rd axis (depth).
    concatenate : Join arrays together.
    hstack : Stack arrays in sequence horizontally (column wise).
    vstack : Stack arrays in sequence vertically (row wise).
    dstack : Stack arrays in sequence depth wise (along third dimension).

    Examples
    --------
    >>> x = np.arange(9.0)
    >>> np.split(x, 3)
    [array([ 0.,  1.,  2.]), array([ 3.,  4.,  5.]), array([ 6.,  7.,  8.])]

    >>> x = np.arange(8.0)
    >>> np.split(x, [3, 5, 6, 10])
    [array([ 0.,  1.,  2.]),
     array([ 3.,  4.]),
     array([ 5.]),
     array([ 6.,  7.]),
     array([], dtype=float64)]

    """
    try:
        len(indices_or_sections)
    except TypeError:
        sections = indices_or_sections
        N = ary.shape[axis]
        if N % sections:
            raise ValueError(
                'array split does not result in an equal division')
    res = array_split(ary, indices_or_sections, axis)
    return res

def hsplit(ary, indices_or_sections):
    """
    Split an array into multiple sub-arrays horizontally (column-wise).

    Please refer to the `split` documentation.  `hsplit` is equivalent
    to `split` with ``axis=1``, the array is always split along the second
    axis regardless of the array dimension.

    See Also
    --------
    split : Split an array into multiple sub-arrays of equal size.

    Examples
    --------
    >>> x = np.arange(16.0).reshape(4, 4)
    >>> x
    array([[  0.,   1.,   2.,   3.],
           [  4.,   5.,   6.,   7.],
           [  8.,   9.,  10.,  11.],
           [ 12.,  13.,  14.,  15.]])
    >>> np.hsplit(x, 2)
    [array([[  0.,   1.],
           [  4.,   5.],
           [  8.,   9.],
           [ 12.,  13.]]),
     array([[  2.,   3.],
           [  6.,   7.],
           [ 10.,  11.],
           [ 14.,  15.]])]
    >>> np.hsplit(x, np.array([3, 6]))
    [array([[  0.,   1.,   2.],
           [  4.,   5.,   6.],
           [  8.,   9.,  10.],
           [ 12.,  13.,  14.]]),
     array([[  3.],
           [  7.],
           [ 11.],
           [ 15.]]),
     array([], dtype=float64)]

    With a higher dimensional array the split is still along the second axis.

    >>> x = np.arange(8.0).reshape(2, 2, 2)
    >>> x
    array([[[ 0.,  1.],
            [ 2.,  3.]],
           [[ 4.,  5.],
            [ 6.,  7.]]])
    >>> np.hsplit(x, 2)
    [array([[[ 0.,  1.]],
           [[ 4.,  5.]]]),
     array([[[ 2.,  3.]],
           [[ 6.,  7.]]])]

    """
    if len(_nx.shape(ary)) == 0:
        raise ValueError('hsplit only works on arrays of 1 or more dimensions')
    if len(ary.shape) > 1:
        return split(ary, indices_or_sections, 1)
    else:
        return split(ary, indices_or_sections, 0)

def vsplit(ary, indices_or_sections):
    """
    Split an array into multiple sub-arrays vertically (row-wise).

    Please refer to the ``split`` documentation.  ``vsplit`` is equivalent
    to ``split`` with `axis=0` (default), the array is always split along the
    first axis regardless of the array dimension.

    See Also
    --------
    split : Split an array into multiple sub-arrays of equal size.

    Examples
    --------
    >>> x = np.arange(16.0).reshape(4, 4)
    >>> x
    array([[  0.,   1.,   2.,   3.],
           [  4.,   5.,   6.,   7.],
           [  8.,   9.,  10.,  11.],
           [ 12.,  13.,  14.,  15.]])
    >>> np.vsplit(x, 2)
    [array([[ 0.,  1.,  2.,  3.],
           [ 4.,  5.,  6.,  7.]]),
     array([[  8.,   9.,  10.,  11.],
           [ 12.,  13.,  14.,  15.]])]
    >>> np.vsplit(x, np.array([3, 6]))
    [array([[  0.,   1.,   2.,   3.],
           [  4.,   5.,   6.,   7.],
           [  8.,   9.,  10.,  11.]]),
     array([[ 12.,  13.,  14.,  15.]]),
     array([], dtype=float64)]

    With a higher dimensional array the split is still along the first axis.

    >>> x = np.arange(8.0).reshape(2, 2, 2)
    >>> x
    array([[[ 0.,  1.],
            [ 2.,  3.]],
           [[ 4.,  5.],
            [ 6.,  7.]]])
    >>> np.vsplit(x, 2)
    [array([[[ 0.,  1.],
            [ 2.,  3.]]]),
     array([[[ 4.,  5.],
            [ 6.,  7.]]])]

    """
    if len(_nx.shape(ary)) < 2:
        raise ValueError('vsplit only works on arrays of 2 or more dimensions')
    return split(ary, indices_or_sections, 0)

def dsplit(ary, indices_or_sections):
    """
    Split array into multiple sub-arrays along the 3rd axis (depth).

    Please refer to the `split` documentation.  `dsplit` is equivalent
    to `split` with ``axis=2``, the array is always split along the third
    axis provided the array dimension is greater than or equal to 3.

    See Also
    --------
    split : Split an array into multiple sub-arrays of equal size.

    Examples
    --------
    >>> x = np.arange(16.0).reshape(2, 2, 4)
    >>> x
    array([[[  0.,   1.,   2.,   3.],
            [  4.,   5.,   6.,   7.]],
           [[  8.,   9.,  10.,  11.],
            [ 12.,  13.,  14.,  15.]]])
    >>> np.dsplit(x, 2)
    [array([[[  0.,   1.],
            [  4.,   5.]],
           [[  8.,   9.],
            [ 12.,  13.]]]),
     array([[[  2.,   3.],
            [  6.,   7.]],
           [[ 10.,  11.],
            [ 14.,  15.]]])]
    >>> np.dsplit(x, np.array([3, 6]))
    [array([[[  0.,   1.,   2.],
            [  4.,   5.,   6.]],
           [[  8.,   9.,  10.],
            [ 12.,  13.,  14.]]]),
     array([[[  3.],
            [  7.]],
           [[ 11.],
            [ 15.]]]),
     array([], dtype=float64)]

    """
    if len(_nx.shape(ary)) < 3:
        raise ValueError('dsplit only works on arrays of 3 or more dimensions')
    return split(ary, indices_or_sections, 2)

def get_array_prepare(*args):
    """Find the wrapper for the array with the highest priority.

    In case of ties, leftmost wins. If no wrapper is found, return None
    """
    wrappers = sorted((getattr(x, '__array_priority__', 0), -i,
                 x.__array_prepare__) for i, x in enumerate(args)
                                   if hasattr(x, '__array_prepare__'))
    if wrappers:
        return wrappers[-1][-1]
    return None

def get_array_wrap(*args):
    """Find the wrapper for the array with the highest priority.

    In case of ties, leftmost wins. If no wrapper is found, return None
    """
    wrappers = sorted((getattr(x, '__array_priority__', 0), -i,
                 x.__array_wrap__) for i, x in enumerate(args)
                                   if hasattr(x, '__array_wrap__'))
    if wrappers:
        return wrappers[-1][-1]
    return None

def kron(a, b):
    """
    Kronecker product of two arrays.

    Computes the Kronecker product, a composite array made of blocks of the
    second array scaled by the first.

    Parameters
    ----------
    a, b : array_like

    Returns
    -------
    out : ndarray

    See Also
    --------
    outer : The outer product

    Notes
    -----
    The function assumes that the number of dimenensions of `a` and `b`
    are the same, if necessary prepending the smallest with ones.
    If `a.shape = (r0,r1,..,rN)` and `b.shape = (s0,s1,...,sN)`,
    the Kronecker product has shape `(r0*s0, r1*s1, ..., rN*SN)`.
    The elements are products of elements from `a` and `b`, organized
    explicitly by::

        kron(a,b)[k0,k1,...,kN] = a[i0,i1,...,iN] * b[j0,j1,...,jN]

    where::

        kt = it * st + jt,  t = 0,...,N

    In the common 2-D case (N=1), the block structure can be visualized::

        [[ a[0,0]*b,   a[0,1]*b,  ... , a[0,-1]*b  ],
         [  ...                              ...   ],
         [ a[-1,0]*b,  a[-1,1]*b, ... , a[-1,-1]*b ]]


    Examples
    --------
    >>> np.kron([1,10,100], [5,6,7])
    array([  5,   6,   7,  50,  60,  70, 500, 600, 700])
    >>> np.kron([5,6,7], [1,10,100])
    array([  5,  50, 500,   6,  60, 600,   7,  70, 700])

    >>> np.kron(np.eye(2), np.ones((2,2)))
    array([[ 1.,  1.,  0.,  0.],
           [ 1.,  1.,  0.,  0.],
           [ 0.,  0.,  1.,  1.],
           [ 0.,  0.,  1.,  1.]])

    >>> a = np.arange(100).reshape((2,5,2,5))
    >>> b = np.arange(24).reshape((2,3,4))
    >>> c = np.kron(a,b)
    >>> c.shape
    (2, 10, 6, 20)
    >>> I = (1,3,0,2)
    >>> J = (0,2,1)
    >>> J1 = (0,) + J             # extend to ndim=4
    >>> S1 = (1,) + b.shape
    >>> K = tuple(np.array(I) * np.array(S1) + np.array(J1))
    >>> c[K] == a[I]*b[J]
    True

    """
    b = asanyarray(b)
    a = array(a, copy=False, subok=True, ndmin=b.ndim)
    ndb, nda = b.ndim, a.ndim
    if (nda == 0 or ndb == 0):
        return _nx.multiply(a, b)
    as_ = a.shape
    bs = b.shape
    if not a.flags.contiguous:
        a = reshape(a, as_)
    if not b.flags.contiguous:
        b = reshape(b, bs)
    nd = ndb
    if (ndb != nda):
        if (ndb > nda):
            as_ = (1,)*(ndb-nda) + as_
        else:
            bs = (1,)*(nda-ndb) + bs
            nd = nda
    result = outer(a, b).reshape(as_+bs)
    axis = nd-1
    for _ in range(nd):
        result = concatenate(result, axis=axis)
    wrapper = get_array_prepare(a, b)
    if wrapper is not None:
        result = wrapper(result)
    wrapper = get_array_wrap(a, b)
    if wrapper is not None:
        result = wrapper(result)
    return result


def tile(A, reps):
    """
    Construct an array by repeating A the number of times given by reps.

    If `reps` has length ``d``, the result will have dimension of
    ``max(d, A.ndim)``.

    If ``A.ndim < d``, `A` is promoted to be d-dimensional by prepending new
    axes. So a shape (3,) array is promoted to (1, 3) for 2-D replication,
    or shape (1, 1, 3) for 3-D replication. If this is not the desired
    behavior, promote `A` to d-dimensions manually before calling this
    function.

    If ``A.ndim > d``, `reps` is promoted to `A`.ndim by pre-pending 1's to it.
    Thus for an `A` of shape (2, 3, 4, 5), a `reps` of (2, 2) is treated as
    (1, 1, 2, 2).

    Parameters
    ----------
    A : array_like
        The input array.
    reps : array_like
        The number of repetitions of `A` along each axis.

    Returns
    -------
    c : ndarray
        The tiled output array.

    See Also
    --------
    repeat : Repeat elements of an array.

    Examples
    --------
    >>> a = np.array([0, 1, 2])
    >>> np.tile(a, 2)
    array([0, 1, 2, 0, 1, 2])
    >>> np.tile(a, (2, 2))
    array([[0, 1, 2, 0, 1, 2],
           [0, 1, 2, 0, 1, 2]])
    >>> np.tile(a, (2, 1, 2))
    array([[[0, 1, 2, 0, 1, 2]],
           [[0, 1, 2, 0, 1, 2]]])

    >>> b = np.array([[1, 2], [3, 4]])
    >>> np.tile(b, 2)
    array([[1, 2, 1, 2],
           [3, 4, 3, 4]])
    >>> np.tile(b, (2, 1))
    array([[1, 2],
           [3, 4],
           [1, 2],
           [3, 4]])

    """
    try:
        tup = tuple(reps)
    except TypeError:
        tup = (reps,)
    d = len(tup)
    c = _nx.array(A, copy=False, subok=True, ndmin=d)
    shape = list(c.shape)
    n = max(c.size, 1)
    if (d < c.ndim):
        tup = (1,)*(c.ndim-d) + tup
    for i, nrep in enumerate(tup):
        if nrep != 1:
            c = c.reshape(-1, n).repeat(nrep, 0)
        dim_in = shape[i]
        dim_out = dim_in*nrep
        shape[i] = dim_out
        n //= max(dim_in, 1)
    return c.reshape(shape)
