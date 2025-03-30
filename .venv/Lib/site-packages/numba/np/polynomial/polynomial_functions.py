"""
Implementation of operations involving polynomials.
"""


import numpy as np
from numpy.polynomial import polynomial as poly
from numpy.polynomial import polyutils as pu

from numba import literal_unroll
from numba.core import types, errors
from numba.core.extending import overload
from numba.np.numpy_support import type_can_asarray, as_dtype, from_dtype


@overload(np.roots)
def roots_impl(p):

    # cast int vectors to float cf. numpy, this is a bit dicey as
    # the roots could be complex which will fail anyway
    ty = getattr(p, 'dtype', p)
    if isinstance(ty, types.Integer):
        cast_t = np.float64
    else:
        cast_t = as_dtype(ty)

    def roots_impl(p):
        # impl based on numpy:
        # https://github.com/numpy/numpy/blob/master/numpy/lib/polynomial.py

        if len(p.shape) != 1:
            raise ValueError("Input must be a 1d array.")

        non_zero = np.nonzero(p)[0]

        if len(non_zero) == 0:
            return np.zeros(0, dtype=cast_t)

        tz = len(p) - non_zero[-1] - 1

        # pull out the coeffs selecting between possible zero pads
        p = p[int(non_zero[0]):int(non_zero[-1]) + 1]

        n = len(p)
        if n > 1:
            # construct companion matrix, ensure fortran order
            # to give to eigvals, write to upper diag and then
            # transpose.
            A = np.diag(np.ones((n - 2,), cast_t), 1).T
            A[0, :] = -p[1:] / p[0]  # normalize
            roots = np.linalg.eigvals(A)
        else:
            roots = np.zeros(0, dtype=cast_t)

        # add in additional zeros on the end if needed
        if tz > 0:
            return np.hstack((roots, np.zeros(tz, dtype=cast_t)))
        else:
            return roots

    return roots_impl


@overload(pu.trimseq)
def polyutils_trimseq(seq):
    if not type_can_asarray(seq):
        msg = 'The argument "seq" must be array-like'
        raise errors.TypingError(msg)

    if isinstance(seq, types.BaseTuple):
        msg = 'Unsupported type %r for argument "seq"'
        raise errors.TypingError(msg % (seq))

    if np.ndim(seq) > 1:
        msg = 'Coefficient array is not 1-d'
        raise errors.NumbaValueError(msg)

    def impl(seq):
        if len(seq) == 0:
            return seq
        else:
            for i in range(len(seq) - 1, -1, -1):
                if seq[i] != 0:
                    break
            return seq[:i + 1]

    return impl


@overload(pu.as_series)
def polyutils_as_series(alist, trim=True):
    if not type_can_asarray(alist):
        msg = 'The argument "alist" must be array-like'
        raise errors.TypingError(msg)

    if not isinstance(trim, (bool, types.Boolean)):
        msg = 'The argument "trim" must be boolean'
        raise errors.TypingError(msg)

    res_dtype = np.float64

    tuple_input = isinstance(alist, types.BaseTuple)
    list_input = isinstance(alist, types.List)
    if tuple_input:
        if np.any(np.array([np.ndim(a) > 1 for a in alist])):
            raise errors.NumbaValueError("Coefficient array is not 1-d")

        res_dtype = _poly_result_dtype(*alist)

    elif list_input:
        dt = as_dtype(_get_list_type(alist))
        res_dtype = np.result_type(dt, np.float64)

    else:
        if np.ndim(alist) <= 2:
            res_dtype = np.result_type(res_dtype, as_dtype(alist.dtype))
        else:
            # If total dimension has ndim > 2, then coeff arrays are not 1D
            raise errors.NumbaValueError("Coefficient array is not 1-d")

    def impl(alist, trim=True):
        if tuple_input:
            arrays = []
            for item in literal_unroll(alist):
                arrays.append(np.atleast_1d(np.asarray(item)).astype(res_dtype))

        elif list_input:
            arrays = [np.atleast_1d(np.asarray(a)).astype(res_dtype)
                      for a in alist]

        else:
            alist_arr = np.asarray(alist)
            arrays = [np.atleast_1d(np.asarray(a)).astype(res_dtype)
                      for a in alist_arr]

        if min([a.size for a in arrays]) == 0:
            raise ValueError("Coefficient array is empty")

        if trim:
            arrays = [pu.trimseq(a) for a in arrays]

        ret = arrays
        return ret

    return impl


def _get_list_type(l):
    # A helper function that takes a list (possibly nested) and returns its
    # dtype. Returns a Numba type.
    dt = l.dtype
    if (not isinstance(dt, types.Number)) and type_can_asarray(dt):
        return _get_list_type(dt)
    else:
        return dt


def _poly_result_dtype(*args):
    # A helper function that takes a tuple of inputs and returns their result
    # dtype. Used for poly functions. Returns a NumPy dtype.
    res_dtype = np.float64
    for item in args:
        if isinstance(item, types.BaseTuple):
            s1 = item.types
        elif isinstance(item, types.List):
            s1 = [_get_list_type(item)]
        elif isinstance(item, types.Number):
            s1 = [item]
        elif isinstance(item, types.Array):
            s1 = [item.dtype]
        else:
            msg = 'Input dtype must be scalar'
            raise errors.TypingError(msg)

        try:
            l = [as_dtype(t) for t in s1]
            l.append(res_dtype)
            res_dtype = (np.result_type(*l))
        except errors.NumbaNotImplementedError:
            msg = 'Input dtype must be scalar.'
            raise errors.TypingError(msg)

    return from_dtype(res_dtype)


@overload(poly.polyadd)
def numpy_polyadd(c1, c2):
    if not type_can_asarray(c1):
        msg = 'The argument "c1" must be array-like'
        raise errors.TypingError(msg)

    if not type_can_asarray(c2):
        msg = 'The argument "c2" must be array-like'
        raise errors.TypingError(msg)

    def impl(c1, c2):
        arr1, arr2 = pu.as_series((c1, c2))
        diff = len(arr2) - len(arr1)
        if diff > 0:
            zr = np.zeros(diff)
            arr1 = np.concatenate((arr1, zr))
        if diff < 0:
            zr = np.zeros(-diff)
            arr2 = np.concatenate((arr2, zr))
        val = arr1 + arr2
        return pu.trimseq(val)

    return impl


@overload(poly.polysub)
def numpy_polysub(c1, c2):
    if not type_can_asarray(c1):
        msg = 'The argument "c1" must be array-like'
        raise errors.TypingError(msg)

    if not type_can_asarray(c2):
        msg = 'The argument "c2" must be array-like'
        raise errors.TypingError(msg)

    def impl(c1, c2):
        arr1, arr2 = pu.as_series((c1, c2))
        diff = len(arr2) - len(arr1)
        if diff > 0:
            zr = np.zeros(diff)
            arr1 = np.concatenate((arr1, zr))
        if diff < 0:
            zr = np.zeros(-diff)
            arr2 = np.concatenate((arr2, zr))
        val = arr1 - arr2
        return pu.trimseq(val)

    return impl


@overload(poly.polymul)
def numpy_polymul(c1, c2):
    if not type_can_asarray(c1):
        msg = 'The argument "c1" must be array-like'
        raise errors.TypingError(msg)

    if not type_can_asarray(c2):
        msg = 'The argument "c2" must be array-like'
        raise errors.TypingError(msg)

    def impl(c1, c2):
        arr1, arr2 = pu.as_series((c1, c2))
        val = np.convolve(arr1, arr2)
        return pu.trimseq(val)

    return impl


@overload(poly.polyval, prefer_literal=True)
def poly_polyval(x, c, tensor=True):
    if not type_can_asarray(x):
        msg = 'The argument "x" must be array-like'
        raise errors.TypingError(msg)

    if not type_can_asarray(c):
        msg = 'The argument "c" must be array-like'
        raise errors.TypingError(msg)

    if not isinstance(tensor, (bool, types.BooleanLiteral)):
        msg = 'The argument "tensor" must be boolean'
        raise errors.RequireLiteralValue(msg)

    res_dtype = _poly_result_dtype(c, x)

    # Simulate new_shape = (1,) * np.ndim(x) in the general case
    # If x is a number, new_shape is not used
    # If x is a tuple or a list, then it's 1d hence new_shape=(1,)
    x_nd_array = not isinstance(x, types.Number)
    new_shape = (1,)
    if isinstance(x, types.Array):
        # If x is a np.array, then take its dimension
        new_shape = (1,) * np.ndim(x)

    if isinstance(tensor, bool):
        tensor_arg = tensor
    else:
        tensor_arg = tensor.literal_value

    def impl(x, c, tensor=True):
        arr = np.asarray(c).astype(res_dtype)
        inputs = np.asarray(x).astype(res_dtype)
        if x_nd_array and tensor_arg:
            arr = arr.reshape(arr.shape + new_shape)

        l = len(arr)
        y = arr[l - 1] + inputs * 0

        for i in range(l - 1, 0, -1):
            y = arr[i - 1] + y * inputs

        return y

    return impl


@overload(poly.polyint)
def poly_polyint(c, m=1):

    if not type_can_asarray(c):
        msg = 'The argument "c" must be array-like'
        raise errors.TypingError(msg)

    if not isinstance(m, (int, types.Integer)):
        msg = 'The argument "m" must be an integer'
        raise errors.TypingError(msg)

    res_dtype = as_dtype(_poly_result_dtype(c))

    if not np.issubdtype(res_dtype, np.number):
        msg = f'Input dtype must be scalar. Found {res_dtype} instead'
        raise errors.TypingError(msg)

    is1D = ((np.ndim(c) == 1) or
            (isinstance(c, (types.List, types.BaseTuple))
             and isinstance(c.dtype, types.Number)))

    def impl(c, m=1):
        c = np.asarray(c).astype(res_dtype)
        cdt = c.dtype
        for i in range(m):
            n = len(c)

            tmp = np.empty((n + 1,) + c.shape[1:], dtype=cdt)
            tmp[0] = c[0] * 0
            tmp[1] = c[0]
            for j in range(1, n):
                tmp[j + 1] = c[j] / (j + 1)
            c = tmp
        if is1D:
            return pu.trimseq(c)
        else:
            return c

    return impl


@overload(poly.polydiv)
def numpy_polydiv(c1, c2):
    if not type_can_asarray(c1):
        msg = 'The argument "c1" must be array-like'
        raise errors.TypingError(msg)

    if not type_can_asarray(c2):
        msg = 'The argument "c2" must be array-like'
        raise errors.TypingError(msg)

    def impl(c1, c2):
        arr1, arr2 = pu.as_series((c1, c2))
        if arr2[-1] == 0:
            raise ZeroDivisionError()

        l1 = len(arr1)
        l2 = len(arr2)
        if l1 < l2:
            return arr1[:1] * 0, arr1
        elif l2 == 1:
            return arr1 / arr2[-1], arr1[:1] * 0
        else:
            dlen = l1 - l2
            scl = arr2[-1]
            arr2 = arr2[:-1] / scl
            i = dlen
            j = l1 - 1
            while i >= 0:
                arr1[i:j] -= arr2 * arr1[j]
                i -= 1
                j -= 1
            return arr1[j + 1:] / scl, pu.trimseq(arr1[:j + 1])

    return impl
