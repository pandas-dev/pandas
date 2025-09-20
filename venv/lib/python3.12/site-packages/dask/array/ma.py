from __future__ import annotations

from functools import wraps

import numpy as np

from dask.array import chunk
from dask.array.core import asanyarray, blockwise, elemwise, map_blocks
from dask.array.reductions import reduction
from dask.array.routines import _average
from dask.array.routines import nonzero as _nonzero
from dask.tokenize import normalize_token
from dask.utils import derived_from


@normalize_token.register(np.ma.masked_array)
def normalize_masked_array(x):
    data = normalize_token(x.data)
    mask = normalize_token(x.mask)
    fill_value = normalize_token(x.fill_value)
    return (data, mask, fill_value)


@derived_from(np.ma)
def filled(a, fill_value=None):
    a = asanyarray(a)
    return a.map_blocks(np.ma.filled, fill_value=fill_value)


def _wrap_masked(f):
    @wraps(f)
    def _(a, value):
        a = asanyarray(a)
        value = asanyarray(value)
        ainds = tuple(range(a.ndim))[::-1]
        vinds = tuple(range(value.ndim))[::-1]
        oinds = max(ainds, vinds, key=len)
        return blockwise(f, oinds, a, ainds, value, vinds, dtype=a.dtype)

    return _


masked_greater = _wrap_masked(np.ma.masked_greater)
masked_greater_equal = _wrap_masked(np.ma.masked_greater_equal)
masked_less = _wrap_masked(np.ma.masked_less)
masked_less_equal = _wrap_masked(np.ma.masked_less_equal)
masked_not_equal = _wrap_masked(np.ma.masked_not_equal)


@derived_from(np.ma)
def masked_equal(a, value):
    a = asanyarray(a)
    if getattr(value, "shape", ()):
        raise ValueError("da.ma.masked_equal doesn't support array `value`s")
    inds = tuple(range(a.ndim))
    return blockwise(np.ma.masked_equal, inds, a, inds, value, (), dtype=a.dtype)


@derived_from(np.ma)
def masked_invalid(a):
    return asanyarray(a).map_blocks(np.ma.masked_invalid)


@derived_from(np.ma)
def masked_inside(x, v1, v2):
    x = asanyarray(x)
    return x.map_blocks(np.ma.masked_inside, v1, v2)


@derived_from(np.ma)
def masked_outside(x, v1, v2):
    x = asanyarray(x)
    return x.map_blocks(np.ma.masked_outside, v1, v2)


@derived_from(np.ma)
def masked_where(condition, a):
    cshape = getattr(condition, "shape", ())
    if cshape and cshape != a.shape:
        raise IndexError(
            "Inconsistent shape between the condition and the "
            "input (got %s and %s)" % (cshape, a.shape)
        )
    condition = asanyarray(condition)
    a = asanyarray(a)
    ainds = tuple(range(a.ndim))
    cinds = tuple(range(condition.ndim))
    return blockwise(
        np.ma.masked_where, ainds, condition, cinds, a, ainds, dtype=a.dtype
    )


@derived_from(np.ma)
def masked_values(x, value, rtol=1e-05, atol=1e-08, shrink=True):
    x = asanyarray(x)
    if getattr(value, "shape", ()):
        raise ValueError("da.ma.masked_values doesn't support array `value`s")
    return map_blocks(
        np.ma.masked_values, x, value, rtol=rtol, atol=atol, shrink=shrink
    )


@derived_from(np.ma)
def fix_invalid(a, fill_value=None):
    a = asanyarray(a)
    return a.map_blocks(np.ma.fix_invalid, fill_value=fill_value)


@derived_from(np.ma)
def getdata(a):
    a = asanyarray(a)
    return a.map_blocks(np.ma.getdata)


@derived_from(np.ma)
def getmaskarray(a):
    a = asanyarray(a)
    return a.map_blocks(np.ma.getmaskarray)


def _masked_array(data, mask=np.ma.nomask, masked_dtype=None, **kwargs):
    if "chunks" in kwargs:
        del kwargs["chunks"]  # A Dask kwarg, not NumPy.
    return np.ma.masked_array(data, mask=mask, dtype=masked_dtype, **kwargs)


@derived_from(np.ma)
def masked_array(data, mask=np.ma.nomask, fill_value=None, **kwargs):
    data = asanyarray(data)
    inds = tuple(range(data.ndim))
    arginds = [inds, data, inds]

    if getattr(fill_value, "shape", ()):
        raise ValueError("non-scalar fill_value not supported")
    kwargs["fill_value"] = fill_value

    if mask is not np.ma.nomask:
        mask = asanyarray(mask)
        if mask.size == 1:
            mask = mask.reshape((1,) * data.ndim)
        elif data.shape != mask.shape:
            raise np.ma.MaskError(
                "Mask and data not compatible: data shape "
                "is %s, and mask shape is "
                "%s." % (repr(data.shape), repr(mask.shape))
            )
        arginds.extend([mask, inds])

    if "dtype" in kwargs:
        kwargs["masked_dtype"] = kwargs["dtype"]
    else:
        kwargs["dtype"] = data.dtype

    return blockwise(_masked_array, *arginds, **kwargs)


def _set_fill_value(x, fill_value):
    if isinstance(x, np.ma.masked_array):
        x = x.copy()
        np.ma.set_fill_value(x, fill_value=fill_value)
    return x


@derived_from(np.ma)
def set_fill_value(a, fill_value):
    a = asanyarray(a)
    if getattr(fill_value, "shape", ()):
        raise ValueError("da.ma.set_fill_value doesn't support array `value`s")
    fill_value = np.ma.core._check_fill_value(fill_value, a.dtype)
    res = a.map_blocks(_set_fill_value, fill_value)
    a.dask = res.dask
    a._name = res.name


@derived_from(np.ma)
def average(a, axis=None, weights=None, returned=False, keepdims=False):
    return _average(a, axis, weights, returned, is_masked=True, keepdims=keepdims)


def _chunk_count(x, axis=None, keepdims=None):
    return np.ma.count(x, axis=axis, keepdims=keepdims)


@derived_from(np.ma)
def count(a, axis=None, keepdims=False, split_every=None):
    return reduction(
        a,
        _chunk_count,
        chunk.sum,
        axis=axis,
        keepdims=keepdims,
        dtype=np.intp,
        split_every=split_every,
        out=None,
    )


@derived_from(np.ma.core)
def ones_like(a, **kwargs):
    a = asanyarray(a)
    return a.map_blocks(np.ma.core.ones_like, **kwargs)


@derived_from(np.ma.core)
def zeros_like(a, **kwargs):
    a = asanyarray(a)
    return a.map_blocks(np.ma.core.zeros_like, **kwargs)


@derived_from(np.ma.core)
def empty_like(a, **kwargs):
    a = asanyarray(a)
    return a.map_blocks(np.ma.core.empty_like, **kwargs)


@derived_from(np.ma.core)
def nonzero(a):
    return _nonzero(getdata(a) * ~getmaskarray(a))


@derived_from(np.ma.core)
def where(condition, x=None, y=None):
    if (x is None) != (y is None):
        raise ValueError("either both or neither of x and y should be given")
    if (x is None) and (y is None):
        return nonzero(condition)
    else:
        return elemwise(np.ma.where, condition, x, y)
