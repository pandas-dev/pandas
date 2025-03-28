from __future__ import annotations

import warnings

from dask._compatibility import import_optional_dependency

import_optional_dependency("numpy")
import numpy as np
from packaging.version import Version

from dask.utils import derived_from

_np_version = Version(np.__version__)
NUMPY_GE_125 = _np_version.release >= (1, 25)
NUMPY_GE_200 = _np_version.release >= (2, 0)
NUMPY_GE_210 = _np_version.release >= (2, 1)


if NUMPY_GE_200:
    from numpy.exceptions import AxisError, ComplexWarning  # noqa: F401
    from numpy.lib.array_utils import normalize_axis_index, normalize_axis_tuple
else:
    from numpy import (  # type: ignore[no-redef, attr-defined] # noqa: F401
        AxisError,
        ComplexWarning,
    )
    from numpy.core.numeric import normalize_axis_index  # type: ignore[no-redef]
    from numpy.core.numeric import normalize_axis_tuple  # type: ignore[no-redef]


# Taken from scikit-learn:
# https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/utils/fixes.py#L84
try:
    with warnings.catch_warnings():
        if (
            not np.allclose(
                np.divide(0.4, 1, casting="unsafe"),
                np.divide(0.4, 1, casting="unsafe", dtype=float),
            )
            or not np.allclose(np.divide(1, 0.5, dtype="i8"), 2)
            or not np.allclose(np.divide(0.4, 1), 0.4)
        ):
            raise TypeError(
                "Divide not working with dtype: "
                "https://github.com/numpy/numpy/issues/3484"
            )
        divide = np.divide
        ma_divide = np.ma.divide

except TypeError:
    # Divide with dtype doesn't work on Python 3
    def divide(x1, x2, out=None, dtype=None):  # type: ignore
        """Implementation of numpy.divide that works with dtype kwarg.

        Temporary compatibility fix for a bug in numpy's version. See
        https://github.com/numpy/numpy/issues/3484 for the relevant issue."""
        x = np.divide(x1, x2, out)
        if dtype is not None:
            x = x.astype(dtype)
        return x

    ma_divide = np.ma.core._DomainedBinaryOperation(  # type: ignore
        divide, np.ma.core._DomainSafeDivide(), 0, 1  # type: ignore
    )


class _Recurser:
    """
    Utility class for recursing over nested iterables
    """

    # This was copied almost verbatim from numpy.core.shape_base._Recurser
    # See numpy license at https://github.com/numpy/numpy/blob/master/LICENSE.txt
    # or NUMPY_LICENSE.txt within this directory

    def __init__(self, recurse_if):
        self.recurse_if = recurse_if

    def map_reduce(
        self,
        x,
        f_map=lambda x, **kwargs: x,
        f_reduce=lambda x, **kwargs: x,
        f_kwargs=lambda **kwargs: kwargs,
        **kwargs,
    ):
        """
        Iterate over the nested list, applying:
        * ``f_map`` (T -> U) to items
        * ``f_reduce`` (Iterable[U] -> U) to mapped items

        For instance, ``map_reduce([[1, 2], 3, 4])`` is::

            f_reduce([
              f_reduce([
                f_map(1),
                f_map(2)
              ]),
              f_map(3),
              f_map(4)
            ]])


        State can be passed down through the calls with `f_kwargs`,
        to iterables of mapped items. When kwargs are passed, as in
        ``map_reduce([[1, 2], 3, 4], **kw)``, this becomes::

            kw1 = f_kwargs(**kw)
            kw2 = f_kwargs(**kw1)
            f_reduce([
              f_reduce([
                f_map(1), **kw2)
                f_map(2,  **kw2)
              ],      **kw1),
              f_map(3, **kw1),
              f_map(4, **kw1)
            ]],     **kw)
        """

        def f(x, **kwargs):
            if not self.recurse_if(x):
                return f_map(x, **kwargs)
            else:
                next_kwargs = f_kwargs(**kwargs)
                return f_reduce((f(xi, **next_kwargs) for xi in x), **kwargs)

        return f(x, **kwargs)

    def walk(self, x, index=()):
        """
        Iterate over x, yielding (index, value, entering), where

        * ``index``: a tuple of indices up to this point
        * ``value``: equal to ``x[index[0]][...][index[-1]]``. On the first iteration, is
                     ``x`` itself
        * ``entering``: bool. The result of ``recurse_if(value)``
        """
        do_recurse = self.recurse_if(x)
        yield index, x, do_recurse

        if not do_recurse:
            return
        for i, xi in enumerate(x):
            # yield from ...
            yield from self.walk(xi, index + (i,))


# Implementation taken directly from numpy:
# https://github.com/numpy/numpy/blob/d9b1e32cb8ef90d6b4a47853241db2a28146a57d/numpy/core/numeric.py#L1336-L1405
@derived_from(np)
def moveaxis(a, source, destination):
    source = normalize_axis_tuple(source, a.ndim, "source")
    destination = normalize_axis_tuple(destination, a.ndim, "destination")
    if len(source) != len(destination):
        raise ValueError(
            "`source` and `destination` arguments must have "
            "the same number of elements"
        )

    order = [n for n in range(a.ndim) if n not in source]

    for dest, src in sorted(zip(destination, source)):
        order.insert(dest, src)

    result = a.transpose(order)
    return result


# Implementation adapted directly from numpy:
# https://github.com/numpy/numpy/blob/v1.17.0/numpy/core/numeric.py#L1107-L1204
def rollaxis(a, axis, start=0):
    n = a.ndim
    axis = normalize_axis_index(axis, n)
    if start < 0:
        start += n
    msg = "'%s' arg requires %d <= %s < %d, but %d was passed in"
    if not (0 <= start < n + 1):
        raise ValueError(msg % ("start", -n, "start", n + 1, start))
    if axis < start:
        # it's been removed
        start -= 1
    if axis == start:
        return a[...]
    axes = list(range(0, n))
    axes.remove(axis)
    axes.insert(start, axis)
    return a.transpose(axes)
