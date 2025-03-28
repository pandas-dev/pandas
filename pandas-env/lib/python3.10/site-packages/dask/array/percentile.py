from __future__ import annotations

import warnings
from collections.abc import Iterator
from functools import wraps
from numbers import Number

import numpy as np
from tlz import merge

from dask.array.core import Array
from dask.base import tokenize
from dask.highlevelgraph import HighLevelGraph
from dask.utils import derived_from


@wraps(np.percentile)
def _percentile(a, q, method="linear"):
    n = len(a)
    if not len(a):
        return None, n
    if isinstance(q, Iterator):
        q = list(q)
    if a.dtype.name == "category":
        result = np.percentile(a.cat.codes, q, method=method)
        import pandas as pd

        return pd.Categorical.from_codes(result, a.dtype.categories, a.dtype.ordered), n
    if type(a.dtype).__name__ == "DatetimeTZDtype":
        import pandas as pd

        if isinstance(a, (pd.Series, pd.Index)):
            a = a.values

    if np.issubdtype(a.dtype, np.datetime64):
        values = a
        if type(a).__name__ in ("Series", "Index"):
            a2 = values.astype("i8")
        else:
            a2 = values.view("i8")
        result = np.percentile(a2, q, method=method).astype(values.dtype)
        if q[0] == 0:
            # https://github.com/dask/dask/issues/6864
            result[0] = min(result[0], values.min())
        return result, n
    if not np.issubdtype(a.dtype, np.number):
        method = "nearest"
    return np.percentile(a, q, method=method), n


def _tdigest_chunk(a):
    from crick import TDigest

    t = TDigest()
    t.update(a)

    return t


def _percentiles_from_tdigest(qs, digests):
    from crick import TDigest

    t = TDigest()
    t.merge(*digests)

    return np.array(t.quantile(qs / 100.0))


def percentile(a, q, method="linear", internal_method="default", **kwargs):
    """Approximate percentile of 1-D array

    Parameters
    ----------
    a : Array
    q : array_like of float
        Percentile or sequence of percentiles to compute, which must be between
        0 and 100 inclusive.
    method : {'linear', 'lower', 'higher', 'midpoint', 'nearest'}, optional
        The interpolation method to use when the desired percentile lies
        between two data points ``i < j``. Only valid for ``internal_method='dask'``.

        - 'linear': ``i + (j - i) * fraction``, where ``fraction``
          is the fractional part of the index surrounded by ``i``
          and ``j``.
        - 'lower': ``i``.
        - 'higher': ``j``.
        - 'nearest': ``i`` or ``j``, whichever is nearest.
        - 'midpoint': ``(i + j) / 2``.

        .. versionchanged:: 2022.1.0
            This argument was previously called "interpolation"

    internal_method : {'default', 'dask', 'tdigest'}, optional
        What internal method to use. By default will use dask's internal custom
        algorithm (``'dask'``).  If set to ``'tdigest'`` will use tdigest for
        floats and ints and fallback to the ``'dask'`` otherwise.

        .. versionchanged:: 2022.1.0
            This argument was previously called “method”.

    interpolation : str, optional
        Deprecated name for the method keyword argument.

        .. deprecated:: 2022.1.0

    See Also
    --------
    numpy.percentile : Numpy's equivalent Percentile function
    """
    from dask.array.dispatch import percentile_lookup as _percentile
    from dask.array.reductions import quantile
    from dask.array.utils import array_safe, meta_from_array

    if a.ndim > 1:
        q = np.true_divide(q, a.dtype.type(100) if a.dtype.kind == "f" else 100)

        return quantile(a, q, method=method, **kwargs)

    allowed_internal_methods = ["default", "dask", "tdigest"]

    if method in allowed_internal_methods:
        warnings.warn(
            "The `method=` argument was renamed to `internal_method=`",
            FutureWarning,
        )
        internal_method = method

    if "interpolation" in kwargs:
        warnings.warn(
            "The `interpolation=` argument to percentile was renamed to " "`method= ` ",
            FutureWarning,
        )
        method = kwargs.pop("interpolation")

    if kwargs:
        raise TypeError(
            f"percentile() got an unexpected keyword argument {kwargs.keys()}"
        )

    if not a.ndim == 1:
        raise NotImplementedError("Percentiles only implemented for 1-d arrays")
    if isinstance(q, Number):
        q = [q]
    q = array_safe(q, like=meta_from_array(a))
    token = tokenize(a, q, method)

    dtype = a.dtype
    if np.issubdtype(dtype, np.integer):
        dtype = (array_safe([], dtype=dtype, like=meta_from_array(a)) / 0.5).dtype
    meta = meta_from_array(a, dtype=dtype)

    if internal_method not in allowed_internal_methods:
        raise ValueError(
            f"`internal_method=` must be one of {allowed_internal_methods}"
        )

    # Allow using t-digest if method is allowed and dtype is of floating or integer type
    if (
        internal_method == "tdigest"
        and method == "linear"
        and (np.issubdtype(dtype, np.floating) or np.issubdtype(dtype, np.integer))
    ):
        from dask.utils import import_required

        import_required(
            "crick", "crick is a required dependency for using the t-digest method."
        )

        name = "percentile_tdigest_chunk-" + token
        dsk = {
            (name, i): (_tdigest_chunk, key) for i, key in enumerate(a.__dask_keys__())
        }

        name2 = "percentile_tdigest-" + token

        dsk2 = {(name2, 0): (_percentiles_from_tdigest, q, sorted(dsk))}

    # Otherwise use the custom percentile algorithm
    else:
        # Add 0 and 100 during calculation for more robust behavior (hopefully)
        calc_q = np.pad(q, 1, mode="constant")
        calc_q[-1] = 100
        name = "percentile_chunk-" + token
        dsk = {
            (name, i): (_percentile, key, calc_q, method)
            for i, key in enumerate(a.__dask_keys__())
        }

        name2 = "percentile-" + token
        dsk2 = {
            (name2, 0): (
                merge_percentiles,
                q,
                [calc_q] * len(a.chunks[0]),
                sorted(dsk),
                method,
            )
        }

    dsk = merge(dsk, dsk2)
    graph = HighLevelGraph.from_collections(name2, dsk, dependencies=[a])
    return Array(graph, name2, chunks=((len(q),),), meta=meta)


def merge_percentiles(finalq, qs, vals, method="lower", Ns=None, raise_on_nan=True):
    """Combine several percentile calculations of different data.

    Parameters
    ----------

    finalq : numpy.array
        Percentiles to compute (must use same scale as ``qs``).
    qs : sequence of :class:`numpy.array`s
        Percentiles calculated on different sets of data.
    vals : sequence of :class:`numpy.array`s
        Resulting values associated with percentiles ``qs``.
    Ns : sequence of integers
        The number of data elements associated with each data set.
    method : {'linear', 'lower', 'higher', 'midpoint', 'nearest'}
        Specify the interpolation method to use to calculate final
        percentiles.  For more information, see :func:`numpy.percentile`.

    Examples
    --------

    >>> finalq = [10, 20, 30, 40, 50, 60, 70, 80]
    >>> qs = [[20, 40, 60, 80], [20, 40, 60, 80]]
    >>> vals = [np.array([1, 2, 3, 4]), np.array([10, 11, 12, 13])]
    >>> Ns = [100, 100]  # Both original arrays had 100 elements

    >>> merge_percentiles(finalq, qs, vals, Ns=Ns)
    array([ 1,  2,  3,  4, 10, 11, 12, 13])
    """
    from dask.array.utils import array_safe

    if isinstance(finalq, Iterator):
        finalq = list(finalq)
    finalq = array_safe(finalq, like=finalq)
    qs = list(map(list, qs))
    vals = list(vals)
    if Ns is None:
        vals, Ns = zip(*vals)
    Ns = list(Ns)

    L = list(zip(*[(q, val, N) for q, val, N in zip(qs, vals, Ns) if N]))
    if not L:
        if raise_on_nan:
            raise ValueError("No non-trivial arrays found")
        return np.full(len(qs[0]) - 2, np.nan)
    qs, vals, Ns = L

    # TODO: Perform this check above in percentile once dtype checking is easy
    #       Here we silently change meaning
    if vals[0].dtype.name == "category":
        result = merge_percentiles(
            finalq, qs, [v.codes for v in vals], method, Ns, raise_on_nan
        )
        import pandas as pd

        return pd.Categorical.from_codes(result, vals[0].categories, vals[0].ordered)
    if not np.issubdtype(vals[0].dtype, np.number):
        method = "nearest"

    if len(vals) != len(qs) or len(Ns) != len(qs):
        raise ValueError("qs, vals, and Ns parameters must be the same length")

    # transform qs and Ns into number of observations between percentiles
    counts = []
    for q, N in zip(qs, Ns):
        count = np.empty_like(finalq, shape=len(q))
        count[1:] = np.diff(array_safe(q, like=q[0]))
        count[0] = q[0]
        count *= N
        counts.append(count)

    # Sort by calculated percentile values, then number of observations.
    combined_vals = np.concatenate(vals)
    combined_counts = array_safe(np.concatenate(counts), like=combined_vals)
    sort_order = np.argsort(combined_vals)
    combined_vals = np.take(combined_vals, sort_order)
    combined_counts = np.take(combined_counts, sort_order)

    # percentile-like, but scaled by total number of observations
    combined_q = np.cumsum(combined_counts)

    # rescale finalq percentiles to match combined_q
    finalq = array_safe(finalq, like=combined_vals)
    desired_q = finalq * sum(Ns)

    # the behavior of different interpolation methods should be
    # investigated further.
    if method == "linear":
        rv = np.interp(desired_q, combined_q, combined_vals)
    else:
        left = np.searchsorted(combined_q, desired_q, side="left")
        right = np.searchsorted(combined_q, desired_q, side="right") - 1
        np.minimum(left, len(combined_vals) - 1, left)  # don't exceed max index
        lower = np.minimum(left, right)
        upper = np.maximum(left, right)
        if method == "lower":
            rv = combined_vals[lower]
        elif method == "higher":
            rv = combined_vals[upper]
        elif method == "midpoint":
            rv = 0.5 * (combined_vals[lower] + combined_vals[upper])
        elif method == "nearest":
            lower_residual = np.abs(combined_q[lower] - desired_q)
            upper_residual = np.abs(combined_q[upper] - desired_q)
            mask = lower_residual > upper_residual
            index = lower  # alias; we no longer need lower
            index[mask] = upper[mask]
            rv = combined_vals[index]
        else:
            raise ValueError(
                "interpolation method can only be 'linear', 'lower', "
                "'higher', 'midpoint', or 'nearest'"
            )
    return rv


@derived_from(np)
def nanpercentile(a, q, **kwargs):
    from dask.array.reductions import nanquantile

    q = np.true_divide(q, a.dtype.type(100) if a.dtype.kind == "f" else 100)

    return nanquantile(a, q, **kwargs)
