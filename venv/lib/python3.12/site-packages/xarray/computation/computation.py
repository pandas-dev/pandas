"""
Functions for applying functions that act on arrays to xarray's labeled data.

NOTE: This module is currently large and contains various computational functionality.
The long-term plan is to break it down into more focused submodules.
"""

from __future__ import annotations

import functools
from collections import Counter
from collections.abc import (
    Callable,
    Hashable,
)
from typing import TYPE_CHECKING, Any, Literal, cast, overload

import numpy as np

from xarray.compat.array_api_compat import to_like_array
from xarray.core import dtypes, duck_array_ops, utils
from xarray.core.common import zeros_like
from xarray.core.duck_array_ops import datetime_to_numeric
from xarray.core.options import OPTIONS, _get_keep_attrs
from xarray.core.types import Dims, T_DataArray
from xarray.core.utils import (
    is_scalar,
    parse_dims_as_set,
)
from xarray.core.variable import Variable
from xarray.namedarray.parallelcompat import get_chunked_array_type
from xarray.namedarray.pycompat import is_chunked_array
from xarray.structure.alignment import align
from xarray.util.deprecation_helpers import deprecate_dims

if TYPE_CHECKING:
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset

    MissingCoreDimOptions = Literal["raise", "copy", "drop"]

_NO_FILL_VALUE = utils.ReprObject("<no-fill-value>")
_JOINS_WITHOUT_FILL_VALUES = frozenset({"inner", "exact"})


def cov(
    da_a: T_DataArray,
    da_b: T_DataArray,
    dim: Dims = None,
    ddof: int = 1,
    weights: T_DataArray | None = None,
) -> T_DataArray:
    """
    Compute covariance between two DataArray objects along a shared dimension.

    Parameters
    ----------
    da_a : DataArray
        Array to compute.
    da_b : DataArray
        Array to compute.
    dim : str, iterable of hashable, "..." or None, optional
        The dimension along which the covariance will be computed
    ddof : int, default: 1
        If ddof=1, covariance is normalized by N-1, giving an unbiased estimate,
        else normalization is by N.
    weights : DataArray, optional
        Array of weights.

    Returns
    -------
    covariance : DataArray

    See Also
    --------
    pandas.Series.cov : corresponding pandas function
    xarray.corr : respective function to calculate correlation

    Examples
    --------
    >>> from xarray import DataArray
    >>> da_a = DataArray(
    ...     np.array([[1, 2, 3], [0.1, 0.2, 0.3], [3.2, 0.6, 1.8]]),
    ...     dims=("space", "time"),
    ...     coords=[
    ...         ("space", ["IA", "IL", "IN"]),
    ...         ("time", pd.date_range("2000-01-01", freq="1D", periods=3)),
    ...     ],
    ... )
    >>> da_a
    <xarray.DataArray (space: 3, time: 3)> Size: 72B
    array([[1. , 2. , 3. ],
           [0.1, 0.2, 0.3],
           [3.2, 0.6, 1.8]])
    Coordinates:
      * space    (space) <U2 24B 'IA' 'IL' 'IN'
      * time     (time) datetime64[ns] 24B 2000-01-01 2000-01-02 2000-01-03
    >>> da_b = DataArray(
    ...     np.array([[0.2, 0.4, 0.6], [15, 10, 5], [3.2, 0.6, 1.8]]),
    ...     dims=("space", "time"),
    ...     coords=[
    ...         ("space", ["IA", "IL", "IN"]),
    ...         ("time", pd.date_range("2000-01-01", freq="1D", periods=3)),
    ...     ],
    ... )
    >>> da_b
    <xarray.DataArray (space: 3, time: 3)> Size: 72B
    array([[ 0.2,  0.4,  0.6],
           [15. , 10. ,  5. ],
           [ 3.2,  0.6,  1.8]])
    Coordinates:
      * space    (space) <U2 24B 'IA' 'IL' 'IN'
      * time     (time) datetime64[ns] 24B 2000-01-01 2000-01-02 2000-01-03
    >>> xr.cov(da_a, da_b)
    <xarray.DataArray ()> Size: 8B
    array(-3.53055556)
    >>> xr.cov(da_a, da_b, dim="time")
    <xarray.DataArray (space: 3)> Size: 24B
    array([ 0.2       , -0.5       ,  1.69333333])
    Coordinates:
      * space    (space) <U2 24B 'IA' 'IL' 'IN'
    >>> weights = DataArray(
    ...     [4, 2, 1],
    ...     dims=("space"),
    ...     coords=[
    ...         ("space", ["IA", "IL", "IN"]),
    ...     ],
    ... )
    >>> weights
    <xarray.DataArray (space: 3)> Size: 24B
    array([4, 2, 1])
    Coordinates:
      * space    (space) <U2 24B 'IA' 'IL' 'IN'
    >>> xr.cov(da_a, da_b, dim="space", weights=weights)
    <xarray.DataArray (time: 3)> Size: 24B
    array([-4.69346939, -4.49632653, -3.37959184])
    Coordinates:
      * time     (time) datetime64[ns] 24B 2000-01-01 2000-01-02 2000-01-03
    """
    from xarray.core.dataarray import DataArray

    if any(not isinstance(arr, DataArray) for arr in [da_a, da_b]):
        raise TypeError(
            "Only xr.DataArray is supported."
            f"Given {[type(arr) for arr in [da_a, da_b]]}."
        )
    if weights is not None and not isinstance(weights, DataArray):
        raise TypeError(f"Only xr.DataArray is supported. Given {type(weights)}.")
    return _cov_corr(da_a, da_b, weights=weights, dim=dim, ddof=ddof, method="cov")


def corr(
    da_a: T_DataArray,
    da_b: T_DataArray,
    dim: Dims = None,
    weights: T_DataArray | None = None,
) -> T_DataArray:
    """
    Compute the Pearson correlation coefficient between
    two DataArray objects along a shared dimension.

    Parameters
    ----------
    da_a : DataArray
        Array to compute.
    da_b : DataArray
        Array to compute.
    dim : str, iterable of hashable, "..." or None, optional
        The dimension along which the correlation will be computed
    weights : DataArray, optional
        Array of weights.

    Returns
    -------
    correlation: DataArray

    See Also
    --------
    pandas.Series.corr : corresponding pandas function
    xarray.cov : underlying covariance function

    Examples
    --------
    >>> from xarray import DataArray
    >>> da_a = DataArray(
    ...     np.array([[1, 2, 3], [0.1, 0.2, 0.3], [3.2, 0.6, 1.8]]),
    ...     dims=("space", "time"),
    ...     coords=[
    ...         ("space", ["IA", "IL", "IN"]),
    ...         ("time", pd.date_range("2000-01-01", freq="1D", periods=3)),
    ...     ],
    ... )
    >>> da_a
    <xarray.DataArray (space: 3, time: 3)> Size: 72B
    array([[1. , 2. , 3. ],
           [0.1, 0.2, 0.3],
           [3.2, 0.6, 1.8]])
    Coordinates:
      * space    (space) <U2 24B 'IA' 'IL' 'IN'
      * time     (time) datetime64[ns] 24B 2000-01-01 2000-01-02 2000-01-03
    >>> da_b = DataArray(
    ...     np.array([[0.2, 0.4, 0.6], [15, 10, 5], [3.2, 0.6, 1.8]]),
    ...     dims=("space", "time"),
    ...     coords=[
    ...         ("space", ["IA", "IL", "IN"]),
    ...         ("time", pd.date_range("2000-01-01", freq="1D", periods=3)),
    ...     ],
    ... )
    >>> da_b
    <xarray.DataArray (space: 3, time: 3)> Size: 72B
    array([[ 0.2,  0.4,  0.6],
           [15. , 10. ,  5. ],
           [ 3.2,  0.6,  1.8]])
    Coordinates:
      * space    (space) <U2 24B 'IA' 'IL' 'IN'
      * time     (time) datetime64[ns] 24B 2000-01-01 2000-01-02 2000-01-03
    >>> xr.corr(da_a, da_b)
    <xarray.DataArray ()> Size: 8B
    array(-0.57087777)
    >>> xr.corr(da_a, da_b, dim="time")
    <xarray.DataArray (space: 3)> Size: 24B
    array([ 1., -1.,  1.])
    Coordinates:
      * space    (space) <U2 24B 'IA' 'IL' 'IN'
    >>> weights = DataArray(
    ...     [4, 2, 1],
    ...     dims=("space"),
    ...     coords=[
    ...         ("space", ["IA", "IL", "IN"]),
    ...     ],
    ... )
    >>> weights
    <xarray.DataArray (space: 3)> Size: 24B
    array([4, 2, 1])
    Coordinates:
      * space    (space) <U2 24B 'IA' 'IL' 'IN'
    >>> xr.corr(da_a, da_b, dim="space", weights=weights)
    <xarray.DataArray (time: 3)> Size: 24B
    array([-0.50240504, -0.83215028, -0.99057446])
    Coordinates:
      * time     (time) datetime64[ns] 24B 2000-01-01 2000-01-02 2000-01-03
    """
    from xarray.core.dataarray import DataArray

    if any(not isinstance(arr, DataArray) for arr in [da_a, da_b]):
        raise TypeError(
            "Only xr.DataArray is supported."
            f"Given {[type(arr) for arr in [da_a, da_b]]}."
        )
    if weights is not None and not isinstance(weights, DataArray):
        raise TypeError(f"Only xr.DataArray is supported. Given {type(weights)}.")
    return _cov_corr(da_a, da_b, weights=weights, dim=dim, method="corr")


def _cov_corr(
    da_a: T_DataArray,
    da_b: T_DataArray,
    weights: T_DataArray | None = None,
    dim: Dims = None,
    ddof: int = 0,
    method: Literal["cov", "corr"] | None = None,
) -> T_DataArray:
    """
    Internal method for xr.cov() and xr.corr() so only have to
    sanitize the input arrays once and we don't repeat code.
    """
    # 1. Broadcast the two arrays
    da_a, da_b = align(da_a, da_b, join="inner", copy=False)

    # 2. Ignore the nans
    valid_values = da_a.notnull() & da_b.notnull()
    da_a = da_a.where(valid_values)
    da_b = da_b.where(valid_values)

    # 3. Detrend along the given dim
    if weights is not None:
        demeaned_da_a = da_a - da_a.weighted(weights).mean(dim=dim)
        demeaned_da_b = da_b - da_b.weighted(weights).mean(dim=dim)
    else:
        demeaned_da_a = da_a - da_a.mean(dim=dim)
        demeaned_da_b = da_b - da_b.mean(dim=dim)

    # 4. Compute covariance along the given dim
    # N.B. `skipna=True` is required or auto-covariance is computed incorrectly. E.g.
    # Try xr.cov(da,da) for da = xr.DataArray([[1, 2], [1, np.nan]], dims=["x", "time"])
    if weights is not None:
        cov = (
            (demeaned_da_a.conj() * demeaned_da_b)
            .weighted(weights)
            .mean(dim=dim, skipna=True)
        )
    else:
        cov = (demeaned_da_a.conj() * demeaned_da_b).mean(dim=dim, skipna=True)

    if method == "cov":
        # Adjust covariance for degrees of freedom
        valid_count = valid_values.sum(dim)
        adjust = valid_count / (valid_count - ddof)
        # I think the cast is required because of `T_DataArray` + `T_Xarray` (would be
        # the same with `T_DatasetOrArray`)
        # https://github.com/pydata/xarray/pull/8384#issuecomment-1784228026
        return cast(T_DataArray, cov * adjust)

    else:
        # Compute std and corr
        if weights is not None:
            da_a_std = da_a.weighted(weights).std(dim=dim)
            da_b_std = da_b.weighted(weights).std(dim=dim)
        else:
            da_a_std = da_a.std(dim=dim)
            da_b_std = da_b.std(dim=dim)
        corr = cov / (da_a_std * da_b_std)
        return cast(T_DataArray, corr)


def cross(
    a: DataArray | Variable, b: DataArray | Variable, *, dim: Hashable
) -> DataArray | Variable:
    """
    Compute the cross product of two (arrays of) vectors.

    The cross product of `a` and `b` in :math:`R^3` is a vector
    perpendicular to both `a` and `b`. The vectors in `a` and `b` are
    defined by the values along the dimension `dim` and can have sizes
    1, 2 or 3. Where the size of either `a` or `b` is
    1 or 2, the remaining components of the input vector is assumed to
    be zero and the cross product calculated accordingly. In cases where
    both input vectors have dimension 2, the z-component of the cross
    product is returned.

    Parameters
    ----------
    a, b : DataArray or Variable
        Components of the first and second vector(s).
    dim : hashable
        The dimension along which the cross product will be computed.
        Must be available in both vectors.

    Examples
    --------
    Vector cross-product with 3 dimensions:

    >>> a = xr.DataArray([1, 2, 3])
    >>> b = xr.DataArray([4, 5, 6])
    >>> xr.cross(a, b, dim="dim_0")
    <xarray.DataArray (dim_0: 3)> Size: 24B
    array([-3,  6, -3])
    Dimensions without coordinates: dim_0

    Vector cross-product with 3 dimensions but zeros at the last axis
    yields the same results as with 2 dimensions:

    >>> a = xr.DataArray([1, 2, 0])
    >>> b = xr.DataArray([4, 5, 0])
    >>> xr.cross(a, b, dim="dim_0")
    <xarray.DataArray (dim_0: 3)> Size: 24B
    array([ 0,  0, -3])
    Dimensions without coordinates: dim_0

    Multiple vector cross-products. Note that the direction of the
    cross product vector is defined by the right-hand rule:

    >>> a = xr.DataArray(
    ...     [[1, 2, 3], [4, 5, 6]],
    ...     dims=("time", "cartesian"),
    ...     coords=dict(
    ...         time=(["time"], [0, 1]),
    ...         cartesian=(["cartesian"], ["x", "y", "z"]),
    ...     ),
    ... )
    >>> b = xr.DataArray(
    ...     [[4, 5, 6], [1, 2, 3]],
    ...     dims=("time", "cartesian"),
    ...     coords=dict(
    ...         time=(["time"], [0, 1]),
    ...         cartesian=(["cartesian"], ["x", "y", "z"]),
    ...     ),
    ... )
    >>> xr.cross(a, b, dim="cartesian")
    <xarray.DataArray (time: 2, cartesian: 3)> Size: 48B
    array([[-3,  6, -3],
           [ 3, -6,  3]])
    Coordinates:
      * time       (time) int64 16B 0 1
      * cartesian  (cartesian) <U1 12B 'x' 'y' 'z'

    Cross can be called on Datasets by converting to DataArrays and later
    back to a Dataset:

    >>> ds_a = xr.Dataset(dict(x=("dim_0", [1]), y=("dim_0", [2]), z=("dim_0", [3])))
    >>> ds_b = xr.Dataset(dict(x=("dim_0", [4]), y=("dim_0", [5]), z=("dim_0", [6])))
    >>> c = xr.cross(
    ...     ds_a.to_dataarray("cartesian"),
    ...     ds_b.to_dataarray("cartesian"),
    ...     dim="cartesian",
    ... )
    >>> c.to_dataset(dim="cartesian")
    <xarray.Dataset> Size: 24B
    Dimensions:  (dim_0: 1)
    Dimensions without coordinates: dim_0
    Data variables:
        x        (dim_0) int64 8B -3
        y        (dim_0) int64 8B 6
        z        (dim_0) int64 8B -3

    See Also
    --------
    numpy.cross : Corresponding numpy function
    """

    if dim not in a.dims:
        raise ValueError(f"Dimension {dim!r} not on a")
    elif dim not in b.dims:
        raise ValueError(f"Dimension {dim!r} not on b")

    if not 1 <= a.sizes[dim] <= 3:
        raise ValueError(
            f"The size of {dim!r} on a must be 1, 2, or 3 to be "
            f"compatible with a cross product but is {a.sizes[dim]}"
        )
    elif not 1 <= b.sizes[dim] <= 3:
        raise ValueError(
            f"The size of {dim!r} on b must be 1, 2, or 3 to be "
            f"compatible with a cross product but is {b.sizes[dim]}"
        )

    all_dims = list(dict.fromkeys(a.dims + b.dims))

    if a.sizes[dim] != b.sizes[dim]:
        # Arrays have different sizes. Append zeros where the smaller
        # array is missing a value, zeros will not affect np.cross:

        if (
            not isinstance(a, Variable)  # Only used to make mypy happy.
            and dim in getattr(a, "coords", {})
            and not isinstance(b, Variable)  # Only used to make mypy happy.
            and dim in getattr(b, "coords", {})
        ):
            # If the arrays have coords we know which indexes to fill
            # with zeros:
            a, b = align(
                a,
                b,
                fill_value=0,
                join="outer",
                exclude=set(all_dims) - {dim},
            )
        elif min(a.sizes[dim], b.sizes[dim]) == 2:
            # If the array doesn't have coords we can only infer
            # that it has composite values if the size is at least 2.
            # Once padded, rechunk the padded array because apply_ufunc
            # requires core dimensions not to be chunked:
            if a.sizes[dim] < b.sizes[dim]:
                a = a.pad({dim: (0, 1)}, constant_values=0)
                # TODO: Should pad or apply_ufunc handle correct chunking?
                a = a.chunk({dim: -1}) if is_chunked_array(a.data) else a
            else:
                b = b.pad({dim: (0, 1)}, constant_values=0)
                # TODO: Should pad or apply_ufunc handle correct chunking?
                b = b.chunk({dim: -1}) if is_chunked_array(b.data) else b
        else:
            raise ValueError(
                f"{dim!r} on {'a' if a.sizes[dim] == 1 else 'b'} is incompatible:"
                " dimensions without coordinates must have have a length of 2 or 3"
            )

    from xarray.computation.apply_ufunc import apply_ufunc

    c = apply_ufunc(
        duck_array_ops.cross,
        a,
        b,
        input_core_dims=[[dim], [dim]],
        output_core_dims=[[dim] if a.sizes[dim] == 3 else []],
        dask="parallelized",
        output_dtypes=[np.result_type(a, b)],
    )
    c = c.transpose(*all_dims, missing_dims="ignore")

    return c


@deprecate_dims
def dot(
    *arrays,
    dim: Dims = None,
    **kwargs: Any,
):
    """Generalized dot product for xarray objects. Like ``np.einsum``, but
    provides a simpler interface based on array dimension names.

    Parameters
    ----------
    *arrays : DataArray or Variable
        Arrays to compute.
    dim : str, iterable of hashable, "..." or None, optional
        Which dimensions to sum over. Ellipsis ('...') sums over all dimensions.
        If not specified, then all the common dimensions are summed over.
    **kwargs : dict
        Additional keyword arguments passed to ``numpy.einsum`` or
        ``dask.array.einsum``

    Returns
    -------
    DataArray

    See Also
    --------
    numpy.einsum
    dask.array.einsum
    opt_einsum.contract

    Notes
    -----
    We recommend installing the optional ``opt_einsum`` package, or alternatively passing ``optimize=True``,
    which is passed through to ``np.einsum``, and works for most array backends.

    Examples
    --------
    >>> da_a = xr.DataArray(np.arange(3 * 2).reshape(3, 2), dims=["a", "b"])
    >>> da_b = xr.DataArray(np.arange(3 * 2 * 2).reshape(3, 2, 2), dims=["a", "b", "c"])
    >>> da_c = xr.DataArray(np.arange(2 * 3).reshape(2, 3), dims=["c", "d"])

    >>> da_a
    <xarray.DataArray (a: 3, b: 2)> Size: 48B
    array([[0, 1],
           [2, 3],
           [4, 5]])
    Dimensions without coordinates: a, b

    >>> da_b
    <xarray.DataArray (a: 3, b: 2, c: 2)> Size: 96B
    array([[[ 0,  1],
            [ 2,  3]],
    <BLANKLINE>
           [[ 4,  5],
            [ 6,  7]],
    <BLANKLINE>
           [[ 8,  9],
            [10, 11]]])
    Dimensions without coordinates: a, b, c

    >>> da_c
    <xarray.DataArray (c: 2, d: 3)> Size: 48B
    array([[0, 1, 2],
           [3, 4, 5]])
    Dimensions without coordinates: c, d

    >>> xr.dot(da_a, da_b, dim=["a", "b"])
    <xarray.DataArray (c: 2)> Size: 16B
    array([110, 125])
    Dimensions without coordinates: c

    >>> xr.dot(da_a, da_b, dim=["a"])
    <xarray.DataArray (b: 2, c: 2)> Size: 32B
    array([[40, 46],
           [70, 79]])
    Dimensions without coordinates: b, c

    >>> xr.dot(da_a, da_b, da_c, dim=["b", "c"])
    <xarray.DataArray (a: 3, d: 3)> Size: 72B
    array([[  9,  14,  19],
           [ 93, 150, 207],
           [273, 446, 619]])
    Dimensions without coordinates: a, d

    >>> xr.dot(da_a, da_b)
    <xarray.DataArray (c: 2)> Size: 16B
    array([110, 125])
    Dimensions without coordinates: c

    >>> xr.dot(da_a, da_b, dim=...)
    <xarray.DataArray ()> Size: 8B
    array(235)
    """
    from xarray.core.dataarray import DataArray

    if any(not isinstance(arr, Variable | DataArray) for arr in arrays):
        raise TypeError(
            "Only xr.DataArray and xr.Variable are supported."
            f"Given {[type(arr) for arr in arrays]}."
        )

    if len(arrays) == 0:
        raise TypeError("At least one array should be given.")

    common_dims: set[Hashable] = set.intersection(*(set(arr.dims) for arr in arrays))
    all_dims = []
    for arr in arrays:
        all_dims += [d for d in arr.dims if d not in all_dims]

    einsum_axes = "abcdefghijklmnopqrstuvwxyz"
    dim_map = {d: einsum_axes[i] for i, d in enumerate(all_dims)}

    dot_dims: set[Hashable]
    if dim is None:
        # find dimensions that occur more than once
        dim_counts: Counter = Counter()
        for arr in arrays:
            dim_counts.update(arr.dims)
        dot_dims = {d for d, c in dim_counts.items() if c > 1}
    else:
        dot_dims = parse_dims_as_set(dim, all_dims=set(all_dims))

    # dimensions to be parallelized
    broadcast_dims = common_dims - dot_dims
    input_core_dims = [
        [d for d in arr.dims if d not in broadcast_dims] for arr in arrays
    ]
    output_core_dims = [
        [d for d in all_dims if d not in dot_dims and d not in broadcast_dims]
    ]

    # construct einsum subscripts, such as '...abc,...ab->...c'
    # Note: input_core_dims are always moved to the last position
    subscripts_list = [
        "..." + "".join(dim_map[d] for d in ds) for ds in input_core_dims
    ]
    subscripts = ",".join(subscripts_list)
    subscripts += "->..." + "".join(dim_map[d] for d in output_core_dims[0])

    join = OPTIONS["arithmetic_join"]
    # using "inner" emulates `(a * b).sum()` for all joins (except "exact")
    if join != "exact":
        join = "inner"

    # subscripts should be passed to np.einsum as arg, not as kwargs. We need
    # to construct a partial function for apply_ufunc to work.
    func = functools.partial(duck_array_ops.einsum, subscripts, **kwargs)
    from xarray.computation.apply_ufunc import apply_ufunc

    result = apply_ufunc(
        func,
        *arrays,
        input_core_dims=input_core_dims,
        output_core_dims=output_core_dims,
        join=join,
        dask="allowed",
    )
    return result.transpose(*all_dims, missing_dims="ignore")


def where(cond, x, y, keep_attrs=None):
    """Return elements from `x` or `y` depending on `cond`.

    Performs xarray-like broadcasting across input arguments.

    All dimension coordinates on `x` and `y`  must be aligned with each
    other and with `cond`.

    Parameters
    ----------
    cond : scalar, array, Variable, DataArray or Dataset
        When True, return values from `x`, otherwise returns values from `y`.
    x : scalar, array, Variable, DataArray or Dataset
        values to choose from where `cond` is True
    y : scalar, array, Variable, DataArray or Dataset
        values to choose from where `cond` is False
    keep_attrs : bool or str or callable, optional
        How to treat attrs. If True, keep the attrs of `x`.

    Returns
    -------
    Dataset, DataArray, Variable or array
        In priority order: Dataset, DataArray, Variable or array, whichever
        type appears as an input argument.

    Examples
    --------
    >>> x = xr.DataArray(
    ...     0.1 * np.arange(10),
    ...     dims=["lat"],
    ...     coords={"lat": np.arange(10)},
    ...     name="sst",
    ... )
    >>> x
    <xarray.DataArray 'sst' (lat: 10)> Size: 80B
    array([0. , 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    Coordinates:
      * lat      (lat) int64 80B 0 1 2 3 4 5 6 7 8 9

    >>> xr.where(x < 0.5, x, x * 100)
    <xarray.DataArray 'sst' (lat: 10)> Size: 80B
    array([ 0. ,  0.1,  0.2,  0.3,  0.4, 50. , 60. , 70. , 80. , 90. ])
    Coordinates:
      * lat      (lat) int64 80B 0 1 2 3 4 5 6 7 8 9

    >>> y = xr.DataArray(
    ...     0.1 * np.arange(9).reshape(3, 3),
    ...     dims=["lat", "lon"],
    ...     coords={"lat": np.arange(3), "lon": 10 + np.arange(3)},
    ...     name="sst",
    ... )
    >>> y
    <xarray.DataArray 'sst' (lat: 3, lon: 3)> Size: 72B
    array([[0. , 0.1, 0.2],
           [0.3, 0.4, 0.5],
           [0.6, 0.7, 0.8]])
    Coordinates:
      * lat      (lat) int64 24B 0 1 2
      * lon      (lon) int64 24B 10 11 12

    >>> xr.where(y.lat < 1, y, -1)
    <xarray.DataArray (lat: 3, lon: 3)> Size: 72B
    array([[ 0. ,  0.1,  0.2],
           [-1. , -1. , -1. ],
           [-1. , -1. , -1. ]])
    Coordinates:
      * lat      (lat) int64 24B 0 1 2
      * lon      (lon) int64 24B 10 11 12

    >>> cond = xr.DataArray([True, False], dims=["x"])
    >>> x = xr.DataArray([1, 2], dims=["y"])
    >>> xr.where(cond, x, 0)
    <xarray.DataArray (x: 2, y: 2)> Size: 32B
    array([[1, 2],
           [0, 0]])
    Dimensions without coordinates: x, y

    See Also
    --------
    numpy.where : corresponding numpy function
    Dataset.where, DataArray.where :
        equivalent methods
    """
    from xarray.core.dataset import Dataset

    if keep_attrs is None:
        keep_attrs = _get_keep_attrs(default=False)

    # alignment for three arguments is complicated, so don't support it yet
    from xarray.computation.apply_ufunc import apply_ufunc

    result = apply_ufunc(
        duck_array_ops.where,
        cond,
        x,
        y,
        join="exact",
        dataset_join="exact",
        dask="allowed",
        keep_attrs=keep_attrs,
    )

    # keep the attributes of x, the second parameter, by default to
    # be consistent with the `where` method of `DataArray` and `Dataset`
    # rebuild the attrs from x at each level of the output, which could be
    # Dataset, DataArray, or Variable, and also handle coords
    if keep_attrs is True and hasattr(result, "attrs"):
        if isinstance(y, Dataset) and not isinstance(x, Dataset):
            # handle special case where x gets promoted to Dataset
            result.attrs = {}
            if getattr(x, "name", None) in result.data_vars:
                result[x.name].attrs = getattr(x, "attrs", {})
        else:
            # otherwise, fill in global attrs and variable attrs (if they exist)
            result.attrs = getattr(x, "attrs", {})
            for v in getattr(result, "data_vars", []):
                result[v].attrs = getattr(getattr(x, v, None), "attrs", {})
        for c in getattr(result, "coords", []):
            # always fill coord attrs of x
            result[c].attrs = getattr(getattr(x, c, None), "attrs", {})

    return result


@overload
def polyval(
    coord: DataArray, coeffs: DataArray, degree_dim: Hashable = "degree"
) -> DataArray: ...


@overload
def polyval(
    coord: DataArray, coeffs: Dataset, degree_dim: Hashable = "degree"
) -> Dataset: ...


@overload
def polyval(
    coord: Dataset, coeffs: DataArray, degree_dim: Hashable = "degree"
) -> Dataset: ...


@overload
def polyval(
    coord: Dataset, coeffs: Dataset, degree_dim: Hashable = "degree"
) -> Dataset: ...


@overload
def polyval(
    coord: Dataset | DataArray,
    coeffs: Dataset | DataArray,
    degree_dim: Hashable = "degree",
) -> Dataset | DataArray: ...


def polyval(
    coord: Dataset | DataArray,
    coeffs: Dataset | DataArray,
    degree_dim: Hashable = "degree",
) -> Dataset | DataArray:
    """Evaluate a polynomial at specific values

    Parameters
    ----------
    coord : DataArray or Dataset
        Values at which to evaluate the polynomial.
    coeffs : DataArray or Dataset
        Coefficients of the polynomial.
    degree_dim : Hashable, default: "degree"
        Name of the polynomial degree dimension in `coeffs`.

    Returns
    -------
    DataArray or Dataset
        Evaluated polynomial.

    See Also
    --------
    xarray.DataArray.polyfit
    numpy.polynomial.polynomial.polyval
    """

    if degree_dim not in coeffs._indexes:
        raise ValueError(
            f"Dimension `{degree_dim}` should be a coordinate variable with labels."
        )
    if not np.issubdtype(coeffs[degree_dim].dtype, np.integer):
        raise ValueError(
            f"Dimension `{degree_dim}` should be of integer dtype. Received {coeffs[degree_dim].dtype} instead."
        )
    max_deg = coeffs[degree_dim].max().item()
    coeffs = coeffs.reindex(
        {degree_dim: np.arange(max_deg + 1)}, fill_value=0, copy=False
    )
    coord = _ensure_numeric(coord)

    # using Horner's method
    # https://en.wikipedia.org/wiki/Horner%27s_method
    res = zeros_like(coord) + coeffs.isel({degree_dim: max_deg}, drop=True)
    for deg in range(max_deg - 1, -1, -1):
        res *= coord
        res += coeffs.isel({degree_dim: deg}, drop=True)

    return res


def _ensure_numeric(data: Dataset | DataArray) -> Dataset | DataArray:
    """Converts all datetime64 variables to float64

    Parameters
    ----------
    data : DataArray or Dataset
        Variables with possible datetime dtypes.

    Returns
    -------
    DataArray or Dataset
        Variables with datetime64 dtypes converted to float64.
    """
    from xarray.core.dataset import Dataset

    def _cfoffset(x: DataArray) -> Any:
        scalar = x.compute().data[0]
        if not is_scalar(scalar):
            # we do not get a scalar back on dask == 2021.04.1
            scalar = scalar.item()
        return type(scalar)(1970, 1, 1)

    def to_floatable(x: DataArray) -> DataArray:
        if x.dtype.kind in "MO":
            # datetimes (CFIndexes are object type)
            offset = (
                np.datetime64("1970-01-01") if x.dtype.kind == "M" else _cfoffset(x)
            )
            return x.copy(
                data=datetime_to_numeric(x.data, offset=offset, datetime_unit="ns"),
            )
        elif x.dtype.kind == "m":
            # timedeltas
            return duck_array_ops.astype(x, dtype=float)
        return x

    if isinstance(data, Dataset):
        return data.map(to_floatable)
    else:
        return to_floatable(data)


def _calc_idxminmax(
    *,
    array,
    func: Callable,
    dim: Hashable | None = None,
    skipna: bool | None = None,
    fill_value: Any = dtypes.NA,
    keep_attrs: bool | None = None,
):
    """Apply common operations for idxmin and idxmax."""
    # This function doesn't make sense for scalars so don't try
    if not array.ndim:
        raise ValueError("This function does not apply for scalars")

    if dim is not None:
        pass  # Use the dim if available
    elif array.ndim == 1:
        # it is okay to guess the dim if there is only 1
        dim = array.dims[0]
    else:
        # The dim is not specified and ambiguous.  Don't guess.
        raise ValueError("Must supply 'dim' argument for multidimensional arrays")

    if dim not in array.dims:
        raise KeyError(
            f"Dimension {dim!r} not found in array dimensions {array.dims!r}"
        )
    if dim not in array.coords:
        raise KeyError(
            f"Dimension {dim!r} is not one of the coordinates {tuple(array.coords.keys())}"
        )

    # These are dtypes with NaN values argmin and argmax can handle
    na_dtypes = "cfO"

    if skipna or (skipna is None and array.dtype.kind in na_dtypes):
        # Need to skip NaN values since argmin and argmax can't handle them
        allna = array.isnull().all(dim)
        array = array.where(~allna, 0)

    # This will run argmin or argmax.
    indx = func(array, dim=dim, axis=None, keep_attrs=keep_attrs, skipna=skipna)

    # Handle chunked arrays (e.g. dask).
    coord = array[dim]._variable.to_base_variable()
    if is_chunked_array(array.data):
        chunkmanager = get_chunked_array_type(array.data)
        coord_array = chunkmanager.from_array(
            array[dim].data, chunks=((array.sizes[dim],),)
        )
        coord = coord.copy(data=coord_array)
    else:
        coord = coord.copy(data=to_like_array(array[dim].data, array.data))

    res = indx._replace(coord[(indx.variable,)]).rename(dim)

    if skipna or (skipna is None and array.dtype.kind in na_dtypes):
        # Put the NaN values back in after removing them
        res = res.where(~allna, fill_value)

    # Copy attributes from argmin/argmax, if any
    res.attrs = indx.attrs

    return res
