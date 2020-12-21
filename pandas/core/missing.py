"""
Routines for filling missing data.
"""
from functools import partial
from typing import TYPE_CHECKING, Any, List, Optional, Set, Union

import numpy as np

from pandas._libs import algos, lib
from pandas._typing import ArrayLike, Axis, DtypeObj
from pandas.compat._optional import import_optional_dependency

from pandas.core.dtypes.cast import infer_dtype_from_array
from pandas.core.dtypes.common import (
    ensure_float64,
    is_integer_dtype,
    is_numeric_v_string_like,
    needs_i8_conversion,
)
from pandas.core.dtypes.missing import isna

if TYPE_CHECKING:
    from pandas import Index


def mask_missing(arr: ArrayLike, values_to_mask) -> np.ndarray:
    """
    Return a masking array of same size/shape as arr
    with entries equaling any member of values_to_mask set to True

    Parameters
    ----------
    arr : ArrayLike
    values_to_mask: list, tuple, or scalar

    Returns
    -------
    np.ndarray[bool]
    """
    # When called from Block.replace/replace_list, values_to_mask is a scalar
    #  known to be holdable by arr.
    # When called from Series._single_replace, values_to_mask is tuple or list
    dtype, values_to_mask = infer_dtype_from_array(values_to_mask)
    values_to_mask = np.array(values_to_mask, dtype=dtype)

    na_mask = isna(values_to_mask)
    nonna = values_to_mask[~na_mask]

    # GH 21977
    mask = np.zeros(arr.shape, dtype=bool)
    for x in nonna:
        if is_numeric_v_string_like(arr, x):
            # GH#29553 prevent numpy deprecation warnings
            pass
        else:
            mask |= arr == x

    if na_mask.any():
        mask |= isna(arr)

    return mask


def clean_fill_method(method, allow_nearest: bool = False):
    # asfreq is compat for resampling
    if method in [None, "asfreq"]:
        return None

    if isinstance(method, str):
        method = method.lower()
        if method == "ffill":
            method = "pad"
        elif method == "bfill":
            method = "backfill"

    valid_methods = ["pad", "backfill"]
    expecting = "pad (ffill) or backfill (bfill)"
    if allow_nearest:
        valid_methods.append("nearest")
        expecting = "pad (ffill), backfill (bfill) or nearest"
    if method not in valid_methods:
        raise ValueError(f"Invalid fill method. Expecting {expecting}. Got {method}")
    return method


# interpolation methods that dispatch to np.interp

NP_METHODS = ["linear", "time", "index", "values"]

# interpolation methods that dispatch to _interpolate_scipy_wrapper

SP_METHODS = [
    "nearest",
    "zero",
    "slinear",
    "quadratic",
    "cubic",
    "barycentric",
    "krogh",
    "spline",
    "polynomial",
    "from_derivatives",
    "piecewise_polynomial",
    "pchip",
    "akima",
    "cubicspline",
]


def clean_interp_method(method: str, **kwargs) -> str:
    order = kwargs.get("order")

    if method in ("spline", "polynomial") and order is None:
        raise ValueError("You must specify the order of the spline or polynomial.")

    valid = NP_METHODS + SP_METHODS
    if method not in valid:
        raise ValueError(f"method must be one of {valid}. Got '{method}' instead.")

    return method


def find_valid_index(values, how: str):
    """
    Retrieves the index of the first valid value.

    Parameters
    ----------
    values : ndarray or ExtensionArray
    how : {'first', 'last'}
        Use this parameter to change between the first or last valid index.

    Returns
    -------
    int or None
    """
    assert how in ["first", "last"]

    if len(values) == 0:  # early stop
        return None

    is_valid = ~isna(values)

    if values.ndim == 2:
        is_valid = is_valid.any(1)  # reduce axis 1

    if how == "first":
        idxpos = is_valid[::].argmax()

    if how == "last":
        idxpos = len(values) - 1 - is_valid[::-1].argmax()

    chk_notna = is_valid[idxpos]

    if not chk_notna:
        return None
    return idxpos


def interpolate_1d(
    xvalues: "Index",
    yvalues: np.ndarray,
    method: Optional[str] = "linear",
    limit: Optional[int] = None,
    limit_direction: str = "forward",
    limit_area: Optional[str] = None,
    fill_value: Optional[Any] = None,
    bounds_error: bool = False,
    order: Optional[int] = None,
    **kwargs,
):
    """
    Logic for the 1-d interpolation.  The result should be 1-d, inputs
    xvalues and yvalues will each be 1-d arrays of the same length.

    Bounds_error is currently hardcoded to False since non-scipy ones don't
    take it as an argument.
    """
    invalid = isna(yvalues)
    valid = ~invalid

    if not valid.any():
        result = np.empty(xvalues.shape, dtype=np.float64)
        result.fill(np.nan)
        return result

    if valid.all():
        return yvalues

    if method == "time":
        if not needs_i8_conversion(xvalues.dtype):
            raise ValueError(
                "time-weighted interpolation only works "
                "on Series or DataFrames with a "
                "DatetimeIndex"
            )
        method = "values"

    valid_limit_directions = ["forward", "backward", "both"]
    limit_direction = limit_direction.lower()
    if limit_direction not in valid_limit_directions:
        raise ValueError(
            "Invalid limit_direction: expecting one of "
            f"{valid_limit_directions}, got '{limit_direction}'."
        )

    if limit_area is not None:
        valid_limit_areas = ["inside", "outside"]
        limit_area = limit_area.lower()
        if limit_area not in valid_limit_areas:
            raise ValueError(
                f"Invalid limit_area: expecting one of {valid_limit_areas}, got "
                f"{limit_area}."
            )

    # default limit is unlimited GH #16282
    limit = algos.validate_limit(nobs=None, limit=limit)

    # These are sets of index pointers to invalid values... i.e. {0, 1, etc...
    all_nans = set(np.flatnonzero(invalid))
    start_nans = set(range(find_valid_index(yvalues, "first")))
    end_nans = set(range(1 + find_valid_index(yvalues, "last"), len(valid)))
    mid_nans = all_nans - start_nans - end_nans

    # Like the sets above, preserve_nans contains indices of invalid values,
    # but in this case, it is the final set of indices that need to be
    # preserved as NaN after the interpolation.

    # For example if limit_direction='forward' then preserve_nans will
    # contain indices of NaNs at the beginning of the series, and NaNs that
    # are more than'limit' away from the prior non-NaN.

    # set preserve_nans based on direction using _interp_limit
    preserve_nans: Union[List, Set]
    if limit_direction == "forward":
        preserve_nans = start_nans | set(_interp_limit(invalid, limit, 0))
    elif limit_direction == "backward":
        preserve_nans = end_nans | set(_interp_limit(invalid, 0, limit))
    else:
        # both directions... just use _interp_limit
        preserve_nans = set(_interp_limit(invalid, limit, limit))

    # if limit_area is set, add either mid or outside indices
    # to preserve_nans GH #16284
    if limit_area == "inside":
        # preserve NaNs on the outside
        preserve_nans |= start_nans | end_nans
    elif limit_area == "outside":
        # preserve NaNs on the inside
        preserve_nans |= mid_nans

    # sort preserve_nans and covert to list
    preserve_nans = sorted(preserve_nans)

    result = yvalues.copy()

    # xarr to pass to NumPy/SciPy
    xarr = xvalues._values
    if needs_i8_conversion(xarr.dtype):
        # GH#1646 for dt64tz
        xarr = xarr.view("i8")

    if method == "linear":
        inds = xarr
    else:
        inds = np.asarray(xarr)

        if method in ("values", "index"):
            if inds.dtype == np.object_:
                inds = lib.maybe_convert_objects(inds)

    if method in NP_METHODS:
        # np.interp requires sorted X values, #21037
        indexer = np.argsort(inds[valid])
        result[invalid] = np.interp(
            inds[invalid], inds[valid][indexer], yvalues[valid][indexer]
        )
    else:
        result[invalid] = _interpolate_scipy_wrapper(
            inds[valid],
            yvalues[valid],
            inds[invalid],
            method=method,
            fill_value=fill_value,
            bounds_error=bounds_error,
            order=order,
            **kwargs,
        )

    result[preserve_nans] = np.nan
    return result


def _interpolate_scipy_wrapper(
    x, y, new_x, method, fill_value=None, bounds_error=False, order=None, **kwargs
):
    """
    Passed off to scipy.interpolate.interp1d. method is scipy's kind.
    Returns an array interpolated at new_x.  Add any new methods to
    the list in _clean_interp_method.
    """
    extra = f"{method} interpolation requires SciPy."
    import_optional_dependency("scipy", extra=extra)
    from scipy import interpolate

    new_x = np.asarray(new_x)

    # ignores some kwargs that could be passed along.
    alt_methods = {
        "barycentric": interpolate.barycentric_interpolate,
        "krogh": interpolate.krogh_interpolate,
        "from_derivatives": _from_derivatives,
        "piecewise_polynomial": _from_derivatives,
    }

    if getattr(x, "_is_all_dates", False):
        # GH 5975, scipy.interp1d can't handle datetime64s
        x, new_x = x._values.astype("i8"), new_x.astype("i8")

    if method == "pchip":
        alt_methods["pchip"] = interpolate.pchip_interpolate
    elif method == "akima":
        alt_methods["akima"] = _akima_interpolate
    elif method == "cubicspline":
        alt_methods["cubicspline"] = _cubicspline_interpolate

    interp1d_methods = [
        "nearest",
        "zero",
        "slinear",
        "quadratic",
        "cubic",
        "polynomial",
    ]
    if method in interp1d_methods:
        if method == "polynomial":
            method = order
        terp = interpolate.interp1d(
            x, y, kind=method, fill_value=fill_value, bounds_error=bounds_error
        )
        new_y = terp(new_x)
    elif method == "spline":
        # GH #10633, #24014
        if isna(order) or (order <= 0):
            raise ValueError(
                f"order needs to be specified and greater than 0; got order: {order}"
            )
        terp = interpolate.UnivariateSpline(x, y, k=order, **kwargs)
        new_y = terp(new_x)
    else:
        # GH 7295: need to be able to write for some reason
        # in some circumstances: check all three
        if not x.flags.writeable:
            x = x.copy()
        if not y.flags.writeable:
            y = y.copy()
        if not new_x.flags.writeable:
            new_x = new_x.copy()
        method = alt_methods[method]
        new_y = method(x, y, new_x, **kwargs)
    return new_y


def _from_derivatives(xi, yi, x, order=None, der=0, extrapolate=False):
    """
    Convenience function for interpolate.BPoly.from_derivatives.

    Construct a piecewise polynomial in the Bernstein basis, compatible
    with the specified values and derivatives at breakpoints.

    Parameters
    ----------
    xi : array_like
        sorted 1D array of x-coordinates
    yi : array_like or list of array-likes
        yi[i][j] is the j-th derivative known at xi[i]
    order: None or int or array_like of ints. Default: None.
        Specifies the degree of local polynomials. If not None, some
        derivatives are ignored.
    der : int or list
        How many derivatives to extract; None for all potentially nonzero
        derivatives (that is a number equal to the number of points), or a
        list of derivatives to extract. This number includes the function
        value as 0th derivative.
     extrapolate : bool, optional
        Whether to extrapolate to ouf-of-bounds points based on first and last
        intervals, or to return NaNs. Default: True.

    See Also
    --------
    scipy.interpolate.BPoly.from_derivatives

    Returns
    -------
    y : scalar or array_like
        The result, of length R or length M or M by R.
    """
    from scipy import interpolate

    # return the method for compat with scipy version & backwards compat
    method = interpolate.BPoly.from_derivatives
    m = method(xi, yi.reshape(-1, 1), orders=order, extrapolate=extrapolate)

    return m(x)


def _akima_interpolate(xi, yi, x, der=0, axis=0):
    """
    Convenience function for akima interpolation.
    xi and yi are arrays of values used to approximate some function f,
    with ``yi = f(xi)``.

    See `Akima1DInterpolator` for details.

    Parameters
    ----------
    xi : array_like
        A sorted list of x-coordinates, of length N.
    yi : array_like
        A 1-D array of real values.  `yi`'s length along the interpolation
        axis must be equal to the length of `xi`. If N-D array, use axis
        parameter to select correct axis.
    x : scalar or array_like
        Of length M.
    der : int, optional
        How many derivatives to extract; None for all potentially
        nonzero derivatives (that is a number equal to the number
        of points), or a list of derivatives to extract. This number
        includes the function value as 0th derivative.
    axis : int, optional
        Axis in the yi array corresponding to the x-coordinate values.

    See Also
    --------
    scipy.interpolate.Akima1DInterpolator

    Returns
    -------
    y : scalar or array_like
        The result, of length R or length M or M by R,

    """
    from scipy import interpolate

    P = interpolate.Akima1DInterpolator(xi, yi, axis=axis)

    return P(x, nu=der)


def _cubicspline_interpolate(xi, yi, x, axis=0, bc_type="not-a-knot", extrapolate=None):
    """
    Convenience function for cubic spline data interpolator.

    See `scipy.interpolate.CubicSpline` for details.

    Parameters
    ----------
    xi : array_like, shape (n,)
        1-d array containing values of the independent variable.
        Values must be real, finite and in strictly increasing order.
    yi : array_like
        Array containing values of the dependent variable. It can have
        arbitrary number of dimensions, but the length along ``axis``
        (see below) must match the length of ``x``. Values must be finite.
    x : scalar or array_like, shape (m,)
    axis : int, optional
        Axis along which `y` is assumed to be varying. Meaning that for
        ``x[i]`` the corresponding values are ``np.take(y, i, axis=axis)``.
        Default is 0.
    bc_type : string or 2-tuple, optional
        Boundary condition type. Two additional equations, given by the
        boundary conditions, are required to determine all coefficients of
        polynomials on each segment [2]_.
        If `bc_type` is a string, then the specified condition will be applied
        at both ends of a spline. Available conditions are:
        * 'not-a-knot' (default): The first and second segment at a curve end
          are the same polynomial. It is a good default when there is no
          information on boundary conditions.
        * 'periodic': The interpolated functions is assumed to be periodic
          of period ``x[-1] - x[0]``. The first and last value of `y` must be
          identical: ``y[0] == y[-1]``. This boundary condition will result in
          ``y'[0] == y'[-1]`` and ``y''[0] == y''[-1]``.
        * 'clamped': The first derivative at curves ends are zero. Assuming
          a 1D `y`, ``bc_type=((1, 0.0), (1, 0.0))`` is the same condition.
        * 'natural': The second derivative at curve ends are zero. Assuming
          a 1D `y`, ``bc_type=((2, 0.0), (2, 0.0))`` is the same condition.
        If `bc_type` is a 2-tuple, the first and the second value will be
        applied at the curve start and end respectively. The tuple values can
        be one of the previously mentioned strings (except 'periodic') or a
        tuple `(order, deriv_values)` allowing to specify arbitrary
        derivatives at curve ends:
        * `order`: the derivative order, 1 or 2.
        * `deriv_value`: array_like containing derivative values, shape must
          be the same as `y`, excluding ``axis`` dimension. For example, if
          `y` is 1D, then `deriv_value` must be a scalar. If `y` is 3D with
          the shape (n0, n1, n2) and axis=2, then `deriv_value` must be 2D
          and have the shape (n0, n1).
    extrapolate : {bool, 'periodic', None}, optional
        If bool, determines whether to extrapolate to out-of-bounds points
        based on first and last intervals, or to return NaNs. If 'periodic',
        periodic extrapolation is used. If None (default), ``extrapolate`` is
        set to 'periodic' for ``bc_type='periodic'`` and to True otherwise.

    See Also
    --------
    scipy.interpolate.CubicHermiteSpline

    Returns
    -------
    y : scalar or array_like
        The result, of shape (m,)

    References
    ----------
    .. [1] `Cubic Spline Interpolation
            <https://en.wikiversity.org/wiki/Cubic_Spline_Interpolation>`_
            on Wikiversity.
    .. [2] Carl de Boor, "A Practical Guide to Splines", Springer-Verlag, 1978.
    """
    from scipy import interpolate

    P = interpolate.CubicSpline(
        xi, yi, axis=axis, bc_type=bc_type, extrapolate=extrapolate
    )

    return P(x)


def _interpolate_with_limit_area(
    values: ArrayLike, method: str, limit: Optional[int], limit_area: Optional[str]
) -> ArrayLike:
    """
    Apply interpolation and limit_area logic to values along a to-be-specified axis.

    Parameters
    ----------
    values: array-like
        Input array.
    method: str
        Interpolation method. Could be "bfill" or "pad"
    limit: int, optional
        Index limit on interpolation.
    limit_area: str
        Limit area for interpolation. Can be "inside" or "outside"

    Returns
    -------
    values: array-like
        Interpolated array.
    """

    invalid = isna(values)

    if not invalid.all():
        first = find_valid_index(values, "first")
        last = find_valid_index(values, "last")

        values = interpolate_2d(
            values,
            method=method,
            limit=limit,
        )

        if limit_area == "inside":
            invalid[first : last + 1] = False
        elif limit_area == "outside":
            invalid[:first] = invalid[last + 1 :] = False

        values[invalid] = np.nan

    return values


def interpolate_2d(
    values,
    method: str = "pad",
    axis: Axis = 0,
    limit: Optional[int] = None,
    limit_area: Optional[str] = None,
):
    """
    Perform an actual interpolation of values, values will be make 2-d if
    needed fills inplace, returns the result.

       Parameters
    ----------
    values: array-like
        Input array.
    method: str, default "pad"
        Interpolation method. Could be "bfill" or "pad"
    axis: 0 or 1
        Interpolation axis
    limit: int, optional
        Index limit on interpolation.
    limit_area: str, optional
        Limit area for interpolation. Can be "inside" or "outside"

    Returns
    -------
    values: array-like
        Interpolated array.
    """
    if limit_area is not None:
        return np.apply_along_axis(
            partial(
                _interpolate_with_limit_area,
                method=method,
                limit=limit,
                limit_area=limit_area,
            ),
            axis,
            values,
        )

    orig_values = values

    transf = (lambda x: x) if axis == 0 else (lambda x: x.T)

    # reshape a 1 dim if needed
    ndim = values.ndim
    if values.ndim == 1:
        if axis != 0:  # pragma: no cover
            raise AssertionError("cannot interpolate on a ndim == 1 with axis != 0")
        values = values.reshape(tuple((1,) + values.shape))

    method = clean_fill_method(method)
    tvalues = transf(values)
    if method == "pad":
        result = _pad_2d(tvalues, limit=limit)
    else:
        result = _backfill_2d(tvalues, limit=limit)

    result = transf(result)
    # reshape back
    if ndim == 1:
        result = result[0]

    if orig_values.dtype.kind in ["m", "M"]:
        # convert float back to datetime64/timedelta64
        result = result.view(orig_values.dtype)

    return result


def _cast_values_for_fillna(values, dtype: DtypeObj, has_mask: bool):
    """
    Cast values to a dtype that algos.pad and algos.backfill can handle.
    """
    # TODO: for int-dtypes we make a copy, but for everything else this
    #  alters the values in-place.  Is this intentional?

    if needs_i8_conversion(dtype):
        values = values.view(np.int64)

    elif is_integer_dtype(values) and not has_mask:
        # NB: this check needs to come after the datetime64 check above
        # has_mask check to avoid casting i8 values that have already
        #  been cast from PeriodDtype
        values = ensure_float64(values)

    return values


def _fillna_prep(values, mask=None):
    # boilerplate for _pad_1d, _backfill_1d, _pad_2d, _backfill_2d
    dtype = values.dtype

    has_mask = mask is not None
    if not has_mask:
        # This needs to occur before datetime/timedeltas are cast to int64
        mask = isna(values)

    values = _cast_values_for_fillna(values, dtype, has_mask)

    mask = mask.view(np.uint8)
    return values, mask


def _pad_1d(values, limit=None, mask=None):
    values, mask = _fillna_prep(values, mask)
    algos.pad_inplace(values, mask, limit=limit)
    return values


def _backfill_1d(values, limit=None, mask=None):
    values, mask = _fillna_prep(values, mask)
    algos.backfill_inplace(values, mask, limit=limit)
    return values


def _pad_2d(values, limit=None, mask=None):
    values, mask = _fillna_prep(values, mask)

    if np.all(values.shape):
        algos.pad_2d_inplace(values, mask, limit=limit)
    else:
        # for test coverage
        pass
    return values


def _backfill_2d(values, limit=None, mask=None):
    values, mask = _fillna_prep(values, mask)

    if np.all(values.shape):
        algos.backfill_2d_inplace(values, mask, limit=limit)
    else:
        # for test coverage
        pass
    return values


_fill_methods = {"pad": _pad_1d, "backfill": _backfill_1d}


def get_fill_func(method):
    method = clean_fill_method(method)
    return _fill_methods[method]


def clean_reindex_fill_method(method):
    return clean_fill_method(method, allow_nearest=True)


def _interp_limit(invalid, fw_limit, bw_limit):
    """
    Get indexers of values that won't be filled
    because they exceed the limits.

    Parameters
    ----------
    invalid : boolean ndarray
    fw_limit : int or None
        forward limit to index
    bw_limit : int or None
        backward limit to index

    Returns
    -------
    set of indexers

    Notes
    -----
    This is equivalent to the more readable, but slower

    .. code-block:: python

        def _interp_limit(invalid, fw_limit, bw_limit):
            for x in np.where(invalid)[0]:
                if invalid[max(0, x - fw_limit):x + bw_limit + 1].all():
                    yield x
    """
    # handle forward first; the backward direction is the same except
    # 1. operate on the reversed array
    # 2. subtract the returned indices from N - 1
    N = len(invalid)
    f_idx = set()
    b_idx = set()

    def inner(invalid, limit):
        limit = min(limit, N)
        windowed = _rolling_window(invalid, limit + 1).all(1)
        idx = set(np.where(windowed)[0] + limit) | set(
            np.where((~invalid[: limit + 1]).cumsum() == 0)[0]
        )
        return idx

    if fw_limit is not None:

        if fw_limit == 0:
            f_idx = set(np.where(invalid)[0])
        else:
            f_idx = inner(invalid, fw_limit)

    if bw_limit is not None:

        if bw_limit == 0:
            # then we don't even need to care about backwards
            # just use forwards
            return f_idx
        else:
            b_idx_inv = list(inner(invalid[::-1], bw_limit))
            b_idx = set(N - 1 - np.asarray(b_idx_inv))
            if fw_limit == 0:
                return b_idx

    return f_idx & b_idx


def _rolling_window(a: np.ndarray, window: int):
    """
    [True, True, False, True, False], 2 ->

    [
        [True,  True],
        [True, False],
        [False, True],
        [True, False],
    ]
    """
    # https://stackoverflow.com/a/6811241
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)
