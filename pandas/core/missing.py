"""
Routines for filling missing data.
"""

from typing import Any, Optional

import numpy as np

from pandas._libs import algos, lib
from pandas.compat._optional import import_optional_dependency

from pandas.core.dtypes.cast import infer_dtype_from_array
from pandas.core.dtypes.common import (
    ensure_float64,
    is_datetime64_dtype,
    is_datetime64tz_dtype,
    is_integer_dtype,
    is_numeric_v_string_like,
    is_scalar,
    is_timedelta64_dtype,
    needs_i8_conversion,
)
from pandas.core.dtypes.missing import isna


def mask_missing(arr, values_to_mask):
    """
    Return a masking array of same size/shape as arr
    with entries equaling any member of values_to_mask set to True
    """
    dtype, values_to_mask = infer_dtype_from_array(values_to_mask)

    try:
        values_to_mask = np.array(values_to_mask, dtype=dtype)

    except Exception:
        values_to_mask = np.array(values_to_mask, dtype=object)

    na_mask = isna(values_to_mask)
    nonna = values_to_mask[~na_mask]

    mask = None
    for x in nonna:
        if mask is None:
            if is_numeric_v_string_like(arr, x):
                # GH#29553 prevent numpy deprecation warnings
                mask = False
            else:
                mask = arr == x

            # if x is a string and arr is not, then we get False and we must
            # expand the mask to size arr.shape
            if is_scalar(mask):
                mask = np.zeros(arr.shape, dtype=bool)
        else:
            if is_numeric_v_string_like(arr, x):
                # GH#29553 prevent numpy deprecation warnings
                mask |= False
            else:
                mask |= arr == x

    if na_mask.any():
        if mask is None:
            mask = isna(arr)
        else:
            mask |= isna(arr)

    # GH 21977
    if mask is None:
        mask = np.zeros(arr.shape, dtype=bool)

    return mask


def clean_fill_method(method, allow_nearest=False):
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


def clean_interp_method(method: str, **kwargs) -> str:
    order = kwargs.get("order")
    valid = [
        "linear",
        "time",
        "index",
        "values",
        "nearest",
        "zero",
        "slinear",
        "quadratic",
        "cubic",
        "barycentric",
        "polynomial",
        "krogh",
        "piecewise_polynomial",
        "pchip",
        "akima",
        "spline",
        "from_derivatives",
        "cubicspline",
    ]
    if method in ("spline", "polynomial") and order is None:
        raise ValueError("You must specify the order of the spline or polynomial.")
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
    xvalues: np.ndarray,
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
    # Treat the original, non-scipy methods first.

    invalid = isna(yvalues)
    valid = ~invalid

    if not valid.any():
        # have to call np.asarray(xvalues) since xvalues could be an Index
        # which can't be mutated
        result = np.empty_like(np.asarray(xvalues), dtype=np.float64)
        result.fill(np.nan)
        return result

    if valid.all():
        return yvalues

    if method == "time":
        if not getattr(xvalues, "is_all_dates", None):
            # if not issubclass(xvalues.dtype.type, np.datetime64):
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
    limit = algos._validate_limit(nobs=None, limit=limit)

    if limit_direction == "forward":
        nans_to_interpolate = _interp_limit(invalid, limit, 0)
    elif limit_direction == "backward":
        nans_to_interpolate = _interp_limit(invalid, 0, limit)
    else:
        # both directions... just use _interp_limit
        nans_to_interpolate = _interp_limit(invalid, limit, limit)

    # if limit_area is set, add either mid or outside indices
    # to preserve_nans GH #16284
    if limit_area:
        first = find_valid_index(yvalues, "first")
        last = find_valid_index(yvalues, "last")
        if limit_area == "inside":
            # preserve NaNs on the outside
            nans_to_interpolate[:first] = False
            nans_to_interpolate[last + 1 :] = False
        else:
            # preserve NaNs on the inside
            nans_to_interpolate[first : last + 1] = False

    xvalues = getattr(xvalues, "values", xvalues)
    yvalues = getattr(yvalues, "values", yvalues)
    result = yvalues.copy()

    if method in ["linear", "time", "index", "values"]:
        if method in ("values", "index"):
            inds = np.asarray(xvalues)
            # hack for DatetimeIndex, #1646
            if needs_i8_conversion(inds.dtype):
                inds = inds.view(np.int64)
            if inds.dtype == np.object_:
                inds = lib.maybe_convert_objects(inds)
        else:
            inds = xvalues
        # np.interp requires sorted X values, #21037
        indexer = np.argsort(inds[valid])
        result[nans_to_interpolate] = np.interp(
            inds[nans_to_interpolate], inds[valid][indexer], yvalues[valid][indexer]
        )
        return result

    sp_methods = [
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

    if method in sp_methods:
        inds = np.asarray(xvalues)
        # hack for DatetimeIndex, #1646
        if issubclass(inds.dtype.type, np.datetime64):
            inds = inds.view(np.int64)
        result[nans_to_interpolate] = _interpolate_scipy_wrapper(
            inds[valid],
            yvalues[valid],
            inds[nans_to_interpolate],
            method=method,
            fill_value=fill_value,
            bounds_error=bounds_error,
            order=order,
            **kwargs,
        )
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

    if getattr(x, "is_all_dates", False):
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


def interpolate_2d(
    values, method="pad", axis=0, limit=None, fill_value=None, dtype=None
):
    """
    Perform an actual interpolation of values, values will be make 2-d if
    needed fills inplace, returns the result.
    """
    orig_values = values

    transf = (lambda x: x) if axis == 0 else (lambda x: x.T)

    # reshape a 1 dim if needed
    ndim = values.ndim
    if values.ndim == 1:
        if axis != 0:  # pragma: no cover
            raise AssertionError("cannot interpolate on a ndim == 1 with axis != 0")
        values = values.reshape(tuple((1,) + values.shape))

    if fill_value is None:
        mask = None
    else:  # todo create faster fill func without masking
        mask = mask_missing(transf(values), fill_value)

    method = clean_fill_method(method)
    if method == "pad":
        values = transf(pad_2d(transf(values), limit=limit, mask=mask, dtype=dtype))
    else:
        values = transf(
            backfill_2d(transf(values), limit=limit, mask=mask, dtype=dtype)
        )

    # reshape back
    if ndim == 1:
        values = values[0]

    if orig_values.dtype.kind == "M":
        # convert float back to datetime64
        values = values.astype(orig_values.dtype)

    return values


def _cast_values_for_fillna(values, dtype):
    """
    Cast values to a dtype that algos.pad and algos.backfill can handle.
    """
    # TODO: for int-dtypes we make a copy, but for everything else this
    #  alters the values in-place.  Is this intentional?

    if (
        is_datetime64_dtype(dtype)
        or is_datetime64tz_dtype(dtype)
        or is_timedelta64_dtype(dtype)
    ):
        values = values.view(np.int64)

    elif is_integer_dtype(values):
        # NB: this check needs to come after the datetime64 check above
        values = ensure_float64(values)

    return values


def _fillna_prep(values, mask=None, dtype=None):
    # boilerplate for pad_1d, backfill_1d, pad_2d, backfill_2d
    if dtype is None:
        dtype = values.dtype

    if mask is None:
        # This needs to occur before datetime/timedeltas are cast to int64
        mask = isna(values)

    values = _cast_values_for_fillna(values, dtype)

    mask = mask.view(np.uint8)
    return values, mask


def pad_1d(values, limit=None, mask=None, dtype=None):
    values, mask = _fillna_prep(values, mask, dtype)
    algos.pad_inplace(values, mask, limit=limit)
    return values


def backfill_1d(values, limit=None, mask=None, dtype=None):
    values, mask = _fillna_prep(values, mask, dtype)
    algos.backfill_inplace(values, mask, limit=limit)
    return values


def pad_2d(values, limit=None, mask=None, dtype=None):
    values, mask = _fillna_prep(values, mask, dtype)

    if np.all(values.shape):
        algos.pad_2d_inplace(values, mask, limit=limit)
    else:
        # for test coverage
        pass
    return values


def backfill_2d(values, limit=None, mask=None, dtype=None):
    values, mask = _fillna_prep(values, mask, dtype)

    if np.all(values.shape):
        algos.backfill_2d_inplace(values, mask, limit=limit)
    else:
        # for test coverage
        pass
    return values


_fill_methods = {"pad": pad_1d, "backfill": backfill_1d}


def get_fill_func(method):
    method = clean_fill_method(method)
    return _fill_methods[method]


def clean_reindex_fill_method(method):
    return clean_fill_method(method, allow_nearest=True)


def _interp_limit(
    invalid: np.ndarray, fw_limit: Optional[int], bw_limit: Optional[int]
) -> np.ndarray:
    """
    Update mask to exclude elements not within limits

    Parameters
    ----------
    invalid : boolean ndarray
    fw_limit : int or None
        forward limit to index
    bw_limit : int or None
        backward limit to index

    Returns
    -------
    boolean ndarray

    Notes
    -----
    There follows a description of the implementation used for creating a mask
    for forward interpolation with a limit. To create a backwards fill, we first
    reverse the array and use the same algorithm.
    To fill in both directions we combine the masks from both forward and backwards
    fills.

    Say we start with the following array

    array([nan, nan,  1.,  3., nan, nan, nan, 11., nan, nan])

    create (or get from masked arrays) a boolean array of missing values

    >>> arr = pd.core.missing.isna(arr)
    >>> arr
    array([ True,  True, False, False,  True,  True,  True, False,  True,
            True])

    we convert the boolean array to integer array for counting the streaks

    >>> arr = arr.astype(int)
    >>> arr
    array([1, 1, 0, 0, 1, 1, 1, 0, 1, 1])

    cumsum will get us off to a good start, we store this as we will need this later

    >>> cumsum = arr.cumsum()
    >>> cumsum
    array([1, 2, 2, 2, 3, 4, 5, 5, 6, 7], dtype=int32)

    multiplying this accumulation with the original array of ones to get non-zero
    values where we originally had ones

    >>> arr = cumsum * arr
    >>> arr
    array([1, 2, 0, 0, 3, 4, 5, 0, 6, 7])

    the previous result is close to what we want, but we want to restart
    each streak at one. start by using the diff method to substract the previous
    value for each element

    >>> arr = np.diff(arr, prepend=0)
    >>> arr
    array([ 1,  1, -2,  0,  3,  1,  1, -5,  6,  1])

    a negative value now represents the end of a streak of missing values
    so let's first select just the negative values

    >>> arr = np.where(arr < 0, arr, 0)
    >>> arr
    array([ 0,  0, -2,  0,  0,  0,  0, -5,  0,  0])

    we will need to propegate the negative values

    >>> arr = np.minimum.accumulate(arr)
    >>> arr
    array([ 0,  0, -2, -2, -2, -2, -2, -5, -5, -5], dtype=int32)

    and then subtract the excess accumlation

    >>> arr = arr + cumsum
    >>> arr
    array([1, 2, 0, 0, 1, 2, 3, 0, 1, 2], dtype=int32)

    remember that positive values represent missing values and zeros represent
    valid values. We have a array with some missing values at the start. For a
    forward fill algorithm, we want to update the mask to leave these missing
    values unchanged.

    >>> arr[: arr.argmin()] = 0
    >>> arr
    array([0, 0, 0, 0, 1, 2, 3, 0, 1, 2], dtype=int32)

    we will now select only values within a set limit, say 2

    >>> arr = np.where(arr > 2, 0, arr)
    >>> arr
    array([0, 0, 0, 0, 1, 2, 0, 0, 1, 2], dtype=int32)

    and finally convert back into a boolean mask

    >>> arr.astype(bool)
    array([ False,  False, False, False,  True,  True, False, False,  True,
            True])
    """

    def inner(arr, limit):
        arr = arr.astype(int)
        arr[: arr.argmin()] = 0
        if limit:
            cumsum = arr.cumsum()
            arr = cumsum * arr
            arr = np.diff(arr, prepend=0)
            arr = np.where(arr < 0, arr, 0)
            arr = np.minimum.accumulate(arr)
            arr = arr + cumsum
            arr = np.where(arr > limit, 0, arr)
        return arr.astype(bool)

    if fw_limit == 0:
        f_idx = invalid
    else:
        f_idx = inner(invalid, fw_limit)

    if bw_limit == 0:
        # then we don't even need to care about backwards
        # just use forwards
        return f_idx
    else:
        b_idx = inner(invalid[::-1], bw_limit)[::-1]
        if fw_limit == 0:
            return b_idx

    return f_idx | b_idx
