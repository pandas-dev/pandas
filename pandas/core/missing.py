"""
Routines for filling missing data.
"""

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


def clean_interp_method(method, **kwargs):
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
    xvalues,
    yvalues,
    method="linear",
    limit=None,
    limit_direction="forward",
    limit_area=None,
    fill_value=None,
    bounds_error=False,
    order=None,
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

    xvalues = getattr(xvalues, "values", xvalues)
    yvalues = getattr(yvalues, "values", yvalues)
    result = yvalues.copy()

    if method in ["linear", "time", "index", "values"]:
        if method in ("values", "index"):
            inds = np.asarray(xvalues)
            # hack for DatetimeIndex, #1646
            if needs_i8_conversion(inds.dtype.type):
                inds = inds.view(np.int64)
            if inds.dtype == np.object_:
                inds = lib.maybe_convert_objects(inds)
        else:
            inds = xvalues
        # np.interp requires sorted X values, #21037
        indexer = np.argsort(inds[valid])
        result[invalid] = np.interp(
            inds[invalid], inds[valid][indexer], yvalues[valid][indexer]
        )
        result[preserve_nans] = np.nan
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
    ]

    if method in sp_methods:
        inds = np.asarray(xvalues)
        # hack for DatetimeIndex, #1646
        if issubclass(inds.dtype.type, np.datetime64):
            inds = inds.view(np.int64)
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

    if getattr(x, "is_all_dates", False):
        # GH 5975, scipy.interp1d can't handle datetime64s
        x, new_x = x._values.astype("i8"), new_x.astype("i8")

    if method == "pchip":
        try:
            alt_methods["pchip"] = interpolate.pchip_interpolate
        except AttributeError as err:
            raise ImportError(
                "Your version of Scipy does not support PCHIP interpolation."
            ) from err
    elif method == "akima":
        alt_methods["akima"] = _akima_interpolate

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
        list of derivatives to extract. This numberincludes the function
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
    der : int or list, optional
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

    if der == 0:
        return P(x)
    elif interpolate._isscalar(der):
        return P(x, der=der)
    else:
        return [P(x, nu) for nu in der]


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
            b_idx = list(inner(invalid[::-1], bw_limit))
            b_idx = set(N - 1 - np.asarray(b_idx))
            if fw_limit == 0:
                return b_idx

    return f_idx & b_idx


def _rolling_window(a, window):
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
