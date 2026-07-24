import math
import numpy as np
from scipy.special import betainc
from scipy._lib._array_api import (
    xp_capabilities,
    xp_ravel,
    array_namespace,
    xp_promote,
    xp_device,
    xp_size,
    _count_nonmasked,
    is_torch,
    is_lazy_array,
)
import scipy._external.array_api_extra as xpx
from scipy.stats._axis_nan_policy import _broadcast_arrays, _contains_nan


def _quantile_iv(x, p, method, axis, nan_policy, keepdims, weights, fun='quantile'):
    xp = array_namespace(x, p, weights)

    if fun == "quantile":
        methods = {'inverted_cdf', 'averaged_inverted_cdf', 'closest_observation',
                   'hazen', 'interpolated_inverted_cdf', 'linear',
                   'median_unbiased', 'normal_unbiased', 'weibull',
                   'harrell-davis', '_lower', '_midpoint', '_higher', '_nearest',
                   'round_outward', 'round_inward', 'round_nearest'}
        allowed_types = 'real floating'
        def mask_fun(p): return (p > 1) | (p < 0) | xp.isnan(p)
        var2_name = 'p'
        var2_type_msg = '`p` must have real floating dtype.'
    else:
        methods = set(_estimated_cdf_methods)
        allowed_types = ('integral', 'real floating')
        mask_fun = xp.isnan
        var2_name = 'y'
        var2_type_msg = '`y` must have real dtype.'

    if not xp.isdtype(xp.asarray(x).dtype, ('integral', 'real floating')):
        raise ValueError("`x` must have real dtype.")

    if not xp.isdtype(xp.asarray(p).dtype, allowed_types):
        raise ValueError(var2_type_msg)

    if not (weights is None
            or xp.isdtype(xp.asarray(weights).dtype, ('integral', 'real floating'))):
        raise ValueError("`weights` must have real dtype.")

    x, p, weights = xp_promote(x, p, weights, force_floating=True, xp=xp)
    p = xp.asarray(p, device=xp_device(x))
    dtype = x.dtype

    axis_none = axis is None
    ndim = max(x.ndim, p.ndim)
    if axis_none:
        x = xp_ravel(x)
        p = xp_ravel(p)
        axis = 0
    elif np.iterable(axis) or int(axis) != axis:
        message = "`axis` must be an integer or None."
        raise ValueError(message)
    elif (axis >= ndim) or (axis < -ndim):
        message = "`axis` is not compatible with the shapes of the inputs."
        raise ValueError(message)
    axis = int(axis)

    if method not in methods:
        message = f"`method` must be one of {methods}"
        raise ValueError(message)

    no_weights = {'_lower', '_midpoint', '_higher', '_nearest', 'harrell-davis',
                  'round_nearest', 'round_inward', 'round_outward'}
    if weights is not None and method in no_weights:
        message = f"`method='{method}'` does not support `weights`."
        raise ValueError(message)

    contains_nans = _contains_nan(x, nan_policy, xp_omit_okay=True, xp=xp)

    if keepdims not in {None, True, False}:
        message = "If specified, `keepdims` must be True or False."
        raise ValueError(message)

    # If data has length zero along `axis`, the result will be an array of NaNs just
    # as if the data had length 1 along axis and were filled with NaNs. This is treated
    # naturally below whether `nan_policy` is `'propagate'` or `'omit'`.
    if x.shape[axis] == 0 and fun == 'quantile':
        shape = list(x.shape)
        shape[axis] = 1
        x = xp.full(shape, xp.nan, dtype=dtype, device=xp_device(x))

    if weights is None:
        y = xp.sort(x, axis=axis, stable=False)
        y, p = _broadcast_arrays((y, p), axis=axis)
        n_zero_weight = None
    else:
        x, weights = xp.broadcast_arrays(x, weights)
        i_zero_weight = (weights == 0)
        n_zero_weight = xp.count_nonzero(i_zero_weight, axis=axis, keepdims=True)
        x = xpx.at(x)[i_zero_weight].set(xp.inf, copy=True)
        i_y = xp.argsort(x, axis=axis, stable=False)
        y = xp.take_along_axis(x, i_y, axis=axis)
        weights = xp.take_along_axis(weights, i_y, axis=axis)
        y, p, weights, i_y, n_zero_weight = _broadcast_arrays(
            (y, p, weights, i_y, n_zero_weight), axis=axis)

    if (keepdims is False) and (p.shape[axis] != 1):
        message = ("`keepdims` may be False only if the length of "
                   f"`{var2_name}` along `axis` is 1.")
        raise ValueError(message)
    keepdims = (p.shape[axis] != 1) if keepdims is None else keepdims

    y = xp.moveaxis(y, axis, -1)
    p = xp.moveaxis(p, axis, -1)
    weights = weights if weights is None else xp.moveaxis(weights, axis, -1)
    n_zero_weight = (n_zero_weight if n_zero_weight is None
                     else xp.moveaxis(n_zero_weight, axis, -1))

    n = _count_nonmasked(y, -1, xp=xp, keepdims=True)
    n = n if n_zero_weight is None else n - n_zero_weight

    nan_out = None
    if is_lazy_array(y) or contains_nans:
        nans = xp.isnan(y)

        # Note that if length along `axis` were 0 to begin with,
        # it is now length 1 and filled with NaNs.
        if nan_policy == 'propagate':
            nan_out = xp.any(nans, axis=-1, keepdims=True)
        else:  # 'omit'
            n_int = n - xp.count_nonzero(nans, axis=-1, keepdims=True)
            n = xp.astype(n_int, dtype)
            # NaNs are produced only if slice is empty after removing NaNs
            nan_out = n == 0
            n = xpx.at(n, nan_out).set(y.shape[-1])  # avoids pytorch/pytorch#146211

        if (is_lazy_array(nans) or xp.any(nans)) and (method == 'harrell-davis'):
            y = xp.asarray(y, copy=True)  # ensure writable
            y = xpx.at(y, nans).set(0)  # any non-nan will prevent NaN from propagating
        if is_lazy_array(nan_out) or xp.any(nan_out):
            y = xp.asarray(y, copy=True)  # ensure writable
            y = xpx.at(y, xp.broadcast_to(nan_out, y.shape)).set(xp.nan)

    n = xp.asarray(n, dtype=dtype, device=xp_device(y))

    # apparently xpx.at is accepting nan_out as a mask even though it doesn't have the
    # same number of dimensions as `y`, yet it still appears to work correctly?
    # should refactor for clarity, and p_mask that gets returned here should probably
    # get a different name - `nan_out` would be more appropriate!
    p_mask = mask_fun(p) if nan_out is None else mask_fun(p) | nan_out
    if is_lazy_array(p_mask) or xp.any(p_mask):
        p = xp.where(p_mask, 0.5, p)  # these get NaN-ed out at the end

    return (y, p, method, axis, nan_policy, keepdims,
            n, axis_none, ndim, p_mask, weights, xp)


@xp_capabilities(skip_backends=[("dask.array", "No take_along_axis yet.")],
                 marray=True,
                 extra_note=("Use of `weights` is incompatible with MArray. "
                             "``method='harrell-davis'`` is CPU-only."))
def quantile(x, p, *, method='linear', axis=0, nan_policy='propagate', keepdims=None,
             weights=None):
    """
    Compute the p-th quantile of the data along the specified axis.

    Parameters
    ----------
    x : array_like of real numbers
        Data array.
    p : array_like of float
        Probability or sequence of probabilities of the quantiles to compute.
        Values must be between 0 and 1 (inclusive).
        While `numpy.quantile` can only compute quantiles according to the Cartesian
        product of the first two arguments, this function enables calculation of
        quantiles at different probabilities for each axis slice by following
        broadcasting rules like those of `scipy.stats` reducing functions.
        See `axis`, `keepdims`, and the examples.
    method : str, default: 'linear'
        The method to use for estimating the quantile.
        The available options, numbered as they appear in [1]_, are:

        1. 'inverted_cdf'
        2. 'averaged_inverted_cdf'
        3. 'closest_observation'
        4. 'interpolated_inverted_cdf'
        5. 'hazen'
        6. 'weibull'
        7. 'linear'  (default)
        8. 'median_unbiased'
        9. 'normal_unbiased'

        'harrell-davis' is also available to compute the quantile estimate
        according to [2]_.

        'round_outward', 'round_inward', and 'round_nearest' are available for use
        in trimming and winsorizing data.

        See Notes for details.
    axis : int or None, default: 0
        Axis along which the quantiles are computed.
        ``None`` ravels both `x` and `p` before performing the calculation,
        without checking whether the original shapes were compatible.
        As in other `scipy.stats` functions, a positive integer `axis` is resolved
        after prepending 1s to the shape of `x` or `p` as needed until the two arrays
        have the same dimensionality. When providing `x` and `p` with different
        dimensionality, consider using negative `axis` integers for clarity.
    nan_policy : str, default: 'propagate'
        Defines how to handle NaNs in the input data `x`.

        - ``propagate``: if a NaN is present in the axis slice (e.g. row) along
          which the  statistic is computed, the corresponding slice of the output
          will contain NaN(s).
        - ``omit``: NaNs will be omitted when performing the calculation.
          If insufficient data remains in the axis slice along which the
          statistic is computed, the corresponding slice of the output will
          contain NaN(s).
        - ``raise``: if a NaN is present, a ``ValueError`` will be raised.

        If NaNs are present in `p`, a ``ValueError`` will be raised.
    keepdims : bool, optional
        Consider the case in which `x` is 1-D and `p` is a scalar: the quantile
        is a reducing statistic, and the default behavior is to return a scalar.
        If `keepdims` is set to True, the axis will not be reduced away, and the
        result will be a 1-D array with one element.

        The general case is more subtle, since multiple quantiles may be
        requested for each axis-slice of `x`. For instance, if both `x` and `p`
        are 1-D and ``p.size > 1``, no axis can be reduced away; there must be an
        axis to contain the number of quantiles given by ``p.size``. Therefore:

        - By default, the axis will be reduced away if possible (i.e. if there is
          exactly one element of `p` per axis-slice of `x`).
        - If `keepdims` is set to True, the axis will not be reduced away.
        - If `keepdims` is set to False, the axis will be reduced away
          if possible, and an error will be raised otherwise.

    weights : array_like of finite, non-negative real numbers
        Frequency weights; e.g., for counting number weights,
        ``quantile(x, p, weights=weights)`` is equivalent to
        ``quantile(np.repeat(x, weights), p)``. Values other than finite counting
        numbers are accepted, but may not have valid statistical interpretations.
        Not compatible with ``method='harrell-davis'`` or those that begin with
        ``'round_'``.

    Returns
    -------
    quantile : scalar or ndarray
        The resulting quantile(s). The dtype is the result dtype of `x` and `p`.

    See Also
    --------
    numpy.quantile
    estimated_cdf
    :ref:`outliers`

    Notes
    -----
    Given a sample `x` from an underlying distribution, `quantile` provides a
    nonparametric estimate of the inverse cumulative distribution function.

    By default, this is done by interpolating between adjacent elements in
    ``y``, a sorted copy of `x`::

        (1-g)*y[j] + g*y[j+1]

    where the index ``j`` and coefficient ``g`` are the integral and
    fractional components of ``p * (n-1)``, and ``n`` is the number of
    elements in the sample.

    This is a special case of Equation 1 of H&F [1]_. More generally,

    - ``j = (p*n + m - 1) // 1``, and
    - ``g = (p*n + m - 1) % 1``,

    where ``m`` may be defined according to several different conventions.
    The preferred convention may be selected using the ``method`` parameter:

    =============================== =============== ===============
    ``method``                      number in H&F   ``m``
    =============================== =============== ===============
    ``interpolated_inverted_cdf``   4               ``0``
    ``hazen``                       5               ``1/2``
    ``weibull``                     6               ``p``
    ``linear`` (default)            7               ``1 - p``
    ``median_unbiased``             8               ``p/3 + 1/3``
    ``normal_unbiased``             9               ``p/4 + 3/8``
    =============================== =============== ===============

    Note that indices ``j`` and ``j + 1`` are clipped to the range ``0`` to
    ``n - 1`` when the results of the formula would be outside the allowed
    range of non-negative indices. When ``j`` is clipped to zero, ``g`` is
    set to zero as well. The ``-1`` in the formulas for ``j`` and ``g``
    accounts for Python's 0-based indexing.

    The table above includes only the estimators from [1]_ that are continuous
    functions of probability `p` (estimators 4-9). SciPy also provides the
    three discontinuous estimators from [1]_ (estimators 1-3), where ``j`` is
    defined as above, ``m`` is defined as follows, and ``g`` is ``0`` when
    ``index = p*n + m - 1`` is less than ``0`` and otherwise is defined below.

    1. ``inverted_cdf``: ``m = 0`` and ``g = int(index - j > 0)``
    2. ``averaged_inverted_cdf``: ``m = 0`` and
       ``g = (1 + int(index - j > 0)) / 2``
    3. ``closest_observation``: ``m = -1/2`` and
       ``g = 1 - int((index == j) & (j%2 == 1))``

    Note that for methods ``inverted_cdf`` and ``averaged_inverted_cdf``, only the
    relative proportions of tied observations (and relative weights) affect the
    results; for all other methods, the total number of observations (and absolute
    weights) matter.

    A different strategy for computing quantiles from [2]_, ``method='harrell-davis'``,
    uses a weighted combination of all elements. The weights are computed as:

    .. math::

        w_{n, i} = I_{i/n}(a, b) - I_{(i - 1)/n}(a, b)

    where :math:`n` is the number of elements in the sample,
    :math:`i` are the indices :math:`1, 2, ..., n-1, n` of the sorted elements,
    :math:`a = p (n + 1)`, :math:`b = (1 - p)(n + 1)`,
    :math:`p` is the probability of the quantile, and
    :math:`I` is the regularized, lower incomplete beta function
    (`scipy.special.betainc`).

    ``method='round_nearest'`` is equivalent to indexing ``y[j]``, where::

        j = int(np.round(p*n) if p < 0.5 else np.round(n*p - 1))

    This is useful when winsorizing data: replacing ``p*n`` of the most extreme
    observations with the next most extreme observation. ``method='round_outward'``
    adjusts the direction of rounding to winsorize fewer elements::

        j = int(np.floor(p*n) if p < 0.5 else np.ceil(n*p - 1))

    and ``method='round_inward'`` rounds to winsorize more elements::

        j = int(np.ceil(p*n) if p < 0.5 else np.floor(n*p - 1))

    These methods are also useful for trimming data: removing ``p*n`` of the most
    extreme observations. See :ref:`outliers` for example applications.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy import stats
    >>> x = np.asarray([[10, 8, 7, 5, 4],
    ...                 [0, 1, 2, 3, 5]])

    Take the median of each row.

    >>> stats.quantile(x, 0.5, axis=-1)
    array([7.,  2.])

    Take a different quantile for each row.

    >>> stats.quantile(x, [[0.25], [0.75]], axis=-1, keepdims=True)
    array([[5.],
           [3.]])

    Take multiple quantiles for each row.

    >>> stats.quantile(x, [0.25, 0.75], axis=-1)
    array([[5., 8.],
           [1., 3.]])

    Take different quantiles for each row.

    >>> p = np.asarray([[0.25, 0.75],
    ...                 [0.5, 1.0]])
    >>> stats.quantile(x, p, axis=-1)
    array([[5., 8.],
           [2., 5.]])

    Take different quantiles for each column.

    >>> stats.quantile(x.T, p.T, axis=0)
    array([[5., 2.],
           [8., 5.]])

    References
    ----------
    .. [1] R. J. Hyndman and Y. Fan,
       "Sample quantiles in statistical packages,"
       The American Statistician, 50(4), pp. 361-365, 1996
    .. [2] Harrell, Frank E., and C. E. Davis.
       "A new distribution-free quantile estimator."
       Biometrika 69.3 (1982): 635-640.

    """
    # Input validation / standardization

    temp = _quantile_iv(x, p, method, axis, nan_policy, keepdims, weights)
    (y, p, method, axis, nan_policy, keepdims,
     n, axis_none, ndim, p_mask, weights, xp) = temp

    if method in {'inverted_cdf', 'averaged_inverted_cdf', 'closest_observation',
                  'hazen', 'interpolated_inverted_cdf', 'linear',
                  'median_unbiased', 'normal_unbiased', 'weibull'}:
        res = _quantile_hf(y, p, n, method, weights, xp)
    elif method in {'harrell-davis'}:
        res = _quantile_hd(y, p, n, xp)
    elif method in {'_lower', '_midpoint', '_higher', '_nearest'}:
        res = _quantile_bc(y, p, n, method, xp)
    else:  # method.startswith('round'):
        res = _quantile_winsor(y, p, n, method, xp)

    return _post_quantile(res, p_mask, axis, axis_none, ndim, keepdims, xp)


def _post_quantile(res, p_mask, axis, axis_none, ndim, keepdims, xp):
    res = xpx.at(res, p_mask).set(xp.nan)

    # Reshape per axis/keepdims
    if axis_none and keepdims:
        shape = (1,)*(ndim - 1) + res.shape
        res = xp.reshape(res, shape)
        axis = -1

    res = xp.moveaxis(res, -1, axis)

    if not keepdims:
        res = xp.squeeze(res, axis=axis)

    return res[()] if res.ndim == 0 else res


def _quantile_hf(y, p, n, method, weights, xp):
    ms = dict(inverted_cdf=0, averaged_inverted_cdf=0, closest_observation=-0.5,
              interpolated_inverted_cdf=0, hazen=0.5, weibull=p, linear=1 - p,
              median_unbiased=p/3 + 1/3, normal_unbiased=p/4 + 3/8)
    m = ms[method]

    if weights is None:
        jg = p * n + m
        jp1 = jg // 1
        j = jp1 - 1
    else:
        cumulative_weights = xp.cumulative_sum(weights, axis=-1)
        n_int = xp.asarray(n, dtype=xp.int64)
        n_int = xp.broadcast_to(n_int, cumulative_weights.shape[:-1] + (1,))
        total_weight = xp.take_along_axis(cumulative_weights, n_int-1, axis=-1)
        jg = p * total_weight + m
        jp1 = _xp_searchsorted(cumulative_weights, jg, side='right')
        j = _xp_searchsorted(cumulative_weights, jg-1, side='right')
        j, jp1 = xp.astype(j, y.dtype), xp.astype(jp1, y.dtype)

    g = jg % 1
    if method == 'inverted_cdf':
        g = xp.astype((g > 0), jg.dtype)
    elif method == 'averaged_inverted_cdf':
        g = (1 + xp.astype((g > 0), jg.dtype)) / 2
    elif method == 'closest_observation':
        g = (1 - xp.astype((g == 0) & (j % 2 == 1), jg.dtype))
    if method in {'inverted_cdf', 'averaged_inverted_cdf', 'closest_observation'}:
        g = xp.asarray(g)
        g = xpx.at(g, jg < 0).set(0)

    g = xpx.at(g)[j < 0].set(0)
    j = xp.clip(j, 0., n - 1)
    jp1 = xp.clip(jp1, 0., n - 1)

    return ((1 - g) * xp.take_along_axis(y, xp.astype(j, xp.int64), axis=-1)
            + g * xp.take_along_axis(y, xp.astype(jp1, xp.int64), axis=-1))


def _quantile_hd(y, p, n, xp):
    # RE axis handling: We need to perform a reducing operation over rows of `y` for
    # each element in the corresponding row of `p` (a la Cartesian product). Strategy:
    # move rows of `p` to an axis at the front that is orthogonal to all the rest,
    # perform the reducing operating over the last axis, then move the front axis back
    # to the end.
    p = xp.moveaxis(p, -1, 0)[..., xp.newaxis]
    a = p * (n + 1)
    b = (1 - p) * (n + 1)
    i = xp.arange(y.shape[-1] + 1, dtype=y.dtype, device=xp_device(y))
    w = betainc(a, b, i / n)
    w = w[..., 1:] - w[..., :-1]
    w = xpx.at(w, xp.isnan(w)).set(0)
    res = xp.vecdot(w, y, axis=-1)
    return xp.moveaxis(res, 0, -1)


def _quantile_winsor(y, p, n, method, xp):
    ops = dict(round_outward=(xp.floor, xp.ceil),
               round_inward=(xp.ceil, xp.floor),
               round_nearest=(xp.round, xp.round))
    op_left, op_right = ops[method]
    j = xp.where(p < 0.5, op_left(p*n), op_right(n*p - 1))
    return xp.take_along_axis(y, xp.astype(j, xp.int64), axis=-1)


def _quantile_bc(y, p, n, method, xp):
    # Methods retained for backward compatibility. NumPy documentation is not
    # quite right about what these methods do: if `p * (n - 1)` is integral,
    # that is used as the index. See numpy/numpy#28910.
    ij = p * (n - 1)
    if method == '_midpoint':
        return (xp.take_along_axis(y, xp.astype(xp.floor(ij), xp.int64), axis=-1)
                + xp.take_along_axis(y, xp.astype(xp.ceil(ij), xp.int64), axis=-1)) / 2
    elif method == '_lower':
        k = xp.floor(ij)
    elif method == '_higher':
        k = xp.ceil(ij)
    elif method == '_nearest':
        k = xp.round(ij)
    return xp.take_along_axis(y, xp.astype(k, xp.int64), axis=-1)


@xp_capabilities(skip_backends=[("dask.array", "No take_along_axis yet.")])
def _xp_searchsorted(x, y, *, side='left', xp=None):
    # Vectorized xp.searchsorted. Search is performed along last axis. The shape of the
    # output is that of `y`, broadcasting the batch dimensions with those of `x` if
    # necessary.
    xp = array_namespace(x, y) if xp is None else xp
    xp_default_int = xp.asarray(1).dtype
    y_0d = xp.asarray(y).ndim == 0
    x, y = _broadcast_arrays((x, y), axis=-1, xp=xp)
    x_1d = x.ndim <= 1

    if x_1d or is_torch(xp):
        y = xp.reshape(y, ()) if (y_0d and x_1d) else y
        out = xp.searchsorted(x, y, side=side)
        out = xp.astype(out, xp_default_int, copy=False)
        return out

    a = xp.full(y.shape, 0, device=xp_device(x))

    if x.shape[-1] == 0:
        return a

    n = xp.count_nonzero(~xp.isnan(x), axis=-1, keepdims=True)
    b = xp.broadcast_to(n, y.shape)

    compare = xp.less_equal if side == 'left' else xp.less

    # while xp.any(b - a > 1):
    # refactored to for loop with ~log2(n) iterations for JAX JIT
    for i in range(int(math.log2(x.shape[-1])) + 1):
        c = (a + b) // 2
        x0 = xp.take_along_axis(x, c, axis=-1)
        j = compare(y, x0)
        b = xp.where(j, c, b)
        a = xp.where(j, a, c)

    out = xp.where(compare(y, xp.min(x, axis=-1, keepdims=True)), 0, b)
    out = xp.where(xp.isnan(y), x.shape[-1], out) if side == 'right' else out
    out = xp.astype(out, xp_default_int, copy=False)
    return out


@xp_capabilities(skip_backends=[("dask.array", "No take_along_axis yet.")],
                 marray=True)
def estimated_cdf(x, y, *, method='linear',
                  axis=0, nan_policy='propagate', keepdims=None):
    """
    Estimate the CDF of the distribution underlying a sample.

    Parameters
    ----------
    x : array_like of real numbers
        Data array.
    y : array_like of real numbers
        Points at which to evaluate the estimated cdf.
        See `axis`, `keepdims`, and the examples for broadcasting behavior.
    method : str, default: 'linear'
        The method to use for estimating the cumulative distribution function.
        The available options, numbered as they appear in [1]_, are:

        1. 'inverted_cdf' (AKA *the* empirical cumulative distribution function)
        2. 'averaged_inverted_cdf'
        3. 'closest_observation'
        4. 'interpolated_inverted_cdf'
        5. 'hazen'
        6. 'weibull'
        7. 'linear'  (default)
        8. 'median_unbiased'
        9. 'normal_unbiased'

        See Notes for details.
    axis : int or None, default: 0
        Axis along which samples in `x` are given in ND case.
        ``None`` ravels both `x` and `y` before performing the calculation,
        without checking whether the original shapes were compatible.
        As in other `scipy.stats` functions, a positive integer `axis` is resolved
        after prepending 1s to the shape of `x` or `y` as needed until the two arrays
        have the same dimensionality. When providing `x` and `y` with different
        dimensionality, consider using negative `axis` integers for clarity.
    nan_policy : str, default: 'propagate'
        Defines how to handle NaNs in the input data `x`.

        - ``propagate``: if a NaN is present in the axis slice (e.g. row) along
          which the  statistic is computed, the corresponding slice of the output
          will contain NaN(s).
        - ``omit``: NaNs will be omitted when performing the calculation.
          If insufficient data remains in the axis slice along which the
          statistic is computed, the corresponding slice of the output will
          contain NaN(s).
        - ``raise``: if a NaN is present, a ``ValueError`` will be raised.

        If NaNs are present in `y`, the corresponding entries in the output will be NaN.
    keepdims : bool, optional
        Consider the case in which `x` is 1-D and `y` is a scalar: the estimated
        cumulative distribution function at a point is a reducing statistic, and the
        default behavior is to return a scalar.
        If `keepdims` is set to True, the axis will not be reduced away, and the
        result will be a 1-D array with one element.

        The general case is more subtle, since multiple values of `y` may be
        specified for each axis-slice of `x`. For instance, if both `x` and `y`
        are 1-D and ``y.size > 1``, no axis can be reduced away; there must be an
        axis to contain the number of values given by ``y.size``. Therefore:

        - By default, the axis will be reduced away if possible (i.e. if there is
          exactly one element of `y` per axis-slice of `x`).
        - If `keepdims` is set to True, the axis will not be reduced away.
        - If `keepdims` is set to False, the axis will be reduced away
          if possible, and an error will be raised otherwise.

    Returns
    -------
    probability : scalar or ndarray
        The resulting probabilities(s). The dtype is the result dtype of `x`, `y`, and
        the Python ``float`` type.

    See Also
    --------
    quantile

    Notes
    -----
    Given a sample `x` from an underlying distribution, `estimated_cdf` provides a
    nonparametric estimate of the cumulative distribution function.

    By default, this is done by computing the "fractional index" ``p`` at which ``y``
    would appear within ``z``, a sorted copy of `x`::

        p = 1 / (n - 1) * (j +  (     y - z[j])
                              / (z[j+1] - z[j]))

    where the index ``j`` is that of the largest element of ``z`` that does not exceed
    ``y``, and ``n`` is the number of elements in the sample. Note that if ``y`` is an
    element of ``z``, then ``j`` is the index such that ``y = z[j]``, and the formula
    simplifies to the intuitive ``j / (n - 1)``. The full formula linearly interpolates
    between ``j / (n - 1)`` and ``(j + 1) / (n - 1)``.

    This is a special case of the more general::

        p = (j + (y - z[j]) / (z[j+1] - z[j] + 1 - m) / n

    where ``m`` may be defined according to several different conventions.
    The preferred convention may be selected using the ``method`` parameter:

    =============================== =============== ===============
    ``method``                      number in H&F   ``m``
    =============================== =============== ===============
    ``interpolated_inverted_cdf``   4               ``0``
    ``hazen``                       5               ``1/2``
    ``weibull``                     6               ``p``
    ``linear`` (default)            7               ``1 - p``
    ``median_unbiased``             8               ``p/3 + 1/3``
    ``normal_unbiased``             9               ``p/4 + 3/8``
    =============================== =============== ===============

    Note that indices ``j`` and ``j + 1`` are clipped to the range ``0`` to
    ``n - 1`` when the results of the formula would be outside the allowed
    range of non-negative indices, and resulting ``p`` is clipped to the range
    ``0`` to ``1``.

    The table above includes only the estimators from [1]_ that are continuous
    functions of probability `p` (estimators 4-9). SciPy also provides the
    three discontinuous estimators from [1]_ (estimators 1-3), where ``j`` is
    defined as above, ``m`` is defined as follows, and ``g`` is ``0`` when
    ``index = p*n + m - 1`` is less than ``0`` and otherwise is defined below.

    1. ``inverted_cdf``: ``m = 0`` and ``g = int(index - j > 0)``
    2. ``averaged_inverted_cdf``: ``m = 0`` and
       ``g = (1 + int(index - j > 0)) / 2``
    3. ``closest_observation``: ``m = -1/2`` and
       ``g = 1 - int((index == j) & (j%2 == 1))``

    When all the data in ``x`` are unique, `estimated_cdf` and `quantile` are are
    inverses of one another within a certain domain.
    Although `quantile` with ``method='linear'`` is invertible over the whole domain
    of ``p`` from ``0`` to ``1``, this is not true of other methods.

    References
    ----------
    .. [1] R. J. Hyndman and Y. Fan,
           "Sample quantiles in statistical packages,"
           The American Statistician, 50(4), pp. 361-365, 1996

    Examples
    --------
    >>> import numpy as np
    >>> from scipy import stats
    >>> x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    Estimate the cumulative distribution function of one sample at a single point.

    >>> stats.estimated_cdf(x, 5, axis=-1)
    np.float64(0.5)

    Estimate the cumulative distribution function of one sample at two points.

    >>> stats.estimated_cdf(x, [2.5, 7.5], axis=-1)
    array([0.25, 0.75])

    Estimate the cumulative distribution function of two samples at different points.

    >>> x = np.stack((np.arange(0, 11), np.arange(10, 21)))
    >>> stats.estimated_cdf(x, [[2.5], [17.5]], axis=-1, keepdims=True)
    array([[0.25],
           [0.75]])

    Estimate the cumulative distribution function at many points for each of two
    samples.

    >>> import matplotlib.pyplot as plt
    >>> rng = np.random.default_rng(6110515095)
    >>> x = stats.Normal(mu=[-1, 1]).sample(10000)
    >>> y = np.linspace(-4, 4, 5000)[:, np.newaxis]
    >>> p = stats.estimated_cdf(x, y, axis=0)
    >>> plt.plot(y, p)
    >>> plt.show()

    Note that the `quantile` and `estimated_cdf` functions are inverses of one another
    within a certain domain.

    >>> p = np.linspace(0, 1, 300)
    >>> x = rng.standard_normal(300)
    >>> y = stats.quantile(x, p)
    >>> p2 = stats.estimated_cdf(x, y)
    >>> np.testing.assert_allclose(p2, p)
    >>> y2 = stats.quantile(x, p2)
    >>> np.testing.assert_allclose(y2, y)

    However, the domain over which `quantile` can be inverted by `estimated_cdf` depends
    on the `method` used. This is most noticeable when there are few observations in the
    sample.

    >>> x = np.asarray([0, 1])
    >>> y_linear = stats.quantile(x, p, method='linear')
    >>> y_weibull = stats.quantile(x, p, method='weibull')
    >>> y_iicdf = stats.quantile(x, p, method='interpolated_inverted_cdf')
    >>> plt.plot(p, y_linear, p, y_weibull, p, y_iicdf)
    >>> plt.legend(['linear', 'weibull', 'iicdf'])
    >>> plt.xlabel('p')
    >>> plt.ylabel('y = quantile(x, p)')
    >>> plt.show()

    For example, in the case above, `quantile` is only invertible from
    ``p = 0.5`` to ``p = 1.0`` with ``method = 'interpolated_inverted_cdf'``. This is a
    fundamental characteristic of the methods, not a shortcoming of `estimated_cdf`.

    """
    temp = _quantile_iv(x, y, method, axis, nan_policy, keepdims, weights=None,
                        fun='estimated_cdf')
    x, y, method, axis, nan_policy, keepdims, n, axis_none, ndim, y_mask, _, xp = temp

    if xp_size(x) == 0:
        res = xp.full_like(y, xp.nan)
    else:
        res = _estimated_cdf_hf(x, y, n, method, xp)

    return _post_quantile(res, y_mask, axis, axis_none, ndim, keepdims, xp)


_estimated_cdf_discontinuous_methods = dict(
    inverted_cdf=0.0,
    averaged_inverted_cdf=0.0,
    closest_observation=0.5,
)


_estimated_cdf_continuous_methods = dict(
    interpolated_inverted_cdf=(0, 1),
    hazen=(0.5, 0.5),
    weibull=(0, 0),
    linear=(1, 1),
    median_unbiased=(1 / 3, 1 / 3),
    normal_unbiased=(3 / 8, 3 / 8),
)


_estimated_cdf_methods = (set(_estimated_cdf_continuous_methods).union(
                          set(_estimated_cdf_discontinuous_methods)))


def _estimated_cdf_hf(x, y, n, method, xp):
    n_int = xp.astype(n, xp.int64)
    j_max = n_int - 1
    j_min = xp.minimum(j_max, xp.asarray(1, dtype=j_max.dtype))
    jp1 = _xp_searchsorted(x, y, side='right')

    if method in _estimated_cdf_discontinuous_methods:
        dp = _estimated_cdf_discontinuous_methods[method]
        p = (xp.astype(jp1, x.dtype)+dp)/n

    else:
        jp1 = xp.clip(jp1, j_min, j_max)
        j = xp.clip(jp1-1, 0)
        xj = xp.take_along_axis(x, j, axis=-1)
        xjp1 = xp.take_along_axis(x, jp1, axis=-1)
        with np.errstate(divide='ignore', invalid='ignore'):  # refactor to apply_where?
            delta = xp.where((xjp1 > xj) & xp.isfinite(xj), (y - xj) / (xjp1 - xj), 1.)

        a, b = _estimated_cdf_continuous_methods[method]
        p = (xp.astype(jp1, x.dtype) + delta - a) / (n + 1 - a - b)

    xmin = x[..., :1]
    xmax = (x[..., -1:] if n.shape == () else xp.take_along_axis(x, j_max))
    p = xpx.at(p)[y < xmin].set(0.)
    p = xpx.at(p)[y > xmax].set(1.)

    return xp.clip(p, 0., 1.)
