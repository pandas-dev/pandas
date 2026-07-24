import numpy as np
import math
from scipy import stats, special
from scipy._external import array_api_extra as xpx
from scipy._lib._array_api import (xp_capabilities, array_namespace, xp_promote,
                                   is_numpy, _share_masks, _count_nonmasked, is_marray)
from scipy.stats._stats_py import (_SimpleNormal, SignificanceResult, _get_pvalue,
                                   _rankdata)
from scipy.stats._axis_nan_policy import _axis_nan_policy_factory
from scipy.stats._stats_mstats_common import (TheilslopesResult, _n_samples_optional_x,
                                              SiegelslopesResult)


__all__ = ['chatterjeexi', 'spearmanrho', 'theilslopes', 'siegelslopes']


def _xi_statistic(x, y, y_continuous, xp):
    # Compute xi correlation statistic

    # `axis=-1` is guaranteed by _axis_nan_policy decorator
    n = x.shape[-1]

    # "Rearrange the data as (X(1), Y(1)), . . . ,(X(n), Y(n))
    # such that X(1) ≤ ··· ≤ X(n)"
    j = xp.argsort(x, axis=-1)
    j, y = xp.broadcast_arrays(j, y)
    y = xp.take_along_axis(y, j, axis=-1)

    # "Let ri be the rank of Y(i), that is, the number of j such that Y(j) ≤ Y(i)"
    r = stats.rankdata(y, method='max', axis=-1)
    # " additionally define li to be the number of j such that Y(j) ≥ Y(i)"
    # Could probably compute this from r, but that can be an enhancement
    l = stats.rankdata(-y, method='max', axis=-1)

    num = xp.sum(xp.abs(xp.diff(r, axis=-1)), axis=-1)
    if y_continuous:  # [1] Eq. 1.1
        statistic = 1 - 3 * num / (n ** 2 - 1)
    else:  # [1] Eq. 1.2
        den = 2 * xp.sum((n - l) * l, axis=-1)
        statistic = 1 - n * num / den

    return statistic, r, l


def _xi_std(r, l, y_continuous, xp):
    # Compute asymptotic standard deviation of xi under null hypothesis of independence

    # `axis=-1` is guaranteed by _axis_nan_policy decorator
    n = r.shape[-1]

    # "Suppose that X and Y are independent and Y is continuous. Then
    # √n·ξn(X, Y) → N(0, 2/5) in distribution as n → ∞"
    if y_continuous:  # [1] Theorem 2.1
        return xp.asarray(math.sqrt(2 / 5) / math.sqrt(n), dtype=r.dtype)

    # "Suppose that X and Y are independent. Then √n·ξn(X, Y)
    # converges to N(0, τ²) in distribution as n → ∞
    # [1] Eq. 2.2 and surrounding math
    i = xp.arange(1, n + 1, dtype=r.dtype)
    u = xp.sort(r, axis=-1)
    v = xp.cumulative_sum(u, axis=-1)
    an = 1 / n**4 * xp.sum((2*n - 2*i + 1) * u**2, axis=-1)
    bn = 1 / n**5 * xp.sum((v + (n - i)*u)**2, axis=-1)
    cn = 1 / n**3 * xp.sum((2*n - 2*i + 1) * u, axis=-1)
    dn = 1 / n**3 * xp.sum((l * (n - l)), axis=-1)
    tau2 = (an - 2*bn + cn**2) / dn**2

    return xp.sqrt(tau2) / math.sqrt(n)


def _chatterjeexi_iv(y_continuous, method):
    # Input validation for `chatterjeexi`
    # x, y, `axis` input validation taken care of by decorator

    if y_continuous not in {True, False}:
        raise ValueError('`y_continuous` must be boolean.')

    if not isinstance(method, stats.PermutationMethod):
        method = method.lower()
        message = "`method` must be 'asymptotic' or a `PermutationMethod` instance."
        if method != 'asymptotic':
            raise ValueError(message)

    return y_continuous, method


def _unpack(res, _):
    return res.statistic, res.pvalue


@xp_capabilities(skip_backends=[('dask.array', 'no take_along_axis')])
@_axis_nan_policy_factory(SignificanceResult, paired=True, n_samples=2,
                          result_to_tuple=_unpack, n_outputs=2, too_small=1)
def chatterjeexi(x, y, *, axis=0, y_continuous=False, method='asymptotic'):
    r"""Compute the xi correlation and perform a test of independence.

    The xi correlation coefficient is a measure of association between two
    variables; the value tends to be close to zero when the variables are
    independent and close to 1 when there is a strong association. Unlike
    other correlation coefficients, the xi correlation is effective even
    when the association is not monotonic.

    Parameters
    ----------
    x, y : array-like
        The samples: corresponding observations of the independent and
        dependent variable. The (N-d) arrays must be broadcastable.
    axis : int, default: 0
        Axis along which to perform the test.
    y_continuous : bool, default: False
        Whether `y` is assumed to be drawn from a continuous distribution.
        If `y` is drawn from a continuous distribution, results are valid
        whether this is assumed or not, but enabling this assumption will
        result in faster computation and typically produce similar results.
    method : 'asymptotic' or `PermutationMethod` instance, optional
        Selects the method used to calculate the *p*-value.
        Default is 'asymptotic'. The following options are available.

        * ``'asymptotic'``: compares the standardized test statistic
          against the normal distribution.
        * `PermutationMethod` instance. In this case, the p-value
          is computed using `permutation_test` with the provided
          configuration options and other appropriate settings.

    Returns
    -------
    res : SignificanceResult
        An object containing attributes:

        statistic : float
            The xi correlation statistic.
        pvalue : float
            The associated *p*-value: the probability of a statistic at least as
            high as the observed value under the null hypothesis of independence.

    See Also
    --------
    scipy.stats.pearsonr, scipy.stats.spearmanr, scipy.stats.kendalltau

    Notes
    -----
    There is currently no special handling of ties in `x`; they are broken arbitrarily
    by the implementation. [1]_ recommends: "if there are ties among the Xi's, then
    choose an increasing rearrangement as above by breaking ties uniformly at random."
    This is easily accomplished by adding a small amount of random noise to `x`; see
    examples.

    [1]_ notes that the statistic is not symmetric in `x` and `y` *by design*:
    "...we may want to understand if :math:`Y` is a function :math:`X`, and not just
    if one of the variables is a function of the other." See [1]_ Remark 1.

    References
    ----------
    .. [1] Chatterjee, Sourav. "A new coefficient of correlation." Journal of
           the American Statistical Association 116.536 (2021): 2009-2022.
           :doi:`10.1080/01621459.2020.1758115`.

    Examples
    --------
    Generate perfectly correlated data, and observe that the xi correlation is
    nearly 1.0.

    >>> import numpy as np
    >>> from scipy import stats
    >>> rng = np.random.default_rng(348932549825235)
    >>> x = rng.uniform(0, 10, size=100)
    >>> y = np.sin(x)
    >>> res = stats.chatterjeexi(x, y)
    >>> res.statistic
    np.float64(0.9012901290129013)

    The probability of observing such a high value of the statistic under the
    null hypothesis of independence is very low.

    >>> res.pvalue
    np.float64(2.2206974648177804e-46)

    As noise is introduced, the correlation coefficient decreases.

    >>> noise = rng.normal(scale=[[0.1], [0.5], [1]], size=(3, 100))
    >>> res = stats.chatterjeexi(x, y + noise, axis=-1)
    >>> res.statistic
    array([0.79507951, 0.41824182, 0.16651665])

    Because the distribution of `y` is continuous, it is valid to pass
    ``y_continuous=True``. The statistic is identical, and the p-value
    (not shown) is only slightly different.

    >>> stats.chatterjeexi(x, y + noise, y_continuous=True, axis=-1).statistic
    array([0.79507951, 0.41824182, 0.16651665])

    Consider a case in which there are ties in `x`.

    >>> x = rng.integers(10, size=1000)
    >>> y = rng.integers(10, size=1000)

    [1]_ recommends breaking the ties uniformly at random.

    >>> d = rng.uniform(1e-5, size=x.size)
    >>> res = stats.chatterjeexi(x + d, y)
    >>> res.statistic
    -0.029919991638798438

    Since this gives a randomized estimate of the statistic, [1]_ also suggests
    considering the average over all possibilities of breaking ties. This is
    computationally infeasible when there are many ties, but a randomized estimate of
    *this* quantity can be obtained by considering many random possibilities of breaking
    ties.

    >>> d = rng.uniform(1e-5, size=(9999, x.size))
    >>> res = stats.chatterjeexi(x + d, y, axis=1)
    >>> np.mean(res.statistic)
    0.001186895213756626

    """
    xp = array_namespace(x, y)

    # x, y, `axis` input validation taken care of by decorator
    # In fact, `axis` is guaranteed to be -1
    y_continuous, method = _chatterjeexi_iv(y_continuous, method)
    x, y = xp_promote(x, y, force_floating=True, xp=xp)

    # A highly negative statistic is possible, e.g.
    # x = np.arange(100.), y = (x % 2 == 0)
    # Unclear whether we should expose `alternative`, though.
    alternative = 'greater'

    if method == 'asymptotic':
        xi, r, l = _xi_statistic(x, y, y_continuous, xp=xp)
        std = _xi_std(r, l, y_continuous, xp=xp)
        norm = _SimpleNormal()
        pvalue = _get_pvalue(xi / std, norm, alternative=alternative, xp=xp)
    elif isinstance(method, stats.PermutationMethod):
        res = stats.permutation_test(
            # Could be faster if we just permuted the ranks; for now, keep it simple.
            data=(y,),
            statistic=lambda y, axis: _xi_statistic(x, y, y_continuous, xp=xp)[0],
            alternative=alternative, permutation_type='pairings', **method._asdict(),
            axis=-1)  # `axis=-1` is guaranteed by _axis_nan_policy decorator

        xi, pvalue = res.statistic, res.pvalue

    xi = xi[()] if xi.ndim == 0 else xi
    pvalue = pvalue[()] if pvalue.ndim == 0 else pvalue
    return SignificanceResult(xi, pvalue)


@xp_capabilities(cpu_only=True, exceptions=['jax.numpy'], marray=True,
    skip_backends=[('dask.array', 'not supported by rankdata (take_along_axis)')],
    extra_note='Only the default `method` is compatible with MArray input.'
)
@_axis_nan_policy_factory(SignificanceResult, paired=True, n_samples=2,
                          result_to_tuple=_unpack, n_outputs=2, too_small=1)
def spearmanrho(x, y, /, *, alternative='two-sided', method=None, axis=0):
    r"""Calculate a Spearman rho correlation coefficient with associated p-value.

    The Spearman rank-order correlation coefficient is a nonparametric measure
    of the monotonicity of the relationship between two datasets.
    Like other correlation coefficients, it varies between -1 and +1 with 0
    implying no correlation. Coefficients of -1 or +1 are associated with an exact
    monotonic relationship.  Positive correlations indicate that as `x` increases,
    so does `y`; negative correlations indicate that as `x` increases, `y` decreases.
    The p-value is the probability of an uncorrelated system producing datasets
    with a Spearman correlation at least as extreme as the one computed from the
    observed dataset.

    Parameters
    ----------
    x, y : array-like
        The samples: corresponding observations of the independent and
        dependent variable. The (N-d) arrays must be broadcastable.
    alternative : {'two-sided', 'less', 'greater'}, optional
        Defines the alternative hypothesis. Default is 'two-sided'.
        The following options are available:

        * 'two-sided': the correlation is nonzero
        * 'less': the correlation is negative (less than zero)
        * 'greater':  the correlation is positive (greater than zero)

    method : ResamplingMethod, optional
        Defines the method used to compute the p-value. If `method` is an
        instance of `PermutationMethod`/`MonteCarloMethod`, the p-value is
        computed using
        `scipy.stats.permutation_test`/`scipy.stats.monte_carlo_test` with the
        provided configuration options and other appropriate settings.
        Otherwise, the p-value is computed using an asymptotic approximation of
        the null distribution.
    axis : int or None, optional
        If axis=0 (default), then each column represents a variable, with
        observations in the rows. If axis=1, the relationship is transposed:
        each row represents a variable, while the columns contain observations.
        If axis=None, then both arrays will be raveled.
        Like other `scipy.stats` functions, `axis` is interpreted after the
        arrays are broadcasted.

    Returns
    -------
    res : SignificanceResult
        An object containing attributes:

        statistic : floating point array or NumPy scalar
            Spearman correlation coefficient
        pvalue : floating point array NumPy scalar
            The p-value - the probabilitiy of realizing such an extreme statistic
            value under the null hypothesis that two samples have no ordinal
            correlation. See `alternative` above for alternative hypotheses.

    Warns
    -----
    `~scipy.stats.ConstantInputWarning`
        Raised if an input is a constant array.  The correlation coefficient
        is not defined in this case, so ``np.nan`` is returned.

    Notes
    -----
    `spearmanrho` was created to make improvements to SciPy's implementation of
    the Spearman correlation test without making backward-incompatible changes
    to `spearmanr`. Advantages of `spearmanrho` over `spearmanr` include:

    - `spearmanrho` follows standard array broadcasting rules.
    - `spearmanrho` is compatible with some non-NumPy arrays.
    - `spearmanrho` can compute exact p-values, even in the presence of ties,
      when an appropriate instance of `PermutationMethod` is provided via the
      `method` argument.

    References
    ----------
    .. [1] Zwillinger, D. and Kokoska, S. (2000). CRC Standard
       Probability and Statistics Tables and Formulae. Chapman & Hall: New
       York. 2000.
       Section  14.7
    .. [2] Kendall, M. G. and Stuart, A. (1973).
       The Advanced Theory of Statistics, Volume 2: Inference and Relationship.
       Griffin. 1973.
       Section 31.18

    Examples
    --------
    Univariate samples, approximate p-value.

    >>> import numpy as np
    >>> from scipy import stats
    >>> x = [1, 2, 3, 4, 5]
    >>> y = [5, 6, 7, 8, 7]
    >>> res = stats.spearmanrho(x, y)
    >>> res.statistic
    np.float64(0.8207826816681233)
    >>> res.pvalue
    np.float64(0.08858700531354405)

    Univariate samples, exact p-value.

    >>> res = stats.spearmanrho(x, y, method=stats.PermutationMethod())
    >>> res.statistic
    np.float64(0.8207826816681233)
    >>> res.pvalue
    np.float64(0.13333333333333333)

    Batch of univariate samples, one vectorized call.

    >>> rng = np.random.default_rng(98145152315484)
    >>> x2 = rng.standard_normal((2, 100))
    >>> y2 = rng.standard_normal((2, 100))
    >>> res = stats.spearmanrho(x2, y2, axis=-1)
    >>> res.statistic
    array([ 0.16585659, -0.12151215])
    >>> res.pvalue
    array([0.0991155 , 0.22846869])

    Bivariate samples using standard broadcasting rules.

    >>> res = stats.spearmanrho(x2[np.newaxis, :], x2[:, np.newaxis], axis=-1)
    >>> res.statistic
    array([[ 1.        , -0.14670267],
           [-0.14670267,  1.        ]])
    >>> res.pvalue
    array([[0.        , 0.14526128],
           [0.14526128, 0.        ]])

    """
    xp = array_namespace(x, y)
    x, y = _share_masks(x, y, xp=xp)
    rx = stats.rankdata(x, axis=axis)
    ry = stats.rankdata(y, axis=axis)
    res = stats.pearsonr(rx, ry, method=method, alternative=alternative, axis=axis)
    return SignificanceResult(res.statistic, res.pvalue)


@xp_capabilities(skip_backends=[("dask.array", "no take_along_axis"),
                                ("jax.numpy", "non-concrete boolean indexing")],
                 marray=True)
@_axis_nan_policy_factory(TheilslopesResult, default_axis=None, n_outputs=4,
                          n_samples=_n_samples_optional_x,
                          result_to_tuple=lambda x, _: tuple(x), paired=True,
                          too_small=1)
def theilslopes(y, x=None, alpha=0.95, method='separate', *, axis=None):
    r"""
    Computes the Theil-Sen estimator for a set of points (x, y).

    `theilslopes` implements a method for robust linear regression.  It
    computes the slope as the median of all slopes between paired values.

    Parameters
    ----------
    y : array_like
        Dependent variable.
    x : array_like or None, optional
        Independent variable. If None, use ``arange(len(y))`` instead.
    alpha : float, optional
        Confidence degree between 0 and 1. Default is 95% confidence.
        Note that `alpha` is symmetric around 0.5, i.e. both 0.1 and 0.9 are
        interpreted as "find the 90% confidence interval".
    method : {'joint', 'separate'}, optional
        Method to be used for computing estimate for intercept.
        Following methods are supported,

        * 'joint': Uses ``np.median(y - slope * x)`` as intercept.
        * 'separate': Uses ``np.median(y) - slope * np.median(x)``
                      as intercept.

        The default is 'separate'.

        .. versionadded:: 1.8.0

    axis : int or tuple of ints, default: None
        If an int or tuple of ints, the axis or axes of the input along which
        to compute the statistic. The statistic of each axis-slice (e.g. row)
        of the input will appear in a corresponding element of the output.
        If ``None``, the input will be raveled before computing the statistic.

    Returns
    -------
    result : ``TheilslopesResult`` instance
        The return value is an object with the following attributes:

        slope : float
            Theil slope.
        intercept : float
            Intercept of the Theil line.
        low_slope : float
            Lower bound of the confidence interval on `slope`.
        high_slope : float
            Upper bound of the confidence interval on `slope`.

    See Also
    --------
    siegelslopes : a similar technique using repeated medians

    Notes
    -----
    The implementation of `theilslopes` follows [1]_. The intercept is
    not defined in [1]_, and here it is defined as ``median(y) -
    slope*median(x)``, which is given in [3]_. Other definitions of
    the intercept exist in the literature such as  ``median(y - slope*x)``
    in [4]_. The approach to compute the intercept can be determined by the
    parameter ``method``. A confidence interval for the intercept is not
    given as this question is not addressed in [1]_.

    For compatibility with older versions of SciPy, the return value acts
    like a ``namedtuple`` of length 4, with fields ``slope``, ``intercept``,
    ``low_slope``, and ``high_slope``, so one can continue to write::

        slope, intercept, low_slope, high_slope = theilslopes(y, x)

    References
    ----------
    .. [1] P.K. Sen, "Estimates of the regression coefficient based on
           Kendall's tau", J. Am. Stat. Assoc., Vol. 63, pp. 1379-1389, 1968.
    .. [2] H. Theil, "A rank-invariant method of linear and polynomial
           regression analysis I, II and III",  Nederl. Akad. Wetensch., Proc.
           53:, pp. 386-392, pp. 521-525, pp. 1397-1412, 1950.
    .. [3] W.L. Conover, "Practical nonparametric statistics", 2nd ed.,
           John Wiley and Sons, New York, pp. 493.
    .. [4] https://en.wikipedia.org/wiki/Theil%E2%80%93Sen_estimator

    Examples
    --------
    >>> import numpy as np
    >>> from scipy import stats
    >>> import matplotlib.pyplot as plt

    >>> x = np.linspace(-5, 5, num=150)
    >>> y = x + np.random.normal(size=x.size)
    >>> y[11:15] += 10  # add outliers
    >>> y[-5:] -= 7

    Compute the slope, intercept and 90% confidence interval.  For comparison,
    also compute the least-squares fit with `linregress`:

    >>> res = stats.theilslopes(y, x, 0.90, method='separate')
    >>> lsq_res = stats.linregress(x, y)

    Plot the results. The Theil-Sen regression line is shown in red, with the
    dashed red lines illustrating the confidence interval of the slope (note
    that the dashed red lines are not the confidence interval of the regression
    as the confidence interval of the intercept is not included). The green
    line shows the least-squares fit for comparison.

    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.plot(x, y, 'b.')
    >>> ax.plot(x, res[1] + res[0] * x, 'r-')
    >>> ax.plot(x, res[1] + res[2] * x, 'r--')
    >>> ax.plot(x, res[1] + res[3] * x, 'r--')
    >>> ax.plot(x, lsq_res[1] + lsq_res[0] * x, 'g-')
    >>> plt.show()

    """
    return _robust_slopes(y, x=x, alpha=alpha, method=method, pfun='theilslopes')


@xp_capabilities(skip_backends=[("dask.array", "no take_along_axis"),
                                ("jax.numpy", "quantile needs lazy nan_policy")],
                 marray=True)
@_axis_nan_policy_factory(SiegelslopesResult, default_axis=None, n_outputs=2,
                          n_samples=_n_samples_optional_x,
                          result_to_tuple=lambda x, _: tuple(x), paired=True,
                          too_small=1)
def siegelslopes(y, x=None, method='hierarchical', *, axis=None):
    r"""
    Computes the Siegel estimator for a set of points (x, y).

    `siegelslopes` implements a method for robust linear regression
    using repeated medians (see [1]_) to fit a line to the points (x, y).
    The method is robust to outliers with an asymptotic breakdown point
    of 50%.

    Parameters
    ----------
    y : array_like
        Dependent variable.
    x : array_like or None, optional
        Independent variable. If None, use ``arange(len(y))`` instead.
    method : {'hierarchical', 'separate'}
        If 'hierarchical', estimate the intercept using the estimated
        slope ``slope`` (default option).
        If 'separate', estimate the intercept independent of the estimated
        slope. See Notes for details.
    axis : int or tuple of ints, default: None
        If an int or tuple of ints, the axis or axes of the input along which
        to compute the statistic. The statistic of each axis-slice (e.g. row)
        of the input will appear in a corresponding element of the output.
        If ``None``, the input will be raveled before computing the statistic.

    Returns
    -------
    result : ``SiegelslopesResult`` instance
        The return value is an object with the following attributes:

        slope : float
            Estimate of the slope of the regression line.
        intercept : float
            Estimate of the intercept of the regression line.

    See Also
    --------
    theilslopes : a similar technique without repeated medians

    Notes
    -----
    With ``n = len(y)``, compute ``m_j`` as the median of
    the slopes of the lines from the point ``(x[j], y[j])`` to all other ``n-1`` points.
    ``slope`` is then the median of all slopes ``m_j``.
    Two ways are given to estimate the intercept in [1]_ which can be chosen
    via the parameter ``method``.
    The hierarchical approach uses the estimated slope ``slope``
    and computes ``intercept`` as the median of ``y - slope*x``.
    The other approach estimates the intercept separately as follows: for
    each point ``(x[j], y[j])``, compute the intercepts of all the ``n-1``
    lines through the remaining points and take the median ``i_j``.
    ``intercept`` is the median of the ``i_j``.

    The implementation computes `n` times the median of a vector of size `n`
    which can be slow for large vectors. There are more efficient algorithms
    (see [2]_) which are not implemented here.

    For compatibility with older versions of SciPy, the return value acts
    like a ``namedtuple`` of length 2, with fields ``slope`` and
    ``intercept``, so one can continue to write::

        slope, intercept = siegelslopes(y, x)

    References
    ----------
    .. [1] A. Siegel, "Robust Regression Using Repeated Medians",
           Biometrika, Vol. 69, pp. 242-244, 1982.

    .. [2] A. Stein and M. Werman, "Finding the repeated median regression
           line", Proceedings of the Third Annual ACM-SIAM Symposium on
           Discrete Algorithms, pp. 409-413, 1992.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy import stats
    >>> import matplotlib.pyplot as plt

    >>> x = np.linspace(-5, 5, num=150)
    >>> y = x + np.random.normal(size=x.size)
    >>> y[11:15] += 10  # add outliers
    >>> y[-5:] -= 7

    Compute the slope and intercept.  For comparison, also compute the
    least-squares fit with `linregress`:

    >>> res = stats.siegelslopes(y, x)
    >>> lsq_res = stats.linregress(x, y)

    Plot the results. The Siegel regression line is shown in red. The green
    line shows the least-squares fit for comparison.

    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.plot(x, y, 'b.')
    >>> ax.plot(x, res[1] + res[0] * x, 'r-')
    >>> ax.plot(x, lsq_res[1] + lsq_res[0] * x, 'g-')
    >>> plt.show()

    """
    return _robust_slopes(y, x=x, method=method, pfun='siegelslopes')


def _robust_slopes(y, *, x, alpha=None, method, pfun):
    other_method = 'joint' if pfun == 'theilslopes' else 'hierarchical'
    if method not in {other_method, 'separate'}:
        raise ValueError(f"method must be either '{other_method}' or 'separate'. "
                         f"'{method}' is invalid.")

    xp = array_namespace(y, x)
    y, x = xp_promote(y, x, force_floating=True, xp=xp)
    x = xp.arange(y.shape[-1], dtype=y.dtype) if x is None else x
    y, x = xp.broadcast_arrays(y, x)
    x, y = _share_masks(x, y, xp=xp)

    if x.shape[-1] < 2 or y.shape[-1] < 2:  # only needed by test_axis_nan_policy
        raise ValueError("`x` and `y` must have length at least 2.")

    # Compute sorted slopes only when deltax > 0
    deltax = x[..., :, xp.newaxis] - x[..., xp.newaxis, :]
    deltay = y[..., :, xp.newaxis] - y[..., xp.newaxis, :]

    if pfun == 'theilslopes':
        i = xp.astype(xp.triu(xp.ones(deltax.shape[-2:]), k=1), xp.bool)
        if is_numpy(xp):
            deltax, deltay = deltax[..., i], deltay[..., i]
        else:
            # With array API, mask must be sole index, so we need to broadcast it.
            # Indexing ravels the array, so we need to reshape the results.
            i = xp.broadcast_to(i, deltax.shape)
            deltax, deltay = deltax[i], deltay[i]
            deltax = xp.reshape(deltax, (*i.shape[:-2], -1))
            deltay = xp.reshape(deltay, (*i.shape[:-2], -1))

    # Use `.data` to avoid indexing with masked indices. If masked elements of deltax=0,
    # they can safely be set to xp.nan, too - they're masked, after all.
    deltax0 = (deltax.data == 0) if is_marray(xp) else (deltax == 0)
    deltax = xpx.at(deltax)[deltax0].set(xp.nan)
    slopes = deltay / deltax

    if is_marray(xp):
        def nanmedian(x, axis):  # use mask to ignore NaNs
            x_nans_masked = xp.asarray(x.data, mask=x.mask | x._xp.isnan(x.data))
            return stats.quantile(x_nans_masked, 0.5, axis=axis)
    else:
        def nanmedian(x, axis):  # use nan_policy to ignore NaNs
            return stats.quantile(x, 0.5, axis=axis, nan_policy='omit')

    def median(x, axis): return stats.quantile(x, 0.5, axis=axis)

    # `theilslopes` is median of all slopes. Indexing with `i` above has already raveled
    # all the slopes, so we only need to take the median along the last axis.
    # `siegelslope` is a median of medians: we take the median of slopes from point i
    # to all other points, then the median of those medians.
    # The slope is NaN wherever the two points are the same, and those don't contribute,
    # hence the first median omits NaNs. NaNs in the input are propagated by the
    # `_axis_nan_policy` decorator.
    medslope = (median(nanmedian(slopes, axis=-1), axis=-1) if pfun == 'siegelslopes'
                else nanmedian(slopes, axis=-1))

    if method in {'joint', 'hierarchical'}:
        medinter = median(y - medslope[..., np.newaxis] * x, axis=-1)
    elif pfun == 'theilslopes':
        medinter = median(y, axis=-1) - medslope * median(x, axis=-1)
    else:
        # Calculate pairwise intercepts given each point (row i) and the slope to each
        # other point (column j). Then calculate the median of (row) medians.
        intercepts = y[..., :, xp.newaxis] - slopes*x[..., :, xp.newaxis]
        medinter = median(nanmedian(intercepts, axis=-1), axis=-1)

    if pfun == 'siegelslopes':
        return SiegelslopesResult(slope=medslope[()], intercept=medinter[()])

    # Now compute confidence intervals
    if alpha > 0.5:
        alpha = 1. - alpha

    z = float(special.ndtri(alpha / 2.))
    # This implements (2.6) from Sen (1968)
    # we don't actually need ranks, so an enhancement could be to have
    # `rankdata` return only the third output. In the meantime, use the
    # least expensive `method`.
    _, _, nxreps = _rankdata(x, method='min', return_ties=True)
    _, _, nyreps = _rankdata(y, method='min', return_ties=True)
    nt = xp.count_nonzero(xp.isfinite(slopes), axis=-1, keepdims=True)  # N in Sen 1968
    nt = xp.asarray(nt, dtype=y.dtype)
    ny = _count_nonmasked(y, keepdims=True, axis=-1)                    # n in Sen 1968
    # Equation 2.6 in Sen (1968):
    sigsq = 1/18. * (
        ny * (ny-1) * (2*ny+5)
        - xp.sum(nxreps * (nxreps-1) * (2*nxreps + 5), axis=-1, keepdims=True)
        - xp.sum(nyreps * (nyreps-1) * (2*nyreps + 5), axis=-1, keepdims=True))
    # Find the confidence interval indices in `slopes`
    sigma = xp.sqrt(xp.maximum(sigsq, xp.zeros_like(sigsq)))
    Ru = xp.minimum(xp.astype(xp.round((nt - z*sigma)/2.), xp.int64),
                    xp.astype(nt, xp.int64)-1)
    Rl = xp.maximum(xp.astype(xp.round((nt + z*sigma)/2.), xp.int64) - 1,
                    xp.asarray(0, dtype=xp.int64))
    R = xp.concat((xpx.atleast_nd(Rl, ndim=1), xpx.atleast_nd(Ru, ndim=1)), axis=-1)
    slopes = xp.sort(slopes, axis=-1)
    delta = xp.take_along_axis(slopes, R, axis=-1)
    i_nan = xp.broadcast_to(sigsq < 0, delta.shape)
    delta = xpx.at(delta)[i_nan].set(xp.nan)

    slope = medslope[()] if medslope.ndim == 0 else medslope
    intercept = medinter[()] if medinter.ndim == 0 else medinter
    low_slope = delta[..., 0]
    high_slope = delta[..., 1]
    low_slope = low_slope[()] if low_slope.ndim == 0 else low_slope
    high_slope = high_slope[()] if high_slope.ndim == 0 else high_slope
    return TheilslopesResult(slope=slope, intercept=intercept,
                             low_slope=low_slope, high_slope=high_slope)
