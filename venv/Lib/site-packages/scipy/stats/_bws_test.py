import numpy as np
from functools import partial

from scipy import stats
from scipy.stats._axis_nan_policy import _broadcast_arrays
from scipy._external import array_api_extra as xpx
from scipy._lib._array_api import (xp_capabilities, array_namespace, xp_ravel,
                                   xp_result_type)


def _bws_input_validation(x, y, alternative, axis, method):
    ''' Input validation and standardization for bws test'''
    xp = array_namespace(x, y)
    x, y = xp.asarray(x), xp.asarray(y)
    x, y = xpx.atleast_nd(x, ndim=1), xpx.atleast_nd(y, ndim=1)

    axis_none = axis is None
    ndim = max(x.ndim, y.ndim)
    if axis_none:
        x = xp_ravel(x, xp=xp)
        y = xp_ravel(y, xp=xp)
        axis = 0
    elif np.iterable(axis) or int(axis) != axis:
        message = "`axis` must be an integer or None."
        raise ValueError(message)
    elif (axis >= ndim) or (axis < -ndim):
        message = "`axis` is not compatible with the shapes of the inputs."
        raise ValueError(message)
    axis = int(axis)

    x, y = _broadcast_arrays((x, y), axis=axis, xp=xp)
    nx, ny = x.shape[-1], y.shape[-1]
    if nx < 2 or ny < 2:
        raise ValueError('`x` and `y` must contain at least two observations each.')

    x, y = xp.moveaxis(x, axis, -1), xp.moveaxis(y, axis, -1)
    z = stats.rankdata(xp.concat((x, y), axis=-1), axis=-1)
    z = xp.astype(z, xp_result_type(x, y, force_floating=True, xp=xp))
    x, y = z[..., :x.shape[-1]], z[..., x.shape[-1]:]
    x, y = xp.moveaxis(x, -1, axis), xp.moveaxis(y, -1, axis)

    alternatives = {'two-sided', 'less', 'greater'}
    alternative = alternative.lower()
    if alternative not in alternatives:
        raise ValueError(f'`alternative` must be one of {alternatives}.')

    method = stats.PermutationMethod() if method is None else method
    if not isinstance(method, stats.PermutationMethod):
        raise ValueError('`method` must be an instance of '
                         '`scipy.stats.PermutationMethod`')

    return x, y, alternative, axis, method, xp


def _bws_statistic(x, y, alternative, axis, xp):
    '''Compute the BWS test statistic for two independent samples'''
    # Public function currently does not accept `axis`, but `permutation_test`
    # uses `axis` to make vectorized call.

    Ri, Hj = xp.sort(x, axis=axis), xp.sort(y, axis=axis)
    n, m = Ri.shape[axis], Hj.shape[axis]
    i, j = xp.arange(1, n+1, dtype=Ri.dtype), xp.arange(1, m+1, dtype=Hj.dtype)

    Bx_num = Ri - (m + n)/n * i
    By_num = Hj - (m + n)/m * j

    if alternative == 'two-sided':
        Bx_num = xpx.at(Bx_num)[...].multiply(Bx_num)
        By_num = xpx.at(By_num)[...].multiply(By_num)
    else:
        Bx_num = xpx.at(Bx_num)[...].multiply(xp.abs(Bx_num))
        By_num = xpx.at(By_num)[...].multiply(xp.abs(By_num))

    Bx_den = i/(n+1) * (1 - i/(n+1)) * m*(m+n)/n
    By_den = j/(m+1) * (1 - j/(m+1)) * n*(m+n)/m

    Bx = 1/n * xp.sum(Bx_num/Bx_den, axis=axis)
    By = 1/m * xp.sum(By_num/By_den, axis=axis)

    B = (Bx + By) / 2 if alternative == 'two-sided' else (Bx - By) / 2

    return B


@xp_capabilities(skip_backends=[('dask.array', 'no rankdata')])
def bws_test(x, y, *, alternative="two-sided", axis=0, method=None):
    r'''Perform the Baumgartner-Weiss-Schindler test on two independent samples.

    The Baumgartner-Weiss-Schindler (BWS) test is a nonparametric test of
    the null hypothesis that the distribution underlying sample `x`
    is the same as the distribution underlying sample `y`. Unlike
    the Kolmogorov-Smirnov, Wilcoxon, and Cramer-Von Mises tests,
    the BWS test weights the integral by the variance of the difference
    in cumulative distribution functions (CDFs), emphasizing the tails of the
    distributions, which increases the power of the test in many applications.

    Parameters
    ----------
    x, y : array-like
        Arrays of samples. Shapes must be broadcastable except along `axis`.
    alternative : {'two-sided', 'less', 'greater'}, optional
        Defines the alternative hypothesis. Default is 'two-sided'.
        Let *F(u)* and *G(u)* be the cumulative distribution functions of the
        distributions underlying `x` and `y`, respectively. Then the following
        alternative hypotheses are available:

        * 'two-sided': the distributions are not equal, i.e. *F(u) ≠ G(u)* for
          at least one *u*.
        * 'less': the distribution underlying `x` is stochastically less than
          the distribution underlying `y`, i.e. *F(u) >= G(u)* for all *u*.
        * 'greater': the distribution underlying `x` is stochastically greater
          than the distribution underlying `y`, i.e. *F(u) <= G(u)* for all
          *u*.

        Under a more restrictive set of assumptions, the alternative hypotheses
        can be expressed in terms of the locations of the distributions;
        see [2] section 5.1.
    axis : int or None, default: 0
        Axis along which the quantiles are computed.
        ``None`` ravels both `x` and `y` before performing the calculation.
    method : PermutationMethod, optional
        Configures the method used to compute the p-value. The default is
        the default `PermutationMethod` object.

    Returns
    -------
    res : PermutationTestResult
        An object with attributes:

        statistic : float
            The observed test statistic of the data.
        pvalue : float
            The p-value for the given alternative.
        null_distribution : ndarray
            The values of the test statistic generated under the null hypothesis.

    See Also
    --------
    scipy.stats.wilcoxon, scipy.stats.mannwhitneyu, scipy.stats.ttest_ind

    Notes
    -----
    When ``alternative=='two-sided'``, the statistic is defined by the
    equations given in [1]_ Section 2. This statistic is not appropriate for
    one-sided alternatives; in that case, the statistic is the *negative* of
    that given by the equations in [1]_ Section 2. Consequently, when the
    distribution of the first sample is stochastically greater than that of the
    second sample, the statistic will tend to be positive.

    References
    ----------
    .. [1] Neuhäuser, M. (2005). Exact Tests Based on the
           Baumgartner-Weiss-Schindler Statistic: A Survey. Statistical Papers,
           46(1), 1-29.
    .. [2] Fay, M. P., & Proschan, M. A. (2010). Wilcoxon-Mann-Whitney or t-test?
           On assumptions for hypothesis tests and multiple interpretations of
           decision rules. Statistics surveys, 4, 1.

    Examples
    --------
    We follow the example of table 3 in [1]_: Fourteen children were divided
    randomly into two groups. Their ranks at performing a specific tests are
    as follows.

    >>> import numpy as np
    >>> x = [1, 2, 3, 4, 6, 7, 8]
    >>> y = [5, 9, 10, 11, 12, 13, 14]

    We use the BWS test to assess whether there is a statistically significant
    difference between the two groups.
    The null hypothesis is that there is no difference in the distributions of
    performance between the two groups. We decide that a significance level of
    1% is required to reject the null hypothesis in favor of the alternative
    that the distributions are different.
    Since the number of samples is very small, we can compare the observed test
    statistic against the *exact* distribution of the test statistic under the
    null hypothesis.

    >>> from scipy.stats import bws_test
    >>> res = bws_test(x, y)
    >>> print(res.statistic)
    5.132167152575315

    This agrees with :math:`B = 5.132` reported in [1]_. The *p*-value produced
    by `bws_test` also agrees with :math:`p = 0.0029` reported in [1]_.

    >>> print(res.pvalue)
    0.002913752913752914

    Because the p-value is below our threshold of 1%, we take this as evidence
    against the null hypothesis in favor of the alternative that there is a
    difference in performance between the two groups.
    '''

    x, y, alternative, axis, method, xp = _bws_input_validation(x, y, alternative,
                                                                axis, method)
    bws_statistic = partial(_bws_statistic, alternative=alternative, xp=xp)

    permutation_alternative = 'less' if alternative == 'less' else 'greater'
    res = stats.permutation_test((x, y), bws_statistic,
                                 alternative=permutation_alternative, axis=axis,
                                 **method._asdict())

    return res
