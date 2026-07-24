import math


def get_err(result, stats):
    """
    Computes an 'error measure' based on the interquartile range of the
    measurement results.

    #### Parameters
    **result** (`any`)
    : The measurement results. Currently unused.

    **stats** (`dict`)
    : A dictionary of statistics computed from the measurement results.
      It should contain the keys "q_25" and "q_75" representing the 25th and
      75th percentiles respectively.

    #### Returns
    **error** (`float`)
    : The error measure, defined as half the interquartile range
    (i.e., (Q3 - Q1) / 2).
    """
    a, b = stats["q_25"], stats["q_75"]
    return (b - a) / 2


def binom_pmf(n, k, p):
    """
    Computes the Probability Mass Function (PMF) for a binomial distribution.

    #### Parameters
    **n** (`int`)
    : The number of trials in the binomial distribution.

    **k** (`int`)
    : The number of successful trials.

    **p** (`float`)
    : The probability of success on each trial.

    #### Returns
    **pmf** (`float`)
    : The binomial PMF computed as (n choose k) * p**k * (1 - p)**(n - k).

    #### Notes
    Handles edge cases where p equals 0 or 1.
    """
    if not (0 <= k <= n):
        return 0
    if p == 0:
        return 1.0 * (k == 0)
    elif p == 1.0:
        return 1.0 * (k == n)

    logp = math.log(p)
    log1mp = math.log(1 - p)
    lg1pn = math.lgamma(1 + n)
    lg1pnmk = math.lgamma(1 + n - k)
    lg1pk = math.lgamma(1 + k)
    return math.exp(lg1pn - lg1pnmk - lg1pk + k * logp + (n - k) * log1mp)


def quantile(x, q):
    """
    Computes a quantile/percentile from a dataset.

    #### Parameters
    **x** (`list` of `float`)
    : The dataset for which the quantile is to be computed.

    **q** (`float`)
    : The quantile to compute. Must be in the range [0, 1].

    #### Returns
    **m** (`float`)
    : The computed quantile from the dataset.

    #### Raises
    **ValueError**
    : If the provided quantile q is not in the range [0, 1].

    #### Notes
    This function sorts the input data and calculates the quantile
    using a linear interpolation method if the desired quantile lies
    between two data points.
    """
    if not 0 <= q <= 1:
        raise ValueError("Invalid quantile")

    y = sorted(x)
    n = len(y)

    z = (n - 1) * q
    j = int(math.floor(z))
    z -= j

    return y[-1] if j == n - 1 else (1 - z) * y[j] + z * y[j + 1]


def quantile_ci(x, q, alpha_min=0.01):
    """
    Compute a quantile and a confidence interval for a given dataset.

    #### Parameters
    **x** (`list` of `float`)
    : The dataset from which the quantile and confidence interval are computed.

    **q** (`float`)
    : The quantile to compute. Must be in the range [0, 1].

    **alpha_min** (`float`, optional)
    : Limit for coverage. The result has coverage >= 1 - alpha_min. Defaults to 0.01.

    #### Returns
    **m** (`float`)
    : The computed quantile from the dataset.

    **ci** (`tuple` of `float`)
    : Confidence interval (a, b), of coverage >= alpha_min.

    #### Notes
    This function assumes independence but is otherwise nonparametric. It sorts
    the input data and calculates the quantile using a linear interpolation
    method if the desired quantile lies between two data points. The confidence
    interval is computed using a known property of the cumulative distribution
    function (CDF) of a binomial distribution. This method calculates the
    smallest range (y[r-1], y[s-1]) for which the coverage is at least
    alpha_min.
    """

    y = sorted(x)
    n = len(y)

    alpha_min = min(alpha_min, 1 - alpha_min)
    pa = alpha_min / 2
    pb = 1 - pa

    a = -math.inf
    b = math.inf

    # It's known that
    #
    # Pr[X_{(r)} < m < X_{(s)}] = Pr[r <= K <= s-1], K ~ Bin(n,p)
    #
    # where cdf(m) = p defines the quantile.
    #
    # Simplest median CI follows by picking r,s such that
    #
    #   F(r;n,q) <= alpha/2
    #   F(s;n,q) >= 1 - alpha/2
    #
    #   F(k;n,q) = sum(binom_pmf(n, j, q) for j in range(k))
    #
    # Then (y[r-1], y[s-1]) is a CI.
    # If no such r or s exists, replace by +-inf.

    F = 0
    for k, yp in enumerate(y):
        F += binom_pmf(n, k, q)
        # F = F(k+1;n,q)

        if F <= pa:
            a = yp

        if F >= pb:
            b = yp
            break

    m = quantile(y, q)
    return m, (a, b)


class LaplacePosterior:
    """
    Class to represent univariate Laplace posterior distribution.

    #### Description
    This class represents the univariate posterior distribution defined as
    `p(beta|y) = N [sum(|y_j - beta|)]**(-nu-1)` where N is the normalization factor.

    #### Parameters
    **y** (`list` of `float`)
    : A list of sample values from the distribution.

    **nu** (`float`, optional)
    : Degrees of freedom. Default is `len(y) - 1`.

    #### Attributes
    **mle** (`float`)
    : The maximum likelihood estimate for beta which is the median of y.

    #### Notes
    This is the posterior distribution in the Bayesian model assuming Laplace
    distributed noise, where `p(y|beta,sigma) = N exp(- sum_j (1/sigma) |y_j -
    beta|)`, `p(sigma) ~ 1/sigma`, and `nu = len(y) - 1`. The MLE for beta is
    `median(y)`. Applying the same approach to a Gaussian model results to
    `p(beta|y) = N T(t, m-1)`, `t = (beta - mean(y)) / (sstd(y) / sqrt(m))`
    where `T(t, nu)` is the Student t-distribution pdf, which gives the standard
    textbook formulas for the mean.
    """

    def __init__(self, y, nu=None):
        """
        Initializes an instance of the LaplacePosterior class.

        #### Parameters
        **y** (`list` of `float`):
        : The samples from the distribution.

        **nu** (`float`, optional):
        : The degrees of freedom. Default is `len(y) - 1`.

        #### Raises
        `ValueError`: If `y` is an empty list.

        #### Notes
        This constructor sorts the input data `y` and calculates the MLE
        (Maximum Likelihood Estimate).  It computes a scale factor, `_y_scale`,
        to prevent overflows when computing unnormalized CDF integrals.  The
        input data `y` is then shifted and scaled according to this computed
        scale.  The method also initializes a memoization dictionary `_cdf_memo`
        for the unnormalized CDF, and a placeholder `_cdf_norm` for the
        normalization constant of the CDF.
        """
        if len(y) == 0:
            raise ValueError("empty input")

        self.nu = len(y) - 1 if nu is None else nu
        # Sort input
        y = sorted(y)

        # Get location and scale so that data is centered at MLE, and
        # the unnormalized PDF at MLE has amplitude ~ 1/nu.
        #
        # Proper scaling of inputs is important to avoid overflows
        # when computing the unnormalized CDF integrals below.
        self.mle = quantile(y, 0.5)
        self._y_scale = sum(abs(yp - self.mle) for yp in y)
        self._y_scale *= self.nu ** (1 / (self.nu + 1))

        # Shift and scale
        if self._y_scale != 0:
            self.y = [(yp - self.mle) / self._y_scale for yp in y]
        else:
            self.y = [0 for _ in y]

        self._cdf_norm = None
        self._cdf_memo = {}

    def _cdf_unnorm(self, beta):
        """
        Computes the unnormalized cumulative distribution function (CDF).

        #### Parameters
        **beta** (`float`):
        : The upper limit of the integration for the CDF.

        #### Returns
        Returns the unnormalized CDF evaluated at `beta`.

        #### Notes
        The method computes the unnormalized CDF as:

            cdf_unnorm(b) = int_{-oo}^{b} 1/(sum_j |y - b'|)**(m+1) db'

        The method integrates piecewise, resolving the absolute values
        separately for each section. The results of these calculations
        are memoized to speed up subsequent computations.

        It also handles special cases, such as when `beta` is not a number
        (returns `beta` as is), or when `beta` is positive infinity
        (memoizes the integral value at the end of the list `y`).
        """
        if beta != beta:
            return beta

        k0 = next((k for k, y in enumerate(self.y) if y > beta), len(self.y))
        cdf = 0

        nu = self.nu

        # Save some work by memoizing intermediate results
        if k0 - 1 in self._cdf_memo:
            k_start = k0
            cdf = self._cdf_memo[k0 - 1]
        else:
            k_start = 0
            cdf = 0

        # Do the integral piecewise, resolving the absolute values
        for k in range(k_start, k0 + 1):
            c = 2 * k - len(self.y)
            y = sum(self.y[k:]) - sum(self.y[:k])

            a = -math.inf if k == 0 else self.y[k - 1]
            b = beta if k == k0 else self.y[k]
            if c == 0:
                term = (b - a) / y ** (nu + 1)
            else:
                term = 1 / (nu * c) * ((a * c + y) ** (-nu) - (b * c + y) ** (-nu))

            cdf += max(0, term)  # avoid rounding error

            if k != k0:
                self._cdf_memo[k] = cdf

        if beta == math.inf:
            self._cdf_memo[len(self.y)] = cdf

        return cdf

    def _ppf_unnorm(self, cdfx):
        """
        Computes the inverse function of `_cdf_unnorm`.

        #### Parameters
        **cdfx** (`float`):
        : The value for which to compute the inverse cumulative distribution
        function (CDF).

        #### Returns
        Returns the unnormalized quantile function evaluated at `cdfx`.

        #### Notes
        This method computes the inverse of `_cdf_unnorm`. It first finds the
        interval within which `cdfx` lies, then performs the inversion on this
        interval.

        Special cases are handled when the interval index `k` is 0 (the
        computation of `beta` involves a check for negative infinity), or when
        the calculated `c` is 0. The result `beta` is clipped at the upper bound
        of the interval, ensuring it does not exceed `self.y[k]`.
        """
        # Find interval
        for k in range(len(self.y) + 1):
            if cdfx <= self._cdf_memo[k]:
                break

        # Invert on interval
        c = 2 * k - len(self.y)
        y = sum(self.y[k:]) - sum(self.y[:k])
        nu = self.nu
        if k == 0:
            term = cdfx
        else:
            a = self.y[k - 1]
            term = cdfx - self._cdf_memo[k - 1]

        if k == 0:
            z = -nu * c * term
            beta = (z ** (-1 / nu) - y) / c if z > 0 else -math.inf
        elif c == 0:
            beta = a + term * y ** (nu + 1)
        else:
            z = (a * c + y) ** (-nu) - nu * c * term
            beta = (z ** (-1 / nu) - y) / c if z > 0 else math.inf
        if k < len(self.y):
            beta = min(beta, self.y[k])

        return beta

    def pdf(self, beta):
        """
        Computes the probability distribution function (PDF).

        #### Parameters
        **beta** (`float`)
        : The point at which to evaluate the PDF.

        #### Returns
        A `float` which is the probability density function evaluated at `beta`.

        #### Notes
        This function computes the PDF by exponentiating the result of
        `self.logpdf(beta)`.  The `logpdf` method should therefore be
        implemented in the class that uses this method.
        """
        return math.exp(self.logpdf(beta))

    def logpdf(self, beta):
        """
        Computes the logarithm of the probability distribution function (log-PDF).

        #### Parameters
        **beta** (`float`)
        : The point at which to evaluate the log-PDF.

        #### Returns
        A `float` which is the logarithm of the probability density function
        evaluated at `beta`.

        #### Notes
        This function computes the log-PDF by first checking if the scale of the
        distribution `_y_scale` is zero.  If so, it returns `math.inf` if `beta`
        equals the maximum likelihood estimate `mle`, otherwise it returns
        `-math.inf`.

        The `beta` value is then transformed by subtracting the maximum
        likelihood estimate `mle` and dividing by `_y_scale`.

        If the unnormalized cumulative distribution function `_cdf_norm` has not
        been computed yet, it is computed by calling `_cdf_unnorm(math.inf)`.

        The function then computes the sum of absolute differences between
        `beta` and all points in `y`, applies the log-PDF formula and returns
        the result.
        """
        if self._y_scale == 0:
            return math.inf if beta == self.mle else -math.inf

        beta = (beta - self.mle) / self._y_scale

        if self._cdf_norm is None:
            self._cdf_norm = self._cdf_unnorm(math.inf)

        ws = sum(abs(yp - beta) for yp in self.y)
        m = self.nu
        return (
            -(m + 1) * math.log(ws) - math.log(self._cdf_norm) - math.log(self._y_scale)
        )

    def cdf(self, beta):
        """
        Computes the cumulative distribution function (CDF).

        #### Parameters
        **beta** (`float`)
        : The point at which to evaluate the CDF.

        #### Returns
        A `float` which is the value of the cumulative distribution function
        evaluated at `beta`.

        #### Notes
        This function computes the CDF by first checking if the scale of the
        distribution `_y_scale` is zero.  If so, it returns 1 if `beta` is
        greater than the maximum likelihood estimate `mle`, and 0 otherwise.

        The `beta` value is then transformed by subtracting the maximum
        likelihood estimate `mle` and dividing by `_y_scale`.

        If the unnormalized cumulative distribution function `_cdf_norm` has not
        been computed yet, it is computed by calling `_cdf_unnorm(math.inf)`.

        The function then computes the unnormalized CDF at `beta` and normalizes
        it by dividing with `_cdf_norm`.
        """
        if self._y_scale == 0:
            return 1.0 * (beta > self.mle)

        beta = (beta - self.mle) / self._y_scale

        if self._cdf_norm is None:
            self._cdf_norm = self._cdf_unnorm(math.inf)
        return self._cdf_unnorm(beta) / self._cdf_norm

    def ppf(self, cdf):
        """
        Computes the percent point function (PPF), also known as the inverse
        cumulative distribution function.

        #### Parameters
        **cdf** (`float`)
        : The cumulative probability for which to compute the inverse CDF. It
        must be between 0 and 1 (inclusive).

        #### Returns
        A `float` which is the value of the percent point function evaluated at
        `cdf`.

        #### Notes
        This function computes the PPF by first checking if `cdf` is not between
        0 and 1. If it is not, it returns `math.nan`.

        If the scale of the distribution `_y_scale` is zero, it returns the
        maximum likelihood estimate `mle`.

        If the unnormalized cumulative distribution function `_cdf_norm` has not
        been computed yet, it is computed by calling `_cdf_unnorm(math.inf)`.

        The function then scales `cdf` by `_cdf_norm` (making sure it does not
        exceed `_cdf_norm`), computes the unnormalized PPF at this scaled value,
        and transforms it back to the original scale.
        """
        if cdf < 0 or cdf > 1.0:
            return math.nan

        if self._y_scale == 0:
            return self.mle

        if self._cdf_norm is None:
            self._cdf_norm = self._cdf_unnorm(math.inf)

        cdfx = min(cdf * self._cdf_norm, self._cdf_norm)
        beta = self._ppf_unnorm(cdfx)
        return beta * self._y_scale + self.mle


def compute_stats(samples, number):
    """
    Performs statistical analysis on the provided samples.

    #### Parameters
    **samples** (`list` of `float`)
    : A list of total times (in seconds) of benchmarks.

    **number** (`int`)
    : The number of times each benchmark was repeated.

    #### Returns
    **beta_hat** (`float`)
    : The estimated time per iteration.

    **stats** (`dict`)
    : A dictionary containing various statistical measures of the estimator. It
    includes:
        - **"ci_99_a"**: The lower bound of the 99% confidence interval.
        - **"ci_99_b"**: The upper bound of the 99% confidence interval.
        - **"q_25"**: The 25th percentile of the sample times.
        - **"q_75"**: The 75th percentile of the sample times.
        - **"repeat"**: The total number of samples.
        - **"number"**: The repeat number for each sample.

    #### Notes
    This function first checks if there are any samples. If there are none, it
    returns `None, None`.

    It then calculates the median and the 25th and 75th percentiles of the
    samples. If the nonparametric confidence interval estimation did not provide
    an estimate, it computes the posterior distribution for the location,
    assuming exponential noise. The Maximum Likelihood Estimate (MLE) is equal
    to the median. The function uses the confidence interval from that
    distribution to extend beyond the sample bounds if necessary.

    Finally, it produces the median as the result and a dictionary of the
    computed statistics.
    """

    if len(samples) < 1:
        return None, None

    Y = list(samples)

    # Median and quantiles
    y_50, ci_50 = quantile_ci(Y, 0.5, alpha_min=0.99)
    y_25 = quantile(Y, 0.25)
    y_75 = quantile(Y, 0.75)

    # If nonparametric CI estimation didn't give an estimate,
    # use the credible interval of a bayesian posterior distribution.
    a, b = ci_50
    if (math.isinf(a) or math.isinf(b)) and len(Y) > 1:
        # Compute posterior distribution for location, assuming
        # exponential noise. The MLE is equal to the median.
        c = LaplacePosterior(Y)

        # Use the CI from that distribution to extend beyond sample
        # bounds
        if math.isinf(a):
            a = min(c.ppf(0.01 / 2), min(Y))
        if math.isinf(b):
            b = max(c.ppf(1 - 0.01 / 2), max(Y))

        ci_50 = (a, b)

    # Produce results
    result = y_50

    stats = {
        "ci_99_a": ci_50[0],
        "ci_99_b": ci_50[1],
        "q_25": y_25,
        "q_75": y_75,
        "repeat": len(Y),
        "number": number,
    }

    return result, stats
