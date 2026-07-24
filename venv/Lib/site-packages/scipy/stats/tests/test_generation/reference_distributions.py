import numpy as np
import mpmath
from mpmath import mp


class ReferenceDistribution:
    """Minimalist distribution infrastructure for generating reference data.

    The purpose is to generate reference values for unit tests of SciPy
    distribution accuracy and robustness.

    Handles array input with standard broadcasting rules, and method
    implementations are easily compared against their mathematical definitions.
    No attempt is made to handle edge cases or be fast, and arbitrary precision
    arithmetic is trusted for accuracy rather than making the method
    implementations "smart".

    Notes
    -----

    In this infrastructure, distributions families are classes, and
    fully-specified distributions (i.e. with definite values of all family
    parameters) are instances of these classes. Typically, the public methods
    accept as input only the argument at which the at which the function is to
    be evaluated. Unlike SciPy distributions, they never accept values of
    distribution family shape, location, or scale parameters. A few
    other parameters are noteworthy:

    - All methods accept `dtype` to control the output data type. The default
      is `np.float64`, but `object` or `mp.mpf` may be
      specified to output the full `mpf`.
    - `ppf`/`isf` accept a `guess` because they use a scalar rootfinder
      to invert the `cdf`/`sf`. This is passed directly into the `x0` method
      of `mpmath.findroot`; see its documentation for details.
    - moment accepts `order`, an integer that specifies the order of the (raw)
      moment, and `center`, which is the value about which the moment is
      taken. The default is to calculate the mean and use it to calculate
      central moments; passing ``0`` results in a noncentral moment. For
      efficiency, the mean can be passed explicitly if it is already known.

    Follow the example of SkewNormal to generate new reference distributions,
    overriding only `__init__` and `_pdf`*. Use the reference distributions to
    generate reference values for unit tests of SciPy distribution method
    precision and robustness (e.g. for extreme arguments). If the a SciPy
    methods implementation is independent and yet the output matches reference
    values generated with this infrastructure, it is unlikely that the SciPy
    and reference values are both inaccurate.

    * If the SciPy output *doesn't* match and the cause appears to be
    inaccuracy of the reference values (e.g. due to numerical issues that
    mpmath's arbitrary precision arithmetic doesn't handle), then it may be
    appropriate to override a method of the reference distribution rather than
    relying on the generic implementation. Otherwise, hesitate to override
    methods: the generic implementations are mathematically correct and easy
    to verify, whereas an override introduces many possibilities of mistakes,
    requires more time to write, and requires more time to review.

    In general, do not create custom unit tests to ensure that
    SciPy distribution methods are *correct* (in the sense of being consistent
    with the rest of the distribution methods); generic tests take care of
    that.
    """

    def __init__(self, **kwargs):
        try:
            if mpmath.dps is not None:
                message = ("`mpmath.dps` has been assigned. This is not "
                           "intended usage; instead, assign the desired "
                           "precision to `mpmath.mp.dps` (e.g. `from mpmath "
                           "as mp; mp.dps = 50.")
                raise RuntimeError(message)
        except AttributeError:
            mpmath.dps = None

        if mp.dps <= 15:
            message = ("`mpmath.mp.dps <= 15`. Set a higher precision (e.g."
                       "`50`) to use this distribution.")
            raise RuntimeError(message)

        self._params = {key:self._make_mpf_array(val)
                        for key, val in kwargs.items()}

    def _make_mpf_array(self, x):
        shape = np.shape(x)
        x = np.asarray(x, dtype=np.float64).ravel()
        return np.asarray([mp.mpf(xi) for xi in x]).reshape(shape)[()]

    def _pdf(self, x):
        raise NotImplementedError("_pdf must be overridden.")

    def _cdf(self, x, **kwargs):
        if ((self._cdf.__func__ is ReferenceDistribution._cdf)
                and (self._sf.__func__ is not ReferenceDistribution._sf)):
            return mp.one - self._sf(x, **kwargs)

        a, b = self._support(**kwargs)
        res = mp.quad(lambda x: self._pdf(x, **kwargs), (a, x))
        res = res if res < 0.5 else mp.one - self._sf(x, **kwargs)
        return res

    def _sf(self, x, **kwargs):
        if ((self._sf.__func__ is ReferenceDistribution._sf)
                and (self._cdf.__func__ is not ReferenceDistribution._cdf)):
            return mp.one - self._cdf(x, **kwargs)

        a, b = self._support(**kwargs)
        res = mp.quad(lambda x: self._pdf(x, **kwargs), (x, b))
        res = res if res < 0.5 else mp.one - self._cdf(x, **kwargs)
        return res

    def _ppf(self, p, guess=0, **kwargs):
        if ((self._ppf.__func__ is ReferenceDistribution._ppf)
                and (self._isf.__func__ is not ReferenceDistribution._isf)):
            return self._isf(mp.one - p, guess, **kwargs)
        def f(x):
            return self._cdf(x, **kwargs) - p
        return mp.findroot(f, guess)

    def _isf(self, p, guess=0, **kwargs):
        if ((self._isf.__func__ is ReferenceDistribution._isf)
                and (self._ppf.__func__ is not ReferenceDistribution._ppf)):
            return self._ppf(mp.one - p, guess, **kwargs)
        def f(x):
            return self._sf(x, **kwargs) - p
        return mp.findroot(f, guess)

    def _logpdf(self, x, **kwargs):
        return mp.log(self._pdf(x, **kwargs))

    def _logcdf(self, x, **kwargs):
        return mp.log(self._cdf(x, **kwargs))

    def _logsf(self, x, **kwargs):
        return mp.log(self._sf(x, **kwargs))

    def _support(self, **kwargs):
        return -mp.inf, mp.inf

    def _entropy(self, **kwargs):
        def integrand(x):
            logpdf = self._logpdf(x, **kwargs)
            pdf = mp.exp(logpdf)
            return  -pdf*logpdf

        a, b = self._support(**kwargs)
        return mp.quad(integrand, (a, b))

    def _mean(self, **kwargs):
        return self._moment(order=1, center=0, **kwargs)

    def _var(self, **kwargs):
        mu = self._mean(**kwargs)
        return self._moment(order=2, center=mu, **kwargs)

    def _skew(self, **kwargs):
        mu = self._mean(**kwargs)
        u2 = self._moment(order=2, center=mu, **kwargs)
        sigma = mp.sqrt(u2)
        u3 = self._moment(order=3, center=mu, **kwargs)
        return u3 / sigma**3

    def _kurtosis(self, **kwargs):
        mu = self._mean(**kwargs)
        u2 = self._moment(order=2, center=mu, **kwargs)
        u4 = self._moment(order=4, center=mu, **kwargs)
        return u4 / u2**2 - 3

    def _moment(self, order, center, **kwargs):
        def integrand(x):
            return self._pdf(x, **kwargs)*(x - center)**order

        if center is None:
            center = self._mean(**kwargs)

        a, b = self._support(**kwargs)
        return mp.quad(integrand, (a, b))

    def pdf(self, x, dtype=np.float64):
        fun = np.vectorize(self._pdf)
        x = self._make_mpf_array(x)
        res = fun(x, **self._params)
        return np.asarray(res, dtype=dtype)[()]

    def cdf(self, x, dtype=np.float64):
        fun = np.vectorize(self._cdf)
        x = self._make_mpf_array(x)
        res = fun(x, **self._params)
        return np.asarray(res, dtype=dtype)[()]

    def sf(self, x, dtype=np.float64):
        fun = np.vectorize(self._sf)
        x = self._make_mpf_array(x)
        res = fun(x, **self._params)
        return np.asarray(res, dtype=dtype)[()]

    def ppf(self, x, guess=0, dtype=np.float64):
        fun = np.vectorize(self._ppf, excluded={1})  # don't vectorize guess
        x = self._make_mpf_array(x)
        res = fun(x, guess, **self._params)
        return np.asarray(res, dtype=dtype)[()]

    def isf(self, x, guess=0, dtype=np.float64):
        fun = np.vectorize(self._isf, excluded={1})  # don't vectorize guess
        x = self._make_mpf_array(x)
        res = fun(x, guess, **self._params)
        return np.asarray(res, dtype=dtype)[()]

    def logpdf(self, x, dtype=np.float64):
        fun = np.vectorize(self._logpdf)
        x = self._make_mpf_array(x)
        res = fun(x, **self._params)
        return np.asarray(res, dtype=dtype)[()]

    def logcdf(self, x, dtype=np.float64):
        fun = np.vectorize(self._logcdf)
        x = self._make_mpf_array(x)
        res = fun(x, **self._params)
        return np.asarray(res, dtype=dtype)[()]

    def logsf(self, x, dtype=np.float64):
        fun = np.vectorize(self._logsf)
        x = self._make_mpf_array(x)
        res = fun(x, **self._params)
        return np.asarray(res, dtype=dtype)[()]

    def support(self, dtype=np.float64):
        fun = np.vectorize(self._support)
        res = fun(**self._params)
        return np.asarray(res, dtype=dtype)[()]

    def entropy(self, dtype=np.float64):
        fun = np.vectorize(self._entropy)
        res = fun(**self._params)
        return np.asarray(res, dtype=dtype)[()]

    def mean(self, dtype=np.float64):
        fun = np.vectorize(self._mean)
        res = fun(**self._params)
        return np.asarray(res, dtype=dtype)[()]

    def var(self, dtype=np.float64):
        fun = np.vectorize(self._var)
        res = fun(**self._params)
        return np.asarray(res, dtype=dtype)[()]

    def skew(self, dtype=np.float64):
        fun = np.vectorize(self._skew)
        res = fun(**self._params)
        return np.asarray(res, dtype=dtype)[()]

    def kurtosis(self, dtype=np.float64):
        fun = np.vectorize(self._kurtosis)
        res = fun(**self._params)
        return np.asarray(res, dtype=dtype)[()]

    def moment(self, order, center=None, dtype=np.float64):
        fun = np.vectorize(self._moment)
        order = self._make_mpf_array(order)
        res = fun(order, **self._params)
        return np.asarray(res, dtype=dtype)[()]


class SkewNormal(ReferenceDistribution):
    """Reference implementation of the SkewNormal distribution.

    Follow the example here to generate new reference distributions.
    Use the reference distributions to generate reference values of
    distributions functions. For now, copy-paste the output into unit
    tests. Full code to generate reference values does not need to be
    included as a comment in the test; just refer to the reference
    distribution used and the settings (e.g. mp.dps=50).
    """

    def __init__(self, *, a):
        # Overriding __init__ is not necessary, but it allows IDEs to hint at
        # shape parameters. All parameters are keyword only to avoid the
        # ambiguity inherent in positional arguments. The infrastructure does
        # not take care of location and scale; nonetheless, assume standard
        # location and scale. Typically, there is no need to test the SciPy
        # distribution infrastructure's treatment of location and scale
        # separately for a specific distribution.
        super().__init__(a=a)

    def _support(self, a):
        # Override _support if the support of the distribution is a subset of
        # the real line
        return -mp.inf, mp.inf

    def _pdf(self, x, a):
        # Write PDFs following a scholarly reference as closely as possible.
        # Trust mpmath for the accuracy, and don't worry about speed. What's
        # important is the ease of verifying the PDF against the reference. If
        # the result is inaccurate, it will not match SciPy's output (whether
        # SciPy is accurate or not). If this occurs, try increasing dps before
        # implementing a numerically favorable (but presumably more complex)
        # implementation.
        return 2 * mp.npdf(x) * mp.ncdf(a * x)

    # Avoid overriding other methods unless the generic implementation is
    # believed to be inaccurate (e.g. due to numerical difficulties) or it is
    # too slow. Why? Less code to write, less code to review, and a guarantee
    # that there is no *mistake* in the implementation (e.g. wrong formula).


class BetaPrime(ReferenceDistribution):

    def __init__(self, *, a, b):
        super().__init__(a=a, b=b)

    def _support(self, **kwargs):
        return mp.zero, mp.inf

    def _logpdf(self, x, a, b):
        return (a - mp.one)*mp.log(x) - (a + b)*mp.log1p(x) - mp.log(mp.beta(a, b))

    def _pdf(self, x, a, b):
        return mp.exp(self._logpdf(x=x, a=a, b=b))

    def _sf(self, x, a, b):
        return 1.0 - mp.betainc(a, b, 0, x/(1+x), regularized=True)


class Burr(ReferenceDistribution):

    def __init__(self, *, c, d):
        super().__init__(c=c, d=d)

    def _support(self, c, d):
        return 0, mp.inf

    def _pdf(self, x, c, d):
        return c * d * x ** (-c - 1) * (1 + x ** (-c)) ** (-d - 1)

    def _ppf(self, p, guess, c, d):
        return (p**(-1.0/d) - 1)**(-1.0/c)


class LogLaplace(ReferenceDistribution):

    def __init__(self, *, c):
        super().__init__(c=c)

    def _support(self, c):
        return 0, mp.inf

    def _pdf(self, x, c):
        if x < mp.one:
            return c / 2 * x**(c - mp.one)
        else:
            return c / 2 * x**(-c - mp.one)

    def _ppf(self, q, guess, c):
        if q < 0.5:
            return (2.0 * q)**(mp.one / c)
        else:
            return (2 * (mp.one - q))**(-mp.one / c)


class LogNormal(ReferenceDistribution):

    def __init__(self, *, s):
        super().__init__(s=s)

    def _support(self, s):
        return 0, mp.inf

    def _pdf(self, x, s):
        return (
            mp.one / (s * x * mp.sqrt(2 * mp.pi))
            * mp.exp(-mp.one / 2 * (mp.log(x) / s)**2)
        )

    def _cdf(self, x, s):
        return mp.ncdf(mp.log(x) / s)


class Normal(ReferenceDistribution):

    def _pdf(self, x):
        return mp.npdf(x)


class NormInvGauss(ReferenceDistribution):

    def __init__(self, *, alpha, beta):
        super().__init__(alpha=alpha, beta=beta)

    def _pdf(self, x, alpha, beta):
        # Implemented as described in https://www.jstor.org/stable/4616433
        # Equations 2.1 - 2.3
        q = mp.sqrt(1 + x**2)
        a = mp.pi**-1 * alpha * mp.exp(mp.sqrt(alpha**2 - beta**2))
        return a * q**-1 * mp.besselk(1, alpha*q) * mp.exp(beta*x)


class Pearson3(ReferenceDistribution):
    def __init__(self, *, skew):
        super().__init__(skew=skew)

    def _pdf(self, x, skew):
        b = 2 / skew
        a = b**2
        c = -b
        res = abs(b)/mp.gamma(a) * (b*(x-c))**(a-1) * mp.exp(-b*(x-c))
        return res if abs(res.real) == res else 0


class StudentT(ReferenceDistribution):

    def __init(self, *, df):
        super().__init__(df=df)

    def _pdf(self, x, df):
        return (mp.gamma((df + mp.one)/2)/(mp.sqrt(df * mp.pi) * mp.gamma(df/2))
                * (mp.one + x*x/df)**(-(df + mp.one)/2))


class TruncExpon(ReferenceDistribution):

    def __init__(self, *, b):
        super().__init__(b=b)

    def _support(self, b):
        return 0, b

    def _pdf(self, x, b):
        return -mp.exp(-x)/mp.expm1(-b)

    def _sf(self, x, b):
        return (mp.exp(-b) - mp.exp(-x))/mp.expm1(-b)
