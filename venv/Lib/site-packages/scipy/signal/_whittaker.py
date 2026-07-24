import warnings
import numpy as np
from scipy._lib._util import _RichResult, _validate_int
from scipy.linalg.lapack import get_lapack_funcs
from scipy.optimize import minimize_scalar
from scipy.special import binom


def _solveh_banded(ab, b, calc_logdet=False):
    """
    Solve the equation ``a @ x = b`` for ``x``,  where ``a`` is the 
    Hermitian positive-definite banded matrix defined by `ab`.

    Same as scipy.linalg.solveh_banded(lower=True, check_finite=False), but:
    - also returns the log of the determinant and info
    - no error is raised if info > 0
    - no input validation
    - only real values, no complex
    - only `lower = True` code path
    - always overwrite_XX = False
    - b only a 1-dim array

    Parameters
    ----------
    ab : (``u`` + 1, M) array_like
        Banded matrix
    b : (M,) array_like
        Right-hand side

    Returns
    -------
    x : (M,) ndarray
        The solution to the system ``a x = b``. Shape of return matches shape of `b`.
    logdet : float
        Logarithm of the determinant of `ab`. Returns 0 if ``calc_logdet=False``.
    info : int
    """
    a1 = ab
    b1 = b
    overwrite_b = False
    overwrite_ab = False
    logdet = 0.0

    if a1.shape[0] == 2:
        method = "ptsv"
        ptsv = get_lapack_funcs(method, (a1, b1), ilp64="preferred")
        # We assume lower=True and real arrays
        d = a1[0, :]
        e = a1[1, :-1]
        # ptsv uses LDL', returnes d=diag(D), du=diag(L, -1)
        d, du, x, info = ptsv(d, e, b1, overwrite_ab, overwrite_ab, overwrite_b)
        if calc_logdet and info == 0:
            logdet = np.log(d).sum()
    else:
        method = "pbsv"
        pbsv = get_lapack_funcs(method, (a1, b1), ilp64="preferred")
        # pbsv uses Cholesky LL', returns c=L in ab-storage format
        c, x, info = pbsv(a1, b1, lower=True, overwrite_ab=overwrite_ab,
                          overwrite_b=overwrite_b)
        if calc_logdet and info == 0:
            logdet = 2 * np.log(c[0, :]).sum()
    if info < 0:
        raise ValueError(f"illegal value in {-info}th argument of internal {method}")
    return x, logdet, info


def whittaker_henderson(signal, *, lamb="reml", order=2, weights=None):
    r"""
    Whittaker-Henderson (WH) smoothing/graduation of a discrete signal.

    This implements WH smoothing with a difference penalty of the specified ``order``
    and penalty strength ``lamb``, see [1]_, [2]_ and [3]_. WH can be seen as a
    penalized B-Spline (P-Spline) of degree zero for equidistant knots at the signal
    positions.

    In econometrics, the WH graduation of order 2 is referred to as the
    Hodrick-Prescott filter [4]_.
    
    Parameters
    ----------
    signal : array_like
        A one-dimensional array of length ``order + 1`` representing equidistant data
        points of a signal, e.g. a time series with constant time lag.
    
    lamb : str or float, optional
        Smoothing or penalty parameter, default is ``"reml"`` which minimizes the
        restricted maximum likelihood (REML) criterion to find the parameter ``lamb``.
        If a number is passed, it must be non-negative and it is used directly.

    order : int, default: 2
        The order of the difference penalty, must be at least 1.

    weights : array_like or None, optional
        A one-dimensional array of case weights with the same length as `signal`.
        ``None`` is equivalent to an array of all ones, ``np.ones_like(signal)``.

    Returns
    -------
     res : _RichResult
        An object similar to an instance of `scipy.optimize.OptimizeResult` with the
        following attributes:

        x : ndarray
            The WH smoothed signal.
        lamb : float
            The penalty parameter. If the input ``lamb`` is a number, this value is
            returned. If the input was ``"reml"``, the penalty parameter that minimized
            the REML criterion is returned.

    Notes
    -----
    For the signal :math:`y = (y_1, y_2, \ldots, y_n)` and weights
    :math:`w = (w_1, w_2, \ldots, w_n)`, WH of order :math:`p` with smoothing or
    penalty parameter :math:`\lambda` solves the following optimization problem for
    :math:`x_i`:

    .. math::

        \operatorname{argmin}_{x_i} \sum_i^n w_i (y_i - x_i)^2
        + \lambda \sum_i^{n-p} (\Delta^p x_i)^2 \,,

    where :math:`\Delta^p` is the forward difference operator of order :math:`p`,
    :math:`\Delta x_i = x_{i+1} - x_i` and
    :math:`\Delta^2 x_i = \Delta(\Delta x_i) = x_{i+2} - 2x_{i+1} + x_i`.
    For every input value :math:`y_i`, it generates a smoothed value :math:`x_i`.

    One of the nice properties of WH is that it automatically performs inter- and
    extrapolation of data. Interpolation means filling in values when signal data is
    only available on the left and on the right. Extrapolation means filling in values
    after the observed signal ends.
    Set ``weights = 0`` for regions where you want to extra- or interpolate. The values
    of `signal` don't matter if ``weights = 0``.
    WH interpolates a polynomial of ``degree = 2 * order - 1`` and it extrapolates a
    polynomial of ``degree = order - 1``. For ``order = 2``, this means cubic
    interpolation and linear extrapolation.

    References
    ----------
    .. [1] Whittaker-Henderson smoothing,
           https://en.wikipedia.org/wiki/Whittaker%E2%80%93Henderson_smoothing
    .. [2] Eilers, P.H.C. (2003).
           "A perfect smoother". Analytical Chem. 75, 3631-3636.
           :doi:`10.1021/AC034173T`
    .. [3] Weinert, Howard L. (2007).
           "Efficient computation for Whittaker-Henderson smoothing".
           Computational Statistics and Data Analysis 52:959-74.
           :doi:`10.1016/j.csda.2006.11.038`
    .. [4] Hodrick, R. J., and Prescott, E. C. (1997).
           "Postwar U.S. Business Cycles: An Empirical Investigation".
           :doi:`10.2307/2953682`

    Examples
    --------
    We use data from https://data.giss.nasa.gov/gistemp/ for global temperature
    anomalies, i.e., deviations from the corresponding 1951-1980 means (Combined
    Land-Surface Air and Sea-Surface Water Temperature Anomalies, Land-Ocean Temperature
    Index, L-OTI).
    We use weights to indicate where to inter- and extrapolate the missing data.

    >>> from pathlib import Path
    >>> import numpy as np
    >>> import scipy
    >>> from scipy.signal import whittaker_henderson
    >>> # For the most recent data, use
    ... # fname="https://data.giss.nasa.gov/gistemp/tabledata_v4/GLB.Ts+dSST.csv"
    ... # Here, we instead use a copy made in 2026.
    ... fname = Path(scipy.signal.__file__).parent / "tests/data/GLB.Ts+dSST.csv"
    >>> data = np.genfromtxt(
    ...     fname=fname,
    ...     delimiter=",", skip_header=2, missing_values="***"
    ... )
    >>> year = data[:, 0]
    >>> temperature = data[:, 1:13].ravel()  # monthly temperature anomalies
    >>> w = np.ones_like(temperature)
    >>> # We might have some nan values.
    ... np.sum(np.isnan(temperature))
    np.int64(10)
    >>> w[np.isnan(temperature)] = 0
    >>> res = whittaker_henderson(temperature, weights=w)
    >>> temperature[:5]
    array([-0.19, -0.25, -0.09, -0.16, -0.1])
    >>> res.x[:5]
    array([-0.18244619, -0.17823282, -0.17409373, -0.17080896, -0.16833158])

    Let us plot measurements and Whittaker-Henderson smoothing.

    >>> import matplotlib.pyplot as plt
    >>> x = year[0] + np.arange(len(temperature)) / 12
    >>> plt.plot(x, temperature, label="measurement")
    >>> plt.plot(x, res.x, label="WH smooth")
    >>> # Above, we set w = 0 for nan values of temperature. WH automatically
    ... # inter- or extrapolates for all data points with w = 0.
    ... plt.plot(x[w==0], res.x[w==0], color="red", label="inter-/extrapolation")
    >>> plt.xlabel("year")
    >>> plt.ylabel("temperature deviation [°C]")
    >>> plt.title("Global Temperature Anomalies (ref. 1951-1980)")
    >>> plt.legend()
    >>> plt.show()

    We can see that extrapolation has occurred at the right end of the signal, meaning
    that NaNs existed in the data for the most recent dates, in particular months of
    2026 that have not yet happened at the time of the data download in March 2026.

    """
    order = _validate_int(order, name="order", minimum=1)

    signal = np.asarray(signal)
    if signal.ndim != 1:
        msg = f"Input array signal must be of shape (n,); got {signal.shape}"
        raise ValueError(msg)

    n = signal.shape[0]
    if n < order + 1:
        msg = f"Input array signal must be at least of shape ({order + 1},); got {n}."
        raise ValueError(msg)

    if weights is not None:
        weights = np.asarray(weights)
        if weights.shape != signal.shape:
            msg = "Input array weights must have the same shape as the signal array."
            raise ValueError(msg)
        
    if weights is None:
        if not np.isfinite(signal).all():
            raise ValueError("Input array signal must be finite.")
    else:
        if not np.isfinite(weights).all():
            raise ValueError("Input array weights must be finite.")
        if (mask := ~np.isfinite(signal)).any():
            # Only weights * y matter in the end and weights must be 0 for all
            # non-finite elements of signal.
            if not (weights[mask] == 0).all():
                raise ValueError("Input array weights must be zero for all non-finite "
                                 "elements of signal.")
            signal = np.nan_to_num(signal)


    msg = f"Parameter lamb must be string 'reml' or a non-negative float; got {lamb=}."
    if isinstance(lamb, str):
        if lamb != "reml":
            raise ValueError(msg)
        def criterion(loglamb):
            return -_reml(lamb=np.exp(loglamb), y=signal, order=order, weights=weights)
        opt = minimize_scalar(criterion, bracket=[-10, 10])
        lamb = np.exp(opt.x)

    if lamb < 0:
        raise ValueError(msg)
    elif lamb == 0.0:
        x = np.asarray(signal).copy()
    else:
        x, _ = _solve_WH_banded(signal, lamb=lamb, order=order, weights=weights)
    return _RichResult(x=x, lamb=lamb)


def _polynomial_fit(y, lamb, order=2, weights=None, calc_logdet=False):
    """Polynomial fit equivalent to WH for lamb -> infinity."""
    n = len(y)
    x_range = np.arange(n)
    poly = np.polynomial.Polynomial.fit(x=x_range, y=y, deg=order - 1, w=weights)
    if calc_logdet:
        # For large lambda, log|W + lambda D'D| ~ log|lambda D'D|
        # (with determinant understood as product of non-zero eigenvalues). 
        logdet_DtD = _logdet_difference_matrix(order=order, n=n)
        logdet = (n - order) * np.log(lamb) + logdet_DtD
    else:
        logdet = 0.0
    return poly(x_range), logdet


def _solve_WH_banded(y, lamb, order=2, weights=None, calc_logdet=False, warn_user=True):
    """
    Solve the WH optimization problem via the normal equations.
    
    A @ x = y
    A = I + lamb * P = I + lamb * D' @ D
    D = difference matrix of order=`order` 

    With weights W = diag(weights):
    A = W + lamb * P
    A @ x = W @ y

    Returns
    -------
    x : ndarray
        The solution.
    logdet : float
        Logarithm of the determinant of matrix A. Returns 0 if ``calc_logdet=False``.
    """
    n = y.shape[0]  # n >= p + 1 was already checked
    p = order  # order of difference penalty
    # Construct penalty matrix P = D'D of shape (n-p, n) as if n = 2p+1 (to save
    # memory).
    if n < 2*p + 1:
        D = np.diff(np.eye(n), n=p, axis=0)  # shape (n-p, n)
    else:
        D = np.diff(np.eye(2*p + 1), n=p, axis=0)  # shape (p+1, 2p+1)
    P_raw = D.T @ D  # shape (2p+1, 2p+1) if n >= 2p+1 else (n, n)

    # Because our matrix A = np.eye(n, dtype=np.float64) + lamb * (D.T @ D) is
    # symmetric and banded with u = l = p, we construct it in the lower "ab"-format
    # for use in solveh_banded, i.e. each row in ab is a subdiagonal of A:
    #   ab[0, :]   = np.diagonal(A, 0)
    #   ab[1, :-1] = np.diagonal(A, 1)
    #   ab[2, :-2] = np.diagonal(A, 2)
    #   ..
    ab = np.zeros((p + 1, min(2*p + 1, n)))
    for i in range(p + 1):
        ab[i, :ab.shape[1] - i] = np.diagonal(P_raw, i)
    ab *= lamb
    if n > 2*p + 1:
        ab = np.concat(
            [
                ab[:, :p+1],
                np.repeat(ab[:, p:p+1], n - (2*p+1), axis=1),
                ab[:, -p:],
            ],
            axis=1,
        )

    if weights is None:
        # Check if lambda is so large that A = I + lambda D'D = lambda D'D. We even add
        # a factor of 8, i.e. A should have at least 4 bits from I (not only D'D).
        # Note that the minimal diagonal element of D'D is always 1.
        if lamb * np.finfo(np.float64).eps > 8:
            # If lambda approaches infinity, WH approaches a polynomial fit.
            x, logdet = _polynomial_fit(
                y, lamb=lamb, order=order, calc_logdet=calc_logdet
            )
            info = 0
        else:
            ab[0, :] += 1.0  # This corresponds to np.eye(n).
            x, logdet, info = _solveh_banded(ab, y, calc_logdet=calc_logdet)
    else:
        if (ab[0, :] == ab[0, :] + weights * 8).all():
            # If lambda approaches infinity, WH approaches a polynomial fit.
            x, logdet = _polynomial_fit(
                y, lamb=lamb, order=order, weights=weights, calc_logdet=calc_logdet
            )
            info = 0
        else:
            ab[0, :] += weights
            x, logdet, info = _solveh_banded(ab, weights * y, calc_logdet=calc_logdet)

    if info > 0:
        # LinAlgError(f"{info}th leading minor not positive definite")
        # For very large values of lamb, we know that
        #   - the linear solver breaks down
        #   - the solution approaches a polynomial least squares fit of degree
        #     order - 1.
        # Note that for a certain large lamb, WH already reaches the polynomial least
        # squares fit almost exactly. For larger lamb, WH starts to deviate from the
        # polynomial (=worse solution due to numerical instability), until the solver
        # breaks down and reports info > 0.
        if warn_user:
            msg_weights = "" if weights is None else " or due to the weights"
            msg = (
                "The linear solver in Whittaker-Henderson smoothing detected a "
                "numerical instability. This is likely due to a very large value of "
                f"{lamb=}"
                + msg_weights + ". "
                "As Whittaker-Henderson approaches a polynomial of degree 'order - 1' "
                "for large lamb, this polynomial (via least squares) is returned."
            )
            warnings.warn(msg, UserWarning, stacklevel=2)
        x, logdet = _polynomial_fit(
            y, lamb, order=order, weights=weights, calc_logdet=calc_logdet
        )
    return x, logdet


def _logdet_difference_matrix(order, n):
    """Logarithm of the determinant of the difference matrix.

    If D is the difference matrix of order=p, then this computes `log det(D @ D.T)`
    which equals the sum of the log of non-zero eigenvalues of `D.T @ D`.
    """
    # product of eigenvalues =
    # prod(binom(n+i-1, 2i-1), i=1..p) / prod(binom(2i, i), i=1..p-1)
    # How to derive this formula? Well, ... some magic.
    p = order
    if order == 1:
        return np.log(n)
    logdet = 0.0
    for i in range(1, p+1):
        logdet += np.log(binom(n + i - 1, 2*i - 1) / binom(2*i, i)) 
    logdet += np.log(binom(2*p, p))
    return logdet


def _reml(lamb, y, order, weights=None):
    """Calculate the restricted maximum likelihood (REML).
    
    Parameters
    ----------
    lamb : penalty
    y : signal
    x : smoothed signal
    order : order of the difference penalty.
    weights : case weights

    Returns
    -------
    reml : REML criterion

    References
    ----------
    - Biessy https://arxiv.org/abs/2306.06932 (version 4)
    - Wood https://doi.org/10.1111/j.1467-9868.2010.00749.x
    """
    n = y.shape[0]
    x, logdet = _solve_WH_banded(
        y=y, lamb=lamb, order=order, weights=weights, calc_logdet=True, warn_user=False
    )
    logdet_DtD = _logdet_difference_matrix(order=order, n=n)
    residual = y - x
    # Eq. 12 of Biessy gives the REML criterion:
    # REML(lambda, sigma) = (log of restriced maximum likelihood)
    #     = -1/2 ((y - theta) W (y - theta) / sigma^2 + lambda theta D'D theta / sigma^2
    #             - log|lambda D'D| + log|(W + lambda D'D)| + (n - p) log(sigma^2)
    #             + const
    #            )
    # where the constant term "const" does not depend on lambda or sigma and p is the
    # order of the difference penalty.
    # Note that Biessy then does not mention to use the profiled REML criterion, i.e.,
    # analytically plug in the optimal sigma^2. This gives us
    #     sigma^2 = r2 / (n - p)
    #     r^2     = (y - theta) W (y - theta) + lambda theta D'D theta
    #     profiled REML(lambda) =
    #         -1/2 (
    #               (n-p) (1 + log(r^2 / (n-p)))
    #               -log|lambda D'D| + log|W + lambda D'D| + const
    #              )
    # This can be compared to Eq. 41 of Bates et al
    # https://doi.org/10.18637/jss.v067.i01.
    # An alternative derivation stems from a mixed model formulation of P-splines, see
    # Currie and Durban https://doi.org/10.1191/1471082x02st039ob or Boer
    # https://doi.org/10.1177/1471082X231178591. One then has 2 variance parameters
    # sigma^2 (from y) and tau^2 (from the random effect) leading to
    #     -1/2 ((y - theta) W (y - theta) / sigma^2 + theta D'D theta / tau^2 + ...
    # One then sets tau^2 = sigma^2 / lambda.
    if weights is None:
        r2 = residual @ residual
    else:
        r2 = residual @ (weights * residual)
    r2 += lamb * np.sum(np.diff(x, n=order)**2)  # + lambda theta D'D theta
    reml = (n - order) * (1 + np.log(r2 / (n - order)))
    reml -= (n - order) * np.log(lamb) + logdet_DtD  # -log|lambda D'D|
    reml += logdet  # +log|W + lambda D'D|
    reml *= -0.5
    return reml
