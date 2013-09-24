# pylint: disable-msg=E1103
# pylint: disable-msg=W0212

from __future__ import division

from pandas.compat import range
import numpy as np
import numpy.linalg as linalg


def rank(X, cond=1.0e-12):
    """
    Return the rank of a matrix X based on its generalized inverse,
    not the SVD.
    """
    X = np.asarray(X)
    if len(X.shape) == 2:
        import scipy.linalg as SL
        D = SL.svdvals(X)
        result = np.add.reduce(np.greater(D / D.max(), cond))
        return int(result.astype(np.int32))
    else:
        return int(not np.alltrue(np.equal(X, 0.)))


def solve(a, b):
    """Returns the solution of A X = B."""
    try:
        return linalg.solve(a, b)
    except linalg.LinAlgError:
        return np.dot(linalg.pinv(a), b)


def inv(a):
    """Returns the inverse of A."""
    try:
        return np.linalg.inv(a)
    except linalg.LinAlgError:
        return np.linalg.pinv(a)


def is_psd(m):
    eigvals = linalg.eigvals(m)
    return np.isreal(eigvals).all() and (eigvals >= 0).all()


def newey_west(m, max_lags, nobs, df, nw_overlap=False):
    """
    Compute Newey-West adjusted covariance matrix, taking into account
    specified number of leads / lags

    Parameters
    ----------
    m : (N x K)
    max_lags : int
    nobs : int
        Number of observations in model
    df : int
        Degrees of freedom in explanatory variables
    nw_overlap : boolean, default False
        Assume data is overlapping

    Returns
    -------
    ndarray (K x K)

    Reference
    ---------
    Newey, W. K. & West, K. D. (1987) A Simple, Positive
    Semi-definite, Heteroskedasticity and Autocorrelation Consistent
    Covariance Matrix, Econometrica, vol. 55(3), 703-708
    """
    Xeps = np.dot(m.T, m)
    for lag in range(1, max_lags + 1):
        auto_cov = np.dot(m[:-lag].T, m[lag:])
        weight = lag / (max_lags + 1)
        if nw_overlap:
            weight = 0
        bb = auto_cov + auto_cov.T
        dd = (1 - weight) * bb
        Xeps += dd

    Xeps *= nobs / (nobs - df)

    if nw_overlap and not is_psd(Xeps):
        new_max_lags = int(np.ceil(max_lags * 1.5))
#         print('nw_overlap is True and newey_west generated a non positive '
#               'semidefinite matrix, so using newey_west with max_lags of %d.'
#               % new_max_lags)
        return newey_west(m, new_max_lags, nobs, df)

    return Xeps


def calc_F(R, r, beta, var_beta, nobs, df):
    """
    Computes the standard F-test statistic for linear restriction
    hypothesis testing

    Parameters
    ----------
    R: ndarray (N x N)
        Restriction matrix
    r: ndarray (N x 1)
        Restriction vector
    beta: ndarray (N x 1)
        Estimated model coefficients
    var_beta: ndarray (N x N)
        Variance covariance matrix of regressors
    nobs: int
        Number of observations in model
    df: int
        Model degrees of freedom

    Returns
    -------
    F value, (q, df_resid), p value
    """
    from scipy.stats import f

    hyp = np.dot(R, beta.reshape(len(beta), 1)) - r
    RSR = np.dot(R, np.dot(var_beta, R.T))

    q = len(r)

    F = np.dot(hyp.T, np.dot(inv(RSR), hyp)).squeeze() / q

    p_value = 1 - f.cdf(F, q, nobs - df)

    return F, (q, nobs - df), p_value
