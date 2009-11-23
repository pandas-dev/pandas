# pylint: disable-msg=E1103
# pylint: disable-msg=W0212

from __future__ import division

from scipy import stats
import numpy as np
import numpy.linalg as linalg

from pandas.stats.common import FULL_SAMPLE, EXPANDING, ROLLING, TIME, ENTITY

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

def calc_xx_with_time_effects(x, y):
    """
    Returns X'X - (X'T) (T'T)^-1 (T'X)
    """
    # X'X
    xx = np.dot(x.values.T, x.values)
    xt = x.sum().values

    # X'X - (T'T)^-1 (T'X)
    count = y.count()
    selector = count > 0

    xt = xt[selector]
    count = count[selector]

    return xx - np.dot(xt.T / count, xt)

def is_psd(m):
    eigvals = linalg.eigvals(m)

    return np.isreal(eigvals).all() and (eigvals >= 0).all()

def newey_west(m, nw_lags, nobs, df, nw_overlap=False):
    """Returns the Newey West values."""
    Xeps = np.dot(m.T, m)
    for lag in xrange(1, nw_lags + 1):
        auto_cov = np.dot(m[:-lag].T, m[lag:])
        weight = lag / (nw_lags + 1)
        if nw_overlap:
            weight = 0
        bb = auto_cov + auto_cov.T
        dd = (1 - weight) * bb
        Xeps += dd

    Xeps *= nobs / (nobs - df)

    if nw_overlap and not is_psd(Xeps):
        new_nw_lags = int(np.ceil(nw_lags * 1.5))
        print ('nw_overlap is True and newey_west generated a non positive '
               'semidefinite matrix, so using newey_west with nw_lags of %d.'
               % new_nw_lags)
        return newey_west(m, new_nw_lags, nobs, df)

    return Xeps

def var_beta_panel(y, x, beta, xx, rmse, cluster_axis,
                   nw_lags, nobs, df, nw_overlap):

    from pandas.core.panel import LongPanel, group_agg

    xx_inv = inv(xx)

    if cluster_axis is None:
        if nw_lags is None:
            return xx_inv * (rmse ** 2)
        else:
            resid = y.values.squeeze() - np.dot(x.values, beta)
            m = (x.values.T * resid).T

            xeps = newey_west(m, nw_lags, nobs, df, nw_overlap)

            return np.dot(xx_inv, np.dot(xeps, xx_inv))
    else:
        Xb = np.dot(x.values, beta).reshape((len(x.values), 1))
        resid = LongPanel(y.values - Xb, ['resid'], y.index)

        if cluster_axis == 1:
            x = x.swapaxes()
            resid = resid.swapaxes()

        m = group_agg(x.values * resid.values, x.index._bounds,
                      lambda x: np.sum(x, axis=0))

        if nw_lags is None:
            nw_lags = 0

        xox = 0
        for i in range(len(x.major_axis)):
            xox += newey_west(m[i : i + 1], nw_lags,
                              nobs, df, nw_overlap)

        return np.dot(xx_inv, np.dot(xox, xx_inv))


def xx_time_effects(x, y):
    """
    Returns X'X - (X'T) (T'T)^-1 (T'X)
    """
    # X'X
    xx = np.dot(x.values.T, x.values)
    xt = x.sum().values

    count = y.count()
    selector = count > 0

    # X'X - (T'T)^-1 (T'X)
    xt = xt[selector]
    count = count[selector]

    return xx - np.dot(xt.T / count, xt)

def calc_var_beta(x, y, nw_lags, rmse, beta, nobs, df, window_type=FULL_SAMPLE,
                  window=None, nw_overlap=False):
    """Returns the covariance of beta.

    For a full-sample regression, this returns the covariance matrix of betas.
    For rolling/expanding regressions, this returns the variances of betas.
    """
    if window_type == FULL_SAMPLE:
        xx = np.dot(x.T, x)
    else:
        cum_xx = []
        cum_xx.append(np.dot(x[0 : 1].T, x[0 : 1]))

        for i in xrange(1, len(y)):
            cum_xx.append(cum_xx[i - 1] + np.dot(x[i : i + 1].T,
                      x[i : i + 1]))

    if window_type == FULL_SAMPLE:
        if nw_lags is None:
            return inv(xx) * (rmse ** 2)
        else:
            resid = y - np.dot(x, beta)
            m = (x.T * resid).T

            xeps = newey_west(m, nw_lags, nobs, df, nw_overlap)

            xx_inv = inv(xx)
            return np.dot(xx_inv, np.dot(xeps, xx_inv))
    else:
        results = []
        start = window - 1
        for i in xrange(start, len(y)):
            if nw_lags is None:
                temp_xx = cum_xx[i]
                if window_type == ROLLING and i >= window:
                    temp_xx = temp_xx - cum_xx[i - window]
                result = inv(temp_xx) * (rmse[i - start] ** 2)
            else:
                temp_xx = cum_xx[i]

                if window_type == EXPANDING:
                    begin = 0
                else:
                    begin = i - start
                    if i >= window:
                        temp_xx = temp_xx - cum_xx[i - window]

                section = slice(begin, i + 1)

                resid = y[section] - np.dot(x[section], beta[i - start])
                m = (x[section].T * resid).T

                window_nobs = i + 1 - begin
                window_df = df[i - start]

                xeps = newey_west(m, nw_lags, window_nobs,
                                  window_df, nw_overlap)

                xx_inv = inv(temp_xx)
                result = np.dot(xx_inv, np.dot(xeps, xx_inv))

            results.append(result)

        return results

def calc_F(R, r, beta, var_beta, nobs, df):
    hyp = np.dot(R, beta.reshape(len(beta), 1)) - r
    RSR = np.dot(R, np.dot(var_beta, R.T))

    q = len(r)

    F = np.dot(hyp.T, np.dot(inv(RSR), hyp)).squeeze() / q

    p_value = 1 - stats.f.cdf(F, q, nobs - df)

    return F, (q, nobs - df), p_value

def calc_f_stat(nw_lags, r2, r2_adj, cols, beta, var_beta, nobs, df,
                window=None, T=None):
    if nw_lags is None:
        F = r2 / (r2 - r2_adj)

        q = len(cols)
        if 'intercept' in cols:
            q -= 1

        if window is None:
            shape = q, nobs - df
            p_value = 1 - stats.f.cdf(F, shape[0], shape[1])
            return F, shape, p_value

        results = []

        start = window - 1
        for i in xrange(start, T):
            shape = q, nobs[i - start] - df[i - start]
            p_value = 1 - stats.f.cdf(F[i - start], shape[0], shape[1])
            result = F[i - start], shape, p_value
            results.append(result)

        return results

    k = len(cols)

    R = np.eye(k)
    r = np.zeros((k, 1))

    intercept = cols.indexMap.get('intercept')

    if intercept is not None:
        R = np.concatenate((R[0 : intercept], R[intercept + 1:]))
        r = np.concatenate((r[0 : intercept], r[intercept + 1:]))

    if window is None:
        return calc_F(R, r, beta, var_beta, nobs, df)

    results = []

    start = window - 1
    for i in xrange(start, T):
        b = beta[i - start]
        vb = var_beta[i - start]
        n = nobs[i - start]
        d = df[i - start]
        result = calc_F(R, r, b, vb, n, d)
        results.append(result)

    return results
