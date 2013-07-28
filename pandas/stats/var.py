from __future__ import division

from pandas.compat import range, lrange, zip, reduce
from pandas import compat
import numpy as np
from pandas.core.base import StringMixin
from pandas.util.decorators import cache_readonly
from pandas.core.frame import DataFrame
from pandas.core.panel import Panel
from pandas.core.series import Series
import pandas.stats.common as common
from pandas.stats.math import inv
from pandas.stats.ols import _combine_rhs


class VAR(StringMixin):
    """
    Estimates VAR(p) regression on multivariate time series data
    presented in pandas data structures.

    Parameters
    ----------
    data : DataFrame or dict of Series
    p : lags to include

    """

    def __init__(self, data, p=1, intercept=True):
        try:
            import statsmodels.tsa.vector_ar.api as sm_var
        except ImportError:
            import scikits.statsmodels.tsa.var as sm_var

        self._data = DataFrame(_combine_rhs(data))
        self._p = p

        self._columns = self._data.columns
        self._index = self._data.index

        self._intercept = intercept

    @cache_readonly
    def aic(self):
        """Returns the Akaike information criterion."""
        return self._ic['aic']

    @cache_readonly
    def bic(self):
        """Returns the Bayesian information criterion."""
        return self._ic['bic']

    @cache_readonly
    def beta(self):
        """
        Returns a DataFrame, where each column x1 contains the betas
        calculated by regressing the x1 column of the VAR input with
        the lagged input.

        Returns
        -------
        DataFrame
        """
        d = dict([(key, value.beta)
                  for (key, value) in compat.iteritems(self.ols_results)])
        return DataFrame(d)

    def forecast(self, h):
        """
        Returns a DataFrame containing the forecasts for 1, 2, ..., n time
        steps.  Each column x1 contains the forecasts of the x1 column.

        Parameters
        ----------
        n: int
            Number of time steps ahead to forecast.

        Returns
        -------
        DataFrame
        """
        forecast = self._forecast_raw(h)[:, 0, :]
        return DataFrame(forecast, index=lrange(1, 1 + h),
                         columns=self._columns)

    def forecast_cov(self, h):
        """
        Returns the covariance of the forecast residuals.

        Returns
        -------
        DataFrame
        """
        return [DataFrame(value, index=self._columns, columns=self._columns)
                for value in self._forecast_cov_raw(h)]

    def forecast_std_err(self, h):
        """
        Returns the standard errors of the forecast residuals.

        Returns
        -------
        DataFrame
        """
        return DataFrame(self._forecast_std_err_raw(h),
                         index=lrange(1, 1 + h), columns=self._columns)

    @cache_readonly
    def granger_causality(self):
        """Returns the f-stats and p-values from the Granger Causality Test.

        If the data consists of columns x1, x2, x3, then we perform the
        following regressions:

        x1 ~ L(x2, x3)
        x1 ~ L(x1, x3)
        x1 ~ L(x1, x2)

        The f-stats of these results are placed in the 'x1' column of the
        returned DataFrame.  We then repeat for x2, x3.

        Returns
        -------
        Dict, where 'f-stat' returns the DataFrame containing the f-stats,
        and 'p-value' returns the DataFrame containing the corresponding
        p-values of the f-stats.
        """
        from pandas.stats.api import ols
        from scipy.stats import f

        d = {}
        for col in self._columns:
            d[col] = {}
            for i in range(1, 1 + self._p):
                lagged_data = self._lagged_data[i].filter(
                    self._columns - [col])

                for key, value in compat.iteritems(lagged_data):
                    d[col][_make_param_name(i, key)] = value

        f_stat_dict = {}
        p_value_dict = {}

        for col, y in compat.iteritems(self._data):
            ssr_full = (self.resid[col] ** 2).sum()

            f_stats = []
            p_values = []

            for col2 in self._columns:
                result = ols(y=y, x=d[col2])

                resid = result.resid
                ssr_reduced = (resid ** 2).sum()

                M = self._p
                N = self._nobs
                K = self._k * self._p + 1
                f_stat = ((ssr_reduced - ssr_full) / M) / (ssr_full / (N - K))
                f_stats.append(f_stat)

                p_value = f.sf(f_stat, M, N - K)
                p_values.append(p_value)

            f_stat_dict[col] = Series(f_stats, self._columns)
            p_value_dict[col] = Series(p_values, self._columns)

        f_stat_mat = DataFrame(f_stat_dict)
        p_value_mat = DataFrame(p_value_dict)

        return {
            'f-stat': f_stat_mat,
            'p-value': p_value_mat,
        }

    @cache_readonly
    def ols_results(self):
        """
        Returns the results of the regressions:
        x_1 ~ L(X)
        x_2 ~ L(X)
        ...
        x_k ~ L(X)

        where X = [x_1, x_2, ..., x_k]
        and L(X) represents the columns of X lagged 1, 2, ..., n lags
        (n is the user-provided number of lags).

        Returns
        -------
        dict
        """
        from pandas.stats.api import ols

        d = {}
        for i in range(1, 1 + self._p):
            for col, series in compat.iteritems(self._lagged_data[i]):
                d[_make_param_name(i, col)] = series

        result = dict([(col, ols(y=y, x=d, intercept=self._intercept))
                       for col, y in compat.iteritems(self._data)])

        return result

    @cache_readonly
    def resid(self):
        """
        Returns the DataFrame containing the residuals of the VAR regressions.
        Each column x1 contains the residuals generated by regressing the x1
        column of the input against the lagged input.

        Returns
        -------
        DataFrame
        """
        d = dict([(col, series.resid)
                  for (col, series) in compat.iteritems(self.ols_results)])
        return DataFrame(d, index=self._index)

    @cache_readonly
    def summary(self):
        template = """
%(banner_top)s

Number of Observations:         %(nobs)d
AIC:                            %(aic).3f
BIC:                            %(bic).3f

%(banner_coef)s
%(coef_table)s
%(banner_end)s
"""
        params = {
            'banner_top': common.banner('Summary of VAR'),
            'banner_coef': common.banner('Summary of Estimated Coefficients'),
            'banner_end': common.banner('End of Summary'),
            'coef_table': self.beta,
            'aic': self.aic,
            'bic': self.bic,
            'nobs': self._nobs,
        }

        return template % params

    @cache_readonly
    def _alpha(self):
        """
        Returns array where the i-th element contains the intercept
        when regressing the i-th column of self._data with the lagged data.
        """
        if self._intercept:
            return self._beta_raw[-1]
        else:
            return np.zeros(self._k)

    @cache_readonly
    def _beta_raw(self):
        return np.array([list(self.beta[col].values()) for col in self._columns]).T

    def _trans_B(self, h):
        """
        Returns 0, 1, ..., (h-1)-th power of transpose of B as defined in
        equation (4) on p. 142 of the Stata 11 Time Series reference book.
        """
        result = [np.eye(1 + self._k * self._p)]

        row1 = np.zeros((1, 1 + self._k * self._p))
        row1[0, 0] = 1

        v = self._alpha.reshape((self._k, 1))
        row2 = np.hstack(tuple([v] + self._lag_betas))

        m = self._k * (self._p - 1)
        row3 = np.hstack((
            np.zeros((m, 1)),
            np.eye(m),
            np.zeros((m, self._k))
        ))

        trans_B = np.vstack((row1, row2, row3)).T

        result.append(trans_B)

        for i in range(2, h):
            result.append(np.dot(trans_B, result[i - 1]))

        return result

    @cache_readonly
    def _x(self):
        values = np.array([
            list(self._lagged_data[i][col].values())
            for i in range(1, 1 + self._p)
            for col in self._columns
        ]).T

        x = np.hstack((np.ones((len(values), 1)), values))[self._p:]

        return x

    @cache_readonly
    def _cov_beta(self):
        cov_resid = self._sigma

        x = self._x

        inv_cov_x = inv(np.dot(x.T, x))

        return np.kron(inv_cov_x, cov_resid)

    def _data_xs(self, i):
        """
        Returns the cross-section of the data at the given timestep.
        """
        return self._data.values[i]

    def _forecast_cov_raw(self, n):
        resid = self._forecast_cov_resid_raw(n)
        # beta = self._forecast_cov_beta_raw(n)

        # return [a + b for a, b in zip(resid, beta)]
        # TODO: ignore the beta forecast std err until it's verified

        return resid

    def _forecast_cov_beta_raw(self, n):
        """
        Returns the covariance of the beta errors for the forecast at
        1, 2, ..., n timesteps.
        """
        p = self._p

        values = self._data.values
        T = len(values) - self._p - 1

        results = []

        for h in range(1, n + 1):
            psi = self._psi(h)
            trans_B = self._trans_B(h)

            sum = 0

            cov_beta = self._cov_beta

            for t in range(T + 1):
                index = t + p
                y = values.take(lrange(index, index - p, -1), axis=0).ravel()
                trans_Z = np.hstack(([1], y))
                trans_Z = trans_Z.reshape(1, len(trans_Z))

                sum2 = 0
                for i in range(h):
                    ZB = np.dot(trans_Z, trans_B[h - 1 - i])

                    prod = np.kron(ZB, psi[i])
                    sum2 = sum2 + prod

                sum = sum + chain_dot(sum2, cov_beta, sum2.T)

            results.append(sum / (T + 1))

        return results

    def _forecast_cov_resid_raw(self, h):
        """
        Returns the covariance of the residual errors for the forecast at
        1, 2, ..., h timesteps.
        """
        psi_values = self._psi(h)
        sum = 0
        result = []
        for i in range(h):
            psi = psi_values[i]
            sum = sum + chain_dot(psi, self._sigma, psi.T)
            result.append(sum)

        return result

    def _forecast_raw(self, h):
        """
        Returns the forecast at 1, 2, ..., h timesteps in the future.
        """
        k = self._k
        result = []
        for i in range(h):
            sum = self._alpha.reshape(1, k)
            for j in range(self._p):
                beta = self._lag_betas[j]
                idx = i - j
                if idx > 0:
                    y = result[idx - 1]
                else:
                    y = self._data_xs(idx - 1)

                sum = sum + np.dot(beta, y.T).T
            result.append(sum)

        return np.array(result)

    def _forecast_std_err_raw(self, h):
        """
        Returns the standard error of the forecasts
        at 1, 2, ..., n timesteps.
        """
        return np.array([np.sqrt(np.diag(value))
                         for value in self._forecast_cov_raw(h)])

    @cache_readonly
    def _ic(self):
        """
        Returns the Akaike/Bayesian information criteria.
        """
        RSS = self._rss
        k = self._p * (self._k * self._p + 1)
        n = self._nobs * self._k

        return {'aic': 2 * k + n * np.log(RSS / n),
                'bic': n * np.log(RSS / n) + k * np.log(n)}

    @cache_readonly
    def _k(self):
        return len(self._columns)

    @cache_readonly
    def _lag_betas(self):
        """
        Returns list of B_i, where B_i represents the (k, k) matrix
        with the j-th row containing the betas of regressing the j-th
        column of self._data with self._data lagged i time steps.
        First element is B_1, second element is B_2, etc.
        """
        k = self._k
        b = self._beta_raw
        return [b[k * i: k * (i + 1)].T for i in range(self._p)]

    @cache_readonly
    def _lagged_data(self):
        return dict([(i, self._data.shift(i))
                     for i in range(1, 1 + self._p)])

    @cache_readonly
    def _nobs(self):
        return len(self._data) - self._p

    def _psi(self, h):
        """
        psi value used for calculating standard error.

        Returns [psi_0, psi_1, ..., psi_(h - 1)]
        """
        k = self._k
        result = [np.eye(k)]
        for i in range(1, h):
            result.append(sum(
                [np.dot(result[i - j], self._lag_betas[j - 1])
                 for j in range(1, 1 + i)
                 if j <= self._p]))

        return result

    @cache_readonly
    def _resid_raw(self):
        resid = np.array([self.ols_results[col]._resid_raw
                          for col in self._columns])
        return resid

    @cache_readonly
    def _rss(self):
        """Returns the sum of the squares of the residuals."""
        return (self._resid_raw ** 2).sum()

    @cache_readonly
    def _sigma(self):
        """Returns covariance of resids."""
        k = self._k
        n = self._nobs

        resid = self._resid_raw

        return np.dot(resid, resid.T) / (n - k)

    def __unicode__(self):
        return self.summary


def lag_select(data, max_lags=5, ic=None):
    """
    Select number of lags based on a variety of information criteria

    Parameters
    ----------
    data : DataFrame-like
    max_lags : int
        Maximum number of lags to evaluate
    ic : {None, 'aic', 'bic', ...}
        Choosing None will just display the results

    Returns
    -------
    None
    """
    pass


class PanelVAR(VAR):
    """
    Performs Vector Autoregression on panel data.

    Parameters
    ----------
    data: Panel or dict of DataFrame
    lags: int
    """
    def __init__(self, data, lags, intercept=True):
        self._data = _prep_panel_data(data)
        self._p = lags
        self._intercept = intercept

        self._columns = self._data.items

    @cache_readonly
    def _nobs(self):
        """Returns the number of observations."""
        _, timesteps, entities = self._data.values.shape
        return (timesteps - self._p) * entities

    @cache_readonly
    def _rss(self):
        """Returns the sum of the squares of the residuals."""
        return (self.resid.values ** 2).sum()

    def forecast(self, h):
        """
        Returns the forecasts at 1, 2, ..., n timesteps in the future.
        """
        forecast = self._forecast_raw(h).T.swapaxes(1, 2)
        index = lrange(1, 1 + h)
        w = Panel(forecast, items=self._data.items, major_axis=index,
                  minor_axis=self._data.minor_axis)
        return w

    @cache_readonly
    def resid(self):
        """
        Returns the DataFrame containing the residuals of the VAR regressions.
        Each column x1 contains the residuals generated by regressing the x1
        column of the input against the lagged input.

        Returns
        -------
        DataFrame
        """
        d = dict([(key, value.resid)
                  for (key, value) in compat.iteritems(self.ols_results)])
        return Panel.fromDict(d)

    def _data_xs(self, i):
        return self._data.values[:, i, :].T

    @cache_readonly
    def _sigma(self):
        """Returns covariance of resids."""
        k = self._k
        resid = _drop_incomplete_rows(self.resid.toLong().values)
        n = len(resid)
        return np.dot(resid.T, resid) / (n - k)


def _prep_panel_data(data):
    """Converts the given data into a Panel."""
    if isinstance(data, Panel):
        return data

    return Panel.fromDict(data)


def _drop_incomplete_rows(array):
    mask = np.isfinite(array).all(1)
    indices = np.arange(len(array))[mask]
    return array.take(indices, 0)


def _make_param_name(lag, name):
    return 'L%d.%s' % (lag, name)


def chain_dot(*matrices):
    """
    Returns the dot product of the given matrices.

    Parameters
    ----------
    matrices: argument list of ndarray
    """
    return reduce(lambda x, y: np.dot(y, x), matrices[::-1])
