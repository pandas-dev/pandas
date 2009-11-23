"""
Simple OLS.
"""
from __future__ import division

from StringIO import StringIO

import numpy as np
from scipy import stats

from pandas.core.api import DataFrame, DataMatrix, Series
from pandas.util.decorators import cache_readonly
import pandas.lib.tseries as tseries
import pandas.stats.common as common
import pandas.stats.math as math

class OLS(object):
    """
    Runs a full sample ordinary least squares regression

    Parameters
    ----------
    y: Series
    x: Series, DataFrame, or dict of Series
    intercept: bool
        True if you want an intercept.
    nw_lags: None or int
        Number of Newey-West lags.
    """
    def __init__(self, y, x, intercept=True, nw_lags=None, nw_overlap=False):
        import scikits.statsmodels as sm

        self._x_orig = x
        self._y_orig = y
        self._intercept = intercept
        self._nw_lags = nw_lags
        self._nw_overlap = nw_overlap

        (self._y, self._x, self._x_filtered,
         self._index, self._time_has_obs) = self._prepare_data()

        self._x_raw = self._x.values
        self._y_raw = self._y.view(np.ndarray)
        self._nobs = len(self._y_raw)

        self.sm_ols = sm.OLS(self._y_raw, self._x_raw).fit()

    def _prepare_data(self):
        """
        Filters the data and sets up an intercept if necessary.

        Returns
        -------
        (DataFrame, Series).
        """
        (y, x, x_filtered,
         union_index, valid) = _filter_data(self._y_orig, self._x_orig)

        if self._intercept:
            x['intercept'] = x_filtered['intercept'] = 1.

        return y, x, x_filtered, union_index, valid

    @property
    def nobs(self):
        return self._nobs

    @property
    def nw_lags(self):
        return self._nw_lags

    @property
    def x(self):
        """Returns the filtered x used in the regression."""
        return self._x

    @property
    def y(self):
        """Returns the filtered y used in the regression."""
        return self._y

    @cache_readonly
    def _beta_raw(self):
        """Runs the regression and returns the beta."""
        return self.sm_ols.params

    @cache_readonly
    def beta(self):
        """Returns the betas in Series form."""
        return Series(self._beta_raw, index=self._x.cols())

    @cache_readonly
    def _df_raw(self):
        """Returns the degrees of freedom."""
        return math.rank(self._x_raw)

    @cache_readonly
    def df(self):
        """Returns the degrees of freedom.

        This equals the rank of the X matrix."""
        return self._df_raw

    @cache_readonly
    def _df_model_raw(self):
        """Returns the raw model degrees of freedom."""
        return self.sm_ols.df_model

    @cache_readonly
    def df_model(self):
        """Returns the degrees of freedom of the model."""
        return self._df_model_raw

    @cache_readonly
    def _df_resid_raw(self):
        """Returns the raw residual degrees of freedom."""
        return self.sm_ols.df_resid

    @cache_readonly
    def df_resid(self):
        """Returns the degrees of freedom of the residuals."""
        return self._df_resid_raw

    @cache_readonly
    def _f_stat_raw(self):
        """Returns the raw f-stat value."""
        return math.calc_f_stat(self._nw_lags, self._r2_raw, self._r2_adj_raw,
                                self._x.columns, self._beta_raw,
                                self._var_beta_raw,self._nobs, self.df)
    @cache_readonly
    def f_stat(self):
        """Returns the f-stat value."""
        return common.f_stat_to_dict(self._f_stat_raw)

    def f_test(self, hypothesis):
        """Runs the F test, given a joint hypothesis.  The hypothesis is
        represented by a collection of equations, in the form

        A*x_1+B*x_2=C

        You must provide the coefficients even if they're 1.  No spaces.

        The equations can be passed as either a single string or a
        list of strings.

        Examples:
        o = ols(...)
        o.f_test('1*x1+2*x2=0,1*x3=0')
        o.f_test(['1*x1+2*x2=0','1*x3=0'])
        """

        x_names = self._x.columns

        R = []
        r = []

        if isinstance(hypothesis, str):
            eqs = hypothesis.split(',')
        elif isinstance(hypothesis, list):
            eqs = hypothesis
        else:
            raise Exception('hypothesis must be either string or list')
        for equation in eqs:
            row = np.zeros(len(x_names))
            lhs, rhs = equation.split('=')
            for s in lhs.split('+'):
                ss = s.split('*')
                coeff = float(ss[0])
                x_name = ss[1]
                idx = x_names.indexMap[x_name]
                row[idx] = coeff
            rhs = float(rhs)

            R.append(row)
            r.append(rhs)

        R = np.array(R)
        q = len(r)
        r = np.array(r).reshape(q, 1)

        result = math.calc_F(R, r, self._beta_raw, self._var_beta_raw,
                             self.nobs, self.df)

        return common.f_stat_to_dict(result)

    @cache_readonly
    def _p_value_raw(self):
        """Returns the raw p values."""
        t_stat = self._t_stat_raw
        p_value = 2 * (1 - stats.t.cdf(np.fabs(t_stat),
            (self._nobs - self._df_raw)))
        return np.array(p_value)

    @cache_readonly
    def p_value(self):
        """Returns the p values."""
        index = self.beta.index
        return Series(self._p_value_raw, index=index)

    @cache_readonly
    def _r2_raw(self):
        """Returns the raw r-squared values."""
        return self.sm_ols.rsquared

    @cache_readonly
    def r2(self):
        """Returns the r-squared values."""
        return self._r2_raw

    @cache_readonly
    def _r2_adj_raw(self):
        """Returns the raw r-squared adjusted values."""
        return self.sm_ols.rsquared_adj

    @cache_readonly
    def r2_adj(self):
        """Returns the r-squared adjusted values."""
        return self._r2_adj_raw

    @cache_readonly
    def _resid_raw(self):
        """Returns the raw residuals."""
        return self.sm_ols.resid

    @cache_readonly
    def resid(self):
        """Returns the residuals."""
        index = self._x.index

        return Series(self._resid_raw, index=index)

    @cache_readonly
    def _rmse_raw(self):
        """Returns the raw rmse values."""
        return np.sqrt(self.sm_ols.mse_resid)

    @cache_readonly
    def rmse(self):
        """Returns the rmse value."""
        return self._rmse_raw

    @cache_readonly
    def _std_err_raw(self):
        """Returns the raw standard err values."""
        return np.nan_to_num(np.sqrt(np.diag(self._var_beta_raw)))

    @cache_readonly
    def std_err(self):
        """Returns the standard err values of the betas."""
        index = self.beta.index
        return Series(self._std_err_raw, index=index)

    @cache_readonly
    def _t_stat_raw(self):
        """Returns the raw t-stat value."""
        return np.nan_to_num(self._beta_raw / self._std_err_raw)

    @cache_readonly
    def t_stat(self):
        """Returns the t-stat values of the betas."""
        return Series(self._t_stat_raw, index=self.beta.index)

    @cache_readonly
    def _var_beta_raw(self):
        """Returns the raw covariance of beta."""
        result = math.calc_var_beta(x=self._x_raw, y=self._y_raw,
                                    nw_lags=self._nw_lags, rmse=self._rmse_raw,
                                    beta=self._beta_raw, nobs=self.nobs,
                                    df=self._df_raw, nw_overlap=self._nw_overlap)
        return np.array(result)

    @cache_readonly
    def var_beta(self):
        """Returns the variance-covariance matrix of beta."""
        return DataMatrix(self._var_beta_raw, index=self.beta.index,
                          columns=self.beta.index)

    @cache_readonly
    def _y_fitted_raw(self):
        """Returns the raw fitted y values."""
        return self.sm_ols.fittedvalues

    @cache_readonly
    def y_fitted(self):
        """Returns the fitted y values.  This equals BX."""
        index = self._x_filtered.index
        return Series(self._y_fitted_raw, index=index)

    @cache_readonly
    def _y_predict_raw(self):
        """Returns the raw predicted y values."""
        return self._y_fitted_raw

    @cache_readonly
    def y_predict(self):
        """Returns the predicted y values.

        For in-sample, this is same as y_fitted."""
        return self.y_fitted

    RESULT_FIELDS = ['r2', 'r2_adj', 'df', 'df_model', 'df_resid', 'rmse',
                     'f_stat', 'beta', 'std_err', 't_stat', 'p_value']

    @cache_readonly
    def _results(self):
        results = {}
        for result in self.RESULT_FIELDS:
            results[result] = getattr(self, result)

        return results

    @cache_readonly
    def _coef_table(self):
        buf = StringIO()

        buf.write('%14s %10s %10s %10s %10s %10s %10s\n' %
                  ('Variable', 'Coef', 'Std Err', 't-stat',
                   'p-value', 'CI 2.5%', 'CI 97.5%'))
        buf.write(common.banner(''))
        coef_template = '\n%14s %10.4f %10.4f %10.2f %10.4f %10.4f %10.4f'

        results = self._results

        beta = results['beta']

        for i, name in enumerate(beta.index):
            if i and not (i % 5):
                buf.write('\n' + common.banner(''))

            std_err = results['std_err'][name]
            CI1 = beta[name] - 1.96 * std_err
            CI2 = beta[name] + 1.96 * std_err

            t_stat = results['t_stat'][name]
            p_value = results['p_value'][name]

            line = coef_template % (name,
                beta[name], std_err, t_stat, p_value, CI1, CI2)

            buf.write(line)

        if self.nw_lags is not None:
            buf.write('\n')
            buf.write('*** The calculations are Newey-West '
                      'adjusted with lags %5d\n' % self.nw_lags)

        return buf.getvalue()

    @cache_readonly
    def summary_as_matrix(self):
        """Returns the formatted results of the OLS as a DataMatrix."""
        results = self._results
        beta = results['beta']
        data = {'beta' : results['beta'],
                't-stat' : results['t_stat'],
                'p-value' : results['p_value'],
                'std err' : results['std_err']}
        return DataFrame(data, beta.index).T

    @cache_readonly
    def summary(self):
        """
        This returns the formatted result of the OLS computation
        """
        template = """
%(bannerTop)s

Formula: Y ~ %(formula)s

Number of Observations:         %(nobs)d
Number of Degrees of Freedom:   %(df)d

R-squared:     %(r2)10.4f
Adj R-squared: %(r2_adj)10.4f

Rmse:          %(rmse)10.4f

F-stat %(f_stat_shape)s: %(f_stat)10.4f, p-value: %(f_stat_p_value)10.4f

Degrees of Freedom: model %(df_model)d, resid %(df_resid)d

%(bannerCoef)s
%(coef_table)s
%(bannerEnd)s
"""
        coef_table = self._coef_table

        results = self._results

        f_stat = results['f_stat']

        bracketed = ['<%s>' % c for c in results['beta'].index]

        formula = StringIO()
        formula.write(bracketed[0])
        tot = len(bracketed[0])
        line = 1
        for coef in bracketed[1:]:
            tot = tot + len(coef) + 3

            if tot // (68 * line):
                formula.write('\n' + ' ' * 12)
                line += 1

            formula.write(' + ' + coef)

        params = {
            'bannerTop' : common.banner('Summary of Regression Analysis'),
            'bannerCoef' : common.banner('Summary of Estimated Coefficients'),
            'bannerEnd' : common.banner('End of Summary'),
            'formula' : formula.getvalue(),
            'r2' : results['r2'],
            'r2_adj' : results['r2_adj'],
            'nobs' : self.nobs,
            'df'  : results['df'],
            'df_model'  : results['df_model'],
            'df_resid'  : results['df_resid'],
            'coef_table' : coef_table,
            'rmse' : results['rmse'],
            'f_stat' : f_stat['f-stat'],
            'f_stat_shape' : '(%d, %d)' % (f_stat['DF X'], f_stat['DF Resid']),
            'f_stat_p_value' : f_stat['p-value'],
        }

        return template % params

    def __repr__(self):
        return self.summary

    @property
    def _total_times(self):
        return self._time_has_obs.sum()


class MovingOLS(OLS):
    """
    Runs a rolling/expanding simple OLS.

    Parameters
    ----------
    y: Series
    x: Series, DataFrame, or dict of Series
    intercept: bool
        True if you want an intercept.
    nw_lags: None or int
        Number of Newey-West lags.
    window_type: int
        FULL_SAMPLE, ROLLING, EXPANDING.  FULL_SAMPLE by default.
    window: int
        size of window (for rolling/expanding OLS)
    """
    def __init__(self, y, x, window_type=common.ROLLING,
                 window=10, min_periods=None, intercept=True,
                 nw_lags=None, nw_overlap=False):

        self._args = dict(intercept=intercept, nw_lags=nw_lags,
                          nw_overlap=nw_overlap)

        OLS.__init__(self, y=y, x=x, **self._args)

        self._window_type = common._get_window_type(window_type)
        self._window = window
        self._min_periods = window if min_periods is None else min_periods

    @property
    def _is_rolling(self):
        return self._window_type == common.ROLLING

    @cache_readonly
    def beta(self):
        """Returns the betas in Series/DataMatrix form."""
        return DataMatrix(self._beta_raw,
                          index=self._result_index,
                          columns=self._x.cols())

    @cache_readonly
    def _beta_raw(self):
        """Runs the regression and returns the beta."""
        beta, indices = self._rolling_ols_call

        return beta[indices]

    @cache_readonly
    def _result_index(self):
        return self._index[self._valid_indices]

    @property
    def _valid_indices(self):
        return self._rolling_ols_call[1]

    @cache_readonly
    def _rolling_ols_call(self):
        return self._calc_betas()

    def _calc_betas(self):
        N = len(self._index)
        K = len(self._x.cols())

        betas = np.empty((N, K), dtype=float)
        betas[:] = np.NaN

        valid = self._time_has_obs
        enough = self._enough_obs
        window = self._window

        # Use transformed (demeaned) Y, X variables
        cum_xx = self._cum_xx(self._x)
        cum_xy = self._cum_xy(self._x, self._y)

        for i in xrange(N):
            # XXX
            if not valid[i] or not enough[i]:
                continue

            xx = cum_xx[i]
            xy = cum_xy[i]
            if self._is_rolling and i >= window:
                xx = xx - cum_xx[i - window]
                xy = xy - cum_xy[i - window]

            betas[i] = math.solve(xx, xy)

        have_betas = np.arange(N)[-np.isnan(betas).any(axis=1)]

        return betas, have_betas

    def _cum_xx(self, x):
        K = len(x.cols())
        valid = self._time_has_obs
        cum_xx = []

        last = np.zeros((K, K))
        for i, date in enumerate(self._index):
            if not valid[i]:
                cum_xx.append(last)
                continue

            xs = x.getXS(date)
            xx = last = last + np.outer(xs, xs)
            cum_xx.append(xx)

        return cum_xx

    def _cum_xy(self, x, y):
        valid = self._time_has_obs
        cum_xy = []
        last = len(x.cols())
        for i, date in enumerate(self._index):
            if not valid[i]:
                cum_xy.append(last)
                continue

            xs = np.asarray(x.getXS(date))
            xy = last = last + xs * y[date]
            cum_xy.append(xy)

        return cum_xy

    @cache_readonly
    def rank(self):
        return Series(self._rank_raw, index=self._result_index)

    @cache_readonly
    def _rank_raw(self):
        rank = self._rolling_rank

        return rank[self._valid_indices]

    @cache_readonly
    def _rolling_rank(self):
        dates = self._index
        enough = self._enough_obs
        window = self._window

        ranks = np.empty(len(dates), dtype=float)
        ranks[:] = np.NaN
        for i, date in enumerate(dates):
            if self._is_rolling and i >= window:
                prior_date = dates[i - window + 1]
            else:
                prior_date = dates[0]

            x_slice = self._x.truncate(before=prior_date, after=date)
            ranks[i] = math.rank(x_slice.values)

        return ranks

    @cache_readonly
    def _df_raw(self):
        """Returns the degrees of freedom."""
        return self._rank_raw

    @cache_readonly
    def df(self):
        """Returns the degrees of freedom."""
        index = self.beta.index
        return Series(self._df_raw, index=index)

    @cache_readonly
    def _df_model_raw(self):
        """Returns the raw model degrees of freedom."""
        return self._df_raw - 1

    @cache_readonly
    def df_model(self):
        """Returns the model degrees of freedom."""
        index = self.beta.index

        return Series(self._df_model_raw, index=index)

    @cache_readonly
    def _df_resid_raw(self):
        """Returns the raw residual degrees of freedom."""
        return self._window_nobs - self._df_raw

    @cache_readonly
    def df_resid(self):
        """Returns the residual degrees of freedom."""
        index = self.beta.index

        return Series(self._df_resid_raw, index=index)

    @cache_readonly
    def _f_stat_raw(self):
        """Returns the raw f-stat value."""
        return math.calc_f_stat(self._nw_lags, self._r2_raw, self._r2_adj_raw,
                                self._x.columns, self._beta_raw,
                                self._var_beta_raw, self._window_nobs, self.df,
                                self._window, self._nobs)

    @cache_readonly
    def f_stat(self):
        """Returns the f-stat value."""
        f_stat_dicts = dict((date, common.f_stat_to_dict(f_stat))
                            for date, f_stat in zip(self.beta.index,
                                                    self._f_stat_raw))

        return DataFrame.fromDict(f_stat_dicts).T

    def f_test(self, hypothesis):
        raise Exception('f_test not supported for rolling/expanding OLS')

    @cache_readonly
    def _p_value_raw(self):
        """Returns the raw p values."""
        get_prob = lambda a, b: 2 * (1 - stats.t.cdf(a, b))

        result = starmap(get_prob,
                         izip(np.fabs(self._t_stat_raw), self._df_resid_raw))

        result = np.array(list(result))

        return result

    @cache_readonly
    def p_value(self):
        """Returns the p values."""
        cols = self.beta.cols()
        rows = self.beta.index
        return DataMatrix(self._p_value_raw, columns=cols, index=rows)

    @cache_readonly
    def _r2_raw(self):
        """Returns the raw r-squared values."""
        window = self._window
        _r2 = []

        X = self._x_raw
        Y = self._y_raw

        for i in xrange(window - 1, self._nobs):
            if self._window_type == common.EXPANDING:
                section = slice(None, i + 1)
            else:
                section = slice(i - window + 1, i + 1)
            SSerr = ((Y[section] - np.dot(X[section],
                      self._beta_raw[i - window + 1])) ** 2).sum()
            SStotal = ((Y[section] - np.mean(Y[section])) ** 2).sum()
            _r2.append(1 - SSerr / SStotal)

        return np.array(_r2)

    @cache_readonly
    def r2(self):
        """Returns the r-squared values."""
        index = self.beta.index

        return Series(self._r2_raw, index=index)

    @cache_readonly
    def _resid_raw(self):
        """Returns the raw residuals."""
        start = self._window - 1
        return self._y_raw[start:] - self._y_fitted_raw

    @cache_readonly
    def resid(self):
        """Returns the residuals."""
        index = self.beta.index
        return Series(self._resid_raw, index=index)

    @cache_readonly
    def _r2_adj_raw(self):
        """Returns the raw r-squared adjusted values."""
        Pa = []
        start = self._window - 1
        for i in xrange(start, self._nobs):
            if self._window_type == common.EXPANDING:
                nobs = i + 1
            else:
                nobs = self._window

            Pa.append((nobs - 1) / (nobs - self._df_raw[i - start]))

        _r2_adj_raw = 1 - (1 - self._r2_raw) * (Pa)

        return np.array(_r2_adj_raw)

    @cache_readonly
    def r2_adj(self):
        """Returns the r-squared adjusted values."""
        index = self.r2.index

        return Series(self._r2_adj_raw, index=index)


    @cache_readonly
    def _rmse_raw(self):
        """Returns the raw rmse values."""
        window = self._window
        X = self._x_raw
        Y = self._y_raw

        results = []
        start = window - 1
        for i in xrange(start, self._nobs):
            if self._window_type == common.EXPANDING:
                section = slice(0, i + 1)
            else:
                section = slice(i - start, i + 1)
            estimate = np.dot(X[section], self._beta_raw[i - start])
            s = ((Y[section] - estimate) ** 2).sum()
            df = self._df_raw[i - start]
            nobs = len(Y[section])
            result = np.sqrt(s / (nobs - df))
            results.append(result)

        return np.array(results)

    @cache_readonly
    def rmse(self):
        """Returns the rmse values."""
        index = self.beta.index

        return Series(self._rmse_raw, index=index)

    @cache_readonly
    def _std_err_raw(self):
        """Returns the raw standard err values."""
        results = []
        for i in xrange(len(self._var_beta_raw)):
            results.append(np.sqrt(np.diag(self._var_beta_raw[i])))

        return np.nan_to_num(np.array(results))

    @cache_readonly
    def std_err(self):
        """Returns the standard err values."""
        index = self.beta.index
        cols = self.beta.cols()
        return DataMatrix(self._std_err_raw, columns=cols, index=index)

    @cache_readonly
    def _t_stat_raw(self):
        """Returns the raw t-stat value."""
        results = []
        start = self._window - 1
        for i in xrange(start, self._nobs):
            results.append(np.nan_to_num(self._beta_raw[i - start] /
                self._std_err_raw[i - start]))

        return np.array(results)

    @cache_readonly
    def t_stat(self):
        """Returns the t-stat value."""
        cols = self.beta.cols()
        rows = self.beta.index
        return DataMatrix(self._t_stat_raw, columns=cols, index=rows)

    @cache_readonly
    def _var_beta_raw(self):
        """Returns the raw covariance of beta."""
        result = math.calc_var_beta(x=self._x_raw, y=self._y_raw,
                                    window_type=self._window_type,
                                    window=self._window, nw_lags=self._nw_lags,
                                    rmse=self._rmse_raw, beta=self._beta_raw,
                                    nobs=self.nobs, df=self._df_raw,
                                    nw_overlap=self._nw_overlap)
        return np.array(result)

    @cache_readonly
    def var_beta(self):
        """Returns the covariance of beta."""
        result = []
        for i in xrange(len(self._var_beta_raw)):
            result.append(DataMatrix(
                self._var_beta_raw[i], columns=self.beta.cols(),
                index=self.beta.cols()))

        return Series(result, index=self.beta.index)

    @cache_readonly
    def _y_fitted_raw(self):
        """Returns the raw fitted y values."""
        start = self._window - 1
        return (self._x_raw[start:] * self._beta_raw).sum(1)

    @cache_readonly
    def y_fitted(self):
        """Returns the fitted y values."""
        index = self.beta.index
        return Series(self._y_fitted_raw, index=index)

    @cache_readonly
    def _y_predict_raw(self):
        """Returns the raw predicted y values."""
        bx = self._beta_raw[: -1] * self._x_raw[self._window :]
        return bx.sum(1)

    @cache_readonly
    def y_predict(self):
        """Returns the predicted y values."""
        index = self.beta.index[1 :]
        return Series(self._y_predict_raw, index=index)

    @cache_readonly
    def _results(self):
        results = {}
        for result in self.RESULT_FIELDS:
            value = getattr(self, result)
            if isinstance(value, Series):
                value = value[self.beta.index[-1]]
            elif isinstance(value, DataFrame):
                value = value.getXS(self.beta.index[-1])
            else:
                raise Exception('Problem retrieving %s' % result)
            results[result] = value

        return results

    @cache_readonly
    def _window_nobs(self):
        results = []
        start = self._window - 1
        for i in xrange(start, self._nobs):
            if self._window_type == common.EXPANDING:
                result = i + 1
            else:
                result = self._window

            results.append(result)

        return results

    def _beta_matrix(self, lag=0):
        assert(lag >= 0)

        labels = self._y_trans.index.major_labels - lag
        indexer = self._valid_indices.searchsorted(labels, side='left')

        beta_matrix = self._beta_raw[indexer]
        beta_matrix[labels < 0] = np.NaN

        return beta_matrix

    @cache_readonly
    def _window_nobs_raw(self):
        if self._is_rolling:
            window = self._window
        else:
            # expanding case
            window = len(self._index)

        result = tseries.rolling_sum(self._time_has_obs, window,
                                     minp=1)

        return result.astype(int)

    @cache_readonly
    def _window_nobs(self):
        return self._window_nobs_raw[self._valid_indices]

    @cache_readonly
    def _enough_obs(self):
        # XXX: what's the best way to determine where to start?
        return self._window_nobs_raw >= max(self._min_periods,
                                            len(self._x.columns) + 1)

def _safe_update(d, other):
    """
    Combine dictionaries with non-overlapping keys
    """
    for k, v in other.iteritems():
        if k in d:
            raise Exception('Duplicate regressor: %s' % k)

        d[k] = v

def _combine_rhs(rhs):
    """
    Glue input X variables together while checking for potential
    duplicates
    """
    series = {}

    if isinstance(rhs, Series):
        series['x'] = rhs
    elif isinstance(rhs, DataFrame):
        _safe_update(series, rhs)
    elif isinstance(rhs, dict):
        for name, value in rhs.iteritems():
            if isinstance(value, Series):
                _safe_update(series, {name : value})
            elif isinstance(value, (dict, DataFrame)):
                _safe_update(series, value)
            else:
                raise Exception('Invalid RHS data type: %s' % type(value))
    else:
        raise Exception('Invalid RHS type: %s' % type(rhs))

    return series

def _filter_data(lhs, rhs):
    """
    Cleans the input for single OLS.

    Parameters
    ----------
    lhs: Series
        Dependent variable in the regression.
    rhs: dict, whose values are Series, DataFrame, or dict
        Explanatory variables of the regression.

    Returns
    -------
    Series, DataFrame
        Cleaned lhs and rhs
    """
    if not isinstance(lhs, Series):
        raise Exception('lhs must be a Series')

    combined_rhs = _combine_rhs(rhs)

    pre_filtered_rhs = DataMatrix.fromDict(combined_rhs).dropIncompleteRows()

    # Union of all indices
    combined_rhs['_y'] = lhs
    full_dataset = DataMatrix.fromDict(combined_rhs)

    index = full_dataset.index

    obs_count = full_dataset.count(axis=1, asarray=True)
    valid = obs_count == len(full_dataset.cols())

    filtered_rhs = full_dataset.reindex(index[valid])
    filtered_lhs = filtered_rhs.pop('_y')

    return filtered_lhs, filtered_rhs, pre_filtered_rhs, index, valid
