from pandas.core.base import StringMixin
from pandas.compat import StringIO, range

import numpy as np

from pandas.core.api import Series, DataFrame
import pandas.stats.common as common
from pandas.util.decorators import cache_readonly

# flake8: noqa

def fama_macbeth(**kwargs):
    """Runs Fama-MacBeth regression.

    Parameters
    ----------
    Takes the same arguments as a panel OLS, in addition to:

    nw_lags_beta: int
       Newey-West adjusts the betas by the given lags
       """
    window_type = kwargs.get('window_type')
    if window_type is None:
        klass = FamaMacBeth
    else:
        klass = MovingFamaMacBeth

    return klass(**kwargs)


class FamaMacBeth(StringMixin):

    def __init__(self, y, x, intercept=True, nw_lags=None,
                 nw_lags_beta=None,
                 entity_effects=False, time_effects=False, x_effects=None,
                 cluster=None, dropped_dummies=None, verbose=False):
        import warnings
        warnings.warn("The pandas.stats.fama_macbeth module is deprecated and will be "
                      "removed in a future version. We refer to external packages "
                      "like statsmodels, see here: "
                      "http://www.statsmodels.org/stable/index.html",
                      FutureWarning, stacklevel=4)

        if dropped_dummies is None:
            dropped_dummies = {}
        self._nw_lags_beta = nw_lags_beta

        from pandas.stats.plm import MovingPanelOLS
        self._ols_result = MovingPanelOLS(
            y=y, x=x, window_type='rolling', window=1,
            intercept=intercept,
            nw_lags=nw_lags, entity_effects=entity_effects,
            time_effects=time_effects, x_effects=x_effects, cluster=cluster,
            dropped_dummies=dropped_dummies, verbose=verbose)

        self._cols = self._ols_result._x.columns

    @cache_readonly
    def _beta_raw(self):
        return self._ols_result._beta_raw

    @cache_readonly
    def _stats(self):
        return _calc_t_stat(self._beta_raw, self._nw_lags_beta)

    @cache_readonly
    def _mean_beta_raw(self):
        return self._stats[0]

    @cache_readonly
    def _std_beta_raw(self):
        return self._stats[1]

    @cache_readonly
    def _t_stat_raw(self):
        return self._stats[2]

    def _make_result(self, result):
        return Series(result, index=self._cols)

    @cache_readonly
    def mean_beta(self):
        return self._make_result(self._mean_beta_raw)

    @cache_readonly
    def std_beta(self):
        return self._make_result(self._std_beta_raw)

    @cache_readonly
    def t_stat(self):
        return self._make_result(self._t_stat_raw)

    @cache_readonly
    def _results(self):
        return {
            'mean_beta': self._mean_beta_raw,
            'std_beta': self._std_beta_raw,
            't_stat': self._t_stat_raw,
        }

    @cache_readonly
    def _coef_table(self):
        buffer = StringIO()
        buffer.write('%13s %13s %13s %13s %13s %13s\n' %
                     ('Variable', 'Beta', 'Std Err', 't-stat', 'CI 2.5%', 'CI 97.5%'))
        template = '%13s %13.4f %13.4f %13.2f %13.4f %13.4f\n'

        for i, name in enumerate(self._cols):
            if i and not (i % 5):
                buffer.write('\n' + common.banner(''))

            mean_beta = self._results['mean_beta'][i]
            std_beta = self._results['std_beta'][i]
            t_stat = self._results['t_stat'][i]
            ci1 = mean_beta - 1.96 * std_beta
            ci2 = mean_beta + 1.96 * std_beta

            values = '(%s)' % name, mean_beta, std_beta, t_stat, ci1, ci2

            buffer.write(template % values)

        if self._nw_lags_beta is not None:
            buffer.write('\n')
            buffer.write('*** The Std Err, t-stat are Newey-West '
                         'adjusted with Lags %5d\n' % self._nw_lags_beta)

        return buffer.getvalue()

    def __unicode__(self):
        return self.summary

    @cache_readonly
    def summary(self):
        template = """
----------------------Summary of Fama-MacBeth Analysis-------------------------

Formula: Y ~ %(formulaRHS)s
# betas : %(nu)3d

----------------------Summary of Estimated Coefficients------------------------
%(coefTable)s
--------------------------------End of Summary---------------------------------
"""
        params = {
            'formulaRHS': ' + '.join(self._cols),
            'nu': len(self._beta_raw),
            'coefTable': self._coef_table,
        }

        return template % params


class MovingFamaMacBeth(FamaMacBeth):

    def __init__(self, y, x, window_type='rolling', window=10,
                 intercept=True, nw_lags=None, nw_lags_beta=None,
                 entity_effects=False, time_effects=False, x_effects=None,
                 cluster=None, dropped_dummies=None, verbose=False):
        if dropped_dummies is None:
            dropped_dummies = {}
        self._window_type = common._get_window_type(window_type)
        self._window = window

        FamaMacBeth.__init__(
            self, y=y, x=x, intercept=intercept,
            nw_lags=nw_lags, nw_lags_beta=nw_lags_beta,
            entity_effects=entity_effects, time_effects=time_effects,
            x_effects=x_effects, cluster=cluster,
            dropped_dummies=dropped_dummies, verbose=verbose)

        self._index = self._ols_result._index
        self._T = len(self._index)

    @property
    def _is_rolling(self):
        return self._window_type == 'rolling'

    def _calc_stats(self):
        mean_betas = []
        std_betas = []
        t_stats = []

        # XXX

        mask = self._ols_result._rolling_ols_call[2]
        obs_total = mask.astype(int).cumsum()

        start = self._window - 1
        betas = self._beta_raw
        for i in range(start, self._T):
            if self._is_rolling:
                begin = i - start
            else:
                begin = 0

            B = betas[max(obs_total[begin] - 1, 0): obs_total[i]]
            mean_beta, std_beta, t_stat = _calc_t_stat(B, self._nw_lags_beta)
            mean_betas.append(mean_beta)
            std_betas.append(std_beta)
            t_stats.append(t_stat)

        return np.array([mean_betas, std_betas, t_stats])

    _stats = cache_readonly(_calc_stats)

    def _make_result(self, result):
        return DataFrame(result, index=self._result_index, columns=self._cols)

    @cache_readonly
    def _result_index(self):
        mask = self._ols_result._rolling_ols_call[2]
        # HACK XXX
        return self._index[mask.cumsum() >= self._window]

    @cache_readonly
    def _results(self):
        return {
            'mean_beta': self._mean_beta_raw[-1],
            'std_beta': self._std_beta_raw[-1],
            't_stat': self._t_stat_raw[-1],
        }


def _calc_t_stat(beta, nw_lags_beta):
    N = len(beta)
    B = beta - beta.mean(0)
    C = np.dot(B.T, B) / N

    if nw_lags_beta is not None:
        for i in range(nw_lags_beta + 1):

            cov = np.dot(B[i:].T, B[:(N - i)]) / N
            weight = i / (nw_lags_beta + 1)
            C += 2 * (1 - weight) * cov

    mean_beta = beta.mean(0)
    std_beta = np.sqrt(np.diag(C)) / np.sqrt(N)
    t_stat = mean_beta / std_beta

    return mean_beta, std_beta, t_stat
