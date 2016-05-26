# flake8: noqa

from __future__ import print_function

import pandas.util.testing as tm

from pandas.compat import range
import nose
import unittest

raise nose.SkipTest('skipping this for now')

try:
    import statsmodels.tsa.var as sm_var
    import statsmodels as sm
except ImportError:
    import scikits.statsmodels.tsa.var as sm_var
    import scikits.statsmodels as sm


import pandas.stats.var as _pvar
reload(_pvar)
from pandas.stats.var import VAR

DECIMAL_6 = 6
DECIMAL_5 = 5
DECIMAL_4 = 4
DECIMAL_3 = 3
DECIMAL_2 = 2


class CheckVAR(object):

    def test_params(self):
        tm.assert_almost_equal(self.res1.params, self.res2.params, DECIMAL_3)

    def test_neqs(self):
        tm.assert_numpy_array_equal(self.res1.neqs, self.res2.neqs)

    def test_nobs(self):
        tm.assert_numpy_array_equal(self.res1.avobs, self.res2.nobs)

    def test_df_eq(self):
        tm.assert_numpy_array_equal(self.res1.df_eq, self.res2.df_eq)

    def test_rmse(self):
        results = self.res1.results
        for i in range(len(results)):
            tm.assert_almost_equal(results[i].mse_resid ** .5,
                                   eval('self.res2.rmse_' + str(i + 1)),
                                   DECIMAL_6)

    def test_rsquared(self):
        results = self.res1.results
        for i in range(len(results)):
            tm.assert_almost_equal(results[i].rsquared,
                                   eval('self.res2.rsquared_' + str(i + 1)),
                                   DECIMAL_3)

    def test_llf(self):
        results = self.res1.results
        tm.assert_almost_equal(self.res1.llf, self.res2.llf, DECIMAL_2)
        for i in range(len(results)):
            tm.assert_almost_equal(results[i].llf,
                                   eval('self.res2.llf_' + str(i + 1)),
                                   DECIMAL_2)

    def test_aic(self):
        tm.assert_almost_equal(self.res1.aic, self.res2.aic)

    def test_bic(self):
        tm.assert_almost_equal(self.res1.bic, self.res2.bic)

    def test_hqic(self):
        tm.assert_almost_equal(self.res1.hqic, self.res2.hqic)

    def test_fpe(self):
        tm.assert_almost_equal(self.res1.fpe, self.res2.fpe)

    def test_detsig(self):
        tm.assert_almost_equal(self.res1.detomega, self.res2.detsig)

    def test_bse(self):
        tm.assert_almost_equal(self.res1.bse, self.res2.bse, DECIMAL_4)


class Foo(object):

    def __init__(self):
        data = sm.datasets.macrodata.load()
        data = data.data[['realinv', 'realgdp', 'realcons']].view((float, 3))
        data = diff(log(data), axis=0)
        self.res1 = VAR2(endog=data).fit(maxlag=2)
        from results import results_var
        self.res2 = results_var.MacrodataResults()


if __name__ == '__main__':
    unittest.main()
