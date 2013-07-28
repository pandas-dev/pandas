from __future__ import print_function
from numpy.testing import run_module_suite, assert_equal, TestCase

from pandas.util.testing import assert_almost_equal

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

try:
    import rpy2.robjects as robj
    from rpy2.robjects import r
    from rpy2.robjects.packages import importr
    import pandas.rpy.common as rpy
    vars = importr('vars')
    urca = importr('urca')
except ImportError:
    pass

DECIMAL_6 = 6
DECIMAL_5 = 5
DECIMAL_4 = 4
DECIMAL_3 = 3
DECIMAL_2 = 2


class CheckVAR(object):
    def test_params(self):
        assert_almost_equal(self.res1.params, self.res2.params, DECIMAL_3)

    def test_neqs(self):
        assert_equal(self.res1.neqs, self.res2.neqs)

    def test_nobs(self):
        assert_equal(self.res1.avobs, self.res2.nobs)

    def test_df_eq(self):
        assert_equal(self.res1.df_eq, self.res2.df_eq)

    def test_rmse(self):
        results = self.res1.results
        for i in range(len(results)):
            assert_almost_equal(results[i].mse_resid ** .5,
                                eval('self.res2.rmse_' + str(i + 1)), DECIMAL_6)

    def test_rsquared(self):
        results = self.res1.results
        for i in range(len(results)):
            assert_almost_equal(results[i].rsquared,
                                eval('self.res2.rsquared_' + str(i + 1)), DECIMAL_3)

    def test_llf(self):
        results = self.res1.results
        assert_almost_equal(self.res1.llf, self.res2.llf, DECIMAL_2)
        for i in range(len(results)):
            assert_almost_equal(results[i].llf,
                                eval('self.res2.llf_' + str(i + 1)), DECIMAL_2)

    def test_aic(self):
        assert_almost_equal(self.res1.aic, self.res2.aic)

    def test_bic(self):
        assert_almost_equal(self.res1.bic, self.res2.bic)

    def test_hqic(self):
        assert_almost_equal(self.res1.hqic, self.res2.hqic)

    def test_fpe(self):
        assert_almost_equal(self.res1.fpe, self.res2.fpe)

    def test_detsig(self):
        assert_almost_equal(self.res1.detomega, self.res2.detsig)

    def test_bse(self):
        assert_almost_equal(self.res1.bse, self.res2.bse, DECIMAL_4)


class Foo(object):
    def __init__(self):
        data = sm.datasets.macrodata.load()
        data = data.data[['realinv', 'realgdp', 'realcons']].view((float, 3))
        data = diff(log(data), axis=0)
        self.res1 = VAR2(endog=data).fit(maxlag=2)
        from results import results_var
        self.res2 = results_var.MacrodataResults()


class RVAR(object):
    """
    Estimates VAR model using R vars package and rpy
    """

    def __init__(self, data, p=1, type='both'):
        self.rdata = data
        self.p = p
        self.type = type

        self.pydata = rpy.convert_robj(data)
        self._estimate = None
        self.estimate()

    @property
    def aic(self):
        pass

    @property
    def bic(self):
        pass

    @property
    def beta(self):
        return rpy.convert_robj(r.coef(self._estimate))

    def summary(self, equation=None):
        print(r.summary(self._estimate, equation=equation))

    def output(self):
        print(self._estimate)

    def estimate(self):
        self._estimate = r.VAR(self.rdata, p=self.p, type=self.type)

    def plot(self, names=None):
        r.plot(model._estimate, names=names)

    def serial_test(self, lags_pt=16, type='PT.asymptotic'):
        f = r['serial.test']

        test = f(self._estimate, **{'lags.pt': lags_pt,
                                    'type': type})

        return test

    def data_summary(self):
        print(r.summary(self.rdata))


class TestVAR(TestCase):

    def setUp(self):
        try:
            import rpy2
        except ImportError:
            raise nose.SkipTest("No rpy2")

        self.rdata = rpy.load_data('Canada', package='vars', convert=False)
        self.data = rpy.load_data('Canada', package='vars', convert=True)

        self.res = VAR(self.data)
        self.ref = RVAR(self.rdata)

    def test_foo(self):
        pass

if __name__ == '__main__':
    # canada = rpy.load_data('Canada', package='vars', convert=False)

    # model = RVAR(canada, p=1)

    # summary(Canada)

    # plot(Canada, nc=2, xlab="")ppp

    # adf1 <- summary(ur.df(Canada[, "prod"], type = "trend", lags = 2))
    # adf1

    # adf2 <- summary(ur.df(diff(Canada[, "prod"]), type = "drift", lags = 1))
    # adf2

    # VARselect(Canada, lag.max = 8, type = "both")

    # Canada <- Canada[, c("prod", "e", "U", "rw")]

    # p1ct <- VAR(Canada, p = 1, type = "both")
    # p1ct

    # coefs <- coef(p1ct)
    # class(coefs)

    # run_module_suite()
    unittest.main()
