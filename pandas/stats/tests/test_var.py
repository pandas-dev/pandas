
import rpy2.robjects as robj
from rpy2.robjects import r
from rpy2.robjects.packages import importr
import pandas.rpy.common as rpy

import scikits.statsmodels.tsa.var as sm_var
import pandas.stats.var as pvar

vars = importr('vars')
urca = importr('urca')

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

    def summary(self, equation=None):
        print r.summary(self._estimate, equation=equation)

    def output(self):
        print self._estimate

    def estimate(self):
        self._estimate = r.VAR(self.rdata, p=self.p, type=self.type)

    def plot(self, names=None):
        r.plot(model._estimate, names=names)

    def serial_test(self, lags_pt=16, type='PT.asymptotic'):
        f = r['serial.test']

        test = f(self._estimate, **{'lags.pt' : lags_pt,
                                    'type' : type})

        return test

    @property
    def beta(self):
        return rpy.convert_robj(r.coef(self._estimate))

    def data_summary(self):
        print r.summary(self.rdata)

if __name__ == '__main__':
    canada = rpy.load_data('Canada', package='vars', convert=False)

    model = RVAR(canada, p=1)

    # summary(Canada)

    # plot(Canada, nc=2, xlab="")

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

