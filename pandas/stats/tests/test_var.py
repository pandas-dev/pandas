from rpy2.robjects import r
from rpy2.robjects.packages import importr

vars = importr('vars')
urca = importr('urca')

class RVAR(object):
    pass

if __name__ == '__main__':
    r.data("Canada")
    can = r['Canada']
    p1ct = r('p1ct <- VAR(Canada, p=1, type="both")')
    coefs = r('coefs <- coef(p1ct)')

    ecoef = coefs.rx2('e')

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

