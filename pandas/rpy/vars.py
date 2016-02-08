import pandas.rpy.util as util


class VAR(object):
    """

    Parameters
    ----------
    y :
    p :
    type : {"const", "trend", "both", "none"}
    season :
    exogen :
    lag_max :
    ic : {"AIC", "HQ", "SC", "FPE"}
        Information criterion to use, if lag_max is not None
    """
    def __init__(y, p=1, type="none", season=None, exogen=None,
                 lag_max=None, ic=None):
        pass
