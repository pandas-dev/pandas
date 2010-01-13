"""
Common namespace of statistical functions
"""

# pylint: disable-msg=W0611,W0614

from pandas.stats.moments import *

try:
    import scipy.stats as _stats

    from pandas.stats.interface import ols
    from pandas.stats.fama_macbeth import fama_macbeth

    del _stats
except ImportError:
    def ols(*args, **kwargs):
        """
        Stub to give error message for missing SciPy
        """
        raise Exception('Must install SciPy to use OLS functionality')
