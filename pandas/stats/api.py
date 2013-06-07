"""
Common namespace of statistical functions
"""

# pylint: disable-msg=W0611,W0614,W0401

from pandas.stats.moments import *
from pandas.stats.interface import ols
from pandas.stats.fama_macbeth import fama_macbeth

__all__ = ['ols', 'fama_macbeth', 'rolling_count', 'rolling_max', 'rolling_min',
        'rolling_sum', 'rolling_mean', 'rolling_std', 'rolling_cov', 'rolling_corr',
        'rolling_var', 'rolling_skew', 'rolling_kurt', 'rolling_quantile',
        'rolling_median', 'rolling_apply', 'rolling_corr_pairwise', 'rolling_window',
        'ewma', 'ewmvar', 'ewmstd', 'ewmvol', 'ewmcorr', 'ewmcov', 'expanding_count',
        'expanding_max', 'expanding_min', 'expanding_sum', 'expanding_mean',
        'expanding_std', 'expanding_cov', 'expanding_corr', 'expanding_var',
        'expanding_skew', 'expanding_kurt', 'expanding_quantile', 'expanding_median',
        'expanding_apply', 'expanding_corr_pairwise']
