"""
Compute common portfolio statistics
"""

import pandas._tseries as _tseries
import pandas.stats.moments as moments

def leverage(weights):
    """
    Parameters
    ----------

    Returns
    -------
    y : TimeSeries
    """
    pass

def turnover(weights):
    """

    Returns
    -------
    y : TimeSeries
    """
    pass

def max_drawdown(returns):
    """
    Parameters
    ----------
    returns : TimeSeries or DataFrame

    """
    pass

def sharpe_ratio(returns):
    """

    Parameters
    ----------
    returns : TimeSeries or DataFrame

    Returns
    -------
    y : Series
    """
    pass

def rolling_sharpe_ratio(returns, window, min_periods=None):
    """

    Parameters
    ----------
    returns : TimeSeries or DataFrame


    Returns
    -------
    y : TimeSeries or DataFrame
    """
    pass

def beta(returns, market_returns):
    """

    Parameters
    ----------
    returns : TimeSeries or DataFrame


    Returns
    -------
    y : TimeSeries or DataFrame
    """
    pass
