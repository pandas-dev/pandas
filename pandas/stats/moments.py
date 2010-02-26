"""
Provides rolling statistical moments and related descriptive
statistics implemented in Cython
"""

from numpy import NaN
from pandas.core.api import (DataFrame, TimeSeries, DataMatrix,
                                 Series, notnull)

import pandas.lib.tseries as tseries
import numpy as np

__all__ = ['rolling_count', 'rolling_sum', 'rolling_mean',
           'rolling_std', 'rolling_cov', 'rolling_corr',
           'rolling_var', 'rolling_skew', 'rolling_kurt',
           'rolling_median', 'ewma', 'ewmvol', 'ewmcorr', 'ewmcov']

#-------------------------------------------------------------------------------
# Rolling statistics

def rolling_count(arg, window, timeRule=None):
    """
    Rolling count of number of observations inside provided window.

    Parameters
    ----------
    arg :  DataFrame or numpy ndarray-like
    window : Number of observations used for calculating statistic

    Note
    ----
    Uses Cython
    """
    types = (DataFrame, DataMatrix, TimeSeries)
    if timeRule is not None and isinstance(arg, types):
        # Conform to whatever frequency needed.
        arg = arg.asfreq(timeRule)

    window = min(window, len(arg))
    if isinstance(arg, DataMatrix):
        arg = arg.copy()
        arg.values = np.isfinite(arg.values).astype(float)
        result = rolling_sum(arg, window, minPeriods=1, timeRule=timeRule)
        result.values[np.isnan(result.values)] = 0
    elif isinstance(arg, DataFrame):
        converter = lambda x: np.isfinite(x).astype(float)
        arg = arg.apply(converter)
        result = rolling_sum(arg, window, minPeriods=1, timeRule=timeRule)
        result = result.fill(value=0)
    else:
        arg = np.isfinite(arg).astype(float)
        result = rolling_sum(arg, window, minPeriods=1, timeRule=timeRule)
        result[np.isnan(result)] = 0
    return result

def rolling_sum(arg, window=0, minPeriods=None, timeRule=None):
    """
    Rolling sum of observations inside provided window.

    Parameters
    ----------
    arg :  DataFrame or numpy ndarray-like
    window : Number of observations used for calculating statistic
    """
    return _rollingMoment(arg, window, tseries.rolling_sum, minp=minPeriods,
                          timeRule=timeRule)

def _rollingMoment(arg, window, func, minp=None, timeRule=None):
    """
    Rolling statistical measure using supplied function. Designed to be
    used with passed-in Cython array-based functions.

    Parameters
    ----------
    arg :  DataFrame or numpy ndarray-like
    window : Number of observations used for calculating statistic
    func : Cython function to compute rolling statistic on raw series
    minp : int
        Minimum number of observations required to have a value
    """
    if minp is None:
        minp = window

    types = (DataFrame, DataMatrix, TimeSeries)
    if timeRule is not None and isinstance(arg, types):
        # Conform to whatever frequency needed.
        arg = arg.asfreq(timeRule)

    if isinstance(arg, DataMatrix):
        T, N = arg.values.shape
        resultMatrix = np.empty((T, N), dtype=arg.values.dtype)
        arg.values[np.isinf(arg.values)] = NaN
        for i in range(N):
            resultMatrix[:,i] = func(arg.values[:, i], window, minp=minp)
        output = DataMatrix(resultMatrix, index=arg.index,
                            columns=arg.columns)

    elif isinstance(arg, DataFrame):
        output = DataFrame(index = arg.index)
        for col, series in arg.iteritems():
            series[np.isinf(series)] = NaN
            output[col] = TimeSeries(func(series,window,minp=minp),
                                     index = series.index)
    elif isinstance(arg, TimeSeries):
        arg[np.isinf(arg)] = NaN
        output = TimeSeries(func(arg,window,minp=minp), index = arg.index)
    else:
        try:
            assert(hasattr(arg, '__iter__'))
        except AssertionError:
            raise AssertionError('Expected DataFrame or array-like argument')
        arg[np.isinf(arg)] = NaN
        output = func(arg, window, minp=minp)
    return output

def rolling_cov(arg1, arg2, window, minPeriods=None, timeRule=None):
    """
    Unbiased Rolling Covariance

    Parameters
    ----------
    arg :  DataFrame or numpy ndarray-like
    window : Number of observations used for calculating statistic
    minPeriods : int
        Minimum number of observations in window required to have a value
    timeRule : {None, 'WEEKDAY', 'EOM', 'W@MON', ...}, default=None
        Name of time rule to conform to before computing statistic
    """
    num1 = rolling_mean(arg1*arg2, window, minPeriods, timeRule) #E(XY)
    num2 = (rolling_mean(arg1, window, minPeriods, timeRule) *
            rolling_mean(arg2, window, minPeriods, timeRule)) #E(X)E(Y)
    return (num1-num2)*window/(window-1)

def rolling_corr(arg1, arg2, window, minPeriods=None, timeRule=None):
    """
    Rolling Correl

    Parameters
    ----------
    arg :  DataFrame or numpy ndarray-like
    window : Number of observations used for calculating statistic
    minPeriods : int
        Minimum number of observations in window required to have a value
    timeRule : {None, 'WEEKDAY', 'EOM', 'W@MON', ...}, default=None
        Name of time rule to conform to before computing statistic
    """
    num = rolling_cov(arg1, arg2, window, minPeriods, timeRule) #Cov(X,Y)
    den  = (rolling_std(arg1, window, minPeriods, timeRule) *
            rolling_std(arg2, window, minPeriods, timeRule)) #STDEV(X)STDEV(Y)
    return (num/den)

def rolling_mean(arg, window, minPeriods=None, timeRule=None):
    """
    Unbiased Rolling Mean

    Parameters
    ----------
    arg :  DataFrame or numpy ndarray-like
    window : Number of observations used for calculating statistic
    minPeriods : int
        Minimum number of observations in window required to have a value
    timeRule : {None, 'WEEKDAY', 'EOM', 'W@MON', ...}, default=None
        Name of time rule to conform to before computing statistic
    """
    return _rollingMoment(arg, window, tseries.rolling_mean,
                          minp=minPeriods, timeRule=timeRule)

def rolling_std(arg, window, minPeriods=None, timeRule=None):
    """
    Unbiased Rolling Volatility

    Parameters
    ----------
    arg :  DataFrame or numpy ndarray-like
    window : Number of observations used for calculating statistic
    minPeriods : int
        Minimum number of observations in window required to have a value
    timeRule : {None, 'WEEKDAY', 'EOM', 'W@MON', ...}, default=None
        Name of time rule to conform to before computing statistic
    """
    return _rollingMoment(arg, window, tseries.rolling_std,
                          minp=minPeriods, timeRule=timeRule)

def rolling_var(arg, window, minPeriods=None, timeRule=None):
    """
    Unbiased Rolling Variance

    Parameters
    ----------
    arg :  DataFrame or numpy ndarray-like
    window : Number of observations used for calculating statistic
    minPeriods : int
        Minimum number of observations in window required to have a value
    timeRule : {None, 'WEEKDAY', 'EOM', 'W@MON', ...}, default=None
        Name of time rule to conform to before computing statistic
    """
    return _rollingMoment(arg, window, tseries.rolling_var,
                          minp=minPeriods, timeRule=timeRule)

def rolling_skew(arg, window, minPeriods=None, timeRule=None):
    """
    Unbiased Rolling Skewness

    Parameters
    ----------
    arg :  DataFrame or numpy ndarray-like
    window : Number of observations used for calculating statistic
    minPeriods : int
        Minimum number of observations in window required to have a value
    timeRule : {None, 'WEEKDAY', 'EOM', 'W@MON', ...}, default=None
        Name of time rule to conform to before computing statistic
    """
    return _rollingMoment(arg, window, tseries.rolling_skew,
                          minp=minPeriods, timeRule=timeRule)

def rolling_kurt(arg, window, minPeriods=None, timeRule=None):
    """
    Unbiased Rolling Kurtosis

    Parameters
    ----------
    arg :  DataFrame or numpy ndarray-like
    window : Number of observations used for calculating statistic
    minPeriods : int
        Minimum number of observations in window required to have a value
    timeRule : {None, 'WEEKDAY', 'EOM', 'W@MON', ...}, default=None
        Name of time rule to conform to before computing statistic
    """
    return _rollingMoment(arg, window, tseries.rolling_kurt,
                          minp=minPeriods, timeRule=timeRule)

def rolling_median(arg, window, minPeriods=None, timeRule=None):
    return _rollingMoment(arg, window, tseries.rolling_median,
                          minp=minPeriods, timeRule=timeRule)

def test_rolling_median():
    arr = np.random.randn(100)
    arr[20:40] = np.NaN

    result = rolling_median(arr, 50)

    assert(np.isnan(result[20]))

    assert(result[-1] == np.median(arr[-50:]))

    result = rolling_median(arr, 49)

    assert(np.isnan(result[20]))

    assert(result[-1] == np.median(arr[-49:]))

#-------------------------------------------------------------------------------
# Exponential moving moments

def _getMinPeriods(minPct, rho):
    if not (0 <= minPct <= 1):
        raise Exception('minPct must be between 0 and 1!')

    if minPct == 0:
        return 0

    return int(np.ceil(np.log(minPct) / np.log(rho)))

def _ewmoment(values, func, minPeriods=None, minPct=0.95,
              biasCorrection=None):
    """
    Generic rolling exponential moment function using blended accumulator
    method.

    Parameters
    ----------
    values : ndarray or Series
    func : function
        taking previous value and next value

    biasCorrection : float
        Optional bias correction

    ** Mutually exclusive options **

    minPct : float
        Minimum percentage of weight "in window" to require to have a value
    minPeriods : int, optional
        require a particular number of periods "in window" to compute statistic
        If provided, overrides the minPct argument

    Returns
    -------
    Same type and length as values argument
    """
    isSeries = isinstance(values, Series)
    okLocs = notnull(values)

    if isSeries:
        index = values.index

    cleanValues = values[okLocs]

    result = np.frompyfunc(func, 2, 1).accumulate(cleanValues)
    result = result.astype(float)

    if minPeriods is not None:
        if minPeriods < 0:
            raise Exception('minPeriods cannot be less than 0!')

        result[:minPeriods] = np.NaN

    output = values.copy()
    output[okLocs] = result

    if biasCorrection is not None:
        if biasCorrection <= 0:
            raise Exception('Bias correction cannot be negative!')

        output *= biasCorrection

    return output

def ewma(arg, com, minCom=0):
    """
    Calculates the rolling exponentially weighted moving average of a series.

    Parameters
    ----------
    arg : Series, DataFrame, or DataMatrix
    com : integer
        Center of Mass for exponentially weighted moving average
        decay = com / (1 + com) maps center of mass to decay parameter

    minCom : int, default 0
        Optionally require that at least a certain number of periods as
        a multiple of the Center of Mass be included in the sample.
    """
    def ewmaFunc(series):
        series[np.isinf(series)] = NaN
        result = tseries.ewma(series, com)

        firstIndex = np.arange(len(series))[notnull(series)][0]

        result[firstIndex : firstIndex + minCom*com] = NaN
        result = Series(result, index=arg.index)
        return result

    if isinstance(arg, TimeSeries):
        output = ewmaFunc(arg)
    elif isinstance(arg, DataFrame):
        output = arg.apply(ewmaFunc)
    else:
        output = tseries.ewma(arg, com)
        firstIndex = np.arange(len(arg))[notnull(arg)][0]
        output[firstIndex : firstIndex + minCom * com] = NaN

    return output

def ewmvar(arg, com, minCom = 0, correctBias = True):
    """
    Calculates the rolling exponentially weighted moving variance of a series.

    Parameters
    ----------
    series : Series
    com : integer
        Center of Mass for exponentially weighted moving average
        decay = com / (1 + com) maps center of mass to decay parameter

    minCom : int, default None
        Optionally require that at least a certain number of periods as
        a multiple of the Center of Mass be included in the sample.

    correctBias : boolean
        Use a standard bias correction
    """

    if correctBias:
        biasCorrection = ( 1.0 + 2.0 * com ) / (2.0 * com)
    else:
        biasCorrection = 1.0

    moment2nd = ewma(arg * arg, com=com, minCom=minCom)
    moment1st = ewma(arg, com=com, minCom=minCom)

    return biasCorrection * (moment2nd - moment1st**2)

def ewmvol(arg, com, minCom = 0, correctBias = True):
    """
    Calculates the rolling exponentially weighted moving variance of a series.

    Parameters
    ----------
    series : Series
    com : integer
        Center of Mass for exponentially weighted moving average
        decay = com / (1 + com) maps center of mass to decay parameter

    minCom : int, default None
        Optionally require that at least a certain number of periods as
        a multiple of the Center of Mass be included in the sample.

    correctBias : boolean
        Use a standard bias correction
    """
    result = ewmvar(arg, com=com, minCom=minCom, correctBias=correctBias)

    if isinstance(result, DataFrame):
        result = result.apply(np.sqrt)
    else:
        result = np.sqrt(result)

    return result

def ewmcov(seriesA, seriesB, com, minCom = 0, correctBias = True):
    """
    Calculates the rolling exponentially weighted moving variance of a
    series.

    Parameters
    ----------
    series : Series
    com : integer
        Center of Mass for exponentially weighted moving average
        decay = com / (1 + com) maps center of mass to decay parameter

    minCom : int, default None
        Optionally require that at least a certain number of periods as
        a multiple of the Center of Mass be included in the sample.

    correctBias : boolean
        Use a standard bias correction
    """

    if correctBias:
        biasCorrection = ( 1.0 + 2.0 * com ) / (2.0 * com)
    else:
        biasCorrection = 1.0

    if not isinstance(seriesB, type(seriesA)):
        raise Exception('Input arrays must be of the same type!')

    if isinstance(seriesA, Series):
        if seriesA.index is not seriesB.index:
            commonIndex = seriesA.index.intersection(seriesB.index)
            seriesA = seriesA.reindex(commonIndex)
            seriesB = seriesB.reindex(commonIndex)

    okLocs = notnull(seriesA) & notnull(seriesB)

    cleanSeriesA = seriesA[okLocs]
    cleanSeriesB = seriesB.reindex(cleanSeriesA.index)

    XY = ewma(cleanSeriesA * cleanSeriesB, com = com, minCom = minCom)
    X  = ewma(cleanSeriesA, com = com, minCom = minCom)
    Y  = ewma(cleanSeriesB, com = com, minCom = minCom)

    return biasCorrection * (XY - X * Y)


def ewmcorr(seriesA, seriesB, com, minCom = 0):
    """
    Calculates a rolling exponentially weighted moving correlation of
    2 series.

    Parameters
    ----------
    seriesA : Series
    seriesB : Series
    com : integer
        Center of Mass for exponentially weighted moving average
        decay = com / (1 + com) maps center of mass to decay parameter

    minCom : int, default None
        Optionally require that at least a certain number of periods as
        a multiple of the Center of Mass be included in the sample.
    """

    if not isinstance(seriesB, type(seriesA)):
        raise Exception('Input arrays must be of the same type!')

    if isinstance(seriesA, Series):
        if seriesA.index is not seriesB.index:
            commonIndex = seriesA.index.intersection(seriesB.index)
            seriesA = seriesA.reindex(commonIndex)
            seriesB = seriesB.reindex(commonIndex)

    okLocs = notnull(seriesA) & notnull(seriesB)

    cleanSeriesA = seriesA[okLocs]
    cleanSeriesB = seriesB.reindex(cleanSeriesA.index)

    XY = ewma(cleanSeriesA * cleanSeriesB, com = com, minCom = minCom)
    X  = ewma(cleanSeriesA, com = com, minCom = minCom)
    Y  = ewma(cleanSeriesB, com = com, minCom = minCom)
    varX = ewmvar(cleanSeriesA, com = com, minCom = minCom, correctBias = False)
    varY = ewmvar(cleanSeriesB, com = com, minCom = minCom, correctBias = False)


    return (XY - X * Y) / np.sqrt( varX * varY )

