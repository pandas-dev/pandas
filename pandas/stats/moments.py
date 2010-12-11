"""
Provides rolling statistical moments and related descriptive
statistics implemented in Cython
"""
from __future__ import division

from functools import wraps

from numpy import NaN
import numpy as np

from pandas.core.api import (DataFrame, DataMatrix, Series, notnull)
import pandas.lib.tseries as tseries

__all__ = ['rolling_count', 'rolling_max', 'rolling_min',
           'rolling_sum', 'rolling_mean', 'rolling_std', 'rolling_cov',
           'rolling_corr', 'rolling_var', 'rolling_skew', 'rolling_kurt',
           'rolling_median', 'ewma', 'ewmvar', 'ewmstd', 'ewmvol',
           'ewmcorr', 'ewmcov']

def rolling_count(arg, window, time_rule=None):
    """
    Rolling count of number of non-NaN observations inside provided window.

    Parameters
    ----------
    arg :  DataFrame or numpy ndarray-like
    window : Number of observations used for calculating statistic
    """
    arg = _conv_timerule(arg, time_rule)
    window = min(window, len(arg))

    return_hook, values = _process_data_structure(arg, kill_inf=False)

    converted = np.isfinite(values).astype(float)
    result = rolling_sum(converted, window, min_periods=1,
                         time_rule=time_rule)

    # putmask here?
    result[np.isnan(result)] = 0

    return return_hook(result)

def rolling_cov(arg1, arg2, window, min_periods=None, time_rule=None):
    X, Y = _prep_binary(arg1, arg2)
    mean = lambda x: rolling_mean(x, window, min_periods, time_rule)
    bias_adj = window / (window - 1)
    return (mean(X * Y) - mean(X) * mean(Y)) * bias_adj

def rolling_corr(arg1, arg2, window, min_periods=None, time_rule=None):
    X, Y = _prep_binary(arg1, arg2)
    num = rolling_cov(X, Y, window, min_periods, time_rule)
    den  = (rolling_std(X, window, min_periods, time_rule) *
            rolling_std(Y, window, min_periods, time_rule))
    return num / den

def _rolling_moment(arg, window, func, minp, axis=0, time_rule=None):
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
    axis : int, default 0
    time_rule : string or DateOffset
        Time rule to conform to before computing result

    Returns
    -------
    y : type of input
    """
    arg = _conv_timerule(arg, time_rule)
    calc = lambda x: func(x, window, minp=minp)
    return_hook, values = _process_data_structure(arg)
    # actually calculate the moment. Faster way to do this?
    result = np.apply_along_axis(calc, axis, values)

    return return_hook(result)

def _process_data_structure(arg, kill_inf=True):
    if isinstance(arg, DataFrame):
        if isinstance(arg, DataMatrix):
            return_hook = lambda v: DataMatrix(v, index=arg.index,
                                               columns=arg.columns,
                                               objects=arg.objects)
        else:
            return_hook = lambda v: DataFrame(v, index=arg.index,
                                              columns=arg.columns)
        values = arg.values
    elif isinstance(arg, Series):
        values = arg.values
        return_hook = lambda v: Series(v, arg.index)
    else:
        return_hook = lambda v: v
        values = arg

    if not issubclass(values.dtype.type, float):
        values = values.astype(float)

    if kill_inf:
        values = values.copy()
        values[np.isinf(values)] = np.NaN

    return return_hook, values

#-------------------------------------------------------------------------------
# Exponential moving moments

def _get_center_of_mass(com, span):
    if span is not None:
        if com is not None:
            raise Exception("com and span are mutually exclusive")

        # convert span to center of mass
        com = (span - 1) / 2.

    elif com is None:
        raise Exception("Must pass either com or span")

    return float(com)


def ewma(arg, com=None, span=None, min_periods=0, time_rule=None):
    com = _get_center_of_mass(com, span)
    arg = _conv_timerule(arg, time_rule)

    def _ewma(v):
        result = tseries.ewma(v, com)
        first_index = _first_valid_index(v)
        result[first_index : first_index + min_periods] = NaN
        return result

    return_hook, values = _process_data_structure(arg)
    output = np.apply_along_axis(_ewma, 0, values)
    return return_hook(output)

def _first_valid_index(arr):
    # argmax scans from left
    return notnull(arr).argmax()

def ewmvar(arg, com=None, span=None, min_periods=0, bias=False,
           time_rule=None):
    com = _get_center_of_mass(com, span)
    arg = _conv_timerule(arg, time_rule)
    moment2nd = ewma(arg * arg, com=com, min_periods=min_periods)
    moment1st = ewma(arg, com=com, min_periods=min_periods)

    result = moment2nd - moment1st ** 2
    if not bias:
        result *= (1.0 + 2.0 * com) / (2.0 * com)

    return result

def ewmstd(arg, com=None, span=None, min_periods=0, bias=False,
           time_rule=None):
    result = ewmvar(arg, com=com, span=span, time_rule=time_rule,
                    min_periods=min_periods, bias=bias)
    return np.sqrt(result)

ewmvol = ewmstd

def ewmcov(arg1, arg2, com=None, span=None, min_periods=0, bias=False,
           time_rule=None):
    X, Y = _prep_binary(arg1, arg2)

    X = _conv_timerule(X, time_rule)
    Y = _conv_timerule(Y, time_rule)

    mean = lambda x: ewma(x, com=com, span=span, min_periods=min_periods)

    result = (mean(X*Y) - mean(X) * mean(Y))

    if not bias:
        result *= (1.0 + 2.0 * com) / (2.0 * com)

    return result

def ewmcorr(arg1, arg2, com=None, span=None, min_periods=0,
            time_rule=None):
    X, Y = _prep_binary(arg1, arg2)

    X = _conv_timerule(X, time_rule)
    Y = _conv_timerule(Y, time_rule)

    mean = lambda x: ewma(x, com=com, span=span, min_periods=min_periods)
    var = lambda x: ewmvar(x, com=com, span=span, min_periods=min_periods,
                           bias=True)
    return (mean(X*Y) - mean(X)*mean(Y)) / np.sqrt(var(X) * var(Y))

def _prep_binary(arg1, arg2):
    if not isinstance(arg2, type(arg1)):
        raise Exception('Input arrays must be of the same type!')

    # mask out values, this also makes a common index...
    X = arg1 + 0 * arg2
    Y = arg2 + 0 * arg1

    return X, Y

#-------------------------------------------------------------------------------
# Docs

_doc_template = """
%s

Parameters
----------
%s
window : Number of observations used for calculating statistic
min_periods : int
    Minimum number of observations in window required to have a value
time_rule : {None, 'WEEKDAY', 'EOM', 'W@MON', ...}, default=None
    Name of time rule to conform to before computing statistic

Returns
-------
y : type of input argument
"""


_ewm_doc = r"""%s

Parameters
----------
%s
com : float. optional
    Center of mass: \alpha = com / (1 + com),
span : float, optional
    Specify decay in terms of span, \alpha = 2 / (span + 1)
min_periods : int, default 0
    Number of observations in sample to require (only affects
    beginning)
time_rule : {None, 'WEEKDAY', 'EOM', 'W@MON', ...}, default None
    Name of time rule to conform to before computing statistic
%s
Notes
-----
Either center of mass or span must be specified

EWMA is sometimes specified using a "span" parameter s, we have have that the
decay parameter \alpha is related to the span as :math:`\alpha = 1 - 2 / (s + 1)
= c / (1 + c)`

where c is the center of mass. Given a span, the associated center of mass is
:math:`c = (s - 1) / 2`

So a "20-day EWMA" would have center 9.5.

Returns
-------
y : type of input argument
"""

_unary_arg = "arg : Series, DataFrame, or DataMatrix"
_binary_arg = """arg1 : Series, DataFrame, or DataMatrix, or ndarray
arg2 : type of arg1"""

_bias_doc = r"""bias : boolean, default False
    Use a standard estimation bias correction
"""

rolling_cov.__doc__ = _doc_template % ("Unbiased moving covariance",
                                       _binary_arg)
rolling_corr.__doc__ = _doc_template % ("Moving sample correlation",
                                        _binary_arg)

ewma.__doc__ = _ewm_doc % ("Exponentially-weighted moving average",
                           _unary_arg, "")
ewmstd.__doc__ = _ewm_doc % ("Exponentially-weighted moving std",
                             _unary_arg, _bias_doc)
ewmvar.__doc__ = _ewm_doc % ("Exponentially-weighted moving variance",
                             _unary_arg, _bias_doc)
ewmcorr.__doc__ = _ewm_doc % ("Exponentially-weighted moving "
                              "correlation", _binary_arg, "")
ewmcov.__doc__ = _ewm_doc % ("Exponentially-weighted moving covariance",
                             _binary_arg, "")

#-------------------------------------------------------------------------------
# Python interface to Cython functions

def _conv_timerule(arg, time_rule):
    types = (DataFrame, Series)
    if time_rule is not None and isinstance(arg, types):
        # Conform to whatever frequency needed.
        arg = arg.asfreq(time_rule)

    return arg

def _two_periods(minp, window):
    if minp is None:
        return window
    else:
        return max(2, minp)

def _use_window(minp, window):
    if minp is None:
        return window
    else:
        return minp

def _rolling_func(func, desc, check_minp=_use_window):
    @wraps(func)
    def f(arg, window, min_periods=None, time_rule=None):
        def call_cython(arg, window, minp):
            minp = check_minp(minp, window)
            return func(arg, window, minp)
        return _rolling_moment(arg, window, call_cython, min_periods,
                               time_rule=time_rule)

    f.__doc__ = _doc_template % (desc, _unary_arg)

    return f

rolling_max = _rolling_func(tseries.roll_max, 'Moving maximum')
rolling_min = _rolling_func(tseries.roll_min, 'Moving minimum')
rolling_sum = _rolling_func(tseries.roll_sum, 'Moving sum')
rolling_mean = _rolling_func(tseries.roll_mean, 'Moving mean')
rolling_median = _rolling_func(tseries.roll_median, 'Moving median')

_ts_std = lambda *a, **kw: np.sqrt(tseries.roll_var(*a, **kw))
rolling_std = _rolling_func(_ts_std, 'Unbiased moving standard deviation',
                            check_minp=_two_periods)
rolling_var = _rolling_func(tseries.roll_var, 'Unbiased moving variance',
                            check_minp=_two_periods)
rolling_skew = _rolling_func(tseries.roll_skew, 'Unbiased moving skewness',
                             check_minp=_two_periods)
rolling_kurt = _rolling_func(tseries.roll_kurt, 'Unbiased moving kurtosis',
                             check_minp=_two_periods)

