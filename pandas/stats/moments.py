"""
Provides rolling statistical moments and related descriptive
statistics implemented in Cython
"""
from __future__ import division

from functools import wraps
from collections import defaultdict

from numpy import NaN
import numpy as np

from pandas.core.api import DataFrame, Series, Panel, notnull
import pandas.algos as algos
import pandas.core.common as pdcom

from pandas.util.decorators import Substitution, Appender

__all__ = ['rolling_count', 'rolling_max', 'rolling_min',
           'rolling_sum', 'rolling_mean', 'rolling_std', 'rolling_cov',
           'rolling_corr', 'rolling_var', 'rolling_skew', 'rolling_kurt',
           'rolling_quantile', 'rolling_median', 'rolling_apply',
           'rolling_corr_pairwise', 'rolling_window',
           'ewma', 'ewmvar', 'ewmstd', 'ewmvol', 'ewmcorr', 'ewmcov',
           'expanding_count', 'expanding_max', 'expanding_min',
           'expanding_sum', 'expanding_mean', 'expanding_std',
           'expanding_cov', 'expanding_corr', 'expanding_var',
           'expanding_skew', 'expanding_kurt', 'expanding_quantile',
           'expanding_median', 'expanding_apply', 'expanding_corr_pairwise']

#------------------------------------------------------------------------------
# Docs

# The order of arguments for the _doc_template is:
# (header, args, kwargs, returns, notes)

_doc_template = """
%s

Parameters
----------
%s%s
Returns
-------
%s
%s
"""

_roll_kw = """window : int
    Size of the moving window. This is the number of observations used for
    calculating the statistic.
min_periods : int, default None
    Minimum number of observations in window required to have a value
    (otherwise result is NA).
freq : string or DateOffset object, optional (default None)
    Frequency to conform the data to before computing the statistic. Specified
    as a frequency string or DateOffset object.
center : boolean, default False
    Set the labels at the center of the window.
how : string, default '%s'
    Method for down- or re-sampling
"""

_roll_notes = r"""
Notes
-----
By default, the result is set to the right edge of the window. This can be
changed to the center of the window by setting ``center=True``.

The `freq` keyword is used to conform time series data to a specified
frequency by resampling the data. This is done with the default parameters
of :meth:`~pandas.Series.resample` (i.e. using the `mean`).
"""


_ewm_kw = r"""com : float. optional
    Center of mass: :math:`\alpha = 1 / (1 + com)`,
span : float, optional
    Specify decay in terms of span, :math:`\alpha = 2 / (span + 1)`
halflife : float, optional
    Specify decay in terms of halflife, :math:`\alpha = 1 - exp(log(0.5) / halflife)`
min_periods : int, default 0
    Number of observations in sample to require (only affects
    beginning)
freq : None or string alias / date offset object, default=None
    Frequency to conform to before computing statistic
adjust : boolean, default True
    Divide by decaying adjustment factor in beginning periods to account for
    imbalance in relative weightings (viewing EWMA as a moving average)
how : string, default 'mean'
    Method for down- or re-sampling
"""

_ewm_notes = r"""
Notes
-----
Either center of mass or span must be specified

EWMA is sometimes specified using a "span" parameter `s`, we have that the
decay parameter :math:`\alpha` is related to the span as
:math:`\alpha = 2 / (s + 1) = 1 / (1 + c)`

where `c` is the center of mass. Given a span, the associated center of mass is
:math:`c = (s - 1) / 2`

So a "20-day EWMA" would have center 9.5.
"""

_expanding_kw = """min_periods : int, default None
    Minimum number of observations in window required to have a value
    (otherwise result is NA).
freq : string or DateOffset object, optional (default None)
    Frequency to conform the data to before computing the statistic. Specified
    as a frequency string or DateOffset object.
"""


_type_of_input_retval = "y : type of input argument"

_flex_retval = """y : type depends on inputs
    DataFrame / DataFrame -> DataFrame (matches on columns) or Panel (pairwise)
    DataFrame / Series -> Computes result for each column
    Series / Series -> Series"""

_pairwise_retval = "y : Panel whose items are df1.index values"

_unary_arg = "arg : Series, DataFrame\n"

_binary_arg_flex = """arg1 : Series, DataFrame, or ndarray
arg2 : Series, DataFrame, or ndarray, optional
    if not supplied then will default to arg1 and produce pairwise output
"""

_binary_arg = """arg1 : Series, DataFrame, or ndarray
arg2 : Series, DataFrame, or ndarray
"""

_pairwise_arg = """df1 : DataFrame
df2 : DataFrame
"""

_pairwise_kw = """pairwise : bool, default False
    If False then only matching columns between arg1 and arg2 will be used and
    the output will be a DataFrame.
    If True then all pairwise combinations will be calculated and the output
    will be a Panel in the case of DataFrame inputs. In the case of missing
    elements, only complete pairwise observations will be used.
"""

_bias_kw = r"""bias : boolean, default False
    Use a standard estimation bias correction
"""


def rolling_count(arg, window, freq=None, center=False, how=None):
    """
    Rolling count of number of non-NaN observations inside provided window.

    Parameters
    ----------
    arg :  DataFrame or numpy ndarray-like
    window : int
        Size of the moving window. This is the number of observations used for
        calculating the statistic.
    freq : string or DateOffset object, optional (default None)
        Frequency to conform the data to before computing the statistic. Specified
        as a frequency string or DateOffset object.
    center : boolean, default False
        Whether the label should correspond with center of window
    how : string, default 'mean'
        Method for down- or re-sampling

    Returns
    -------
    rolling_count : type of caller

    Notes
    -----
    The `freq` keyword is used to conform time series data to a specified
    frequency by resampling the data. This is done with the default parameters
    of :meth:`~pandas.Series.resample` (i.e. using the `mean`).
    """
    arg = _conv_timerule(arg, freq, how)
    window = min(window, len(arg))

    return_hook, values = _process_data_structure(arg, kill_inf=False)

    converted = np.isfinite(values).astype(float)
    result = rolling_sum(converted, window, min_periods=1,
                         center=center)  # already converted

    # putmask here?
    result[np.isnan(result)] = 0

    return return_hook(result)


@Substitution("Unbiased moving covariance.", _binary_arg_flex,
              _roll_kw%'None'+_pairwise_kw, _flex_retval, _roll_notes)
@Appender(_doc_template)
def rolling_cov(arg1, arg2=None, window=None, min_periods=None, freq=None,
                center=False, pairwise=None, how=None):
    if window is None and isinstance(arg2, (int, float)):
        window = arg2
        arg2 = arg1
        pairwise = True if pairwise is None else pairwise  # only default unset
    elif arg2 is None:
        arg2 = arg1
        pairwise = True if pairwise is None else pairwise  # only default unset
    arg1 = _conv_timerule(arg1, freq, how)
    arg2 = _conv_timerule(arg2, freq, how)
    window = min(window, len(arg1), len(arg2))

    def _get_cov(X, Y):
        mean = lambda x: rolling_mean(x, window, min_periods, center=center)
        count = rolling_count(X + Y, window, center=center)
        bias_adj = count / (count - 1)
        return (mean(X * Y) - mean(X) * mean(Y)) * bias_adj
    rs = _flex_binary_moment(arg1, arg2, _get_cov, pairwise=bool(pairwise))
    return rs


@Substitution("Moving sample correlation.", _binary_arg_flex,
              _roll_kw%'None'+_pairwise_kw, _flex_retval, _roll_notes)
@Appender(_doc_template)
def rolling_corr(arg1, arg2=None, window=None, min_periods=None, freq=None,
                 center=False, pairwise=None, how=None):
    if window is None and isinstance(arg2, (int, float)):
        window = arg2
        arg2 = arg1
        pairwise = True if pairwise is None else pairwise  # only default unset
    elif arg2 is None:
        arg2 = arg1
        pairwise = True if pairwise is None else pairwise  # only default unset
    arg1 = _conv_timerule(arg1, freq, how)
    arg2 = _conv_timerule(arg2, freq, how)
    window = min(window, len(arg1), len(arg2))

    def _get_corr(a, b):
        num = rolling_cov(a, b, window, min_periods, freq=freq,
                          center=center)
        den = (rolling_std(a, window, min_periods, freq=freq,
                           center=center) *
               rolling_std(b, window, min_periods, freq=freq,
                           center=center))
        return num / den
    return _flex_binary_moment(arg1, arg2, _get_corr, pairwise=bool(pairwise))


def _flex_binary_moment(arg1, arg2, f, pairwise=False):
    if not (isinstance(arg1,(np.ndarray, Series, DataFrame)) and
            isinstance(arg2,(np.ndarray, Series, DataFrame))):
        raise TypeError("arguments to moment function must be of type "
                         "np.ndarray/Series/DataFrame")

    if isinstance(arg1, (np.ndarray, Series)) and \
            isinstance(arg2, (np.ndarray,Series)):
        X, Y = _prep_binary(arg1, arg2)
        return f(X, Y)
    elif isinstance(arg1, DataFrame):
        results = {}
        if isinstance(arg2, DataFrame):
            X, Y = arg1.align(arg2, join='outer')
            X = X + 0 * Y
            Y = Y + 0 * X
            if pairwise is False:
                res_columns = arg1.columns.union(arg2.columns)
                for col in res_columns:
                    if col in X and col in Y:
                        results[col] = f(X[col], Y[col])
            elif pairwise is True:
                results = defaultdict(dict)
                for i, k1 in enumerate(arg1.columns):
                    for j, k2 in enumerate(arg2.columns):
                        if j<i and arg2 is arg1:
                            # Symmetric case
                            results[k1][k2] = results[k2][k1]
                        else:
                            results[k1][k2] = f(arg1[k1], arg2[k2])
                return Panel.from_dict(results).swapaxes('items', 'major')
            else:
                raise ValueError("'pairwise' is not True/False")
        else:
            res_columns = arg1.columns
            X, Y = arg1.align(arg2, axis=0, join='outer')
            results = {}

            for col in res_columns:
                results[col] = f(X[col], Y)

        return DataFrame(results, index=X.index, columns=res_columns)
    else:
        return _flex_binary_moment(arg2, arg1, f)


@Substitution("Deprecated. Use rolling_corr(..., pairwise=True) instead.\n\n"
              "Pairwise moving sample correlation", _pairwise_arg,
              _roll_kw%'None', _pairwise_retval, _roll_notes)
@Appender(_doc_template)
def rolling_corr_pairwise(df1, df2=None, window=None, min_periods=None,
                          freq=None, center=False):
    import warnings
    warnings.warn("rolling_corr_pairwise is deprecated, use rolling_corr(..., pairwise=True)", FutureWarning)
    return rolling_corr(df1, df2, window=window, min_periods=min_periods,
                        freq=freq, center=center,
                        pairwise=True)


def _rolling_moment(arg, window, func, minp, axis=0, freq=None, center=False,
                    how=None, args=(), kwargs={}, **kwds):
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
    freq : None or string alias / date offset object, default=None
        Frequency to conform to before computing statistic
    center : boolean, default False
        Whether the label should correspond with center of window
    how : string, default 'mean'
        Method for down- or re-sampling
    args : tuple
        Passed on to func
    kwargs : dict
        Passed on to func

    Returns
    -------
    y : type of input
    """
    arg = _conv_timerule(arg, freq, how)
    calc = lambda x: func(x, window, minp=minp, args=args, kwargs=kwargs,
                          **kwds)
    return_hook, values = _process_data_structure(arg)
    # actually calculate the moment. Faster way to do this?
    if values.ndim > 1:
        result = np.apply_along_axis(calc, axis, values)
    else:
        result = calc(values)

    rs = return_hook(result)
    if center:
        rs = _center_window(rs, window, axis)
    return rs


def _center_window(rs, window, axis):
    if axis > rs.ndim-1:
        raise ValueError("Requested axis is larger then no. of argument "
                         "dimensions")

    offset = int((window - 1) / 2.)
    if isinstance(rs, (Series, DataFrame, Panel)):
        rs = rs.shift(-offset, axis=axis)
    else:
        rs_indexer = [slice(None)] * rs.ndim
        rs_indexer[axis] = slice(None, -offset)

        lead_indexer = [slice(None)] * rs.ndim
        lead_indexer[axis] = slice(offset, None)

        na_indexer = [slice(None)] * rs.ndim
        na_indexer[axis] = slice(-offset, None)

        rs[tuple(rs_indexer)] = np.copy(rs[tuple(lead_indexer)])
        rs[tuple(na_indexer)] = np.nan
    return rs


def _process_data_structure(arg, kill_inf=True):
    if isinstance(arg, DataFrame):
        return_hook = lambda v: type(arg)(v, index=arg.index,
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

#------------------------------------------------------------------------------
# Exponential moving moments


def _get_center_of_mass(com, span, halflife):
    valid_count = len([x for x in [com, span, halflife] if x is not None])
    if valid_count > 1:
        raise Exception("com, span, and halflife are mutually exclusive")

    if span is not None:
        # convert span to center of mass
        com = (span - 1) / 2.
    elif halflife is not None:
        # convert halflife to center of mass
        decay = 1 - np.exp(np.log(0.5) / halflife)
        com = 1 / decay - 1
    elif com is None:
        raise Exception("Must pass one of com, span, or halflife")

    return float(com)


@Substitution("Exponentially-weighted moving average", _unary_arg, _ewm_kw,
              _type_of_input_retval, _ewm_notes)
@Appender(_doc_template)
def ewma(arg, com=None, span=None, halflife=None, min_periods=0, freq=None,
         adjust=True, how=None):
    com = _get_center_of_mass(com, span, halflife)
    arg = _conv_timerule(arg, freq, how)

    def _ewma(v):
        result = algos.ewma(v, com, int(adjust))
        first_index = _first_valid_index(v)
        result[first_index: first_index + min_periods] = NaN
        return result

    return_hook, values = _process_data_structure(arg)
    output = np.apply_along_axis(_ewma, 0, values)
    return return_hook(output)


def _first_valid_index(arr):
    # argmax scans from left
    return notnull(arr).argmax() if len(arr) else 0


@Substitution("Exponentially-weighted moving variance", _unary_arg,
              _ewm_kw+_bias_kw, _type_of_input_retval, _ewm_notes)
@Appender(_doc_template)
def ewmvar(arg, com=None, span=None, halflife=None, min_periods=0, bias=False,
           freq=None, how=None):
    com = _get_center_of_mass(com, span, halflife)
    arg = _conv_timerule(arg, freq, how)
    moment2nd = ewma(arg * arg, com=com, min_periods=min_periods)
    moment1st = ewma(arg, com=com, min_periods=min_periods)

    result = moment2nd - moment1st ** 2
    if not bias:
        result *= (1.0 + 2.0 * com) / (2.0 * com)

    return result


@Substitution("Exponentially-weighted moving std", _unary_arg,
              _ewm_kw+_bias_kw, _type_of_input_retval, _ewm_notes)
@Appender(_doc_template)
def ewmstd(arg, com=None, span=None, halflife=None, min_periods=0, bias=False):
    result = ewmvar(arg, com=com, span=span, halflife=halflife,
                    min_periods=min_periods, bias=bias)
    return _zsqrt(result)

ewmvol = ewmstd


@Substitution("Exponentially-weighted moving covariance", _binary_arg_flex,
              _ewm_kw+_pairwise_kw, _type_of_input_retval, _ewm_notes)
@Appender(_doc_template)
def ewmcov(arg1, arg2=None, com=None, span=None, halflife=None, min_periods=0,
           bias=False, freq=None, pairwise=None, how=None):
    if arg2 is None:
        arg2 = arg1
        pairwise = True if pairwise is None else pairwise
    elif isinstance(arg2, (int, float)) and com is None:
        com = arg2
        arg2 = arg1
        pairwise = True if pairwise is None else pairwise
    arg1 = _conv_timerule(arg1, freq, how)
    arg2 = _conv_timerule(arg2, freq, how)

    def _get_ewmcov(X, Y):
        mean = lambda x: ewma(x, com=com, span=span, halflife=halflife, min_periods=min_periods)
        return (mean(X * Y) - mean(X) * mean(Y))
    result = _flex_binary_moment(arg1, arg2, _get_ewmcov,
                                 pairwise=bool(pairwise))
    if not bias:
        com = _get_center_of_mass(com, span, halflife)
        result *= (1.0 + 2.0 * com) / (2.0 * com)

    return result


@Substitution("Exponentially-weighted moving correlation", _binary_arg_flex,
              _ewm_kw+_pairwise_kw, _type_of_input_retval, _ewm_notes)
@Appender(_doc_template)
def ewmcorr(arg1, arg2=None, com=None, span=None, halflife=None, min_periods=0,
            freq=None, pairwise=None, how=None):
    if arg2 is None:
        arg2 = arg1
        pairwise = True if pairwise is None else pairwise
    elif isinstance(arg2, (int, float)) and com is None:
        com = arg2
        arg2 = arg1
        pairwise = True if pairwise is None else pairwise
    arg1 = _conv_timerule(arg1, freq, how)
    arg2 = _conv_timerule(arg2, freq, how)

    def _get_ewmcorr(X, Y):
        mean = lambda x: ewma(x, com=com, span=span, halflife=halflife, min_periods=min_periods)
        var = lambda x: ewmvar(x, com=com, span=span, halflife=halflife, min_periods=min_periods,
                               bias=True)
        return (mean(X * Y) - mean(X) * mean(Y)) / _zsqrt(var(X) * var(Y))
    result = _flex_binary_moment(arg1, arg2, _get_ewmcorr,
                                 pairwise=bool(pairwise))
    return result


def _zsqrt(x):
    result = np.sqrt(x)
    mask = x < 0

    if isinstance(x, DataFrame):
        if mask.values.any():
            result[mask] = 0
    else:
        if mask.any():
            result[mask] = 0

    return result


def _prep_binary(arg1, arg2):
    if not isinstance(arg2, type(arg1)):
        raise Exception('Input arrays must be of the same type!')

    # mask out values, this also makes a common index...
    X = arg1 + 0 * arg2
    Y = arg2 + 0 * arg1

    return X, Y

#----------------------------------------------------------------------
# Python interface to Cython functions


def _conv_timerule(arg, freq, how):

    types = (DataFrame, Series)
    if freq is not None and isinstance(arg, types):
        # Conform to whatever frequency needed.
        arg = arg.resample(freq, how=how)

    return arg


def _require_min_periods(p):
    def _check_func(minp, window):
        if minp is None:
            return window
        else:
            return max(p, minp)
    return _check_func


def _use_window(minp, window):
    if minp is None:
        return window
    else:
        return minp


def _rolling_func(func, desc, check_minp=_use_window, how=None):
    if how is None:
        how_arg_str = 'None'
    else:
        how_arg_str = "'%s"%how

    @Substitution(desc, _unary_arg, _roll_kw%how_arg_str, _type_of_input_retval,
                  _roll_notes)
    @Appender(_doc_template)
    @wraps(func)
    def f(arg, window, min_periods=None, freq=None, center=False, how=how,
          **kwargs):
        def call_cython(arg, window, minp, args=(), kwargs={}, **kwds):
            minp = check_minp(minp, window)
            return func(arg, window, minp, **kwds)
        return _rolling_moment(arg, window, call_cython, min_periods, freq=freq,
                               center=center, how=how, **kwargs)

    return f

rolling_max = _rolling_func(algos.roll_max2, 'Moving maximum.', how='max')
rolling_min = _rolling_func(algos.roll_min2, 'Moving minimum.', how='min')
rolling_sum = _rolling_func(algos.roll_sum, 'Moving sum.')
rolling_mean = _rolling_func(algos.roll_mean, 'Moving mean.')
rolling_median = _rolling_func(algos.roll_median_cython, 'Moving median.',
                               how='median')

_ts_std = lambda *a, **kw: _zsqrt(algos.roll_var(*a, **kw))
rolling_std = _rolling_func(_ts_std, 'Unbiased moving standard deviation.',
                            check_minp=_require_min_periods(1))
rolling_var = _rolling_func(algos.roll_var, 'Unbiased moving variance.',
                            check_minp=_require_min_periods(1))
rolling_skew = _rolling_func(algos.roll_skew, 'Unbiased moving skewness.',
                             check_minp=_require_min_periods(3))
rolling_kurt = _rolling_func(algos.roll_kurt, 'Unbiased moving kurtosis.',
                             check_minp=_require_min_periods(4))


def rolling_quantile(arg, window, quantile, min_periods=None, freq=None,
                     center=False):
    """Moving quantile.

    Parameters
    ----------
    arg : Series, DataFrame
    window : int
        Size of the moving window. This is the number of observations used for
        calculating the statistic.
    quantile : float
        0 <= quantile <= 1
    min_periods : int, default None
        Minimum number of observations in window required to have a value
        (otherwise result is NA).
    freq : string or DateOffset object, optional (default None)
        Frequency to conform the data to before computing the statistic. Specified
        as a frequency string or DateOffset object.
    center : boolean, default False
        Whether the label should correspond with center of window

    Returns
    -------
    y : type of input argument

    Notes
    -----
    By default, the result is set to the right edge of the window. This can be
    changed to the center of the window by setting ``center=True``.

    The `freq` keyword is used to conform time series data to a specified
    frequency by resampling the data. This is done with the default parameters
    of :meth:`~pandas.Series.resample` (i.e. using the `mean`).
    """

    def call_cython(arg, window, minp, args=(), kwargs={}):
        minp = _use_window(minp, window)
        return algos.roll_quantile(arg, window, minp, quantile)
    return _rolling_moment(arg, window, call_cython, min_periods, freq=freq,
                           center=center)


def rolling_apply(arg, window, func, min_periods=None, freq=None,
                  center=False, args=(), kwargs={}):
    """Generic moving function application.

    Parameters
    ----------
    arg : Series, DataFrame
    window : int
        Size of the moving window. This is the number of observations used for
        calculating the statistic.
    func : function
        Must produce a single value from an ndarray input
    min_periods : int, default None
        Minimum number of observations in window required to have a value
        (otherwise result is NA).
    freq : string or DateOffset object, optional (default None)
        Frequency to conform the data to before computing the statistic. Specified
        as a frequency string or DateOffset object.
    center : boolean, default False
        Whether the label should correspond with center of window
    args : tuple
        Passed on to func
    kwargs : dict
        Passed on to func

    Returns
    -------
    y : type of input argument

    Notes
    -----
    By default, the result is set to the right edge of the window. This can be
    changed to the center of the window by setting ``center=True``.

    The `freq` keyword is used to conform time series data to a specified
    frequency by resampling the data. This is done with the default parameters
    of :meth:`~pandas.Series.resample` (i.e. using the `mean`).
    """
    def call_cython(arg, window, minp, args, kwargs):
        minp = _use_window(minp, window)
        return algos.roll_generic(arg, window, minp, func, args, kwargs)
    return _rolling_moment(arg, window, call_cython, min_periods, freq=freq,
                           center=center, args=args, kwargs=kwargs)


def rolling_window(arg, window=None, win_type=None, min_periods=None,
                   freq=None, center=False, mean=True,
                   axis=0, how=None, **kwargs):
    """
    Applies a moving window of type ``window_type`` and size ``window``
    on the data.

    Parameters
    ----------
    arg : Series, DataFrame
    window : int or ndarray
        Weighting window specification. If the window is an integer, then it is
        treated as the window length and win_type is required
    win_type : str, default None
        Window type (see Notes)
    min_periods : int, default None
        Minimum number of observations in window required to have a value
        (otherwise result is NA).
    freq : string or DateOffset object, optional (default None)
        Frequency to conform the data to before computing the statistic. Specified
        as a frequency string or DateOffset object.
    center : boolean, default False
        Whether the label should correspond with center of window
    mean : boolean, default True
        If True computes weighted mean, else weighted sum
    axis : {0, 1}, default 0
    how : string, default 'mean'
        Method for down- or re-sampling

    Returns
    -------
    y : type of input argument

    Notes
    -----
    The recognized window types are:

    * ``boxcar``
    * ``triang``
    * ``blackman``
    * ``hamming``
    * ``bartlett``
    * ``parzen``
    * ``bohman``
    * ``blackmanharris``
    * ``nuttall``
    * ``barthann``
    * ``kaiser`` (needs beta)
    * ``gaussian`` (needs std)
    * ``general_gaussian`` (needs power, width)
    * ``slepian`` (needs width).

    By default, the result is set to the right edge of the window. This can be
    changed to the center of the window by setting ``center=True``.

    The `freq` keyword is used to conform time series data to a specified
    frequency by resampling the data. This is done with the default parameters
    of :meth:`~pandas.Series.resample` (i.e. using the `mean`).
    """
    if isinstance(window, (list, tuple, np.ndarray)):
        if win_type is not None:
            raise ValueError(('Do not specify window type if using custom '
                              'weights'))
        window = pdcom._asarray_tuplesafe(window).astype(float)
    elif pdcom.is_integer(window):  # window size
        if win_type is None:
            raise ValueError('Must specify window type')
        try:
            import scipy.signal as sig
        except ImportError:
            raise ImportError('Please install scipy to generate window weight')
        win_type = _validate_win_type(win_type, kwargs)  # may pop from kwargs
        window = sig.get_window(win_type, window).astype(float)
    else:
        raise ValueError('Invalid window %s' % str(window))

    minp = _use_window(min_periods, len(window))

    arg = _conv_timerule(arg, freq, how)
    return_hook, values = _process_data_structure(arg)

    f = lambda x: algos.roll_window(x, window, minp, avg=mean)
    result = np.apply_along_axis(f, axis, values)

    rs = return_hook(result)
    if center:
        rs = _center_window(rs, len(window), axis)
    return rs


def _validate_win_type(win_type, kwargs):
    # may pop from kwargs
    arg_map = {'kaiser': ['beta'],
               'gaussian': ['std'],
               'general_gaussian': ['power', 'width'],
               'slepian': ['width']}
    if win_type in arg_map:
        return tuple([win_type] +
                     _pop_args(win_type, arg_map[win_type], kwargs))
    return win_type


def _pop_args(win_type, arg_names, kwargs):
    msg = '%s window requires %%s' % win_type
    all_args = []
    for n in arg_names:
        if n not in kwargs:
            raise ValueError(msg % n)
        all_args.append(kwargs.pop(n))
    return all_args


def _expanding_func(func, desc, check_minp=_use_window):
    @Substitution(desc, _unary_arg, _expanding_kw, _type_of_input_retval, "")
    @Appender(_doc_template)
    @wraps(func)
    def f(arg, min_periods=1, freq=None, center=False, **kwargs):
        window = len(arg)

        def call_cython(arg, window, minp, args=(), kwargs={}, **kwds):
            minp = check_minp(minp, window)
            return func(arg, window, minp, **kwds)
        return _rolling_moment(arg, window, call_cython, min_periods, freq=freq,
                               center=center, **kwargs)

    return f

expanding_max = _expanding_func(algos.roll_max2, 'Expanding maximum.')
expanding_min = _expanding_func(algos.roll_min2, 'Expanding minimum.')
expanding_sum = _expanding_func(algos.roll_sum, 'Expanding sum.')
expanding_mean = _expanding_func(algos.roll_mean, 'Expanding mean.')
expanding_median = _expanding_func(
    algos.roll_median_cython, 'Expanding median.')

expanding_std = _expanding_func(_ts_std,
                                'Unbiased expanding standard deviation.',
                                check_minp=_require_min_periods(2))
expanding_var = _expanding_func(algos.roll_var, 'Unbiased expanding variance.',
                                check_minp=_require_min_periods(2))
expanding_skew = _expanding_func(
    algos.roll_skew, 'Unbiased expanding skewness.',
    check_minp=_require_min_periods(3))
expanding_kurt = _expanding_func(
    algos.roll_kurt, 'Unbiased expanding kurtosis.',
    check_minp=_require_min_periods(4))


def expanding_count(arg, freq=None, center=False):
    """
    Expanding count of number of non-NaN observations.

    Parameters
    ----------
    arg :  DataFrame or numpy ndarray-like
    freq : string or DateOffset object, optional (default None)
        Frequency to conform the data to before computing the statistic. Specified
        as a frequency string or DateOffset object.
    center : boolean, default False
        Whether the label should correspond with center of window.

    Returns
    -------
    expanding_count : type of caller

    Notes
    -----
    The `freq` keyword is used to conform time series data to a specified
    frequency by resampling the data. This is done with the default parameters
    of :meth:`~pandas.Series.resample` (i.e. using the `mean`).
    """
    return rolling_count(arg, len(arg), freq=freq, center=center)


def expanding_quantile(arg, quantile, min_periods=1, freq=None,
                       center=False):
    """Expanding quantile.

    Parameters
    ----------
    arg : Series, DataFrame
    quantile : float
        0 <= quantile <= 1
    min_periods : int, default None
        Minimum number of observations in window required to have a value
        (otherwise result is NA).
    freq : string or DateOffset object, optional (default None)
        Frequency to conform the data to before computing the statistic. Specified
        as a frequency string or DateOffset object.
    center : boolean, default False
        Whether the label should correspond with center of window.

    Returns
    -------
    y : type of input argument

    Notes
    -----
    The `freq` keyword is used to conform time series data to a specified
    frequency by resampling the data. This is done with the default parameters
    of :meth:`~pandas.Series.resample` (i.e. using the `mean`).
    """
    return rolling_quantile(arg, len(arg), quantile, min_periods=min_periods,
                            freq=freq, center=center)


@Substitution("Unbiased expanding covariance.", _binary_arg_flex,
              _expanding_kw+_pairwise_kw, _flex_retval, "")
@Appender(_doc_template)
def expanding_cov(arg1, arg2=None, min_periods=1, freq=None, center=False,
                  pairwise=None):
    if arg2 is None:
        arg2 = arg1
        pairwise = True if pairwise is None else pairwise
    elif isinstance(arg2, (int, float)) and min_periods is None:
        min_periods = arg2
        arg2 = arg1
        pairwise = True if pairwise is None else pairwise
    window = max(len(arg1), len(arg2))
    return rolling_cov(arg1, arg2, window,
                       min_periods=min_periods, freq=freq,
                       center=center, pairwise=pairwise)


@Substitution("Expanding sample correlation.", _binary_arg_flex,
              _expanding_kw+_pairwise_kw, _flex_retval, "")
@Appender(_doc_template)
def expanding_corr(arg1, arg2=None, min_periods=1, freq=None, center=False,
                   pairwise=None):
    if arg2 is None:
        arg2 = arg1
        pairwise = True if pairwise is None else pairwise
    elif isinstance(arg2, (int, float)) and min_periods is None:
        min_periods = arg2
        arg2 = arg1
        pairwise = True if pairwise is None else pairwise
    window = max(len(arg1), len(arg2))
    return rolling_corr(arg1, arg2, window,
                        min_periods=min_periods,
                        freq=freq, center=center, pairwise=pairwise)


@Substitution("Deprecated. Use expanding_corr(..., pairwise=True) instead.\n\n"
              "Pairwise expanding sample correlation", _pairwise_arg,
              _expanding_kw, _pairwise_retval, "")
@Appender(_doc_template)
def expanding_corr_pairwise(df1, df2=None, min_periods=1, freq=None,
                            center=False):
    import warnings
    warnings.warn("expanding_corr_pairwise is deprecated, use expanding_corr(..., pairwise=True)", FutureWarning)
    return expanding_corr(df1, df2, min_periods=min_periods,
                          freq=freq, center=center, pairwise=True)


def expanding_apply(arg, func, min_periods=1, freq=None, center=False,
                    args=(), kwargs={}):
    """Generic expanding function application.

    Parameters
    ----------
    arg : Series, DataFrame
    func : function
        Must produce a single value from an ndarray input
    min_periods : int, default None
        Minimum number of observations in window required to have a value
        (otherwise result is NA).
    freq : string or DateOffset object, optional (default None)
        Frequency to conform the data to before computing the statistic. Specified
        as a frequency string or DateOffset object.
    center : boolean, default False
        Whether the label should correspond with center of window.
    args : tuple
        Passed on to func
    kwargs : dict
        Passed on to func

    Returns
    -------
    y : type of input argument

    Notes
    -----
    The `freq` keyword is used to conform time series data to a specified
    frequency by resampling the data. This is done with the default parameters
    of :meth:`~pandas.Series.resample` (i.e. using the `mean`).
    """
    window = len(arg)
    return rolling_apply(arg, window, func, min_periods=min_periods, freq=freq,
                         center=center, args=args, kwargs=kwargs)
