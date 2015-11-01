"""

provide a generic structure to support window functions,
similar to how we have a Groupby object


"""
from __future__ import division

import numpy as np
from functools import wraps
from collections import defaultdict

import pandas as pd
from pandas.core.base import PandasObject, SelectionMixin, AbstractMethodError
import pandas.core.common as com
import pandas.algos as algos
from pandas import compat
from pandas.util.decorators import Substitution, Appender

class _Window(PandasObject, SelectionMixin):
    _attributes = ['window','min_periods','freq','center','how','win_type','axis']
    exclusions = set()

    def __init__(self, obj, window=None, min_periods=None, freq=None, center=False,
                 how=None, win_type=None, axis=0):
        self.blocks = []
        self.obj = obj
        self.window = window
        self.min_periods = min_periods
        self.freq = freq
        self.center = center
        self.how = how
        self.win_type = win_type
        self.axis = axis
        self._convert_freq()
        self._setup()

    @property
    def _constructor(self):
        return Window

    def _setup(self):
        pass

    def _create_blocks(self):
        """ split data into blocks """
        return self._selected_obj.as_blocks(copy=False).values()

    def _gotitem(self, key, ndim, subset=None):
        """
        sub-classes to define
        return a sliced object

        Parameters
        ----------
        key : string / list of selections
        ndim : 1,2
            requested ndim of result
        subset : object, default None
            subset to act on
        """

        # create a new object to prevent aliasing
        if subset is None:
            subset = self.obj
        new_self = self._shallow_copy(subset)
        if ndim==2 and key in subset:
            new_self._selection = key
        new_self._reset_cache()
        return new_self

    def __getattr__(self, attr):
        if attr in self._internal_names_set:
            return object.__getattribute__(self, attr)
        if attr in self.obj:
            return self[attr]

        raise AttributeError("%r object has no attribute %r" %
                             (type(self).__name__, attr))

    def _dir_additions(self):
        return self.obj._dir_additions()

    def _get_window(self, other=None):
        return self.window

    def __unicode__(self):
        """ provide a nice str repr of our rolling object """

        attrs = [ "{k}->{v}".format(k=k,v=getattr(self,k)) \
                  for k in self._attributes if getattr(self,k,None) is not None ]
        return "{klass} [{attrs}]".format(klass=self.__class__.__name__,
                                          attrs=','.join(attrs))

    def _shallow_copy(self, obj=None, **kwargs):
        """ return a new object with the replacement attributes """
        if obj is None:
            obj = self._selected_obj.copy()
        if isinstance(obj, self.__class__):
            obj = obj.obj
        for attr in self._attributes:
            if attr not in kwargs:
                kwargs[attr] = getattr(self,attr)
        return self._constructor(obj, **kwargs)

    def _prep_values(self, values=None, kill_inf=True):

        if values is None:
            values = getattr(self._selected_obj,'values',self._selected_obj)

        # coerce dtypes as appropriate
        if com.is_float_dtype(values.dtype):
            pass
        elif com.is_integer_dtype(values.dtype):
            values = values.astype(float)
        elif com.is_timedelta64_dtype(values.dtype):
            values = values.view('i8').astype(float)
        else:
            try:
                values = values.astype(float)
            except (ValueError, TypeError):
                raise TypeError("cannot handle this type -> {0}".format(values.dtype))

        if kill_inf:
            values = values.copy()
            values[np.isinf(values)] = np.NaN

        return values

    def _wrap_result(self, result, block=None):
        """ wrap a single result """

        obj = self._selected_obj
        if isinstance(result, np.ndarray):

            # coerce if necessary
            if block is not None:
                if com.is_timedelta64_dtype(block.values.dtype):
                    result = pd.to_timedelta(result.ravel(),unit='ns').values.reshape(result.shape)

            if result.ndim == 1:
                from pandas import Series
                return Series(result, obj.index, name=obj.name)

            return type(obj)(result,
                             index=obj.index,
                             columns=block.columns)
        return result

    def _wrap_results(self, results, blocks):
        """ wrap lists of results, blocks """

        obj = self._selected_obj
        final = []
        for result, block in zip(results, blocks):

            result = self._wrap_result(result, block)
            if result.ndim == 1:
                return result
            final.append(result)

        if not len(final):
            return obj.astype('float64')
        return pd.concat(final,axis=1).reindex(columns=obj.columns)

    def _center_window(self, result, window):
        """ center the result in the window """
        if self.axis > result.ndim-1:
            raise ValueError("Requested axis is larger then no. of argument "
                             "dimensions")

        from pandas import Series, DataFrame
        offset = _offset(window, True)
        if offset > 0:
            if isinstance(result, (Series, DataFrame)):
                result = result.slice_shift(-offset, axis=self.axis)
            else:
                lead_indexer = [slice(None)] * result.ndim
                lead_indexer[self.axis] = slice(offset, None)
                result = np.copy(result[tuple(lead_indexer)])
        return result

    def _convert_freq(self):
        """ conform to our freq """

        from pandas import Series, DataFrame
        if self.freq is not None and isinstance(self.obj, (Series, DataFrame)):
            self.obj = self.obj.resample(self.freq, how=self.how)

    @Appender(SelectionMixin._agg_doc)
    def aggregate(self, arg, *args, **kwargs):
        result, how = self._aggregate(arg, *args, **kwargs)
        if result is None:
            import pdb; pdb.set_trace()
        return result

class Window(_Window):

    def _prep_window(self, **kwargs):
        """ provide validation for our window type, return the window """
        window = self._get_window()

        if isinstance(window, (list, tuple, np.ndarray)):
            return com._asarray_tuplesafe(window).astype(float)
        elif com.is_integer(window):
            try:
                import scipy.signal as sig
            except ImportError:
                raise ImportError('Please install scipy to generate window weight')
            win_type = _validate_win_type(self.win_type, kwargs)  # may pop from kwargs
            return sig.get_window(win_type, window).astype(float)

        raise ValueError('Invalid window %s' % str(window))

    def _apply_window(self, mean=True, **kwargs):
        """
        Applies a moving window of type ``window_type`` on the data.

        Parameters
        ----------
        mean : boolean, default True
            If True computes weighted mean, else weighted sum

        Returns
        -------
        y : type of input argument

        """
        window = self._prep_window(**kwargs)
        center = self.center

        results, blocks = [], self._create_blocks()
        for b in blocks:
            try:
                values = self._prep_values(b.values)
            except TypeError:
                results.append(b.values.copy())
                continue

            if values.size == 0:
                results.append(values.copy())
                continue

            offset = _offset(window, center)
            additional_nans = np.array([np.NaN] * offset)
            def f(arg, *args, **kwargs):
                minp = _use_window(self.min_periods, len(window))
                return algos.roll_window(np.concatenate((arg, additional_nans)) if center else arg,
                                         window, minp, avg=mean)

            result = np.apply_along_axis(f, self.axis, values)

            if center:
                result = self._center_window(result, window)
            results.append(result)

        return self._wrap_results(results, blocks)

    def sum(self, **kwargs):
        return self._apply_window(mean=False, **kwargs)

    def mean(self, **kwargs):
        return self._apply_window(mean=True, **kwargs)

class _Rolling(_Window):

    @property
    def _constructor(self):
        return Rolling

    def _apply(self, func, window=None, center=None, check_minp=None, how=None, **kwargs):
        """
        Rolling statistical measure using supplied function. Designed to be
        used with passed-in Cython array-based functions.

        Parameters
        ----------
        func : string/callable to apply
        window : int/array, default to _get_window()
        center : boolean, default to self.center
        check_minp : function, default to _use_window
        how : string, default to None

        Returns
        -------
        y : type of input
        """

        if center is None:
            center = self.center
        if window is None:
            window = self._get_window()

        if check_minp is None:
            check_minp = _use_window

        results, blocks = [], self._create_blocks()
        for b in blocks:
            try:
                values = self._prep_values(b.values)
            except TypeError:
                results.append(b.values.copy())
                continue

            if values.size == 0:
                results.append(values.copy())
                continue

            # if we have a string function name, wrap it
            if isinstance(func, compat.string_types):
                if not hasattr(algos, func):
                    raise ValueError("we do not support this function algos.{0}".format(func))

                cfunc = getattr(algos, func)
                def func(arg, window, min_periods=None):
                    minp = check_minp(min_periods, window)
                    return cfunc(arg, window, minp, **kwargs)

            # calculation function
            if center:
                offset = _offset(window, center)
                additional_nans = np.array([np.NaN] * offset)
                def calc(x):
                    return func(np.concatenate((x, additional_nans)),
                                window, min_periods=self.min_periods)
            else:
                def calc(x):
                    return func(x,window, min_periods=self.min_periods)

            if values.ndim > 1:
                result = np.apply_along_axis(calc, self.axis, values)
            else:
                result = calc(values)

            if center:
                result = self._center_window(result, window)

            results.append(result)

        return self._wrap_results(results, blocks)

class Rolling(_Rolling):

    def count(self):
        """
        Rolling count of number of non-NaN observations inside provided window.

        Returns
        -------
        same type as input
        """

        obj = self._selected_obj
        window = self._get_window()
        window = min(window, len(obj)) if not self.center else window
        try:
            converted = np.isfinite(obj).astype(float)
        except TypeError:
            converted = np.isfinite(obj.astype(float)).astype(float)
        result = self._constructor(converted,
                                   window=window,
                                   min_periods=0,
                                   center=self.center).sum()

        result[result.isnull()] = 0
        return result

    def apply(self, func, args=(), kwargs={}):
        """
        Moving function apply

        Parameters
        ----------
        func : function
            Must produce a single value from an ndarray input
        *args and **kwargs are passed to the function
        """
        window = self._get_window()
        offset = _offset(window, self.center)
        def f(arg, window, min_periods):
            minp = _use_window(min_periods, window)
            return algos.roll_generic(arg, window, minp, offset, func, args, kwargs)

        return self._apply(f, center=False)

    def sum(self):
        """
        Moving sum
        """
        return self._apply('roll_sum')

    def max(self, how='max'):
        """
        Moving max

        Parameters
        ----------
        how : string, default max
          Method for down- or re-sampling
        """
        return self._apply('roll_max', how=how)

    def min(self, how='min'):
        """
        Moving min

        Parameters
        ----------
        how : string, default min
          Method for down- or re-sampling
        """
        return self._apply('roll_min', how=how)

    def mean(self):
        """
        Moving mean
        """
        return self._apply('roll_mean')

    def median(self, how='median'):
        """
        Moving median

        Parameters
        ----------
        how : string, default median
          Method for down- or re-sampling
        """

        return self._apply('roll_median_c', how=how)

    def std(self, ddof=1):
        """
        Moving standard deviation

        Parameters
        ----------
        ddof : int, default 1
           Delta Degrees of Freedom.  The divisor used in calculations
           is ``N - ddof``, where ``N`` represents the number of elements.
        """
        window = self._get_window()
        def f(arg, *args, **kwargs):
            minp = _require_min_periods(1)(self.min_periods, window)
            return _zsqrt(algos.roll_var(arg, window, minp, ddof))

        return self._apply(f, check_minp=_require_min_periods(1))

    def var(self, ddof=1):
        """
        Moving variance

        Parameters
        ----------
        ddof : int, default 1
           Delta Degrees of Freedom.  The divisor used in calculations
           is ``N - ddof``, where ``N`` represents the number of elements.
        """
        return self._apply('roll_var',
                           check_minp=_require_min_periods(1),
                           ddof=ddof)

    def skew(self):
        """
        Unbiased moving skewness
        """
        return self._apply('roll_skew',
                           check_minp=_require_min_periods(3))

    def kurt(self):
        """
        Unbiased moving kurtosis
        """
        return self._apply('roll_kurt',
                           check_minp=_require_min_periods(4))

    def quantile(self, quantile):
        """
        Rolling quantile

        Parameters
        ----------
        quantile : float
            0 <= quantile <= 1
        """
        window = self._get_window()
        def f(arg, *args, **kwargs):
            minp = _use_window(self.min_periods, window)
            return algos.roll_quantile(arg, window, minp, quantile)

        return self._apply(f)

    def cov(self, other=None, pairwise=False, ddof=1):
        """
        Moving sample covariance

        Parameters
        ----------
        other : Series, DataFrame, or ndarray, optional
            if not supplied then will default to self and produce pairwise output
        pairwise : bool, default False
            If False then only matching columns between self and other will be used and
            the output will be a DataFrame.
            If True then all pairwise combinations will be calculated and the output
            will be a Panel in the case of DataFrame inputs. In the case of missing
            elements, only complete pairwise observations will be used.
        ddof : int, default 1
            Delta Degrees of Freedom.  The divisor used in calculations
           is ``N - ddof``, where ``N`` represents the number of elements.
        """
        if other is None:
            other = self._selected_obj
            pairwise = True
        other = self._shallow_copy(other)
        window = self._get_window(other)

        def _get_cov(X, Y):
            mean = lambda x: x.rolling(window, self.min_periods, center=self.center).mean()
            count = (X+Y).rolling(window=window, center=self.center).count()
            bias_adj = count / (count - ddof)
            return (mean(X * Y) - mean(X) * mean(Y)) * bias_adj
        return _flex_binary_moment(self._selected_obj, other._selected_obj, _get_cov, pairwise=bool(pairwise))

    def corr(self, other=None, pairwise=False):
        """
        Moving sample correlation

        Parameters
        ----------
        other : Series, DataFrame, or ndarray, optional
            if not supplied then will default to self and produce pairwise output
        pairwise : bool, default False
            If False then only matching columns between self and other will be used and
            the output will be a DataFrame.
            If True then all pairwise combinations will be calculated and the output
            will be a Panel in the case of DataFrame inputs. In the case of missing
            elements, only complete pairwise observations will be used.
        """

        if other is None:
            other = self._selected_obj
            pairwise = True
        other = self._shallow_copy(other)
        window = self._get_window(other)

        def _get_corr(a, b):
            a = a.rolling(window=window,
                          min_periods=self.min_periods,
                          freq=self.freq,
                          center=self.center)
            b = b.rolling(window=window,
                          min_periods=self.min_periods,
                          freq=self.freq,
                          center=self.center)

            return a.cov(b) / (a.std() * b.std())
        return _flex_binary_moment(self._selected_obj, other._selected_obj, _get_corr, pairwise=bool(pairwise))

class Expanding(Rolling):
    _attributes = ['min_periods','freq','center','how','axis']

    @property
    def _constructor(self):
        return Expanding

    def _get_window(self, other=None):
        obj = self._selected_obj
        if other is None:
            return max(len(obj), self.min_periods) if self.min_periods else len(obj)
        return max((len(obj) + len(obj)), self.min_periods) if self.min_periods else (len(obj) + len(obj))

class EWM(_Rolling):
    _attributes = ['com','min_periods','freq','adjust','how','ignore_na','axis']

    def __init__(self, obj, com=None, span=None, halflife=None, min_periods=0, freq=None,
                 adjust=True, how=None, ignore_na=False, axis=0):
        self.obj = obj
        self.com = _get_center_of_mass(com, span, halflife)
        self.min_periods = min_periods
        self.freq = freq
        self.adjust = adjust
        self.how = how
        self.ignore_na = ignore_na
        self.axis = axis
        self._convert_freq()

    @property
    def _constructor(self):
        return EWM

    def _apply(self, func, **kwargs):
        """
        Rolling statistical measure using supplied function. Designed to be
        used with passed-in Cython array-based functions.

        Parameters
        ----------
        func : string/callable to apply

        Returns
        -------
        y : type of input argument

        """
        results, blocks = [], self._create_blocks()
        for b in blocks:
            try:
                values = self._prep_values(b.values)
            except TypeError:
                results.append(b.values.copy())
                continue

            if values.size == 0:
                results.append(values.copy())
                continue

            # if we have a string function name, wrap it
            if isinstance(func, compat.string_types):
                if not hasattr(algos, func):
                    raise ValueError("we do not support this function algos.{0}".format(func))

                cfunc = getattr(algos, func)
                def func(arg):
                    return cfunc(arg, self.com, int(self.adjust), int(self.ignore_na), int(self.min_periods))

            results.append(np.apply_along_axis(func, self.axis, values))

        return self._wrap_results(results, blocks)

    def mean(self):
        """
        exponential weighted moving average
        """
        return self._apply('ewma')

    def std(self, bias=False):
        """
        exponential weighted moving stddev

        Parameters
        ----------
        bias : boolean, default False
           Use a standard estimation bias correction
        """
        return _zsqrt(self.var(bias=bias))
    vol=std

    def var(self, bias=False):
        """
        exponential weighted moving average

        Parameters
        ----------
        bias : boolean, default False
           Use a standard estimation bias correction
        """
        def f(arg):
            return algos.ewmcov(arg,
                                arg,
                                self.com,
                                int(self.adjust),
                                int(self.ignore_na),
                                int(self.min_periods),
                                int(bias))

        return self._apply(f)

    def cov(self, other=None, pairwise=False, bias=False):
        """
        exponential weighted sample covariance

        Parameters
        ----------
        other : Series, DataFrame, or ndarray, optional
            if not supplied then will default to self and produce pairwise output
        pairwise : bool, default False
            If False then only matching columns between self and other will be used and
            the output will be a DataFrame.
            If True then all pairwise combinations will be calculated and the output
            will be a Panel in the case of DataFrame inputs. In the case of missing
            elements, only complete pairwise observations will be used.
        bias : boolean, default False
           Use a standard estimation bias correction
        """
        if other is None:
            other = self._selected_obj
            pairwise = True
        other = self._shallow_copy(other)

        def _get_cov(X, Y):
            X = self._shallow_copy(X)
            Y = self._shallow_copy(Y)
            cov = algos.ewmcov(X._prep_values(),
                               Y._prep_values(),
                               self.com,
                               int(self.adjust),
                               int(self.ignore_na),
                               int(self.min_periods),
                               int(bias))
            return X._wrap_result(cov)

        return _flex_binary_moment(self._selected_obj, other._selected_obj, _get_cov, pairwise=bool(pairwise))

    def corr(self, other=None, pairwise=False):
        """
        exponential weighted sample correlation

        Parameters
        ----------
        other : Series, DataFrame, or ndarray, optional
            if not supplied then will default to self and produce pairwise output
        pairwise : bool, default False
            If False then only matching columns between self and other will be used and
            the output will be a DataFrame.
            If True then all pairwise combinations will be calculated and the output
            will be a Panel in the case of DataFrame inputs. In the case of missing
            elements, only complete pairwise observations will be used.
        """
        if other is None:
            other = self._selected_obj
            pairwise = True
        other = self._shallow_copy(other)

        def _get_corr(X, Y):
            X = self._shallow_copy(X)
            Y = self._shallow_copy(Y)
            def _cov(x, y):
                return algos.ewmcov(x, y, self.com, int(self.adjust), int(self.ignore_na), int(self.min_periods), 1)

            x_values = X._prep_values()
            y_values = Y._prep_values()
            cov = _cov(x_values, y_values)
            x_var = _cov(x_values, x_values)
            y_var = _cov(y_values, y_values)
            corr = cov / _zsqrt(x_var * y_var)
            return X._wrap_result(corr)

        return _flex_binary_moment(self._selected_obj, other._selected_obj, _get_corr, pairwise=bool(pairwise))

########################
##### Helper Funcs #####
########################

def _flex_binary_moment(arg1, arg2, f, pairwise=False):
    from pandas import Series, DataFrame, Panel
    if not (isinstance(arg1,(np.ndarray, Series, DataFrame)) and
            isinstance(arg2,(np.ndarray, Series, DataFrame))):
        raise TypeError("arguments to moment function must be of type "
                         "np.ndarray/Series/DataFrame")

    if isinstance(arg1, (np.ndarray, Series)) and \
            isinstance(arg2, (np.ndarray,Series)):
        X, Y = _prep_binary(arg1, arg2)
        return f(X, Y)

    elif isinstance(arg1, DataFrame):
        def dataframe_from_int_dict(data, frame_template):
            result = DataFrame(data, index=frame_template.index)
            if len(result.columns) > 0:
                result.columns = frame_template.columns[result.columns]
            return result

        results = {}
        if isinstance(arg2, DataFrame):
            if pairwise is False:
                if arg1 is arg2:
                    # special case in order to handle duplicate column names
                    for i, col in enumerate(arg1.columns):
                        results[i] = f(arg1.iloc[:, i], arg2.iloc[:, i])
                    return dataframe_from_int_dict(results, arg1)
                else:
                    if not arg1.columns.is_unique:
                        raise ValueError("'arg1' columns are not unique")
                    if not arg2.columns.is_unique:
                        raise ValueError("'arg2' columns are not unique")
                    X, Y = arg1.align(arg2, join='outer')
                    X = X + 0 * Y
                    Y = Y + 0 * X
                    res_columns = arg1.columns.union(arg2.columns)
                    for col in res_columns:
                        if col in X and col in Y:
                            results[col] = f(X[col], Y[col])
                    return DataFrame(results, index=X.index, columns=res_columns)
            elif pairwise is True:
                results = defaultdict(dict)
                for i, k1 in enumerate(arg1.columns):
                    for j, k2 in enumerate(arg2.columns):
                        if j<i and arg2 is arg1:
                            # Symmetric case
                            results[i][j] = results[j][i]
                        else:
                            results[i][j] = f(*_prep_binary(arg1.iloc[:, i], arg2.iloc[:, j]))
                p = Panel.from_dict(results).swapaxes('items', 'major')
                if len(p.major_axis) > 0:
                    p.major_axis = arg1.columns[p.major_axis]
                if len(p.minor_axis) > 0:
                    p.minor_axis = arg2.columns[p.minor_axis]
                return p
            else:
                raise ValueError("'pairwise' is not True/False")
        else:
            results = {}
            for i, col in enumerate(arg1.columns):
                results[i] = f(*_prep_binary(arg1.iloc[:, i], arg2))
            return dataframe_from_int_dict(results, arg1)

    else:
        return _flex_binary_moment(arg2, arg1, f)

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

def _offset(window, center):
    if not com.is_integer(window):
        window = len(window)
    offset = (window - 1) / 2. if center else 0
    try:
        return int(offset)
    except:
        return offset.astype(int)

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

def _zsqrt(x):
    result = np.sqrt(x)
    mask = x < 0

    from pandas import DataFrame
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

#############################
##### top-level exports #####
#############################

def rolling(obj, win_type=None, **kwds):
    """
    Provides rolling transformations.

    .. versionadded:: 0.18.0

    Parameters
    ----------
    window : int
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
    how : string, default None
        Method for down- or re-sampling
    win_type : string, default None
        prove a window type, see the notes below
    axis : int, default 0

    Returns
    -------
    a Window sub-classed for the particular operation

    Notes
    -----
    By default, the result is set to the right edge of the window. This can be
    changed to the center of the window by setting ``center=True``.

    The `freq` keyword is used to conform time series data to a specified
    frequency by resampling the data. This is done with the default parameters
    of :meth:`~pandas.Series.resample` (i.e. using the `mean`).

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
    """
    from pandas import Series, DataFrame
    if not isinstance(obj, (Series, DataFrame)):
        raise TypeError('invalid type: %s' % type(obj))

    if win_type is not None:
        return Window(obj, win_type=win_type, **kwds)

    return Rolling(obj, **kwds)

def expanding(obj, **kwds):
    """
    Provides expanding transformations.

    .. versionadded:: 0.18.0

    Parameters
    ----------
    min_periods : int, default None
        Minimum number of observations in window required to have a value
        (otherwise result is NA).
    freq : string or DateOffset object, optional (default None)
        Frequency to conform the data to before computing the statistic. Specified
        as a frequency string or DateOffset object.
    center : boolean, default False
        Set the labels at the center of the window.
    how : string, default None
        Method for down- or re-sampling
    axis : int, default 0

    Returns
    -------
    a Window sub-classed for the particular operation

    Notes
    -----
    By default, the result is set to the right edge of the window. This can be
    changed to the center of the window by setting ``center=True``.

    The `freq` keyword is used to conform time series data to a specified
    frequency by resampling the data. This is done with the default parameters
    of :meth:`~pandas.Series.resample` (i.e. using the `mean`).
    """

    from pandas import Series, DataFrame
    if not isinstance(obj, (Series, DataFrame)):
        raise TypeError('invalid type: %s' % type(obj))

    return Expanding(obj, **kwds)

def ewm(obj, **kwds):
    """
    .. versionadded:: 0.18.0

    Provides exponential weighted functions

    Parameters
    ----------
    com : float. optional
        Center of mass: :math:`\alpha = 1 / (1 + com)`,
    span : float, optional
        Specify decay in terms of span, :math:`\alpha = 2 / (span + 1)`
    halflife : float, optional
        Specify decay in terms of halflife, :math:`\alpha = 1 - exp(log(0.5) / halflife)`
    min_periods : int, default 0
        Minimum number of observations in window required to have a value
        (otherwise result is NA).
    freq : None or string alias / date offset object, default=None
        Frequency to conform to before computing statistic
    adjust : boolean, default True
        Divide by decaying adjustment factor in beginning periods to account for
        imbalance in relative weightings (viewing EWMA as a moving average)
    how : string, default 'mean'
        Method for down- or re-sampling
    ignore_na : boolean, default False
        Ignore missing values when calculating weights;
        specify True to reproduce pre-0.15.0 behavior

    Returns
    -------
    a Window sub-classed for the particular operation

    Notes
    -----
    Either center of mass, span or halflife must be specified

    EWMA is sometimes specified using a "span" parameter `s`, we have that the
    decay parameter :math:`\alpha` is related to the span as
    :math:`\alpha = 2 / (s + 1) = 1 / (1 + c)`

    where `c` is the center of mass. Given a span, the associated center of mass is
    :math:`c = (s - 1) / 2`

    So a "20-day EWMA" would have center 9.5.

    The `freq` keyword is used to conform time series data to a specified
    frequency by resampling the data. This is done with the default parameters
    of :meth:`~pandas.Series.resample` (i.e. using the `mean`).

    When adjust is True (default), weighted averages are calculated using weights
        (1-alpha)**(n-1), (1-alpha)**(n-2), ..., 1-alpha, 1.

    When adjust is False, weighted averages are calculated recursively as:
       weighted_average[0] = arg[0];
       weighted_average[i] = (1-alpha)*weighted_average[i-1] + alpha*arg[i].

    When ignore_na is False (default), weights are based on absolute positions.
    For example, the weights of x and y used in calculating the final weighted
    average of [x, None, y] are (1-alpha)**2 and 1 (if adjust is True), and
    (1-alpha)**2 and alpha (if adjust is False).

    When ignore_na is True (reproducing pre-0.15.0 behavior), weights are based on
    relative positions. For example, the weights of x and y used in calculating
    the final weighted average of [x, None, y] are 1-alpha and 1 (if adjust is
    True), and 1-alpha and alpha (if adjust is False).

    More details can be found at
    http://pandas.pydata.org/pandas-docs/stable/computation.html#exponentially-weighted-moment-functions
    """
    from pandas import Series, DataFrame
    if not isinstance(obj, (Series, DataFrame)):
        raise TypeError('invalid type: %s' % type(obj))

    return EWM(obj, **kwds)
