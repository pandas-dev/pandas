"""

provide a generic structure to support window functions,
similar to how we have a Groupby object


"""
from __future__ import division

import warnings
import numpy as np
from collections import defaultdict

import pandas as pd
from pandas.lib import isscalar
from pandas.core.base import (PandasObject, SelectionMixin,
                              GroupByMixin)
import pandas.core.common as com
import pandas.algos as algos
from pandas import compat
from pandas.util.decorators import Substitution, Appender
from textwrap import dedent

_shared_docs = dict()
_doc_template = """

Returns
-------
same type as input

See also
--------
pandas.Series.%(name)s
pandas.DataFrame.%(name)s
"""


class _Window(PandasObject, SelectionMixin):
    _attributes = ['window', 'min_periods', 'freq', 'center', 'win_type',
                   'axis']
    exclusions = set()

    def __init__(self, obj, window=None, min_periods=None, freq=None,
                 center=False, win_type=None, axis=0, **kwargs):

        if freq is not None:
            warnings.warn("The freq kw is deprecated and will be removed in a "
                          "future version. You can resample prior to passing "
                          "to a window function", FutureWarning, stacklevel=3)

        self.blocks = []
        self.obj = obj
        self.window = window
        self.min_periods = min_periods
        self.freq = freq
        self.center = center
        self.win_type = win_type
        self.axis = obj._get_axis_number(axis) if axis is not None else None
        self.validate()

    @property
    def _constructor(self):
        return Window

    def validate(self):
        if self.center is not None and not com.is_bool(self.center):
            raise ValueError("center must be a boolean")
        if self.min_periods is not None and not \
           com.is_integer(self.min_periods):
            raise ValueError("min_periods must be an integer")

    def _convert_freq(self, how=None):
        """ resample according to the how, return a new object """

        obj = self._selected_obj
        if (self.freq is not None and
                isinstance(obj, (com.ABCSeries, com.ABCDataFrame))):
            if how is not None:
                warnings.warn("The how kw argument is deprecated and removed "
                              "in a future version. You can resample prior "
                              "to passing to a window function", FutureWarning,
                              stacklevel=6)

            obj = obj.resample(self.freq).aggregate(how or 'asfreq')
        return obj

    def _create_blocks(self, how):
        """ split data into blocks & return conformed data """

        obj = self._convert_freq(how)
        return obj.as_blocks(copy=False).values(), obj

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
        self = self._shallow_copy(subset)
        self._reset_cache()
        if subset.ndim == 2:
            if isscalar(key) and key in subset or com.is_list_like(key):
                self._selection = key
        return self

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

    @property
    def _window_type(self):
        return self.__class__.__name__

    def __unicode__(self):
        """ provide a nice str repr of our rolling object """

        attrs = ["{k}={v}".format(k=k, v=getattr(self, k))
                 for k in self._attributes
                 if getattr(self, k, None) is not None]
        return "{klass} [{attrs}]".format(klass=self._window_type,
                                          attrs=','.join(attrs))

    def _prep_values(self, values=None, kill_inf=True, how=None):

        if values is None:
            values = getattr(self._selected_obj, 'values', self._selected_obj)

        # GH #12373 : rolling functions error on float32 data
        # make sure the data is coerced to float64
        if com.is_float_dtype(values.dtype):
            values = com._ensure_float64(values)
        elif com.is_integer_dtype(values.dtype):
            values = com._ensure_float64(values)
        elif com.needs_i8_conversion(values.dtype):
            raise NotImplementedError("ops for {action} for this "
                                      "dtype {dtype} are not "
                                      "implemented".format(
                                          action=self._window_type,
                                          dtype=values.dtype))
        else:
            try:
                values = com._ensure_float64(values)
            except (ValueError, TypeError):
                raise TypeError("cannot handle this type -> {0}"
                                "".format(values.dtype))

        if kill_inf:
            values = values.copy()
            values[np.isinf(values)] = np.NaN

        return values

    def _wrap_result(self, result, block=None, obj=None):
        """ wrap a single result """

        if obj is None:
            obj = self._selected_obj

        index = obj.index
        if isinstance(result, np.ndarray):

            # coerce if necessary
            if block is not None:
                if com.is_timedelta64_dtype(block.values.dtype):
                    result = pd.to_timedelta(
                        result.ravel(), unit='ns').values.reshape(result.shape)

            if result.ndim == 1:
                from pandas import Series
                return Series(result, index, name=obj.name)

            return type(obj)(result, index=index, columns=block.columns)
        return result

    def _wrap_results(self, results, blocks, obj):
        """
        wrap the results

        Paramters
        ---------
        results : list of ndarrays
        blocks : list of blocks
        obj : conformed data (may be resampled)
        """

        final = []
        for result, block in zip(results, blocks):

            result = self._wrap_result(result, block=block, obj=obj)
            if result.ndim == 1:
                return result
            final.append(result)

        if not len(final):
            return obj.astype('float64')
        return pd.concat(final, axis=1).reindex(columns=obj.columns)

    def _center_window(self, result, window):
        """ center the result in the window """
        if self.axis > result.ndim - 1:
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

    def aggregate(self, arg, *args, **kwargs):
        result, how = self._aggregate(arg, *args, **kwargs)
        if result is None:
            return self.apply(arg, args=args, kwargs=kwargs)
        return result

    agg = aggregate

    _shared_docs['sum'] = dedent("""
    %(name)s sum

    Parameters
    ----------
    how : string, default None (DEPRECATED)
        Method for down- or re-sampling""")

    _shared_docs['mean'] = dedent("""
    %(name)s mean

    Parameters
    ----------
    how : string, default None (DEPRECATED)
        Method for down- or re-sampling""")


class Window(_Window):
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
    freq : string or DateOffset object, optional (default None) (DEPRECATED)
        Frequency to conform the data to before computing the statistic.
        Specified as a frequency string or DateOffset object.
    center : boolean, default False
        Set the labels at the center of the window.
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

    def validate(self):
        super(Window, self).validate()

        window = self.window
        if isinstance(window, (list, tuple, np.ndarray)):
            pass
        elif com.is_integer(window):
            try:
                import scipy.signal as sig
            except ImportError:
                raise ImportError('Please install scipy to generate window '
                                  'weight')

            if not isinstance(self.win_type, compat.string_types):
                raise ValueError('Invalid win_type {0}'.format(self.win_type))
            if getattr(sig, self.win_type, None) is None:
                raise ValueError('Invalid win_type {0}'.format(self.win_type))
        else:
            raise ValueError('Invalid window {0}'.format(window))

    def _prep_window(self, **kwargs):
        """
        provide validation for our window type, return the window
        we have already been validated
        """

        window = self._get_window()
        if isinstance(window, (list, tuple, np.ndarray)):
            return com._asarray_tuplesafe(window).astype(float)
        elif com.is_integer(window):
            import scipy.signal as sig

            # the below may pop from kwargs
            def _validate_win_type(win_type, kwargs):
                arg_map = {'kaiser': ['beta'],
                           'gaussian': ['std'],
                           'general_gaussian': ['power', 'width'],
                           'slepian': ['width']}
                if win_type in arg_map:
                    return tuple([win_type] + _pop_args(win_type,
                                                        arg_map[win_type],
                                                        kwargs))
                return win_type

            def _pop_args(win_type, arg_names, kwargs):
                msg = '%s window requires %%s' % win_type
                all_args = []
                for n in arg_names:
                    if n not in kwargs:
                        raise ValueError(msg % n)
                    all_args.append(kwargs.pop(n))
                return all_args

            win_type = _validate_win_type(self.win_type, kwargs)
            return sig.get_window(win_type, window).astype(float)

    def _apply_window(self, mean=True, how=None, **kwargs):
        """
        Applies a moving window of type ``window_type`` on the data.

        Parameters
        ----------
        mean : boolean, default True
            If True computes weighted mean, else weighted sum
        how : string, default to None (DEPRECATED)
            how to resample

        Returns
        -------
        y : type of input argument

        """
        window = self._prep_window(**kwargs)
        center = self.center

        blocks, obj = self._create_blocks(how=how)
        results = []
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
                return algos.roll_window(np.concatenate((arg, additional_nans))
                                         if center else arg, window, minp,
                                         avg=mean)

            result = np.apply_along_axis(f, self.axis, values)

            if center:
                result = self._center_window(result, window)
            results.append(result)

        return self._wrap_results(results, blocks, obj)

    @Substitution(name='rolling')
    @Appender(SelectionMixin._see_also_template)
    @Appender(SelectionMixin._agg_doc)
    def aggregate(self, arg, *args, **kwargs):
        result, how = self._aggregate(arg, *args, **kwargs)
        if result is None:

            # these must apply directly
            result = arg(self)

        return result

    agg = aggregate

    @Substitution(name='window')
    @Appender(_doc_template)
    @Appender(_shared_docs['sum'])
    def sum(self, **kwargs):
        return self._apply_window(mean=False, **kwargs)

    @Substitution(name='window')
    @Appender(_doc_template)
    @Appender(_shared_docs['mean'])
    def mean(self, **kwargs):
        return self._apply_window(mean=True, **kwargs)


class _GroupByMixin(GroupByMixin):
    """ provide the groupby facilities """

    def __init__(self, obj, *args, **kwargs):
        parent = kwargs.pop('parent', None)  # noqa
        groupby = kwargs.pop('groupby', None)
        if groupby is None:
            groupby, obj = obj, obj.obj
        self._groupby = groupby
        self._groupby.mutated = True
        self._groupby.grouper.mutated = True
        super(GroupByMixin, self).__init__(obj, *args, **kwargs)

    count = GroupByMixin._dispatch('count')
    corr = GroupByMixin._dispatch('corr', other=None, pairwise=None)
    cov = GroupByMixin._dispatch('cov', other=None, pairwise=None)

    def _apply(self, func, name, window=None, center=None,
               check_minp=None, how=None, **kwargs):
        """
        dispatch to apply; we are stripping all of the _apply kwargs and
        performing the original function call on the grouped object
        """

        def f(x, name=name, *args):
            x = self._shallow_copy(x)

            if isinstance(name, compat.string_types):
                return getattr(x, name)(*args, **kwargs)

            return x.apply(name, *args, **kwargs)

        return self._groupby.apply(f)


class _Rolling(_Window):
    @property
    def _constructor(self):
        return Rolling

    def _apply(self, func, name=None, window=None, center=None,
               check_minp=None, how=None, **kwargs):
        """
        Rolling statistical measure using supplied function. Designed to be
        used with passed-in Cython array-based functions.

        Parameters
        ----------
        func : string/callable to apply
        name : string, optional
           name of this function
        window : int/array, default to _get_window()
        center : boolean, default to self.center
        check_minp : function, default to _use_window
        how : string, default to None (DEPRECATED)
            how to resample

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

        blocks, obj = self._create_blocks(how=how)
        results = []
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
                    raise ValueError("we do not support this function "
                                     "algos.{0}".format(func))

                cfunc = getattr(algos, func)

                def func(arg, window, min_periods=None):
                    minp = check_minp(min_periods, window)
                    # GH #12373: rolling functions error on float32 data
                    return cfunc(com._ensure_float64(arg),
                                 window, minp, **kwargs)

            # calculation function
            if center:
                offset = _offset(window, center)
                additional_nans = np.array([np.NaN] * offset)

                def calc(x):
                    return func(np.concatenate((x, additional_nans)),
                                window, min_periods=self.min_periods)
            else:

                def calc(x):
                    return func(x, window, min_periods=self.min_periods)

            if values.ndim > 1:
                result = np.apply_along_axis(calc, self.axis, values)
            else:
                result = calc(values)

            if center:
                result = self._center_window(result, window)

            results.append(result)

        return self._wrap_results(results, blocks, obj)


class _Rolling_and_Expanding(_Rolling):

    _shared_docs['count'] = """%(name)s count of number of non-NaN
    observations inside provided window."""

    def count(self):
        obj = self._convert_freq()
        window = self._get_window()
        window = min(window, len(obj)) if not self.center else window

        blocks, obj = self._create_blocks(how=None)
        results = []
        for b in blocks:

            if com.needs_i8_conversion(b.values):
                result = b.notnull().astype(int)
            else:
                try:
                    result = np.isfinite(b).astype(float)
                except TypeError:
                    result = np.isfinite(b.astype(float)).astype(float)

                result[pd.isnull(result)] = 0

            result = self._constructor(result, window=window, min_periods=0,
                                       center=self.center).sum()
            results.append(result)

        return self._wrap_results(results, blocks, obj)

    _shared_docs['apply'] = dedent("""
    %(name)s function apply

    Parameters
    ----------
    func : function
        Must produce a single value from an ndarray input
        \*args and \*\*kwargs are passed to the function""")

    def apply(self, func, args=(), kwargs={}):
        # TODO: _level is unused?
        _level = kwargs.pop('_level', None)  # noqa
        window = self._get_window()
        offset = _offset(window, self.center)

        def f(arg, window, min_periods):
            minp = _use_window(min_periods, window)
            return algos.roll_generic(arg, window, minp, offset, func, args,
                                      kwargs)

        return self._apply(f, func, args=args, kwargs=kwargs,
                           center=False)

    def sum(self, **kwargs):
        return self._apply('roll_sum', 'sum', **kwargs)

    _shared_docs['max'] = dedent("""
    %(name)s maximum

    Parameters
    ----------
    how : string, default 'max' (DEPRECATED)
        Method for down- or re-sampling""")

    def max(self, how=None, **kwargs):
        if self.freq is not None and how is None:
            how = 'max'
        return self._apply('roll_max', 'max', how=how, **kwargs)

    _shared_docs['min'] = dedent("""
    %(name)s minimum

    Parameters
    ----------
    how : string, default 'min' (DEPRECATED)
        Method for down- or re-sampling""")

    def min(self, how=None, **kwargs):
        if self.freq is not None and how is None:
            how = 'min'
        return self._apply('roll_min', 'min', how=how, **kwargs)

    def mean(self, **kwargs):
        return self._apply('roll_mean', 'mean', **kwargs)

    _shared_docs['median'] = dedent("""
    %(name)s median

    Parameters
    ----------
    how : string, default 'median' (DEPRECATED)
        Method for down- or re-sampling""")

    def median(self, how=None, **kwargs):
        if self.freq is not None and how is None:
            how = 'median'
        return self._apply('roll_median_c', 'median', how=how, **kwargs)

    _shared_docs['std'] = dedent("""
    %(name)s standard deviation

    Parameters
    ----------
    ddof : int, default 1
        Delta Degrees of Freedom.  The divisor used in calculations
        is ``N - ddof``, where ``N`` represents the number of elements.""")

    def std(self, ddof=1, **kwargs):
        window = self._get_window()

        def f(arg, *args, **kwargs):
            minp = _require_min_periods(1)(self.min_periods, window)
            return _zsqrt(algos.roll_var(arg, window, minp, ddof))

        return self._apply(f, 'std', check_minp=_require_min_periods(1),
                           ddof=ddof, **kwargs)

    _shared_docs['var'] = dedent("""
    %(name)s variance

    Parameters
    ----------
    ddof : int, default 1
        Delta Degrees of Freedom.  The divisor used in calculations
        is ``N - ddof``, where ``N`` represents the number of elements.""")

    def var(self, ddof=1, **kwargs):
        return self._apply('roll_var', 'var',
                           check_minp=_require_min_periods(1), ddof=ddof,
                           **kwargs)

    _shared_docs['skew'] = """Unbiased %(name)s skewness"""

    def skew(self, **kwargs):
        return self._apply('roll_skew', 'skew',
                           check_minp=_require_min_periods(3), **kwargs)

    _shared_docs['kurt'] = """Unbiased %(name)s kurtosis"""

    def kurt(self, **kwargs):
        return self._apply('roll_kurt', 'kurt',
                           check_minp=_require_min_periods(4), **kwargs)

    _shared_docs['quantile'] = dedent("""
    %(name)s quantile

    Parameters
    ----------
    quantile : float
        0 <= quantile <= 1""")

    def quantile(self, quantile, **kwargs):
        window = self._get_window()

        def f(arg, *args, **kwargs):
            minp = _use_window(self.min_periods, window)
            return algos.roll_quantile(arg, window, minp, quantile)

        return self._apply(f, 'quantile', quantile=quantile,
                           **kwargs)

    _shared_docs['cov'] = dedent("""
    %(name)s sample covariance

    Parameters
    ----------
    other : Series, DataFrame, or ndarray, optional
        if not supplied then will default to self and produce pairwise output
    pairwise : bool, default None
        If False then only matching columns between self and other will be used
        and the output will be a DataFrame.
        If True then all pairwise combinations will be calculated and the
        output will be a Panel in the case of DataFrame inputs. In the case of
        missing elements, only complete pairwise observations will be used.
    ddof : int, default 1
        Delta Degrees of Freedom.  The divisor used in calculations
        is ``N - ddof``, where ``N`` represents the number of elements.""")

    def cov(self, other=None, pairwise=None, ddof=1, **kwargs):
        if other is None:
            other = self._selected_obj
            # only default unset
            pairwise = True if pairwise is None else pairwise
        other = self._shallow_copy(other)
        window = self._get_window(other)

        def _get_cov(X, Y):
            # GH #12373 : rolling functions error on float32 data
            # to avoid potential overflow, cast the data to float64
            X = X.astype('float64')
            Y = Y.astype('float64')
            mean = lambda x: x.rolling(window, self.min_periods,
                                       center=self.center).mean(**kwargs)
            count = (X + Y).rolling(window=window,
                                    center=self.center).count(**kwargs)
            bias_adj = count / (count - ddof)
            return (mean(X * Y) - mean(X) * mean(Y)) * bias_adj

        return _flex_binary_moment(self._selected_obj, other._selected_obj,
                                   _get_cov, pairwise=bool(pairwise))

    _shared_docs['corr'] = dedent("""
    %(name)s sample correlation

    Parameters
    ----------
    other : Series, DataFrame, or ndarray, optional
        if not supplied then will default to self and produce pairwise output
    pairwise : bool, default None
        If False then only matching columns between self and other will be used
        and the output will be a DataFrame.
        If True then all pairwise combinations will be calculated and the
        output will be a Panel in the case of DataFrame inputs. In the case of
        missing elements, only complete pairwise observations will be used.""")

    def corr(self, other=None, pairwise=None, **kwargs):
        if other is None:
            other = self._selected_obj
            # only default unset
            pairwise = True if pairwise is None else pairwise
        other = self._shallow_copy(other)
        window = self._get_window(other)

        def _get_corr(a, b):
            a = a.rolling(window=window, min_periods=self.min_periods,
                          freq=self.freq, center=self.center)
            b = b.rolling(window=window, min_periods=self.min_periods,
                          freq=self.freq, center=self.center)

            return a.cov(b, **kwargs) / (a.std(**kwargs) * b.std(**kwargs))

        return _flex_binary_moment(self._selected_obj, other._selected_obj,
                                   _get_corr, pairwise=bool(pairwise))


class Rolling(_Rolling_and_Expanding):
    """
    Provides rolling window calculcations.

    .. versionadded:: 0.18.0

    Parameters
    ----------
    window : int
        Size of the moving window. This is the number of observations used for
        calculating the statistic.
    min_periods : int, default None
        Minimum number of observations in window required to have a value
        (otherwise result is NA).
    freq : string or DateOffset object, optional (default None) (DEPRECATED)
        Frequency to conform the data to before computing the statistic.
        Specified as a frequency string or DateOffset object.
    center : boolean, default False
        Set the labels at the center of the window.
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

    def validate(self):
        super(Rolling, self).validate()
        if not com.is_integer(self.window):
            raise ValueError("window must be an integer")

    @Substitution(name='rolling')
    @Appender(SelectionMixin._see_also_template)
    @Appender(SelectionMixin._agg_doc)
    def aggregate(self, arg, *args, **kwargs):
        return super(Rolling, self).aggregate(arg, *args, **kwargs)

    agg = aggregate

    @Substitution(name='rolling')
    @Appender(_doc_template)
    @Appender(_shared_docs['count'])
    def count(self):
        return super(Rolling, self).count()

    @Substitution(name='rolling')
    @Appender(_doc_template)
    @Appender(_shared_docs['apply'])
    def apply(self, func, args=(), kwargs={}):
        return super(Rolling, self).apply(func, args=args, kwargs=kwargs)

    @Substitution(name='rolling')
    @Appender(_doc_template)
    @Appender(_shared_docs['sum'])
    def sum(self, **kwargs):
        return super(Rolling, self).sum(**kwargs)

    @Substitution(name='rolling')
    @Appender(_doc_template)
    @Appender(_shared_docs['max'])
    def max(self, **kwargs):
        return super(Rolling, self).max(**kwargs)

    @Substitution(name='rolling')
    @Appender(_doc_template)
    @Appender(_shared_docs['min'])
    def min(self, **kwargs):
        return super(Rolling, self).min(**kwargs)

    @Substitution(name='rolling')
    @Appender(_doc_template)
    @Appender(_shared_docs['mean'])
    def mean(self, **kwargs):
        return super(Rolling, self).mean(**kwargs)

    @Substitution(name='rolling')
    @Appender(_doc_template)
    @Appender(_shared_docs['median'])
    def median(self, **kwargs):
        return super(Rolling, self).median(**kwargs)

    @Substitution(name='rolling')
    @Appender(_doc_template)
    @Appender(_shared_docs['std'])
    def std(self, ddof=1, **kwargs):
        return super(Rolling, self).std(ddof=ddof, **kwargs)

    @Substitution(name='rolling')
    @Appender(_doc_template)
    @Appender(_shared_docs['var'])
    def var(self, ddof=1, **kwargs):
        return super(Rolling, self).var(ddof=ddof, **kwargs)

    @Substitution(name='rolling')
    @Appender(_doc_template)
    @Appender(_shared_docs['skew'])
    def skew(self, **kwargs):
        return super(Rolling, self).skew(**kwargs)

    @Substitution(name='rolling')
    @Appender(_doc_template)
    @Appender(_shared_docs['kurt'])
    def kurt(self, **kwargs):
        return super(Rolling, self).kurt(**kwargs)

    @Substitution(name='rolling')
    @Appender(_doc_template)
    @Appender(_shared_docs['quantile'])
    def quantile(self, quantile, **kwargs):
        return super(Rolling, self).quantile(quantile=quantile, **kwargs)

    @Substitution(name='rolling')
    @Appender(_doc_template)
    @Appender(_shared_docs['cov'])
    def cov(self, other=None, pairwise=None, ddof=1, **kwargs):
        return super(Rolling, self).cov(other=other, pairwise=pairwise,
                                        ddof=ddof, **kwargs)

    @Substitution(name='rolling')
    @Appender(_doc_template)
    @Appender(_shared_docs['corr'])
    def corr(self, other=None, pairwise=None, **kwargs):
        return super(Rolling, self).corr(other=other, pairwise=pairwise,
                                         **kwargs)


class RollingGroupby(_GroupByMixin, Rolling):
    """
    Provides a rolling groupby implementation

    .. versionadded:: 0.18.1

    """
    @property
    def _constructor(self):
        return Rolling


class Expanding(_Rolling_and_Expanding):
    """
    Provides expanding transformations.

    .. versionadded:: 0.18.0

    Parameters
    ----------
    min_periods : int, default None
        Minimum number of observations in window required to have a value
        (otherwise result is NA).
    freq : string or DateOffset object, optional (default None) (DEPRECATED)
        Frequency to conform the data to before computing the statistic.
        Specified as a frequency string or DateOffset object.
    center : boolean, default False
        Set the labels at the center of the window.
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

    _attributes = ['min_periods', 'freq', 'center', 'axis']

    def __init__(self, obj, min_periods=1, freq=None, center=False, axis=0,
                 **kwargs):
        return super(Expanding, self).__init__(obj=obj,
                                               min_periods=min_periods,
                                               freq=freq, center=center,
                                               axis=axis)

    @property
    def _constructor(self):
        return Expanding

    def _get_window(self, other=None):
        obj = self._selected_obj
        if other is None:
            return (max(len(obj), self.min_periods) if self.min_periods
                    else len(obj))
        return (max((len(obj) + len(obj)), self.min_periods)
                if self.min_periods else (len(obj) + len(obj)))

    @Substitution(name='expanding')
    @Appender(SelectionMixin._see_also_template)
    @Appender(SelectionMixin._agg_doc)
    def aggregate(self, arg, *args, **kwargs):
        return super(Expanding, self).aggregate(arg, *args, **kwargs)

    agg = aggregate

    @Substitution(name='expanding')
    @Appender(_doc_template)
    @Appender(_shared_docs['count'])
    def count(self, **kwargs):
        return super(Expanding, self).count(**kwargs)

    @Substitution(name='expanding')
    @Appender(_doc_template)
    @Appender(_shared_docs['apply'])
    def apply(self, func, args=(), kwargs={}):
        return super(Expanding, self).apply(func, args=args, kwargs=kwargs)

    @Substitution(name='expanding')
    @Appender(_doc_template)
    @Appender(_shared_docs['sum'])
    def sum(self, **kwargs):
        return super(Expanding, self).sum(**kwargs)

    @Substitution(name='expanding')
    @Appender(_doc_template)
    @Appender(_shared_docs['max'])
    def max(self, **kwargs):
        return super(Expanding, self).max(**kwargs)

    @Substitution(name='expanding')
    @Appender(_doc_template)
    @Appender(_shared_docs['min'])
    def min(self, **kwargs):
        return super(Expanding, self).min(**kwargs)

    @Substitution(name='expanding')
    @Appender(_doc_template)
    @Appender(_shared_docs['mean'])
    def mean(self, **kwargs):
        return super(Expanding, self).mean(**kwargs)

    @Substitution(name='expanding')
    @Appender(_doc_template)
    @Appender(_shared_docs['median'])
    def median(self, **kwargs):
        return super(Expanding, self).median(**kwargs)

    @Substitution(name='expanding')
    @Appender(_doc_template)
    @Appender(_shared_docs['std'])
    def std(self, ddof=1, **kwargs):
        return super(Expanding, self).std(ddof=ddof, **kwargs)

    @Substitution(name='expanding')
    @Appender(_doc_template)
    @Appender(_shared_docs['var'])
    def var(self, ddof=1, **kwargs):
        return super(Expanding, self).var(ddof=ddof, **kwargs)

    @Substitution(name='expanding')
    @Appender(_doc_template)
    @Appender(_shared_docs['skew'])
    def skew(self, **kwargs):
        return super(Expanding, self).skew(**kwargs)

    @Substitution(name='expanding')
    @Appender(_doc_template)
    @Appender(_shared_docs['kurt'])
    def kurt(self, **kwargs):
        return super(Expanding, self).kurt(**kwargs)

    @Substitution(name='expanding')
    @Appender(_doc_template)
    @Appender(_shared_docs['quantile'])
    def quantile(self, quantile, **kwargs):
        return super(Expanding, self).quantile(quantile=quantile, **kwargs)

    @Substitution(name='expanding')
    @Appender(_doc_template)
    @Appender(_shared_docs['cov'])
    def cov(self, other=None, pairwise=None, ddof=1, **kwargs):
        return super(Expanding, self).cov(other=other, pairwise=pairwise,
                                          ddof=ddof, **kwargs)

    @Substitution(name='expanding')
    @Appender(_doc_template)
    @Appender(_shared_docs['corr'])
    def corr(self, other=None, pairwise=None, **kwargs):
        return super(Expanding, self).corr(other=other, pairwise=pairwise,
                                           **kwargs)


class ExpandingGroupby(_GroupByMixin, Expanding):
    """
    Provides a expanding groupby implementation

    .. versionadded:: 0.18.1

    """
    @property
    def _constructor(self):
        return Expanding


_bias_template = """

Parameters
----------
bias : boolean, default False
    Use a standard estimation bias correction
"""

_pairwise_template = """

Parameters
----------
other : Series, DataFrame, or ndarray, optional
    if not supplied then will default to self and produce pairwise output
pairwise : bool, default None
    If False then only matching columns between self and other will be used and
    the output will be a DataFrame.
    If True then all pairwise combinations will be calculated and the output
    will be a Panel in the case of DataFrame inputs. In the case of missing
    elements, only complete pairwise observations will be used.
bias : boolean, default False
   Use a standard estimation bias correction
"""


class EWM(_Rolling):
    r"""
    Provides exponential weighted functions

    .. versionadded:: 0.18.0

    Parameters
    ----------
    com : float, optional
        Specify decay in terms of center of mass,
        :math:`\alpha = 1 / (1 + com),\text{ for } com \geq 0`
    span : float, optional
        Specify decay in terms of span,
        :math:`\alpha = 2 / (span + 1),\text{ for } span \geq 1`
    halflife : float, optional
        Specify decay in terms of half-life,
        :math:`\alpha = 1 - exp(log(0.5) / halflife),\text{ for } halflife > 0`
    alpha : float, optional
        Specify smoothing factor :math:`\alpha` directly,
        :math:`0 < \alpha \leq 1`

        .. versionadded:: 0.18.0

    min_periods : int, default 0
        Minimum number of observations in window required to have a value
        (otherwise result is NA).
    freq : None or string alias / date offset object, default=None (DEPRECATED)
        Frequency to conform to before computing statistic
    adjust : boolean, default True
        Divide by decaying adjustment factor in beginning periods to account
        for imbalance in relative weightings (viewing EWMA as a moving average)
    ignore_na : boolean, default False
        Ignore missing values when calculating weights;
        specify True to reproduce pre-0.15.0 behavior

    Returns
    -------
    a Window sub-classed for the particular operation

    Notes
    -----
    Exactly one of center of mass, span, half-life, and alpha must be provided.
    Allowed values and relationship between the parameters are specified in the
    parameter descriptions above; see the link at the end of this section for
    a detailed explanation.

    The `freq` keyword is used to conform time series data to a specified
    frequency by resampling the data. This is done with the default parameters
    of :meth:`~pandas.Series.resample` (i.e. using the `mean`).

    When adjust is True (default), weighted averages are calculated using
    weights (1-alpha)**(n-1), (1-alpha)**(n-2), ..., 1-alpha, 1.

    When adjust is False, weighted averages are calculated recursively as:
       weighted_average[0] = arg[0];
       weighted_average[i] = (1-alpha)*weighted_average[i-1] + alpha*arg[i].

    When ignore_na is False (default), weights are based on absolute positions.
    For example, the weights of x and y used in calculating the final weighted
    average of [x, None, y] are (1-alpha)**2 and 1 (if adjust is True), and
    (1-alpha)**2 and alpha (if adjust is False).

    When ignore_na is True (reproducing pre-0.15.0 behavior), weights are based
    on relative positions. For example, the weights of x and y used in
    calculating the final weighted average of [x, None, y] are 1-alpha and 1
    (if adjust is True), and 1-alpha and alpha (if adjust is False).

    More details can be found at
    http://pandas.pydata.org/pandas-docs/stable/computation.html#exponentially-weighted-windows
    """
    _attributes = ['com', 'min_periods', 'freq', 'adjust', 'ignore_na', 'axis']

    def __init__(self, obj, com=None, span=None, halflife=None, alpha=None,
                 min_periods=0, freq=None, adjust=True, ignore_na=False,
                 axis=0):
        self.obj = obj
        self.com = _get_center_of_mass(com, span, halflife, alpha)
        self.min_periods = min_periods
        self.freq = freq
        self.adjust = adjust
        self.ignore_na = ignore_na
        self.axis = axis

    @property
    def _constructor(self):
        return EWM

    @Substitution(name='ewm')
    @Appender(SelectionMixin._see_also_template)
    @Appender(SelectionMixin._agg_doc)
    def aggregate(self, arg, *args, **kwargs):
        return super(EWM, self).aggregate(arg, *args, **kwargs)

    agg = aggregate

    def _apply(self, func, how=None, **kwargs):
        """Rolling statistical measure using supplied function. Designed to be
        used with passed-in Cython array-based functions.

        Parameters
        ----------
        func : string/callable to apply
        how : string, default to None (DEPRECATED)
            how to resample

        Returns
        -------
        y : type of input argument

        """
        blocks, obj = self._create_blocks(how=how)
        results = []
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
                    raise ValueError("we do not support this function "
                                     "algos.{0}".format(func))

                cfunc = getattr(algos, func)

                def func(arg):
                    return cfunc(arg, self.com, int(self.adjust),
                                 int(self.ignore_na), int(self.min_periods))

            results.append(np.apply_along_axis(func, self.axis, values))

        return self._wrap_results(results, blocks, obj)

    @Substitution(name='ewm')
    @Appender(_doc_template)
    def mean(self, **kwargs):
        """exponential weighted moving average"""
        return self._apply('ewma', **kwargs)

    @Substitution(name='ewm')
    @Appender(_doc_template)
    @Appender(_bias_template)
    def std(self, bias=False, **kwargs):
        """exponential weighted moving stddev"""
        return _zsqrt(self.var(bias=bias, **kwargs))

    vol = std

    @Substitution(name='ewm')
    @Appender(_doc_template)
    @Appender(_bias_template)
    def var(self, bias=False, **kwargs):
        """exponential weighted moving variance"""

        def f(arg):
            return algos.ewmcov(arg, arg, self.com, int(self.adjust),
                                int(self.ignore_na), int(self.min_periods),
                                int(bias))

        return self._apply(f, **kwargs)

    @Substitution(name='ewm')
    @Appender(_doc_template)
    @Appender(_pairwise_template)
    def cov(self, other=None, pairwise=None, bias=False, **kwargs):
        """exponential weighted sample covariance"""
        if other is None:
            other = self._selected_obj
            # only default unset
            pairwise = True if pairwise is None else pairwise
        other = self._shallow_copy(other)

        def _get_cov(X, Y):
            X = self._shallow_copy(X)
            Y = self._shallow_copy(Y)
            cov = algos.ewmcov(X._prep_values(), Y._prep_values(), self.com,
                               int(self.adjust), int(self.ignore_na),
                               int(self.min_periods), int(bias))
            return X._wrap_result(cov)

        return _flex_binary_moment(self._selected_obj, other._selected_obj,
                                   _get_cov, pairwise=bool(pairwise))

    @Substitution(name='ewm')
    @Appender(_doc_template)
    @Appender(_pairwise_template)
    def corr(self, other=None, pairwise=None, **kwargs):
        """exponential weighted sample correlation"""
        if other is None:
            other = self._selected_obj
            # only default unset
            pairwise = True if pairwise is None else pairwise
        other = self._shallow_copy(other)

        def _get_corr(X, Y):
            X = self._shallow_copy(X)
            Y = self._shallow_copy(Y)

            def _cov(x, y):
                return algos.ewmcov(x, y, self.com, int(self.adjust),
                                    int(self.ignore_na), int(self.min_periods),
                                    1)

            x_values = X._prep_values()
            y_values = Y._prep_values()
            cov = _cov(x_values, y_values)
            x_var = _cov(x_values, x_values)
            y_var = _cov(y_values, y_values)
            corr = cov / _zsqrt(x_var * y_var)
            return X._wrap_result(corr)

        return _flex_binary_moment(self._selected_obj, other._selected_obj,
                                   _get_corr, pairwise=bool(pairwise))

# Helper Funcs


def _flex_binary_moment(arg1, arg2, f, pairwise=False):
    from pandas import Series, DataFrame, Panel
    if not (isinstance(arg1, (np.ndarray, Series, DataFrame)) and
            isinstance(arg2, (np.ndarray, Series, DataFrame))):
        raise TypeError("arguments to moment function must be of type "
                        "np.ndarray/Series/DataFrame")

    if (isinstance(arg1, (np.ndarray, Series)) and
            isinstance(arg2, (np.ndarray, Series))):
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
                    return DataFrame(results, index=X.index,
                                     columns=res_columns)
            elif pairwise is True:
                results = defaultdict(dict)
                for i, k1 in enumerate(arg1.columns):
                    for j, k2 in enumerate(arg2.columns):
                        if j < i and arg2 is arg1:
                            # Symmetric case
                            results[i][j] = results[j][i]
                        else:
                            results[i][j] = f(*_prep_binary(arg1.iloc[:, i],
                                                            arg2.iloc[:, j]))
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


def _get_center_of_mass(com, span, halflife, alpha):
    valid_count = len([x for x in [com, span, halflife, alpha]
                      if x is not None])
    if valid_count > 1:
        raise ValueError("com, span, halflife, and alpha "
                         "are mutually exclusive")

    # Convert to center of mass; domain checks ensure 0 < alpha <= 1
    if com is not None:
        if com < 0:
            raise ValueError("com must satisfy: com >= 0")
    elif span is not None:
        if span < 1:
            raise ValueError("span must satisfy: span >= 1")
        com = (span - 1) / 2.
    elif halflife is not None:
        if halflife <= 0:
            raise ValueError("halflife must satisfy: halflife > 0")
        decay = 1 - np.exp(np.log(0.5) / halflife)
        com = 1 / decay - 1
    elif alpha is not None:
        if alpha <= 0 or alpha > 1:
            raise ValueError("alpha must satisfy: 0 < alpha <= 1")
        com = (1.0 - alpha) / alpha
    else:
        raise ValueError("Must pass one of com, span, halflife, or alpha")

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


# Top-level exports


def rolling(obj, win_type=None, **kwds):
    from pandas import Series, DataFrame
    if not isinstance(obj, (Series, DataFrame)):
        raise TypeError('invalid type: %s' % type(obj))

    if win_type is not None:
        return Window(obj, win_type=win_type, **kwds)

    return Rolling(obj, **kwds)


rolling.__doc__ = Window.__doc__


def expanding(obj, **kwds):
    from pandas import Series, DataFrame
    if not isinstance(obj, (Series, DataFrame)):
        raise TypeError('invalid type: %s' % type(obj))

    return Expanding(obj, **kwds)


expanding.__doc__ = Expanding.__doc__


def ewm(obj, **kwds):
    from pandas import Series, DataFrame
    if not isinstance(obj, (Series, DataFrame)):
        raise TypeError('invalid type: %s' % type(obj))

    return EWM(obj, **kwds)


ewm.__doc__ = EWM.__doc__
