# -*- coding: utf-8 -*-
"""
Base and utility classes for tseries type pandas objects.
"""
import operator
import warnings

import numpy as np

from pandas._libs import NaT, iNaT, lib
from pandas.compat.numpy import function as nv
from pandas.errors import AbstractMethodError
from pandas.util._decorators import Appender, cache_readonly

from pandas.core.dtypes.common import (
    ensure_int64, is_bool_dtype, is_categorical_dtype,
    is_datetime_or_timedelta_dtype, is_dtype_equal, is_float, is_float_dtype,
    is_integer, is_integer_dtype, is_list_like, is_object_dtype,
    is_period_dtype, is_scalar, is_string_dtype)
from pandas.core.dtypes.generic import ABCIndex, ABCIndexClass, ABCSeries

from pandas.core import algorithms, ops
from pandas.core.accessor import PandasDelegate
from pandas.core.arrays.datetimelike import (
    DatetimeLikeArrayMixin, _ensure_datetimelike_to_i8)
import pandas.core.indexes.base as ibase
from pandas.core.indexes.base import Index, _index_shared_docs
from pandas.core.tools.timedeltas import to_timedelta

import pandas.io.formats.printing as printing

_index_doc_kwargs = dict(ibase._index_doc_kwargs)


class DatetimeIndexOpsMixin(DatetimeLikeArrayMixin):
    """
    common ops mixin to support a unified interface datetimelike Index
    """

    # override DatetimeLikeArrayMixin method
    copy = Index.copy
    unique = Index.unique
    take = Index.take

    # DatetimeLikeArrayMixin assumes subclasses are mutable, so these are
    # properties there.  They can be made into cache_readonly for Index
    # subclasses bc they are immutable
    inferred_freq = cache_readonly(DatetimeLikeArrayMixin.inferred_freq.fget)
    _isnan = cache_readonly(DatetimeLikeArrayMixin._isnan.fget)
    hasnans = cache_readonly(DatetimeLikeArrayMixin.hasnans.fget)
    _resolution = cache_readonly(DatetimeLikeArrayMixin._resolution.fget)
    resolution = cache_readonly(DatetimeLikeArrayMixin.resolution.fget)

    def equals(self, other):
        """
        Determines if two Index objects contain the same elements.
        """
        if self.is_(other):
            return True

        if not isinstance(other, ABCIndexClass):
            return False
        elif not isinstance(other, type(self)):
            try:
                other = type(self)(other)
            except Exception:
                return False

        if not is_dtype_equal(self.dtype, other.dtype):
            # have different timezone
            return False

        elif is_period_dtype(self):
            if not is_period_dtype(other):
                return False
            if self.freq != other.freq:
                return False

        return np.array_equal(self.asi8, other.asi8)

    @staticmethod
    def _join_i8_wrapper(joinf, dtype, with_indexers=True):
        """
        Create the join wrapper methods.
        """

        @staticmethod
        def wrapper(left, right):
            if isinstance(left, (np.ndarray, ABCIndex, ABCSeries)):
                left = left.view('i8')
            if isinstance(right, (np.ndarray, ABCIndex, ABCSeries)):
                right = right.view('i8')
            results = joinf(left, right)
            if with_indexers:
                join_index, left_indexer, right_indexer = results
                join_index = join_index.view(dtype)
                return join_index, left_indexer, right_indexer
            return results

        return wrapper

    @Appender(DatetimeLikeArrayMixin._evaluate_compare.__doc__)
    def _evaluate_compare(self, other, op):
        result = DatetimeLikeArrayMixin._evaluate_compare(self, other, op)
        if is_bool_dtype(result):
            return result
        try:
            return Index(result)
        except TypeError:
            return result

    def _ensure_localized(self, arg, ambiguous='raise', nonexistent='raise',
                          from_utc=False):
        """
        Ensure that we are re-localized.

        This is for compat as we can then call this on all datetimelike
        indexes generally (ignored for Period/Timedelta)

        Parameters
        ----------
        arg : DatetimeIndex / i8 ndarray
        ambiguous : str, bool, or bool-ndarray, default 'raise'
        nonexistent : str, default 'raise'
        from_utc : bool, default False
            If True, localize the i8 ndarray to UTC first before converting to
            the appropriate tz. If False, localize directly to the tz.

        Returns
        -------
        localized DTI
        """

        # reconvert to local tz
        if getattr(self, 'tz', None) is not None:
            if not isinstance(arg, ABCIndexClass):
                arg = self._simple_new(arg)
            if from_utc:
                arg = arg.tz_localize('UTC').tz_convert(self.tz)
            else:
                arg = arg.tz_localize(
                    self.tz, ambiguous=ambiguous, nonexistent=nonexistent
                )
        return arg

    def _box_values_as_index(self):
        """
        Return object Index which contains boxed values.
        """
        from pandas.core.index import Index
        return Index(self._box_values(self.asi8), name=self.name, dtype=object)

    @Appender(_index_shared_docs['contains'] % _index_doc_kwargs)
    def __contains__(self, key):
        try:
            res = self.get_loc(key)
            return (is_scalar(res) or isinstance(res, slice) or
                    (is_list_like(res) and len(res)))
        except (KeyError, TypeError, ValueError):
            return False

    contains = __contains__

    # Try to run function on index first, and then on elements of index
    # Especially important for group-by functionality
    def map(self, f):
        try:
            result = f(self)

            # Try to use this result if we can
            if isinstance(result, np.ndarray):
                result = Index(result)

            if not isinstance(result, Index):
                raise TypeError('The map function must return an Index object')
            return result
        except Exception:
            return self.astype(object).map(f)

    def sort_values(self, return_indexer=False, ascending=True):
        """
        Return sorted copy of Index.
        """
        if return_indexer:
            _as = self.argsort()
            if not ascending:
                _as = _as[::-1]
            sorted_index = self.take(_as)
            return sorted_index, _as
        else:
            sorted_values = np.sort(self._ndarray_values)
            attribs = self._get_attributes_dict()
            freq = attribs['freq']

            if freq is not None and not is_period_dtype(self):
                if freq.n > 0 and not ascending:
                    freq = freq * -1
                elif freq.n < 0 and ascending:
                    freq = freq * -1
            attribs['freq'] = freq

            if not ascending:
                sorted_values = sorted_values[::-1]

            sorted_values = self._maybe_box_as_values(sorted_values,
                                                      **attribs)

            return self._simple_new(sorted_values, **attribs)

    @Appender(_index_shared_docs['take'] % _index_doc_kwargs)
    def take(self, indices, axis=0, allow_fill=True,
             fill_value=None, **kwargs):
        nv.validate_take(tuple(), kwargs)
        indices = ensure_int64(indices)

        maybe_slice = lib.maybe_indices_to_slice(indices, len(self))
        if isinstance(maybe_slice, slice):
            return self[maybe_slice]

        taken = self._assert_take_fillable(self.asi8, indices,
                                           allow_fill=allow_fill,
                                           fill_value=fill_value,
                                           na_value=iNaT)

        # keep freq in PeriodArray/Index, reset otherwise
        freq = self.freq if is_period_dtype(self) else None
        return self._shallow_copy(taken, freq=freq)

    _can_hold_na = True

    _na_value = NaT
    """The expected NA value to use with this index."""

    @property
    def asobject(self):
        """
        Return object Index which contains boxed values.

        .. deprecated:: 0.23.0
            Use ``astype(object)`` instead.

        *this is an internal non-public method*
        """
        warnings.warn("'asobject' is deprecated. Use 'astype(object)'"
                      " instead", FutureWarning, stacklevel=2)
        return self.astype(object)

    def _convert_tolerance(self, tolerance, target):
        tolerance = np.asarray(to_timedelta(tolerance, box=False))
        if target.size != tolerance.size and tolerance.size > 1:
            raise ValueError('list-like tolerance size must match '
                             'target index size')
        return tolerance

    def tolist(self):
        """
        Return a list of the underlying data.
        """
        return list(self.astype(object))

    def min(self, axis=None, *args, **kwargs):
        """
        Return the minimum value of the Index or minimum along
        an axis.

        See Also
        --------
        numpy.ndarray.min
        """
        nv.validate_min(args, kwargs)
        nv.validate_minmax_axis(axis)

        try:
            i8 = self.asi8

            # quick check
            if len(i8) and self.is_monotonic:
                if i8[0] != iNaT:
                    return self._box_func(i8[0])

            if self.hasnans:
                min_stamp = self[~self._isnan].asi8.min()
            else:
                min_stamp = i8.min()
            return self._box_func(min_stamp)
        except ValueError:
            return self._na_value

    def argmin(self, axis=None, *args, **kwargs):
        """
        Returns the indices of the minimum values along an axis.

        See `numpy.ndarray.argmin` for more information on the
        `axis` parameter.

        See Also
        --------
        numpy.ndarray.argmin
        """
        nv.validate_argmin(args, kwargs)
        nv.validate_minmax_axis(axis)

        i8 = self.asi8
        if self.hasnans:
            mask = self._isnan
            if mask.all():
                return -1
            i8 = i8.copy()
            i8[mask] = np.iinfo('int64').max
        return i8.argmin()

    def max(self, axis=None, *args, **kwargs):
        """
        Return the maximum value of the Index or maximum along
        an axis.

        See Also
        --------
        numpy.ndarray.max
        """
        nv.validate_max(args, kwargs)
        nv.validate_minmax_axis(axis)

        try:
            i8 = self.asi8

            # quick check
            if len(i8) and self.is_monotonic:
                if i8[-1] != iNaT:
                    return self._box_func(i8[-1])

            if self.hasnans:
                max_stamp = self[~self._isnan].asi8.max()
            else:
                max_stamp = i8.max()
            return self._box_func(max_stamp)
        except ValueError:
            return self._na_value

    def argmax(self, axis=None, *args, **kwargs):
        """
        Returns the indices of the maximum values along an axis.

        See `numpy.ndarray.argmax` for more information on the
        `axis` parameter.

        See Also
        --------
        numpy.ndarray.argmax
        """
        nv.validate_argmax(args, kwargs)
        nv.validate_minmax_axis(axis)

        i8 = self.asi8
        if self.hasnans:
            mask = self._isnan
            if mask.all():
                return -1
            i8 = i8.copy()
            i8[mask] = 0
        return i8.argmax()

    # --------------------------------------------------------------------
    # Rendering Methods

    def _format_with_header(self, header, **kwargs):
        return header + list(self._format_native_types(**kwargs))

    @property
    def _formatter_func(self):
        raise AbstractMethodError(self)

    def _format_attrs(self):
        """
        Return a list of tuples of the (attr,formatted_value).
        """
        attrs = super(DatetimeIndexOpsMixin, self)._format_attrs()
        for attrib in self._attributes:
            if attrib == 'freq':
                freq = self.freqstr
                if freq is not None:
                    freq = "'%s'" % freq
                attrs.append(('freq', freq))
        return attrs

    # --------------------------------------------------------------------

    def _convert_scalar_indexer(self, key, kind=None):
        """
        We don't allow integer or float indexing on datetime-like when using
        loc.

        Parameters
        ----------
        key : label of the slice bound
        kind : {'ix', 'loc', 'getitem', 'iloc'} or None
        """

        assert kind in ['ix', 'loc', 'getitem', 'iloc', None]

        # we don't allow integer/float indexing for loc
        # we don't allow float indexing for ix/getitem
        if is_scalar(key):
            is_int = is_integer(key)
            is_flt = is_float(key)
            if kind in ['loc'] and (is_int or is_flt):
                self._invalid_indexer('index', key)
            elif kind in ['ix', 'getitem'] and is_flt:
                self._invalid_indexer('index', key)

        return (super(DatetimeIndexOpsMixin, self)
                ._convert_scalar_indexer(key, kind=kind))

    @classmethod
    def _add_datetimelike_methods(cls):
        """
        Add in the datetimelike methods (as we may have to override the
        superclass).
        """

        def __add__(self, other):
            # dispatch to ExtensionArray implementation
            result = super(cls, self).__add__(other)
            return wrap_arithmetic_op(self, other, result)

        cls.__add__ = __add__

        def __radd__(self, other):
            # alias for __add__
            return self.__add__(other)
        cls.__radd__ = __radd__

        def __sub__(self, other):
            # dispatch to ExtensionArray implementation
            result = super(cls, self).__sub__(other)
            return wrap_arithmetic_op(self, other, result)

        cls.__sub__ = __sub__

        def __rsub__(self, other):
            result = super(cls, self).__rsub__(other)
            return wrap_arithmetic_op(self, other, result)

        cls.__rsub__ = __rsub__

    def isin(self, values):
        """
        Compute boolean array of whether each index value is found in the
        passed set of values.

        Parameters
        ----------
        values : set or sequence of values

        Returns
        -------
        is_contained : ndarray (boolean dtype)
        """
        if not isinstance(values, type(self)):
            try:
                values = type(self)(values)
            except ValueError:
                return self.astype(object).isin(values)

        return algorithms.isin(self.asi8, values.asi8)

    def repeat(self, repeats, *args, **kwargs):
        """
        Analogous to ndarray.repeat.
        """
        nv.validate_repeat(args, kwargs)
        if is_period_dtype(self):
            freq = self.freq
        else:
            freq = None
        return self._shallow_copy(self.asi8.repeat(repeats),
                                  freq=freq)

    @Appender(_index_shared_docs['where'] % _index_doc_kwargs)
    def where(self, cond, other=None):
        other = _ensure_datetimelike_to_i8(other, to_utc=True)
        values = _ensure_datetimelike_to_i8(self, to_utc=True)
        result = np.where(cond, values, other).astype('i8')

        result = self._ensure_localized(result, from_utc=True)
        return self._shallow_copy(result)

    def _summary(self, name=None):
        """
        Return a summarized representation.

        Parameters
        ----------
        name : str
            name to use in the summary representation

        Returns
        -------
        String with a summarized representation of the index
        """
        formatter = self._formatter_func
        if len(self) > 0:
            index_summary = ', %s to %s' % (formatter(self[0]),
                                            formatter(self[-1]))
        else:
            index_summary = ''

        if name is None:
            name = type(self).__name__
        result = '%s: %s entries%s' % (printing.pprint_thing(name),
                                       len(self), index_summary)
        if self.freq:
            result += '\nFreq: %s' % self.freqstr

        # display as values, not quoted
        result = result.replace("'", "")
        return result

    def _concat_same_dtype(self, to_concat, name):
        """
        Concatenate to_concat which has the same class.
        """
        attribs = self._get_attributes_dict()
        attribs['name'] = name
        # do not pass tz to set because tzlocal cannot be hashed
        if len({str(x.dtype) for x in to_concat}) != 1:
            raise ValueError('to_concat must have the same tz')

        if not is_period_dtype(self):
            # reset freq
            attribs['freq'] = None
            # TODO(DatetimeArray)
            # - remove the .asi8 here
            # - remove the _maybe_box_as_values
            # - combine with the `else` block
            new_data = self._concat_same_type(to_concat).asi8
        else:
            new_data = type(self._values)._concat_same_type(to_concat)

        return self._simple_new(new_data, **attribs)

    def _maybe_box_as_values(self, values, **attribs):
        # TODO(DatetimeArray): remove
        # This is a temporary shim while PeriodArray is an ExtensoinArray,
        # but others are not. When everyone is an ExtensionArray, this can
        # be removed. Currently used in
        # - sort_values
        return values

    def astype(self, dtype, copy=True):
        if is_object_dtype(dtype):
            return self._box_values_as_index()
        elif is_string_dtype(dtype) and not is_categorical_dtype(dtype):
            return Index(self.format(), name=self.name, dtype=object)
        elif is_integer_dtype(dtype):
            # TODO(DatetimeArray): use self._values here.
            # Can't use ._values currently, because that returns a
            # DatetimeIndex, which throws us in an infinite loop.
            return Index(self.values.astype('i8', copy=copy), name=self.name,
                         dtype='i8')
        elif (is_datetime_or_timedelta_dtype(dtype) and
              not is_dtype_equal(self.dtype, dtype)) or is_float_dtype(dtype):
            # disallow conversion between datetime/timedelta,
            # and conversions for any datetimelike to float
            msg = 'Cannot cast {name} to dtype {dtype}'
            raise TypeError(msg.format(name=type(self).__name__, dtype=dtype))
        return super(DatetimeIndexOpsMixin, self).astype(dtype, copy=copy)

    @Appender(DatetimeLikeArrayMixin._time_shift.__doc__)
    def _time_shift(self, periods, freq=None):
        result = DatetimeLikeArrayMixin._time_shift(self, periods, freq=freq)
        result.name = self.name
        return result


def wrap_arithmetic_op(self, other, result):
    if result is NotImplemented:
        return NotImplemented

    if isinstance(result, tuple):
        # divmod, rdivmod
        assert len(result) == 2
        return (wrap_arithmetic_op(self, other, result[0]),
                wrap_arithmetic_op(self, other, result[1]))

    if not isinstance(result, Index):
        # Index.__new__ will choose appropriate subclass for dtype
        result = Index(result)

    res_name = ops.get_op_result_name(self, other)
    result.name = res_name
    return result


def wrap_array_method(method, pin_name=False):
    """
    Wrap a DatetimeArray/TimedeltaArray/PeriodArray method so that the
    returned object is an Index subclass instead of ndarray or ExtensionArray
    subclass.

    Parameters
    ----------
    method : method of Datetime/Timedelta/Period Array class
    pin_name : bool
        Whether to set name=self.name on the output Index

    Returns
    -------
    method
    """
    def index_method(self, *args, **kwargs):
        result = method(self, *args, **kwargs)

        # Index.__new__ will choose the appropriate subclass to return
        result = Index(result)
        if pin_name:
            result.name = self.name
        return result

    index_method.__name__ = method.__name__
    index_method.__doc__ = method.__doc__
    return index_method


def wrap_field_accessor(prop):
    """
    Wrap a DatetimeArray/TimedeltaArray/PeriodArray array-returning property
    to return an Index subclass instead of ndarray or ExtensionArray subclass.

    Parameters
    ----------
    prop : property

    Returns
    -------
    new_prop : property
    """
    fget = prop.fget

    def f(self):
        result = fget(self)
        if is_bool_dtype(result):
            # return numpy array b/c there is no BoolIndex
            return result
        return Index(result, name=self.name)

    f.__name__ = fget.__name__
    f.__doc__ = fget.__doc__
    return property(f)


class DatetimelikeDelegateMixin(PandasDelegate):
    """
    Delegation mechanism, specific for Datetime, Timedelta, and Period types.

    Functionality is delegated from the Index class to an Array class. A
    few things can be customized

    * _delegate_class : type
        The class being delegated to.
    * _delegated_methods, delegated_properties : List
        The list of property / method names being delagated.
    * raw_methods : Set
        The set of methods whose results should should *not* be
        boxed in an index, after being returned from the array
    * raw_properties : Set
        The set of properties whose results should should *not* be
        boxed in an index, after being returned from the array
    """
    # raw_methods : dispatch methods that shouldn't be boxed in an Index
    _raw_methods = set()
    # raw_properties : dispatch properties that shouldn't be boxed in an Index
    _raw_properties = set()
    name = None
    _data = None

    @property
    def _delegate_class(self):
        raise AbstractMethodError

    def _delegate_property_get(self, name, *args, **kwargs):
        result = getattr(self._data, name)
        if name not in self._raw_properties:
            result = Index(result, name=self.name)
        return result

    def _delegate_property_set(self, name, value, *args, **kwargs):
        setattr(self._data, name, value)

    def _delegate_method(self, name, *args, **kwargs):
        result = operator.methodcaller(name, *args, **kwargs)(self._data)
        if name not in self._raw_methods:
            result = Index(result, name=self.name)
        return result
