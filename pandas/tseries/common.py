## datetimelike delegation ##

import numpy as np
from pandas.core.base import PandasDelegate
from pandas.core import common as com
from pandas.tseries.index import DatetimeIndex
from pandas.tseries.period import PeriodIndex
from pandas.tseries.tdi import TimedeltaIndex
from pandas import tslib
from pandas.core.common import (_NS_DTYPE, _TD_DTYPE, is_period_arraylike,
                                is_datetime_arraylike, is_integer_dtype, is_list_like,
                                is_datetime64_dtype, is_datetime64tz_dtype,
                                is_timedelta64_dtype,
                                get_dtype_kinds)

def is_datetimelike(data):
    """ return a boolean if we can be successfully converted to a datetimelike """
    try:
        maybe_to_datetimelike(data)
        return True
    except (Exception):
        pass
    return False

def maybe_to_datetimelike(data, copy=False):
    """
    return a DelegatedClass of a Series that is datetimelike
      (e.g. datetime64[ns],timedelta64[ns] dtype or a Series of Periods)
    raise TypeError if this is not possible.

    Parameters
    ----------
    data : Series
    copy : boolean, default False
           copy the input data

    Returns
    -------
    DelegatedClass

    """
    from pandas import Series

    if not isinstance(data, Series):
        raise TypeError("cannot convert an object of type {0} to a datetimelike index".format(type(data)))

    index = data.index
    if is_datetime64_dtype(data.dtype) or is_datetime64tz_dtype(data.dtype):
        return DatetimeProperties(DatetimeIndex(data, copy=copy, freq='infer'), index, name=data.name)
    elif is_timedelta64_dtype(data.dtype):
        return TimedeltaProperties(TimedeltaIndex(data, copy=copy, freq='infer'), index, name=data.name)
    else:
        if is_period_arraylike(data):
            return PeriodProperties(PeriodIndex(data, copy=copy), index, name=data.name)
        if is_datetime_arraylike(data):
            return DatetimeProperties(DatetimeIndex(data, copy=copy, freq='infer'), index, name=data.name)

    raise TypeError("cannot convert an object of type {0} to a datetimelike index".format(type(data)))

class Properties(PandasDelegate):

    def __init__(self, values, index, name):
        self.values = values
        self.index = index
        self.name = name

    def _delegate_property_get(self, name):
        from pandas import Series

        result = getattr(self.values,name)

        # maybe need to upcast (ints)
        if isinstance(result, np.ndarray):
            if is_integer_dtype(result):
                result = result.astype('int64')
        elif not is_list_like(result):
            return result

        # return the result as a Series, which is by definition a copy
        result = Series(result, index=self.index, name=self.name)

        # setting this object will show a SettingWithCopyWarning/Error
        result.is_copy = ("modifications to a property of a datetimelike object are not "
                          "supported and are discarded. Change values on the original.")

        return result

    def _delegate_property_set(self, name, value, *args, **kwargs):
        raise ValueError("modifications to a property of a datetimelike object are not "
                         "supported. Change values on the original.")

    def _delegate_method(self, name, *args, **kwargs):
        from pandas import Series

        method = getattr(self.values, name)
        result = method(*args, **kwargs)

        if not com.is_list_like(result):
            return result

        result = Series(result, index=self.index, name=self.name)

        # setting this object will show a SettingWithCopyWarning/Error
        result.is_copy = ("modifications to a method of a datetimelike object are not "
                          "supported and are discarded. Change values on the original.")

        return result


class DatetimeProperties(Properties):
    """
    Accessor object for datetimelike properties of the Series values.

    Examples
    --------
    >>> s.dt.hour
    >>> s.dt.second
    >>> s.dt.quarter

    Returns a Series indexed like the original Series.
    Raises TypeError if the Series does not contain datetimelike values.
    """

    def to_pydatetime(self):
        return self.values.to_pydatetime()

DatetimeProperties._add_delegate_accessors(delegate=DatetimeIndex,
                                           accessors=DatetimeIndex._datetimelike_ops,
                                           typ='property')
DatetimeProperties._add_delegate_accessors(delegate=DatetimeIndex,
                                           accessors=["to_period","tz_localize","tz_convert","normalize","strftime"],
                                           typ='method')

class TimedeltaProperties(Properties):
    """
    Accessor object for datetimelike properties of the Series values.

    Examples
    --------
    >>> s.dt.hours
    >>> s.dt.seconds

    Returns a Series indexed like the original Series.
    Raises TypeError if the Series does not contain datetimelike values.
    """

    def to_pytimedelta(self):
        return self.values.to_pytimedelta()

    @property
    def components(self):
        """
        Return a dataframe of the components (days, hours, minutes,
        seconds, milliseconds, microseconds, nanoseconds) of the Timedeltas.

        Returns
        -------
        a DataFrame

        """
        return self.values.components.set_index(self.index)

TimedeltaProperties._add_delegate_accessors(delegate=TimedeltaIndex,
                                            accessors=TimedeltaIndex._datetimelike_ops,
                                            typ='property')
TimedeltaProperties._add_delegate_accessors(delegate=TimedeltaIndex,
                                            accessors=["to_pytimedelta", "total_seconds"],
                                            typ='method')

class PeriodProperties(Properties):
    """
    Accessor object for datetimelike properties of the Series values.

    Examples
    --------
    >>> s.dt.hour
    >>> s.dt.second
    >>> s.dt.quarter

    Returns a Series indexed like the original Series.
    Raises TypeError if the Series does not contain datetimelike values.
    """

PeriodProperties._add_delegate_accessors(delegate=PeriodIndex,
                                         accessors=PeriodIndex._datetimelike_ops,
                                         typ='property')
PeriodProperties._add_delegate_accessors(delegate=PeriodIndex,
                                         accessors=["strftime"],
                                         typ='method')


class CombinedDatetimelikeProperties(DatetimeProperties, TimedeltaProperties):
    # This class is never instantiated, and exists solely for the benefit of
    # the Series.dt class property. For Series objects, .dt will always be one
    # of the more specific classes above.
    __doc__ = DatetimeProperties.__doc__


def _concat_compat(to_concat, axis=0):
    """
    provide concatenation of an datetimelike array of arrays each of which is a single
    M8[ns], datetimet64[ns, tz] or m8[ns] dtype

    Parameters
    ----------
    to_concat : array of arrays
    axis : axis to provide concatenation

    Returns
    -------
    a single array, preserving the combined dtypes
    """

    def convert_to_pydatetime(x, axis):
        # coerce to an object dtype
        if x.dtype == _NS_DTYPE:

            if hasattr(x, 'tz'):
                x = x.asobject

            shape = x.shape
            x = tslib.ints_to_pydatetime(x.view(np.int64).ravel())
            x = x.reshape(shape)
        elif x.dtype == _TD_DTYPE:
            shape = x.shape
            x = tslib.ints_to_pytimedelta(x.view(np.int64).ravel())
            x = x.reshape(shape)

        return x

    typs = get_dtype_kinds(to_concat)

    # datetimetz
    if 'datetimetz' in typs:

        # we require ALL of the same tz for datetimetz
        tzs = set([ getattr(x,'tz',None) for x in to_concat ])-set([None])
        if len(tzs) == 1:
            return DatetimeIndex(np.concatenate([ x.tz_localize(None).asi8 for x in to_concat ]), tz=list(tzs)[0])

    # single dtype
    if len(typs) == 1:

        if not len(typs-set(['datetime'])):
            new_values = np.concatenate([x.view(np.int64) for x in to_concat],
                                        axis=axis)
            return new_values.view(_NS_DTYPE)

        elif not len(typs-set(['timedelta'])):
            new_values = np.concatenate([x.view(np.int64) for x in to_concat],
                                        axis=axis)
            return new_values.view(_TD_DTYPE)

    # need to coerce to object
    to_concat = [convert_to_pydatetime(x, axis) for x in to_concat]

    return np.concatenate(to_concat,axis=axis)
