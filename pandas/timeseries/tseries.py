"""
The :class:`TimeSeries` class provides a base for the definition of time series.
A time series is defined here as the combination of two arrays:

- an array storing the time information
  (as a :class:`~scikits.timeseries.tdates.DateArray` instance);
- an array storing the data (as a :class:`MaskedArray` instance.)

These two classes were liberally adapted from :class:`MaskedArray` class.


:author: Pierre GF Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu - mattknox_ca_at_hotmail_dot_com
"""

#!!!: * Allow different lengths for data and dates to handle 2D data more easily
#!!!:    In that case, just make sure that the data is (n,rows,cols) where n is the nb of dates
#!!!: * Add some kind of marker telling whether we are 1D or nD:
#!!!:    That could be done by checking the ratio series.size/series._dates.size
#!!!: * Disable some of the tests on date compatibility if we are nD
#!!!: * Adapt reshaping to preserve the first dimension: that goes for squeeze

__author__ = "Pierre GF Gerard-Marchant & Matt Knox"
__revision__ = "$Revision$"
__date__ = '$Date$'

import sys
import warnings

import numpy as np
from numpy import bool_, complex_, float_, int_, object_, dtype, \
    ndarray, recarray
import numpy.core.umath as umath
from numpy.core.records import fromarrays as recfromarrays

from numpy import ma
from numpy.ma import MaskedArray, MAError, masked, nomask, \
    filled, getmask, getmaskarray, hsplit, make_mask_none, mask_or, make_mask, \
    masked_array

import tdates
from tdates import \
    DateError, FrequencyDateError, InsufficientDateError, Date, DateArray, \
    date_array, now, check_freq, check_freq_str, nodates

import const as _c
import pandas._skts

__all__ = ['TimeSeries', 'TimeSeriesCompatibilityError', 'TimeSeriesError',
           'adjust_endpoints', 'align_series', 'align_with', 'aligned',
           'asrecords',
           'compressed', 'concatenate', 'convert',
           'day', 'day_of_year',
           'empty_like',
           'fill_missing_dates', 'find_duplicated_dates', 'first_unmasked_val',
           'flatten',
           'hour',
           'last_unmasked_val',
           'minute', 'month',
           'pct', 'pct_log', 'pct_symmetric',
           'quarter',
           'remove_duplicated_dates',
           'second', 'split', 'stack',
           'time_series', 'tofile', 'tshift', 'masked', 'nomask',
           'week', 'weekday',
           'year',
           ]


def _unmasked_val(a, kind, axis=None):
    "helper function for first_unmasked_val and last_unmasked_val"
    if axis is None or a.ndim == 1:
        a = a.ravel()
        m = getmask(a)
        if m is nomask or not np.any(m):
            if kind == 0:
                indx = 0
            else:
                indx = -1
        else:
            indx = np.flatnonzero(~m)[[0, -1]][kind]
    else:
        m = ma.getmaskarray(a)
        indx = ma.array(np.indices(a.shape), mask=np.asarray([m] * a.ndim))
        if kind == 0:
            indx = tuple([indx[i].min(axis=axis).filled(0)
                          for i in range(a.ndim)])
        else:
            indx = tuple([indx[i].max(axis=axis).filled(0)
                          for i in range(a.ndim)])
    return a[indx]

def first_unmasked_val(a, axis=None):
    """
    Retrieve the first unmasked value along the given axis in a MaskedArray.

    Parameters
    ----------
    a : MaskedArray
        Input MaskedArray (or a subclass of).
    axis : int, optional
        Axis along which to perform the operation.
        If None, applies to a flattened version of the array.

    Returns
    -------
    val : {singleton of type marray.dtype}
        First unmasked value in a.
        If all values in a are masked, returns the numpy.ma.masked constant.
    """
    return _unmasked_val(a, 0, axis=axis)


def last_unmasked_val(a, axis=None):
    """
    Retrieve the last unmasked value along the given axis in a MaskedArray.

    Parameters
    ----------
    a : MaskedArray
        Input MaskedArray (or a subclass of).
    axis : int, optional
        Axis along which to perform the operation.
        If None, applies to a flattened version of the array.

    Returns
    -------
    val : {singleton of type marray.dtype}
        Last unmasked value in a.
        If all values in a are masked, returns the numpy.ma.masked constant.
    """
    return _unmasked_val(a, 1, axis=axis)

#### -------------------------------------------------------------------------
#--- ... TimeSeriesError class ...
#### -------------------------------------------------------------------------
class TimeSeriesError(Exception):
    "Class for TS related errors."
    def __init__ (self, value=None):
        "Creates an exception."
        self.value = value
    def __str__(self):
        "Calculates the string representation."
        return str(self.value)
    __repr__ = __str__


class TimeSeriesCompatibilityError(TimeSeriesError):
    """
    Defines the exception raised when series are incompatible.
    Incompatibility can arise from:

    * Inconsistent frequency;
    * Inconsistent starting dates;
    * Inconsistent size and/or shape.

    """
    def __init__(self, mode, first, second):
        if mode == 'freq':
            msg = "Incompatible time steps! (%s <> %s)"
        elif mode == 'start_date':
            msg = "Incompatible starting dates! (%s <> %s)"
        elif mode in ('size', 'shape'):
            msg = "Incompatible sizes! (%s <> %s)"
        elif mode == 'order':
            msg = "The series must be sorted in chronological order !"
        else:
            msg = "Incompatibility !  (%s <> %s)"
        msg = msg % (first, second)
        TimeSeriesError.__init__(self, msg)

#???: Should we go crazy and add some new exceptions ?
#???: TimeSeriesShapeCompatibilityError
#???: TimeSeriesStepCompatibilityError


def _timeseriescompat(a, b, raise_error=True):
    """
    Checks the date compatibility of two TimeSeries object.
    Returns True if everything's fine, or raises an exception.
    """
    #!!!: We need to use _varshape to simplify the analysis
    # Check the frequency ..............
    (afreq, bfreq) = (getattr(a, 'freq', None), getattr(b, 'freq', None))
    if afreq != bfreq:
        if raise_error:
            raise TimeSeriesCompatibilityError('freq', afreq, bfreq)
        return False
    # Make sure a.freq is not None
    if afreq is None:
        return True
    # Make sure that the series are sorted in chronological order
    if (not a.is_chronological()) or (not b.is_chronological()):
        if raise_error:
            raise TimeSeriesCompatiblityError('sort', None, None)
        return False
    # Check the starting dates ..........
    (astart, bstart) = (getattr(a, 'start_date'), getattr(b, 'start_date'))
    if astart != bstart:
        if raise_error:
            raise TimeSeriesCompatibilityError('start_date', astart, bstart)
        return False
    # Check the time steps ..............
    asteps = getattr(a, '_dates', a).get_steps()
    bsteps = getattr(b, '_dates', b).get_steps()
    step_diff = (asteps != bsteps)
    if (step_diff is True) or \
       (hasattr(step_diff, "any") and step_diff.any()):
        if raise_error:
            raise TimeSeriesCompatibilityError('time_steps', asteps, bsteps)
        return False
    elif a.shape != b.shape:
        if raise_error:
            raise TimeSeriesCompatibilityError('size', "1: %s" % str(a.shape),
                                                       "2: %s" % str(b.shape))
        return False
    return True

def _timeseriescompat_multiple(*series):
    """
    Checks the date compatibility of multiple TimeSeries objects.
    Returns True if everything's fine, or raises an exception. Unlike
    the binary version, all items must be TimeSeries objects.
    """

    defsteps = series[0]._dates.get_steps()

    def _check_steps(ser):
        _defsteps = s._dates.get_steps()
        _ds_comp = (_defsteps != defsteps)
        if not hasattr(_ds_comp, "any") or _ds_comp.any():
            return True
        else:
            return False

    (freqs, start_dates, steps, shapes) = \
                                zip(*[(s.freq,
                                       s.start_date,
                                       _check_steps(s),
                                       s.shape) for s in series])
    # Check the frequencies ................
    freqset = set(freqs)
    if len(set(freqs)) > 1:
        err_items = tuple(freqset)
        raise TimeSeriesCompatibilityError('freq', err_items[0], err_items[1])
    # Check the strting dates ..............
    startset = set(start_dates)
    if len(startset) > 1:
        err_items = tuple(startset)
        raise TimeSeriesCompatibilityError('start_dates',
                                           err_items[0], err_items[1])
    # Check the shapes .....................
    shapeset = set(shapes)
    if len(shapeset) > 1:
        err_items = tuple(shapeset)
        raise TimeSeriesCompatibilityError('size',
                                           "1: %s" % str(err_items[0]),
                                           "2: %s" % str(err_items[1]))
    # Check the steps ......................
    if max(steps) == True:
        bad_index = [x for (x, val) in enumerate(steps) if val][0]
        raise TimeSeriesCompatibilityError('time_steps',
                                           defsteps,
                                           series[bad_index]._dates.get_steps())
    return True


def get_varshape(data, dates):
    """
    Checks the compatibility of dates and data.

    Parameters
    ----------
    data : array-like
        Array of data
    dates : Date, DateArray
        Sequence of dates

    Returns
    -------
    varshape : tuple
        A tuple indicating the shape of the data at any date.

    Raises
    ------
    A :exc:`TimeSeriesCompatibilityError` exception is raised if something goes
    wrong.

    """

    dshape = data.shape
    dates = np.array(dates, copy=False, ndmin=1)
    tshape = dates.shape
    err_args = ('shape', "data: %s" % str(dshape), "dates: %s" % str(tshape))
    # Same size: all is well
    #???: The (not dates.size) is introduced to deal with masked
    if (not dates.size):
        return ()
    if (dates.size == data.size):
        if (dates.ndim > 1) or (data.ndim < 2):
            return ()
#    if (dates.size == data.size):
#        if (dates.ndim > 1) or (data.ndim == 1):
#            return ()
    # More dates than data: not good
    if (dates.size > data.size) or (data.ndim == 1):
        raise TimeSeriesCompatibilityError(*err_args)
    #....................
    dcumulshape = np.cumprod(dshape).tolist()
    try:
        k = dcumulshape.index(dates.size)
    except ValueError:
        raise TimeSeriesCompatibilityError(*err_args)
    else:
        return dshape[k + 1:]


def _getdatalength(data):
    "Estimates the length of a series (size/nb of variables)."
    if np.ndim(data) >= 2:
        return np.asarray(np.shape(data))[:-1].prod()
    else:
        return np.size(data)

def _compare_frequencies(*series):
    """Compares the frequencies of a sequence of series.

Returns the common frequency, or raises an exception if series have different
frequencies.
"""
    unique_freqs = np.unique([x.freqstr for x in series])
    try:
        common_freq = unique_freqs.item()
    except ValueError:
        raise TimeSeriesError, \
            "All series must have same frequency! (got %s instead)" % \
            unique_freqs
    return common_freq



##### ------------------------------------------------------------------------
##--- ... Time Series ...
##### ------------------------------------------------------------------------
_print_templates = dict(desc="""\
timeseries(
 %(data)s,
    dates =
 %(time)s,
    freq  = %(freq)s)
""",
                        desc_short="""\
timeseries(%(data)s,
   dates = %(time)s,
   freq  = %(freq)s)
""",
                        desc_flx="""\
timeseries(
 %(data)s,
   dtype = %(dtype)s,
   dates =
 %(time)s,
   freq  = %(freq)s)
""",
                        desc_flx_short="""\
timeseries(%(data)s,
   dtype = %(dtype)s,
   dates = %(time)s,
   freq  = %(freq)s)
"""
)



class _tsmathmethod(object):
    """
    Defines a wrapper for arithmetic array methods (add, mul...).
    When called, returns a new TimeSeries object, with the new series the result
    of the method applied on the original series. The `_dates` part remains
    unchanged.
    """
    def __init__ (self, methodname):
        self.__name__ = methodname
        self.__doc__ = getattr(MaskedArray, methodname).__doc__
        self.obj = None

    def __get__(self, obj, objtype=None):
        "Gets the calling object."
        self.obj = obj
        return self

    def __call__ (self, other, *args):
        "Execute the call behavior."
        instance = self.obj
        if isinstance(other, TimeSeries):
            compat = _timeseriescompat(instance, other, raise_error=False)
        else:
            compat = True
        func = getattr(super(TimeSeries, instance), self.__name__)
        if compat:
            result = np.array(func(other, *args), subok=True).view(type(instance))
            result._dates = instance._dates
        else:
            other_ = getattr(other, '_series', other)
            result_ = func(other_, *args)
            result = getattr(result_, '_series', result_)
        return result


class _tsarraymethod(object):
    """
    Defines a wrapper for basic array methods.
    When called, returns a new TimeSeries object, with the new series the result
    of the method applied on the original series.
    If `ondates` is True, the same operation is performed on the `_dates`.
    If `ondates` is False, the `_dates` part remains unchanged.
    """
    def __init__ (self, methodname, ondates=False):
        self.__name__ = methodname
        self.__doc__ = getattr(MaskedArray, methodname).__doc__
        self._ondates = ondates
        self.obj = None

    def __get__(self, obj, objtype=None):
        self.obj = obj
        return self

    def __call__ (self, *args, **kwargs):
        "Execute the call behavior."
        _name = self.__name__
        instance = self.obj
        # Fallback: if the instance has not been set, use the first argument
        if instance is None:
            args = list(args)
            instance = args.pop(0)
        _series = ndarray.__getattribute__(instance, '_series')
        _dates = ndarray.__getattribute__(instance, '_dates')
        func_series = getattr(_series, _name)
        result = func_series(*args, **kwargs).view(type(instance))
        if self._ondates:
            result._dates = getattr(_dates, _name)(*args, **kwargs)
        else:
            result._dates = _dates
        return result


class _tsaxismethod(object):
    """
    Defines a wrapper for array methods working on an axis (mean...).

    When called, returns a ndarray, as the result of the method applied on the
    series.
    """
    def __init__ (self, methodname):
        """abfunc(fillx, filly) must be defined.
           abinop(x, filly) = x for all x to enable reduce.
        """
        self.__name__ = methodname
        self.__doc__ = getattr(MaskedArray, methodname).__doc__
        self.obj = None

    def __get__(self, obj, objtype=None):
        self.obj = obj
        return self

    def __call__ (self, *args, **params):
        "Execute the call behavior."
        instance = self.obj
        if instance is None:
            args = list(args)
            instance = args.pop(0)
        (_dates, _series) = (instance._dates, instance._series)
        func = getattr(_series, self.__name__)
        result = func(*args, **params)
        if _dates.size != _series.size:
            axis = params.get('axis', None)
            if axis is None and len(args):
                axis = args[0]
            if axis in [-1, _series.ndim - 1]:
                result = result.view(type(instance))
                result._dates = _dates
        return result






class TimeSeries(MaskedArray, object):
    """
    Base class for the definition of time series.

    Parameters
    ----------
    data : {array_like}
        Data portion of the array.
        Any data that is valid for constructing a MaskedArray can be used here.
    dates : {DateArray}
        A `DateArray` instance.
    **optional_parameters:
        All the parameters recognized by `MaskedArray` are also recognized by
        TimeSeries.

    Notes
    -----
    It is recommended to use the :func:`time_series` function for construction,
    as it is more flexible and convenient.

    See Also
    --------
    numpy.ma.MaskedArray
        ndarray with support for missing data.
    scikits.timeseries.DateArray
    """

    __array_priority__ = 20

    def __new__(cls, data, dates, mask=nomask, dtype=None, copy=False,
                fill_value=None, subok=True, keep_mask=True, hard_mask=False,
                autosort=True, **options):

        maparms = dict(copy=copy, dtype=dtype, fill_value=fill_value,
                       subok=subok, keep_mask=keep_mask, hard_mask=hard_mask)
        _data = MaskedArray.__new__(cls, data, mask=mask, **maparms)

        # Get the data .......................................................
        if not subok or not isinstance(_data, TimeSeries):
            _data = _data.view(cls)
        if _data is masked:
            assert(np.size(dates) == 1)
            return _data.view(cls)
        # Check that the dates and data are compatible in shape.
        _data._varshape = get_varshape(_data, dates)
        # Set the dates
        _data._dates = dates
        if autosort:
            _data.sort_chronologically()
        return _data

    def __array_finalize__(self, obj):
        self._varshape = getattr(obj, '_varshape', ())
        MaskedArray.__array_finalize__(self, obj)

    def _update_from(self, obj):
        _dates = getattr(self, '_dates', nodates)
        newdates = getattr(obj, '_dates', nodates)
        # Only update the dates if we don't have any
        if not getattr(_dates, 'size', 0):
            self.__setdates__(newdates)
        MaskedArray._update_from(self, obj)


    def view(self, dtype=None, type=None):
        try:
            output = super(TimeSeries, self).view(dtype=dtype, type=type)
        except ValueError:
            output = super(TimeSeries, self).view(dtype)
        if isinstance(output, TimeSeries):
            if output.dtype.fields:
                output._varshape = ()
            else:
                fields = self.dtype.fields
                if fields:
                    output._varshape = (len(fields),)
        return output


    def _get_series(self):
        """
    Returns a view of the instance as a regular masked array.

    This attribute is read-only.
        """
        _mask = self._mask
        if _mask.ndim == 0 and _mask:
            return masked
        return self.view(MaskedArray)
    series = _series = property(fget=_get_series)

    @property
    def varshape(self):
        """
    Returns the shape of the underlying variables.
        """
        return self._varshape


    def _index_checker(self, indx):
        """
    Private function to process the index.
        """
        # Basic index ............
        if isinstance(indx, int):
            return (indx, indx, False)
        _dates = self._dates
        # String index : field name or date ?
        if isinstance(indx, basestring):
            if indx in (self.dtype.names or ()):
                return (indx, slice(None, None, None), False)
            try:
                indx = _dates.date_to_index(Date(_dates.freq, string=indx))
            except IndexError:
                # Trap the exception: we need the traceback
                exc_info = sys.exc_info()
                msg = "Invalid field or date '%s'" % indx
                raise IndexError(msg), None, exc_info[2]
            return (indx, indx, False)
        # Date or DateArray index ..........
        if isinstance(indx, (Date, DateArray)):
            indx = _dates.date_to_index(indx)
            return (indx, indx, False)
        # Slice ............................
        if isinstance(indx, slice):
            indx = slice(self._slicebound_checker(indx.start),
                         self._slicebound_checker(indx.stop),
                         indx.step)
            return (indx, indx, False)
        # Tuple index ......................
        if isinstance(indx, tuple):
            if not self._varshape:
                return (indx, indx, False)
            else:
                return (indx, indx[0], False)
        return (indx, indx, True)


    def _slicebound_checker(self, bound):
        "Private functions to check the bounds of a slice"
        # Integer bound (or None) ..........
        if bound is None or isinstance(bound, int):
            return bound
        # The bound is a date (string or Date)
        _dates = self._dates
        if isinstance(bound, str):
            bound = Date(_dates.freq, string=bound)
        if not isinstance(bound, Date):
            raise ValueError(
                "invalid object used in slice: %s" % repr(bound))
        if bound.freq != _dates.freq:
            raise TimeSeriesCompatibilityError('freq',
                                               _dates.freq, bound.freq)
        # this allows for slicing with dates outside the end points of the
        # series and slicing on series with missing dates
        return np.sum(self._dates < bound)


    def __getitem__(self, indx):
        """x.__getitem__(y) <==> x[y]

    Returns the item described by i. Not a copy.
        """
        (sindx, dindx, recheck) = self._index_checker(indx)
        _data = ndarray.__getattribute__(self, '_data')
        _mask = ndarray.__getattribute__(self, '_mask')
        _dates = ndarray.__getattribute__(self, '_dates')
        try:
            output = _data.__getitem__(sindx)
        except IndexError:
            # We don't need to recheck the index: just raise an exception
            if not recheck:
                raise
            # Maybe the index is a list of Dates ?
            try:
                indx = _dates.date_to_index(indx)
            except (IndexError, ValueError):
                # Mmh, is it a list of dates as strings ?
                try:
                    indx = _dates.date_to_index(date_array(indx,
                                                           freq=_dates.freq))
                except (IndexError, ValueError, DateError):
                    exc_info = sys.exc_info()
                    msg = "Invalid index or date '%s'" % indx
                    raise IndexError(msg), None, exc_info[2]
                else:
                    output = _data.__getitem__(indx)
                    sindx = dindx = indx
            else:
                output = _data.__getitem__(indx)
                sindx = dindx = indx
        # Don't find the date if it's not needed......
        if not getattr(output, 'ndim', False):
            # A record ................
            if isinstance(output, np.void):
                mask = _mask[sindx]
                if mask.view((bool, len(mask.dtype))).any():
                    output = masked_array(output, mask=mask)
                else:
                    return output
            elif _mask is not nomask and _mask[sindx]:
                return masked
            return output
        # Get the date................................
        newdates = _dates.__getitem__(dindx)
        if not getattr(newdates, 'shape', 0):
            # No dates ? Output a MaskedArray
            newseries = output.view(MaskedArray)
        else:
            # Some dates: output a TimeSeries
            newseries = output.view(type(self))
            newseries._dates = newdates
        # Update some info from self (fill_value, _basedict...)
        MaskedArray._update_from(newseries, self)
        # Fix the fill_value if we were accessing a field of a flexible array
        if isinstance(sindx, basestring):
            _fv = self._fill_value
            if _fv is not None and not np.isscalar(_fv):
                 newseries._fill_value = _fv[sindx]
            newseries._isfield = True
        # Fix the mask
        if _mask is not nomask:
            newseries._mask = _mask[sindx]
            newseries._sharedmask = True
        return newseries


    def __setitem__(self, indx, value):
        """x.__setitem__(i, y) <==> x[i]=y

    Sets item described by index. If value is masked, masks those locations.
        """
        (sindx, dindx, recheck) = self._index_checker(indx)
        _dates = ndarray.__getattribute__(self, '_dates')
        try:
            MaskedArray.__setitem__(self, sindx, value)
        except IndexError:
            # We don't need to recheck the index: just raise an exception
            if not recheck:
                raise
            # Maybe the index is a list of Dates ?
            try:
                indx = _dates.date_to_index(indx)
            except (IndexError, ValueError):
                # Mmh, is it a list of dates as strings ?
                try:
                    indx = _dates.date_to_index(date_array(indx,
                                                           freq=_dates.freq))
                except (IndexError, ValueError, DateError):
                    exc_info = sys.exc_info()
                    msg = "Invalid index or date '%s'" % indx
                    raise IndexError(msg), None, exc_info[2]
                else:
                    MaskedArray.__setitem__(self, indx, value)
            else:
                MaskedArray.__setitem__(self, indx, value)

    def __setattr__(self, attr, value):
        if attr in ['_dates', 'dates']:
            return self.__setdates__(value)
        elif attr == 'shape':
            if self._varshape:
                err_msg = "Reshaping a nV/nD series is not implemented yet !"
                raise NotImplementedError(err_msg)
            else:
                self._dates.shape = value
        return ndarray.__setattr__(self, attr, value)


    def __setdates__(self, value):
        """
    Sets the dates to `value`.
        """
        # Make sure it's a DateArray
        if not isinstance(value, DateArray):
            err_msg = "The input dates should be a valid "\
                      "DateArray object (got %s instead)" % type(value)
            raise TypeError(err_msg)
        # Skip if dates is nodates (or empty)\
        if value is nodates or not getattr(value, 'size', 0):
            return super(TimeSeries, self).__setattr__('_dates', value)
        # Make sure it has the proper size
        tsize = getattr(value, 'size', 1)
        # Check the _varshape
        varshape = self._varshape
        if not varshape:
            # We may be using the default: retry
            varshape = self._varshape = get_varshape(self, value)
        # Get the data length (independently of the nb of variables)
        dsize = self.size // int(np.prod(varshape))
        if tsize != dsize:
            raise TimeSeriesCompatibilityError("size",
                                               "data: %s" % dsize,
                                               "dates: %s" % tsize)
#        # Check whether the dates are already sorted
#        if not value.is_chronological():
#            _cached = value._cachedinfo
#            idx = _cached['chronidx']
#            _series = self._series
#            if not varshape:
#                if self.ndim > 1:
#                    flatseries = _series.flat
#                    flatseries[:] = flatseries[idx]
#                else:
#                    _series[:] = _series[idx]
#            else:
#                inishape = self.shape
#                _series.shape = tuple([-1,]+list(varshape))
#                _series[:] = _series[idx]
#                _series.shape = inishape
#            _cached['chronidx'] = np.array([], dtype=int)
#            _cached['ischrono'] = True
        #
        if not varshape and (value.shape != self.shape):
            # The data is 1D
            value = value.reshape(self.shape)
        super(TimeSeries, self).__setattr__('_dates', value)
        return

    dates = property(fget=lambda self:self._dates,
                     fset=__setdates__)


    #
    def sort_chronologically(self):
        """
    Sort the series by chronological order (in place).

    Notes
    -----
    This method sorts the series **in place**.
    To sort the a copy of the series, use the :func:`sort_chronologically`
    function.

    See Also
    --------
    sort_chronologically
        Equivalent function.
        """
        _dates = self._dates
        _series = self._series
        if not _dates.is_chronological():
            _cached = _dates._cachedinfo
            idx = _cached['chronidx']
            if not self._varshape:
                flatseries = _series.flat
                flatseries[:] = flatseries[idx]
            else:
                inishape = self.shape
                _series.shape = tuple([-1, ] + list(self._varshape))
                _series[:] = _series[idx]
                _series.shape = inishape
            # Sort the dates and reset the cache
            flatdates = _dates.ravel()
            flatdates[:] = flatdates[idx]
            _cached['chronidx'] = np.array([], dtype=int)
            _cached['ischrono'] = True



    def __str__(self):
        """Returns a string representation of self (w/o the dates...)"""
        return str(self._series)


    def __repr__(self):
        """
    Calculates the repr representation, using masked for fill if it is
    enabled. Otherwise fill with fill value.
    """
        _dates = self._dates
        if np.size(self._dates) > 2 and self.is_valid():
            timestr = "[%s ... %s]" % (str(_dates[0]), str(_dates[-1]))
        else:
            timestr = str(_dates)
        kwargs = {'data': str(self._series), 'time': timestr,
                  'freq': self.freqstr, 'dtype': self.dtype}
        names = kwargs['dtype'].names
        if self.ndim <= 1:
            if names:
                return _print_templates['desc_flx_short'] % kwargs
            return _print_templates['desc_short'] % kwargs
        if names:
            return _print_templates['desc_flx'] % kwargs
        return _print_templates['desc'] % kwargs

    #............................................
    __add__ = _tsmathmethod('__add__')
    __radd__ = _tsmathmethod('__add__')
    __sub__ = _tsmathmethod('__sub__')
    __rsub__ = _tsmathmethod('__rsub__')
    __pow__ = _tsmathmethod('__pow__')
    __mul__ = _tsmathmethod('__mul__')
    __rmul__ = _tsmathmethod('__mul__')
    __div__ = _tsmathmethod('__div__')
    __rdiv__ = _tsmathmethod('__rdiv__')
    __truediv__ = _tsmathmethod('__truediv__')
    __rtruediv__ = _tsmathmethod('__rtruediv__')
    __floordiv__ = _tsmathmethod('__floordiv__')
    __rfloordiv__ = _tsmathmethod('__rfloordiv__')
    __eq__ = _tsmathmethod('__eq__')
    __ne__ = _tsmathmethod('__ne__')
    __lt__ = _tsmathmethod('__lt__')
    __le__ = _tsmathmethod('__le__')
    __gt__ = _tsmathmethod('__gt__')
    __ge__ = _tsmathmethod('__ge__')

    copy = _tsarraymethod('copy', ondates=True)
    compress = _tsarraymethod('compress', ondates=True)
    cumsum = _tsarraymethod('cumsum', ondates=False)
    cumprod = _tsarraymethod('cumprod', ondates=False)
    anom = _tsarraymethod('anom', ondates=False)

    sum = _tsaxismethod('sum')
    prod = _tsaxismethod('prod')
    mean = _tsaxismethod('mean')
    var = _tsaxismethod('var')
    std = _tsaxismethod('std')
    all = _tsaxismethod('all')
    any = _tsaxismethod('any')


    def ravel(self):
        """
    Returns a ravelled view of the instance.

    If the instance corresponds to one variable (e.g., ``self.varshape == ()``),
    the result is the ravelled view of the input.
    Otherwise, the result is actually a TimeSeries where the `series` attribute
    is a reshaped view of the input and where the `dates` attribute is
    a ravelled view of the input `dates` attribute.

    Examples
    --------
    >>> start_date = ts.Date('M', '2001-01')
    >>> dates = ts.date_array(start_date=start_date, length=4)
    >>> series = ts.time_series([[1, 2], [3, 4]], dates=dates)
    >>> series
    timeseries(
     [[1 2]
     [3 4]],
        dates =
     [[Jan-2001 Feb-2001] ... [Mar-2001 Apr-2001]],
        freq  = M)
    >>> series.ravel()
    timeseries([1 2 3 4],
       dates = [Jan-2001 ... Apr-2001],
       freq  = M)
    >>> series = ts.time_series([[1, 2], [3, 4]], start_date=start_date)
    >>> series
    timeseries(
     [[1 2]
     [3 4]],
        dates =
     [Jan-2001 Feb-2001],
        freq  = M)
    >>> series.ravel()
    timeseries(
     [[1 2]
     [3 4]],
        dates =
     [Jan-2001 Feb-2001],
        freq  = M)
        """
        _varshape = self._varshape
        if _varshape:
            newshape = tuple([-1] + [np.prod(_varshape), ])
            result = MaskedArray.reshape(self, *newshape)
        else:
            result = MaskedArray.ravel(self)
        result._dates = self._dates.ravel()
        return result


    def reshape(self, *newshape, **kwargs):
        """
    Returns a time series containing the data of a, but with a new shape.

    The result is a view to the original array; if this is not possible,
    a ValueError is raised.

    Parameters
    ----------
    shape : shape tuple or int
       The new shape should be compatible with the original shape. If an
       integer, then the result will be a 1D array of that length.
    order : {'C', 'F'}, optional
        Determines whether the array data should be viewed as in C
        (row-major) order or FORTRAN (column-major) order.

    Returns
    -------
    reshaped_array : array
        A new view to the timeseries.

    Warnings
    --------
    The `._dates` part is reshaped, but the order is NOT ensured.

        """
        kwargs.update(order=kwargs.get('order', 'C'))
        if self._varshape:
            try:
                bkdtype = (self.dtype, self._varshape)
                _series = self._series.view([('', bkdtype)])
                result = _series.reshape(*newshape, **kwargs)
                result = result.view(dtype=bkdtype, type=type(self))
                result._dates = ndarray.reshape(self._dates, *newshape, **kwargs)
            except:
                err_msg = "Reshaping a nV/nD series is not implemented yet !"
                raise NotImplementedError(err_msg)
        # 1D series : reshape the dates as well
        else:
            result = MaskedArray.reshape(self, *newshape, **kwargs)
            result._dates = ndarray.reshape(self._dates, *newshape, **kwargs)
            result._varshape = ()
        return result

    #.........................................................................
    def ids (self):
        """Return the ids of the data, dates and mask areas"""
        return (id(self._series), id(self.dates),)

    #.........................................................................
    @property
    def freq(self):
        """Returns the corresponding frequency (as an integer)."""
        return self._dates.freq
    @property
    def freqstr(self):
        """Returns the corresponding frequency (as a string)."""
        return self._dates.freqstr

    @property
    def year(self):
        """Returns the year for each date of the instance."""
        return self._dates.year
    years = year
    @property
    def qyear(self):
        """
    For quarterly frequencies, returns the equivalent of the 'fiscal' year
    for each date of the instance.
    For non-quarterly frequencies, returns the year.
        """
        return self._dates.qyear
    @property
    def quarter(self):
        """Returns the quarter for each date of the instance."""
        return self._dates.quarter
    quarters = quarter
    @property
    def month(self):
        """Returns the month for each date of the instance."""
        return self._dates.month
    months = month
    @property
    def week(self):
        """Returns the week for each date in self._dates."""
        return self._dates.week
    weeks = week
    @property
    def day(self):
        """Returns the day of month for each date of the instance."""
        return self._dates.day
    days = day
    @property
    def day_of_week(self):
        """Returns the day of week for each date of the instance."""
        return self._dates.weekday
    weekdays = weekday = day_of_week
    @property
    def day_of_year(self):
        """Returns the day of year for each date of the instance."""
        return self._dates.day_of_year
    yeardays = day_of_year
    @property
    def hour(self):
        """Returns the hour for each date in self._dates."""
        return self._dates.hour
    hours = hour
    @property
    def minute(self):
        """Returns the minute for each date in self._dates."""
        return self._dates.minute
    minutes = minute
    @property
    def second(self):
        """Returns the second for each date in self._dates."""
        return self._dates.second
    seconds = second

    @property
    def start_date(self):
        """Returns the first date of the series."""
        _dates = self._dates
        dsize = _dates.size
        if dsize == 0:
            return None
        elif dsize == 1:
            return _dates[0]
        else:
            return Date(self.freq, _dates.flat[0])

    @property
    def end_date(self):
        """Returns the last date of the series."""
        _dates = self._dates
        dsize = _dates.size
        if dsize == 0:
            return None
        elif dsize == 1:
            return _dates[-1]
        else:
            return Date(self.freq, _dates.flat[-1])

    def is_valid(self):
        """Returns whether the series has no duplicate/missing dates."""
        return self._dates.is_valid()

    def isvalid(self):
        """Deprecated name: use '.is_valid' instead."""
        return self._dates.isvalid()

    def has_missing_dates(self):
        """Returns whether there's a date gap in the series."""
        return self._dates.has_missing_dates()

    def is_full(self):
        """Returns whether there's no date gap in the series."""
        return self._dates.is_full()

    def isfull(self):
        """Deprecated name: use '.is_full' instead."""
        return self._dates.isfull()

    def has_duplicated_dates(self):
        """Returns whether there are duplicated dates in the series."""
        return self._dates.has_duplicated_dates()

    def get_steps(self):
        """
    Returns the time steps between consecutive dates, in the same unit as
    the frequency of the instance.
        """
        return self._dates.get_steps()

    def is_chronological(self):
        """Returns whether the series is in chronological order."""
        return self._dates.is_chronological()

    def date_to_index(self, date):
        """Returns the index corresponding to a given date, as an integer."""
        return self._dates.date_to_index(date)

    #.....................................................
    def asfreq(self, freq, relation="END"):
        """
    Converts the dates portion of the TimeSeries to another frequency.

    The resulting TimeSeries will have the same shape and dimensions
    as the original series (unlike the :meth:`convert` method).

    Parameters
    ----------
    freq : {freq_spec}
    relation : {'END', 'START'} (optional)

    Returns
    -------
    A new TimeSeries with the :attr:`.dates` :class:`DateArray` at the
    specified frequency (the :meth`.asfreq` method of the :attr:`.dates`
    property will be called).
    The data in the resulting series will be a VIEW of the original series.

    Notes
    -----
    The parameters are the exact same as for
    :meth:`~scikit.timeseries.DateArray.asfreq`. Please see the docstring for
    that method for details on the parameters and how the actual conversion is
    performed.

    """
        if freq is None: return self

        return TimeSeries(self._series,
                          dates=self._dates.asfreq(freq, relation=relation))
    #.....................................................
    def transpose(self, *axes):
        if self._dates.size == self.size:
            result = MaskedArray.transpose(self, *axes)
            result._dates = self._dates.transpose(*axes)
        else:
            errmsg = "Operation not permitted on multi-variable series"
            if (len(axes) == 0) or axes[0] != 0:
                raise TimeSeriesError, errmsg
            else:
                result = MaskedArray.transpose(self, *axes)
                result._dates = self._dates
        return result
    transpose.__doc__ = np.transpose.__doc__

    def split(self):
        """Split a multi-dimensional series into individual columns."""
        if self.ndim == 1:
            return [self]
        else:
            n = self.shape[1]
            arr = hsplit(self, n)[0]
            return [self.__class__(np.squeeze(a),
                                   self._dates,
                                   **_attrib_dict(self)) for a in arr]

    def filled(self, fill_value=None):
        """
    Returns an array of the same class as `_data`,  with masked values
    filled with `fill_value`. Subclassing is preserved.

    Parameters
    ----------
    fill_value : {None, singleton of type self.dtype}, optional
        The value to fill in masked values with.
        If `fill_value` is None, uses ``self.fill_value``.

    """
        result = self._series.filled(fill_value=fill_value).view(type(self))
        result._dates = self._dates
        return result


    def tolist(self):
        """
    Returns the dates and data portion of the TimeSeries "zipped" up in
    a list of standard python objects (eg. datetime, int, etc...).

        """
        if self.ndim > 0:
            return zip(self.dates.tolist(), self.series.tolist())
        else:
            return self.series.tolist()

    def torecords(self):
        """
    Transforms a TimeSeries into a structured array with three fields:.

    * the ``_dates`` field stores the date information;
    * the ``_data`` field stores the ``_data`` part of the series;
    * the ``_mask`` field stores the mask.


    Returns
    -------
    record : ndarray
        A new flexible-type ndarray with three fields: the first element
        contains the date (as an integer), the second element contains the
        corresponding value and the third the corresponding mask boolean.
        The returned record shape matches the shape of the instance.

        """
        _varshape = self._varshape
        if not _varshape:
            desctype = [('_dates', int),
                        ('_data', self.dtype),
                        ('_mask', self.mask.dtype)]
        else:
            desctype = [('_dates', int),
                        ('_data', (self.dtype, _varshape)),
                        ('_mask', (self.mask.dtype, _varshape))]
        flat = self.ravel()
        _dates = np.asarray(flat._dates)
        if flat.size > 0:
            result = np.empty(len(_dates), dtype=desctype)
            result['_dates'] = _dates
            result['_data'] = flat._data
            result['_mask'] = flat._mask
            return result
        else:
            return np.array([[], [], []], dtype=desctype,)

    # for backwards compatibility
    toflex = torecords
    #......................................................
    # Pickling
    def __getstate__(self):
        """

    Returns the internal state of the TimeSeries, for pickling purposes.
        """
    #    raise NotImplementedError,"Please use timeseries.archive/unarchive instead."""
        state = (1,
                 self.shape,
                 self.dtype,
                 self.flags.fnc,
                 self._data.tostring(),
                 getmaskarray(self).tostring(),
                 self._fill_value,
                 self._dates.shape,
                 self._dates.__array__().tostring(),
                 self.freq,
                 self._optinfo,
                 )
        return state
    #
    def __setstate__(self, state):
        """

    Restores the internal state of the TimeSeries, for pickling purposes.
    `state` is typically the output of the ``__getstate__`` output, and is a 5-tuple:

        - class name
        - a tuple giving the shape of the data
        - a typecode for the data
        - a binary string for the data
        - a binary string for the mask.
        """
        (ver, shp, typ, isf, raw, msk, flv, dsh, dtm, frq, infodict) = state
        MaskedArray.__setstate__(self, (ver, shp, typ, isf, raw, msk, flv))
        _dates = self._dates
        _dates.__setstate__((ver, dsh, dtype(int_), isf, dtm, frq))
        _dates.freq = frq
        _dates._cachedinfo.update(dict(full=None, hasdups=None, steps=None,
                                       toobj=None, toord=None, tostr=None))
        # Update the _optinfo dictionary
        self._optinfo.update(infodict)
#
    def __reduce__(self):
        """Returns a 3-tuple for pickling a MaskedArray."""
        return (_tsreconstruct,
                (self.__class__, self._baseclass,
                 self.shape, self._dates.shape, self.dtype, self._fill_value),
                self.__getstate__())

def _tsreconstruct(genclass, baseclass, baseshape, dateshape, basetype, fill_value):
    """Internal function that builds a new TimeSeries from the information stored
    in a pickle."""
    #    raise NotImplementedError,"Please use timeseries.archive/unarchive instead."""
    _series = ndarray.__new__(baseclass, baseshape, basetype)
    _dates = ndarray.__new__(DateArray, dateshape, int_)
    _mask = ndarray.__new__(ndarray, baseshape, bool_)
    return genclass.__new__(genclass, _series, dates=_dates, mask=_mask,
                            dtype=basetype, fill_value=fill_value)


def _attrib_dict(series, exclude=[]):
    """this function is used for passing through attributes of one
time series to a new one being created"""
    result = {'fill_value':series.fill_value}
    return dict(filter(lambda x: x[0] not in exclude, result.iteritems()))


##### --------------------------------------------------------------------------
##--- ... Additional methods ...
##### --------------------------------------------------------------------------

def _extrema(self, method, axis=None, fill_value=None):
    "Private function used by max/min"
    (_series, _dates) = (self._series, self._dates)
    func = getattr(_series, method)
    idx = func(axis, fill_value)

    # 1D series .......................
    if (_dates.size == _series.size):
        if axis is None:
            return self.ravel()[idx]
        else:
            return self[idx]
    # nD series .......................
    else:
        if axis is None:
            idces = np.unravel_index(idx, _series.shape)
            result = time_series(_series[idces], dates=_dates[idces[0]])
        else:
            _shape = _series.shape
            _dates = np.repeat(_dates, np.prod(_shape[1:])).reshape(_shape)
            _s = ma.choose(idx, np.rollaxis(_series, axis, 0))
            _d = ma.choose(idx, np.rollaxis(_dates, axis, 0))
            result = time_series(_s, dates=_d, freq=_dates.freq)
        return result

def _max(self, axis=None, fill_value=None):
    """Return the maximum of self along the given axis.
    Masked values are filled with fill_value.

    Parameters
    ----------
    axis : int, optional
        Axis along which to perform the operation.
        If None, applies to a flattened view of the array.
    fill_value : {var}, optional
        Value used to fill in the masked values.
        If None, use the the output of maximum_fill_value().
    """
    return _extrema(self, 'argmax', axis, fill_value)
TimeSeries.max = _max

def _min(self, axis=None, fill_value=None):
    """Return the minimum of self along the given axis.
    Masked values are filled with fill_value.

    Parameters
    ----------
    axis : int, optional
        Axis along which to perform the operation.
        If None, applies to a flattened view of the array.
    fill_value : {var}, optional
        Value used to fill in the masked values.
        If None, use the the output of minimum_fill_value().
    """
    return _extrema(self, 'argmin', axis, fill_value)
TimeSeries.min = _min


#.......................................


class _tsblockedmethods(object):
    """Defines a wrapper for array methods that should be temporarily disabled.
    """
    def __init__ (self, methodname):
        """abfunc(fillx, filly) must be defined.
           abinop(x, filly) = x for all x to enable reduce.
        """
        self._name = methodname
        self.obj = None
    #
    def __get__(self, obj, objtype=None):
        self.obj = obj
        return self
    #
    def __call__ (self, *args, **params):
        raise NotImplementedError

TimeSeries.swapaxes = _tsarraymethod('swapaxes', ondates=True)

#####---------------------------------------------------------------------------
#---- --- Definition of functions from the corresponding methods ---
#####---------------------------------------------------------------------------
class _frommethod(object):
    """Defines functions from existing MaskedArray methods.
:ivar _methodname (String): Name of the method to transform.
    """
    def __init__(self, methodname):
        self.__name__ = methodname
        self.__doc__ = self.getdoc()
    def getdoc(self):
        "Returns the doc of the function (from the doc of the method)."
        try:
            return getattr(TimeSeries, self.__name__).__doc__
        except:
            return "???"
    #
    def __call__ (self, caller, *args, **params):
        if hasattr(caller, self.__name__):
            method = getattr(caller, self.__name__)
            # If method is not callable, it's a property, and don't call it
            if hasattr(method, '__call__'):
                return method.__call__(*args, **params)
            return method
        method = getattr(np.asarray(caller), self.__name__)
        try:
            return method(*args, **params)
        except SystemError:
            return getattr(np, self.__name__).__call__(caller, *args, **params)
#............................
weekday = _frommethod('weekday')
day_of_year = _frommethod('day_of_year')
week = _frommethod('week')
year = _frommethod('year')
quarter = _frommethod('quarter')
month = _frommethod('month')
day = _frommethod('day')
hour = _frommethod('hour')
minute = _frommethod('minute')
second = _frommethod('second')

split = _frommethod('split')

torecords = toflex = _frommethod('toflex')

#
##### ---------------------------------------------------------------------------
#---- ... Additional methods ...
##### ---------------------------------------------------------------------------
def tofile(series, fileobject, format=None,
           separator=" ", linesep='\n', precision=5,
           suppress_small=False, keep_open=False):
    """
    Writes the TimeSeries to a file. The series should be 2D at most.

    Parameters
    ----------
    series : TimeSeries
        The array to write.
    fileobject
        An open file object or a string to a valid filename.
    format : {None, string}, optional
        Format string for the date.
        If None, uses the default date format.
    separator : {string}, optional
        Separator to write between elements of the array.
    linesep : {string}, optional
        Separator to write between rows of array.
    precision : {integer}, optional
        Number of digits after the decimal place to write.
    suppress_small : {boolean}, optional
        Whether on-zero to round small numbers down to 0.0
    keep_open : {boolean}, optional
        Whether to close the file or to return the open file.

    Returns
    -------
    file : file object
        The open file (if ``keep_open`` is non-zero).
    """

    try:
        import scipy.io
    except ImportError:
        raise ImportError("scipy is required for the tofile function/method")

    (_dates, _data) = (series._dates, series._series)
    optpars = dict(separator=separator, linesep=linesep, precision=precision,
                   suppress_small=suppress_small, keep_open=keep_open)
    if _dates.size == _data.size:
        # 1D version
        tmpfiller = ma.empty((_dates.size, 2), dtype=np.object_)
        _data = _data.reshape(-1)
        tmpfiller[:, 1:] = ma.atleast_2d(_data).T
    else:
        sshape = list(_data.shape)
        sshape[-1] += 1
        tmpfiller = ma.empty(sshape, dtype=np.object_)
        tmpfiller[:, 1:] = _data
    #
    if format is None:
        tmpfiller[:, 0] = _dates.ravel().tostring()
    else:
        tmpfiller[:, 0] = [_.strftime(format) for _ in _dates.ravel()]
    return scipy.io.write_array(fileobject, tmpfiller, **optpars)


TimeSeries.tofile = tofile


def flatten(series):
    """
    Flattens a (multi-) time series to 1D series.

    """
    shp_ini = series.shape
    # Already flat time series....
    if len(shp_ini) == 1:
        return series
    # Folded single time series ..
    newdates = series._dates.ravel()
    if series._dates.size == series._series.size:
        newshape = (series._series.size,)
    else:
        newshape = (np.asarray(shp_ini[:-1]).prod(), shp_ini[-1])
    newseries = series._series.reshape(newshape)
    return time_series(newseries, newdates)

TimeSeries.flatten = flatten


def compressed(series):
    """
    Suppresses missing values from a time series.

    Returns a :class:`TimeSeries` object.

    """
    if series._mask is nomask:
        return series
    if series.ndim == 1:
        keeper = ~(series._mask)
    elif series.ndim == 2:
        _dates = series._dates
        _series = series._series
        # Both dates and data are 2D: ravel first
        if _dates.ndim == 2:
            series = series.ravel()
            keeper = ~(series._mask)
        # 2D series w/ only one date : return a new series ....
        elif _dates.size == 1:
            result = _series.compressed().view(type(series))
            result._dates = series.dates
            return result
        # a 2D series: suppress the rows (dates are in columns)
        else:
            keeper = ~(series._mask.any(-1))
    else:
        raise NotImplementedError
    return series[keeper]
TimeSeries.compressed = compressed


##### -------------------------------------------------------------------------
#---- --- TimeSeries constructor ---
##### -------------------------------------------------------------------------

def time_series(data, dates=None, start_date=None, length=None, freq=None,
                mask=nomask, dtype=None, copy=False, fill_value=None,
                keep_mask=True, hard_mask=False, autosort=True):
    """
    Creates a TimeSeries object.

    The ``data`` parameter can be a valid :class:`TimeSeries` object.
    In that case, the ``dates``, ``start_date`` or ``freq`` parameters are
    optional: if none of them is given, the dates of the result are the dates of
    ``data``.

    If ``data`` is not a :class:`TimeSeries`, then ``dates`` must be either
    ``None`` or an object recognized by the :func:`date_array` function (used
    internally):

        * an existing :class:`DateArray` object;
        * a sequence of :class:`Date` objects with the same frequency;
        * a sequence of :class:`datetime.datetime` objects;
        * a sequence of dates in string format;
        * a sequence of integers corresponding to the representation of
          :class:`Date` objects.

    In any of the last four possibilities, the ``freq`` parameter is mandatory.

    If ``dates`` is ``None``, a continuous :class:`DateArray` is automatically
    constructed as an array of size ``len(data)`` starting at ``start_date`` and
    with a frequency ``freq``.


    Parameters
    ----------
    data : array_like
        Data portion of the array. Any data that is valid for constructing a
        :class:`~numpy.ma.MaskedArray` can be used here.
        :keyword:`data` can also be a :class:`TimeSeries` object.
    dates : {None, var}, optional
        A sequence of dates corresponding to each entry.
    start_date : {Date}, optional
        Date corresponding to the first entry of the data (index 0).
        This parameter must be a valid :class:`Date` object, and is mandatory
        if ``dates`` is None and if ``data`` has a length greater or equal to 1.
    length : {integer}, optional
        Length of the dates.
    freq : {freq_spec}, optional
        A valid frequency specification, as a string or an integer.
        This parameter is mandatory if ``dates`` is None.
        Otherwise, the frequency of the series is set to the frequency
        of the ``dates`` input.

    Notes
    -----
    * All other parameters recognized by the :func:`numpy.ma.array` constructor
      are also recognized by the function.
    * If ``data`` is zero-sized, only the ``freq`` parameter is mandatory.

    See Also
    --------
    numpy.ma.masked_array
        Constructor for the :class:`~numpy.ma.MaskedArray` class.
    scikits.timeseries.date_array
        Constructor for the :class:`DateArray` class.

    """
    freq = check_freq(freq)

    if dates is None:
        _dates = getattr(data, '_dates', None)
    elif isinstance(dates, (Date, DateArray)):
        if copy:
            _dates = date_array(dates, autosort=False).copy()
        else:
            _dates = date_array(dates, autosort=False)
    elif isinstance(dates, (tuple, list, ndarray)):
        _dates = date_array(dlist=dates, freq=freq, autosort=False)
    else:
        _dates = date_array([], freq=freq)

    if _dates is not None:
        # Make sure _dates has the proper frequency
        if (freq != _c.FR_UND) and (_dates.freq != freq):
            _dates = _dates.asfreq(freq)
    else:
        dshape = np.shape(data)
        if len(dshape) > 0:
            length = length or dshape[0]
            _dates = date_array(start_date=start_date, freq=freq, length=length)
        else:
            _dates = date_array([], freq=freq)

    return TimeSeries(data=data, mask=mask, dates=_dates,
                      copy=copy, dtype=dtype, subok=True,
                      fill_value=fill_value, keep_mask=keep_mask,
                      hard_mask=hard_mask, autosort=autosort)




##### --------------------------------------------------------------------------
#---- ... Additional functions ...
##### --------------------------------------------------------------------------


def sort_chronologically(series):
    """
    Returns a copy of series, sorted chronologically.
    """
    series = series.copy()
    series.sort_chronologically()
    return series


def asrecords(series):
    "Deprecated version of torecords"
    warnings.warn("Deprecated function: use torecords instead")
    return torecords(series)


def adjust_endpoints(a, start_date=None, end_date=None, copy=False):
    """
    Returns a TimeSeries going from `start_date` to `end_date`.

    Parameters
    ----------
    a : TimeSeries
        TimeSeries object whose dates must be adjusted
    start_date : Date, optional
        New starting date. If not specified, the current starting date is used.
    end_date : Date, optional
        New ending date. If not specified, the current ending date is used.
    copy : {False, True}, optional
        Whether to return a copy of the initial array (:const:`True`)
        or a reference to the array (:const:`False`), in the case where both
        the `start_date` and `end_date` both fall into the initial range
        of dates.

    """
    # Series validity tests .....................
    if not isinstance(a, TimeSeries):
        raise TypeError, "Argument should be a valid TimeSeries object!"
    if a.freq == 'U':
        errmsg = "Cannot adjust a series with 'Undefined' frequency."
        raise TimeSeriesError(errmsg,)

    if not a._dates.is_valid():
        errmsg = "Cannot adjust a series with missing or duplicated dates."
        raise TimeSeriesError(errmsg,)
    # Flatten the series if needed ..............
    a = a.flatten()
    shp_flat = a.shape
    # Dates validity checks .,...................
    msg = "%s should be a valid Date object! (got %s instead)"
    if a.dates.size >= 1:
        (dstart, dend) = a.dates[[0, -1]]
    else:
        (dstart, dend) = (None, None)
    # Skip the empty series case
    if dstart is None and (start_date is None or end_date is None):
        errmsg = "Both start_date and end_date must be specified"\
                 " to adjust endpoints of a zero length series!"
        raise TimeSeriesError(errmsg,)
    #....
    if start_date is None:
        start_date = dstart
        start_lag = 0
    else:
        if isinstance(start_date, basestring):
            start_date = Date(a.freq, string=start_date)
        elif not isinstance(start_date, Date):
            raise TypeError, msg % ('start_date', type(start_date))
        if dstart is not None:
            start_lag = start_date - dstart
        else:
            start_lag = start_date
    #....
    if end_date is None:
        end_date = dend
        end_lag = 0
    else:
        if isinstance(end_date, basestring):
            end_date = Date(a.freq, string=end_date)
        elif not isinstance(end_date, Date):
            raise TypeError, msg % ('end_date', type(end_date))
        if dend is not None:
            end_lag = end_date - dend
        else:
            end_lag = end_date
    # Check if the new range is included in the old one
    if start_lag >= 0:
        if end_lag == 0:
            if not copy:
                return a[start_lag:]
            else:
                return a[start_lag:].copy()
        elif end_lag < 0:
            if not copy:
                return a[start_lag:end_lag]
            else:
                return a[start_lag:end_lag].copy()
    # Create a new series .......................
    newdates = date_array(start_date=start_date, end_date=end_date)

    newshape = list(shp_flat)
    newshape[0] = len(newdates)
    newshape = tuple(newshape)

    newseries = np.empty(newshape, dtype=a.dtype).view(type(a))
    #!!!: Here, we may wanna use something else than MaskType
    newseries.__setmask__(np.ones(newseries.shape, dtype=bool_))
    newseries._dates = newdates
    newseries._update_from(a)
    if dstart is not None:
        start_date = max(start_date, dstart)
        end_date = min(end_date, dend) + 1
        if not copy:
            newseries[start_date:end_date] = a[start_date:end_date]
        else:
            newseries[start_date:end_date] = a[start_date:end_date].copy()
    return newseries
TimeSeries.adjust_endpoints = adjust_endpoints



def align_series(*series, **kwargs):
    """
    Aligns several TimeSeries, so that their starting and ending dates match.

    Series are resized and filled with masked values accordingly.

    The resulting series have no missing dates (ie. ``series.is_valid() == True``
    for each of the resulting series).

    The function accepts two extras parameters:
    - `start_date` forces the series to start at that given date,
    - `end_date` forces the series to end at that given date.

    By default, `start_date` and `end_date` are set respectively to the smallest
    and largest dates of the series.
    """
    if len(series) < 2:
        return series
    unique_freqs = np.unique([x.freqstr for x in series])
    common_freq = _compare_frequencies(*series)

    # if any of the series have missing dates, fill them in first
    filled_series = []
    for ser in series:
        if ser.is_valid():
            filled_series.append(ser)
        else:
            filled_series.append(ser.fill_missing_dates())

    start_date = kwargs.pop('start_date',
                            min([x.start_date for x in filled_series
                                     if x.start_date is not None]))
    if isinstance(start_date, str):
        start_date = Date(common_freq, string=start_date)
    end_date = kwargs.pop('end_date',
                          max([x.end_date for x in filled_series
                                   if x.end_date is not None]))
    if isinstance(end_date, str):
        end_date = Date(common_freq, string=end_date)

    return [adjust_endpoints(x, start_date, end_date) for x in filled_series]
aligned = align_series



def align_with(*series):
    """
    Aligns several TimeSeries to the first of the list, so that their
    starting and ending dates match.

    The series are resized and padded with masked values accordingly.
    """
    if len(series) < 2:
        return series
    dates = series[0]._dates[[0, -1]]
    if len(series) == 2:
        return adjust_endpoints(series[-1], dates[0], dates[-1])
    return [adjust_endpoints(x, dates[0], dates[-1]) for x in series[1:]]


#....................................................................
def _convert1d(series, freq, func, position, *args, **kwargs):
    "helper function for `convert` function"
    # Check the frequencies ..........................
    to_freq = check_freq(freq)
    from_freq = series.freq
    # Don't do anything if not needed
    if from_freq == to_freq:
        return series
    if from_freq == _c.FR_UND:
        err_msg = "Cannot convert a series with UNDEFINED frequency."
        raise TimeSeriesError(err_msg)
    if to_freq == _c.FR_UND:
        err_msg = "Cannot convert a series to UNDEFINED frequency."
        raise TimeSeriesError(err_msg)
    # Check the validity of the series .....
    if not series.is_valid():
        err_msg = "Cannot adjust a series with missing or duplicated dates."
        raise TimeSeriesError(err_msg)

    # Check the position parameter..........
    position = position.upper()
    if position not in ('END', 'START'):
        err_msg = "Invalid value for position argument: (%s). "\
                  "Should be in ['END','START']," % str(position)
        raise ValueError(err_msg)

    start_date = series._dates[0]

    if series.size == 0:
        return TimeSeries(series, freq=to_freq,
                          start_date=start_date.asfreq(to_freq))

    data_ = series._series.filled()
    mask_ = getmaskarray(series)

    if (data_.size // series._dates.size) > 1:
        raise TimeSeriesError("convert works with 1D data only !")

    cdictresult = pandas._skts.TS_convert(data_, from_freq, to_freq, position,
                                     int(start_date), mask_)
    start_date = Date(freq=to_freq, value=cdictresult['startindex'])
    data_ = masked_array(cdictresult['values'], mask=cdictresult['mask'])

    if data_.ndim == 2:
        if func is None:
            newvarshape = data_.shape[1:]
        else:
            # Try to use an axis argument
            try:
                data_ = func(data_, axis= -1, *args, **kwargs)
            # Fall back to apply_along_axis (slower)
            except TypeError:
                data_ = ma.apply_along_axis(func, -1, data_, *args, **kwargs)
            newvarshape = ()
    elif data_.ndim == 1:
        newvarshape = ()

    newdates = DateArray(np.arange(len(data_)) + start_date, freq=to_freq)

    newseries = data_.view(type(series))
    newseries._varshape = newvarshape
    newseries._dates = newdates
    newseries._update_from(series)
    return newseries



def convert(series, freq, func=None, position='END', *args, **kwargs):
    """
    Converts a series from one frequency to another, by manipulating both the
    `data` and `dates` attributes.

    If the input series has any missing dates, it will first be filled in with
    masked values prior to doing the conversion.

    Parameters
    ----------
    series : TimeSeries
        Series to convert. Skip this parameter if you are calling this as
        a method of the TimeSeries object instead of the module function.
    freq : freq_spec
        Frequency to convert the TimeSeries to. Accepts any valid frequency
        specification (string or integer)
    func : function, optional
        When converting a series to a lower frequency, the :keyword:`func`
        parameter to perform a calculation on each period of values
        to aggregate results.
        For example, when converting a daily series to a monthly series, use
        :func:`numpy.ma.mean` to get a series of monthly averages.
        If the first or last value from a period, the functions
        :func:`~scikits.timeseries.first_unmasked_val` and
        :func:`~scikits.timeseries.last_unmasked_val` should be used instead.
        If :keyword:`func` is not given, the output series group the points
        of the initial series that share the same new date. Thus, if the
        initial series has a daily frequency and is 1D, the output series is
        2D.
    position : {'END', 'START'}, optional
        When converting a series to a higher frequency, use this parameter to
        determine where the points should fall in the new period.
        For example, when converting a monthly series to daily, using
        position='START' will cause the values to fall on the first day of
        each month (with all other values being masked).
    *args : {extra arguments for func parameter}, optional
        Mandatory parameters of the :keyword:`func` function.
    **kwargs : {extra keyword arguments for func parameter}, optional
        Optional keyword parameters of the :keyword:`func` function.

    Returns
    -------
    converted_series
        A new :class:`TimeSeries` at the given frequency, without any missing
        nor duplicated dates

    """
    #!!!: Raise some kind of proper exception if the underlying dtype will mess things up
    #!!!: For example, mean on string array...

    if series.ndim > 2 or series.ndim == 0:
        raise ValueError(
            "only series with ndim == 1 or ndim == 2 may be converted")

    if series.has_duplicated_dates():
        raise TimeSeriesError("The input series must not have duplicated dates!")

    if series.has_missing_dates():
        # can only convert continuous time series, so fill in missing dates
        series = fill_missing_dates(series)

    if series.ndim == 1:
        obj = _convert1d(series, freq, func, position, *args, **kwargs)
    elif series.ndim == 2:
        base = _convert1d(series[:, 0], freq, func, position, *args, **kwargs)
        obj = ma.column_stack([_convert1d(m, freq, func, position,
                                          *args, **kwargs)._series
                               for m in series.split()]).view(type(series))
        obj._dates = base._dates
        if func is None:
            shp = obj.shape
            ncols = base.shape[-1]
            obj.shape = (shp[0], shp[-1] // ncols, ncols)
            obj = np.swapaxes(obj, 1, 2)

    return obj
TimeSeries.convert = convert



def tshift(series, nper, copy=True):
    """
    Returns a series of the same size as `series`, with the same `start_date`
    and `end_date`, but values shifted by `nper`.

    Parameters
    ----------
    series : TimeSeries
        TimeSeries object to shift. Ignore this parameter if calling this as a
        method.
    nper : int
        Number of periods to shift. Negative numbers shift values to the right,
        positive to the left.
    copy : {True, False}, optional
        copies the data if True, returns a view if False.

    Examples
    --------
    >>> series = time_series([0,1,2,3], start_date=Date(freq='A', year=2005))
    >>> series
    timeseries(data  = [0 1 2 3],
               dates = [2005 ... 2008],
               freq  = A-DEC)
    >>> tshift(series, -1)
    timeseries(data  = [-- 0 1 2],
               dates = [2005 ... 2008],
               freq  = A-DEC)
    >>> pct_change = 100 * (series/series.tshift(-1, copy=False) - 1)

    """
    newdata = masked_array(np.empty(series.shape, dtype=series.dtype),
                           mask=True)
    if copy:
        inidata = series._series.copy()
    else:
        inidata = series._series
    if nper < 0:
        nper = max(-len(series), nper)
        newdata[-nper:] = inidata[:nper]
    elif nper > 0:
        nper = min(len(series), nper)
        newdata[:-nper] = inidata[nper:]
    else:
        newdata = inidata
    newseries = newdata.view(type(series))
    newseries._dates = series._dates
    newseries._update_from(series)
    return newseries
TimeSeries.tshift = tshift

#...............................................................................
def _get_type_num_double(dtype):
    """
    Private used to force dtypes upcasting in certain functions
    (eg. int -> float in pct function).
    Adapted from function of the same name in the C source code.
    """
    if dtype.num < np.dtype('f').num:
        return np.dtype('d')
    return dtype

def _pct_generic(series, nper, pct_func):
    "helper function for the pct_* functions"
    _dtype = _get_type_num_double(series.dtype)
    if _dtype != series.dtype:
        series = series.astype(_dtype)
    newdata = masked_array(np.empty(series.shape, dtype=series.dtype),
                           mask=True)
    if nper < newdata.size:
        mseries = series.view(MaskedArray)
        newdata[nper:] = pct_func(mseries, nper)
    newseries = newdata.view(type(series))
    newseries._dates = series._dates
    newseries._update_from(series)
    return newseries

def pct(series, nper=1):
    """
    Returns the rolling percentage change of the series.

    Parameters
    ----------
    series : {TimeSeries}
        TimeSeries object to to calculate percentage chage for. Ignore this
        parameter if calling this as a method.
    nper : {int}
        Number of periods for percentage change.

    Notes
    -----
    Series of integer types will be upcast
    1.0 == 100% in result

    Examples
    --------
    >>> series = ts.time_series(
    ...     [2.,1.,2.,3.], start_date=ts.Date(freq='A', year=2005))
    >>> series.pct()
    timeseries([-- -0.5 1.0 0.5],
               dates = [2005 ... 2008],
               freq  = A-DEC)
    >>> series.pct(2)
    timeseries([-- -- 0.0 2.0],
               dates = [2005 ... 2008],
               freq  = A-DEC)

    """
    def pct_func(series, nper):
        return series[nper:] / series[:-nper] - 1
    return _pct_generic(series, nper, pct_func)
TimeSeries.pct = pct

def pct_log(series, nper=1):
    """
    Returns the rolling log percentage change of the series. This is defined as
    the log of the ratio of series[T]/series[T-nper]

    Parameters
    ----------
    series : {TimeSeries}
        TimeSeries object to to calculate log percentage chage for. Ignore this
        parameter if calling this as a method.
    nper : {int}
        Number of periods for percentage change.

    Notes
    -----
    Series of integer types will be upcast
    1.0 == 100% in result

    Examples
    --------
    >>> series = ts.time_series(
    ...     [2.,1.,2.,3.], start_date=ts.Date(freq='A', year=2005))
    >>> series.pct_log()
    timeseries([-- -0.69314718056 0.69314718056 0.405465108108],
               dates = [2005 ... 2008],
               freq  = A-DEC)
    >>> series.pct_log(2)
    timeseries([-- -- 0.0 1.09861228867],
               dates = [2005 ... 2008],
               freq  = A-DEC)

    """
    def pct_func(series, nper):
        return ma.log(series[nper:] / series[:-nper])
    return _pct_generic(series, nper, pct_func)
TimeSeries.pct_log = pct_log

def pct_symmetric(series, nper=1):
    """
    Returns the rolling symmetric percentage change of the series. This is
    defined as 2*(series[T] - series[T-nper])/(series[T] - series[T-nper])

    Parameters
    ----------
    series : {TimeSeries}
        TimeSeries object to to calculate symmetric percentage chage for. Ignore
        this parameter if calling this as a method.
    nper : {int}
        Number of periods for percentage change.

    Notes
    -----
    Series of integer types will be upcast
    1.0 == 100% in result

    Examples
    --------
    >>> series = ts.time_series(
    ...     [2.,1.,2.,3.], start_date=ts.Date(freq='A', year=2005))
    >>> series.pct_symmetric()
    timeseries([-- -0.666666666667 0.666666666667 0.4],
               dates = [2005 ... 2008],
               freq  = A-DEC)
    >>> series.pct_symmetric(2)
    timeseries([-- -- 0.0 1.0],
               dates = [2005 ... 2008],
               freq  = A-DEC)

    """
    def pct_func(series, nper):
        return \
            2 * (series[nper:] - series[:-nper]) / \
                (series[nper:] + series[:-nper])
    return _pct_generic(series, nper, pct_func)
TimeSeries.pct_symmetric = pct_symmetric



def fill_missing_dates(data, dates=None, freq=None, fill_value=None):
    """
    Finds and fills the missing dates in a time series. The data
    corresponding to the initially missing dates are masked, or filled to
    `fill_value`.

    Parameters
    ----------
    data : {TimeSeries, ndarray}
        Initial array of data.
    dates : {DateArray} (optional)
        Initial array of dates. Specify this if you are passing a plain ndarray
        for the data instead of a :class:`TimeSeries`.
    freq : {freq_spec} (optional)
        Frequency of result. If not specified, the initial frequency is used.
    fill_value : {scalar of type data.dtype} (optional)
        Default value for missing data. If Not specified, the data are just
        masked.

    """
    # Check the frequency ........
    orig_freq = freq
    freq = check_freq(freq)
    if orig_freq is not None and freq == _c.FR_UND:
        freqstr = check_freq_str(freq)
        raise ValueError, \
              "Unable to define a proper date resolution (found %s)." % freqstr
    # Check the dates .............
    if dates is None:
        if not isinstance(data, TimeSeries):
            raise InsufficientDateError
        dates = data._dates
    else:
        if not isinstance(dates, DateArray):
            dates = DateArray(dates, freq)
    dflat = dates.asfreq(freq).ravel()
    if not dflat.has_missing_dates():
        if isinstance(data, TimeSeries):
            return data
        data = data.view(TimeSeries)
        data._dates = dflat
        return data
    # Check the data ..............
    if isinstance(data, MaskedArray):
        datad = data._data
        datam = data._mask
        if isinstance(data, TimeSeries):
            datat = type(data)
            datas = data._varshape
        else:
            datat = TimeSeries
            datas = ()
    else:
        datad = np.asarray(data)
        datam = nomask
        datat = TimeSeries
    # Check whether we need to flatten the data
    if data.ndim > 1:
        if (not datas):
            datad.shape = -1
        elif dflat.size != len(datad):
            err_msg = "fill_missing_dates is not yet implemented for nD series!"
            raise NotImplementedError(err_msg)
    # ...and now, fill it ! ......
    (tstart, tend) = dflat[[0, -1]]
    newdates = date_array(start_date=tstart, end_date=tend)
    (osize, nsize) = (dflat.size, newdates.size)
    #.............................
    # Get the steps between consecutive data.
    delta = dflat.get_steps() - 1
    gap = delta.nonzero()
    slcid = np.concatenate(([0, ], np.arange(1, osize)[gap], [osize, ]))
    oldslc = np.array([slice(i, e)
                       for (i, e) in np.broadcast(slcid[:-1], slcid[1:])])
    addidx = delta[gap].astype(int).cumsum()
    newslc = np.concatenate(([oldslc[0]],
                             [slice(i + d, e + d) for (i, e, d) in \
                              np.broadcast(slcid[1:-1], slcid[2:], addidx)]
                             ))
    #.............................
    # Just a quick check
    vdflat = np.asarray(dflat)
    vnewdates = np.asarray(newdates)
    for (osl, nsl) in zip(oldslc, newslc):
        assert np.equal(vdflat[osl], vnewdates[nsl]).all(), \
            "Slicing mishap ! Please check %s (old) and %s (new)" % (osl, nsl)
    #.............................
    newshape = list(datad.shape)
    newshape[0] = nsize
    newdatad = np.empty(newshape, dtype=data.dtype)
    newdatam = np.ones(newshape, dtype=ma.make_mask_descr(datad.dtype))
    #....
    if datam is nomask:
        for (new, old) in zip(newslc, oldslc):
            newdatad[new] = datad[old]
            newdatam[new] = False
    else:
        for (new, old) in zip(newslc, oldslc):
            newdatad[new] = datad[old]
            newdatam[new] = datam[old]
    if fill_value is None:
        fill_value = getattr(data, '_fill_value', None)
    newdata = ma.masked_array(newdatad, mask=newdatam, fill_value=fill_value)
    _data = newdata.view(datat)
    _data._dates = newdates
    return _data
TimeSeries.fill_missing_dates = fill_missing_dates



def find_duplicated_dates(series):
    """
    Return a dictionary (duplicated dates <> indices) for the input series.

    The indices are given as a tuple of ndarrays, a la :meth:`nonzero`.

    Parameters
    ----------
    series : TimeSeries, DateArray
        A valid :class:`TimeSeries` or :class:`DateArray` object.

    Examples
    --------
    >>> series = time_series(np.arange(10),
                            dates=[2000, 2001, 2002, 2003, 2003,
                                   2003, 2004, 2005, 2005, 2006], freq='A')
    >>> test = find_duplicated_dates(series)
     {<A-DEC : 2003>: (array([3, 4, 5]),), <A-DEC : 2005>: (array([7, 8]),)}
    """
    dates = getattr(series, '_dates', series)
    steps = dates.get_steps()
    duplicated_dates = tuple(set(dates[steps == 0]))
    indices = {}
    for d in duplicated_dates:
        indices[d] = (dates == d).nonzero()
    return indices



def remove_duplicated_dates(series):
    """
    Remove the entries of `series` corresponding to duplicated dates.

    The series is first sorted in chronological order.
    Only the first occurence of a date is then kept, the others are discarded.

    Parameters
    ----------
    series : TimeSeries
        Time series to process
    """
    dates = getattr(series, '_dates', series)
    steps = np.concatenate(([1, ], dates.get_steps()))
    if not dates.is_chronological():
        series = series.copy()
        series.sort_chronologically()
        dates = series._dates
    return series[steps.nonzero()]




def stack(*series):
    """
    Performs a column_stack on the data from each series, and the
    resulting series has the same dates as each individual series. The series
    must have compatible dates (same starting and ending dates, same frequency).

    Parameters
    ----------
    series : the series to be stacked
    """
    _timeseriescompat_multiple(*series)
    return time_series(ma.column_stack(series), series[0]._dates,
                       **_attrib_dict(series[0]))



def concatenate(series, axis=0, remove_duplicates=True, fill_missing=False):
    """
    Joins series together.

    The series are joined in chronological order.
    Duplicated dates are handled with the `remove_duplicates` parameter.
    If `remove_duplicate` is False, duplicated dates are saved.
    Otherwise, only the first occurence of the date is conserved.


    Parameters
    ----------
    series : {sequence}
        Sequence of time series to join
    axis : {0, None, int}, optional
        Axis along which to join
    remove_duplicates : {False, True}, optional
        Whether to remove duplicated dates.
    fill_missing : {False, True}, optional
        Whether to fill the missing dates with missing values.

    Examples
    --------
    >>> a = time_series([1,2,3], start_date=now('D'))
    >>> b = time_series([10,20,30], start_date=now('D')+1)
    >>> c = concatenate((a,b))
    >>> c._series
    masked_array(data = [ 1  2  3 30],
          mask = False,
          fill_value=999999)

    """
    # Get the common frequency, raise an error if incompatibility
    common_f = _compare_frequencies(*series)
    # Concatenate the order of series
    sidx = np.concatenate([np.repeat(i, len(s))
                           for (i, s) in enumerate(series)], axis=axis)
    # Concatenate the dates and data
    ndates = np.concatenate([s._dates for s in series], axis=axis)
    ndata = ma.concatenate([s._series for s in series], axis=axis)
    # Resort the data chronologically
    norder = ndates.argsort(kind='mergesort')
    ndates = ndates[norder]
    ndata = ndata[norder]
    sidx = sidx[norder]
    #
    if not remove_duplicates:
        ndates = date_array(ndates, freq=common_f)
        result = time_series(ndata, dates=ndates)
    else:
        # Find the original dates
        orig = np.concatenate([[True], (np.diff(ndates) != 0)])
        result = time_series(ndata.compress(orig, axis=axis),
                             dates=ndates.compress(orig, axis=axis),
                             freq=common_f)
    if fill_missing:
        result = fill_missing_dates(result)
    return result



def empty_like(series):
    """
    Returns an empty series with the same dtype, mask and dates as series.
    """
    result = np.empty_like(series).view(type(series))
    result._dates = series._dates
    result._mask = series._mask.copy()
    return result

################################################################################
