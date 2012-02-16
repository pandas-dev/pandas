# pylint: disable-msg=W0201, W0212
"""
Support for multi-variable time series, through masked record arrays.

Individual fields can be accessed as keys or attributes.

:author: Pierre GF Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu - mattknox_ca_at_hotmail_dot_com
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author$)"
__revision__ = "$Revision$"
__date__ = '$Date$'


import sys

import numpy as np
from numpy import bool_, complex_, float_, int_, str_, object_, \
    ndarray, chararray, recarray
import numpy.core.numerictypes as ntypes
import numpy.core.umath as umath
from numpy.core.records import find_duplicate, format_parser, record, \
    fromarrays as recfromarrays

import numpy.ma as ma
from numpy.ma import MaskedArray, MAError, \
     default_fill_value, masked_print_option, masked, nomask, \
     getmask, getmaskarray, make_mask, make_mask_none, mask_or, filled

import numpy.ma.mrecords
from numpy.ma.mrecords import _checknames, \
     _guessvartypes, openfile, MaskedRecords, mrecarray, addfield, \
     fromrecords as mrecfromrecords, fromarrays as mrecfromarrays

from tseries import TimeSeries, TimeSeriesCompatibilityError, \
    time_series, _getdatalength, nodates, get_varshape
from tdates import Date, DateArray, date_array
from extras import tsfromtxt

_byteorderconv = numpy.core.records._byteorderconv
_typestr = ntypes._typestr

reserved_fields = numpy.ma.mrecords.reserved_fields + ['_dates']

import warnings

__all__ = [
'TimeSeriesRecords', 'time_records',
'fromarrays', 'fromrecords', 'fromtextfile',
]

def _getformats(data):
    """
    Returns the formats of each array of arraylist as a comma-separated string.
    """
    if isinstance(data, record):
        return ",".join([desc[1] for desc in data.dtype.descr])

    formats = ''
    for obj in data:
        obj = np.asarray(obj)
        formats += _typestr[obj.dtype.type]
        if issubclass(obj.dtype.type, ntypes.flexible):
            formats += `obj.itemsize`
        formats += ','
    return formats[:-1]


def _getdates(dates=None, newdates=None, length=None, freq=None,
              start_date=None):
    """
    Determines new dates (private function not meant to be used).
    """
    if dates is None:
        if newdates is not None:
            if not hasattr(newdates, 'freq'):
                newdates = date_array(dlist=newdates, freq=freq)
        else:
            newdates = date_array(start_date=start_date, length=length,
                                  freq=freq)
    elif not hasattr(dates, 'freq'):
        newdates = date_array(dlist=dates, freq=freq)
    else:
        newdates = dates
    return newdates


class TimeSeriesRecords(TimeSeries, MaskedRecords, object):
    """
    MaskedRecords with support for time-indexing.

    Fields can be retrieved either as indices (using the indexing scheme
    based on `__getitem__`) or as attributes (using `__getattribute__`).

    The type of the output of `__getitem__` is variable:
    
    field
       returns a :class:`TimeSeries` object.
    single record, no masked fields
        returns a ``numpy.void`` object
    single record with at least one masked field
        returns a :class:`MaskedRecords` object.
    slice
        return a :class:`TimeSeriesRecords`.
    """
    def __new__(cls, shape, dtype=None, buf=None, offset=0, strides=None,
                formats=None, names=None, titles=None,
                byteorder=None, aligned=False,
                mask=nomask, hard_mask=False, fill_value=None, keep_mask=True,
                copy=False,
                dates=None, freq='U', start_date=None, observed=None,
                **options):
        _data = mrecarray.__new__(cls, shape, dtype=dtype, buf=buf, offset=offset,
                                  strides=strides, formats=formats,
                                  byteorder=byteorder, aligned=aligned,
                                  mask=mask, hard_mask=hard_mask, copy=copy,
                                  keep_mask=keep_mask, fill_value=fill_value,
                                  )
        #
        newdates = _getdates(dates, length=len(_data),
                             start_date=start_date, freq=freq)
        _data._dates = newdates
        _data._observed = observed
        #
        return _data

    def __array_finalize__(self, obj):
        self.__dict__.update(_varshape=getattr(obj, '_varshape', ()),
                             _dates=getattr(obj, '_dates', DateArray([])),
                             _observed=getattr(obj, '_observed', None),
                             _optinfo=getattr(obj, '_optinfo', {}))
        MaskedRecords.__array_finalize__(self, obj)
        return


    def _getdata(self):
        "Returns the data as a recarray."
        return ndarray.view(self, recarray)
    _data = property(fget=_getdata)

    def _getseries(self):
        "Returns the data as a MaskedRecord array."
        return MaskedArray.view(self, mrecarray)
    _series = property(fget=_getseries)


    def __getattribute__(self, attr):
        getattribute = MaskedRecords.__getattribute__
        _dict = getattribute(self, '__dict__')
        if attr == '_dict':
            return _dict
        _names = ndarray.__getattribute__(self, 'dtype').names
        if attr in (_names or []):
            obj = getattribute(self, attr).view(TimeSeries)
            obj._dates = _dict['_dates']
            return obj
        return getattribute(self, attr)


    def __setattr__(self, attr, value):
        if attr in ['_dates', 'dates']:
            self.__setdates__(value)
        elif attr == 'shape':
            if self._varshape:
                err_msg = "Reshaping a nV/nD series is not implemented yet !"
                raise NotImplementedError(err_msg)
        return MaskedRecords.__setattr__(self, attr, value)

    #......................................................
    def __getitem__(self, indx):
        """Returns all the fields sharing the same fieldname base.
    The fieldname base is either `_data` or `_mask`."""
        _localdict = self.__dict__
        # We want a field ........
        if indx in ndarray.__getattribute__(self, 'dtype').names:
            obj = self._data[indx].view(TimeSeries)
            obj._dates = _localdict['_dates']
            obj._mask = make_mask(_localdict['_mask'][indx])
            return obj
        # We want some elements ..
        obj = TimeSeries.__getitem__(self, indx)
        if isinstance(obj, MaskedArray) and not isinstance(obj, TimeSeries):
            obj = ndarray.view(obj, MaskedRecords)
        return obj


    def __setslice__(self, i, j, value):
        """Sets the slice described by [i,j] to `value`."""
        MaskedRecords.__setitem__(self, slice(i, j), value)
        return

    #......................................................
    def __str__(self):
        """x.__str__() <==> str(x)
    Calculates the string representation, using masked for fill if it is enabled.
    Otherwise, fills with fill value.
        """
        if self.size > 1:
            mstr = ["(%s)" % ",".join([str(i) for i in s])
                    for s in zip(*[getattr(self, f)._series
                                   for f in self.dtype.names])]
            return "[%s]" % ", ".join(mstr)
        else:
            mstr = ["%s" % ",".join([str(i) for i in s])
                    for s in zip([getattr(self, f)._series
                                  for f in self.dtype.names])]
            return "(%s)" % ", ".join(mstr)

    def __repr__(self):
        """x.__repr__() <==> repr(x)
    Calculates the repr representation, using masked for fill if it is enabled.
    Otherwise fill with fill value.
        """
        _names = self.dtype.names
        _dates = self._dates
        if np.size(_dates) > 2 and self._dates.is_valid():
            timestr = "[%s ... %s]" % (str(_dates[0]), str(_dates[-1]))
        else:
            timestr = str(_dates)
        fmt = "%%%is : %%s" % (max([len(n) for n in _names]) + 4,)
        reprstr = [fmt % (f, getattr(self, f)) for f in self.dtype.names]
        reprstr.insert(0, 'TimeSeriesRecords(')
        reprstr.extend([fmt % ('dates', timestr),
                        fmt % ('    fill_value', self.fill_value),
                         '               )'])
        return str("\n".join(reprstr))


    def copy(self):
        "Returns a copy of the argument."
        copied = MaskedRecords.copy(self)
        copied._dates = self._dates.copy()
        return copied


    def convert(self, freq, func=None, position='END', *args, **kwargs):
        """
    Converts a series to another frequency.

    Parameters
    ----------
    series : TimeSeries
        the series to convert. Skip this parameter if you are calling this as
        a method of the TimeSeries object instead of the module function.
    freq : freq_spec
        Frequency to convert the TimeSeries to. Accepts any valid frequency
        specification (string or integer)
    func : {None,function}, optional
        When converting to a lower frequency, `func` is a function that acts on
        one date's worth of data. `func` should handle masked values appropriately.
        If `func` is None, then each entry of the resulting series is the group
        of data points that fall into the date at the lower frequency.
        For example, if converting from monthly to daily and you wanted each
        data point in the resulting series to be the average value for each
        month, you could specify numpy.ma.average for the 'func' parameter.
    position : {'END', 'START'}, optional
        When converting to a higher frequency, position is 'START' or 'END'
        and determines where the data point is in each period. For example, if
        going from monthly to daily, and position is 'END', then each data
        point is placed at the end of the month.
    *args : {extra arguments for func parameter}, optional
        if a func is specified that requires additional parameters, specify
        them here.
    **kwargs : {extra keyword arguments for func parameter}, optional
        if a func is specified that requires additional keyword parameters,
        specify them here.

        """
        kwargs.update(func=func, position=position)
        field_names = self.dtype.names
        by_field = [self[f].convert(freq, **kwargs) for f in field_names]
        output = fromarrays(by_field,
                            dates=by_field[0].dates,
                            names=field_names)
        output.fill_value = self._fill_value
        return output
trecarray = TimeSeriesRecords


#####---------------------------------------------------------------------------
#---- --- Constructors ---
#####---------------------------------------------------------------------------

def time_records(data, dates=None, start_date=None, freq=None, mask=nomask,
                dtype=None, copy=False, fill_value=None, keep_mask=True,
                hard_mask=False):
    """
    Creates a TimeSeriesRecords object.

    Parameters
    ----------
    data : array_like
        Data portion of the array. Any data that is valid for constructing a
        MaskedArray can be used here. May also be a TimeSeries object.
    dates : {None, DateArray}, optional
        A sequence of dates corresponding to each entry.
        If None, the dates will be constructed as a DateArray with the same
        length as ``data``, starting at ``start_date`` with frequency ``freq``.
    start_date : {Date}, optional
        Date corresponding to the first entry of the data (index 0).
        This parameter must be a valid Date object, and is mandatory if ``dates``
        is None and if ``data`` has a length greater or equal to 1.
    freq : {freq_spec}, optional
        A valid frequency specification, as a string or an integer.
        This parameter is mandatory if ``dates`` is None.
    mask : {nomask, sequence}, optional
        Mask.  Must be convertible to an array of booleans with
        the same shape as data: True indicates a masked (eg.,
        invalid) data.
    dtype : {dtype}, optional
        Data type of the output.
        If dtype is None, the type of the data argument (`data.dtype`) is used.
        If dtype is not None and different from `data.dtype`, a copy is performed.
    copy : {False, True}, optional
        Whether to copy the input data (True), or to use a reference instead.
        Note: data are NOT copied by default.
    fill_value : {var}, optional
        Value used to fill in the masked values when necessary.
        If None, a default based on the datatype is used.
    keep_mask : {True, boolean}, optional
        Whether to combine mask with the mask of the input data,
        if any (True), or to use only mask for the output (False).
    hard_mask : {False, boolean}, optional
        Whether to use a hard mask or not.
        With a hard mask, masked values cannot be unmasked.

    Notes
    -----
    * All other parameters that are accepted by the :func:`numpy.ma.array`
      function in the :mod:`numpy.ma` module are also accepted by this function.
    * The date portion of the time series must be specified in one of the
      following ways:

       * specify a TimeSeries object for the ``data`` parameter.
       * pass a DateArray for the ``dates`` parameter.
       * specify a start_date (a continuous DateArray will be automatically
         constructed for the dates portion).
       * specify just a frequency (for TimeSeries of size zero).

    """
    series = time_series(data, dates=dates, start_date=start_date, freq=freq,
                          mask=mask, dtype=dtype, copy=copy,
                          fill_value=fill_value, keep_mask=keep_mask,
                          hard_mask=hard_mask)
    return series.view(TimeSeriesRecords)

#!!!: * The docstrings of the following functions need some serious work ;)
#!!!: * We should try to have a list of TimeSeries sufficient to build a record...
#!!!:   without having to precise a list of dates...
#!!!:   > check the compatibility of dates
#!!!:   > try to adjust endpoints if needed
#!!!:   > if one of the series is not a TimeSeries, keep going.

def fromarrays(arraylist, dates=None, start_date=None, freq='U',
               fill_value=None, autosort=True,
               dtype=None, shape=None, formats=None,
               names=None, titles=None, aligned=False, byteorder=None,):
    """
    Creates a mrecarray from a (flat) list of masked arrays.

    Parameters
    ----------
    arraylist : array_like
        A list of (masked) arrays. Each element of the sequence is first converted
        to a masked array if needed. If a 2D array is passed as argument, it is
        processed line by line
    dates : {DateArray}, optional
        Array of dates corresponding to each entry.
        If None, a DateArray is constructed from `start_date` and the length
        of the arrays in the input list.
    start_date : {Date}, optional
        First date of the output.
        This parameter is inly needed if `dates` is None.
    freq : {var}, optional
        Frequency of the DateArray
    fill_value : {var}, optional
        Value used to fill in the masked values when necessary.
        If None, a default based on the datatype is used.
    autosort : {True, False}, optional
        Whether the records should be sorted chronologically.

    See Also
    --------
    numpy.core.records.fromarrays : equivalent function for ndarrays
        The docstring of this function describes the additional optional
        input parameters.
    

    Notes
    -----
    * Lists of tuples should be preferred over lists of lists as inputs 
      for faster processing.
    """
    _array = mrecfromarrays(arraylist, dtype=dtype, shape=shape, formats=formats,
                            names=names, titles=titles, aligned=aligned,
                            byteorder=byteorder, fill_value=fill_value)
    _dates = _getdates(dates, length=len(_array), start_date=start_date,
                       freq=freq)
#    if _dates._unsorted is not None:
#        idx = _dates._unsorted
#        _array = _array[idx]
#        _dates._unsorted = None
    result = _array.view(trecarray)
    result._dates = _dates
    if autosort:
        result.sort_chronologically()
    return result


#..............................................................................
def fromrecords(reclist, dates=None, freq=None, start_date=None,
                fill_value=None, mask=nomask, autosort=True,
                dtype=None, shape=None, formats=None, names=None,
                titles=None, aligned=False, byteorder=None):
    """
    Creates a TimeSeriesRecords from a list of records.

    The data in the same field can be heterogeneous, they will be promoted
    to the highest data type.  This method is intended for creating
    smaller record arrays.  If used to create large array without formats
    defined, it can be slow.

    If formats is None, then this will auto-detect formats. Use a list of
    tuples rather than a list of lists for faster processing.


    Parameters
    ----------
    reclist : array_like
        A list of records. Each element of the sequence is first converted
        to a masked array if needed. If a 2D array is passed as argument, it is
        processed line by line
    dates : {DateArray}, optional
        Array of dates corresponding to each entry.
        If None, a DateArray is constructed from `start_date` and the length
        of the arrays in the input list.
    freq : {var}, optional
        Frequency of the DateArray
    start_date : {Date}, optional
        First date of the output.
        This parameter is inly needed if `dates` is None.
    fill_value : {var}, optional
        Value used to fill in the masked values when necessary.
        If None, a default based on the datatype is used.
    autosort : {True, False}, optional
        Whether the records should be sorted chronologically.
        

    See Also
    --------
    numpy.core.records.fromrecords : equivalent function for ndarrays


    """
    _data = mrecfromrecords(reclist, dtype=dtype, shape=shape, formats=formats,
                            names=names, titles=titles, aligned=aligned,
                            byteorder=byteorder, mask=mask)
    _dtype = _data.dtype
    # Check the names for a '_dates' .................
    newdates = None
    _names = list(_dtype.names)
    reserved = [n for n in _names if n.lower() in ['dates', '_dates']]
    if len(reserved) > 0:
        newdates = _data[reserved[-1]]
        [_names.remove(n) for n in reserved]
        _dtype = np.dtype([t for t in _dtype.descr \
                                    if t[0] not in reserved ])
        _data = mrecfromarrays([_data[n] for n in _names], dtype=_dtype)
    #
    if dates is None:
        dates = getattr(reclist, '_dates', None)
    _dates = _getdates(dates=dates, newdates=newdates, length=len(_data),
                       freq=freq, start_date=start_date)
    #
    result = _data.view(trecarray)
    result._dates = _dates
    if autosort:
        result.sort_chronologically()
    return result



def fromtextfile(fname, delimitor=None, commentchar='#', missingchar='',
                 dates_column=None, varnames=None, vartypes=None,
                 dates=None, freq=None, skiprows=0):
    """
    Deprecated function: please use tsfromtxt instead.
    """
    msg = "This function is deprecated.\nPlease use `tsfromtxt` instead."
    warnings.warn(msg, DeprecationWarning)
    # Split the names by comma first, then per character
    if isinstance(varnames, basestring):
        varnames = varnames.split(",")
        if len(varnames) == 1:
            varnames = zip(*varnames[0])[0]
    return tsfromtxt(fname, dtype=vartypes, names=(varnames or True), freq=freq,
                     datecols=dates_column, skiprows=skiprows,
                     delimiter=delimitor, comments=commentchar,
                     missing=missingchar, asrecarray=True)

