"""
Extras functions for time series.

:author: Pierre GF Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu - mattknox_ca_at_hotmail_dot_com
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author$)"
__revision__ = "$Revision$"
__date__ = '$Date$'


import numpy as np
import numpy.ma as ma
from numpy.ma import masked

import const as _c
from tdates import Date, date_array, DateArray
from tseries import TimeSeries, time_series
from pandas._skts import DateCalc_Error

from _preview import genfromtxt, easy_dtype

__all__ = ['accept_atmost_missing',
           'convert_to_annual', 'count_missing',
           'guess_freq',
           'isleapyear',
           'tsfromtxt']

#..............................................................................
def isleapyear(year):
    """
    Returns true if year is a leap year.

    Parameters
    ----------
    year : integer / sequence
        A given (list of) year(s).
    """
    year = np.asarray(year)
    return np.logical_or(year % 400 == 0,
                         np.logical_and(year % 4 == 0, year % 100 > 0))

#..............................................................................
def count_missing(series):
    """
    Returns the number of missing data per period.

    Notes
    -----
    This function is designed to return the actual number of missing values when
    a series has been converted from one frequency to a smaller frequency.

    For example, converting a 12-month-long daily series to months will yield
    a (12x31) array, with missing values in February, April, June...
    count_missing will discard these extra missing values.
    """
    if not isinstance(series, TimeSeries):
        raise TypeError, "The input data should be a valid TimeSeries object! "\
                         "(got %s instead)" % type(series)
    if series.ndim == 1:
        return len(series) - series.count()
    elif series.ndim != 2:
        raise NotImplementedError
    #
    missing = series.shape[-1] - series.count(axis= -1)
    period = series.shape[-1]
    freq = series.freq
    if (period == 366) and (freq // _c.FR_ANN == 1):
        # row: years, cols: days
        missing -= ~isleapyear(series.year)
    elif period == 31 and (freq // _c.FR_MTH == 1):
        months = series.months
        # row: months, cols: days
        missing[np.array([m in [4, 6, 9, 11] for m in months])] -= 1
        isfeb = (months == 2)
        missing[isfeb] -= 2
        missing[isfeb & ~isleapyear(series.year)] -= 1
    elif period == 92 and (freq // _c.FR_QTR == 1):
        # row: quarters, cold:days
        months = series.months
        if freq in (_c.FR_QTREJAN, _c.FR_QTRSJAN, _c.FR_QTREAPR, _c.FR_QTRSAPR,
                    _c.FR_QTREOCT, _c.FR_QTRSOCT, _c.FR_QTREOCT, _c.FR_QTRSOCT):
            isfeb = (months == 4)
            missing[isfeb] -= 2
        elif freq in (_c.FR_QTREFEB, _c.FR_QTRSFEB, _c.FR_QTREMAY, _c.FR_QTRSMAY,
                      _c.FR_QTREAUG, _c.FR_QTRSAUG, _c.FR_QTRENOV, _c.FR_QTRSNOV):
            missing[np.array([m in [2, 11] for m in months])] -= 1
            isfeb = (months == 2)
        elif freq in (_c.FR_QTREMAR, _c.FR_QTRSMAR, _c.FR_QTREJUN, _c.FR_QTRSJUN,
                      _c.FR_QTRESEP, _c.FR_QTRSSEP, _c.FR_QTREDEC, _c.FR_QTRSDEC):
            missing[np.array([m in [3, 6] for m in months])] -= 1
            isfeb = (months == 3)
        missing[isfeb & ~isleapyear(series.year)] -= 1
    elif period not in (12, 7):
        raise NotImplementedError, "Not yet implemented for that frequency..."
    return missing



def convert_to_annual(series):
    """
    Group a series by years, taking leap years into account.

    The output has as many rows as distinct years in the original series,
    and as many columns as the length of a leap year in the units corresponding
    to the original frequency (366 for daily frequency, 366*24 for hourly...).
    The fist column of the output corresponds to Jan. 1st, 00:00:00,
    while the last column corresponds to Dec, 31st, 23:59:59.
    Entries corresponding to Feb. 29th are masked for non-leap years.

    For example, if the initial series has a daily frequency, the 59th column
    of the output always corresponds to Feb. 28th, the 61st column to Mar. 1st,
    and the 60th column is masked for non-leap years.
    With a hourly initial frequency, the (59*24)th column of the output always
    correspond to Feb. 28th 23:00, the (61*24)th column to Mar. 1st, 00:00, and
    the 24 columns between (59*24) and (61*24) are masked.

    If the original frequency is less than daily, the output is equivalent to
    ``series.convert('A', func=None)``.


    Parameters
    ----------
    series : TimeSeries
        A valid :class:`~scikits.timeseries.TimeSeries` object.

    Returns
    -------
    aseries : TimeSeries
        A 2D  :class:`~scikits.timeseries.TimeSeries` object with annual ('A')
        frequency.

    """
    freq = series._dates.freq
    if freq < _c.FR_DAY:
        return series.convert('A')
    baseidx = np.array((59, 60), dtype=int)
    if (freq == _c.FR_DAY):
        (idx0228, idx0301) = baseidx
    elif (freq == _c.FR_HR):
        (idx0228, idx0301) = baseidx * 24
    elif (freq == _c.FR_MIN):
        (idx0228, idx0301) = baseidx * 24 * 60
    elif (freq == _c.FR_SEC):
        (idx0228, idx0301) = baseidx * 24 * 3600
    aseries = series.convert('A')
    leapcondition = isleapyear(aseries.dates.years)
    leapidx = np.arange(len(aseries), dtype=int)[~leapcondition]
    aseries[leapidx, idx0301:] = aseries[leapidx, idx0228:idx0228 - idx0301]
    aseries[leapidx, idx0228:idx0301] = ma.masked
    return aseries



#.............................................................................
def accept_atmost_missing(series, max_missing, strict=False):
    """
    Masks the rows of `series` that contain more than `max_missing` missing data.
    Returns a new masked series.

    Parameters
    ----------
    series : TimeSeries
        Input time series.
    max_missing : float
        Number of maximum acceptable missing values per row (if larger than 1),
        or maximum acceptable percentage of missing values (if lower than 1).
    strict : boolean *[False]*
        Whether the number of missing values should be strictly greater than
        `max_missing` or not.

    Returns
    -------
    output : TimeSeries
        A new TimeSeries object
    """
    series = np.array(series, copy=True, subok=True)
    if not isinstance(series, TimeSeries):
        raise TypeError, "The input data should be a valid TimeSeries object! "\
                         "(got %s instead)" % type(series)
    # Find the number of missing values ....
    missing = count_missing(series)
    # Transform an acceptable percentage in a number
    if max_missing < 1:
        max_missing = np.round(max_missing * series.shape[-1], 0)
    #
    series.unshare_mask()
    if strict:
        series[missing > max_missing] = masked
    else:
        series[missing >= max_missing] = masked
    return series


def guess_freq(dates):
    """
    Return an estimate of the frequency from a list of dates.

    The frequency is estimated from the difference of dates (in days or seconds)
    after chronological sorting.


    Parameters
    ----------
    dates : var
        Sequence of dates

    Notes
    -----
    * In practice, the list of dates is first transformed into a list of
      :class:`datetime.datetime` objects.
    """
    if isinstance(dates, DateArray):
        dates = dates.copy()

    try:
        dates = date_array(dates, freq='S', autosort=True)
    except c_dates.DateCalc_Error:
        # contains dates prior to 1979, assume lower frequency
        dates = date_array(dates, freq='D', autosort=True)

    ddif = np.diff(dates)
    mind = np.min(ddif)

    if dates.freq == _c.FR_SEC and mind < 86400:

        # hourly, minutely, or secondly frequency
        if (mind > 3599) and not np.all(ddif % 3600 > 0):
            freq = _c.FR_HR
        elif mind < 59:
            freq = _c.FR_SEC
        else:
            freq = _c.FR_MIN
        return freq

    # daily or lower frequency
    if dates.freq == _c.FR_SEC:
        dates = dates.asfreq('D')
        ddif = np.diff(dates)
        mind = np.min(ddif)

    if mind > 360:
        return _c.FR_ANN

    if mind > 88:
        qincs = [89, 90, 91, 92, 273, 274, 275, 276, 277]
        if np.all([i in qincs for i in (ddif % 365)]):
            freq = _c.FR_QTR
        else:
            freq = _c.FR_MTH
        return freq

    if (mind > 27):
        return _c.FR_MTH

    dow = dates.day_of_week
    if (mind % 7 == 0) and np.all((ddif % 7) == 0):

        mdow = np.min(dow)
        freq = _c.FR_WKSUN + ((mdow + 1) % 7)
        return freq
    else:
        if np.any((dow == 5) | (dow == 6)):
            freq = _c.FR_DAY
        else:
            # no weekends, assume business frequency
            freq = _c.FR_BUS
        return freq




def tsfromtxt(fname, dtype=None, freq='U', comments='#', delimiter=None,
              skip_header=0, skip_footer=0, skiprows=0,
              converters=None, dateconverter=None,
              missing='', missing_values=None, filling_values=None,
              usecols=None, datecols=None,
              names=None, excludelist=None, deletechars=None, autostrip=True,
              case_sensitive=True, defaultfmt="f%i", unpack=None, loose=True,
              asrecarray=False, invalid_raise=True):
    """
    Load a TimeSeries from a text file.

    Each line of the input after the first `skiprows` ones is split at
    `delimiter`. Characters occuring after `comments` are discarded.

    If a column is named ``'dates'`` (case insensitive), it is used to define
    the dates. The ``freq`` parameter should be set to the expected frequency of
    the output series.
    If the date information spans several columns (for example, year in col #1,
    month in col #2...), a specific conversion function must be defined with
    the ``dateconverter`` parameter. This function should accept as many inputs
    as date columns, and return a valid :class:`Date` object.

    Parameters
    ----------
    fname : file or string
        File or filename to read.
        If the file extension is ``.gz`` or ``.bz2``, the file is first
        decompressed.
    dtype : data-type, optional
        Data type of the resulting array.
        If it is a structured data-type, the resulting array is 1-dimensional,
        and each row is interpreted as an element of the array. In this case,
        the number of columns used must match the number of fields in the dtype
        and the names of each field are set by the corresponding name of the dtype.
        If None, the dtypes will be determined by the contents of each
        column, individually.
    comments : {string}, optional
        The character used to indicate the start of a comment.
        All the characters occurring on a line after a comment are discarded.
    delimiter : {string}, optional
        The string used to separate values.  By default, any consecutive
        whitespace act as delimiter.
    skip_header : int, optional
        The numbers of lines to skip at the beginning of the file.
    skip_footer : int, optional
        The numbers of lines to skip at the end of the file
    converters : variable or None, optional
        The set of functions that convert the data of a column to a value.
        The converters can also be used to provide a default value
        for missing data: ``converters = {3: lambda s: float(s or 0)}``.
    dateconverter : {function}, optional
        The function to convert the date information to a :class:`Date` object.
        This function requires as many parameters as number of ``datecols``.
        This parameter is mandatory if ``dtype=None``.
    missing_values : variable or None, optional
        The set of strings corresponding to missing data.
    filling_values : variable or None, optional
        The set of values to be used as default when the data are missing.
    usecols : sequence or None, optional
        Which columns to read, with 0 being the first.  For example,
        ``usecols = (1, 4, 5)`` will extract the 2nd, 5th and 6th columns.
    datecols : {None, int, sequence}, optional
        Which columns store the date information.
    names : {None, True, str, sequence}, optional
        If `names` is True, the field names are read from the first valid line
        after the first `skiprows` lines.
        If `names` is a sequence or a single-string of comma-separated names,
        the names will be used to define the field names in a structured dtype.
        If `names` is None, the names of the dtype fields will be used, if any.
    excludelist : sequence, optional
        A list of names to exclude. This list is appended to the default list
        ['return','file','print']. Excluded names are appended an underscore:
        for example, `file` would become `file_`.
    deletechars : str, optional
        A string combining invalid characters that must be deleted from the
        names.
    defaultfmt : str, optional
        A format used to define default field names, such as "f%i" or "f_%02i".
    autostrip : bool, optional
        Whether to automatically strip white spaces from the variables.
    case_sensitive : {True, False, 'upper', 'lower'}, optional
        If True, field names are case sensitive.
        If False or 'upper', field names are converted to upper case.
        If 'lower', field names are converted to lower case.
    unpack : bool, optional
        If True, the returned array is transposed, so that arguments may be
        unpacked using ``x, y, z = loadtxt(...)``
    usemask : bool, optional
        If True, return a masked array.
        If False, return a regular array.
    asrecarray : {False, True}, optional
        Whether to return a TimeSeriesRecords or a series with a structured
        dtype.
    invalid_raise : bool, optional
        If True, an exception is raised if an inconsistency is detected in the
        number of columns.
        If False, a warning is emitted and the offending lines are skipped.


    Returns
    -------
    out : MaskedArray
        Data read from the text file.

    See Also
    --------
    numpy.lib.io.genfromtxt
        Equivalent function for standard arrays

    Notes
    -----
    * When spaces are used as delimiters, or when no delimiter has been given
      as input, there should not be any missing data between two fields.
    * When the variable are named (either by a flexible dtype or with `names`,
      there must not be any header in the file (else a :exc:`ValueError`
      exception is raised).
    * If ``names`` is True or a sequence of strings, these names overwrite
      the fields names of a structured array.
    * The sequence of names must NOT take the date columns into account.
    * If the datatype is not given explicitly (``dtype=None``),
      a :keyword:`dateconverter` must be given explicitly.
    * If the ``dtype`` is given explicitly,
      it must NOT refer to the date columns.
    * By default, the types of variables is defined from the values encountered
      in the file (``dtype=None``). This is *NOT* the default for np.genfromtxt.

    Examples
    --------
    >>> data = "year, month, a, b\\n 2001, 01, 0.0, 10.\\n 2001, 02, 1.1, 11."
    >>> dateconverter = lambda y, m: Date('M', year=int(y), month=int(m))
    >>> series = tsfromtxt(StringIO.StringIO(data), delimiter=',', names=True,
    ...                    datecols=(0,1), dateconverter=dateconverter,)
    >>> series
    timeseries([(0.0, 10.0) (1.1, 11.0)],
       dtype = [('a', '<f8'), ('b', '<f8')],
       dates = [Jan-2001 Feb-2001],
       freq  = M)
    >>> series = tsfromtxt(StringIO.StringIO(data), delimiter=",",
    ...                    datecols=(0, 1), dateconverter=dateconverter,
    ...                    names="A, B", skip_header=1)
    timeseries([(0.0, 10.0) (1.1000000000000001, 11.0)],
       dtype = [('A', '<f8'), ('B', '<f8')],
       dates = [Jan-2001 Feb-2001],
       freq  = M)

    """
    # Update the date converter ...........................
    converters = converters or {}
    dateconv = dateconverter or None
    if dateconv is None:
        dateconv = lambda s: Date(freq, string=s)
    if 'dates' in converters:
        dateconv = converters['dates']
        del(converters['dates'])

    # Make sure `datecols` is a sequence ..................
    if datecols is not None:
        try:
            datecols = [_.strip() for _ in datecols.split(",")]
        except AttributeError:
            try:
                datecols = list(datecols)
            except TypeError:
                datecols = [datecols, ]
        # ... and update the converters
        converters.update((i, str) for i in datecols)

    # Save the initial names and dtypes ...................
    idtype = dtype
    if isinstance(names, basestring):
        names = names.split(",")
    inames = names

    # Update the dtype (if needed) ........................
    if (dtype is not None):
        # Crash if we can't find the datecols
        if datecols is None:
            raise TypeError("No column selected for the dates!")
        # Make sure dtype is a valid np.dtype and make a copy
        dtype = easy_dtype(dtype, names=names)
        idtype = dtype
        inames = dtype.names
        if inames is not None:
            nbfields = len(inames) + len(datecols)
            # Create a new dtype description and a set of names
            dtype = [''] * nbfields
            names = [''] * nbfields
            idx = range(nbfields)
            for i in datecols:
                if i < 0:
                    i += nbfields
                del idx[idx.index(i)]
                # Set the default dtype for date columns, as np.object
                # (we can't use string as we don't know the final size)
                dtype[i] = ('', np.object)
            convdict = {'b': bool, 'i': int, 'l':int, 'u': int,
                        'f': float, 'd': float, 'g': float,
                        'c': complex, 'D': complex,
                        'S': str, 'U': str, 'a': str}
            converter_update = []
            for (i, name) in zip(idx, inames):
                field = idtype[name]
                dtype[i] = (name, field)
                converter_update.append((i, convdict[field.char]))
                names[i] = name
            converters.update(converter_update)
    elif names not in (True, None):
        # Make sure that we saved the names as a list
        inames = list(inames)
        # Get the list of columns to use
        if usecols is None:
            nbcols = len(datecols) + len(inames)
            names = [''] * nbcols
            ucols = range(nbcols)
        else:
            names = [''] * (max(usecols) + 1)
            ucols = usecols
        # Fill the list of names:
        for i in ucols:
            if i in datecols:
                names[i] = "__%i" % i
            else:
                names[i] = inames.pop(0)
    #
    # Update the optional arguments ...
    kwargs = dict(dtype=dtype, comments=comments, delimiter=delimiter,
                  skiprows=skiprows, converters=converters,
                  skip_header=skip_header, skip_footer=skip_footer,
                  missing=missing, missing_values=missing_values,
                  filling_values=filling_values,
                  usecols=usecols, unpack=unpack, names=names,
                  excludelist=excludelist, deletechars=deletechars,
                  case_sensitive=case_sensitive, defaultfmt=defaultfmt,
                  autostrip=autostrip, loose=loose, invalid_raise=invalid_raise,
                  usemask=True)
    # Get the raw data ................
    mrec = genfromtxt(fname, **kwargs)
    if not mrec.shape:
        mrec.shape = -1
    names = mrec.dtype.names
    # Revert to the original dtype.....
    dtype = idtype
    # Get the date columns ................................
    if datecols is None:
        import re
        datespattern = re.compile("'?_?dates?'?", re.IGNORECASE)
        datecols = [i for (i, name) in enumerate(names or ())
                     if datespattern.search(name)]
        if not datecols:
            raise TypeError("No column selected for the dates!")
    else:
        # We have `datecols` already, make sure the indices are positive
        # (the nb of fields might still be undefined)
        nbfields = len(names)
        for (i, v) in enumerate(datecols):
            if (v < 0):
                datecols[i] = v + nbfields
    # Fix the date columns if usecols was given
    if usecols is not None:
        datecols = tuple([list(usecols).index(d) for d in datecols])
    # Get the date info ...............
    if names:
        _dates = [mrec[names[i]] for i in datecols]
    else:
        _dates = [mrec[:, i] for i in datecols]
    # Convert the date columns to a date_array
    if len(_dates) == 1:
        _dates = np.array(_dates[0], copy=False, ndmin=1)
        dates = date_array([dateconv(args) for args in _dates],
                           freq=freq, autosort=False)
    else:
        dates = date_array([dateconv(*args) for args in zip(*_dates)],
                           freq=freq, autosort=False)
    # Resort the array according to the dates
    sortidx = dates.argsort()
    dates = dates[sortidx]
    mrec = mrec[sortidx]
    # Get the dtype from the named columns (if any), or just use the initial one
    mdtype = mrec.dtype
    if mdtype.names:
        newdescr = [descr for (i, descr) in enumerate(mdtype.descr)
                    if i not in datecols]
        output = time_series(ma.empty((len(mrec),), dtype=newdescr),
                             dates=dates)
        for name in output.dtype.names:
            output[name] = mrec[name]
        if (idtype is not None):
            if (idtype.names is None):
                dtype = (idtype, len(output.dtype.names))
            else:
                dtype = idtype
            output = output.view(dtype)
    else:
        dataidx = [i for i in range(mrec.shape[-1]) if i not in datecols]
        if len(dataidx) == 1:
            dataidx = dataidx[0]
        output = time_series(mrec[:, dataidx], dates=dates)
    #
    if asrecarray:
        from trecords import TimeSeriesRecords
        return output.view(TimeSeriesRecords)
    return output
