from datetime import datetime, timedelta
import re
import sys

import numpy as np

import pandas.lib as lib
import pandas.tslib as tslib
import pandas.core.common as com
from pandas.compat import StringIO, callable
from pandas.core.common import ABCIndexClass
import pandas.compat as compat
from pandas.util.decorators import deprecate_kwarg

try:
    import dateutil
    # raise exception if dateutil 2.0 install on 2.x platform
    if (sys.version_info[0] == 2 and
            dateutil.__version__ == '2.0'):  # pragma: no cover
        raise Exception('dateutil 2.0 incompatible with Python 2.x, you must '
                        'install version 1.5 or 2.1+!')
except ImportError:  # pragma: no cover
    print('Please install python-dateutil via easy_install or some method!')
    raise  # otherwise a 2nd import won't show the message

_DATEUTIL_LEXER_SPLIT = None
try:
    # Since these are private methods from dateutil, it is safely imported
    # here so in case this interface changes, pandas will just fallback
    # to not using the functionality
    from dateutil.parser import _timelex

    if hasattr(_timelex, 'split'):
        def _lexer_split_from_str(dt_str):
            # The StringIO(str(_)) is for dateutil 2.2 compatibility
            return _timelex.split(StringIO(str(dt_str)))

        _DATEUTIL_LEXER_SPLIT = _lexer_split_from_str
except (ImportError, AttributeError):
    pass

def _infer_tzinfo(start, end):
    def _infer(a, b):
        tz = a.tzinfo
        if b and b.tzinfo:
            if not (tslib.get_timezone(tz) == tslib.get_timezone(b.tzinfo)):
                raise AssertionError('Inputs must both have the same timezone,'
                                     ' {0} != {1}'.format(tz, b.tzinfo))
        return tz
    tz = None
    if start is not None:
        tz = _infer(start, end)
    elif end is not None:
        tz = _infer(end, start)
    return tz


def _guess_datetime_format(dt_str, dayfirst=False,
                           dt_str_parse=compat.parse_date,
                           dt_str_split=_DATEUTIL_LEXER_SPLIT):
    """
    Guess the datetime format of a given datetime string.

    Parameters
    ----------
    dt_str : string, datetime string to guess the format of
    dayfirst : boolean, default False
        If True parses dates with the day first, eg 20/01/2005
        Warning: dayfirst=True is not strict, but will prefer to parse
        with day first (this is a known bug).
    dt_str_parse : function, defaults to `compate.parse_date` (dateutil)
        This function should take in a datetime string and return
        a `datetime.datetime` guess that the datetime string represents
    dt_str_split : function, defaults to `_DATEUTIL_LEXER_SPLIT` (dateutil)
        This function should take in a datetime string and return
        a list of strings, the guess of the various specific parts
        e.g. '2011/12/30' -> ['2011', '/', '12', '/', '30']

    Returns
    -------
    ret : datetime formatt string (for `strftime` or `strptime`)
    """
    if dt_str_parse is None or dt_str_split is None:
        return None

    if not isinstance(dt_str, compat.string_types):
        return None

    day_attribute_and_format = (('day',), '%d', 2)

    # attr name, format, padding (if any)
    datetime_attrs_to_format = [
        (('year', 'month', 'day'), '%Y%m%d', 0),
        (('year',), '%Y', 0),
        (('month',), '%B', 0),
        (('month',), '%b', 0),
        (('month',), '%m', 2),
        day_attribute_and_format,
        (('hour',), '%H', 2),
        (('minute',), '%M', 2),
        (('second',), '%S', 2),
        (('microsecond',), '%f', 6),
        (('second', 'microsecond'), '%S.%f', 0),
    ]

    if dayfirst:
        datetime_attrs_to_format.remove(day_attribute_and_format)
        datetime_attrs_to_format.insert(0, day_attribute_and_format)

    try:
        parsed_datetime = dt_str_parse(dt_str, dayfirst=dayfirst)
    except:
        # In case the datetime can't be parsed, its format cannot be guessed
        return None

    if parsed_datetime is None:
        return None

    try:
        tokens = dt_str_split(dt_str)
    except:
        # In case the datetime string can't be split, its format cannot
        # be guessed
        return None

    format_guess = [None] * len(tokens)
    found_attrs = set()

    for attrs, attr_format, padding in datetime_attrs_to_format:
        # If a given attribute has been placed in the format string, skip
        # over other formats for that same underlying attribute (IE, month
        # can be represented in multiple different ways)
        if set(attrs) & found_attrs:
            continue

        if all(getattr(parsed_datetime, attr) is not None for attr in attrs):
            for i, token_format in enumerate(format_guess):
                token_filled = tokens[i].zfill(padding)
                if (token_format is None and
                        token_filled == parsed_datetime.strftime(attr_format)):
                    format_guess[i] = attr_format
                    tokens[i] = token_filled
                    found_attrs.update(attrs)
                    break

    # Only consider it a valid guess if we have a year, month and day
    if len(set(['year', 'month', 'day']) & found_attrs) != 3:
        return None

    output_format = []
    for i, guess in enumerate(format_guess):
        if guess is not None:
            # Either fill in the format placeholder (like %Y)
            output_format.append(guess)
        else:
            # Or just the token separate (IE, the dashes in "01-01-2013")
            try:
                # If the token is numeric, then we likely didn't parse it
                # properly, so our guess is wrong
                float(tokens[i])
                return None
            except ValueError:
                pass

            output_format.append(tokens[i])

    guessed_format = ''.join(output_format)

    # rebuild string, capturing any inferred padding
    dt_str = ''.join(tokens)
    if parsed_datetime.strftime(guessed_format) == dt_str:
        return guessed_format

def _guess_datetime_format_for_array(arr, **kwargs):
    # Try to guess the format based on the first non-NaN element
    non_nan_elements = com.notnull(arr).nonzero()[0]
    if len(non_nan_elements):
        return _guess_datetime_format(arr[non_nan_elements[0]], **kwargs)


@deprecate_kwarg(old_arg_name='coerce', new_arg_name='errors',
                 mapping={True: 'coerce', False: 'raise'})
def to_datetime(arg, errors='raise', dayfirst=False, yearfirst=False,
                utc=None, box=True, format=None, exact=True, coerce=None,
                unit='ns', infer_datetime_format=False):
    """
    Convert argument to datetime.

    Parameters
    ----------
    arg : string, datetime, array of strings (with possible NAs)
    errors : {'ignore', 'raise', 'coerce'}, default 'raise'
        - If 'raise', then invalid parsing will raise an exception
        - If 'coerce', then invalid parsing will be set as NaT
        - If 'ignore', then invalid parsing will return the input
    dayfirst : boolean, default False
        Specify a date parse order if `arg` is str or its list-likes.
        If True, parses dates with the day first, eg 10/11/12 is parsed as 2012-11-10.
        Warning: dayfirst=True is not strict, but will prefer to parse
        with day first (this is a known bug, based on dateutil behavior).
    yearfirst : boolean, default False
        Specify a date parse order if `arg` is str or its list-likes.
        - If True parses dates with the year first, eg 10/11/12 is parsed as 2010-11-12.
        - If both dayfirst and yearfirst are True, yearfirst is preceded (same as dateutil).
        Warning: yearfirst=True is not strict, but will prefer to parse
        with year first (this is a known bug, based on dateutil beahavior).

        .. versionadded: 0.16.1

    utc : boolean, default None
        Return UTC DatetimeIndex if True (converting any tz-aware
        datetime.datetime objects as well).
    box : boolean, default True
        - If True returns a DatetimeIndex
        - If False returns ndarray of values.
    format : string, default None
        strftime to parse time, eg "%d/%m/%Y", note that "%f" will parse
        all the way up to nanoseconds.
    exact : boolean, True by default
        - If True, require an exact format match.
        - If False, allow the format to match anywhere in the target string.
    unit : unit of the arg (D,s,ms,us,ns) denote the unit in epoch
        (e.g. a unix timestamp), which is an integer/float number.
    infer_datetime_format : boolean, default False
        If no `format` is given, try to infer the format based on the first
        datetime string. Provides a large speed-up in many cases.

    Returns
    -------
    ret : datetime if parsing succeeded.
        Return type depends on input:

        - list-like: DatetimeIndex
        - Series: Series of datetime64 dtype
        - scalar: Timestamp

        In case when it is not possible to return designated types (e.g. when
        any element of input is before Timestamp.min or after Timestamp.max)
        return will have datetime.datetime type (or correspoding array/Series).

    Examples
    --------
    Take separate series and convert to datetime

    >>> import pandas as pd
    >>> i = pd.date_range('20000101',periods=100)
    >>> df = pd.DataFrame(dict(year = i.year, month = i.month, day = i.day))
    >>> pd.to_datetime(df.year*10000 + df.month*100 + df.day, format='%Y%m%d')
    0    2000-01-01
    1    2000-01-02
    ...
    98   2000-04-08
    99   2000-04-09
    Length: 100, dtype: datetime64[ns]

    Or from strings

    >>> df = df.astype(str)
    >>> pd.to_datetime(df.day + df.month + df.year, format="%d%m%Y")
    0    2000-01-01
    1    2000-01-02
    ...
    98   2000-04-08
    99   2000-04-09
    Length: 100, dtype: datetime64[ns]

    Date that does not meet timestamp limitations:

    >>> pd.to_datetime('13000101', format='%Y%m%d')
    datetime.datetime(1300, 1, 1, 0, 0)
    >>> pd.to_datetime('13000101', format='%Y%m%d', errors='coerce')
    NaT
    """
    return _to_datetime(arg, errors=errors, dayfirst=dayfirst, yearfirst=yearfirst,
                        utc=utc, box=box, format=format, exact=exact,
                        unit=unit, infer_datetime_format=infer_datetime_format)


def _to_datetime(arg, errors='raise', dayfirst=False, yearfirst=False,
                 utc=None, box=True, format=None, exact=True,
                 unit='ns', freq=None, infer_datetime_format=False):
    """
    Same as to_datetime, but accept freq for
    DatetimeIndex internal construction
    """
    from pandas.core.series import Series
    from pandas.tseries.index import DatetimeIndex

    def _convert_listlike(arg, box, format, name=None):

        if isinstance(arg, (list,tuple)):
            arg = np.array(arg, dtype='O')

        # these are shortcutable
        if com.is_datetime64_ns_dtype(arg):
            if box and not isinstance(arg, DatetimeIndex):
                try:
                    return DatetimeIndex(arg, tz='utc' if utc else None, name=name)
                except ValueError:
                    pass

            return arg

        elif com.is_datetime64tz_dtype(arg):
            if not isinstance(arg, DatetimeIndex):
                return DatetimeIndex(arg, tz='utc' if utc else None)
            if utc:
                arg = arg.tz_convert(None)
            return arg

        elif format is None and com.is_integer_dtype(arg) and unit=='ns':
            result = arg.astype('datetime64[ns]')
            if box:
                return DatetimeIndex(result, tz='utc' if utc else None, name=name)

            return result

        arg = com._ensure_object(arg)
        require_iso8601 = False

        if infer_datetime_format and format is None:
            format = _guess_datetime_format_for_array(arg, dayfirst=dayfirst)

        if format is not None:
            # There is a special fast-path for iso8601 formatted
            # datetime strings, so in those cases don't use the inferred
            # format because this path makes process slower in this
            # special case
            format_is_iso8601 = (
                ('%Y-%m-%dT%H:%M:%S.%f'.startswith(format) or
                '%Y-%m-%d %H:%M:%S.%f'.startswith(format)) and
				format != '%Y'
            )
            if format_is_iso8601:
                require_iso8601 = not infer_datetime_format
                format = None

        try:
            result = None

            if format is not None:
                # shortcut formatting here
                if format == '%Y%m%d':
                    try:
                        result = _attempt_YYYYMMDD(arg, errors=errors)
                    except:
                        raise ValueError("cannot convert the input to '%Y%m%d' date format")

                # fallback
                if result is None:
                    try:
                        result = tslib.array_strptime(
                            arg, format, exact=exact, errors=errors)
                    except (tslib.OutOfBoundsDatetime):
                        if errors == 'raise':
                            raise
                        result = arg
                    except ValueError:
                        # if format was inferred, try falling back
                        # to array_to_datetime - terminate here
                        # for specified formats
                        if not infer_datetime_format:
                            if errors == 'raise':
                                raise
                            result = arg

            if result is None and (format is None or infer_datetime_format):
                result = tslib.array_to_datetime(arg, errors=errors,
                                                 utc=utc, dayfirst=dayfirst,
                                                 yearfirst=yearfirst, freq=freq,
                                                 unit=unit,
                                                 require_iso8601=require_iso8601)

            if com.is_datetime64_dtype(result) and box:
                result = DatetimeIndex(result, tz='utc' if utc else None, name=name)
            return result

        except ValueError as e:
            try:
                values, tz = tslib.datetime_to_datetime64(arg)
                return DatetimeIndex._simple_new(values, name=name, tz=tz)
            except (ValueError, TypeError):
                raise e

    if arg is None:
        return arg
    elif isinstance(arg, tslib.Timestamp):
        return arg
    elif isinstance(arg, Series):
        values = _convert_listlike(arg._values, False, format)
        return Series(values, index=arg.index, name=arg.name)
    elif isinstance(arg, ABCIndexClass):
        return _convert_listlike(arg, box, format, name=arg.name)
    elif com.is_list_like(arg):
        return _convert_listlike(arg, box, format)

    return _convert_listlike(np.array([ arg ]), box, format)[0]


def _attempt_YYYYMMDD(arg, errors):
    """ try to parse the YYYYMMDD/%Y%m%d format, try to deal with NaT-like,
        arg is a passed in as an object dtype, but could really be ints/strings with nan-like/or floats (e.g. with nan)

        Parameters
        ----------
        arg : passed value
        errors : 'raise','ignore','coerce'
        """

    def calc(carg):
        # calculate the actual result
        carg = carg.astype(object)
        return tslib.array_to_datetime(lib.try_parse_year_month_day(carg/10000,carg/100 % 100, carg % 100), errors=errors)

    def calc_with_mask(carg,mask):
        result = np.empty(carg.shape, dtype='M8[ns]')
        iresult = result.view('i8')
        iresult[~mask] = tslib.iNaT
        result[mask] = calc(carg[mask].astype(np.float64).astype(np.int64)).astype('M8[ns]')
        return result

    # try intlike / strings that are ints
    try:
        return calc(arg.astype(np.int64))
    except:
        pass

    # a float with actual np.nan
    try:
        carg = arg.astype(np.float64)
        return calc_with_mask(carg,com.notnull(carg))
    except:
        pass

    # string with NaN-like
    try:
        mask = ~lib.ismember(arg, tslib._nat_strings)
        return calc_with_mask(arg,mask)
    except:
        pass

    return None


def parse_time_string(arg, freq=None, dayfirst=None, yearfirst=None):
    """
    Try hard to parse datetime string, leveraging dateutil plus some extra
    goodies like quarter recognition.

    Parameters
    ----------
    arg : compat.string_types
    freq : str or DateOffset, default None
        Helps with interpreting time string if supplied
    dayfirst : bool, default None
        If None uses default from print_config
    yearfirst : bool, default None
        If None uses default from print_config

    Returns
    -------
    datetime, datetime/dateutil.parser._result, str
    """
    from pandas.core.config import get_option
    if not isinstance(arg, compat.string_types):
        return arg

    from pandas.tseries.offsets import DateOffset
    if isinstance(freq, DateOffset):
        freq = freq.rule_code

    if dayfirst is None:
        dayfirst = get_option("display.date_dayfirst")
    if yearfirst is None:
        yearfirst = get_option("display.date_yearfirst")

    return tslib.parse_datetime_string_with_reso(arg, freq=freq, dayfirst=dayfirst,
                                                 yearfirst=yearfirst)


DateParseError = tslib.DateParseError
normalize_date = tslib.normalize_date


def format(dt):
    """Returns date in YYYYMMDD format."""
    return dt.strftime('%Y%m%d')

OLE_TIME_ZERO = datetime(1899, 12, 30, 0, 0, 0)


def ole2datetime(oledt):
    """function for converting excel date to normal date format"""
    val = float(oledt)

    # Excel has a bug where it thinks the date 2/29/1900 exists
    # we just reject any date before 3/1/1900.
    if val < 61:
        raise ValueError("Value is outside of acceptable range: %s " % val)

    return OLE_TIME_ZERO + timedelta(days=val)
