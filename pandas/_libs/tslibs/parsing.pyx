# -*- coding: utf-8 -*-
"""
Parsing functions for datetime and datetime-like strings.
"""
import sys
import re

import cython
from cython import Py_ssize_t


from cpython.datetime cimport datetime
import time

import numpy as np

# Avoid import from outside _libs
if sys.version_info.major == 2:
    from StringIO import StringIO
else:
    from io import StringIO


# dateutil compat
from dateutil.tz import (tzoffset,
                         tzlocal as _dateutil_tzlocal,
                         tzutc as _dateutil_tzutc,
                         tzstr as _dateutil_tzstr)
from dateutil.relativedelta import relativedelta
from dateutil.parser import DEFAULTPARSER
from dateutil.parser import parse as du_parse

from ccalendar import MONTH_NUMBERS
from nattype import nat_strings, NaT

# ----------------------------------------------------------------------
# Constants


class DateParseError(ValueError):
    pass


_DEFAULT_DATETIME = datetime(1, 1, 1).replace(hour=0, minute=0,
                                              second=0, microsecond=0)

cdef object _TIMEPAT = re.compile(r'^([01]?[0-9]|2[0-3]):([0-5][0-9])')

cdef set _not_datelike_strings = {'a', 'A', 'm', 'M', 'p', 'P', 't', 'T'}

# ----------------------------------------------------------------------


def parse_datetime_string(date_string, freq=None, dayfirst=False,
                          yearfirst=False, **kwargs):
    """parse datetime string, only returns datetime.
    Also cares special handling matching time patterns.

    Returns
    -------
    datetime
    """

    cdef:
        object dt

    if not _does_string_look_like_datetime(date_string):
        raise ValueError('Given date string not likely a datetime.')

    if _TIMEPAT.match(date_string):
        # use current datetime as default, not pass _DEFAULT_DATETIME
        dt = du_parse(date_string, dayfirst=dayfirst,
                      yearfirst=yearfirst, **kwargs)
        return dt

    try:
        dt, _, _ = _parse_dateabbr_string(date_string, _DEFAULT_DATETIME, freq)
        return dt
    except DateParseError:
        raise
    except ValueError:
        pass

    try:
        dt = du_parse(date_string, default=_DEFAULT_DATETIME,
                      dayfirst=dayfirst, yearfirst=yearfirst, **kwargs)
    except TypeError:
        # following may be raised from dateutil
        # TypeError: 'NoneType' object is not iterable
        raise ValueError('Given date string not likely a datetime.')

    return dt


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
    if not isinstance(arg, (str, unicode)):
        # Note: cython recognizes `unicode` in both py2/py3, optimizes
        # this check into a C call.
        return arg

    if getattr(freq, "_typ", None) == "dateoffset":
        freq = freq.rule_code

    if dayfirst is None:
        from pandas.core.config import get_option
        dayfirst = get_option("display.date_dayfirst")
    if yearfirst is None:
        from pandas.core.config import get_option
        yearfirst = get_option("display.date_yearfirst")

    res = parse_datetime_string_with_reso(arg, freq=freq,
                                          dayfirst=dayfirst,
                                          yearfirst=yearfirst)
    return res


cdef parse_datetime_string_with_reso(date_string, freq=None, dayfirst=False,
                                     yearfirst=False):
    """parse datetime string, only returns datetime

    Returns
    -------
    parsed : datetime
    parsed2 : datetime/dateutil.parser._result
    reso : str
        inferred resolution

    Raises
    ------
    ValueError : preliminary check suggests string is not datetime
    DateParseError : error within dateutil
    """
    cdef:
        object parsed, reso

    if not _does_string_look_like_datetime(date_string):
        raise ValueError('Given date string not likely a datetime.')

    try:
        return _parse_dateabbr_string(date_string, _DEFAULT_DATETIME, freq)
    except DateParseError:
        raise
    except ValueError:
        pass

    try:
        parsed, reso = dateutil_parse(date_string, _DEFAULT_DATETIME,
                                      dayfirst=dayfirst, yearfirst=yearfirst,
                                      ignoretz=False, tzinfos=None)
    except Exception as e:
        # TODO: allow raise of errors within instead
        raise DateParseError(e)
    if parsed is None:
        raise DateParseError("Could not parse {dstr}".format(dstr=date_string))
    return parsed, parsed, reso


cpdef bint _does_string_look_like_datetime(object date_string):
    if date_string.startswith('0'):
        # Strings starting with 0 are more consistent with a
        # date-like string than a number
        return True

    try:
        if float(date_string) < 1000:
            return False
    except ValueError:
        pass

    if date_string in _not_datelike_strings:
        return False

    return True


cdef inline object _parse_dateabbr_string(object date_string, object default,
                                          object freq):
    cdef:
        object ret
        int year, quarter = -1, month, mnum, date_len

    # special handling for possibilities eg, 2Q2005, 2Q05, 2005Q1, 05Q1
    assert isinstance(date_string, (str, unicode))

    # len(date_string) == 0
    # should be NaT???

    if date_string in nat_strings:
        return NaT, NaT, ''

    date_string = date_string.upper()
    date_len = len(date_string)

    if date_len == 4:
        # parse year only like 2000
        try:
            ret = default.replace(year=int(date_string))
            return ret, ret, 'year'
        except ValueError:
            pass

    try:
        if 4 <= date_len <= 7:
            i = date_string.index('Q', 1, 6)
            if i == 1:
                quarter = int(date_string[0])
                if date_len == 4 or (date_len == 5
                                     and date_string[i + 1] == '-'):
                    # r'(\d)Q-?(\d\d)')
                    year = 2000 + int(date_string[-2:])
                elif date_len == 6 or (date_len == 7
                                       and date_string[i + 1] == '-'):
                    # r'(\d)Q-?(\d\d\d\d)')
                    year = int(date_string[-4:])
                else:
                    raise ValueError
            elif i == 2 or i == 3:
                # r'(\d\d)-?Q(\d)'
                if date_len == 4 or (date_len == 5
                                     and date_string[i - 1] == '-'):
                    quarter = int(date_string[-1])
                    year = 2000 + int(date_string[:2])
                else:
                    raise ValueError
            elif i == 4 or i == 5:
                if date_len == 6 or (date_len == 7
                                     and date_string[i - 1] == '-'):
                    # r'(\d\d\d\d)-?Q(\d)'
                    quarter = int(date_string[-1])
                    year = int(date_string[:4])
                else:
                    raise ValueError

            if not (1 <= quarter <= 4):
                msg = ('Incorrect quarterly string is given, quarter must be '
                       'between 1 and 4: {dstr}')
                raise DateParseError(msg.format(dstr=date_string))

            if freq is not None:
                # hack attack, #1228
                try:
                    mnum = MONTH_NUMBERS[_get_rule_month(freq)] + 1
                except (KeyError, ValueError):
                    msg = ('Unable to retrieve month information from given '
                           'freq: {freq}'.format(freq=freq))
                    raise DateParseError(msg)

                month = (mnum + (quarter - 1) * 3) % 12 + 1
                if month > mnum:
                    year -= 1
            else:
                month = (quarter - 1) * 3 + 1

            ret = default.replace(year=year, month=month)
            return ret, ret, 'quarter'

    except DateParseError:
        raise
    except ValueError:
        pass

    if date_len == 6 and (freq == 'M' or
                          getattr(freq, 'rule_code', None) == 'M'):
        year = int(date_string[:4])
        month = int(date_string[4:6])
        try:
            ret = default.replace(year=year, month=month)
            return ret, ret, 'month'
        except ValueError:
            pass

    for pat in ['%Y-%m', '%m-%Y', '%b %Y', '%b-%Y']:
        try:
            ret = datetime.strptime(date_string, pat)
            return ret, ret, 'month'
        except ValueError:
            pass

    raise ValueError('Unable to parse {0}'.format(date_string))


cdef dateutil_parse(object timestr, object default, ignoretz=False,
                    tzinfos=None, dayfirst=None, yearfirst=None):
    """ lifted from dateutil to get resolution"""

    cdef:
        object fobj, res, attr, ret, tzdata
        object reso = None
        dict repl = {}

    fobj = StringIO(str(timestr))
    res = DEFAULTPARSER._parse(fobj, dayfirst=dayfirst, yearfirst=yearfirst)

    # dateutil 2.2 compat
    if isinstance(res, tuple):  # PyTuple_Check
        res, _ = res

    if res is None:
        msg = "Unknown datetime string format, unable to parse: {timestr}"
        raise ValueError(msg.format(timestr=timestr))

    for attr in ["year", "month", "day", "hour",
                 "minute", "second", "microsecond"]:
        value = getattr(res, attr)
        if value is not None:
            repl[attr] = value
            reso = attr

    if reso is None:
        msg = "Unable to parse datetime string: {timestr}"
        raise ValueError(msg.format(timestr=timestr))

    if reso == 'microsecond':
        if repl['microsecond'] == 0:
            reso = 'second'
        elif repl['microsecond'] % 1000 == 0:
            reso = 'millisecond'

    ret = default.replace(**repl)
    if res.weekday is not None and not res.day:
        ret = ret + relativedelta.relativedelta(weekday=res.weekday)
    if not ignoretz:
        if callable(tzinfos) or tzinfos and res.tzname in tzinfos:
            if callable(tzinfos):
                tzdata = tzinfos(res.tzname, res.tzoffset)
            else:
                tzdata = tzinfos.get(res.tzname)
            if isinstance(tzdata, datetime.tzinfo):
                tzinfo = tzdata
            elif isinstance(tzdata, (str, unicode)):
                tzinfo = _dateutil_tzstr(tzdata)
            elif isinstance(tzdata, int):
                tzinfo = tzoffset(res.tzname, tzdata)
            else:
                raise ValueError("offset must be tzinfo subclass, "
                                 "tz string, or int offset")
            ret = ret.replace(tzinfo=tzinfo)
        elif res.tzname and res.tzname in time.tzname:
            ret = ret.replace(tzinfo=_dateutil_tzlocal())
        elif res.tzoffset == 0:
            ret = ret.replace(tzinfo=_dateutil_tzutc())
        elif res.tzoffset:
            ret = ret.replace(tzinfo=tzoffset(res.tzname, res.tzoffset))
    return ret, reso


cpdef object _get_rule_month(object source, object default='DEC'):
    """
    Return starting month of given freq, default is December.

    Example
    -------
    >>> _get_rule_month('D')
    'DEC'

    >>> _get_rule_month('A-JAN')
    'JAN'
    """
    if hasattr(source, 'freqstr'):
        source = source.freqstr
    source = source.upper()
    if '-' not in source:
        return default
    else:
        return source.split('-')[1]


# ----------------------------------------------------------------------
# Parsing for type-inference


def try_parse_dates(object[:] values, parser=None,
                    dayfirst=False, default=None):
    cdef:
        Py_ssize_t i, n
        object[:] result

    n = len(values)
    result = np.empty(n, dtype='O')

    if parser is None:
        if default is None:  # GH2618
            date = datetime.now()
            default = datetime(date.year, date.month, 1)

        parse_date = lambda x: du_parse(x, dayfirst=dayfirst, default=default)

        # EAFP here
        try:
            for i in range(n):
                if values[i] == '':
                    result[i] = np.nan
                else:
                    result[i] = parse_date(values[i])
        except Exception:
            # failed
            return values
    else:
        parse_date = parser

        try:
            for i in range(n):
                if values[i] == '':
                    result[i] = np.nan
                else:
                    result[i] = parse_date(values[i])
        except Exception:
            # raise if passed parser and it failed
            raise

    return result.base  # .base to access underlying ndarray


def try_parse_date_and_time(object[:] dates, object[:] times,
                            date_parser=None, time_parser=None,
                            dayfirst=False, default=None):
    cdef:
        Py_ssize_t i, n
        object[:] result

    n = len(dates)
    if len(times) != n:
        raise ValueError('Length of dates and times must be equal')
    result = np.empty(n, dtype='O')

    if date_parser is None:
        if default is None:  # GH2618
            date = datetime.now()
            default = datetime(date.year, date.month, 1)

        parse_date = lambda x: du_parse(x, dayfirst=dayfirst, default=default)

    else:
        parse_date = date_parser

    if time_parser is None:
        parse_time = lambda x: du_parse(x)

    else:
        parse_time = time_parser

    for i in range(n):
        d = parse_date(str(dates[i]))
        t = parse_time(str(times[i]))
        result[i] = datetime(d.year, d.month, d.day,
                             t.hour, t.minute, t.second)

    return result.base  # .base to access underlying ndarray


def try_parse_year_month_day(object[:] years, object[:] months,
                             object[:] days):
    cdef:
        Py_ssize_t i, n
        object[:] result

    n = len(years)
    if len(months) != n or len(days) != n:
        raise ValueError('Length of years/months/days must all be equal')
    result = np.empty(n, dtype='O')

    for i in range(n):
        result[i] = datetime(int(years[i]), int(months[i]), int(days[i]))

    return result.base  # .base to access underlying ndarray


def try_parse_datetime_components(object[:] years,
                                  object[:] months,
                                  object[:] days,
                                  object[:] hours,
                                  object[:] minutes,
                                  object[:] seconds):

    cdef:
        Py_ssize_t i, n
        object[:] result
        int secs
        double float_secs
        double micros

    n = len(years)
    if (len(months) != n or len(days) != n or len(hours) != n or
            len(minutes) != n or len(seconds) != n):
        raise ValueError('Length of all datetime components must be equal')
    result = np.empty(n, dtype='O')

    for i in range(n):
        float_secs = float(seconds[i])
        secs = int(float_secs)

        micros = float_secs - secs
        if micros > 0:
            micros = micros * 1000000

        result[i] = datetime(int(years[i]), int(months[i]), int(days[i]),
                             int(hours[i]), int(minutes[i]), secs,
                             int(micros))

    return result.base  # .base to access underlying ndarray


# ----------------------------------------------------------------------
# Miscellaneous

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


def _format_is_iso(f):
    """
    Does format match the iso8601 set that can be handled by the C parser?
    Generally of form YYYY-MM-DDTHH:MM:SS - date separator can be different
    but must be consistent.  Leading 0s in dates and times are optional.
    """
    iso_template = '%Y{date_sep}%m{date_sep}%d{time_sep}%H:%M:%S.%f'.format
    excluded_formats = ['%Y%m%d', '%Y%m', '%Y']

    for date_sep in [' ', '/', '\\', '-', '.', '']:
        for time_sep in [' ', 'T']:
            if (iso_template(date_sep=date_sep,
                             time_sep=time_sep
                             ).startswith(f) and f not in excluded_formats):
                return True
    return False


def _guess_datetime_format(dt_str, dayfirst=False, dt_str_parse=du_parse,
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
    dt_str_parse : function, defaults to `compat.parse_date` (dateutil)
        This function should take in a datetime string and return
        a `datetime.datetime` guess that the datetime string represents
    dt_str_split : function, defaults to `_DATEUTIL_LEXER_SPLIT` (dateutil)
        This function should take in a datetime string and return
        a list of strings, the guess of the various specific parts
        e.g. '2011/12/30' -> ['2011', '/', '12', '/', '30']

    Returns
    -------
    ret : datetime format string (for `strftime` or `strptime`)
    """
    if dt_str_parse is None or dt_str_split is None:
        return None

    if not isinstance(dt_str, (str, unicode)):
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
    if len({'year', 'month', 'day'} & found_attrs) != 3:
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
    else:
        return None
