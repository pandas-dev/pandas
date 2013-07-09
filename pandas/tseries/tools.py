from datetime import datetime, timedelta
import re
import sys

import numpy as np

import pandas.lib as lib
import pandas.tslib as tslib
import pandas.core.common as com
from pandas.util.py3compat import StringIO

try:
    import dateutil
    from dateutil.parser import parse, DEFAULTPARSER
    from dateutil.relativedelta import relativedelta

    # raise exception if dateutil 2.0 install on 2.x platform
    if (sys.version_info[0] == 2 and
            dateutil.__version__ == '2.0'):  # pragma: no cover
        raise Exception('dateutil 2.0 incompatible with Python 2.x, you must '
                        'install version 1.5 or 2.1+!')
except ImportError:  # pragma: no cover
    print ('Please install python-dateutil via easy_install or some method!')
    raise  # otherwise a 2nd import won't show the message


def _infer_tzinfo(start, end):
    def _infer(a, b):
        tz = a.tzinfo
        if b and b.tzinfo:
            if not (tslib.get_timezone(tz) == tslib.get_timezone(b.tzinfo)):
                raise AssertionError()
        return tz
    tz = None
    if start is not None:
        tz = _infer(start, end)
    elif end is not None:
        tz = _infer(end, start)
    return tz


def _maybe_get_tz(tz):
    if isinstance(tz, basestring):
        import pytz
        tz = pytz.timezone(tz)
    if com.is_integer(tz):
        import pytz
        tz = pytz.FixedOffset(tz / 60)
    return tz


def to_datetime(arg, errors='ignore', dayfirst=False, utc=None, box=True,
                format=None, coerce=False, unit='ns'):
    """
    Convert argument to datetime

    Parameters
    ----------
    arg : string, datetime, array of strings (with possible NAs)
    errors : {'ignore', 'raise'}, default 'ignore'
        Errors are ignored by default (values left untouched)
    dayfirst : boolean, default False
        If True parses dates with the day first, eg 20/01/2005
        Warning: dayfirst=True is not strict, but will prefer to parse
        with day first (this is a known bug).
    utc : boolean, default None
        Return UTC DatetimeIndex if True (converting any tz-aware
        datetime.datetime objects as well)
    box : boolean, default True
        If True returns a DatetimeIndex, if False returns ndarray of values
    format : string, default None
        strftime to parse time, eg "%d/%m/%Y"
    coerce : force errors to NaT (False by default)
    unit : unit of the arg (D,s,ms,us,ns) denote the unit in epoch
        (e.g. a unix timestamp), which is an integer/float number

    Returns
    -------
    ret : datetime if parsing succeeded
    """
    from pandas import Timestamp
    from pandas.core.series import Series
    from pandas.tseries.index import DatetimeIndex

    def _convert_listlike(arg, box):

        if isinstance(arg, (list,tuple)):
            arg = np.array(arg, dtype='O')

        if com.is_datetime64_dtype(arg):
            if box and not isinstance(arg, DatetimeIndex):
                try:
                    return DatetimeIndex(arg, tz='utc' if utc else None)
                except ValueError, e:
                    values, tz = tslib.datetime_to_datetime64(arg)
                    return DatetimeIndex._simple_new(values, None, tz=tz)

            return arg

        arg = com._ensure_object(arg)
        try:
            if format is not None:
                result = tslib.array_strptime(arg, format)
            else:
                result = tslib.array_to_datetime(arg, raise_=errors == 'raise',
                                                 utc=utc, dayfirst=dayfirst,
                                                 coerce=coerce, unit=unit)
            if com.is_datetime64_dtype(result) and box:
                result = DatetimeIndex(result, tz='utc' if utc else None)
            return result

        except ValueError, e:
            try:
                values, tz = tslib.datetime_to_datetime64(arg)
                return DatetimeIndex._simple_new(values, None, tz=tz)
            except (ValueError, TypeError):
                raise e

    if arg is None:
        return arg
    elif isinstance(arg, Timestamp):
        return arg
    elif isinstance(arg, Series):
        values = _convert_listlike(arg.values, box=False)
        return Series(values, index=arg.index, name=arg.name)
    elif com.is_list_like(arg):
        return _convert_listlike(arg, box=box)

    return _convert_listlike(np.array([ arg ]), box=box)[0]

class DateParseError(ValueError):
    pass


# patterns for quarters like '4Q2005', '05Q1'
qpat1full = re.compile(r'(\d)Q(\d\d\d\d)')
qpat2full = re.compile(r'(\d\d\d\d)Q(\d)')
qpat1 = re.compile(r'(\d)Q(\d\d)')
qpat2 = re.compile(r'(\d\d)Q(\d)')
ypat = re.compile(r'(\d\d\d\d)$')
has_time = re.compile('(.+)([\s]|T)+(.+)')


def parse_time_string(arg, freq=None, dayfirst=None, yearfirst=None):
    """
    Try hard to parse datetime string, leveraging dateutil plus some extra
    goodies like quarter recognition.

    Parameters
    ----------
    arg : basestring
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
    from pandas.tseries.offsets import DateOffset
    from pandas.tseries.frequencies import (_get_rule_month, _month_numbers,
                                            _get_freq_str)

    if not isinstance(arg, basestring):
        return arg

    arg = arg.upper()

    default = datetime(1, 1, 1).replace(hour=0, minute=0,
                                        second=0, microsecond=0)

    # special handling for possibilities eg, 2Q2005, 2Q05, 2005Q1, 05Q1
    if len(arg) in [4, 6]:
        m = ypat.match(arg)
        if m:
            ret = default.replace(year=int(m.group(1)))
            return ret, ret, 'year'

        add_century = False
        if len(arg) == 4:
            add_century = True
            qpats = [(qpat1, 1), (qpat2, 0)]
        else:
            qpats = [(qpat1full, 1), (qpat2full, 0)]

        for pat, yfirst in qpats:
            qparse = pat.match(arg)
            if qparse is not None:
                if yfirst:
                    yi, qi = 1, 2
                else:
                    yi, qi = 2, 1
                q = int(qparse.group(yi))
                y_str = qparse.group(qi)
                y = int(y_str)
                if add_century:
                    y += 2000

                if freq is not None:
                    # hack attack, #1228
                    mnum = _month_numbers[_get_rule_month(freq)] + 1
                    month = (mnum + (q - 1) * 3) % 12 + 1
                    if month > mnum:
                        y -= 1
                else:
                    month = (q - 1) * 3 + 1

                ret = default.replace(year=y, month=month)
                return ret, ret, 'quarter'

        is_mo_str = freq is not None and freq == 'M'
        is_mo_off = getattr(freq, 'rule_code', None) == 'M'
        is_monthly = is_mo_str or is_mo_off
        if len(arg) == 6 and is_monthly:
            try:
                ret = _try_parse_monthly(arg)
                if ret is not None:
                    return ret, ret, 'month'
            except Exception:
                pass

    # montly f7u12
    mresult = _attempt_monthly(arg)
    if mresult:
        return mresult

    if dayfirst is None:
        dayfirst = get_option("display.date_dayfirst")
    if yearfirst is None:
        yearfirst = get_option("display.date_yearfirst")

    try:
        parsed, reso = dateutil_parse(arg, default, dayfirst=dayfirst,
                                      yearfirst=yearfirst)
    except Exception, e:
        raise DateParseError(e)

    if parsed is None:
        raise DateParseError("Could not parse %s" % arg)

    return parsed, parsed, reso  # datetime, resolution


def dateutil_parse(timestr, default,
                   ignoretz=False, tzinfos=None,
                   **kwargs):
    """ lifted from dateutil to get resolution"""
    from dateutil import tz
    import time

    res = DEFAULTPARSER._parse(StringIO(timestr), **kwargs)

    if res is None:
        raise ValueError("unknown string format")

    repl = {}
    for attr in ["year", "month", "day", "hour",
                 "minute", "second", "microsecond"]:
        value = getattr(res, attr)
        if value is not None:
            repl[attr] = value
            reso = attr
    if reso == 'microsecond' and repl['microsecond'] == 0:
        reso = 'second'

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
            elif isinstance(tzdata, basestring):
                tzinfo = tz.tzstr(tzdata)
            elif isinstance(tzdata, int):
                tzinfo = tz.tzoffset(res.tzname, tzdata)
            else:
                raise ValueError("offset must be tzinfo subclass, "
                                 "tz string, or int offset")
            ret = ret.replace(tzinfo=tzinfo)
        elif res.tzname and res.tzname in time.tzname:
            ret = ret.replace(tzinfo=tz.tzlocal())
        elif res.tzoffset == 0:
            ret = ret.replace(tzinfo=tz.tzutc())
        elif res.tzoffset:
            ret = ret.replace(tzinfo=tz.tzoffset(res.tzname, res.tzoffset))
    return ret, reso


def _attempt_monthly(val):
    pats = ['%Y-%m', '%m-%Y', '%b %Y', '%b-%Y']
    for pat in pats:
        try:
            ret = datetime.strptime(val, pat)
            return ret, ret, 'month'
        except Exception:
            pass


def _try_parse_monthly(arg):
    base = 2000
    add_base = False
    default = datetime(1, 1, 1).replace(hour=0, minute=0, second=0,
                                        microsecond=0)

    if len(arg) == 4:
        add_base = True
        y = int(arg[:2])
        m = int(arg[2:4])
    elif len(arg) >= 6:  # 201201
        y = int(arg[:4])
        m = int(arg[4:6])
    if add_base:
        y += base
    ret = default.replace(year=y, month=m)
    return ret


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
