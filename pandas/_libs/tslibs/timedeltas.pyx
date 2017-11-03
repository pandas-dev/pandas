# -*- coding: utf-8 -*-
# cython: profile=False
import collections

import sys
cdef bint PY3 = (sys.version_info[0] >= 3)

from cpython cimport PyUnicode_Check

import numpy as np
cimport numpy as np
from numpy cimport int64_t
np.import_array()

from cpython.datetime cimport (datetime, timedelta,
                               PyDelta_Check, PyDateTime_IMPORT)
PyDateTime_IMPORT


cimport util
from util cimport is_timedelta64_object

from nattype import nat_strings

# ----------------------------------------------------------------------
# Constants

cdef int64_t NPY_NAT = util.get_nat()

cdef int64_t DAY_NS = 86400000000000LL

# components named tuple
Components = collections.namedtuple('Components', [
    'days', 'hours', 'minutes', 'seconds',
    'milliseconds', 'microseconds', 'nanoseconds'])

cdef dict timedelta_abbrevs = { 'D': 'd',
                                'd': 'd',
                                'days': 'd',
                                'day': 'd',
                                'hours': 'h',
                                'hour': 'h',
                                'hr': 'h',
                                'h': 'h',
                                'm': 'm',
                                'minute': 'm',
                                'min': 'm',
                                'minutes': 'm',
                                's': 's',
                                'seconds': 's',
                                'sec': 's',
                                'second': 's',
                                'ms': 'ms',
                                'milliseconds': 'ms',
                                'millisecond': 'ms',
                                'milli': 'ms',
                                'millis': 'ms',
                                'us': 'us',
                                'microseconds': 'us',
                                'microsecond': 'us',
                                'micro': 'us',
                                'micros': 'us',
                                'ns': 'ns',
                                'nanoseconds': 'ns',
                                'nano': 'ns',
                                'nanos': 'ns',
                                'nanosecond': 'ns'}

# ----------------------------------------------------------------------


cpdef inline int64_t cast_from_unit(object ts, object unit) except? -1:
    """ return a casting of the unit represented to nanoseconds
        round the fractional part of a float to our precision, p """
    cdef:
        int64_t m
        int p

    if unit == 'D' or unit == 'd':
        m = 1000000000L * 86400
        p = 6
    elif unit == 'h':
        m = 1000000000L * 3600
        p = 6
    elif unit == 'm':
        m = 1000000000L * 60
        p = 6
    elif unit == 's':
        m = 1000000000L
        p = 6
    elif unit == 'ms':
        m = 1000000L
        p = 3
    elif unit == 'us':
        m = 1000L
        p = 0
    elif unit == 'ns' or unit is None:
        m = 1L
        p = 0
    else:
        raise ValueError("cannot cast unit {0}".format(unit))

    # just give me the unit back
    if ts is None:
        return m

    # cast the unit, multiply base/frace separately
    # to avoid precision issues from float -> int
    base = <int64_t> ts
    frac = ts -base
    if p:
        frac = round(frac, p)
    return <int64_t> (base *m) + <int64_t> (frac *m)


cdef inline parse_timedelta_string(object ts):
    """
    Parse a regular format timedelta string. Return an int64_t (in ns)
    or raise a ValueError on an invalid parse.
    """

    cdef:
        unicode c
        bint neg=0, have_dot=0, have_value=0, have_hhmmss=0
        object current_unit=None
        int64_t result=0, m=0, r
        list number=[], frac=[], unit=[]

    # neg : tracks if we have a leading negative for the value
    # have_dot : tracks if we are processing a dot (either post hhmmss or
    #            inside an expression)
    # have_value : track if we have at least 1 leading unit
    # have_hhmmss : tracks if we have a regular format hh:mm:ss

    if len(ts) == 0 or ts in nat_strings:
        return NPY_NAT

    # decode ts if necessary
    if not PyUnicode_Check(ts) and not PY3:
        ts = str(ts).decode('utf-8')

    for c in ts:

        # skip whitespace / commas
        if c == ' ' or c == ',':
            pass

        # positive signs are ignored
        elif c == '+':
            pass

        # neg
        elif c == '-':

            if neg or have_value or have_hhmmss:
                raise ValueError("only leading negative signs are allowed")

            neg = 1

        # number (ascii codes)
        elif ord(c) >= 48 and ord(c) <= 57:

            if have_dot:

                # we found a dot, but now its just a fraction
                if len(unit):
                    number.append(c)
                    have_dot = 0
                else:
                    frac.append(c)

            elif not len(unit):
                number.append(c)

            else:
                r = timedelta_from_spec(number, frac, unit)
                unit, number, frac = [], [c], []

                result += timedelta_as_neg(r, neg)

        # hh:mm:ss.
        elif c == ':':

            # we flip this off if we have a leading value
            if have_value:
                neg = 0

            # we are in the pattern hh:mm:ss pattern
            if len(number):
                if current_unit is None:
                    current_unit = 'h'
                    m = 1000000000L * 3600
                elif current_unit == 'h':
                    current_unit = 'm'
                    m = 1000000000L * 60
                elif current_unit == 'm':
                    current_unit = 's'
                    m = 1000000000L
                r = <int64_t> int(''.join(number)) * m
                result += timedelta_as_neg(r, neg)
                have_hhmmss = 1
            else:
                raise ValueError("expecting hh:mm:ss format, "
                                 "received: {0}".format(ts))

            unit, number = [], []

        # after the decimal point
        elif c == '.':

            if len(number) and current_unit is not None:

                # by definition we had something like
                # so we need to evaluate the final field from a
                # hh:mm:ss (so current_unit is 'm')
                if current_unit != 'm':
                    raise ValueError("expected hh:mm:ss format before .")
                m = 1000000000L
                r = <int64_t> int(''.join(number)) * m
                result += timedelta_as_neg(r, neg)
                have_value = 1
                unit, number, frac = [], [], []

            have_dot = 1

        # unit
        else:
            unit.append(c)
            have_value = 1
            have_dot = 0

    # we had a dot, but we have a fractional
    # value since we have an unit
    if have_dot and len(unit):
        r = timedelta_from_spec(number, frac, unit)
        result += timedelta_as_neg(r, neg)

    # we have a dot as part of a regular format
    # e.g. hh:mm:ss.fffffff
    elif have_dot:

        if ((len(number) or len(frac)) and not len(unit)
                and current_unit is None):
            raise ValueError("no units specified")

        if len(frac) > 0 and len(frac) <= 3:
            m = 10**(3 -len(frac)) * 1000L * 1000L
        elif len(frac) > 3 and len(frac) <= 6:
            m = 10**(6 -len(frac)) * 1000L
        else:
            m = 10**(9 -len(frac))

        r = <int64_t> int(''.join(frac)) * m
        result += timedelta_as_neg(r, neg)

    # we have a regular format
    # we must have seconds at this point (hence the unit is still 'm')
    elif current_unit is not None:
        if current_unit != 'm':
            raise ValueError("expected hh:mm:ss format")
        m = 1000000000L
        r = <int64_t> int(''.join(number)) * m
        result += timedelta_as_neg(r, neg)

    # we have a last abbreviation
    elif len(unit):
        if len(number):
            r = timedelta_from_spec(number, frac, unit)
            result += timedelta_as_neg(r, neg)
        else:
            raise ValueError("unit abbreviation w/o a number")

    # treat as nanoseconds
    # but only if we don't have anything else
    else:
        if have_value:
            raise ValueError("have leftover units")
        if len(number):
            r = timedelta_from_spec(number, frac, 'ns')
            result += timedelta_as_neg(r, neg)

    return result


cdef inline int64_t timedelta_as_neg(int64_t value, bint neg):
    """

    Parameters
    ----------
    value : int64_t of the timedelta value
    neg : boolean if the a negative value
    """
    if neg:
        return -value
    return value


cdef inline timedelta_from_spec(object number, object frac, object unit):
    """

    Parameters
    ----------
    number : a list of number digits
    frac : a list of frac digits
    unit : a list of unit characters
    """
    cdef object n

    try:
        unit = ''.join(unit)
        unit = timedelta_abbrevs[unit.lower()]
    except KeyError:
        raise ValueError("invalid abbreviation: {0}".format(unit))

    n = ''.join(number) + '.' + ''.join(frac)
    return cast_from_unit(float(n), unit)

# ----------------------------------------------------------------------
# Timedelta Construction

# Similar to Timestamp/datetime, this is a construction requirement for
# timedeltas that we need to do object instantiation in python. This will
# serve as a C extension type that shadows the Python class, where we do any
# heavy lifting.
cdef class _Timedelta(timedelta):
    # cdef readonly:
    #     int64_t value     # nanoseconds
    #     object freq       # frequency reference
    #     bint is_populated # are my components populated
    #     int64_t _sign, _d, _h, _m, _s, _ms, _us, _ns

    # higher than np.ndarray and np.matrix
    __array_priority__ = 100

    cpdef bint _has_ns(self):
        return self.value % 1000 != 0

    def _ensure_components(_Timedelta self):
        """
        compute the components
        """
        cdef int64_t sfrac, ifrac, frac, ivalue = self.value

        if self.is_populated:
            return

        # put frac in seconds
        frac = ivalue / (1000 * 1000 * 1000)
        if frac < 0:
            self._sign = -1

            # even fraction
            if (-frac % 86400) != 0:
                self._d = -frac / 86400 + 1
                frac += 86400 * self._d
            else:
                frac = -frac
        else:
            self._sign = 1
            self._d = 0

        if frac >= 86400:
            self._d += frac / 86400
            frac -= self._d * 86400

        if frac >= 3600:
            self._h = frac / 3600
            frac -= self._h * 3600
        else:
            self._h = 0

        if frac >= 60:
            self._m = frac / 60
            frac -= self._m * 60
        else:
            self._m = 0

        if frac >= 0:
            self._s = frac
            frac -= self._s
        else:
            self._s = 0

        sfrac = (self._h * 3600 + self._m * 60
                 + self._s) * (1000 * 1000 * 1000)
        if self._sign < 0:
            ifrac = ivalue + self._d * DAY_NS - sfrac
        else:
            ifrac = ivalue - (self._d * DAY_NS + sfrac)

        if ifrac != 0:
            self._ms = ifrac / (1000 * 1000)
            ifrac -= self._ms * 1000 * 1000
            self._us = ifrac / 1000
            ifrac -= self._us * 1000
            self._ns = ifrac
        else:
            self._ms = 0
            self._us = 0
            self._ns = 0

        self.is_populated = 1

    cpdef timedelta to_pytimedelta(_Timedelta self):
        """
        return an actual datetime.timedelta object
        note: we lose nanosecond resolution if any
        """
        return timedelta(microseconds=int(self.value) / 1000)

    def to_timedelta64(self):
        """ Returns a numpy.timedelta64 object with 'ns' precision """
        return np.timedelta64(self.value, 'ns')

    def total_seconds(self):
        """
        Total duration of timedelta in seconds (to ns precision)
        """
        return 1e-9 * self.value

    def view(self, dtype):
        """ array view compat """
        return np.timedelta64(self.value).view(dtype)

    @property
    def components(self):
        """ Return a Components NamedTuple-like """
        self._ensure_components()
        if self._sign < 0:
            return Components(-self._d, self._h, self._m, self._s,
                              self._ms, self._us, self._ns)

        # return the named tuple
        return Components(self._d, self._h, self._m, self._s,
                          self._ms, self._us, self._ns)

    @property
    def delta(self):
        """ return out delta in ns (for internal compat) """
        return self.value

    @property
    def asm8(self):
        """ return a numpy timedelta64 array view of myself """
        return np.int64(self.value).view('m8[ns]')

    @property
    def resolution(self):
        """ return a string representing the lowest resolution that we have """

        self._ensure_components()
        if self._ns:
            return "N"
        elif self._us:
            return "U"
        elif self._ms:
            return "L"
        elif self._s:
            return "S"
        elif self._m:
            return "T"
        elif self._h:
            return "H"
        else:
            return "D"

    @property
    def days(self):
        """
        Number of Days

        .components will return the shown components
        """
        self._ensure_components()
        if self._sign < 0:
            return -1 * self._d
        return self._d

    @property
    def seconds(self):
        """
        Number of seconds (>= 0 and less than 1 day).

        .components will return the shown components
        """
        self._ensure_components()
        return self._h * 3600 + self._m * 60 + self._s

    @property
    def microseconds(self):
        """
        Number of microseconds (>= 0 and less than 1 second).

        .components will return the shown components
        """
        self._ensure_components()
        return self._ms * 1000 + self._us

    @property
    def nanoseconds(self):
        """
        Number of nanoseconds (>= 0 and less than 1 microsecond).

        .components will return the shown components
        """
        self._ensure_components()
        return self._ns

    def _repr_base(self, format=None):
        """

        Parameters
        ----------
        format : None|all|even_day|sub_day|long

        Returns
        -------
        converted : string of a Timedelta

        """
        cdef object sign_pretty, sign2_pretty, seconds_pretty, subs

        self._ensure_components()

        if self._sign < 0:
            sign_pretty = "-"
            sign2_pretty = " +"
        else:
            sign_pretty = ""
            sign2_pretty = " "

        # show everything
        if format == 'all':
            seconds_pretty = "%02d.%03d%03d%03d" % (
                self._s, self._ms, self._us, self._ns)
            return "%s%d days%s%02d:%02d:%s" % (sign_pretty, self._d,
                                                sign2_pretty, self._h,
                                                self._m, seconds_pretty)

        # by default not showing nano
        if self._ms or self._us or self._ns:
            seconds_pretty = "%02d.%03d%03d" % (self._s, self._ms, self._us)
        else:
            seconds_pretty = "%02d" % self._s

        # if we have a partial day
        subs = (self._h or self._m or self._s or
                self._ms or self._us or self._ns)

        if format == 'even_day':
            if not subs:
                return "%s%d days" % (sign_pretty, self._d)

        elif format == 'sub_day':
            if not self._d:

                # degenerate, don't need the extra space
                if self._sign > 0:
                    sign2_pretty = ""
                return "%s%s%02d:%02d:%s" % (sign_pretty, sign2_pretty,
                                             self._h, self._m, seconds_pretty)

        if subs or format=='long':
            return "%s%d days%s%02d:%02d:%s" % (sign_pretty, self._d,
                                                sign2_pretty, self._h,
                                                self._m, seconds_pretty)
        return "%s%d days" % (sign_pretty, self._d)

    def __repr__(self):
        return "Timedelta('{0}')".format(self._repr_base(format='long'))

    def __str__(self):
        return self._repr_base(format='long')

    def isoformat(self):
        """
        Format Timedelta as ISO 8601 Duration like
        `P[n]Y[n]M[n]DT[n]H[n]M[n]S`, where the `[n]`s are replaced by the
        values. See https://en.wikipedia.org/wiki/ISO_8601#Durations

        .. versionadded:: 0.20.0

        Returns
        -------
        formatted : str

        Notes
        -----
        The longest component is days, whose value may be larger than
        365.
        Every component is always included, even if its value is 0.
        Pandas uses nanosecond precision, so up to 9 decimal places may
        be included in the seconds component.
        Trailing 0's are removed from the seconds component after the decimal.
        We do not 0 pad components, so it's `...T5H...`, not `...T05H...`

        Examples
        --------
        >>> td = pd.Timedelta(days=6, minutes=50, seconds=3,
        ...                   milliseconds=10, microseconds=10, nanoseconds=12)
        >>> td.isoformat()
        'P6DT0H50M3.010010012S'
        >>> pd.Timedelta(hours=1, seconds=10).isoformat()
        'P0DT0H0M10S'
        >>> pd.Timedelta(hours=1, seconds=10).isoformat()
        'P0DT0H0M10S'
        >>> pd.Timedelta(days=500.5).isoformat()
        'P500DT12H0MS'

        See Also
        --------
        Timestamp.isoformat
        """
        components = self.components
        seconds = '{}.{:0>3}{:0>3}{:0>3}'.format(components.seconds,
                                                 components.milliseconds,
                                                 components.microseconds,
                                                 components.nanoseconds)
        # Trim unnecessary 0s, 1.000000000 -> 1
        seconds = seconds.rstrip('0').rstrip('.')
        tpl = 'P{td.days}DT{td.hours}H{td.minutes}M{seconds}S'.format(
            td=components, seconds=seconds)
        return tpl
