# -*- coding: utf-8 -*-
# cython: profile=False
import sys
cdef bint PY3 = (sys.version_info[0] >= 3)

from cpython cimport PyUnicode_Check

from numpy cimport int64_t

cimport util

# ----------------------------------------------------------------------
# Constants

# TODO: Get this from tslibs.nattype once available
_nat_strings = set(['NaT', 'nat', 'NAT', 'nan', 'NaN', 'NAN'])

cdef int64_t NPY_NAT = util.get_nat()

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

    if len(ts) == 0 or ts in _nat_strings:
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
