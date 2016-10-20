import sys
cimport util
from tslib import NaT, get_timezone
from datetime import datetime, timedelta
iNaT = util.get_nat()

cdef bint PY2 = sys.version_info[0] == 2

from util cimport (UINT8_MAX, UINT16_MAX, UINT32_MAX, UINT64_MAX,
                   INT8_MIN, INT8_MAX, INT16_MIN, INT16_MAX,
                   INT32_MAX, INT32_MIN, INT64_MAX, INT64_MIN)

# core.common import for fast inference checks


def is_float(object obj):
    return util.is_float_object(obj)


def is_integer(object obj):
    return util.is_integer_object(obj)


def is_bool(object obj):
    return util.is_bool_object(obj)


def is_complex(object obj):
    return util.is_complex_object(obj)

cpdef bint is_period(object val):
    """ Return a boolean if this is a Period object """
    return util.is_period_object(val)

_TYPE_MAP = {
    'categorical': 'categorical',
    'category': 'categorical',
    'int8': 'integer',
    'int16': 'integer',
    'int32': 'integer',
    'int64': 'integer',
    'i': 'integer',
    'uint8': 'integer',
    'uint16': 'integer',
    'uint32': 'integer',
    'uint64': 'integer',
    'u': 'integer',
    'float32': 'floating',
    'float64': 'floating',
    'f': 'floating',
    'complex128': 'complex',
    'c': 'complex',
    'string': 'string' if PY2 else 'bytes',
    'S': 'string' if PY2 else 'bytes',
    'unicode': 'unicode' if PY2 else 'string',
    'U': 'unicode' if PY2 else 'string',
    'bool': 'boolean',
    'b': 'boolean',
    'datetime64[ns]': 'datetime64',
    'M': 'datetime64',
    'timedelta64[ns]': 'timedelta64',
    'm': 'timedelta64',
}

# types only exist on certain platform
try:
    np.float128
    _TYPE_MAP['float128'] = 'floating'
except AttributeError:
    pass
try:
    np.complex256
    _TYPE_MAP['complex256'] = 'complex'
except AttributeError:
    pass
try:
    np.float16
    _TYPE_MAP['float16'] = 'floating'
except AttributeError:
    pass

cdef _try_infer_map(v):
    """ if its in our map, just return the dtype """
    cdef:
        object attr, val
    for attr in ['name', 'kind', 'base']:
        val = getattr(v.dtype, attr)
        if val in _TYPE_MAP:
            return _TYPE_MAP[val]
    return None


def infer_dtype(object _values):
    """
    we are coercing to an ndarray here
    """

    cdef:
        Py_ssize_t i, n
        object val
        ndarray values
        bint seen_pdnat = False, seen_val = False

    if isinstance(_values, np.ndarray):
        values = _values
    elif hasattr(_values, 'dtype'):

        # this will handle ndarray-like
        # e.g. categoricals
        try:
            values = getattr(_values, '_values', getattr(
                _values, 'values', _values))
        except:
            val = _try_infer_map(_values)
            if val is not None:
                return val

            # its ndarray like but we can't handle
            raise ValueError("cannot infer type for {0}".format(type(_values)))

    else:
        if not isinstance(_values, list):
            _values = list(_values)
        values = list_to_object_array(_values)

    values = getattr(values, 'values', values)
    val = _try_infer_map(values)
    if val is not None:
        return val

    if values.dtype != np.object_:
        values = values.astype('O')

    n = len(values)
    if n == 0:
        return 'empty'

    # make contiguous
    values = values.ravel()

    # try to use a valid value
    for i from 0 <= i < n:
        val = util.get_value_1d(values, i)

        # do not use is_nul_datetimelike to keep
        # np.datetime64('nat') and np.timedelta64('nat')
        if util._checknull(val):
            pass
        elif val is NaT:
            seen_pdnat = True
        else:
            seen_val = True
            break

    # if all values are nan/NaT
    if seen_val is False and seen_pdnat is True:
        return 'datetime'
        # float/object nan is handled in latter logic

    if util.is_datetime64_object(val):
        if is_datetime64_array(values):
            return 'datetime64'
        elif is_timedelta_or_timedelta64_array(values):
            return 'timedelta'

    elif is_timedelta(val):
        if is_timedelta_or_timedelta64_array(values):
            return 'timedelta'

    elif util.is_integer_object(val):
        # a timedelta will show true here as well
        if is_timedelta(val):
            if is_timedelta_or_timedelta64_array(values):
                return 'timedelta'

        if is_integer_array(values):
            return 'integer'
        elif is_integer_float_array(values):
            return 'mixed-integer-float'
        elif is_timedelta_or_timedelta64_array(values):
            return 'timedelta'
        return 'mixed-integer'

    elif is_datetime(val):
        if is_datetime_array(values):
            return 'datetime'

    elif is_date(val):
        if is_date_array(values):
            return 'date'

    elif is_time(val):
        if is_time_array(values):
            return 'time'

    elif util.is_float_object(val):
        if is_float_array(values):
            return 'floating'
        elif is_integer_float_array(values):
            return 'mixed-integer-float'

    elif util.is_bool_object(val):
        if is_bool_array(values):
            return 'boolean'

    elif PyString_Check(val):
        if is_string_array(values):
            return 'string'

    elif PyUnicode_Check(val):
        if is_unicode_array(values):
            return 'unicode'

    elif PyBytes_Check(val):
        if is_bytes_array(values):
            return 'bytes'

    elif is_period(val):
        if is_period_array(values):
            return 'period'

    for i in range(n):
        val = util.get_value_1d(values, i)
        if (util.is_integer_object(val) and
            not util.is_timedelta64_object(val) and
            not util.is_datetime64_object(val)):
            return 'mixed-integer'

    return 'mixed'


def is_possible_datetimelike_array(object arr):
    # determine if we have a possible datetimelike (or null-like) array
    cdef:
        Py_ssize_t i, n = len(arr)
        bint seen_timedelta = 0, seen_datetime = 0
        object v

    for i in range(n):
        v = arr[i]
        if util.is_string_object(v):
            continue
        elif util._checknull(v):
            continue
        elif is_datetime(v):
            seen_datetime=1
        elif is_timedelta(v):
            seen_timedelta=1
        else:
            return False
    return seen_datetime or seen_timedelta


cdef inline bint is_null_datetimelike(v):
    # determine if we have a null for a timedelta/datetime (or integer
    # versions)x
    if util._checknull(v):
        return True
    elif v is NaT:
        return True
    elif util.is_timedelta64_object(v):
        return v.view('int64') == iNaT
    elif util.is_datetime64_object(v):
        return v.view('int64') == iNaT
    elif util.is_integer_object(v):
        return v == iNaT
    return False


cdef inline bint is_null_datetime64(v):
    # determine if we have a null for a datetime (or integer versions),
    # excluding np.timedelta64('nat')
    if util._checknull(v):
        return True
    elif v is NaT:
        return True
    elif util.is_datetime64_object(v):
        return v.view('int64') == iNaT
    return False


cdef inline bint is_null_timedelta64(v):
    # determine if we have a null for a timedelta (or integer versions),
    # excluding np.datetime64('nat')
    if util._checknull(v):
        return True
    elif v is NaT:
        return True
    elif util.is_timedelta64_object(v):
        return v.view('int64') == iNaT
    return False


cdef inline bint is_null_period(v):
    # determine if we have a null for a Period (or integer versions),
    # excluding np.datetime64('nat') and np.timedelta64('nat')
    if util._checknull(v):
        return True
    elif v is NaT:
        return True
    return False


cdef inline bint is_datetime(object o):
    return PyDateTime_Check(o)

cdef inline bint is_date(object o):
    return PyDate_Check(o)

cdef inline bint is_time(object o):
    return PyTime_Check(o)

cdef inline bint is_timedelta(object o):
    return PyDelta_Check(o) or util.is_timedelta64_object(o)


def is_bool_array(ndarray values):
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[object] objbuf

    if issubclass(values.dtype.type, np.bool_):
        return True
    elif values.dtype == np.object_:
        objbuf = values

        if n == 0:
            return False

        for i in range(n):
            if not util.is_bool_object(objbuf[i]):
                return False
        return True
    else:
        return False


def is_integer(object o):
    return util.is_integer_object(o)


def is_integer_array(ndarray values):
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[object] objbuf

    if issubclass(values.dtype.type, np.integer):
        return True
    elif values.dtype == np.object_:
        objbuf = values

        if n == 0:
            return False

        for i in range(n):
            if not util.is_integer_object(objbuf[i]):
                return False
        return True
    else:
        return False


def is_integer_float_array(ndarray values):
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[object] objbuf

    if issubclass(values.dtype.type, np.integer):
        return True
    elif values.dtype == np.object_:
        objbuf = values

        if n == 0:
            return False

        for i in range(n):
            if not (util.is_integer_object(objbuf[i]) or
                    util.is_float_object(objbuf[i])):

                return False
        return True
    else:
        return False


def is_float_array(ndarray values):
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[object] objbuf

    if issubclass(values.dtype.type, np.floating):
        return True
    elif values.dtype == np.object_:
        objbuf = values

        if n == 0:
            return False

        for i in range(n):
            if not util.is_float_object(objbuf[i]):
                return False
        return True
    else:
        return False


def is_string_array(ndarray values):
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[object] objbuf

    if ((PY2 and issubclass(values.dtype.type, np.string_)) or
        not PY2 and issubclass(values.dtype.type, np.unicode_)):
        return True
    elif values.dtype == np.object_:
        objbuf = values

        if n == 0:
            return False

        for i in range(n):
            if not PyString_Check(objbuf[i]):
                return False
        return True
    else:
        return False


def is_unicode_array(ndarray values):
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[object] objbuf

    if issubclass(values.dtype.type, np.unicode_):
        return True
    elif values.dtype == np.object_:
        objbuf = values

        if n == 0:
            return False

        for i in range(n):
            if not PyUnicode_Check(objbuf[i]):
                return False
        return True
    else:
        return False


def is_bytes_array(ndarray values):
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[object] objbuf

    if issubclass(values.dtype.type, np.bytes_):
        return True
    elif values.dtype == np.object_:
        objbuf = values

        if n == 0:
            return False

        for i in range(n):
            if not PyBytes_Check(objbuf[i]):
                return False
        return True
    else:
        return False


def is_datetime_array(ndarray[object] values):
    cdef Py_ssize_t i, null_count = 0, n = len(values)
    cdef object v
    if n == 0:
        return False

    # return False for all nulls
    for i in range(n):
        v = values[i]
        if is_null_datetime64(v):
            # we are a regular null
            if util._checknull(v):
                null_count += 1
        elif not is_datetime(v):
            return False
    return null_count != n


def is_datetime64_array(ndarray values):
    cdef Py_ssize_t i, null_count = 0, n = len(values)
    cdef object v
    if n == 0:
        return False

    # return False for all nulls
    for i in range(n):
        v = values[i]
        if is_null_datetime64(v):
            # we are a regular null
            if util._checknull(v):
                null_count += 1
        elif not util.is_datetime64_object(v):
            return False
    return null_count != n


cpdef is_datetime_with_singletz_array(ndarray[object] values):
    """
    Check values have the same tzinfo attribute.
    Doesn't check values are datetime-like types.
    """

    cdef Py_ssize_t i, j, n = len(values)
    cdef object base_val, base_tz, val, tz

    if n == 0:
        return False

    for i in range(n):
        base_val = values[i]
        if base_val is not NaT:
            base_tz = get_timezone(getattr(base_val, 'tzinfo', None))

            for j in range(i, n):
                val = values[j]
                if val is not NaT:
                    tz = getattr(val, 'tzinfo', None)
                    if base_tz != tz and base_tz != get_timezone(tz):
                        return False
            break

    return True


def is_timedelta_array(ndarray values):
    cdef Py_ssize_t i, null_count = 0, n = len(values)
    cdef object v
    if n == 0:
        return False
    for i in range(n):
        v = values[i]
        if is_null_timedelta64(v):
            # we are a regular null
            if util._checknull(v):
                null_count += 1
        elif not PyDelta_Check(v):
            return False
    return null_count != n


def is_timedelta64_array(ndarray values):
    cdef Py_ssize_t i, null_count = 0, n = len(values)
    cdef object v
    if n == 0:
        return False
    for i in range(n):
        v = values[i]
        if is_null_timedelta64(v):
            # we are a regular null
            if util._checknull(v):
                null_count += 1
        elif not util.is_timedelta64_object(v):
            return False
    return null_count != n


def is_timedelta_or_timedelta64_array(ndarray values):
    """ infer with timedeltas and/or nat/none """
    cdef Py_ssize_t i, null_count = 0, n = len(values)
    cdef object v
    if n == 0:
        return False
    for i in range(n):
        v = values[i]
        if is_null_timedelta64(v):
            # we are a regular null
            if util._checknull(v):
                null_count += 1
        elif not is_timedelta(v):
            return False
    return null_count != n


def is_date_array(ndarray[object] values):
    cdef Py_ssize_t i, n = len(values)
    if n == 0:
        return False
    for i in range(n):
        if not is_date(values[i]):
            return False
    return True


def is_time_array(ndarray[object] values):
    cdef Py_ssize_t i, n = len(values)
    if n == 0:
        return False
    for i in range(n):
        if not is_time(values[i]):
            return False
    return True


def is_period_array(ndarray[object] values):
    cdef Py_ssize_t i, null_count = 0, n = len(values)
    cdef object v
    if n == 0:
        return False

    # return False for all nulls
    for i in range(n):
        v = values[i]
        if is_null_period(v):
            # we are a regular null
            if util._checknull(v):
                null_count += 1
        elif not is_period(v):
            return False
    return null_count != n


cdef extern from "parse_helper.h":
    inline int floatify(object, double *result, int *maybe_int) except -1

cdef int64_t iINT64_MAX = <int64_t> INT64_MAX
cdef int64_t iINT64_MIN = <int64_t> INT64_MIN


def maybe_convert_numeric(object[:] values, set na_values,
                          bint convert_empty=True, bint coerce_numeric=False):
    """
    Type inference function-- convert strings to numeric (potentially) and
    convert to proper dtype array
    """
    cdef:
        int status, maybe_int
        Py_ssize_t i, n = values.size
        ndarray[float64_t] floats = np.empty(n, dtype='f8')
        ndarray[complex128_t] complexes = np.empty(n, dtype='c16')
        ndarray[int64_t] ints = np.empty(n, dtype='i8')
        ndarray[uint8_t] bools = np.empty(n, dtype='u1')
        bint seen_float = False
        bint seen_complex = False
        bint seen_int = False
        bint seen_bool = False
        object val
        float64_t fval

    for i in range(n):
        val = values[i]

        if val.__hash__ is not None and val in na_values:
            floats[i] = complexes[i] = nan
            seen_float = True
        elif util.is_float_object(val):
            floats[i] = complexes[i] = val
            seen_float = True
        elif util.is_integer_object(val):
            floats[i] = ints[i] = val
            seen_int = True
        elif util.is_bool_object(val):
            floats[i] = ints[i] = bools[i] = val
            seen_bool = True
        elif val is None:
            floats[i] = complexes[i] = nan
            seen_float = True
        elif hasattr(val, '__len__') and len(val) == 0:
            if convert_empty or coerce_numeric:
                floats[i] = complexes[i] = nan
                seen_float = True
            else:
                raise ValueError('Empty string encountered')
        elif util.is_complex_object(val):
            complexes[i] = val
            seen_complex = True
        else:
            try:
                status = floatify(val, &fval, &maybe_int)

                if fval in na_values:
                    floats[i] = complexes[i] = nan
                    seen_float = True
                else:
                    floats[i] = fval

                if not seen_float:
                    if maybe_int:
                        as_int = int(val)

                        if as_int <= iINT64_MAX and as_int >= iINT64_MIN:
                            ints[i] = as_int
                        else:
                            raise ValueError('integer out of range')
                    else:
                        seen_float = True
            except (TypeError, ValueError) as e:
                if not coerce_numeric:
                    raise type(e)(str(e) + ' at position {}'.format(i))

                floats[i] = nan
                seen_float = True

    if seen_complex:
        return complexes
    elif seen_float:
        return floats
    elif seen_int:
        return ints
    elif seen_bool:
        return bools.view(np.bool_)
    return ints


def maybe_convert_objects(ndarray[object] objects, bint try_float=0,
                          bint safe=0, bint convert_datetime=0,
                          bint convert_timedelta=0):
    """
    Type inference function-- convert object array to proper dtype
    """
    cdef:
        Py_ssize_t i, n
        ndarray[float64_t] floats
        ndarray[complex128_t] complexes
        ndarray[int64_t] ints
        ndarray[uint8_t] bools
        ndarray[int64_t] idatetimes
        ndarray[int64_t] itimedeltas
        bint seen_float = 0
        bint seen_complex = 0
        bint seen_datetime = 0
        bint seen_datetimetz = 0
        bint seen_timedelta = 0
        bint seen_int = 0
        bint seen_bool = 0
        bint seen_object = 0
        bint seen_null = 0
        bint seen_numeric = 0
        object val, onan
        float64_t fval, fnan

    n = len(objects)

    floats = np.empty(n, dtype='f8')
    complexes = np.empty(n, dtype='c16')
    ints = np.empty(n, dtype='i8')
    bools = np.empty(n, dtype=np.uint8)

    if convert_datetime:
        datetimes = np.empty(n, dtype='M8[ns]')
        idatetimes = datetimes.view(np.int64)

    if convert_timedelta:
        timedeltas = np.empty(n, dtype='m8[ns]')
        itimedeltas = timedeltas.view(np.int64)

    onan = np.nan
    fnan = np.nan

    for i from 0 <= i < n:
        val = objects[i]

        if val is None:
            seen_null = 1
            floats[i] = complexes[i] = fnan
        elif val is NaT:
            if convert_datetime:
                idatetimes[i] = iNaT
                seen_datetime = 1
            if convert_timedelta:
                itimedeltas[i] = iNaT
                seen_timedelta = 1
            if not (convert_datetime or convert_timedelta):
                seen_object = 1
        elif util.is_bool_object(val):
            seen_bool = 1
            bools[i] = val
        elif util.is_float_object(val):
            floats[i] = complexes[i] = val
            seen_float = 1
        elif util.is_datetime64_object(val):
            if convert_datetime:
                idatetimes[i] = convert_to_tsobject(
                    val, None, None, 0, 0).value
                seen_datetime = 1
            else:
                seen_object = 1
                # objects[i] = val.astype('O')
                break
        elif is_timedelta(val):
            if convert_timedelta:
                itimedeltas[i] = convert_to_timedelta64(val, 'ns')
                seen_timedelta = 1
            else:
                seen_object = 1
                break
        elif util.is_integer_object(val):
            seen_int = 1
            floats[i] = <float64_t> val
            complexes[i] = <double complex> val
            if not seen_null:
                try:
                    ints[i] = val
                except OverflowError:
                    seen_object = 1
                    break
        elif util.is_complex_object(val):
            complexes[i] = val
            seen_complex = 1
        elif PyDateTime_Check(val) or util.is_datetime64_object(val):

            # if we have an tz's attached then return the objects
            if convert_datetime:
                if getattr(val, 'tzinfo', None) is not None:
                    seen_datetimetz = 1
                    break
                else:
                    seen_datetime = 1
                    idatetimes[i] = convert_to_tsobject(
                        val, None, None, 0, 0).value
            else:
                seen_object = 1
                break
        elif try_float and not util.is_string_object(val):
            # this will convert Decimal objects
            try:
                floats[i] = float(val)
                complexes[i] = complex(val)
                seen_float = 1
            except Exception:
                seen_object = 1
                break
        else:
            seen_object = 1
            break

    seen_numeric = seen_complex or seen_float or seen_int

    # we try to coerce datetime w/tz but must all have the same tz
    if seen_datetimetz:
        if len(set([ getattr(val, 'tz', None) for val in objects ])) == 1:
            from pandas import DatetimeIndex
            return DatetimeIndex(objects)
        seen_object = 1

    if not seen_object:

        if not safe:
            if seen_null:
                if not seen_bool and not seen_datetime and not seen_timedelta:
                    if seen_complex:
                        return complexes
                    elif seen_float or seen_int:
                        return floats
            else:
                if not seen_bool:
                    if seen_datetime:
                        if not seen_numeric:
                            return datetimes
                    elif seen_timedelta:
                        if not seen_numeric:
                            return timedeltas
                    else:
                        if seen_complex:
                            return complexes
                        elif seen_float:
                            return floats
                        elif seen_int:
                            return ints
                elif (not seen_datetime and not seen_numeric
                      and not seen_timedelta):
                    return bools.view(np.bool_)

        else:
            # don't cast int to float, etc.
            if seen_null:
                if not seen_bool and not seen_datetime and not seen_timedelta:
                    if seen_complex:
                        if not seen_int:
                            return complexes
                    elif seen_float:
                        if not seen_int:
                            return floats
            else:
                if not seen_bool:
                    if seen_datetime:
                        if not seen_numeric:
                            return datetimes
                    elif seen_timedelta:
                        if not seen_numeric:
                            return timedeltas
                    else:
                        if seen_complex:
                            if not seen_int:
                                return complexes
                        elif seen_float:
                            if not seen_int:
                                return floats
                        elif seen_int:
                            return ints
                elif (not seen_datetime and not seen_numeric
                      and not seen_timedelta):
                    return bools.view(np.bool_)

    return objects


def convert_sql_column(x):
    return maybe_convert_objects(x, try_float=1)


def try_parse_dates(ndarray[object] values, parser=None,
                    dayfirst=False, default=None):
    cdef:
        Py_ssize_t i, n
        ndarray[object] result

    n = len(values)
    result = np.empty(n, dtype='O')

    if parser is None:
        if default is None: # GH2618
            date=datetime.now()
            default=datetime(date.year, date.month, 1)

        try:
            from dateutil.parser import parse
            parse_date = lambda x: parse(x, dayfirst=dayfirst, default=default)
        except ImportError: # pragma: no cover
            def parse_date(s):
                try:
                    return datetime.strptime(s, '%m/%d/%Y')
                except Exception:
                    return s
        # EAFP here
        try:
            for i from 0 <= i < n:
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
            for i from 0 <= i < n:
                if values[i] == '':
                    result[i] = np.nan
                else:
                    result[i] = parse_date(values[i])
        except Exception:
            # raise if passed parser and it failed
            raise

    return result


def try_parse_date_and_time(ndarray[object] dates, ndarray[object] times,
                            date_parser=None, time_parser=None,
                            dayfirst=False, default=None):
    cdef:
        Py_ssize_t i, n
        ndarray[object] result

    from datetime import date, time, datetime, timedelta

    n = len(dates)
    if len(times) != n:
        raise ValueError('Length of dates and times must be equal')
    result = np.empty(n, dtype='O')

    if date_parser is None:
        if default is None: # GH2618
            date=datetime.now()
            default=datetime(date.year, date.month, 1)

        try:
            from dateutil.parser import parse
            parse_date = lambda x: parse(x, dayfirst=dayfirst, default=default)
        except ImportError: # pragma: no cover
            def parse_date(s):
                try:
                    return date.strptime(s, '%m/%d/%Y')
                except Exception:
                    return s
    else:
        parse_date = date_parser

    if time_parser is None:
        try:
            from dateutil.parser import parse
            parse_time = lambda x: parse(x)
        except ImportError: # pragma: no cover
            def parse_time(s):
                try:
                    return time.strptime(s, '%H:%M:%S')
                except Exception:
                    return s

    else:
        parse_time = time_parser

    for i from 0 <= i < n:
        d = parse_date(str(dates[i]))
        t = parse_time(str(times[i]))
        result[i] = datetime(d.year, d.month, d.day,
                             t.hour, t.minute, t.second)

    return result


def try_parse_year_month_day(ndarray[object] years, ndarray[object] months,
                             ndarray[object] days):
    cdef:
        Py_ssize_t i, n
        ndarray[object] result

    from datetime import datetime

    n = len(years)
    if len(months) != n or len(days) != n:
        raise ValueError('Length of years/months/days must all be equal')
    result = np.empty(n, dtype='O')

    for i from 0 <= i < n:
        result[i] = datetime(int(years[i]), int(months[i]), int(days[i]))

    return result


def try_parse_datetime_components(ndarray[object] years,
                                  ndarray[object] months,
                                  ndarray[object] days,
                                  ndarray[object] hours,
                                  ndarray[object] minutes,
                                  ndarray[object] seconds):

    cdef:
        Py_ssize_t i, n
        ndarray[object] result
        int secs
        double float_secs
        double micros

    from datetime import datetime

    n = len(years)
    if (len(months) != n or len(days) != n or len(hours) != n or
        len(minutes) != n or len(seconds) != n):
        raise ValueError('Length of all datetime components must be equal')
    result = np.empty(n, dtype='O')

    for i from 0 <= i < n:
        float_secs = float(seconds[i])
        secs = int(float_secs)

        micros = float_secs - secs
        if micros > 0:
            micros = micros * 1000000

        result[i] = datetime(int(years[i]), int(months[i]), int(days[i]),
                             int(hours[i]), int(minutes[i]), secs,
                             int(micros))

    return result


def sanitize_objects(ndarray[object] values, set na_values,
                     convert_empty=True):
    cdef:
        Py_ssize_t i, n
        object val, onan
        Py_ssize_t na_count = 0
        dict memo = {}

    n = len(values)
    onan = np.nan

    for i from 0 <= i < n:
        val = values[i]
        if (convert_empty and val == '') or (val in na_values):
            values[i] = onan
            na_count += 1
        elif val in memo:
            values[i] = memo[val]
        else:
            memo[val] = val

    return na_count


def maybe_convert_bool(ndarray[object] arr,
                       true_values=None, false_values=None):
    cdef:
        Py_ssize_t i, n
        ndarray[uint8_t] result
        object val
        set true_vals, false_vals
        int na_count = 0

    n = len(arr)
    result = np.empty(n, dtype=np.uint8)

    # the defaults
    true_vals = set(('True', 'TRUE', 'true'))
    false_vals = set(('False', 'FALSE', 'false'))

    if true_values is not None:
        true_vals = true_vals | set(true_values)

    if false_values is not None:
        false_vals = false_vals | set(false_values)

    for i from 0 <= i < n:
        val = arr[i]

        if cpython.PyBool_Check(val):
            if val is True:
                result[i] = 1
            else:
                result[i] = 0
        elif val in true_vals:
            result[i] = 1
        elif val in false_vals:
            result[i] = 0
        elif PyFloat_Check(val):
            result[i] = UINT8_MAX
            na_count += 1
        else:
            return arr

    if na_count > 0:
        mask = result == UINT8_MAX
        arr = result.view(np.bool_).astype(object)
        np.putmask(arr, mask, np.nan)
        return arr
    else:
        return result.view(np.bool_)


def map_infer_mask(ndarray arr, object f, ndarray[uint8_t] mask,
                   bint convert=1):
    """
    Substitute for np.vectorize with pandas-friendly dtype inference

    Parameters
    ----------
    arr : ndarray
    f : function

    Returns
    -------
    mapped : ndarray
    """
    cdef:
        Py_ssize_t i, n
        ndarray[object] result
        object val

    n = len(arr)
    result = np.empty(n, dtype=object)
    for i in range(n):
        if mask[i]:
            val = util.get_value_at(arr, i)
        else:
            val = f(util.get_value_at(arr, i))

            # unbox 0-dim arrays, GH #690
            if is_array(val) and PyArray_NDIM(val) == 0:
                # is there a faster way to unbox?
                val = val.item()

        result[i] = val

    if convert:
        return maybe_convert_objects(result,
                                     try_float=0,
                                     convert_datetime=0,
                                     convert_timedelta=0)

    return result


def map_infer(ndarray arr, object f, bint convert=1):
    """
    Substitute for np.vectorize with pandas-friendly dtype inference

    Parameters
    ----------
    arr : ndarray
    f : function

    Returns
    -------
    mapped : ndarray
    """
    cdef:
        Py_ssize_t i, n
        ndarray[object] result
        object val

    n = len(arr)
    result = np.empty(n, dtype=object)
    for i in range(n):
        val = f(util.get_value_at(arr, i))

        # unbox 0-dim arrays, GH #690
        if is_array(val) and PyArray_NDIM(val) == 0:
            # is there a faster way to unbox?
            val = val.item()

        result[i] = val

    if convert:
        return maybe_convert_objects(result,
                                     try_float=0,
                                     convert_datetime=0,
                                     convert_timedelta=0)

    return result


def to_object_array(list rows, int min_width=0):
    """
    Convert a list of lists into an object array.

    Parameters
    ----------
    rows : 2-d array (N, K)
        A list of lists to be converted into an array
    min_width : int
        The minimum width of the object array. If a list
        in `rows` contains fewer than `width` elements,
        the remaining elements in the corresponding row
        will all be `NaN`.

    Returns
    -------
    obj_array : numpy array of the object dtype
    """
    cdef:
        Py_ssize_t i, j, n, k, tmp
        ndarray[object, ndim=2] result
        list row

    n = len(rows)

    k = min_width
    for i from 0 <= i < n:
        tmp = len(rows[i])
        if tmp > k:
            k = tmp

    result = np.empty((n, k), dtype=object)

    for i from 0 <= i < n:
        row = rows[i]

        for j from 0 <= j < len(row):
            result[i, j] = row[j]

    return result


def tuples_to_object_array(ndarray[object] tuples):
    cdef:
        Py_ssize_t i, j, n, k, tmp
        ndarray[object, ndim=2] result
        tuple tup

    n = len(tuples)
    k = len(tuples[0])
    result = np.empty((n, k), dtype=object)
    for i in range(n):
        tup = tuples[i]
        for j in range(k):
            result[i, j] = tup[j]

    return result


def to_object_array_tuples(list rows):
    cdef:
        Py_ssize_t i, j, n, k, tmp
        ndarray[object, ndim=2] result
        tuple row

    n = len(rows)

    k = 0
    for i from 0 <= i < n:
        tmp = len(rows[i])
        if tmp > k:
            k = tmp

    result = np.empty((n, k), dtype=object)

    try:
        for i in range(n):
            row = rows[i]
            for j from 0 <= j < len(row):
                result[i, j] = row[j]
    except Exception:
        # upcast any subclasses to tuple
        for i in range(n):
            row = tuple(rows[i])
            for j from 0 <= j < len(row):
                result[i, j] = row[j]

    return result


def fast_multiget(dict mapping, ndarray keys, default=np.nan):
    cdef:
        Py_ssize_t i, n = len(keys)
        object val
        ndarray[object] output = np.empty(n, dtype='O')

    if n == 0:
        # kludge, for Series
        return np.empty(0, dtype='f8')

    keys = getattr(keys, 'values', keys)

    for i in range(n):
        val = util.get_value_1d(keys, i)
        if val in mapping:
            output[i] = mapping[val]
        else:
            output[i] = default

    return maybe_convert_objects(output)


def downcast_int64(ndarray[int64_t] arr, object na_values,
                   bint use_unsigned=0):
    cdef:
        Py_ssize_t i, n = len(arr)
        int64_t mx = INT64_MIN + 1, mn = INT64_MAX
        int64_t NA = na_values[np.int64]
        int64_t val
        ndarray[uint8_t] mask
        int na_count = 0

    _mask = np.empty(n, dtype=bool)
    mask = _mask.view(np.uint8)

    for i in range(n):
        val = arr[i]

        if val == NA:
            mask[i] = 1
            na_count += 1
            continue

        # not NA
        mask[i] = 0

        if val > mx:
            mx = val

        if val < mn:
            mn = val

    if mn >= 0 and use_unsigned:
        if mx <= UINT8_MAX - 1:
            result = arr.astype(np.uint8)
            if na_count:
                np.putmask(result, _mask, na_values[np.uint8])
            return result

        if mx <= UINT16_MAX - 1:
            result = arr.astype(np.uint16)
            if na_count:
                np.putmask(result, _mask, na_values[np.uint16])
            return result

        if mx <= UINT32_MAX - 1:
            result = arr.astype(np.uint32)
            if na_count:
                np.putmask(result, _mask, na_values[np.uint32])
            return result

    else:
        if mn >= INT8_MIN + 1 and mx <= INT8_MAX:
            result = arr.astype(np.int8)
            if na_count:
                np.putmask(result, _mask, na_values[np.int8])
            return result

        if mn >= INT16_MIN + 1 and mx <= INT16_MAX:
            result = arr.astype(np.int16)
            if na_count:
                np.putmask(result, _mask, na_values[np.int16])
            return result

        if mn >= INT32_MIN + 1 and mx <= INT32_MAX:
            result = arr.astype(np.int32)
            if na_count:
                np.putmask(result, _mask, na_values[np.int32])
            return result

    return arr
