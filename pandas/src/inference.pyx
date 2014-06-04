cimport util
from tslib import NaT
from datetime import datetime, timedelta
iNaT = util.get_nat()

# core.common import for fast inference checks
def is_float(object obj):
    return util.is_float_object(obj)

def is_integer(object obj):
    return util.is_integer_object(obj)

def is_bool(object obj):
    return util.is_bool_object(obj)

def is_complex(object obj):
    return util.is_complex_object(obj)

_TYPE_MAP = {
    'int8': 'integer',
    'int16': 'integer',
    'int32': 'integer',
    'int64': 'integer',
    'i' : 'integer',
    'uint8': 'integer',
    'uint16': 'integer',
    'uint32': 'integer',
    'uint64': 'integer',
    'u' : 'integer',
    'float32': 'floating',
    'float64': 'floating',
    'f' : 'floating',
    'complex128': 'complex',
    'c' : 'complex',
    'string': 'string',
    'S' : 'string',
    'unicode': 'unicode',
    'U' : 'unicode',
    'bool': 'boolean',
    'b' : 'boolean',
    'datetime64[ns]' : 'datetime64',
    'M' : 'datetime64',
    'timedelta64[ns]' : 'timedelta64',
    'm' : 'timedelta64',
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

def infer_dtype(object _values):
    cdef:
        Py_ssize_t i, n
        object val
        ndarray values

    if isinstance(_values, np.ndarray):
        values = _values
    elif hasattr(_values,'values'):
        values = _values.values
    else:
        if not isinstance(_values, list):
            _values = list(_values)
        values = list_to_object_array(_values)

    values = getattr(values, 'values', values)

    val_name = values.dtype.name
    if val_name in _TYPE_MAP:
        return _TYPE_MAP[val_name]
    val_kind = values.dtype.kind
    if val_kind in _TYPE_MAP:
        return _TYPE_MAP[val_kind]

    if values.dtype != np.object_:
        values = values.astype('O')

    n = len(values)
    if n == 0:
        return 'empty'

    # make contiguous
    values = values.ravel()

    # try to use a valid value
    for i in range(n):
       val = util.get_value_1d(values, i)
       if not is_null_datetimelike(val):
           break

    if util.is_datetime64_object(val) or val is NaT:
        if is_datetime64_array(values):
            return 'datetime64'
        elif is_timedelta_or_timedelta64_array(values):
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

    elif is_timedelta(val):
        if is_timedelta_or_timedelta64_array(values):
            return 'timedelta'

    elif is_period(val):
        if is_period_array(values):
            return 'period'

    for i in range(n):
        val = util.get_value_1d(values, i)
        if util.is_integer_object(val):
            return 'mixed-integer'

    return 'mixed'

def infer_dtype_list(list values):
    cdef:
        Py_ssize_t i, n = len(values)
    pass


cdef inline bint is_null_datetimelike(v):
    # determine if we have a null for a timedelta/datetime (or integer versions)x
    if util._checknull(v):
        return True
    elif util.is_timedelta64_object(v):
        return v.view('int64') == iNaT
    elif util.is_datetime64_object(v):
        return v.view('int64') == iNaT
    elif util.is_integer_object(v):
        return v == iNaT
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
        object obj

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
        object obj

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
        object obj

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
        object obj

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
        object obj

    if issubclass(values.dtype.type, (np.string_, np.unicode_)):
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
        object obj

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


def is_datetime_array(ndarray[object] values):
    cdef int i, n = len(values)
    cdef object v
    if n == 0:
        return False
    for i in range(n):
        v = values[i]
        if not (is_datetime(v) or is_null_datetimelike(v)):
            return False
    return True


def is_datetime64_array(ndarray values):
    cdef int i, n = len(values)
    cdef object v
    if n == 0:
        return False
    for i in range(n):
        v = values[i]
        if not (util.is_datetime64_object(v) or is_null_datetimelike(v)):
            return False
    return True

def is_timedelta_array(ndarray values):
    cdef int i, n = len(values)
    cdef object v
    if n == 0:
        return False
    for i in range(n):
        v = values[i]
        if not (PyDelta_Check(v) or is_null_datetimelike(v)):
            return False
    return True

def is_timedelta64_array(ndarray values):
    cdef int i, n = len(values)
    cdef object v
    if n == 0:
        return False
    for i in range(n):
        v = values[i]
        if not (util.is_timedelta64_object(v) or is_null_datetimelike(v)):
            return False
    return True

def is_timedelta_or_timedelta64_array(ndarray values):
    """ infer with timedeltas and/or nat/none """
    cdef int i, n = len(values)
    cdef object v
    if n == 0:
        return False
    for i in range(n):
        v = values[i]
        if not (is_timedelta(v) or is_null_datetimelike(v)):
            return False
    return True

def is_date_array(ndarray[object] values):
    cdef int i, n = len(values)
    if n == 0:
        return False
    for i in range(n):
        if not is_date(values[i]):
            return False
    return True

def is_time_array(ndarray[object] values):
    cdef int i, n = len(values)
    if n == 0:
        return False
    for i in range(n):
        if not is_time(values[i]):
            return False
    return True

def is_period(object o):
    from pandas import Period
    return isinstance(o,Period)

def is_period_array(ndarray[object] values):
    cdef int i, n = len(values)
    from pandas.tseries.period import Period

    if n == 0:
        return False
    for i in range(n):
        if not isinstance(values[i], Period):
            return False
    return True


cdef extern from "parse_helper.h":
    inline int floatify(object, double *result) except -1

cdef double fINT64_MAX = <double> INT64_MAX
cdef double fINT64_MIN = <double> INT64_MIN

def maybe_convert_numeric(ndarray[object] values, set na_values,
                          convert_empty=True, coerce_numeric=False):
    '''
    Type inference function-- convert strings to numeric (potentially) and
    convert to proper dtype array
    '''
    cdef:
        int status
        Py_ssize_t i, n
        ndarray[float64_t] floats
        ndarray[complex128_t] complexes
        ndarray[int64_t] ints
        bint seen_float = 0
        bint seen_complex = 0
        object val
        float64_t fval

    n = len(values)

    floats = np.empty(n, dtype='f8')
    complexes = np.empty(n, dtype='c16')
    ints = np.empty(n, dtype='i8')

    for i from 0 <= i < n:
        val = values[i]

        if val in na_values:
            floats[i] = complexes[i] = nan
            seen_float = 1
        elif util.is_float_object(val):
            floats[i] = complexes[i] = val
            seen_float = 1
        elif util.is_integer_object(val):
            floats[i] = ints[i] = val
            seen_int = 1
        elif val is None:
            floats[i] = complexes[i] = nan
            seen_float = 1
        elif hasattr(val,'__len__') and len(val) == 0:
            if convert_empty or coerce_numeric:
                floats[i] = complexes[i] = nan
                seen_float = 1
            else:
                raise ValueError('Empty string encountered')
        elif util.is_complex_object(val):
            complexes[i] = val
            seen_complex = 1
        else:
            try:
                status = floatify(val, &fval)
                floats[i] = fval
                if not seen_float:
                    if '.' in val or fval == INF or fval == NEGINF:
                        seen_float = 1
                    elif 'inf' in val:  # special case to handle +/-inf
                        seen_float = 1
                    elif fval < fINT64_MAX and fval > fINT64_MIN:
                        try:
                            ints[i] = int(val)
                        except ValueError:
                            ints[i] = <int64_t> fval
                    else:
                        seen_float = 1
            except:
                if not coerce_numeric:
                    raise

                floats[i] = nan
                seen_float = 1


    if seen_complex:
        return complexes
    elif seen_float:
        return floats
    else:
        return ints

def maybe_convert_objects(ndarray[object] objects, bint try_float=0,
                          bint safe=0, bint convert_datetime=0, bint convert_timedelta=0):
    '''
    Type inference function-- convert object array to proper dtype
    '''
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
        elif util.is_bool_object(val):
            seen_bool = 1
            bools[i] = val
        elif util.is_float_object(val):
            floats[i] = complexes[i] = val
            seen_float = 1
        elif util.is_datetime64_object(val):
            if convert_datetime:
                idatetimes[i] = convert_to_tsobject(val, None, None).value
                seen_datetime = 1
            else:
                seen_object = 1
                # objects[i] = val.astype('O')
                break
        elif is_timedelta(val):
            if convert_timedelta:
                itimedeltas[i] = convert_to_timedelta64(val, 'ns', False)
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
            if convert_datetime:
                seen_datetime = 1
                idatetimes[i] = convert_to_tsobject(val, None, None).value
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
                elif not seen_datetime and not seen_numeric and not seen_timedelta:
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
                elif not seen_datetime and not seen_numeric and not seen_timedelta:
                    return bools.view(np.bool_)

    return objects


def convert_sql_column(x):
    return maybe_convert_objects(x, try_float=1)

def try_parse_dates(ndarray[object] values, parser=None,
                    dayfirst=False,default=None):
    cdef:
        Py_ssize_t i, n
        ndarray[object] result

    n = len(values)
    result = np.empty(n, dtype='O')

    if parser is None:
        if default is None: # GH2618
           date=datetime.now()
           default=datetime(date.year,date.month,1)

        try:
            from dateutil.parser import parse
            parse_date = lambda x: parse(x, dayfirst=dayfirst,default=default)
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
                            dayfirst=False,default=None):
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
           default=datetime(date.year,date.month,1)

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
    '''
    Substitute for np.vectorize with pandas-friendly dtype inference

    Parameters
    ----------
    arr : ndarray
    f : function

    Returns
    -------
    mapped : ndarray
    '''
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
    '''
    Substitute for np.vectorize with pandas-friendly dtype inference

    Parameters
    ----------
    arr : ndarray
    f : function

    Returns
    -------
    mapped : ndarray
    '''
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


def to_object_array(list rows):
    cdef:
        Py_ssize_t i, j, n, k, tmp
        ndarray[object, ndim=2] result
        list row

    n = len(rows)

    k = 0
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
