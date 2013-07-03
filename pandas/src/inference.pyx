cimport util

_TYPE_MAP = {
    np.int8: 'integer',
    np.int16: 'integer',
    np.int32: 'integer',
    np.int64: 'integer',
    np.uint8: 'integer',
    np.uint16: 'integer',
    np.uint32: 'integer',
    np.uint64: 'integer',
    np.float32: 'floating',
    np.float64: 'floating',
    np.complex128: 'complex',
    np.complex128: 'complex',
    np.string_: 'string',
    np.unicode_: 'unicode',
    np.bool_: 'boolean',
    np.datetime64 : 'datetime64'
}

try:
    _TYPE_MAP[np.float128] = 'floating'
    _TYPE_MAP[np.complex256] = 'complex'
    _TYPE_MAP[np.float16] = 'floating'
    _TYPE_MAP[np.datetime64] = 'datetime64'
    _TYPE_MAP[np.timedelta64] = 'timedelta64'
except AttributeError:
    pass

def infer_dtype(object _values):
    cdef:
        Py_ssize_t i, n
        object val
        ndarray values

    if isinstance(_values, np.ndarray):
        values = _values
    else:
        if not isinstance(_values, list):
            _values = list(_values)
        values = list_to_object_array(_values)

    n = len(values)
    if n == 0:
        return 'empty'

    val_kind = values.dtype.type
    if val_kind in _TYPE_MAP:
        return _TYPE_MAP[val_kind]

    if values.dtype != np.object_:
        values = values.astype('O')

    val = util.get_value_1d(values, 0)

    if util.is_datetime64_object(val):
        if is_datetime64_array(values):
            return 'datetime64'
    elif util.is_integer_object(val):
        if is_integer_array(values):
            return 'integer'
        elif is_integer_float_array(values):
            return 'mixed-integer-float'
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


cdef inline bint is_datetime(object o):
    return PyDateTime_Check(o)

cdef inline bint is_date(object o):
    return PyDate_Check(o)

cdef inline bint is_time(object o):
    return PyTime_Check(o)

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
    if n == 0:
        return False
    for i in range(n):
        if not is_datetime(values[i]):
            return False
    return True


def is_datetime64_array(ndarray values):
    cdef int i, n = len(values)
    if n == 0:
        return False
    for i in range(n):
        if not util.is_datetime64_object(values[i]):
            return False
    return True

def is_timedelta(object o):
    import datetime
    return isinstance(o,datetime.timedelta) or isinstance(o,np.timedelta64)

def is_timedelta_array(ndarray values):
    import datetime
    cdef int i, n = len(values)
    if n == 0:
        return False
    for i in range(n):
        if not isinstance(values[i],datetime.timedelta):
            return False
    return True

def is_timedelta64_array(ndarray values):
    cdef int i, n = len(values)
    if n == 0:
        return False
    for i in range(n):
        if not isinstance(values[i],np.timedelta64):
            return False
    return True

def is_timedelta_or_timedelta64_array(ndarray values):
    import datetime
    cdef int i, n = len(values)
    if n == 0:
        return False
    for i in range(n):
        if not (isinstance(values[i],datetime.timedelta) or isinstance(values[i],np.timedelta64)):
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
    from pandas import Period

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
                          bint safe=0, bint convert_datetime=0):
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
        bint seen_float = 0
        bint seen_complex = 0
        bint seen_datetime = 0
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
    datetimes = np.empty(n, dtype='M8[ns]')
    idatetimes = datetimes.view(np.int64)

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
                if not seen_bool and not seen_datetime:
                    if seen_complex:
                        return complexes
                    elif seen_float or seen_int:
                        return floats
            else:
                if not seen_bool:
                    if seen_datetime:
                        if not seen_numeric:
                            return datetimes
                    else:
                        if seen_complex:
                            return complexes
                        elif seen_float:
                            return floats
                        elif seen_int:
                            return ints
                elif not seen_datetime and not seen_numeric:
                    return bools.view(np.bool_)

        else:
            # don't cast int to float, etc.
            if seen_null:
                if not seen_bool and not seen_datetime:
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
                    else:
                        if seen_complex:
                            if not seen_int:
                                return complexes
                        elif seen_float:
                            if not seen_int:
                                return floats
                        elif seen_int:
                            return ints
                elif not seen_datetime and not seen_numeric:
                    return bools.view(np.bool_)

    return objects


def convert_sql_column(x):
    return maybe_convert_objects(x, try_float=1)

def try_parse_dates(ndarray[object] values, parser=None,
                    dayfirst=False,default=None):
    cdef:
        Py_ssize_t i, n
        ndarray[object] result

    from datetime import datetime, timedelta

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
        double micros

    from datetime import datetime

    n = len(years)
    if (len(months) != n and len(days) != n and len(hours) != n and
        len(minutes) != n and len(seconds) != n):
        raise ValueError('Length of all datetime components must be equal')
    result = np.empty(n, dtype='O')

    for i from 0 <= i < n:
        secs = int(seconds[i])

        micros = seconds[i] - secs
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
        return maybe_convert_objects(result, try_float=0,
                                     convert_datetime=0)

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
        return maybe_convert_objects(result, try_float=0,
                                     convert_datetime=0)

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

    for i in range(n):
        val = util.get_value_1d(keys, i)
        if val in mapping:
            output[i] = mapping[val]
        else:
            output[i] = default

    return maybe_convert_objects(output)
