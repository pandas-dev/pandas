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
    np.complex64: 'complex',
    np.complex128: 'complex',
    np.string_: 'string',
    np.unicode_: 'unicode',
    np.bool_: 'boolean'
}

try:
    _TYPE_MAP[np.float128] = 'floating'
    _TYPE_MAP[np.complex256] = 'complex'
    _TYPE_MAP[np.float16] = 'floating'
    _TYPE_MAP[np.datetime64] = 'datetime64'
except AttributeError:
    pass

def infer_dtype(object _values):
    cdef:
        Py_ssize_t i, n
        object test_val
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

    test_val = util.get_value_1d(values, 0)

    if util.is_integer_object(test_val):
        if is_integer_array(values):
            return 'integer'

    elif is_datetime(test_val):
        if is_datetime_array(values):
            return 'datetime'

    elif util.is_float_object(test_val):
        if is_float_array(values):

            return 'floating'

    elif util.is_bool_object(test_val):
        if is_bool_array(values):
            return 'boolean'

    elif util.is_string_object(test_val):
        if is_string_array(values):
            return 'string'

    return 'mixed'

def infer_dtype_list(list values):
    cdef:
        Py_ssize_t i, n = len(values)
    pass

cdef inline bint is_datetime(object o):
    return PyDateTime_Check(o)

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
            if not util.is_string_object(objbuf[i]):
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


def maybe_convert_numeric(ndarray[object] values, set na_values):
    '''
    Type inference function-- convert strings to numeric (potentially) and
    convert to proper dtype array
    '''
    cdef:
        Py_ssize_t i, n
        ndarray[float64_t] floats
        ndarray[int64_t] ints
        bint seen_float = 0
        object val
        float64_t fval

    n = len(values)

    floats = np.empty(n, dtype='f8')
    ints = np.empty(n, dtype='i8')

    for i from 0 <= i < n:
        val = values[i]

        if util.is_float_object(val):
            floats[i] = val
            seen_float = 1
        elif val in na_values:
            floats[i] = nan
            seen_float = 1
        elif val is None:
            floats[i] = nan
            seen_float = 1
        elif len(val) == 0:
            floats[i] = nan
            seen_float = 1
        else:
            fval = util.floatify(val)
            floats[i] = fval
            if not seen_float:
                if '.' in val:
                    seen_float = 1
                else:
                    ints[i] = <int64_t> fval

    if seen_float:
        return floats
    else:
        return ints

def maybe_convert_objects(ndarray[object] objects, bint try_float=0):
    '''
    Type inference function-- convert object array to proper dtype
    '''
    cdef:
        Py_ssize_t i, n
        ndarray[float64_t] floats
        ndarray[int64_t] ints
        ndarray[uint8_t] bools
        bint seen_float = 0
        bint seen_int = 0
        bint seen_bool = 0
        bint seen_object = 0
        bint seen_null = 0
        object val, onan
        float64_t fval, fnan

    n = len(objects)

    floats = np.empty(n, dtype='f8')
    ints = np.empty(n, dtype='i8')
    bools = np.empty(n, dtype=np.uint8)

    onan = np.nan
    fnan = np.nan

    for i from 0 <= i < n:
        val = objects[i]

        if val is None:
            seen_null = 1
            floats[i] = fnan
        elif util.is_bool_object(val):
            seen_bool = 1
            bools[i] = val
        elif util.is_integer_object(val):
            seen_int = 1
            floats[i] = <float64_t> val
            if not seen_null:
                ints[i] = val
        elif util.is_float_object(val):
            floats[i] = val
            seen_float = 1
        elif try_float and not util.is_string_object(val):
            # this will convert Decimal objects
            try:
                floats[i] = float(val)
                seen_float = 1
            except Exception:
                seen_object = 1
        else:
            seen_object = 1

    if seen_null:
        if (seen_float or seen_int) and not seen_object:
            return floats
        else:
            return objects
    else:
        if seen_object:
            return objects
        elif not seen_bool:
            if seen_float:
                return floats
            elif seen_int:
                return ints
        else:
            if not seen_float and not seen_int:
                return bools.view(np.bool_)

        return objects

def convert_sql_column(x):
    return maybe_convert_objects(x, try_float=1)

def try_parse_dates(ndarray[object] values, parser=None):
    cdef:
        Py_ssize_t i, n
        ndarray[object] result

    from datetime import datetime

    n = len(values)
    result = np.empty(n, dtype='O')

    if parser is None:
        try:
            from dateutil import parser
            parse_date = parser.parse
        except ImportError: # pragma: no cover
            def parse_date(s):
                try:
                    return datetime.strptime(s, '%m/%d/%Y')
                except Exception:
                    return s
    else:
        parse_date = parser

    # EAFP
    try:
        for i from 0 <= i < n:
            result[i] = parse_date(values[i])
    except Exception:
        # failed
        return values

    return result

def sanitize_objects(ndarray[object] values, set na_values):
    cdef:
        Py_ssize_t i, n
        object val, onan
        Py_ssize_t na_count = 0
        dict memo = {}

    n = len(values)
    onan = np.nan

    for i from 0 <= i < n:
        val = values[i]
        if val == '' or val in na_values:
            values[i] = onan
            na_count += 1
        elif val in memo:
            values[i] = memo[val]
        else:
            memo[val] = val

    return na_count

def maybe_convert_bool(ndarray[object] arr):
    cdef:
        Py_ssize_t i, n
        ndarray[uint8_t] result
        object val

    n = len(arr)
    result = np.empty(n, dtype=np.uint8)

    for i from 0 <= i < n:
        val = arr[i]

        if val == 'True':
            result[i] = 1
        elif val == 'False':
            result[i] = 0
        else:
            return arr

    return result.view(np.bool_)


def map_infer(ndarray arr, object f):
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
        flatiter it
        ndarray[object] result
        object val

    it = <flatiter> PyArray_IterNew(arr)
    n = len(arr)
    result = np.empty(n, dtype=object)
    for i in range(n):
        val = f(PyArray_GETITEM(arr, PyArray_ITER_DATA(it)))

        # unbox 0-dim arrays, GH #690
        if is_array(val) and PyArray_NDIM(val) == 0:
            # is there a faster way to unbox?
            val = val.item()

        result[i] = val


        PyArray_ITER_NEXT(it)

    return maybe_convert_objects(result, try_float=0)

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

