cimport cpython

cdef extern from "math.h":
    double fabs(double)

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

def maybe_convert_numeric(ndarray[object] values, set na_values):
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

        if cpython.PyFloat_Check(val):
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
            fval = float(val)
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

def sanitize_objects(ndarray[object] values):
    cdef:
        Py_ssize_t i, n
        object val, onan

    n = len(values)
    onan = np.nan

    for i from 0 <= i < n:
        val = values[i]
        if val == '':
            values[i] = onan

def maybe_convert_bool(ndarray[object] arr):
    cdef:
        Py_ssize_t i, n
        ndarray[uint8_t, cast=True] result
        object val

    n = len(arr)
    result = np.empty(n, dtype=bool)

    for i from 0 <= i < n:
        val = arr[i]

        if val == 'True':
            result[i] = 1
        elif val == 'False':
            result[i] = 0
        else:
            return arr

    return result
