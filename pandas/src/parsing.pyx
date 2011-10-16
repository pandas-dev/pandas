cimport cpython

cdef extern from "math.h":
    double fabs(double)

def maybe_convert_float_list(tuple values):
    cdef:
        Py_ssize_t i, n
        ndarray[float64_t] result
        object val

    n = len(values)
    result = np.empty(n, dtype='f8')

    for i from 0 <= i < n:
        val = values[i]
        result[i] = <float64_t> val

    return val

def maybe_convert_numeric(tuple values, set na_values):
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

def maybe_convert_numeric_list(list values, set na_values):
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

def string_to_ndarray_tuple(tuple values):
    cdef:
        Py_ssize_t i, n
        ndarray[object] result
        object val, onan

    n = len(values)
    result = np.empty(n, dtype=object)
    onan = np.nan

    for i from 0 <= i < n:
        val = values[i]

        if val == '':
            result[i] = onan
        else:
            result[i] = val

    return result

def string_to_ndarray_list(list values):
    cdef:
        Py_ssize_t i, n
        ndarray[object] result
        object val, onan

    n = len(values)
    result = np.empty(n, dtype=object)
    onan = np.nan

    for i from 0 <= i < n:
        val = values[i]

        if val == '':
            result[i] = onan
        else:
            result[i] = val

    return result

def maybe_convert_bool_object(ndarray[object] arr):
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

cdef float64_t FP_ERR = 1e-10

def maybe_convert_int(ndarray[float64_t] arr):
    cdef:
        Py_ssize_t i, n
        ndarray[int64_t] result
        float64_t val

    n = len(arr)
    result = np.empty(n, dtype='i8')
    for i from 0 <= i < n:
        val = arr[i]
        result[i] = <int64_t> val

        # NA
        if val != val:
            return arr

        if fabs(result[i] - val) > FP_ERR:
            return arr

    return result

def maybe_convert_int_object(ndarray[object] arr):
    cdef:
        Py_ssize_t i, n
        ndarray[int64_t] result
        object val

    n = len(arr)
    result = np.empty(n, dtype='i8')
    for i from 0 <= i < n:
        val = arr[i]
        result[i] = <int64_t> val

        # NA
        if val != val:
            return arr

        if fabs(result[i] - val) > FP_ERR:
            return arr

    return result

