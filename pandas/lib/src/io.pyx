cdef int _EPOCH_ORD = 719163

from datetime import date as pydate

cdef inline int64_t gmtime(object date):
    cdef int y, m, d, h, mn, s, ms, days

    y = PyDateTime_GET_YEAR(date)
    m = PyDateTime_GET_MONTH(date)
    d = PyDateTime_GET_DAY(date)
    h = PyDateTime_DATE_GET_HOUR(date)
    mn = PyDateTime_DATE_GET_MINUTE(date)
    s = PyDateTime_DATE_GET_SECOND(date)
    ms = PyDateTime_DATE_GET_MICROSECOND(date) / 1000

    days = pydate(y, m, 1).toordinal() - _EPOCH_ORD + d - 1
    return ((<int64_t> (((days * 24 + h) * 60 + mn))) * 60 + s) * 1000

cpdef object to_datetime(int64_t timestamp):
    return pydatetime.utcfromtimestamp(timestamp / 1000.0)

cpdef object to_timestamp(object dt):
    return gmtime(dt)

def array_to_timestamp(ndarray[object, ndim=1] arr):
    cdef int i, n
    cdef ndarray[int64_t, ndim=1] result

    n = len(arr)
    result = np.empty(n, dtype=np.int64)

    for i from 0 <= i < n:
        result[i] = gmtime(arr[i])

    return result

def array_to_datetime(ndarray[int64_t, ndim=1] arr):
    cdef int i, n
    cdef ndarray[object, ndim=1] result

    n = len(arr)
    result = np.empty(n, dtype=object)

    for i from 0 <= i < n:
        result[i] = to_datetime(arr[i])

    return result
