cpdef is_array(object o):
    return np.PyArray_Check(o)


def is_bool_array(ndarray[object] values):
    cdef Py_ssize_t i, n = len(values)
    for i in range(n):
        if not util.is_bool_object(values[i]):
            return False
    return True



def isAllDates(ndarray[object, ndim=1] arr):
    cdef int i, size = len(arr)
    cdef object date

    if size == 0:
        return False

    for i from 0 <= i < size:
        date = arr[i]

        if not PyDateTime_Check(date):
            return False

    return True
