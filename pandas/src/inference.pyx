def is_bool_array(ndarray[object] values):
    cdef Py_ssize_t i, n = len(values)
    for i in range(n):
        if not util.is_bool_object(values[i]):
            return False
    return True
