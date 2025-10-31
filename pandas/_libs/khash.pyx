from cpython.exc cimport (
    PyErr_Fetch,
    PyErr_Occurred,
)
from cpython.object cimport PyObject
from cpython.ref cimport Py_XDECREF


cdef inline raise_if_errors():
    cdef:
        object exc_type
        object exc_value
        PyObject *type
        PyObject *value
        PyObject *traceback

    if PyErr_Occurred():
        PyErr_Fetch(&type, &value, &traceback)
        Py_XDECREF(traceback)
        if value != NULL:
            exc_value = <object>value
            if isinstance(exc_value, str):
                if type != NULL:
                    exc_type = <object>type
                else:
                    exc_type = RuntimeError
                Py_XDECREF(type)
                raise exc_type(exc_value)
            else:
                Py_XDECREF(type)
                raise exc_value


cdef kh_pymap_t* kh_init_pymap_checked():
    cdef kh_pymap_t* table = kh_init_pymap()
    if PyErr_Occurred():
        kh_destroy_pymap(table)
        table = NULL
    raise_if_errors()
    return table


cdef void kh_destroy_pymap_checked(kh_pymap_t* table):
    kh_destroy_pymap(table)
    raise_if_errors()


cdef void kh_clear_pymap_checked(kh_pymap_t* table):
    kh_clear_pymap(table)
    raise_if_errors()


cdef khuint_t kh_get_pymap_checked(kh_pymap_t* table, PyObject* key):
    cdef khuint_t k = kh_get_pymap(table, key)
    raise_if_errors()
    return k


cdef void kh_resize_pymap_checked(kh_pymap_t* table, khuint_t new_n_buckets):
    kh_resize_pymap(table, new_n_buckets)
    raise_if_errors()


cdef khuint_t kh_put_pymap_checked(kh_pymap_t* table, PyObject* key, int* ret):
    cdef khuint_t result = kh_put_pymap(table, key, ret)
    raise_if_errors()
    return result


cdef void kh_del_pymap_checked(kh_pymap_t* table, khuint_t k):
    kh_del_pymap(table, k)
    raise_if_errors()


cdef bint kh_exist_pymap_checked(kh_pymap_t* table, khiter_t k):
    cdef bint res = kh_exist_pymap(table, k)
    raise_if_errors()
    return res
