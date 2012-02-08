from numpy cimport ndarray
cimport numpy as cnp
cimport cpython

cnp.import_array()
cnp.import_ufunc()

cimport util

cdef extern from "Python.h":
    int PySlice_Check(object)

def get_value_at(ndarray arr, object loc):
    return util.get_value_at(arr, loc)

def set_value_at(ndarray arr, object loc, object val):
    return util.set_value_at(arr, loc, val)

cdef class IndexEngine:

    cpdef get_value(self, ndarray arr, object key):
        '''
        arr : 1-dimensional ndarray
        '''
        cdef:
            Py_ssize_t loc
            void* data_ptr

        loc = self.get_loc(key)
        return util.get_value_at(arr, loc)

    cpdef set_value(self, ndarray arr, object key, object value):
        '''
        arr : 1-dimensional ndarray
        '''
        cdef:
            Py_ssize_t loc
            void* data_ptr

        loc = self.get_loc(key)
        util.set_value_at(arr, loc, value)

cdef class DictIndexEngine(IndexEngine):
    '''
    For accelerating low-level internal details of indexing
    '''

    cdef readonly:
        object index_weakref
        dict mapping
        object mapfun

    cdef:
        bint initialized, integrity

    def __init__(self, index_weakref, object mapfun):
        self.index_weakref = index_weakref
        self.initialized = 0
        self.integrity = 0
        self.mapfun = mapfun

    def __contains__(self, object val):
        self._ensure_initialized()
        return val in self.mapping

    cpdef get_mapping(self, bint check_integrity):
        self._ensure_initialized()
        if check_integrity and self.integrity == 0:
            raise Exception('Index cannot contain duplicate values!')

        return self.mapping

    def clear_mapping(self):
        self.mapping = None
        self.initialized = 0
        self.integrity = 0

    cdef inline _ensure_initialized(self):
        if not self.initialized:
            self.initialize()

    property mapping_prop:

        def __get__(self):
            self._ensure_initialized()
            return self.mapfun

    property has_integrity:

        def __get__(self):
            self._ensure_initialized()
            return self.integrity == 1

    cdef initialize(self):
        values = self.index_weakref().values
        self.mapping = self.mapfun(values)
        if len(self.mapping) == len(values):
            self.integrity = 1
        self.initialized = 1

    cpdef get_loc(self, object val):
        if is_definitely_invalid_key(val):
            raise TypeError

        self._ensure_initialized()
        if not self.integrity:
            raise Exception('Index values are not unique')
        return self.mapping[val]


cdef inline is_definitely_invalid_key(object val):
    return PySlice_Check(val) or cnp.PyArray_Check(val)
