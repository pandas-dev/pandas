from numpy cimport ndarray
cimport numpy as cnp
cimport cpython

cnp.import_array()

cimport util

cpdef inline object get_value_at(ndarray arr, object loc):
    cdef:
        Py_ssize_t i
        void* data_ptr
    if util.is_float_object(loc):
        casted = int(loc)
        if casted == loc:
            loc = casted
    i = <Py_ssize_t> loc
    if i < 0:
        i += cnp.PyArray_SIZE(arr)
    data_ptr = cnp.PyArray_GETPTR1(arr, i)
    return cnp.PyArray_GETITEM(arr, data_ptr)

cpdef inline set_value_at(ndarray arr, object loc, object value):
    cdef:
        Py_ssize_t i
    if util.is_float_object(loc):
        casted = int(loc)
        if casted == loc:
            loc = casted
    i = <Py_ssize_t> loc
    if i < 0:
        i += cnp.PyArray_SIZE(arr)

    util.assign_value_1d(arr, i, value)


cdef class IndexEngine:

    cpdef get_value(self, ndarray arr, object key):
        '''
        arr : 1-dimensional ndarray
        '''
        cdef:
            Py_ssize_t loc
            void* data_ptr

        loc = self.get_loc(key)
        return get_value_at(arr, loc)

    cpdef set_value(self, ndarray arr, object key, object value):
        '''
        arr : 1-dimensional ndarray
        '''
        cdef:
            Py_ssize_t loc
            void* data_ptr

        loc = self.get_loc(key)
        set_value_at(arr, loc, value)

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

    def get_mapping(self, bint check_integrity):
        self._ensure_initialized()
        if check_integrity and self.integrity == 0:
            raise Exception('Index cannot contain duplicate values!')

        return self.mapping

    def clear_mapping(self):
        self.mapping = None
        self.initialized = 0
        self.integrity = 0

    cdef _ensure_initialized(self):
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
        self._ensure_initialized()
        return self.mapping[val]


