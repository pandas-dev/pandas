from numpy cimport ndarray

cdef class IndexEngine:

    pass

cdef class DictIndexEngine(IndexEngine):
    '''
    For accelerating low-level internal details of indexing
    '''

    cdef readonly:
        ndarray values
        dict mapping
        object mapfun

    cdef:
        bint initialized, integrity

    def __init__(self, ndarray values, object mapfun):
        self.values = values
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
        self.mapping = self.mapfun(self.values)
        if len(self.mapping) == len(self.values):
            self.integrity = 1
        self.initialized = 1

    cpdef get_loc(self, object val):
        self._ensure_initialized()
        return self.mapping[val]

