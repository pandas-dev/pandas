from cpython cimport PyDict_Contains, PyDict_GetItem, PyDict_GetItem


cdef class cache_readonly(object):

    cdef readonly:
        object func, name, allow_setting

    def __init__(self, func=None, allow_setting=False):
        if func is not None:
            self.func = func
            self.name = func.__name__
        self.allow_setting = allow_setting

    def __call__(self, func, doc=None):
        self.func = func
        self.name = func.__name__
        return self

    def __get__(self, obj, typ):
        # Get the cache or set a default one if needed

        cache = getattr(obj, '_cache', None)
        if cache is None:
            try:
                cache = obj._cache = {}
            except (AttributeError):
                return

        if PyDict_Contains(cache, self.name):
            # not necessary to Py_INCREF
            val = <object> PyDict_GetItem(cache, self.name)
        else:
            val = self.func(obj)
            PyDict_SetItem(cache, self.name, val)
        return val

    def __set__(self, obj, value):

        if not self.allow_setting:
            raise Exception("cannot set values for [%s]" % self.name)

        # Get the cache or set a default one if needed
        cache = getattr(obj, '_cache', None)
        if cache is None:
            try:
                cache = obj._cache = {}
            except (AttributeError):
                return

        PyDict_SetItem(cache, self.name, value)

cdef class AxisProperty(object):
    cdef:
        Py_ssize_t axis

    def __init__(self, axis=0):
        self.axis = axis

    def __get__(self, obj, type):
        cdef list axes = obj._data.axes
        return axes[self.axis]

    def __set__(self, obj, value):
        obj._set_axis(self.axis, value)
