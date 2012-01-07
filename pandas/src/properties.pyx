from cpython cimport PyDict_Contains, PyDict_GetItem, PyDict_GetItem

cdef class cache_readonly(object):

    cdef readonly:
        object fget, name

    def __init__(self, func):
        self.fget = func
        self.name = func.__name__

    def __get__(self, obj, type):
        if obj is None:
            return self.fget

        # Get the cache or set a default one if needed

        cache = getattr(obj, '_cache', None)
        if cache is None:
            cache = obj._cache = {}

        if PyDict_Contains(cache, self.name):
            # not necessary to Py_INCREF
            val = <object> PyDict_GetItem(cache, self.name)
            return val
        else:
            val = self.fget(obj)
            PyDict_SetItem(cache, self.name, val)
            return val

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

cdef class SeriesIndex(object):
    cdef:
        object _check_type

    def __init__(self):
        from pandas.core.index import _ensure_index
        self._check_type = _ensure_index

    def __get__(self, obj, type):
        return obj._index

    def __set__(self, obj, value):
        if len(obj) != len(value):
            raise AssertionError('Index length did not match values')
        obj._index = self._check_type(value)

cdef class ValuesProperty(object):

    def __get__(self, obj, type):
        cdef:
            ndarray arr = obj
            object base

        base = np.get_array_base(arr)
        if base is None or not np.PyArray_CheckExact(base):
            arr = arr.view(np.ndarray)
        else:
            arr = base
        return arr
