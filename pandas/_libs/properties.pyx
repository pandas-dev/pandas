# -*- coding: utf-8 -*-

from cython cimport Py_ssize_t

from cpython cimport (
    PyDict_Contains, PyDict_GetItem, PyDict_SetItem)


cdef class cache_readonly(object):

    cdef readonly:
        object func, name, allow_setting, __doc__

    def __init__(self, func=None, allow_setting=False):
        if func is not None:
            self.func = func
            self.name = func.__name__
            self.__doc__ = func.__doc__
        self.allow_setting = allow_setting

    def __call__(self, func, doc=None):
        self.func = func
        self.name = func.__name__
        if doc is not None:
            self.__doc__ = doc
        else:
            self.__doc__ = func.__doc__
        return self

    def __get__(self, obj, typ):
        if obj is None:
            # accessed on the class, not the instance
            return self

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


cdef class maybe_cache(object):
    """
    A configurable property-like class that acts like a cache_readonly
    if an "_immutable" flag is True and a like property otherwise
    """

    cdef readonly:
        object func, name, __doc__

    def __init__(self, func=None):
        if func is not None:
            self.func = func
            self.name = func.__name__
            self.__doc__ = func.__doc__

    def __call__(self, func, doc=None):
        self.func = func
        self.name = func.__name__
        if doc is not None:
            self.__doc__ = doc
        else:
            self.__doc__ = func.__doc__
        return self

    def __get__(self, obj, typ):
        cdef:
            bint immutable

        if obj is None:
            # accessed on the class, not the instance
            return self

        # Default to non-caching
        immutable = getattr(typ, '_immutable', False)
        if not immutable:
            # behave like a property
            val = self.func(obj)
            return val

        # Get the cache or set a default one if needed
        cache = getattr(obj, '_cache', None)
        if cache is None:
            try:
                cache = obj._cache = {}
            except AttributeError:
                return

        if PyDict_Contains(cache, self.name):
            # not necessary to Py_INCREF
            val = <object> PyDict_GetItem(cache, self.name)
        else:
            val = self.func(obj)
            PyDict_SetItem(cache, self.name, val)
        return val

    def __set__(self, obj, value):
        raise Exception("cannot set values for [%s]" % self.name)


cdef class AxisProperty(object):
    cdef:
        Py_ssize_t axis

    def __init__(self, axis=0):
        self.axis = axis

    def __get__(self, obj, type):
        cdef:
            list axes

        if obj is None:
            # Only instances have _data, not classes
            return None
        else:
            axes = obj._data.axes
        return axes[self.axis]

    def __set__(self, obj, value):
        obj._set_axis(self.axis, value)
