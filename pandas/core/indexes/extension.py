"""
Shared methods for Index subclasses backed by ExtensionArray.
"""
from pandas.util._decorators import cache_readonly

from .base import Index

from pandas.core.arrays import ExtensionArray


def inherit_from_data(name, delegate, cache=False):
    """
    Make an alias for a method of the underlying ExtensionArray.

    Parameters
    ----------
    name : str
    delegate : class
    cache : bool, default False
        Whether to convert wrapped properties into cache_readonly

    Returns
    -------
    method, property, or cache_readonly
    """
    attr = getattr(delegate, name)

    if isinstance(attr, property):
        # TODO: are we getting the right name/doc here?
        if cache:
            method = cache_readonly(attr.fget)

        else:
            @property
            def method(self):
                return getattr(self._data, name)

            @method.setter
            def method(self, value):
                setattr(self._data, name, value)

    elif not callable(attr):
        # just a normal attribute, no wrapping
        method = attr

    else:
        def method(self, *args, **kwargs):
            result = attr(self._data, *args, **kwargs)
            return result

        method.__name__ = name
        method.__doc__ = attr.__doc__
    return method


def inherit_names(names, delegate, cache=False):
    def wrapper(cls):
        for name in names:
            meth = inherit_from_data(name, delegate, cache=cache)
            setattr(cls, name, meth)

        return cls

    return wrapper
