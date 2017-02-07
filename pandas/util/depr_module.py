"""
This module houses a utility class for mocking deprecated modules.
It is for internal use only and should not be used beyond this purpose.
"""

import warnings
import importlib


class _DeprecatedModule(object):
    """ Class for mocking deprecated modules.

    Parameters
    ----------
    deprmod : name of module to be deprecated.
    removals : objects or methods in module that will no longer be
               accessible once module is removed.
    """

    def __init__(self, deprmod, removals=None):
        self.deprmod = deprmod
        self.removals = removals
        if self.removals is not None:
            self.removals = frozenset(self.removals)

        # For introspection purposes.
        self.self_dir = frozenset(dir(self.__class__))

    def __dir__(self):
        deprmodule = self._import_deprmod()
        return dir(deprmodule)

    def __repr__(self):
        deprmodule = self._import_deprmod()
        return repr(deprmodule)

    __str__ = __repr__

    def __getattr__(self, name):
        if name in self.self_dir:
            return object.__getattribute__(self, name)

        deprmodule = self._import_deprmod()
        obj = getattr(deprmodule, name)

        if self.removals is not None and name in self.removals:
            warnings.warn(
                "{deprmod}.{name} is deprecated and will be removed in "
                "a future version.".format(deprmod=self.deprmod, name=name),
                FutureWarning, stacklevel=2)
        else:
            # The object is actually located in another module.
            warnings.warn(
                "{deprmod}.{name} is deprecated. Please use "
                "{modname}.{name} instead.".format(
                    deprmod=self.deprmod, modname=obj.__module__, name=name),
                FutureWarning, stacklevel=2)

        return obj

    def _import_deprmod(self):
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=FutureWarning)
            deprmodule = importlib.import_module(self.deprmod)
            return deprmodule
