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
    alts : alternative modules to be used to access objects or methods
           available in module.
    removals : objects or methods in module that will no longer be
               accessible once module is removed.
    """
    def __init__(self, deprmod, alts=None, removals=None):
        self.deprmod = deprmod

        self.alts = alts
        if self.alts is not None:
            self.alts = frozenset(self.alts)

        self.removals = removals
        if self.removals is not None:
            self.removals = frozenset(self.removals)

        # For introspection purposes.
        self.self_dir = frozenset(dir(self.__class__))

    def __dir__(self):
        _dir = object.__dir__(self)

        if self.removals is not None:
            _dir.extend(list(self.removals))

        if self.alts is not None:
            for modname in self.alts:
                module = importlib.import_module(modname)
                _dir.extend(dir(module))

        return _dir

    def __getattr__(self, name):
        if name in self.self_dir:
            return object.__getattribute__(self, name)

        if self.removals is not None and name in self.removals:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', category=FutureWarning)
                module = importlib.import_module(self.deprmod)

            warnings.warn(
                "{deprmod}.{name} is deprecated and will be removed in "
                "a future version.".format(deprmod=self.deprmod, name=name),
                FutureWarning, stacklevel=2)

            return object.__getattribute__(module, name)

        if self.alts is not None:
            for modname in self.alts:
                module = importlib.import_module(modname)

                if hasattr(module, name):
                    warnings.warn(
                        "{deprmod}.{name} is deprecated. Please use "
                        "{modname}.{name} instead.".format(
                            deprmod=self.deprmod, modname=modname, name=name),
                        FutureWarning, stacklevel=2)

                    return getattr(module, name)

        raise AttributeError("module '{deprmod}' has no attribute "
                             "'{name}'".format(deprmod=self.deprmod,
                                               name=name))
