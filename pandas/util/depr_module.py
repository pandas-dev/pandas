"""
This module houses utilities mocking deprecated enties.
They are for internal use only and should not be used beyond this purpose.
"""

import sys
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


def _add_proxies(old_mod_name, new_mod_name, entities):
    """ Mock entities moved between modules

    Parameters
    ----------
    old_mod_name : module that used to contain the entity implementation
    new_mod_name : module contains the implementations of the entities
    entities : iterable of the names of the mocked entities

    The mechanics are as follows:

    1. Physically move the entity from 'old_mod_name' to 'new_mod_name'
    2. Add the name of the above entity to the 'entities' iterable
    3. Repeat the (1-2) for each entity you want to move

    Invoking the moved entity from 'old_mod_name' will act as proxy to the
    actual entity in 'new_mod_name'. If warnings are enabled a deprecation
    warning will be issued.
    """

    def create_proxy(entity):

        def wrapper(*args, **kwargs):
            warnings.warn("{old}.{entity} has been deprecated. Use "
                          "{new}.{entity} instead.".format(entity=entity,
                                                           old=old_mod_name,
                                                           new=new_mod_name),
                          DeprecationWarning, stacklevel=2)

            return getattr(new_mod, entity)(*args, **kwargs)

        return wrapper

    importlib.import_module(new_mod_name)
    old_mod = sys.modules[old_mod_name]
    new_mod = sys.modules[new_mod_name]

    for entity in entities:
        setattr(old_mod, entity, create_proxy(entity))
