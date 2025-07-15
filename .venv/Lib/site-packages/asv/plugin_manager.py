# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys
import pkgutil
import importlib

from . import commands, plugins
from .console import log


class PluginManager:
    """
    A class to load and manage plugins.

    By default in asv, plugins are searched for in the `asv.plugins`
    namespace package and in the `asv.commands` package.

    Then, any modules specified in the ``plugins`` entry in the
    ``asv.conf.json`` file are loaded.
    """

    def __init__(self):
        self._plugins = []

    def load_plugins(self, package):
        prefix = package.__name__ + "."
        for module_finder, name, ispkg in pkgutil.iter_modules(package.__path__, prefix):
            try:
                mod = importlib.import_module(name)
                self.init_plugin(mod)
                self._plugins.append(mod)
            except ModuleNotFoundError as err:
                if any(keyword in name for keyword in [".mamba", ".virtualenv", ".conda"]):
                    continue  # Fine to not have these
                else:
                    log.error(f"Couldn't load {name} because\n{err}")

    def _load_plugin_by_name(self, name):
        prefix = plugins.__name__ + "."
        for module_finder, module_name, ispkg in pkgutil.iter_modules(plugins.__path__, prefix):
            if name in module_name:
                mod = importlib.import_module(module_name)
                return mod
        return None

    def import_plugin(self, name):
        extended = False
        if name.startswith("."):
            extended = True
            sys.path.insert(0, ".")
            name = name[1:]
        try:
            if extended:
                mod = importlib.import_module(name)
            else:
                mod = self._load_plugin_by_name(name)
            if mod:
                self.init_plugin(mod)
                self._plugins.append(mod)
        finally:
            if extended:
                del sys.path[0]

    def init_plugin(self, mod):
        if hasattr(mod, "setup"):
            mod.setup()

    def run_hook(self, hook_name, args, kwargs):
        for plugin in self._plugins:
            if hasattr(plugin, hook_name):
                getattr(plugin, hook_name)(*args, **kwargs)


plugin_manager = PluginManager()
plugin_manager.load_plugins(commands)
plugin_manager.load_plugins(plugins)

commands.__doc__ = commands._make_docstring()
