import sys
import importlib


def entry_points_impl():
    # entry_points is available since 3.8 but "horrible inefficient"
    if sys.version_info < (3, 10):
        from importlib_metadata import entry_points
    else:
        from importlib.metadata import entry_points
    return entry_points


class PluginEntry(object):
    """A Plugin entry point

    The entry-point is not loaded/imported on creation.
    Use the method `get()` to import the module and get the attribute.
    """

    class Sentinel(object):
        pass

    # indicate the entry-point object is not loaded yet
    NOT_LOADED = Sentinel()

    def __init__(self, category, name, location):
        """
        :param category str: plugin category name
        :param name str: plugin name (as used by doit)
        :param location str: python object location as <module>:<attr>
        """
        self.obj = self.NOT_LOADED
        self.category = category
        self.name = name
        self.location = location

    def __repr__(self):
        return "PluginEntry('{}', '{}', '{}')".format(
            self.category, self.name, self.location)

    def get(self):
        """return obj, get from cache or load"""
        if self.obj is self.NOT_LOADED:
            self.obj = self.load()
        return self.obj

    def load(self):
        """load/import reference to obj from named module/obj"""
        module_name, obj_name = self.location.split(':')
        try:
            module = importlib.import_module(module_name)
        except ImportError:
            raise Exception('Plugin {} module `{}` not found.'.format(
                self.category, module_name))
        try:
            obj = getattr(module, obj_name)
        except AttributeError:
            raise Exception('Plugin {}:{} module `{}` has no {}.'.format(
                self.category, self.name, module_name, obj_name))
        return obj


class PluginDict(dict):
    """Item values *might* be a PluginEntry or a direct reference to class/obj.

    Values should not be accessed directly, use `get_plugin()`
    to make sure the plugin is loaded.

    Typically, one dict is created for each kind of plugin.
    doit supports 4 categories:
     - COMMAND
     - LOADER
     - BACKEND
     - REPORTER
    """

    entry_point_prefix = 'doit'

    def add_plugins(self, cfg_data, category):
        """read all items from a ConfigParser section containing plugins & entry-points.

        Plugins from entry-point have higher priority
        """
        self.update(self._from_ini(cfg_data, category))
        self.update(self._from_entry_points(category))

    def _from_ini(self, cfg_data, category):
        """plugins from INI file

        INI `section` names map exactly to plugin `category`.
        """
        result = {}
        if category in cfg_data:
            for name, location in cfg_data[category].items():
                result[name] = PluginEntry(category, name, location)
        return result

    def _from_entry_points(self, category):
        """get all plugins from setuptools entry_points"""
        result = {}
        group = f"{self.entry_point_prefix}.{category}"
        entry_points = entry_points_impl()
        for point in entry_points(group=group):
            name = point.name
            location = "{}:{}".format(point.module, point.attr)
            result[name] = PluginEntry(category, name, location)
        return result


    def get_plugin(self, key):
        """load and return a single plugin"""
        val = self[key]
        if isinstance(val, PluginEntry):
            val.name = key  # overwrite obj name attribute
            return val.get()
        else:
            return val

    def to_dict(self):
        """return a standard dict with all plugins values loaded (no PluginEntry)"""
        return {k: self.get_plugin(k) for k in self.keys()}
