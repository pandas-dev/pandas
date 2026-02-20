"""
This module contains a Jupyter Server extension that attempts to
make classic server and notebook extensions work in the new server.

Unfortunately, you'll notice that requires some major monkey-patching.
The goal is that this extension will only be used as a temporary
patch to transition extension authors from classic notebook server to jupyter_server.
"""
import os
import types
import inspect
from functools import wraps
from jupyter_core.paths import jupyter_config_path
from traitlets.traitlets import is_trait


from jupyter_server.services.config.manager import ConfigManager
from .traits import NotebookAppTraits


class ClassProxyError(Exception):
    pass


def proxy(obj1, obj2, name, overwrite=False):
    """Redirects a method, property, or trait from object 1 to object 2."""
    if hasattr(obj1, name) and overwrite is False:
        raise ClassProxyError(
            "Cannot proxy the attribute '{name}' from {cls2} because "
            "{cls1} already has this attribute.".format(
                name=name,
                cls1=obj1.__class__,
                cls2=obj2.__class__
            )
        )
    attr = getattr(obj2, name)

    # First check if this thing is a trait (see traitlets)
    cls_attr = getattr(obj2.__class__, name)
    if is_trait(cls_attr) or type(attr) == property:
        thing = property(lambda self: getattr(obj2, name))

    elif isinstance(attr, types.MethodType):
        @wraps(attr)
        def thing(self, *args, **kwargs):
            return attr(*args, **kwargs)

    # Anything else appended on the class is just an attribute of the class.
    else:
        thing = attr

    setattr(obj1.__class__, name, thing)


def public_members(obj):
    members = inspect.getmembers(obj)
    return [m for m, _ in members if not m.startswith('_')]


def diff_members(obj1, obj2):
    """Return all attribute names found in obj2 but not obj1"""
    m1 = public_members(obj1)
    m2 = public_members(obj2)
    return set(m2).difference(m1)


def get_nbserver_extensions(config_dirs):
    cm = ConfigManager(read_config_path=config_dirs)
    section = cm.get("jupyter_notebook_config")
    extensions = section.get('NotebookApp', {}).get('nbserver_extensions', {})
    return extensions


def _link_jupyter_server_extension(serverapp):
    # Get the extension manager from the server
    manager = serverapp.extension_manager
    logger = serverapp.log

    # Hack that patches the enabled extensions list, prioritizing
    # jupyter nbclassic. In the future, it would be much better
    # to incorporate a dependency injection system in the
    # Extension manager that allows extensions to list
    # their dependency tree and sort that way.
    def sorted_extensions(self):
        """Dictionary with extension package names as keys
        and an ExtensionPackage objects as values.
        """
        # Sort the keys and
        keys = sorted(self.extensions.keys())
        keys.remove("notebook_shim")
        keys = ["notebook_shim"] + keys
        return {key: self.extensions[key] for key in keys}

    manager.__class__.sorted_extensions = property(sorted_extensions)

    # Look to see if nbclassic is enabled. if so,
    # link the nbclassic extension here to load
    # its config. Then, port its config to the serverapp
    # for backwards compatibility.
    try:
        pkg = manager.extensions["notebook_shim"]
        pkg.link_point("notebook_shim", serverapp)
        point = pkg.extension_points["notebook_shim"]
        nbapp = point.app
    except Exception:
        nbapp = NotebookAppTraits()

    # Proxy NotebookApp traits through serverapp to notebookapp.
    members = diff_members(serverapp, nbapp)
    for m in members:
        proxy(serverapp, nbapp, m)

    # Find jupyter server extensions listed as notebook server extensions.
    jupyter_paths = jupyter_config_path()
    config_dirs = jupyter_paths + [serverapp.config_dir]
    nbserver_extensions = get_nbserver_extensions(config_dirs)

    # Link all extensions found in the old locations for
    # notebook server extensions.
    for name, enabled in nbserver_extensions.items():
        # If the extension is already enabled in the manager, i.e.
        # because it was discovered already by Jupyter Server
        # through its jupyter_server_config, then don't re-enable here.
        if name not in manager.extensions:
            successful = manager.add_extension(name, enabled=enabled)
            if successful:
                logger.info(
                    "{name} | extension was found and enabled by notebook_shim. "
                    "Consider moving the extension to Jupyter Server's "
                    "extension paths.".format(name=name)
                )
                manager.link_extension(name)

def _load_jupyter_server_extension(serverapp):
    # Patch the config service manager to find the
    # proper path for old notebook frontend extensions
    config_manager = serverapp.config_manager
    read_config_path = config_manager.read_config_path
    read_config_path += [os.path.join(p, 'nbconfig')
                         for p in jupyter_config_path()]
    config_manager.read_config_path = read_config_path
