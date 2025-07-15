"""The extension manager."""

from __future__ import annotations

import importlib
from itertools import starmap

from tornado.gen import multi
from traitlets import Any, Bool, Dict, HasTraits, Instance, List, Unicode, default, observe
from traitlets import validate as validate_trait
from traitlets.config import LoggingConfigurable

from .config import ExtensionConfigManager
from .utils import ExtensionMetadataError, ExtensionModuleNotFound, get_loader, get_metadata


class ExtensionPoint(HasTraits):
    """A simple API for connecting to a Jupyter Server extension
    point defined by metadata and importable from a Python package.
    """

    _linked = Bool(False)
    _app = Any(None, allow_none=True)

    metadata = Dict()

    @validate_trait("metadata")
    def _valid_metadata(self, proposed):
        """Validate metadata."""
        metadata = proposed["value"]
        # Verify that the metadata has a "name" key.
        try:
            self._module_name = metadata["module"]
        except KeyError:
            msg = "There is no 'module' key in the extension's metadata packet."
            raise ExtensionMetadataError(msg) from None

        try:
            self._module = importlib.import_module(self._module_name)
        except ImportError:
            msg = (
                f"The submodule '{self._module_name}' could not be found. Are you "
                "sure the extension is installed?"
            )
            raise ExtensionModuleNotFound(msg) from None
        # If the metadata includes an ExtensionApp, create an instance.
        if "app" in metadata:
            self._app = metadata["app"]()
        return metadata

    @property
    def linked(self):
        """Has this extension point been linked to the server.

        Will pull from ExtensionApp's trait, if this point
        is an instance of ExtensionApp.
        """
        if self.app:
            return self.app._linked
        return self._linked

    @property
    def app(self):
        """If the metadata includes an `app` field"""
        return self._app

    @property
    def config(self):
        """Return any configuration provided by this extension point."""
        if self.app:
            return self.app._jupyter_server_config()
        # At some point, we might want to add logic to load config from
        # disk when extensions don't use ExtensionApp.
        else:
            return {}

    @property
    def module_name(self):
        """Name of the Python package module where the extension's
        _load_jupyter_server_extension can be found.
        """
        return self._module_name

    @property
    def name(self):
        """Name of the extension.

        If it's not provided in the metadata, `name` is set
        to the extensions' module name.
        """
        if self.app:
            return self.app.name
        return self.metadata.get("name", self.module_name)

    @property
    def module(self):
        """The imported module (using importlib.import_module)"""
        return self._module

    def _get_linker(self):
        """Get a linker."""
        if self.app:
            linker = self.app._link_jupyter_server_extension
        else:
            linker = getattr(
                self.module,
                # Search for a _link_jupyter_extension
                "_link_jupyter_server_extension",
                # Otherwise return a dummy function.
                lambda serverapp: None,
            )
        return linker

    def _get_loader(self):
        """Get a loader."""
        loc = self.app
        if not loc:
            loc = self.module
        loader = get_loader(loc)
        return loader

    def _get_starter(self):
        """Get a starter function."""
        if self.app:
            linker = self.app._start_jupyter_server_extension
        else:

            async def _noop_start(serverapp):
                return

            linker = getattr(
                self.module,
                # Search for a _start_jupyter_extension
                "_start_jupyter_server_extension",
                # Otherwise return a no-op function.
                _noop_start,
            )
        return linker

    def validate(self):
        """Check that both a linker and loader exists."""
        try:
            self._get_linker()
            self._get_loader()
        except Exception:
            return False
        else:
            return True

    def link(self, serverapp):
        """Link the extension to a Jupyter ServerApp object.

        This looks for a `_link_jupyter_server_extension` function
        in the extension's module or ExtensionApp class.
        """
        if not self.linked:
            linker = self._get_linker()
            linker(serverapp)
            # Store this extension as already linked.
            self._linked = True

    def load(self, serverapp):
        """Load the extension in a Jupyter ServerApp object.

        This looks for a `_load_jupyter_server_extension` function
        in the extension's module or ExtensionApp class.
        """
        loader = self._get_loader()
        return loader(serverapp)

    async def start(self, serverapp):
        """Call's the extensions 'start' hook where it can
        start (possibly async) tasks _after_ the event loop is running.
        """
        starter = self._get_starter()
        return await starter(serverapp)


class ExtensionPackage(LoggingConfigurable):
    """An API for interfacing with a Jupyter Server extension package.

    Usage:

    ext_name = "my_extensions"
    extpkg = ExtensionPackage(name=ext_name)
    """

    name = Unicode(help="Name of the an importable Python package.")
    enabled = Bool(False, help="Whether the extension package is enabled.")

    _linked_points = Dict()
    extension_points = Dict()
    module = Any(allow_none=True, help="The module for this extension package. None if not enabled")
    metadata = List(Dict(), help="Extension metadata loaded from the extension package.")
    version = Unicode(
        help="""
            The version of this extension package, if it can be found.
            Otherwise, an empty string.
            """,
    )

    @default("version")
    def _load_version(self):
        if not self.enabled:
            return ""
        return getattr(self.module, "__version__", "")

    def __init__(self, **kwargs):
        """Initialize an extension package."""
        super().__init__(**kwargs)
        if self.enabled:
            self._load_metadata()

    def _load_metadata(self):
        """Import package and load metadata

        Only used if extension package is enabled
        """
        name = self.name
        try:
            self.module, self.metadata = get_metadata(name, logger=self.log)
        except ImportError as e:
            msg = (
                f"The module '{name}' could not be found ({e}). Are you "
                "sure the extension is installed?"
            )
            raise ExtensionModuleNotFound(msg) from None
        # Create extension point interfaces for each extension path.
        for m in self.metadata:
            point = ExtensionPoint(metadata=m)
            self.extension_points[point.name] = point
        return name

    def validate(self):
        """Validate all extension points in this package."""
        return all(extension.validate() for extension in self.extension_points.values())

    def link_point(self, point_name, serverapp):
        """Link an extension point."""
        linked = self._linked_points.get(point_name, False)
        if not linked:
            point = self.extension_points[point_name]
            point.link(serverapp)

    def load_point(self, point_name, serverapp):
        """Load an extension point."""
        point = self.extension_points[point_name]
        return point.load(serverapp)

    async def start_point(self, point_name, serverapp):
        """Load an extension point."""
        point = self.extension_points[point_name]
        return await point.start(serverapp)

    def link_all_points(self, serverapp):
        """Link all extension points."""
        for point_name in self.extension_points:
            self.link_point(point_name, serverapp)

    def load_all_points(self, serverapp):
        """Load all extension points."""
        return [self.load_point(point_name, serverapp) for point_name in self.extension_points]

    async def start_all_points(self, serverapp):
        """Load all extension points."""
        for point_name in self.extension_points:
            await self.start_point(point_name, serverapp)


class ExtensionManager(LoggingConfigurable):
    """High level interface for finding, validating,
    linking, loading, and managing Jupyter Server extensions.

    Usage:
    m = ExtensionManager(config_manager=...)
    """

    config_manager = Instance(ExtensionConfigManager, allow_none=True)

    serverapp = Any()  # Use Any to avoid circular import of Instance(ServerApp)

    @default("config_manager")
    def _load_default_config_manager(self):
        config_manager = ExtensionConfigManager()
        self._load_config_manager(config_manager)
        return config_manager

    @observe("config_manager")
    def _config_manager_changed(self, change):
        if change.new:
            self._load_config_manager(change.new)

    # The `extensions` attribute provides a dictionary
    # with extension (package) names mapped to their ExtensionPackage interface
    # (see above). This manager simplifies the interaction between the
    # ServerApp and the extensions being appended.
    extensions = Dict(
        help="""
        Dictionary with extension package names as keys
        and ExtensionPackage objects as values.
        """
    )

    @property
    def sorted_extensions(self):
        """Returns an extensions dictionary, sorted alphabetically."""
        return dict(sorted(self.extensions.items()))

    # The `_linked_extensions` attribute tracks when each extension
    # has been successfully linked to a ServerApp. This helps prevent
    # extensions from being re-linked recursively unintentionally if another
    # extension attempts to link extensions again.
    linked_extensions = Dict(
        help="""
        Dictionary with extension names as keys

        values are True if the extension is linked, False if not.
        """
    )

    @property
    def extension_apps(self):
        """Return mapping of extension names and sets of ExtensionApp objects."""
        return {
            name: {point.app for point in extension.extension_points.values() if point.app}
            for name, extension in self.extensions.items()
        }

    @property
    def extension_points(self):
        """Return mapping of extension point names and ExtensionPoint objects."""
        return {
            name: point
            for value in self.extensions.values()
            for name, point in value.extension_points.items()
        }

    def from_config_manager(self, config_manager):
        """Add extensions found by an ExtensionConfigManager"""
        # load triggered via config_manager trait observer
        self.config_manager = config_manager

    def _load_config_manager(self, config_manager):
        """Actually load our config manager"""
        jpserver_extensions = config_manager.get_jpserver_extensions()
        self.from_jpserver_extensions(jpserver_extensions)

    def from_jpserver_extensions(self, jpserver_extensions):
        """Add extensions from 'jpserver_extensions'-like dictionary."""
        for name, enabled in jpserver_extensions.items():
            self.add_extension(name, enabled=enabled)

    def add_extension(self, extension_name, enabled=False):
        """Try to add extension to manager, return True if successful.
        Otherwise, return False.
        """
        try:
            extpkg = ExtensionPackage(name=extension_name, enabled=enabled)
            self.extensions[extension_name] = extpkg
            return True
        # Raise a warning if the extension cannot be loaded.
        except Exception as e:
            if self.serverapp and self.serverapp.reraise_server_extension_failures:
                raise
            self.log.warning(
                "%s | error adding extension (enabled: %s): %s",
                extension_name,
                enabled,
                e,
                exc_info=True,
            )
        return False

    def link_extension(self, name):
        """Link an extension by name."""
        linked = self.linked_extensions.get(name, False)
        extension = self.extensions[name]
        if not linked and extension.enabled:
            try:
                # Link extension and store links
                extension.link_all_points(self.serverapp)
                self.linked_extensions[name] = True
                self.log.info("%s | extension was successfully linked.", name)
            except Exception as e:
                if self.serverapp and self.serverapp.reraise_server_extension_failures:
                    raise
                self.log.warning("%s | error linking extension: %s", name, e, exc_info=True)

    def load_extension(self, name):
        """Load an extension by name."""
        extension = self.extensions.get(name)

        if extension and extension.enabled:
            try:
                extension.load_all_points(self.serverapp)
            except Exception as e:
                if self.serverapp and self.serverapp.reraise_server_extension_failures:
                    raise
                self.log.warning(
                    "%s | extension failed loading with message: %r", name, e, exc_info=True
                )
            else:
                self.log.info("%s | extension was successfully loaded.", name)

    async def start_extension(self, name):
        """Start an extension by name."""
        extension = self.extensions.get(name)

        if extension and extension.enabled:
            try:
                await extension.start_all_points(self.serverapp)
            except Exception as e:
                if self.serverapp and self.serverapp.reraise_server_extension_failures:
                    raise
                self.log.warning(
                    "%s | extension failed starting with message: %r", name, e, exc_info=True
                )
            else:
                self.log.debug("%s | extension was successfully started.", name)

    async def stop_extension(self, name, apps):
        """Call the shutdown hooks in the specified apps."""
        for app in apps:
            self.log.debug("%s | extension app %r stopping", name, app.name)
            await app.stop_extension()
            self.log.debug("%s | extension app %r stopped", name, app.name)

    def link_all_extensions(self):
        """Link all enabled extensions
        to an instance of ServerApp
        """
        # Sort the extension names to enforce deterministic linking
        # order.
        for name in self.sorted_extensions:
            self.link_extension(name)

    def load_all_extensions(self):
        """Load all enabled extensions and append them to
        the parent ServerApp.
        """
        # Sort the extension names to enforce deterministic loading
        # order.
        for name in self.sorted_extensions:
            self.load_extension(name)

    async def start_all_extensions(self):
        """Start all enabled extensions."""
        # Sort the extension names to enforce deterministic loading
        # order.
        await multi([self.start_extension(name) for name in self.sorted_extensions])

    async def stop_all_extensions(self):
        """Call the shutdown hooks in all extensions."""
        await multi(list(starmap(self.stop_extension, sorted(dict(self.extension_apps).items()))))

    def any_activity(self):
        """Check for any activity currently happening across all extension applications."""
        for _, apps in sorted(dict(self.extension_apps).items()):
            for app in apps:
                if app.current_activity():
                    return True
