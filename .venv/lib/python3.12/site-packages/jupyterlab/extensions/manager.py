"""Base classes for the extension manager."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import json
import re
from dataclasses import dataclass, field, fields, replace
from pathlib import Path
from typing import Optional, Union

import tornado
from jupyterlab_server.translation_utils import translator
from traitlets import Enum
from traitlets.config import Configurable, LoggingConfigurable

from jupyterlab.commands import (
    _AppHandler,
    _ensure_options,
    disable_extension,
    enable_extension,
    get_app_info,
)

PYTHON_TO_SEMVER = {"a": "-alpha.", "b": "-beta.", "rc": "-rc."}


def _ensure_compat_errors(info, app_options):
    """Ensure that the app info has compat_errors field"""
    handler = _AppHandler(app_options)
    info["compat_errors"] = handler._get_extension_compat()


_message_map = {
    "install": re.compile(r"(?P<name>.*) needs to be included in build"),
    "uninstall": re.compile(r"(?P<name>.*) needs to be removed from build"),
    "update": re.compile(r"(?P<name>.*) changed from (?P<oldver>.*) to (?P<newver>.*)"),
}


def _build_check_info(app_options):
    """Get info about packages scheduled for (un)install/update"""
    handler = _AppHandler(app_options)
    messages = handler.build_check(fast=True)
    # Decode the messages into a dict:
    status = {"install": [], "uninstall": [], "update": []}
    for msg in messages:
        for key, pattern in _message_map.items():
            match = pattern.match(msg)
            if match:
                status[key].append(match.group("name"))
    return status


@dataclass(frozen=True)
class ExtensionPackage:
    """Extension package entry.

    Attributes:
        name: Package name
        description: Package description
        homepage_url: Package home page
        pkg_type: Type of package - ["prebuilt", "source"]
        allowed: [optional] Whether this extension is allowed or not - default True
        approved: [optional] Whether the package is approved by your administrators - default False
        companion: [optional] Type of companion for the frontend extension - [None, "kernel", "server"]; default None
        core: [optional] Whether the package is a core package or not - default False
        enabled: [optional] Whether the package is enabled or not - default False
        install: [optional] Extension package installation instructions - default None
        installed: [optional] Whether the extension is currently installed - default None
        installed_version: [optional] Installed version - default ""
        latest_version: [optional] Latest available version - default ""
        status: [optional] Package status - ["ok", "warning", "error"]; default "ok"
        author: [optional] Package author - default None
        license: [optional] Package license - default None
        bug_tracker_url: [optional] Package bug tracker URL - default None
        documentation_url: [optional] Package documentation URL - default None
        package_manager_url: Package home page in the package manager - default None
        repository_url: [optional] Package code repository URL - default None
    """

    name: str
    description: str
    homepage_url: str
    pkg_type: str
    allowed: bool = True
    approved: bool = False
    companion: Optional[str] = None
    core: bool = False
    enabled: bool = False
    install: Optional[dict] = None
    installed: Optional[bool] = None
    installed_version: str = ""
    latest_version: str = ""
    status: str = "ok"
    author: Optional[str] = None
    license: Optional[str] = None
    bug_tracker_url: Optional[str] = None
    documentation_url: Optional[str] = None
    package_manager_url: Optional[str] = None
    repository_url: Optional[str] = None


@dataclass(frozen=True)
class ActionResult:
    """Action result

    Attributes:
        status: Action status - ["ok", "warning", "error"]
        message: Action status explanation
        needs_restart: Required action follow-up - Valid follow-up are "frontend", "kernel" and "server"
    """

    # Note: no simple way to use Enum in dataclass - https://stackoverflow.com/questions/72859557/typing-dataclass-that-can-only-take-enum-values
    #       keeping str for simplicity
    status: str
    message: Optional[str] = None
    needs_restart: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class PluginManagerOptions:
    """Plugin manager options.

    Attributes:
        lock_all: Whether to lock (prevent enabling/disabling) all plugins.
        lock_rules: A list of plugins or extensions that cannot be toggled.
            If extension name is provided, all its plugins will be disabled.
            The plugin names need to follow colon-separated format of `extension:plugin`.
    """

    lock_rules: frozenset[str] = field(default_factory=frozenset)
    lock_all: bool = False


@dataclass(frozen=True)
class ExtensionManagerOptions(PluginManagerOptions):
    """Extension manager options.

    Attributes:
        allowed_extensions_uris: A list of comma-separated URIs to get the allowed extensions list
        blocked_extensions_uris: A list of comma-separated URIs to get the blocked extensions list
        listings_refresh_seconds: The interval delay in seconds to refresh the lists
        listings_tornado_options: The optional kwargs to use for the listings HTTP requests as described on https://www.tornadoweb.org/en/stable/httpclient.html#tornado.httpclient.HTTPRequest
    """

    allowed_extensions_uris: set[str] = field(default_factory=set)
    blocked_extensions_uris: set[str] = field(default_factory=set)
    listings_refresh_seconds: int = 60 * 60
    listings_tornado_options: dict = field(default_factory=dict)


@dataclass(frozen=True)
class ExtensionManagerMetadata:
    """Extension manager metadata.

    Attributes:
        name: Extension manager name to be displayed
        can_install: Whether the extension manager can un-/install packages (default False)
        install_path: Installation path for the extensions (default None); e.g. environment path
    """

    name: str
    can_install: bool = False
    install_path: Optional[str] = None


@dataclass
class ExtensionsCache:
    """Extensions cache

    Attributes:
        cache: Extension list per page
        last_page: Last available page result
    """

    cache: dict[int, Optional[dict[str, ExtensionPackage]]] = field(default_factory=dict)
    last_page: int = 1


class PluginManager(LoggingConfigurable):
    """Plugin manager enables or disables plugins unless locked.

    It can also disable/enable all plugins in an extension.

    Args:
        app_options: Application options
        ext_options: Plugin manager (subset of extension manager) options
        parent: Configurable parent

    Attributes:
        app_options: Application options
        options: Plugin manager options
    """

    level = Enum(
        values=["sys_prefix", "user", "system"],
        default_value="sys_prefix",
        help="Level at which to manage plugins: sys_prefix, user, system",
    ).tag(config=True)

    def __init__(
        self,
        app_options: Optional[dict] = None,
        ext_options: Optional[dict] = None,
        parent: Optional[Configurable] = None,
    ) -> None:
        super().__init__(parent=parent)
        self.log.debug(
            f"Plugins in {self.__class__.__name__} will managed on the {self.level} level"
        )
        self.app_options = _ensure_options(app_options)
        plugin_options_field = {f.name for f in fields(PluginManagerOptions)}
        plugin_options = {
            option: value
            for option, value in (ext_options or {}).items()
            if option in plugin_options_field
        }
        self.options = PluginManagerOptions(**plugin_options)

    async def plugin_locks(self) -> dict:
        """Get information about locks on plugin enabling/disabling"""
        return {
            "lockRules": list(self.options.lock_rules),
            "allLocked": self.options.lock_all,
        }

    def _find_locked(self, plugins_or_extensions: list[str]) -> frozenset[str]:
        """Find a subset of plugins (or extensions) which are locked"""
        if self.options.lock_all:
            return set(plugins_or_extensions)
        locked_subset = set()
        extensions_with_locked_plugins = {
            plugin.split(":")[0] for plugin in self.options.lock_rules
        }
        for plugin in plugins_or_extensions:
            if ":" in plugin:
                # check directly if this is a plugin identifier (has colon)
                if plugin in self.options.lock_rules:
                    locked_subset.add(plugin)
            elif plugin in extensions_with_locked_plugins:
                # this is an extension - we need to check for >any< plugin
                # belonging to said extension
                locked_subset.add(plugin)
        return locked_subset

    async def disable(self, plugins: Union[str, list[str]]) -> ActionResult:
        """Disable a set of plugins (or an extension).

        Args:
            plugins: The list of plugins to disable
        Returns:
            The action result
        """
        plugins = plugins if isinstance(plugins, list) else [plugins]
        locked = self._find_locked(plugins)
        trans = translator.load("jupyterlab")
        if locked:
            return ActionResult(
                status="error",
                message=trans.gettext(
                    "The following plugins cannot be disabled as they are locked: "
                )
                + ", ".join(locked),
            )
        try:
            for plugin in plugins:
                disable_extension(plugin, app_options=self.app_options, level=self.level)
            return ActionResult(status="ok", needs_restart=["frontend"])
        except Exception as err:
            return ActionResult(status="error", message=repr(err))

    async def enable(self, plugins: Union[str, list[str]]) -> ActionResult:
        """Enable a set of plugins (or an extension).

        Args:
            plugins: The list of plugins to enable
        Returns:
            The action result
        """
        plugins = plugins if isinstance(plugins, list) else [plugins]
        locked = self._find_locked(plugins)
        trans = translator.load("jupyterlab")
        if locked:
            return ActionResult(
                status="error",
                message=trans.gettext(
                    "The following plugins cannot be enabled as they are locked: "
                )
                + ", ".join(locked),
            )
        try:
            for plugin in plugins:
                enable_extension(plugin, app_options=self.app_options, level=self.level)
            return ActionResult(status="ok", needs_restart=["frontend"])
        except Exception as err:
            return ActionResult(status="error", message=repr(err))


class ExtensionManager(PluginManager):
    """Base abstract extension manager.

    Note:
        Any concrete implementation will need to implement the five
        following abstract methods:
        - :ref:`metadata`
        - :ref:`get_latest_version`
        - :ref:`list_packages`
        - :ref:`install`
        - :ref:`uninstall`

        It could be interesting to override the :ref:`get_normalized_name`
        method too.

    Args:
        app_options: Application options
        ext_options: Extension manager options
        parent: Configurable parent

    Attributes:
        log: Logger
        app_dir: Application directory
        core_config: Core configuration
        app_options: Application options
        options: Extension manager options
    """

    def __init__(
        self,
        app_options: Optional[dict] = None,
        ext_options: Optional[dict] = None,
        parent: Optional[Configurable] = None,
    ) -> None:
        super().__init__(app_options=app_options, ext_options=ext_options, parent=parent)
        self.log = self.app_options.logger
        self.app_dir = Path(self.app_options.app_dir)
        self.core_config = self.app_options.core_config
        self.options = ExtensionManagerOptions(**(ext_options or {}))
        self._extensions_cache: dict[Optional[str], ExtensionsCache] = {}
        self._listings_cache: Optional[dict] = None
        self._listings_block_mode = True
        self._listing_fetch: Optional[tornado.ioloop.PeriodicCallback] = None

        if len(self.options.allowed_extensions_uris) or len(self.options.blocked_extensions_uris):
            self._listings_block_mode = len(self.options.allowed_extensions_uris) == 0
            if not self._listings_block_mode and len(self.options.blocked_extensions_uris) > 0:
                self.log.warning(
                    "You have define simultaneously blocked and allowed extensions listings. The allowed listing will take precedence."
                )

            self._listing_fetch = tornado.ioloop.PeriodicCallback(
                self._fetch_listings,
                callback_time=self.options.listings_refresh_seconds * 1000,
                jitter=0.1,
            )
            self._listing_fetch.start()

    def __del__(self):
        if self._listing_fetch is not None:
            self._listing_fetch.stop()

    @property
    def metadata(self) -> ExtensionManagerMetadata:
        """Extension manager metadata."""
        raise NotImplementedError()

    async def get_latest_version(self, extension: str) -> Optional[str]:
        """Return the latest available version for a given extension.

        Args:
            pkg: The extension name
        Returns:
            The latest available version
        """
        raise NotImplementedError()

    async def list_packages(
        self, query: str, page: int, per_page: int
    ) -> tuple[dict[str, ExtensionPackage], Optional[int]]:
        """List the available extensions.

        Args:
            query: The search extension query
            page: The result page
            per_page: The number of results per page
        Returns:
            The available extensions in a mapping {name: metadata}
            The results last page; None if the manager does not support pagination
        """
        raise NotImplementedError()

    async def install(self, extension: str, version: Optional[str] = None) -> ActionResult:
        """Install the required extension.

        Note:
            If the user must be notified with a message (like asking to restart the
            server), the result should be
            {"status": "warning", "message": "<explanation for the user>"}

        Args:
            extension: The extension name
            version: The version to install; default None (i.e. the latest possible)
        Returns:
            The action result
        """
        raise NotImplementedError()

    async def uninstall(self, extension: str) -> ActionResult:
        """Uninstall the required extension.

        Note:
            If the user must be notified with a message (like asking to restart the
            server), the result should be
            {"status": "warning", "message": "<explanation for the user>"}

        Args:
            extension: The extension name
        Returns:
            The action result
        """
        raise NotImplementedError()

    @staticmethod
    def get_semver_version(version: str) -> str:
        """Convert a Python version to Semver version.

        It:

        - drops ``.devN`` and ``.postN``
        - converts ``aN``, ``bN`` and ``rcN`` to ``-alpha.N``, ``-beta.N``, ``-rc.N`` respectively

        Args:
            version: Version to convert
        Returns
            Semver compatible version
        """
        return re.sub(
            r"(a|b|rc)(\d+)$",
            lambda m: f"{PYTHON_TO_SEMVER[m.group(1)]}{m.group(2)}",
            re.subn(r"\.(dev|post)\d+", "", version)[0],
        )

    def get_normalized_name(self, extension: ExtensionPackage) -> str:
        """Normalize extension name.

        Extension have multiple parts, npm package, Python package,...
        Sub-classes may override this method to ensure the name of
        an extension from the service provider and the local installed
        listing is matching.

        Args:
            extension: The extension metadata
        Returns:
            The normalized name
        """
        return extension.name

    async def list_extensions(
        self, query: Optional[str] = None, page: int = 1, per_page: int = 30
    ) -> tuple[list[ExtensionPackage], Optional[int]]:
        """List extensions for a given ``query`` search term.

        This will return the extensions installed (if ``query`` is None) or
        available if allowed by the listing settings.

        Args:
            query: [optional] Query search term.

        Returns:
            The extensions
            Last page of results
        """
        if query not in self._extensions_cache or page not in self._extensions_cache[query].cache:
            await self.refresh(query, page, per_page)

        # filter using listings settings
        if self._listings_cache is None and self._listing_fetch is not None:
            await self._listing_fetch.callback()

        cache = self._extensions_cache[query].cache[page]
        if cache is None:
            cache = {}
        extensions = list(cache.values())
        if query is not None and self._listings_cache is not None:
            listing = list(self._listings_cache)
            extensions = []
            if self._listings_block_mode:
                for name, ext in cache.items():
                    if name not in listing:
                        extensions.append(replace(ext, allowed=True))
                    elif ext.installed_version:
                        self.log.warning(f"Blocked extension '{name}' is installed.")
                        extensions.append(replace(ext, allowed=False))
            else:
                for name, ext in cache.items():
                    if name in listing:
                        extensions.append(replace(ext, allowed=True))
                    elif ext.installed_version:
                        self.log.warning(f"Not allowed extension '{name}' is installed.")
                        extensions.append(replace(ext, allowed=False))

        return extensions, self._extensions_cache[query].last_page

    async def refresh(self, query: Optional[str], page: int, per_page: int) -> None:
        """Refresh the list of extensions."""
        if query in self._extensions_cache:
            self._extensions_cache[query].cache[page] = None
        await self._update_extensions_list(query, page, per_page)

    async def _fetch_listings(self) -> None:
        """Fetch the listings for the extension manager."""
        rules = []
        client = tornado.httpclient.AsyncHTTPClient()
        if self._listings_block_mode:
            if len(self.options.blocked_extensions_uris):
                self.log.info(
                    f"Fetching blocked extensions from {self.options.blocked_extensions_uris}"
                )
                for blocked_extensions_uri in self.options.blocked_extensions_uris:
                    r = await client.fetch(
                        blocked_extensions_uri,
                        **self.options.listings_tornado_options,
                    )
                    j = json.loads(r.body)
                    rules.extend(j.get("blocked_extensions", []))
        elif len(self.options.allowed_extensions_uris):
            self.log.info(
                f"Fetching allowed extensions from {self.options.allowed_extensions_uris}"
            )
            for allowed_extensions_uri in self.options.allowed_extensions_uris:
                r = await client.fetch(
                    allowed_extensions_uri,
                    **self.options.listings_tornado_options,
                )
                j = json.loads(r.body)
                rules.extend(j.get("allowed_extensions", []))

        self._listings_cache = {r["name"]: r for r in rules}

    async def _get_installed_extensions(
        self, get_latest_version=True
    ) -> dict[str, ExtensionPackage]:
        """Get the installed extensions.

        Args:
            get_latest_version: Whether to fetch the latest extension version or not.
        Returns:
            The installed extensions as a mapping {name: metadata}
        """
        app_options = self.app_options
        info = get_app_info(app_options=app_options)
        build_check_info = _build_check_info(app_options)
        _ensure_compat_errors(info, app_options)
        extensions = {}

        # TODO: the three for-loops below can be run concurrently
        for name, data in info["federated_extensions"].items():
            status = "ok"
            pkg_info = data
            if info["compat_errors"].get(name, None):
                status = "error"

            normalized_name = self._normalize_name(name)
            pkg = ExtensionPackage(
                name=normalized_name,
                description=pkg_info.get("description", ""),
                homepage_url=data.get("url", ""),
                enabled=(name not in info["disabled"]),
                core=False,
                latest_version=ExtensionManager.get_semver_version(data["version"]),
                installed=True,
                installed_version=ExtensionManager.get_semver_version(data["version"]),
                status=status,
                install=data.get("install", {}),
                pkg_type="prebuilt",
                companion=self._get_companion(data),
                author=data.get("author", {}).get("name", data.get("author")),
                license=data.get("license"),
                bug_tracker_url=data.get("bugs", {}).get("url"),
                repository_url=data.get("repository", {}).get("url", data.get("repository")),
            )

            if get_latest_version:
                pkg = replace(pkg, latest_version=await self.get_latest_version(pkg.name))

            extensions[normalized_name] = pkg

        for name, data in info["extensions"].items():
            if name in info["shadowed_exts"]:
                continue
            status = "ok"

            if info["compat_errors"].get(name, None):
                status = "error"
            else:
                for packages in build_check_info.values():
                    if name in packages:
                        status = "warning"

            normalized_name = self._normalize_name(name)
            pkg = ExtensionPackage(
                name=normalized_name,
                description=data.get("description", ""),
                homepage_url=data["url"],
                enabled=(name not in info["disabled"]),
                core=False,
                latest_version=ExtensionManager.get_semver_version(data["version"]),
                installed=True,
                installed_version=ExtensionManager.get_semver_version(data["version"]),
                status=status,
                pkg_type="source",
                companion=self._get_companion(data),
                author=data.get("author", {}).get("name", data.get("author")),
                license=data.get("license"),
                bug_tracker_url=data.get("bugs", {}).get("url"),
                repository_url=data.get("repository", {}).get("url", data.get("repository")),
            )
            if get_latest_version:
                pkg = replace(pkg, latest_version=await self.get_latest_version(pkg.name))
            extensions[normalized_name] = pkg

        for name in build_check_info["uninstall"]:
            data = self._get_scheduled_uninstall_info(name)
            if data is not None:
                normalized_name = self._normalize_name(name)
                pkg = ExtensionPackage(
                    name=normalized_name,
                    description=data.get("description", ""),
                    homepage_url=data.get("homepage", ""),
                    installed=False,
                    enabled=False,
                    core=False,
                    latest_version=ExtensionManager.get_semver_version(data["version"]),
                    installed_version=ExtensionManager.get_semver_version(data["version"]),
                    status="warning",
                    pkg_type="prebuilt",
                    author=data.get("author", {}).get("name", data.get("author")),
                    license=data.get("license"),
                    bug_tracker_url=data.get("bugs", {}).get("url"),
                    repository_url=data.get("repository", {}).get("url", data.get("repository")),
                )
                extensions[normalized_name] = pkg

        return extensions

    def _get_companion(self, data: dict) -> Optional[str]:
        companion = None
        if "discovery" in data["jupyterlab"]:
            if "server" in data["jupyterlab"]["discovery"]:
                companion = "server"
            elif "kernel" in data["jupyterlab"]["discovery"]:
                companion = "kernel"
        return companion

    def _get_scheduled_uninstall_info(self, name) -> Optional[dict]:
        """Get information about a package that is scheduled for uninstallation"""
        target = self.app_dir / "staging" / "node_modules" / name / "package.json"
        if target.exists():
            with target.open() as fid:
                return json.load(fid)
        else:
            return None

    def _normalize_name(self, name: str) -> str:
        """Normalize extension name; by default does nothing.

        Args:
            name: Extension name
        Returns:
            Normalized name
        """
        return name

    async def _update_extensions_list(
        self, query: Optional[str] = None, page: int = 1, per_page: int = 30
    ) -> None:
        """Update the list of extensions"""
        last_page = None
        if query is not None:
            # Get the available extensions
            extensions, last_page = await self.list_packages(query, page, per_page)
        else:
            # Get the installed extensions
            extensions = await self._get_installed_extensions()

        if query in self._extensions_cache:
            self._extensions_cache[query].cache[page] = extensions
            self._extensions_cache[query].last_page = last_page or 1
        else:
            self._extensions_cache[query] = ExtensionsCache({page: extensions}, last_page or 1)
