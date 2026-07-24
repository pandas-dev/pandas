"""JupyterLab Server config"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import json
import os.path as osp
from glob import iglob
from itertools import chain
from logging import Logger
from os.path import join as pjoin
from typing import Any

import json5
from jupyter_core.paths import SYSTEM_CONFIG_PATH, jupyter_config_dir, jupyter_path
from jupyter_server.services.config.manager import ConfigManager, recursive_update
from jupyter_server.utils import url_path_join as ujoin
from traitlets import Bool, HasTraits, List, Unicode, default

# -----------------------------------------------------------------------------
# Module globals
# -----------------------------------------------------------------------------

DEFAULT_TEMPLATE_PATH = osp.join(osp.dirname(__file__), "templates")


def get_package_url(data: dict[str, Any]) -> str:
    """Get the url from the extension data"""
    # homepage, repository  are optional
    if "homepage" in data:
        url = data["homepage"]
    elif "repository" in data and isinstance(data["repository"], dict):
        url = data["repository"].get("url", "")
    else:
        url = ""
    return url


def get_federated_extensions(labextensions_path: list[str]) -> dict[str, Any]:
    """Get the metadata about federated extensions"""
    federated_extensions = {}
    for ext_dir in labextensions_path:
        # extensions are either top-level directories, or two-deep in @org directories
        for ext_path in chain(
            iglob(pjoin(ext_dir, "[!@]*", "package.json")),
            iglob(pjoin(ext_dir, "@*", "*", "package.json")),
        ):
            with open(ext_path, encoding="utf-8") as fid:
                pkgdata = json.load(fid)
            if pkgdata["name"] not in federated_extensions:
                data = dict(
                    name=pkgdata["name"],
                    version=pkgdata["version"],
                    description=pkgdata.get("description", ""),
                    url=get_package_url(pkgdata),
                    ext_dir=ext_dir,
                    ext_path=osp.dirname(ext_path),
                    is_local=False,
                    dependencies=pkgdata.get("dependencies", dict()),
                    jupyterlab=pkgdata.get("jupyterlab", dict()),
                )

                # Add repository info if available
                if "repository" in pkgdata and "url" in pkgdata.get("repository", {}):
                    data["repository"] = dict(url=pkgdata.get("repository").get("url"))

                install_path = osp.join(osp.dirname(ext_path), "install.json")
                if osp.exists(install_path):
                    with open(install_path, encoding="utf-8") as fid:
                        data["install"] = json.load(fid)
                federated_extensions[data["name"]] = data
    return federated_extensions


def get_static_page_config(
    app_settings_dir: str | None = None,  # noqa: ARG001
    logger: Logger | None = None,  # noqa: ARG001
    level: str = "all",
    include_higher_levels: bool = False,
) -> dict[str, Any]:
    """Get the static page config for JupyterLab

    Parameters
    ----------
    logger: logger, optional
        An optional logging object
    level: string, optional ['all']
        The level at which to get config: can be 'all', 'user', 'sys_prefix', or 'system'
    """
    cm = _get_config_manager(level, include_higher_levels)
    return cm.get("page_config")  # type:ignore[no-untyped-call]


def load_config(path: str) -> Any:
    """Load either a json5 or a json config file.

    Parameters
    ----------
    path : str
        Path to the file to be loaded

    Returns
    -------
    Dict[Any, Any]
        Dictionary of json or json5 data
    """
    with open(path, encoding="utf-8") as fid:
        if path.endswith(".json5"):
            return json5.load(fid)
        return json.load(fid)


def get_page_config(
    labextensions_path: list[str], app_settings_dir: str | None = None, logger: Logger | None = None
) -> dict[str, Any]:
    """Get the page config for the application handler"""
    # Build up the full page config
    page_config: dict = {}

    disabled_key = "disabledExtensions"

    # Start with the app_settings_dir as lowest priority
    if app_settings_dir:
        config_paths = [
            pjoin(app_settings_dir, "page_config.json5"),
            pjoin(app_settings_dir, "page_config.json"),
        ]
        for path in config_paths:
            if osp.exists(path) and osp.getsize(path):
                data = load_config(path)
                # Convert lists to dicts
                for key in [disabled_key, "deferredExtensions"]:
                    if key in data:
                        data[key] = {key: True for key in data[key]}

                recursive_update(page_config, data)
                break

    # Get the traitlets config
    static_page_config = get_static_page_config(logger=logger, level="all")
    recursive_update(page_config, static_page_config)

    # Handle federated extensions that disable other extensions
    disabled_by_extensions_all = {}
    extensions = page_config["federated_extensions"] = []

    federated_exts = get_federated_extensions(labextensions_path)

    # Ensure there is a disabled key
    page_config.setdefault(disabled_key, {})

    for _, ext_data in federated_exts.items():
        if "_build" not in ext_data["jupyterlab"]:
            if logger:
                logger.warning("%s is not a valid extension", ext_data["name"])
            continue
        extbuild = ext_data["jupyterlab"]["_build"]
        extension = {"name": ext_data["name"], "load": extbuild["load"]}

        if "extension" in extbuild:
            extension["extension"] = extbuild["extension"]
        if "mimeExtension" in extbuild:
            extension["mimeExtension"] = extbuild["mimeExtension"]
        if "style" in extbuild:
            extension["style"] = extbuild["style"]
        # FIXME @experimental for plugin with no-code entrypoints.
        extension["entrypoints"] = extbuild.get("entrypoints")
        extensions.append(extension)

        # If there is disabledExtensions metadata, consume it.
        name = ext_data["name"]

        if ext_data["jupyterlab"].get(disabled_key):
            disabled_by_extensions_all[ext_data["name"]] = ext_data["jupyterlab"][disabled_key]

    # Handle source extensions that disable other extensions
    # Check for `jupyterlab`:`extensionMetadata` in the built application directory's package.json
    if app_settings_dir:
        app_dir = osp.dirname(app_settings_dir)
        package_data_file = pjoin(app_dir, "static", "package.json")
        if osp.exists(package_data_file):
            with open(package_data_file, encoding="utf-8") as fid:
                app_data = json.load(fid)
            all_ext_data = app_data["jupyterlab"].get("extensionMetadata", {})
            for ext, ext_data in all_ext_data.items():
                if ext in disabled_by_extensions_all:
                    continue
                if ext_data.get(disabled_key):
                    disabled_by_extensions_all[ext] = ext_data[disabled_key]

    disabled_by_extensions = {}
    for name in sorted(disabled_by_extensions_all):
        # skip if the extension itself is disabled by other config
        if page_config[disabled_key].get(name) is True:
            continue

        disabled_list = disabled_by_extensions_all[name]
        for item in disabled_list:
            disabled_by_extensions[item] = True

    rollup_disabled = disabled_by_extensions
    rollup_disabled.update(page_config.get(disabled_key, []))
    page_config[disabled_key] = rollup_disabled

    # Convert dictionaries to lists to give to the front end
    for key, value in page_config.items():
        if isinstance(value, dict):
            page_config[key] = [subkey for subkey in value if value[subkey]]

    return page_config


def write_page_config(page_config: dict[str, Any], level: str = "all") -> None:
    """Write page config to disk"""
    cm = _get_config_manager(level)
    cm.set("page_config", page_config)  # type:ignore[no-untyped-call]


class LabConfig(HasTraits):
    """The lab application configuration object."""

    app_name = Unicode("", help="The name of the application.").tag(config=True)

    app_version = Unicode("", help="The version of the application.").tag(config=True)

    app_namespace = Unicode("", help="The namespace of the application.").tag(config=True)

    app_url = Unicode("/lab", help="The url path for the application.").tag(config=True)

    app_settings_dir = Unicode("", help="The application settings directory.").tag(config=True)

    extra_labextensions_path = List(
        Unicode(), help="""Extra paths to look for federated JupyterLab extensions"""
    ).tag(config=True)

    labextensions_path = List(
        Unicode(), help="The standard paths to look in for federated JupyterLab extensions"
    ).tag(config=True)

    templates_dir = Unicode("", help="The application templates directory.").tag(config=True)

    static_dir = Unicode(
        "",
        help=(
            "The optional location of local static files. "
            "If given, a static file handler will be "
            "added."
        ),
    ).tag(config=True)

    labextensions_url = Unicode("", help="The url for federated JupyterLab extensions").tag(
        config=True
    )

    settings_url = Unicode(help="The url path of the settings handler.").tag(config=True)

    user_settings_dir = Unicode(
        "", help=("The optional location of the user settings directory.")
    ).tag(config=True)

    schemas_dir = Unicode(
        "",
        help=(
            "The optional location of the settings "
            "schemas directory. If given, a handler will "
            "be added for settings."
        ),
    ).tag(config=True)

    workspaces_api_url = Unicode(help="The url path of the workspaces API.").tag(config=True)

    workspaces_dir = Unicode(
        "",
        help=(
            "The optional location of the saved "
            "workspaces directory. If given, a handler "
            "will be added for workspaces."
        ),
    ).tag(config=True)

    listings_url = Unicode(help="The listings url.").tag(config=True)

    themes_url = Unicode(help="The theme url.").tag(config=True)

    licenses_url = Unicode(help="The third-party licenses url.")

    themes_dir = Unicode(
        "",
        help=(
            "The optional location of the themes "
            "directory. If given, a handler will be added "
            "for themes."
        ),
    ).tag(config=True)

    translations_api_url = Unicode(help="The url path of the translations handler.").tag(
        config=True
    )

    tree_url = Unicode(help="The url path of the tree handler.").tag(config=True)

    cache_files = Bool(
        True,
        help=("Whether to cache files on the server. This should be `True` except in dev mode."),
    ).tag(config=True)

    notebook_starts_kernel = Bool(
        True, help="Whether a notebook should start a kernel automatically."
    ).tag(config=True)

    copy_absolute_path = Bool(
        False,
        help="Whether getting a relative (False) or absolute (True) path when copying a path.",
    ).tag(config=True)

    @default("templates_dir")
    def _default_templates_dir(self) -> str:
        return DEFAULT_TEMPLATE_PATH

    @default("labextensions_url")
    def _default_labextensions_url(self) -> str:
        return ujoin(self.app_url, "extensions/")

    @default("labextensions_path")
    def _default_labextensions_path(self) -> list[str]:
        return jupyter_path("labextensions")

    @default("workspaces_url")
    def _default_workspaces_url(self) -> str:
        return ujoin(self.app_url, "workspaces/")

    @default("workspaces_api_url")
    def _default_workspaces_api_url(self) -> str:
        return ujoin(self.app_url, "api", "workspaces/")

    @default("settings_url")
    def _default_settings_url(self) -> str:
        return ujoin(self.app_url, "api", "settings/")

    @default("listings_url")
    def _default_listings_url(self) -> str:
        return ujoin(self.app_url, "api", "listings/")

    @default("themes_url")
    def _default_themes_url(self) -> str:
        return ujoin(self.app_url, "api", "themes/")

    @default("licenses_url")
    def _default_licenses_url(self) -> str:
        return ujoin(self.app_url, "api", "licenses/")

    @default("tree_url")
    def _default_tree_url(self) -> str:
        return ujoin(self.app_url, "tree/")

    @default("translations_api_url")
    def _default_translations_api_url(self) -> str:
        return ujoin(self.app_url, "api", "translations/")


def get_allowed_levels() -> list[str]:
    """
    Returns the levels where configs can be stored.
    """
    return ["all", "user", "sys_prefix", "system", "app", "extension"]


def _get_config_manager(level: str, include_higher_levels: bool = False) -> ConfigManager:
    """Get the location of config files for the current context
    Returns the string to the environment
    """
    # Delayed import since this gets monkey-patched in tests
    from jupyter_core.paths import ENV_CONFIG_PATH

    allowed = get_allowed_levels()
    if level not in allowed:
        msg = f"Page config level must be one of: {allowed}"
        raise ValueError(msg)

    config_name = "labconfig"

    if level == "all":
        return ConfigManager(config_dir_name=config_name)

    paths: dict[str, list] = {
        "app": [],
        "system": SYSTEM_CONFIG_PATH,
        "sys_prefix": [ENV_CONFIG_PATH[0]],
        "user": [jupyter_config_dir()],
        "extension": [],
    }

    levels = allowed[allowed.index(level) :] if include_higher_levels else [level]

    read_config_paths, write_config_dir = [], None

    for _level in levels:
        for p in paths[_level]:
            read_config_paths.append(osp.join(p, config_name))
        if write_config_dir is None and paths[_level]:  # type: ignore[redundant-expr]
            write_config_dir = osp.join(paths[_level][0], config_name)

    return ConfigManager(read_config_path=read_config_paths, write_config_dir=write_config_dir)
