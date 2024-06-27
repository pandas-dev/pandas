"""Frontend config storage helpers."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import json
import os
from glob import glob
from typing import Any

import json5
from jsonschema import Draft7Validator as Validator
from jsonschema import ValidationError
from jupyter_server import _tz as tz
from jupyter_server.base.handlers import APIHandler
from jupyter_server.services.config.manager import ConfigManager, recursive_update
from tornado import web

from .translation_utils import DEFAULT_LOCALE, L10N_SCHEMA_NAME, SYS_LOCALE, is_valid_locale

# The JupyterLab settings file extension.
SETTINGS_EXTENSION = ".jupyterlab-settings"


def _get_schema(
    schemas_dir: str,
    schema_name: str,
    overrides: dict[str, Any],
    labextensions_path: list[str] | None,
) -> tuple[dict[str, Any], str]:
    """Returns a dict containing a parsed and validated JSON schema."""
    notfound_error = "Schema not found: %s"
    parse_error = "Failed parsing schema (%s): %s"
    validation_error = "Failed validating schema (%s): %s"

    path = None

    # Look for the setting in all of the labextension paths first
    # Use the first one
    if labextensions_path is not None:
        ext_name, _, plugin_name = schema_name.partition(":")
        for ext_path in labextensions_path:
            target = os.path.join(ext_path, ext_name, "schemas", ext_name, plugin_name + ".json")
            if os.path.exists(target):
                schemas_dir = os.path.join(ext_path, ext_name, "schemas")
                path = target
                break

    # Fall back on the default location
    if path is None:
        path = _path(schemas_dir, schema_name)

    if not os.path.exists(path):
        raise web.HTTPError(404, notfound_error % path)

    with open(path, encoding="utf-8") as fid:
        # Attempt to load the schema file.
        try:
            schema = json.load(fid)
        except Exception as e:
            name = schema_name
            raise web.HTTPError(500, parse_error % (name, str(e))) from None

    schema = _override(schema_name, schema, overrides)

    # Validate the schema.
    try:
        Validator.check_schema(schema)
    except Exception as e:
        name = schema_name
        raise web.HTTPError(500, validation_error % (name, str(e))) from None

    version = _get_version(schemas_dir, schema_name)

    return schema, version


def _get_user_settings(settings_dir: str, schema_name: str, schema: Any) -> dict[str, Any]:
    """
    Returns a dictionary containing the raw user settings, the parsed user
    settings, a validation warning for a schema, and file times.
    """
    path = _path(settings_dir, schema_name, False, SETTINGS_EXTENSION)
    raw = "{}"
    settings = {}
    warning = None
    validation_warning = "Failed validating settings (%s): %s"
    parse_error = "Failed loading settings (%s): %s"
    last_modified = None
    created = None

    if os.path.exists(path):
        stat = os.stat(path)
        last_modified = tz.utcfromtimestamp(stat.st_mtime).isoformat()
        created = tz.utcfromtimestamp(stat.st_ctime).isoformat()
        with open(path, encoding="utf-8") as fid:
            try:  # to load and parse the settings file.
                raw = fid.read() or raw
                settings = json5.loads(raw)
            except Exception as e:
                raise web.HTTPError(500, parse_error % (schema_name, str(e))) from None

    # Validate the parsed data against the schema.
    if len(settings):
        validator = Validator(schema)
        try:
            validator.validate(settings)
        except ValidationError as e:
            warning = validation_warning % (schema_name, str(e))
            raw = "{}"
            settings = {}

    return dict(
        raw=raw, settings=settings, warning=warning, last_modified=last_modified, created=created
    )


def _get_version(schemas_dir: str, schema_name: str) -> str:
    """Returns the package version for a given schema or 'N/A' if not found."""

    path = _path(schemas_dir, schema_name)
    package_path = os.path.join(os.path.split(path)[0], "package.json.orig")

    try:  # to load and parse the package.json.orig file.
        with open(package_path, encoding="utf-8") as fid:
            package = json.load(fid)
            return package["version"]
    except Exception:
        return "N/A"


def _list_settings(
    schemas_dir: str,
    settings_dir: str,
    overrides: dict[str, Any],
    extension: str = ".json",
    labextensions_path: list[str] | None = None,
    translator: Any = None,
    ids_only: bool = False,
) -> tuple[list[Any], list[Any]]:
    """
    Returns a tuple containing:
     - the list of plugins, schemas, and their settings,
       respecting any defaults that may have been overridden if `ids_only=False`,
       otherwise a list of dict containing only the ids of plugins.
     - the list of warnings that were generated when
       validating the user overrides against the schemas.
    """

    settings: dict[str, Any] = {}
    federated_settings: dict[str, Any] = {}
    warnings = []

    if not os.path.exists(schemas_dir):
        warnings = ["Settings directory does not exist at %s" % schemas_dir]
        return ([], warnings)

    schema_pattern = schemas_dir + "/**/*" + extension
    schema_paths = [path for path in glob(schema_pattern, recursive=True)]  # noqa: C416
    schema_paths.sort()

    for schema_path in schema_paths:
        # Generate the schema_name used to request individual settings.
        rel_path = os.path.relpath(schema_path, schemas_dir)
        rel_schema_dir, schema_base = os.path.split(rel_path)
        _id = schema_name = ":".join(
            [rel_schema_dir, schema_base[: -len(extension)]]  # Remove file extension.
        ).replace("\\", "/")  # Normalize slashes.

        if ids_only:
            settings[_id] = dict(id=_id)
        else:
            schema, version = _get_schema(schemas_dir, schema_name, overrides, None)
            if translator is not None:
                schema = translator(schema)
            user_settings = _get_user_settings(settings_dir, schema_name, schema)

            if user_settings["warning"]:
                warnings.append(user_settings.pop("warning"))

            # Add the plugin to the list of settings.
            settings[_id] = dict(id=_id, schema=schema, version=version, **user_settings)

    if labextensions_path is not None:
        schema_paths = []
        for ext_dir in labextensions_path:
            schema_pattern = ext_dir + "/**/schemas/**/*" + extension
            schema_paths.extend(path for path in glob(schema_pattern, recursive=True))

        schema_paths.sort()

        for schema_path_ in schema_paths:
            schema_path = schema_path_.replace(os.sep, "/")

            base_dir, rel_path = schema_path.split("schemas/")

            # Generate the schema_name used to request individual settings.
            rel_schema_dir, schema_base = os.path.split(rel_path)
            _id = schema_name = ":".join(
                [rel_schema_dir, schema_base[: -len(extension)]]  # Remove file extension.
            ).replace("\\", "/")  # Normalize slashes.

            # bail if we've already handled the highest federated setting
            if _id in federated_settings:
                continue

            if ids_only:
                federated_settings[_id] = dict(id=_id)
            else:
                schema, version = _get_schema(
                    schemas_dir, schema_name, overrides, labextensions_path=labextensions_path
                )
                user_settings = _get_user_settings(settings_dir, schema_name, schema)

                if user_settings["warning"]:
                    warnings.append(user_settings.pop("warning"))

                # Add the plugin to the list of settings.
                federated_settings[_id] = dict(
                    id=_id, schema=schema, version=version, **user_settings
                )

    settings.update(federated_settings)
    settings_list = [settings[key] for key in sorted(settings.keys(), reverse=True)]

    return (settings_list, warnings)


def _override(
    schema_name: str, schema: dict[str, Any], overrides: dict[str, Any]
) -> dict[str, Any]:
    """Override default values in the schema if necessary."""
    if schema_name in overrides:
        defaults = overrides[schema_name]
        for key in defaults:
            if key in schema["properties"]:
                new_defaults = schema["properties"][key]["default"]
                # If values for defaults are dicts do a recursive update
                if isinstance(new_defaults, dict):
                    recursive_update(new_defaults, defaults[key])
                else:
                    new_defaults = defaults[key]

                schema["properties"][key]["default"] = new_defaults
            else:
                schema["properties"][key] = dict(default=defaults[key])

    return schema


def _path(
    root_dir: str, schema_name: str, make_dirs: bool = False, extension: str = ".json"
) -> str:
    """
    Returns the local file system path for a schema name in the given root
    directory. This function can be used to filed user overrides in addition to
    schema files. If the `make_dirs` flag is set to `True` it will create the
    parent directory for the calculated path if it does not exist.
    """

    notfound_error = "Settings not found (%s)"
    write_error = "Failed writing settings (%s): %s"

    try:  # to parse path, e.g. @jupyterlab/apputils-extension:themes.
        package_dir, plugin = schema_name.split(":")
        parent_dir = os.path.join(root_dir, package_dir)
        path = os.path.join(parent_dir, plugin + extension)
    except Exception:
        raise web.HTTPError(404, notfound_error % schema_name) from None

    if make_dirs and not os.path.exists(parent_dir):
        try:
            os.makedirs(parent_dir)
        except Exception as e:
            raise web.HTTPError(500, write_error % (schema_name, str(e))) from None

    return path


def _get_overrides(app_settings_dir: str) -> tuple[dict[str, Any], str]:
    """Get overrides settings from `app_settings_dir`.

    The ordering of paths is:
    - {app_settings_dir}/overrides.d/*.{json,json5} (many, namespaced by package)
    - {app_settings_dir}/overrides.{json,json5} (singleton, owned by the user)
    """
    overrides: dict[str, Any]
    error: str
    overrides, error = {}, ""

    overrides_d = os.path.join(app_settings_dir, "overrides.d")

    # find (and sort) the conf.d overrides files
    all_override_paths = sorted(
        [
            *(glob(os.path.join(overrides_d, "*.json"))),
            *(glob(os.path.join(overrides_d, "*.json5"))),
        ]
    )

    all_override_paths += [
        os.path.join(app_settings_dir, "overrides.json"),
        os.path.join(app_settings_dir, "overrides.json5"),
    ]

    for overrides_path in all_override_paths:
        if not os.path.exists(overrides_path):
            continue

        with open(overrides_path, encoding="utf-8") as fid:
            try:
                if overrides_path.endswith(".json5"):
                    path_overrides = json5.load(fid)
                else:
                    path_overrides = json.load(fid)
                for plugin_id, config in path_overrides.items():
                    recursive_update(overrides.setdefault(plugin_id, {}), config)
            except Exception as e:
                error = e  # type:ignore[assignment]

    # Allow `default_settings_overrides.json` files in <jupyter_config>/labconfig dirs
    # to allow layering of defaults
    cm = ConfigManager(config_dir_name="labconfig")

    for plugin_id, config in cm.get("default_setting_overrides").items():  # type:ignore[no-untyped-call]
        recursive_update(overrides.setdefault(plugin_id, {}), config)

    return overrides, error


def get_settings(
    app_settings_dir: str,
    schemas_dir: str,
    settings_dir: str,
    schema_name: str = "",
    overrides: dict[str, Any] | None = None,
    labextensions_path: list[str] | None = None,
    translator: Any = None,
    ids_only: bool = False,
) -> tuple[dict[str, Any], list[Any]]:
    """
    Get settings.

    Parameters
    ----------
    app_settings_dir:
        Path to applications settings.
    schemas_dir: str
        Path to schemas.
    settings_dir:
        Path to settings.
    schema_name str, optional
        Schema name. Default is "".
    overrides: dict, optional
        Settings overrides. If not provided, the overrides will be loaded
        from the `app_settings_dir`. Default is None.
    labextensions_path: list, optional
        List of paths to federated labextensions containing their own schema files.
    translator: Callable[[Dict], Dict] or None, optional
        Translate a schema. It requires the schema dictionary and returns its translation

    Returns
    -------
    tuple
        The first item is a dictionary with a list of setting if no `schema_name`
        was provided (only the ids if `ids_only=True`), otherwise it is a dictionary
        with id, raw, scheme, settings and version keys.
        The second item is a list of warnings. Warnings will either be a list of
        i) strings with the warning messages or ii) `None`.
    """
    result = {}
    warnings = []

    if overrides is None:
        overrides, _error = _get_overrides(app_settings_dir)

    if schema_name:
        schema, version = _get_schema(schemas_dir, schema_name, overrides, labextensions_path)
        if translator is not None:
            schema = translator(schema)
        user_settings = _get_user_settings(settings_dir, schema_name, schema)
        warnings = [user_settings.pop("warning")]
        result = {"id": schema_name, "schema": schema, "version": version, **user_settings}
    else:
        settings_list, warnings = _list_settings(
            schemas_dir,
            settings_dir,
            overrides,
            labextensions_path=labextensions_path,
            translator=translator,
            ids_only=ids_only,
        )
        result = {
            "settings": settings_list,
        }

    return result, warnings


def save_settings(
    schemas_dir: str,
    settings_dir: str,
    schema_name: str,
    raw_settings: str,
    overrides: dict[str, Any],
    labextensions_path: list[str] | None = None,
) -> None:
    """
    Save ``raw_settings`` settings for ``schema_name``.

    Parameters
    ----------
    schemas_dir: str
        Path to schemas.
    settings_dir: str
        Path to settings.
    schema_name str
        Schema name.
    raw_settings: str
        Raw serialized settings dictionary
    overrides: dict
        Settings overrides.
    labextensions_path: list, optional
        List of paths to federated labextensions containing their own schema files.
    """
    payload = json5.loads(raw_settings)

    # Validate the data against the schema.
    schema, _ = _get_schema(
        schemas_dir, schema_name, overrides, labextensions_path=labextensions_path
    )
    validator = Validator(schema)
    validator.validate(payload)

    # Write the raw data (comments included) to a file.
    path = _path(settings_dir, schema_name, True, SETTINGS_EXTENSION)
    with open(path, "w", encoding="utf-8") as fid:
        fid.write(raw_settings)


class SchemaHandler(APIHandler):
    """Base handler for handler requiring access to settings."""

    def initialize(
        self,
        app_settings_dir: str,
        schemas_dir: str,
        settings_dir: str,
        labextensions_path: list[str] | None,
        overrides: dict[str, Any] | None = None,
        **kwargs: Any,
    ) -> None:
        """Initialize the handler."""
        super().initialize(**kwargs)
        error = None
        if not overrides:
            overrides, error = _get_overrides(app_settings_dir)
        self.overrides = overrides
        self.app_settings_dir = app_settings_dir
        self.schemas_dir = schemas_dir
        self.settings_dir = settings_dir
        self.labextensions_path = labextensions_path

        if error:
            overrides_warning = "Failed loading overrides: %s"
            self.log.warning(overrides_warning, error)

    def get_current_locale(self) -> str:
        """
        Get the current locale as specified in the translation-extension settings.

        Returns
        -------
        str
            The current locale string.

        Notes
        -----
        If the locale setting is not available or not valid, it will default to jupyterlab_server.translation_utils.DEFAULT_LOCALE.
        """
        try:
            settings, _ = get_settings(
                self.app_settings_dir,
                self.schemas_dir,
                self.settings_dir,
                schema_name=L10N_SCHEMA_NAME,
                overrides=self.overrides,
                labextensions_path=self.labextensions_path,
            )
        except web.HTTPError as e:
            schema_warning = "Missing or misshapen translation settings schema:\n%s"
            self.log.warning(schema_warning, e)

            settings = {}

        current_locale = settings.get("settings", {}).get("locale") or SYS_LOCALE
        if current_locale == "default":
            current_locale = SYS_LOCALE
        if not is_valid_locale(current_locale):
            current_locale = DEFAULT_LOCALE

        return current_locale
