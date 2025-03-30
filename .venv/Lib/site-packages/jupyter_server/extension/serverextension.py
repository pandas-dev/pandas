"""Utilities for installing extensions"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import logging
import os
import sys
import typing as t

from jupyter_core.application import JupyterApp
from jupyter_core.paths import ENV_CONFIG_PATH, SYSTEM_CONFIG_PATH, jupyter_config_dir
from tornado.log import LogFormatter
from traitlets import Bool

from jupyter_server._version import __version__
from jupyter_server.extension.config import ExtensionConfigManager
from jupyter_server.extension.manager import ExtensionManager, ExtensionPackage


def _get_config_dir(user: bool = False, sys_prefix: bool = False) -> str:
    """Get the location of config files for the current context

    Returns the string to the environment

    Parameters
    ----------
    user : bool [default: False]
        Get the user's .jupyter config directory
    sys_prefix : bool [default: False]
        Get sys.prefix, i.e. ~/.envs/my-env/etc/jupyter
    """
    if user and sys_prefix:
        sys_prefix = False
    if user:
        extdir = jupyter_config_dir()
    elif sys_prefix:
        extdir = ENV_CONFIG_PATH[0]
    else:
        extdir = SYSTEM_CONFIG_PATH[0]
    return extdir


def _get_extmanager_for_context(
    write_dir: str = "jupyter_server_config.d", user: bool = False, sys_prefix: bool = False
) -> tuple[str, ExtensionManager]:
    """Get an extension manager pointing at the current context

    Returns the path to the current context and an ExtensionManager object.

    Parameters
    ----------
    write_dir : str [default: 'jupyter_server_config.d']
        Name of config directory to write extension config.
    user : bool [default: False]
        Get the user's .jupyter config directory
    sys_prefix : bool [default: False]
        Get sys.prefix, i.e. ~/.envs/my-env/etc/jupyter
    """
    config_dir = _get_config_dir(user=user, sys_prefix=sys_prefix)
    config_manager = ExtensionConfigManager(
        read_config_path=[config_dir],
        write_config_dir=os.path.join(config_dir, write_dir),
    )
    extension_manager = ExtensionManager(
        config_manager=config_manager,
    )
    return config_dir, extension_manager


class ArgumentConflict(ValueError):
    pass


_base_flags: dict[str, t.Any] = {}
_base_flags.update(JupyterApp.flags)
_base_flags.pop("y", None)
_base_flags.pop("generate-config", None)
_base_flags.update(
    {
        "user": (
            {
                "BaseExtensionApp": {
                    "user": True,
                }
            },
            "Apply the operation only for the given user",
        ),
        "system": (
            {
                "BaseExtensionApp": {
                    "user": False,
                    "sys_prefix": False,
                }
            },
            "Apply the operation system-wide",
        ),
        "sys-prefix": (
            {
                "BaseExtensionApp": {
                    "sys_prefix": True,
                }
            },
            "Use sys.prefix as the prefix for installing extensions (for environments, packaging)",
        ),
        "py": (
            {
                "BaseExtensionApp": {
                    "python": True,
                }
            },
            "Install from a Python package",
        ),
    }
)
_base_flags["python"] = _base_flags["py"]

_base_aliases: dict[str, t.Any] = {}
_base_aliases.update(JupyterApp.aliases)


class BaseExtensionApp(JupyterApp):
    """Base extension installer app"""

    _log_formatter_cls = LogFormatter  # type:ignore[assignment]
    flags = _base_flags
    aliases = _base_aliases
    version = __version__

    user = Bool(False, config=True, help="Whether to do a user install")
    sys_prefix = Bool(True, config=True, help="Use the sys.prefix as the prefix")
    python = Bool(False, config=True, help="Install from a Python package")

    def _log_format_default(self) -> str:
        """A default format for messages"""
        return "%(message)s"

    @property
    def config_dir(self) -> str:  # type:ignore[override]
        return _get_config_dir(user=self.user, sys_prefix=self.sys_prefix)


# Constants for pretty print extension listing function.
# Window doesn't support coloring in the commandline
GREEN_ENABLED = "\033[32menabled\033[0m" if os.name != "nt" else "enabled"
RED_DISABLED = "\033[31mdisabled\033[0m" if os.name != "nt" else "disabled"
GREEN_OK = "\033[32mOK\033[0m" if os.name != "nt" else "ok"
RED_X = "\033[31m X\033[0m" if os.name != "nt" else " X"

# ------------------------------------------------------------------------------
# Public API
# ------------------------------------------------------------------------------


def toggle_server_extension_python(
    import_name: str,
    enabled: bool | None = None,
    parent: t.Any = None,
    user: bool = False,
    sys_prefix: bool = True,
) -> None:
    """Toggle the boolean setting for a given server extension
    in a Jupyter config file.
    """
    sys_prefix = False if user else sys_prefix
    config_dir = _get_config_dir(user=user, sys_prefix=sys_prefix)
    manager = ExtensionConfigManager(
        read_config_path=[config_dir],
        write_config_dir=os.path.join(config_dir, "jupyter_server_config.d"),
    )
    if enabled:
        manager.enable(import_name)
    else:
        manager.disable(import_name)


# ----------------------------------------------------------------------
# Applications
# ----------------------------------------------------------------------

flags = {}
flags.update(BaseExtensionApp.flags)
flags.pop("y", None)
flags.pop("generate-config", None)
flags.update(
    {
        "user": (
            {
                "ToggleServerExtensionApp": {
                    "user": True,
                }
            },
            "Perform the operation for the current user",
        ),
        "system": (
            {
                "ToggleServerExtensionApp": {
                    "user": False,
                    "sys_prefix": False,
                }
            },
            "Perform the operation system-wide",
        ),
        "sys-prefix": (
            {
                "ToggleServerExtensionApp": {
                    "sys_prefix": True,
                }
            },
            "Use sys.prefix as the prefix for installing server extensions",
        ),
        "py": (
            {
                "ToggleServerExtensionApp": {
                    "python": True,
                }
            },
            "Install from a Python package",
        ),
    }
)
flags["python"] = flags["py"]


_desc = "Enable/disable a server extension using frontend configuration files."


class ToggleServerExtensionApp(BaseExtensionApp):
    """A base class for enabling/disabling extensions"""

    name = "jupyter server extension enable/disable"
    description = _desc

    flags = flags

    _toggle_value = Bool()
    _toggle_pre_message = ""
    _toggle_post_message = ""

    def toggle_server_extension(self, import_name: str) -> None:
        """Change the status of a named server extension.

        Uses the value of `self._toggle_value`.

        Parameters
        ---------

        import_name : str
            Importable Python module (dotted-notation) exposing the magic-named
            `load_jupyter_server_extension` function
        """
        # Create an extension manager for this instance.
        config_dir, extension_manager = _get_extmanager_for_context(
            user=self.user, sys_prefix=self.sys_prefix
        )
        try:
            self.log.info(f"{self._toggle_pre_message.capitalize()}: {import_name}")
            self.log.info(f"- Writing config: {config_dir}")
            # Validate the server extension.
            self.log.info(f"    - Validating {import_name}...")
            # Interface with the Extension Package and validate.
            extpkg = ExtensionPackage(name=import_name)
            extpkg.validate()
            version = extpkg.version
            self.log.info(f"      {import_name} {version} {GREEN_OK}")

            # Toggle extension config.
            config = extension_manager.config_manager
            if config:
                if self._toggle_value is True:
                    config.enable(import_name)
                else:
                    config.disable(import_name)

            # If successful, let's log.
            self.log.info(f"    - Extension successfully {self._toggle_post_message}.")
        except Exception as err:
            self.log.error(f"     {RED_X} Validation failed: {err}")

    def start(self) -> None:
        """Perform the App's actions as configured"""
        if not self.extra_args:
            sys.exit("Please specify a server extension/package to enable or disable")
        for arg in self.extra_args:
            self.toggle_server_extension(arg)


class EnableServerExtensionApp(ToggleServerExtensionApp):
    """An App that enables (and validates) Server Extensions"""

    name = "jupyter server extension enable"
    description = """
    Enable a server extension in configuration.

    Usage
        jupyter server extension enable [--system|--sys-prefix]
    """
    _toggle_value = True  # type:ignore[assignment]
    _toggle_pre_message = "enabling"
    _toggle_post_message = "enabled"


class DisableServerExtensionApp(ToggleServerExtensionApp):
    """An App that disables Server Extensions"""

    name = "jupyter server extension disable"
    description = """
    Disable a server extension in configuration.

    Usage
        jupyter server extension disable [--system|--sys-prefix]
    """
    _toggle_value = False  # type:ignore[assignment]
    _toggle_pre_message = "disabling"
    _toggle_post_message = "disabled"


class ListServerExtensionsApp(BaseExtensionApp):
    """An App that lists (and validates) Server Extensions"""

    name = "jupyter server extension list"
    version = __version__
    description = "List all server extensions known by the configuration system"

    def list_server_extensions(self) -> None:
        """List all enabled and disabled server extensions, by config path

        Enabled extensions are validated, potentially generating warnings.
        """
        configurations = (
            {"user": True, "sys_prefix": False},
            {"user": False, "sys_prefix": True},
            {"user": False, "sys_prefix": False},
        )

        for option in configurations:
            config_dir = _get_config_dir(**option)
            print(f"Config dir: {config_dir}")
            write_dir = "jupyter_server_config.d"
            config_manager = ExtensionConfigManager(
                read_config_path=[config_dir],
                write_config_dir=os.path.join(config_dir, write_dir),
            )
            jpserver_extensions = config_manager.get_jpserver_extensions()
            for name, enabled in jpserver_extensions.items():
                # Attempt to get extension metadata
                print(f"    {name} {GREEN_ENABLED if enabled else RED_DISABLED}")
                try:
                    print(f"    - Validating {name}...")
                    extension = ExtensionPackage(name=name, enabled=enabled)
                    if not extension.validate():
                        msg = "validation failed"
                        raise ValueError(msg)
                    version = extension.version
                    print(f"      {name} {version} {GREEN_OK}")
                except Exception as err:
                    self.log.debug("", exc_info=True)
                    print(f"      {RED_X} {err}")
            # Add a blank line between paths.
            self.log.info("")

    def start(self) -> None:
        """Perform the App's actions as configured"""
        self.list_server_extensions()


_examples = """
jupyter server extension list                        # list all configured server extensions
jupyter server extension enable --py <packagename>   # enable all server extensions in a Python package
jupyter server extension disable --py <packagename>  # disable all server extensions in a Python package
"""


class ServerExtensionApp(BaseExtensionApp):
    """Root level server extension app"""

    name = "jupyter server extension"
    version = __version__
    description: str = "Work with Jupyter server extensions"
    examples = _examples

    subcommands: dict[str, t.Any] = {
        "enable": (EnableServerExtensionApp, "Enable a server extension"),
        "disable": (DisableServerExtensionApp, "Disable a server extension"),
        "list": (ListServerExtensionsApp, "List server extensions"),
    }

    def start(self) -> None:
        """Perform the App's actions as configured"""
        super().start()

        # The above should have called a subcommand and raised NoStart; if we
        # get here, it didn't, so we should self.log.info a message.
        subcmds = ", ".join(sorted(self.subcommands))
        sys.exit("Please supply at least one subcommand: %s" % subcmds)


main = ServerExtensionApp.launch_instance


if __name__ == "__main__":
    main()
