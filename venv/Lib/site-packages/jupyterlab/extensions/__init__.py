"""Extension manager for JupyterLab."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from importlib.metadata import entry_points

from traitlets.config import Configurable

from .manager import ActionResult, ExtensionManager, ExtensionPackage  # noqa: F401
from .pypi import PyPIExtensionManager
from .readonly import ReadOnlyExtensionManager

# Supported third-party services
MANAGERS = {}

for entry in entry_points(group="jupyterlab.extension_manager_v1"):
    MANAGERS[entry.name] = entry


# Entry points


def get_readonly_manager(
    app_options: dict | None = None,
    ext_options: dict | None = None,
    parent: Configurable | None = None,
) -> ExtensionManager:
    """Read-Only Extension Manager factory"""
    return ReadOnlyExtensionManager(app_options, ext_options, parent)


def get_pypi_manager(
    app_options: dict | None = None,
    ext_options: dict | None = None,
    parent: Configurable | None = None,
) -> ExtensionManager:
    """PyPi Extension Manager factory"""
    return PyPIExtensionManager(app_options, ext_options, parent)
