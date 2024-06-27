"""Extension manager for JupyterLab."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import sys
from typing import Optional

from traitlets.config import Configurable

from .manager import ActionResult, ExtensionManager, ExtensionPackage  # noqa: F401
from .pypi import PyPIExtensionManager
from .readonly import ReadOnlyExtensionManager

# See compatibility note on `group` keyword in https://docs.python.org/3/library/importlib.metadata.html#entry-points
if sys.version_info < (3, 10):
    from importlib_metadata import entry_points
else:
    from importlib.metadata import entry_points

# Supported third-party services
MANAGERS = {}

for entry in entry_points(group="jupyterlab.extension_manager_v1"):
    MANAGERS[entry.name] = entry


# Entry points


def get_readonly_manager(
    app_options: Optional[dict] = None,
    ext_options: Optional[dict] = None,
    parent: Optional[Configurable] = None,
) -> ExtensionManager:
    """Read-Only Extension Manager factory"""
    return ReadOnlyExtensionManager(app_options, ext_options, parent)


def get_pypi_manager(
    app_options: Optional[dict] = None,
    ext_options: Optional[dict] = None,
    parent: Optional[Configurable] = None,
) -> ExtensionManager:
    """PyPi Extension Manager factory"""
    return PyPIExtensionManager(app_options, ext_options, parent)
