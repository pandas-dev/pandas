"""Jupyter LabExtension Entry Points."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from __future__ import annotations

from copy import copy
from pathlib import Path

from jupyter_core.application import JupyterApp, base_aliases, base_flags
from jupyter_core.paths import jupyter_path
from traitlets import List, Unicode, default

from .debug_log_file_mixin import DebugLogFileMixin

HERE = str(Path(__file__).resolve().parent)

flags = dict(base_flags)

develop_flags = copy(flags)
develop_flags["overwrite"] = (
    {"DevelopLabExtensionApp": {"overwrite": True}},
    "Overwrite files",
)

aliases = dict(base_aliases)
aliases["debug-log-path"] = "DebugLogFileMixin.debug_log_path"

VERSION = "1"


class BaseExtensionApp(JupyterApp, DebugLogFileMixin):
    """Base application class for JupyterLab extension CLI commands."""

    version = VERSION
    flags = flags
    aliases = aliases
    name = "lab"

    labextensions_path = List(
        Unicode(),
        help="The standard paths to look in for prebuilt JupyterLab extensions",
    )

    @default("labextensions_path")
    def _default_labextensions_path(self) -> list[str]:
        return jupyter_path("labextensions")

    def start(self) -> None:
        """Start the extension app and run the configured task."""
        with self.debug_logging():
            self.run_task()

    def run_task(self) -> None:
        """Execute the app's primary task; override in subclasses."""

    def _log_format_default(self) -> str:
        """Return the default log format string."""
        return "%(message)s"
