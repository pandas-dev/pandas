# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""A workspace management CLI"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path
from typing import Any

from jupyter_core.application import JupyterApp
from traitlets import Bool, Unicode

from ._version import __version__
from .config import LabConfig
from .workspaces_handler import WorkspacesManager

# Default workspace ID
#  Needs to match PageConfig.defaultWorkspace define in packages/coreutils/src/pageconfig.ts
DEFAULT_WORKSPACE = "default"


class WorkspaceListApp(JupyterApp, LabConfig):
    """An app to list workspaces."""

    version = __version__
    description = """
    Print all the workspaces available

    If '--json' flag is passed in, a single 'json' object is printed.
    If '--jsonlines' flag is passed in, 'json' object of each workspace separated by a new line is printed.
    If nothing is passed in, workspace ids list is printed.
    """
    flags = dict(
        jsonlines=(
            {"WorkspaceListApp": {"jsonlines": True}},
            ("Produce machine-readable JSON Lines output."),
        ),
        json=(
            {"WorkspaceListApp": {"json": True}},
            ("Produce machine-readable JSON object output."),
        ),
    )

    jsonlines = Bool(
        False,
        config=True,
        help=(
            "If True, the output will be a newline-delimited JSON (see https://jsonlines.org/) of objects, "
            "one per JupyterLab workspace, each with the details of the relevant workspace"
        ),
    )
    json = Bool(
        False,
        config=True,
        help=(
            "If True, each line of output will be a JSON object with the "
            "details of the workspace."
        ),
    )

    def initialize(self, *args: Any, **kwargs: Any) -> None:
        """Initialize the app."""
        super().initialize(*args, **kwargs)
        self.manager = WorkspacesManager(self.workspaces_dir)

    def start(self) -> None:
        """Start the app."""
        workspaces = self.manager.list_workspaces()
        if self.jsonlines:
            for workspace in workspaces:
                print(json.dumps(workspace))
        elif self.json:
            print(json.dumps(workspaces))
        else:
            for workspace in workspaces:
                print(workspace["metadata"]["id"])


class WorkspaceExportApp(JupyterApp, LabConfig):
    """A workspace export app."""

    version = __version__
    description = """
    Export a JupyterLab workspace

    If no arguments are passed in, this command will export the default
        workspace.
    If a workspace name is passed in, this command will export that workspace.
    If no workspace is found, this command will export an empty workspace.
    """

    def initialize(self, *args: Any, **kwargs: Any) -> None:
        """Initialize the app."""
        super().initialize(*args, **kwargs)
        self.manager = WorkspacesManager(self.workspaces_dir)

    def start(self) -> None:
        """Start the app."""
        if len(self.extra_args) > 1:  # pragma: no cover
            warnings.warn("Too many arguments were provided for workspace export.")
            self.exit(1)

        raw = DEFAULT_WORKSPACE if not self.extra_args else self.extra_args[0]
        try:
            workspace = self.manager.load(raw)
            print(json.dumps(workspace))
        except Exception:  # pragma: no cover
            self.log.error(json.dumps(dict(data=dict(), metadata=dict(id=raw))))


class WorkspaceImportApp(JupyterApp, LabConfig):
    """A workspace import app."""

    version = __version__
    description = """
    Import a JupyterLab workspace

    This command will import a workspace from a JSON file. The format of the
        file must be the same as what the export functionality emits.
    """
    workspace_name = Unicode(
        None,
        config=True,
        allow_none=True,
        help="""
        Workspace name. If given, the workspace ID in the imported
        file will be replaced with a new ID pointing to this
        workspace name.
        """,
    )

    aliases = {"name": "WorkspaceImportApp.workspace_name"}

    def initialize(self, *args: Any, **kwargs: Any) -> None:
        """Initialize the app."""
        super().initialize(*args, **kwargs)
        self.manager = WorkspacesManager(self.workspaces_dir)

    def start(self) -> None:
        """Start the app."""
        if len(self.extra_args) != 1:  # pragma: no cover
            self.log.info("One argument is required for workspace import.")
            self.exit(1)

        with self._smart_open() as fid:
            try:  # to load, parse, and validate the workspace file.
                workspace = self._validate(fid)
            except Exception as e:  # pragma: no cover
                self.log.info("%s is not a valid workspace:\n%s", fid.name, e)
                self.exit(1)

        try:
            workspace_path = self.manager.save(workspace["metadata"]["id"], json.dumps(workspace))
        except Exception as e:  # pragma: no cover
            self.log.info("Workspace could not be exported:\n%s", e)
            self.exit(1)

        self.log.info("Saved workspace: %s", workspace_path)

    def _smart_open(self) -> Any:
        file_name = self.extra_args[0]

        if file_name == "-":  # pragma: no cover
            return sys.stdin

        file_path = Path(file_name).resolve()

        if not file_path.exists():  # pragma: no cover
            self.log.info("%s does not exist.", file_name)
            self.exit(1)

        return file_path.open(encoding="utf-8")

    def _validate(self, data: Any) -> Any:
        workspace = json.load(data)

        if "data" not in workspace:
            msg = "The `data` field is missing."
            raise Exception(msg)

        # If workspace_name is set in config, inject the
        # name into the workspace metadata.
        if self.workspace_name is not None and self.workspace_name:
            workspace["metadata"] = {"id": self.workspace_name}
        elif "id" not in workspace["metadata"]:
            msg = "The `id` field is missing in `metadata`."
            raise Exception(msg)

        return workspace
