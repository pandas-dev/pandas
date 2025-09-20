"""Tornado handlers for frontend config storage."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import hashlib
import json
import re
import unicodedata
import urllib
from pathlib import Path
from typing import Any

from jupyter_server import _tz as tz
from jupyter_server.base.handlers import APIHandler
from jupyter_server.extension.handler import ExtensionHandlerJinjaMixin, ExtensionHandlerMixin
from jupyter_server.utils import url_path_join as ujoin
from tornado import web
from traitlets.config import LoggingConfigurable

# The JupyterLab workspace file extension.
WORKSPACE_EXTENSION = ".jupyterlab-workspace"


def _list_workspaces(directory: Path, prefix: str) -> list[dict[str, Any]]:
    """
    Return the list of workspaces in a given directory beginning with the
    given prefix.
    """
    workspaces: list = []
    if not directory.exists():
        return workspaces

    items = [
        item
        for item in directory.iterdir()
        if item.name.startswith(prefix) and item.name.endswith(WORKSPACE_EXTENSION)
    ]
    items.sort()

    for slug in items:
        workspace_path: Path = directory / slug
        if workspace_path.exists():
            workspace = _load_with_file_times(workspace_path)
            workspaces.append(workspace)

    return workspaces


def _load_with_file_times(workspace_path: Path) -> dict:
    """
    Load workspace JSON from disk, overwriting the `created` and `last_modified`
    metadata with current file stat information
    """
    stat = workspace_path.stat()
    with workspace_path.open(encoding="utf-8") as fid:
        workspace = json.load(fid)
        workspace["metadata"].update(
            last_modified=tz.utcfromtimestamp(stat.st_mtime).isoformat(),
            created=tz.utcfromtimestamp(stat.st_ctime).isoformat(),
        )
    return workspace


def slugify(
    raw: str, base: str = "", sign: bool = True, max_length: int = 128 - len(WORKSPACE_EXTENSION)
) -> str:
    """
    Use the common superset of raw and base values to build a slug shorter
    than max_length. By default, base value is an empty string.
    Convert spaces to hyphens. Remove characters that aren't alphanumerics
    underscores, or hyphens. Convert to lowercase. Strip leading and trailing
    whitespace.
    Add an optional short signature suffix to prevent collisions.
    Modified from Django utils:
    https://github.com/django/django/blob/master/django/utils/text.py
    """
    raw = raw if raw.startswith("/") else "/" + raw
    signature = ""
    if sign:
        data = raw[1:]  # Remove initial slash that always exists for digest.
        signature = "-" + hashlib.sha256(data.encode("utf-8")).hexdigest()[:4]
    base = (base if base.startswith("/") else "/" + base).lower()
    raw = raw.lower()
    common = 0
    limit = min(len(base), len(raw))
    while common < limit and base[common] == raw[common]:
        common += 1
    value = ujoin(base[common:], raw)
    value = urllib.parse.unquote(value)
    value = unicodedata.normalize("NFKC", value).encode("ascii", "ignore").decode("ascii")
    value = re.sub(r"[^\w\s-]", "", value).strip()
    value = re.sub(r"[-\s]+", "-", value)
    return value[: max_length - len(signature)] + signature


class WorkspacesManager(LoggingConfigurable):
    """A manager for workspaces."""

    def __init__(self, path: str) -> None:
        """Initialize a workspaces manager with content in ``path``."""
        super()
        if not path:
            msg = "Workspaces directory is not set"
            raise ValueError(msg)
        self.workspaces_dir = Path(path)

    def delete(self, space_name: str) -> None:
        """Remove a workspace ``space_name``."""
        slug = slugify(space_name)
        workspace_path = self.workspaces_dir / (slug + WORKSPACE_EXTENSION)

        if not workspace_path.exists():
            msg = f"Workspace {space_name!r} ({slug!r}) not found"
            raise FileNotFoundError(msg)

        # to delete the workspace file.
        workspace_path.unlink()

    def list_workspaces(self) -> list:
        """List all available workspaces."""
        prefix = slugify("", sign=False)
        return _list_workspaces(self.workspaces_dir, prefix)

    def load(self, space_name: str) -> dict:
        """Load the workspace ``space_name``."""
        slug = slugify(space_name)
        workspace_path = self.workspaces_dir / (slug + WORKSPACE_EXTENSION)

        if workspace_path.exists():
            # to load and parse the workspace file.
            return _load_with_file_times(workspace_path)
        _id = space_name if space_name.startswith("/") else "/" + space_name
        return dict(data=dict(), metadata=dict(id=_id))

    def save(self, space_name: str, raw: str) -> Path:
        """Save the ``raw`` data as workspace ``space_name``."""
        if not self.workspaces_dir.exists():
            self.workspaces_dir.mkdir(parents=True)

        workspace = {}

        # Make sure the data is valid JSON.
        try:
            decoder = json.JSONDecoder()
            workspace = decoder.decode(raw)
        except Exception as e:
            raise ValueError(str(e)) from e

        # Make sure metadata ID matches the workspace name.
        # Transparently support an optional initial root `/`.
        metadata_id = workspace["metadata"]["id"]
        metadata_id = metadata_id if metadata_id.startswith("/") else "/" + metadata_id
        metadata_id = urllib.parse.unquote(metadata_id)
        if metadata_id != "/" + space_name:
            message = f"Workspace metadata ID mismatch: expected {space_name!r} got {metadata_id!r}"
            raise ValueError(message)

        slug = slugify(space_name)
        workspace_path = self.workspaces_dir / (slug + WORKSPACE_EXTENSION)

        # Write the workspace data to a file.
        workspace_path.write_text(raw, encoding="utf-8")

        return workspace_path


class WorkspacesHandler(ExtensionHandlerMixin, ExtensionHandlerJinjaMixin, APIHandler):
    """A workspaces API handler."""

    def initialize(self, name: str, manager: WorkspacesManager, **kwargs: Any) -> None:  # noqa: ARG002
        """Initialize the handler."""
        super().initialize(name)
        self.manager = manager

    @web.authenticated
    def delete(self, space_name: str) -> None:
        """Remove a workspace"""
        if not space_name:
            raise web.HTTPError(400, "Workspace name is required for DELETE")

        try:
            self.manager.delete(space_name)
            return self.set_status(204)
        except FileNotFoundError as e:
            raise web.HTTPError(404, str(e)) from e
        except Exception as e:  # pragma: no cover
            raise web.HTTPError(500, str(e)) from e

    @web.authenticated
    async def get(self, space_name: str = "") -> Any:
        """Get workspace(s) data"""

        try:
            if not space_name:
                workspaces = self.manager.list_workspaces()
                ids = []
                values = []
                for workspace in workspaces:
                    ids.append(workspace["metadata"]["id"])
                    values.append(workspace)
                return self.finish(json.dumps({"workspaces": {"ids": ids, "values": values}}))

            workspace = self.manager.load(space_name)
            return self.finish(json.dumps(workspace))
        except Exception as e:  # pragma: no cover
            raise web.HTTPError(500, str(e)) from e

    @web.authenticated
    def put(self, space_name: str = "") -> None:
        """Update workspace data"""
        if not space_name:
            raise web.HTTPError(400, "Workspace name is required for PUT.")

        raw = self.request.body.strip().decode("utf-8")

        # Make sure the data is valid JSON.
        try:
            self.manager.save(space_name, raw)
        except ValueError as e:
            raise web.HTTPError(400, str(e)) from e
        except Exception as e:  # pragma: no cover
            raise web.HTTPError(500, str(e)) from e

        self.set_status(204)
