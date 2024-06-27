"""API handlers for terminals."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from jupyter_server.auth.decorator import authorized
from jupyter_server.base.handlers import APIHandler
from tornado import web

from .base import TerminalsMixin

AUTH_RESOURCE = "terminals"


class TerminalAPIHandler(APIHandler):
    """The base terminal handler."""

    auth_resource = AUTH_RESOURCE


class TerminalRootHandler(TerminalsMixin, TerminalAPIHandler):
    """The root termanal API handler."""

    @web.authenticated
    @authorized
    def get(self) -> None:
        """Get the list of terminals."""
        models = self.terminal_manager.list()
        self.finish(json.dumps(models))

    @web.authenticated
    @authorized
    def post(self) -> None:
        """POST /terminals creates a new terminal and redirects to it"""
        data = self.get_json_body() or {}

        # if cwd is a relative path, it should be relative to the root_dir,
        # but if we pass it as relative, it will we be considered as relative to
        # the path jupyter_server was started in
        if "cwd" in data:
            cwd: Path | None = Path(data["cwd"])
            assert cwd is not None
            if not cwd.resolve().exists():
                cwd = Path(self.settings["server_root_dir"]).expanduser() / cwd
                if not cwd.resolve().exists():
                    cwd = None

            if cwd is None:
                server_root_dir = self.settings["server_root_dir"]
                self.log.debug(
                    "Failed to find requested terminal cwd: %s\n"
                    "  It was not found within the server root neither: %s.",
                    data.get("cwd"),
                    server_root_dir,
                )
                del data["cwd"]
            else:
                self.log.debug("Opening terminal in: %s", cwd.resolve())
                data["cwd"] = str(cwd.resolve())

        model = self.terminal_manager.create(**data)
        self.finish(json.dumps(model))


class TerminalHandler(TerminalsMixin, TerminalAPIHandler):
    """A handler for a specific terminal."""

    SUPPORTED_METHODS = ("GET", "DELETE", "OPTIONS")  # type:ignore[assignment]

    @web.authenticated
    @authorized
    def get(self, name: str) -> None:
        """Get a terminal by name."""
        model = self.terminal_manager.get(name)
        self.finish(json.dumps(model))

    @web.authenticated
    @authorized
    async def delete(self, name: str) -> None:
        """Remove a terminal by name."""
        await self.terminal_manager.terminate(name, force=True)
        self.set_status(204)
        self.finish()


default_handlers: list[tuple[str, type[Any]]] = [
    (r"/api/terminals", TerminalRootHandler),
    (r"/api/terminals/(\w+)", TerminalHandler),
]
