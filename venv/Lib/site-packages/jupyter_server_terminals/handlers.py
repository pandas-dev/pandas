"""Tornado handlers for the terminal emulator."""
# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import typing as t

from jupyter_core.utils import ensure_async
from jupyter_server._tz import utcnow
from jupyter_server.auth.utils import warn_disabled_authorization
from jupyter_server.base.handlers import JupyterHandler
from jupyter_server.base.websocket import WebSocketMixin
from terminado.management import NamedTermManager
from terminado.websocket import TermSocket as BaseTermSocket
from tornado import web

from .base import TerminalsMixin

AUTH_RESOURCE = "terminals"


class TermSocket(TerminalsMixin, WebSocketMixin, JupyterHandler, BaseTermSocket):
    """A terminal websocket."""

    auth_resource = AUTH_RESOURCE

    def initialize(  # type:ignore[override]
        self, name: str, term_manager: NamedTermManager, **kwargs: t.Any
    ) -> None:
        """Initialize the socket."""
        BaseTermSocket.initialize(self, term_manager, **kwargs)
        TerminalsMixin.initialize(self, name)

    def origin_check(self, origin: t.Any = None) -> bool:
        """Terminado adds redundant origin_check
        Tornado already calls check_origin, so don't do anything here.
        """
        return True

    async def get(self, *args: t.Any, **kwargs: t.Any) -> None:
        """Get the terminal socket."""
        user = self.current_user

        if not user:
            raise web.HTTPError(403)

        # authorize the user.
        if self.authorizer is None:
            # Warn if an authorizer is unavailable.
            warn_disabled_authorization()  # type:ignore[unreachable]
        elif not await ensure_async(
            self.authorizer.is_authorized(self, user, "execute", self.auth_resource)
        ):
            raise web.HTTPError(403)

        if args[0] not in self.term_manager.terminals:  # type:ignore[attr-defined]
            raise web.HTTPError(404)
        resp = super().get(*args, **kwargs)
        if resp is not None:
            await ensure_async(resp)  # type:ignore[arg-type]

    async def on_message(self, message: t.Any) -> None:  # type:ignore[override]
        """Handle a socket message."""
        await ensure_async(super().on_message(message))  # type:ignore[arg-type]
        self._update_activity()

    def write_message(self, message: t.Any, binary: bool = False) -> None:  # type:ignore[override]
        """Write a message to the socket."""
        super().write_message(message, binary=binary)
        self._update_activity()

    def _update_activity(self) -> None:
        self.application.settings["terminal_last_activity"] = utcnow()
        # terminal may not be around on deletion/cull
        if self.term_name in self.terminal_manager.terminals:
            self.terminal_manager.terminals[self.term_name].last_activity = utcnow()  # type:ignore[attr-defined]
