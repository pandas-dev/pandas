"""A Tornado UI module for a terminal backed by terminado.

See the Tornado docs for information on UI modules:
http://www.tornadoweb.org/en/stable/guide/templates.html#ui-modules
"""
# Copyright (c) Jupyter Development Team
# Copyright (c) 2014, Ramalingam Saravanan <sarava@sarava.net>
# Distributed under the terms of the Simplified BSD License.
from __future__ import annotations

from pathlib import Path

import tornado.web


class Terminal(tornado.web.UIModule):
    """A terminal UI module."""

    def render(self, ws_url: str, cols: int = 80, rows: int = 25) -> str:
        """Render the module."""
        return (
            '<div class="terminado-container" '
            f'data-ws-url="{ws_url}" '
            f'data-rows="{rows}" data-cols="{cols}"/>'
        )

    def javascript_files(self) -> list[str]:
        """Get the list of JS files to include."""
        # TODO: Can we calculate these dynamically?
        return ["/xstatic/termjs/term.js", "/static/terminado.js"]

    def embedded_javascript(self) -> str:
        """Get the embedded JS content as a string."""
        file = Path(__file__).parent / "uimod_embed.js"
        with file.open() as f:
            return f.read()
