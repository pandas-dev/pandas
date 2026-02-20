"""Serve files directly from the ContentsManager."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import mimetypes
from base64 import decodebytes
from typing import TYPE_CHECKING

from jupyter_core.utils import ensure_async
from tornado import web

from jupyter_server.auth.decorator import authorized
from jupyter_server.base.handlers import JupyterHandler

if TYPE_CHECKING:
    from collections.abc import Awaitable

AUTH_RESOURCE = "contents"


class FilesHandler(JupyterHandler, web.StaticFileHandler):
    """serve files via ContentsManager

    Normally used when ContentsManager is not a FileContentsManager.

    FileContentsManager subclasses use AuthenticatedFilesHandler by default,
    a subclass of StaticFileHandler.
    """

    auth_resource = AUTH_RESOURCE

    @property
    def content_security_policy(self):
        """The content security policy."""
        # In case we're serving HTML/SVG, confine any Javascript to a unique
        # origin so it can't interact with the notebook server.
        return super().content_security_policy + "; sandbox allow-scripts"

    @web.authenticated
    @authorized
    def head(self, path: str) -> Awaitable[None] | None:  # type:ignore[override]
        """The head response."""
        self.get(path, include_body=False)
        self.check_xsrf_cookie()
        return self.get(path, include_body=False)

    @web.authenticated
    @authorized
    async def get(self, path, include_body=True):
        """Get a file by path."""
        # /files/ requests must originate from the same site
        self.check_xsrf_cookie()
        cm = self.contents_manager

        if not cm.allow_hidden and await ensure_async(cm.is_hidden(path)):
            self.log.info("Refusing to serve hidden file, via 404 Error")
            raise web.HTTPError(404)

        path = path.strip("/")
        if "/" in path:
            _, name = path.rsplit("/", 1)
        else:
            name = path

        model = await ensure_async(cm.get(path, type="file", content=include_body))

        if self.get_argument("download", None):
            self.set_attachment_header(name)

        # get mimetype from filename
        if name.lower().endswith(".ipynb"):
            self.set_header("Content-Type", "application/x-ipynb+json")
        else:
            cur_mime, encoding = mimetypes.guess_type(name)
            if cur_mime == "text/plain":
                self.set_header("Content-Type", "text/plain; charset=UTF-8")
            # RFC 6713
            if encoding == "gzip":
                self.set_header("Content-Type", "application/gzip")
            elif encoding is not None:
                self.set_header("Content-Type", "application/octet-stream")
            elif cur_mime is not None:
                self.set_header("Content-Type", cur_mime)
            elif model["format"] == "base64":
                self.set_header("Content-Type", "application/octet-stream")
            else:
                self.set_header("Content-Type", "text/plain; charset=UTF-8")

        if include_body:
            if model["format"] == "base64":
                b64_bytes = model["content"].encode("ascii")
                self.write(decodebytes(b64_bytes))
            else:
                self.write(model["content"])
            self.flush()


default_handlers: list[JupyterHandler] = []
