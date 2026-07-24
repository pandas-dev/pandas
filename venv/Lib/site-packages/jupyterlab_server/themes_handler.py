"""Tornado handlers for dynamic theme loading."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import os
import re
from collections.abc import Generator
from glob import glob
from typing import Any
from urllib.parse import urlparse

from jupyter_server.base.handlers import FileFindHandler
from jupyter_server.utils import url_path_join as ujoin


class ThemesHandler(FileFindHandler):
    """A file handler that mangles local urls in CSS files."""

    def initialize(
        self,
        path: str | list[str],
        default_filename: str | None = None,
        no_cache_paths: list[str] | None = None,
        themes_url: str | None = None,
        labextensions_path: list[str] | None = None,
        **kwargs: Any,  # noqa: ARG002
    ) -> None:
        """Initialize the handler."""
        # Get all of the available theme paths in order
        labextensions_path = labextensions_path or []
        ext_paths: list[str] = []
        for ext_dir in labextensions_path:
            theme_pattern = ext_dir + "/**/themes"
            ext_paths.extend(path for path in glob(theme_pattern, recursive=True))

        # Add the core theme path last
        if not isinstance(path, list):
            path = [path]
        path = ext_paths + path

        FileFindHandler.initialize(
            self, path, default_filename=default_filename, no_cache_paths=no_cache_paths
        )
        self.themes_url = themes_url

    def get_content(  # type:ignore[override]
        self, abspath: str, start: int | None = None, end: int | None = None
    ) -> bytes | Generator[bytes, None, None]:
        """Retrieve the content of the requested resource which is located
        at the given absolute path.

        This method should either return a byte string or an iterator
        of byte strings.
        """
        base, ext = os.path.splitext(abspath)
        if ext != ".css":
            return FileFindHandler.get_content(abspath, start, end)

        return self._get_css()

    def get_content_size(self) -> int:
        """Retrieve the total size of the resource at the given path."""
        assert self.absolute_path is not None
        base, ext = os.path.splitext(self.absolute_path)
        if ext != ".css":
            return FileFindHandler.get_content_size(self)
        return len(self._get_css())

    def _get_css(self) -> bytes:
        """Get the mangled css file contents."""
        assert self.absolute_path is not None
        with open(self.absolute_path, "rb") as fid:
            data = fid.read().decode("utf-8")

        if not self.themes_url:
            return b""

        basedir = os.path.dirname(self.path).replace(os.sep, "/")
        basepath = ujoin(self.themes_url, basedir)

        # Replace local paths with mangled paths.
        # We only match strings that are local urls,
        # e.g. `url('../foo.css')`, `url('images/foo.png')`
        pattern = r"url\('(.*)'\)|url\('(.*)'\)"

        def replacer(m: Any) -> Any:
            """Replace the matched relative url with the mangled url."""
            group = m.group()
            # Get the part that matched
            part = next(g for g in m.groups() if g)

            # Ignore urls that start with `/` or have a protocol like `http`.
            parsed = urlparse(part)
            if part.startswith("/") or parsed.scheme:
                return group

            return group.replace(part, ujoin(basepath, part))

        return re.sub(pattern, replacer, data).encode("utf-8")
