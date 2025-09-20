"""Tornado handlers for the contents web service.

Preliminary documentation at https://github.com/ipython/ipython/wiki/IPEP-27%3A-Contents-Service
"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
import json
from http import HTTPStatus
from typing import Any

try:
    from jupyter_client.jsonutil import json_default
except ImportError:
    from jupyter_client.jsonutil import date_default as json_default

from jupyter_core.utils import ensure_async
from tornado import web

from jupyter_server.auth.decorator import allow_unauthenticated, authorized
from jupyter_server.base.handlers import APIHandler, JupyterHandler, path_regex
from jupyter_server.utils import url_escape, url_path_join

AUTH_RESOURCE = "contents"


def _validate_keys(expect_defined: bool, model: dict[str, Any], keys: list[str]):
    """
    Validate that the keys are defined (i.e. not None) or not (i.e. None)
    """

    if expect_defined:
        errors = [key for key in keys if model[key] is None]
        if errors:
            raise web.HTTPError(
                500,
                f"Keys unexpectedly None: {errors}",
            )
    else:
        errors = {key: model[key] for key in keys if model[key] is not None}  # type: ignore[assignment]
        if errors:
            raise web.HTTPError(
                500,
                f"Keys unexpectedly not None: {errors}",
            )


def validate_model(model, expect_content=False, expect_hash=False):
    """
    Validate a model returned by a ContentsManager method.

    If expect_content is True, then we expect non-null entries for 'content'
    and 'format'.

    If expect_hash is True, then we expect non-null entries for 'hash' and 'hash_algorithm'.
    """
    required_keys = {
        "name",
        "path",
        "type",
        "writable",
        "created",
        "last_modified",
        "mimetype",
        "content",
        "format",
    }
    if expect_hash:
        required_keys.update(["hash", "hash_algorithm"])
    missing = required_keys - set(model.keys())
    if missing:
        raise web.HTTPError(
            500,
            f"Missing Model Keys: {missing}",
        )

    content_keys = ["content", "format"]
    _validate_keys(expect_content, model, content_keys)
    if expect_hash:
        _validate_keys(expect_hash, model, ["hash", "hash_algorithm"])


class ContentsAPIHandler(APIHandler):
    """A contents API handler."""

    auth_resource = AUTH_RESOURCE


class ContentsHandler(ContentsAPIHandler):
    """A contents handler."""

    def location_url(self, path):
        """Return the full URL location of a file.

        Parameters
        ----------
        path : unicode
            The API path of the file, such as "foo/bar.txt".
        """
        return url_path_join(self.base_url, "api", "contents", url_escape(path))

    def _finish_model(self, model, location=True):
        """Finish a JSON request with a model, setting relevant headers, etc."""
        if location:
            location = self.location_url(model["path"])
            self.set_header("Location", location)
        self.set_header("Last-Modified", model["last_modified"])
        self.set_header("Content-Type", "application/json")
        self.finish(json.dumps(model, default=json_default))

    async def _finish_error(self, code, message):
        """Finish a JSON request with an error code and descriptive message"""
        self.set_status(code)
        self.write(message)
        await self.finish()

    @web.authenticated
    @authorized
    async def get(self, path=""):
        """Return a model for a file or directory.

        A directory model contains a list of models (without content)
        of the files and directories it contains.
        """
        path = path or ""
        cm = self.contents_manager

        type = self.get_query_argument("type", default=None)
        if type not in {None, "directory", "file", "notebook"}:
            # fall back to file if unknown type
            type = "file"

        format = self.get_query_argument("format", default=None)
        if format not in {None, "text", "base64"}:
            raise web.HTTPError(400, "Format %r is invalid" % format)
        content_str = self.get_query_argument("content", default="1")
        if content_str not in {"0", "1"}:
            raise web.HTTPError(400, "Content %r is invalid" % content_str)
        content = int(content_str or "")

        hash_str = self.get_query_argument("hash", default="0")
        if hash_str not in {"0", "1"}:
            raise web.HTTPError(
                400, f"Hash argument {hash_str!r} is invalid. It must be '0' or '1'."
            )
        require_hash = int(hash_str)

        if not cm.allow_hidden and await ensure_async(cm.is_hidden(path)):
            await self._finish_error(
                HTTPStatus.NOT_FOUND, f"file or directory {path!r} does not exist"
            )

        try:
            expect_hash = require_hash
            try:
                model = await ensure_async(
                    self.contents_manager.get(
                        path=path,
                        type=type,
                        format=format,
                        content=content,
                        require_hash=require_hash,
                    )
                )
            except TypeError:
                # Fallback for ContentsManager not handling the require_hash argument
                # introduced in 2.11
                expect_hash = False
                model = await ensure_async(
                    self.contents_manager.get(
                        path=path,
                        type=type,
                        format=format,
                        content=content,
                    )
                )
            validate_model(model, expect_content=content, expect_hash=expect_hash)
            self._finish_model(model, location=False)
        except web.HTTPError as exc:
            # 404 is okay in this context, catch exception and return 404 code to prevent stack trace on client
            if exc.status_code == HTTPStatus.NOT_FOUND:
                await self._finish_error(
                    HTTPStatus.NOT_FOUND, f"file or directory {path!r} does not exist"
                )
            raise

    @web.authenticated
    @authorized
    async def patch(self, path=""):
        """PATCH renames a file or directory without re-uploading content."""
        cm = self.contents_manager
        model = self.get_json_body()
        if model is None:
            raise web.HTTPError(400, "JSON body missing")

        old_path = model.get("path")
        if (
            old_path
            and not cm.allow_hidden
            and (
                await ensure_async(cm.is_hidden(path)) or await ensure_async(cm.is_hidden(old_path))
            )
        ):
            raise web.HTTPError(400, f"Cannot rename file or directory {path!r}")

        model = await ensure_async(cm.update(model, path))
        validate_model(model)
        self._finish_model(model)

    async def _copy(self, copy_from, copy_to=None):
        """Copy a file, optionally specifying a target directory."""
        self.log.info(
            "Copying %r to %r",
            copy_from,
            copy_to or "",
        )
        model = await ensure_async(self.contents_manager.copy(copy_from, copy_to))
        self.set_status(201)
        validate_model(model)
        self._finish_model(model)

    async def _upload(self, model, path):
        """Handle upload of a new file to path"""
        self.log.info("Uploading file to %s", path)
        model = await ensure_async(self.contents_manager.new(model, path))
        self.set_status(201)
        validate_model(model)
        self._finish_model(model)

    async def _new_untitled(self, path, type="", ext=""):
        """Create a new, empty untitled entity"""
        self.log.info("Creating new %s in %s", type or "file", path)
        model = await ensure_async(
            self.contents_manager.new_untitled(path=path, type=type, ext=ext)
        )
        self.set_status(201)
        validate_model(model)
        self._finish_model(model)

    async def _save(self, model, path):
        """Save an existing file."""
        chunk = model.get("chunk", None)
        if not chunk or chunk == -1:  # Avoid tedious log information
            self.log.info("Saving file at %s", path)
        model = await ensure_async(self.contents_manager.save(model, path))
        validate_model(model)
        self._finish_model(model)

    @web.authenticated
    @authorized
    async def post(self, path=""):
        """Create a new file in the specified path.

        POST creates new files. The server always decides on the name.

        POST /api/contents/path
          New untitled, empty file or directory.
        POST /api/contents/path
          with body {"copy_from" : "/path/to/OtherNotebook.ipynb"}
          New copy of OtherNotebook in path
        """

        cm = self.contents_manager

        file_exists = await ensure_async(cm.file_exists(path))
        if file_exists:
            raise web.HTTPError(400, "Cannot POST to files, use PUT instead.")

        model = self.get_json_body()
        if model:
            copy_from = model.get("copy_from")
            if copy_from:
                if not cm.allow_hidden and (
                    await ensure_async(cm.is_hidden(path))
                    or await ensure_async(cm.is_hidden(copy_from))
                ):
                    raise web.HTTPError(400, f"Cannot copy file or directory {path!r}")
                else:
                    await self._copy(copy_from, path)
            else:
                ext = model.get("ext", "")
                type = model.get("type", "")
                if type not in {None, "", "directory", "file", "notebook"}:
                    # fall back to file if unknown type
                    type = "file"
                await self._new_untitled(path, type=type, ext=ext)
        else:
            await self._new_untitled(path)

    @web.authenticated
    @authorized
    async def put(self, path=""):
        """Saves the file in the location specified by name and path.

        PUT is very similar to POST, but the requester specifies the name,
        whereas with POST, the server picks the name.

        PUT /api/contents/path/Name.ipynb
          Save notebook at ``path/Name.ipynb``. Notebook structure is specified
          in `content` key of JSON request body. If content is not specified,
          create a new empty notebook.
        """
        model = self.get_json_body()
        cm = self.contents_manager

        if model:
            if model.get("copy_from"):
                raise web.HTTPError(400, "Cannot copy with PUT, only POST")
            if not cm.allow_hidden and (
                (model.get("path") and await ensure_async(cm.is_hidden(model.get("path"))))
                or await ensure_async(cm.is_hidden(path))
            ):
                raise web.HTTPError(400, f"Cannot create file or directory {path!r}")

            exists = await ensure_async(self.contents_manager.file_exists(path))
            if model.get("type", "") not in {None, "", "directory", "file", "notebook"}:
                # fall back to file if unknown type
                model["type"] = "file"
            if exists:
                await self._save(model, path)
            else:
                await self._upload(model, path)
        else:
            await self._new_untitled(path)

    @web.authenticated
    @authorized
    async def delete(self, path=""):
        """delete a file in the given path"""
        cm = self.contents_manager

        if not cm.allow_hidden and await ensure_async(cm.is_hidden(path)):
            raise web.HTTPError(400, f"Cannot delete file or directory {path!r}")

        self.log.warning("delete %s", path)
        await ensure_async(cm.delete(path))
        self.set_status(204)
        self.finish()


class CheckpointsHandler(ContentsAPIHandler):
    """A checkpoints API handler."""

    @web.authenticated
    @authorized
    async def get(self, path=""):
        """get lists checkpoints for a file"""
        cm = self.contents_manager
        checkpoints = await ensure_async(cm.list_checkpoints(path))
        data = json.dumps(checkpoints, default=json_default)
        self.finish(data)

    @web.authenticated
    @authorized
    async def post(self, path=""):
        """post creates a new checkpoint"""
        cm = self.contents_manager
        checkpoint = await ensure_async(cm.create_checkpoint(path))
        data = json.dumps(checkpoint, default=json_default)
        location = url_path_join(
            self.base_url,
            "api/contents",
            url_escape(path),
            "checkpoints",
            url_escape(checkpoint["id"]),
        )
        self.set_header("Location", location)
        self.set_status(201)
        self.finish(data)


class ModifyCheckpointsHandler(ContentsAPIHandler):
    """A checkpoints modification handler."""

    @web.authenticated
    @authorized
    async def post(self, path, checkpoint_id):
        """post restores a file from a checkpoint"""
        cm = self.contents_manager
        await ensure_async(cm.restore_checkpoint(checkpoint_id, path))
        self.set_status(204)
        self.finish()

    @web.authenticated
    @authorized
    async def delete(self, path, checkpoint_id):
        """delete clears a checkpoint for a given file"""
        cm = self.contents_manager
        await ensure_async(cm.delete_checkpoint(checkpoint_id, path))
        self.set_status(204)
        self.finish()


class NotebooksRedirectHandler(JupyterHandler):
    """Redirect /api/notebooks to /api/contents"""

    SUPPORTED_METHODS = (
        "GET",
        "PUT",
        "PATCH",
        "POST",
        "DELETE",
    )

    @allow_unauthenticated
    def get(self, path):
        """Handle a notebooks redirect."""
        self.log.warning("/api/notebooks is deprecated, use /api/contents")
        self.redirect(url_path_join(self.base_url, "api/contents", url_escape(path)))

    put = patch = post = delete = get


class TrustNotebooksHandler(JupyterHandler):
    """Handles trust/signing of notebooks"""

    @web.authenticated  # type:ignore[misc]
    @authorized(resource=AUTH_RESOURCE)
    async def post(self, path=""):
        """Trust a notebook by path."""
        cm = self.contents_manager
        await ensure_async(cm.trust_notebook(path))
        self.set_status(201)
        self.finish()


# -----------------------------------------------------------------------------
# URL to handler mappings
# -----------------------------------------------------------------------------


_checkpoint_id_regex = r"(?P<checkpoint_id>[\w-]+)"


default_handlers = [
    (r"/api/contents%s/checkpoints" % path_regex, CheckpointsHandler),
    (
        rf"/api/contents{path_regex}/checkpoints/{_checkpoint_id_regex}",
        ModifyCheckpointsHandler,
    ),
    (r"/api/contents%s/trust" % path_regex, TrustNotebooksHandler),
    (r"/api/contents%s" % path_regex, ContentsHandler),
    (r"/api/notebooks/?(.*)", NotebooksRedirectHandler),
]
