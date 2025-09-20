import base64
import os

from anyio.to_thread import run_sync
from tornado import web

from jupyter_server.services.contents.filemanager import (
    AsyncFileContentsManager,
    FileContentsManager,
)


class LargeFileManager(FileContentsManager):
    """Handle large file upload."""

    def save(self, model, path=""):
        """Save the file model and return the model with no content."""
        chunk = model.get("chunk", None)
        if chunk is not None:
            path = path.strip("/")

            if chunk == 1:
                self.run_pre_save_hooks(model=model, path=path)

            if "type" not in model:
                raise web.HTTPError(400, "No file type provided")
            if model["type"] != "file":
                raise web.HTTPError(
                    400,
                    'File type "{}" is not supported for large file transfer'.format(model["type"]),
                )
            if "content" not in model and model["type"] != "directory":
                raise web.HTTPError(400, "No file content provided")

            os_path = self._get_os_path(path)
            if chunk == -1:
                self.log.debug(f"Saving last chunk of file {os_path}")
            else:
                self.log.debug(f"Saving chunk {chunk} of file {os_path}")

            try:
                if chunk == 1:
                    super()._save_file(os_path, model["content"], model.get("format"))
                else:
                    self._save_large_file(os_path, model["content"], model.get("format"))
            except web.HTTPError:
                raise
            except Exception as e:
                self.log.error("Error while saving file: %s %s", path, e, exc_info=True)
                raise web.HTTPError(500, f"Unexpected error while saving file: {path} {e}") from e

            model = self.get(path, content=False)

            # Last chunk
            if chunk == -1:
                self.run_post_save_hooks(model=model, os_path=os_path)
            self.emit(data={"action": "save", "path": path})
            return model
        else:
            return super().save(model, path)

    def _save_large_file(self, os_path, content, format):
        """Save content of a generic file."""
        if format not in {"text", "base64"}:
            raise web.HTTPError(
                400,
                "Must specify format of file contents as 'text' or 'base64'",
            )
        try:
            if format == "text":
                bcontent = content.encode("utf8")
            else:
                b64_bytes = content.encode("ascii")
                bcontent = base64.b64decode(b64_bytes)
        except Exception as e:
            raise web.HTTPError(400, f"Encoding error saving {os_path}: {e}") from e

        with self.perm_to_403(os_path):
            if os.path.islink(os_path):
                os_path = os.path.join(os.path.dirname(os_path), os.readlink(os_path))
            with open(os_path, "ab") as f:
                f.write(bcontent)


class AsyncLargeFileManager(AsyncFileContentsManager):
    """Handle large file upload asynchronously"""

    async def save(self, model, path=""):
        """Save the file model and return the model with no content."""
        chunk = model.get("chunk", None)
        if chunk is not None:
            path = path.strip("/")

            if chunk == 1:
                self.run_pre_save_hooks(model=model, path=path)

            if "type" not in model:
                raise web.HTTPError(400, "No file type provided")
            if model["type"] != "file":
                raise web.HTTPError(
                    400,
                    'File type "{}" is not supported for large file transfer'.format(model["type"]),
                )
            if "content" not in model and model["type"] != "directory":
                raise web.HTTPError(400, "No file content provided")

            os_path = self._get_os_path(path)
            if chunk == -1:
                self.log.debug(f"Saving last chunk of file {os_path}")
            else:
                self.log.debug(f"Saving chunk {chunk} of file {os_path}")

            try:
                if chunk == 1:
                    await super()._save_file(os_path, model["content"], model.get("format"))
                else:
                    await self._save_large_file(os_path, model["content"], model.get("format"))
            except web.HTTPError:
                raise
            except Exception as e:
                self.log.error("Error while saving file: %s %s", path, e, exc_info=True)
                raise web.HTTPError(500, f"Unexpected error while saving file: {path} {e}") from e

            model = await self.get(path, content=False)

            # Last chunk
            if chunk == -1:
                self.run_post_save_hooks(model=model, os_path=os_path)

            self.emit(data={"action": "save", "path": path})
            return model
        else:
            return await super().save(model, path)

    async def _save_large_file(self, os_path, content, format):
        """Save content of a generic file."""
        if format not in {"text", "base64"}:
            raise web.HTTPError(
                400,
                "Must specify format of file contents as 'text' or 'base64'",
            )
        try:
            if format == "text":
                bcontent = content.encode("utf8")
            else:
                b64_bytes = content.encode("ascii")
                bcontent = base64.b64decode(b64_bytes)
        except Exception as e:
            raise web.HTTPError(400, f"Encoding error saving {os_path}: {e}") from e

        with self.perm_to_403(os_path):
            if os.path.islink(os_path):
                os_path = os.path.join(os.path.dirname(os_path), os.readlink(os_path))
            with open(os_path, "ab") as f:  # noqa: ASYNC101
                await run_sync(f.write, bcontent)
