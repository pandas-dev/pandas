"""
Utilities for file-based Contents/Checkpoints managers.
"""
# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from __future__ import annotations

import errno
import hashlib
import os
import shutil
from base64 import decodebytes, encodebytes
from contextlib import contextmanager
from functools import partial

import nbformat
from anyio.to_thread import run_sync
from tornado.web import HTTPError
from traitlets import Bool, Enum
from traitlets.config import Configurable
from traitlets.config.configurable import LoggingConfigurable

from jupyter_server.utils import ApiPath, to_api_path, to_os_path


def replace_file(src, dst):
    """replace dst with src"""
    os.replace(src, dst)


async def async_replace_file(src, dst):
    """replace dst with src asynchronously"""
    await run_sync(os.replace, src, dst)


def copy2_safe(src, dst, log=None):
    """copy src to dst

    like shutil.copy2, but log errors in copystat instead of raising
    """
    shutil.copyfile(src, dst)
    try:
        shutil.copystat(src, dst)
    except OSError:
        if log:
            log.debug("copystat on %s failed", dst, exc_info=True)


async def async_copy2_safe(src, dst, log=None):
    """copy src to dst asynchronously

    like shutil.copy2, but log errors in copystat instead of raising
    """
    await run_sync(shutil.copyfile, src, dst)
    try:
        await run_sync(shutil.copystat, src, dst)
    except OSError:
        if log:
            log.debug("copystat on %s failed", dst, exc_info=True)


def path_to_intermediate(path):
    """Name of the intermediate file used in atomic writes.

    The .~ prefix will make Dropbox ignore the temporary file."""
    dirname, basename = os.path.split(path)
    return os.path.join(dirname, ".~" + basename)


def path_to_invalid(path):
    """Name of invalid file after a failed atomic write and subsequent read."""
    dirname, basename = os.path.split(path)
    return os.path.join(dirname, basename + ".invalid")


@contextmanager
def atomic_writing(path, text=True, encoding="utf-8", log=None, **kwargs):
    """Context manager to write to a file only if the entire write is successful.

    This works by copying the previous file contents to a temporary file in the
    same directory, and renaming that file back to the target if the context
    exits with an error. If the context is successful, the new data is synced to
    disk and the temporary file is removed.

    Parameters
    ----------
    path : str
        The target file to write to.
    text : bool, optional
        Whether to open the file in text mode (i.e. to write unicode). Default is
        True.
    encoding : str, optional
        The encoding to use for files opened in text mode. Default is UTF-8.
    **kwargs
        Passed to :func:`io.open`.
    """
    # realpath doesn't work on Windows: https://bugs.python.org/issue9949
    # Luckily, we only need to resolve the file itself being a symlink, not
    # any of its directories, so this will suffice:
    if os.path.islink(path):
        path = os.path.join(os.path.dirname(path), os.readlink(path))

    tmp_path = path_to_intermediate(path)

    if os.path.isfile(path):
        copy2_safe(path, tmp_path, log=log)

    if text:
        # Make sure that text files have Unix linefeeds by default
        kwargs.setdefault("newline", "\n")
        fileobj = open(path, "w", encoding=encoding, **kwargs)  # noqa: SIM115
    else:
        fileobj = open(path, "wb", **kwargs)  # noqa: SIM115

    try:
        yield fileobj
    except BaseException:
        # Failed! Move the backup file back to the real path to avoid corruption
        fileobj.close()
        replace_file(tmp_path, path)
        raise

    # Flush to disk
    fileobj.flush()
    os.fsync(fileobj.fileno())
    fileobj.close()

    # Written successfully, now remove the backup copy
    if os.path.isfile(tmp_path):
        os.remove(tmp_path)


@contextmanager
def _simple_writing(path, text=True, encoding="utf-8", log=None, **kwargs):
    """Context manager to write file without doing atomic writing
    (for weird filesystem eg: nfs).

    Parameters
    ----------
    path : str
        The target file to write to.
    text : bool, optional
        Whether to open the file in text mode (i.e. to write unicode). Default is
        True.
    encoding : str, optional
        The encoding to use for files opened in text mode. Default is UTF-8.
    **kwargs
        Passed to :func:`io.open`.
    """
    # realpath doesn't work on Windows: https://bugs.python.org/issue9949
    # Luckily, we only need to resolve the file itself being a symlink, not
    # any of its directories, so this will suffice:
    if os.path.islink(path):
        path = os.path.join(os.path.dirname(path), os.readlink(path))

    if text:
        # Make sure that text files have Unix linefeeds by default
        kwargs.setdefault("newline", "\n")
        fileobj = open(path, "w", encoding=encoding, **kwargs)  # noqa: SIM115
    else:
        fileobj = open(path, "wb", **kwargs)  # noqa: SIM115

    try:
        yield fileobj
    except BaseException:
        fileobj.close()
        raise

    fileobj.close()


class FileManagerMixin(LoggingConfigurable, Configurable):
    """
    Mixin for ContentsAPI classes that interact with the filesystem.

    Provides facilities for reading, writing, and copying files.

    Shared by FileContentsManager and FileCheckpoints.

    Note
    ----
    Classes using this mixin must provide the following attributes:

    root_dir : unicode
        A directory against against which API-style paths are to be resolved.

    log : logging.Logger
    """

    use_atomic_writing = Bool(
        True,
        config=True,
        help="""By default notebooks are saved on disk on a temporary file and then if succefully written, it replaces the old ones.
      This procedure, namely 'atomic_writing', causes some bugs on file system without operation order enforcement (like some networked fs).
      If set to False, the new notebook is written directly on the old one which could fail (eg: full filesystem or quota )""",
    )

    hash_algorithm = Enum(  # type: ignore[call-overload]
        hashlib.algorithms_available,
        default_value="sha256",
        config=True,
        help="Hash algorithm to use for file content, support by hashlib",
    )

    @contextmanager
    def open(self, os_path, *args, **kwargs):
        """wrapper around io.open that turns permission errors into 403"""
        with self.perm_to_403(os_path), open(os_path, *args, **kwargs) as f:
            yield f

    @contextmanager
    def atomic_writing(self, os_path, *args, **kwargs):
        """wrapper around atomic_writing that turns permission errors to 403.
        Depending on flag 'use_atomic_writing', the wrapper perform an actual atomic writing or
        simply writes the file (whatever an old exists or not)"""
        with self.perm_to_403(os_path):
            kwargs["log"] = self.log
            if self.use_atomic_writing:
                with atomic_writing(os_path, *args, **kwargs) as f:
                    yield f
            else:
                with _simple_writing(os_path, *args, **kwargs) as f:
                    yield f

    @contextmanager
    def perm_to_403(self, os_path=""):
        """context manager for turning permission errors into 403."""
        try:
            yield
        except OSError as e:
            if e.errno in {errno.EPERM, errno.EACCES}:
                # make 403 error message without root prefix
                # this may not work perfectly on unicode paths on Python 2,
                # but nobody should be doing that anyway.
                if not os_path:
                    os_path = e.filename or "unknown file"
                path = to_api_path(os_path, root=self.root_dir)  # type:ignore[attr-defined]
                raise HTTPError(403, "Permission denied: %s" % path) from e
            else:
                raise

    def _copy(self, src, dest):
        """copy src to dest

        like shutil.copy2, but log errors in copystat
        """
        copy2_safe(src, dest, log=self.log)

    def _get_os_path(self, path):
        """Given an API path, return its file system path.

        Parameters
        ----------
        path : str
            The relative API path to the named file.

        Returns
        -------
        path : str
            Native, absolute OS path to for a file.

        Raises
        ------
        404: if path is outside root
        """
        # This statement can cause excessive logging, uncomment if necessary when troubleshooting.
        # self.log.debug("Reading path from disk: %s", path)
        root = os.path.abspath(self.root_dir)  # type:ignore[attr-defined]
        # to_os_path is not safe if path starts with a drive, since os.path.join discards first part
        if os.path.splitdrive(path)[0]:
            raise HTTPError(404, "%s is not a relative API path" % path)
        os_path = to_os_path(ApiPath(path), root)
        # validate os path
        # e.g. "foo\0" raises ValueError: embedded null byte
        try:
            os.lstat(os_path)
        except OSError:
            # OSError could be FileNotFound, PermissionError, etc.
            # those should raise (or not) elsewhere
            pass
        except ValueError:
            raise HTTPError(404, f"{path} is not a valid path") from None

        if not (os.path.abspath(os_path) + os.path.sep).startswith(root):
            raise HTTPError(404, "%s is outside root contents directory" % path)
        return os_path

    def _read_notebook(
        self, os_path, as_version=4, capture_validation_error=None, raw: bool = False
    ):
        """Read a notebook from an os path."""
        answer = self._read_file(os_path, "text", raw=raw)

        try:
            nb = nbformat.reads(
                answer[0],
                as_version=as_version,
                capture_validation_error=capture_validation_error,
            )

            return (nb, answer[2]) if raw else nb  # type:ignore[misc]
        except Exception as e:
            e_orig = e

        # If use_atomic_writing is enabled, we'll guess that it was also
        # enabled when this notebook was written and look for a valid
        # atomic intermediate.
        tmp_path = path_to_intermediate(os_path)

        if not self.use_atomic_writing or not os.path.exists(tmp_path):
            raise HTTPError(
                400,
                f"Unreadable Notebook: {os_path} {e_orig!r}",
            )

        # Move the bad file aside, restore the intermediate, and try again.
        invalid_file = path_to_invalid(os_path)
        replace_file(os_path, invalid_file)
        replace_file(tmp_path, os_path)
        return self._read_notebook(
            os_path, as_version, capture_validation_error=capture_validation_error, raw=raw
        )

    def _save_notebook(self, os_path, nb, capture_validation_error=None):
        """Save a notebook to an os_path."""
        with self.atomic_writing(os_path, encoding="utf-8") as f:
            nbformat.write(
                nb,
                f,
                version=nbformat.NO_CONVERT,
                capture_validation_error=capture_validation_error,
            )

    def _get_hash(self, byte_content: bytes) -> dict[str, str]:
        """Compute the hash hexdigest for the provided bytes.

        The hash algorithm is provided by the `hash_algorithm` attribute.

        Parameters
        ----------
        byte_content : bytes
            The bytes to hash

        Returns
        -------
        A dictionary to be appended to a model {"hash": str, "hash_algorithm": str}.
        """
        algorithm = self.hash_algorithm
        h = hashlib.new(algorithm)
        h.update(byte_content)
        return {"hash": h.hexdigest(), "hash_algorithm": algorithm}

    def _read_file(
        self, os_path: str, format: str | None, raw: bool = False
    ) -> tuple[str | bytes, str] | tuple[str | bytes, str, bytes]:
        """Read a non-notebook file.

        Parameters
        ----------
        os_path: str
            The path to be read.
        format: str
            If 'text', the contents will be decoded as UTF-8.
            If 'base64', the raw bytes contents will be encoded as base64.
            If 'byte', the raw bytes contents will be returned.
            If not specified, try to decode as UTF-8, and fall back to base64
        raw: bool
            [Optional] If True, will return as third argument the raw bytes content

        Returns
        -------
        (content, format, byte_content) It returns the content in the given format
        as well as the raw byte content.
        """
        if not os.path.isfile(os_path):
            raise HTTPError(400, "Cannot read non-file %s" % os_path)

        with self.open(os_path, "rb") as f:
            bcontent = f.read()

        if format == "byte":
            # Not for http response but internal use
            return (bcontent, "byte", bcontent) if raw else (bcontent, "byte")

        if format is None or format == "text":
            # Try to interpret as unicode if format is unknown or if unicode
            # was explicitly requested.
            try:
                return (
                    (bcontent.decode("utf8"), "text", bcontent)
                    if raw
                    else (
                        bcontent.decode("utf8"),
                        "text",
                    )
                )
            except UnicodeError as e:
                if format == "text":
                    raise HTTPError(
                        400,
                        "%s is not UTF-8 encoded" % os_path,
                        reason="bad format",
                    ) from e
        return (
            (encodebytes(bcontent).decode("ascii"), "base64", bcontent)
            if raw
            else (
                encodebytes(bcontent).decode("ascii"),
                "base64",
            )
        )

    def _save_file(self, os_path, content, format):
        """Save content of a generic file."""
        if format not in {"text", "base64"}:
            raise HTTPError(
                400,
                "Must specify format of file contents as 'text' or 'base64'",
            )
        try:
            if format == "text":
                bcontent = content.encode("utf8")
            else:
                b64_bytes = content.encode("ascii")
                bcontent = decodebytes(b64_bytes)
        except Exception as e:
            raise HTTPError(400, f"Encoding error saving {os_path}: {e}") from e

        with self.atomic_writing(os_path, text=False) as f:
            f.write(bcontent)


class AsyncFileManagerMixin(FileManagerMixin):
    """
    Mixin for ContentsAPI classes that interact with the filesystem asynchronously.
    """

    async def _copy(self, src, dest):
        """copy src to dest

        like shutil.copy2, but log errors in copystat
        """
        await async_copy2_safe(src, dest, log=self.log)

    async def _read_notebook(
        self, os_path, as_version=4, capture_validation_error=None, raw: bool = False
    ):
        """Read a notebook from an os path."""
        answer = await self._read_file(os_path, "text", raw)

        try:
            nb = await run_sync(
                partial(
                    nbformat.reads,
                    as_version=as_version,
                    capture_validation_error=capture_validation_error,
                ),
                answer[0],
            )
            return (nb, answer[2]) if raw else nb  # type:ignore[misc]
        except Exception as e:
            e_orig = e

        # If use_atomic_writing is enabled, we'll guess that it was also
        # enabled when this notebook was written and look for a valid
        # atomic intermediate.
        tmp_path = path_to_intermediate(os_path)

        if not self.use_atomic_writing or not os.path.exists(tmp_path):
            raise HTTPError(
                400,
                f"Unreadable Notebook: {os_path} {e_orig!r}",
            )

        # Move the bad file aside, restore the intermediate, and try again.
        invalid_file = path_to_invalid(os_path)
        await async_replace_file(os_path, invalid_file)
        await async_replace_file(tmp_path, os_path)
        answer = await self._read_notebook(
            os_path, as_version, capture_validation_error=capture_validation_error, raw=raw
        )

        return answer

    async def _save_notebook(self, os_path, nb, capture_validation_error=None):
        """Save a notebook to an os_path."""
        with self.atomic_writing(os_path, encoding="utf-8") as f:
            await run_sync(
                partial(
                    nbformat.write,
                    version=nbformat.NO_CONVERT,
                    capture_validation_error=capture_validation_error,
                ),
                nb,
                f,
            )

    async def _read_file(  # type: ignore[override]
        self, os_path: str, format: str | None, raw: bool = False
    ) -> tuple[str | bytes, str] | tuple[str | bytes, str, bytes]:
        """Read a non-notebook file.

        Parameters
        ----------
        os_path: str
            The path to be read.
        format: str
            If 'text', the contents will be decoded as UTF-8.
            If 'base64', the raw bytes contents will be encoded as base64.
            If 'byte', the raw bytes contents will be returned.
            If not specified, try to decode as UTF-8, and fall back to base64
        raw: bool
            [Optional] If True, will return as third argument the raw bytes content

        Returns
        -------
        (content, format, byte_content) It returns the content in the given format
        as well as the raw byte content.
        """
        if not os.path.isfile(os_path):
            raise HTTPError(400, "Cannot read non-file %s" % os_path)

        with self.open(os_path, "rb") as f:
            bcontent = await run_sync(f.read)

        if format == "byte":
            # Not for http response but internal use
            return (bcontent, "byte", bcontent) if raw else (bcontent, "byte")

        if format is None or format == "text":
            # Try to interpret as unicode if format is unknown or if unicode
            # was explicitly requested.
            try:
                return (
                    (bcontent.decode("utf8"), "text", bcontent)
                    if raw
                    else (
                        bcontent.decode("utf8"),
                        "text",
                    )
                )
            except UnicodeError as e:
                if format == "text":
                    raise HTTPError(
                        400,
                        "%s is not UTF-8 encoded" % os_path,
                        reason="bad format",
                    ) from e
        return (
            (encodebytes(bcontent).decode("ascii"), "base64", bcontent)
            if raw
            else (encodebytes(bcontent).decode("ascii"), "base64")
        )

    async def _save_file(self, os_path, content, format):
        """Save content of a generic file."""
        if format not in {"text", "base64"}:
            raise HTTPError(
                400,
                "Must specify format of file contents as 'text' or 'base64'",
            )
        try:
            if format == "text":
                bcontent = content.encode("utf8")
            else:
                b64_bytes = content.encode("ascii")
                bcontent = decodebytes(b64_bytes)
        except Exception as e:
            raise HTTPError(400, f"Encoding error saving {os_path}: {e}") from e

        with self.atomic_writing(os_path, text=False) as f:
            await run_sync(f.write, bcontent)
