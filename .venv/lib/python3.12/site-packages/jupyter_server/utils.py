"""Notebook related utilities"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import errno
import importlib.util
import os
import socket
import sys
import warnings
from _frozen_importlib_external import _NamespacePath
from contextlib import contextmanager
from pathlib import Path
from typing import TYPE_CHECKING, Any, NewType
from urllib.parse import (
    SplitResult,
    quote,
    unquote,
    urlparse,
    urlsplit,
    urlunsplit,
)
from urllib.parse import (
    urljoin as _urljoin,
)
from urllib.request import pathname2url as _pathname2url

from jupyter_core.utils import ensure_async as _ensure_async
from packaging.version import Version
from tornado.httpclient import AsyncHTTPClient, HTTPClient, HTTPRequest, HTTPResponse
from tornado.netutil import Resolver

if TYPE_CHECKING:
    from collections.abc import Generator, Sequence

ApiPath = NewType("ApiPath", str)

# Re-export
urljoin = _urljoin
pathname2url = _pathname2url
ensure_async = _ensure_async


def url_path_join(*pieces: str) -> str:
    """Join components of url into a relative url

    Use to prevent double slash when joining subpath. This will leave the
    initial and final / in place
    """
    initial = pieces[0].startswith("/")
    final = pieces[-1].endswith("/")
    stripped = [s.strip("/") for s in pieces]
    result = "/".join(s for s in stripped if s)
    if initial:
        result = "/" + result
    if final:
        result = result + "/"
    if result == "//":
        result = "/"
    return result


def url_is_absolute(url: str) -> bool:
    """Determine whether a given URL is absolute"""
    return urlparse(url).path.startswith("/")


def path2url(path: str) -> str:
    """Convert a local file path to a URL"""
    pieces = [quote(p) for p in path.split(os.sep)]
    # preserve trailing /
    if pieces[-1] == "":
        pieces[-1] = "/"
    url = url_path_join(*pieces)
    return url


def url2path(url: str) -> str:
    """Convert a URL to a local file path"""
    pieces = [unquote(p) for p in url.split("/")]
    path = os.path.join(*pieces)
    return path


def url_escape(path: str) -> str:
    """Escape special characters in a URL path

    Turns '/foo bar/' into '/foo%20bar/'
    """
    parts = path.split("/")
    return "/".join([quote(p) for p in parts])


def url_unescape(path: str) -> str:
    """Unescape special characters in a URL path

    Turns '/foo%20bar/' into '/foo bar/'
    """
    return "/".join([unquote(p) for p in path.split("/")])


def samefile_simple(path: str, other_path: str) -> bool:
    """
    Fill in for os.path.samefile when it is unavailable (Windows+py2).

    Do a case-insensitive string comparison in this case
    plus comparing the full stat result (including times)
    because Windows + py2 doesn't support the stat fields
    needed for identifying if it's the same file (st_ino, st_dev).

    Only to be used if os.path.samefile is not available.

    Parameters
    ----------
    path : str
        representing a path to a file
    other_path : str
        representing a path to another file

    Returns
    -------
    same:   Boolean that is True if both path and other path are the same
    """
    path_stat = os.stat(path)
    other_path_stat = os.stat(other_path)
    return path.lower() == other_path.lower() and path_stat == other_path_stat


def to_os_path(path: ApiPath, root: str = "") -> str:
    """Convert an API path to a filesystem path

    If given, root will be prepended to the path.
    root must be a filesystem path already.
    """
    parts = str(path).strip("/").split("/")
    parts = [p for p in parts if p != ""]  #  remove duplicate splits
    path_ = os.path.join(root, *parts)
    return os.path.normpath(path_)


def to_api_path(os_path: str, root: str = "") -> ApiPath:
    """Convert a filesystem path to an API path

    If given, root will be removed from the path.
    root must be a filesystem path already.
    """
    if os_path.startswith(root):
        os_path = os_path[len(root) :]
    parts = os_path.strip(os.path.sep).split(os.path.sep)
    parts = [p for p in parts if p != ""]  # remove duplicate splits
    path = "/".join(parts)
    return ApiPath(path)


def check_version(v: str, check: str) -> bool:
    """check version string v >= check

    If dev/prerelease tags result in TypeError for string-number comparison,
    it is assumed that the dependency is satisfied.
    Users on dev branches are responsible for keeping their own packages up to date.
    """
    try:
        return bool(Version(v) >= Version(check))
    except TypeError:
        return True


# Copy of IPython.utils.process.check_pid:


def _check_pid_win32(pid: int) -> bool:
    import ctypes

    # OpenProcess returns 0 if no such process (of ours) exists
    # positive int otherwise
    return bool(ctypes.windll.kernel32.OpenProcess(1, 0, pid))  # type:ignore[attr-defined]


def _check_pid_posix(pid: int) -> bool:
    """Copy of IPython.utils.process.check_pid"""
    try:
        os.kill(pid, 0)
    except OSError as err:
        if err.errno == errno.ESRCH:
            return False
        elif err.errno == errno.EPERM:
            # Don't have permission to signal the process - probably means it exists
            return True
        raise
    else:
        return True


if sys.platform == "win32":
    check_pid = _check_pid_win32
else:
    check_pid = _check_pid_posix


async def run_sync_in_loop(maybe_async):
    """**DEPRECATED**: Use ``ensure_async`` from jupyter_core instead."""
    warnings.warn(
        "run_sync_in_loop is deprecated since Jupyter Server 2.0, use 'ensure_async' from jupyter_core instead",
        DeprecationWarning,
        stacklevel=2,
    )
    return ensure_async(maybe_async)


def urlencode_unix_socket_path(socket_path: str) -> str:
    """Encodes a UNIX socket path string from a socket path for the `http+unix` URI form."""
    return socket_path.replace("/", "%2F")


def urldecode_unix_socket_path(socket_path: str) -> str:
    """Decodes a UNIX sock path string from an encoded sock path for the `http+unix` URI form."""
    return socket_path.replace("%2F", "/")


def urlencode_unix_socket(socket_path: str) -> str:
    """Encodes a UNIX socket URL from a socket path for the `http+unix` URI form."""
    return "http+unix://%s" % urlencode_unix_socket_path(socket_path)


def unix_socket_in_use(socket_path: str) -> bool:
    """Checks whether a UNIX socket path on disk is in use by attempting to connect to it."""
    if not os.path.exists(socket_path):
        return False

    try:
        sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        sock.connect(socket_path)
    except OSError:
        return False
    else:
        return True
    finally:
        sock.close()


@contextmanager
def _request_for_tornado_client(
    urlstring: str, method: str = "GET", body: Any = None, headers: Any = None
) -> Generator[HTTPRequest, None, None]:
    """A utility that provides a context that handles
    HTTP, HTTPS, and HTTP+UNIX request.
    Creates a tornado HTTPRequest object with a URL
    that tornado's HTTPClients can accept.
    If the request is made to a unix socket, temporarily
    configure the AsyncHTTPClient to resolve the URL
    and connect to the proper socket.
    """
    parts = urlsplit(urlstring)
    if parts.scheme in ["http", "https"]:
        pass
    elif parts.scheme == "http+unix":
        # If unix socket, mimic HTTP.
        parts = SplitResult(
            scheme="http",
            netloc=parts.netloc,
            path=parts.path,
            query=parts.query,
            fragment=parts.fragment,
        )

        class UnixSocketResolver(Resolver):
            """A resolver that routes HTTP requests to unix sockets
            in tornado HTTP clients.
            Due to constraints in Tornados' API, the scheme of the
            must be `http` (not `http+unix`). Applications should replace
            the scheme in URLS before making a request to the HTTP client.
            """

            def initialize(self, resolver):
                self.resolver = resolver

            def close(self):
                self.resolver.close()

            async def resolve(self, host, port, *args, **kwargs):
                return [(socket.AF_UNIX, urldecode_unix_socket_path(host))]

        resolver = UnixSocketResolver(resolver=Resolver())
        AsyncHTTPClient.configure(None, resolver=resolver)
    else:
        msg = "Unknown URL scheme."
        raise Exception(msg)

    # Yield the request for the given client.
    url = urlunsplit(parts)
    request = HTTPRequest(url, method=method, body=body, headers=headers, validate_cert=False)
    yield request


def fetch(
    urlstring: str, method: str = "GET", body: Any = None, headers: Any = None
) -> HTTPResponse:
    """
    Send a HTTP, HTTPS, or HTTP+UNIX request
    to a Tornado Web Server. Returns a tornado HTTPResponse.
    """
    with _request_for_tornado_client(
        urlstring, method=method, body=body, headers=headers
    ) as request:
        response = HTTPClient(AsyncHTTPClient).fetch(request)
    return response


async def async_fetch(
    urlstring: str, method: str = "GET", body: Any = None, headers: Any = None, io_loop: Any = None
) -> HTTPResponse:
    """
    Send an asynchronous HTTP, HTTPS, or HTTP+UNIX request
    to a Tornado Web Server. Returns a tornado HTTPResponse.
    """
    with _request_for_tornado_client(
        urlstring, method=method, body=body, headers=headers
    ) as request:
        response = await AsyncHTTPClient(io_loop).fetch(request)
    return response


def is_namespace_package(namespace: str) -> bool | None:
    """Is the provided namespace a Python Namespace Package (PEP420).

    https://www.python.org/dev/peps/pep-0420/#specification

    Returns `None` if module is not importable.

    """
    # NOTE: using submodule_search_locations because the loader can be None
    try:
        spec = importlib.util.find_spec(namespace)
    except ValueError:  # spec is not set - see https://docs.python.org/3/library/importlib.html#importlib.util.find_spec
        return None

    if not spec:
        # e.g. module not installed
        return None
    return isinstance(spec.submodule_search_locations, _NamespacePath)


def filefind(filename: str, path_dirs: Sequence[str]) -> str:
    """Find a file by looking through a sequence of paths.

    For use in FileFindHandler.

    Iterates through a sequence of paths looking for a file and returns
    the full, absolute path of the first occurrence of the file.

    Absolute paths are not accepted for inputs.

    This function does not automatically try any paths,
    such as the cwd or the user's home directory.

    Parameters
    ----------
    filename : str
        The filename to look for. Must be a relative path.
    path_dirs : sequence of str
        The sequence of paths to look in for the file.
        Walk through each element and join with ``filename``.
        Only after ensuring the path resolves within the directory is it checked for existence.

    Returns
    -------
    Raises :exc:`OSError` or returns absolute path to file.
    """
    file_path = Path(filename)

    # If the input is an absolute path, reject it
    if file_path.is_absolute():
        msg = f"{filename} is absolute, filefind only accepts relative paths."
        raise OSError(msg)

    for path_str in path_dirs:
        path = Path(path_str).absolute()
        test_path = path / file_path
        # os.path.abspath resolves '..', but Path.absolute() doesn't
        # Path.resolve() does, but traverses symlinks, which we don't want
        test_path = Path(os.path.abspath(test_path))
        if not test_path.is_relative_to(path):
            # points outside root, e.g. via `filename='../foo'`
            continue
        # make sure we don't call is_file before we know it's a file within a prefix
        # GHSA-hrw6-wg82-cm62 - can leak password hash on windows.
        if test_path.is_file():
            return os.path.abspath(test_path)

    msg = f"File {filename!r} does not exist in any of the search paths: {path_dirs!r}"
    raise OSError(msg)


def import_item(name: str) -> Any:
    """Import and return ``bar`` given the string ``foo.bar``.
    Calling ``bar = import_item("foo.bar")`` is the functional equivalent of
    executing the code ``from foo import bar``.
    Parameters
    ----------
    name : str
      The fully qualified name of the module/package being imported.
    Returns
    -------
    mod : module object
       The module that was imported.
    """

    parts = name.rsplit(".", 1)
    if len(parts) == 2:
        # called with 'foo.bar....'
        package, obj = parts
        module = __import__(package, fromlist=[obj])
        try:
            pak = getattr(module, obj)
        except AttributeError as e:
            raise ImportError("No module named %s" % obj) from e
        return pak
    else:
        # called with un-dotted string
        return __import__(parts[0])


class JupyterServerAuthWarning(RuntimeWarning):
    """Emitted when authentication configuration issue is detected.

    Intended for filtering out expected warnings in tests, including
    downstream tests, rather than for users to silence this warning.
    """
