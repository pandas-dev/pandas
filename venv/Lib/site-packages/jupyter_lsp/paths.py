import os
import re
from pathlib import Path
from typing import Union
from urllib.parse import unquote, urlparse

RE_PATH_ANCHOR = r"^file://([^/]+|/[A-Z]:)"


def normalized_uri(root_dir):
    """Attempt to make an LSP rootUri from a ContentsManager root_dir

    Special care must be taken around windows paths: the canonical form of
    windows drives and UNC paths is lower case
    """
    root_uri = Path(root_dir).expanduser().resolve().as_uri()
    root_uri = re.sub(
        RE_PATH_ANCHOR, lambda m: "file://{}".format(m.group(1).lower()), root_uri
    )
    return root_uri


def file_uri_to_path(file_uri):
    """Return a path string for give file:/// URI.

    Respect the different path convention on Windows.
    Based on https://stackoverflow.com/a/57463161/6646912, BSD 0
    """
    windows_path = os.name == "nt"
    file_uri_parsed = urlparse(file_uri)
    file_uri_path_unquoted = unquote(file_uri_parsed.path)
    if windows_path and file_uri_path_unquoted.startswith("/"):
        result = file_uri_path_unquoted[1:]  # pragma: no cover
    else:
        result = file_uri_path_unquoted  # pragma: no cover
    return result


def is_relative(root: Union[str, Path], path: Union[str, Path]) -> bool:
    """Return if path is relative to root"""
    try:
        Path(path).resolve().relative_to(Path(root).resolve())
        return True
    except ValueError:
        return False
