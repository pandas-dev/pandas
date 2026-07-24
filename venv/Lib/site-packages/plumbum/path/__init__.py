from __future__ import annotations

__lazy_modules__ = {
    "plumbum.path.base",
    "plumbum.path.local",
    "plumbum.path.remote",
    "plumbum.path.utils",
}

from plumbum.path.base import FSUser, Path, RelativePath
from plumbum.path.local import LocalPath, LocalWorkdir
from plumbum.path.remote import RemotePath, RemoteWorkdir
from plumbum.path.utils import copy, delete, move

__all__ = (
    "FSUser",
    "LocalPath",
    "LocalWorkdir",
    "Path",
    "RelativePath",
    "RemotePath",
    "RemoteWorkdir",
    "copy",
    "delete",
    "move",
)
