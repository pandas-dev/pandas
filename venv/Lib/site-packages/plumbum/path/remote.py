from __future__ import annotations

__lazy_modules__ = {
    "contextlib",
    "copy",
    "plumbum.commands",
    "urllib",
    "urllib.request",
}

import copy
import errno
import os
import typing
import urllib.request as urllib
from contextlib import contextmanager

from plumbum.commands import ProcessExecutionError, shquote
from plumbum.path.base import FSUser, Path

if typing.TYPE_CHECKING:
    import builtins
    from collections.abc import Generator, Iterable

    from plumbum._compat.typing import Self
    from plumbum.machines.remote import BaseRemoteMachine
else:
    BaseRemoteMachine = typing.Any
    Self = typing.Any


class StatRes:
    """POSIX-like stat result"""

    __slots__ = ("_tup",)
    _tup: tuple[int, int, int, int, int, int, int, float, float, float]

    def __init__(
        self, tup: tuple[int, int, int, int, int, int, int, float, float, float]
    ):
        self._tup = copy.copy(tup)

    def __getitem__(self, index: int) -> int | float:
        return self._tup[index]

    @property
    def st_mode(self) -> int:
        return self._tup[0]

    @property
    def st_ino(self) -> int:
        return self._tup[1]

    @property
    def st_dev(self) -> int:
        return self._tup[2]

    @property
    def st_nlink(self) -> int:
        return self._tup[3]

    @property
    def st_uid(self) -> int:
        return self._tup[4]

    @property
    def st_gid(self) -> int:
        return self._tup[5]

    @property
    def st_size(self) -> int:
        return self._tup[6]

    @property
    def st_atime(self) -> float:
        return self._tup[7]

    @property
    def st_mtime(self) -> float:
        return self._tup[8]

    @property
    def st_ctime(self) -> float:
        return self._tup[9]

    # Aliases for backward compatibility / convenience
    mode = st_mode
    ino = st_ino
    dev = st_dev
    nlink = st_nlink
    uid = st_uid
    gid = st_gid
    size = st_size
    atime = st_atime
    mtime = st_mtime
    ctime = st_ctime


class RemoteStatRes(StatRes):
    """Remote POSIX-like stat result"""

    __slots__ = ("text_mode",)

    text_mode: str  # e.g., "directory", "regular file", etc.


class RemotePath(Path):
    """The class implementing remote-machine paths"""

    __slots__ = ("CASE_SENSITIVE", "remote")
    remote: BaseRemoteMachine

    def __new__(cls, remote: BaseRemoteMachine, *parts_str: str | Path) -> Self:
        if not parts_str:
            raise TypeError("At least one path part is required (none given)")
        windows = remote.uname.lower() == "windows"
        normed: list[str] = []

        parts: tuple[str, ...] = tuple(
            map(str, parts_str)
        )  # force the paths into string, so subscription works properly
        # Simple skip if path is absolute
        if parts[0] and parts[0][0] not in ("/", "\\"):
            cwd = (
                remote._cwd
                if hasattr(remote, "_cwd")
                else remote._session.run("pwd")[1].strip()
            )
            parts = (cwd, *parts)

        for p in parts:
            if windows:
                plist = str(p).replace("\\", "/").split("/")
            else:
                plist = str(p).split("/")
            if not plist[0]:
                plist.pop(0)
                del normed[:]
            for item in plist:
                if item in {"", "."}:
                    continue
                if item == "..":
                    if normed:
                        normed.pop(-1)
                else:
                    normed.append(item)
        if windows:
            self = super().__new__(cls, "\\".join(normed))
            self.CASE_SENSITIVE = False  # On this object only
        else:
            self = super().__new__(cls, "/" + "/".join(normed))
            self.CASE_SENSITIVE = True

        self.remote = remote
        return self

    def _form(self, *parts: str | Path) -> RemotePath:
        return RemotePath(self.remote, *parts)

    @property
    def _path(self) -> str:
        return str(self)

    @property
    def name(self) -> str:
        if "/" not in str(self):
            return str(self)
        return str(self).rsplit("/", 1)[1]

    @property
    def dirname(self) -> RemotePath | str:  # type: ignore[override]
        if "/" not in str(self):
            return str(self)
        return self.__class__(self.remote, str(self).rsplit("/", 1)[0])

    @property
    def suffix(self) -> str:
        return os.path.splitext(self.name)[1]

    @property
    def suffixes(self) -> list[str]:
        exts = []
        name = self.name
        while True:
            name, ext = os.path.splitext(name)
            if ext:
                exts.append(ext)
            else:
                return list(reversed(exts))

    @property
    def uid(self) -> FSUser:
        uid, name = self.remote._path_getuid(self)
        return FSUser(int(uid), name)

    @property
    def gid(self) -> FSUser:
        gid, name = self.remote._path_getgid(self)
        return FSUser(int(gid), name)

    def _get_info(self) -> tuple[BaseRemoteMachine, str]:  # type: ignore[override]
        return (self.remote, self._path)

    def join(self, *parts: str) -> RemotePath:  # type: ignore[override]
        return RemotePath(self.remote, self, *parts)

    def list(self) -> list[RemotePath]:
        if not self.is_dir():
            return []
        return [self.join(fn) for fn in self.remote._path_listdir(self)]

    def iterdir(self) -> Iterable[RemotePath]:
        if not self.is_dir():
            return ()
        return (self.join(fn) for fn in self.remote._path_listdir(self))

    def is_dir(self) -> bool:
        res = self.remote._path_stat(self)
        if not res:
            return False
        return res.text_mode == "directory"

    def is_file(self) -> bool:
        res = self.remote._path_stat(self)
        if not res:
            return False
        return res.text_mode in ("regular file", "regular empty file")

    def is_symlink(self) -> bool:
        res = self.remote._path_stat(self)
        if not res:
            return False
        return res.text_mode == "symbolic link"

    def exists(self) -> bool:
        return self.remote._path_stat(self) is not None

    def stat(self) -> RemoteStatRes:  # type: ignore[override]
        res = self.remote._path_stat(self)
        if res is None:
            raise OSError(errno.ENOENT, os.strerror(errno.ENOENT), "")
        return res

    def with_name(self, name: str) -> RemotePath:
        return self.__class__(self.remote, self.dirname) / name

    def with_suffix(self, suffix: str, depth: int | None = 1) -> RemotePath:
        if (suffix and not suffix.startswith(".")) or suffix == ".":
            raise ValueError(f"Invalid suffix {suffix!r}")
        name = self.name
        depth = len(self.suffixes) if depth is None else min(depth, len(self.suffixes))
        for _ in range(depth):
            name, _ = name.rsplit(".", 1)
        return self.__class__(self.remote, self.dirname) / (name + suffix)

    def glob(self, pattern: str) -> builtins.list[RemotePath]:
        return self._glob(
            pattern,
            lambda pat: [
                RemotePath(self.remote, m) for m in self.remote._path_glob(self, pat)
            ],
        )

    def delete(self) -> None:
        if not self.exists():
            return
        self.remote._path_delete(self)

    def unlink(self) -> None:
        """Removes this file or symbolic link.

        Unlike :meth:`delete`, this refuses to recursively remove a directory
        (matching :func:`os.unlink` / :meth:`LocalPath.unlink
        <plumbum.path.local.LocalPath.unlink>`); call :meth:`delete` for that.
        """
        # Use lstat so a symlink to a directory is still unlinkable and a
        # dangling symlink is still found; only a real directory is refused.
        res = self.remote._path_lstat(self)
        if res is None:
            return
        if res.text_mode == "directory":
            raise OSError(errno.EISDIR, os.strerror(errno.EISDIR), str(self))
        self.remote._path_unlink(self)

    def move(self, dst: RemotePath | str) -> RemotePath:
        if isinstance(dst, RemotePath):
            if dst.remote is not self.remote:
                raise TypeError("dst points to a different remote machine")
        elif not isinstance(dst, str):
            raise TypeError(
                f"dst must be a string or a RemotePath (to the same remote machine), got {dst!r}"
            )
        return self.remote._path_move(self, dst)

    def copy(self, dst: RemotePath | str, override: bool | None = False) -> RemotePath:
        if isinstance(dst, RemotePath):
            if dst.remote is not self.remote:
                raise TypeError("dst points to a different remote machine")
        elif not isinstance(dst, str):
            raise TypeError(
                f"dst must be a string or a RemotePath (to the same remote machine), got {dst!r}"
            )
        if override:
            if isinstance(dst, str):
                dst = RemotePath(self.remote, dst)
            dst.delete()
        else:
            if isinstance(dst, str):
                dst = RemotePath(self.remote, dst)
            if dst.exists():
                raise TypeError("Override not specified and dst exists")

        return self.remote._path_copy(self, dst)

    def mkdir(
        self, mode: int | None = None, parents: bool = True, exist_ok: bool = True
    ) -> None:
        if parents and exist_ok:
            self.remote._path_mkdir(self, mode=mode, minus_p=True)
        else:
            if parents and len(self.parts) > 1:
                self.remote._path_mkdir(self.parent, mode=mode, minus_p=True)
            try:
                self.remote._path_mkdir(self, mode=mode, minus_p=False)
            except ProcessExecutionError as ex:
                if "File exists" not in ex.stderr:
                    raise

                if not exist_ok:
                    raise OSError(
                        errno.EEXIST, "File exists (on remote end)", str(self)
                    ) from None

    def read(self, encoding: str | None = None) -> str | bytes:
        data = self.remote._path_read(self)
        if encoding:
            return data.decode(encoding)
        return data

    def write(self, data: str | bytes, encoding: str | None = None) -> None:
        if encoding:
            assert isinstance(data, str)
            data = data.encode(encoding)
        self.remote._path_write(self, data)

    def touch(self) -> None:
        self.remote._path_touch(str(self))

    def chown(
        self,
        owner: int | str | None = None,
        group: int | str | None = None,
        recursive: bool | None = None,
    ) -> None:
        self.remote._path_chown(
            self, owner, group, self.is_dir() if recursive is None else recursive
        )

    def chmod(self, mode: int) -> None:
        self.remote._path_chmod(mode, self)

    def access(self, mode: int | str = 0) -> bool:
        """Test file existence or permission bits.

        This is a best-effort check: it inspects the file's permission bits but
        does not compare the caller's uid/gid against the file's owner/group, so
        the owner, group and other bits are all considered; every requested
        bit must be granted by some class. A bare existence
        check (``mode=0`` / ``os.F_OK`` / ``"f"``) simply reports whether the
        path exists.
        """
        mode = self._access_mode_to_flags(mode)
        res = self.remote._path_stat(self)
        if res is None:
            return False
        if mode == 0:  # os.F_OK: existence check only
            return True
        mask = res.st_mode & 0x1FF
        merged = (mask >> 6) | (mask >> 3) | mask
        return (merged & mode) == mode

    def link(self, dst: RemotePath | str) -> None:
        if isinstance(dst, RemotePath):
            if dst.remote is not self.remote:
                raise TypeError("dst points to a different remote machine")
        elif not isinstance(dst, str):
            raise TypeError(
                f"dst must be a string or a RemotePath (to the same remote machine), got {dst!r}"
            )
        self.remote._path_link(self, dst, False)

    def symlink(self, dst: RemotePath | str) -> None:
        if isinstance(dst, RemotePath):
            if dst.remote is not self.remote:
                raise TypeError("dst points to a different remote machine")
        elif not isinstance(dst, str):
            raise TypeError(
                f"dst must be a string or a RemotePath (to the same remote machine), got {dst!r}"
            )
        self.remote._path_link(self, dst, True)

    def open(
        self, mode: str = "r", bufsize: int = -1, *, encoding: str | None = None
    ) -> typing.IO[str] | typing.IO[bytes]:
        """
        Opens this path as a file.

        Only works for ParamikoMachine-associated paths for now.
        """
        if encoding is not None:
            raise NotImplementedError(
                "encoding not supported for ParamikoMachine paths"
            )

        if hasattr(self.remote, "sftp") and hasattr(self.remote.sftp, "open"):
            return self.remote.sftp.open(self, mode, bufsize)  # type: ignore[no-any-return]

        raise NotImplementedError(
            "RemotePath.open only works for ParamikoMachine-associated paths for now"
        )

    def as_uri(self, scheme: str = "ssh") -> str:
        suffix = urllib.pathname2url(str(self))
        # TODO: BaseRemoteMachine doesn't have this, but maybe just because it's not typed yet
        return f"{scheme}://{self.remote._fqhost}{suffix}"  # type: ignore[attr-defined]

    @property
    def stem(self) -> str:
        return os.path.splitext(self.name)[0]

    @property
    def root(self) -> str:
        return "/"

    @property
    def drive(self) -> str:
        return ""


class RemoteWorkdir(RemotePath):
    """Remote working directory manipulator"""

    __slots__ = ()

    def __new__(cls, remote: BaseRemoteMachine) -> Self:
        return super().__new__(cls, remote, remote._session.run("pwd")[1].strip())

    __hash__ = None  # type: ignore[assignment]

    def chdir(self, newdir: RemotePath | str) -> Self:
        """Changes the current working directory to the given one"""
        self.remote._session.run(f"cd {shquote(newdir)}")
        if hasattr(self.remote, "_cwd"):
            del self.remote._cwd
        return self.__class__(self.remote)

    def getpath(self) -> RemotePath:
        """Returns the current working directory as a
        `remote path <plumbum.path.remote.RemotePath>` object"""
        return RemotePath(self.remote, self)

    @contextmanager
    def __call__(self, newdir: RemotePath | str) -> Generator[RemotePath, None, None]:
        """A context manager used to ``chdir`` into a directory and then ``chdir`` back to
        the previous location; much like ``pushd``/``popd``.

        :param newdir: The destination director (a string or a
                       :class:`RemotePath <plumbum.path.remote.RemotePath>`)
        """
        prev = self._path
        changed_dir = self.chdir(newdir)
        try:
            yield changed_dir
        finally:
            self.chdir(prev)


__all__ = [
    "RemotePath",
    "RemoteStatRes",
    "RemoteWorkdir",
    "StatRes",
]


def __dir__() -> list[str]:
    return list(__all__)
