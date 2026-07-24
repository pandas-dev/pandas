from __future__ import annotations

__lazy_modules__ = {"fnmatch", "functools", "itertools", "operator"}

import fnmatch
import itertools
import operator
import os
import typing
from abc import ABC, abstractmethod
from functools import reduce
from typing import IO, SupportsIndex, TypeVar

if typing.TYPE_CHECKING:
    import builtins
    from collections.abc import Callable, Generator, Iterable, Iterator

    from plumbum._compat.typing import Self

FLAGS = {"f": os.F_OK, "w": os.W_OK, "r": os.R_OK, "x": os.X_OK}

P = TypeVar("P", bound="Path")


class FSUser(int):
    """A special object that represents a file-system user. It derives from ``int``, so it behaves
    just like a number (``uid``/``gid``), but also have a ``.name`` attribute that holds the
    string-name of the user, if given (otherwise ``None``)
    """

    name: str | None

    def __new__(cls, val: int, name: str | None = None) -> Self:
        self = int.__new__(cls, val)
        self.name = name
        return self


class Path(str, ABC):
    """An abstraction over file system paths. This class is abstract, and the two implementations
    are :class:`LocalPath <plumbum.path.local.LocalPath>` and
    :class:`RemotePath <plumbum.path.remote.RemotePath>`.
    """

    __slots__ = ()
    CASE_SENSITIVE = True

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__} {self}>"

    def __truediv__(self, other: str | Path) -> Self:
        """Joins two paths"""
        return self.join(other)

    def __rtruediv__(self, other: str | Path) -> Self:
        """Joins two paths when this path appears on the right-hand side."""
        other = self._form(other)
        return other / self

    @typing.overload
    def __getitem__(self, key: str | Path) -> Self: ...

    @typing.overload
    def __getitem__(self, key: SupportsIndex | slice) -> str: ...

    def __getitem__(self, key: str | Path | SupportsIndex | slice) -> Self | str:
        if isinstance(key, (str, Path)):
            return self / key
        return str(self)[key]

    def __floordiv__(self, expr: str) -> list[Self]:
        """Returns a (possibly empty) list of paths that matched the glob-pattern under this path"""
        return self.glob(expr)

    def __iter__(self) -> Iterator[Self]:
        """Iterate over the files in this directory"""
        return iter(self.list())

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Path):
            return self._get_info() == other._get_info()
        if isinstance(other, str):
            if self.CASE_SENSITIVE:
                return str(self) == other

            return str(self).lower() == other.lower()

        return NotImplemented

    def __ne__(self, other: object) -> bool:
        return not self == other

    def __gt__(self, other: object) -> bool:
        return str(self) > str(other)

    def __ge__(self, other: object) -> bool:
        return str(self) >= str(other)

    def __lt__(self, other: object) -> bool:
        return str(self) < str(other)

    def __le__(self, other: object) -> bool:
        return str(self) <= str(other)

    def __hash__(self) -> int:
        return hash(str(self)) if self.CASE_SENSITIVE else hash(str(self).lower())

    def __bool__(self) -> bool:
        return bool(str(self))

    def __fspath__(self) -> str:
        """Added for Python 3.6 support"""
        return str(self)

    def __contains__(self, item: object) -> bool:
        """Paths should support checking to see if an file or folder is in them."""
        try:
            return (self / item.name).exists()  # type: ignore[attr-defined, no-any-return]
        except AttributeError:
            return (self / item).exists()  # type: ignore[operator]

    @abstractmethod
    def _form(self, *parts: str) -> Self:
        pass

    def up(self, count: int = 1) -> Self:
        """Go up in ``count`` directories (the default is 1)"""
        return self.join("../" * count)

    # Currently, typed as just "str" to match .join, though Path should mostly work too.
    def joinpath(self, *others: str) -> Self:
        """Pathlib-compatible alias of :meth:`join`.

        .. versionadded:: 2.0
        """
        return self.join(*others)

    def walk(
        self,
        filter: Callable[[Self], bool] = lambda _: True,  # pylint: disable=redefined-builtin
        dir_filter: Callable[[Self], bool] = lambda _: True,
    ) -> Generator[Self, None, None]:
        """traverse all (recursive) sub-elements under this directory, that match the given filter.
        By default, the filter accepts everything; you can provide a custom filter function that
        takes a path as an argument and returns a boolean

        :param filter: the filter (predicate function) for matching results. Only paths matching
                       this predicate are returned. Defaults to everything.
        :param dir_filter: the filter (predicate function) for matching directories. Only directories
                           matching this predicate are recursed into. Defaults to everything.
        """
        for p in self.list():
            if filter(p):
                yield p
            if p.is_dir() and dir_filter(p):
                yield from p.walk(filter, dir_filter)

    @property
    @abstractmethod
    def name(self) -> str:
        """The basename component of this path"""

    @property
    @abstractmethod
    def stem(self) -> str:
        """The name without an extension, or the last component of the path"""

    @property
    @abstractmethod
    def dirname(self) -> Self:
        """The dirname component of this path"""

    @property
    @abstractmethod
    def root(self) -> str:
        """The root of the file tree (`/` on Unix)"""

    @property
    @abstractmethod
    def drive(self) -> str:
        """The drive letter (on Windows)"""

    @property
    @abstractmethod
    def suffix(self) -> str:
        """The suffix of this file"""

    @property
    @abstractmethod
    def suffixes(self) -> list[str]:
        """This is a list of all suffixes"""

    @property
    @abstractmethod
    def uid(self) -> FSUser:
        """The user that owns this path. The returned value is a :class:`FSUser <plumbum.path.base.FSUser>`
        object which behaves like an ``int`` (as expected from ``uid``), but it also has a ``.name``
        attribute that holds the string-name of the user"""

    @property
    @abstractmethod
    def gid(self) -> FSUser:
        """The group that owns this path. The returned value is a :class:`FSUser <plumbum.path.base.FSUser>`
        object which behaves like an ``int`` (as expected from ``gid``), but it also has a ``.name``
        attribute that holds the string-name of the group"""

    @abstractmethod
    def as_uri(self, scheme: str = ...) -> str:
        """Returns a universal resource identifier. Use ``scheme`` to force a scheme."""

    @abstractmethod
    def _get_info(self) -> str | tuple[str, str]:
        pass

    @abstractmethod
    def join(self, *parts: str) -> Self:  # type: ignore[override]
        """Joins this path with any number of paths"""

    @abstractmethod
    def list(self) -> builtins.list[Self]:
        """Returns the files in this directory"""

    @abstractmethod
    def iterdir(self) -> Iterable[Self]:
        """Returns an iterator over the directory. Might be slightly faster on Python 3.5 than .list()"""

    @abstractmethod
    def is_dir(self) -> bool:
        """Returns ``True`` if this path is a directory, ``False`` otherwise"""

    @abstractmethod
    def is_file(self) -> bool:
        """Returns ``True`` if this path is a regular file, ``False`` otherwise"""

    @abstractmethod
    def is_symlink(self) -> bool:
        """Returns ``True`` if this path is a symbolic link, ``False`` otherwise"""

    @abstractmethod
    def exists(self) -> bool:
        """Returns ``True`` if this path exists, ``False`` otherwise"""

    @abstractmethod
    def stat(self) -> os.stat_result:
        """Returns the os.stats for a file"""

    @abstractmethod
    def with_name(self, name: str) -> Self:
        """Returns a path with the name replaced"""

    @abstractmethod
    def with_suffix(self, suffix: str, depth: int | None = 1) -> Self:
        """Returns a path with the suffix replaced. Up to last ``depth`` suffixes will be
        replaced. None will replace all suffixes. If there are less than ``depth`` suffixes,
        this will replace all suffixes. ``.tar.gz`` is an example where ``depth=2`` or
        ``depth=None`` is useful"""

    def preferred_suffix(self, suffix: str) -> Self:
        """Adds a suffix if one does not currently exist (otherwise, no change). Useful
        for loading files with a default suffix"""
        return self if len(self.suffixes) > 0 else self.with_suffix(suffix)

    def with_stem(self, stem: str) -> Self:
        """Returns a path with the stem replaced.

        .. versionadded:: 2.0
        """
        return self.with_name(stem + "".join(self.suffixes))

    @abstractmethod
    def glob(self, pattern: str) -> builtins.list[Self]:
        """Returns a (possibly empty) list of paths that matched the glob-pattern under this path"""

    @abstractmethod
    def delete(self) -> None:
        """Deletes this path (recursively, if a directory)"""

    @abstractmethod
    def move(self, dst: Self | str) -> Self:
        """Moves this path to a different location"""

    def rename(self, newname: str) -> Self:
        """Renames this path to the ``new name`` (only the basename is changed)"""
        return self.move(self.up() / newname)

    def rmdir(self) -> None:
        """Removes this directory if it is empty.

        .. versionadded:: 2.0
        """
        if not self.is_dir():
            raise NotADirectoryError(str(self))
        if any(True for _ in self.iterdir()):
            raise OSError(f"Directory not empty: {self}")
        self.delete()

    @abstractmethod
    def copy(self, dst: Self | str, override: bool | None = None) -> Self:
        """Copies this path (recursively, if a directory) to the destination path "dst".
        Raises TypeError if dst exists and override is False.
        Will overwrite if override is True.
        Will silently fail to copy if override is None (the default)."""

    @abstractmethod
    def mkdir(
        self, mode: int = 0o777, parents: bool = True, exist_ok: bool = True
    ) -> None:
        """
        Creates a directory at this path.

        :param mode: **Currently only implemented for local paths!** Numeric mode to use for directory
                     creation, which may be ignored on some systems. The current implementation
                     reproduces the behavior of ``os.mkdir`` (i.e., the current umask is first masked
                     out), but this may change for remote paths. As with ``os.mkdir``, it is recommended
                     to call :func:`chmod` explicitly if you need to be sure.
        :param parents: If this is true (the default), the directory's parents will also be created if
                        necessary.
        :param exist_ok: If this is true (the default), no exception will be raised if the directory
                         already exists (otherwise ``OSError``).

        Note that the defaults for ``parents`` and ``exist_ok`` are the opposite of what they are in
        Python's own ``pathlib`` - this is to maintain backwards-compatibility with Plumbum's behaviour
        from before they were implemented.
        """

    @abstractmethod
    def open(
        self, mode: str = "r", *, encoding: str | None = None
    ) -> IO[str] | IO[bytes]:
        """opens this path as a file"""

    def read_bytes(self) -> bytes:
        """Returns the contents of this file as ``bytes``.

        .. versionadded:: 2.0
        """
        try:
            with self.open("rb") as f:
                data = f.read()
            assert isinstance(data, bytes)
        except NotImplementedError:
            data = self.read()
            if isinstance(data, str):
                return data.encode("utf-8")
            assert isinstance(data, bytes)
            return data
        return data

    def read_text(self, encoding: str | None = None, errors: str | None = None) -> str:
        """Returns the contents of this file as ``str``.

        .. versionadded:: 2.0
        """
        if errors not in {None, "strict"}:
            raise NotImplementedError("Only errors='strict' is currently supported")
        data = self.read(encoding=encoding or "utf-8")
        assert isinstance(data, str)
        return data

    @abstractmethod
    def read(self, encoding: str | None = None) -> str | bytes:
        """returns the contents of this file as a ``str``. By default the data is read
        as text, but you can specify the encoding, e.g., ``'latin1'`` or ``'utf8'``"""

    @abstractmethod
    def write(self, data: str | bytes, encoding: str | None = None) -> None:
        """writes the given data to this file. By default the data is written as-is
        (either text or binary), but you can specify the encoding, e.g., ``'latin1'``
        or ``'utf8'``"""

    def write_bytes(self, data: bytes) -> int:
        """Writes bytes to this file and returns the number of bytes written.

        .. versionadded:: 2.0
        """
        if not isinstance(data, bytes):
            raise TypeError("data must be bytes")
        self.write(data)
        return len(data)

    def write_text(
        self,
        data: str,
        encoding: str | None = None,
        errors: str | None = None,
        newline: str | None = None,
    ) -> int:
        """Writes text to this file and returns the number of characters written.

        .. versionadded:: 2.0
        """
        if errors not in {None, "strict"}:
            raise NotImplementedError("Only errors='strict' is currently supported")
        if newline is not None:
            data = data.replace("\n", newline)
        self.write(data, encoding=encoding or "utf-8")
        return len(data)

    @abstractmethod
    def touch(self) -> None:
        """Update the access time. Creates an empty file if none exists."""

    @abstractmethod
    def chown(
        self,
        owner: int | str | None = None,
        group: int | str | None = None,
        recursive: bool | None = None,
    ) -> None:
        """Change ownership of this path.

        :param owner: The owner to set (either ``uid`` or ``username``), optional
        :param group: The group to set (either ``gid`` or ``groupname``), optional
        :param recursive: whether to change ownership of all contained files and subdirectories.
                          Only meaningful when ``self`` is a directory. If ``None``, the value
                          will default to ``True`` if ``self`` is a directory, ``False`` otherwise.
        """

    @abstractmethod
    def chmod(self, mode: int) -> None:
        """Change the mode of path to the numeric mode.

        :param mode: file mode as for os.chmod
        """

    @staticmethod
    def _access_mode_to_flags(
        mode: int | str, flags: dict[str, int] | None = None
    ) -> int:
        if flags is None:
            flags = FLAGS

        if isinstance(mode, str):
            return reduce(operator.or_, [flags[m] for m in mode.lower()], 0)

        return mode

    @abstractmethod
    def access(self, mode: int | str = 0) -> bool:
        """Test file existence or permission bits

        :param mode: a bitwise-or of access bits, or a string-representation thereof:
                     ``'f'``, ``'x'``, ``'r'``, ``'w'`` for ``os.F_OK``, ``os.X_OK``,
                     ``os.R_OK``, ``os.W_OK``
        """

    @abstractmethod
    def link(self, dst: Self | str) -> None:
        """Creates a hard link from ``self`` to ``dst``

        :param dst: the destination path
        """

    @abstractmethod
    def symlink(self, dst: Self | str) -> None:
        """Creates a symbolic link from ``self`` to ``dst``

        :param dst: the destination path
        """

    def symlink_to(self, target: Self | str) -> None:
        """Pathlib-compatible symbolic-link creation.

        Creates a symbolic link at ``self`` that points to ``target``.

        .. versionadded:: 2.0
        """
        target = self._form(target)
        target.symlink(self)

    def hardlink_to(self, target: Self | str) -> None:
        """Pathlib-compatible hard-link creation.

        Creates a hard link at ``self`` that points to ``target``.

        .. versionadded:: 2.0
        """
        target = self._form(target)
        target.link(self)

    @abstractmethod
    def unlink(self) -> None:
        """Deletes a symbolic link"""

    def lstat(self) -> os.stat_result:
        """Like :meth:`stat`. Backends may override for symlink-specific behavior.

        .. versionadded:: 2.0
        """
        return self.stat()

    def split(self, *_args: object, **_kargs: object) -> builtins.list[str]:
        """Splits the path on directory separators, yielding a list of directories, e.g,
        ``"/var/log/messages"`` will yield ``['var', 'log', 'messages']``.
        """
        parts: list[str] = []
        path = self
        while path != path.dirname:
            parts.append(path.name)
            path = path.dirname
        return parts[::-1]

    @property
    def parts(self) -> tuple[str, ...]:
        """Splits the directory into parts, including the base directory, returns a tuple"""
        return (self.drive + self.root, *self.split())

    @property
    def anchor(self) -> str:
        """Pathlib-compatible anchor (drive + root).

        .. versionadded:: 2.0
        """
        return self.drive + self.root

    def as_posix(self) -> str:
        """Returns the path string with POSIX separators.

        .. versionadded:: 2.0
        """
        return str(self).replace("\\", "/")

    def absolute(self) -> Self:
        """Returns an absolute version of this path.

        .. versionadded:: 2.0
        """
        return self.resolve()

    def is_absolute(self) -> bool:
        """Returns ``True`` if this path is absolute.

        .. versionadded:: 2.0
        """
        path = str(self)
        return path.startswith(("/", "\\")) or (
            len(path) >= 3 and path[1] == ":" and path[2] in ("/", "\\")
        )

    def is_relative_to(self, other: Self | str) -> bool:
        """Returns ``True`` when this path is within ``other``.

        .. versionadded:: 2.0
        """
        other = self._form(other)
        if self.anchor != other.anchor:
            return False

        if self.CASE_SENSITIVE:
            parts = self.split()
            base_parts = other.split()
        else:
            parts = [part.lower() for part in self.split()]
            base_parts = [part.lower() for part in other.split()]
        return len(parts) >= len(base_parts) and parts[: len(base_parts)] == base_parts

    def samefile(self, other: Self | str) -> bool:
        """Returns ``True`` if this path and ``other`` point to the same file.

        .. versionadded:: 2.0
        """
        other = self._form(other)
        st = self.stat()
        other_st = other.stat()
        return (st.st_dev, st.st_ino) == (other_st.st_dev, other_st.st_ino)

    def match(self, pattern: str) -> bool:
        """Tests if this path matches a glob-style pattern.

        .. versionadded:: 2.0
        """
        path = self.as_posix()
        pat = pattern.replace("\\", "/")
        if not self.CASE_SENSITIVE:
            path = path.lower()
            pat = pat.lower()

        if pat.startswith("/") or (len(pat) >= 3 and pat[1] == ":" and pat[2] == "/"):
            return fnmatch.fnmatchcase(path, pat)

        if "/" not in pat:
            return fnmatch.fnmatchcase(path.rsplit("/", 1)[-1], pat)

        pat_parts = [segment for segment in pat.split("/") if segment]
        path_parts = [segment for segment in path.split("/") if segment]
        if len(pat_parts) > len(path_parts):
            return False
        right_side = "/".join(path_parts[-len(pat_parts) :])
        return fnmatch.fnmatchcase(right_side, "/".join(pat_parts))

    def rglob(self, pattern: str) -> builtins.list[Self]:
        """Recursively glob under this path.

        .. versionadded:: 2.0
        """
        return [path for path in self.walk() if path.match(pattern)]

    def relative_to(self, source: Self | str) -> RelativePath:
        """Computes the "relative path" require to get from ``source`` to ``self``. They satisfy the invariant
        ``source_path + (target_path - source_path) == target_path``. For example::

            /var/log/messages - /var/log/messages = []
            /var/log/messages - /var              = [log, messages]
            /var/log/messages - /                 = [var, log, messages]
            /var/log/messages - /var/tmp          = [.., log, messages]
            /var/log/messages - /opt              = [.., var, log, messages]
            /var/log/messages - /opt/lib          = [.., .., var, log, messages]
        """
        source = self._form(source)
        parts = self.split()
        baseparts = source.split()
        ancestors = len(
            list(itertools.takewhile(lambda p: p[0] == p[1], zip(parts, baseparts)))
        )
        return RelativePath([".."] * (len(baseparts) - ancestors) + parts[ancestors:])

    def __sub__(self, other: Self | str) -> RelativePath:
        """Same as ``self.relative_to(other)``"""
        return self.relative_to(other)

    @classmethod
    def _glob(
        cls, pattern: str | Iterable[str], fn: Callable[[str], builtins.list[Self]]
    ) -> builtins.list[Self]:
        """Applies a glob string or list/tuple/iterable to the current path, using ``fn``"""
        if isinstance(pattern, str):
            return fn(pattern)

        results = {value for single_pattern in pattern for value in fn(single_pattern)}
        return sorted(results)

    def resolve(self, strict: bool = False) -> Self:  # noqa: ARG002
        """Added to allow pathlib like syntax. Does nothing since
        Plumbum paths are always absolute. Does not (currently) resolve
        symlinks."""
        # TODO: Resolve symlinks here
        return self

    @property
    def parents(self) -> tuple[Self, ...]:
        """Pathlib like sequence of ancestors"""
        # reduce() infers ``str`` here because Path subclasses str, but each
        # element is really a freshly-formed Self (from ``_form`` / ``/``).
        as_list = (
            typing.cast(
                "Self",
                reduce(lambda x, y: self._form(x) / y, self.parts[:i], self.parts[0]),
            )
            for i in range(len(self.parts) - 1, 0, -1)
        )
        return tuple(as_list)

    @property
    def parent(self) -> Self:
        """Pathlib like parent of the path.

        When called on the root path (e.g. ``/``), returns itself, matching
        the behaviour of :class:`pathlib.Path`.
        """
        parents = self.parents
        return parents[0] if parents else self


class RelativePath:
    """
    Relative paths are the "delta" required to get from one path to another.
    Note that relative path do not point at anything, and thus are not paths.
    Therefore they are system agnostic (but closed under addition)
    Paths are always absolute and point at "something", whether existent or not.

    Relative paths are created by subtracting paths (``Path.relative_to``)
    """

    __slots__ = ("parts",)

    def __init__(self, parts: list[str]):
        self.parts = parts

    def __str__(self) -> str:
        return "/".join(self.parts)

    def __iter__(self) -> Iterator[str]:
        return iter(self.parts)

    def __len__(self) -> int:
        return len(self.parts)

    @typing.overload
    def __getitem__(self, index: SupportsIndex) -> str: ...

    @typing.overload
    def __getitem__(self, index: slice) -> list[str]: ...

    def __getitem__(self, index: SupportsIndex | slice) -> list[str] | str:
        return self.parts[index]

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.parts!r})"

    def __eq__(self, other: object) -> bool:
        return str(self) == str(other)

    def __ne__(self, other: object) -> bool:
        return not self == other

    def __gt__(self, other: object) -> bool:
        return str(self) > str(other)

    def __ge__(self, other: object) -> bool:
        return str(self) >= str(other)

    def __lt__(self, other: object) -> bool:
        return str(self) < str(other)

    def __le__(self, other: object) -> bool:
        return str(self) <= str(other)

    def __hash__(self) -> int:
        return hash(str(self))

    def __bool__(self) -> bool:
        return bool(str(self))

    def up(self, count: int = 1) -> Self:
        return self.__class__(self.parts[:-count])

    def __radd__(self, path: P) -> P:
        return path.join(*self.parts)


__all__ = [
    "FLAGS",
    "FSUser",
    "Path",
    "RelativePath",
]


def __dir__() -> list[str]:
    return list(__all__)
