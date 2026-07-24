from __future__ import annotations

__lazy_modules__ = {"contextlib", "copy"}

import copy
import os
import typing
from contextlib import contextmanager

if typing.TYPE_CHECKING:
    from collections.abc import (
        Generator,
        ItemsView,
        Iterable,
        Iterator,
        KeysView,
        MutableMapping,
        ValuesView,
    )
    from typing import Any, Callable, SupportsIndex

    from plumbum.path.base import Path

AnyPath = typing.TypeVar("AnyPath", bound="Path")
AnyPath_co = typing.TypeVar("AnyPath_co", bound="Path", covariant=True)


class EnvPathList(list[AnyPath], typing.Generic[AnyPath]):
    __slots__ = ("__weakref__", "_path_factory", "_pathsep")

    def __init__(self, path_factory: Callable[[str], AnyPath], pathsep: str) -> None:
        super().__init__()
        self._path_factory = path_factory
        self._pathsep: str = pathsep

    def append(self, path: str) -> None:
        list.append(self, self._path_factory(path))

    def extend(self, paths: Iterable[str]) -> None:
        list.extend(self, (self._path_factory(p) for p in paths))

    def insert(self, index: SupportsIndex, path: str) -> None:
        list.insert(self, index, self._path_factory(path))

    def index(self, path: str) -> int:  # type: ignore[override]
        return list.index(self, self._path_factory(path))

    def __contains__(self, path: object) -> bool:
        return list.__contains__(self, self._path_factory(path))  # type: ignore[arg-type]

    def remove(self, path: str) -> None:
        list.remove(self, self._path_factory(path))

    def update(self, text: str) -> None:
        self[:] = [self._path_factory(p) for p in text.split(self._pathsep)]

    def join(self) -> str:
        return self._pathsep.join(str(p) for p in self)


class _PathFactory(typing.Protocol[AnyPath_co]):
    def __call__(self, *args: str) -> AnyPath_co: ...


class BaseEnv(typing.Generic[AnyPath]):
    """The base class of LocalEnv and RemoteEnv"""

    __slots__ = ("__weakref__", "_curr", "_path", "_path_factory")
    CASE_SENSITIVE = True

    def __init__(
        self,
        path_factory: _PathFactory[AnyPath],
        pathsep: str,
        *,
        _curr: MutableMapping[str, str],
    ) -> None:
        self._curr = _curr
        self._path_factory = path_factory
        self._path = EnvPathList[AnyPath](path_factory, pathsep)
        self._update_path()

    def _update_path(self) -> None:
        self._path.update(self.get("PATH", ""))

    @contextmanager
    def __call__(self, *args: Any, **kwargs: Any) -> Generator[None, None, None]:
        """A context manager that can be used for temporal modifications of the environment.
        Any time you enter the context, a copy of the old environment is stored, and then restored,
        when the context exits.

        :param args: Any positional arguments for ``update()``
        :param kwargs: Any keyword arguments for ``update()``
        """
        prev = copy.copy(self._curr)
        self.update(*args, **kwargs)
        try:
            yield
        finally:
            self._curr = prev
            self._update_path()

    def __iter__(self) -> Iterator[tuple[str, Any]]:
        """Returns an iterator over the items ``(key, value)`` of current environment
        (like dict.items)"""
        return iter(self._curr.items())

    __hash__ = None  # type: ignore[assignment]

    def __len__(self) -> int:
        """Returns the number of elements of the current environment"""
        return len(self._curr)

    def __contains__(self, name: str) -> bool:
        """Tests whether an environment variable exists in the current environment"""
        return (name if self.CASE_SENSITIVE else name.upper()) in self._curr

    def __getitem__(self, name: str) -> str:
        """Returns the value of the given environment variable from current environment,
        raising a ``KeyError`` if it does not exist"""
        return self._curr[name if self.CASE_SENSITIVE else name.upper()]

    def keys(self) -> KeysView[str]:
        """Returns the keys of the current environment (like dict.keys)"""
        return self._curr.keys()

    def items(self) -> ItemsView[str, str]:
        """Returns the items of the current environment (like dict.items)"""
        return self._curr.items()

    def values(self) -> ValuesView[str]:
        """Returns the values of the current environment (like dict.values)"""
        return self._curr.values()

    @typing.overload
    def get(self, name: str, default: None = ...) -> str | None: ...

    @typing.overload
    def get(self, name: str, default: str) -> str: ...

    def get(self, name: str, default: str | None = None) -> str | None:
        """Returns the keys of the current environment (like dict.keys)"""
        return self._curr.get((name if self.CASE_SENSITIVE else name.upper()), default)

    def __delitem__(self, name: str) -> None:
        """Deletes an environment variable from the current environment"""
        name = name if self.CASE_SENSITIVE else name.upper()
        del self._curr[name]
        if name == "PATH":
            self._update_path()

    def __setitem__(self, name: str, value: str) -> None:
        """Sets/replaces an environment variable's value in the current environment"""
        name = name if self.CASE_SENSITIVE else name.upper()
        self._curr[name] = value
        if name == "PATH":
            self._update_path()

    def pop(self, name: str, *default: str) -> str | None:
        """Pops an element from the current environment (like dict.pop)"""
        name = name if self.CASE_SENSITIVE else name.upper()
        res = self._curr.pop(name, *default)
        if name == "PATH":
            self._update_path()
        return res

    def clear(self) -> None:
        """Clears the current environment (like dict.clear)"""
        self._curr.clear()
        self._update_path()

    def update(self, *args: Any, **kwargs: Any) -> None:
        """Updates the current environment (like dict.update)"""
        self._curr.update(*args, **kwargs)
        if not self.CASE_SENSITIVE:
            for k, v in list(self._curr.items()):
                self._curr[k.upper()] = v
        self._update_path()

    def getdict(self) -> dict[str, str]:
        """Returns the environment as a real dictionary"""
        self._curr["PATH"] = self.path.join()
        return {k: str(v) for k, v in self._curr.items()}

    @property
    def path(self) -> EnvPathList[AnyPath]:
        """The system's ``PATH`` (as an easy-to-manipulate list)"""
        return self._path

    @property
    def home(self) -> Path | None:
        """Get or set the home path"""
        if "HOME" in self:
            return self._path_factory(self["HOME"])
        if "USERPROFILE" in self:  # pragma: no cover
            return self._path_factory(self["USERPROFILE"])
        if "HOMEPATH" in self:  # pragma: no cover
            return self._path_factory(self.get("HOMEDRIVE", ""), self["HOMEPATH"])
        return None

    @home.setter
    def home(self, p: Path) -> None:
        if "HOME" in self:
            self["HOME"] = str(p)
        elif "USERPROFILE" in self:  # pragma: no cover
            self["USERPROFILE"] = str(p)
        elif "HOMEPATH" in self:  # pragma: no cover
            self["HOMEPATH"] = str(p)
        else:  # pragma: no cover
            self["HOME"] = str(p)

    @property
    def user(self) -> str | None:
        """Return the user name, or ``None`` if it is not set"""
        # adapted from getpass.getuser()
        for name in ("LOGNAME", "USER", "LNAME", "USERNAME"):  # pragma: no branch
            if name in self:
                return self[name]
        try:
            # POSIX only
            import pwd
        except ImportError:
            return None
        return pwd.getpwuid(os.getuid())[0]  # @UndefinedVariable


__all__ = [
    "BaseEnv",
    "EnvPathList",
]


def __dir__() -> list[str]:
    return list(__all__)
