from __future__ import annotations

__lazy_modules__ = {"configparser", "contextlib"}

import contextlib
from abc import ABC, abstractmethod
from configparser import ConfigParser, NoOptionError, NoSectionError
from typing import TYPE_CHECKING, Any

from plumbum import local

if TYPE_CHECKING:
    from plumbum._compat.typing import Self
    from plumbum.path.local import LocalPath


class ConfigBase(ABC):
    """Base class for Config parsers.

    :param filename: The file to use

    The ``with`` statement can be used to automatically try to read on entering and write if changed on exiting. Otherwise, use ``.read`` and ``.write`` as needed. Set and get the options using ``[]`` syntax.

    Usage:

        with Config("~/.myprog_rc") as conf:
            value = conf.get("option", "default")
            value2 = conf["option"] # shortcut for default=None

    """

    __slots__ = ["changed", "filename"]

    def __init__(self, filename: str):
        self.filename: LocalPath = local.path(filename)
        self.changed = False

    def __enter__(self) -> Self:
        with contextlib.suppress(FileNotFoundError):
            self.read()
        return self

    def __exit__(self, exc_type: object, exc_val: object, exc_tb: object) -> None:
        if self.changed:
            self.write()

    @abstractmethod
    def read(self) -> None:
        """Read in the linked file"""

    @abstractmethod
    def write(self) -> None:
        """Write out the linked file"""
        self.changed = False

    @abstractmethod
    def _get(self, option: str) -> Any:
        """Internal get function for subclasses"""

    @abstractmethod
    def _set(self, option: str, value: Any) -> str:
        """Internal set function for subclasses. Must return the value that was set."""

    def get(self, option: str, default: Any = None) -> Any:
        "Get an item from the store, returns default if fails"
        try:
            return self._get(option)
        except KeyError:
            self.changed = True
            return self._set(option, default)

    def set(self, option: str, value: Any) -> None:
        """Set an item, mark this object as changed"""
        self.changed = True
        self._set(option, value)

    def __getitem__(self, option: str) -> Any:
        return self._get(option)

    def __setitem__(self, option: str, value: Any) -> None:
        return self.set(option, value)


class ConfigINI(ConfigBase):
    DEFAULT_SECTION = "DEFAULT"
    __slots__ = ("parser",)

    def __init__(self, filename: str):
        super().__init__(filename)
        self.parser = ConfigParser()

    def read(self) -> None:
        self.parser.read(self.filename)

    def write(self) -> None:
        with open(self.filename, "w", encoding="utf-8") as f:
            self.parser.write(f)
        super().write()

    @classmethod
    def _sec_opt(cls, option: str) -> tuple[str, str]:
        if "." not in option:
            sec = cls.DEFAULT_SECTION
        else:
            sec, option = option.split(".", 1)
        return sec, option

    def _get(self, option: str) -> Any:
        sec, option = self._sec_opt(option)

        try:
            return self.parser.get(sec, option)
        except (NoSectionError, NoOptionError):
            raise KeyError(f"{sec}:{option}") from None

    def _set(self, option: str, value: Any) -> str:
        sec, option = self._sec_opt(option)
        try:
            self.parser.set(sec, option, str(value))
        except NoSectionError:
            self.parser.add_section(sec)
            self.parser.set(sec, option, str(value))
        return str(value)


Config = ConfigINI

__all__ = [
    "Config",
    "ConfigBase",
    "ConfigINI",
]


def __dir__() -> list[str]:
    return list(__all__)
