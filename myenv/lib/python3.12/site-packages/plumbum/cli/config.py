import contextlib
from abc import ABC, abstractmethod
from configparser import ConfigParser, NoOptionError, NoSectionError

from plumbum import local


class ConfigBase(ABC):
    """Base class for Config parsers.

    :param filename: The file to use

    The ``with`` statement can be used to automatically try to read on entering and write if changed on exiting. Otherwise, use ``.read`` and ``.write`` as needed. Set and get the options using ``[]`` syntax.

    Usage:

        with Config("~/.myprog_rc") as conf:
            value = conf.get("option", "default")
            value2 = conf["option"] # shortcut for default=None

    """

    __slots__ = "filename changed".split()

    def __init__(self, filename):
        self.filename = local.path(filename)
        self.changed = False

    def __enter__(self):
        with contextlib.suppress(FileNotFoundError):
            self.read()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.changed:
            self.write()

    @abstractmethod
    def read(self):
        """Read in the linked file"""

    @abstractmethod
    def write(self):
        """Write out the linked file"""
        self.changed = False

    @abstractmethod
    def _get(self, option):
        """Internal get function for subclasses"""

    @abstractmethod
    def _set(self, option, value):
        """Internal set function for subclasses. Must return the value that was set."""

    def get(self, option, default=None):
        "Get an item from the store, returns default if fails"
        try:
            return self._get(option)
        except KeyError:
            self.changed = True
            return self._set(option, default)

    def set(self, option, value):
        """Set an item, mark this object as changed"""
        self.changed = True
        self._set(option, value)

    def __getitem__(self, option):
        return self._get(option)

    def __setitem__(self, option, value):
        return self.set(option, value)


class ConfigINI(ConfigBase):
    DEFAULT_SECTION = "DEFAULT"
    slots = "parser".split()

    def __init__(self, filename):
        super().__init__(filename)
        self.parser = ConfigParser()

    def read(self):
        self.parser.read(self.filename)
        super().read()

    def write(self):
        with open(self.filename, "w", encoding="utf-8") as f:
            self.parser.write(f)
        super().write()

    @classmethod
    def _sec_opt(cls, option):
        if "." not in option:
            sec = cls.DEFAULT_SECTION
        else:
            sec, option = option.split(".", 1)
        return sec, option

    def _get(self, option):
        sec, option = self._sec_opt(option)

        try:
            return self.parser.get(sec, option)
        except (NoSectionError, NoOptionError):
            raise KeyError(f"{sec}:{option}") from None

    def _set(self, option, value):
        sec, option = self._sec_opt(option)
        try:
            self.parser.set(sec, option, str(value))
        except NoSectionError:
            self.parser.add_section(sec)
            self.parser.set(sec, option, str(value))
        return str(value)


Config = ConfigINI
