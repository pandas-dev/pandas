from __future__ import annotations

import os
from contextlib import contextmanager


class EnvPathList(list):
    __slots__ = ["_path_factory", "_pathsep", "__weakref__"]

    def __init__(self, path_factory, pathsep):
        super().__init__()
        self._path_factory = path_factory
        self._pathsep = pathsep

    def append(self, path):
        list.append(self, self._path_factory(path))

    def extend(self, paths):
        list.extend(self, (self._path_factory(p) for p in paths))

    def insert(self, index, path):
        list.insert(self, index, self._path_factory(path))

    def index(self, path):
        list.index(self, self._path_factory(path))

    def __contains__(self, path):
        return list.__contains__(self, self._path_factory(path))

    def remove(self, path):
        list.remove(self, self._path_factory(path))

    def update(self, text):
        self[:] = [self._path_factory(p) for p in text.split(self._pathsep)]

    def join(self):
        return self._pathsep.join(str(p) for p in self)


class BaseEnv:
    """The base class of LocalEnv and RemoteEnv"""

    __slots__ = ["_curr", "_path", "_path_factory", "__weakref__"]
    CASE_SENSITIVE = True

    def __init__(self, path_factory, pathsep, *, _curr):
        self._curr = _curr
        self._path_factory = path_factory
        self._path = EnvPathList(path_factory, pathsep)
        self._update_path()

    def _update_path(self):
        self._path.update(self.get("PATH", ""))

    @contextmanager
    def __call__(self, *args, **kwargs):
        """A context manager that can be used for temporal modifications of the environment.
        Any time you enter the context, a copy of the old environment is stored, and then restored,
        when the context exits.

        :param args: Any positional arguments for ``update()``
        :param kwargs: Any keyword arguments for ``update()``
        """
        prev = self._curr.copy()
        self.update(**kwargs)
        try:
            yield
        finally:
            self._curr = prev
            self._update_path()

    def __iter__(self):
        """Returns an iterator over the items ``(key, value)`` of current environment
        (like dict.items)"""
        return iter(self._curr.items())

    def __hash__(self):
        raise TypeError("unhashable type")

    def __len__(self):
        """Returns the number of elements of the current environment"""
        return len(self._curr)

    def __contains__(self, name):
        """Tests whether an environment variable exists in the current environment"""
        return (name if self.CASE_SENSITIVE else name.upper()) in self._curr

    def __getitem__(self, name):
        """Returns the value of the given environment variable from current environment,
        raising a ``KeyError`` if it does not exist"""
        return self._curr[name if self.CASE_SENSITIVE else name.upper()]

    def keys(self):
        """Returns the keys of the current environment (like dict.keys)"""
        return self._curr.keys()

    def items(self):
        """Returns the items of the current environment (like dict.items)"""
        return self._curr.items()

    def values(self):
        """Returns the values of the current environment (like dict.values)"""
        return self._curr.values()

    def get(self, name, *default):
        """Returns the keys of the current environment (like dict.keys)"""
        return self._curr.get((name if self.CASE_SENSITIVE else name.upper()), *default)

    def __delitem__(self, name):
        """Deletes an environment variable from the current environment"""
        name = name if self.CASE_SENSITIVE else name.upper()
        del self._curr[name]
        if name == "PATH":
            self._update_path()

    def __setitem__(self, name, value):
        """Sets/replaces an environment variable's value in the current environment"""
        name = name if self.CASE_SENSITIVE else name.upper()
        self._curr[name] = value
        if name == "PATH":
            self._update_path()

    def pop(self, name, *default):
        """Pops an element from the current environment (like dict.pop)"""
        name = name if self.CASE_SENSITIVE else name.upper()
        res = self._curr.pop(name, *default)
        if name == "PATH":
            self._update_path()
        return res

    def clear(self):
        """Clears the current environment (like dict.clear)"""
        self._curr.clear()
        self._update_path()

    def update(self, *args, **kwargs):
        """Updates the current environment (like dict.update)"""
        self._curr.update(*args, **kwargs)
        if not self.CASE_SENSITIVE:
            for k, v in list(self._curr.items()):
                self._curr[k.upper()] = v
        self._update_path()

    def getdict(self):
        """Returns the environment as a real dictionary"""
        self._curr["PATH"] = self.path.join()
        return {k: str(v) for k, v in self._curr.items()}

    @property
    def path(self):
        """The system's ``PATH`` (as an easy-to-manipulate list)"""
        return self._path

    def _get_home(self):
        if "HOME" in self:
            return self._path_factory(self["HOME"])
        if "USERPROFILE" in self:  # pragma: no cover
            return self._path_factory(self["USERPROFILE"])
        if "HOMEPATH" in self:  # pragma: no cover
            return self._path_factory(self.get("HOMEDRIVE", ""), self["HOMEPATH"])
        return None

    def _set_home(self, p):
        if "HOME" in self:
            self["HOME"] = str(p)
        elif "USERPROFILE" in self:  # pragma: no cover
            self["USERPROFILE"] = str(p)
        elif "HOMEPATH" in self:  # pragma: no cover
            self["HOMEPATH"] = str(p)
        else:  # pragma: no cover
            self["HOME"] = str(p)

    home = property(_get_home, _set_home)
    """Get or set the home path"""

    @property
    def user(self):
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
