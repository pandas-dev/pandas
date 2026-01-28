"""
(c) 2008-2013, holger krekel
"""

from __future__ import annotations


class XSpec:
    """Execution Specification: key1=value1//key2=value2 ...

    * Keys need to be unique within the specification scope
    * Neither key nor value are allowed to contain "//"
    * Keys are not allowed to contain "="
    * Keys are not allowed to start with underscore
    * If no "=value" is given, assume a boolean True value
    """

    # XXX allow customization, for only allow specific key names
    chdir: str | None = None
    dont_write_bytecode: bool | None = None
    execmodel: str | None = None
    id: str | None = None
    installvia: str | None = None
    nice: str | None = None
    popen: bool | None = None
    python: str | None = None
    socket: str | None = None
    ssh: str | None = None
    ssh_config: str | None = None
    vagrant_ssh: str | None = None
    via: str | None = None

    def __init__(self, string: str) -> None:
        self._spec = string
        self.env = {}
        for keyvalue in string.split("//"):
            i = keyvalue.find("=")
            value: str | bool
            if i == -1:
                key, value = keyvalue, True
            else:
                key, value = keyvalue[:i], keyvalue[i + 1 :]
            if key[0] == "_":
                raise AttributeError("%r not a valid XSpec key" % key)
            if key in self.__dict__:
                raise ValueError(f"duplicate key: {key!r} in {string!r}")
            if key.startswith("env:"):
                self.env[key[4:]] = value
            else:
                setattr(self, key, value)

    def __getattr__(self, name: str) -> None | bool | str:
        if name[0] == "_":
            raise AttributeError(name)
        return None

    def __repr__(self) -> str:
        return f"<XSpec {self._spec!r}>"

    def __str__(self) -> str:
        return self._spec

    def __hash__(self) -> int:
        return hash(self._spec)

    def __eq__(self, other: object) -> bool:
        return self._spec == getattr(other, "_spec", None)

    def __ne__(self, other: object) -> bool:
        return self._spec != getattr(other, "_spec", None)

    def _samefilesystem(self) -> bool:
        return self.popen is not None and self.chdir is None
