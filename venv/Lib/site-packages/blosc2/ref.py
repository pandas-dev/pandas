#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True, slots=True)
class Ref:
    """A durable reference to a Blosc2 object.

    ``Ref`` can describe:

    - a persistent local Blosc2 object reopenable from ``urlpath``
    - a member inside a :class:`blosc2.DictStore`
    - a remote :class:`blosc2.C2Array`

    Instances can be created directly, from dictionaries via :meth:`from_dict`,
    or from supported objects via :meth:`from_object`. Use :meth:`open` to
    resolve the reference back into a live Blosc2 object.
    """

    kind: str
    urlpath: str | None = None
    key: str | None = None
    path: str | None = None
    urlbase: str | None = None

    def __post_init__(self) -> None:
        if self.kind == "urlpath":
            if not isinstance(self.urlpath, str):
                raise TypeError("Ref(kind='urlpath') requires a string 'urlpath'")
            if self.key is not None or self.path is not None or self.urlbase is not None:
                raise ValueError("Ref(kind='urlpath') only supports the 'urlpath' field")
            return
        if self.kind == "dictstore_key":
            if not isinstance(self.urlpath, str):
                raise TypeError("Ref(kind='dictstore_key') requires a string 'urlpath'")
            if not isinstance(self.key, str):
                raise TypeError("Ref(kind='dictstore_key') requires a string 'key'")
            if self.path is not None or self.urlbase is not None:
                raise ValueError("Ref(kind='dictstore_key') only supports 'urlpath' and 'key'")
            return
        if self.kind == "c2array":
            if not isinstance(self.path, str):
                raise TypeError("Ref(kind='c2array') requires a string 'path'")
            if self.urlbase is not None and not isinstance(self.urlbase, str):
                raise TypeError("Ref(kind='c2array') requires 'urlbase' to be a string or None")
            if self.urlpath is not None or self.key is not None:
                raise ValueError("Ref(kind='c2array') only supports 'path' and 'urlbase'")
            return
        raise ValueError(f"Unsupported Ref kind: {self.kind!r}")

    @classmethod
    def urlpath_ref(cls, urlpath: str) -> Ref:
        return cls(kind="urlpath", urlpath=urlpath)

    @classmethod
    def dictstore_key(cls, urlpath: str, key: str) -> Ref:
        return cls(kind="dictstore_key", urlpath=urlpath, key=key)

    @classmethod
    def c2array_ref(cls, path: str, urlbase: str | None = None) -> Ref:
        return cls(kind="c2array", path=path, urlbase=urlbase)

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> Ref:
        if not isinstance(payload, dict):
            raise TypeError("Ref payload must be a mapping")
        version = payload.get("version")
        if version != 1:
            raise ValueError(f"Unsupported Ref payload version: {version!r}")
        return cls(
            kind=payload.get("kind"),
            urlpath=payload.get("urlpath"),
            key=payload.get("key"),
            path=payload.get("path"),
            urlbase=payload.get("urlbase"),
        )

    @classmethod
    def from_object(cls, obj: Any) -> Ref:
        import blosc2

        if isinstance(obj, blosc2.C2Array):
            return cls.c2array_ref(obj.path, obj.urlbase)
        if isinstance(obj, blosc2.Proxy):
            obj = obj._cache
        ref = getattr(obj, "_blosc2_ref", None)
        if isinstance(ref, cls):
            return ref
        if hasattr(obj, "schunk"):
            urlpath = obj.schunk.urlpath
            if urlpath is None:
                raise ValueError("Durable Blosc2 references require operands to be stored on disk/network")
            return cls.urlpath_ref(urlpath)
        raise TypeError("Durable Blosc2 references require NDArray, C2Array, or Proxy operands")

    def to_dict(self) -> dict[str, Any]:
        payload = {"kind": self.kind, "version": 1}
        if self.kind == "urlpath":
            payload["urlpath"] = self.urlpath
        elif self.kind == "dictstore_key":
            payload["urlpath"] = self.urlpath
            payload["key"] = self.key
        elif self.kind == "c2array":
            payload["path"] = self.path
            payload["urlbase"] = self.urlbase
        return payload

    def open(self):
        import blosc2

        if self.kind == "urlpath":
            # Structured refs are used to reopen operands for persisted recipes.
            # Read-only access avoids allocating unnecessary writable state.
            return blosc2.open(self.urlpath, mode="r")
        if self.kind == "dictstore_key":
            return blosc2.DictStore(self.urlpath, mode="r")[self.key]
        if self.kind == "c2array":
            return blosc2.C2Array(self.path, urlbase=self.urlbase)
        raise ValueError(f"Unsupported Ref kind: {self.kind!r}")
