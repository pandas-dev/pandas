"""Yaml utilities."""
from __future__ import annotations

from pathlib import Path, PurePath
from typing import Any

from yaml import dump as ydump
from yaml import load as yload

try:
    from yaml import CSafeDumper as SafeDumper
    from yaml import CSafeLoader as SafeLoader
except ImportError:  # pragma: no cover
    from yaml import SafeDumper, SafeLoader  # type:ignore[assignment]


def loads(stream: Any) -> Any:
    """Load yaml from a stream."""
    return yload(stream, Loader=SafeLoader)


def dumps(stream: Any) -> str:
    """Parse the first YAML document in a stream as an object."""
    return ydump(stream, Dumper=SafeDumper)


def load(fpath: str | PurePath) -> Any:
    """Load yaml from a file."""
    # coerce PurePath into Path, then read its contents
    data = Path(str(fpath)).read_text(encoding="utf-8")
    return loads(data)


def dump(data: Any, outpath: str | PurePath) -> None:
    """Parse the a YAML document in a file as an object."""
    Path(outpath).write_text(dumps(data), encoding="utf-8")
