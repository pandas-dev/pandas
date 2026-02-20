from __future__ import annotations

import sys

from pathlib import Path
from typing import TYPE_CHECKING
from typing import Any
from typing import Callable
from typing import Dict
from typing import TypedDict
from typing import cast

if sys.version_info >= (3, 11):
    from tomllib import loads as load_toml
else:
    from tomli import loads as load_toml

if TYPE_CHECKING:
    if sys.version_info >= (3, 10):
        from typing import TypeAlias
    else:
        from typing_extensions import TypeAlias

from .. import _log

log = _log.log.getChild("toml")

TOML_RESULT: TypeAlias = Dict[str, Any]
TOML_LOADER: TypeAlias = Callable[[str], TOML_RESULT]


class InvalidTomlError(ValueError):
    """Raised when TOML data cannot be parsed."""


def read_toml_content(path: Path, default: TOML_RESULT | None = None) -> TOML_RESULT:
    try:
        data = path.read_text(encoding="utf-8")
    except FileNotFoundError:
        if default is None:
            raise
        else:
            log.debug("%s missing, presuming default %r", path, default)
            return default
    else:
        try:
            return load_toml(data)
        except Exception as e:  # tomllib/tomli raise different decode errors
            raise InvalidTomlError(f"Invalid TOML in {path}") from e


class _CheatTomlData(TypedDict):
    cheat: dict[str, Any]


def load_toml_or_inline_map(data: str | None) -> dict[str, Any]:
    """
    load toml data - with a special hack if only a inline map is given
    """
    if not data:
        return {}
    try:
        if data[0] == "{":
            data = "cheat=" + data
            loaded: _CheatTomlData = cast(_CheatTomlData, load_toml(data))
            return loaded["cheat"]
        return load_toml(data)
    except Exception as e:  # tomllib/tomli raise different decode errors
        raise InvalidTomlError("Invalid TOML content") from e
