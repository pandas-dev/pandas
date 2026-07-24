from __future__ import annotations

from contextlib import suppress
from typing import TYPE_CHECKING

from .convert import convert

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any

    from .convert import TypeData


def get_env_var(key: str, as_type: TypeData, env: Mapping[str, str]) -> tuple[Any, str] | None:
    """Get the environment variable option.

    :param key: the config key requested
    :param as_type: the type we would like to convert it to
    :param env: environment variables to use

    :returns: the converted value and source, or None if not set

    """
    environ_key = f"VIRTUALENV_{key.upper()}"
    if env.get(environ_key):
        value = env[environ_key]

        with suppress(Exception):  # note the converter already logs a warning when failures happen
            source = f"env var {environ_key}"
            as_type = convert(value, as_type, source)
            return as_type, source
    return None


__all__ = [
    "get_env_var",
]
