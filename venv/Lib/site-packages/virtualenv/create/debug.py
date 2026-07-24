"""Inspect a target Python interpreter virtual environment wise."""

from __future__ import annotations

import sys  # built-in


def encode_path(value: object) -> str | None:
    if value is None:
        return None
    if not isinstance(value, (str, bytes)):
        value = repr(value) if isinstance(value, type) else repr(type(value))
    if isinstance(value, bytes):
        value = value.decode(sys.getfilesystemencoding())
    return value


def encode_list_path(value: list[object]) -> list[str | None]:
    return [encode_path(i) for i in value]


def run() -> None:
    """Print debug data about the virtual environment."""
    sys_info: dict[str, str | list[str | None] | None] = {}
    result: dict[str, str | dict[str, str | list[str | None] | None] | None] = {"sys": sys_info}
    path_keys = (
        "executable",
        "_base_executable",
        "prefix",
        "base_prefix",
        "exec_prefix",
        "base_exec_prefix",
        "path",
        "meta_path",
    )
    for key in path_keys:
        value = getattr(sys, key, None)
        value = encode_list_path(value) if isinstance(value, list) else encode_path(value)
        sys_info[key] = value
    sys_info["fs_encoding"] = sys.getfilesystemencoding()
    sys_info["io_encoding"] = getattr(sys.stdout, "encoding", None)
    result["version"] = sys.version

    try:
        import sysconfig  # ruff:ignore[import-outside-top-level]

        result["makefile_filename"] = encode_path(sysconfig.get_makefile_filename())
    except ImportError:
        pass

    import os  # landmark  # ruff:ignore[import-outside-top-level]

    result["os"] = repr(os)

    try:
        import site  # site  # ruff:ignore[import-outside-top-level]

        result["site"] = repr(site)
    except ImportError as exception:  # pragma: no cover
        result["site"] = repr(exception)  # pragma: no cover

    try:
        import datetime  # site  # ruff:ignore[import-outside-top-level]

        result["datetime"] = repr(datetime)
    except ImportError as exception:  # pragma: no cover
        result["datetime"] = repr(exception)  # pragma: no cover

    try:
        import math  # site  # ruff:ignore[import-outside-top-level]

        result["math"] = repr(math)
    except ImportError as exception:  # pragma: no cover
        result["math"] = repr(exception)  # pragma: no cover

    # try to print out, this will validate if other core modules are available (json in this case)
    try:
        import json  # ruff:ignore[import-outside-top-level]

        result["json"] = repr(json)
    except ImportError as exception:
        result["json"] = repr(exception)
    else:
        try:
            content = json.dumps(result, indent=2)
            sys.stdout.write(content)
        except (ValueError, TypeError) as exception:  # pragma: no cover
            sys.stderr.write(repr(exception))
            sys.stdout.write(repr(result))  # pragma: no cover
            raise SystemExit(1)  # ruff:ignore[raise-without-from-inside-except]  # pragma: no cover


if __name__ == "__main__":
    run()
