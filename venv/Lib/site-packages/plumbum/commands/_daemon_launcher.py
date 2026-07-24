from __future__ import annotations

__lazy_modules__ = {"pickle", "traceback"}

import os
import pickle
import sys
import traceback
import typing

if typing.TYPE_CHECKING:
    from plumbum.commands.base import BaseCommand


def _main() -> int:
    if len(sys.argv) != 2:
        return 2

    response_fd = int(sys.argv[1])
    payload = pickle.loads(sys.stdin.buffer.read())

    command: BaseCommand = payload["command"]
    cwd: str = payload["cwd"]
    stdout_path: str = payload["stdout"]
    stderr_path: str = payload["stderr"]
    append: bool = payload["append"]

    mode = "a" if append else "w"
    rc = 0
    try:
        with (
            open(os.devnull, encoding="utf-8") as stdin_file,
            open(stdout_path, mode, encoding="utf-8") as stdout_file,
            open(stderr_path, mode, encoding="utf-8") as stderr_file,
        ):
            proc = command.popen(
                cwd=cwd,
                close_fds=True,
                stdin=stdin_file.fileno(),
                stdout=stdout_file.fileno(),
                stderr=stderr_file.fileno(),
                start_new_session=True,
            )
            os.write(response_fd, str(proc.pid).encode("utf8"))
    except Exception:
        rc = 1
        tbtext = "".join(traceback.format_exception(*sys.exc_info()))[-16384:]
        os.write(response_fd, tbtext.encode("utf8"))
    finally:
        os.close(response_fd)

    return rc


__all__: list[str] = ["_main"]


def __dir__() -> list[str]:
    return list(__all__)


if __name__ == "__main__":
    raise SystemExit(_main())
