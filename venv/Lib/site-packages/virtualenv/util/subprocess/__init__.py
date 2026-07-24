from __future__ import annotations

import subprocess
from shlex import quote
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Mapping

CREATE_NO_WINDOW = 0x80000000


class LogCmd:
    def __init__(self, cmd: list[str], env: Mapping[str, str] | None = None) -> None:
        self.cmd = cmd
        self.env = env

    def __repr__(self) -> str:
        cmd_repr = " ".join(quote(str(c)) for c in self.cmd)
        if self.env is not None:
            cmd_repr = f"{cmd_repr} env of {self.env!r}"
        return cmd_repr


def run_cmd(cmd: list[str]) -> tuple[int, str, str]:
    try:
        process = subprocess.Popen(
            cmd,
            universal_newlines=True,
            stdin=subprocess.PIPE,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            encoding="utf-8",
        )
        out, err = process.communicate()  # input disabled
        code = process.returncode
    except OSError as error:
        code, out, err = error.errno, "", error.strerror
        if code == 2 and err is not None and "file" in err:  # ruff:ignore[magic-value-comparison]
            err = str(error)  # FileNotFoundError in Python >= 3.3
    return code, out, err  # ty: ignore[invalid-return-type]


__all__ = (
    "CREATE_NO_WINDOW",
    "LogCmd",
    "run_cmd",
)
