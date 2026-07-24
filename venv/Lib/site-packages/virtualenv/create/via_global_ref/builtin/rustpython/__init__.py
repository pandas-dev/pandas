from __future__ import annotations

from abc import ABC
from pathlib import Path
from typing import TYPE_CHECKING

from virtualenv.create.describe import PosixSupports, Python3Supports, WindowsSupports
from virtualenv.create.via_global_ref.builtin.ref import RefMust, RefWhen
from virtualenv.create.via_global_ref.builtin.via_global_self_do import ViaGlobalRefVirtualenvBuiltin

if TYPE_CHECKING:
    from collections.abc import Generator

    from python_discovery import PythonInfo


class RustPython(ViaGlobalRefVirtualenvBuiltin, Python3Supports, ABC):
    @classmethod
    def can_describe(cls, interpreter: PythonInfo) -> bool:
        return interpreter.implementation == "RustPython" and super().can_describe(interpreter)

    @classmethod
    def exe_stem(cls) -> str:
        return "rustpython"

    @classmethod
    def exe_names(cls, interpreter: PythonInfo) -> set[str]:
        return {
            cls.exe_stem(),
            "python",
            f"python{interpreter.version_info.major}",
            f"python{interpreter.version_info.major}.{interpreter.version_info.minor}",
        }

    @classmethod
    def _executables(cls, interpreter: PythonInfo) -> Generator[tuple[Path, list[str], RefMust, RefWhen], None, None]:  # ty: ignore[invalid-method-override]
        host = Path(interpreter.system_executable)  # ty: ignore[invalid-argument-type]
        targets = sorted(f"{name}{cls.suffix}" for name in cls.exe_names(interpreter))
        yield host, targets, RefMust.NA, RefWhen.ANY


class RustPythonPosix(RustPython, PosixSupports):
    """RustPython on POSIX."""


class RustPythonWindows(RustPython, WindowsSupports):
    """RustPython on Windows."""


__all__ = [
    "RustPythonPosix",
    "RustPythonWindows",
]
