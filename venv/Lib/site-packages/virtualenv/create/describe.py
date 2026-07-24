from __future__ import annotations

from abc import ABC
from collections import OrderedDict
from pathlib import Path
from typing import TYPE_CHECKING

from virtualenv.info import IS_WIN

if TYPE_CHECKING:
    from typing import Any

    from python_discovery import PythonInfo


class Describe:
    """Given a host interpreter tell us information about what the created interpreter might look like."""

    suffix = ".exe" if IS_WIN else ""

    def __init__(self, dest: Path, interpreter: PythonInfo) -> None:
        self.interpreter = interpreter
        self.dest = dest
        self._stdlib = None
        self._stdlib_platform = None
        self._system_stdlib = None
        self._conf_vars = None

    @property
    def bin_dir(self) -> Path:
        return self.script_dir

    @property
    def script_dir(self) -> Path:
        return self.dest / self.interpreter.install_path("scripts")

    @property
    def purelib(self) -> Path:
        return self.dest / self.interpreter.install_path("purelib")

    @property
    def platlib(self) -> Path:
        return self.dest / self.interpreter.install_path("platlib")

    @property
    def libs(self) -> list[Path]:
        return list(OrderedDict(((self.platlib, None), (self.purelib, None))).keys())

    @property
    def stdlib(self) -> Path:
        if self._stdlib is None:
            self._stdlib = Path(self.interpreter.sysconfig_path("stdlib", config_var=self._config_vars))
        return self._stdlib

    @property
    def stdlib_platform(self) -> Path:
        if self._stdlib_platform is None:
            self._stdlib_platform = Path(self.interpreter.sysconfig_path("platstdlib", config_var=self._config_vars))
        return self._stdlib_platform

    @property
    def _config_vars(self) -> dict[str, Any]:
        if self._conf_vars is None:
            self._conf_vars = self._calc_config_vars(self.dest)
        return self._conf_vars

    def _calc_config_vars(self, to: Path) -> dict[str, Any]:
        sys_vars = self.interpreter.sysconfig_vars
        return {
            k: (to if isinstance(v, str) and v.startswith(self.interpreter.prefix) else v) for k, v in sys_vars.items()
        }

    @classmethod
    def can_describe(cls, interpreter: PythonInfo) -> bool:  # ruff:ignore[unused-class-method-argument]
        """Knows means it knows how the output will look."""
        return True

    @property
    def env_name(self) -> str:
        return self.dest.parts[-1]

    @property
    def exe(self) -> Path:
        return self.bin_dir / f"{self.exe_stem()}{self.suffix}"

    @classmethod
    def exe_stem(cls) -> str:
        """Executable name without suffix - there seems to be no standard way to get this without creating it."""
        raise NotImplementedError

    def script(self, name: str) -> Path:
        return self.script_dir / f"{name}{self.suffix}"


class Python3Supports(Describe, ABC):
    pass


class PosixSupports(Describe, ABC):
    @classmethod
    def can_describe(cls, interpreter: PythonInfo) -> bool:
        return interpreter.os == "posix" and super().can_describe(interpreter)


class WindowsSupports(Describe, ABC):
    @classmethod
    def can_describe(cls, interpreter: PythonInfo) -> bool:
        return interpreter.os == "nt" and super().can_describe(interpreter)


__all__ = [
    "Describe",
    "PosixSupports",
    "Python3Supports",
    "WindowsSupports",
]
