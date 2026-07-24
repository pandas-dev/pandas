from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING

from virtualenv.create.describe import PosixSupports, WindowsSupports
from virtualenv.create.via_global_ref.builtin.ref import PathRefToDest, RefMust, RefWhen
from virtualenv.create.via_global_ref.builtin.via_global_self_do import ViaGlobalRefVirtualenvBuiltin

if TYPE_CHECKING:
    from collections.abc import Generator, Iterator

    from python_discovery import PythonInfo


class GraalPy(ViaGlobalRefVirtualenvBuiltin, ABC):
    @classmethod
    @abstractmethod
    def _native_lib(cls, lib_dir: Path, platform: str) -> Path:
        """Return the path to the native library for this platform."""
        raise NotImplementedError

    @classmethod
    def can_describe(cls, interpreter: PythonInfo) -> bool:
        return interpreter.implementation == "GraalPy" and super().can_describe(interpreter)

    @classmethod
    def exe_stem(cls) -> str:
        return "graalpy"

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

    @classmethod
    def sources(cls, interpreter: PythonInfo) -> Generator[PathRefToDest]:  # ty: ignore[invalid-method-override]
        yield from super().sources(interpreter)
        python_dir = Path(interpreter.system_executable).resolve().parent  # ty: ignore[invalid-argument-type]
        if python_dir.name in {"bin", "Scripts"}:
            python_dir = python_dir.parent

        native_lib = cls._native_lib(python_dir / "lib", interpreter.platform)
        if native_lib.exists():
            yield PathRefToDest(native_lib, dest=lambda self, s: self.bin_dir.parent / "lib" / s.name)

        for jvm_dir_name in ("jvm", "jvmlibs", "modules"):
            jvm_dir = python_dir / jvm_dir_name
            if jvm_dir.exists():
                yield PathRefToDest(jvm_dir, dest=lambda self, s: self.bin_dir.parent / s.name)

    @classmethod
    def _shared_libs(cls, python_dir: Path) -> Iterator[Path]:
        raise NotImplementedError

    def set_pyenv_cfg(self) -> None:
        super().set_pyenv_cfg()
        # GraalPy 24.0 and older had home without the bin
        version = self.interpreter.version_info
        if version.minor <= 10:  # ruff:ignore[magic-value-comparison]
            home = Path(self.pyenv_cfg["home"])
            if home.name == "bin":
                self.pyenv_cfg["home"] = str(home.parent)


class GraalPyPosix(GraalPy, PosixSupports):
    @classmethod
    def _native_lib(cls, lib_dir: Path, platform: str) -> Path:
        if platform == "darwin":
            return lib_dir / "libpythonvm.dylib"
        return lib_dir / "libpythonvm.so"


class GraalPyWindows(GraalPy, WindowsSupports):
    @classmethod
    def _native_lib(cls, lib_dir: Path, platform: str) -> Path:  # ruff:ignore[unused-class-method-argument]
        return lib_dir / "pythonvm.dll"

    def set_pyenv_cfg(self) -> None:
        # GraalPy needs an additional entry in pyvenv.cfg on Windows
        super().set_pyenv_cfg()
        self.pyenv_cfg["venvlauncher_command"] = self.interpreter.system_executable  # ty: ignore[invalid-assignment]


__all__ = [
    "GraalPyPosix",
    "GraalPyWindows",
]
