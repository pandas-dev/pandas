from __future__ import annotations

from abc import ABC
from typing import TYPE_CHECKING

from virtualenv.create.via_global_ref.api import ViaGlobalRefApi, ViaGlobalRefMeta
from virtualenv.create.via_global_ref.builtin.ref import (
    ExePathRefToDest,
    RefMust,
    RefWhen,
)
from virtualenv.util.path import ensure_dir

from .builtin_way import VirtualenvBuiltin

if TYPE_CHECKING:
    from collections.abc import Generator
    from pathlib import Path

    from python_discovery import PythonInfo

    from virtualenv.config.cli.parser import VirtualEnvOptions
    from virtualenv.create.via_global_ref.builtin.ref import PathRef
    from virtualenv.create.via_global_ref.venv import Venv


class BuiltinViaGlobalRefMeta(ViaGlobalRefMeta):
    def __init__(self) -> None:
        super().__init__()
        self.sources: list[PathRef] = []


class ViaGlobalRefVirtualenvBuiltin(ViaGlobalRefApi, VirtualenvBuiltin, ABC):
    def __init__(self, options: VirtualEnvOptions, interpreter: PythonInfo) -> None:
        super().__init__(options, interpreter)
        self._sources: list[PathRef] = (
            getattr(options.meta, "sources", None) or []
        )  # if created as a describer this might be missing

    @classmethod
    def can_create(cls, interpreter: PythonInfo) -> BuiltinViaGlobalRefMeta | None:
        """By default, all built-in methods assume that if we can describe it we can create it."""
        # first we must be able to describe it
        if not cls.can_describe(interpreter):
            return None
        meta = cls.setup_meta(interpreter)
        if meta is not None and meta:
            cls._sources_can_be_applied(interpreter, meta)
        return meta

    @classmethod
    def _sources_can_be_applied(cls, interpreter: PythonInfo, meta: BuiltinViaGlobalRefMeta) -> None:
        for src in cls.sources(interpreter):
            if src.exists:
                if meta.can_copy and not src.can_copy:
                    meta.copy_error = f"cannot copy {src}"
                if meta.can_symlink and not src.can_symlink:
                    meta.symlink_error = f"cannot symlink {src}"
            else:
                msg = f"missing required file {src}"
                if src.when == RefMust.NA:
                    meta.error = msg
                elif src.when == RefMust.COPY:
                    meta.copy_error = msg
                elif src.when == RefMust.SYMLINK:
                    meta.symlink_error = msg
            if not meta.can_copy and not meta.can_symlink:
                meta.error = f"neither copy or symlink supported, copy: {meta.copy_error} symlink: {meta.symlink_error}"
            if meta.error:
                break
            meta.sources.append(src)

    @classmethod
    def setup_meta(cls, interpreter: PythonInfo) -> BuiltinViaGlobalRefMeta:  # ruff:ignore[unused-class-method-argument]
        return BuiltinViaGlobalRefMeta()

    @classmethod
    def sources(cls, interpreter: PythonInfo) -> Generator[ExePathRefToDest]:
        for host_exe, targets, must, when in cls._executables(interpreter):
            yield ExePathRefToDest(host_exe, dest=cls.to_bin, targets=targets, must=must, when=when)

    def to_bin(self, src: Path) -> Path:
        return self.bin_dir / src.name

    @classmethod
    def _executables(cls, interpreter: PythonInfo) -> Generator[tuple[Path, list[str], str, str]]:
        raise NotImplementedError

    def create(self) -> None:
        dirs = self.ensure_directories()
        for directory in list(dirs):
            if any(i for i in dirs if i is not directory and directory.parts == i.parts[: len(directory.parts)]):
                dirs.remove(directory)
        for directory in sorted(dirs):
            ensure_dir(directory)

        self.set_pyenv_cfg()
        self.pyenv_cfg.write()
        true_system_site = self.enable_system_site_package
        try:
            self.enable_system_site_package = False
            for src in self._sources:
                if (
                    src.when == RefWhen.ANY
                    or (src.when == RefWhen.SYMLINK and self.symlinks is True)
                    or (src.when == RefWhen.COPY and self.symlinks is False)
                ):
                    src.run(self, self.symlinks)
        finally:
            if true_system_site != self.enable_system_site_package:
                self.enable_system_site_package = true_system_site
        super().create()

    @property
    def include_dir(self) -> Path:
        return self.dest / ("Include" if self.interpreter.os == "nt" else "include")

    def install_venv_shared_libs(self, venv_creator: Venv) -> None:
        pass

    def ensure_directories(self) -> set[Path]:
        return {self.dest, self.bin_dir, self.script_dir, self.stdlib, self.include_dir} | set(self.libs)

    def set_pyenv_cfg(self) -> None:
        """We directly inject the base prefix and base exec prefix to avoid site.py needing to discover these from home (which usually is done within the interpreter itself)."""
        super().set_pyenv_cfg()
        self.pyenv_cfg["base-prefix"] = self.interpreter.system_prefix
        self.pyenv_cfg["base-exec-prefix"] = self.interpreter.system_exec_prefix
        self.pyenv_cfg["base-executable"] = self.interpreter.system_executable  # ty: ignore[invalid-assignment]


__all__ = [
    "BuiltinViaGlobalRefMeta",
    "ViaGlobalRefVirtualenvBuiltin",
]
