from __future__ import annotations

import logging
from abc import ABC
from pathlib import Path
from typing import TYPE_CHECKING

from virtualenv.create.creator import Creator, CreatorMeta
from virtualenv.info import fs_supports_symlink

if TYPE_CHECKING:
    from argparse import ArgumentParser
    from typing import Any

    from python_discovery import PythonInfo

    from virtualenv.app_data.base import AppData
    from virtualenv.config.cli.parser import VirtualEnvOptions

LOGGER = logging.getLogger(__name__)


class ViaGlobalRefMeta(CreatorMeta):
    def __init__(self) -> None:
        super().__init__()
        self.copy_error = None
        self.symlink_error = None
        if not fs_supports_symlink():
            self.symlink_error = "the filesystem does not supports symlink"

    @property
    def can_copy(self) -> bool:
        return not self.copy_error

    @property
    def can_symlink(self) -> bool:
        return not self.symlink_error


class ViaGlobalRefApi(Creator, ABC):
    def __init__(self, options: VirtualEnvOptions, interpreter: PythonInfo) -> None:
        super().__init__(options, interpreter)
        self.symlinks = self._should_symlink(options)
        self.enable_system_site_package = options.system_site

    if TYPE_CHECKING:

        @property
        def purelib(self) -> Path: ...

        @property
        def script_dir(self) -> Path: ...

    @staticmethod
    def _should_symlink(options: VirtualEnvOptions) -> bool:
        # Priority of where the option is set to follow the order: CLI, env var, file, hardcoded.
        # If both set at same level prefers copy over symlink.
        copies, symlinks = getattr(options, "copies", False), getattr(options, "symlinks", False)
        copy_src, sym_src = options.get_source("copies"), options.get_source("symlinks")
        for level in ["cli", "env var", "file", "default"]:
            s_opt = symlinks if sym_src == level else None
            c_opt = copies if copy_src == level else None
            if s_opt is True and c_opt is True:
                return False
            if s_opt is True:
                return True
            if c_opt is True:
                return False
        return False  # fallback to copy

    @classmethod
    def add_parser_arguments(
        cls, parser: ArgumentParser, interpreter: PythonInfo, meta: ViaGlobalRefMeta, app_data: AppData
    ) -> None:  # ty: ignore[invalid-method-override]
        super().add_parser_arguments(parser, interpreter, meta, app_data)
        parser.add_argument(
            "--system-site-packages",
            default=False,
            action="store_true",
            dest="system_site",
            help="give the virtual environment access to the system site-packages dir",
        )
        if not meta.can_symlink and not meta.can_copy:
            errors = []
            if meta.symlink_error:
                errors.append(f"symlink: {meta.symlink_error}")
            if meta.copy_error:
                errors.append(f"copy: {meta.copy_error}")
            msg = f"neither symlink or copy method supported: {', '.join(errors)}"
            raise RuntimeError(msg)
        group = parser.add_mutually_exclusive_group()
        if meta.can_symlink:
            group.add_argument(
                "--symlinks",
                default=True,
                action="store_true",
                dest="symlinks",
                help="try to use symlinks rather than copies, when symlinks are not the default for the platform",
            )
        if meta.can_copy:
            group.add_argument(
                "--copies",
                "--always-copy",
                default=not meta.can_symlink,
                action="store_true",
                dest="copies",
                help="try to use copies rather than symlinks, even when symlinks are the default for the platform",
            )

    def create(self) -> None:
        self.install_patch()

    def install_patch(self) -> None:
        # Python 3.10+ ignores the distutils install config keys this patch guards against (pip does so by default,
        # setuptools and CPython's vendored distutils drop them too), so the runtime import hook only earns its
        # startup cost on 3.9. See https://github.com/pypa/virtualenv/issues/3181.
        if self.interpreter.version_info >= (3, 10):
            return
        text = self.env_patch_text()
        if text:
            pth = self.purelib / "_virtualenv.pth"
            LOGGER.debug("create virtualenv import hook file %s", pth)
            pth.write_text("import _virtualenv", encoding="utf-8")
            dest_path = self.purelib / "_virtualenv.py"
            LOGGER.debug("create %s", dest_path)
            dest_path.write_text(text, encoding="utf-8")

    def env_patch_text(self) -> str:
        """Patch the distutils package to not be derailed by its configuration files."""
        with self.app_data.ensure_extracted(Path(__file__).parent / "_virtualenv.py") as resolved_path:
            return resolved_path.read_text(encoding="utf-8")

    def _args(self) -> list[tuple[str, Any]]:
        return [*super()._args(), ("global", self.enable_system_site_package)]

    def set_pyenv_cfg(self) -> None:
        super().set_pyenv_cfg()
        self.pyenv_cfg["include-system-site-packages"] = "true" if self.enable_system_site_package else "false"


__all__ = [
    "ViaGlobalRefApi",
    "ViaGlobalRefMeta",
]
