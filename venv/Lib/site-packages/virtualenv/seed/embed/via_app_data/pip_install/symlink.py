from __future__ import annotations

import os
from stat import S_IREAD, S_IRGRP, S_IROTH
from subprocess import PIPE, Popen
from typing import TYPE_CHECKING

from virtualenv.util.path import safe_delete, set_tree

from .base import PipInstall

if TYPE_CHECKING:
    from pathlib import Path


class SymlinkPipInstall(PipInstall):
    def _sync(self, src: Path, dst: Path) -> None:
        os.symlink(str(src), str(dst))

    def _generate_new_files(self) -> set[Path]:
        # create the pyc files, as the build image will be R/O
        cmd = [str(self._creator.exe), "-m", "compileall", str(self._image_dir)]
        process = Popen(cmd, stdout=PIPE, stderr=PIPE)
        process.communicate()
        # the root pyc is shared, so we'll not symlink that - but still add the pyc files to the RECORD for close
        root_py_cache = self._image_dir / "__pycache__"
        new_files = set()
        if root_py_cache.exists():
            new_files.update(root_py_cache.iterdir())
            new_files.add(root_py_cache)
            safe_delete(root_py_cache)
        core_new_files = super()._generate_new_files()
        # remove files that are within the image folder deeper than one level (as these will be not linked directly)
        for file in core_new_files:
            try:
                rel = file.relative_to(self._image_dir)
                if len(rel.parts) > 1:
                    continue
            except ValueError:
                pass
            new_files.add(file)
        return new_files

    def _fix_records(self, extra_record_data: set[Path]) -> None:
        extra_record_data.update(i for i in self._image_dir.iterdir())
        extra_record_data_str = self._records_text(sorted(extra_record_data, key=str))  # ty: ignore[invalid-argument-type]
        (self._dist_info / "RECORD").write_text(extra_record_data_str, encoding="utf-8")  # ty: ignore[unsupported-operator]

    def build_image(self) -> None:
        super().build_image()
        # protect the image by making it read only
        set_tree(self._image_dir, S_IREAD | S_IRGRP | S_IROTH)

    def clear(self) -> None:
        if self._image_dir.exists():
            safe_delete(self._image_dir)
        super().clear()


__all__ = [
    "SymlinkPipInstall",
]
