from __future__ import annotations

import os
from collections import OrderedDict
from typing import TYPE_CHECKING

from virtualenv.activation.via_template import ViaTemplateActivator

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

    from virtualenv.create.creator import Creator


class PythonActivator(ViaTemplateActivator):
    def templates(self) -> Iterator[str]:
        yield "activate_this.py"

    @staticmethod
    def quote(string: str) -> str:
        return repr(string)

    def replacements(self, creator: Creator, dest_folder: Path) -> dict[str, str]:
        replacements = super().replacements(creator, dest_folder)
        lib_folders = OrderedDict((os.path.relpath(str(i), str(dest_folder)), None) for i in creator.libs)
        lib_folders = os.pathsep.join(lib_folders.keys())
        replacements.update(
            {
                "__LIB_FOLDERS__": lib_folders,
                "__DECODE_PATH__": "",
            },
        )
        return replacements


__all__ = [
    "PythonActivator",
]
