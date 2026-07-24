from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from virtualenv.activation.via_template import ViaTemplateActivator

if TYPE_CHECKING:
    from collections.abc import Iterator

    from virtualenv.create.creator import Creator


class BashActivator(ViaTemplateActivator):
    def templates(self) -> Iterator[str]:
        yield "activate.sh"

    def as_name(self, template: str) -> str:
        return Path(template).stem

    def replacements(self, creator: Creator, dest_folder: Path) -> dict[str, str]:
        data = super().replacements(creator, dest_folder)
        data.update({
            "__TCL_LIBRARY__": getattr(creator.interpreter, "tcl_lib", None) or "",
            "__TK_LIBRARY__": getattr(creator.interpreter, "tk_lib", None) or "",
        })
        return data


__all__ = [
    "BashActivator",
]
