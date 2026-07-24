from __future__ import annotations

import os
from typing import TYPE_CHECKING

from virtualenv.activation.via_template import ViaTemplateActivator

if TYPE_CHECKING:
    from collections.abc import Iterator

    from python_discovery import PythonInfo

    from virtualenv.create.creator import Creator


class BatchActivator(ViaTemplateActivator):
    @classmethod
    def supports(cls, interpreter: PythonInfo) -> bool:
        return interpreter.os == "nt"

    def templates(self) -> Iterator[str]:
        yield "activate.bat"
        yield "deactivate.bat"
        yield "pydoc.bat"

    @staticmethod
    def quote(string: str) -> str:
        return string

    def instantiate_template(self, replacements: dict[str, str], template: str, creator: Creator) -> str:
        # ensure the text has all newlines as \r\n - required by batch
        base = super().instantiate_template(replacements, template, creator)
        return base.replace(os.linesep, "\n").replace("\n", os.linesep)


__all__ = [
    "BatchActivator",
]
