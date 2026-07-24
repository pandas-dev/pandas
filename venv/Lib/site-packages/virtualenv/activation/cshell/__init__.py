from __future__ import annotations

from typing import TYPE_CHECKING

from virtualenv.activation.via_template import ViaTemplateActivator

if TYPE_CHECKING:
    from collections.abc import Iterator

    from python_discovery import PythonInfo


class CShellActivator(ViaTemplateActivator):
    @classmethod
    def supports(cls, interpreter: PythonInfo) -> bool:
        return interpreter.os != "nt"

    def templates(self) -> Iterator[str]:
        yield "activate.csh"


__all__ = [
    "CShellActivator",
]
