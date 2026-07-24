from __future__ import annotations

from abc import ABC
from typing import TYPE_CHECKING

from virtualenv.create.creator import Creator
from virtualenv.create.describe import Describe

if TYPE_CHECKING:
    from python_discovery import PythonInfo

    from virtualenv.config.cli.parser import VirtualEnvOptions


class VirtualenvBuiltin(Creator, Describe, ABC):
    """A creator that does operations itself without delegation, if we can create it we can also describe it."""

    def __init__(self, options: VirtualEnvOptions, interpreter: PythonInfo) -> None:
        Creator.__init__(self, options, interpreter)
        Describe.__init__(self, self.dest, interpreter)


__all__ = [
    "VirtualenvBuiltin",
]
