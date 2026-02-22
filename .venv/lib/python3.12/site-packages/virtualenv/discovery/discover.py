from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from argparse import ArgumentParser

    from virtualenv.config.cli.parser import VirtualEnvOptions
    from virtualenv.discovery.py_info import PythonInfo


class Discover(ABC):
    """Discover and provide the requested Python interpreter."""

    @classmethod
    def add_parser_arguments(cls, parser: ArgumentParser) -> None:
        """Add CLI arguments for this discovery mechanisms.

        :param parser: the CLI parser

        """
        raise NotImplementedError

    def __init__(self, options: VirtualEnvOptions) -> None:
        """Create a new discovery mechanism.

        :param options: the parsed options as defined within :meth:`add_parser_arguments`

        """
        self._has_run = False
        self._interpreter: PythonInfo | None = None
        self._env = options.env

    @abstractmethod
    def run(self) -> PythonInfo | None:
        """Discovers an interpreter.

        :returns: the interpreter ready to use for virtual environment creation

        """
        raise NotImplementedError

    @property
    def interpreter(self) -> PythonInfo | None:
        """:returns: the interpreter as returned by :meth:`run`, cached"""
        if self._has_run is False:
            self._interpreter = self.run()
            self._has_run = True
        return self._interpreter


__all__ = [
    "Discover",
]
