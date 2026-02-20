from __future__ import annotations

import os
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from argparse import ArgumentParser
    from pathlib import Path

    from virtualenv.config.cli.parser import VirtualEnvOptions
    from virtualenv.create.creator import Creator
    from virtualenv.discovery.py_info import PythonInfo


class Activator(ABC):
    """Generates activate script for the virtual environment."""

    def __init__(self, options: VirtualEnvOptions) -> None:
        """Create a new activator generator.

        :param options: the parsed options as defined within :meth:`add_parser_arguments`

        """
        self.flag_prompt = os.path.basename(os.getcwd()) if options.prompt == "." else options.prompt

    @classmethod
    def supports(cls, interpreter: PythonInfo) -> bool:  # noqa: ARG003
        """Check if the activation script is supported in the given interpreter.

        :param interpreter: the interpreter we need to support

        :returns: ``True`` if supported, ``False`` otherwise

        """
        return True

    @classmethod  # noqa: B027
    def add_parser_arguments(cls, parser: ArgumentParser, interpreter: PythonInfo) -> None:
        """Add CLI arguments for this activation script.

        :param parser: the CLI parser
        :param interpreter: the interpreter this virtual environment is based of

        """

    @abstractmethod
    def generate(self, creator: Creator) -> list[Path]:
        """Generate activate script for the given creator.

        :param creator: the creator (based of :class:`virtualenv.create.creator.Creator`) we used to create this virtual
            environment

        """
        raise NotImplementedError


__all__ = [
    "Activator",
]
