from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from argparse import ArgumentParser

    from python_discovery import PythonInfo

    from virtualenv.app_data.base import AppData
    from virtualenv.config.cli.parser import VirtualEnvOptions
    from virtualenv.create.creator import Creator


class Seeder(ABC):
    """A seeder will install some seed packages into a virtual environment."""

    def __init__(self, options: VirtualEnvOptions, enabled: bool) -> None:
        """Create.

        :param options: the parsed options as defined within :meth:`add_parser_arguments`
        :param enabled: a flag weather the seeder is enabled or not

        """
        self.enabled = enabled
        self.env = options.env

    @classmethod
    def cannot_seed(cls, interpreter: PythonInfo) -> str | None:  # ruff:ignore[unused-class-method-argument]
        """Explain why this seeder cannot install seed packages for the given interpreter.

        :param interpreter: the interpreter the environment is based on

        :returns: ``None`` when the seeder supports the interpreter, otherwise a message describing why it cannot;
            selection rejects a seeder that returns a message and surfaces it to the user

        """
        return None

    @classmethod
    def add_parser_arguments(cls, parser: ArgumentParser, interpreter: PythonInfo, app_data: AppData) -> None:
        """Add CLI arguments for this seed mechanisms.

        :param parser: the CLI parser
        :param app_data: the CLI parser
        :param interpreter: the interpreter this virtual environment is based of

        """
        raise NotImplementedError

    @abstractmethod
    def run(self, creator: Creator) -> None:
        """Perform the seed operation.

        :param creator: the creator (based of :class:`virtualenv.create.creator.Creator`) we used to create this virtual
            environment

        """
        raise NotImplementedError


__all__ = [
    "Seeder",
]
