from __future__ import annotations

import json
import logging
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from types import TracebackType

    from virtualenv.activation.activator import Activator
    from virtualenv.app_data.base import AppData
    from virtualenv.create.creator import Creator
    from virtualenv.discovery.py_info import PythonInfo
    from virtualenv.seed.seeder import Seeder

if sys.version_info >= (3, 11):
    from typing import Self
else:
    from typing_extensions import Self

LOGGER = logging.getLogger(__name__)


class Session:
    """Represents a virtual environment creation session."""

    def __init__(  # noqa: PLR0913
        self,
        verbosity: int,
        app_data: AppData,
        interpreter: PythonInfo,
        creator: Creator,
        seeder: Seeder,
        activators: list[Activator],
    ) -> None:
        self._verbosity = verbosity
        self._app_data = app_data
        self._interpreter = interpreter
        self._creator = creator
        self._seeder = seeder
        self._activators = activators

    @property
    def verbosity(self) -> int:
        """The verbosity of the run."""
        return self._verbosity

    @property
    def interpreter(self) -> PythonInfo:
        """Create a virtual environment based on this reference interpreter."""
        return self._interpreter

    @property
    def creator(self) -> Creator:
        """The creator used to build the virtual environment (must be compatible with the interpreter)."""
        return self._creator

    @property
    def seeder(self) -> Seeder:
        """The mechanism used to provide the seed packages (pip, setuptools, wheel)."""
        return self._seeder

    @property
    def activators(self) -> list[Activator]:
        """Activators used to generate activations scripts."""
        return self._activators

    def run(self) -> None:
        self._create()
        self._seed()
        self._activate()
        self.creator.pyenv_cfg.write()

    def _create(self) -> None:
        LOGGER.info("create virtual environment via %s", self.creator)
        self.creator.run()
        LOGGER.debug(_DEBUG_MARKER)
        LOGGER.debug("%s", _Debug(self.creator))

    def _seed(self) -> None:
        if self.seeder is not None and self.seeder.enabled:
            LOGGER.info("add seed packages via %s", self.seeder)
            self.seeder.run(self.creator)

    def _activate(self) -> None:
        if self.activators:
            active = ", ".join(type(i).__name__.replace("Activator", "") for i in self.activators)
            LOGGER.info("add activators for %s", active)
            for activator in self.activators:
                activator.generate(self.creator)

    def __enter__(self) -> Self:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: TracebackType | None,
    ) -> None:
        self._app_data.close()


_DEBUG_MARKER = "=" * 30 + " target debug " + "=" * 30


class _Debug:
    """lazily populate debug."""

    def __init__(self, creator) -> None:
        self.creator = creator

    def __repr__(self) -> str:
        return json.dumps(self.creator.debug, indent=2)


__all__ = [
    "Session",
]
