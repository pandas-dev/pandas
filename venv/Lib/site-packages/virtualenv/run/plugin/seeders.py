from __future__ import annotations

from typing import TYPE_CHECKING

from .base import ComponentBuilder

if TYPE_CHECKING:
    from collections.abc import Sequence

    from python_discovery import PythonInfo

    from virtualenv.config.cli.parser import VirtualEnvConfigParser, VirtualEnvOptions
    from virtualenv.seed.seeder import Seeder


class SeederSelector(ComponentBuilder):
    def __init__(self, interpreter: PythonInfo, parser: VirtualEnvConfigParser) -> None:
        possible = self.options("virtualenv.seed")
        super().__init__(interpreter, parser, "seeder", possible)

    def add_selector_arg_parse(self, name: str, choices: Sequence[str]) -> None:
        self.parser.add_argument(
            f"--{name}",
            choices=choices,
            default=self._get_default(),
            required=False,
            help="seed packages install method",
        )
        self.parser.add_argument(
            "--no-seed",
            "--without-pip",
            help="do not install seed packages",
            action="store_true",
            dest="no_seed",
        )

    @staticmethod
    def _get_default() -> str:
        return "app-data"

    def handle_selected_arg_parse(self, options: VirtualEnvOptions) -> str:
        return super().handle_selected_arg_parse(options)

    def create(self, options: VirtualEnvOptions) -> Seeder:
        assert self._impl_class is not None  # ruff:ignore[assert]  # Set by handle_selected_arg_parse
        seeder = self._impl_class(options)
        if seeder.enabled and (reason := seeder.cannot_seed(self.interpreter)) is not None:
            raise RuntimeError(reason)
        return seeder


__all__ = [
    "SeederSelector",
]
