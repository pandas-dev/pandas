from argparse import Namespace
from ast import AST
from collections.abc import Generator
from logging import Logger
from typing import Any

from pyflakes.checker import Checker

from ..options.manager import OptionManager

LOG: Logger
FLAKE8_PYFLAKES_CODES: dict[str, str]

class FlakesChecker(Checker):
    with_doctest: bool
    def __init__(self, tree: AST, filename: str) -> None: ...
    @classmethod
    def add_options(cls, parser: OptionManager) -> None: ...
    @classmethod
    def parse_options(cls, options: Namespace) -> None: ...
    def run(self) -> Generator[tuple[int, int, str, type[Any]]]: ...
