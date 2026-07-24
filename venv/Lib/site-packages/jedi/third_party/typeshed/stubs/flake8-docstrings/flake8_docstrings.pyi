import argparse
import ast
from collections.abc import Generator, Iterable
from typing import Any, ClassVar, Final, Literal
from typing_extensions import Self

from flake8.options.manager import OptionManager

__version__: Final[str]
__all__ = ("pep257Checker",)

class pep257Checker:
    name: ClassVar[str]
    version: ClassVar[str]
    tree: ast.AST
    filename: str
    checker: Any  # actual type: pep257.ConventionChecker
    source: str
    def __init__(self, tree: ast.AST, filename: str, lines: Iterable[str]) -> None: ...
    @classmethod
    def add_options(cls, parser: OptionManager) -> None: ...
    @classmethod
    def parse_options(cls, options: argparse.Namespace) -> None: ...
    def run(self) -> Generator[tuple[int, Literal[0], str, type[Self]]]: ...
