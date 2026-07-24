import ast
from argparse import Namespace
from collections.abc import Container, Generator
from typing import Any

rst_prefix: str
rst_fail_load: int
rst_fail_lint: int
code_mapping_info: dict[str, int]
code_mapping_warning: dict[str, int]
code_mapping_error: dict[str, int]
code_mapping_severe: dict[str, int]
code_mappings_by_level: dict[int, dict[str, int]]

def code_mapping(
    level: int,
    msg: str,
    extra_directives: Container[str],
    extra_roles: Container[str],
    extra_substitutions: Container[str],
    default: int = ...,
) -> int: ...

class reStructuredTextChecker:
    name: str
    version: str
    tree: ast.AST
    filename: str
    def __init__(self, tree: ast.AST, filename: str = ...) -> None: ...
    @classmethod
    def add_options(cls, parser: Any) -> None: ...
    @classmethod
    def parse_options(cls, options: Namespace) -> None: ...
    def run(self) -> Generator[tuple[int, int, str, type[reStructuredTextChecker]]]: ...
