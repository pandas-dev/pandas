from argparse import ArgumentParser
from collections.abc import Iterable, Mapping, Sequence
from typing import Any, overload
from typing_extensions import TypeAlias

from xmldiff.actions import (
    DeleteAttrib,
    DeleteNamespace,
    DeleteNode,
    InsertAttrib,
    InsertComment,
    InsertNamespace,
    InsertNode,
    MoveNode,
    RenameAttrib,
    RenameNode,
    UpdateAttrib,
    UpdateTextAfter,
    UpdateTextIn,
)
from xmldiff.formatting import BaseFormatter

__version__: str
FORMATTERS: Mapping[str, BaseFormatter]

_ACTIONS: TypeAlias = (
    DeleteNode
    | InsertNode
    | RenameNode
    | MoveNode
    | UpdateTextIn
    | UpdateTextAfter
    | UpdateAttrib
    | DeleteAttrib
    | InsertAttrib
    | RenameAttrib
    | InsertComment
    | InsertNamespace
    | DeleteNamespace
)
_ET: TypeAlias = Any  # lxml.etree._ElementTree

@overload
def diff_trees(left: _ET, right: _ET, *, diff_options: dict[str, Any] | None = None, formatter: BaseFormatter) -> str: ...
@overload
def diff_trees(
    left: _ET, right: _ET, diff_options: dict[str, Any] | None = None, formatter: None = None
) -> Iterable[_ACTIONS]: ...
@overload
def diff_texts(
    left: str | bytes, right: str | bytes, *, diff_options: dict[str, Any] | None = None, formatter: BaseFormatter
) -> str: ...
@overload
def diff_texts(
    left: str | bytes, right: str | bytes, diff_options: dict[str, Any] | None = None, formatter: None = None
) -> Iterable[_ACTIONS]: ...
@overload
def diff_files(left: str, right: str, *, diff_options: dict[str, Any] | None = None, formatter: BaseFormatter) -> str: ...
@overload
def diff_files(
    left: str, right: str, diff_options: dict[str, Any] | None = None, formatter: None = None
) -> Iterable[_ACTIONS]: ...
def validate_F(arg: float | str) -> float: ...
def make_diff_parser() -> ArgumentParser: ...
def diff_command(args: Sequence[str] | None = None) -> int | None: ...
def patch_tree(actions, tree): ...
def patch_text(actions, tree): ...
def patch_file(actions, tree, diff_encoding=None): ...
def make_patch_parser(): ...
def patch_command(args=None) -> None: ...
