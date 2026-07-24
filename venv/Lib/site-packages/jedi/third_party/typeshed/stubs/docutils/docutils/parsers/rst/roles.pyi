from collections.abc import Callable, Mapping, Sequence
from typing import Any, Final
from typing_extensions import TypeAlias, deprecated

import docutils.parsers.rst.states
from docutils import nodes
from docutils.languages import _LanguageModule
from docutils.nodes import Node, system_message
from docutils.parsers.rst.states import Inliner
from docutils.utils import Reporter

__docformat__: Final = "reStructuredText"
DEFAULT_INTERPRETED_ROLE: Final = "title-reference"

_RoleFn: TypeAlias = Callable[
    [str, str, str, int, docutils.parsers.rst.states.Inliner, Mapping[str, Any], Sequence[str]],
    tuple[Sequence[nodes.reference], Sequence[nodes.reference]],
]

def register_canonical_role(name: str, role_fn: _RoleFn) -> None: ...
def register_local_role(name: str, role_fn: _RoleFn) -> None: ...
def role(
    role_name: str, language_module: _LanguageModule, lineno: int, reporter: Reporter
) -> tuple[_RoleFn | None, list[system_message]]: ...
def set_implicit_options(role_fn: _RoleFn) -> None: ...
def register_generic_role(canonical_name: str, node_class: type[Node]) -> None: ...

class GenericRole:
    name: str
    node_class: type[Node]
    def __init__(self, role_name: str, node_class: type[Node]) -> None: ...
    def __call__(
        self,
        role: str,
        rawtext: str,
        text: str,
        lineno: int,
        inliner: Inliner,
        options: Mapping[str, Any] | None = None,
        content: Sequence[str] | None = None,
    ) -> tuple[list[Node], list[system_message]]: ...

class CustomRole:
    name: str
    base_role: _RoleFn | CustomRole
    options: Mapping[str, Any]
    content: Sequence[str]
    supplied_options: Mapping[str, Any]
    supplied_content: Sequence[str]
    def __init__(
        self,
        role_name: str,
        base_role: _RoleFn | CustomRole,
        options: Mapping[str, Any] | None = None,
        content: Sequence[str] | None = None,
    ) -> None: ...
    def __call__(
        self,
        role: str,
        rawtext: str,
        text: str,
        lineno: int,
        inliner: Inliner,
        options: Mapping[str, Any] | None = None,
        content: Sequence[str] | None = None,
    ) -> tuple[list[Node], list[system_message]]: ...

def generic_custom_role(
    role: str,
    rawtext: str,
    text: str,
    lineno: int,
    inliner: Inliner,
    options: Mapping[str, Any] | None = None,
    content: Sequence[str] | None = None,
) -> tuple[list[Node], list[system_message]]: ...
def pep_reference_role(
    role: str,
    rawtext: str,
    text: str,
    lineno: int,
    inliner: Inliner,
    options: Mapping[str, Any] | None = None,
    content: Sequence[str] | None = None,
) -> tuple[list[Node], list[system_message]]: ...
def rfc_reference_role(
    role: str,
    rawtext: str,
    text: str,
    lineno: int,
    inliner: Inliner,
    options: Mapping[str, Any] | None = None,
    content: Sequence[str] | None = None,
) -> tuple[list[Node], list[system_message]]: ...
def raw_role(
    role: str,
    rawtext: str,
    text: str,
    lineno: int,
    inliner: Inliner,
    options: Mapping[str, Any] | None = None,
    content: Sequence[str] | None = None,
) -> tuple[list[Node], list[system_message]]: ...
def code_role(
    role: str,
    rawtext: str,
    text: str,
    lineno: int,
    inliner: Inliner,
    options: Mapping[str, Any] | None = None,
    content: Sequence[str] | None = None,
) -> tuple[list[Node], list[system_message]]: ...
def math_role(
    role: str,
    rawtext: str,
    text: str,
    lineno: int,
    inliner: Inliner,
    options: Mapping[str, Any] | None = None,
    content: Sequence[str] | None = None,
) -> tuple[list[Node], list[system_message]]: ...
def unimplemented_role(
    role: str,
    rawtext: str,
    text: str,
    lineno: int,
    inliner: Inliner,
    options: Mapping[str, Any] | None = None,
    content: Sequence[str] | None = None,
) -> tuple[list[Node], list[system_message]]: ...
@deprecated("Deprecated and will be removed in Docutils 2.0, Use `roles.normalize_options()` instead.")
def set_classes(options: dict[str, str]) -> None: ...
@deprecated("Deprecated and will be removed in Docutils 2.0, Use `roles.normalize_options()` instead.")
def normalized_role_options(options: Mapping[str, Any] | None) -> dict[str, Any]: ...
def normalize_options(options: Mapping[str, Any] | None) -> dict[str, Any]: ...
