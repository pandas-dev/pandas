from _typeshed import Incomplete
from collections.abc import Callable, Iterator, Sequence
from re import Match
from typing import Any, Generic, TypeVar

from parsimonious.exceptions import VisitationError as VisitationError
from parsimonious.expressions import Expression
from parsimonious.grammar import Grammar

class Node:
    __slots__ = ["expr", "full_text", "start", "end", "children"]
    expr: Expression
    full_text: str
    start: int
    end: int
    children: Sequence[Node]
    def __init__(
        self, expr: Expression, full_text: str, start: int, end: int, children: Sequence[Node] | None = None
    ) -> None: ...
    @property
    def expr_name(self) -> str: ...
    def __iter__(self) -> Iterator[Node]: ...
    @property
    def text(self) -> str: ...
    def prettily(self, error: Node | None = None) -> str: ...
    def __repr__(self, top_level: bool = True) -> str: ...

class RegexNode(Node):
    __slots__ = ["match"]
    match: Match[str]

class RuleDecoratorMeta(type): ...

_VisitResultT = TypeVar("_VisitResultT")
_ChildT = TypeVar("_ChildT")

class NodeVisitor(Generic[_VisitResultT], metaclass=RuleDecoratorMeta):
    grammar: Grammar | Incomplete
    unwrapped_exceptions: tuple[type[BaseException], ...]
    def visit(self, node: Node) -> _VisitResultT: ...
    def generic_visit(self, node: Node, visited_children: Sequence[Any]): ...
    def parse(self, text: str, pos: int = 0) -> _VisitResultT: ...
    def match(self, text: str, pos: int = 0) -> _VisitResultT: ...
    def lift_child(self, node: Node, children: Sequence[_ChildT]) -> _ChildT: ...

_CallableT = TypeVar("_CallableT", bound=Callable[..., Any])

def rule(rule_string: str) -> Callable[[_CallableT], _CallableT]: ...
