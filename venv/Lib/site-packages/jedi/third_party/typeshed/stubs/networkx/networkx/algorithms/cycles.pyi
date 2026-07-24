from collections.abc import Generator
from typing import Literal

from networkx.classes.digraph import DiGraph
from networkx.classes.graph import Graph, _NBunch, _Node
from networkx.utils.backends import _dispatchable

__all__ = [
    "cycle_basis",
    "simple_cycles",
    "recursive_simple_cycles",
    "find_cycle",
    "minimum_cycle_basis",
    "chordless_cycles",
    "girth",
]

@_dispatchable
def cycle_basis(G: Graph[_Node], root: _Node | None = None) -> list[list[_Node]]: ...
@_dispatchable
def simple_cycles(G: Graph[_Node], length_bound: int | None = None) -> Generator[list[_Node]]: ...
@_dispatchable
def chordless_cycles(G: DiGraph[_Node], length_bound: int | None = None) -> Generator[list[_Node]]: ...
@_dispatchable
def recursive_simple_cycles(G: DiGraph[_Node]) -> list[list[_Node]]: ...
@_dispatchable
def find_cycle(
    G: Graph[_Node], source: _NBunch[_Node] = None, orientation: Literal["original", "reverse", "ignore"] | None = None
): ...
@_dispatchable
def minimum_cycle_basis(G: Graph[_Node], weight: str | None = None) -> list[list[_Node]]: ...
@_dispatchable
def girth(G: Graph[_Node]) -> float | int: ...  # accepts any graph type
