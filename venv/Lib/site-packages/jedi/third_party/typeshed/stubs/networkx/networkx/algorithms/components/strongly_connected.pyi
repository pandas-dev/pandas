from collections.abc import Generator, Iterable

from networkx.classes.digraph import DiGraph
from networkx.classes.graph import _Node
from networkx.utils.backends import _dispatchable

__all__ = [
    "number_strongly_connected_components",
    "strongly_connected_components",
    "is_strongly_connected",
    "kosaraju_strongly_connected_components",
    "condensation",
]

@_dispatchable
def strongly_connected_components(G: DiGraph[_Node]) -> Generator[set[_Node]]: ...
@_dispatchable
def kosaraju_strongly_connected_components(G: DiGraph[_Node], source: _Node | None = None) -> Generator[set[_Node]]: ...
@_dispatchable
def number_strongly_connected_components(G: DiGraph[_Node]) -> int: ...
@_dispatchable
def is_strongly_connected(G: DiGraph[_Node]) -> bool: ...
@_dispatchable
def condensation(G: DiGraph[_Node], scc: Iterable[Iterable[_Node]] | None = None) -> DiGraph[int]: ...
