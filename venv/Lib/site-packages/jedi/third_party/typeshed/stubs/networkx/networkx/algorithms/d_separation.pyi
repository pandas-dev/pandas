from _typeshed import Incomplete

from networkx.classes.digraph import DiGraph
from networkx.classes.graph import _Node
from networkx.utils.backends import _dispatchable

__all__ = ["is_d_separator", "is_minimal_d_separator", "find_minimal_d_separator"]

@_dispatchable
def is_d_separator(G: DiGraph[_Node], x: _Node | set[_Node], y: _Node | set[_Node], z: _Node | set[_Node]) -> bool: ...
@_dispatchable
def find_minimal_d_separator(G: DiGraph[_Node], x, y, *, included=None, restricted=None) -> set[Incomplete] | None: ...
@_dispatchable
def is_minimal_d_separator(
    G: DiGraph[_Node],
    x: _Node | set[_Node],
    y: _Node | set[_Node],
    z: _Node | set[_Node],
    *,
    included: _Node | set[_Node] | None = None,
    restricted: _Node | set[_Node] | None = None,
) -> bool: ...
