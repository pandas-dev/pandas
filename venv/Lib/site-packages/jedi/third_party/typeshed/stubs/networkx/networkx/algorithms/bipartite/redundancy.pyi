from _typeshed import Incomplete
from collections.abc import Iterable

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["node_redundancy"]

@_dispatchable
def node_redundancy(G: Graph[_Node], nodes: Iterable[Incomplete] | None = None) -> dict[Incomplete, float]: ...
