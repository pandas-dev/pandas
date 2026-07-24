from _typeshed import Incomplete
from collections.abc import Iterable

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["reciprocity", "overall_reciprocity"]

@_dispatchable
def reciprocity(G: Graph[_Node], nodes: Iterable[_Node] | None = None) -> float | dict[Incomplete, float | None]: ...
@_dispatchable
def overall_reciprocity(G: Graph[_Node]): ...
