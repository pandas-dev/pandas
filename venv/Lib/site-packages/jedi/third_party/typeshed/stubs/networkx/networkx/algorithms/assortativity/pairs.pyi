from _typeshed import Incomplete
from collections.abc import Generator, Iterable

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["node_attribute_xy", "node_degree_xy"]

@_dispatchable
def node_attribute_xy(G: Graph[_Node], attribute, nodes: Iterable[Incomplete] | None = None) -> Generator[Incomplete]: ...
@_dispatchable
def node_degree_xy(
    G: Graph[_Node], x: str = "out", y: str = "in", weight: str | None = None, nodes: Iterable[Incomplete] | None = None
) -> Generator[Incomplete]: ...
