from collections.abc import Container, Iterable

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["is_partition"]

@_dispatchable
def is_partition(G: Graph[_Node], communities: Iterable[Container[_Node]]) -> bool: ...
