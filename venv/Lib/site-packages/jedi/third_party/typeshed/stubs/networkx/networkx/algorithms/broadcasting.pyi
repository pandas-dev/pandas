import networkx as nx
from networkx.classes.graph import Graph, _Node

__all__ = ["tree_broadcast_center", "tree_broadcast_time"]

@nx._dispatchable
def tree_broadcast_center(G: Graph[_Node]) -> tuple[int, set[_Node]]: ...
@nx._dispatchable
def tree_broadcast_time(G: Graph[_Node], node: int | None = None) -> int: ...
