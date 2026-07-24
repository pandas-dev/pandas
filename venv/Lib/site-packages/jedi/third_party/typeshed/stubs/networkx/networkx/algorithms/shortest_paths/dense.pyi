from _typeshed import Incomplete, SupportsGetItem
from collections import defaultdict
from collections.abc import Collection

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["floyd_warshall", "floyd_warshall_predecessor_and_distance", "reconstruct_path", "floyd_warshall_numpy"]

@_dispatchable
def floyd_warshall_numpy(G: Graph[_Node], nodelist: Collection[_Node] | None = None, weight: str | None = "weight"): ...
@_dispatchable
def floyd_warshall_predecessor_and_distance(
    G: Graph[_Node], weight: str | None = "weight"
) -> tuple[dict[Incomplete, dict[Incomplete, Incomplete]], dict[Incomplete, dict[Incomplete, float]]]: ...
@_dispatchable
def reconstruct_path(source: _Node, target: _Node, predecessors: SupportsGetItem[Incomplete, Incomplete]) -> list[Incomplete]: ...
@_dispatchable
def floyd_warshall(G: Graph[_Node], weight: str | None = "weight") -> dict[Incomplete, defaultdict[Incomplete, float]]: ...
