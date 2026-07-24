from collections.abc import Callable

from networkx.algorithms.shortest_paths.weighted import _WeightFunc
from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["astar_path", "astar_path_length"]

@_dispatchable
def astar_path(
    G: Graph[_Node],
    source: _Node,
    target: _Node,
    heuristic: Callable[[_Node, _Node], float] | None = None,
    weight: str | _WeightFunc[_Node] | None = "weight",
    *,
    cutoff: float | None = None,
) -> list[_Node]: ...
@_dispatchable
def astar_path_length(
    G: Graph[_Node],
    source: _Node,
    target: _Node,
    heuristic: Callable[[_Node, _Node], float] | None = None,
    weight: str | _WeightFunc[_Node] | None = "weight",
    *,
    cutoff: float | None = None,
) -> float: ...
