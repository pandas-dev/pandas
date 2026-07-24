from collections.abc import Callable, Collection, Generator
from typing import Any
from typing_extensions import TypeAlias

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = [
    "dijkstra_path",
    "dijkstra_path_length",
    "bidirectional_dijkstra",
    "single_source_dijkstra",
    "single_source_dijkstra_path",
    "single_source_dijkstra_path_length",
    "multi_source_dijkstra",
    "multi_source_dijkstra_path",
    "multi_source_dijkstra_path_length",
    "all_pairs_dijkstra",
    "all_pairs_dijkstra_path",
    "all_pairs_dijkstra_path_length",
    "dijkstra_predecessor_and_distance",
    "bellman_ford_path",
    "bellman_ford_path_length",
    "single_source_bellman_ford",
    "single_source_bellman_ford_path",
    "single_source_bellman_ford_path_length",
    "all_pairs_bellman_ford_path",
    "all_pairs_bellman_ford_path_length",
    "bellman_ford_predecessor_and_distance",
    "negative_edge_cycle",
    "find_negative_cycle",
    "goldberg_radzik",
    "johnson",
]

_WeightFunc: TypeAlias = Callable[
    [_Node, _Node, dict[str, Any]],  # Any: type of edge data cannot be known statically
    float | None,  # the weight or None to indicate a hidden edge
]

@_dispatchable
def dijkstra_path(
    G: Graph[_Node], source: _Node, target: _Node, weight: str | _WeightFunc[_Node] | None = "weight"
) -> list[_Node]: ...
@_dispatchable
def dijkstra_path_length(
    G: Graph[_Node], source: _Node, target: _Node, weight: str | _WeightFunc[_Node] | None = "weight"
) -> float: ...
@_dispatchable
def single_source_dijkstra_path(
    G: Graph[_Node], source: _Node, cutoff: float | None = None, weight: str | _WeightFunc[_Node] | None = "weight"
) -> dict[_Node, list[_Node]]: ...
@_dispatchable
def single_source_dijkstra_path_length(
    G: Graph[_Node], source: _Node, cutoff: float | None = None, weight: str | _WeightFunc[_Node] | None = "weight"
) -> dict[_Node, float]: ...
@_dispatchable
def single_source_dijkstra(
    G: Graph[_Node],
    source: _Node,
    target: _Node | None = None,
    cutoff: float | None = None,
    weight: str | _WeightFunc[_Node] | None = "weight",
) -> tuple[dict[_Node, float], dict[_Node, list[_Node]]] | tuple[float, list[_Node]]: ...  # TODO: overload on target
@_dispatchable
def multi_source_dijkstra_path(
    G: Graph[_Node], sources: Collection[_Node], cutoff: float | None = None, weight: str | _WeightFunc[_Node] | None = "weight"
) -> dict[_Node, list[_Node]]: ...
@_dispatchable
def multi_source_dijkstra_path_length(
    G: Graph[_Node], sources: Collection[_Node], cutoff: float | None = None, weight: str | _WeightFunc[_Node] | None = "weight"
) -> dict[_Node, float]: ...
@_dispatchable
def multi_source_dijkstra(
    G: Graph[_Node],
    sources: Collection[_Node],
    target: _Node | None = None,
    cutoff: float | None = None,
    weight: str | _WeightFunc[_Node] | None = "weight",
) -> tuple[dict[_Node, float], dict[_Node, list[_Node]]] | tuple[float, list[_Node]]: ...  # TODO: overload on target
@_dispatchable
def dijkstra_predecessor_and_distance(
    G: Graph[_Node], source: _Node, cutoff: float | None = None, weight: str | _WeightFunc[_Node] | None = "weight"
) -> tuple[dict[_Node, list[_Node]], dict[_Node, float]]: ...
@_dispatchable
def all_pairs_dijkstra(
    G: Graph[_Node], cutoff: float | None = None, weight: str | _WeightFunc[_Node] | None = "weight"
) -> Generator[tuple[_Node, tuple[dict[_Node, float], dict[_Node, list[_Node]]]]]: ...
@_dispatchable
def all_pairs_dijkstra_path_length(
    G: Graph[_Node], cutoff: float | None = None, weight: str | _WeightFunc[_Node] | None = "weight"
) -> Generator[tuple[_Node, dict[_Node, float]]]: ...
@_dispatchable
def all_pairs_dijkstra_path(
    G: Graph[_Node], cutoff: float | None = None, weight: str | _WeightFunc[_Node] | None = "weight"
) -> Generator[tuple[_Node, dict[_Node, list[_Node]]]]: ...
@_dispatchable
def bellman_ford_predecessor_and_distance(
    G: Graph[_Node],
    source: _Node,
    target: _Node | None = None,
    weight: str | _WeightFunc[_Node] | None = "weight",
    heuristic: bool = False,
) -> tuple[dict[_Node, list[_Node]], dict[_Node, float]]: ...
@_dispatchable
def bellman_ford_path(
    G: Graph[_Node], source: _Node, target: _Node, weight: str | _WeightFunc[_Node] | None = "weight"
) -> list[_Node]: ...
@_dispatchable
def bellman_ford_path_length(
    G: Graph[_Node], source: _Node, target: _Node, weight: str | _WeightFunc[_Node] | None = "weight"
) -> float: ...
@_dispatchable
def single_source_bellman_ford_path(
    G: Graph[_Node], source: _Node, weight: str | _WeightFunc[_Node] | None = "weight"
) -> dict[_Node, list[_Node]]: ...
@_dispatchable
def single_source_bellman_ford_path_length(
    G: Graph[_Node], source: _Node, weight: str | _WeightFunc[_Node] | None = "weight"
) -> dict[_Node, float]: ...
@_dispatchable
def single_source_bellman_ford(
    G: Graph[_Node], source: _Node, target: _Node | None = None, weight: str | _WeightFunc[_Node] | None = "weight"
) -> tuple[dict[_Node, float], dict[_Node, list[_Node]]] | tuple[float, list[_Node]]: ...  # TODO: overload on target
@_dispatchable
def all_pairs_bellman_ford_path_length(
    G: Graph[_Node], weight: str | _WeightFunc[_Node] | None = "weight"
) -> Generator[tuple[_Node, dict[_Node, float]]]: ...
@_dispatchable
def all_pairs_bellman_ford_path(
    G: Graph[_Node], weight: str | _WeightFunc[_Node] | None = "weight"
) -> Generator[tuple[_Node, dict[_Node, list[_Node]]]]: ...
@_dispatchable
def goldberg_radzik(
    G: Graph[_Node], source: _Node, weight: str | _WeightFunc[_Node] | None = "weight"
) -> tuple[dict[_Node, _Node | None], dict[_Node, float]]: ...
@_dispatchable
def negative_edge_cycle(G: Graph[_Node], weight: str | _WeightFunc[_Node] | None = "weight", heuristic: bool = True) -> bool: ...
@_dispatchable
def find_negative_cycle(G: Graph[_Node], source: _Node, weight: str | _WeightFunc[_Node] | None = "weight") -> list[_Node]: ...
@_dispatchable
def bidirectional_dijkstra(
    G: Graph[_Node], source: _Node, target: _Node, weight: str | _WeightFunc[_Node] | None = "weight"
) -> tuple[float, list[_Node]]: ...
@_dispatchable
def johnson(G: Graph[_Node], weight: str | _WeightFunc[_Node] | None = "weight") -> dict[_Node, dict[_Node, list[_Node]]]: ...
