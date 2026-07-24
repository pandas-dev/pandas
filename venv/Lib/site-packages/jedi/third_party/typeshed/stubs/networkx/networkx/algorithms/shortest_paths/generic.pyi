from collections.abc import Generator
from typing import overload

from networkx.algorithms.shortest_paths.weighted import _WeightFunc
from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = [
    "shortest_path",
    "all_shortest_paths",
    "single_source_all_shortest_paths",
    "all_pairs_all_shortest_paths",
    "shortest_path_length",
    "average_shortest_path_length",
    "has_path",
]

@_dispatchable
def has_path(G: Graph[_Node], source: _Node, target: _Node) -> bool: ...
@overload  # both source and target are specified => (s -> t)
def shortest_path(
    G: Graph[_Node],
    source: _Node,
    target: _Node,
    weight: str | _WeightFunc[_Node] | None = None,
    method: str | None = "dijkstra",
    *,
    backend: str | None = None,
    **backend_kwargs,
) -> list[_Node]: ...
@overload  # only source is specified => {t1: (s -> t), t2: (s -> t), ...}
def shortest_path(
    G: Graph[_Node],
    source: _Node,
    target: None = None,
    weight: str | _WeightFunc[_Node] | None = None,
    method: str | None = "dijkstra",
    *,
    backend: str | None = None,
    **backend_kwargs,
) -> dict[_Node, list[_Node]]: ...
@overload  # only target is specified (positional) => {s1: (s1 -> t), s2: (s2 -> t), ...}
def shortest_path(
    G: Graph[_Node],
    source: None,
    target: _Node,
    weight: str | _WeightFunc[_Node] | None = None,
    method: str | None = "dijkstra",
    *,
    backend: str | None = None,
    **backend_kwargs,
) -> dict[_Node, list[_Node]]: ...
@overload  # only target is specified (keyword) => {s1: (s1 -> t), s2: (s2 -> t), ...}
def shortest_path(
    G: Graph[_Node],
    source: None = None,
    *,
    target: _Node,
    weight: str | _WeightFunc[_Node] | None = None,
    method: str | None = "dijkstra",
    backend: str | None = None,
    **backend_kwargs,
) -> dict[_Node, list[_Node]]: ...
@overload
def shortest_path(  # source and target are not specified => generator of (t, {s1: (s1 -> t), s2: (s2 -> t), ...})
    G: Graph[_Node],
    source: None = None,
    target: None = None,
    weight: str | _WeightFunc[_Node] | None = None,
    method: str | None = "dijkstra",
    *,
    backend: str | None = None,
    **backend_kwargs,
) -> Generator[tuple[_Node, dict[str, list[_Node]]]]: ...
@overload  # both source and target are specified => len(s -> t)
def shortest_path_length(
    G: Graph[_Node],
    source: _Node,
    target: _Node,
    weight: str | _WeightFunc[_Node] | None = None,
    method: str | None = "dijkstra",
    *,
    backend: str | None = None,
    **backend_kwargs,
) -> float: ...
@overload  # only source is specified => {t1: len(s -> t1), t2: len(s -> t2), ...}
def shortest_path_length(
    G: Graph[_Node],
    source: _Node,
    target: None = None,
    weight: str | _WeightFunc[_Node] | None = None,
    method: str | None = "dijkstra",
    *,
    backend: str | None = None,
    **backend_kwargs,
) -> dict[_Node, float]: ...
@overload  # only target is specified (positional) => {s1: len(s1 -> t), s2: len(s2 -> t), ...}
def shortest_path_length(
    G: Graph[_Node],
    source: None,
    target: _Node,
    weight: str | _WeightFunc[_Node] | None = None,
    method: str | None = "dijkstra",
    *,
    backend: str | None = None,
    **backend_kwargs,
) -> dict[_Node, float]: ...
@overload  # only target is specified (keyword) => {s1: len(s1 -> t), s2: len(s2 -> t), ...}
def shortest_path_length(
    G: Graph[_Node],
    source: None = None,
    *,
    target: _Node,
    weight: str | _WeightFunc[_Node] | None = None,
    method: str | None = "dijkstra",
    backend: str | None = None,
    **backend_kwargs,
) -> dict[_Node, float]: ...
@overload
def shortest_path_length(  # source and target are not specified => generator of (t, {s1: len(s1 -> t), s2: len(s2 -> t), ...})
    G: Graph[_Node],
    source: None = None,
    target: None = None,
    weight: str | _WeightFunc[_Node] | None = None,
    method: str | None = "dijkstra",
    *,
    backend: str | None = None,
    **backend_kwargs,
) -> Generator[tuple[_Node, dict[_Node, float]]]: ...
@_dispatchable
def average_shortest_path_length(
    G: Graph[_Node], weight: str | _WeightFunc[_Node] | None = None, method: str | None = None
) -> float: ...
@_dispatchable
def all_shortest_paths(
    G: Graph[_Node], source: _Node, target: _Node, weight: str | _WeightFunc[_Node] | None = None, method: str | None = "dijkstra"
) -> Generator[list[_Node]]: ...
@_dispatchable
def single_source_all_shortest_paths(
    G: Graph[_Node], source: _Node, weight: str | _WeightFunc[_Node] | None = None, method: str | None = "dijkstra"
) -> Generator[tuple[_Node, list[list[_Node]]]]: ...
@_dispatchable
def all_pairs_all_shortest_paths(
    G: Graph[_Node], weight: str | _WeightFunc[_Node] | None = None, method: str | None = "dijkstra"
) -> Generator[tuple[_Node, dict[_Node, list[list[_Node]]]]]: ...
