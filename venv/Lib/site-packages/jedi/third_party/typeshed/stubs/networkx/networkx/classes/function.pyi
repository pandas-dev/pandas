from _typeshed import Incomplete, SupportsItems, SupportsKeysAndGetItem, Unused
from collections.abc import Callable, Generator, Hashable, Iterable, Iterator
from typing import Literal, TypeVar, overload

from networkx import _dispatchable
from networkx.algorithms.planarity import PlanarEmbedding
from networkx.classes.digraph import DiGraph
from networkx.classes.graph import Graph, _NBunch, _Node
from networkx.classes.multigraph import MultiGraph

__all__ = [
    "nodes",
    "edges",
    "degree",
    "degree_histogram",
    "neighbors",
    "number_of_nodes",
    "number_of_edges",
    "density",
    "is_directed",
    "freeze",
    "is_frozen",
    "subgraph",
    "induced_subgraph",
    "edge_subgraph",
    "restricted_view",
    "to_directed",
    "to_undirected",
    "add_star",
    "add_path",
    "add_cycle",
    "create_empty_copy",
    "set_node_attributes",
    "get_node_attributes",
    "remove_node_attributes",
    "set_edge_attributes",
    "get_edge_attributes",
    "remove_edge_attributes",
    "all_neighbors",
    "non_neighbors",
    "non_edges",
    "common_neighbors",
    "is_weighted",
    "is_negatively_weighted",
    "is_empty",
    "selfloop_edges",
    "nodes_with_selfloops",
    "number_of_selfloops",
    "path_weight",
    "is_path",
    "describe",
]

_U = TypeVar("_U")

def nodes(G: Graph[_Node]): ...
def edges(G: Graph[_Node], nbunch=None): ...
def degree(G: Graph[_Node], nbunch=None, weight=None): ...
def neighbors(G: Graph[_Node], n): ...
def number_of_nodes(G: Graph[_Node]): ...
def number_of_edges(G: Graph[_Node]): ...
def density(G: Graph[_Node]): ...
def degree_histogram(G: Graph[_Node]) -> list[int]: ...
@overload
def is_directed(G: PlanarEmbedding[Hashable]) -> Literal[False]: ...  # type: ignore[misc] # Incompatible return types
@overload
def is_directed(G: DiGraph[Hashable]) -> Literal[True]: ...  # type: ignore[misc] # Incompatible return types
@overload
def is_directed(G: Graph[Hashable]) -> Literal[False]: ...
def freeze(G: Graph[_Node]): ...
def is_frozen(G: Graph[Incomplete]) -> bool: ...
def add_star(G_to_add_to, nodes_for_star, **attr) -> None: ...
def add_path(G_to_add_to, nodes_for_path, **attr) -> None: ...
def add_cycle(G_to_add_to, nodes_for_cycle, **attr) -> None: ...
def subgraph(G: Graph[_Node], nbunch): ...
def induced_subgraph(G: Graph[_Node], nbunch: _NBunch[_Node]) -> Graph[_Node]: ...
def edge_subgraph(G: Graph[_Node], edges): ...
def restricted_view(G: Graph[_Node], nodes, edges): ...
def to_directed(graph): ...
def to_undirected(graph): ...
def create_empty_copy(G: Graph[_Node], with_data: bool = True): ...

# incomplete: Can "Any scalar value" be enforced?
@overload
def set_node_attributes(
    G: Graph[Hashable],
    values: SupportsItems[_Node, Unused],
    name: str,
    *,
    backend=None,  # @_dispatchable adds these arguments, but we can't use this decorator with @overload
    **backend_kwargs,
) -> None: ...
@overload
def set_node_attributes(
    G: Graph[_Node],
    values: SupportsItems[_Node, SupportsKeysAndGetItem[Incomplete, Incomplete] | Iterable[tuple[Incomplete, Incomplete]]],
    name: None = None,
    *,
    backend=None,
    **backend_kwargs,
) -> None: ...
@_dispatchable
def get_node_attributes(G: Graph[_Node], name: str, default=None) -> dict[_Node, Incomplete]: ...
@_dispatchable
def remove_node_attributes(G: Graph[_Node], *attr_names, nbunch=None) -> None: ...
@overload
def set_edge_attributes(
    G: Graph[_Node],
    values: SupportsItems[tuple[_Node, _Node], Incomplete],
    name: str,
    *,
    backend: str | None = None,  # @_dispatchable adds these arguments, but we can't use this decorator with @overload
    **backend_kwargs,
) -> None: ...
@overload
def set_edge_attributes(
    G: MultiGraph[_Node],
    values: dict[tuple[_Node, _Node, Incomplete], Incomplete],
    name: str,
    *,
    backend: str | None = None,
    **backend_kwargs,
) -> None: ...
@overload
def set_edge_attributes(
    G: Graph[Hashable], values, name: None = None, *, backend: str | None = None, **backend_kwargs
) -> None: ...
@_dispatchable
def get_edge_attributes(G: Graph[_Node], name: str, default=None) -> dict[tuple[_Node, _Node], Incomplete]: ...
@_dispatchable
def remove_edge_attributes(G: Graph[_Node], *attr_names, ebunch=None) -> None: ...
def all_neighbors(graph: Graph[_Node], node: _Node) -> Iterator[_Node]: ...
def non_neighbors(graph: Graph[_Node], node: _Node) -> Generator[_Node]: ...
def non_edges(graph: Graph[_Node]) -> Generator[tuple[_Node, _Node]]: ...
def common_neighbors(G: Graph[_Node], u: _Node, v: _Node) -> Generator[_Node]: ...
@_dispatchable
def is_weighted(G: Graph[_Node], edge: tuple[_Node, _Node] | None = None, weight: str = "weight") -> bool: ...
@_dispatchable
def is_negatively_weighted(G: Graph[_Node], edge: tuple[_Node, _Node] | None = None, weight: str = "weight") -> bool: ...
@_dispatchable
def is_empty(G: Graph[Hashable]) -> bool: ...
def nodes_with_selfloops(G: Graph[_Node]) -> Generator[_Node]: ...
@overload
def selfloop_edges(
    G: Graph[_Node], data: Literal[False] = False, keys: Literal[False] = False, default=None
) -> Generator[tuple[_Node, _Node]]: ...
@overload
def selfloop_edges(
    G: Graph[_Node], data: Literal[True], keys: Literal[False] = False, default=None
) -> Generator[tuple[_Node, _Node, dict[str, Incomplete]]]: ...
@overload
def selfloop_edges(
    G: Graph[_Node], data: str, keys: Literal[False] = False, default: _U | None = None
) -> Generator[tuple[_Node, _Node, _U]]: ...
@overload
def selfloop_edges(
    G: Graph[_Node], data: Literal[False], keys: Literal[True], default=None
) -> Generator[tuple[_Node, _Node, int]]: ...
@overload
def selfloop_edges(
    G: Graph[_Node], data: Literal[False] = False, *, keys: Literal[True], default=None
) -> Generator[tuple[_Node, _Node, int]]: ...
@overload
def selfloop_edges(
    G: Graph[_Node], data: Literal[True], keys: Literal[True], default=None
) -> Generator[tuple[_Node, _Node, int, dict[str, Incomplete]]]: ...
@overload
def selfloop_edges(
    G: Graph[_Node], data: str, keys: Literal[True], default: _U | None = None
) -> Generator[tuple[_Node, _Node, int, _U]]: ...
@_dispatchable
def number_of_selfloops(G: Graph[Hashable]) -> int: ...
def is_path(G: Graph[_Node], path: Iterable[Incomplete]) -> bool: ...
def path_weight(G: Graph[_Node], path, weight) -> int: ...
def describe(G: Graph[_Node], describe_hook: Callable[[Graph[_Node]], dict[str, Incomplete]] | None = None) -> None: ...
