from _typeshed import Incomplete
from collections.abc import Callable, Generator, Iterable

from networkx.classes.digraph import DiGraph
from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = [
    "descendants",
    "ancestors",
    "topological_sort",
    "lexicographical_topological_sort",
    "all_topological_sorts",
    "topological_generations",
    "is_directed_acyclic_graph",
    "is_aperiodic",
    "transitive_closure",
    "transitive_closure_dag",
    "transitive_reduction",
    "antichains",
    "dag_longest_path",
    "dag_longest_path_length",
    "dag_to_branching",
]

@_dispatchable
def descendants(G: Graph[_Node], source) -> set[_Node]: ...
@_dispatchable
def ancestors(G: Graph[_Node], source) -> set[_Node]: ...
@_dispatchable
def is_directed_acyclic_graph(G: Graph[_Node]) -> bool: ...
@_dispatchable
def topological_generations(G: DiGraph[_Node]) -> Generator[list[_Node]]: ...
@_dispatchable
def topological_sort(G: DiGraph[_Node]) -> Generator[_Node]: ...
@_dispatchable
def lexicographical_topological_sort(G: DiGraph[_Node], key: Callable[..., Incomplete] | None = None) -> Generator[_Node]: ...
@_dispatchable
def all_topological_sorts(G: DiGraph[_Node]) -> Generator[list[_Node]]: ...
@_dispatchable
def is_aperiodic(G: DiGraph[_Node]) -> bool: ...
@_dispatchable
def transitive_closure(G: Graph[_Node], reflexive=False) -> Graph[_Node]: ...
@_dispatchable
def transitive_closure_dag(G: DiGraph[_Node], topo_order: Iterable[Incomplete] | None = None) -> DiGraph[_Node]: ...
@_dispatchable
def transitive_reduction(G: DiGraph[_Node]) -> DiGraph[_Node]: ...
@_dispatchable
def antichains(G: DiGraph[_Node], topo_order: Iterable[Incomplete] | None = None) -> Generator[list[_Node]]: ...
@_dispatchable
def dag_longest_path(
    G: DiGraph[_Node],
    weight: str | None = "weight",
    default_weight: int | None = 1,
    topo_order: Iterable[Incomplete] | None = None,
) -> list[_Node]: ...
@_dispatchable
def dag_longest_path_length(G: DiGraph[_Node], weight: str | None = "weight", default_weight: int | None = 1) -> int: ...
@_dispatchable
def dag_to_branching(G: DiGraph[_Node]) -> DiGraph[_Node]: ...
