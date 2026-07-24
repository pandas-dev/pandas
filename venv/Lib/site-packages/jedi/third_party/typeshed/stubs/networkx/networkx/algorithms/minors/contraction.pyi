from _typeshed import Incomplete
from collections.abc import Callable, Iterable

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["contracted_edge", "contracted_nodes", "equivalence_classes", "identified_nodes", "quotient_graph"]

@_dispatchable
def equivalence_classes(iterable: Iterable[_Node], relation: Callable[[_Node, _Node], bool]) -> set[frozenset[_Node]]: ...
@_dispatchable
def quotient_graph(
    G: Graph[_Node],
    partition,
    edge_relation=None,
    node_data: Callable[..., Incomplete] | None = None,
    edge_data: Callable[..., Incomplete] | None = None,
    weight: str | None = "weight",
    relabel: bool = False,
    create_using: Graph[_Node] | None = None,
): ...
@_dispatchable
def contracted_nodes(
    G: Graph[_Node], u, v, self_loops: bool = True, copy: bool = True, *, store_contraction_as: str | None = "contraction"
): ...

identified_nodes = contracted_nodes

@_dispatchable
def contracted_edge(
    G: Graph[_Node],
    edge: tuple[Incomplete],
    self_loops: bool = True,
    copy: bool = True,
    *,
    store_contraction_as: str | None = "contraction",
): ...
