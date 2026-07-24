from _typeshed import Incomplete
from collections.abc import Iterable

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable
from numpy.random import RandomState

__all__ = [
    "configuration_model",
    "havel_hakimi_graph",
    "reverse_havel_hakimi_graph",
    "alternating_havel_hakimi_graph",
    "preferential_attachment_graph",
    "random_graph",
    "gnmk_random_graph",
    "complete_bipartite_graph",
]

@_dispatchable
def complete_bipartite_graph(n1, n2, create_using: Graph[_Node] | None = None): ...
@_dispatchable
def configuration_model(
    aseq: Iterable[Incomplete],
    bseq: Iterable[Incomplete],
    create_using: Graph[_Node] | None = None,
    seed: int | RandomState | None = None,
): ...
@_dispatchable
def havel_hakimi_graph(aseq: Iterable[Incomplete], bseq: Iterable[Incomplete], create_using: Graph[_Node] | None = None): ...
@_dispatchable
def reverse_havel_hakimi_graph(
    aseq: Iterable[Incomplete], bseq: Iterable[Incomplete], create_using: Graph[_Node] | None = None
): ...
@_dispatchable
def alternating_havel_hakimi_graph(
    aseq: Iterable[Incomplete], bseq: Iterable[Incomplete], create_using: Graph[_Node] | None = None
): ...
@_dispatchable
def preferential_attachment_graph(
    aseq: Iterable[Incomplete], p: float, create_using: Graph[_Node] | None = None, seed: int | RandomState | None = None
): ...
@_dispatchable
def random_graph(n: int, m: int, p: float, seed: int | RandomState | None = None, directed: bool | None = False): ...
@_dispatchable
def gnmk_random_graph(n: int, m: int, k: int, seed: int | RandomState | None = None, directed: bool | None = False): ...
