from _typeshed import Incomplete

from networkx.classes.digraph import DiGraph
from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable
from numpy.random import RandomState

__all__ = [
    "hamiltonian_path",
    "is_reachable",
    "is_strongly_connected",
    "is_tournament",
    "random_tournament",
    "score_sequence",
    "tournament_matrix",
]

@_dispatchable
def is_tournament(G: Graph[_Node]) -> bool: ...
@_dispatchable
def hamiltonian_path(G: Graph[_Node]) -> list[Incomplete]: ...
@_dispatchable
def random_tournament(n: int, seed: int | RandomState | None = None) -> DiGraph[Incomplete]: ...
@_dispatchable
def tournament_matrix(G: Graph[_Node]): ...
@_dispatchable
def score_sequence(G: Graph[_Node]) -> list[Incomplete]: ...
@_dispatchable
def is_reachable(G: Graph[_Node], s: _Node, t: _Node) -> bool: ...
@_dispatchable
def is_strongly_connected(G: Graph[_Node]) -> bool: ...
