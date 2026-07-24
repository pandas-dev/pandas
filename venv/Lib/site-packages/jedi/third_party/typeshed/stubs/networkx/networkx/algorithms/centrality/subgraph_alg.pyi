from _typeshed import Incomplete

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["subgraph_centrality_exp", "subgraph_centrality", "communicability_betweenness_centrality", "estrada_index"]

@_dispatchable
def subgraph_centrality_exp(G: Graph[_Node], *, normalized: bool = False) -> dict[Incomplete, float]: ...
@_dispatchable
def subgraph_centrality(G: Graph[_Node], *, normalized: bool = False) -> dict[Incomplete, float]: ...
@_dispatchable
def communicability_betweenness_centrality(G: Graph[_Node]) -> dict[Incomplete, Incomplete]: ...
@_dispatchable
def estrada_index(G: Graph[_Node]) -> float: ...
