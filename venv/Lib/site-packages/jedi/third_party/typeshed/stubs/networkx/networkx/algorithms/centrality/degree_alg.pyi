from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["degree_centrality", "in_degree_centrality", "out_degree_centrality"]

@_dispatchable
def degree_centrality(G: Graph[_Node]) -> dict[_Node, float]: ...
@_dispatchable
def in_degree_centrality(G: Graph[_Node]) -> dict[_Node, float]: ...
@_dispatchable
def out_degree_centrality(G: Graph[_Node]) -> dict[_Node, float]: ...
