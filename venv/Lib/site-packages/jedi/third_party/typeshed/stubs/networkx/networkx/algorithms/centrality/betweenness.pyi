from networkx.classes.graph import Graph, _Edge, _Node
from networkx.utils.backends import _dispatchable
from numpy.random import RandomState

__all__ = ["betweenness_centrality", "edge_betweenness_centrality"]

@_dispatchable
def betweenness_centrality(
    G: Graph[_Node],
    k: int | None = None,
    normalized: bool | None = True,
    weight: str | None = None,
    endpoints: bool | None = False,
    seed: int | RandomState | None = None,
) -> dict[_Node, float]: ...
@_dispatchable
def edge_betweenness_centrality(
    G: Graph[_Node],
    k: int | None = None,
    normalized: bool | None = True,
    weight: str | None = None,
    seed: int | RandomState | None = None,
) -> dict[_Edge[_Node], float]: ...
