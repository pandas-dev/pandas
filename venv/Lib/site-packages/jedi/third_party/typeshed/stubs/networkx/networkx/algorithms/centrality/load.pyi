from _typeshed import Incomplete

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["load_centrality", "edge_load_centrality"]

@_dispatchable
def newman_betweenness_centrality(
    G: Graph[_Node], v=None, cutoff: bool | None = None, normalized: bool | None = True, weight: str | None = None
) -> float | dict[Incomplete, float]: ...

load_centrality = newman_betweenness_centrality

@_dispatchable
def edge_load_centrality(G: Graph[_Node], cutoff: bool | None = False): ...
