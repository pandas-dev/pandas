from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = [
    "resource_allocation_index",
    "jaccard_coefficient",
    "adamic_adar_index",
    "preferential_attachment",
    "cn_soundarajan_hopcroft",
    "ra_index_soundarajan_hopcroft",
    "within_inter_cluster",
    "common_neighbor_centrality",
]

@_dispatchable
def resource_allocation_index(G: Graph[_Node], ebunch=None): ...
@_dispatchable
def jaccard_coefficient(G: Graph[_Node], ebunch=None): ...
@_dispatchable
def adamic_adar_index(G: Graph[_Node], ebunch=None): ...
@_dispatchable
def common_neighbor_centrality(G: Graph[_Node], ebunch=None, alpha=0.8): ...
@_dispatchable
def preferential_attachment(G: Graph[_Node], ebunch=None): ...
@_dispatchable
def cn_soundarajan_hopcroft(G: Graph[_Node], ebunch=None, community: str | None = "community"): ...
@_dispatchable
def ra_index_soundarajan_hopcroft(G: Graph[_Node], ebunch=None, community: str | None = "community"): ...
@_dispatchable
def within_inter_cluster(G: Graph[_Node], ebunch=None, delta: float | None = 0.001, community: str | None = "community"): ...
