from _typeshed import Incomplete

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["current_flow_closeness_centrality", "information_centrality"]

@_dispatchable
def current_flow_closeness_centrality(
    G: Graph[_Node], weight: str | None = None, dtype: type = ..., solver: str = "lu"
) -> dict[Incomplete, float]: ...

information_centrality = current_flow_closeness_centrality
