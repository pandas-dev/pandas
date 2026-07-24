from collections.abc import Collection

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["laplacian_centrality"]

@_dispatchable
def laplacian_centrality(
    G: Graph[_Node],
    normalized: bool = True,
    nodelist: Collection[_Node] | None = None,
    weight: str | None = "weight",
    walk_type: str | None = None,
    alpha: float = 0.95,
): ...
