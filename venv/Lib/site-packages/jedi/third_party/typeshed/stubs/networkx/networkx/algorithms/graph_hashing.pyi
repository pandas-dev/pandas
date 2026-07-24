from _typeshed import Incomplete

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["weisfeiler_lehman_graph_hash", "weisfeiler_lehman_subgraph_hashes"]

@_dispatchable
def weisfeiler_lehman_graph_hash(
    G: Graph[_Node],
    edge_attr: str | None = None,
    node_attr: str | None = None,
    iterations: int | None = 3,
    digest_size: int | None = 16,
) -> str: ...
@_dispatchable
def weisfeiler_lehman_subgraph_hashes(
    G: Graph[_Node],
    edge_attr: str | None = None,
    node_attr: str | None = None,
    iterations: int | None = 3,
    digest_size: int | None = 16,
    include_initial_labels: bool | None = False,
) -> dict[Incomplete, list[str]]: ...
