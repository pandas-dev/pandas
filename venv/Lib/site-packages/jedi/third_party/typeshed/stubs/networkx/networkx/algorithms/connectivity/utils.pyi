from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["build_auxiliary_node_connectivity", "build_auxiliary_edge_connectivity"]

@_dispatchable
def build_auxiliary_node_connectivity(G: Graph[_Node]): ...
@_dispatchable
def build_auxiliary_edge_connectivity(G: Graph[_Node]): ...
