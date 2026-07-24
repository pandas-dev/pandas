from networkx.classes.digraph import DiGraph
from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["is_arborescence", "is_branching", "is_forest", "is_tree"]

@_dispatchable
def is_arborescence(G: Graph[_Node]) -> bool: ...
@_dispatchable
def is_branching(G: DiGraph[_Node]) -> bool: ...
@_dispatchable
def is_forest(G: Graph[_Node]) -> bool: ...
@_dispatchable
def is_tree(G: Graph[_Node]) -> bool: ...
