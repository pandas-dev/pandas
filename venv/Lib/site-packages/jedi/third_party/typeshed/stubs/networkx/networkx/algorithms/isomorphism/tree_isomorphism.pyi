from _typeshed import Incomplete

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["rooted_tree_isomorphism", "tree_isomorphism"]

@_dispatchable
def root_trees(t1, root1, t2, root2): ...
@_dispatchable
def rooted_tree_isomorphism(t1, root1, t2, root2) -> list[tuple[Incomplete, Incomplete]]: ...
@_dispatchable
def tree_isomorphism(t1: Graph[_Node], t2: Graph[_Node]) -> list[tuple[Incomplete, Incomplete]]: ...
