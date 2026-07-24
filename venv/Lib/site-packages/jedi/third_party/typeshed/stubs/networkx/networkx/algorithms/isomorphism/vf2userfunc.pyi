from _typeshed import Incomplete

from . import isomorphvf2 as vf2

__all__ = ["GraphMatcher", "DiGraphMatcher", "MultiGraphMatcher", "MultiDiGraphMatcher"]

class GraphMatcher(vf2.GraphMatcher):
    node_match: Incomplete
    edge_match: Incomplete
    G1_adj: Incomplete
    G2_adj: Incomplete

    def __init__(self, G1, G2, node_match=None, edge_match=None) -> None: ...
    semantic_feasibility: Incomplete

class DiGraphMatcher(vf2.DiGraphMatcher):
    node_match: Incomplete
    edge_match: Incomplete
    G1_adj: Incomplete
    G2_adj: Incomplete

    def __init__(self, G1, G2, node_match=None, edge_match=None) -> None: ...
    def semantic_feasibility(self, G1_node, G2_node): ...

class MultiGraphMatcher(GraphMatcher): ...
class MultiDiGraphMatcher(DiGraphMatcher): ...
