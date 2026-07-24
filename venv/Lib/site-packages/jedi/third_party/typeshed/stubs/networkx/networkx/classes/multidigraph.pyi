from functools import cached_property
from typing import Any

from networkx.classes.coreviews import MultiAdjacencyView
from networkx.classes.digraph import DiGraph
from networkx.classes.graph import _Node
from networkx.classes.multigraph import MultiGraph
from networkx.classes.reportviews import (
    DiMultiDegreeView,
    InMultiDegreeView,
    InMultiEdgeView,
    OutMultiDegreeView,
    OutMultiEdgeView,
)

__all__ = ["MultiDiGraph"]

# NOTE: Graph subclasses relationships are so complex
# we're only overriding methods that differ in signature from the base classes
# to use inheritance to our advantage and reduce complexity
class MultiDiGraph(MultiGraph[_Node], DiGraph[_Node]):
    @cached_property
    def succ(self) -> MultiAdjacencyView[_Node, _Node, dict[str, Any]]: ...
    @cached_property
    def pred(self) -> MultiAdjacencyView[_Node, _Node, dict[str, Any]]: ...
    @cached_property
    def edges(self) -> OutMultiEdgeView[_Node]: ...
    @cached_property
    def out_edges(self) -> OutMultiEdgeView[_Node]: ...
    @cached_property
    def in_edges(self) -> InMultiEdgeView[_Node]: ...
    @cached_property
    def degree(self) -> DiMultiDegreeView[_Node]: ...
    @cached_property
    def in_degree(self) -> InMultiDegreeView[_Node]: ...
    @cached_property
    def out_degree(self) -> OutMultiDegreeView[_Node]: ...
    def to_undirected(self, reciprocal: bool = False, as_view: bool = False) -> MultiGraph[_Node]: ...  # type: ignore[override] # Has an additional `reciprocal` keyword argument
