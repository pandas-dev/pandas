from collections.abc import Hashable
from functools import cached_property
from typing import Any, ClassVar, overload
from typing_extensions import Self, TypeAlias, TypeVar

from networkx.classes.coreviews import MultiAdjacencyView
from networkx.classes.graph import Graph, _MapFactory, _Node
from networkx.classes.multidigraph import MultiDiGraph
from networkx.classes.reportviews import DiMultiDegreeView, MultiDegreeView, MultiEdgeView, OutMultiEdgeView

_MultiEdge: TypeAlias = tuple[_Node, _Node, int]  # noqa: Y047

_DefaultT = TypeVar("_DefaultT")
_KeyT = TypeVar("_KeyT", bound=Hashable)

__all__ = ["MultiGraph"]

# NOTE: Graph subclasses relationships are so complex
# we're only overriding methods that differ in signature from the base classes
# to use inheritance to our advantage and reduce complexity
class MultiGraph(Graph[_Node]):
    edge_key_dict_factory: ClassVar[_MapFactory]
    def to_directed_class(self) -> type[MultiDiGraph[_Node]]: ...
    def to_undirected_class(self) -> type[MultiGraph[_Node]]: ...
    # @_dispatchable adds `backend` argument, but this decorated is unsupported constructor type here
    # and __init__() ignores this argument
    def __new__(cls, *args, backend=None, **kwargs) -> Self: ...
    def __init__(self, incoming_graph_data=None, multigraph_input: bool | None = None, **attr: Any) -> None: ...
    @cached_property
    def adj(self) -> MultiAdjacencyView[_Node, _Node, dict[str, Any]]: ...  # data can be any type
    def new_edge_key(self, u: _Node, v: _Node) -> int: ...
    # key : hashable identifier, optional (default=lowest unused integer)
    @overload  # type: ignore[override] # More complex overload
    def add_edge(self, u_for_edge: _Node, v_for_edge: _Node, key: int | None = None, **attr: Any) -> int: ...
    @overload
    def add_edge(self, u_for_edge: _Node, v_for_edge: _Node, key: _KeyT, **attr: Any) -> _KeyT: ...
    def remove_edge(self, u: _Node, v: _Node, key: Hashable | None = None) -> None: ...
    def has_edge(self, u: _Node, v: _Node, key: Hashable | None = None) -> bool: ...
    @cached_property
    # Including subtypes' possible return types for LSP
    def edges(self) -> MultiEdgeView[_Node] | OutMultiEdgeView[_Node]: ...
    # key : hashable identifier, optional (default=None).
    # default : any Python object (default=None). Value to return if the specific edge (u, v, key) is not found.
    # Returns: The edge attribute dictionary.
    @overload  # type: ignore[override]
    def get_edge_data(
        self, u: _Node, v: _Node, key: Hashable, default: _DefaultT | None = None
    ) -> dict[str, Any] | _DefaultT: ...
    # default : any Python object (default=None). Value to return if there are no edges between u and v and no key is specified.
    # Returns: A dictionary mapping edge keys to attribute dictionaries for each of those edges if no specific key is provided.
    @overload
    def get_edge_data(
        self, u: _Node, v: _Node, key: None = None, default: _DefaultT | None = None
    ) -> dict[Hashable, dict[str, Any] | _DefaultT]: ...
    def copy(self, as_view: bool = False) -> Self: ...
    @cached_property
    # Including subtypes' possible return types for LSP
    def degree(self) -> MultiDegreeView[_Node] | DiMultiDegreeView[_Node]: ...
    def to_directed(self, as_view: bool = False) -> MultiDiGraph[_Node]: ...
    def to_undirected(self, as_view: bool = False) -> MultiGraph[_Node]: ...
