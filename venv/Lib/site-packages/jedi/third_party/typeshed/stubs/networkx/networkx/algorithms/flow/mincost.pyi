from _typeshed import Incomplete, SupportsGetItem

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["min_cost_flow_cost", "min_cost_flow", "cost_of_flow", "max_flow_min_cost"]

@_dispatchable
def min_cost_flow_cost(
    G: Graph[_Node], demand: str = "demand", capacity: str = "capacity", weight: str = "weight"
) -> int | float: ...
@_dispatchable
def min_cost_flow(
    G: Graph[_Node], demand: str = "demand", capacity: str = "capacity", weight: str = "weight"
) -> dict[Incomplete, dict[Incomplete, Incomplete]]: ...
@_dispatchable
def cost_of_flow(G: Graph[_Node], flowDict: SupportsGetItem[Incomplete, Incomplete], weight: str = "weight") -> int | float: ...
@_dispatchable
def max_flow_min_cost(
    G: Graph[_Node], s: str, t: str, capacity: str = "capacity", weight: str = "weight"
) -> dict[Incomplete, dict[Incomplete, Incomplete]]: ...
