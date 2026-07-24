from _typeshed import Incomplete
from collections.abc import Callable

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

from .preflowpush import preflow_push

__all__ = ["maximum_flow", "maximum_flow_value", "minimum_cut", "minimum_cut_value"]
default_flow_func = preflow_push

@_dispatchable
def maximum_flow(
    flowG: Graph[_Node],
    _s: _Node,
    _t: _Node,
    capacity: str = "capacity",
    flow_func: Callable[..., Incomplete] | None = None,
    **kwargs,
): ...
@_dispatchable
def maximum_flow_value(
    flowG: Graph[_Node],
    _s: _Node,
    _t: _Node,
    capacity: str = "capacity",
    flow_func: Callable[..., Incomplete] | None = None,
    **kwargs,
): ...
@_dispatchable
def minimum_cut(
    flowG: Graph[_Node],
    _s: _Node,
    _t: _Node,
    capacity: str = "capacity",
    flow_func: Callable[..., Incomplete] | None = None,
    **kwargs,
): ...
@_dispatchable
def minimum_cut_value(
    flowG: Graph[_Node],
    _s: _Node,
    _t: _Node,
    capacity: str = "capacity",
    flow_func: Callable[..., Incomplete] | None = None,
    **kwargs,
): ...
