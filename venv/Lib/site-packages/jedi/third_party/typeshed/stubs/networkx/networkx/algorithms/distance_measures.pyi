from collections.abc import Callable, Mapping
from typing_extensions import TypeAlias

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

_WeightFunction: TypeAlias = Callable[..., int]

__all__ = [
    "eccentricity",
    "diameter",
    "harmonic_diameter",
    "radius",
    "periphery",
    "center",
    "barycenter",
    "resistance_distance",
    "kemeny_constant",
    "effective_graph_resistance",
]

@_dispatchable
def eccentricity(
    G: Graph[_Node],
    v: _Node | None = None,
    sp: Mapping[_Node, Mapping[_Node, int]] | None = None,
    weight: str | _WeightFunction | None = None,
) -> int | dict[_Node, int]: ...  # TODO: overload on v: dict if v is None else int
@_dispatchable
def diameter(
    G: Graph[_Node], e: Mapping[_Node, int] | None = None, usebounds: bool = False, weight: str | _WeightFunction | None = None
) -> int: ...
@_dispatchable
def harmonic_diameter(
    G: Graph[_Node], sp: Mapping[_Node, Mapping[_Node, int]] | None = None, *, weight: str | _WeightFunction | None = None
) -> float: ...
@_dispatchable
def periphery(
    G: Graph[_Node], e: Mapping[_Node, int] | None = None, usebounds: bool = False, weight: str | _WeightFunction | None = None
) -> list[_Node]: ...
@_dispatchable
def radius(
    G: Graph[_Node], e: Mapping[_Node, int] | None = None, usebounds: bool = False, weight: str | _WeightFunction | None = None
) -> int: ...
@_dispatchable
def center(
    G: Graph[_Node], e: Mapping[_Node, int] | None = None, usebounds: bool = False, weight: str | _WeightFunction | None = None
) -> list[_Node]: ...
@_dispatchable
def barycenter(
    G: Graph[_Node],
    weight: str | _WeightFunction | None = None,
    attr: str | None = None,
    sp: Mapping[_Node, Mapping[_Node, int]] | None = None,
) -> list[_Node]: ...
@_dispatchable
def resistance_distance(
    G: Graph[_Node], nodeA: _Node | None = None, nodeB: _Node | None = None, weight: str | None = None, invert_weight: bool = True
) -> float | dict[_Node, float]: ...  # TODO: overload on the nodes: float if both are specified else dict
@_dispatchable
def effective_graph_resistance(G: Graph[_Node], weight: str | None = None, invert_weight: bool = True) -> float: ...
@_dispatchable
def kemeny_constant(G: Graph[_Node], *, weight: str | None = None) -> float: ...
