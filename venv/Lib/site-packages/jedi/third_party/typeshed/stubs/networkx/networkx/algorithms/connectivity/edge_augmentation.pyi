from _typeshed import Incomplete, SupportsGetItem
from collections.abc import Generator
from typing import NamedTuple

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["k_edge_augmentation", "is_k_edge_connected", "is_locally_k_edge_connected"]

@_dispatchable
def is_k_edge_connected(G: Graph[_Node], k: int) -> bool: ...
@_dispatchable
def is_locally_k_edge_connected(G: Graph[_Node], s: _Node, t: _Node, k: int) -> bool: ...
@_dispatchable
def k_edge_augmentation(
    G: Graph[_Node],
    k: int,
    avail: set[tuple[int, int]] | set[tuple[int, int, float]] | SupportsGetItem[tuple[int, int], float] | None = None,
    weight: str | None = None,
    partial: bool = False,
) -> Generator[tuple[_Node, _Node]]: ...
@_dispatchable
def partial_k_edge_augmentation(G: Graph[_Node], k, avail, weight: str | None = None): ...
@_dispatchable
def one_edge_augmentation(G: Graph[_Node], avail=None, weight: str | None = None, partial: bool = False): ...
@_dispatchable
def bridge_augmentation(G: Graph[_Node], avail=None, weight: str | None = None): ...

class MetaEdge(NamedTuple):
    meta_uv: Incomplete
    uv: Incomplete
    w: Incomplete

@_dispatchable
def unconstrained_one_edge_augmentation(G: Graph[_Node]): ...
@_dispatchable
def weighted_one_edge_augmentation(G: Graph[_Node], avail, weight: str | None = None, partial: bool = False): ...
@_dispatchable
def unconstrained_bridge_augmentation(G: Graph[_Node]): ...
@_dispatchable
def weighted_bridge_augmentation(G: Graph[_Node], avail, weight: str | None = None): ...
@_dispatchable
def collapse(G: Graph[_Node], grouped_nodes): ...
@_dispatchable
def complement_edges(G: Graph[_Node]): ...
@_dispatchable
def greedy_k_edge_augmentation(G: Graph[_Node], k, avail=None, weight: str | None = None, seed=None): ...
