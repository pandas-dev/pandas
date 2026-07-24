from _typeshed import Incomplete
from collections.abc import Generator
from typing import NamedTuple

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["vf2pp_isomorphism", "vf2pp_is_isomorphic", "vf2pp_all_isomorphisms"]

class _GraphParameters(NamedTuple):
    G1: Incomplete
    G2: Incomplete
    G1_labels: Incomplete
    G2_labels: Incomplete
    nodes_of_G1Labels: Incomplete
    nodes_of_G2Labels: Incomplete
    G2_nodes_of_degree: Incomplete

class _StateParameters(NamedTuple):
    mapping: Incomplete
    reverse_mapping: Incomplete
    T1: Incomplete
    T1_in: Incomplete
    T1_tilde: Incomplete
    T1_tilde_in: Incomplete
    T2: Incomplete
    T2_in: Incomplete
    T2_tilde: Incomplete
    T2_tilde_in: Incomplete

@_dispatchable
def vf2pp_isomorphism(G1: Graph[_Node], G2: Graph[_Node], node_label: str | None = None, default_label: float | None = None): ...
@_dispatchable
def vf2pp_is_isomorphic(
    G1: Graph[_Node], G2: Graph[_Node], node_label: str | None = None, default_label: float | None = None
): ...
@_dispatchable
def vf2pp_all_isomorphisms(
    G1: Graph[_Node], G2: Graph[_Node], node_label: str | None = None, default_label: float | None = None
) -> Generator[Incomplete, None, Incomplete]: ...
