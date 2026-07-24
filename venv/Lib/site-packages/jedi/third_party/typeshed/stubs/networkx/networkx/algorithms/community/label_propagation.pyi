from _collections_abc import dict_values
from _typeshed import Incomplete
from collections.abc import Generator

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable
from numpy.random import RandomState

__all__ = ["label_propagation_communities", "asyn_lpa_communities", "fast_label_propagation_communities"]

@_dispatchable
def fast_label_propagation_communities(G: Graph[_Node], *, weight=None, seed=None) -> Generator[Incomplete]: ...
@_dispatchable
def asyn_lpa_communities(
    G: Graph[_Node], weight: str | None = None, seed: int | RandomState | None = None
) -> Generator[Incomplete, Incomplete]: ...
@_dispatchable
def label_propagation_communities(G: Graph[_Node]) -> dict_values[Incomplete, set[Incomplete]]: ...
