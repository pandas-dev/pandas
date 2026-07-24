from _typeshed import Incomplete
from collections.abc import Mapping, Sequence

from networkx.utils.backends import _dispatchable
from numpy.random import RandomState

from ..classes.graph import Graph

__all__ = ["is_valid_joint_degree", "is_valid_directed_joint_degree", "joint_degree_graph", "directed_joint_degree_graph"]

@_dispatchable
def is_valid_joint_degree(joint_degrees: Mapping[int, Mapping[int, int]]) -> bool: ...
@_dispatchable
def joint_degree_graph(
    joint_degrees: Mapping[int, Mapping[int, int]], seed: int | RandomState | None = None
) -> Graph[Incomplete]: ...
@_dispatchable
def is_valid_directed_joint_degree(
    in_degrees: Sequence[int], out_degrees: Sequence[int], nkk: Mapping[int, Mapping[int, int]]
) -> bool: ...
@_dispatchable
def directed_joint_degree_graph(
    in_degrees: Sequence[int],
    out_degrees: Sequence[int],
    nkk: Mapping[int, Mapping[int, int]],
    seed: int | RandomState | None = None,
) -> Graph[Incomplete]: ...
