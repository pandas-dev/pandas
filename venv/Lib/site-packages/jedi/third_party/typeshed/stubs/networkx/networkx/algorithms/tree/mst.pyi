from _typeshed import Incomplete
from collections.abc import Callable, Generator
from dataclasses import dataclass
from enum import Enum
from typing import Final, Literal
from typing_extensions import Self

from networkx.classes.digraph import DiGraph
from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable
from numpy.random import RandomState

__all__ = [
    "minimum_spanning_edges",
    "maximum_spanning_edges",
    "minimum_spanning_tree",
    "maximum_spanning_tree",
    "number_of_spanning_trees",
    "random_spanning_tree",
    "partition_spanning_tree",
    "EdgePartition",
    "SpanningTreeIterator",
]

class EdgePartition(Enum):
    OPEN = 0
    INCLUDED = 1
    EXCLUDED = 2

@_dispatchable
def boruvka_mst_edges(G: Graph[_Node], minimum=True, weight="weight", keys=False, data=True, ignore_nan=False): ...
@_dispatchable
def kruskal_mst_edges(G: Graph[_Node], minimum, weight="weight", keys=True, data=True, ignore_nan=False, partition=None): ...
@_dispatchable
def prim_mst_edges(G: Graph[_Node], minimum, weight="weight", keys=True, data=True, ignore_nan=False): ...

ALGORITHMS: Final[dict[str, Callable[..., Generator[Incomplete, Incomplete, Incomplete]]]]

@_dispatchable
def minimum_spanning_edges(
    G: Graph[_Node],
    algorithm: str = "kruskal",
    weight: str = "weight",
    keys: bool = True,
    data: bool | None = True,
    ignore_nan: bool = False,
): ...
@_dispatchable
def maximum_spanning_edges(
    G: Graph[_Node],
    algorithm: str = "kruskal",
    weight: str = "weight",
    keys: bool = True,
    data: bool | None = True,
    ignore_nan: bool = False,
): ...
@_dispatchable
def minimum_spanning_tree(G: Graph[_Node], weight: str = "weight", algorithm: str = "kruskal", ignore_nan: bool = False): ...
@_dispatchable
def partition_spanning_tree(
    G: Graph[_Node], minimum: bool = True, weight: str = "weight", partition: str = "partition", ignore_nan: bool = False
): ...
@_dispatchable
def maximum_spanning_tree(G: Graph[_Node], weight: str = "weight", algorithm: str = "kruskal", ignore_nan: bool = False): ...
@_dispatchable
def random_spanning_tree(
    G: Graph[_Node], weight: str | None = None, *, multiplicative=True, seed: int | RandomState | None = None
): ...

class SpanningTreeIterator:
    @dataclass(order=True)
    class Partition:
        mst_weight: float
        partition_dict: dict[Incomplete, Incomplete]
        def __copy__(self) -> SpanningTreeIterator.Partition: ...

    G: Incomplete
    weight: Incomplete
    minimum: Incomplete
    ignore_nan: Incomplete
    partition_key: str

    def __init__(self, G: DiGraph[_Node], weight: str = "weight", minimum: bool = True, ignore_nan: bool = False) -> None: ...
    partition_queue: Incomplete

    def __iter__(self) -> Self: ...
    def __next__(self): ...

@_dispatchable
def number_of_spanning_trees(G: Graph[_Node], *, root=None, weight=None) -> float | Literal[0]: ...
