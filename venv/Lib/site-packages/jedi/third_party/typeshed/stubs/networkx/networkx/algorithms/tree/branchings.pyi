from _typeshed import Incomplete
from dataclasses import dataclass
from typing import Final
from typing_extensions import Self

from networkx.classes.digraph import DiGraph
from networkx.classes.graph import _Node
from networkx.utils.backends import _dispatchable
from numpy.random import RandomState

__all__ = [
    "branching_weight",
    "greedy_branching",
    "maximum_branching",
    "minimum_branching",
    "minimal_branching",
    "maximum_spanning_arborescence",
    "minimum_spanning_arborescence",
    "ArborescenceIterator",
]

KINDS: Final[set[str]]
STYLES: Final[dict[str, str]]
INF: Final[float]

def random_string(L=15, seed=None): ...
@_dispatchable
def branching_weight(G: DiGraph[_Node], attr: str = "weight", default: float = 1): ...
@_dispatchable
def greedy_branching(
    G: DiGraph[_Node], attr: str = "weight", default: float = 1, kind: str = "max", seed: int | RandomState | None = None
): ...
@_dispatchable
def maximum_branching(
    G: DiGraph[_Node], attr: str = "weight", default: float = 1, preserve_attrs: bool = False, partition: str | None = None
): ...
@_dispatchable
def minimum_branching(
    G: DiGraph[_Node], attr: str = "weight", default: float = 1, preserve_attrs: bool = False, partition: str | None = None
): ...
@_dispatchable
def minimal_branching(G: DiGraph[_Node], /, *, attr="weight", default=1, preserve_attrs=False, partition=None): ...
@_dispatchable
def maximum_spanning_arborescence(
    G: DiGraph[_Node], attr: str = "weight", default: float = 1, preserve_attrs: bool = False, partition: str | None = None
): ...
@_dispatchable
def minimum_spanning_arborescence(
    G: DiGraph[_Node], attr: str = "weight", default: float = 1, preserve_attrs: bool = False, partition: str | None = None
): ...

class ArborescenceIterator:
    @dataclass(order=True)
    class Partition:
        mst_weight: float
        partition_dict: dict[Incomplete, Incomplete]
        def __copy__(self) -> ArborescenceIterator.Partition: ...

    G: Incomplete
    weight: Incomplete
    minimum: Incomplete
    method: Incomplete
    partition_key: str
    init_partition: Incomplete

    def __init__(self, G: DiGraph[_Node], weight: str = "weight", minimum: bool = True, init_partition=None) -> None: ...
    partition_queue: Incomplete

    def __iter__(self) -> Self: ...
    def __next__(self): ...
