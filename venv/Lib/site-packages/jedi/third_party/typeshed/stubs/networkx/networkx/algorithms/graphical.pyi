from _typeshed import Incomplete
from collections.abc import Iterable
from typing import Literal

from networkx.utils.backends import _dispatchable

__all__ = [
    "is_graphical",
    "is_multigraphical",
    "is_pseudographical",
    "is_digraphical",
    "is_valid_degree_sequence_erdos_gallai",
    "is_valid_degree_sequence_havel_hakimi",
]

@_dispatchable
def is_graphical(sequence: Iterable[Incomplete], method: Literal["eg", "hh"] = "eg") -> bool: ...
@_dispatchable
def is_valid_degree_sequence_havel_hakimi(deg_sequence: Iterable[Incomplete]) -> bool: ...
@_dispatchable
def is_valid_degree_sequence_erdos_gallai(deg_sequence: Iterable[Incomplete]) -> bool: ...
@_dispatchable
def is_multigraphical(sequence: Iterable[Incomplete]) -> bool: ...
@_dispatchable
def is_pseudographical(sequence: Iterable[Incomplete]) -> bool: ...
@_dispatchable
def is_digraphical(in_sequence: Iterable[Incomplete], out_sequence: Iterable[Incomplete]) -> bool: ...
