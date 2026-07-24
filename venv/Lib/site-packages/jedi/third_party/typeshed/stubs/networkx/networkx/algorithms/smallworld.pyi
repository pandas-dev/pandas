from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable
from numpy.random import RandomState

__all__ = ["random_reference", "lattice_reference", "sigma", "omega"]

@_dispatchable
def random_reference(G: Graph[_Node], niter: int = 1, connectivity: bool = True, seed: int | RandomState | None = None): ...
@_dispatchable
def lattice_reference(
    G: Graph[_Node], niter: int = 5, D=None, connectivity: bool = True, seed: int | RandomState | None = None
): ...
@_dispatchable
def sigma(G: Graph[_Node], niter: int = 100, nrand: int = 10, seed: int | RandomState | None = None) -> float: ...
@_dispatchable
def omega(G: Graph[_Node], niter: int = 5, nrand: int = 10, seed: int | RandomState | None = None) -> float: ...
