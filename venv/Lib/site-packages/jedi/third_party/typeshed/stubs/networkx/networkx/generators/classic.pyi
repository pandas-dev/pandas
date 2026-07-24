from _typeshed import Incomplete

from networkx.classes import Graph
from networkx.utils.backends import _dispatchable

__all__ = [
    "balanced_tree",
    "barbell_graph",
    "binomial_tree",
    "complete_graph",
    "complete_multipartite_graph",
    "circular_ladder_graph",
    "circulant_graph",
    "cycle_graph",
    "dorogovtsev_goltsev_mendes_graph",
    "empty_graph",
    "full_rary_tree",
    "kneser_graph",
    "ladder_graph",
    "lollipop_graph",
    "null_graph",
    "path_graph",
    "star_graph",
    "tadpole_graph",
    "trivial_graph",
    "turan_graph",
    "wheel_graph",
]

@_dispatchable
def full_rary_tree(r, n, create_using=None): ...
@_dispatchable
def kneser_graph(n, k) -> Graph[Incomplete]: ...
@_dispatchable
def balanced_tree(r, h, create_using=None): ...
@_dispatchable
def barbell_graph(m1, m2, create_using=None): ...
@_dispatchable
def binomial_tree(n, create_using=None): ...
@_dispatchable
def complete_graph(n, create_using=None): ...
@_dispatchable
def circular_ladder_graph(n, create_using=None): ...
@_dispatchable
def circulant_graph(n, offsets, create_using=None): ...
@_dispatchable
def cycle_graph(n, create_using=None): ...
@_dispatchable
def dorogovtsev_goltsev_mendes_graph(n, create_using=None): ...
@_dispatchable
def empty_graph(n: Incomplete | int = 0, create_using=None, default=...): ...
@_dispatchable
def ladder_graph(n, create_using=None): ...
@_dispatchable
def lollipop_graph(m, n, create_using=None): ...
@_dispatchable
def null_graph(create_using=None): ...
@_dispatchable
def path_graph(n, create_using=None): ...
@_dispatchable
def star_graph(n, create_using=None): ...
@_dispatchable
def tadpole_graph(m, n, create_using=None) -> Graph[Incomplete] | Incomplete: ...
@_dispatchable
def trivial_graph(create_using=None): ...
@_dispatchable
def turan_graph(n, r): ...
@_dispatchable
def wheel_graph(n, create_using=None): ...
@_dispatchable
def complete_multipartite_graph(*subset_sizes): ...
