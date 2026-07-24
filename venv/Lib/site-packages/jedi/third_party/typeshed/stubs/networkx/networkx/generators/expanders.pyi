from typing_extensions import deprecated

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = [
    "margulis_gabber_galil_graph",
    "chordal_cycle_graph",
    "paley_graph",
    "maybe_regular_expander",
    "maybe_regular_expander_graph",
    "is_regular_expander",
    "random_regular_expander_graph",
]

@_dispatchable
def margulis_gabber_galil_graph(n, create_using=None): ...
@_dispatchable
def chordal_cycle_graph(p, create_using=None): ...
@_dispatchable
def paley_graph(p, create_using=None): ...
@_dispatchable
def maybe_regular_expander_graph(n, d, *, create_using=None, max_tries: int = 100, seed=None): ...
@deprecated(
    "`maybe_regular_expander` is a deprecated alias for `maybe_regular_expander_graph`. "
    "Use `maybe_regular_expander_graph` instead."
)
def maybe_regular_expander(n, d, *, create_using=None, max_tries: int = 100, seed=None): ...
@_dispatchable
def is_regular_expander(G: Graph[_Node], *, epsilon: float = 0) -> bool: ...
@_dispatchable
def random_regular_expander_graph(n, d, *, epsilon=0, create_using=None, max_tries=100, seed=None): ...
