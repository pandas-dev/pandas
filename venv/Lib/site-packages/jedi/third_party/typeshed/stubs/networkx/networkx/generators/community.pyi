from _typeshed import Incomplete
from collections.abc import Collection

from networkx.utils.backends import _dispatchable

__all__ = [
    "caveman_graph",
    "connected_caveman_graph",
    "relaxed_caveman_graph",
    "random_partition_graph",
    "planted_partition_graph",
    "gaussian_random_partition_graph",
    "ring_of_cliques",
    "windmill_graph",
    "stochastic_block_model",
    "LFR_benchmark_graph",
]

@_dispatchable
def caveman_graph(l, k): ...
@_dispatchable
def connected_caveman_graph(l, k): ...
@_dispatchable
def relaxed_caveman_graph(l, k, p, seed=None): ...
@_dispatchable
def random_partition_graph(sizes, p_in, p_out, seed=None, directed: bool = False): ...
@_dispatchable
def planted_partition_graph(l, k, p_in, p_out, seed=None, directed: bool = False): ...
@_dispatchable
def gaussian_random_partition_graph(n, s, v, p_in, p_out, directed: bool = False, seed=None): ...
@_dispatchable
def ring_of_cliques(num_cliques, clique_size): ...
@_dispatchable
def windmill_graph(n, k): ...
@_dispatchable
def stochastic_block_model(
    sizes,
    p,
    nodelist: Collection[Incomplete] | None = None,
    seed=None,
    directed: bool = False,
    selfloops: bool = False,
    sparse: bool = True,
): ...
@_dispatchable
def LFR_benchmark_graph(
    n,
    tau1,
    tau2,
    mu,
    average_degree=None,
    min_degree=None,
    max_degree=None,
    min_community=None,
    max_community=None,
    tol: float = 1e-07,
    max_iters: int = 500,
    seed=None,
): ...
