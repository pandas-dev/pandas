from typing import Final

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["lukes_partitioning"]

D_EDGE_W: Final = "weight"
D_EDGE_VALUE: Final[float]
D_NODE_W: Final = "weight"
D_NODE_VALUE: Final = 1
PKEY: Final = "partitions"
CLUSTER_EVAL_CACHE_SIZE: Final = 2048

@_dispatchable
def lukes_partitioning(G: Graph[_Node], max_size: int, node_weight=None, edge_weight=None): ...
