from networkx.classes.digraph import DiGraph
from networkx.classes.graph import _Node
from networkx.utils.backends import _dispatchable

__all__ = ["stochastic_graph"]

@_dispatchable
def stochastic_graph(G: DiGraph[_Node], copy: bool = True, weight: str = "weight"): ...
