from networkx.utils.backends import _dispatchable

__all__ = ["sudoku_graph"]

@_dispatchable
def sudoku_graph(n: int = 3): ...
