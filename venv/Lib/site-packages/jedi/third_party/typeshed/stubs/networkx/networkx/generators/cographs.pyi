from networkx.utils.backends import _dispatchable

__all__ = ["random_cograph"]

@_dispatchable
def random_cograph(n, seed=None): ...
