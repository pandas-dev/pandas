from abc import ABC
from typing_extensions import TypeIs

__all__ = ["SparseABC", "issparse"]

class SparseABC(ABC): ...

def issparse(x: object) -> TypeIs[SparseABC]: ...
