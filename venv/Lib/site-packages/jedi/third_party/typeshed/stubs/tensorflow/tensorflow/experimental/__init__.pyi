from _typeshed import Incomplete
from abc import ABC
from typing import Generic, TypeVar

_T_co = TypeVar("_T_co", covariant=True)

class Optional(ABC, Generic[_T_co]):
    def __getattr__(self, name: str) -> Incomplete: ...

def __getattr__(name: str): ...  # incomplete module
