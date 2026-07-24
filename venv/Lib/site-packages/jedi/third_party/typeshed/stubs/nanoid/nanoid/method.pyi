from collections.abc import Callable, Sequence
from typing_extensions import TypeAlias

_Algorithm: TypeAlias = Callable[[int], Sequence[int]]

def method(algorithm: _Algorithm, alphabet: str, size: int) -> str: ...
