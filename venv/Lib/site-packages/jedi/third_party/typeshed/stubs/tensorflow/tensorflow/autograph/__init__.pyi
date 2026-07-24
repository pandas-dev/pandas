from collections.abc import Callable
from typing import Any, TypeVar

from tensorflow.autograph.experimental import Feature

_Type = TypeVar("_Type")

def set_verbosity(level: int, alsologtostdout: bool = False) -> None: ...
def to_code(
    entity: Callable[..., Any],
    recursive: bool = True,
    experimental_optional_features: None | Feature | tuple[Feature, ...] = None,
) -> str: ...
def to_graph(
    entity: _Type, recursive: bool = True, experimental_optional_features: None | Feature | tuple[Feature, ...] = None
) -> _Type: ...
def trace(*args: Any) -> None: ...
