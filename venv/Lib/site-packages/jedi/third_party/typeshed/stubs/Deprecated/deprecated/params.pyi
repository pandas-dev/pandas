from collections.abc import Callable, Iterable
from inspect import Signature
from typing import Any, TypeVar
from typing_extensions import ParamSpec

_P = ParamSpec("_P")
_R = TypeVar("_R")

class DeprecatedParams:
    messages: dict[str, str]
    category: type[Warning]
    def __init__(self, param: str | dict[str, str], reason: str = "", category: type[Warning] = ...) -> None: ...
    def populate_messages(self, param: str | dict[str, str], reason: str = "") -> None: ...
    def check_params(
        self, signature: Signature, *args: Any, **kwargs: Any  # args and kwargs passing to Signature.bind method
    ) -> list[str]: ...
    def warn_messages(self, messages: Iterable[str]) -> None: ...
    def __call__(self, f: Callable[_P, _R]) -> Callable[_P, _R]: ...

deprecated_params = DeprecatedParams
