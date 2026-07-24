from _typeshed import Incomplete
from collections.abc import Callable, Mapping
from typing import Any, Final, Generic, TypeVar, overload
from typing_extensions import ParamSpec, Self

_P = ParamSpec("_P")
_R = TypeVar("_R")
__all__ = ["_dispatchable"]
FAILED_TO_CONVERT: Final[str]

class _dispatchable(Generic[_P, _R]):
    __defaults__: Incomplete
    __kwdefaults__: Incomplete
    __module__: Incomplete
    __qualname__: Incomplete
    __wrapped__: Incomplete
    orig_func: Callable[_P, _R] | None
    name: str
    edge_attrs: dict[str, Any] | None
    node_attrs: dict[str, Any] | None
    preserve_edge_attrs: bool
    preserve_node_attrs: bool
    preserve_graph_attrs: bool
    mutates_input: bool
    optional_graphs: Incomplete
    list_graphs: Incomplete
    graphs: dict[str, int]
    backends: dict[str, Incomplete]
    # Incomplete: Ignoring the case where func=None returns a partial,
    # we only care about `_dispatchable` used as a static-typing decorator
    def __new__(
        cls,
        func: Callable[_P, _R] | None = None,
        *,
        name: str | None = None,
        graphs: str | None | Mapping[str, int] = "G",
        edge_attrs: str | dict[str, Any] | None = None,
        node_attrs: str | dict[str, Any] | None = None,
        preserve_edge_attrs: bool = False,
        preserve_node_attrs: bool = False,
        preserve_graph_attrs: bool = False,
        preserve_all_attrs: bool = False,
        mutates_input: bool = False,
        returns_graph: bool = False,
        implemented_by_nx: bool = True,
    ) -> Self: ...
    @property
    def __doc__(self): ...
    @__doc__.setter
    def __doc__(self, val) -> None: ...
    @property
    def __signature__(self): ...
    # Type system limitations doesn't allow us to define this as it truly should.
    # But specifying backend with backend_kwargs isn't a common usecase anyway
    # and specifying backend as explicitly None is possible but not intended.
    # If this ever changes, update stubs/networkx/@tests/test_cases/check_dispatch_decorator.py
    @overload
    def __call__(self, *args: _P.args, **kwargs: _P.kwargs) -> _R: ...
    @overload
    def __call__(self, *args: Any, backend: str, **backend_kwargs: Any) -> _R: ...
    # @overload
    # def __call__(self, *args: _P.args, backend: None = None, **kwargs: _P.kwargs) -> _R: ...
    # @overload
    # def __call__(self, *args: _P.args, backend: str, **kwargs: _P.kwargs, **backend_kwargs: Any) -> _R: ...
    def __reduce__(self): ...
