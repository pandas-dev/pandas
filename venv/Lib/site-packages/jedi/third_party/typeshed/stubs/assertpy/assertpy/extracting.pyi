from _typeshed import SupportsRichComparison
from collections.abc import Callable, Iterable, Mapping
from typing import Any
from typing_extensions import Self

__tracebackhide__: bool

class ExtractingMixin:
    def extracting(
        self,
        *names: str,
        # The callable must accept the type of the items in the self.val collection.
        filter: str | Mapping[str, Any] | Callable[[Any], bool] = ...,
        sort: str | Iterable[str] | Callable[[Any], SupportsRichComparison] = ...,
    ) -> Self: ...
