from typing import Any
from typing_extensions import Self

__tracebackhide__: bool

class ExceptionMixin:
    def raises(self, ex: type[BaseException] | BaseException) -> Self: ...
    # The types of some_args and some_kwargs must equal the types of the called function.
    def when_called_with(self, *some_args: Any, **some_kwargs: Any) -> Self: ...
