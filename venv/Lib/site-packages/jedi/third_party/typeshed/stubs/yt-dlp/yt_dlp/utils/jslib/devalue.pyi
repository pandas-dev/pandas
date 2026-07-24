from collections.abc import Callable, Generator
from typing import Any

def parse_iter(
    parsed: Any, /, *, revivers: dict[str, Callable[[list[Any]], Any]] | None = None  # parsed: Any from the signature.
) -> Generator[TypeError | IndexError | ValueError, Any, float | None]: ...
def parse(
    parsed: Any, /, *, revivers: dict[str, Callable[[Any], Any]] | None = None  # parsed: Any from the signature.
) -> Any: ...  # returns StopIteration.value if it is raised.
