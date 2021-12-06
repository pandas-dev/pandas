from contextlib import contextmanager
from typing import Optional, Tuple, Iterator

# These are global mutable state. Don't add anything here unless there's a very
# good reason.

# Value varies by file being processed
strict_optional = False
find_occurrences = None  # type: Optional[Tuple[str, str]]


@contextmanager
def strict_optional_set(value: bool) -> Iterator[None]:
    global strict_optional
    saved = strict_optional
    strict_optional = value
    yield
    strict_optional = saved
